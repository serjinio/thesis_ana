
#include "TChain.h"
#include "TGraph.h"
#include "TH2.h"
#include "TH1.h"

#include <naiarraymodel.hpp>

#include "init_ds.hpp"
#include "macroincludes.hpp"
#include "drawcuts.h"

#include "msgoutput.hpp"


using namespace s13::ana;
using namespace s13::sim;


TChain *g_local_chain = init_4he_partial_dataset(20);

using Interp = ROOT::Math::Interpolator;

bool is_good_event(ScatteringEvent& evt, Espri espri) {
  if (evt.triggers[5] == false) {
    return false;
  }

  double xpos, ypos, e_eff;
  double he_phi, he_theta;
  he_phi = evt.he_phi;
  he_theta = evt.he_theta;
  if (espri == Espri::left) {
    xpos = evt.esl_xpos;
    ypos = evt.esl_ypos;
    e_eff = evt.esl_e_eff;
  } else if (espri == Espri::right) {
    xpos = evt.esr_xpos;
    ypos = evt.esr_ypos;
    e_eff = evt.esr_e_eff;
  } else {
    throw "invalid ESPRI selector!";
  }

  // espri selection cut
  if (espri == Espri::left && evt.he_phi < 0) {
    return false;
  } else if (espri == Espri::right && evt.he_phi > 0) {
    return false;
  }

  // scattered fragment angles cut
  if (std::abs(he_phi) > 7.5) {
    return false;
  }
  if (he_theta < 3 || he_theta > 12) {
    return false;
  }

  // welldef-proton_e cut
  if (isnan(e_eff) || e_eff < 0.45) {
    return false;
  }

  // pde cuts
  if (evt.p_de < 2 || evt.p_de > 6) {
    return false;
  }

  return true;
}

void write_conv_vals(std::vector<double>& energy_means,
                     std::vector<double>& thetas_sorted,
                     Espri espri) {
  assert(espri != Espri::both);
  auto fname = tstrfmt("data/%s_rdc_e_to_theta.csv",
                       espri == Espri::left ? "esl" : "esr");
  std::ofstream ofs(fname);
  for (int i = 0; i < energy_means.size(); i++) {
    ofs << energy_means[i] << "," << thetas_sorted[i] << std::endl;
  }
}

std::shared_ptr<Interp> make_energy_to_rdcxy_interp(TChain& chain, Espri espri) {
  ScatteringYield thetas(7, 55, 70);
  double e_sums[thetas.GetBinsNumber()];
  double e_means[thetas.GetBinsNumber()];
  double e_bin_counts[thetas.GetBinsNumber()];
  int e_indices[thetas.GetBinsNumber()];

  for (int i = 0; i < thetas.GetBinsNumber(); i++) {
    e_sums[i] = 0;
    e_means[i] = 0;
    e_bin_counts[i] = 0;
    e_indices[i] = i;
  }

  ScatteringTreeWalker tree_walker;
  tree_walker.Walk(chain, [&e_bin_counts, &e_sums, &thetas, &espri]
                   (ScatteringEvent& evt) {
      if (!is_good_event(evt, espri)) {
        return;
      }

      double xpos, ypos;
      if (espri == Espri::left) {
        xpos = evt.esl_xpos;
        ypos = evt.esl_ypos;
      } else if (espri == Espri::right) {
        xpos = evt.esr_xpos;
        ypos = evt.esr_ypos;
      } else {
        throw "invalid ESPRI selector!";
      }

      if (evt.p_theta_eff > 55 && evt.p_theta_eff < 70) {
        // if (espri == Espri::left && evt.p_theta_eff > 67.5) {
        //   std::cout << evt.p_theta_eff << ": " << evt.p_e << std::endl;
        // }
        int bin_idx = thetas.GetBinNo(evt.p_theta_eff);
        e_sums[bin_idx] += evt.p_e;
        e_bin_counts[bin_idx] += 1;
      }
    });

  std::cout << "Computing interpolation values for E to Xpos conversion..."
            << std::endl;
  for (int i = 0; i < thetas.GetBinsNumber(); i++) {
    if (e_bin_counts[i] < 1) {
      // throw "Some bins are empty in active range!";
      e_means[i] = 0;
    } else {
      e_means[i] = e_sums[i] / e_bin_counts[i];
    }
  }
  std::sort(e_indices, e_indices + thetas.GetBinsNumber(),
            [&e_means] (int i, int j) {
              return e_means[i] < e_means[j];
            });

  std::vector<double> energy_means;
  std::vector<double> thetas_sorted;
  for (int i = 0; i < thetas.GetBinsNumber(); i++) {
    double e_avg = e_means[e_indices[i]];
    double theta = thetas.xs.at(e_indices[i]);
    energy_means.push_back(e_avg);
    thetas_sorted.push_back(theta);
    std::cout << "\tE:   " << e_avg << " -> theta:   " << theta << std::endl;
  }

  write_conv_vals(energy_means, thetas_sorted, espri);

  return std::shared_ptr<Interp>(new Interp(energy_means, thetas_sorted));
}

class RdcEff {

public:

  struct NaiHitEvt {
    double xpos;
    double ypos;
    double aang;
    double bang;
    double* naie_cal;
  };

  RdcEff(TString name) :
    name_{name}, thresh_e_{0.3}, log_{"RdcEff"} {
    }

  RdcEff(const RdcEff&) = delete;
  RdcEff(const RdcEff&&) = delete;
  RdcEff& operator=(const RdcEff&) = delete;
  RdcEff& operator=(const RdcEff&&) = delete;

  void compute(TChain& chain, s13::ana::Espri espri_sel) {
    espri_ = espri_sel;
    reset();
    init_hists();

    std::string fname = tstrfmt("data/%s_rdc_e_to_theta.csv",
                                espri_sel == Espri::left ? "esl" : "esr").Data();
    e_to_theta_interp_ = LoadRdcEtoThetaConv(fname);

    ScatteringTreeWalker tree_walker;
    tree_walker.Walk(chain, *this);

    hist_rdc_eff_ = static_cast<TH1*>(hist_rdc_yields_->Clone());
    hist_rdc_eff_->Divide(hist_nai_yields_);
    hist_rdc_eff_->SetMinimum(0);
    hist_rdc_eff_->SetMaximum(1);
    hist_rdc_eff_->SetTitle(tstrfmt("%s: RDC efficiency", name_.Data()));

    hist_rdc_eeff_ = static_cast<TH1*>(hist_rdc_eyields_->Clone());
    hist_rdc_eeff_->Divide(hist_nai_eyields_);
    hist_rdc_eeff_->SetMinimum(0);
    hist_rdc_eeff_->SetMaximum(1);
    hist_rdc_eeff_->SetTitle(tstrfmt("%s: RDC efficiency", name_.Data()));

    write_summary();
  }

  void reset() {
    total_evts_ = 0;
    good_e_evts_ = 0;
    rdc_evts_ = 0;
  }

  void write_summary() {
    log_.info("Total events: %.0f", total_evts_);
    log_.info("Good E events: %.0f (%.2f%% relative to total)",
              good_e_evts_, good_e_evts_ / total_evts_ * 100);
    log_.info("\tRDC events: %.2f%% (relative to good E)",
              rdc_evts_ / good_e_evts_ * 100);
  }

  void draw_all(RootScript& script) {
    script.NewPage(1,1).cd();
    hist_nai_yields_->Draw();
    script.NewPage(1,1).cd();
    hist_rdc_yields_->Draw();
    script.NewPage(1,1).cd();
    hist_rdc_eff_->Draw();

    script.NewPage(1,1).cd();
    hist_nai_eyields_->Draw();
    script.NewPage(1,1).cd();
    hist_rdc_eyields_->Draw();
    script.NewPage(1,1).cd();
    hist_rdc_eeff_->Draw();
  }

  void operator() (ScatteringEvent& evt) {
    if (!is_good_event(evt, espri_)) {
      return;
    }

    double xpos, ypos, proton_e;
    double he_phi, he_theta;
    he_phi = evt.he_phi;
    he_theta = evt.he_theta;
    if (espri_ == Espri::left) {
      xpos = evt.esl_xpos;
      ypos = evt.esl_ypos;
      proton_e = evt.esl_e_cal;
    } else if (espri_ == Espri::right) {
      xpos = evt.esr_xpos;
      ypos = evt.esr_ypos;
      proton_e = evt.esr_e_cal;
    } else {
      throw "invalid ESPRI selector!";
    }

    double rdc_theta = e_to_theta_interp_->Eval(evt.p_e);
    log_.info("for E: %.2f found theta: %.2f", evt.p_e, rdc_theta);
    if (isnan(rdc_theta)) {
      return;
    }

    good_e_evts_ += 1;
    hist_nai_yields_->Fill(rdc_theta);
    hist_nai_eyields_->Fill(evt.p_e);

    // if (evt.p_theta > 55 && evt.p_theta < 70) {
    if (abs(xpos) < 500) {
      rdc_evts_ += 1;
      hist_rdc_yields_->Fill(rdc_theta);
      hist_rdc_eyields_->Fill(evt.p_e);
    }
  }

private:

  TVector3 find_extrp_nai_hit_pos(NaiHitEvt& nai_evt) {
    TVector3 position(nai_evt.xpos, nai_evt.ypos, 0);
    TVector3 direction(nai_evt.aang, nai_evt.bang, 1);
    TVector3 hit_position = nais_array_.find_hit_pos(position, direction);
    return hit_position;
  }

  double find_nearest_nai_id(NaiHitEvt& nai_evt) {
    TVector3 position(nai_evt.xpos, nai_evt.ypos, 0);
    TVector3 direction(nai_evt.aang, nai_evt.bang, 1);
    int id = nais_array_.find_nearest_nai_id(position, direction);
    return id;
  }

  void init_hists() {
    hist_rdc_eff_ =
      new TH1F("esl_rdc_eff",
               TString::Format("%s: RDC efficiency vs. RDC X", name_.Data()),
               15, 55, 70);
    hist_rdc_eff_->GetYaxis()->SetTitle("Counts");
    hist_rdc_eff_->GetXaxis()->SetTitle("#theta_p [lab. deg.]");

    hist_rdc_yields_ =
      new TH1F("rdc_yields",
               TString::Format("%s: RDC yields vs. RDC X", name_.Data()),
               15, 55, 70);
    hist_rdc_yields_->GetYaxis()->SetTitle("Counts");
    hist_rdc_yields_->GetXaxis()->SetTitle("#theta_p [lab. deg.]");

    hist_nai_yields_ =
      new TH1F("nai_yields",
               TString::Format("%s: NaIs yields vs. RDC X", name_.Data()),
               15, 55, 70);
    hist_nai_yields_->GetYaxis()->SetTitle("Counts");
    hist_nai_yields_->GetXaxis()->SetTitle("#theta_p [lab. deg.]");

    hist_nai_eyields_ =
      new TH1F("nai_eyields",
               TString::Format("%s: NaIs E yields", name_.Data()),
               15, 0, 200);
    hist_nai_eyields_->GetYaxis()->SetTitle("Counts");
    hist_nai_eyields_->GetXaxis()->SetTitle("E [MeV]");

    hist_rdc_eyields_ =
      new TH1F("rdc_eyields",
               TString::Format("%s: RDC E yields", name_.Data()),
               15, 0, 200);
    hist_rdc_eyields_->GetYaxis()->SetTitle("Counts");
    hist_rdc_eyields_->GetXaxis()->SetTitle("E [MeV]");

    hist_rdc_eeff_ =
      make_th1(15, 0, 200,
               tstrfmt("%s: RDC efficiency", name_.Data()), "E [MeV]");

 }

  // misc stuff
  TString name_;
  ElasticScatteringCuts cuts_;
  s13::ana::Espri espri_;
  NaiArrayBox nais_array_;
  std::shared_ptr<Interp> e_to_theta_interp_;

  // histograms
  TH1* hist_rdc_eff_;
  TH1* hist_rdc_yields_;
  TH1* hist_nai_yields_;

  TH1* hist_nai_eyields_;
  TH1* hist_rdc_eyields_;
  TH1* hist_rdc_eeff_;

  // settings
  double thresh_e_;

  // stats
  double total_evts_;
  double good_e_evts_;
  double rdc_evts_;

  // utilities
  s13::misc::MessageLogger log_;

};

void rdc_eff2() {
  gStyle->SetOptStat(00000000);
  RootScript script("rdc_eff2", g_local_chain);

  RdcEff esl_rdc_eff("ESL");
  esl_rdc_eff.compute(*g_local_chain, s13::ana::Espri::left);
  esl_rdc_eff.draw_all(script);

  RdcEff esr_rdc_eff("ESR");
  esr_rdc_eff.compute(*g_local_chain, s13::ana::Espri::right);
  esr_rdc_eff.draw_all(script);

  // auto esl_e_to_x_interp_ = make_energy_to_rdcxy_interp(*g_local_chain, Espri::left);
  // auto esr_e_to_x_interp_ = make_energy_to_rdcxy_interp(*g_local_chain, Espri::right);
}

int main() {
  rdc_eff2();
}
