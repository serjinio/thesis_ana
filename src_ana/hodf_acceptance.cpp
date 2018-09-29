/**
 *   \file hodf_acceptance.cpp
 *   \brief Test of HODF acceptance
 *
 */

#include "cli.hpp"
#include "init_ds.hpp"
#include "consts.hpp"
#include "csalg.hpp"
#include "drawing.hpp"
#include "scattyield.hpp"
#include "cuts_conf.hpp"
#include "rootscript.hpp"


static s13::misc::CliOptions cli_opts;


using ElaCuts = s13::ana::ElasticScatteringCuts;
using Evt = s13::ana::ScatteringEvent;


class HodfAcceptanceAlg : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  HodfAcceptanceAlg(ElaCuts cuts) :
    cuts_{cuts}, logger_{"HodfAcceptanceAlg"} {
    init_hists();
  }

  virtual void setup() {
    events_num_ = 0;
    for (int i = 0; i < 5; i++) {
      multi_hit_events_num[i] = 0;
    }
    hist_pla_spectra_->Reset();
    hist_pla_spectra_zoom_->Reset();
    hist_mhit_evts_->Reset();
    hist_double_hit_adj_evts_->Reset();
    hist_acceptance_->Reset();
    hist_fdc2_dist_->Reset();
    hist_theta_corr_1d_->Reset();
    hist_theta_corr_2d_->Reset();
  }

  virtual void process_event(Evt& evt) {
    if (cuts_.IsElasticEvent(evt)) {
      events_num_ += 1;
      double max_q = max_hodf_q(evt);
      int pla_idx = hodf_max_q_idx(evt);
      std::vector<int> fired_plastics = hodf_fired_pla_idx(evt);

      // HODF delta-E cut (to suppress tritons)
      // if (max_q < 7) {
      //   return;
      // }

      // FDC2 A/Z cut
      // if ((cuts_.espri_selector == s13::ana::Espri::left)
      //     && !(evt.fdc2_xpos > -300 && evt.fdc2_xpos < 200)) {
      //   return;
      // } else if ((cuts_.espri_selector == s13::ana::Espri::right)
      //            && !(evt.fdc2_xpos > 500 && evt.fdc2_xpos < 1000)) {
      //   return;
      // } else if ( (cuts_.espri_selector == s13::ana::Espri::both)
      //             && !( (evt.fdc2_xpos > -300 && evt.fdc2_xpos < 200) ||
      //                   (evt.fdc2_xpos > 500 && evt.fdc2_xpos < 1000) ) ) {
      //   return;
      // }

      // threshold for bang to avoid acceptance problems with FDC2
      // double thresh_bang = 0.3;
      // double bang_deg = evt.s1dc_bang * 57.2;
      // if (std::abs(bang_deg) > thresh_bang) {
      //   logger_.debug("cutting event due to high bang value: %.2f deg.", bang_deg);
      //   return;
      // }

      if (pla_idx > -1) {
        hodf_events_num_ += 1;
      }

      if (std::abs(evt.fdc2_xpos) < 2500 && std::abs(evt.fdc2_ypos) < 500) {
        fdc2_events_num_ += 1;
      }

      int multiplicity = event_multiplicity(evt);
      hist_mhit_evts_->Fill(multiplicity);

      // logger_.debug("fired_plastics.size(): %d", fired_plastics.size());
      if (fired_plastics.size() == 2) {
        int fired_pla_diff = std::abs(fired_plastics[0] - fired_plastics[1]);
        // logger_.debug("double hit event fired plastics difference: %i", fired_pla_diff);
        hist_double_hit_adj_evts_->Fill(fired_pla_diff);
      }

      hist_theta_corr_2d_->Fill(evt.he_theta, evt.p_theta);
      hist_theta_corr_1d_->Fill(evt.he_theta - evt.he_theta_theor);

      hist_fdc2_dist_->Fill(evt.fdc2_xpos, evt.fdc2_ypos);
      if (pla_idx != -1) {
        hist_acceptance_->Fill(evt.fdc2_xpos, evt.fdc2_ypos);
      }

      if (( (cuts_.espri_selector == s13::ana::Espri::right)
          && (pla_idx > -1 && pla_idx < 5) ) ||
          ( (cuts_.espri_selector == s13::ana::Espri::left)
          && (pla_idx >= 11 && pla_idx <= 15) ) ) {
        hist_pla_spectra_->Fill(max_q);
        hist_pla_spectra_zoom_->Fill(max_q);
      }

      if ((cuts_.espri_selector == s13::ana::Espri::left)
          && (evt.fdc2_xpos > -300 && evt.fdc2_xpos < 200)) {
        hist_pla_spectra_fdc2_->Fill(max_q);
        if (cuts_.IsPhiCut(evt)) {
          hist_pla_spectra_copl_->Fill(max_q);
        }
      } else if ((cuts_.espri_selector == s13::ana::Espri::right)
                 && (evt.fdc2_xpos > 500 && evt.fdc2_xpos < 1000)) {
        hist_pla_spectra_fdc2_->Fill(max_q);
        if (cuts_.IsPhiCut(evt)) {
          hist_pla_spectra_copl_->Fill(max_q);
        }
      } else if ( (cuts_.espri_selector == s13::ana::Espri::both)
                  && ( (evt.fdc2_xpos > -300 && evt.fdc2_xpos < 200) ||
                        (evt.fdc2_xpos > 500 && evt.fdc2_xpos < 1000) ) ) {
        hist_pla_spectra_fdc2_->Fill(max_q);
        if (cuts_.IsPhiCut(evt)) {
          hist_pla_spectra_copl_->Fill(max_q);
        }
      }
    }
  }

  virtual void finalize() {
    logger_.info("Algorithm run end. Total events: %d.", events_num_);
    logger_.info("FDC2 valid events num.: %u", fdc2_events_num_);
    logger_.info("HODF valid events num.: %u", hodf_events_num_);
  }

  virtual void merge(const TTreeAlgorithmBase& other) {
    throw std::logic_error("Not implemented");
  }

  virtual SimpleTTreeAlgorithmBase* clone() const {
    throw std::logic_error("Not implemented");
  }

  virtual void serialize(std::ostream& os) {
    throw std::logic_error("Not implemented");
  }

  virtual void deserialize(std::istream& is) {
    throw std::logic_error("Not implemented");
  }

  TH1* hist_pla_spectra() {
    return hist_pla_spectra_;
  }

  TH1* hist_pla_spectra_copl() {
    return hist_pla_spectra_copl_;
  }

  TH1* hist_pla_spectra_fdc2() {
    return hist_pla_spectra_fdc2_;
  }

  TH1* hist_pla_spectra_zoom() {
    return hist_pla_spectra_zoom_;
  }

  TH1* hist_multi_hit_events() {
    return hist_mhit_evts_;
  }

  TH1* hist_double_hit_adj_evts() {
    return hist_double_hit_adj_evts_;
  }

  TH2* hist_fdc2_dist() {
    return hist_fdc2_dist_;
  }

  TH2* hist_acceptance() {
    return hist_acceptance_;
  }

  TH2* hist_theta_corr_2d() {
    return hist_theta_corr_2d_;
  }

  TH1* hist_theta_corr_1d() {
    return hist_theta_corr_1d_;
  }

private:

  void init_hists() {
    hist_pla_spectra_zoom_ =
      s13::ana::make_th1(100, 0, 2,
                         "HODF Q", "Q [MeV]", "Counts");
    hist_pla_spectra_ =
      s13::ana::make_th1(100, 0, 28,
                         "HODF Q", "Q [MeV]", "Counts");
    hist_pla_spectra_fdc2_ =
      s13::ana::make_th1(100, 0, 28,
                         "HODF Q (FDC2 cut)",
                         "Q [MeV]", "Counts");
    hist_pla_spectra_copl_ =
      s13::ana::make_th1(100, 0, 28,
                         "HODF Q (FDC2 & copl. cuts)",
                         "Q [MeV]", "Counts");
    hist_fdc2_dist_ =
      s13::ana::make_th2(400, -1200, 1200, 100, -500, 500,
                         "FDC2 X&Y", "X [mm]", "Y [mm]");
    hist_acceptance_ =
      s13::ana::make_th2(400, -1200, 1200, 100, -500, 500,
                         "FDC2 X&Y", "X [mm]", "Y [mm]");
    hist_mhit_evts_ =
      s13::ana::make_th1(5, 1, 6, "Events multiplicity map",
                         "Multiplicity", "Count");
    hist_double_hit_adj_evts_ =
      s13::ana::make_th1(23, 1, 24, "Fired plastics distance (for double-hit events)",
                         "Difference [# of plastics]", "Count");
    hist_theta_corr_2d_ =
      s13::ana::make_th2(100, 0, 10, 100, 55, 70,
                         "Kinematical correlation",
                         "He angle [lab. deg.]", "Proton angle [lab. deg.]");
    hist_theta_corr_1d_ =
      s13::ana::make_th1(100, -8, 8,
                         "Kinematical correlation (scattered He)",
                         "He_theor - He_exp [lab. deg.]", "Counts");

  }

  int event_multiplicity(Evt& evt) {
    int mult = 0;
    double thresh = 0.15;
    for (size_t i = 0; i < 24; i++) {
      if (evt.hodf_q[i] > thresh) {
        mult += 1;
      }
    }
    return mult;
  }

  double max_hodf_q(Evt& evt) {
    double q = 0.15;
    int pla_idx = -1;
    for (int i = 0; i < 24; i++) {
      if (q < evt.hodf_q[i]) {
        q = evt.hodf_q[i];
        pla_idx = i;
      }
    }

    if (pla_idx == -1) {
      return 0;
    } else {
      return q;
    }
  }

  double hodf_max_q_idx(Evt& evt) {
    double q = 0.15;
    int pla_idx = -1;
    for (int i = 0; i < 24; i++) {
      if (q < evt.hodf_q[i]) {
        q = evt.hodf_q[i];
        pla_idx = i;
      }
    }
    return pla_idx;
  }

  std::vector<int> hodf_fired_pla_idx(Evt& evt) {
    std::vector<int> res;
    double thresh_q = 0.15;
    for (int i = 0; i < 24; i++) {
      if (evt.hodf_q[i] > thresh_q) {
        res.push_back(i);
      }
    }
    return res;
  }

  ElaCuts cuts_;

  unsigned events_num_;
  unsigned hodf_events_num_ = 0;
  unsigned fdc2_events_num_ = 0;
  unsigned multi_hit_events_num[5];

  s13::misc::MessageLogger logger_;

  TH1* hist_pla_spectra_;
  TH1* hist_pla_spectra_fdc2_;
  TH1* hist_pla_spectra_copl_;
  TH1* hist_pla_spectra_zoom_;
  TH1* hist_mhit_evts_;
  TH1* hist_double_hit_adj_evts_;
  TH2* hist_acceptance_;
  TH2* hist_fdc2_dist_;
  TH2* hist_theta_corr_2d_;
  TH1* hist_theta_corr_1d_;
};


using PtrTreeWalker = std::shared_ptr<s13::ana::ScatteringTreeWalker>;
using PtrHodfAccAlg = std::shared_ptr<HodfAcceptanceAlg>;


PtrHodfAccAlg check_gaps(ElaCuts cuts) {
  auto alg = PtrHodfAccAlg(new HodfAcceptanceAlg(cuts));
  auto tree_walker = PtrTreeWalker(new s13::ana::ScatteringTreeWalker());
  tree_walker->add(alg);

  auto ds = init_dataset_from_cli(cli_opts);
  tree_walker->walk(*ds[0]);

  return alg;
}

void draw_hodf_acc(s13::ana::RootScript& script, PtrHodfAccAlg alg) {
  s13::misc::MessageLogger logger("draw_hodf_acc()");

  script.NewPage(1).cd();
  alg->hist_pla_spectra()->Draw();
  alg->hist_pla_spectra_fdc2()->SetLineColor(kRed);
  alg->hist_pla_spectra_fdc2()->Draw("SAME");
  gPad->SetLogy();

  script.NewPage(1).cd();
  alg->hist_pla_spectra()->GetYaxis()->SetRangeUser(0, 125);
  alg->hist_pla_spectra()->Draw();
  alg->hist_pla_spectra_fdc2()->SetLineColor(kRed);
  alg->hist_pla_spectra_fdc2()->Draw("SAME");
  alg->hist_pla_spectra_copl()->SetLineColor(kGreen);
  alg->hist_pla_spectra_copl()->Draw("SAME");
  gPad->BuildLegend();

  script.NewPage(1).cd();
  alg->hist_pla_spectra()->Draw();
  gPad->SetLogy();

  script.NewPage(1).cd();
  alg->hist_pla_spectra_zoom()->Draw();
  gPad->SetLogy();

  script.NewPage(1).cd();
  alg->hist_multi_hit_events()->Draw();

  script.NewPage(1).cd();
  alg->hist_multi_hit_events()->Draw();
  gPad->SetLogy();

  script.NewPage(1).cd();
  alg->hist_double_hit_adj_evts()->Draw();

  script.NewPage(1).cd();
  alg->hist_fdc2_dist()->Draw();

  script.NewPage(1).cd();
  alg->hist_acceptance()->Draw();

  TH1* hist_acc_projx = alg->hist_acceptance()->ProjectionX();
  TH1* hist_acc_projy = alg->hist_acceptance()->ProjectionY();
  hist_acc_projx->SetTitle(TString::Format("FDC2 X"));
  hist_acc_projy->SetTitle(TString::Format("FDC2 Y"));
  // hist_acc_projx->GetYaxis()->SetRangeUser(0, 31);
  script.NewPage(1).cd();
  hist_acc_projx->Draw();
  script.NewPage(1).cd();
  hist_acc_projy->Draw();

  double hodf_eff = alg->hist_acceptance()->Integral()
    / alg->hist_fdc2_dist()->Integral() * 100;
  logger.info("HODF efficiency: %.2f%%", hodf_eff);

  script.NewPage(1).cd();
  alg->hist_theta_corr_2d()->Draw("colz");

  script.NewPage(1).cd();
  alg->hist_theta_corr_1d()->GetYaxis()->SetRangeUser(0, 166);
  alg->hist_theta_corr_1d()->Draw();
}

void hodf_acceptance() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  auto alg = check_gaps(cuts);

  s13::ana::RootScript script("hodf_acceptance");
  draw_hodf_acc(script, alg);
}

int main(int argc, char** argv) {
  s13::misc::MessageLogger log("main()");
  TThread::Initialize();
  try {
    cli_opts.parse_options(argc, argv);
  } catch (std::exception& ex) {
    log.error("Invalid argument provided: %s", ex.what());
    return -1;
  }

  hodf_acceptance();
}
