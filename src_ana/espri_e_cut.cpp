

#include "TROOT.h"
#include "TPDF.h"
#include "TPaveStats.h"

#include "cli.hpp"
#include "init_ds.hpp"
#include "consts.hpp"
#include "csalg.hpp"
#include "drawing.hpp"
#include "scattyield.hpp"
#include "protonede.hpp"
#include "treewalk.hpp"
#include "cuts_conf.hpp"
#include "rootscript.hpp"


static s13::misc::CliOptions cli_opts;


using ElaCuts = s13::ana::ElasticScatteringCuts;
using Evt = s13::ana::ScatteringEvent;


class EspriCutAlg : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  static constexpr double k_min_theta = 55;
  static constexpr double k_max_theta = 70;

  EspriCutAlg(const ElaCuts& cuts) :
    logger_{"EspriCutAlg"}, cuts_{cuts} {
      init_hists();
    }

  EspriCutAlg(const EspriCutAlg&) = delete;
  EspriCutAlg(EspriCutAlg&&) = delete;
  EspriCutAlg operator=(const EspriCutAlg&) = delete;
  EspriCutAlg operator=(EspriCutAlg&&) = delete;

  virtual void setup() {
    hist_e_de->Reset();
  }

  virtual void process_event(Evt& evt) {
    if (cuts_.IsSignalEvent(evt)) {
      compute_e_de(evt);
      compute_e_cut_stats(evt);
      compute_e_vs_theta(evt);
    }
  }

  virtual void finalize() {
  }

  virtual void merge(const TTreeAlgorithmBase& other) {
    auto other_alg = static_cast<const EspriCutAlg*>(&other);
    hist_e_de->Add(other_alg->hist_e_de);
    for (int nai_id = 0; nai_id < 7; nai_id++) {
      hists_nai_e_de.at(nai_id)->Add(other_alg->hists_nai_e_de.at(nai_id));
    }
  }

  virtual EspriCutAlg* clone() const {
    EspriCutAlg* alg = new EspriCutAlg(cuts_);
    alg->merge(*this);
    return alg;
  }

  ElaCuts cuts() {
    return cuts_;
  }

  TH2* hist_e_de;
  TH2* hist_theta_e;
  std::vector<TH2*> hists_nai_e_de;

  std::vector<double> e_cut_widths = {0.5, 1, 2, 6};
  std::vector<TH1*> e_cut_dists;
  std::vector<TH1*> e_diff_dists;
  std::vector<TH1*> e_cut_total_kin_corrs;
  std::vector<TH1*> e_cut_in_kin_corrs;
  std::vector<TH1*> e_cut_out_kin_corrs;

private:

  TH2* init_nai_hist(int nai_id) {
    auto th2 = &s13::ana::make_th2;
    auto tit = s13::ana::tstrfmt("NaI %i dE-E", nai_id);
    TH2* hist = th2(240, 0, 12, 240, 0, 160,
                    tit, "{#Delta}E [MeV]", "E [MeV]");
    return hist;
  }

  TH1* init_e_cut_dist_hist(double e_cut_width) {
    auto th1 = &s13::ana::make_th1;
    auto tit = s13::ana::tstrfmt("E cut dist. ( (E_{exp}-E_{sim})/E_{sim}= +/-%.1f)", e_cut_width/2);
    auto hist = th1(50, -1.5, 1.5, tit, "(E_{exp} - E_{sim}) / E_{sim}", "Counts");
    return hist;
  }

  TH1* init_e_diff_dist_hist() {
    auto th1 = &s13::ana::make_th1;
    auto tit = s13::ana::tstrfmt("E diff. dist.");
    auto hist = th1(120, -160, 160, tit, "E_{exp} - E_{sim}", "Counts");
    return hist;
  }

  TH1* init_kin_corr_hist(TString title) {
    auto th1 = &s13::ana::make_th1;
    auto hist = th1(100, -8, 8, title, "#theta_{exp} - #theta_{theor}", "Counts");
    return hist;
  }

  void init_hists() {
    // auto th1 = &s13::ana::make_th1;
    auto th2 = &s13::ana::make_th2;

    hist_e_de = th2(240, 0, 10, 240, 0, 160,
                    "dE-E", "#DeltaE [MeV]", "E [MeV]");
    for (int nai_id = 0; nai_id < 7; nai_id++) {
      hists_nai_e_de.push_back(init_nai_hist(nai_id));
    }

    for (auto el : e_cut_widths) {
      e_cut_dists.push_back(init_e_cut_dist_hist(el));
      e_cut_in_kin_corrs.push_back(init_kin_corr_hist("In cut"));
      e_cut_out_kin_corrs.push_back(init_kin_corr_hist("Out cut"));
      e_cut_total_kin_corrs.push_back(init_kin_corr_hist("Total"));
      e_diff_dists.push_back(init_e_diff_dist_hist());
    }

    hist_theta_e = th2(100, 50, 75, 100, 0, 160,
                       "E vs. proton #theta", "#theta_{p}", "E_{p}");
  }

  int get_hit_nai_id(Evt& evt) {
    double* nais_e_cal = cuts_.espri_selector == s13::ana::Espri::right
      ? evt.esr_naie_cal : evt.esl_naie_cal;
    int hit_nai_id = -1;
    for (int nai_id = 0; nai_id < 7; nai_id++) {
      if (nais_e_cal[nai_id] > 5) {
        hit_nai_id = nai_id;
        break;
      }
    }
    return hit_nai_id;
  }

  void compute_e_de(Evt& evt) {
    double* nais_e_cal = cuts_.espri_selector == s13::ana::Espri::right
      ? evt.esr_naie_cal : evt.esl_naie_cal;

    auto proton_e = evt.p_e_eff;
    auto proton_de = evt.p_de;

    if (proton_e < 15) {
      return; // bad E
    }

    int hit_nai_id = get_hit_nai_id(evt);
    if (hit_nai_id < 4 && hit_nai_id > 1) {
      hist_e_de->Fill(proton_de, proton_e);
    }

    for (int nai_id = 0; nai_id < 7; nai_id++) {
      if (nais_e_cal[nai_id] < 5) continue;
      // logger_.debug("NaI %i E: %.1f", nai_id, nais_e_cal[nai_id]);
      hists_nai_e_de[nai_id]->Fill(proton_de, proton_e);
    }
  }

  void compute_e_vs_theta(Evt& evt) {
    if (!(is_well_defined(evt.p_e))) {
      return;
    }
    if (!(is_well_defined(evt.he_theta))) {
      return;
    }

    hist_theta_e->Fill(evt.p_theta_eff, evt.p_e_eff);
  }

  void compute_e_cut_stats(Evt& evt) {
    if (!(is_well_defined(evt.p_e))) {
      return;
    }
    if (!(is_well_defined(evt.p_e_sim))) {
      return;
    }
    if (!(is_well_defined(evt.he_theta))) {
      return;
    }
    if (evt.p_e - evt.p_de < 0.25) {
      return;
    }

    auto e_cuts = cuts_;
    for (size_t i = 0; i < e_cut_widths.size(); i++) {
      double e_diff = evt.p_e_eff - evt.p_e_sim;
      double e_diff_rel = (evt.p_e_eff - evt.p_e_sim) / evt.p_e_sim;
      double theta_diff = evt.he_theta - evt.he_theta_theor;

      e_cut_dists[i]->Fill(e_diff_rel);
      e_diff_dists[i]->Fill(e_diff);
      e_cut_total_kin_corrs[i]->Fill(theta_diff);

      e_cuts.espri_e_width_stdev = e_cut_widths[i];
      if (e_cuts.IsEspriECut(evt)) {
        e_cut_in_kin_corrs[i]->Fill(theta_diff);
      } else {
        e_cut_out_kin_corrs[i]->Fill(theta_diff);
      }
    }
  }

  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;

};


void draw_nais_e_de(s13::ana::RootScript& script, std::shared_ptr<EspriCutAlg> alg) {
  auto de_e_graph = s13::ana::make_he6_proton_de_e_graph();
  script.NewPage(4,2).cd(s13::ana::PadSeq::row);
  for (int nai_id = 0; nai_id < 7; nai_id++) {
    script.cd();
    alg->hists_nai_e_de[nai_id]->Draw("COL");
    de_e_graph->Draw("SAME");
  }
}

void draw_espri_e_cut(s13::ana::RootScript& script, std::shared_ptr<EspriCutAlg> alg) {
  script.NewPage(1).cd();
  alg->hist_e_de->Draw("COL");
  auto de_e_graph = s13::ana::make_he6_proton_de_e_graph();
  de_e_graph->Draw("SAME");
  gPad->SetLogz();
  // s13::ana::draw_text(alg->cuts().GetTotalCutTitle());

  // draw_nais_e_de(script, alg);
}

void draw_e_cut_stats(s13::ana::RootScript& script, std::shared_ptr<EspriCutAlg> alg) {
  for (size_t i = 0; i < alg->e_cut_widths.size(); i++) {
    if (i % 4 == 0) {
      script.NewPage(4, 3).cd(s13::ana::PadSeq::column);
    }

    script.cd();
    alg->e_cut_dists[i]->Draw();
    alg->e_cut_dists[i]->SetMaximum(alg->e_cut_dists[i]->GetMaximum() * 1.5);
    s13::ana::draw_marker_line_ysym(alg->e_cut_widths[i] / 2, 0,
                                    alg->e_cut_dists[i]->GetMaximum());
    s13::ana::draw_text(alg->cuts().GetTotalCutTitle(), 0.15, 0.75, 0.7, 0.9);
    // gPad->SetLogy();

    // script.cd();
    // alg->e_diff_dists[i]->Draw();

    script.cd();
    auto stack = s13::ana::draw_thstack({alg->e_cut_total_kin_corrs[i],
          alg->e_cut_in_kin_corrs[i], alg->e_cut_out_kin_corrs[i]},
      "1-d kin. corrs.", "NOSTACK");
    stack->SetMaximum(stack->GetMaximum() * 0.8);
    gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);
    script.GetCanvas().Update();

    script.cd();
    alg->e_cut_in_kin_corrs[i]->Draw();
  }
  gPad->Update();
}

void draw_e_vs_theta(s13::ana::RootScript& script, std::shared_ptr<EspriCutAlg> alg) {
  script.NewPage(1,1).cd();
  alg->hist_theta_e->Draw("COLZ");
}

void draw_polt_qfs_e(s13::ana::RootScript& script,
                     std::shared_ptr<EspriCutAlg> alg_polt,
                     std::shared_ptr<EspriCutAlg> alg_carbont) {
  script.NewPage(1,1);
  auto polt_e_dist = alg_polt->e_cut_dists[0];
  polt_e_dist->SetTitle("Pol. target");
  auto carbont_e_dist = alg_carbont->e_cut_dists[0];
  carbont_e_dist->Scale(2.65);
  carbont_e_dist->SetTitle("Carbon target (scaled)");
  auto polt_no_carbon = static_cast<TH2*>(polt_e_dist->Clone());
  polt_no_carbon->SetTitle("Pol. - carbon target");
  polt_no_carbon->Add(carbont_e_dist, -1);

  s13::ana::draw_thstack({polt_no_carbon, polt_e_dist, carbont_e_dist},
                         "Pol. target vs. QFS E distributions", "NOSTACK");
  gPad->BuildLegend();
}

void espri_e_cut() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  auto alg =
    s13::ana::walk_alg<EspriCutAlg>(cuts,
                                    init_dataset_from_cli(cli_opts),
                                    cli_opts.use_mt());
  // auto alg_carbon =
  //   s13::ana::walk_alg<EspriCutAlg>(cuts,
  //                                   init_carbon_dataset_from_cli(cli_opts),
  //                                   cli_opts.use_mt());

  s13::ana::RootScript script(cli_opts.output_file_basename());
  gStyle->SetOptStat(1111111);
  // gStyle->SetOptFit();
  draw_espri_e_cut(script, alg);
  draw_e_cut_stats(script, alg);
  draw_e_vs_theta(script, alg);
  // draw_e_vs_theta(script, alg_carbon);
  // draw_polt_qfs_e(script, alg, alg_carbon);
}


int main(int argc, char** argv) {
  s13::misc::MessageLogger log("main()");
  // s13::ana::SetRootStyle(s13::ana::RootStyle::publication);

  TThread::Initialize();
  cli_opts.require_cuts_conf_opt(true);
  try {
    cli_opts.parse_options(argc, argv);
  } catch (std::exception& ex) {
    log.error("Invalid argument provided: %s", ex.what());
    return -1;
  }

  espri_e_cut();
  std::cout << "exiting main\n" << std::flush;
}
