

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
#include "vertzalg.hpp"


static s13::misc::CliOptions cli_opts;


using ElaCuts = s13::ana::ElasticScatteringCuts;
using Evt = s13::ana::ScatteringEvent;


class VertzCutAlg : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  static constexpr double k_min_theta = 55;
  static constexpr double k_max_theta = 70;

  VertzCutAlg(const ElaCuts& cuts) :
    logger_{"VertzCutAlg"}, cuts_{cuts} {
      init_hists();
    }

  VertzCutAlg(const VertzCutAlg&) = delete;
  VertzCutAlg(VertzCutAlg&&) = delete;
  VertzCutAlg operator=(const VertzCutAlg&) = delete;
  VertzCutAlg operator=(VertzCutAlg&&) = delete;

  virtual void setup() {
  }

  virtual void process_event(Evt& evt) {
    if (cuts_.IsSignalEvent(evt)) {
      compute_vertz_cut_stats(evt);
    }
    compute_xpos_vs_aang_stats(evt);
  }

  virtual void finalize() {
  }

  ElaCuts cuts() {
    return cuts_;
  }

  ElaCuts cuts_no_aang() {
    auto tmp = cuts_;
    tmp.apply_espri_aang_cut = false;
    return tmp;
  }

  std::vector<double> vertz_widths = {90, 132, 264, 720};
  std::vector<TH1*> vertz_cut_dists;
  std::vector<TH1*> vertz_aang_dists;
  std::vector<TH1*> vertz_in_cut;
  std::vector<TH1*> vertz_out_cut;
  std::vector<TH1*> vertz_total;
  std::vector<TH1*> zpos_total;
  std::vector<TH1*> zpos_in_cut;
  std::vector<TH1*> zpos_out_cut;
  std::vector<TH2*> kin_corrs_total_2d;
  std::vector<TH2*> kin_corrs_in_cut_2d;

  TH2* xpos_vs_aang_total;
  TH2* xpos_vs_aang_in;
  TH2* xpos_vs_aang_kin_corr;
  TH1* xpos_vs_aang_sn_total;
  TH1* xpos_vs_aang_sn_in;
  TH1* xpos_vs_aang_sn_out;
  TH1* xpos_vs_aaang_1d_in;
  TH1* xpos_vs_aaang_1d_total;
  TH1* xpos_vs_aang_e_corr_in;

private:

  TH1* init_vertz_cut_dist_hist(double width, double xlim = 350) {
    auto th1 = &s13::ana::make_th1;
    auto tit = s13::ana::tstrfmt("RDC aang difference (aang = +/- %.1f)", width/2);
    auto hist = th1(100, -xlim, xlim, tit, "RDC (aang - aang_{tgt}) [mrad]", "Counts");
    return hist;
  }

  TH1* init_vertz_dist_hist() {
    auto th1 = &s13::ana::make_th1;
    auto tit = s13::ana::tstrfmt("RDC aang distribution");
    auto hist = th1(100, -2000, 2000, tit, "RDC aang [mrad]", "Counts");
    return hist;
  }

  TH1* init_zpos_dist_hist() {
    auto th1 = &s13::ana::make_th1;
    auto tit = s13::ana::tstrfmt("RDC vertZ distribution");
    auto hist = th1(100, -500, 500, tit, "Vert Z [mm]", "Counts");
    return hist;
  }

  TH1* init_kin_corr_hist(TString title) {
    auto th1 = &s13::ana::make_th1;
    auto hist = th1(100, -8, 8, title, "#theta_{exp} - #theta_{theor}", "Counts");
    return hist;
  }

  TH2* init_2d_kin_corr_hist(TString title) {
    auto th2 = &s13::ana::make_th2;

    auto hist = th2(200, 0, 12, 200, 50, 75, title, "#theta_{He}", "#theta_{p}");
    return hist;
  }

  TH2* init_xpos_vs_aang_hist(TString title) {
    auto th2 = &s13::ana::make_th2;
    auto hist = th2(300, -250, 250, 300, -350, 350, title,
                    "X [mm]", "aang [mrad]");
    return hist;
  }

  TH1* init_e_cut_dist_hist(double e_cut_width_stdev) {
    auto th1 = &s13::ana::make_th1;
    auto tit = s13::ana::tstrfmt("E cut dist. ( (E_{exp}-E_{sim})/E_{sim}="
                                 " +/-%.1f (sigma))", e_cut_width_stdev);
    auto hist = th1(50, -1.5, 1.5, tit, "(E_{exp} - E_{sim}) / E_{sim}", "Counts");
    return hist;
  }

  void init_hists() {
    for (size_t i = 0; i < vertz_widths.size(); i++) {
      if (i < 2) {
        vertz_cut_dists.push_back(init_vertz_cut_dist_hist(vertz_widths[i], 150));
      } else if (i == 2) {
        vertz_cut_dists.push_back(init_vertz_cut_dist_hist(vertz_widths[i], 250));
      } else {
        vertz_cut_dists.push_back(init_vertz_cut_dist_hist(vertz_widths[i], 600));
      }
      vertz_aang_dists.push_back(init_vertz_dist_hist());
      vertz_in_cut.push_back(init_kin_corr_hist("In cut"));
      vertz_out_cut.push_back(init_kin_corr_hist("Out cut"));
      vertz_total.push_back(init_kin_corr_hist("Total"));
      kin_corrs_total_2d.push_back(init_2d_kin_corr_hist("Total"));
      kin_corrs_in_cut_2d.push_back(init_2d_kin_corr_hist("In cut"));
    }

    xpos_vs_aang_total = init_xpos_vs_aang_hist("X vs. aang - total");
    xpos_vs_aang_in = init_xpos_vs_aang_hist("X vs. aang - in cut");
    xpos_vs_aang_kin_corr = init_2d_kin_corr_hist("S/N (X vs. aang)");
    xpos_vs_aaang_1d_in = init_vertz_cut_dist_hist(cuts_.espri_aang_width);
    xpos_vs_aaang_1d_in->SetTitle(s13::ana::tstrfmt("aang diff. in cut (+/-%.1f)",
                                                    cuts_.espri_aang_width));
    xpos_vs_aaang_1d_total = init_vertz_cut_dist_hist(cuts_.espri_aang_width);
    xpos_vs_aaang_1d_total->SetTitle("aang diff. total");
    xpos_vs_aang_sn_total = init_kin_corr_hist("Total");
    xpos_vs_aang_sn_in = init_kin_corr_hist("In cut");
    xpos_vs_aang_sn_out = init_kin_corr_hist("Out cut");
    xpos_vs_aang_e_corr_in = init_e_cut_dist_hist(cuts_.espri_e_width_stdev);
  }

  void compute_vertz_cut_stats(Evt& evt) {
    if (!is_well_defined(evt.p_theta)) {
      return;
    }
    if (!(is_well_defined(evt.he_theta))) {
      return;
    }

    double aang = cuts_.IsLeftEspriEvent(evt) ? evt.esl_aang : evt.esr_aang;
    double xpos = cuts_.IsLeftEspriEvent(evt) ? evt.esl_xpos : evt.esr_xpos;
    // double nhitplanes = cuts_.IsLeftEspriEvent(evt) ? evt.esl_nhitwires :
    //   evt.esr_nhitwires;
    double theta_diff = evt.he_theta - evt.he_theta_theor;

    if (std::abs(xpos) > 9000) {
      return;
    }

    // if ((aang > -225 && aang < -175)) {
    //   return;
    // }
    // double aang_diff = aang - (xpos/s13::gk_dist_rdc_target * 1e3) - ElaCuts::k_esr_aang_center_offset;
    // bool is_pos_peak = aang_diff > 260 && aang_diff < 280;
    // if (!is_pos_peak) {
    //   return;
    // }

    // if (std::abs((aang - (xpos/s13::gk_dist_rdc_target * 1e3)) - 220) < 10) {
    //   // logger_.debug("aang: %.2f; xpos: %.2f", aang, xpos);
    //   // logger_.debug("nhitwires: %.2f", nhitplanes);

    //   return;
    // }

    auto vz_cuts = cuts_;
    for (size_t i = 0; i < vertz_widths.size(); i++) {
      vz_cuts.espri_aang_width = vertz_widths[i];
      vertz_cut_dists[i]->Fill(aang - (xpos/s13::gk_dist_rdc_target * 1e3));
      vertz_aang_dists[i]->Fill(evt.p_aang);
      vertz_total[i]->Fill(theta_diff);
      kin_corrs_total_2d[i]->Fill(evt.he_theta, evt.p_theta_eff);

      if (vz_cuts.IsEspriAangCut(evt)) {
        vertz_in_cut[i]->Fill(theta_diff);
        kin_corrs_in_cut_2d[i]->Fill(evt.he_theta, evt.p_theta_eff);
      } else {
        vertz_out_cut[i]->Fill(theta_diff);
      }
    }
  }

  bool is_hodf_z2(Evt& evt) {
    double hodf_q = cuts_.get_hodf_q(evt);
    return hodf_q > 7;
  }

  bool is_hodf_he6(Evt& evt) {
    return cuts_.IsHodfHe6Cut(evt);
  }

  bool is_hodf_he4(Evt& evt) {
    return cuts_.get_hodf_hit_id(evt) >= 16 && cuts_.get_hodf_q(evt);
  }

  void compute_xpos_vs_aang_stats(Evt& evt) {
    auto local_cuts = cuts_;
    local_cuts.apply_espri_aang_cut = false;
    if (!local_cuts.IsSignalEvent(evt)) {
      return;
    }
    if (!is_well_defined(evt.p_theta)) {
      return;
    }
    if (!(is_well_defined(evt.he_theta))) {
      return;
    }
    if (!(evt.IsLeftEspriEvent() || evt.IsRightEspriEvent())) {
      return;
    }

    double xpos = evt.IsLeftEspriEvent() ? evt.esl_xpos : evt.esr_xpos;
    double aang = evt.IsLeftEspriEvent() ? evt.esl_aang : evt.esr_aang;

    // use HODF only for he6 data
    // if (xpos < -160 || !is_hodf_z2(evt)) {
    //   return;
    // }

    double aang_diff = aang - (xpos/s13::gk_dist_rdc_target * 1e3);
    bool is_in_cut = aang_diff > 40 && aang_diff < 110;
    double e_corr = (evt.p_e_eff - evt.p_e_sim) / evt.p_e_sim;

    xpos_vs_aang_sn_total->Fill(evt.he_theta - evt.he_theta_theor);
    xpos_vs_aaang_1d_total->Fill(aang_diff);
    xpos_vs_aang_total->Fill(xpos, aang);
    if (is_in_cut) {
      xpos_vs_aang_in->Fill(xpos, aang);
      xpos_vs_aaang_1d_in->Fill(aang_diff);
      xpos_vs_aang_kin_corr->Fill(evt.he_theta, evt.p_theta_eff);
      xpos_vs_aang_sn_in->Fill(evt.he_theta - evt.he_theta_theor);
      xpos_vs_aang_e_corr_in->Fill(e_corr);
    } else {
      xpos_vs_aang_sn_out->Fill(evt.he_theta - evt.he_theta_theor);
    }
  }

  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;

};

void draw_aang_cut_marker_lines(double cut_width, double offset, double ymax) {
  s13::ana::draw_marker_line(cut_width / 2 + offset, 0,
                             cut_width / 2 + offset,
                             ymax);
  s13::ana::draw_marker_line(-cut_width / 2 + offset, 0,
                             -cut_width / 2 + offset,
                             ymax);
}

void draw_e_cut_stats(s13::ana::RootScript& script, std::shared_ptr<VertzCutAlg> alg) {
  double aang_offset = alg->cuts().espri_selector == s13::ana::Espri::right ?
    ElaCuts::k_esr_aang_center_offset : ElaCuts::k_esl_aang_center_offset;
  for (size_t i = 0; i < alg->vertz_widths.size(); i++) {
    if (i % 4 == 0) {
      script.NewPage(4, 3).cd(s13::ana::PadSeq::column);
    }

    script.cd();
    alg->vertz_cut_dists[i]->Draw();
    alg->vertz_cut_dists[i]->Fit("gaus");
    alg->vertz_cut_dists[i]->SetMaximum(alg->vertz_cut_dists[i]->GetMaximum() * 1.5);
    s13::ana::draw_marker_line(alg->vertz_widths[i] / 2 + aang_offset, 0,
                               alg->vertz_widths[i] / 2 + aang_offset,
                               alg->vertz_cut_dists[i]->GetMaximum());
    s13::ana::draw_marker_line(-alg->vertz_widths[i] / 2 + aang_offset, 0,
                               -alg->vertz_widths[i] / 2 + aang_offset,
                               alg->vertz_cut_dists[i]->GetMaximum());
    // s13::ana::draw_text(alg->cuts().GetTotalCutTitle(), 0.15, 0.75, 0.7, 0.9);
    // gPad->SetLogy();

    // script.cd();
    // alg->vertz_aang_dists[i]->Draw();
    // gPad->SetLogy();

    script.cd();
    auto s2 = s13::ana::draw_thstack({alg->vertz_total[i],
          alg->vertz_in_cut[i], alg->vertz_out_cut[i]},
      "1-d kin. corrs.", "NOSTACK");
    s2->SetMaximum(s2->GetMaximum() * 0.8);
    gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);

    // script.cd();
    // alg->vertz_in_cut[i]->Draw();

    auto kin_graph = s13::ana::BuildHe6ProtonKinematicsGraph();
    kin_graph->SetLineWidth(2);
    // script.cd();
    // alg->kin_corrs_total_2d[i]->Draw("col");
    // kin_graph->Draw("SAME");
    script.cd();
    alg->kin_corrs_in_cut_2d[i]->Draw("col");
    kin_graph->Draw("SAME");
  }
  gPad->Update();
}

void draw_xpos_vs_aang_stats(s13::ana::RootScript& script, std::shared_ptr<VertzCutAlg> alg) {
  script.NewPage(2,2).cd(s13::ana::PadSeq::column);
  script.cd();
  s13::ana::draw_thstack({alg->xpos_vs_aang_in, alg->xpos_vs_aang_total},
                         "aang vs. X (total & in cut)", "COL");
  s13::ana::draw_text(alg->cuts_no_aang().GetTotalCutTitle());
  s13::ana::draw_text(s13::ana::tstrfmt("Total: %.1f", alg->xpos_vs_aang_total->Integral()),
                      0.65, 0.25, 0.8, 0.35);
  s13::ana::draw_text(s13::ana::tstrfmt("In cut: %.1f", alg->xpos_vs_aang_in->Integral()),
                      0.65, 0.15, 0.8, 0.25);
  // alg->xpos_vs_aang_total->Draw("COL");
  // s13::ana::draw_text(alg->cuts_no_aang().GetTotalCutTitle());
  gPad->SetLogz();

  script.cd();
  // auto s1 = s13::ana::draw_thstack({alg->xpos_vs_aaang_1d_total, alg->xpos_vs_aaang_1d_in},
  //                                  "aang. diff. distributions", "NOSTACK");
  s13::ana::draw_thstack({alg->xpos_vs_aaang_1d_total},
                         "aang. diff. distributions", "NOSTACK");
  s13::ana::draw_text(s13::ana::tstrfmt("Total: %.1f", alg->xpos_vs_aaang_1d_total->Integral()),
                      0.65, 0.7, 0.8, 0.8);
  s13::ana::draw_text(s13::ana::tstrfmt("In cut: %.1f", alg->xpos_vs_aaang_1d_in->Integral()),
                      0.65, 0.6, 0.8, 0.7);
  gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);
  // draw_aang_cut_marker_lines(alg->cuts().espri_aang_width, ElaCuts::k_esr_aang_center_offset,
  //                            s1->GetMaximum());
  gPad->Update();

  script.cd();
  alg->xpos_vs_aang_kin_corr->Draw("COL");
  auto kin_graph = s13::ana::BuildHe6ProtonKinematicsGraph();
  kin_graph->Draw("SAME");
  gPad->SetLogz();
  s13::ana::draw_text(alg->cuts().GetTotalCutTitle(), 0.15, 0.15, 0.45, 0.25);

  script.cd();
  s13::ana::draw_thstack({alg->xpos_vs_aang_sn_total,
        alg->xpos_vs_aang_sn_in, alg->xpos_vs_aang_sn_out},
    "1-d kin. corrs.", "NOSTACK");
  gPad->BuildLegend(0.65, 0.65, 0.90, 0.85);
  gPad->Update();

  script.NewPage(2,2).cd(s13::ana::PadSeq::column);
  script.cd();
  alg->xpos_vs_aang_e_corr_in->Draw();
}

void vertz_cut() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  auto alg =
    s13::ana::walk_alg<VertzCutAlg>(cuts,
                                    init_dataset_from_cli(cli_opts),
                                    cli_opts.use_mt());

  s13::ana::RootScript script(cli_opts.output_file_basename());
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit();
  // draw_vertz_cut(script, alg);
  draw_e_cut_stats(script, alg);
  draw_xpos_vs_aang_stats(script, alg);
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

  vertz_cut();
  std::cout << "exiting main\n" << std::flush;
}
