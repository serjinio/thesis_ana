//
// Created by serj on 16/10/27.
//

#include <fstream>
#include <future>

#include "TChain.h"
#include "TGraph.h"
#include "TFrame.h"

#include "init_4he_ds.hpp"
#include "macroincludes.hpp"


using namespace s13::ana;


TChain *g_local_chain = init_4he_partial_dataset(5);


void show_sn(RootScript &script, ElasticScatteringCuts &cuts);

void draw_espri_e_de_1d(ElasticScatteringCuts &cuts, Double_t espri_e_cut_width_marker) {
  assert(cuts.espri_selector != Espri::both);

  static int hist_num = 1;
  TString title = TString::Format("dE-E (%s)",
                                  cuts.GetTotalCutTitle().Data());
  TString hist_name = TString::Format("de_e_1d%i", hist_num);
  TString draw_cmd = TString::Format("(p_e_eff - p_e_sim) / p_e_sim");
  TCut user_cut = (TString::Format("p_de > 1 && p_e_eff > 5")).Data();
  //user_cut = "";

  TH1* hist = hdraw(*g_local_chain, hist_name, draw_cmd, "(100,-2.0,2.0)",
                    cuts.GetTotalCut() && user_cut, title, "(p_e_eff - p_e_sim)/p_e_sim (x100%)",
                    "Counts", "");

  if (espri_e_cut_width_marker > 0) {
    Double_t cut_line_x = espri_e_cut_width_marker / 2;
    Double_t offset = cuts.espri_selector == Espri::left ?
                      ElasticScatteringCuts::k_esl_e_dist_center_offset :
                      ElasticScatteringCuts::k_esr_e_dist_center_offset;
    gPad->Update();
    draw_marker_line(-cut_line_x + offset, 0, -cut_line_x + offset, gPad->GetUymax());
    draw_marker_line(cut_line_x + offset, 0, cut_line_x + offset, gPad->GetUymax());
  }

  ++hist_num;
}

void draw_kin_corr(RootScript &script, ElasticScatteringCuts &cuts,
                   bool do_cut, bool in_cut) {
  cuts.apply_espri_e_cut = false;
  TCut espri_e_cut = "";
  TString cut_type = "no_cut";
  if (do_cut) {
    espri_e_cut = cuts.GetEspriECut();
    cut_type = "in_cut";
  }
  if (do_cut && !in_cut) {
    espri_e_cut = !espri_e_cut;
    cut_type = "out_cut";
  }
  TString title = TString::Format("%s %s %s",
                                  cuts.GetTotalCutTitle().Data(),
                                  cut_type.Data(),
                                  cuts.GetEspriECutTitle().Data());

  script.cd();
  draw_theta_corr_2d(*g_local_chain,
                     cuts.GetTotalCut() && espri_e_cut,
                     title);
  script.cd();
  draw_theta_corr_1d(*g_local_chain,
                     cuts.GetTotalCut() && espri_e_cut,
                     title);
}

void draw_kin_corr_sn(RootScript &script, ElasticScatteringCuts &cuts) {
  cuts.apply_espri_e_cut = false;
  TString title = TString::Format("%s", cuts.GetTotalCutTitle().Data());
  TH1 *hist_no_cut = draw_theta_corr_1d(*g_local_chain,
                                        cuts.GetTotalCut(),
                                        title);
  title = TString::Format("%s %s",
                          "IN:",
                          cuts.GetEspriECutTitle().Data());
  TH1 *hist_in_cut =
      draw_theta_corr_1d(*g_local_chain,
                         cuts.GetTotalCut() && cuts.GetEspriECut(),
                         title.Data());
  title = TString::Format("%s %s",
                          "OUT:",
                          cuts.GetEspriECutTitle().Data());
  TH1 *hist_out_cut =
      draw_theta_corr_1d(*g_local_chain,
                         cuts.GetTotalCut() && !cuts.GetEspriECut(),
                         title.Data());

  gStyle->SetOptStat(1111100);
  hist_no_cut->Draw();
  hist_in_cut->SetLineColor(kGreen);
  hist_in_cut->Draw("SAME");
  hist_out_cut->SetLineColor(kRed);
  hist_out_cut->Draw("SAME");
  gPad->BuildLegend(0.6, 0.3, 0.98, 0.6);
}

void draw_kin_corr_in_cut(RootScript &script, ElasticScatteringCuts &cuts) {
  draw_kin_corr(script, cuts, true, true);
}

void draw_kin_corr_out_cut(RootScript &script, ElasticScatteringCuts &cuts) {
  draw_kin_corr(script, cuts, true, false);
}

void draw_kin_corr_no_cut(RootScript &script, ElasticScatteringCuts &cuts) {
  draw_kin_corr(script, cuts, false, false);
}

void draw_cut_range(RootScript &script, ElasticScatteringCuts &cuts) {
    Double_t widths[] = {0.35, 0.7, 0.95};

    for (int i = 1; i <= 3; i++) {
      cuts.espri_e_width = widths[i-1];
      script.cd();
      draw_espri_e_de_1d(cuts, cuts.espri_e_width);
      script.cd();
      draw_kin_corr_sn(script, cuts);
    }
}

void show_minE_bad_events(RootScript &script, ElasticScatteringCuts cuts) {
  cuts.apply_espri_min_e_de_cut = false;
  cuts.apply_espri_e_cut = false;
  cuts.apply_user_cut = true;

  script.NewPage(2, 2).cd(PadSequence::column);

  cuts.espri_selector = Espri::left;
  script.cd();
  cuts.user_cut = "(p_e_eff - p_e_sim) / p_e_sim < -0.98";
  draw_theta_corr_2d(*g_local_chain, cuts.GetTotalCut(), cuts.GetTotalCutTitle());
  script.cd();
  TH1* hist_out_cut = draw_theta_corr_1d(*g_local_chain, cuts.GetTotalCut(),
                                         cuts.GetTotalCutTitle());
  cuts.user_cut = "(p_e_eff - p_e_sim) / p_e_sim > -0.98";
  TH1* hist_in_cut = draw_theta_corr_1d(*g_local_chain, cuts.GetTotalCut(),
                                        cuts.GetTotalCutTitle());
  hist_in_cut->SetLineColor(kGreen);
  hist_out_cut->SetLineColor(kRed);
  hist_in_cut->Draw();
  hist_out_cut->Draw("SAME");

  cuts.espri_selector = Espri::right;
  script.cd();
  cuts.user_cut = "(p_e_eff - p_e_sim) / p_e_sim < -0.98";
  draw_theta_corr_2d(*g_local_chain, cuts.GetTotalCut(), cuts.GetTotalCutTitle());
  script.cd();
  hist_out_cut = draw_theta_corr_1d(*g_local_chain, cuts.GetTotalCut(),
                                         cuts.GetTotalCutTitle());
  cuts.user_cut = "(p_e_eff - p_e_sim) / p_e_sim > -0.98";
  hist_in_cut = draw_theta_corr_1d(*g_local_chain, cuts.GetTotalCut(),
                                        cuts.GetTotalCutTitle());
  hist_in_cut->SetLineColor(kGreen);
  hist_out_cut->SetLineColor(kRed);
  hist_in_cut->Draw();
  hist_out_cut->Draw("SAME");
}

void show_minE_bad_events_2(RootScript &script, ElasticScatteringCuts cuts) {
  static int hist_num = 0;
  cuts.apply_espri_min_e_de_cut = false;
  cuts.apply_espri_e_cut = false;
  cuts.apply_user_cut = true;
  cuts.user_cut = "(p_e_eff - p_e_sim) / p_e_sim > -0.98";

  script.NewPage(1, 2).cd(PadSequence::column);

  cuts.espri_selector = Espri::left;
  script.cd();
  TH1* hist = hdraw(*g_local_chain,
                    TString::Format("minEb%i", hist_num),
                    "p_theta_eff*57.3",
                    "(100, 50, 80)",
                    cuts.GetTotalCut(),
                    "p_theta (" + cuts.GetTotalCutTitle() + ")",
                    "Proton theta [deg.]",
                    "Counts", "colz");
  hist_num += 1;

  cuts.espri_selector = Espri::right;
  script.cd();
  hist = hdraw(*g_local_chain,
               TString::Format("minEb%i", hist_num),
               "p_theta_eff*57.3",
               "(100, 50, 80)",
               cuts.GetTotalCut(),
               "p_theta (" + cuts.GetTotalCutTitle() + ")",
               "Proton theta [deg.]",
               "Counts", "colz");
  hist_num += 1;
}

void espri_e_cut_check() {
  RootScript script("espri_e_cut_check");
  ElasticScatteringCuts cuts;
  cuts.apply_hodf_he6_cut = false;
  cuts.apply_phi_cut = true;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_espri_aang_cut = true;
  cuts.apply_espri_min_e_de_cut = true;
  cuts.apply_espri_e_cut = false;

  // cuts.espri_selector = Espri::left;
  // script.NewPage(3, 2).cd(PadSequence::column);
  // draw_cut_range(script, cuts);

  // cuts.espri_selector = Espri::right;
  // script.NewPage(3, 2).cd(PadSequence::column);
  // draw_cut_range(script, cuts);

  show_minE_bad_events_2(script, cuts);
  
  //show_sn(script, cuts);
}

void show_sn(RootScript &script, ElasticScatteringCuts &cuts) {
  script.NewPage(2, 2).cd(PadSequence::column);
  cuts.apply_espri_e_cut = false;
  cuts.espri_selector = Espri::left;
  script.cd();
  draw_theta_corr_2d(*g_local_chain, cuts.GetTotalCut(), cuts.GetTotalCutTitle());
  script.cd();
  draw_theta_corr_1d(*g_local_chain, cuts.GetTotalCut(), cuts.GetTotalCutTitle());
  cuts.espri_selector = Espri::right;
  script.cd();
  draw_theta_corr_2d(*g_local_chain, cuts.GetTotalCut(), cuts.GetTotalCutTitle());
  script.cd();
  draw_theta_corr_1d(*g_local_chain, cuts.GetTotalCut(), cuts.GetTotalCutTitle());

  script.NewPage(2, 2).cd(PadSequence::column);
  cuts.apply_espri_e_cut = true;
  cuts.espri_e_width = 0.7;
  cuts.espri_selector = Espri::left;
  script.cd();
  draw_theta_corr_2d(*g_local_chain, cuts.GetTotalCut(), cuts.GetTotalCutTitle());
  script.cd();
  draw_theta_corr_1d(*g_local_chain, cuts.GetTotalCut(), cuts.GetTotalCutTitle());
  cuts.espri_selector = Espri::right;
  script.cd();
  draw_theta_corr_2d(*g_local_chain, cuts.GetTotalCut(), cuts.GetTotalCutTitle());
  script.cd();
  draw_theta_corr_1d(*g_local_chain, cuts.GetTotalCut(), cuts.GetTotalCutTitle());
}
