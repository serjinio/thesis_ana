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


void draw_kin_corr(RootScript &script, ElasticScatteringCuts &cuts,
                   bool do_cut, bool in_cut) {
  TCut abang_cut = "";
  TString cut_type = "no_cut";
  if (do_cut) {
    abang_cut = cuts.GetEspriAangCut();
    cut_type = "in_cut";
  }
  if (do_cut && !in_cut) {
    abang_cut = !abang_cut;
    cut_type = "out_cut";
  }
  TString title = TString::Format("%s %s %s",
                                  cuts.GetTotalCutTitle().Data(),
                                  cut_type.Data(),
                                  cuts.GetEspriAangCutTitle().Data());

  script.cd();
  draw_theta_corr_2d(*g_local_chain,
                     cuts.GetTotalCut() && abang_cut,
                     title);
  script.cd();
  draw_theta_corr_1d(*g_local_chain,
                     cuts.GetTotalCut() && abang_cut,
                     title);
}

void draw_kin_corr_sn(RootScript &script, ElasticScatteringCuts &cuts) {
  cuts.apply_espri_aang_cut = false;
  TString title = TString::Format("%s", cuts.GetTotalCutTitle().Data());
  TH1 *hist_no_cut = draw_theta_corr_1d(*g_local_chain,
                                        cuts.GetTotalCut(),
                                        title);
  title = TString::Format("%s %s",
                          "IN:",
                          cuts.GetEspriAangCutTitle().Data());
  TH1 *hist_in_cut =
    draw_theta_corr_1d(*g_local_chain,
                       cuts.GetTotalCut() && cuts.GetEspriAangCut(),
                       title.Data());
  title = TString::Format("%s %s",
                          "OUT:",
                          cuts.GetEspriAangCutTitle().Data());
  TH1 *hist_out_cut =
    draw_theta_corr_1d(*g_local_chain,
                       cuts.GetTotalCut() && !cuts.GetEspriAangCut(),
                       title.Data());

  gStyle->SetOptStat(1111100);
  hist_no_cut->Draw();
  hist_in_cut->SetLineColor(kGreen);
  hist_in_cut->Draw("SAME");
  hist_out_cut->SetLineColor(kRed);
  hist_out_cut->Draw("SAME");
  gPad->BuildLegend(0.6, 0.3, 0.95, 0.6);
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

void draw_abang_cut(RootScript &script, ElasticScatteringCuts &cuts) {
  static int hist_num = 1;
  TString prefix = cuts.espri_selector == Espri::left ? "esl" : "esr";
  TString draw_cmd_aang = TString::Format("%s_aang - %s_xpos/1004.75 * 1e3",
                                          prefix.Data(), prefix.Data());
  TString draw_cmd_bang = TString::Format("%s_bang - %s_ypos/1004.75 * 1e3",
                                          prefix.Data(), prefix.Data());
  TString hist_aang_name = TString::Format("%s_xvsa_%i", prefix.Data(), hist_num);
  TString hist_bang_name = TString::Format("%s_yvsb_%i", prefix.Data(), hist_num);
  const char *binning = "(200, -1000, 1000)";
  TString title = TString::Format("abs(a/b_ang - a/b_ang_ref) (%s)",
                                  cuts.GetTotalCutTitle().Data());

  script.cd();
  hdraw(*g_local_chain, hist_aang_name,draw_cmd_aang, binning, cuts.GetTotalCut(),
        title, "aang - aang_ref [mrad]", "Counts");
  gPad->Update();
  Double_t aang_offset = cuts.espri_selector == Espri::left ?
                         ElasticScatteringCuts::k_esl_aang_center_offset :
                         ElasticScatteringCuts::k_esr_aang_center_offset;
  draw_marker_line(-cuts.espri_aang_width/2. + aang_offset, 0,
                   -cuts.espri_aang_width/2. + aang_offset, gPad->GetFrame()->GetY2());
  draw_marker_line(cuts.espri_aang_width/2. + aang_offset, 0,
                   cuts.espri_aang_width/2. + aang_offset, gPad->GetFrame()->GetY2());
  gPad->SetLogy();

  // script.cd();
  // hdraw(*g_local_chain, hist_bang_name, draw_cmd_bang, binning, cuts.GetTotalCut(),
  //       title, "bang - bang_ref [mrad]", "Counts");
  // gPad->Update();
  // Double_t bang_offset = cuts.espri_selector == Espri::left ?
  //                        ElasticScatteringCuts::k_esl_bang_center_offset :
  //                        ElasticScatteringCuts::k_esr_bang_center_offset;
  // draw_marker_line(-cuts.espri_aang_width/2. + bang_offset, 0,
  //                  -cuts.espri_aang_width/2. + bang_offset, gPad->GetFrame()->GetY2());
  // draw_marker_line(cuts.espri_aang_width/2. + bang_offset, 0,
  //                  cuts.espri_aang_width/2. + bang_offset, gPad->GetFrame()->GetY2());

  hist_num += 1;
}

void draw_rdcab_range(RootScript &script, ElasticScatteringCuts &cuts) {
  Double_t widths[] = {90, 360, 2000};

  for (int i = 1; i <= 3; i++) {
    cuts.espri_aang_width = widths[i-1];
    draw_abang_cut(script, cuts);
    script.cd();
    draw_kin_corr_sn(script, cuts);
  }
}

void draw_rdc_pos_vs_angle(RootScript &script, ElasticScatteringCuts &cuts) {
  assert(cuts.espri_selector != Espri::both);
  static int hist_num = 1;
  TString prefix = cuts.espri_selector == Espri::left ? "esl" : "esr";
  TString draw_cmd_aang = TString::Format("%s_aang:%s_xpos", prefix.Data(), prefix.Data());
  TString draw_cmd_bang = TString::Format("%s_bang:%s_ypos", prefix.Data(), prefix.Data());
  TString hist_aang_name = TString::Format("%s_xvsa_%i", prefix.Data(), hist_num);
  TString hist_bang_name = TString::Format("%s_yvsb_%i", prefix.Data(), hist_num);
  const char *binning = "(200, -300, 300, 200, -300, 300)";

  script.cd();
  hdraw(*g_local_chain, hist_aang_name, draw_cmd_aang, binning,
        cuts.GetTotalCut(), cuts.GetTotalCutTitle(),
        "xpos [mm]", "aang [mrad]");

  script.cd();
  hdraw(*g_local_chain, hist_bang_name, draw_cmd_bang, binning,
        cuts.GetTotalCut(), cuts.GetTotalCutTitle(),
        "ypos [mm]", "bang [mrad]");

  hist_num += 1;
}

void draw_sn(RootScript &script, ElasticScatteringCuts &cuts) {
  assert(cuts.espri_selector != Espri::both);
  TString prefix = cuts.espri_selector == Espri::left ? "esl" : "esr";
  TString title = TString::Format("%s (%s)", "4He",
                                  cuts.GetTotalCutTitle().Data());

  script.cd();
  TH1 *hist_2d = draw_theta_corr_2d(*g_local_chain, cuts.GetTotalCut(), title);
  script.cd();
  TH1 *hist_1d = draw_theta_corr_1d(*g_local_chain, cuts.GetTotalCut(), title);
}

void draw_cut_and_kin_corr(RootScript &script, ElasticScatteringCuts &cuts, Double_t abang_width) {
  cuts.espri_aang_width = abang_width;

  cuts.apply_espri_aang_cut = true;
  script.NewPage(2, 2).cd(PadSequence::column);
  cuts.espri_selector = Espri::left;
  draw_rdc_pos_vs_angle(script, cuts);
  cuts.espri_selector = Espri::right;
  draw_rdc_pos_vs_angle(script, cuts);

  cuts.apply_espri_aang_cut = false;
  cuts.espri_selector = Espri::both;
  script.NewPage(3,2).cd(PadSequence::column);
  draw_kin_corr_no_cut(script, cuts);
  draw_kin_corr_in_cut(script, cuts);
  draw_kin_corr_out_cut(script, cuts);
}

void rdc_abang_check() {
  RootScript script("rdc_abang_check");
  ElasticScatteringCuts cuts;
  cuts.apply_hodf_he6_cut = false;
  cuts.apply_phi_cut = true;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_espri_aang_cut = false;

  cuts.espri_selector = Espri::left;
  script.NewPage(3,2).cd(PadSequence::column);
  draw_rdcab_range(script, cuts);

  cuts.espri_selector = Espri::right;
  script.NewPage(3,2).cd(PadSequence::column);
  draw_rdcab_range(script, cuts);

  // draw_cut_and_kin_corr(script, cuts, 180);
  // draw_cut_and_kin_corr(script, cuts, 90);
  // draw_cut_and_kin_corr(script, cuts, 30);

  // cuts.apply_espri_aang_cut = false;

  // script.NewPage(2, 2).cd(PadSeq::column);
  // cuts.espri_selector = Espri::left;
  // draw_rdc_pos_vs_angle(script, cuts);
  // cuts.espri_selector = Espri::right;
  // draw_rdc_pos_vs_angle(script, cuts);

  // script.NewPage(2, 2).cd(PadSeq::column);
  // cuts.espri_selector = Espri::left;
  // draw_sn(script, cuts);
  // cuts.espri_selector = Espri::right;
  // draw_sn(script, cuts);
}

