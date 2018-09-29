//
// Created by serj on 16/10/25.
//

#include <fstream>
#include <future>

#include "TChain.h"
#include "TGraph.h"
#include "TFrame.h"
#include "Math/Interpolator.h"

#include "init_6he_ds.hpp"
#include "macroincludes.hpp"


using namespace s13::ana;

TChain *g_local_chain = init_6he_partial_dataset(20);

TH1 *draw_theta_corr_in_range(
    TChain &chain, TCut &&cuts, TString title = "",
    TString draw_opts = "colz", int bin_num = 100) {
  static int hist_num = 1;
  TCut p_angle_acceptance_cut = "p_theta_eff*57.3>50 && p_theta_eff*57.3<75";
  TString hist_name = TString::Format("thcorr_%i", hist_num);
  if (title == "") {
    title = cuts;
  }

  TString limits = TString::Format("(%i,0,15,%i,50,75)", bin_num, bin_num);
  TH1 *hist = hdraw(
      chain, hist_name,
      "p_theta_eff*57.3:s1dc_theta*57.3",
      limits,
      cuts,
      title,
      "Fragment #theta [lab. deg.]",
      "Proton #theta [lab. deg.]", draw_opts);

  hist_num += 1;
  return hist;
}

void draw_kin_corr(RootScript &script, ElasticScatteringCuts &cuts,
                   bool do_cut, bool in_cut) {
  TCut vertz_cut = "";
  TString cut_type = "no_cut";
  if (do_cut) {
    vertz_cut = cuts.GetVertexZCut();
    cut_type = "in_cut";
  }
  if (do_cut && !in_cut) {
    vertz_cut = !vertz_cut;
    cut_type = "out_cut";
  }
  TString title = TString::Format("%s %s %s",
                                  cuts.GetTotalCutTitle().Data(),
                                  cut_type.Data(), 
                                  cuts.GetVertexZCutTitle().Data());

  script.cd();
  draw_theta_corr_in_range(*g_local_chain,
                           cuts.GetTotalCut() && vertz_cut,
                           title);
  script.cd();
  draw_theta_corr_1d(*g_local_chain,
                     cuts.GetTotalCut() && vertz_cut,
                     title);
}

void draw_kin_corr_sn(RootScript &script, ElasticScatteringCuts &cuts) {
  TString title = TString::Format("%s", cuts.GetTotalCutTitle().Data());
  TH1 *hist_no_cut = draw_theta_corr_1d(*g_local_chain,
                                        cuts.GetTotalCut(),
                                        title);
  title = TString::Format("%s %s",
                          "in cut: ",
                          cuts.GetVertexZCutTitle().Data());
  TH1 *hist_in_cut =
    draw_theta_corr_1d(*g_local_chain,
                       cuts.GetTotalCut() && cuts.GetVertexZCut(),
                       title.Data());
  title = TString::Format("%s %s",
                          "out cut",
                          cuts.GetVertexZCutTitle().Data());
  TH1 *hist_out_cut =
    draw_theta_corr_1d(*g_local_chain,
                       cuts.GetTotalCut() && !cuts.GetVertexZCut(),
                       title.Data());
  gStyle->SetOptStat(0000000);
  hist_no_cut->Draw();
  hist_in_cut->SetLineColor(kGreen);
  hist_in_cut->Draw("SAME");
  hist_out_cut->SetLineColor(kRed);
  hist_out_cut->Draw("SAME");
  gPad->BuildLegend(0.6, 0.7, 0.9, 0.9);
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

void draw_vertz_cut(RootScript &script, ElasticScatteringCuts &cuts) {
  static int hist_num = 1;
  TString hist_name = TString::Format("vZ_%i", hist_num);
  TString espri_prefix = cuts.espri_selector == Espri::left ? "esl" : "esr";
  TString draw_cmd = TString::Format("%s_vertex_zpos", espri_prefix.Data());
  TString title = TString::Format("vert-Z dist. (%s)",
                                  cuts.GetTotalCutTitle().Data());

  hdraw(*g_local_chain, hist_name,
        draw_cmd,
        "(200,-200,200)",
        cuts.GetTotalCut(),
        title,
        "Z [mm]", "Counts");
  gPad->Update();

  Double_t zpos_offset = cuts.espri_selector == Espri::left ? \
    ElasticScatteringCuts::k_esl_zpos_offset :
    ElasticScatteringCuts::k_esr_zpos_offset;
  draw_marker_line(-cuts.vertexZ_width/2. + zpos_offset, 0,
                   -cuts.vertexZ_width/2. + zpos_offset, gPad->GetFrame()->GetY2());
  draw_marker_line(cuts.vertexZ_width/2. + zpos_offset, 0,
                   cuts.vertexZ_width/2. + zpos_offset, gPad->GetFrame()->GetY2());

  hist_num += 1;
}

void draw_rdcab_range(RootScript &script, ElasticScatteringCuts &cuts) {
  Double_t widths[] = {10, 40, 300};

  for (int i = 1; i <= 3; i++) {
    cuts.vertexZ_width = widths[i-1];
    script.cd();
    draw_vertz_cut(script, cuts);
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
  const char *binning = "(100, -300, 300, 100, -300, 300";

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
  TString title = TString::Format("%s (%s)", "6He",
                                  cuts.GetTotalCutTitle().Data());

  script.cd();
  TH1 *hist_2d = draw_theta_corr_2d(*g_local_chain, cuts.GetTotalCut(), title);
  script.cd();
  TH1 *hist_1d = draw_theta_corr_1d(*g_local_chain, cuts.GetTotalCut(), title);
}

void vertz_check() {
  RootScript script("vertz_cut_check");
  ElasticScatteringCuts cuts;
  cuts.apply_hodf_he6_cut = true;
  cuts.apply_phi_cut = true;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_espri_aang_cut = false;

  script.NewPage(3,2).cd(PadSeq::column);
  cuts.vertexZ_width = 40;
  draw_kin_corr_no_cut(script, cuts);
  draw_kin_corr_in_cut(script, cuts);
  draw_kin_corr_out_cut(script, cuts);

  cuts.espri_selector = Espri::left;
  script.NewPage(3,2).cd(PadSequence::column);
  draw_rdcab_range(script, cuts);

  cuts.espri_selector = Espri::right;
  script.NewPage(3,2).cd(PadSequence::column);
  draw_rdcab_range(script, cuts);

//  script.NewPage(2, 2).cd(PadSeq::column);
//  cuts.espri_selector = Espri::left;
//  draw_rdc_pos_vs_angle(script, cuts);
//  cuts.espri_selector = Espri::right;
//  draw_rdc_pos_vs_angle(script, cuts);
//
//  script.NewPage(2, 2).cd(PadSeq::column);
//  cuts.espri_selector = Espri::left;
//  draw_sn(script, cuts);
//  cuts.espri_selector = Espri::right;
//  draw_sn(script, cuts);
//
//  cuts.apply_espri_aang_cut = true;
//
//  script.NewPage(2, 2).cd(PadSeq::column);
//  cuts.espri_selector = Espri::left;
//  draw_rdc_pos_vs_angle(script, cuts);
//  cuts.espri_selector = Espri::right;
//  draw_rdc_pos_vs_angle(script, cuts);
//
//  script.NewPage(2, 2).cd(PadSeq::column);
//  cuts.espri_selector = Espri::left;
//  draw_sn(script, cuts);
//  cuts.espri_selector = Espri::right;
//  draw_sn(script, cuts);
}
