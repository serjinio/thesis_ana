
#include "TChain.h"
#include "TGraph.h"

#include "init_4he_ds.hpp"
#include "macroincludes.hpp"
#include "drawcuts.h"

using namespace s13::ana;
using namespace s13::ana::algs;


TChain *g_local_chain = init_4he_partial_dataset(25);
RootScript* script = new RootScript("s1dc_pos_dep_eff_4", g_local_chain);


TH1* draw_rel_diff(ElasticScatteringCuts cuts,
                   const char *user_title,
                   TCut fdc0_cut, TCut s1dc_cut,
                   TString binning = "(200, 0, 9)") {
  static int hist_num = 1;

  script->NewPage(2,1).cd(PadSeq::column);
  script->cd();
  TH1* fdc0_hist = hdraw(*g_local_chain, "fdc0_cts", "fdc0_theta*57.3",
                         binning, cuts.GetTotalCut() && fdc0_cut,
                         TString::Format("%s FDC0 #theta (%s)",
                                         user_title,
                                         cuts.GetTotalCutTitle().Data()),
                         "#theta [deg.]", "Counts");
  script->cd();
  TH1* s1dc_hist = hdraw(*g_local_chain, "s1dc_cts", "s1dc_theta*57.3",
                         binning, cuts.GetTotalCut() && s1dc_cut,
                         TString::Format("%s S1DC #theta (%s)",
                                         user_title,
                                         cuts.GetTotalCutTitle().Data()),
                         "#theta [deg.]", "Counts");

  script->NewPage(1,1).cd(PadSeq::column);
  script->cd();
  TH1* rel_diff_hist = static_cast<TH1*>(s1dc_hist->Clone());
  rel_diff_hist->
    SetNameTitle(TString::Format("rel_diff_%i", hist_num),
                 TString::Format("%s; (S1DC - FDC0)/FDC0", user_title).Data());
  rel_diff_hist->GetYaxis()->SetTitle("Rel diff x100 [%]");
  rel_diff_hist->Add(fdc0_hist, -1);
  rel_diff_hist->Divide(fdc0_hist);
  rel_diff_hist->Draw();

  hist_num += 1;
  return rel_diff_hist;
}

void s1dc_pos_dep_eff() {
  gStyle->SetOptFit();
  ElasticScatteringCuts cuts;
  TCut lr_scat_tgt_cut = "(pow(tgt_up_xpos, 2) + pow(tgt_up_ypos, 2) < 9.00)";
  cuts.apply_vertexXY_cut = true;
  cuts.apply_phi_cut = false;
  cuts.vertexXY_radius = 6;
  cuts.apply_user_cut = false;
  // cuts.user_cut = ("abs(fdc0_xpos) < 3 && abs(fdc0_ypos) < 3 "
  //                  "&& abs(tgt_up_aang)*57.3 < 0.3"
  //                  "&& abs(tgt_up_bang)*57.3 < 0.3");

  TString binning = "(200, 0, 13)";


  /// 2D theta distribution
  // script->NewPage(2).cd();
  // cuts.apply_user_cut = true;
  // cuts.user_cut = ("s1dc_xpos > 0 && fdc0_xpos > 0 && s1dc_phi < 0"
  //                  "&& s1dc_phi*57.3 > -200");
  // hdraw(*g_local_chain, "s1dc_fdc0_corr1", "s1dc_theta*57.3:fdc0_theta*57.3",
  //       "(100, 0, 10, 100, 0, 10)", cuts.GetTotalCut(),
  //       TString::Format("S1DC vs. FDC0 #theta (left scat.)"),
  //       "FDC0 #theta [deg]", "S1DC #theta [deg]");
  // script->cd();
  // cuts.user_cut = "s1dc_xpos < 0 && fdc0_xpos < 0 && s1dc_phi > 0";
  // hdraw(*g_local_chain, "s1dc_fdc0_corr2", "s1dc_theta*57.3:fdc0_theta*57.3",
  //       "(100, 0, 10, 100, 0, 10)", cuts.GetTotalCut(),
  //       TString::Format("S1DC vs. FDC0 #theta (right scat.)"),
  //       "FDC0 #theta [deg]", "S1DC #theta [deg]");

  cuts.apply_user_cut = false;
  auto hist_left_diff =
    draw_rel_diff(cuts, "Left scat.",
                  "fdc0_xpos > 0 && fdc0_phi < 0 && fdc0_phi*57.3 > -200",
                  "s1dc_xpos > 0 && s1dc_phi < 0 && s1dc_phi*57.3 > -200");
  auto hist_right_diff =
    draw_rel_diff(cuts, "Right scat.",
                  "fdc0_xpos < 0 && fdc0_phi > 0",
                  "s1dc_xpos < 0 && s1dc_phi > 0");

  script->NewPage(2,1);
  script->cd();
  hist_left_diff->Draw();
  script->cd();
  hist_right_diff->Draw();

  script->NewPage(1).cd();
  hist_left_diff->Add(hist_right_diff, -1);
  hist_left_diff->SetNameTitle("diff_lr", "Diff_left - diff_right");
  hist_left_diff->Draw();


  // script->NewPage(1).cd();
  // TH1* align_hist =
  //   hdraw(*g_local_chain, "s1dc_left_scat1", "s1dc_xpos",
  //         "(600, -100, 100)", cuts.GetTotalCut() &&
  //         "abs(fdc0_xpos) < 2 && abs(tgt_up_aang)*57.3 < 0.3",
  //         TString::Format("S1DC center vs. FDC0 (abs(fdc0_xpos) < 2)"),
  //         "X [mm]", "Counts");
  // align_hist->Fit("gaus", "", "", -30, 30);

  /// phi angles dist
  // script->NewPage(2,1).cd(PadSeq::row);
  // cuts.apply_user_cut = false;
  // script->cd();
  // hdraw(*g_local_chain, "s1dc_phi1", "s1dc_phi*57.3",
  //       "(400, -200, 200)", cuts.GetTotalCut(),
  //       TString::Format("S1DC #phi dist."),
  //       "#phi [deg.]", "Counts");
  // gPad->SetLogy();
  // script->cd();
  // hdraw(*g_local_chain, "fdc0_phi1", "fdc0_phi*57.3",
  //       "(400, -200, 200)", cuts.GetTotalCut(),
  //       TString::Format("FDC0 #phi dist."),
  //       "#phi [deg.]", "Counts");
  // gPad->SetLogy();

  // script->NewPage(2,2).cd(PadSeq::row);
  // script->cd();
  // hdraw(*g_local_chain, "s1dc_left_scat1", "s1dc_xpos",
  //       "(100, -100, 100)", cuts.GetTotalCut() &&
  //       "s1dc_phi < 0 && s1dc_phi*57.3 > -200 && s1dc_xpos > 0",
  //       TString::Format("S1DC left scat. (s1dc_phi < 0)"),
  //       "X [mm]", "Counts");
  // gPad->SetLogy();
  // script->cd();
  // hdraw(*g_local_chain, "s1dc_right_scat1", "s1dc_xpos",
  //       "(100, -100, 100)", cuts.GetTotalCut() &&
  //       "s1dc_phi > 0 && s1dc_xpos < 0",
  //       TString::Format("S1DC right scat. (s1dc_phi > 0)"),
  //       "X [mm]", "Counts");
  // gPad->SetLogy();

  // script->cd();
  // hdraw(*g_local_chain, "s1dc_left_scat2", "s1dc_xpos",
  //       "(100, -100, 100)", cuts.GetTotalCut() &&
  //       "fdc0_phi < 0 && fdc0_phi*57.3 > -200 && fdc0_xpos > 0",
  //       TString::Format("S1DC left scat. (fdc0_phi < 0)"),
  //       "X [mm]", "Counts");
  // gPad->SetLogy();
  // script->cd();
  // hdraw(*g_local_chain, "s1dc_right_scat2", "s1dc_xpos",
  //       "(100, -100, 100)", cuts.GetTotalCut() &&
  //       "fdc0_phi > 0 && fdc0_xpos < 0",
  //       TString::Format("S1DC right scat. (fdc0_phi > 0)"),
  //       "X [mm]", "Counts");
  // gPad->SetLogy();

  // script->NewPage(2,2).cd(PadSeq::row);
  // script->cd();
  // hdraw(*g_local_chain, "s1dc_left_scat1", "s1dc_xpos",
  //       "(100, -100, 100)", cuts.GetTotalCut() &&
  //       "s1dc_xpos > 0",
  //       TString::Format("S1DC left scat. (s1dc_xpos < 0)"),
  //       "X [mm]", "Counts");
  // gPad->SetLogy();
  // script->cd();
  // hdraw(*g_local_chain, "s1dc_right_scat1", "s1dc_xpos",
  //       "(100, -100, 100)", cuts.GetTotalCut() &&
  //       "s1dc_xpos < 0",
  //       TString::Format("S1DC right scat. (s1dc_xpos > 0)"),
  //       "X [mm]", "Counts");
  // gPad->SetLogy();

  // script->cd();
  // hdraw(*g_local_chain, "s1dc_left_scat2", "s1dc_xpos",
  //       "(100, -100, 100)", cuts.GetTotalCut() &&
  //       "fdc0_xpos > 0",
  //       TString::Format("S1DC left scat. (fdc0_xpos < 0)"),
  //       "X [mm]", "Counts");
  // gPad->SetLogy();
  // script->cd();
  // hdraw(*g_local_chain, "s1dc_right_scat2", "s1dc_xpos",
  //       "(100, -100, 100)", cuts.GetTotalCut() &&
  //       "fdc0_xpos < 0",
  //       TString::Format("S1DC right scat. (fdc0_xpos > 0)"),
  //       "X [mm]", "Counts");
  // gPad->SetLogy();

  //hist->Fit("gaus");

  delete script;
}
