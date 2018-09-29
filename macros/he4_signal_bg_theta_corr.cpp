
#include "TChain.h"
#include "TGraph.h"

#include "init_ds.hpp"
#include "macroincludes.hpp"
#include "drawcuts.h"


using namespace s13::ana;

TChain *g_local_chain = init_6he_partial_dataset(60);
//TChain *g_local_chain = &g_chain_total;


// phi cut defined as 3 * sigma_phi_corr
const Double_t gk_phi_cut_width = 37.8;
// center of the phi_corr: p_theta - He_theta [deg]
// ideally should be 180 [deg]
const Double_t gk_phi_corr_center = 177.5;


void draw_theta_corr_range(RootScript& script, Double_t theta_p_start,
                           Double_t theta_p_step = 2.) {
  static int hist_num = 1;
  int pad_num = 1;
  auto p_interp = MakeBgNormFactorInterpolator();

  ElasticScatteringCuts cuts;
  cuts.espri_selector = Espri::both;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_phi_cut = true;
  cuts.apply_espri_e_cut = false;
  cuts.apply_espri_aang_cut = true;
  cuts.apply_user_cut = true;
  script.NewPage(4,2);
  for (Double_t p_theta = theta_p_start;
       p_theta <= theta_p_start + theta_p_step * 3;
       p_theta += theta_p_step) {
    Double_t p_theta_mean = (p_theta + theta_p_step + p_theta) / 2.0;
    TH1 *hist;
    cuts.user_cut = TString::Format("p_theta*57.3 > %.1f && p_theta*57.3<  %.1f",
                                    p_theta, p_theta + theta_p_step).Data();

    script.cd(pad_num);
    TString title = TString::Format("Kin. corr @ p_{#theta} = %.1f", p_theta_mean);
    draw_theta_corr_2d(*g_local_chain, cuts.GetTotalCut(), title);

    script.cd(pad_num + 4);
    hist = draw_theta_corr_1d(*g_local_chain, cuts.GetTotalCut(),
                              "IN & OUT #phi cut");
    hist->SetTitle("IN #phi");

    hist = draw_theta_corr_1d(*g_local_chain, cuts.GetTotalPhiBgCut(),
                              "OUT #phi", "SAME");
    //hist->SetTitle("BG unscaled");
    hist->SetLineColor(kRed);

    TH1 *hist_bg_scaled = (TH1*)hist->Clone();
    hist_bg_scaled->SetTitle("OUT #phi scaled");
    std::cout << "Scaling BG data by: R - 1 = " <<
      p_interp->Eval(p_theta_mean) - 1 << std::endl;
    hist_bg_scaled->Scale(p_interp->Eval(p_theta_mean) - 1);
    hist_bg_scaled->SetLineColor(kGreen);
    hist_bg_scaled->Draw("SAME");

    gPad->BuildLegend(0.15, 0.67, 0.45, 0.88);
    gPad->Update();

    //hist->GetYaxis()->SetRangeUser(0, 500);

    pad_num += 1;
    hist_num += 1;
  }
}

void he4_signal_bg_theta_corr() {
  RootScript script("he6_signal_bg_theta_corr_check", g_local_chain);
  ElasticScatteringCuts cuts;
  cuts.espri_selector = Espri::both;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_phi_cut = true;
  cuts.apply_espri_aang_cut = true;
  cuts.apply_hodf_he6_cut = true;
  cuts.apply_theta_cut = false;

  cuts.phi_width = gk_phi_cut_width;

  //draw_theta_corr_range(script, 54, 4);

  script.NewPage(4, 2).cd(PadSeq::column);
  s13::cuts::draw_he4_bg_range(script, cuts, 4, 55, 63);
  script.NewPage(4, 2).cd(PadSeq::column);
  s13::cuts::draw_he4_bg_range(script, cuts, 4, 63, 71);
  // script.NewPage(5, 2).cd(PadSeq::column);
  // s13::cuts::draw_he4_bg_range(script, cuts, 5, 65, 70);
}
