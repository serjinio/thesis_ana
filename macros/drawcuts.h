//
// Created by serj on 16/11/07.
//

#ifndef ANA_DRAWCUTS_H
#define ANA_DRAWCUTS_H

#include <TH1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TChain.h>
#include <TFrame.h>
#include "common.hpp"
#include "elacuts.hpp"

namespace s13 {
namespace cuts {

  constexpr int gk_default_bin_num = 48;

  // utility functions

  TCut make_proton_theta_eff_cut(Double_t theta_min_deg, Double_t theta_max_deg) {
    return TString::Format("p_theta_eff > %.5f && p_theta_eff < %.5f",
                           theta_min_deg * s13::ana::D2R,
                           theta_max_deg * s13::ana::D2R).Data();
  }

  template <typename FDist, typename FSn>
  void draw_cut(s13::ana::RootScript& script, s13::ana::ElasticScatteringCuts cuts, TString title,
                FDist fn_cut_dist, FSn fn_sn) {
    script.cd();
    fn_cut_dist(script.GetChain(), cuts, title);
    script.cd();
    fn_sn(script.GetChain(), cuts);
  }

  template <typename FDrCut>
  void draw_cut_range(s13::ana::RootScript& script, s13::ana::ElasticScatteringCuts cuts,
                      int bins_num, FDrCut fn_draw_cut) {
    Double_t bin_width = (70. - 55.) / bins_num;
    TCut old_user_cut = cuts.user_cut;

    for (int i = 0; i < bins_num; ++i) {
      Double_t theta_min = 55 + (i * bin_width);
      Double_t theta_max = theta_min + bin_width;
      TString title = TString::Format("%.2f<p_theta<%.2f", theta_min, theta_max);
      cuts.apply_user_cut = true;
      cuts.user_cut = old_user_cut && make_proton_theta_eff_cut(theta_min, theta_max);

      fn_draw_cut(script, cuts, title);
    }

    cuts.user_cut = old_user_cut;
  }

  // RDC aang cut drawing

  // TODO: modify scattree to include aang - common for both ESPRIs, then update this cut
  TH1* draw_rdc_aang_dist(TChain *ttree, s13::ana::ElasticScatteringCuts cuts,
                          TString title = "") {
    // remove the cut which dist we want to show
    cuts.apply_espri_aang_cut = false;
    // remove theta cut to show the same stats as in sn function
    cuts.apply_theta_cut = false;

    static int hist_num = 1;

    TString prefix = cuts.espri_selector == s13::ana::Espri::left ? "esl" : "esr";
    TString draw_cmd_aang = TString::Format("%s_aang - %s_xpos/1004.75 * 1e3",
                                            prefix.Data(), prefix.Data());
    TString hist_aang_name = TString::Format("%s_xvsa_%i",
                                             prefix.Data(), hist_num);
    const char *binning = "(200, -200, 200)";
    title = TString::Format("abs(a_ang - a_ang_ref) %s (%s)",
                            title.Data(), cuts.GetTotalCutTitle().Data());

    TH1* hist = s13::ana::hdraw(*ttree, hist_aang_name, draw_cmd_aang,
                                binning, cuts.GetTotalCut(),
                                title, "aang - aang_ref [mrad]", "Counts");
    gPad->Update();
    Double_t aang_offset = cuts.espri_selector == s13::ana::Espri::left ?
                           s13::ana::ElasticScatteringCuts::k_esl_aang_center_offset :
                           s13::ana::ElasticScatteringCuts::k_esr_aang_center_offset;
    s13::ana::draw_marker_line(-cuts.espri_aang_width/2. + aang_offset, 0,
                               -cuts.espri_aang_width/2. + aang_offset,
                               gPad->GetFrame()->GetY2());
    s13::ana::draw_marker_line(cuts.espri_aang_width/2. + aang_offset, 0,
                               cuts.espri_aang_width/2. + aang_offset,
                               gPad->GetFrame()->GetY2());

    hist_num += 1;
    return hist;
  }

  TH1* draw_rdc_aang_sn(TChain *ttree, s13::ana::ElasticScatteringCuts cuts) {
    // remove the cut which s/n we want to show
    cuts.apply_espri_aang_cut = false;
    // remove theta cut as otherwise hist will not be informative
    cuts.apply_theta_cut = false;

    TString title = TString::Format("%s", "NO");
    TH1 *hist_no_cut = s13::ana::draw_theta_corr_1d(*ttree,
                                                    cuts.GetTotalCut(),
                                                    title);
    title = TString::Format("%s %s",
                            "IN:",
                            cuts.GetEspriAangCutTitle().Data());
    TH1 *hist_in_cut =
        s13::ana::draw_theta_corr_1d(*ttree,
                                     cuts.GetTotalCut() && cuts.GetEspriAangCut(),
                                     title.Data());
    title = TString::Format("%s %s",
                            "OUT:",
                            cuts.GetEspriAangCutTitle().Data());
    TH1 *hist_out_cut =
        s13::ana::draw_theta_corr_1d(*ttree,
                                     cuts.GetTotalCut() && !cuts.GetEspriAangCut(),
                                     title.Data());

    gStyle->SetOptStat(1111100);
    hist_no_cut->Draw();
    hist_in_cut->SetLineColor(kGreen);
    hist_in_cut->Draw("SAME");
    hist_out_cut->SetLineColor(kRed);
    hist_out_cut->Draw("SAME");
    gPad->BuildLegend(0.6, 0.3, 0.98, 0.6);

    return hist_in_cut;
  }

  void draw_rdc_aang_cut(s13::ana::RootScript& script, s13::ana::ElasticScatteringCuts cuts,
                         TString title = "") {
    draw_cut(script, cuts, title, draw_rdc_aang_dist, draw_rdc_aang_sn);
  }

  void draw_rdc_aang_cut_range(s13::ana::RootScript& script, s13::ana::ElasticScatteringCuts cuts,
                               int bins_num) {
    draw_cut_range(script, cuts, bins_num, draw_rdc_aang_cut);
  }

  // ESPRI E cut drawing

  TH1* draw_espri_e_dist(TChain *ttree, s13::ana::ElasticScatteringCuts cuts,
                         TString title = "") {
    // TODO: check if this is correct to run on both espris
    //assert(cuts.espri_selector != s13::ana::Espri::both);

    // remove the cut which dist we want to show
    cuts.apply_espri_e_cut = false;
    // remove theta cut as to show the same stat as in sn function
    cuts.apply_theta_cut = false;

    static int hist_num = 1;

    title = TString::Format("dE-E %s (%s)",
                            title.Data(), cuts.GetTotalCutTitle().Data());
    TString hist_name = TString::Format("de_e_1d%i", hist_num);
    TString draw_cmd = TString::Format("(p_e_eff - p_e_sim) / p_e_sim");
    // TODO: check its effect
    TCut user_cut = (TString::Format("p_de > 1 && p_e_eff > 5")).Data();
    //user_cut = "";

    TH1* hist = s13::ana::hdraw(*ttree, hist_name, draw_cmd, "(100,-2.0,2.0)",
                                cuts.GetTotalCut() && user_cut, title,
                                "(p_e_eff - p_e_sim)/p_e_sim (x100%)",
                                "Counts", "");

      Double_t cut_line_x = cuts.espri_e_width / 2;
      Double_t offset = cuts.espri_selector == s13::ana::Espri::left ?
                        s13::ana::ElasticScatteringCuts::k_esl_e_dist_center_offset :
                        s13::ana::ElasticScatteringCuts::k_esr_e_dist_center_offset;
      gPad->Update();
      s13::ana::draw_marker_line(-cut_line_x + offset, 0,
                                 -cut_line_x + offset, gPad->GetUymax());
      s13::ana::draw_marker_line(cut_line_x + offset, 0,
                                 cut_line_x + offset, gPad->GetUymax());

    ++hist_num;
    return hist;
  }

  TH1* draw_espri_e_cut_sn(TChain *ttree, s13::ana::ElasticScatteringCuts cuts) {
    // remove the cut which s/n we want to show
    cuts.apply_espri_e_cut = false;
    // remove theta cut as otherwise hist will not be informative
    cuts.apply_theta_cut = false;

    TString title = TString::Format("%s", "NO");
    TH1 *hist_no_cut = s13::ana::draw_theta_corr_1d(*ttree,
                                          cuts.GetTotalCut(),
                                          title);
    title = TString::Format("%s %s",
                            "IN:",
                            cuts.GetEspriECutTitle().Data());
    TH1 *hist_in_cut =
        s13::ana::draw_theta_corr_1d(*ttree,
                           cuts.GetTotalCut() && cuts.GetEspriECut(),
                           title.Data());
    title = TString::Format("%s %s",
                            "OUT:",
                            cuts.GetEspriECutTitle().Data());
    TH1 *hist_out_cut =
        s13::ana::draw_theta_corr_1d(*ttree,
                           cuts.GetTotalCut() && !cuts.GetEspriECut(),
                           title.Data());

    gStyle->SetOptStat(1111100);
    hist_no_cut->Draw();
    hist_in_cut->SetLineColor(kGreen);
    hist_in_cut->Draw("SAME");
    hist_out_cut->SetLineColor(kRed);
    hist_out_cut->Draw("SAME");
    gPad->BuildLegend(0.6, 0.3, 0.98, 0.6);

    return hist_in_cut;
  }

  void draw_espri_e_cut(s13::ana::RootScript& script,
                        s13::ana::ElasticScatteringCuts cuts,
                        TString title = "") {
    draw_cut(script, cuts, title, draw_espri_e_dist, draw_espri_e_cut_sn);
  }

  void draw_espri_e_cut_range(s13::ana::RootScript& script,
                              s13::ana::ElasticScatteringCuts cuts,
                               int bins_num) {
    draw_cut_range(script, cuts, bins_num, draw_espri_e_cut);
  }

  // phi cut drawing

  TH1* draw_phi_dist(TChain *ttree, s13::ana::ElasticScatteringCuts cuts,
                     TString title = "") {
    // remove the cut which dist we want to show
    cuts.apply_phi_cut = false;
    // remove theta cut to show the same stats as in sn function
    cuts.apply_theta_cut = false;

    static int hist_num = 0;

    TString hist_name = TString::Format("phi_dist_%i", hist_num);
    TH1* hist = s13::ana::hdraw(
        g_chain_up, hist_name,
        "abs(p_phi*57.3-s1dc_phi*57.3)",
        "(300,150,210)",
        cuts.GetTotalCut(),
        "abs(p_{#phi} - He_{#phi}) " + title + "(" + cuts.GetTotalCutTitle() + ")",
        "p_{#phi} - He_{#phi} [deg]",
        "Counts");
    gPad->Update();

    Double_t phi_corr_center = s13::ana::ElasticScatteringCuts::k_phi_dist_center_deg;
    s13::ana::draw_marker_line(phi_corr_center - -cuts.phi_width/2.,
                               gPad->GetFrame()->GetY1(),
                               phi_corr_center - -cuts.phi_width/2.,
                               gPad->GetFrame()->GetY2());
    s13::ana::draw_marker_line(phi_corr_center - cuts.phi_width/2.,
                               gPad->GetFrame()->GetY1(),
                               phi_corr_center - cuts.phi_width/2.,
                               gPad->GetFrame()->GetY2());

    hist_num += 1;
    return hist;
  }

  TH1* draw_phi_cut_sn(TChain *ttree, s13::ana::ElasticScatteringCuts cuts) {
    // remove the cut which S/N we want to show
    cuts.apply_phi_cut = false;
    // remove theta cut as otherwise the hist will not be informative
    cuts.apply_theta_cut = false;

    TString title = TString::Format("%s", "NO");
    TH1 *hist_no_cut = s13::ana::draw_theta_corr_1d(*ttree,
                                                    cuts.GetTotalCut(),
                                                    title);
    title = TString::Format("%s %s",
                            "IN:",
                            cuts.GetPhiCutTitle().Data());
    TH1 *hist_in_cut =
        s13::ana::draw_theta_corr_1d(*ttree,
                                     cuts.GetTotalCut() && cuts.GetPhiCut(),
                                     title.Data());
    title = TString::Format("%s %s",
                            "OUT:",
                            cuts.GetPhiCutTitle().Data());
    TH1 *hist_out_cut =
        s13::ana::draw_theta_corr_1d(*ttree,
                                     cuts.GetTotalCut() && !cuts.GetPhiCut(),
                                     title.Data());

    gStyle->SetOptStat(1111100);
    hist_no_cut->Draw();
    hist_in_cut->SetLineColor(kGreen);
    hist_in_cut->Draw("SAME");
    hist_out_cut->SetLineColor(kRed);
    hist_out_cut->Draw("SAME");
    gPad->BuildLegend(0.6, 0.3, 0.98, 0.6);

    return hist_in_cut;
  }

  void draw_phi_cut(s13::ana::RootScript& script,
                    s13::ana::ElasticScatteringCuts cuts, TString title = "") {
    draw_cut(script, cuts, title, draw_phi_dist, draw_phi_cut_sn);
  }

  void draw_phi_cut_range(s13::ana::RootScript& script,
                          s13::ana::ElasticScatteringCuts cuts, int bins_num) {
    draw_cut_range(script, cuts, bins_num, draw_phi_cut);
  }

  // BG extraction / BG cut

  TH1* make_he4_bg_hist(TChain *ttree, s13::ana::ElasticScatteringCuts cuts,
                        Double_t p_theta_mean) {
    std::cout << "Computing BG for angle: " << p_theta_mean << std::endl;
    TH1* hist_out_cut = s13::ana::
      draw_theta_corr_1d(*ttree, cuts.GetTotalPhiBgCut(),
                         "OUT #phi", "col", gk_default_bin_num);

    std::cout << "[DEBUG:make_he4_bg_hist] BG yield (OUT phi cut): "
              << hist_out_cut->Integral() << std::endl;
    TH1 *hist_bg_scaled = (TH1*)hist_out_cut->Clone();

    auto p_interp = s13::ana::MakeBgNormFactorInterpolator();
    Double_t norm_fact = p_interp->Eval(p_theta_mean) - 1;
    std::cout << "[DEBUG:make_he4_bg_hist] Scaling BG data by: R - 1 = "
              << norm_fact << std::endl;
    hist_bg_scaled->Scale(norm_fact);
    std::cout << "[DEBUG:make_he4_bg_hist] BG yield (IN phi cut): "
              << hist_bg_scaled->Integral() << std::endl;
    hist_bg_scaled->SetTitle(TString::Format("He4 BG @%.2f deg.", p_theta_mean));

    return hist_bg_scaled;
    //return hist_out_cut;
  }

  TH1* draw_he4_bg_cut_sn(TChain *ttree, s13::ana::ElasticScatteringCuts cuts,
                          Double_t p_theta_mean) {
    TH1* hist_in_cut = s13::ana::
      draw_theta_corr_1d(*ttree, cuts.GetTotalCut(),
                         "IN & OUT #phi cut", "col",
                         gk_default_bin_num);

    hist_in_cut->SetTitle(TString::Format(
        "He4 IN #phi @%.2f deg.", p_theta_mean));
    TH1* hist_bg_scaled = make_he4_bg_hist(ttree, cuts, p_theta_mean);
    hist_bg_scaled->SetLineColor(kRed);

    hist_in_cut->Draw();
    hist_bg_scaled->Draw("SAME");

    gPad->BuildLegend(0.15, 0.67, 0.45, 0.88);
    gPad->Update();
    return hist_bg_scaled;
  }

  TH1* draw_he4_sn_wo_bg(TChain *ttree, s13::ana::ElasticScatteringCuts cuts,
                         Double_t p_theta_mean) {
    TH1* hist_signal = s13::ana::
      draw_theta_corr_1d(*ttree, cuts.GetTotalCut(),
                         "He4 IN #phi total", "col", gk_default_bin_num);
    TH1* hist_bg = make_he4_bg_hist(ttree, cuts, p_theta_mean);

    std::cout <<  "[DEBUG:draw_he4_sn_wo_bg] Total IN phi yield: "
              << hist_signal->Integral() << std::endl;
    if (!hist_signal->Add(hist_bg, -1)) {
      throw std::domain_error("hist subtraction failed!");
    }
    std::cout <<  "[DEBUG:draw_he4_sn_wo_bg] IN phi yield - BG: "
              << hist_signal->Integral() << std::endl;

    hist_signal->SetTitle("He4 S/N w/o BG");
    hist_signal->SetMinimum(0);
    hist_signal->SetLineColor(kGreen);
    hist_signal->Draw();

    gPad->Update();

    return hist_signal;
  }

  void draw_he4_bg(s13::ana::RootScript& script,
                   s13::ana::ElasticScatteringCuts cuts,
                   Double_t p_theta_mean) {
    script.cd();
    draw_he4_bg_cut_sn(script.GetChain(), cuts, p_theta_mean);
    script.cd();
    draw_he4_sn_wo_bg(script.GetChain(), cuts, p_theta_mean);
  }

  void draw_he4_bg_range(s13::ana::RootScript& script,
                         s13::ana::ElasticScatteringCuts cuts,
                         int bins_num = 4,
                         Double_t start_angle = 55, Double_t finish_angle = 70) {
    Double_t bin_width = (finish_angle - start_angle) / bins_num;
    TCut old_user_cut = cuts.user_cut;
    cuts.apply_user_cut = true;

    for (int i = 0; i < bins_num; ++i) {
      Double_t theta_min = start_angle + (i * bin_width);
      Double_t theta_max = theta_min + bin_width;
      Double_t theta_mean = start_angle + (bin_width * i) + (bin_width/2);
      cuts.user_cut = old_user_cut && make_proton_theta_eff_cut(theta_min, theta_max);

      draw_he4_bg(script, cuts, theta_mean);
    }

    cuts.user_cut = old_user_cut;
  }
}
}

#endif //ANA_DRAWCUTS_H
