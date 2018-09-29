//
// Created by serj on 16/10/20.
//


#include <TROOT.h>
#include <TStyle.h>

#include "init_4he_ds.hpp"
#include "macroincludes.hpp"


using namespace s13::ana;
using namespace s13::ana::algs;

TH1* draw_theta_diff(ElasticScatteringCuts& cuts) {
  static int hist_num = 1;
  TString title = TString::Format("p_theta - p_theta_eff (%s)",
                                  cuts.GetTotalCutTitle().Data());
  TString hist_name = TString::Format("theta_diff_%i", hist_num);
  TString draw_cmd = TString::Format("p_theta*57.3 - p_theta_eff*57.3");
  TString user_cut = "p_theta_eff*57.3 >= 55 && p_theta_eff*57.3 <= 70";
  TH1* hist = hdraw(g_chain_total, hist_name, draw_cmd, "(100,-2,2)",
                    cuts.GetTotalCut() && user_cut, title,
                    "p_theta - p_theta_eff (deg.)", "Yield");

  ++hist_num;
  return hist;
}

void p_theta_vs_p_theta_eff() {
  init_dataset();
  RootScript script("p_theta_vs_p_theta_eff");
  gStyle->SetOptStat(1111100);
  ElasticScatteringCuts cuts;
  cuts.espri_selector = Espri::both;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_phi_cut = true;
  cuts.phi_width = 16;
  cuts.apply_espri_e_cut = false;
  cuts.espri_e_width = 30;
  cuts.apply_theta_cut = true;

  script.NewPage(1, 1);
  auto hist = draw_theta_diff(cuts);
  //hist->Draw("AC");
}