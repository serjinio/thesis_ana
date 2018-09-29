
#include "TChain.h"
#include "TGraph.h"

#include "init_4he_ds.hpp"
#include "macroincludes.hpp"
#include "drawcuts.h"

using namespace s13::ana;
using namespace s13::ana::algs;


TChain *g_local_chain = init_4he_partial_dataset(15);
RootScript* script = new RootScript("rdc_pos_dep_eff_2", g_local_chain);

constexpr double proton_min_angle = 55;
constexpr double proton_max_angle = 70;
constexpr double he4_min_angle = 4;
constexpr double he4_max_angle = 12;


void draw_eff_text(TH1* hist) {
  gPad->Update();
  double eff = hist->Integral() / hist->GetEntries();
  auto eff_str = TString:: Format("Eff.: %.2f %%", eff * 100);
  auto eff_text =
    new TText(47, gPad->GetUymax() * 0.95,
              eff_str);
  eff_text->SetTextSize(0.08);
  eff_text->Draw();
}

void draw_cut_text(TH1* hist, TCut cut) {
  //auto eff_str = TString:: Format("Eff.: %.2f %%", eff * 100);
  auto text = new TText(47, gPad->GetUymax() * 0.15, cut);
  text->SetTextSize(0.02);
  text->Draw();
}

TH1* draw_rdc_eff(TCut cuts, TString title) {
  std::cout << "[DEBUG:draw_rdc_eff] total cut: " << cuts << std::endl;
  static int hist_num = 0;
  TString hist_name = TString::Format("rdc_eff%i", hist_num);
  // TString draw_cmd = cuts.espri_selector == Espri::left ?
  //   "esl_xpos" : "esr_xpos";
  TString draw_cmd = "p_theta*57.3";
  TH1* hist =
    hdraw(*g_local_chain, hist_name, draw_cmd,
          "(300, 45, 80)", cuts,
          title, "#theta [lab. deg.]", "Counts");

  ++hist_num;
  return hist;
}

void draw_rdc_eff_range(TCut cuts,
                        Espri espri_sel,
                        int bins_num) {
  assert(espri_sel != Espri::both);

  Double_t bin_width = (he4_max_angle - he4_min_angle) / bins_num;
  TString prefix = espri_sel  == Espri::left ? "esl" : "esr";

  for (int i = 0; i < bins_num; i++) {
    if (i % 3 == 0) {
      script->NewPage(3,1);
    }
    Double_t bin_min_angle = he4_min_angle + (i * bin_width);
    Double_t bin_max_angle = he4_min_angle + ((i + 1) * bin_width);
    TCut he_theta_cut =
      TString::Format("s1dc_theta*57.3 > %.2f && s1dc_theta*57.3 < %.2f",
                      bin_min_angle, bin_max_angle).Data();
    TCut he_phi_cut = "abs(s1dc_ypos) < 60";
    TCut proton_de_cut = TString::
      Format("%s_de > 230", prefix.Data()).Data();
    TCut proton_e_cut = TString::
      Format("p_e_eff > 0.45").Data();

    TCut user_cut = he_theta_cut && he_phi_cut
      && proton_de_cut && proton_e_cut;
    TString title = TString::
      Format("%s RDC X dist (%.1f<#theta<%.1f)",
             espri_sel == Espri::left ? "ESL" : "ESR",
             bin_min_angle, bin_max_angle);
    script->cd();
    auto hist = draw_rdc_eff(cuts && user_cut, title);
    draw_eff_text(hist);
    //draw_cut_text(hist, cuts && user_cut);
  }
}

void rdc_pos_dep_eff() {
  gStyle->SetOptFit();
  ElasticScatteringCuts cuts;
  cuts.vertexXY_radius = 6;

  draw_rdc_eff_range(cuts.GetVertexXYCut(), Espri::left, 6);
  draw_rdc_eff_range(cuts.GetVertexXYCut(), Espri::right, 6);

  delete script;
}
