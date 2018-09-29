
#include <fstream>

#include "TChain.h"
#include "TGraph.h"

#include "init_4he_ds.hpp"
#include "macroincludes.hpp"


using namespace s13::ana;
using namespace s13::ana::algs;


auto g_moss_ay = LoadHe4MossAy("data/he4_200mev_moss_cs_ay.csv");

// phi cut defined as 3 * sigma_phi_corr
const Double_t gk_phi_cut_width = 16;


void show_sn(RootScript &script, ElasticScatteringCuts cuts);
void draw_asymmetry(RootScript &script, const ElasticScatteringCuts &cuts);

Double_t ay_fit_fn(Double_t *v, Double_t *par) {
//   std::cout << "Evaluating Moss Ay at: " << v[0] <<
//     " with param value: " << par[0] << " fn value: " <<
//     g_moss_ay->Eval(v[0] * D2R) << std::endl;
  return g_moss_ay->Eval(v[0] * D2R) * par[0];
}

Double_t fit_moss_to_asymmetry_hist(TGraph *asymmetry_graph) {
  //create a function with 1 parameter in the range [0, 1e6]
  TF1 *func = new TF1("ay_fit", ay_fit_fn, 55, 70, 1);
  func->SetParLimits(0, 1e-3, 1);
  func->SetParNames("Target_pol");
  func->SetParameter("Target_pol", 0.9);
  asymmetry_graph->Fit("ay_fit", "R");
  const Double_t *fit_params_raw = asymmetry_graph->GetFunction("ay_fit")->GetParameters();
  const Double_t *fit_errors_raw = asymmetry_graph->GetFunction("ay_fit")->GetParErrors();
  std::cout << "Asymmetry fit parameter (target polarization): " << fit_params_raw[0] << std::endl;
  return fit_params_raw[0];
}

ScatteringYield compute_asymmetry(ElasticScatteringCuts& cuts) {
  // Yields
  int bin_num = 5;
  ScatteringYield left_total{bin_num, 55, 70};
  ScatteringYield right_total{bin_num, 55, 70};
  ScatteringYield left_bg_out{bin_num, 55, 70};
  ScatteringYield right_bg_out{bin_num, 55,  70};
  ScatteringYield left_up{bin_num, 55, 70};
  ScatteringYield right_up{bin_num, 55, 70};
  ScatteringYield left_down{bin_num, 55, 70};
  ScatteringYield right_down{bin_num, 55, 70};
  int evt_counter = 0;
  int total_evt_counter = 0;

  auto compute_asym = [&left_total, &right_total, &left_bg_out, &right_bg_out,
                       &cuts, &evt_counter, &total_evt_counter]
    (ScatteringEvent& evt) -> void {

    // common cuts - trigger condition && proton angular acceptance
    if (!cuts.IsTriggerCut(evt) || !cuts.IsProtonInAcceptanceCut(evt)) {
      return;
    }
    // common cuts - beam on target cut
    if (!cuts.IsVertexXYCut(evt)) {
      return;
    }
    if (!cuts.IsEspriECut(evt)) {
      return;
    }
    if (!cuts.IsThetaCut(evt)) {
      return;
    }

    // BG counts - out of phi cut (to subtract BG)
    if (cuts.IsPhiBgCut(evt)) {
      if (evt.IsLeftScattering()) {
        left_bg_out.AddEvent(evt.p_theta);
      } else {
        right_bg_out.AddEvent(evt.p_theta);
      }
    }

    // elastic scattering counts
    if (cuts.IsPhiCut(evt)) {
      if (evt.IsLeftScattering()) {
        left_total.AddEvent(evt.p_theta);
      } else {
        right_total.AddEvent(evt.p_theta);
      }

      ++evt_counter;
      ++total_evt_counter;
    }
  };

  auto bg_norm_fact_interp = MakeBgNormFactorInterpolator();
  ScatteringYield left(bin_num, 55, 70), right(bin_num, 55, 70);
  ScatteringTreeWalker tree_walker;

  tree_walker.Walk(g_chain_up, compute_asym);
  left = He4SubtractBg(left_total, left_bg_out, *bg_norm_fact_interp);
  right = He4SubtractBg(right_total, right_bg_out, *bg_norm_fact_interp);

  std::cout << "Left up total: " << left_total << std::endl;
  std::cout << "Right up total: " << right_total << std::endl;
  std::cout << "Left up: " << left << std::endl;
  std::cout << "Right up: " << right << std::endl;
  std::cout << "Pol. up total events # " << evt_counter << std::endl;

  // reset event counter
  evt_counter = 0;

  left_up = left; right_up = right;
  left_total.Reset(); right_total.Reset();
  left.Reset(); right.Reset();
  left_bg_out.Reset(); right_bg_out.Reset();

  tree_walker.Walk(g_chain_down, compute_asym);
  left = He4SubtractBg(left_total, left_bg_out, *bg_norm_fact_interp);
  right = He4SubtractBg(right_total, right_bg_out, *bg_norm_fact_interp);

  std::cout << "Left down total: " << left_total << std::endl;
  std::cout << "Right down total: " << right_total << std::endl;
  std::cout << "Left down: " << left << std::endl;
  std::cout << "Right down: " << right << std::endl;
  std::cout << "Pol. down total events # " << evt_counter << std::endl;

  left_down = left; right_down = right;
  left_total.Reset(); right_total.Reset();
  left.Reset(); right.Reset();
  left_bg_out.Reset(); right_bg_out.Reset();

  auto asymmetry_num = (left_up * right_down).sqrt()
    - (right_up * left_down).sqrt();
  auto asymmetry_den = (left_up * right_down).sqrt()
    + (right_up * left_down).sqrt();
  auto asymmetry = asymmetry_num / asymmetry_den;
  std::cout << "Asymmetry: " << asymmetry << std::endl;
  std::cout << "Total events # " << total_evt_counter << std::endl;

  return asymmetry;
}

void asymmetry_he4_wo_bg() {
  init_dataset();

  RootScript script("asymmetry_he4_wo_bg");
  gStyle->SetOptFit();
  gStyle->SetOptStat(1111100);
  ElasticScatteringCuts cuts;
  cuts.theta_width = 4;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_theta_cut = true;

  cuts.apply_phi_cut = true;
  cuts.apply_espri_min_e_de_cut = true;
  cuts.apply_espri_e_cut = true;
  cuts.apply_espri_aang_cut = true;

  cuts.phi_width = 8;
  cuts.espri_aang_width = 30;
  cuts.espri_e_width = 0.4;
  show_sn(script, cuts);
  draw_asymmetry(script, cuts);


//  for (Double_t radius = 12; radius >= 2; radius -= 2) {
//    cuts.vertexXY_radius = radius;
//    show_sn(script, cuts);
//    draw_asymmetry(script, cuts);
//  }
//  cuts.vertexXY_radius = 10;
//
//  for (Double_t phi = 32; phi >= 4; phi -= 4) {
//    cuts.phi_width = phi;
//    show_sn(script, cuts);
//    draw_asymmetry(script, cuts);
//  }
//  cuts.phi_width = 16;
//
//  for (Double_t e_width = 0.9; e_width >= 0.1; e_width -= 0.2) {
//    cuts.espri_e_width = e_width;
//    show_sn(script, cuts);
//    draw_asymmetry(script, cuts);
//  }
//  cuts.espri_e_width = 0.7;
//
//  for (Double_t ab_width = 120; ab_width >= 30; ab_width -= 30) {
//    cuts.espri_aang_width = ab_width;
//    show_sn(script, cuts);
//    draw_asymmetry(script, cuts);
//  }
//  cuts.espri_aang_width = 90;

}

void draw_asymmetry(RootScript &script, const ElasticScatteringCuts &cuts) {
  // apparently pol. up & down definition is actually inversed
  // due to this exchange chain up with chain down
  auto asym = He4ScatteringAsymmetry(g_chain_down, g_chain_up, cuts, 3);
  std::cout << "\n\n\ncomputed asymmetry: " << asym << std::endl;
  ofstream asymmetry_he4_ofs;
  asymmetry_he4_ofs.open("data/he4_asym.csv");
  asym.Serialize(asymmetry_he4_ofs);
  asymmetry_he4_ofs.close();

  script.NewPage(1, 1);
  script.cd();
  auto asymmetry_graph =
    BuildAsymmetryGraphErrors(asym, cuts.GetTotalCutTitle());
  asymmetry_graph->Draw("ACP");
  fit_moss_to_asymmetry_hist(asymmetry_graph);
}

void show_sn(RootScript &script, ElasticScatteringCuts cuts) {
  cuts.apply_theta_cut = false;

  script.NewPage(2, 1).cd(PadSeq::row);
  script.cd();
  draw_theta_corr_2d(g_chain_total, cuts.GetTotalCut(), cuts.GetTotalCutTitle());
  script.cd();
  draw_theta_corr_1d(g_chain_total, cuts.GetTotalCut(), cuts.GetTotalCutTitle());
  gPad->Update();

  draw_marker_line(-cuts.theta_width/2, 0, -cuts.theta_width/2, gPad->GetUymax());
  draw_marker_line(cuts.theta_width/2, 0, cuts.theta_width/2, gPad->GetUymax());
}


int main() {
  asymmetry_he4_wo_bg();
}
