
#include "TChain.h"
#include "TGraph.h"

#include "macroincludes.hpp"


using namespace s13::ana;


/*
  Function returns a graph of Ay from Moss paper
*/
TGraph* build_graph(std::vector<Double_t>& xs, std::vector<Double_t>& ys) {
  assert(xs.size() == ys.size());
  static const int num_points = xs.size();

  TGraph* graph = new TGraph(num_points, &xs[0], &ys[0]);
  graph->SetLineColor(3);
  graph->SetLineWidth(3);
  graph->SetMarkerColor(kRed);
  graph->SetMarkerStyle(22);
  graph->SetMarkerSize(3);
  graph->SetTitle("Proton deflection angle vs. B0");
  graph->GetXaxis()->SetTitle("B0 [T]");
  graph->GetYaxis()->SetTitle("Angle [deg.]");
  // graph->SetMinimum(-1.2);
  // graph->SetMaximum(1.2);

  return graph;
};

Double_t compute_deflection_angle(Double_t B0_field, Double_t proton_E_kin) {
  constexpr Double_t proton_charge = 1.6e-19;
  constexpr Double_t proton_mass = 1.67e-27;
  constexpr Double_t proton_mass_mev = 938;
  constexpr Double_t magnet_pole_radius = 10e-2;
  constexpr Double_t c = 3.0e8;
  constexpr Double_t D2R = 1./57.3;

  Double_t k = proton_charge * magnet_pole_radius / proton_mass;
  cout << "k: " << k << std::endl;
  Double_t gamma = (proton_E_kin / proton_mass_mev) + 1;
  cout << "gamma: " << gamma << std::endl;
  Double_t proton_velocity = sqrt(1 - (1/pow(gamma, 2))) * c;
  cout << "proton velocity: " << proton_velocity << std::endl;

  return atan(k * B0_field / proton_velocity) * 1/D2R;
}

void proton_deflection_angle_in_mag_field() {
  RootScript script("proton_deflection_angle_in_mag_field");
  gROOT->SetStyle("Video");

  script.NewPage(1);
  std::vector<Double_t> xs;
  std::vector<Double_t> ys_40mev;
  std::vector<Double_t> ys_80mev;
  for (Double_t b0 = 0.05; b0 < 0.35; b0 += 0.05) {
    xs.push_back(b0);
    Double_t def_angle = compute_deflection_angle(b0, 40);
    ys_40mev.push_back(def_angle);
    ys_80mev.push_back(compute_deflection_angle(b0, 120));
    std::cout << "Deflection angle for B0=" << b0 << ": " << def_angle << std::endl;
  }
  auto graph = build_graph(xs, ys_40mev);
  auto graph2 = build_graph(xs, ys_80mev);
  graph->Draw("ACP");
  graph2->SetMarkerColor(kGreen);
  graph2->Draw("CP SAME");
  draw_marker_line(0.03,1,0.32,1,3);
}


