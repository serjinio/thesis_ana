
#include "TChain.h"
#include "TGraph.h"

#include "init_6he_ds.hpp"
#include "macroincludes.hpp"


using namespace s13::ana;


ScatteringYield compute_asymmetry(ElasticScatteringCuts& cuts) {
  // Yields
  int bin_num = 5;
  ScatteringYield left{bin_num, 55, 70};
  ScatteringYield right{bin_num, 55, 70};
  ScatteringYield left_up{bin_num, 55, 70};
  ScatteringYield right_up{bin_num, 55, 70};
  ScatteringYield left_down{bin_num, 55, 70};
  ScatteringYield right_down{bin_num, 55, 70};
  int evt_counter = 0;
  int total_evt_counter = 0;

  auto compute_asym = [&left, &right, &cuts, &evt_counter, &total_evt_counter]
    (ScatteringEvent& evt) -> void {
    if (!cuts(evt)) { // not elastic scattering
      return;
    }
    if (evt.p_theta < 55 || evt.p_theta > 70) {
      return;
    }

    if (evt.IsLeftScattering()) {
      left.AddEvent(evt.p_theta);
    } else {
      right.AddEvent(evt.p_theta);
    }

    ++evt_counter;
    ++total_evt_counter;
  };

  ScatteringTreeWalker tree_walker;
  tree_walker.Walk(g_chain_up, compute_asym);
  std::cout << "Left up: " << left << std::endl;
  std::cout << "Right up: " << right << std::endl;
  std::cout << "Pol. up total events # " << evt_counter << std::endl;
  evt_counter = 0;

  left_up = left; right_up = right;
  left.Reset(); right.Reset();

  tree_walker.Walk(g_chain_down, compute_asym);
  std::cout << "Left down: " << left << std::endl;
  std::cout << "Right down: " << right << std::endl;
  std::cout << "Pol. down total events # " << evt_counter << std::endl;

  left_down = left; right_down = right;
  left.Reset(); right.Reset();

  auto asymmetry_num = sqrt(left_up * right_down) - sqrt(right_up * left_down);
  auto asymmetry_den = sqrt(left_up * right_down) + sqrt(right_up * left_down);
  auto asymmetry = asymmetry_num / asymmetry_den;
  std::cout << "Asymmetry: " << asymmetry << std::endl;
  std::cout << "Total events # " << total_evt_counter << std::endl;

  return asymmetry;
}

void test_ay_4he() {
  init_dataset();
  RootScript script("test_asymmetry");
  ElasticScatteringCuts cuts;
  cuts.phi_width = 8;
  cuts.vertexXY_radius = 9;x
  cuts.theta_width = 2;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_phi_cut = true;
  cuts.apply_theta_cut = false;
  cuts.apply_hodf_he6_cut = true;

  // show S/N
  script.NewPage(2, 1);
  draw_theta_corr_2d(g_chain_total, cuts);
  script.cd(2);
  draw_theta_corr_1d(g_chain_total, cuts);
  draw_marker_line(-cuts.theta_width/2, 0, -cuts.theta_width/2, 650);
  draw_marker_line(cuts.theta_width/2, 0, cuts.theta_width/2, 650);

  // compute asymmetry
  cuts.apply_theta_cut = true;
  auto asym = compute_asymmetry(cuts) * -4;

  script.NewPage(1, 1);
  TString title = "Test asymmetry";
  std::vector<TGraph*> graphs = {Build4HeAyGraph(),
                                 BuildAyGraph(asym, title)};
  auto multi_graph = PrepareAyMultiGraph(title, graphs);

  multi_graph->Draw("ACP");
  //script.GetCanvas().BuildLegend(0.15, 0.67, 0.5, 0.88);
}
