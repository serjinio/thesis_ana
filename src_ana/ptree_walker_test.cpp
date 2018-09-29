/**
 *   \file ptree_walker_test.cpp
 *   \brief Test of ParallelTreeWalker
 *
 */

#include <TChain.h>
#include <TGraph.h>
#include <TMultiGraph.h>

#include "init_ds.hpp"
#include "consts.hpp"
#include "csalg.hpp"
#include "drawing.hpp"


s13::ana::ScatteringYield
compute_yields(s13::ana::ElasticScatteringCuts& cuts,
               int bin_num) {
  auto yields = s13::ana::
    make_elastic_scattering_yield_alg(bin_num, 55, 70, cuts);
  auto chains = init_6he_segmented_dataset(4);
  auto tree_walker =
    std::shared_ptr<s13::ana::ScatteringTreeWalker>(new s13::ana::ScatteringTreeWalker());

  tree_walker->add(yields);
  s13::ana::ParallelTreeWalker ptree_walker(tree_walker, chains);
  ptree_walker.walk();

  std::cout << "Computation finished total yields collected:" << std::endl;
  std::cout << yields->total_yields();

  return yields->total_yields();
}

s13::ana::ElasticScatteringCuts init_he6_cuts(s13::ana::Espri espri_selector) {
  s13::ana::ElasticScatteringCuts cuts;

  cuts.espri_selector = espri_selector;
  cuts.apply_vertexXY_cut = true;
  cuts.vertexXY_radius = 6;

  cuts.apply_phi_cut = true;
  cuts.phi_width = 37.8;

  cuts.apply_theta_cut = true;
  cuts.theta_width = 3.2;

  cuts.apply_hodf_he6_cut = true;

  cuts.apply_espri_aang_cut = false;
  cuts.espri_aang_width = 180;

  cuts.apply_espri_min_e_de_cut = false;
  cuts.apply_espri_e_cut = false;

  cuts.apply_espri_y_cut = true;
  cuts.espri_y_width = s13::gk_cs_espri_y_acceptance;

  return cuts;
}

void ptree_walker_test() {
  std::string out_filename = "ptree_walker_test";
  // s13::ana::RootScript script(out_filename, chain);
  s13::ana::ElasticScatteringCuts cuts = init_he6_cuts(s13::ana::Espri::both);

  compute_yields(cuts, 15);
}

int main(int argc, char** argv) {
  ptree_walker_test();
}
