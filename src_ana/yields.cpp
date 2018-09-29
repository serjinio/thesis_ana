/**
 *   \file cs_he.cpp

 *   \brief Computes yields for He & carbon
 *
 */

#include "cli.hpp"
#include "init_ds.hpp"
#include "consts.hpp"
#include "csalg.hpp"
#include "drawing.hpp"
#include "scattyield.hpp"
#include "cuts_conf.hpp"
#include "rootscript.hpp"


using Cuts = s13::ana::ElasticScatteringCuts;
using DsType = s13::ana::DatasetType;
using Espri = s13::ana::Espri;


static s13::misc::CliOptions cli_opts;


std::vector<TChain*> init_dataset(s13::ana::DatasetType ds_type) {
  return init_dataset_from_cli(cli_opts);
}

s13::ana::ElasticScatteringCuts
init_cuts(std::string cuts_config_fname) {
  std::fstream fs(cuts_config_fname);
  return s13::io::parse_cuts_config(fs);
}

std::string
output_fname(s13::ana::Espri espri_sel, s13::ana::PolarizationDirection pol_dir) {
  std::string out_filename = "yields";
  out_filename += "_" + cli_opts.dataset_type_str()
    + "_" + s13::ana::espri_selector_to_str(espri_sel);
  if (pol_dir == s13::ana::PolarizationDirection::up) {
    out_filename += "_pol_up";
  } else if (pol_dir == s13::ana::PolarizationDirection::down) {
    out_filename += "_pol_down";
  }
  return out_filename;
}

void
compute_yields(s13::ana::ElasticScatteringCuts& cuts,
               int bin_num,
               s13::ana::DatasetType ds_type) {
  auto cuts_left = cuts;
  auto cuts_right = cuts;
  cuts_left.espri_selector = Espri::left;
  cuts_right.espri_selector = Espri::right;
  auto yields_left_alg = s13::ana::
    make_elastic_scattering_yield_alg(bin_num, 55, 70, cuts_left);
  auto yields_right_alg = s13::ana::
    make_elastic_scattering_yield_alg(bin_num, 55, 70, cuts_right);

  auto mon_left_cuts = cuts;
  auto mon_right_cuts = cuts;
  mon_left_cuts.espri_selector = Espri::left;
  mon_right_cuts.espri_selector = Espri::right;
  mon_left_cuts.apply_theta_cut = false;
  mon_right_cuts.apply_theta_cut = false;
  auto signal_left_mon_alg = s13::ana::
    make_signal_bg_monitor_alg(bin_num, 55, 70, mon_left_cuts);
  auto signal_right_mon_alg = s13::ana::
    make_signal_bg_monitor_alg(bin_num, 55, 70, mon_right_cuts);

  std::shared_ptr<s13::ana::ScatteringTreeWalker>
    tree_walker(new s13::ana::ScatteringTreeWalker());
  tree_walker->add(yields_left_alg);
  tree_walker->add(yields_right_alg);
  tree_walker->add(signal_left_mon_alg);
  tree_walker->add(signal_right_mon_alg);
  std::vector<TChain*> chains = init_dataset(ds_type);
  s13::ana::ParallelTreeWalker ptree_walker(tree_walker, chains);
  ptree_walker.walk();

  std::cout << "Left yields: " << yields_left_alg->total_yields();
  std::cout << "Right yields: " << yields_right_alg->total_yields();
  auto outl_fname = output_fname(Espri::left, cli_opts.pol_direction());
  auto outr_fname = output_fname(Espri::right, cli_opts.pol_direction());
  std::ofstream os_left(s13::ana::tstrfmt("cs_data/%s.csv", outl_fname.c_str()));
  std::ofstream os_right(s13::ana::tstrfmt("cs_data/%s.csv", outr_fname.c_str()));
  yields_left_alg->serialize(os_left);
  signal_left_mon_alg->serialize(os_left);
  yields_right_alg->serialize(os_right);
  signal_right_mon_alg->serialize(os_right);
}

void yields() {
  s13::misc::MessageLogger log("yields()");
  s13::ana::ElasticScatteringCuts
    cuts = init_cuts(cli_opts.cuts_config());
  if (cli_opts.dataset_type() == DsType::he4) {
    cuts.PrepareVarThetaCut(cli_opts.dataset_type());
  }

  compute_yields(cuts, s13::gk_cs_bin_number, cli_opts.dataset_type());
  log.info("Yields computation finished. Exiting.");
}

int main(int argc, char** argv) {
  s13::misc::MessageLogger log("main()");
  TThread::Initialize();
  try {
    cli_opts.parse_options(argc, argv);
  } catch (std::exception& ex) {
    log.error("Invalid argument provided: %s", ex.what());
    return -1;
  }

  yields();
}
