/**
 *   \file cs_he.cpp

 *   \brief Computes yields_range for He & carbon
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

std::vector<s13::ana::ElasticScatteringCuts>
init_cuts(std::string cuts_config_fname) {
  std::fstream fs(cuts_config_fname);
  return s13::io::parse_cuts_range_config(fs);
}

std::string
output_fname(s13::ana::PolarizationDirection pol_dir, Cuts& cuts) {
  std::string out_filename = "yields";
  out_filename += "_" + cli_opts.dataset_type_str()
    + "_" + cuts.GetTotalCutName();
  if (pol_dir == s13::ana::PolarizationDirection::up) {
    out_filename += "_pol_up";
  } else if (pol_dir == s13::ana::PolarizationDirection::down) {
    out_filename += "_pol_down";
  }
  return out_filename;
}

void
compute_yields(std::vector<std::shared_ptr<s13::ana::TTreeAlgorithmBase> > algs,
               int bin_num,
               s13::ana::DatasetType ds_type) {
  std::shared_ptr<s13::ana::ScatteringTreeWalker>
    tree_walker(new s13::ana::ScatteringTreeWalker());
  for (auto a : algs) {
    tree_walker->add(a);
  }
  std::vector<TChain*> chains = init_dataset(ds_type);
  s13::ana::ParallelTreeWalker ptree_walker(tree_walker, chains);
  ptree_walker.walk();
}

void save_yields(std::vector<std::shared_ptr<s13::ana::TTreeAlgorithmBase> > algs) {
  s13::misc::MessageLogger log("save_yields()");
  for (auto a : algs) {
    std::string prefix = cli_opts.dataset_type_str() + "_" +
      cli_opts.pol_direction_str();
    a->serialize(prefix);
  }
}

void compute_yields_range(std::vector<Cuts>& cuts,
                          int bin_num,
                          s13::ana::DatasetType ds_type) {
  s13::misc::MessageLogger log("compute_yields_range()");
  std::vector<std::shared_ptr<s13::ana::TTreeAlgorithmBase> > algs;
  for (auto& c : cuts) {
    auto cuts_left = c;
    auto cuts_right = c;
    cuts_left.espri_selector = Espri::left;
    cuts_right.espri_selector = Espri::right;
    auto yields_left_alg = s13::ana::
      make_elastic_scattering_yield_alg(bin_num, s13::gk_cs_range_start, s13::gk_cs_range_end, cuts_left);
    auto yields_right_alg = s13::ana::
      make_elastic_scattering_yield_alg(bin_num, s13::gk_cs_range_start, s13::gk_cs_range_end, cuts_right);
    algs.push_back(yields_left_alg);
    algs.push_back(yields_right_alg);

    auto mon_left_cuts = c;
    auto mon_right_cuts = c;
    mon_left_cuts.espri_selector = Espri::left;
    mon_right_cuts.espri_selector = Espri::right;
    mon_left_cuts.apply_theta_cut = false;
    mon_right_cuts.apply_theta_cut = false;
    auto signal_left_mon_alg = s13::ana::
      make_signal_bg_monitor_alg(bin_num, s13::gk_cs_range_start, s13::gk_cs_range_end, mon_left_cuts);
    auto signal_right_mon_alg = s13::ana::
      make_signal_bg_monitor_alg(bin_num, s13::gk_cs_range_start, s13::gk_cs_range_end, mon_right_cuts);
    algs.push_back(signal_left_mon_alg);
    algs.push_back(signal_right_mon_alg);
  }

  compute_yields(algs, bin_num, ds_type);
  save_yields(algs);
  log.debug("ESL elastic scattering yields alg. total yields: ");
  std::cout << dynamic_cast<s13::ana::ElasticScatteringYieldsAlg*>(algs[0].get())->total_yields();
  log.debug("ELS signal/BG monitor alg. total yields: ");
  std::cout << dynamic_cast<s13::ana::SignalBgMonitorAlg*>(algs[2].get())->total_yields();
}

void yields_range() {
  s13::misc::MessageLogger log("yields_range()");
  std::vector<s13::ana::ElasticScatteringCuts>
    cuts = init_cuts(cli_opts.cuts_config());
  // KLUDGE: currently ignore var theta cut setup
  // if (cli_opts.dataset_type() == DsType::he4) {
  //   for (auto& c : cuts) {
  //     c.PrepareVarThetaCut(cli_opts.dataset_type());
  //   }
  // }

  compute_yields_range(cuts, s13::gk_cs_bin_number, cli_opts.dataset_type());
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

  yields_range();
}
