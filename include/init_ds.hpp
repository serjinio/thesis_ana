
#pragma once

#include <iostream>
#include <TChain.h>
#include <TProof.h>

#include "consts.hpp"
#include "common.hpp"
#include "msgoutput.hpp"
#include "cli.hpp"


extern TChain g_chain_up;
extern TChain g_chain_down;
extern TChain g_chain_total;
extern TChain g_chain_carbon;


/*
  Helper function to make dataset from several segments.
*/
template<typename T>
std::vector<TChain*>
make_segmented_dataset(int num_of_segments,
                       int first_run,
                       int num_of_runs,
                       T segment_ctor_fn) {
  s13::misc::MessageLogger log("make_segmented_dataset()");
  log.debug("Will make %d segments with total runs number: %d...",
            num_of_segments, num_of_runs);
  std::vector<TChain*> segs;
  const int runs_per_segment = num_of_runs / num_of_segments;
  if (runs_per_segment == 0) {
    throw std::invalid_argument("Too many segments requested!");
  }

  const int left_over = num_of_runs - (runs_per_segment * num_of_segments);
  int segm_no = 1, consumed_left_overs = 0;
  while (segm_no <= num_of_segments) {
    int seg_start_run = first_run + (runs_per_segment * (segm_no - 1))
      + consumed_left_overs;
    int seg_runs_num = runs_per_segment - 1;
    if (consumed_left_overs < left_over) {
      seg_runs_num += 1;
      consumed_left_overs += 1;
    }

    segs.push_back(segment_ctor_fn(seg_runs_num + 1, seg_start_run));
    log.debug("Made segment #%d with %d runs.", segm_no, seg_runs_num + 1);
    segm_no += 1;
  }
  log.debug("%d segments were made.", num_of_segments);

  return segs;
}


////////////////////////////////////////////////////////////
// Dataset construction functions
////////////////////////////////////////////////////////////


/*
  Makes partial dataset with specified number of runs included.
*/
TChain* init_partial_dataset(int num_of_runs,
                             int start_run_no = s13::gk_he6_start_run);

/*
  Makes specified number of segments of 4he dataset.
  For use in multi-threaded computations.
*/
std::vector<TChain*>
init_4he_segmented_dataset(int num_of_segments,
                           int num_of_runs = s13::gk_he4_num_of_runs,
                           int first_run = s13::gk_he4_start_run);

/*
  Makes specified number of segments of 6he dataset.
  For use in multithreaded computations.
*/
std::vector<TChain*>
init_6he_segmented_dataset(int num_of_segments,
                           int num_of_runs = s13::gk_he6_num_of_runs);
/*
  Makes full 4he dataset
*/
TChain* init_4he_dataset();

/*
  Makes 4he dataset with specified number of runs
*/
TChain* init_4he_dataset(int num_of_runs);

/*
  Makes full 6he dataset
*/
TChain* init_6he_dataset();


TChain* init_partial_carbon_dataset(int num_of_runs, int start_run_no = 324);
TChain* init_carbon_dataset();

/*
  Makes specified number of segments of carbon dataset.
  For use in multithreaded computations.
*/
std::vector<TChain*>
init_carbon_segmented_dataset(int num_of_segments,
                              int num_of_runs = s13::gk_carbon_num_of_runs);

/*
  Following functions initialize a dataset basing on passed command line options.
*/
std::vector<TChain*>
init_4he_dataset_from_cli(const s13::misc::CliOptions& cli_opts);
std::vector<TChain*>
init_6he_dataset_from_cli(const s13::misc::CliOptions& cli_opts);
std::vector<TChain*>
init_carbon_dataset_from_cli(const s13::misc::CliOptions& cli_opts);

std::vector<TChain*>
init_dataset_from_cli(const s13::misc::CliOptions& cli_opts);

TChain*
tchain_from_file(std::string filename, std::string tree_name);
