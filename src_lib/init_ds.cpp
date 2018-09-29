
#include "init_ds.hpp"


TChain g_chain_up{"scattree"};
TChain g_chain_down{"scattree"};
TChain g_chain_total{"scattree"};
TChain g_chain_carbon{"scattree"};


TChain* init_partial_dataset(int num_of_runs,
                             int start_run_no) {
  TChain *chain = new TChain{"scattree"};

  for (int i = start_run_no; i < start_run_no + num_of_runs; i++) {
    TString name;
    name.Form("scattree_rootfiles/scattree_run%i.root", i);
    std::cout << "Adding file to chain: " << name << std::endl;
    chain->Add(name);
  }

  return chain;
}

std::vector<TChain*>
init_4he_segmented_dataset(int num_of_segments,
                           int num_of_runs,
                           int start_run) {
  return make_segmented_dataset(num_of_segments,
                                start_run,
                                num_of_runs,
                                init_partial_dataset);
}

std::vector<TChain*>
init_6he_segmented_dataset(int num_of_segments,
                           int num_of_runs) {
  return make_segmented_dataset(num_of_segments,
                                s13::gk_he6_start_run,
                                num_of_runs,
                                init_partial_dataset);
}

TChain* init_4he_dataset() {
  //for (int i = 272; i < 297; i++){ // full
  for (int i = 272; i < 291; i++){ // full for asym
  //for (int i = 272; i < 277; i++){ // sample
    TString name;
    name.Form("scattree_rootfiles/scattree_run%i.root", i);
    std::cout << "Adding file to chain: " << name << std::endl;
    g_chain_up.Add(name);
    g_chain_total.Add(name);
  }
  //for (int i = 297; i < 316; i++){ // full
  for (int i = 297; i < 316; i++) { // full for asym
  //for (int i = 297; i < 302; i++) { // sample
    TString name;
    name.Form("scattree_rootfiles/scattree_run%i.root", i);
    std::cout << "Adding file to chain_down: " << name << std::endl;
    g_chain_down.Add(name);
    g_chain_total.Add(name);
  }

  return &g_chain_total;
}

TChain* init_4he_dataset(int num_of_runs) {
  return init_partial_dataset(num_of_runs, s13::gk_he4_start_run);
}

TChain* init_6he_dataset() {
  //for (int i = 133; i <= 190; i++){ // he6 runs
  for (int i = 133; i <= 269; i++) { // sample
  //for (int i = 133; i <= 150; i++){ // sample2
    TString name;
    name.Form("scattree_rootfiles/scattree_run%i.root", i);
    std::cout << "Adding file to chain_up: " << name << std::endl;
    g_chain_up.Add(name);
    g_chain_total.Add(name);
  }

  //for (int i = 191; i <= 269; i++){ // full dataset
    //for (int i = 194; i <= 252; i++){ // full dataset - skip "bad pol." runs
  for (int i = 194; i <= 193; i++){ // sample
    TString name;
    name.Form("scattree_rootfiles/scattree_run%i.root", i);
    std::cout << "Adding file to chain_down: " << name << std::endl;
    g_chain_down.Add(name);
    g_chain_total.Add(name);
  }

  TProof::Open("");
  g_chain_total.SetProof();
  return &g_chain_total;
}

TChain* init_partial_carbon_dataset(int num_of_runs, int start_run_no) {
  TChain *chain = new TChain{"scattree"};
  for (int i = start_run_no; i < start_run_no + num_of_runs; i++) { // carbon runs
    TString name;
    name.Form("scattree_rootfiles/scattree_run%i.root", i);
    std::cout << "Adding file to carbon chain: " << name << std::endl;
    chain->Add(name);
  }
  return chain;
}

TChain* init_carbon_dataset() {
  //for (int i = 324; i <= 350; i++) { // carbon runs
  for (int i = 324; i <= 350; i++) { // carbon runs
    TString name;
    name.Form("scattree_rootfiles/scattree_run%i.root", i);
    std::cout << "Adding file to carbon chain: " << name << std::endl;
    g_chain_carbon.Add(name);
  }

  return &g_chain_carbon;
}

std::vector<TChain*>
init_carbon_segmented_dataset(int num_of_segments,
                              int num_of_runs) {
  return make_segmented_dataset(num_of_segments,
                                s13::gk_carbon_start_run,
                                num_of_runs,
                                init_partial_dataset);
}

std::vector<TChain*>
init_6he_dataset_from_cli(const s13::misc::CliOptions& cli_opts) {
  int num_of_runs = s13::gk_he6_num_of_runs;
  int start_run = s13::gk_he6_start_run;

  if (cli_opts.pol_direction() == s13::ana::PolarizationDirection::up) {
    start_run = s13::gk_he6_pol_up_start_run;
    num_of_runs = s13::gk_he6_pol_up_finish_run - start_run;
  } else if (cli_opts.pol_direction() == s13::ana::PolarizationDirection::down) {
    start_run = s13::gk_he6_pol_down_start_run;
    num_of_runs = s13::gk_he6_pol_down_finish_run - start_run;
  } else {
    start_run = s13::gk_he6_start_run;
    num_of_runs = s13::gk_he6_finish_run - start_run + 1;
  }

  if (cli_opts.dataset_size() != -1) {
    if (cli_opts.dataset_size() > s13::gk_he6_num_of_runs) {
      throw std::invalid_argument("Too many runs for this dataset!");
    }
    num_of_runs = cli_opts.dataset_size();
  }

  return init_6he_segmented_dataset(cli_opts.threads_num(),
                                    num_of_runs);
}

std::vector<TChain*>
init_4he_dataset_from_cli(const s13::misc::CliOptions& cli_opts) {
  int num_of_runs = s13::gk_he4_num_of_runs;
  int start_run = s13::gk_he4_start_run;

  if (cli_opts.pol_direction() == s13::ana::PolarizationDirection::up) {
    start_run = s13::gk_he4_pol_up_start_run;
    num_of_runs = s13::gk_he4_pol_up_finish_run - start_run;
  } else if (cli_opts.pol_direction() == s13::ana::PolarizationDirection::down) {
    start_run = s13::gk_he4_pol_down_start_run;
    num_of_runs = s13::gk_he4_pol_down_finish_run - start_run;
  } else {
    start_run = s13::gk_he4_start_run;
    num_of_runs = s13::gk_he4_finish_run - start_run + 1;
  }

  if (cli_opts.dataset_size() != -1) {
    if (start_run + cli_opts.dataset_size() > s13::gk_he4_finish_run) {
      throw std::invalid_argument("Too many runs for this dataset!");
    }
    num_of_runs = cli_opts.dataset_size();
  }

  return init_4he_segmented_dataset(cli_opts.threads_num(),
                                    num_of_runs,
                                    start_run);
}

std::vector<TChain*>
init_carbon_dataset_from_cli(const s13::misc::CliOptions& cli_opts) {
  int num_of_runs = s13::gk_carbon_num_of_runs;
  if (cli_opts.dataset_size() != -1) {
    if (cli_opts.dataset_size() > s13::gk_carbon_num_of_runs) {
      throw std::invalid_argument("Too many runs for this dataset!");
    }
    num_of_runs = cli_opts.dataset_size();
  }
  return init_carbon_segmented_dataset(cli_opts.threads_num(),
                                       num_of_runs);
}

std::vector<TChain*>
init_dataset_from_cli(const s13::misc::CliOptions& cli_opts) {
  if (cli_opts.dataset_type() == s13::ana::DatasetType::he4) {
    return init_4he_dataset_from_cli(cli_opts);
  } else if (cli_opts.dataset_type() == s13::ana::DatasetType::he6) {
    return init_6he_dataset_from_cli(cli_opts);
  } else if (cli_opts.dataset_type() == s13::ana::DatasetType::carbon) {
    return init_carbon_dataset_from_cli(cli_opts);
  } else {
    throw std::invalid_argument("Invalid dataset type specified!");
  }
}

TChain* tchain_from_file(std::string filename, std::string tree_name) {
  TChain *chain = new TChain{tree_name.data()};
  chain->Add(filename.data());
  return chain;
}
