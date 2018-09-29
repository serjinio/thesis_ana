
#pragma once

#include <iostream>
#include <TProof.h>


TChain g_chain_up{"scattree"};
TChain g_chain_down{"scattree"};
TChain g_chain_total{"scattree"};
// carbon runs dataset
TChain g_c_chain{"scattree"};


void init_6he_dataset() {
  //TProof::Open("");

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

  //g_chain_total.SetProof();
}

TChain* init_6he_partial_dataset(int num_of_runs, int start_run_no = 133) {
  TChain *chain = new TChain{"scattree"};

  for (int i = start_run_no; i < start_run_no + num_of_runs; i++) {
    TString name;
    name.Form("scattree_rootfiles/scattree_run%i.root", i);
    std::cout << "Adding file to chain: " << name << std::endl;
    chain->Add(name);
  }

  // TProof::Open("");
  // chain->SetProof();
  return chain;
}


void init_carbon_dataset() {
  //for (int i = 324; i <= 350; i++) { // carbon runs
  for (int i = 324; i <= 350; i++) { // carbon runs
    TString name;
    name.Form("scattree_rootfiles/scattree_run%i.root", i);
    std::cout << "Adding file to carbon chain: " << name << std::endl;
    g_c_chain.Add(name);
  }
}

