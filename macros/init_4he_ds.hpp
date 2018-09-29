
#pragma once

#include <iostream>
#include <TProof.h>


TChain g_chain_up{"scattree"};
TChain g_chain_down{"scattree"};
TChain g_chain_total{"scattree"};


void init_4he_dataset() {
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
}

TChain* init_4he_partial_dataset(int num_of_runs) {
  TChain *chain = new TChain{"scattree"};

  int he4_start_run = 272;
  for (int i = he4_start_run; i <= he4_start_run + num_of_runs; i++) {
    TString name;
    name.Form("scattree_rootfiles/scattree_run%i.root", i);
    std::cout << "Adding file to chain: " << name << std::endl;
    chain->Add(name);
  }

//  TProof::Open("");
//  chain->SetProof();
  return chain;
}



