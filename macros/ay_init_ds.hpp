
#pragma once

#include <iostream>


TChain g_chain_up{"scattree"};
TChain g_chain_down{"scattree"};


void init_dataset() {
  //for (int i = 272; i < 297; i++){
  for (int i = 272; i < 282; i++){
    TString name;
    name.Form("scattree_rootfiles/scattree_run%i.root", i);
    std::cout << "Adding file to chain_up: " << name << std::endl;
    g_chain_up.Add(name);
  }
  // for (int i = 297; i < 316; i++){
  for (int i = 297; i < 307; i++){
    TString name;
    name.Form("scattree_rootfiles/scattree_run%i.root", i);
    std::cout << "Adding file to chain_down: " << name << std::endl;
    g_chain_down.Add(name);
  }
}
