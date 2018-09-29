//
// Created by serj on 16/10/25.
//

#include <fstream>
#include <future>

#include <TROOT.h>
#include <TRandom.h>
#include "TChain.h"
#include "TGraph.h"

#include "init_6he_ds.hpp"
#include "macroincludes.hpp"

//#include "ProofTest.h"


using namespace s13::ana;

TH1 *compute_hist(TString name, TString title) {
  TH1 *hist = new TH1F(name, title, 100, 55,70);
  TTree *ttree = build_6he_total_dataset();
  TRandom random;
  ElasticScatteringCuts cuts;

//  for (int i = 0; i < 10000; ++i) {
//    hist->Fill(random.Gaus(500,100));
//  }
  //draw_theta_corr_2d(g_chain_total, cuts.GetTotalCut());
  ScatteringTreeWalker treeWalker;
  treeWalker.Walk(*ttree, [&] (ScatteringEvent& evt) {
    if (cuts(evt)) {
      std::cout << "got elastic scattering event" << std::endl;
      hist->Fill(evt.p_theta_eff);
    }
  });

  return hist;
}


void async_test() {
  //init_dataset();
  RootScript script("async_test");
  ElasticScatteringCuts cuts;

  std::vector<std::future<TH1*> > futs;
  for (int i = 0; i < 2; ++i) {
    TString name = TString::Format("gaus_%i", i);
    TString title = TString::Format("Gaussian %i", i);
    futs.push_back(std::async(std::launch::async, compute_hist, name, title));
  }

  script.NewPage(3,2);
  for (int j = 0; j < 2; ++j) {
    futs.at(j).wait();
    TH1 *hist = futs[j].get();
    script.cd(j+1);
    hist->Draw();
  }

}

int main() {
  async_test();
}