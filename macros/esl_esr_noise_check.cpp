//
// Created by serj on 16/10/25.
//

#include <fstream>
#include <future>

#include "TChain.h"
#include "TGraph.h"
#include "Math/Interpolator.h"

#include "init_6he_ds.hpp"
#include "macroincludes.hpp"


using namespace s13::ana;


RootScript g_script{"esl_esr_noise_check"};
ElasticScatteringCuts g_cuts;
TChain *g_local_chain = init_6he_partial_dataset(40);


void draw_sn() {
  TString title = TString::Format("%s (%s)", "6He",
                                  g_cuts.GetTotalCutTitle().Data());

  g_script.cd();
  TH1 *hist_2d = draw_theta_corr_2d(*g_local_chain, g_cuts.GetTotalCut(), title);
  g_script.cd();
  TH1 *hist_1d = draw_theta_corr_1d(*g_local_chain, g_cuts.GetTotalCut(), title);
}

void esl_esr_noise_check() {
  g_cuts.apply_hodf_he6_cut = true;
  g_cuts.apply_phi_cut = true;
  //g_cuts.phi_width = 8;
  g_cuts.apply_vertexXY_cut = true;
  //g_cuts.vertexXY_radius = 5;
  g_cuts.apply_espri_e_cut = false;
  g_cuts.apply_vertexZ_cut = false;

  g_script.NewPage(2,2).cd(PadSeq::column);
  g_cuts.espri_selector = Espri::left;
  draw_sn();
  g_cuts.espri_selector = Espri::right;
  draw_sn();

  g_script.Close();
}
