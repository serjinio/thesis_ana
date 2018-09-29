//
// Created by serj on 16/10/27.
//

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


RootScript g_script{"hodf_cut_check"};
ElasticScatteringCuts g_cuts;
TChain *g_local_chain = init_6he_partial_dataset(40);


void draw_hodf_cut(int hodf_bar_idx, TString user_cut = "") {
  TString title = TString::Format("%s %s %i (%s)", "6He", "HODF", hodf_bar_idx,
                                  g_cuts.GetTotalCutTitle().Data());
  TString hodf_bar_cut = TString::Format("hodf_q[%i]>7", hodf_bar_idx);
  const TCut total_cut = g_cuts.GetTotalCut() && hodf_bar_cut && user_cut;

  std::cout << "Total cut: " << total_cut << std::endl;

  g_script.cd();
  TH1 *hist_2d = draw_theta_corr_2d(*g_local_chain,
                                    total_cut, title);
  gPad->SetLogz();
  g_script.cd();
  TH1 *hist_1d = draw_theta_corr_1d(*g_local_chain,
                                    total_cut, title);
}


void draw_espri_sn(Espri espri_sel, TString user_cut = "") {
  g_cuts.espri_selector = espri_sel;
  for (int i = 0; i < 24; ++i) {
    if (i % 4 == 0) {
      g_script.NewPage(4,2).cd(PadSequence::column);
    }
    draw_hodf_cut(i, user_cut);
  }
}

void hodf_cut_check() {
  g_cuts.apply_hodf_he6_cut = false;
  g_cuts.apply_phi_cut = true;
  //g_cuts.phi_width = 8;
  g_cuts.apply_vertexXY_cut = true;
//  g_cuts.vertexXY_radius = 5;
  g_cuts.apply_espri_e_cut = false;
  g_cuts.apply_vertexZ_cut = false;

  draw_espri_sn(Espri::left, "s1dc_xpos < 0 && p_e > 15");
  draw_espri_sn(Espri::right, "s1dc_xpos > 0 && p_e > 15");

  g_script.Close();
}