
#include "TChain.h"
#include "TGraph.h"

#include "init_4he_ds.hpp"
#include "macroincludes.hpp"
#include "drawcuts.h"


using namespace s13::ana;
using namespace s13::cuts;


TChain* g_local_chain = init_4he_partial_dataset(5);
TChain& lch = *g_local_chain;
RootScript g_script("tmp_mac", g_local_chain);
RootScript& sc = g_script;
ElasticScatteringCuts g_cuts;
ElasticScatteringCuts& cuts = g_cuts;


void tmp_mac() {
  cuts.apply_vertexXY_cut = true;
  cuts.apply_phi_cut = true;
  cuts.espri_selector = Espri::left;

  sc.NewPage(2, 1);
  sc.cd();
  lch.Draw("esl_ypos >>esl_ypos1(100, -300,300)",
                       cuts.GetTotalCut());
  sc.cd();
  lch.Draw("esl_ypos >>esl_ypos2(100, -300,300)",
           cuts.GetTotalCut() && cuts.GetEspriMinEdECut());
}

