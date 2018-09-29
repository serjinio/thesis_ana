
#include "TChain.h"
#include "TGraph.h"

#include "init_4he_ds.hpp"
#include "macroincludes.hpp"


using namespace s13::ana;


TH1* plot_tof_de(ElasticScatteringCuts& cuts) {
  static int hist_num = 1;
  TString title = TString::Format("TOF-dE (%s)",
                                  cuts.GetTotalCutTitle().Data());
  TString hist_name = TString::Format("tof_de%i", hist_num);
  TString draw_cmd = TString::Format("p_de:p_tof");
  //TCut user_cut = (TString::Format("p_de > 1 && p_e_eff > 5")).Data();

  TH1* hist = hdraw(g_chain_total, hist_name, draw_cmd, "(100,7000,8800,100,0,14)",
                    cuts.GetTotalCut(), title, "TOF (a.u.)", "dE (MeV)", "colz");
  gPad->SetLogz();

  ++hist_num;
  return hist;
}

void he4_add_tof_de_espri_cut() {
  init_dataset();
  RootScript script("he4_add_espri_cuts");
  ElasticScatteringCuts cuts;

  cuts.apply_vertexXY_cut = false;
  cuts.apply_phi_cut = false;
  cuts.apply_theta_cut = false;
  cuts.apply_espri_e_cut = false;

  script.NewPage(2,2);

  script.cd(1);
  cuts.espri_selector = Espri::left;
  plot_tof_de(cuts);
  script.cd(2);
  cuts.espri_selector = Espri::right;
  plot_tof_de(cuts);

  cuts.apply_vertexXY_cut = true;
  cuts.apply_phi_cut = true;
  cuts.apply_theta_cut = true;
  cuts.apply_espri_e_cut = true;
  cuts.espri_e_width = 10;

  script.cd(3);
  cuts.espri_selector = Espri::left;
  plot_tof_de(cuts);
  script.cd(4);
  cuts.espri_selector = Espri::right;
  plot_tof_de(cuts);


}
