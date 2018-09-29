#include <iostream>

#include "init_6he_ds.hpp"
#include "TChain.h"


void hdraw(TTree& tree, TString name, TString draw_cmd,
           TString binning, TCut cuts = "", TString title = "",
           TString xaxis_title = "", TString yaxis_title = "",
           TString draw_opts = "colz") {
    TString hstr = TString::Format("%s >>%s%s", draw_cmd.Data(), name.Data(),
                                   binning.Data());
    tree.Draw(hstr, cuts, draw_opts);
    TH1 *hist = static_cast<TH1*>(gPad->GetPrimitive(name));
    if (yaxis_title != "") {
        hist->GetYaxis()->SetTitle(yaxis_title);
    }
    if (xaxis_title != "") {
        hist->GetXaxis()->SetTitle(xaxis_title);
    }
    if (title != "") {
        hist->SetTitle(title);
    }
}


void phi_correlation() {
  init_dataset();

    // TCutG* esr_phi_corr_gcut = (TCutG*)
    //     (new TFile("cuts/esr_phi_corr_gcut.root"))->Get("CUTG");
    // esr_phi_corr_gcut->SetName("esr_phi_corr_gcut");
    // TCutG* esl_phi_corr_gcut = (TCutG*)
    //     (new TFile("cuts/esl_phi_corr_gcut.root"))->Get("CUTG");
    // esr_phi_corr_gcut->SetName("esl_phi_corr_gcut");

  TCut p_angle_sn_cut = "p_theta*57.3>65 && p_theta*57.3<68";
  TCut p_angle_acceptance_cut = "p_theta*57.3>55 && p_theta*57.3<70";

    TCut esl_phi_corr_cut = "abs(-1.3 * s1dc_phi*57.3 + 25 - esl_p_phi*57.3) < 2.5";
    TCut esr_phi_corr_cut = "abs(-1.35 * s1dc_phi*57.3 - 33 - esr_p_phi*57.3) < 2.5";

    TCut phi_corr_cut_1d = "abs(abs(p_phi*57.3-s1dc_phi*57.3) - 180) < 8";

    TCut tgt_up_dwn_xpos_cut = "abs((tgt_up_xpos - tgt_dwn_xpos) - -3.32) < 7";
    TCut tgt_up_dwn_ypos_cut = "abs((tgt_up_ypos - tgt_dwn_ypos) - 1.13e-2) < 7";

    TCut target_cut{"tgt_up_xpos*tgt_up_xpos + tgt_up_ypos*tgt_up_ypos < 64"};
    TCut target_cut2{tgt_up_dwn_xpos_cut && tgt_up_dwn_ypos_cut};

    TCut vertex_zpos_cut{"es_vertex_zpos > -30 && es_vertex_zpos < 30"};

    // distribution on HODF plastics
    // pla IDs [7,9]: beam / beam stopper
    // pla ID [0,4]: left scattering 6He
    // pla ID [11,14]: right scattering 6He
    // pla ID [17,23]: 4He (from inelastic shannel)
    std::string hodf_cut_str = "";
    for (int i = 0; i < 5; i++) {
      hodf_cut_str += TString::Format("(hodf_q[%i]>14 && hodf_q[%i]<17) ", i, i);
      if (i != 4) hodf_cut_str += "|| ";
    }
    hodf_cut_str += "&& ";
    for (int i = 11; i < 15; i++) {
      hodf_cut_str += TString::Format("(hodf_q[%i]>14 && hodf_q[%i]<17) ", i, i);
      if (i != 14) hodf_cut_str += "|| ";
    }
    cout << "HODF cut:\n" << hodf_cut_str << std::endl;
    TCut hodf_cut{hodf_cut_str.c_str()};

    TFile hist_out("out/phi_angle_correlation.root", "RECREATE");
    TCanvas c1("c1");
    gStyle->SetOptStat(1111111);

    // proton phi
    c1.Clear();
    c1.Divide(2);
    c1.cd(1);
    hdraw(g_chain_up, "esl_p_phi", "esl_p_phi*57.3", "(360,-180,180)",
          "triggers[5]==1" && target_cut && hodf_cut,
          "ESPRI left, proton phi angle",
          "Proton #phi angle [lab. deg.]", "Counts");
    //gPad->SetLogy();
    c1.cd(2);
    hdraw(g_chain_up, "esr_p_phi", "esr_p_phi*57.3", "(360,-180,180)",
          "triggers[5]==1" && target_cut && hodf_cut,
          "ESPRI right, proton phi angle",
          "Proton #phi angle [lab. deg.]", "Counts");
    //gPad->SetLogy();
    c1.Print("out/phi_angle_correlation.pdf(", "pdf");

    // He phi
    c1.Clear();
    c1.Divide(1);
    c1.cd(1);
    hdraw(g_chain_up, "fragment_phi_s1dc", "s1dc_phi*57.3", "(360,-180,180)",
          "triggers[5]==1" && target_cut && hodf_cut,
          "Fragment azimuthal angle",
          "Fragment #phi angle [lab. deg.]", "Counts");
    //gPad->SetLogy();
    c1.Print("out/phi_angle_correlation.pdf", "pdf");

    // phi 2d corr ESPRI L&R
    c1.Clear();
    c1.Divide(2);
    c1.cd(1);
    hdraw(g_chain_up, "esl_phi_corr", "esl_p_phi*57.3:s1dc_phi*57.3",
          "(100,60,120,100,-120,-60)",
          "triggers[5]==1" && target_cut && hodf_cut,
          "ESPRI left, phi correlation",
          "Fragment #phi angle [lab. deg.]",
          "Proton #phi angle [lab. deg.]");
    c1.cd(2);
    hdraw(g_chain_up, "esr_phi_corr", "esr_p_phi*57.3:s1dc_phi*57.3",
          "(100,-120,-60,100,60,120)",
          "triggers[5]==1" && target_cut && hodf_cut,
          "ESPRI right, phi correlation",
          "Fragment #phi angle [lab. deg.]",
          "Proton #phi angle [lab. deg.]");
    c1.Print("out/phi_angle_correlation.pdf", "pdf");

    // projection of locus, angle determination
    // c1.Clear();
    // c1.Divide(2,2);
    // c1.cd(1);
    // hdraw(g_chain_up, "esl_phi_corr_cut1", "s1dc_phi*57.3",
    //       "(100,60,120)",
    //       "triggers[5]==1" && target_cut && hodf_cut &&
    //       "esl_p_phi*57.3>-84 && esl_p_phi*57.3<-80",
    //       "ESPRI left, azimuthal losus X proj (p phi: -84 - -80)",
    //       "Fragment #phi angle [lab. deg.]",
    //       "Counts");
    // c1.cd(3);
    // hdraw(g_chain_up, "esl_phi_corr_cut2", "s1dc_phi*57.3",
    //       "(100,60,120)",
    //       "triggers[5]==1" && target_cut && hodf_cut &&
    //       "esl_p_phi*57.3>-100 && esl_p_phi*57.3<-96",
    //       "ESPRI left, azimuthal losus X proj (p phi: -100 - -96)",
    //       "Fragment #phi angle [lab. deg.]",
    //       "Counts");
    // c1.cd(2);
    // hdraw(g_chain_up, "esr_phi_corr_cut3", "s1dc_phi*57.3",
    //       "(100,-120,-60)",
    //       "triggers[5]==1" && target_cut && hodf_cut &&
    //       "esr_p_phi*57.3>96 && esr_p_phi*57.3<100",
    //       "ESPRI right, azimuthal losus X proj (p phi: 96 - 100)",
    //       "Fragment #phi angle [lab. deg.]",
    //       "Counts");
    // c1.cd(4);
    // hdraw(g_chain_up, "esr_phi_corr_cu4", "s1dc_phi*57.3",
    //       "(100,-120,-60)",
    //       "triggers[5]==1" && target_cut && hodf_cut &&
    //       "esr_p_phi*57.3>80 && esr_p_phi*57.3<84",
    //       "ESPRI right, azimuthal losust X proj (p phi: 80 - 84)",
    //       "Fragment #phi angle [lab. deg.]",
    //       "Counts");
    // c1.Print("out/phi_angle_correlation.pdf", "pdf");

    // ESPRI L&R phi corr 1d
    c1.Clear();
    c1.Divide(1,2);
    c1.cd(1);
    hdraw(g_chain_up, "esl_phi_corr_1d", "esl_p_phi*57.3-s1dc_phi*57.3",
          "(300,-350,-10)",
          "triggers[5]==1" && target_cut && hodf_cut,
          "ESPRI left: P_{#phi} - 4He_{#phi}",
          "P_{#phi} - 4He_{#phi} [deg]",
          "Counts");
    //gPad->SetLogy();
    c1.cd(2);
    hdraw(g_chain_up, "esr_phi_corr_1d", "esr_p_phi*57.3-s1dc_phi*57.3",
          "(300,10,350)",
          "triggers[5]==1" && target_cut && hodf_cut,
          "ESPRI right: P_{#phi} - 4He_{#phi}",
          "P_{#phi} - 4He_{#phi} [deg]",
          "Counts");
    //gPad->SetLogy();
    c1.Print("out/phi_angle_correlation.pdf", "pdf");

    // ESPRI total phi cut
    c1.Clear();
    c1.Divide(1,2);
    c1.cd(1);
    hdraw(g_chain_up, "tot_phi_corr_1d", "abs(p_phi*57.3-s1dc_phi*57.3)",
          "(300,10,350)",
          "triggers[5]==1" && target_cut && hodf_cut,
          "ESPRI left & right: P_{#phi} - 4He_{#phi} (no #phi cut)",
          "P_{#phi} - 4He_{#phi} [deg.]",
          "Counts");
    //gPad->SetLogy();
    c1.cd(2);
    hdraw(g_chain_up, "tot_phi_corr_1d_cut", "abs(p_phi*57.3-s1dc_phi*57.3)",
          "(300,10,350)",
          "triggers[5]==1" && target_cut && hodf_cut &&
          phi_corr_cut_1d,
          "ESPRI left & right: P_{#phi} - 4He_{#phi} (#phi cut)",
          "P_{#phi} - 4He_{#phi} [deg.]",
          "Counts");
    //gPad->SetLogy();
    c1.Print("out/phi_angle_correlation.pdf", "pdf");

    // c1.Clear();
    // c1.Divide(1);
    // c1.cd(1);
    // hdraw(g_chain_up, "polar_tgt", "p_theta*57.3:s1dc_theta*57.3",
    //       "(200,0,20,200,50,75)",
    //       "triggers[5]==1" && target_cut,
    //       "Polar angle correlation - tgt cut",
    //       "Fragment #theta angle [lab. deg.]",
    //       "Proton #theta angle [lab. deg.]", "colz");
    // c1.Print("out/phi_angle_correlation.pdf", "pdf");

    // polar corr with diff. cuts

    // HODF cut
    c1.Clear();
    c1.Divide(1);
    c1.cd(1);
    hdraw(g_chain_up, "polar_phi_hodf", "p_theta*57.3:s1dc_theta*57.3",
          "(200,0.1,20,200,50,75)",
          "triggers[5]==1" && target_cut && hodf_cut,
          "Polar angle correlation - tgt, hodf cuts",
          "Fragment #theta angle [lab. deg.]",
          "Proton #theta angle [lab. deg.]", "colz");
    c1.Print("out/phi_angle_correlation.pdf", "pdf");

    c1.Clear();
    c1.Divide(1);
    c1.cd(1);
    hdraw(g_chain_up, "polar_tgt_phi1d", "p_theta*57.3:s1dc_theta*57.3",
          "(200,0,20,200,50,75)",
          "triggers[5]==1" && target_cut && phi_corr_cut_1d,
          "Polar angle correlation - tgt, #phi cuts",
          "Fragment #theta angle [lab. deg.]",
          "Proton #theta angle [lab. deg.]", "colz");
    c1.Print("out/phi_angle_correlation.pdf", "pdf");

    // HODF cut
    c1.Clear();
    c1.Divide(1);
    c1.cd(1);
    hdraw(g_chain_up, "polar_tgt_phi1d_hodf", "p_theta*57.3:s1dc_theta*57.3",
          "(200,0.1,20,200,50,75)",
          "triggers[5]==1" && target_cut &&
          phi_corr_cut_1d && hodf_cut,
          "Polar angle correlation - tgt, hodf, #phi cuts",
          "Fragment #theta angle [lab. deg.]",
          "Proton #theta angle [lab. deg.]", "colz");
    c1.Print("out/phi_angle_correlation.pdf", "pdf");

    // 1d hist theta corr S/N
    c1.Clear();
    c1.Divide(1);
    c1.cd(1);
    hdraw(g_chain_up, "polar_tgt_phi1d_hodf_1d", "s1dc_theta*57.3",
          "(200,0,20)",
          "triggers[5]==1" && target_cut &&
          hodf_cut && phi_corr_cut_1d && p_angle_sn_cut,
          "Polar angle cut S/N - tgt, #phi cuts",
          "Fragment #theta angle [lab. deg.]",
          "Counts", "colz");
    c1.Print("out/phi_angle_correlation.pdf", "pdf");

    // theta exp - theta_theor 2d corr
    c1.Clear();
    c1.Divide(1);
    c1.cd(1);
    hdraw(g_chain_up, "polar_t_tgt_hodf_phi",
          "p_theta*57.3:s1dc_theta*57.3 - he_theta_theor*57.3",
          "(200,-10,10,200,50,75)",
          "triggers[5]==1" && target_cut &&
          hodf_cut && phi_corr_cut_1d,
          "Polar angle correlation (exp - theor)",
          "Fragment: #theta_{exp} - #theta_{theor} [lab. deg.]",
          "Proton #theta angle [lab. deg.]", "colz");
    c1.Print("out/phi_angle_correlation.pdf", "pdf");

    // theta exp - theta_theor 1d corr
    c1.Clear();
    c1.Divide(1);
    c1.cd(1);
    hdraw(g_chain_up, "polar10_1d",
          "s1dc_theta*57.3 - he_theta_theor*57.3",
          "(200,-10,10)",
          "triggers[5]==1" && target_cut &&
          hodf_cut && phi_corr_cut_1d && p_angle_acceptance_cut,
          "Polar angle correlation (exp - theor)",
          "Fragment: #theta_{exp} - #theta_{theor} [lab. deg.]",
          "Counts", "colz");
    c1.Print("out/phi_angle_correlation.pdf)", "pdf");

    hist_out.Write();
    hist_out.Close();
}
