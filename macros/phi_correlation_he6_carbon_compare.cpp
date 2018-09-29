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


void phi_correlation_he6_carbon_compare() {
  init_dataset();
  init_carbon_dataset();

  std::cout << "Total number of events for he6: " <<  
    g_chain_total.GetEntries() << std::endl;
  std::cout << "Total number of events for carbon: " <<
    g_c_chain.GetEntries() << std::endl;

  TCut p_angle_sn_cut = "p_theta*57.3>65 && p_theta*57.3<68";
  TCut p_angle_acceptance_cut = "p_theta*57.3>55 && p_theta*57.3<70";

  TCut phi_corr_cut_1d = "abs(abs(p_phi*57.3-s1dc_phi*57.3) - 180) < 8";

  TCut target_cut{"tgt_up_xpos*tgt_up_xpos + tgt_up_ypos*tgt_up_ypos < 81"};

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
  hodf_cut_str += "|| ";
  for (int i = 11; i < 15; i++) {
    hodf_cut_str += TString::Format("(hodf_q[%i]>14 && hodf_q[%i]<17) ", i, i);
    if (i != 14) hodf_cut_str += "|| ";
  }
  cout << "HODF cut:\n" << hodf_cut_str << std::endl;
  TCut hodf_cut{hodf_cut_str.c_str()};

  TFile hist_out("out/phi_correlation_he6_carbon_compare.root", "RECREATE");
  TCanvas c1("c1");
  gStyle->SetOptStat(1111111);

  // // proton phi
  // c1.Clear();
  // c1.Divide(2);
  // c1.cd(1);
  // hdraw(g_chain_total, "esl_p_phi", "esl_p_phi*57.3", "(360,-180,180)",
  //       "triggers[5]==1" && target_cut && hodf_cut,
  //       "ESPRI left, proton phi angle",
  //       "Proton #phi angle [lab. deg.]", "Counts");
  // //gPad->SetLogy();
  // c1.cd(2);
  // hdraw(g_chain_total, "esr_p_phi", "esr_p_phi*57.3", "(360,-180,180)",
  //       "triggers[5]==1" && target_cut && hodf_cut,
  //       "ESPRI right, proton phi angle",
  //       "Proton #phi angle [lab. deg.]", "Counts");
  // //gPad->SetLogy();
  // c1.Print("out/phi_correlation_he6_carbon_compare.pdf(", "pdf");

  // // He phi
  // c1.Clear();
  // c1.Divide(1);
  // c1.cd(1);
  // hdraw(g_chain_total, "fragment_phi_s1dc", "s1dc_phi*57.3", "(360,-180,180)",
  //       "triggers[5]==1" && target_cut && hodf_cut,
  //       "Fragment azimuthal angle",
  //       "Fragment #phi angle [lab. deg.]", "Counts");
  // //gPad->SetLogy();
  // c1.Print("out/phi_correlation_he6_carbon_compare.pdf", "pdf");

  // // phi 2d corr ESPRI L&R
  // c1.Clear();
  // c1.Divide(2);
  // c1.cd(1);
  // hdraw(g_chain_total, "esl_phi_corr", "esl_p_phi*57.3:s1dc_phi*57.3",
  //       "(100,60,120,100,-120,-60)",
  //       "triggers[5]==1" && target_cut && hodf_cut,
  //       "ESPRI left, phi correlation",
  //       "Fragment #phi angle [lab. deg.]",
  //       "Proton #phi angle [lab. deg.]");
  // c1.cd(2);
  // hdraw(g_chain_total, "esr_phi_corr", "esr_p_phi*57.3:s1dc_phi*57.3",
  //       "(100,-120,-60,100,60,120)",
  //       "triggers[5]==1" && target_cut && hodf_cut,
  //       "ESPRI right, phi correlation",
  //       "Fragment #phi angle [lab. deg.]",
  //       "Proton #phi angle [lab. deg.]");
  // c1.Print("out/phi_correlation_he6_carbon_compare.pdf", "pdf");

  // ESPRI L&R phi corr 1d
  // c1.Clear();
  // c1.Divide(1,2);
  // c1.cd(1);
  // hdraw(g_chain_total, "esl_phi_corr_1d", "esl_p_phi*57.3-s1dc_phi*57.3",
  //       "(300,-350,-10)",
  //       "triggers[5]==1" && target_cut && hodf_cut,
  //       "ESPRI left: P_{#phi} - He_{#phi}",
  //       "P_{#phi} - He_{#phi} [deg]",
  //       "Counts");
  // //gPad->SetLogy();
  // c1.cd(2);
  // hdraw(g_chain_total, "esr_phi_corr_1d", "esr_p_phi*57.3-s1dc_phi*57.3",
  //       "(300,10,350)",
  //       "triggers[5]==1" && target_cut && hodf_cut,
  //       "ESPRI right: P_{#phi} - He_{#phi}",
  //       "P_{#phi} - He_{#phi} [deg]",
  //       "Counts");
  // //gPad->SetLogy();
  // c1.Print("out/phi_correlation_he6_carbon_compare.pdf", "pdf");

  // // ESPRI total -  phi, hodf cuts
  // c1.Clear();
  // c1.Divide(1,2);
  // c1.cd(1);
  // hdraw(g_chain_total, "tot_phi_corr_1d", "abs(p_phi*57.3-s1dc_phi*57.3)",
  //       "(300,10,350)",
  //       "triggers[5]==1" && target_cut && hodf_cut,
  //       "ESPRI left & right: P_{#phi} - He_{#phi} (no #phi cut)",
  //       "P_{#phi} - He_{#phi} [deg.]",
  //       "Counts");
  // //gPad->SetLogy();
  // c1.cd(2);
  // hdraw(g_chain_total, "tot_phi_corr_1d_cut", "abs(p_phi*57.3-s1dc_phi*57.3)",
  //       "(300,10,350)",
  //       "triggers[5]==1" && target_cut && hodf_cut &&
  //       phi_corr_cut_1d,
  //       "ESPRI left & right: P_{#phi} - He_{#phi} (#phi cut)",
  //       "P_{#phi} - He_{#phi} [deg.]",
  //       "Counts");
  // //gPad->SetLogy();
  // c1.Print("out/phi_correlation_he6_carbon_compare.pdf", "pdf");

  // polar corr with diff. cuts

  // // HODF cut - He6 vs. carbon
  // c1.Clear();
  // c1.Divide(2,1);
  // c1.cd(1);
  // hdraw(g_chain_total, "polar_phi_hodf", "p_theta*57.3:s1dc_theta*57.3",
  //       "(200,0.1,10,200,50,75)",
  //       "triggers[5]==1" && target_cut && hodf_cut,
  //       "#theta angle correlation for 6he-p (tgt, hodf cuts)",
  //       "Fragment #theta angle [lab. deg.]",
  //       "Proton #theta angle [lab. deg.]", "colz");
  // c1.cd(2);
  // hdraw(g_c_chain, "polar_phi_hodf2", "p_theta*57.3:s1dc_theta*57.3",
  //       "(200,0.1,10,200,50,75)",
  //       "triggers[5]==1" && target_cut && hodf_cut,
  //       "#theta angle correlation for 6he-C (tgt, hodf cuts)",
  //       "Fragment #theta angle [lab. deg.]",
  //       "Proton #theta angle [lab. deg.]", "colz");
  // c1.Print("out/phi_correlation_he6_carbon_compare.pdf(", "pdf");

  // // tgt,phi,hodf cut - He6 vs. carbon
  // c1.Clear();
  // c1.Divide(2,1);
  // c1.cd(1);
  // hdraw(g_chain_total, "polar_tgt_phi1d_hodf", "p_theta*57.3:s1dc_theta*57.3",
  //       "(200,0.1,10,200,50,75)",
  //       "triggers[5]==1" && target_cut &&
  //       phi_corr_cut_1d && hodf_cut,
  //       "#theta angle corr.: 6he-p (tgt, hodf, phi cuts)",
  //       "Fragment #theta angle [lab. deg.]",
  //       "Proton #theta angle [lab. deg.]", "colz");
  // c1.cd(2);
  // hdraw(g_c_chain, "polar_tgt_phi1d_hodf2", "p_theta*57.3:s1dc_theta*57.3",
  //       "(200,0.1,10,200,50,75)",
  //       "triggers[5]==1" && target_cut &&
  //       phi_corr_cut_1d && hodf_cut,
  //       "#theta angle corr.: 6he-C (tgt, hodf, phi cuts)",
  //       "Fragment #theta angle [lab. deg.]",
  //       "Proton #theta angle [lab. deg.]", "colz");
  // c1.Print("out/phi_correlation_he6_carbon_compare.pdf", "pdf");

  // // 1d hist theta corr S/N
  // c1.Clear();
  // c1.Divide(2,1);
  // c1.cd(1);
  // hdraw(g_chain_total, "polar_tgt_phi1d_hodf_1d", "s1dc_theta*57.3",
  //       "(200,0,10)",
  //       "triggers[5]==1" && target_cut &&
  //       hodf_cut && phi_corr_cut_1d && p_angle_sn_cut,
  //       "#theta angle corr.:  6he-p (tgt, hodf, phi cuts)",
  //       "Fragment #theta angle [lab. deg.]",
  //       "Counts", "colz");
  // c1.cd(2);
  // hdraw(g_c_chain, "polar_tgt_phi1d_hodf_1d2", "s1dc_theta*57.3",
  //       "(200,0,10)",
  //       "triggers[5]==1" && target_cut &&
  //       hodf_cut && phi_corr_cut_1d && p_angle_sn_cut,
  //       "#theta angle corr.:  6he-C (tgt, hodf, phi cuts)",
  //       "Fragment #theta angle [lab. deg.]",
  //       "Counts", "colz");
  // c1.Print("out/phi_correlation_he6_carbon_compare.pdf", "pdf");

  // // theta exp - theta_theor 2d corr
  c1.Clear();
  c1.Divide(2,1);
  c1.cd(1);
  hdraw(g_chain_total, "polar_t_tgt_hodf_phi",
        "p_theta*57.3:s1dc_theta*57.3 - he_theta_theor*57.3",
        "(200,-8,8,200,50,75)",
        "triggers[5]==1" && target_cut &&
        hodf_cut && phi_corr_cut_1d,
        "#theta angle corr.: #theta_{exp}-#theta_{theor} for 6he-p (tgt, hodf cuts)",
        "Fragment: #theta_{exp} - #theta_{theor} [lab. deg.]",
        "Proton #theta angle [lab. deg.]", "colz");
  c1.cd(2);
  hdraw(g_c_chain, "polar_t_tgt_hodf_phi2",
        "p_theta*57.3:s1dc_theta*57.3 - he_theta_theor*57.3",
        "(200,-8,8,200,50,75)",
        "triggers[5]==1" && target_cut &&
        hodf_cut && phi_corr_cut_1d,
        "#theta angle corr.: #theta_{exp}-#theta_{theor} for 6he-C (tgt, hodf  cuts)",
        "Fragment: #theta_{exp} - #theta_{theor} [lab. deg.]",
        "Proton #theta angle [lab. deg.]", "colz");
  c1.Print("out/phi_correlation_he6_carbon_compare.pdf(", "pdf");

  // // theta exp - theta_theor 1d corr
  // c1.Clear();
  // c1.Divide(2,1);
  // c1.cd(1);
  // hdraw(g_chain_total, "polar_t_1d",
  //       "s1dc_theta*57.3 - he_theta_theor*57.3",
  //       "(200,-8,8)",
  //       "triggers[5]==1" && target_cut &&
  //       hodf_cut && p_angle_acceptance_cut,
  //       "#theta angle corr.: #theta_{exp}-#theta_{theor} for 6he-p (tgt, hodf cuts)",
  //       "Fragment: #theta_{exp} - #theta_{theor} [lab. deg.]",
  //       "Counts", "colz");
  // gPad->SetLogy();
  // c1.cd(2);
  // hdraw(g_c_chain, "polar_t_1d2",
  //       "s1dc_theta*57.3 - he_theta_theor*57.3",
  //       "(200,-8,8)",
  //       "triggers[5]==1" && target_cut &&
  //       hodf_cut && p_angle_acceptance_cut,
  //       "#theta angle corr.: #theta_{exp}-#theta_{theor} for 6he-C (tgt, hodf cuts)",
  //       "Fragment: #theta_{exp} - #theta_{theor} [lab. deg.]",
  //       "Counts", "colz");
  // gPad->SetLogy();
  // c1.Print("out/phi_correlation_he6_carbon_compare.pdf", "pdf");
  
  // // theta exp - theta_theor 1d corr + carbon BG
  // c1.Clear();
  // c1.Divide(1);
  // c1.cd(1);
  // g_chain_total.SetLineColor(kBlue);
  // hdraw(g_chain_total, "polar_t_1d_hc1",
  //       "s1dc_theta*57.3 - he_theta_theor*57.3",
  //       "(200,-8,8)",
  //       "triggers[5]==1" && target_cut &&
  //       hodf_cut && p_angle_acceptance_cut,
  //       "#theta angle corr.: #theta_{exp}-#theta_{theor} for 6he-p and 6he-C (tgt, hodf cuts)",
  //       "Fragment: #theta_{exp} - #theta_{theor} [lab. deg.]",
  //       "Counts");
  // g_c_chain.SetLineColor(kRed);
  // hdraw(g_c_chain, "polar_t_1d_hc2",
  //       "s1dc_theta*57.3 - he_theta_theor*57.3",
  //       "(200,-8,8)",
  //       "triggers[5]==1" && target_cut &&
  //       hodf_cut && p_angle_acceptance_cut,
  //       "#theta angle corr.: #theta_{exp}-#theta_{theor} for 6he-p and 6he-C (tgt, hodf cuts)",
  //       "Fragment: #theta_{exp} - #theta_{theor} [lab. deg.]",
  //       "Counts", "SAME");
  // gPad->SetLogy();
  // c1.Print("out/phi_correlation_he6_carbon_compare.pdf", "pdf");

// theta exp - theta_theor 1d corr + phi cut
  c1.Clear();
  c1.Divide(2,1);
  c1.cd(1);
  hdraw(g_chain_total, "polar_t_1d_phi",
        "s1dc_theta*57.3 - he_theta_theor*57.3",
        "(200,-8,8)",
        "triggers[5]==1" && target_cut &&
        hodf_cut && phi_corr_cut_1d && p_angle_acceptance_cut,
        "#theta angle corr.: #theta_{exp}-#theta_{theor} for 6he-p (tgt, hodf, phi cuts)",
        "Fragment: #theta_{exp} - #theta_{theor} [lab. deg.]",
        "Counts", "colz");
  gPad->SetLogy();
  c1.cd(2);
  g_c_chain.SetLineColor(kBlue);
  hdraw(g_c_chain, "polar_t_1d2_phi",
        "s1dc_theta*57.3 - he_theta_theor*57.3",
        "(200,-8,8)",
        "triggers[5]==1" && target_cut &&
        hodf_cut && phi_corr_cut_1d && p_angle_acceptance_cut,
        "#theta angle corr.: #theta_{exp}-#theta_{theor} for 6he-C (tgt, hodf, phi cuts)",
        "Fragment: #theta_{exp} - #theta_{theor} [lab. deg.]",
        "Counts", "colz");
  gPad->SetLogy();
  c1.Print("out/phi_correlation_he6_carbon_compare.pdf", "pdf");
  
  // theta exp - theta_theor 1d corr + phi cut + carbon BG 
  c1.Clear();
  c1.Divide(1);
  c1.cd(1);
  g_chain_total.SetLineColor(kBlue);
  hdraw(g_chain_total, "polar_t_1d_hc1_phi",
        "s1dc_theta*57.3 - he_theta_theor*57.3",
        "(200,-8,8)",
        "triggers[5]==1" && target_cut &&
        hodf_cut && phi_corr_cut_1d && p_angle_acceptance_cut,
        "#theta angle corr.: #theta_{exp}-#theta_{theor} for 6he-p and 6he-C (tgt, hodf, phi cuts)",
        "Fragment: #theta_{exp} - #theta_{theor} [lab. deg.]",
        "Counts");
  g_c_chain.SetLineColor(kRed);
  hdraw(g_c_chain, "polar_t_1d_hc2_phi",
        "s1dc_theta*57.3 - he_theta_theor*57.3",
        "(200,-8,8)",
        "triggers[5]==1" && target_cut &&
        hodf_cut && phi_corr_cut_1d && p_angle_acceptance_cut,
        "#theta angle corr.: #theta_{exp}-#theta_{theor} for 6he-p and 6he-C (tgt, hodf, phi cuts)",
        "Fragment: #theta_{exp} - #theta_{theor} [lab. deg.]",
        "Counts", "SAME");
  TH1* he6_polar = static_cast<TH1*>(gPad->GetPrimitive("polar_t_1d_hc1_phi"));
  TH1* carbon_polar = static_cast<TH1*>(gPad->GetPrimitive("polar_t_1d_hc2_phi"));
  carbon_polar->Scale(3);
  he6_polar->Draw();
  carbon_polar->SetLineColor(kRed);
  carbon_polar->Draw("SAME");
  gPad->SetLogy();
  c1.Print("out/phi_correlation_he6_carbon_compare.pdf", "pdf");

  // theta exp - theta_theor 1d corr - carbon BG 
  c1.Clear();
  c1.Divide(1);
  c1.cd(1);
  g_chain_total.SetLineColor(kBlue);
  g_c_chain.SetLineColor(kBlue);
  TH1* he6_polar_min_bg = (TH1*)he6_polar->Clone();
  he6_polar_min_bg->Add(carbon_polar, -1);
  he6_polar_min_bg->SetTitle("#theta angle corr.: #theta_{exp}-#theta_{theor} for "
                             "6he-p - 6he-C (tgt, hodf, phi cuts)");
  he6_polar_min_bg->Draw();
  gPad->SetLogy();
  c1.Print("out/phi_correlation_he6_carbon_compare.pdf", "pdf");

  c1.Clear();
  c1.Divide(1);
  c1.cd(1);
  he6_polar->Draw();
  carbon_polar->SetLineColor(kRed);
  carbon_polar->Draw("SAME");
  //gPad->SetLogy();
  c1.Print("out/phi_correlation_he6_carbon_compare.pdf", "pdf");
  
  c1.Clear();
  c1.Divide(1);
  c1.cd(1);
  g_chain_total.SetLineColor(kBlue);
  g_c_chain.SetLineColor(kBlue);
  he6_polar_min_bg->Draw();
  //gPad->SetLogy();
  c1.Print("out/phi_correlation_he6_carbon_compare.pdf)", "pdf");

  
  hist_out.Write();
  hist_out.Close();
}
