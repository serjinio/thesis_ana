#include <assert.h>

#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <memory>
#include <cmath>

#include <TMath.h>
#include <TVector3.h>

#include "dataset.h"
#include "protontrajsim.hpp"
#include "scattree.hpp"


using namespace s13::ana;
using namespace s13::common;


////////////////////////////////////////////////////////////
// RawTreeConverter implementation
////////////////////////////////////////////////////////////

bool RawTreeConverter::IsHe6Run() {
  if (m_run_no >= 133 && m_run_no <= 269) { // he6 runs
    return true;
  }

  return false;
}

bool RawTreeConverter::IsCarbonRun() {
  if (m_run_no >= 324 && m_run_no <= 350) { // he6 runs
    return true;
  }

  return false;
}

bool RawTreeConverter::IsHe4Run() {
  if (m_run_no >= 272 && m_run_no <= 316) { // he4 runs
    return true;
  }

  return false;
}

bool RawTreeConverter::IsEspriLeftEvent() {
  //return m_proton_phi < 0;
  return m_espri_left_de > 230;
}

bool RawTreeConverter::IsEspriRightEvent() {
  //return m_proton_phi > 0;
  return m_espri_right_de > 230;
}

bool RawTreeConverter::IsWellDefinedValue(double value) {
  if (value == INIT_VAL) {
    return false;
  } else {
    return true;
  }
}

void RawTreeConverter::Initialize() {
  SetupOutputTree();
  SetupInputTree();

  // load kinematics data
  m_he4_kinematics = LoadKinematicsData("data/he4_kin.csv");
  m_he6_kinematics = LoadKinematicsData("data/he6_kin.csv");

  // load stopping powers
  m_pla_stp_power = LoadStoppingPower("data/p_in_BC400.csv");
  m_edeg_stp_power = LoadStoppingPower("data/p_in_brass.csv");

  // load calibration constants
  m_espri_left_nai_cal_consts = LoadEspriNaiCalConsts("data/esl_nai_calib.csv");
  m_espri_right_nai_cal_consts = LoadEspriNaiCalConsts("data/esr_nai_calib.csv");
  m_espri_left_nai_pedestals = LoadEspriNaiCalConsts("data/esl_nai_calib.csv", true);
  m_espri_right_nai_pedestals = LoadEspriNaiCalConsts("data/esr_nai_calib.csv", true);

  std::cout << "Left ESPRI NaIs calibration constants (pedestal [ch]; cal. fact [MeV/ch]):"
            << std::endl;
  for (int i = 0; i < ESPRI_NUM_NAIS; ++i) {
    std::cout << "\t" << m_espri_left_nai_pedestals[i] << ";\t"
              << m_espri_left_nai_cal_consts[i] << std::endl;
  }
  std::cout << "Right ESPRI NaIs calibration constants (pedestal [ch]; cal. fact [MeV/ch]):"
            << std::endl;
  for (int i = 0; i < ESPRI_NUM_NAIS; ++i) {
    std::cout << "\t" << m_espri_right_nai_pedestals[i] << ";\t"
              << m_espri_right_nai_cal_consts[i] << std::endl;
  }

  auto pde_consts = LoadEspriPdeCalConsts("data/es_pde_calib.csv");
  auto pde_pedestals = LoadEspriPdeCalConsts("data/es_pde_calib.csv", true);
  m_espri_left_pde_pedestal = pde_pedestals[0];
  m_espri_right_pde_pedestal = pde_pedestals[1];
  m_espri_left_pde_cal_const = pde_consts[0];
  m_espri_right_pde_cal_const = pde_consts[1];

  std::cout << "ESPRI pdE calibration constants (pedestal [ch]; cal. fact [MeV/ch]):"
            << std::endl;
  std::cout << "\tLeft pdE:  " << m_espri_left_pde_pedestal << ";\t"
             << m_espri_left_pde_cal_const << std::endl;
  std::cout << "\tRight pdE: " << m_espri_right_pde_pedestal << ";\t"
            << m_espri_right_pde_cal_const << std::endl;

  // helper, simulated values - proton E & dE deposits for p-6He case
  auto interps_he6 = ComputeRefHeDeE("data/he6_kin_E.csv");
  m_he6_sim_p_de = std::get<0>(interps_he6);
  m_he6_sim_p_e = std::get<1>(interps_he6);

  // helper, simulated values - proton E & dE deposits for p-4He case
  auto interps_he4 = ComputeRefHeDeE("data/he4_kin_E.csv");
  m_he4_sim_p_de = std::get<0>(interps_he4);
  m_he4_sim_p_e = std::get<1>(interps_he4);

  // proton angle correction data (due to target's magnetic field)
  m_p_he4_esl_theta_eff_vs_meas =
    LoadProtonAngleCorrections("data/p_4he_esl_theta_eff_vs_theta_meas.csv");
  m_p_he4_esr_theta_eff_vs_meas =
    LoadProtonAngleCorrections("data/p_4he_esr_theta_eff_vs_theta_meas.csv");
  m_p_he6_esl_theta_eff_vs_meas =
    LoadProtonAngleCorrections("data/p_6he_esl_theta_eff_vs_theta_meas.csv");
  m_p_he6_esr_theta_eff_vs_meas =
    LoadProtonAngleCorrections("data/p_6he_esr_theta_eff_vs_theta_meas.csv");

  m_p_he4_esl_theta_eff_vs_mtrls_e_loss =
    LoadProtonMtrlsEloss("data/p_4he_esl_theta_eff_vs_theta_meas.csv");
  m_p_he4_esr_theta_eff_vs_mtrls_e_loss =
    LoadProtonMtrlsEloss("data/p_4he_esr_theta_eff_vs_theta_meas.csv");
  m_p_he6_esl_theta_eff_vs_mtrls_e_loss =
    LoadProtonMtrlsEloss("data/p_6he_esl_theta_eff_vs_theta_meas.csv");
  m_p_he6_esr_theta_eff_vs_mtrls_e_loss =
    LoadProtonMtrlsEloss("data/p_6he_esr_theta_eff_vs_theta_meas.csv");

  std::cout << "Raw tree converter initialized." << std::endl;
}

TTree* RawTreeConverter::Run() {
  assert(mp_output_tree != nullptr);
  assert(mp_raw_tree != nullptr);
  if (m_converted) {
    return mp_output_tree;
  }

  std::cout << "Starting conversion of a raw tree..." << std::endl;
  std::cout << "Number of events in the input tree: " << \
    mp_raw_tree->GetEntries() << std::endl;
  for (int i_entry = 0; i_entry < mp_raw_tree->GetEntries(); i_entry++) {
    ClearTreeVars();
    mp_raw_tree->GetEntry(i_entry);

    ComputeDerivedValues();
    ComputeTheoreticalValues();

    // Add new event into output tree
    mp_output_tree->Fill();

    if (i_entry % 50000 == 0) {
      std::cout << "Converted: " << i_entry << " events." << std::endl;
    }
  }

  mp_output_tree->GetCurrentFile()->Write();

  if (mp_raw_tree->GetEntries() != mp_output_tree->GetEntries()) {
    std::cout << "ERROR: Output tree event number is different " << \
      "from input!!!" << std::endl;
    throw "Output tree event number is different!";
  }

  std::cout << "Tree conversion finished. Output filename: " << \
    mp_output_file->GetName() << std::endl;
  std::cout << "Total event number converted: " << \
    mp_output_tree->GetEntries() << std::endl;
  std::cout << "Number of events with good proton data: " << \
    m_good_proton_events_num << std::endl;
  std::cout << "Number of events with bad proton data: " << \
    m_bad_proton_events_num << std::endl;
  std::cout << "Number of events with invalid NaI data: " << \
    m_invalid_nai_data_events_num << std::endl;
  std::cout << "Number of invalid NaI PMT values: " <<
    m_invalid_nai_pmt_values_num << std::endl;

  m_converted = true;
  mp_output_file->Close();

  return mp_output_tree;
}

void RawTreeConverter::SetupOutputTree() {
  mp_output_file = new TFile(
                             Form("scattree_rootfiles/scattree_run%i.root", m_run_no),
                             "recreate");
  mp_output_tree = new TTree("scattree", "scattree");

  mp_output_tree->Branch("event_no", &m_event_no,
                         "event_no/I");
  mp_output_tree->Branch("triggers", &m_triggers,
                         "triggers[8]/O");

  // beamline plastics
  mp_output_tree->Branch("f7_q", &m_f7_q, "f7_q/D");
  mp_output_tree->Branch("f7_t", &m_f7_t, "f7_t/D");
  mp_output_tree->Branch("sbt1_q", &m_sbt1_q, "sbt1_q/D");
  mp_output_tree->Branch("sbt1_t", &m_sbt1_t, "sbt1_t/D");
  mp_output_tree->Branch("sbt2_q", &m_sbt2_q, "sbt2_q/D");
  mp_output_tree->Branch("sbt2_t", &m_sbt2_t, "sbt2_t/D");
  mp_output_tree->Branch("sbv_q", &m_sbv_q, "sbv_q/D");
  mp_output_tree->Branch("sbv_t", &m_sbv_t, "sbv_t/D");

  // BDCs
  mp_output_tree->Branch("bdc1_xpos", &m_bdc1_xpos, "bdc1_xpos/D");
  mp_output_tree->Branch("bdc1_ypos", &m_bdc1_ypos, "bdc1_ypos/D");
  mp_output_tree->Branch("bdc2_xpos", &m_bdc2_xpos, "bdc2_xpos/D");
  mp_output_tree->Branch("bdc2_ypos", &m_bdc2_ypos, "bdc2_ypos/D");
  mp_output_tree->Branch("tgt_up_aang", &m_tgt_up_aang, "tgt_up_aang/D");
  mp_output_tree->Branch("tgt_up_bang", &m_tgt_up_bang, "tgt_up_bang/D");
  mp_output_tree->Branch("tgt_up_theta", &m_tgt_up_theta, "tgt_up_theta/D");
  mp_output_tree->Branch("tgt_up_phi", &m_tgt_up_phi, "tgt_up_phi/D");
  mp_output_tree->Branch("tgt_up_xpos", &m_tgt_up_xpos, "tgt_up_xpos/D");
  mp_output_tree->Branch("tgt_up_ypos", &m_tgt_up_ypos, "tgt_up_ypos/D");

  // DSDCs
  mp_output_tree->Branch("fdc0_xpos", &m_fdc0_xpos, "fdc0_xpos/D");
  mp_output_tree->Branch("fdc0_ypos", &m_fdc0_ypos, "fdc0_ypos/D");
  mp_output_tree->Branch("fdc0_aang", &m_fdc0_aang, "fdc0_aang/D");
  mp_output_tree->Branch("fdc0_bang", &m_fdc0_bang, "fdc0_bang/D");
  mp_output_tree->Branch("fdc0_theta", &m_fdc0_theta, "fdc0_theta/D");
  mp_output_tree->Branch("fdc0_phi", &m_fdc0_phi, "fdc0_phi/D");
  mp_output_tree->Branch("s1dc_xpos", &m_s1dc_xpos, "s1dc_xpos/D");
  mp_output_tree->Branch("s1dc_ypos", &m_s1dc_ypos, "s1dc_ypos/D");
  mp_output_tree->Branch("s1dc_aang", &m_s1dc_aang, "s1dc_aang/D");
  mp_output_tree->Branch("s1dc_bang", &m_s1dc_bang, "s1dc_bang/D");
  mp_output_tree->Branch("s1dc_theta", &m_s1dc_theta, "s1dc_theta/D");
  mp_output_tree->Branch("s1dc_phi", &m_s1dc_phi, "s1dc_phi/D");
  mp_output_tree->Branch("tgt_dwn_aang", &m_tgt_dwn_aang, "tgt_dwn_aang/D");
  mp_output_tree->Branch("tgt_dwn_bang", &m_tgt_dwn_bang, "tgt_dwn_bang/D");
  mp_output_tree->Branch("tgt_dwn_theta", &m_tgt_dwn_theta, "tgt_dwn_theta/D");
  mp_output_tree->Branch("tgt_dwn_phi", &m_tgt_dwn_phi, "tgt_dwn_phi/D");
  mp_output_tree->Branch("tgt_dwn_xpos", &m_tgt_dwn_xpos, "tgt_dwn_xpos/D");
  mp_output_tree->Branch("tgt_dwn_ypos", &m_tgt_dwn_ypos, "tgt_dwn_ypos/D");

  // ESPRI
  mp_output_tree->Branch("esl_xpos", &m_espri_left_xpos, "esl_xpos/D");
  mp_output_tree->Branch("esl_ypos", &m_espri_left_ypos, "esl_ypos/D");

  // NOTE: RDCs report angles in mrads, better to convert them to rads
  mp_output_tree->Branch("esl_aang", &m_espri_left_aang, "esl_truoaang/D");
  mp_output_tree->Branch("esl_bang", &m_espri_left_bang, "esl_bang/D");

  mp_output_tree->Branch("esl_aang_tgt", &m_espri_left_aang_tgt,
                         "esl_aang_tgt/D");
  mp_output_tree->Branch("esl_bang_tgt", &m_espri_left_bang_tgt,
                         "esl_bang_tgt/D");
  mp_output_tree->Branch("esl_de", &m_espri_left_de, "esl_de/D");
  mp_output_tree->Branch("esl_de_cal", &m_espri_left_de_cal, "esl_de_cal/D");
  mp_output_tree->Branch("esl_tof", &m_espri_left_tof, "esl_tof/D");
  mp_output_tree->Branch("esl_naie", &m_espri_left_nai_e, "esl_naie[7]/D");
  mp_output_tree->Branch("esl_naie_cal", &m_espri_left_nai_e_cal, "esl_naie_cal[7]/D");
  mp_output_tree->Branch("esl_e", &m_espri_left_e, "esl_e/D");
  mp_output_tree->Branch("esl_e_cal", &m_espri_left_e_cal, "esl_e_cal/D");
  mp_output_tree->Branch("esl_e_eff", &m_espri_left_e_eff, "esl_e_eff/D");
  mp_output_tree->Branch("esl_vertex_zpos", &m_espri_left_vertex_zpos,
                         "esl_vertex_zpos/D");
  mp_output_tree->Branch("esl_p_theta", &m_espri_left_proton_theta,
                         "esl_p_theta/D");
  mp_output_tree->Branch("esl_p_phi", &m_espri_left_proton_phi,
                         "esl_p_phi/D");
  mp_output_tree->Branch("esl_rdc_nhitwires", &m_espri_left_rdc_nhitwires,
                         "esl_rdc_nhitwires/I");
  mp_output_tree->Branch("esl_rdc_xnhitwires", &m_espri_left_rdc_xnhitwires,
                         "esl_rdc_xnhitwires/I");
  mp_output_tree->Branch("esl_rdc_ynhitwires", &m_espri_left_rdc_ynhitwires,
                         "esl_rdc_ynhitwires/I");

  mp_output_tree->Branch("esr_xpos", &m_espri_right_xpos, "esr_xpos/D");
  mp_output_tree->Branch("esr_ypos", &m_espri_right_ypos, "esr_ypos/D");

  // NOTE: RDCs report angles in mrads, better to convert them to rads
  mp_output_tree->Branch("esr_aang", &m_espri_right_aang, "esr_aang/D");
  mp_output_tree->Branch("esr_bang", &m_espri_right_bang, "esr_bang/D");

  mp_output_tree->Branch("esr_aang_tgt", &m_espri_right_aang_tgt,
                         "esr_aang_tgt/D");
  mp_output_tree->Branch("esr_bang_tgt", &m_espri_right_bang_tgt,
                         "esr_bang_tgt/D");
  mp_output_tree->Branch("esr_de", &m_espri_right_de, "esr_de/D");
  mp_output_tree->Branch("esr_de_cal", &m_espri_right_de_cal, "esr_de_cal/D");
  mp_output_tree->Branch("esr_tof", &m_espri_right_tof, "esr_tof/D");
  mp_output_tree->Branch("esr_naie", &m_espri_right_nai_e, "esr_naie[7]/D");
  mp_output_tree->Branch("esr_naie_cal", &m_espri_right_nai_e_cal, "esr_naie_cal[7]/D");
  mp_output_tree->Branch("esr_e", &m_espri_right_e, "esr_e/D");
  mp_output_tree->Branch("esr_e_cal", &m_espri_right_e_cal, "esr_e_cal/D");
  mp_output_tree->Branch("esr_e_eff", &m_espri_right_e_eff, "esr_e_eff/D");
  mp_output_tree->Branch("esr_vertex_zpos", &m_espri_right_vertex_zpos,
                         "esr_vertex_zpos/D");
  mp_output_tree->Branch("esr_p_theta", &m_espri_right_proton_theta,
                         "esr_p_theta/D");
  mp_output_tree->Branch("esr_p_phi", &m_espri_right_proton_phi,
                         "esr_p_phi/D");
  mp_output_tree->Branch("esr_rdc_nhitwires", &m_espri_right_rdc_nhitwires,
                         "esr_rdc_nhitwires/I");
  mp_output_tree->Branch("esr_rdc_xnhitwires", &m_espri_right_rdc_xnhitwires,
                         "esr_rdc_xnhitwires/I");
  mp_output_tree->Branch("esr_rdc_ynhitwires", &m_espri_right_rdc_ynhitwires,
                         "esr_rdc_ynhitwires/I");

  mp_output_tree->Branch("es_vertex_zpos", &m_espri_vertex_zpos,
                         "es_vertex_zpos/D");
  mp_output_tree->Branch("p_theta", &m_proton_theta, "p_theta/D");
  mp_output_tree->Branch("p_theta_eff", &m_proton_theta_eff, "p_theta_eff/D");
  mp_output_tree->Branch("p_phi", &m_proton_phi, "p_phi/D");
  mp_output_tree->Branch("nai_e", &m_nai_e, "nai_e/D");
  mp_output_tree->Branch("p_e", &m_proton_e, "p_e/D");
  mp_output_tree->Branch("p_e_eff", &m_proton_e_eff, "p_e/D");
  mp_output_tree->Branch("p_e_sim", &m_proton_e_sim, "p_e_sim/D");
  mp_output_tree->Branch("p_de", &m_proton_de, "p_de/D");
  // redundant value, not used, might be removed
  // mp_output_tree->Branch("p_de_eff", &m_proton_de_eff, "p_de_eff/D");
  mp_output_tree->Branch("p_de_sim", &m_proton_de_sim, "p_de_sim/D");
  mp_output_tree->Branch("p_tof", &m_proton_tof, "p_tof/D");

  // FDC2
  mp_output_tree->Branch("fdc2_xpos", &m_fdc2_xpos, "fdc2_xpos/D");
  mp_output_tree->Branch("fdc2_ypos", &m_fdc2_ypos, "fdc2_ypos/D");
  mp_output_tree->Branch("fdc2_aang", &m_fdc2_aang, "fdc2_aang/D");
  mp_output_tree->Branch("fdc2_bang", &m_fdc2_bang, "fdc2_bang/D");

  //HODF
  mp_output_tree->Branch("hodf_q", &m_hodf_q, "hodf_q[24]/D");
  mp_output_tree->Branch("hodf_t", &m_hodf_t, "hodf_t[24]/D");

  // theoretical values
  mp_output_tree->Branch("he_theta_theor", &m_he_theta_theor,
                         "he_theta_theor/D");
}

void RawTreeConverter::SetupInputTree() {
  assert(mp_raw_tree != nullptr);

  mp_raw_tree->SetBranchAddress("BRIPSPla_EventNumber", &m_event_no);
  mp_raw_tree->SetBranchAddress("TriggerBit", &m_triggers);

  mp_raw_tree->SetBranchAddress("F7Pla_QCal", &m_f7_q);
  mp_raw_tree->SetBranchAddress("F7Pla_TCal", &m_f7_t);
  mp_raw_tree->SetBranchAddress("SBT1Pla_QCal", &m_sbt1_q);
  mp_raw_tree->SetBranchAddress("SBT1Pla_TCal", &m_sbt1_t);
  mp_raw_tree->SetBranchAddress("SBT2Pla_QCal", &m_sbt2_q);
  mp_raw_tree->SetBranchAddress("SBT2Pla_TCal", &m_sbt2_t);
  mp_raw_tree->SetBranchAddress("SBVPla_QCal", &m_sbv_q);
  mp_raw_tree->SetBranchAddress("SBVPla_TCal", &m_sbv_t);

  mp_raw_tree->SetBranchAddress("BDC1_XPos", &m_bdc1_xpos);
  mp_raw_tree->SetBranchAddress("BDC1_YPos", &m_bdc1_ypos);
  mp_raw_tree->SetBranchAddress("BDC2_XPos", &m_bdc2_xpos);
  mp_raw_tree->SetBranchAddress("BDC2_YPos", &m_bdc2_ypos);

  mp_raw_tree->SetBranchAddress("FDC0_xPos", &m_fdc0_xpos);
  mp_raw_tree->SetBranchAddress("FDC0_yPos", &m_fdc0_ypos);
  mp_raw_tree->SetBranchAddress("S1MWDC_xPos", &m_s1dc_xpos);
  mp_raw_tree->SetBranchAddress("S1MWDC_yPos", &m_s1dc_ypos);

  mp_raw_tree->SetBranchAddress("ESPRI_RDCL_XPos", &m_espri_left_xpos);
  mp_raw_tree->SetBranchAddress("ESPRI_RDCL_YPos", &m_espri_left_ypos);
  mp_raw_tree->SetBranchAddress("ESPRI_RDCL_aAng", &m_espri_left_aang);
  mp_raw_tree->SetBranchAddress("ESPRI_RDCL_bAng", &m_espri_left_bang);
  mp_raw_tree->SetBranchAddress("ESPRI_RDCL_XNHitPlanes",
                                &m_espri_left_rdc_xnhitwires);
  mp_raw_tree->SetBranchAddress("ESPRI_RDCL_YNHitPlanes",
                                &m_espri_left_rdc_ynhitwires);
  mp_raw_tree->SetBranchAddress("ESPRIPlaL_TCal", &m_espri_left_tof);
  mp_raw_tree->SetBranchAddress("ESPRIPlaL_QCal", &m_espri_left_de);

  mp_raw_tree->SetBranchAddress("ESPRI_RDCR_XPos", &m_espri_right_xpos);
  mp_raw_tree->SetBranchAddress("ESPRI_RDCR_YPos", &m_espri_right_ypos);
  mp_raw_tree->SetBranchAddress("ESPRI_RDCR_aAng", &m_espri_right_aang);
  mp_raw_tree->SetBranchAddress("ESPRI_RDCR_bAng", &m_espri_right_bang);
  mp_raw_tree->SetBranchAddress("ESPRI_RDCR_XNHitPlanes",
                                &m_espri_right_rdc_xnhitwires);
  mp_raw_tree->SetBranchAddress("ESPRI_RDCR_YNHitPlanes",
                                &m_espri_right_rdc_ynhitwires);
  mp_raw_tree->SetBranchAddress("ESPRIPlaR_TCal", &m_espri_right_tof);
  mp_raw_tree->SetBranchAddress("ESPRIPlaR_QCal", &m_espri_right_de);

  mp_raw_tree->SetBranchAddress("ESPRI_NaIID", &m_raw_espri_nai_id);
  mp_raw_tree->SetBranchAddress("ESPRI_NaIPla_QCal", &m_raw_espri_nai_qcal);

  mp_raw_tree->SetBranchAddress("FDC2_XPos", &m_fdc2_xpos);
  mp_raw_tree->SetBranchAddress("FDC2_YPos", &m_fdc2_ypos);
  mp_raw_tree->SetBranchAddress("FDC2_aAng", &m_fdc2_aang);
  mp_raw_tree->SetBranchAddress("FDC2_bAng", &m_fdc2_bang);

  mp_raw_tree->SetBranchAddress("HodoPla_TCal", &m_hodf_t);
  mp_raw_tree->SetBranchAddress("HodoPla_QCal", &m_hodf_q);
}

void RawTreeConverter::ClearTreeVars() {
  m_event_no = INIT_VAL;
  for (int i = 0; i < 8; i++) {
    m_triggers[i] = false;
  }

  m_bdc1_xpos = INIT_VAL;
  m_bdc1_ypos = INIT_VAL;
  m_bdc2_xpos = INIT_VAL;
  m_bdc2_ypos = INIT_VAL;
  m_tgt_up_aang = INIT_VAL;
  m_tgt_up_bang = INIT_VAL;
  m_tgt_up_theta = INIT_VAL;
  m_tgt_up_phi = INIT_VAL;
  m_tgt_up_xpos = INIT_VAL;
  m_tgt_up_ypos = INIT_VAL;

  m_fdc0_xpos = INIT_VAL;
  m_fdc0_ypos = INIT_VAL;
  m_fdc0_aang = INIT_VAL;
  m_fdc0_bang = INIT_VAL;
  m_fdc0_theta = INIT_VAL;
  m_fdc0_phi = INIT_VAL;
  m_s1dc_xpos = INIT_VAL;
  m_s1dc_ypos = INIT_VAL;
  m_s1dc_aang = INIT_VAL;
  m_s1dc_bang = INIT_VAL;
  m_s1dc_theta = INIT_VAL;
  m_s1dc_phi = INIT_VAL;

  m_tgt_dwn_aang = INIT_VAL;
  m_tgt_dwn_bang = INIT_VAL;
  m_tgt_dwn_theta = INIT_VAL;
  m_tgt_dwn_phi = INIT_VAL;
  m_tgt_dwn_xpos = INIT_VAL;
  m_tgt_dwn_ypos = INIT_VAL;

  m_fdc2_xpos = INIT_VAL;
  m_fdc2_ypos = INIT_VAL;
  m_fdc2_aang = INIT_VAL;
  m_fdc2_bang = INIT_VAL;

  m_espri_left_xpos = INIT_VAL;
  m_espri_left_ypos = INIT_VAL;
  m_espri_left_aang = INIT_VAL;
  m_espri_left_bang = INIT_VAL;
  m_espri_left_aang_tgt = INIT_VAL;
  m_espri_left_bang_tgt = INIT_VAL;
  m_espri_left_de = INIT_VAL;
  m_espri_left_de_cal = INIT_VAL;
  m_espri_left_tof = INIT_VAL;
  for (int i = 0; i < ESPRI_NUM_NAIS; i++) {
    m_espri_left_nai_e[i] = INIT_VAL;
    m_espri_left_nai_e_cal[i] = INIT_VAL;
  }
  m_espri_left_e = INIT_VAL;
  m_espri_left_e_cal = INIT_VAL;
  m_espri_left_e_eff = INIT_VAL;
  m_espri_left_vertex_zpos = INIT_VAL;
  m_espri_left_proton_theta = INIT_VAL;
  m_espri_left_proton_phi = INIT_VAL;
  m_espri_left_rdc_nhitwires = INIT_VAL;
  m_espri_left_rdc_xnhitwires = INIT_VAL;
  m_espri_left_rdc_ynhitwires = INIT_VAL;

  m_espri_right_xpos = INIT_VAL;
  m_espri_right_ypos = INIT_VAL;
  m_espri_right_aang = INIT_VAL;
  m_espri_right_bang = INIT_VAL;
  m_espri_right_aang_tgt = INIT_VAL;
  m_espri_right_bang_tgt = INIT_VAL;
  m_espri_right_de = INIT_VAL;
  m_espri_right_de_cal = INIT_VAL;
  m_espri_right_tof = INIT_VAL;
  for (int i = 0; i < ESPRI_NUM_NAIS; i++) {
    m_espri_right_nai_e[i] = INIT_VAL;
    m_espri_right_nai_e_cal[i] = INIT_VAL;
  }
  m_espri_right_e = INIT_VAL;
  m_espri_right_e_cal = INIT_VAL;
  m_espri_right_e_eff = INIT_VAL;
  m_espri_right_vertex_zpos = INIT_VAL;
  m_espri_right_proton_theta = INIT_VAL;
  m_espri_right_proton_phi = INIT_VAL;
  m_espri_right_rdc_nhitwires = INIT_VAL;
  m_espri_right_rdc_xnhitwires = INIT_VAL;
  m_espri_right_rdc_ynhitwires = INIT_VAL;

  for (int i = 0; i < NAI_PMTS_NUM; i++) {
    m_raw_espri_nai_qcal[i] = INIT_VAL;
    m_raw_espri_nai_id[i] = INIT_VAL;
  }

  m_espri_vertex_zpos = INIT_VAL;
  m_proton_theta = INIT_VAL;
  m_proton_theta_eff = INIT_VAL;
  m_proton_phi = INIT_VAL;
  m_nai_e = INIT_VAL;
  m_proton_e = INIT_VAL;
  m_proton_e_eff = INIT_VAL;
  m_proton_de = INIT_VAL;
  m_proton_de_eff = INIT_VAL;
  m_proton_de_sim = INIT_VAL;
  m_proton_tof = INIT_VAL;

  for (int i = 0; i < 24; i++) {
    m_hodf_q[i] = INIT_VAL;
    m_hodf_t[i] = INIT_VAL;
  }

  m_he_theta_theor = INIT_VAL;
}

void RawTreeConverter::ComputeDerivedValues() {
  Double_t nai_pmts_e[NAI_PMTS_NUM];
  for (int i = 0; i < NAI_PMTS_NUM; i++) {
    int nai_index = m_raw_espri_nai_id[i] - 1;
    if (nai_index < 0) {
      // use counter instead of error message
      m_invalid_nai_data_events_num += 1;
      // std::cerr << "ERROR: Invalid NaI PMT index: " << nai_index << \
      //   ". Skipping this NaI PMT." << std::endl;
      continue;
    }
    Double_t nai_pmt_e = m_raw_espri_nai_qcal[i];
    nai_pmts_e[nai_index] = CheckAndFixNaIPmtValue(nai_pmt_e);
  }

  Double_t nais_e[NAIS_NUM];
  for (int i = 0; i < NAIS_NUM; i++) {
    // geometrical mean, does not depend on hit position
    nais_e[i] = std::sqrt(nai_pmts_e[2 * i] * nai_pmts_e[2 * i + 1]);
  }

  // NaIs 0 - 7: left ESPRI; 7 - 14: right ESPRI
  for (int i = 0; i < ESPRI_NUM_NAIS; i++) {
    m_espri_left_nai_e[i] = nais_e[i];
    m_espri_right_nai_e[i] = nais_e[i + 7];
  }

  // e of whole detector - max deposit in any crystal
  m_espri_left_e = *std::max_element(m_espri_left_nai_e,
                                     m_espri_left_nai_e + NAIS_NUM);
  m_espri_right_e = *std::max_element(m_espri_right_nai_e,
                                      m_espri_right_nai_e + NAIS_NUM);

  // correct value definitions
  // data of Y axes of DSDCs were inverted by some reason
  m_s1dc_ypos = -m_s1dc_ypos;
  m_fdc0_ypos = -m_fdc0_ypos;

  m_espri_left_rdc_nhitwires = std::min(m_espri_left_rdc_xnhitwires,
                                        m_espri_left_rdc_ynhitwires);
  m_espri_right_rdc_nhitwires = std::min(m_espri_right_rdc_xnhitwires,
                                         m_espri_right_rdc_ynhitwires);

  ComputeBdcsAngles();
  ComputeProtonAngles();
  ComputeProtonEffAngles();
  ComputeVertexFromRdc('l', m_espri_left_xpos, m_espri_left_ypos,
                       m_espri_left_aang, m_espri_left_bang);
  ComputeVertexFromRdc('r', m_espri_right_xpos, m_espri_right_ypos,
                       m_espri_right_aang, m_espri_right_bang);
  ComputeCalEnergies();
  ComputeSimEnergies();

  // KLUDGE: Used to get proton energies into separate .csv file.
  // ComputeProtonEWithMtrlELossesCsv();
}

double RawTreeConverter::CheckAndFixNaIPmtValue(double nai_pmt_value) {
  if (nai_pmt_value < 0) {
    m_invalid_nai_pmt_values_num += 1;
    return 0;
  } else if (nai_pmt_value != nai_pmt_value) {
    std::cout << "[WARNING] NaI PMT value is NaN!";
    return 0;
  } else {
    return nai_pmt_value;
  }
}

/*
  NOTE: phi angle defined from Y axis:
  phi = -(vector.Phi() + 90); vector.Phi() - as usual from X axis
*/
void RawTreeConverter::ComputeBdcsAngles() {
  TVector3 fragment_vec;
  TVector3 beam_vec;
  TVector3 zaxis_vec(0,0,1);

  // if x or y pos in range - compute angle
  if (std::abs(m_bdc1_xpos) < BDC_WIDTH/2 &&
      std::abs(m_bdc1_ypos) < BDC_HEIGHT/2 &&
      std::abs(m_bdc2_xpos) < BDC_WIDTH/2 &&
      std::abs(m_bdc2_ypos) < BDC_HEIGHT/2) {
    m_tgt_up_aang = (m_bdc2_xpos - m_bdc1_xpos) / DIST_BDCs;
    m_tgt_up_bang = (m_bdc2_ypos - m_bdc1_ypos) / DIST_BDCs;
    beam_vec.SetXYZ(m_tgt_up_aang, m_tgt_up_bang, 1);
    m_tgt_up_theta = beam_vec.Angle(zaxis_vec);
    m_tgt_up_phi = beam_vec.Phi() - (90 * D2R);
    if (m_tgt_up_phi < -180 * D2R) {
      m_tgt_up_phi += 360 * D2R;
    }

    m_tgt_up_xpos = m_bdc2_xpos + \
      m_tgt_up_aang * (DIST_BDC2_TARGET + TARGET_Z / 2.);
    m_tgt_up_ypos = m_bdc2_ypos + \
      m_tgt_up_bang * (DIST_BDC2_TARGET + TARGET_Z / 2.);
  }

  // if x or y pos in range - compute angle
  if (std::abs(m_fdc0_xpos) < FDC0_WIDTH/2 &&
      std::abs(m_fdc0_ypos) < FDC0_HEIGHT/2 &&
      std::abs(m_s1dc_xpos) < S1DC_WIDTH/2 &&
      std::abs(m_s1dc_ypos) < S1DC_HEIGHT/2) {
    m_tgt_dwn_aang = (m_s1dc_xpos - m_fdc0_xpos) / DIST_S1MWDC_FDC0;
    m_tgt_dwn_bang = (m_s1dc_ypos - m_fdc0_ypos) / DIST_S1MWDC_FDC0;
    fragment_vec.SetXYZ(m_tgt_dwn_aang, m_tgt_dwn_bang, 1);
    m_tgt_dwn_theta = fragment_vec.Angle(beam_vec);
    m_tgt_dwn_phi = fragment_vec.Phi() - (90 * D2R);
    if (m_tgt_dwn_phi < -180 * D2R) {
      m_tgt_dwn_phi += 360 * D2R;
    }

    m_tgt_dwn_xpos = m_fdc0_xpos - \
      m_tgt_dwn_aang * (DIST_FDC0_TARGET + TARGET_Z / 2.);
    m_tgt_dwn_ypos = m_fdc0_ypos - \
      m_tgt_dwn_bang * (DIST_FDC0_TARGET + TARGET_Z / 2.);
  }

  // if x or y pos in range - compute angle
  if (std::abs(m_fdc0_xpos) < FDC0_WIDTH/2 &&
      std::abs(m_fdc0_ypos) < FDC0_HEIGHT/2 &&
      std::abs(m_tgt_up_xpos) < 9e3 &&
      std::abs(m_tgt_up_ypos) < 9e3) {
    m_fdc0_aang = (m_fdc0_xpos - m_tgt_up_xpos) / DIST_FDC0_TARGET;
    m_fdc0_bang = (m_fdc0_ypos - m_tgt_up_ypos) / DIST_FDC0_TARGET;
    fragment_vec.SetXYZ(m_fdc0_aang, m_fdc0_bang, 1);
    m_fdc0_theta = fragment_vec.Angle(beam_vec);
    m_fdc0_phi = fragment_vec.Phi() - (90 * D2R);
    if (m_fdc0_phi < -180 * D2R) {
      m_fdc0_phi += 360 * D2R;
    }
  }

  // if x or y pos in range compute angle
  if (std::abs(m_s1dc_xpos) < S1DC_WIDTH/2 &&
      std::abs(m_s1dc_ypos) < S1DC_HEIGHT/2 &&
      std::abs(m_tgt_up_xpos) < 9e3 &&
      std::abs(m_tgt_up_ypos) < 9e3) {
    m_s1dc_aang = (m_s1dc_xpos - m_tgt_up_xpos) / \
      (DIST_FDC0_TARGET + DIST_S1MWDC_FDC0);
    m_s1dc_bang = (m_s1dc_ypos - m_tgt_up_ypos) / \
      (DIST_FDC0_TARGET + DIST_S1MWDC_FDC0);
    fragment_vec.SetXYZ(m_s1dc_aang, m_s1dc_bang, 1);
    m_s1dc_theta = fragment_vec.Angle(beam_vec);
    m_s1dc_phi = fragment_vec.Phi() - (90 * D2R);
    if (m_s1dc_phi < -180 * D2R) {
      m_s1dc_phi += 360 * D2R;
    }
  }

  // proton a & b angles
  m_espri_left_aang_tgt = (m_espri_left_xpos - m_tgt_up_xpos) / \
    DIST_RDC_TARGET;
  m_espri_left_bang_tgt = (m_espri_left_ypos - m_tgt_up_ypos) / \
    DIST_RDC_TARGET;
  m_espri_right_aang_tgt = (m_espri_right_xpos - m_tgt_up_xpos) / \
    DIST_RDC_TARGET;
  m_espri_right_bang_tgt = (m_espri_right_ypos - m_tgt_up_ypos) / \
    DIST_RDC_TARGET;
}

/*
  NOTE: phi angle defined non-standard for polar coordinates
  from Y axis - CCW rotation: positive; CW: negative
*/
void RawTreeConverter::ComputeProtonAngles() {
  TVector3 beam(m_tgt_up_aang, m_tgt_up_bang, 1.);
  TVector3 recoil_p;

  // check if BDC data is well defined as without it
  // we cannot determine proton angles
  if (std::abs(m_tgt_up_xpos) > 9e3 ||
      std::abs(m_tgt_up_ypos) > 9e3) {
    return;
  }

  bool proton_data_error = false;
  if (std::abs(m_espri_left_xpos) < ESPRI_RDC_WIDTH/2 &&
      std::abs(m_espri_left_ypos) < ESPRI_RDC_WIDTH/2 &&
      std::abs(m_espri_right_xpos) < ESPRI_RDC_WIDTH/2 &&
      std::abs(m_espri_right_ypos) < ESPRI_RDC_WIDTH/2) {
    // std::cerr << "ERROR: Proton angular data in both RDCs!" << std::endl;
    proton_data_error = true;
  } else if (std::abs(m_espri_left_xpos) < ESPRI_RDC_WIDTH/2 &&
             std::abs(m_espri_left_ypos) < ESPRI_RDC_WIDTH/2) {
    recoil_p.SetXYZ(m_espri_left_xpos - m_tgt_up_xpos,
                    m_espri_left_ypos - m_tgt_up_ypos,
                    0);
    recoil_p.RotateY(180 * D2R);
    recoil_p.SetZ(DIST_RDC_TARGET);
    recoil_p.RotateY(RDC_CENTER_ANGLE * D2R);
    m_espri_left_proton_theta = recoil_p.Angle(beam);
    m_espri_left_proton_phi = recoil_p.Phi() - (90 * D2R);
    if (m_espri_left_proton_phi < -180 * D2R) {
      m_espri_left_proton_phi += 360 * D2R;
    }

    m_good_proton_events_num += 1;
    m_proton_theta = m_espri_left_proton_theta;
    m_proton_phi = m_espri_left_proton_phi;
  } else if (std::abs(m_espri_right_xpos) < ESPRI_RDC_WIDTH/2 &&
             std::abs(m_espri_right_ypos) < ESPRI_RDC_WIDTH/2) {
    recoil_p.SetXYZ(m_espri_right_xpos - m_tgt_up_xpos,
                    m_espri_right_ypos - m_tgt_up_ypos,
                    0);
    recoil_p.RotateY(180 * D2R);
    recoil_p.SetZ(DIST_RDC_TARGET);
    recoil_p.RotateY(-RDC_CENTER_ANGLE * D2R);
    m_espri_right_proton_theta = recoil_p.Angle(beam);
    m_espri_right_proton_phi = recoil_p.Phi() - (90 * D2R);
    if (m_espri_right_proton_phi < -180 * D2R) {
      m_espri_right_proton_phi += 360 * D2R;
    }

    m_good_proton_events_num += 1;
    m_proton_theta = m_espri_right_proton_theta;
    m_proton_phi = m_espri_right_proton_phi;
  } else {
    // std::cerr << "ERROR: No proton angular data for the event." <<  \
    //     std::endl;
    proton_data_error = true;
  }

  if (proton_data_error) {
    m_bad_proton_events_num += 1;
  }
}

void
RawTreeConverter::ComputeVertexFromRdc(char espri_side, Double_t xpos,
                                       Double_t ypos, Double_t aang,
                                       Double_t bang) {
  assert(espri_side == 'l' || espri_side == 'r');
  if (std::abs(xpos) > ESPRI_ACCEPTANCE_WIDTH / 2. ||
      std::abs(ypos) > ESPRI_ACCEPTANCE_WIDTH / 2.) {
    return;
  }
  // the angle is too big to converge on target - it will miss
  if (std::abs(aang) > 700) {
    return;
  }

  // DEBUG
  // if (xpos < 0) {
  //     return;
  // }
  // espri_side = 'l';
  // xpos = -10;
  // aang = -10;

  Double_t frame_rotation_angle = (180 + RDC_CENTER_ANGLE) * D2R;
  TVector3 rdc_dir_vec(0.0, 0.0, 1.0);
  TVector3 rdc_hit_pos(xpos, 0.0, 0.0);

  // std::cout << "\nRDC hit position: " << \
  //     rdc_hit_pos.X() << "; " << rdc_hit_pos.Y() << \
  //     "; " << rdc_hit_pos.Z() << std::endl;
  // std::cout << "RDC a angle: " << \
  //     aang << "; " << std::endl;

  // set proton hit direction vector from RDC data
  // angle direction inverted as ESPRI uses LHS roordinates
  // see analysis notes for definition
  rdc_dir_vec.RotateY(-aang * 1e-3);

  // std::cout << "RDC dir vec: " << \
  //     rdc_dir_vec.X() << "; " << rdc_dir_vec.Y() << \
  //     "; " << rdc_dir_vec.Z() << std::endl;

  // rotate RDC frame to match target frame axes
  if (espri_side == 'r') {
    frame_rotation_angle *= -1;
  }
  rdc_hit_pos.RotateY(frame_rotation_angle);
  rdc_dir_vec.RotateY(frame_rotation_angle);

  // translate RDC frame to match target frame origin
  Double_t dx = DIST_RDC_TARGET * TMath::Sin(RDC_CENTER_ANGLE * D2R);
  if (espri_side == 'r') {
    dx *= -1;
  }
  Double_t dz = DIST_RDC_TARGET * TMath::Cos(RDC_CENTER_ANGLE * D2R);
  TVector3 coords_translation(dx, 0.0, dz);
  rdc_hit_pos += coords_translation;
  Double_t hit_ray_length = std::abs(rdc_hit_pos.X() / rdc_dir_vec.X());
  TVector3 tgt_hit_pos = rdc_hit_pos + hit_ray_length * rdc_dir_vec;

  // std::cout << "RDC hit position (aft. translation): " << \
  //     rdc_hit_pos.X() << "; " << rdc_hit_pos.Y() << \
  //     "; " << rdc_hit_pos.Z() << std::endl;
  // std::cout << "RDC dir vec (aft. rotation): " << \
  //     rdc_dir_vec.X() << "; " << rdc_dir_vec.Y() << \
  //     "; " << rdc_dir_vec.Z() << std::endl;
  // std::cout << "Hit ray length: " << hit_ray_length << std::endl;
  // std::cout << "Target hit coorinates: " << \
  //     tgt_hit_pos.X() << "; " << tgt_hit_pos.Y() << \
  //     "; " << tgt_hit_pos.Z() << std::endl;

  // bad event if X coordinate not on the beamline
  if (std::abs(tgt_hit_pos.X()) > 3) {
    std::cout <<
      "ERROR: X coodrinate not on the beamline - " <<
      "reconstruction failed." << std::endl;
    return;
  }

  // assign values if extrapolation was successfull
  if (espri_side == 'l') {
    m_espri_left_vertex_zpos = tgt_hit_pos.Z();
  } else {
    m_espri_right_vertex_zpos = tgt_hit_pos.Z();
  }
  m_espri_vertex_zpos = tgt_hit_pos.Z();

  // DEBUG: check if Z is in near target body
  // if (std::abs(tgt_hit_pos.Z()) < 2) {
  //     std::cout << "GOOD EVENT!!!" << std::endl;
  // }
}

void RawTreeConverter::ComputeCalEnergies() {
  for (int i = 0; i < ESPRI_NUM_NAIS; ++i) {
    m_espri_left_nai_e_cal[i] = m_espri_left_nai_cal_consts[i] *
        (m_espri_left_nai_e[i] - m_espri_left_nai_pedestals[i]);
    m_espri_right_nai_e_cal[i] = m_espri_right_nai_cal_consts[i] *
        (m_espri_right_nai_e[i] - m_espri_right_nai_pedestals[i]);
  }

  m_espri_left_de_cal =  m_espri_left_pde_cal_const *
      (m_espri_left_de - m_espri_left_pde_pedestal);
  m_espri_right_de_cal = m_espri_right_pde_cal_const *
      (m_espri_right_de - m_espri_right_pde_pedestal);

  FindEspriCalEnergies2();

  // make sure we know to which of ESPRIs this event belongs
  // if not, do not compute values derived from ESPRI L/R
  if (IsEspriLeftEvent() || IsEspriRightEvent()) {

    if (IsEspriLeftEvent()) {
      m_nai_e = m_espri_left_e_cal;
      m_proton_de = m_espri_left_de_cal;
      m_proton_tof = m_espri_left_tof;
    } else if (IsEspriRightEvent()) {
      m_nai_e = m_espri_right_e_cal;
      m_proton_de = m_espri_right_de_cal;
      m_proton_tof = m_espri_right_tof;
    }
    if (IsWellDefinedValue(m_nai_e) && IsWellDefinedValue(m_proton_de)) {
      m_proton_e = m_nai_e + m_proton_de;
    }
  }

  // std::cout << "Cal. NaI E: " << m_nai_e << std::endl;
  if (std::isnan(m_nai_e)) {
    std::cout << "WARNING: Proton E has NaN value!" << std::endl;
  }
  if (std::isnan(m_proton_de)) {
    std::cout << "WARNING: Proton delta-E has NaN value!" << std::endl;
  }
  // if proton energy is "well_defined" then compute effective energies
  ComputeEffectiveEnergies();
}

void RawTreeConverter::ComputeEffectiveEnergies() {
  // Do not compute anything if proton E is not well-defined
  if (!IsWellDefinedValue(m_proton_e)) {
    return;
  }

  Double_t edeg_thickness = ComputeEdegThickness(m_proton_theta);
  Double_t proton_edeg_de_computed =
    ComputeProtonDe(m_proton_e, m_edeg_stp_power, edeg_thickness, true);

  // restore e-loss in degrader
  m_proton_e_eff = m_proton_e + proton_edeg_de_computed;
  if (IsEspriLeftEvent()) {
    m_espri_left_e_eff = m_proton_e_eff;
  } else if (IsEspriRightEvent()){
    m_espri_right_e_eff = m_proton_e_eff;
  }
}

void RawTreeConverter::FindEspriCalEnergies1() {
  auto espri_left_e_cal = std::max_element(m_espri_left_nai_e_cal,
                                           m_espri_left_nai_e_cal + ESPRI_NUM_NAIS);
  auto espri_right_e_cal = std::max_element(m_espri_right_nai_e_cal,
                                            m_espri_right_nai_e_cal + ESPRI_NUM_NAIS);
  m_espri_left_e_cal = *espri_left_e_cal;
  m_espri_right_e_cal = *espri_right_e_cal;
}

int RawTreeConverter::FindMaxENaiId() {
  double* naie_cal = nullptr;
  if (IsEspriLeftEvent()) {
    naie_cal = m_espri_left_nai_e_cal;
  } else if (IsEspriRightEvent()) {
    naie_cal = m_espri_right_nai_e_cal;
  } else {
    // if we cannot determine if it was ESL or ESR event, then, we cannot
    // if we cannot
    return 0;
  }
  int max_e_nai_id = -1;
  double max_e_cal = 0;
  for (int i = 0; i < 7; i++) {
    if (naie_cal[i] > max_e_cal) {
      max_e_cal = naie_cal[i];
      max_e_nai_id = i;
    }
  }
  return max_e_nai_id;
}
void RawTreeConverter::FindEspriCalEnergies2() {
  TVector3 position, direction;
  if (IsEspriLeftEvent()) {
    position = TVector3(m_espri_left_xpos, m_espri_left_ypos, 0);
    direction = TVector3(m_espri_left_aang / 1e3, m_espri_left_bang / 1e3, 1);
  } else if (IsEspriRightEvent()) {
    position = TVector3(m_espri_right_xpos, m_espri_right_ypos, 0);
    direction = TVector3(m_espri_right_aang / 1e3, m_espri_right_bang / 1e3, 1);
  } else {
    // if we cannot determine if it was ESL or ESR event, then, we cannot
    // determine energy for this event
    return;
  }

  TVector3 nai_hit_pos = m_nais_array.find_hit_pos(position, direction);
  int nai_id = m_nais_array.find_hit_nai_id(nai_hit_pos);
  int max_e_nai_id = FindMaxENaiId();

  if (nai_id == -1) {
    // double nearest_nai_dist = m_nais_array.find_nearest_nai_dist(position, direction);
    nai_id = m_nais_array.find_nearest_nai_id(position, direction);
  }
  if (nai_id != -1 && max_e_nai_id != nai_id) {
    int nais_diff = std::abs(nai_id - max_e_nai_id);
    if (nais_diff < 2) {
      nai_id = max_e_nai_id;
    }
  }

  if (nai_id != -1) {
    if (IsEspriLeftEvent()) {
      m_espri_left_e_cal = m_espri_left_nai_e_cal[nai_id];
    } else if (IsEspriRightEvent()) {
      m_espri_right_e_cal = m_espri_right_nai_e_cal[nai_id];
    }
  }
}

void RawTreeConverter::ComputeProtonEffAngles() {
  if (m_proton_theta * 1/D2R < 50.0 || m_proton_theta * 1/D2R > 73.0) {
    return;
  }

  // std::cout << "Computing effective proton angle for measured: "
  //           << m_proton_theta * 1/D2R << "..." << std::endl;

  if (IsHe4Run()) {
    if (IsEspriLeftEvent()) {
      m_proton_theta_eff = m_p_he4_esl_theta_eff_vs_meas->Eval(m_proton_theta);
    } else if (IsEspriRightEvent()) {
      m_proton_theta_eff = m_p_he4_esr_theta_eff_vs_meas->Eval(m_proton_theta);
    }
  }

  if (IsHe6Run() || IsCarbonRun()) {
    if (IsEspriLeftEvent()) {
      m_proton_theta_eff = m_p_he6_esl_theta_eff_vs_meas->Eval(m_proton_theta);
    } else if (IsEspriRightEvent()) {
      m_proton_theta_eff = m_p_he6_esr_theta_eff_vs_meas->Eval(m_proton_theta);
    }
  }

  assert(!std::isnan(m_proton_theta_eff) && "Effective proton angle interpolation failed!");
}

Double_t
RawTreeConverter::ComputeProtonMtrlEloss() {
  if (m_proton_theta_eff * 1/D2R < 51.7 || m_proton_theta_eff * 1/D2R > 73.15) {
    throw std::invalid_argument("Proton theta is out of range for materials e-loss computation!");
  }

  // std::cout << "Computing proton e-loss in materials for measured theta: "
  //           << m_proton_theta * 1/D2R << "..." << std::endl;

  Double_t e_loss = 0;
  if (IsHe4Run()) {
    if (IsEspriLeftEvent()) {
      e_loss = m_p_he4_esl_theta_eff_vs_mtrls_e_loss->Eval(m_proton_theta_eff);
    } else if (IsEspriRightEvent()) {
      e_loss = m_p_he4_esr_theta_eff_vs_mtrls_e_loss->Eval(m_proton_theta_eff);
    }
  }

  if (IsHe6Run() || IsCarbonRun()) {
    if (IsEspriLeftEvent()) {
      e_loss = m_p_he6_esl_theta_eff_vs_mtrls_e_loss->Eval(m_proton_theta_eff);
    } else if (IsEspriRightEvent()) {
      e_loss = m_p_he6_esr_theta_eff_vs_mtrls_e_loss->Eval(m_proton_theta_eff);
    }
  }

  return e_loss;
}

void RawTreeConverter::ComputeProtonEWithMtrlELossesCsv() {
  std::ofstream ofs_esr("out/p_he4_esr_p_e_with_e_losses.csv");
  ofs_esr << "#p_theta_rdc,p_e (data for NaIs cal.)" << std::endl;
  for (m_proton_theta = 50.0 * D2R; m_proton_theta < 73 * D2R;
       m_proton_theta += 0.1 * D2R) {
    m_espri_right_de = 240;
    ComputeProtonEffAngles();
    // std::cout << "proton theta & theta_eff: " << m_proton_theta / D2R
    //           << "; " << m_proton_theta_eff / D2R << std::endl;
    ComputeSimEnergies();
    ofs_esr << m_proton_theta / D2R << "," << m_proton_e_sim << std::endl;
  }

  std::ofstream ofs_esl("out/p_he4_esl_p_e_with_e_losses.csv");
  ofs_esl << "#p_theta_rdc,p_e (data for NaIs cal.)" << std::endl;
  for (m_proton_theta = 50.0 * D2R; m_proton_theta < 73 * D2R;
       m_proton_theta += 0.1 * D2R) {
    m_espri_left_de = 240;
    ComputeProtonEffAngles();
    // std::cout << "proton theta & theta_eff: " << m_proton_theta / D2R
    //           << "; " << m_proton_theta_eff / D2R << std::endl;
    ComputeSimEnergies();
    ofs_esl << m_proton_theta / D2R << "," << m_proton_e_sim << std::endl;
  }

  std::exit(0);
}

void RawTreeConverter::ComputeSimEnergies() {
  if (!IsWellDefinedValue(m_proton_theta_eff)) {
    return;
  }
  if (m_proton_theta_eff * 1/D2R < 51.7 || m_proton_theta_eff * 1/D2R > 73.15) {
    // std::cout << "WARN: Out of range proton angle: " << m_proton_theta_eff * 1/D2R
    //           << "!" << std::endl;
    return;
  }

  Double_t proton_e_kin = 0;
  if (IsHe4Run()) {
    proton_e_kin = m_he4_sim_p_e->Eval(m_proton_theta_eff);
    m_proton_de_sim = m_he4_sim_p_de->Eval(m_proton_theta_eff);
  } else if (IsHe6Run() || IsCarbonRun()) {
    proton_e_kin = m_he6_sim_p_e->Eval(m_proton_theta_eff);
    m_proton_de_sim = m_he6_sim_p_de->Eval(m_proton_theta_eff);
  } else {
    throw std::invalid_argument("Invalid run type for proton energy computation!");
  }

  Double_t edeg_thickness = ComputeEdegThickness(m_proton_theta);
  Double_t mtrls_e_loss = ComputeProtonMtrlEloss();
  Double_t edeg_e_loss =
    ComputeProtonDe(proton_e_kin - mtrls_e_loss, m_edeg_stp_power,
                    edeg_thickness, false);
  m_proton_e_sim = proton_e_kin - mtrls_e_loss - edeg_e_loss;

  // std::cout << "Computed simulated proton E for measured theta: "
  //           << m_proton_theta * 1./D2R << std::endl;
  // std::cout << "\tProton E_kin: " << proton_e_kin << std::endl;
  // std::cout << "\tProton mtrls e_loss: " << mtrls_e_loss
  //           << " MeV." << std::endl;
  // std::cout << "\tE-deg thickness: " << edeg_thickness << " mm;"
  //           << " E-loss: " << edeg_e_loss << " MeV. " << std::endl;
  // std::cout << "\tProton E_sim: " << m_proton_e_sim << std::endl;
  // std::cout << m_proton_theta_eff/D2R << "," << m_proton_e_sim << std::endl;
}

void RawTreeConverter::ComputeTheoreticalValues() {
  if (m_proton_theta_eff == INIT_VAL) {
    return;
  }
  if (m_proton_theta_eff < 45 * D2R || m_proton_theta_eff > 80 * D2R) {
    // std::cout << "WARNING: Unexpected proton angle: " <<  \
    //   m_proton_theta * 1/D2R << std::endl;
    return;
  }

  // std::cerr << "computing theoretical he4/6 angle given proton angle: " << \
  //   m_proton_theta * 1/D2R << std::endl;
  if (IsHe4Run()) {
    m_he_theta_theor = m_he4_kinematics.GetFragmentAngle(m_proton_theta_eff);
  } else if (IsHe6Run() || IsCarbonRun()) {
    m_he_theta_theor = m_he6_kinematics.GetFragmentAngle(m_proton_theta_eff);
  }
}

ScatteringKinematics
RawTreeConverter::LoadKinematicsData(std::string filename) {
  std::cout << "Loading kinematics data file: " << filename << "..." << std::endl;
  std::ifstream ifile(filename);

  ScatteringKinematics kin_data;
  for (CSVIterator iter(ifile, 1); iter != CSVIterator(); iter++) {
    auto csv_row = *iter;
    assert(csv_row.size() == 2);

    Double_t p_angle = std::stod(csv_row[0]) * D2R;
    Double_t he_angle = std::stod(csv_row[1]) * D2R;
    kin_data.AddPoint(he_angle, p_angle);
  }

  return kin_data;
}

std::vector<Double_t>
RawTreeConverter::LoadEspriNaiCalConsts(std::string filename, bool load_pedestals) {
  std::ifstream ifile(filename);

  std::vector<Double_t> cal_consts;
  std::cout << "Loading calib. consts for NaIs from : " << filename << "..." << std::endl;
  for (CSVIterator iter(ifile, 1); iter != CSVIterator(); iter++) {
    auto csv_row = *iter;
    assert(csv_row.size() == 3 && "Invalid format of input data");

    // Double_t nai_id = std::stod(csv_row[0]);
    Double_t nai_ped = std::stod(csv_row[1]);
    Double_t nai_const = std::stod(csv_row[2]);
    //std::cout << "NaI #" << nai_id << " : " << nai_const << std::endl;
    cal_consts.push_back(load_pedestals ? nai_ped : nai_const);
  }

  assert(cal_consts.size() == 7 && "Invalid number of NaI consts. Should be 7.");
  return cal_consts;
}

std::vector<Double_t>
RawTreeConverter::LoadEspriPdeCalConsts(std::string filename, bool load_pedestals) {
  std::ifstream ifile(filename);

  std::vector<Double_t> cal_consts;
  std::cout << "Loading calib. consts for ESPRI pdEs from: " << filename << "..." << std::endl;
  for (CSVIterator iter(ifile, 0); iter != CSVIterator(); iter++) {
    auto csv_row = *iter;
    assert(csv_row.size() == 3 && "Invalid format of input data");

    std::string pde_id = csv_row[0];
    Double_t pde_pedestal = std::stod(csv_row[1]);
    Double_t pde_const = std::stod(csv_row[2]);
    //std::cout << "pdE id: " << pde_id << " : " << pde_const << std::endl;
    cal_consts.push_back(load_pedestals ? pde_pedestal : pde_const);
  }

  assert(cal_consts.size() == 2 && "Invalid number of pdE consts. Should be 2 (2 peds or 2 cal. facts)");
  return cal_consts;
}

std::unique_ptr<ROOT::Math::Interpolator>
RawTreeConverter::LoadProtonKinE(std::string filename) {
  std::ifstream ifile(filename);
  std::vector<Double_t> p_angles, p_energies;

  for (CSVIterator iter(ifile, 2); iter != CSVIterator(); iter++) {
    auto csv_row = *iter;
    assert(csv_row.size() == 4 && "Invalid input file format");

    Double_t p_a = std::stod(csv_row[2]) * D2R;
    Double_t p_e = std::stod(csv_row[3]);

    p_angles.push_back(p_a);
    p_energies.push_back(p_e);
  }

  auto interp = new ROOT::Math::Interpolator{
    static_cast<unsigned int>(p_angles.size()),
    ROOT::Math::Interpolation::kLINEAR};
  interp->SetData(p_angles, p_energies);

  using interp_ptr_t = std::unique_ptr<ROOT::Math::Interpolator>;
  return interp_ptr_t(interp);
}

std::unique_ptr<ROOT::Math::Interpolator>
RawTreeConverter::LoadStoppingPower(std::string filename) {
  std::cout << "Loading stopping power file: " << filename << "..." << std::endl;
  std::ifstream ifile(filename);
  std::vector<Double_t> p_energies, stopping_powers;

  for (CSVIterator iter(ifile, 2); iter != CSVIterator(); iter++) {
    auto csv_row = *iter;
    assert(csv_row.size() == 4 && "Invalid input file format");

    Double_t p_e = std::stod(csv_row[0]);
    Double_t stp = std::stod(csv_row[3]);

    p_energies.push_back(p_e);
    stopping_powers.push_back(stp);
  }

  auto interp = new ROOT::Math::Interpolator{
    static_cast<unsigned int>(p_energies.size()),
    ROOT::Math::Interpolation::kLINEAR};
  interp->SetData(p_energies, stopping_powers);

  using interp_ptr_t = std::unique_ptr<ROOT::Math::Interpolator>;
  return interp_ptr_t(interp);
}

Double_t
RawTreeConverter::ComputeProtonDe(Double_t proton_e,
                                  std::unique_ptr<ROOT::Math::Interpolator>& stp_power,
                                  Double_t material_thickness, bool inverse) {
  static const Double_t dx_step = 0.05;

  // check for invalid number
  if (std::isnan(proton_e)) {
    throw std::invalid_argument("Input energy value passed to ComputeProtonDe() is NaN!");
  }
  // non-sensical value cut off
  if (proton_e > 399.9) {
    std::cout << "WARNING: Unusually high input E in ComputeProtonDe(): "
              << proton_e << std::endl;
    return 0;
  }
  // very low e - just return loss of its energy
  if (proton_e < 1.1) {
    return proton_e;
  }

  Double_t p_de = 0;
  Double_t p_e_cur = proton_e;
  for (Double_t z = 0; z < material_thickness; z += dx_step) {
    // if energy is lower than 1 MeV, consider particle stopped in the material
    if (p_e_cur < 1.1) {
      return proton_e;
    }

    assert((p_e_cur > 1.1 && p_e_cur < 399.9) && "Domain error!");
    // if (!(p_e_cur > 1.1 && p_e_cur < 399.9)) {
    //   m_logger.debug("p_e_cur: %.1f", p_e_cur);
    // }
    // std::cout << "Evaluating stopping power with p_e_cur: "
    //           << p_e_cur << std::endl;
    Double_t de = dx_step * stp_power->Eval(p_e_cur);
    //Double_t de = 0;
    p_de += de;
    if (!inverse) { // direct - computing e-loss
      p_e_cur -= de;
    } else { // inverse - restoring e_loss
      p_e_cur += de;
    }
  }

  return p_de;
}

/*
  Returns two interpolator objects to obtain dE & E as a function of proton angle.

  Params:
    he_kin_file - string specifying kin. file for 6 or 4He
*/
std::pair<
  std::shared_ptr<ROOT::Math::Interpolator>,
  std::shared_ptr<ROOT::Math::Interpolator> >
RawTreeConverter::ComputeRefHeDeE(std::string proton_kin_file) {
  std::cout << "Loading proton E kinematics from file: " << proton_kin_file <<
    "..." << std::endl;
  auto p_energies = LoadProtonKinE(proton_kin_file);

  std::vector<Double_t> angles, des, es;
  for (Double_t p_angle = 50; p_angle < 75; p_angle += 0.1) {
    Double_t angle_rad = p_angle * D2R;
    Double_t p_e_kin = p_energies->Eval(angle_rad);
    Double_t p_de = ComputeProtonDe(p_e_kin, m_pla_stp_power, ESPRI_PLASTIC_THICKNESS);
    angles.push_back(angle_rad);
    des.push_back(p_de);
    es.push_back(p_e_kin);
    // std::cout << "angle: " << p_angle << " dE: " << p_de <<
    //   " E: " << p_e << std::endl;
  }

  using interp_ptr_t = std::shared_ptr<ROOT::Math::Interpolator>;

  auto interp_de = new ROOT::Math::Interpolator{
    static_cast<unsigned int>(angles.size()),
    ROOT::Math::Interpolation::kLINEAR};
  interp_de->SetData(angles, des);
  auto p_interp_de = interp_ptr_t(interp_de);

  auto interp_e = new ROOT::Math::Interpolator{
    static_cast<unsigned int>(angles.size()),
    ROOT::Math::Interpolation::kLINEAR};
  interp_e->SetData(angles, es);
  auto p_interp_e = interp_ptr_t(interp_e);

  return make_pair(p_interp_de, p_interp_e);
}

Double_t RawTreeConverter::ComputeEdegThickness(Double_t proton_theta) {
  static const Double_t edeg_dist_from_target = 1123;
  static const Double_t edeg_base_angle = 53 * D2R;
  // static const Double_t edeg_finish_angle = 65 * D2R;
  static const Double_t edeg_thickness = 29;
  static const Double_t edeg_length = 255;
  static const Double_t edeg_min_thickness = 2;
  static const Double_t edeg_angle_tan = edeg_thickness / edeg_length;

  if (proton_theta < edeg_base_angle) {
    return 0;
  }

  Double_t dist_from_base = 2 * edeg_dist_from_target *
    TMath::Tan((proton_theta - edeg_base_angle) / 2);
  Double_t length2 = edeg_length - dist_from_base;
  Double_t thickness = length2 * edeg_angle_tan;
  if (thickness < edeg_min_thickness) {
    thickness = 0;
  }
  // std::cout << "Proton angle: " << proton_theta / D2R <<
  //   "; edeg thickness: " << thickness << std::endl;

  return thickness;
}

std::unique_ptr<ROOT::Math::Interpolator>
RawTreeConverter::LoadProtonAngleCorrections(std::string filename) {
  std::cout << "Loading proton angular corrections file: " <<
            filename << "..." << std::endl;
  std::ifstream ifile(filename);
  std::vector<Double_t> thetas_meas, thetas_eff;

  for (CSVIterator iter(ifile, 1); iter != CSVIterator(); iter++) {
    auto csv_row = *iter;
    assert(csv_row.size() == 6 && "Invalid input file format");

    Double_t theta_meas = std::stod(csv_row[2]) * D2R;
    Double_t theta_eff = std::stod(csv_row[1]) * D2R;

    thetas_meas.push_back(theta_meas);
    thetas_eff.push_back(theta_eff);
  }

  auto interp = new ROOT::Math::Interpolator{
      static_cast<unsigned int>(thetas_meas.size()),
      ROOT::Math::Interpolation::kLINEAR};
  interp->SetData(thetas_meas, thetas_eff);

  using interp_ptr_t = std::unique_ptr<ROOT::Math::Interpolator>;
  return interp_ptr_t(interp);
}

std::unique_ptr<ROOT::Math::Interpolator>
RawTreeConverter::LoadProtonMtrlsEloss(std::string filename) {
  std::cout << "Loading proton angular corrections file: " <<
    filename << "..." << std::endl;
  std::ifstream ifile(filename);
  std::vector<Double_t> thetas_meas, mtrls_e_loss;

  for (CSVIterator iter(ifile, 1); iter != CSVIterator(); iter++) {
    auto csv_row = *iter;
    assert(csv_row.size() == 6 && "Invalid input file format");

    Double_t theta_meas = std::stod(csv_row[1]) * D2R;
    //sum up e-loss in target, chamber & air
    Double_t e_loss = std::stod(csv_row[3])
      + std::stod(csv_row[4]) + std::stod(csv_row[5]);

    thetas_meas.push_back(theta_meas);
    mtrls_e_loss.push_back(e_loss);
  }

  auto interp = new ROOT::Math::Interpolator{
    static_cast<unsigned int>(thetas_meas.size()),
    ROOT::Math::Interpolation::kLINEAR};
  interp->SetData(thetas_meas, mtrls_e_loss);

  using interp_ptr_t = std::unique_ptr<ROOT::Math::Interpolator>;
  return interp_ptr_t(interp);
}


////////////////////////////////////////////////////////////
// Global functions definition
////////////////////////////////////////////////////////////


TTree* MakeScatteringTree(TTree *raw_tree, int run_number) {
  RawTreeConverter raw_converter{raw_tree, run_number};
  return raw_converter.Run();
}

////////////////////////////////////////////////////////////
// main proc
////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "ERROR: This program reqiures run number as first argument!" << std::endl;
    exit(1);
  }

  // everything should be relative to our home directory
  auto ana_root_dir{GetAnaRootDir()};
  std::cout << "Working directory will be set to ana root: " << ana_root_dir \
            << std::endl;
  if (chdir(ana_root_dir.c_str()) != 0) {
    std::cerr << "ERROR: Unable to set program working directory!" << \
      std::endl;
    exit(1);
  }

  int run_no = atoi(argv[1]);
  std::cout << "Will make scattering tree for run. no.: " << run_no << \
    std::endl;
  auto dataset = MakeDataset(std::vector<int> {run_no});
  MakeScatteringTree(dataset.GetTree(), run_no);
  return 0;
}
