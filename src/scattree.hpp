
#pragma once

#include <string>
#include <array>
#include <vector>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
//#include "Math/Polynomial.h"
#include "Math/Interpolator.h"

#include "srccommon.hpp"
#include "naiarraymodel.hpp"
#include "msgoutput.hpp"


namespace s13 {
  namespace ana {

    /*
      Converter for raw ROOT tree, obtains higher-level tree, a "scattering"
      tree. Output tree contains no "raw" data, like signals from individual
      PMTs of a plastic, but averaged and calibrated values for each
      detector. Along with it, this tree contains angular data derived from
      DCs information.
    */
    class RawTreeConverter {

    public:

      RawTreeConverter(TTree *raw_tree, int run_no) :
        m_logger{"RawTreeConverter"},
        mp_raw_tree{raw_tree},
        m_converted{false},
        m_run_no{run_no},
        m_good_proton_events_num{0},
        m_bad_proton_events_num{0},
        m_invalid_nai_data_events_num{0},
        m_invalid_nai_pmt_values_num{0} {
          assert(raw_tree != nullptr);

          Initialize();
        }

      /*
        Main loop
      */
      TTree* Run();

      // constants
      static constexpr Double_t INIT_VAL = -99999.0;
      static constexpr Double_t D2R = 0.0174533;
      static constexpr Double_t DIST_BDCs = 750.;
      static constexpr Double_t DIST_BDC2_TARGET = 1130;
      // thickness of the target [mm]
      static constexpr Double_t TARGET_Z = 2.6;

      // S1 dimensions
      static constexpr Double_t S1DC_WIDTH = 250 * 2;
      static constexpr Double_t S1DC_HEIGHT = 130 * 2;

      // FDC0 dimensions
      static constexpr Double_t FDC0_WIDTH = 80 * 2;
      static constexpr Double_t FDC0_HEIGHT = 80 * 2;

      // BDC dimensions
      static constexpr Double_t BDC_WIDTH = 80 * 2;
      static constexpr Double_t BDC_HEIGHT = 80 * 2;

      static constexpr Double_t DIST_TARGET_FDC0 = 399.5;
      static constexpr Double_t DIST_FDC0_S1DC1 = 272.;
      static constexpr Double_t DIST_FDC0_S1DC2 = 400.;
      static constexpr Double_t DIST_S1MWDC_FDC0 = 336;
      static constexpr Double_t DIST_FDC0_TARGET = 399.5;

      // ESPRIs
      static constexpr int ESPRI_NUM_NAIS = 7;
      static constexpr int NUM_ESPRIS = 2;
      static constexpr int NAIS_NUM = ESPRI_NUM_NAIS * NUM_ESPRIS;
      static constexpr int NAI_PMTS_NUM = NAIS_NUM * 2;
      static constexpr Double_t RDC_CENTER_ANGLE = 62.5;
      static constexpr Double_t DIST_RDC_TARGET = 995.; // 1004.75?
      static constexpr Double_t ESPRI_ACCEPTANCE_WIDTH = 145 * 2;
      static constexpr Double_t ESPRI_RDC_WIDTH = 220 * 2;
      static constexpr Double_t ESPRI_PLASTIC_THICKNESS = 4;

    protected:
      bool IsHe6Run();
      bool IsHe4Run();
      bool IsCarbonRun();
      bool IsEspriLeftEvent();
      bool IsEspriRightEvent();
      bool IsWellDefinedValue(double value);

      /* Various initialization - called during construction of the converter */
      void Initialize();
      void SetupOutputTree();
      void SetupInputTree();
      /* Computes higher level values than in raw tree */
      void ComputeDerivedValues();

      /*
        Checks if NaI PMT value is in range (non-negative) and fixes to zero otherwise.
      */
      double CheckAndFixNaIPmtValue(double nai_pmt_value);


      /* Computes a & b angles for BDCs */
      void ComputeBdcsAngles();
      /* Computes proton recoil angles */
      void ComputeProtonAngles();
      /* Theoretical values useful to build various correlations */
      void ComputeTheoreticalValues();
      /* Computes vertex-Z coordinate using RDCs hit position and a&b angles */
      void ComputeVertexFromRdc(char espri_side, Double_t xpos,
                                Double_t ypos, Double_t aang,
                                Double_t bang);
      /*
        Computes physical energy values from raw detector data
        using calibration cofficients
      */
      void ComputeCalEnergies();

      /*
        Computes effective energy values from calibrated detector data
      */
      void ComputeEffectiveEnergies();

      /*
        Finds energy deposit for an event using calibrated energies
        from NaI array. Method 1 uses max-E as energy deposit.
       */
      void FindEspriCalEnergies1();

      /*
        The same purpose as FindEspriCalEnergies1() but relies on
        extrapolation method to find energy deposit.
       */
      void FindEspriCalEnergies2();

      /*
        Computes NaI ID which corresponds to the highest E deposit.
      */
      int FindMaxENaiId();

      /*
        Computes energy lost by protons in a material by given
        material thickness and its stopping power. E before penetrating material
        will be: E_before = proton_e + ComputeProtonDe(...)
       */
      Double_t ComputeProtonDe(Double_t proton_e,
                               std::unique_ptr<ROOT::Math::Interpolator>& stopping_power,
                               Double_t material_thickness, bool inverse = false);

      /*
        Utility to produce CSV file with proton energies
        including material E-loss & E-deg
      */
      void ComputeProtonEWithMtrlELossesCsv();

      /*
        Computes energy values from kinematics data -
        a simulated values to use as a reference.
      */
      void ComputeSimEnergies();

      /*
        Computes proton e-loss in materials up to RDC front
      */
      Double_t
      ComputeProtonMtrlEloss();

      /*
       * Computes proton angle with correction for magnetic field effect.
       */
      void ComputeProtonEffAngles();

      /*
        Cleans all data vars for trees.
      */
      void ClearTreeVars();

      /*
        Returns He theta, proton theta angles - from reaction
        kinematics (loaded from external CSV file).
      */
      common::ScatteringKinematics LoadKinematicsData(std::string filename);

      /*
        Loads NaI calibration constants to convert ch -> MeV
      */
      std::vector<Double_t> LoadEspriNaiCalConsts(std::string filename, bool load_pedestals = false);

      /*
        Loads dE plastic calibration constants to convert ch -> MeV
      */
      std::vector<Double_t> LoadEspriPdeCalConsts(std::string filename, bool load_pedestals = false);

      /*
        From reaction kinematics file loads energies of recoil proton.

        Params:
          filename: should be a name of a file with kin. data for p+4/6He data
      */
      std::unique_ptr<ROOT::Math::Interpolator> LoadProtonKinE(std::string filename);

      /*
        Loads a stopping power from specified CSV file.
      */
      std::unique_ptr<ROOT::Math::Interpolator> LoadStoppingPower(std::string filename);

      /*
        Computes energy loss of a proton inside ESPRI plastic
        given its stopping power function.
      */
      Double_t ComputeProtonPlaDe(Double_t proton_e,
                                  std::unique_ptr<ROOT::Math::Interpolator>& stp_powers);

      /*
        Computes thickness of the degrader in mm by proton angle.
      */
      Double_t ComputeEdegThickness(Double_t proton_theta);

      /*
       * Loads angular correction to proton data from specified file.
       */
      std::unique_ptr<ROOT::Math::Interpolator> LoadProtonAngleCorrections(std::string filename);

      /*
        Loads proton energy loss data of recoil protons in materials up to
        RDC front as a function of measured proton theta by RDC.
       */
      std::unique_ptr<ROOT::Math::Interpolator>
      LoadProtonMtrlsEloss(std::string filename);

      /*
        Computes a simulated E & dE for a recoil proton to use as a reference
        for dE-E cut.
      */
      std::pair<
        std::shared_ptr<ROOT::Math::Interpolator>,
        std::shared_ptr<ROOT::Math::Interpolator> >
      ComputeRefHeDeE(std::string he_kin_file);

      /// input/output tree vars

      int m_event_no;
      Bool_t m_triggers[8];

      // beamline plastics
      Double_t m_f7_q;
      Double_t m_f7_t;
      Double_t m_sbt1_q;
      Double_t m_sbt1_t;
      Double_t m_sbt2_q;
      Double_t m_sbt2_t;
      Double_t m_sbv_q;
      Double_t m_sbv_t;

      // BDCs
      Double_t m_bdc1_xpos;
      Double_t m_bdc1_ypos;
      Double_t m_bdc2_xpos;
      Double_t m_bdc2_ypos;
      Double_t m_tgt_up_aang;
      Double_t m_tgt_up_bang;
      Double_t m_tgt_up_theta;
      Double_t m_tgt_up_phi;
      Double_t m_tgt_up_xpos;
      Double_t m_tgt_up_ypos;

      // DSDCs
      Double_t m_fdc0_xpos;
      Double_t m_fdc0_ypos;
      Double_t m_fdc0_aang;
      Double_t m_fdc0_bang;
      Double_t m_fdc0_theta;
      Double_t m_fdc0_phi;
      Double_t m_s1dc_xpos;
      Double_t m_s1dc_ypos;
      Double_t m_s1dc_aang;
      Double_t m_s1dc_bang;
      Double_t m_s1dc_theta;
      Double_t m_s1dc_phi;

      Double_t m_tgt_dwn_aang;
      Double_t m_tgt_dwn_bang;
      Double_t m_tgt_dwn_theta;
      Double_t m_tgt_dwn_phi;
      Double_t m_tgt_dwn_xpos;
      Double_t m_tgt_dwn_ypos;

      // ESPRI left
      Double_t m_espri_left_xpos;
      Double_t m_espri_left_ypos;
      Double_t m_espri_left_aang;
      Double_t m_espri_left_bang;
      Double_t m_espri_left_aang_tgt;
      Double_t m_espri_left_bang_tgt;
      Double_t m_espri_left_de;
      Double_t m_espri_left_de_cal;
      Double_t m_espri_left_tof;
      Double_t m_espri_left_nai_e[ESPRI_NUM_NAIS];
      Double_t m_espri_left_nai_e_cal[ESPRI_NUM_NAIS];
      Double_t m_espri_left_e;
      Double_t m_espri_left_e_cal;
      Double_t m_espri_left_e_eff;
      Double_t m_espri_left_vertex_zpos;
      Double_t m_espri_left_proton_theta;
      Double_t m_espri_left_proton_phi;
      Int_t m_espri_left_rdc_nhitwires;
      Int_t m_espri_left_rdc_xnhitwires;
      Int_t m_espri_left_rdc_ynhitwires;

      // ESPRI right
      Double_t m_espri_right_xpos;
      Double_t m_espri_right_ypos;
      Double_t m_espri_right_aang;
      Double_t m_espri_right_bang;
      Double_t m_espri_right_aang_tgt;
      Double_t m_espri_right_bang_tgt;
      Double_t m_espri_right_de;
      Double_t m_espri_right_de_cal;
      Double_t m_espri_right_tof;
      Double_t m_espri_right_nai_e[ESPRI_NUM_NAIS];
      Double_t m_espri_right_nai_e_cal[ESPRI_NUM_NAIS];
      Double_t m_espri_right_e;
      Double_t m_espri_right_e_cal;
      Double_t m_espri_right_e_eff;
      Double_t m_espri_right_vertex_zpos;
      Double_t m_espri_right_proton_theta;
      Double_t m_espri_right_proton_phi;
      Int_t m_espri_right_rdc_nhitwires;
      Int_t m_espri_right_rdc_xnhitwires;
      Int_t m_espri_right_rdc_ynhitwires;

      Double_t m_espri_vertex_zpos;
      Double_t m_proton_theta;
      Double_t m_proton_theta_eff;
      Double_t m_proton_phi;
      Double_t m_proton_e;
      Double_t m_nai_e;
      // effective energy of a proton - e-loss in energy degrader restored
      Double_t m_proton_e_eff;
      Double_t m_proton_e_sim;
      Double_t m_proton_de;
      Double_t m_proton_de_sim;
      Double_t m_proton_de_eff;

      Double_t m_proton_tof;

      // FDC2
      Double_t m_fdc2_xpos;
      Double_t m_fdc2_ypos;
      Double_t m_fdc2_aang;
      Double_t m_fdc2_bang;

      // HODF
      Double_t m_hodf_q[24];
      Double_t m_hodf_t[24];

      // theoretical values
      Double_t m_he_theta_theor;

      /// input tree vars
      Double_t m_raw_espri_nai_qcal[NAI_PMTS_NUM];
      Int_t m_raw_espri_nai_id[NAI_PMTS_NUM];

      // theoretical kinematics
      common::ScatteringKinematics m_he4_kinematics;
      common::ScatteringKinematics m_he6_kinematics;
      // recoil proton E & dE as a function of angle
      std::shared_ptr<ROOT::Math::Interpolator> m_he6_sim_p_e;
      std::shared_ptr<ROOT::Math::Interpolator> m_he6_sim_p_de;
      std::shared_ptr<ROOT::Math::Interpolator> m_he4_sim_p_e;
      std::shared_ptr<ROOT::Math::Interpolator> m_he4_sim_p_de;
      // angular corrections for proton data
      std::shared_ptr<ROOT::Math::Interpolator> m_p_he4_esl_theta_eff_vs_meas;
      std::shared_ptr<ROOT::Math::Interpolator> m_p_he4_esr_theta_eff_vs_meas;
      std::shared_ptr<ROOT::Math::Interpolator> m_p_he6_esl_theta_eff_vs_meas;
      std::shared_ptr<ROOT::Math::Interpolator> m_p_he6_esr_theta_eff_vs_meas;
      // e-loss in materials for recoil protons
      std::shared_ptr<ROOT::Math::Interpolator> m_p_he4_esl_theta_eff_vs_mtrls_e_loss;
      std::shared_ptr<ROOT::Math::Interpolator> m_p_he4_esr_theta_eff_vs_mtrls_e_loss;
      std::shared_ptr<ROOT::Math::Interpolator> m_p_he6_esl_theta_eff_vs_mtrls_e_loss;
      std::shared_ptr<ROOT::Math::Interpolator> m_p_he6_esr_theta_eff_vs_mtrls_e_loss;
      // calibration consts for e & de detectors
      std::vector<Double_t> m_espri_left_nai_cal_consts;
      std::vector<Double_t> m_espri_left_nai_pedestals;
      std::vector<Double_t> m_espri_right_nai_cal_consts;
      std::vector<Double_t> m_espri_right_nai_pedestals;
      Double_t m_espri_left_pde_cal_const;
      Double_t m_espri_left_pde_pedestal;
      Double_t m_espri_right_pde_cal_const;
      Double_t m_espri_right_pde_pedestal;
      // stopping powers
      std::unique_ptr<ROOT::Math::Interpolator> m_pla_stp_power;
      std::unique_ptr<ROOT::Math::Interpolator> m_edeg_stp_power;

      // NaI array model for extrapolation of proton trajectories
      // onto the crystals (for E deposit deduction)
      s13::sim::NaiArrayBox m_nais_array;

    private:
      s13::misc::MessageLogger m_logger;
      TTree *mp_raw_tree;
      bool m_converted;
      int m_run_no;
      TFile *mp_output_file;
      TTree *mp_output_tree;

      // stats of conversion
      int m_good_proton_events_num;
      int m_bad_proton_events_num;
      int m_invalid_nai_data_events_num;
      int m_invalid_nai_pmt_values_num;

    };


    TTree *MakeScatteringTree(TTree *raw_tree);

  }
}
