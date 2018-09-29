/**
 *   \file scattevt.hpp
 *   \brief Definition of ScatteringEvent - used in TTree traversal
 *
 */


#pragma once

#include <assert.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>


namespace s13 {
  namespace ana {

    /*
      Class holds scattering event data plus some convenience methods.
    */
    class ScatteringEvent {

    public:

      ScatteringEvent() {
        Double_t k_inf = -99999;

        event_no = -99999;

        for (int i = 0; i < 8; i++) {
          triggers[i] = false;
        }

        p_phi = k_inf;
        p_theta = k_inf;
        p_theta_eff = k_inf;
        p_aang = k_inf;
        p_bang = k_inf;
        p_e_sim = k_inf;
        p_e_eff = k_inf;
        p_e = k_inf;
        p_de = k_inf;
        p_e_raw = k_inf;
        p_de_raw = k_inf;
        he_phi = k_inf;
        he_theta = k_inf;
        he_theta_theor = k_inf;

        bdc1_xpos = k_inf;
        bdc1_ypos = k_inf;
        bdc2_xpos = k_inf;
        bdc2_ypos = k_inf;

        vert_xpos = k_inf;
        vert_ypos = k_inf;
        tgt_up_aang = k_inf;
        tgt_up_bang = k_inf;

        esl_vert_zpos = k_inf;
        esr_vert_zpos = k_inf;
        es_vert_zpos = k_inf;

        esl_aang = k_inf;
        esl_bang = k_inf;
        esr_aang = k_inf;
        esr_bang = k_inf;

        esl_xpos = k_inf;
        esl_ypos = k_inf;
        esr_xpos = k_inf;
        esr_ypos = k_inf;

        esl_e_cal = k_inf;
        esr_e_cal = k_inf;
        esl_e_eff = k_inf;
        esr_e_eff = k_inf;
        esl_e_raw = k_inf;
        esr_e_raw = k_inf;

        esl_de_raw = k_inf;
        esr_de_raw = k_inf;

        esl_nhitwires = k_inf;
        esl_xnhitwires = k_inf;
        esl_ynhitwires = k_inf;
        esr_nhitwires = k_inf;
        esr_xnhitwires = k_inf;
        esr_ynhitwires = k_inf;

        s1dc_phi = k_inf;
        s1dc_theta = k_inf;
        s1dc_xpos = k_inf;
        s1dc_ypos = k_inf;
        s1dc_aang = k_inf;
        s1dc_bang = k_inf;

        fdc0_phi = k_inf;
        fdc0_theta = k_inf;
        fdc0_xpos = k_inf;
        fdc0_ypos = k_inf;
        fdc0_aang = k_inf;
        fdc0_bang = k_inf;

        tgt_dwn_vert_xpos = k_inf;
        tgt_dwn_vert_ypos = k_inf;
        tgt_dwn_aang = k_inf;
        tgt_dwn_bang = k_inf;

        fdc2_xpos = k_inf;
        fdc2_ypos = k_inf;
        fdc2_aang = k_inf;
        fdc2_bang = k_inf;

        for (int i = 0; i < 7; i++) {
          esl_naie_cal[i] = k_inf;
          esr_naie_cal[i] = k_inf;
          esl_naie_raw[i] = k_inf;
          esr_naie_raw[i] = k_inf;
        }

        for (int i = 0; i < 24; i++) {
          hodf_q[i] = k_inf;
          hodf_t[i] = k_inf;
        }
      }

      bool IsLeftEspriEvent() const {
        // return evt.p_phi < 0;
        return esl_de_raw > 240;
      }

      bool IsRightEspriEvent() const {
        return esr_de_raw > 240;
      }

      bool IsEspriEvent() const {
        return IsLeftEspriEvent() || IsRightEspriEvent();
      }

      long event_no;

      Bool_t triggers[8];

      Double_t p_phi;
      Double_t p_theta;
      Double_t p_theta_eff;
      Double_t p_aang;
      Double_t p_bang;
      Double_t p_e_sim;
      Double_t p_e_eff;
      Double_t p_e;
      Double_t p_de;
      Double_t p_e_raw;
      Double_t p_de_raw;
      Double_t he_phi;
      Double_t he_theta;
      Double_t he_theta_theor;

      Double_t vert_xpos;
      Double_t vert_ypos;
      Double_t bdc1_xpos;
      Double_t bdc1_ypos;
      Double_t bdc2_xpos;
      Double_t bdc2_ypos;
      Double_t tgt_up_aang;
      Double_t tgt_up_bang;

      Double_t esl_vert_zpos;
      Double_t esr_vert_zpos;
      Double_t es_vert_zpos;

      Double_t esl_aang;
      Double_t esl_bang;
      Double_t esr_aang;
      Double_t esr_bang;

      Double_t esl_xpos;
      Double_t esl_ypos;
      Double_t esr_xpos;
      Double_t esr_ypos;

      Double_t esl_e_cal;
      Double_t esr_e_cal;
      Double_t esl_e_eff;
      Double_t esr_e_eff;
      Double_t esl_e_raw;
      Double_t esr_e_raw;

      Double_t esl_de_raw;
      Double_t esr_de_raw;

      Double_t esl_naie_cal[7];
      Double_t esr_naie_cal[7];
      Double_t esl_naie_raw[7];
      Double_t esr_naie_raw[7];

      Int_t esr_nhitwires;
      Int_t esr_xnhitwires;
      Int_t esr_ynhitwires;
      Int_t esl_nhitwires;
      Int_t esl_xnhitwires;
      Int_t esl_ynhitwires;

      Double_t hodf_q[24];
      Double_t hodf_t[24];

      Double_t fdc0_phi;
      Double_t fdc0_theta;
      Double_t fdc0_xpos;
      Double_t fdc0_ypos;
      Double_t fdc0_aang;
      Double_t fdc0_bang;

      Double_t s1dc_phi;
      Double_t s1dc_theta;
      Double_t s1dc_xpos;
      Double_t s1dc_ypos;
      Double_t s1dc_aang;
      Double_t s1dc_bang;

      Double_t tgt_dwn_vert_xpos;
      Double_t tgt_dwn_vert_ypos;
      Double_t tgt_dwn_aang;
      Double_t tgt_dwn_bang;

      Double_t fdc2_xpos;
      Double_t fdc2_ypos;
      Double_t fdc2_aang;
      Double_t fdc2_bang;

    };

  }
}
