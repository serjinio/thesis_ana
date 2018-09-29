/**
 *   \file elacuts.hpp
 *   \brief Cuts for elastic scattering event selection.
 *
 */

#pragma once

#include "TChain.h"
#include "TGraph.h"
#include "TCut.h"

#include "common.hpp"
#include "consts.hpp"
#include "scattevt.hpp"
#include "msgoutput.hpp"
#include "varthetacut.hpp"
#include "vertzalg.hpp"
#include "ecututils.hpp"


namespace s13 {
  namespace ana {

    class ElasticScatteringCuts {

    public:
      ElasticScatteringCuts() :
        apply_reaction_trigger_cut{true},
        apply_beam_trigger_cut{false},
        apply_beam_abang_cut{false},
        apply_vertexXY_cut{true},
        apply_vertexXY_elliptical_cut{false},
        apply_phi_cut{false},
        use_var_phi_cut{true},
        apply_theta_cut{false},
        use_var_theta_cut{false},
        can_use_var_theta_cut{false},
        apply_hodf_he6_cut{false},
        apply_hodf_he4_cut{false},
        subtract_he6_carbon_bg{true},
        apply_espri_aang_cut{false},
        use_espri_aang_corr{true},
        apply_espri_vertz_cut{false},
        apply_espri_min_e_de_cut{false},
        apply_espri_e_cut{false},
        espri_e_cut_use_he6_rel_e_corr{false},
        espri_e_cut_use_he4_rel_e_corr{false},
        apply_espri_y_cut{false},
        apply_rdc_xnhitplanes_cut{false},
        apply_rdc_ynhitplanes_cut{false},
        apply_user_cut{false},
        vertexXY_radius{10},
        vertexXY_ellipse_xR{9},
        vertexXY_ellipse_yR{3},
        beam_aang_width{2.0},
        beam_bang_width{2.0},
        vertexXY_Xdisp{-2},
        vertexXY_Ydisp{-2},
        phi_width{3},
        theta_width{4},
        espri_e_width_stdev{0.9},
        espri_aang_width{92},
        espri_vertz_width{92},
        espri_y_width{132},
        espri_selector{Espri::both},
        dsdc_selector{Dsdcs::s1dc},
        rdc_xnhitplanes{3},
        rdc_ynhitplanes{3},
        user_cut{""},
        lin_phi_cut_fdc0{gk_lin_phi_cut_fdc0_a, gk_lin_phi_cut_fdc0_b},
        lin_phi_cut_s1dc{gk_lin_phi_cut_s1dc_a, gk_lin_phi_cut_s1dc_b},
        esl_aang_mean{gk_espri_left_aang_a, gk_espri_left_aang_b},
        esr_aang_mean{gk_espri_right_aang_a, gk_espri_right_aang_b},
        espri_computed_aang_mean{gk_espri_computed_aang_a, gk_espri_computed_aang_b},
        he6_e_dist_sigma_interp{make_he6_e_dist_sigma_interp()},
        he4_e_dist_sigma_interp{make_he4_e_dist_sigma_interp()},
        he6_esl_e_dist_mean_interp{make_he6_e_dist_mean_interp(Espri::left)},
        he6_esr_e_dist_mean_interp{make_he6_e_dist_mean_interp(Espri::right)},
        he4_esl_e_dist_mean_interp{make_he4_e_dist_mean_interp(Espri::left)},
        he4_esr_e_dist_mean_interp{make_he4_e_dist_mean_interp(Espri::right)},
        logger_{"ElasticScatteringCuts"}
      {}


      // constants
      static constexpr Double_t k_phi_dist_center_deg = 177.5;
      static constexpr Double_t k_phi_bg_min_angle = 120;
      static constexpr Double_t k_phi_bg_max_angle = 240;
      static constexpr Double_t k_esl_zpos_offset = -5.6;
      static constexpr Double_t k_esr_zpos_offset = -13.65;
      static constexpr Double_t k_esl_aang_center_offset = -4.73;
      static constexpr Double_t k_esr_aang_center_offset = -6.22;
      static constexpr Double_t k_dist_to_rdc = 1004.75;
      static constexpr Double_t k_esr_e_dist_center_offset = -0.05;
      static constexpr Double_t k_esl_e_dist_center_offset = 0.1;

      // flags to apply cuts or not
      bool apply_reaction_trigger_cut;
      bool apply_beam_trigger_cut;
      bool apply_beam_abang_cut;
      bool apply_vertexXY_cut;
      bool apply_vertexXY_elliptical_cut;
      bool apply_phi_cut;
      bool use_var_phi_cut;
      bool apply_theta_cut;
      bool use_var_theta_cut;
      bool can_use_var_theta_cut;
      bool apply_hodf_he6_cut;
      bool apply_hodf_he4_cut;
      bool subtract_he6_carbon_bg;
      bool apply_espri_aang_cut;
      bool use_espri_aang_corr;
      bool apply_espri_vertz_cut;
      bool apply_espri_min_e_de_cut;
      bool apply_espri_e_cut;
      bool espri_e_cut_use_he6_rel_e_corr;
      bool espri_e_cut_use_he4_rel_e_corr;
      bool apply_espri_y_cut;
      bool apply_rdc_xnhitplanes_cut;
      bool apply_rdc_ynhitplanes_cut;

      // user-defined cut
      bool apply_user_cut;

      // cut parameters
      Double_t vertexXY_radius;
      Double_t vertexXY_ellipse_xR;
      Double_t vertexXY_ellipse_yR;
      Double_t beam_aang_width;
      Double_t beam_bang_width;
      Double_t vertexXY_Xdisp;
      Double_t vertexXY_Ydisp;
      Double_t phi_width;
      // width of proton theta angle spread from kinematics
      Double_t theta_width;
      // width of E spread from reference kinematics
      Double_t espri_e_width_stdev;
      // width of angular window in mrad
      Double_t espri_aang_width;
      // width of vertZ distribution
      Double_t espri_vertz_width;
      // width of Y cut (Y-acceptance) - used for CS computation
      Double_t espri_y_width;
      // which ESPRI to use
      Espri espri_selector;
      // which of DSDCs to use
      Dsdcs dsdc_selector;
      // how many RDC planes should be hit
      Int_t rdc_xnhitplanes;
      Int_t rdc_ynhitplanes;

      //user cut
      TCut user_cut;

      /*
        Returns target cut by radius.
        NOTE: center of the beam for 4He is (-2,-2)
      */
      TCut GetVertexXYCut() const {
        return GetCut("pow(tgt_up_xpos + 2, 2) + pow(tgt_up_ypos + 2, 2) < %.2f",
                      std::pow(vertexXY_radius, 2));
      }

      TCut GetThetaCut() const {
        return GetCut("abs(s1dc_theta*57.3 - he_theta_theor*57.3) < %.2f",
                      theta_width / 2.);
      }

      TCut GetPhiCut() const {
        return GetCut("abs(abs(p_phi*57.3-s1dc_phi*57.3) - %.2f) < %.2f",
                      k_phi_dist_center_deg, phi_width / 2.);
      }

      TCut GetPhiBgCut() const {
        return !GetPhiCut() &&
          GetCut("abs(p_phi*57.3-s1dc_phi*57.3) > %.2f && "
                 "abs(p_phi*57.3-s1dc_phi*57.3) < %.2f",
                 k_phi_bg_min_angle, k_phi_bg_max_angle);
      }

      TCut GetHe6HodfCut() const {
        std::string hodf_cut_str = "";
        std::string hodf_cut_left, hodf_cut_right;

        for (int i = 0; i <= 5; i++) {
          //hodf_cut_left += TString::Format("(hodf_q[%i]>14 && hodf_q[%i]<17) ", i, i);
          hodf_cut_left += TString::Format("hodf_q[%i]>7", i);
          if (i != 5) hodf_cut_left += " || ";
        }
        hodf_cut_left = TString::Format("( (%s) && s1dc_xpos > 0 )", hodf_cut_left.data());

        for (int i = 11; i <= 15; i++) {
          //hodf_cut_str += TString::Format("(hodf_q[%i]>14 && hodf_q[%i]<17) ", i, i);
          hodf_cut_right += TString::Format("hodf_q[%i]>7", i);
          if (i != 15) hodf_cut_right += " || ";
        }
        hodf_cut_right = TString::Format("( (%s) && s1dc_xpos < 0 )", hodf_cut_right.data());

        hodf_cut_str = TString::Format("%s || %s", hodf_cut_left.data(), hodf_cut_right.data());
        //std::cout << "HODF cut: " << hodf_cut_str << std::endl;
        return TCut{hodf_cut_str.data()};
      }

      /*
       * Cut without adjustable parameters - cuts off garbage events in terms
       * of measured E & dE.
       * BUG: NaI E data is not well-defined for about 50% of the events measured.
       */
      TCut GetEspriMinEdECut() const {
        //return "p_de > 1 && p_e_eff > 5";
        return "p_e_eff > 0.43";
      }

      TCut GetEslECut() const {
        return GetCut("abs((esl_e_eff - p_e_sim) / p_e_sim - %.2f) < %.2f",
                      k_esl_e_dist_center_offset, espri_e_width_stdev / 2);
      }

      TCut GetEsrECut() const {
        return GetCut("abs((esr_e_eff - p_e_sim) / p_e_sim - %.2f) < %.2f",
                      k_esr_e_dist_center_offset, espri_e_width_stdev / 2);
      }

      TCut GetEspriECut() const {
        switch (espri_selector) {
        case Espri::left:
          return GetEslECut();
        case Espri::right:
          return GetEsrECut();
        case Espri::both:
          return (GetEslECut() && "p_phi < 0") || (GetEsrECut() && "p_phi > 0");
        default:
          throw std::logic_error("Invalid ESPRI selector value!");
        }
      }

      TCut GetEslAangCut() const {
        return GetCut("abs(esl_aang - (esl_xpos/%.2f * 1e3) - %.2f) < %.2f",
                      k_dist_to_rdc, k_esl_aang_center_offset,
                      espri_aang_width / 2);
      }

      TCut GetEsrAangCut() const {
        return GetCut("abs(esr_aang - (esr_xpos/%.2f * 1e3) - %.2f) < %.2f",
                      k_dist_to_rdc, k_esr_aang_center_offset,
                      espri_aang_width / 2);
      }

      TCut GetEspriAangCut() const {
        switch (espri_selector) {
        case Espri::left:
          return GetEslAangCut();
        case Espri::right:
          return GetEsrAangCut();
        case Espri::both:
          return (GetEslAangCut() && "p_phi < 0") ||
            (GetEsrAangCut() && "p_phi > 0");
        default:
          throw std::logic_error("Invalid ESPRI selector value!");
        }
      }

      TCut GetEspriSelectionCut() const {
        if (espri_selector == Espri::both) return TCut{};
        else if (espri_selector == Espri::left) return "p_phi < 0";
        else return "p_phi > 0";
      }

      TCut GetEspriYCut() const {
        TCut left = GetCut("p_phi < 0 && abs(esl_ypos) < %.2f", espri_y_width);
        TCut right = GetCut("p_phi > 0 && abs(esr_ypos) < %.2f", espri_y_width);
        if (espri_selector == Espri::both) {
          return left || right;
        } else if (espri_selector == Espri::left) {
          return left;
        } else {
          return right;
        }
      }

      // DEPRECATED
      TCut GetTotalCut() const {
        throw std::invalid_argument("not implemented!");

        // beam && pdE trigger
        TCut reaction_trigger_cut = "triggers[5]==1";
        TCut beam_trigger_cut = "triggers[1]==1";
        // beam only / 2000
        //TCut reaction_trigger_cut = "triggers[1]==1";

        TCut total_cut = GetEspriSelectionCut();
        if (apply_reaction_trigger_cut) total_cut = total_cut && reaction_trigger_cut;
        if (apply_beam_trigger_cut) total_cut = total_cut && beam_trigger_cut;
        if (apply_vertexXY_cut) total_cut = total_cut && GetVertexXYCut();
        if (apply_phi_cut) total_cut = total_cut && GetPhiCut();
        if (apply_theta_cut) total_cut = total_cut && GetThetaCut();
        if (apply_espri_min_e_de_cut) total_cut = total_cut && GetEspriMinEdECut();
        if (apply_espri_e_cut) total_cut = total_cut && GetEspriECut();
        if (apply_espri_aang_cut) total_cut = total_cut && GetEspriAangCut();
        if (apply_hodf_he6_cut) total_cut = total_cut && GetHe6HodfCut();
        //if (apply_hodf_he4_cut) total_cut = total_cut && GetHe4HodfCut();
        if (apply_espri_y_cut) total_cut = total_cut && GetEspriYCut();
        if (apply_user_cut) total_cut = total_cut && user_cut;

        std::cout << "[DEBUG:elacuts] Total cut: " << total_cut << std::endl;
        return total_cut;
      }

      // DEPRECATED
      TCut GetTotalPhiBgCut() const {
        throw std::invalid_argument("not implemented!");

        TCut reaction_trigger_cut = "triggers[5]==1";
        TCut beam_trigger_cut = "triggers[1]==1";
        TCut total_cut = GetEspriSelectionCut();
        if (apply_reaction_trigger_cut) total_cut = total_cut && reaction_trigger_cut;
        if (apply_beam_trigger_cut) total_cut = total_cut && beam_trigger_cut;
        if (apply_vertexXY_cut) total_cut = total_cut && GetVertexXYCut();
        // unconditionally apply phi bg cut. if using BG extraction, then total cut
        // must include phi cut (apply_phi_cut = true)
        total_cut = total_cut && GetPhiBgCut();
        if (apply_theta_cut) total_cut = total_cut && GetThetaCut();
        if (apply_espri_min_e_de_cut) total_cut = total_cut && GetEspriMinEdECut();
        if (apply_espri_e_cut) total_cut = total_cut && GetEspriECut();
        if (apply_espri_aang_cut) total_cut = total_cut && GetEspriAangCut();
        if (apply_hodf_he6_cut) total_cut = total_cut && GetHe6HodfCut();
        if (apply_espri_y_cut) total_cut = total_cut && GetEspriYCut();
        if (apply_user_cut) total_cut = total_cut && user_cut;

        std::cout << "[DEBUG] total phi BG cut: " << total_cut << std::endl;
        return total_cut;
      }

      TString GetVertexXYCutTitle() const {
        return TString::Format("v-XY(R<%.1f)", vertexXY_radius);
      }

      TString GetVertexXYEllipticalCutTitle() const {
        return TString::Format("v-XY(xR<%.1f;yR<%.1f)",
                               vertexXY_ellipse_xR, vertexXY_ellipse_yR);
      }

      TString GetThetaCutTitle() const {
        return TString::Format("#theta(+/-%.1f)", theta_width/2.);
      }

      TString GetEspriMinEdECutTitle() const {
        return "minEdE";
      }

      TString GetEspriECutTitle() const {
        return TString::Format("p_E(+/-%.1f)", espri_e_width_stdev/2.);
      }

      TString GetEspriAangCutTitle() const {
        return TString::Format("rdc_a(+/-%.1f)", espri_aang_width/2);
      }

      TString GetEspriVertzCutTitle() const {
        return TString::Format("rdc_vertz(+/-%.1f)", espri_vertz_width/2);
      }

      TString GetPhiCutTitle() const {
        return TString::Format("#phi(+/-%.1f)", phi_width/2.);
      }

      TString GetHe6HodfCutTitle() const {
        return TString{"he6@hodf"};
      }

      TString GetHe4HodfCutTitle() const {
        return TString{"he4@hodf"};
      }

      TString GetEspriYCutTitle() const {
        return TString::Format("Y(+/-%.0f)", espri_y_width);
      }

      TString GetUserCutTitle() const {
        return TString(user_cut);
      }

      TString GetTotalCutTitle() const {
        TString title = "";
        if (espri_selector == Espri::both) title += "";
        else if (espri_selector == Espri::left) title += "ESL";
        else title += "ESR";
        if (apply_vertexXY_cut) title += " " + GetVertexXYCutTitle();
        if (apply_vertexXY_elliptical_cut) title += " " + GetVertexXYEllipticalCutTitle();
        if (apply_phi_cut) title += " " + GetPhiCutTitle();
        if (apply_theta_cut) title += " " + GetThetaCutTitle();
        if (apply_espri_min_e_de_cut) title += " " + GetEspriMinEdECutTitle();
        if (apply_espri_e_cut) title += " " + GetEspriECutTitle();
        if (apply_espri_aang_cut) title += " " + GetEspriAangCutTitle();
        if (apply_espri_vertz_cut) title += " " + GetEspriVertzCutTitle();
        if (apply_hodf_he6_cut) title += " " + GetHe6HodfCutTitle();
        if (apply_hodf_he4_cut) title += " " + GetHe4HodfCutTitle();
        if (apply_espri_y_cut) title += " " + GetEspriYCutTitle();
        //if (apply_user_cut) title += " " + GetUserCutTitle();

        return title;
      }

      TString GetTotalCutName() const {
        TString name = "";

        if (espri_selector == Espri::both) name += "eslr_";
        else if (espri_selector == Espri::left) name += "esl_";
        else name += "esr_";

        if (apply_vertexXY_cut)
          name += TString::Format("vR%i", (int)vertexXY_radius);
        if (apply_vertexXY_elliptical_cut)
          name += TString::Format("xR%iyR%i",
                                  (int)vertexXY_ellipse_xR,
                                  (int)vertexXY_ellipse_yR);
        if (apply_beam_abang_cut)
          name += TString::Format("Bab%.1f;%.1f", beam_aang_width, beam_bang_width);
        if (apply_espri_vertz_cut)
          name += TString::Format("vZ%i", (int)espri_vertz_width);
        if (apply_espri_aang_cut)
          name += TString::Format("aang%i", (int)espri_aang_width);
        if (apply_phi_cut) {
          name += TString::Format("phi%.2f", phi_width);
        }
        if (apply_theta_cut)
          name += TString::Format("th%.2f", theta_width);
        if (apply_espri_e_cut)
          name += TString::Format("E%.2f", espri_e_width_stdev);
        if (apply_espri_y_cut)
          name += TString::Format("esy%i", (int)espri_y_width);
        if (apply_hodf_he6_cut)
          name += TString::Format("he6hodf");
        if (apply_hodf_he4_cut)
          name += TString::Format("he4hodf");

        return name;
      }

      bool IsVertexXYCut(const ScatteringEvent& evt) const {
        if (std::pow(evt.vert_xpos - vertexXY_Xdisp, 2) +
            std::pow(evt.vert_ypos - vertexXY_Ydisp, 2) > std::pow(vertexXY_radius, 2)) {
          return false;
        }

        // cuts for high pol. region selection
        // if (evt.vert_xpos > -3 && evt.vert_ypos > 3) {
        //   return false;
        // }
        // if (evt.vert_xpos < -3 && std::abs(evt.vert_ypos) < 3) {
        //   return false;
        // }
        // if (std::abs(evt.vert_xpos) < 3 && evt.vert_ypos < -3) {
        //   return false;
        // }

        return true;
      }

      bool IsVertexXYEllipticalCut(const ScatteringEvent& evt) const {
        double vX = evt.vert_xpos;
        double vY = evt.vert_ypos;
        double xdisp = vertexXY_Xdisp;
        double ydisp = vertexXY_Ydisp;
        double xR = vertexXY_ellipse_xR;
        double yR = vertexXY_ellipse_yR;

        double lhs = std::pow((vX - xdisp) / xR, 2) + std::pow((vY - ydisp) / yR, 2);
        // if (lhs <= 1) {
        //   std::cout << "Evaluating event with (X,Y): " << vX << "; " << vY
        //             << ". lhs: " << lhs << std::endl;
        // }
        return lhs <= 1;
      }

      bool IsBeamAbangCut(const ScatteringEvent& evt) const {
        double beam_aang = evt.tgt_up_aang * 1/s13::gk_d2r;
        double beam_bang = evt.tgt_up_bang * 1/s13::gk_d2r;
        if ((std::abs(beam_aang) > beam_aang_width / 2) ||
            (std::abs(beam_bang) > beam_bang_width / 2)) {
          return false;
        }

        return true;
      }

      // DEPRECATED: Only variable cut is used
      bool IsConstPhiCut(const ScatteringEvent& evt) const {
        Double_t cut_halfwidth = phi_width / 2.;
        if (std::abs(std::abs(evt.p_phi - GetHePhi(evt)) - k_phi_dist_center_deg)
            > cut_halfwidth) {
          return false;
        }

        return true;
      }

      Double_t GetVarPhiCutSigma(Double_t theta) const {
        Double_t cut_sigma = dsdc_selector == Dsdcs::fdc0
          ? lin_phi_cut_fdc0.eval(theta) : lin_phi_cut_s1dc.eval(theta);
        return cut_sigma;
      }

      Double_t GetPhiCutWidth(Double_t theta) const {
        if (use_var_phi_cut) {
          return 2 * GetVarPhiCutSigma(theta) * phi_width;
        } else {
          throw std::logic_error("Use of const phi cut is deprecated!");
          return phi_width;
        }
      }

      bool IsVarPhiCut(const ScatteringEvent& evt) const {
        Double_t cut_sigma = GetVarPhiCutSigma(evt.p_theta_eff);
        double min_phi_angle = k_phi_dist_center_deg - phi_width * cut_sigma;
        double max_phi_angle = k_phi_dist_center_deg + phi_width * cut_sigma;
        double he_phi_diff = std::abs(evt.p_phi - GetHePhi(evt));

        if (he_phi_diff > min_phi_angle && he_phi_diff < max_phi_angle) {
          return true;
        } else {
          return false;
        }
      }

      bool IsPhiCut(const ScatteringEvent& evt) const {
        if (use_var_phi_cut) {
          return IsVarPhiCut(evt);
        } else {
          throw std::logic_error("Const phi cut is deprecated and should not be used!");
          return IsConstPhiCut(evt);
        }
      }

      /*
        Selects BG events which fall out of phi cut.
        Also limits region of phi angle 120 - 240 deg. It is done to be intact
        with definition of R_norm - on the same range.
      */
      bool IsPhiBgCut(const ScatteringEvent& evt) const {
        Double_t phi_abs_diff = std::abs(evt.p_phi - GetHePhi(evt));
        // event should NOT fall into phi gate
        if (IsPhiCut(evt)) {
          return false;
        }

        // event should stay in same region where we defined R_norm
        if (phi_abs_diff < k_phi_bg_min_angle || phi_abs_diff > k_phi_bg_max_angle) {
          return false;
        }

        return true;
      }

      bool IsConstThetaCut(const ScatteringEvent& evt) const {
        Double_t theta_cut_halfwidth = theta_width / 2.;
        if (std::abs(GetHeTheta(evt) - evt.he_theta_theor)
            > theta_cut_halfwidth) {
          return false;
        }

        return true;
      }

      bool IsVarThetaCut(const ScatteringEvent& evt) const {
        if (!can_use_var_theta_cut) {
          throw std::logic_error("Variable theta cut resources were not prepared!");
        }
        // We cannot compute theta cut parameters if proton angle
        // is out of RDC acceptance.
        if (evt.p_theta_eff < 51 || evt.p_theta_eff > 75) {
          return false;
        }
        // Also we should be able to determine to which
        // ESPRI given event belongs
        if (!IsLeftEspriEvent(evt) && !IsRightEspriEvent(evt)) {
          return false;
        }
        double mean = GetThetaMeanInterp(evt)->Eval(evt.p_theta_eff);
        double sigma = GetThetaSigmaInterp(evt)->Eval(evt.p_theta_eff);
        double min_theta_angle = mean - 3 * sigma;
        double max_theta_angle = mean + 3 * sigma;

        // logger_.debug("Evaluating event at p_theta_eff: %.2f, "
        //               "theta cut (min, max): %.2f, %.2f.",
        //               evt.p_theta_eff, min_theta_angle, max_theta_angle);
        double he_theta_diff = GetHeTheta(evt) - evt.he_theta_theor;
        if (he_theta_diff > min_theta_angle && he_theta_diff < max_theta_angle) {
          // logger_.debug("good event!");
          return true;
        } else {
          return false;
        }
      }

      bool IsThetaCut(const ScatteringEvent& evt) const {
        if (use_var_theta_cut) {
          return IsVarThetaCut(evt);
        } else {
          return IsConstThetaCut(evt);
        }
      }

      bool IsEspriMinEdECut(const ScatteringEvent& evt) const {
        return evt.p_de > 1 && evt.p_e - evt.p_de > 0.25;
      }

      bool IsEspriECut(const ScatteringEvent& evt) const {
        if (!IsWellDefined(evt.p_theta_eff) || !evt.IsEspriEvent()) {
          return false;
        }
        Double_t cut_hw = (*he6_e_dist_sigma_interp)(evt.p_theta_eff) *
          espri_e_width_stdev;
        Double_t cut_dist_mean = 0.0;
        if (espri_e_cut_use_he6_rel_e_corr) {
          auto& interp = IsLeftEspriEvent(evt) ?
            *he6_esl_e_dist_mean_interp : *he6_esr_e_dist_mean_interp;
          // NaI mean correction
          cut_dist_mean += he6_nai_e_rel_interp->eval(evt);
          // Overall E-dist mean correction
          cut_dist_mean += interp(evt.p_theta_eff);
        }
        if (espri_e_cut_use_he4_rel_e_corr) {
          auto& interp = IsLeftEspriEvent(evt) ?
            *he4_esl_e_dist_mean_interp : *he4_esr_e_dist_mean_interp;
          // NaI mean correction
          cut_dist_mean += he4_nai_e_rel_interp->eval(evt);
          // Overall E-dist mean correction
          cut_dist_mean += interp(evt.p_theta_eff);
        }

        if (evt.p_theta_eff < 49.7 || evt.p_theta_eff > 75) {
          logger_.debug("WARN: out of range proton angle! "
                        "Evaluated E cut HW: %.2f and mean: %.2f "
                        "at an angle %.2f deg.",
                        cut_hw, cut_dist_mean, evt.p_theta_eff);
        }
        Double_t esl_e_diff =
          ((evt.p_e - evt.p_e_sim) / evt.p_e_sim) - cut_dist_mean;
        Double_t esr_e_diff =
          ((evt.p_e - evt.p_e_sim) / evt.p_e_sim) - cut_dist_mean;
        bool is_esl = std::abs(esl_e_diff) < cut_hw;
        bool is_esr = std::abs(esr_e_diff) < cut_hw;

        switch (espri_selector) {
        case Espri::left:
          return IsLeftEspriEvent(evt) && is_esl;
        case Espri::right:
          return IsRightEspriEvent(evt) && is_esr;
        case Espri::both:
          return (IsLeftEspriEvent(evt) && is_esl) ||
            (IsRightEspriEvent(evt) && is_esr);
        default:
          throw std::logic_error("Invalid ESPRI selector value!");
        }
      }

      bool IsEspriAangCut(const ScatteringEvent &evt) const {
        Double_t aang_halfwidth = espri_aang_width / 2;
        Double_t esl_aang = GetCorrectedEspriAang(evt, Espri::left);
        Double_t esr_aang = GetCorrectedEspriAang(evt, Espri::right);

        Double_t esl_aang_diff = esl_aang
          - (evt.esl_xpos/s13::gk_dist_rdc_target * 1e3);
        Double_t esr_aang_diff = esr_aang
          - (evt.esr_xpos/s13::gk_dist_rdc_target * 1e3);

        bool is_esl = std::abs(esl_aang_diff) < aang_halfwidth;
        bool is_esr = std::abs(esr_aang_diff) < aang_halfwidth;

        switch (espri_selector) {
        case Espri::left:
          return IsLeftEspriEvent(evt) && is_esl;
        case Espri::right:
          return IsRightEspriEvent(evt) && is_esr;
        case Espri::both:
          return (IsLeftEspriEvent(evt) && is_esl) ||
            (IsRightEspriEvent(evt) && is_esr);
        default:
          throw std::logic_error("Invalid ESPRI selector value!");
        }
      }

      bool IsEspriVertzCut(const ScatteringEvent& evt) const {
        Double_t vertz_halfwidth = espri_vertz_width / 2;
        Double_t esl_aang = GetCorrectedEspriAang(evt, Espri::left);
        Double_t esr_aang = GetCorrectedEspriAang(evt, Espri::right);
        Double_t esl_vertz = compute_vertz(evt.esl_xpos, esl_aang, Espri::left);
        Double_t esr_vertz = compute_vertz(evt.esr_xpos, esr_aang, Espri::right);

        bool is_esl = std::abs(esl_vertz) < vertz_halfwidth;
        bool is_esr = std::abs(esr_vertz) < vertz_halfwidth;
        switch (espri_selector) {
        case Espri::left:
          return IsLeftEspriEvent(evt) && is_esl;
        case Espri::right:
          return IsRightEspriEvent(evt) && is_esr;
        case Espri::both:
          return (IsLeftEspriEvent(evt) && is_esl) ||
            (IsRightEspriEvent(evt) && is_esr);
        default:
          throw std::logic_error("Invalid ESPRI selector value!");
        }
      }

      bool IsHodfHe6Cut(const ScatteringEvent& evt) const {
        bool is_left_scattered_he6 = false;
        bool is_right_scattered_he6 = false;
        double de = get_hodf_q(evt);
        double id = get_hodf_hit_id(evt);

        if (id >= 0 && id <= 4) {
          if (de > 7) {
            is_left_scattered_he6 = true;
          }
        }

        if (id >= 7 && id <= 15) {
          if (de > 7) {
            is_right_scattered_he6 = true;
          }
        }

        if (espri_selector == Espri::both) {
          return is_right_scattered_he6 || is_left_scattered_he6;
        } else if (espri_selector == Espri::left) {
          return is_right_scattered_he6;
        } else if (espri_selector == Espri::right) {
          return is_left_scattered_he6;
        } else {
          throw std::logic_error("Invalid ESPRI selector.");
        }
      }

      bool IsHodfHe4Cut(const ScatteringEvent& evt) const {
        bool is_he4 = false;
        double de = get_hodf_q(evt);
        double id = get_hodf_hit_id(evt);

        if (id >= 17) {
          if (de > 7) {
            is_he4 = true;
          }
        }

        return is_he4;
      }

      bool IsEspriYCut(const ScatteringEvent& evt) const {
        bool is_left = std::abs(evt.esl_ypos) < espri_y_width;
        bool is_right = std::abs(evt.esr_ypos) < espri_y_width;
        if (espri_selector == Espri::both) {
          return is_left || is_right;
        } else if (espri_selector == Espri::left) {
          return is_left;
        } else {
          return is_right;
        }
      }

      bool IsRdcXNHitPlanesCut(const ScatteringEvent& evt) const {
        if (IsLeftEspriEvent(evt) && evt.esl_xnhitwires >= rdc_xnhitplanes) {
          return true;
        } else if (IsRightEspriEvent(evt) && evt.esr_xnhitwires >= rdc_xnhitplanes) {
          return true;
        } else {
          return false;
        }
      }

      bool IsRdcYNHitPlanesCut(const ScatteringEvent& evt) const {
        if (IsLeftEspriEvent(evt) && evt.esl_ynhitwires >= rdc_ynhitplanes) {
          return true;
        } else if (IsRightEspriEvent(evt) && evt.esr_ynhitwires >= rdc_ynhitplanes) {
          return true;
        } else {
          return false;
        }
      }

      bool IsBeamTriggerCut(const ScatteringEvent& evt) const {
        if (evt.triggers[1] != true) {
          return false;
        }

        return true;
      }

      bool IsReactionTriggerCut(const ScatteringEvent& evt) const {
        if (evt.triggers[5] != true) {
          return false;
        }

        return true;
      }

      bool IsProtonInAcceptanceCut(const ScatteringEvent& evt) const {
        if (evt.p_theta_eff < 55 || evt.p_theta_eff > 70) {
          return false;
        }

        return true;
      }

      bool IsLeftEspriEvent(const ScatteringEvent& evt) const {
        return evt.esl_de_raw > 240;
      }

      bool IsRightEspriEvent(const ScatteringEvent& evt) const {
        return evt.esr_de_raw > 240;
      }

      bool IsEspriEvent(const ScatteringEvent& evt) const {
        return IsLeftEspriEvent(evt) || IsRightEspriEvent(evt);
      }

      bool CheckCommonCuts(const ScatteringEvent& evt) const {

        // select by ESPRI left/right
        switch (espri_selector) {
        case Espri::both:
          break;
        case Espri::left:
          if (!IsLeftEspriEvent(evt)) {
            return false;
          }
          break;
        case Espri::right:
          if (!IsRightEspriEvent(evt)) {
            return false;
          }
          break;
        }

        // trigger conditions
        if (apply_reaction_trigger_cut && !IsReactionTriggerCut(evt)) {
          return false;
        }
        if (apply_beam_trigger_cut && !IsBeamTriggerCut(evt)) {
          return false;
        }

        // angular region check
        // if (!IsProtonInAcceptanceCut(evt)) {
        //   return false;
        // }

        // beam ab angles cut
        if (apply_beam_abang_cut && !IsBeamAbangCut(evt)) {
          return false;
        }

        // target cut
        if (apply_vertexXY_cut && !IsVertexXYCut(evt)) {
          return false;
        }
        // target cut - elliptical
        if (apply_vertexXY_elliptical_cut && !IsVertexXYEllipticalCut(evt)) {
          return false;
        }

        // theta angle correlation cut
        if (apply_theta_cut && !IsThetaCut(evt)) {
          return false;
        }

        // min E & dE cuts
        if (apply_espri_min_e_de_cut && !IsEspriMinEdECut(evt)) {
          return false;
        }

        // E cut using ESPRI NaIs calibrated data
        if (apply_espri_e_cut && !IsEspriECut(evt)) {
          return false;
        }

        // ab angle cut from RDC data
        if (apply_espri_aang_cut && !IsEspriAangCut(evt)) {
          return false;
        }

        // ab angle cut from RDC data
        if (apply_espri_vertz_cut && !IsEspriVertzCut(evt)) {
          return false;
        }

        // hodf cuts
        if (apply_hodf_he6_cut && !IsHodfHe6Cut(evt)) {
          return false;
        }
        if (apply_hodf_he4_cut && !IsHodfHe4Cut(evt)) {
          return false;
        }

        // espri Y cut
        if (apply_espri_y_cut && !IsEspriYCut(evt)) {
          return false;
        }

        // RDC nhit planes cut
        if (apply_rdc_xnhitplanes_cut && !IsRdcXNHitPlanesCut(evt)) {
          return false;
        }
        if (apply_rdc_ynhitplanes_cut && !IsRdcYNHitPlanesCut(evt)) {
          return false;
        }

        return true;
      }

      bool IsElasticEvent(const ScatteringEvent& evt) const {

        if (!CheckCommonCuts(evt)) {
          return false;
        }

        // phi angle cut
        if (apply_phi_cut && !IsPhiCut(evt)) {
          return false;
        }

        return true;
      }

      bool IsSignalEvent(const ScatteringEvent& evt) const {
        return IsElasticEvent(evt);
      }

      bool IsBackgroundEvent(const ScatteringEvent& evt) const {

        if (!CheckCommonCuts(evt)) {
          return false;
        }

        // phi angle cut
        if (!IsPhiBgCut(evt)) {
          return false;
        }

        return true;
      }

      bool operator() (const ScatteringEvent& evt) const {
        return IsElasticEvent(evt);
      }

      int get_hodf_hit_id(const ScatteringEvent& evt) const {
        std::vector<int> ids = get_hodf_hit_ids(evt);
        if (ids.size() > 2 || ids.size() == 0) {
          return -1;
        }

        if (is_hodf_double_hit_evt(evt)) {
          return evt.hodf_q[ids[0]] > evt.hodf_q[ids[1]] ?
                                             ids[0] : ids[1];
        } else {
          return ids[0];
        }
      }

      double get_hodf_q(const ScatteringEvent& evt) const {
        std::vector<int> ids = get_hodf_hit_ids(evt);
        if (ids.size() > 2 || ids.size() == 0) {
          return 0;
        }

        if (is_hodf_double_hit_evt(evt)) {
          return evt.hodf_q[ids[0]] + evt.hodf_q[ids[1]];
        } else {
          return evt.hodf_q[ids[0]];
        }
      }

      void PrepareVarThetaCut(DatasetType ds_type) {
        if (can_use_var_theta_cut) {
          // if already prepared no need to repeat
          return;
        }
        theta_sigma_interp_s1dc_esl = make_theta_sigma_interp(ds_type, Dsdcs::s1dc, Espri::left);
        theta_sigma_interp_s1dc_esr = make_theta_sigma_interp(ds_type, Dsdcs::s1dc, Espri::right);
        theta_mean_interp_s1dc_esl = make_theta_mean_interp(ds_type, Dsdcs::s1dc, Espri::left);
        theta_mean_interp_s1dc_esr = make_theta_mean_interp(ds_type, Dsdcs::s1dc, Espri::right);
        theta_sigma_interp_fdc0_esl = make_theta_sigma_interp(ds_type, Dsdcs::fdc0, Espri::left);
        theta_sigma_interp_fdc0_esr = make_theta_sigma_interp(ds_type, Dsdcs::fdc0, Espri::right);
        theta_mean_interp_fdc0_esl = make_theta_mean_interp(ds_type, Dsdcs::fdc0, Espri::left);
        theta_mean_interp_fdc0_esr = make_theta_mean_interp(ds_type, Dsdcs::fdc0, Espri::right);
        logger_.info("Variable theta cut was prepared for use with %s dataset",
                     ds_type == DatasetType::he4 ? "he4" : "he6");
        can_use_var_theta_cut = true;
      }

      double GetHeTheta(const ScatteringEvent& evt) const {
        switch (dsdc_selector) {
        case Dsdcs::s1dc:
          return evt.s1dc_theta;
          break;
        case Dsdcs::fdc0:
          return evt.fdc0_theta;
          break;
        default:
          throw std::logic_error("invalid DSDC selector value!");
        }
      }

      double GetHePhi(const ScatteringEvent& evt) const {
        switch (dsdc_selector) {
        case Dsdcs::s1dc:
          return evt.s1dc_phi;
          break;
        case Dsdcs::fdc0:
          return evt.fdc0_phi;
          break;
        default:
          throw std::logic_error("invalid DSDC selector value!");
        }
      }

      const EDistParamsInterp&
      GetHe6EDistSigmaInterp() const {
        return *he6_e_dist_sigma_interp;
      }

      const EDistParamsInterp&
      GetHe4EDistSigmaInterp() const {
        return *he4_e_dist_sigma_interp;
      }

      const EDistParamsInterp&
      GetEDistSigmaInterp(DatasetType ds_type) const {
        if (ds_type == DatasetType::he4) {
          return GetHe4EDistSigmaInterp();
        } else if (ds_type == DatasetType::he6) {
          return GetHe6EDistSigmaInterp();
        } else {
          throw std::invalid_argument("Invalid dataset type");
        }
      }

      const EDistParamsInterp&
      GetHe6EDistMeanInterp(Espri espri) const {
        if (espri == Espri::left) {
          return *he6_esl_e_dist_mean_interp;
        } else if (espri == Espri::right) {
          return *he6_esr_e_dist_mean_interp;
        } else {
          throw std::invalid_argument("BOTH is invalid selector here!");
        }
      }

      const EDistParamsInterp&
      GetHe4EDistMeanInterp(Espri espri) const {
        if (espri == Espri::left) {
          return *he4_esl_e_dist_mean_interp;
        } else if (espri == Espri::right) {
          return *he4_esr_e_dist_mean_interp;
        } else {
          throw std::invalid_argument("BOTH is invalid selector here!");
        }
      }

      const EDistParamsInterp&
      GetEDistMeanInterp(DatasetType ds_type, Espri espri) const {
        if (ds_type == DatasetType::he4) {
          return GetHe4EDistMeanInterp(espri);
        } else if (ds_type == DatasetType::he6) {
          return GetHe6EDistMeanInterp(espri);
        } else {
          throw std::invalid_argument("Invalid dataset type");
        }
      }

      const std::shared_ptr<ERelMeanInterp>
      GetHe6NaIMeanInterp() const {
        return he6_nai_e_rel_interp;
      }

      double ComputeCorrectedERel(const ScatteringEvent& evt) {
        Double_t cut_dist_mean = 0.0;
        if (espri_e_cut_use_he6_rel_e_corr) {
          auto& interp = IsLeftEspriEvent(evt) ?
            *he6_esl_e_dist_mean_interp : *he6_esr_e_dist_mean_interp;
          // NaI mean correction
          cut_dist_mean += he6_nai_e_rel_interp->eval(evt);
          // Overall E-dist mean correction
          cut_dist_mean += interp(evt.p_theta_eff);
        }
        double e_rel = ((evt.p_e - evt.p_e_sim) / evt.p_e_sim) - cut_dist_mean;
        return e_rel;
      }

      double GetCorrectedEspriAang(const ScatteringEvent& evt, const Espri espri) const {
        assert(espri != Espri::both);
        Double_t aang = espri == Espri::left ? evt.esl_aang : evt.esr_aang;
        if (use_espri_aang_corr) {
          auto& espri_aang_mean = espri == Espri::left ? esl_aang_mean : esr_aang_mean;
          aang -= (espri_aang_mean.eval(evt.p_theta_eff)
                   - espri_computed_aang_mean.eval(evt.p_theta_eff));
        }
        return aang;
      }

    private:

      std::vector<int> get_hodf_hit_ids(const ScatteringEvent& evt) const {
        std::vector<int> hit_ids;
        for (int i = 0; i <= 23; i++) {
          if (evt.hodf_q[i] > gk_hodf_fire_threshold_q) {
            hit_ids.push_back(i);
          }
        }
        return hit_ids;
      }

      bool is_hodf_double_hit_evt(const ScatteringEvent& evt) const {
        std::vector<int> hit_ids = get_hodf_hit_ids(evt);
        if (hit_ids.size() == 2) {
          if (std::abs(hit_ids[0] - hit_ids[1]) == 1) {
            return true;
          }
        }
        return false;
      }

      const std::shared_ptr<ROOT::Math::Interpolator>&
      GetThetaMeanInterp(const ScatteringEvent& evt) const {
        if (IsLeftEspriEvent(evt)) {
          if (dsdc_selector == Dsdcs::s1dc) {
            return theta_mean_interp_s1dc_esl;
          } else if (dsdc_selector == Dsdcs::fdc0) {
            return theta_mean_interp_fdc0_esl;
          }
        } else if (IsRightEspriEvent(evt)) {
          if (dsdc_selector == Dsdcs::s1dc) {
            return theta_mean_interp_s1dc_esr;
          } else if (dsdc_selector == Dsdcs::fdc0) {
            return theta_mean_interp_fdc0_esr;
          }
        }

        throw std::logic_error("No suitable return value!");
      }

      const std::shared_ptr<ROOT::Math::Interpolator>&
      GetThetaSigmaInterp(const ScatteringEvent& evt) const {
        if (IsLeftEspriEvent(evt)) {
          if (dsdc_selector == Dsdcs::s1dc) {
            return theta_sigma_interp_s1dc_esl;
          } else if (dsdc_selector == Dsdcs::fdc0) {
            return theta_sigma_interp_fdc0_esl;
          }
        } else if (IsRightEspriEvent(evt)) {
          if (dsdc_selector == Dsdcs::s1dc) {
            return theta_sigma_interp_s1dc_esr;
          } else if (dsdc_selector == Dsdcs::fdc0) {
            return theta_sigma_interp_fdc0_esr;
          }
        }

        throw std::logic_error("No suitable return value!");
      }

      template<typename... FmtArgTypes>
      TCut GetCut(const char* format_string, FmtArgTypes&&... fmt_args) const {
        // std::tuple<FmtArgTypes...> store( fmt_args... );
        // std::cout << "format parameter 1: " << std::get<0>(store) << std::endl;
        TString cut_str{
          TString::Format(format_string, std::forward<FmtArgTypes>(fmt_args)...)};
        return TCut{cut_str.Data()};
      }

      bool IsWellDefined(double value) const {
        if (static_cast<int>(value) == -99999) {
          return false;
        } else {
          return true;
        }
      }

      // variable phi cuts helper objects
      LinearInterp lin_phi_cut_fdc0;
      LinearInterp lin_phi_cut_s1dc;
      LinearInterp esl_aang_mean;
      LinearInterp esr_aang_mean;
      LinearInterp espri_computed_aang_mean;

      // variable theta cut helper objects
      std::shared_ptr<ROOT::Math::Interpolator> theta_sigma_interp_s1dc_esl;
      std::shared_ptr<ROOT::Math::Interpolator> theta_sigma_interp_s1dc_esr;
      std::shared_ptr<ROOT::Math::Interpolator> theta_mean_interp_s1dc_esl;
      std::shared_ptr<ROOT::Math::Interpolator> theta_mean_interp_s1dc_esr;
      std::shared_ptr<ROOT::Math::Interpolator> theta_sigma_interp_fdc0_esl;
      std::shared_ptr<ROOT::Math::Interpolator> theta_sigma_interp_fdc0_esr;
      std::shared_ptr<ROOT::Math::Interpolator> theta_mean_interp_fdc0_esl;
      std::shared_ptr<ROOT::Math::Interpolator> theta_mean_interp_fdc0_esr;

      std::shared_ptr<EDistParamsInterp> he6_e_dist_sigma_interp;
      std::shared_ptr<EDistParamsInterp> he4_e_dist_sigma_interp;
      std::shared_ptr<EDistParamsInterp> he6_esl_e_dist_mean_interp;
      std::shared_ptr<EDistParamsInterp> he6_esr_e_dist_mean_interp;
      std::shared_ptr<EDistParamsInterp> he4_esl_e_dist_mean_interp;
      std::shared_ptr<EDistParamsInterp> he4_esr_e_dist_mean_interp;
      std::shared_ptr<ERelMeanInterp> he6_nai_e_rel_interp = make_he6_nai_erel_mean_interp();
      std::shared_ptr<ERelMeanInterp> he4_nai_e_rel_interp = make_he4_nai_erel_mean_interp();

      s13::misc::MessageLogger logger_;
    };

  }
}
