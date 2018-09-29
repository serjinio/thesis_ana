
#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TCut.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TVector3.h>

//#include "anapow.hpp"
#include "treewalker.hpp"
#include "scattyield.hpp"
#include "scattevt.hpp"

#include "elacuts.hpp"
#include "common.hpp"


namespace s13 {
  namespace ana {
    namespace algs {

      using asym_yields = std::pair<ScatteringYield, ScatteringYield>;

      static auto empty_evt_handler_fn = [] (ScatteringEvent&) {};

      // class ScatteringAsymmetry {
      // public:
      //   ScatteringAsymmetry(ScatteringYield yields,
      //                       ScatteringYield rel_stat_error) :
      //     m_target_pol{0}, m_yields{yields},
      //     m_rel_stat_error{rel_stat_error} {}

      //   ScatteringYield yields() {
      //     return m_yields;
      //   }

      //   ScatteringYield relative_stat_error() {
      //     return rel_stat_error;
      //   }

      // private:
      //   asym_yields m_yields;
      //   asym_yields m_rel_stat_error;
      //   Double_t m_target_pol;
      // };

      /*
        Computes left & right scattering yields for p-4He with BG subtracted.
      */
      template <typename EvtHandlerFn = decltype(empty_evt_handler_fn)>
      asym_yields
      He4LeftRightScatteringYields(TTree& ttree,
                                   const ElasticScatteringCuts& cuts,
                                   int bin_num = 5,
                                   EvtHandlerFn bg_evt_fn = empty_evt_handler_fn,
                                   EvtHandlerFn signal_evt_fn = empty_evt_handler_fn) {
        static const Double_t angular_region[] {55, 68};

        ScatteringYield left_yield{bin_num,
            angular_region[0], angular_region[1]};
        ScatteringYield right_yield{bin_num,
            angular_region[0], angular_region[1]};
        ScatteringYield left_bg_yield{bin_num,
            angular_region[0], angular_region[1]};
        ScatteringYield right_bg_yield{bin_num,
            angular_region[0], angular_region[1]};
        ScatteringTreeWalker tree_walker;
        int evt_counter = 0;
        int bg_evt_counter = 0;

        tree_walker.Walk
          (ttree,
           [&] (ScatteringEvent& evt) {
            if (cuts.IsBackgroundEvent(evt)) {
              bg_evt_fn(evt);

              if (evt.IsLeftScattering()) {
                left_bg_yield.AddEvent(evt.p_theta_eff);
              } else {
                right_bg_yield.AddEvent(evt.p_theta_eff);
              }

              ++bg_evt_counter;
            }

            if (cuts.IsElasticEvent(evt)) {
              signal_evt_fn(evt);

              if (evt.IsLeftScattering()) {
                left_yield.AddEvent(evt.p_theta_eff);
              } else {
                right_yield.AddEvent(evt.p_theta_eff);
              }

              ++evt_counter;
            }
          });

        auto bg_norm_fact_interp = MakeBgNormFactorInterpolator();
        auto left = He4SubtractBg(left_yield, left_bg_yield, *bg_norm_fact_interp);
        auto right = He4SubtractBg(right_yield, right_bg_yield, *bg_norm_fact_interp);

        return std::make_pair(left, right);
      }

      /*
        Computes relative statistical error for scattering asymmetry.
        Currently not used.
       */
      ScatteringYield ScatteringAsymmetryRelStatError(asym_yields pol_up_yields,
                                                      asym_yields pol_down_yields) {
        auto left_up = std::get<0>(pol_up_yields);
        auto right_up = std::get<1>(pol_up_yields);
        auto left_down = std::get<0>(pol_down_yields);
        auto right_down = std::get<1>(pol_down_yields);

        auto fact1_num = (left_up * right_down * right_up * left_down).sqrt();
        auto fact1_den = left_up * right_down - right_up * left_down;
        auto fact1 = fact1_num / fact1_den;

        // compute 1/Y_l_up + 1/Y_r_dwn + 1/Y_r_up + 1/Y_l_dwn
        auto fact2 = (left_up.invert() + right_down.invert()
                      + right_up.invert() + left_down.invert()).sqrt();

        return fact1 * fact2;
      }

      /*
        Computes scattering asymmetry.
      */
      template <typename EvtHandlerFn = decltype(empty_evt_handler_fn)>
      ScatteringYield
      He4ScatteringAsymmetry(TTree& ttree_pol_up, TTree& ttree_pol_down,
                             const ElasticScatteringCuts& cuts,
                             int bin_num = 5,
                             EvtHandlerFn bg_evt_fn = empty_evt_handler_fn,
                             EvtHandlerFn signal_evt_fn = empty_evt_handler_fn) {
        auto left_right_up =
          He4LeftRightScatteringYields(ttree_pol_up, cuts,
                                       bin_num, bg_evt_fn,
                                       signal_evt_fn);
        auto left_right_down =
          He4LeftRightScatteringYields(ttree_pol_down, cuts,
                                       bin_num, bg_evt_fn,
                                       signal_evt_fn);
        auto left_up = std::get<0>(left_right_up);
        auto right_up = std::get<1>(left_right_up);
        auto left_down = std::get<0>(left_right_down);
        auto right_down = std::get<1>(left_right_down);

        std::cout << "[DEBUG] Pol. up yields: " << std::endl;
        std::cout << "\t Left: " << left_up << std::endl;
        std::cout << "\t Right: " << right_up << std::endl;
        std::cout << "[DEBUG] Pol. down yields: " << std::endl;
        std::cout << "\t Left: " << left_down << std::endl;
        std::cout << "\t Right: " << right_down << std::endl;
        std::cout << "\nPol up vs. down yields difference: " << std::endl;
        std::cout << "\tLeft up - left down: "
                  << left_up.TotalEventCount() - left_down.TotalEventCount()
                  << std::endl;
        std::cout << "\tRight up - right down: "
                  << right_up.TotalEventCount() - right_down.TotalEventCount()
                  << std::endl;

        auto asymmetry_num = (left_up * right_down).sqrt()
          - (right_up * left_down).sqrt();
        auto asymmetry_den = (left_up * right_down).sqrt()
          + (right_up * left_down).sqrt();
        auto asymmetry = asymmetry_num / asymmetry_den;
        // auto rel_stat_error = ScatteringAsymmetryRelStatError(left_right_up,
        //                                                       left_right_down);
        //ScatteringAsymmetry asym(asymmetry, rel_stat_error);
        return asymmetry;
      }


      class CrossSectionUtils {
      public:

        using range = std::pair<Double_t, Double_t>;

        static constexpr Double_t k_rdc_center_angle_deg = 62.5;
        static constexpr Double_t k_dist_to_rdc = 1004.75;
        static constexpr Double_t k_tgt_areal_number_density = 1.046e22;

      public:
        CrossSectionUtils(Double_t total_num_incident, ElasticScatteringCuts cuts,
                          int bin_num = 32) :
            m_total_num_incident{total_num_incident},
            m_bin_num{bin_num},
            m_min_theta{55}, m_max_theta{70},
            m_rdc_min_y{-132}, m_rdc_max_y{132}, m_cuts{cuts} {}

        CrossSectionUtils() : CrossSectionUtils{1, ElasticScatteringCuts()} {}

        int bin_num() const {
          return m_bin_num;
        }

        void bin_num(Double_t new_val) {
          m_bin_num = new_val;
        }

        range theta_range_deg() const {
          return range(m_min_theta, m_max_theta);
        }

        void theta_range_deg(range&& new_range) {
          set_from_range(&m_min_theta, &m_max_theta, std::move(new_range));
        }

        range theta_range_deg_cm() const {
          return range(lab_theta_angle_to_cm(m_max_theta * D2R) / D2R,
                       lab_theta_angle_to_cm(m_min_theta * D2R) / D2R);
        }

        Double_t bin_width_deg() const {
          return std::abs(m_max_theta - m_min_theta) / m_bin_num;
        }

        Double_t bin_width_rad() const {
          return bin_width_deg() * D2R;
        }

        /*
          Returns bin width in CM frame.
        */
        Double_t bin_width_rad_cm() const {
          Double_t min_theta_cm = lab_theta_angle_to_cm(m_min_theta);
          Double_t max_theta_cm =
            lab_theta_angle_to_cm(m_min_theta + bin_width_deg());
          return std::abs(max_theta_cm - min_theta_cm) * D2R;
        }

        range ypos_range_mm() const {
          return range(m_rdc_min_y, m_rdc_max_y);
        }

        void ypos_range_mm(range&& new_range) {
          set_from_range(&m_rdc_min_y, &m_rdc_max_y, std::move(new_range));
        }

        ElasticScatteringCuts& cuts() {
          return m_cuts;
        }

        void cuts(ElasticScatteringCuts cuts) {
          m_cuts = cuts;
        }

        Double_t lab_theta_angle_to_cm(Double_t theta_lab_rad) const {
          return TMath::Pi() - 2 * theta_lab_rad;
        }

        Double_t cm_theta_angle_to_lab(Double_t theta_cm_rad) const {
          return (TMath::Pi() - theta_cm_rad) / 2;
        }

        bool is_event(const ScatteringEvent& evt) const {
          Double_t ypos = evt.IsLeftScattering() ? evt.esl_ypos : evt.esr_ypos;
          if (ypos < ypos_range_mm().first ||
              ypos > ypos_range_mm().second) {
            return false;
          }

          return m_cuts.IsElasticEvent(evt);
        }

        bool is_bg_event(const ScatteringEvent& evt) const {
          Double_t ypos = evt.IsLeftScattering() ? evt.esl_ypos : evt.esr_ypos;
          if (ypos < ypos_range_mm().first ||
              ypos > ypos_range_mm().second) {
            return false;
          }

          return m_cuts.IsBackgroundEvent(evt);
        }

        Double_t compute_bin_solid_angle(Double_t bin_theta_center_angle_deg) const {
          Double_t ksi_angle = abs(k_rdc_center_angle_deg - bin_theta_center_angle_deg);
          Double_t dist_to_bin_center = k_dist_to_rdc / TMath::Cos(ksi_angle * D2R);
          Double_t rdc_segment_area = compute_rdc_segment_area(bin_theta_center_angle_deg);
          Double_t solid_angle = rdc_segment_area * TMath::Cos(ksi_angle * D2R)
                                 / std::pow(dist_to_bin_center, 2);
          // account for two ESPRIs left & right
          if (m_cuts.espri_selector == Espri::both) {
            solid_angle *= 2;
          }

          std::cout << "[DEBUG] ksi_angle: " << ksi_angle
                    << "; dist to bin center [mm]: " << dist_to_bin_center
                    << "; RDC segment area [mm2]: " << rdc_segment_area
                    << "; solid angle [sr]: " << solid_angle << std::endl;
          return solid_angle;
        }

        Double_t compute_bin_solid_angle_2(Double_t bin_theta_center_angle_deg) const {
          std::cout << "computing bin solid angle (2) for theta: "
                    << bin_theta_center_angle_deg << std::endl;

          TVector3 phi1_vec(0, m_rdc_min_y, k_dist_to_rdc);
          phi1_vec.RotateY(bin_theta_center_angle_deg * D2R);
          TVector3 phi2_vec(0, m_rdc_max_y, k_dist_to_rdc);
          phi2_vec.RotateY(bin_theta_center_angle_deg * D2R);

          // we use phi defined from Y axis, use -90 deg.
          // to convert to our definition
          Double_t delta_phi = std::abs((phi1_vec.Phi() - 90*D2R) -
                                        (phi2_vec.Phi() - 90*D2R));
          Double_t solid_angle = TMath::Sin(bin_theta_center_angle_deg * D2R)
            * bin_width_rad() * delta_phi;

          // account for two ESPRIs left & right
          if (m_cuts.espri_selector == Espri::both) {
            solid_angle *= 2;
          }

          std::cout << "[DEBUG] theta_width [rad]: " << bin_width_rad()
                    << "; phi width [rad]: " << delta_phi << std::endl;

          std::cout << "[DEBUG] phi1 angle [deg]: "
                    << (phi2_vec.Phi() - 90*D2R) * 1/D2R
                    << "; phi2 angle [deg]: "
                    << (phi1_vec.Phi() - 90*D2R) * 1/D2R
                    << "; solid angle (2) [sr]: "
                    << solid_angle << std::endl;
          return solid_angle;
        }

        Double_t compute_bin_solid_angle_cm(Double_t bin_theta_center_angle_cm_deg) const {
          std::cout << "computing bin solid angle (CM) for CM theta: "
                    << bin_theta_center_angle_cm_deg << std::endl;
          // convert to lab angles to define radius vectors for the bin
          Double_t bin_theta_ctrang_lab_deg =
            cm_theta_angle_to_lab(bin_theta_center_angle_cm_deg * D2R) / D2R;
          TVector3 phi1_vec(0, m_rdc_min_y, k_dist_to_rdc);
          phi1_vec.RotateY(bin_theta_ctrang_lab_deg * D2R);
          TVector3 phi2_vec(0, m_rdc_max_y, k_dist_to_rdc);
          phi2_vec.RotateY(bin_theta_ctrang_lab_deg * D2R);

          // we use phi defined from Y axis, use -90 deg.
          // to convert to our definition
          Double_t delta_phi = std::abs((phi1_vec.Phi() - 90*D2R) -
                                        (phi2_vec.Phi() - 90*D2R));
          Double_t solid_angle = TMath::Sin(bin_theta_center_angle_cm_deg * D2R)
            * bin_width_rad_cm() * delta_phi;

          // account for two ESPRIs left & right
          if (m_cuts.espri_selector == Espri::both) {
            solid_angle *= 2;
          }

          std::cout << "[DEBUG] CM theta_width [rad]: " << bin_width_rad_cm()
                    << "; phi width [rad]: " << delta_phi << std::endl;

          std::cout << "[DEBUG] phi1 angle [deg]: "
                    << (phi2_vec.Phi() - 90*D2R) * 1/D2R
                    << "; phi2 angle [deg]: "
                    << (phi1_vec.Phi() - 90*D2R) * 1/D2R
                    << "; solid angle (CM) [sr]: "
                    << solid_angle << std::endl;
          return solid_angle;
        }

        Double_t compute_rdc_segment_area(Double_t bin_theta_center_angle_deg) const {
          Double_t ksi_angle = std::abs(k_rdc_center_angle_deg - bin_theta_center_angle_deg);
          Double_t bin_min_angle = ksi_angle - bin_width_deg() / 2;
          Double_t bin_max_angle = ksi_angle + bin_width_deg() / 2;
          Double_t rdc_min_x = k_dist_to_rdc * TMath::Tan(bin_min_angle * D2R);
          Double_t rdc_max_x = k_dist_to_rdc * TMath::Tan(bin_max_angle * D2R);
          std::cout << "rdc_min_x: " << rdc_min_x << "; rdc_max_x :"
                    << rdc_max_x << std::endl;
          // returns area of the RDC detector segment in mm2
          return (rdc_max_x - rdc_min_x) * (m_rdc_max_y - m_rdc_min_y);
        }

        Double_t bin_norm_factor(Double_t bin_theta_center_angle_deg) const {
          Double_t solid_angle = compute_bin_solid_angle(bin_theta_center_angle_deg);
          Double_t solid_angle_2 = compute_bin_solid_angle_2(bin_theta_center_angle_deg);
          std::cout << "difference: solid_angle - solid_angle2: "
                    << solid_angle - solid_angle_2 << std::endl;
          Double_t norm_factor = 1 / (m_total_num_incident * k_tgt_areal_number_density
                                      * solid_angle);
          // from cm2 to barns
          norm_factor *= 1e24;
          // now into millibarns
          norm_factor *= 1e3;

          std::cout << "[DEBUG] Norm factor for bin at theta [deg]: "
                    << bin_theta_center_angle_deg << ": "
                    << norm_factor << "." << std::endl;
          return norm_factor;
        }

        Double_t cm_bin_norm_factor(Double_t bin_theta_center_angle_cm_deg) const {
          Double_t solid_angle =
            compute_bin_solid_angle_cm(bin_theta_center_angle_cm_deg);
          Double_t norm_factor = 1 / (m_total_num_incident * k_tgt_areal_number_density
                                      * solid_angle);
          // from cm2 to barns
          norm_factor *= 1e24;
          // now into millibarns
          norm_factor *= 1e3;

          std::cout << "[DEBUG] Norm factor for bin at theta [CM deg]: "
                    << bin_theta_center_angle_cm_deg << ": "
                    << norm_factor << "." << std::endl;
          return norm_factor;
        }

        void normalize_yields(ScatteringYield& yields) const {
          for (int i = 0; i < yields.bins_num; i++) {
            Double_t norm_factor = bin_norm_factor(yields.GetBinArg(i));
            yields.ScaleBin(i, norm_factor);
          }
        }

        void normalize_cm_yields(ScatteringYield& yields) const {
          for (int i = 0; i < yields.bins_num; i++) {
            Double_t norm_factor = cm_bin_norm_factor(yields.GetBinArg(i));
            yields.ScaleBin(i, norm_factor);
          }
        }

        ScatteringYield lab_cs_to_cm(const ScatteringYield& lab_cs) {
          const double lab_theta_min =
            (lab_cs.GetBinArg(0) - lab_cs.GetBinWidth() / 2) * D2R;
          const double lab_theta_max =
            (lab_cs.GetBinArg(lab_cs.GetBinsNumber() - 1)
             + lab_cs.GetBinWidth() / 2) * D2R;
          const double cm_theta_max = lab_theta_angle_to_cm(lab_theta_min);
          const double cm_theta_min = lab_theta_angle_to_cm(lab_theta_max);

          std::cout << "CM angular range: [" << cm_theta_min * 1/D2R
                    << "; " << cm_theta_max * 1/D2R << "]" << std::endl;

          ScatteringYield cm_cs(lab_cs.bins_num,
                                cm_theta_min * 1/D2R, cm_theta_max * 1/D2R);

          for (int i = cm_cs.bins_num - 1; i >= 0; --i) {
            int cm_cs_idx = cm_cs.bins_num - 1 - i;
            Double_t cm_angle = lab_theta_angle_to_cm(lab_cs.GetBinArg(i) * D2R);
            std::cout << "Computing CM CS for angle: " << cm_angle * 1/D2R
                      << "; actual angle in object: " << cm_cs.GetBinArg(cm_cs_idx)
                      << std::endl;

            Double_t solid_angle_factor = 1 / (4 * TMath::Sin(cm_angle / 2));
            Double_t cm_cs_value = lab_cs.GetBinValue(i) * solid_angle_factor;
            std::cout << "CS value in lab frame: " << lab_cs.GetBinValue(i)
                      << "; in CM frame: " << cm_cs_value << std::endl;

            cm_cs.SetBinValue(cm_angle * 1/D2R, cm_cs_value);
            cm_cs.SetBinError(cm_angle * 1/D2R, lab_cs.yserr[i] * solid_angle_factor);
          }

          return cm_cs;
        }

      private:

        void set_from_range(Double_t* min, Double_t* max, range&& from_range) {
          *min = std::get<0>(from_range);
          *max = std::get<1>(from_range);
        }

        Double_t m_total_num_incident;
        int m_bin_num;
        Double_t m_min_theta;
        Double_t m_max_theta;
        Double_t m_rdc_min_y;
        Double_t m_rdc_max_y;
        ElasticScatteringCuts m_cuts;
      };

      struct CrossSectionRetVal {

        CrossSectionRetVal(Double_t theta_min, Double_t theta_max, int bin_num) :
            yield_signal_total{bin_num, theta_min, theta_max},
            yield_bg_total{bin_num, theta_min, theta_max},
            yield_signal_wo_bg{bin_num, theta_min, theta_max},
            lab_cs{bin_num, theta_min, theta_max},
            cm_cs{bin_num, theta_min, theta_max},
            name{""}
        {}

        std::vector<TH1*> hist_signals;
        std::vector<TH1*> hist_bgs;
        ScatteringYield yield_signal_total;
        ScatteringYield yield_bg_total;
        ScatteringYield yield_signal_wo_bg;
        ScatteringYield lab_cs;
        ScatteringYield cm_cs;
        CrossSectionUtils utils;
        TString name;
      };

      template <typename EvtHandlerFn = decltype(empty_evt_handler_fn)>
      CrossSectionRetVal
      CrossSection(TTree& ttree,
                   CrossSectionUtils& cs_utils,
                   EvtHandlerFn bg_evt_fn = empty_evt_handler_fn,
                   EvtHandlerFn signal_evt_fn = empty_evt_handler_fn) {
        ScatteringYield yields{cs_utils.bin_num(),
                               std::get<0>(cs_utils.theta_range_deg()),
                               std::get<1>(cs_utils.theta_range_deg())};
        ScatteringYield bg_yields{yields};

        // for debugging purposes - to monitor yields distribution into each bin
        std::vector<TH1*> hist_signals;
        std::vector<TH1*> hist_bgs;
        for (int i = 0; i < cs_utils.bin_num(); ++i) {
          Double_t bin_mean_theta = yields.GetBinArg(i);
          TString signal_hist_name = TString::Format("sig%i", i);
          TString signal_hist_title = TString::Format("He4 IN #phi @%.2f lab. deg.",
                                                      bin_mean_theta);
          TString bg_hist_name = TString::Format("bg%i", i);
          TString bg_hist_title = TString::Format("He4 BG @%.2f lab. deg.",
                                                  bin_mean_theta);

          hist_signals.push_back(new TH1F(signal_hist_name.Data(),
                                          signal_hist_title.Data(),
                                          48, -8, 8));
          hist_bgs.push_back(new TH1F(bg_hist_name.Data(),
                                      bg_hist_title.Data(),
                                      48, -8, 8));
        }

        ScatteringTreeWalker tree_walker;
        tree_walker.Walk(ttree, [&] (ScatteringEvent& evt) {
          if (cs_utils.is_bg_event(evt)) {
            bg_evt_fn(evt);
            bg_yields.AddEvent(evt.p_theta_eff);
            int bin_no = bg_yields.GetBinNo(evt.p_theta_eff);
            hist_bgs[bin_no]->Fill(evt.he_theta - evt.he_theta_theor);
          }

          if (cs_utils.is_event(evt)) {
            signal_evt_fn(evt);
            yields.AddEvent(evt.p_theta_eff);
            int bin_no = yields.GetBinNo(evt.p_theta_eff);
            hist_signals[bin_no]->Fill(evt.he_theta - evt.he_theta_theor);
          }
        });

        auto bg_norm_fact_interp = MakeBgNormFactorInterpolator();
        auto yields_wo_bg = He4SubtractBg(yields, bg_yields, *bg_norm_fact_interp);
        //auto yields_wo_bg = yields;
        std::cout << "Yields with BG: " << yields;
        std::cout << "Total BG yields: " << bg_yields;
        std::cout << "Yields w/o BG: " << yields_wo_bg;

        CrossSectionRetVal retval(std::get<0>(cs_utils.theta_range_deg()),
                                  std::get<1>(cs_utils.theta_range_deg()),
                                  cs_utils.bin_num());
        retval.utils = cs_utils;
        retval.hist_bgs = hist_bgs;
        retval.hist_signals = hist_signals;
        retval.yield_signal_wo_bg = yields_wo_bg;
        retval.yield_signal_total = yields;
        retval.yield_bg_total = bg_yields;

        cs_utils.normalize_yields(yields_wo_bg);
        retval.lab_cs = yields_wo_bg;
        retval.cm_cs = cs_utils.lab_cs_to_cm(yields_wo_bg);

        return retval;
      }

      template <typename EvtHandlerFn = decltype(empty_evt_handler_fn)>
      CrossSectionRetVal
      CmCrossSection(TTree& ttree,
                     CrossSectionUtils& cs_utils,
                     EvtHandlerFn bg_evt_fn = empty_evt_handler_fn,
                     EvtHandlerFn signal_evt_fn = empty_evt_handler_fn) {
        ScatteringYield yields{cs_utils.bin_num(),
                               std::get<0>(cs_utils.theta_range_deg_cm()),
                               std::get<1>(cs_utils.theta_range_deg_cm())};
        ScatteringYield bg_yields{yields};

        // for debugging purposes - to monitor yields distribution into each bin
        std::vector<TH1*> hist_signals;
        std::vector<TH1*> hist_bgs;
        for (int i = 0; i < cs_utils.bin_num(); ++i) {
          Double_t bin_mean_theta = yields.GetBinArg(i);
          TString signal_hist_name = TString::Format("sig%i", i);
          TString signal_hist_title = TString::Format("He4 IN #phi @%.2f lab. deg.",
                                                      bin_mean_theta);
          TString bg_hist_name = TString::Format("bg%i", i);
          TString bg_hist_title = TString::Format("He4 BG @%.2f lab. deg.",
                                                  bin_mean_theta);

          hist_signals.push_back(new TH1F(signal_hist_name.Data(),
                                          signal_hist_title.Data(),
                                          48, -8, 8));
          hist_bgs.push_back(new TH1F(bg_hist_name.Data(),
                                      bg_hist_title.Data(),
                                      48, -8, 8));
        }

        ScatteringTreeWalker tree_walker;
        tree_walker.Walk(ttree, [&] (ScatteringEvent& evt) {
            Double_t p_theta_eff_cm =
              cs_utils.lab_theta_angle_to_cm(evt.p_theta_eff * D2R) / D2R;
            if (cs_utils.is_bg_event(evt)) {
              bg_evt_fn(evt);
              bg_yields.AddEvent(p_theta_eff_cm);
              int bin_no =
                bg_yields.GetBinNo(p_theta_eff_cm);
              hist_bgs[bin_no]->Fill(evt.he_theta - evt.he_theta_theor);
            }

            if (cs_utils.is_event(evt)) {
              signal_evt_fn(evt);
              yields.AddEvent(p_theta_eff_cm);
              int bin_no =
                yields.GetBinNo(p_theta_eff_cm);
              hist_signals[bin_no]->Fill(evt.he_theta - evt.he_theta_theor);
            }
          });

        auto bg_norm_fact_interp = MakeBgNormFactorInterpolator();
        auto yields_wo_bg = He4SubtractBgCm(yields, bg_yields, *bg_norm_fact_interp);
        //auto yields_wo_bg = yields;
        std::cout << "Yields with BG: " << yields;
        std::cout << "Total BG yields: " << bg_yields;
        std::cout << "Yields w/o BG: " << yields_wo_bg;

        CrossSectionRetVal retval(std::get<0>(cs_utils.theta_range_deg()),
                                  std::get<1>(cs_utils.theta_range_deg()),
                                  cs_utils.bin_num());
        retval.utils = cs_utils;
        retval.hist_bgs = hist_bgs;
        retval.hist_signals = hist_signals;
        retval.yield_signal_wo_bg = yields_wo_bg;
        retval.yield_signal_total = yields;
        retval.yield_bg_total = bg_yields;

        cs_utils.normalize_cm_yields(yields_wo_bg);
        //retval.lab_cs = yields_wo_bg;
        retval.cm_cs = yields_wo_bg;

        return retval;
      }

    }
  }
}
