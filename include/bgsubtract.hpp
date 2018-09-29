/**
 *   \file bgsubtract.hpp
 *   \brief BG subtraction algorithms
 *
 */

#pragma once

#include <vector>

#include "Math/Interpolator.h"

#include "common.hpp"
#include "commonalg.hpp"
#include "msgoutput.hpp"
#include "scattyield.hpp"
#include "elacuts.hpp"
#include "init_ds.hpp"


namespace s13 {
  namespace ana {

    using InterpPtr = std::unique_ptr<ROOT::Math::Interpolator>;

    enum class BgExtractionMethod { out_phi, he6_carbon_data };


    class OutPhiBgFactor {
    private:
      std::array<double, 2> a_consts_ = { {1, 1} };
      std::array<double, 2> b_consts_ = { {1, 1} };
      double c_const_ = 1;

    public:
      OutPhiBgFactor(std::array<double, 2> a_consts,
                     std::array<double, 2> b_consts,
                     double c_const) :
        a_consts_{a_consts}, b_consts_{b_consts}, c_const_{c_const} {}

      double eval(double proton_theta, double phi_width) {
        double r = c_const_;
        r += std::pow(phi_width, 2) * b_consts_[1] + phi_width * b_consts_[0];
        r += std::pow(proton_theta, 2) * a_consts_[1] +
          proton_theta * a_consts_[0];
        return r;
      }
    };

    class BgExtractor {

    public:
      BgExtractor() {};
      virtual ~BgExtractor() {};

      virtual ScatteringYield
      subtract_bg(const ScatteringYield& total_yield,
                  const ScatteringYield& yield_bg_out) = 0;
    };

    class OutPhiBgExtractor : public BgExtractor {

    public:

      OutPhiBgExtractor(InterpPtr norm_fact_interpolator) :
        norm_fact_interpolator_{std::move(norm_fact_interpolator)},
        log_{"OutPhiBgExtractor"} {}

      virtual ScatteringYield
      subtract_bg(const ScatteringYield& total_yield,
                  const ScatteringYield& yield_bg_out) {
        ScatteringYield yield_bg_in(total_yield.GetBinsNumber(),
                                    total_yield.RangeStart(), total_yield.RangeEnd());

        log_.debug("Extracting BG contribution...");
        for (int i = 0; i < total_yield.GetBinsNumber(); i++) {
          Double_t p_theta = total_yield.GetBinArg(i);
          Double_t bin_yield_bg_out = yield_bg_out.GetBinValue(i);
          Double_t norm_fact = norm_fact_interpolator_->Eval(p_theta);
          Double_t bin_yield_bg_in = bin_yield_bg_out * (norm_fact - 1);

          log_.debug("\tProton theta: %.2f", p_theta);
          log_.debug("\tNormalization factor: %.2f", norm_fact);
          log_.debug("\tbin_yield_bg_in: %.2f", bin_yield_bg_in);

          yield_bg_in.SetBinValue(p_theta, bin_yield_bg_in);
        }

        return total_yield - yield_bg_in;
      }

    private:
      InterpPtr norm_fact_interpolator_;

      misc::MessageLogger log_;
    };

    OutPhiBgFactor
    make_bg_norm_factor_interp2(DatasetType ds_type,
                                Dsdcs dsdc_selector);

    /*
      Converts given bg_yields (out phi cut) to in phi using conversion factor.
    */
    ScatteringYield
    convert_to_in_phi_bg(DatasetType ds_type,
                         const ElasticScatteringCuts& cuts,
                         ScatteringYield bg_yields);

    /*
      Computes carbon yields angular distributions.
      As result provides an array of 1d kin. corrs.
      To use in computation of final spectra with subtracted
      carbon yields.
    */
    class CarbonYields : public s13::ana::SimpleTTreeAlgorithmBase {

    public:

      // type aliases
      using Evt = ScatteringEvent;
      using ElaCuts = ElasticScatteringCuts;

      // cut proc alias
      using CutProc = bool (ElaCuts::*)(const Evt& evt) const;
      // in & out components analysis on ang dependence
      using YieldData = std::vector< std::vector<TH1*> >;

      YieldData yields_in;
      YieldData yields_out;
      YieldData yields_tot;

      bool scale_carbon_bg = true;

    private:

      s13::misc::MessageLogger logger_;

      std::vector<ElaCuts> cuts_range_;
      CutProc cut_proc_;

      size_t bin_num_ = gk_cs_bin_number;
      double min_theta_ = gk_cs_range_start;
      double max_theta_ = gk_cs_range_end;

    public:

      CarbonYields(const std::vector<ElaCuts>& cuts_range,
                   CutProc cut_proc,
                   size_t bin_num = gk_cs_bin_number,
                   double range_start = gk_cs_range_start,
                   double range_end = gk_cs_range_end) :
        logger_{"CarbonYields"}, cuts_range_{cuts_range}, cut_proc_{cut_proc},
        bin_num_{bin_num}, min_theta_{range_start}, max_theta_{range_end} {
          init_hists();
        }

      CarbonYields(const CarbonYields&) = delete;
      CarbonYields(CarbonYields&&) = delete;
      CarbonYields operator=(const CarbonYields&) = delete;
      CarbonYields operator=(CarbonYields&&) = delete;

      virtual void setup() {
      }

      virtual void process_event(Evt& evt) {
        assert(cuts_range_.size() > 0);
        if (cuts_range_[0].IsSignalEvent(evt)) {
          fill_carbon_yield_hists(evt);
        }
      }

      virtual void finalize() {
        if (scale_carbon_bg) {
          for (size_t i = 0; i < yields_in.size(); i++) {
            auto& v_in = yields_in[i];
            auto& v_out = yields_out[i];
            auto& v_tot = yields_tot[i];
            for (size_t c = 0; c < v_in.size(); c++) {
              v_in[c]->Scale(s13::gk_carbon_bg_norm_factor);
              v_out[c]->Scale(s13::gk_carbon_bg_norm_factor);
              v_tot[c]->Scale(s13::gk_carbon_bg_norm_factor);
            }
          }
        }
      }

      const std::vector<ElaCuts>& cuts_range() {
        return cuts_range_;
      }

      size_t bin_num() {
        return bin_num_;
      }

      double bin_width() {
        return (max_theta_ - min_theta_) / bin_num_;
      }

      double bin_min_theta(int idx) {
        double minth = min_theta_ + bin_width() * idx;
        return minth;
      }

      double bin_max_theta(int idx) {
        double maxth = min_theta_ + bin_width() * (idx + 1);
        return maxth;
      }

      double bin_center_theta(int idx) {
        return (bin_min_theta(idx) + bin_max_theta(idx)) / 2;
      }

      int bin_idx(Evt& evt) {
        if (evt.p_theta_eff - min_theta_ < 0) {
          return -1;
        }

        size_t idx = (evt.p_theta_eff - min_theta_) / bin_width();
        if (idx < bin_num_) {
          return idx;
        } else {
          return -1;
        }
      }

    private:

      TH1* init_kin_corr_hist(TString title) {
        auto th1 = &s13::ana::make_th1;
        auto hist = th1(80, -8, 8, title, "#theta_{exp} - #theta_{theor}", "Counts");
        return hist;
      }

      void init_hists() {
        for (size_t i = 0; i < bin_num_; i++) {
          yields_in.push_back(std::vector<TH1*>());
          yields_out.push_back(std::vector<TH1*>());
          yields_tot.push_back(std::vector<TH1*>());
          for (size_t c = 0; c < cuts_range_.size(); c++) {
            yields_in[i].push_back(init_kin_corr_hist("IN"));
            yields_out[i].push_back(init_kin_corr_hist("OUT"));
            yields_tot[i].push_back(init_kin_corr_hist("TOT"));
          }
        }
      }

      void fill_carbon_yield_hists(Evt& evt) {
        Double_t he_theta = cuts_range()[0].dsdc_selector == Dsdcs::s1dc
          ? evt.s1dc_theta : evt.fdc0_theta;
        double p_theta = evt.p_theta_eff;
        double he_theta_diff = he_theta - evt.he_theta_theor;

        if (!is_well_defined(p_theta)
            || !is_well_defined(he_theta) ||
            !evt.IsEspriEvent()) {
          return;
        }

        int evt_bin_idx =  bin_idx(evt);
        if (evt_bin_idx == -1) {
          return;
        }

        for (size_t c = 0; c < cuts_range_.size(); c++) {
          auto lcuts = cuts_range_[c];
          yields_tot[evt_bin_idx][c]->Fill(he_theta_diff);
          if ((lcuts.*cut_proc_)(evt)) {
            yields_in[evt_bin_idx][c]->Fill(he_theta_diff);
          } else {
            yields_out[evt_bin_idx][c]->Fill(he_theta_diff);
          }
        }
      }
    };

    using UseAbsVal = s13::ana::NamedType<bool, struct draw_abs_bg2>;
    using Degrees = s13::ana::NamedType<double, struct degrees_type_tag>;

    /*
      Computes S/N for 1d kin. corr. using -4 - -theta_width/2 region
      as a measure of BG contribution and +/-theta_width/2 as signal measure.
    */
    double
    draw_kin_corr_hist_sn(double angle, TH1* bin_hist,
                          Dsdcs dsdc_selector, double theta_width,
                          UseAbsVal use_abs = UseAbsVal(false),
                          Degrees bg_min_theta_angle = Degrees(-4));

    /*
      Class to compute angular dependent characteristics for
      a given cut.
    */
    class CutAngDep : public s13::ana::SimpleTTreeAlgorithmBase {

    public:

      using Yield = s13::ana::ScatteringYield;
      using DsType = s13::ana::DatasetType;
      using Dsdcs = s13::ana::Dsdcs;
      using ElaCuts = s13::ana::ElasticScatteringCuts;
      using Evt = s13::ana::ScatteringEvent;
      using Espri = s13::ana::Espri;

      static constexpr double k_min_theta = s13::gk_cs_range_start;
      static constexpr double k_max_theta = s13::gk_cs_range_end;

      using CutProc = bool (ElaCuts::*)(const Evt& evt) const;
      using HistInitProc = std::function<TH1*(double, double)>;
      using HistFillProc = std::function<double(const Evt& evt)>;
      using AngularSnData = std::vector< std::vector<TH1*> >;

      std::vector<TH1*> cut_dists;
      AngularSnData sn_in;
      AngularSnData sn_out;
      AngularSnData sn_tot;
      AngularSnData outphi;

    private:

      s13::misc::MessageLogger logger_;

      std::vector<ElaCuts> cuts_range_;

      size_t bin_num_;

      CutProc cut_proc_;
      HistInitProc cut_hist_init_proc_;
      HistFillProc hist_fill_proc_;

    public:

      CutAngDep(const std::vector<ElaCuts>& cuts,
                CutProc cut_proc,
                HistInitProc hist_init_proc,
                HistFillProc hist_fill_proc,
                size_t bin_num = s13::gk_cs_bin_number) :
        logger_{"CutAngDep"}, cuts_range_{cuts}, bin_num_{bin_num},
        cut_proc_{cut_proc}, cut_hist_init_proc_{hist_init_proc},
        hist_fill_proc_{hist_fill_proc} {
          init_hists();
        }

      CutAngDep(const CutAngDep&) = delete;
      CutAngDep(CutAngDep&&) = delete;
      CutAngDep operator=(const CutAngDep&) = delete;
      CutAngDep operator=(CutAngDep&&) = delete;

      virtual void setup() {
      }

      virtual void process_event(Evt& evt) {
        assert(!cuts_range_.empty());
        fill_ang_dep_hists(evt);
      }

      virtual void finalize() {
        auto bg_norm_fact_interp =
          s13::ana::make_bg_norm_factor_interp2(DsType::he4,
                                               cuts_range_[0].dsdc_selector);
        for (size_t i = 0; i < outphi.size(); i++) {
          auto& v = outphi.at(i);
          for (auto& h : v) {
            double r =
              bg_norm_fact_interp.eval(bin_center_theta(i),
                                       cuts_range_[0].GetPhiCutWidth(bin_center_theta(i)));
            h->Scale(r - 1);
          }
        }
      }

      const std::vector<ElaCuts>& cuts_range() {
        return cuts_range_;
      }

      size_t bin_num() {
        return bin_num_;
      }

      double bin_width() {
        return (k_max_theta - k_min_theta) / bin_num_;
      }

      double bin_min_theta(int idx) {
        double minth = k_min_theta + bin_width() * idx;
        return minth;
      }

      double bin_max_theta(int idx) {
        double maxth = k_min_theta + bin_width() * (idx + 1);
        return maxth;
      }

      double bin_center_theta(int idx) {
        return (bin_min_theta(idx) + bin_max_theta(idx)) / 2;
      }

    private:

      TH1* init_kin_corr_hist(TString title) {
        auto th1 = &s13::ana::make_th1;
        auto hist = th1(80, -8, 8, title, "#theta_{exp} - #theta_{theor}", "Counts");
        return hist;
      }

      TH2* init_2d_kin_corr_hist(TString title) {
        auto th2 = &s13::ana::make_th2;
        auto hist = th2(200, 0, 12, 200, 50, 75, title, "#theta_{He}", "#theta_{p}");
        return hist;
      }

      int bin_idx(Evt& evt) {
        if (evt.p_theta_eff - k_min_theta < 0) {
          return -1;
        }

        size_t idx = (evt.p_theta_eff - k_min_theta) / bin_width();
        if (idx < bin_num_) {
          return idx;
        } else {
          return -1;
        }
      }

      void init_hists() {
        for (size_t i = 0; i < bin_num_; i++) {
          auto e_hist = cut_hist_init_proc_(bin_min_theta(i),
                                            bin_max_theta(i));
          cut_dists.push_back(e_hist);

          sn_in.push_back(std::vector<TH1*>());
          sn_out.push_back(std::vector<TH1*>());
          sn_tot.push_back(std::vector<TH1*>());
          outphi.push_back(std::vector<TH1*>());
          for (size_t c = 0; c < cuts_range_.size(); c++) {
            sn_in[i].push_back(init_kin_corr_hist("IN"));
            sn_out[i].push_back(init_kin_corr_hist("OUT"));
            sn_tot[i].push_back(init_kin_corr_hist("TOT"));
            outphi[i].push_back(init_kin_corr_hist("OUT_phi"));
          }
        }
      }

      void fill_ang_dep_hists(Evt& evt) {
        Double_t he_theta = cuts_range()[0].dsdc_selector == Dsdcs::s1dc
          ? evt.s1dc_theta : evt.fdc0_theta;
        double p_theta = evt.p_theta_eff;
        double he_theta_diff = he_theta - evt.he_theta_theor;

        if (!is_well_defined(p_theta)
            || !is_well_defined(he_theta) ||
            !evt.IsEspriEvent()) {
          return;
        }

        int evt_bin_idx =  bin_idx(evt);

        if (evt_bin_idx == -1) {
          return;
        }

        double cut_dist = hist_fill_proc_(evt);

        for (size_t c = 0; c < cuts_range_.size(); c++) {
          auto loc_cut = cuts_range_[c];
          if (loc_cut.IsSignalEvent(evt)) {
            cut_dists[evt_bin_idx]->Fill(cut_dist);
            sn_tot[evt_bin_idx][c]->Fill(he_theta_diff);
            if ((loc_cut.*cut_proc_)(evt)) {
              sn_in[evt_bin_idx][c]->Fill(he_theta_diff);
            } else {
              sn_out[evt_bin_idx][c]->Fill(he_theta_diff);
            }
          }

          if (loc_cut.IsBackgroundEvent(evt) && (loc_cut.*cut_proc_)(evt)) {
            outphi[evt_bin_idx][c]->Fill(he_theta_diff);
          }
        }
      }

    };

  }
}
