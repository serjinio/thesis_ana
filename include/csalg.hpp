/**
 *   \file csalg.hpp
 *   \brief Algorighms for CS computation
 *
 */

#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TMath.h>
#include <TVector3.h>

#include "common.hpp"
#include "elacuts.hpp"
#include "commonalg.hpp"
#include "treewalker.hpp"
#include "scattevt.hpp"
#include "scattyield.hpp"
#include "bgsubtract.hpp"
#include "csutil.hpp"


namespace s13 {
  namespace ana {

    /*
      TTree algorithm - to monitor BG/signal bins conditions.
    */
    class SignalBgMonitorAlg : public TTreeAlgorithmBase {

    public:

      static constexpr int k_mon_hists_bin_num = 80;

      static int hist_counter;

      /*
        Constructs algorithm with given bin numbers and angular range
        and cuts.
      */
      SignalBgMonitorAlg(int bin_num,
                         double theta_range_start,
                         double theta_range_end,
                         ElasticScatteringCuts& cuts_) :
        bin_num_{bin_num},
        range_start_{theta_range_start},
        range_end_{theta_range_end},
        yields_range_{bin_num, theta_range_start, theta_range_end},
        cuts_{cuts_},
        log_{"SignalBgMonitorAlg"} {
          init_hists();
          init_mon_yields();
          init_kin_interp();
        }

      /*
        This ctor is used to default-construct the algorithm.
        Instances constructed this way are used to deserialize
        existing results from file.
      */
      explicit SignalBgMonitorAlg() :
        bin_num_{s13::gk_cs_bin_number},
        range_start_{s13::gk_cs_range_start},
        range_end_{s13::gk_cs_range_end},
        yields_range_{bin_num_, range_start_, range_end_},
        cuts_{},
        log_{"SignalBgMonitorAlg"} {
          init_hists();
          init_mon_yields();
          init_kin_interp();
        }

      SignalBgMonitorAlg(const SignalBgMonitorAlg& other) :
        bin_num_{other.bin_num_},
        range_start_{other.range_start_},
        range_end_{other.range_end_},
        yields_range_{other.yields_range_},
        yield_signals_{other.yield_signals_},
        yield_bgs_{other.yield_bgs_},
        cuts_{other.cuts_},
        he4_kin_interp_{other.he4_kin_interp_},
        log_{"SignalBgMonitorAlg"} {
            log_.trace("Copy-constructing algorithm...");
            log_.trace("other alg signals hists num.: %d",
                       other.hist_signals_.size());
            for (size_t i = 0; i < other.hist_signals_.size(); i++) {
              TString new_name = TString::Format("sig%i", hist_counter);
              auto signal_clone =
                static_cast<TH1*>(other.hist_signals_[i]->Clone(new_name.Data()));
              log_.trace("cloned from hist: %s; hist with new name: %s",
                         other.hist_signals_[i]->GetName(), signal_clone->GetName());
              hist_signals_.push_back(signal_clone);

              new_name = TString::Format("bg%i", hist_counter);
              auto bg_clone =
                static_cast<TH1*>(other.hist_bgs_[i]->Clone(new_name.Data()));
              log_.trace("cloned from hist: %s; hist with new name: %s",
                         other.hist_bgs_[i]->GetName(), bg_clone->GetName());
              hist_bgs_.push_back(bg_clone);

              ++hist_counter;
            }
          }

      virtual void setup() {
        reset();
      }

      virtual void process_event(ScatteringEvent& evt) {
        if (!yields_range_.IsInRange(evt.p_theta_eff)) {
          return;
        }

        Double_t he_theta = cuts_.dsdc_selector == Dsdcs::s1dc
          ? evt.s1dc_theta : evt.fdc0_theta;
        Double_t he_theta_theor = cuts_.apply_hodf_he4_cut ?
          he4_kin_interp_->Eval(evt.p_theta_eff) : evt.he_theta_theor;
        Double_t he_theta_diff = he_theta - he_theta_theor;

        int bin_no = yields_range_.GetBinNo(evt.p_theta_eff);
        if (cuts_.IsBackgroundEvent(evt)) {
          hist_bgs_[bin_no]->Fill(he_theta_diff);
          yield_bgs_[bin_no].AddEvent(he_theta_diff);
        }

        if (cuts_.IsElasticEvent(evt)) {
          yields_range_.AddEvent(evt.p_theta_eff);
          hist_signals_[bin_no]->Fill(he_theta_diff);
          yield_signals_[bin_no].AddEvent(he_theta_diff);
        }
      }

      virtual void finalize() {}

      const std::vector<TH1*>& signal_hists() {
        return hist_signals_;
      }

      const std::vector<TH1*>& bg_hists() {
        return hist_bgs_;
      }

      const std::vector<ScatteringYield>& signal_yields() {
        return yield_signals_;
      }

      const std::vector<ScatteringYield>& bgs_yields() {
        return yield_bgs_;
      }

      const ScatteringYield& total_yields() {
        return yields_range_;
      }

      virtual SignalBgMonitorAlg* clone() const {
        log_.trace("Cloning self...");
        auto* myclone = new SignalBgMonitorAlg(*this);
        return myclone;
      }

      virtual void merge(const TTreeAlgorithmBase& other) {
        const SignalBgMonitorAlg& other_alg =
          static_cast<const SignalBgMonitorAlg&>(other);
        if (bin_num_ != other_alg.bin_num_
            || range_start_ != other_alg.range_start_
            || range_end_ != other_alg.range_end_) {
          throw std::logic_error("Other_Alg algorithm settings are incompatible "
                                 "with this one to perform a merge!");
        }

        log_.trace("merge(): Merging histograms...");
        for (size_t i = 0; i < hist_signals_.size(); i++) {
          hist_signals_[i]->Add(other_alg.hist_signals_[i]);
          log_.trace("merge(): Merged hist.: %s", hist_signals_[i]->GetName());
          hist_bgs_[i]->Add(other_alg.hist_bgs_[i]);

          yield_signals_[i] += other_alg.yield_signals_[i];
          yield_bgs_[i] += other_alg.yield_bgs_[i];
        }

        yields_range_ += other_alg.yields_range_;
      }

      virtual void serialize(std::ostream& os) {
        yields_range_.Serialize(os);
        for (int i = 0; i < bin_num_; i++) {
          yield_signals_.at(i).Serialize(os);
          yield_bgs_.at(i).Serialize(os);
        }
        return;
      }

      virtual void serialize(std::string prefix, std::string dir) {
        std::ofstream os(s13::ana::tstrfmt("%s/%s_%s_snmon.csv",
                                           dir.data(),
                                           prefix.data(),
                                           cuts_.GetTotalCutName().Data()));
        return serialize(os);
      }

      virtual void deserialize(std::istream& is) {
        yields_range_.Deserialize(is);
        bin_num_ = yields_range_.GetBinsNumber();
        range_start_ = yields_range_.RangeStart();
        range_end_ = yields_range_.RangeEnd();
        reset_mon_yields();
        init_mon_yields();
        for (int i = 0; i < bin_num_; i++) {
          log_.trace("Deserializing yield for theta: %.1f...", yields_range_.GetBinArg(i));
          yield_signals_.at(i).Deserialize(is);
          // log_.trace("Deserialized signal yield:");
          // std::cout << yield_signals_.at(i);
          yield_bgs_.at(i).Deserialize(is);
          // log_.trace("Deserialized BG yield:");
          // std::cout << yield_bgs_.at(i);
        }
      }

    private:

      void reset() {
        yields_range_.Reset();
        for (size_t i = 0; i < hist_signals_.size(); i++) {
          hist_signals_[i]->Reset();
          hist_bgs_[i]->Reset();
          yield_signals_[i].Reset();
          yield_bgs_[i].Reset();
        }
      }

      void init_hists() {
        for (int i = 0; i < bin_num_; ++i) {
          Double_t bin_mean_theta = yields_range_.GetBinArg(i);
          TString signal_hist_name = TString::Format("sig%i", hist_counter);
          TString signal_hist_title = TString::Format("IN #phi @%.2f lab. deg.",
                                                      bin_mean_theta);
          TString bg_hist_name = TString::Format("bg%i", hist_counter);
          TString bg_hist_title = TString::Format("BG @%.2f lab. deg.",
                                                  bin_mean_theta);

          hist_signals_.push_back(new TH1F(signal_hist_name.Data(),
                                           signal_hist_title.Data(),
                                           k_mon_hists_bin_num, -8, 8));
          hist_bgs_.push_back(new TH1F(bg_hist_name.Data(),
                                       bg_hist_title.Data(),
                                       k_mon_hists_bin_num, -8, 8));
          ++hist_counter;
        }
      }

      void init_kin_interp() {
        he4_kin_interp_ = s13::ana::LoadKinematicsData("data/he4_kin.csv");
      }

      void reset_mon_yields() {
        yields_range_.Reset();
        yield_signals_.clear();
        yield_bgs_.clear();
      }

      void init_mon_yields() {
        yields_range_.Reset();
        for (int i = 0; i < bin_num_; ++i) {
          yield_signals_.push_back(s13::ana::ScatteringYield(k_mon_hists_bin_num, -8, 8));
          yield_bgs_.push_back(s13::ana::ScatteringYield(k_mon_hists_bin_num, -8, 8));
        }
      }

      int bin_num_;
      double range_start_;
      double range_end_;
      ScatteringYield yields_range_;

      std::vector<s13::ana::ScatteringYield> yield_signals_;
      std::vector<s13::ana::ScatteringYield> yield_bgs_;

      ElasticScatteringCuts cuts_;

      std::vector<TH1*> hist_signals_;
      std::vector<TH1*> hist_bgs_;

      s13::ana::TInterpPtr he4_kin_interp_;

      s13::misc::MessageLogger log_;
    };

    /*
      Factory for SignalBgMonitorAlg
    */
    std::shared_ptr<SignalBgMonitorAlg>
    make_signal_bg_monitor_alg(int bin_num,
                               double theta_range_start,
                               double theta_range_end,
                               ElasticScatteringCuts& cuts);

    /*
      Factory for SignalBgMonitorAlg. Returns alg. with results
      deserialized from file.
    */
    std::shared_ptr<SignalBgMonitorAlg>
    make_signal_bg_monitor_alg_from_file(std::istream& is);

    /*
      Obtains distribution of elastic scattering events.
    */
    class ElasticScatteringYieldsAlg : public TTreeAlgorithmBase {

    public:

      ElasticScatteringYieldsAlg(int bin_num,
                                 double range_start, double range_end,
                                 ElasticScatteringCuts cuts) :
        bin_num_{bin_num}, range_start_{range_start}, range_end_{range_end},
        cuts_{cuts},
        total_yields_{bin_num, range_start, range_end},
        bg_yields_{bin_num, range_start, range_end},
        total_events_{0},
        log_{"ElasticScatteringYieldsAlg"} {
        }

      ElasticScatteringYieldsAlg() :
        bin_num_{s13::gk_cs_bin_number},
        range_start_{s13::gk_cs_range_start},
        range_end_{s13::gk_cs_range_end},
        cuts_{},
        total_yields_{bin_num_, range_start_, range_end_},
        bg_yields_{bin_num_, range_start_, range_end_},
        total_events_{0},
        log_{"ElasticScatteringYieldsAlg"} {
        }

      virtual void setup() {
        total_yields_.Reset();
        bg_yields_.Reset();
      }

      virtual void process_event(ScatteringEvent& evt) {
        total_events_ += 1;

        if (cuts_.IsBackgroundEvent(evt)) {
          bg_yields_.AddEvent(evt.p_theta_eff);
        }

        if (cuts_.IsElasticEvent(evt)) {
          total_yields_.AddEvent(evt.p_theta_eff);
        }
      }

      virtual void finalize() {
        // log_.debug("Total events after run: %d", total_events_);
      }

      virtual ElasticScatteringYieldsAlg* clone() const {
        auto* myclone =
          new ElasticScatteringYieldsAlg(*this);
        return myclone;
      }

      virtual void merge(const TTreeAlgorithmBase& other) {
        const ElasticScatteringYieldsAlg& other_alg =
          static_cast<const ElasticScatteringYieldsAlg&>(other);
        if (bin_num_ != other_alg.bin_num_
            || range_start_ != other_alg.range_start_
            || range_end_ != other_alg.range_end_) {
          throw std::logic_error("Other_Alg algorithm settings are incompatible "
                                 "with this one to perform a merge!");
        }

        total_yields_ += other_alg.total_yields();
        bg_yields_ += other_alg.bg_yields();
      }

      virtual void serialize(std::ostream& os) {
        total_yields_.Serialize(os);
        bg_yields_.Serialize(os);
      }

      virtual void serialize(std::string prefix, std::string dir) {
        TString fname = s13::ana::tstrfmt("%s/%s_%s_yield.csv",
                                          dir.data(),
                                          prefix.data(),
                                          cuts_.GetTotalCutName().Data());
        log_.trace("Serializing yields into: %s...", fname.Data());
        std::ofstream os(fname);
        return serialize(os);
      }

      virtual void deserialize(std::istream& is) {
        total_yields_.Deserialize(is);
        bg_yields_.Deserialize(is);
      }

      const ScatteringYield total_yields() const {
        return total_yields_;
      }

      const ScatteringYield bg_yields() const {
        return bg_yields_;
      }

    private:

      int bin_num_;
      double range_start_;
      double range_end_;

      ElasticScatteringCuts cuts_;

      ScatteringYield total_yields_;
      ScatteringYield bg_yields_;

      int total_events_;

      s13::ana::TInterpPtr he4_kin_interp_;

      misc::MessageLogger log_;
    };

    std::shared_ptr<ElasticScatteringYieldsAlg>
    make_elastic_scattering_yield_alg(int bin_num,
                                      double theta_range_start,
                                      double theta_range_end,
                                      ElasticScatteringCuts& cuts);

    std::shared_ptr<ElasticScatteringYieldsAlg>
    make_elastic_scattering_yield_alg_from_file(std::istream& is);
  }
}
