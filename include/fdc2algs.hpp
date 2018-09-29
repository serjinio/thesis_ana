
/**
 *   \file fdc2algs.hpp
 *   \brief Algorithms for FDC2 data
 *
 */


#pragma once

#include <utility>
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


namespace s13 {
  namespace ana {

    /*
      Algorithm to compute events distribution on FDC2.
    */
    class Fdc2EventsDistributionAlg : public SimpleTTreeAlgorithmBase {

    public:

      static constexpr double abs_he_max_angle = 9;
      static constexpr double hist_he_angle_step = 3;
      static constexpr int he_aang_hists_num_regions =
        (abs_he_max_angle * 2) / hist_he_angle_step;

      Fdc2EventsDistributionAlg(const ElasticScatteringCuts& cuts) :
        cuts_{cuts}, histvec_fdc2_ypos_s1dc_bang_{},
        histvec_fdc2_ypos_s1dc_bang_raw_{},
        log_{"Fdc2EventsDistributionAlg"} {
          init_hists();
        }

      Fdc2EventsDistributionAlg(const Fdc2EventsDistributionAlg& other) :
        cuts_{other.cuts_}, elastic_events_num_{other.elastic_events_num_},
        histvec_fdc2_ypos_s1dc_bang_{} {
          hist_signal_xy_dist_ =
            static_cast<TH2*>(other.hist_signal_xy_dist_->Clone());
          hist_bg_xy_dist_ =
            static_cast<TH2*>(other.hist_bg_xy_dist_->Clone());
          hist_fdc2_xpos_s1dc_aang_ =
            static_cast<TH2*>(other.hist_fdc2_xpos_s1dc_aang_->Clone());
          hist_fdc2_ypos_s1dc_bang_ =
            static_cast<TH2*>(other.hist_fdc2_ypos_s1dc_bang_->Clone());
          for (auto el : other.histvec_fdc2_ypos_s1dc_bang_) {
            histvec_fdc2_ypos_s1dc_bang_.push_back(static_cast<TH2*>(el->Clone()));
          }
          for (auto el : other.histvec_fdc2_ypos_s1dc_bang_raw_) {
            histvec_fdc2_ypos_s1dc_bang_raw_.push_back(static_cast<TH2*>(el->Clone()));
          }
          hist_s1dc_abang_ =
            static_cast<TH2*>(other.hist_s1dc_abang_->Clone());
      }

      Fdc2EventsDistributionAlg
      operator=(const Fdc2EventsDistributionAlg& other) {
        Fdc2EventsDistributionAlg tmp(other);
        swap(tmp);
        return *this;
      }

      void swap(Fdc2EventsDistributionAlg& other) {
        std::swap(cuts_, other.cuts_);
        std::swap(elastic_events_num_, other.elastic_events_num_);
        std::swap(hist_signal_xy_dist_, other.hist_signal_xy_dist_);
        std::swap(hist_bg_xy_dist_, other.hist_bg_xy_dist_);
        std::swap(hist_fdc2_xpos_s1dc_aang_, other.hist_fdc2_xpos_s1dc_aang_);
        std::swap(hist_fdc2_ypos_s1dc_bang_, other.hist_fdc2_ypos_s1dc_bang_);
        std::swap(histvec_fdc2_ypos_s1dc_bang_,
                  other.histvec_fdc2_ypos_s1dc_bang_);
        std::swap(histvec_fdc2_ypos_s1dc_bang_raw_,
                  other.histvec_fdc2_ypos_s1dc_bang_raw_);
        std::swap(hist_s1dc_abang_, other.hist_s1dc_abang_);
      }

      virtual Fdc2EventsDistributionAlg* clone() const {
        auto* myclone = new Fdc2EventsDistributionAlg(*this);
        return myclone;
      }

      virtual void merge(const TTreeAlgorithmBase& other) {
        const Fdc2EventsDistributionAlg& other_alg =
          static_cast<const Fdc2EventsDistributionAlg&>(other);
        assert(other_alg.histvec_fdc2_ypos_s1dc_bang_.size()
               == histvec_fdc2_ypos_s1dc_bang_.size());
        log_.debug("merge(): self.elastic_events_num_: %d; "
                   "other.elastic_events_num_: %d",
                   elastic_events_num_, other_alg.elastic_events_num_);
        elastic_events_num_ += other_alg.elastic_events_num_;
        hist_signal_xy_dist_->Add(other_alg.hist_signal_xy_dist_);
        hist_bg_xy_dist_->Add(other_alg.hist_bg_xy_dist_);
        hist_fdc2_ypos_s1dc_bang_->Add(other_alg.hist_fdc2_ypos_s1dc_bang_);
        hist_fdc2_xpos_s1dc_aang_->Add(other_alg.hist_fdc2_xpos_s1dc_aang_);
        for (size_t i = 0; i < histvec_fdc2_ypos_s1dc_bang_.size(); i++) {
          histvec_fdc2_ypos_s1dc_bang_[i]->
            Add(other_alg.histvec_fdc2_ypos_s1dc_bang_[i]);
        }
        for (size_t i = 0; i < histvec_fdc2_ypos_s1dc_bang_raw_.size(); i++) {
          histvec_fdc2_ypos_s1dc_bang_raw_[i]->
            Add(other_alg.histvec_fdc2_ypos_s1dc_bang_raw_[i]);
        }
        hist_s1dc_abang_->Add(other_alg.hist_s1dc_abang_);
      }

      void init_hists() {
        hist_signal_xy_dist_ =
          s13::ana::make_th2(200, -2500, 2500, 200, -500, 500,
                             "Signals XY dist", "X [mm]", "Y [mm]");
        hist_bg_xy_dist_ =
          s13::ana::make_th2(200, -2500, 2500, 200, -500, 500,
              "BG XY dist", "X [mm]", "Y [mm]");
        hist_fdc2_xpos_s1dc_aang_ =
          s13::ana::make_th2(200, -2500, 2500, 80, -10, 10,
              "FDC2 X vs. S1DC aang", "X [mm]", "S1DC aang [deg.]");
        hist_fdc2_ypos_s1dc_bang_ =
          s13::ana::make_th2(200, -500, 500, 80, -4, 4,
              "FDC2 Y vs. S1DC bang", "Y [mm]", "S1DC bang [deg.]");
        for (int i = 0; i < 6; i++) {
          double max_angle = abs_he_max_angle - (i * hist_he_angle_step);
          double min_angle = abs_he_max_angle - ((i+1) * hist_he_angle_step);
          TH2* hist_clone = static_cast<TH2*>(hist_fdc2_ypos_s1dc_bang_->Clone());
          TString title = s13::ana::tstrfmt("FDC2 Y vs. S1DC bang (%.1f<He(aang)<%.1f))",
                                            min_angle, max_angle);
          hist_clone->SetTitle(title);
          TH2* hist_clone2 = static_cast<TH2*>(hist_clone->Clone());
          log_.debug("init of hist: %s", title.Data());
          histvec_fdc2_ypos_s1dc_bang_.push_back(hist_clone);
          histvec_fdc2_ypos_s1dc_bang_raw_.push_back(hist_clone2);
        }
        hist_s1dc_abang_ =
          s13::ana::make_th2(200, -100, 100, 200, -4, 4,
                             "S1DC Y vs. S1DC bang", "Y [mm]", "bang [deg.]");
      }

      virtual void setup() {
        reset();
      }

      virtual void process_event(ScatteringEvent& evt) {
        if (cuts_.IsBackgroundEvent(evt)) {
          hist_bg_xy_dist_->Fill(evt.fdc2_xpos, evt.fdc2_ypos);
        }

        if (cuts_.IsReactionTriggerCut(evt) && cuts_.IsVertexXYCut(evt)) {
          size_t hist_idx = (abs_he_max_angle - (evt.s1dc_aang * 57.2))
            / hist_he_angle_step;
          if (hist_idx < histvec_fdc2_ypos_s1dc_bang_.size() && hist_idx >= 0) {
            histvec_fdc2_ypos_s1dc_bang_raw_[hist_idx]->Fill(evt.fdc2_ypos,
                                                             evt.s1dc_bang * 57.2);
          }
        }

        if (cuts_.IsElasticEvent(evt)) {
          elastic_events_num_ += 1;
          hist_signal_xy_dist_->Fill(evt.fdc2_xpos, evt.fdc2_ypos);
          hist_fdc2_xpos_s1dc_aang_->Fill(evt.fdc2_xpos, evt.s1dc_aang * 57.2);
          hist_fdc2_ypos_s1dc_bang_ ->Fill(evt.fdc2_ypos, evt.s1dc_bang * 57.2);
          size_t hist_idx = (abs_he_max_angle - (evt.s1dc_aang * 57.2))
            / hist_he_angle_step;
          if (hist_idx < histvec_fdc2_ypos_s1dc_bang_.size() && hist_idx >= 0) {
            histvec_fdc2_ypos_s1dc_bang_[hist_idx]->Fill(evt.fdc2_ypos,
                                                         evt.s1dc_bang * 57.2);
          }
          hist_s1dc_abang_->Fill(evt.s1dc_ypos, evt.s1dc_bang * 57.2);
        }
      }

      virtual void finalize() {
        log_.debug("Total amount of elastic events: %d", elastic_events_num_);
      }

      TH2* hist_signal_xy_dist() {
        return hist_signal_xy_dist_;
      }

      TH2* hist_bg_xy_dist() {
        return hist_bg_xy_dist_;
      }

      TH2* hist_fdc2_xpos_s1dc_aang() {
        return hist_fdc2_xpos_s1dc_aang_;
      }

      TH2* hist_fdc2_ypos_s1dc_bang() {
        return hist_fdc2_ypos_s1dc_bang_;
      }

      std::vector<TH2*>& histvec_fdc2_ypos_s1dc_bang() {
        return histvec_fdc2_ypos_s1dc_bang_;
      }

      std::vector<TH2*>& histvec_fdc2_ypos_s1dc_bang_raw() {
        return histvec_fdc2_ypos_s1dc_bang_raw_;
      }

      int elastic_events_num() {
        return elastic_events_num_;
      }

      TH2* hist_s1dc_abang() {
        return hist_s1dc_abang_;
      }

      virtual void serialize(std::ostream& is) {
        throw std::logic_error("not implemented");
      }

      virtual void deserialize(std::istream& is) {
        throw std::logic_error("not implemented");
      }

    private:

      void reset() {
        elastic_events_num_ = 0;
        hist_bg_xy_dist_->Reset();
        hist_signal_xy_dist_->Reset();
        hist_fdc2_xpos_s1dc_aang_->Reset();
        hist_fdc2_ypos_s1dc_bang_->Reset();
      }

      ElasticScatteringCuts cuts_;

      int elastic_events_num_ = 0;

      TH2* hist_signal_xy_dist_;
      TH2* hist_bg_xy_dist_;
      TH2* hist_fdc2_xpos_s1dc_aang_;
      TH2* hist_fdc2_ypos_s1dc_bang_;
      TH2* hist_s1dc_abang_;
      std::vector<TH2*> histvec_fdc2_ypos_s1dc_bang_;
      std::vector<TH2*> histvec_fdc2_ypos_s1dc_bang_raw_;

      s13::misc::MessageLogger log_;
    };

    /// type alias
    using Fdc2EvtsDistAlgPtr =
      std::shared_ptr<s13::ana::Fdc2EventsDistributionAlg>;

    /*
      Factory function for the Fdc2EventsDistributionAlg

      Returns: Fdc2EvtsDistAlgPtr - a pointer to constructed algorithm.
    */
    Fdc2EvtsDistAlgPtr
    make_fdc2_events_dist_alg(const ElasticScatteringCuts& cuts) {
      Fdc2EvtsDistAlgPtr palg(new Fdc2EventsDistributionAlg(cuts));
      return palg;
    }

  }
}
