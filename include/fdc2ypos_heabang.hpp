
#include "consts.hpp"
#include "commonalg.hpp"
#include "drawing.hpp"
#include "scattyield.hpp"
#include "elacuts.hpp"


namespace s13 {
  namespace ana {

    class Fdc2YposHeAbangAlg : public SimpleTTreeAlgorithmBase {

      using ElaCuts = s13::ana::ElasticScatteringCuts;
      using Evt = s13::ana::ScatteringEvent;

      static double constexpr k_fdc2_heabang_rel_npoints = 10;
      static double constexpr k_fdc2_heabang_rel_min_aang = -10;
      static double constexpr k_fdc2_heabang_rel_max_aang = 10;

    public:

      Fdc2YposHeAbangAlg(ElaCuts cuts) :
        logger_{"Fdc2YposHeAbangAlg"}, cuts_{cuts} {
          init_hists();
        }

      // no copy or move
      Fdc2YposHeAbangAlg(const Fdc2YposHeAbangAlg&) = delete;
      Fdc2YposHeAbangAlg(Fdc2YposHeAbangAlg&&) = delete;
      Fdc2YposHeAbangAlg operator=(const Fdc2YposHeAbangAlg&) = delete;
      Fdc2YposHeAbangAlg operator=(const Fdc2YposHeAbangAlg&&) = delete;

      virtual void setup() {
        for (auto el : ypos_hebang_hists_) {
          el->Reset();
        }
      }

      virtual void process_event(Evt& evt) {
        if (cuts_.IsSignalEvent(evt)) {
          ypos_heabang_relation_process_event(evt);
        }
      }

      virtual void finalize() {
      }

      virtual void merge(const TTreeAlgorithmBase& other) {
        auto* other_alg = static_cast<const Fdc2YposHeAbangAlg*>(&other);
        assert(other_alg->ypos_hebang_hists_.size()
               == ypos_hebang_hists_.size());
        for (size_t i = 0; i < ypos_hebang_hists_.size(); i++) {
          ypos_hebang_hists_.at(i)->Add(other_alg->ypos_hebang_hists_.at(i));
        }
      }

      virtual Fdc2YposHeAbangAlg* clone() const {
        Fdc2YposHeAbangAlg* my_clone = new Fdc2YposHeAbangAlg(cuts_);
        my_clone->merge(*this);
        return my_clone;
      }

      virtual void serialize(std::ostream& os) {
        throw std::logic_error("Not implemented");
      }

      virtual void deserialize(std::istream& is) {
        throw std::logic_error("Not implemented");
      }

      std::vector<TH2*>& ypos_heabang_hists() {
        return ypos_hebang_hists_;
      }

      /*
        Returns center angle of ypos_heabang histogram by its index.
      */
      double ypos_heabang_hist_aang(size_t idx) {
        assert(idx >= 0 && idx < ypos_hebang_hists_.size());
        double bin_width = (k_fdc2_heabang_rel_max_aang - k_fdc2_heabang_rel_min_aang)
          / k_fdc2_heabang_rel_npoints;
        return idx * bin_width + k_fdc2_heabang_rel_min_aang + bin_width/2;
      }

      void fit_ypos_heabang_hists(TString fit_fn) {
        for (auto el : ypos_hebang_hists_) {
          el->Fit(fit_fn);
        }
      }

    private:

      void init_hists() {
        double bin_width = (k_fdc2_heabang_rel_max_aang - k_fdc2_heabang_rel_min_aang)
          / k_fdc2_heabang_rel_npoints;
        for (double bin_mid_point = k_fdc2_heabang_rel_min_aang + bin_width/2;
             bin_mid_point < k_fdc2_heabang_rel_max_aang;
             bin_mid_point += bin_width) {
          auto hist_title =
            s13::ana::tstrfmt("S1DC_bang vs. FDC2_y (%.1f<S1DC_aang<%.1f)",
                              bin_mid_point - bin_width/2, bin_mid_point + bin_width/2);
          logger_.debug("init_hists(): Initializing histogram: %s", hist_title.Data());
          auto hist = s13::ana::make_th2(100, -6, 6, 100, -600, 600,
                                         hist_title,
                                         "S1DC_bang [deg]",
                                         "FDC2 Y [mm]");
          ypos_hebang_hists_.push_back(hist);
        }
      }

      int fdc2ypos_heabang_hist_idx(double he_aang_rad) {
        double bin_width = (k_fdc2_heabang_rel_max_aang - k_fdc2_heabang_rel_min_aang)
          / k_fdc2_heabang_rel_npoints;
        // convert mrads into degrees
        double aang_deg = he_aang_rad / gk_d2r;

        if (aang_deg > k_fdc2_heabang_rel_min_aang
            && aang_deg < k_fdc2_heabang_rel_max_aang) { // in aang range
          size_t hist_idx = (aang_deg - k_fdc2_heabang_rel_min_aang) / bin_width;
          assert(hist_idx < ypos_hebang_hists_.size());
          return hist_idx;
        } else { // not in range
          return -1;
        }
      }

      void ypos_heabang_relation_process_event(Evt& evt) {
        // only accept events with good FDC2 data
        if (evt.fdc2_ypos < -9000) {
          return;
        }
        // limit he bang to the range of interest
        double abs_bang_deg = std::abs(evt.s1dc_bang) / gk_d2r;
        if (abs_bang_deg > 6) {
          return;
        }

        int hist_idx = fdc2ypos_heabang_hist_idx(evt.s1dc_aang);
        if (hist_idx != -1) {
          double bang_deg = evt.s1dc_bang / gk_d2r;
          ypos_hebang_hists_.at(hist_idx)->Fill(bang_deg, evt.fdc2_ypos);
        }
      }

      misc::MessageLogger logger_;
      ElaCuts cuts_;

      std::vector<TH2*> ypos_hebang_hists_;
    };


    std::shared_ptr<Fdc2YposHeAbangAlg>
    make_fdc2_heabang_alg(ElasticScatteringCuts& cuts);

  }
}
