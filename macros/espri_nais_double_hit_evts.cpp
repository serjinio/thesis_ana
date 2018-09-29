
#include "TChain.h"
#include "TGraph.h"
#include "TH2.h"
#include "TH1.h"

#include <naiarraymodel.hpp>

#include "init_ds.hpp"
#include "macroincludes.hpp"
#include "drawcuts.h"

#include "msgoutput.hpp"


using namespace s13::ana;
using namespace s13::sim;


TChain *g_local_chain = init_4he_partial_dataset(10);


class NaiDoubleHits {

public:

  struct NaiHitEvt {
    double xpos;
    double ypos;
    double aang;
    double bang;
    double* naie_cal;
  };

  /*
    param: name A name to include in output histograms
  */
  NaiDoubleHits(TString name = "") : name_{name}, thresh_e_{0.3} {};

  NaiDoubleHits(const NaiDoubleHits&) = delete;
  NaiDoubleHits(const NaiDoubleHits&&) = delete;
  NaiDoubleHits& operator=(const NaiDoubleHits&) = delete;
  NaiDoubleHits& operator=(const NaiDoubleHits&&) = delete;


  void compute(TChain& chain, const ElasticScatteringCuts& cuts) {
    cuts_ = cuts;
    assert(cuts_.espri_selector != Espri::both);

    reset();
    init_hists();

    ScatteringTreeWalker tree_walker;
    tree_walker.Walk(chain, *this);

    write_summary();
  }

  void operator() (ScatteringEvent& evt) {
    if (!cuts_.IsElasticEvent(evt)) {
      return;
    }
    total_evts_ += 1;

    auto nai_evt = get_nai_event_data(evt);

    int max_e_nai_id = find_max_e_nai_id(nai_evt);
    double max_e_cal = find_max_e(nai_evt);
    TVector3 nai_hit_pos = find_extrp_nai_hit_pos(nai_evt);
    int extrp_nai_id = nais_array_.find_hit_nai_id(nai_hit_pos);

    std::vector<int> nai_ids;
    for (int i = 0; i < 7; i++) {
      if (nai_evt.naie_cal[i] > thresh_e_) {
        nai_ids.push_back(i);
      }
    }

    if (extrp_nai_id != -1 && max_e_nai_id != -1
        && max_e_cal > thresh_e_) {
      // log_.info("Max-E NaI ID: %d. Extrp. NaI ID: %d. ",
      //           max_e_nai_id, extrp_nai_id);
      int nais_diff = std::abs(extrp_nai_id - max_e_nai_id);
      double multiplicity = 7 - nais_diff;
      hist_extrp_diff_->Fill(nais_diff, 1/multiplicity);
      if (nais_diff > 0) {
        hist_extrp_diff_2_->Fill(nais_diff, 1/multiplicity);
      }
    }

    if (nai_ids.size() > 1) {
      multi_hit_evts_ += 1;
      if (nai_ids.size() == 2) {
        int nais_diff = std::abs(nai_ids[1] - nai_ids[0]);
        double multiplicity = 7 - nais_diff;

        hist_nais_diff_->Fill(nais_diff, 1/multiplicity);
      }

      int nais_diff = std::abs(extrp_nai_id - max_e_nai_id);
      int min_diff = 8;
      std::for_each(nai_ids.begin(), nai_ids.end(), [&min_diff, &extrp_nai_id] (int el) {
          if (min_diff > std::abs(extrp_nai_id - el)) {
            min_diff = std::abs(extrp_nai_id - el);
          }
        });

      if (min_diff == nais_diff) {
        extrp_min_diff_evts_ += 1;
      } else {
        extrp_no_min_diff_evts_ += 1;
      }

      // log_.info("Multi-hit event. Max-E NaI ID: %d. Extrp. NaI ID: %d. "
      //           "Event NaI IDs: ", max_e_nai_id, extrp_nai_id);
      std::cout << "\t";
      for (auto el : nai_ids) {
        std::cout << el << ", ";
      }
      std::cout << std::endl;
    }

   if (extrp_nai_id == -1 && max_e_cal > thresh_e_) {
     int nearest_nai_id = find_nearest_nai_id(nai_evt);
     double nearest_nai_dist = find_nearest_nai_dist(nai_evt);
     log_.info("Failed extrp. event. Nearest NaI dist: %.2f."
               " Nearest NaI ID: %d. max-E NaI ID: %d",
               nearest_nai_dist, nearest_nai_id, max_e_nai_id);
     if (nearest_nai_id != -1) {
       int nearest_nai_diff = std::abs(max_e_nai_id - nearest_nai_id);
       int multiplicity = 7 - nearest_nai_diff;
       hist_nearest_nai_diff_->Fill(nearest_nai_diff);
     }
   }
  }

  void draw_all(RootScript& script) {
    script.NewPage(1,1).cd();
    hist_nais_diff_->Draw();
    script.NewPage(1,1).cd();
    hist_extrp_diff_->Draw();
    script.NewPage(1,1).cd();
    hist_extrp_diff_2_->Draw();
    script.NewPage(1,1).cd();
    hist_nearest_nai_diff_->Draw();
  }

private:

  void reset() {
    total_evts_ = 0;
    multi_hit_evts_ = 0;
    extrp_min_diff_evts_ = 0;
    extrp_no_min_diff_evts_ = 0;
  }

  void init_hists() {
    hist_nais_diff_ =
      new TH1F("naisdiff",
               TString::Format("%s Double hit NaI events. NaI IDs difference", name_.Data()),
               7, -0.5, 6.5);
    hist_nais_diff_->GetXaxis()->SetTitle("ID(NaI_1) - ID(NaI_2)");
    hist_nais_diff_->GetYaxis()->SetTitle("Counts");

    hist_extrp_diff_ =
      new TH1F("extrpdiff",
               TString::Format("%s All NaI events. Extrp. & max-E NaI IDs diff.",
                               name_.Data()),
               7, -0.5, 6.5);
    hist_extrp_diff_->GetXaxis()->SetTitle("abs(ID(Extrp_NaI) - ID(max-E NaI))");
    hist_extrp_diff_->GetYaxis()->SetTitle("Counts");

    hist_extrp_diff_2_ =
      new TH1F("extrpdiff2",
               TString::Format("%s All NaI events. Extrp. & max-E NaI IDs diff. (no zero)",
                               name_.Data()),
               6, 0.5, 6.5);
    hist_extrp_diff_2_->GetXaxis()->SetTitle("abs(ID(Extrp_NaI) - ID(max-E NaI))");
    hist_extrp_diff_2_->GetYaxis()->SetTitle("Counts");

    hist_nearest_nai_diff_ =
      new TH1F("nnaidiff",
               TString::Format("%s: Bad type-3 events. Nearest NaI vs. max-E NaI diff.",
                               name_.Data()),
               7, -0.5, 6.5);
    hist_nearest_nai_diff_->GetXaxis()->SetTitle("abs(ID(Extrp_NaI) - ID(max-E NaI))");
    hist_nearest_nai_diff_->GetYaxis()->SetTitle("Counts");

}

  void write_summary() {
    log_.info("Summary for %s...", name_.Data());
    log_.info("Total events: %d", total_evts_);
    log_.info("Multi hit events: %d", multi_hit_evts_);
    double ratio = (float)multi_hit_evts_ / total_evts_ * 100;
    log_.info("Ratio of multi hit events to total events: %.2f%%", ratio);
    double ratio_min_diff = (float)extrp_min_diff_evts_
      / (extrp_min_diff_evts_ + extrp_no_min_diff_evts_) * 100;
    log_.info("%% of events where max-E NaI ID is closest"
              " to extrp. NaI ID: %.2f%%", ratio_min_diff);
  }

  NaiHitEvt get_nai_event_data(ScatteringEvent& evt) {
    NaiHitEvt nai_evt;
    if (cuts_.espri_selector == Espri::left) {
      nai_evt.xpos = evt.esl_xpos;
      nai_evt.ypos = evt.esl_ypos;
      nai_evt.aang = evt.esl_aang * 1e-3;
      nai_evt.bang = evt.esl_bang * 1e-3;
      // nai_evt.aang = evt.esr_xpos / 995;
      // nai_evt.bang = evt.esr_ypos / 995;
      nai_evt.naie_cal = evt.esl_naie_cal;
    } else {
      nai_evt.xpos = evt.esr_xpos;
      nai_evt.ypos = evt.esr_ypos;
      nai_evt.aang = evt.esr_aang * 1e-3;
      nai_evt.bang = evt.esr_bang * 1e-3;
      // nai_evt.aang = evt.esr_xpos / 995;
      // nai_evt.bang = evt.esr_ypos / 995;
      nai_evt.naie_cal = evt.esr_naie_cal;
    }
    return nai_evt;
  }

  TVector3 find_extrp_nai_hit_pos(NaiHitEvt& nai_evt) {
    TVector3 position(nai_evt.xpos, nai_evt.ypos, 0);
    TVector3 direction(nai_evt.aang, nai_evt.bang, 1);
    TVector3 hit_position = nais_array_.find_hit_pos(position, direction);
    return hit_position;
  }

  int find_max_e(NaiHitEvt& nai_evt) {
    return *std::max_element(nai_evt.naie_cal, nai_evt.naie_cal + 7);
  }

  int find_max_e_nai_id(NaiHitEvt& nai_evt) {
    int max_e_nai_id = -1;
    double max_e_cal = 0;
    for (int i = 0; i < 7; i++) {
      if (nai_evt.naie_cal[i] > max_e_cal) {
        max_e_cal = nai_evt.naie_cal[i];
        max_e_nai_id = i;
      }
    }
    return max_e_nai_id;
  }

  double find_nearest_nai_id(NaiHitEvt& nai_evt) {
    TVector3 position(nai_evt.xpos, nai_evt.ypos, 0);
    TVector3 direction(nai_evt.aang, nai_evt.bang, 1);
    int id = nais_array_.find_nearest_nai_id(position, direction);
    return id;
  }

  double find_nearest_nai_dist(NaiHitEvt& nai_evt) {
    TVector3 position(nai_evt.xpos, nai_evt.ypos, 0);
    TVector3 direction(nai_evt.aang, nai_evt.bang, 1);
    double dist = nais_array_.find_nearest_nai_dist(position, direction);
    return dist;
  }


  TString name_;
  ElasticScatteringCuts cuts_;
  NaiArrayBox nais_array_;

  TH1* hist_nais_diff_;
  TH1* hist_extrp_diff_;
  TH1* hist_extrp_diff_2_;
  TH1* hist_nearest_nai_diff_;

  int total_evts_;
  int multi_hit_evts_;
  int extrp_min_diff_evts_;
  int extrp_no_min_diff_evts_;
  double thresh_e_;

  // utilities
  s13::misc::MessageLogger log_;
};

void espri_nais_double_hit_evts() {
  RootScript script("espri_nais_double_hit_evts", g_local_chain);
  gStyle->SetOptFit();
  ElasticScatteringCuts cuts;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_phi_cut = true;
  cuts.apply_theta_cut = true;
  cuts.apply_espri_aang_cut = true;
  cuts.apply_espri_min_e_de_cut = false;
  cuts.apply_espri_e_cut = false;

  cuts.espri_selector = Espri::left;
  NaiDoubleHits doubleHits("ESL");
  doubleHits.compute(*g_local_chain, cuts);

  doubleHits.draw_all(script);
}

int main() {
  espri_nais_double_hit_evts();
}
