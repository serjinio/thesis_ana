
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


TChain *g_local_chain = init_4he_partial_dataset(15);


std::ostream& operator<<(ostream& os, const TVector3& vec) {
  os << vec.X() << ", " << vec.Y() << ", " << vec.Z();
  return os;
}

class NaisHitPattern {

public:

  struct NaiHitEvt {
    double xpos;
    double ypos;
    double aang;
    double bang;
    double* naie_cal;
  };

  /*
    Constructs class for hit pattern computation.

    param: name Gives common name to histograms produced
                by this class.
  */
  NaisHitPattern(TString name) :
    name_{name}, thresh_e_{0.3}, log_{"NaisHitPattern"} {
      init_hists();
    }

  NaisHitPattern(const NaisHitPattern&) = delete;
  NaisHitPattern(const NaisHitPattern&&) = delete;
  NaisHitPattern& operator=(const NaisHitPattern&) = delete;
  NaisHitPattern& operator=(const NaisHitPattern&&) = delete;

  void compute(TChain& chain, const ElasticScatteringCuts& cuts) {
    cuts_ = cuts;
    reset();

    ScatteringTreeWalker tree_walker;
    tree_walker.Walk(chain, *this);

    write_summary();
  }

  void reset() {

    total_evts_ = 0;
    good_evts_ = 0;
    bad_type1_evts_ = 0;
    bad_type1_pm1_evts_ = 0;
    bad_type2_evts_ = 0;
    bad_type3_evts_ = 0;
    bad_type4_evts_ = 0;
    extrp_evts_ = 0;
    welldef_e_evts_ = 0;
    bad_ext_evts_ = 0;
    bad_pm1_ext_evts_ = 0;
    no_extrp_naiid_evts_ = 0;
  }

  void write_summary() {
    log_.info("NOTE: all percentages are relative to total events #\n");
    log_.info("Total events (inside of elastic cuts): %.0f", total_evts_);
    log_.info("Well-defined E events: %.0f (%.2f%% from total events)" ,
              welldef_e_evts_, welldef_e_evts_ / total_evts_ * 100);
    log_.info("Good events: %.0f (%.2f%%)",
              good_evts_, good_evts_ / total_evts_ * 100);
    log_.info("'bad' type-1 events: %.0f (%.2f%%)",
              bad_type1_evts_, bad_type1_evts_ / total_evts_ * 100);
    log_.info("\t# events where extrp. NaI ID +/-1"
              " from well-defined E NaI ID: %.0f"
              " (%.2f%% from 'bad' type-1 events); ",
              bad_pm1_ext_evts_, bad_pm1_ext_evts_ / bad_ext_evts_ * 100);

    log_.info("'bad' type-2 events: %.0f (%.2f%%)",
              bad_type2_evts_, bad_type2_evts_ / total_evts_ * 100);
    log_.info("'bad' type-3 events: %.0f (%.2f%%)",
              bad_type3_evts_, bad_type3_evts_ / total_evts_ * 100);
    log_.info("'bad' type-4 events: %.0f (%.2f%%)",
              bad_type4_evts_, bad_type4_evts_ / total_evts_ * 100);

    // log_.info("# events with failed extrp. (no NaI ID): %.0f "
    //           "(%.2f%% from all well-defined E events)",
    //           no_extrp_naiid_evts_, no_extrp_naiid_evts_ / welldef_e_evts_ * 100);
    // log_.info("# events with diff. NaI ID: %.0f "
    //           "(%.2f%% from all well-defined E events)",
    //           bad_ext_evts_, bad_ext_evts_ / welldef_e_evts_ * 100);
    // log_.info("\t# events where extrp. NaI ID +/-1"
    //           " from well-defined E NaI ID: %.0f"
    //           " (%.2f%% from all events with diff NaI ID); ",
    //           bad_pm1_ext_evts_, bad_pm1_ext_evts_ / bad_ext_evts_ * 100);
  }

  void operator() (ScatteringEvent& evt) {
    // we want only good proton events
    if (!cuts_.IsElasticEvent(evt)) {
      return;
    }

    total_evts_ += 1;
    auto nai_evt = get_nai_event_data(evt);

    int max_e_nai_id = find_max_e_nai_id(nai_evt);
    double max_e_cal = find_max_e(nai_evt);
    TVector3 nai_hit_pos = find_extrp_nai_hit_pos(nai_evt);
    int extrp_nai_id = nais_array_.find_hit_nai_id(nai_hit_pos);

    if (extrp_nai_id == -1) {
      extrp_nai_id = find_nearest_nai_id(nai_evt);
    }
    if (extrp_nai_id != -1 && max_e_nai_id != extrp_nai_id) {
      int nais_diff = std::abs(extrp_nai_id - max_e_nai_id);
      if (nais_diff < 2) {
        extrp_nai_id = max_e_nai_id;
      }
    }

    double extrp_e_cal = 0;
    if (extrp_nai_id != -1) {
      extrp_e_cal = nai_evt.naie_cal[extrp_nai_id];
      if (extrp_e_cal > thresh_e_) {
        good_evts_ += 1;
        double nai_hit_aang = compute_nai_hit_aang(nai_hit_pos);
        double nai_hit_bang = compute_nai_hit_bang(nai_hit_pos);
        double rdc_hit_aang = compute_rdc_hit_aang(nai_evt.xpos);
        double rdc_hit_bang = compute_rdc_hit_bang(nai_evt.ypos);

        hist_rdc_dist_->Fill(rdc_hit_aang, rdc_hit_bang);
        hist_rdc_bang_dist_->Fill(rdc_hit_bang);
        hist_nais_dist_->Fill(nai_hit_aang, nai_hit_bang);
        hist_nais_bang_dist_->Fill(nai_hit_bang);
        hist_rdc_nais_aang_->Fill(nai_hit_aang, rdc_hit_aang);
      } else if (extrp_e_cal < thresh_e_ && max_e_cal > thresh_e_) {
        bad_type1_evts_ +=1;
        if (std::abs(extrp_nai_id - max_e_nai_id) == 1) {
          bad_type1_pm1_evts_ += 1;
        }
      } else if (extrp_e_cal < thresh_e_ && max_e_cal < thresh_e_) {
        bad_type2_evts_ += 1;
        hist_type2_theta_->Fill(evt.p_theta_eff);
      }
    } else {
      if (max_e_cal > thresh_e_) {
        bad_type3_evts_ += 1;
      } else {
        bad_type4_evts_ += 1;
      }
    }

    if (max_e_cal > thresh_e_) {
      welldef_e_evts_ += 1;
    }

   if (extrp_nai_id != max_e_nai_id && max_e_cal > thresh_e_) {
      log_.info("extrp NaI ID & E: %i; %.2f; max-E NaI ID & E: %i; %.2f",
                extrp_nai_id, extrp_e_cal, max_e_nai_id, max_e_cal);
      log_.info("RDC hit pos. (X, Y): %.2f, %.2f; "
                "extrp. hit pos. (X, Y, Z): %.2f, %.2f, %.2f; "
                "a & b angles: %.2f, %.2f",
                nai_evt.xpos, nai_evt.ypos,
                nai_hit_pos.X(), nai_hit_pos.Y(), nai_hit_pos.Z(),
                nai_evt.aang, nai_evt.bang);
    }

    if (extrp_nai_id != -1
        && extrp_nai_id != max_e_nai_id
        && max_e_cal > thresh_e_) {
      bad_ext_evts_ += 1;
      if (std::abs(extrp_nai_id - max_e_nai_id) == 1) {
        bad_pm1_ext_evts_ += 1;
      }
    }

    if (extrp_nai_id == -1 && max_e_cal > thresh_e_) {
      no_extrp_naiid_evts_ += 1;
    }

    if (max_e_cal > thresh_e_ && extrp_nai_id != -1) {
      hist_ediff_->Fill(max_e_cal - extrp_e_cal);
    }
  }

  void draw_ediff() {
    hist_ediff_->Draw();
    //gPad->SetLogy();
  }

  void draw_nais_dist() {
    hist_nais_dist_->Draw("colz");
    gPad->SetLogz();
  }

  void draw_nais_bang_dist() {
    hist_nais_bang_dist_->Draw();
  }

  void draw_all(RootScript& script) {
    script.NewPage(1,1).cd();
    hist_rdc_dist_->Draw("colz");
    gPad->SetLogz();
    script.NewPage(1,1).cd();
    draw_nais_dist();

    script.NewPage(1,1).cd();
    hist_rdc_bang_dist_->Draw();
    script.NewPage(1,1).cd();
    draw_nais_bang_dist();

    script.NewPage(1,1).cd();
    hist_rdc_nais_aang_->Draw("colz");
    gPad->SetLogz();

    script.NewPage(1,1).cd();
    draw_ediff();

    script.NewPage(1,1).cd();
    hist_type2_theta_->Draw();
  }

private:

  void init_hists() {
    hist_nais_dist_ =
      new TH2F("nais_ext",
               TString::Format("%s - NaIs hit pattern", name_.Data()),
               400, -250, 250, 200, -250, 250);
    hist_nais_dist_->GetYaxis()->SetTitle("Vertex b angle [mrad]");
    hist_nais_dist_->GetXaxis()->SetTitle("Vertex a angle [mrad]");

    hist_rdc_dist_ =
      new TH2F("rdc_ext",
               TString::Format("%s - RDC hit pattern", name_.Data()),
               400, -250, 250, 200, -250, 250);
    hist_rdc_dist_->GetYaxis()->SetTitle("Vertex b angle [mrad]");
    hist_rdc_dist_->GetXaxis()->SetTitle("Vertex a angle [mrad]");

    hist_ediff_ =
      new TH1F("ediff",
               TString::Format("%s scattree(E) - nai_ext(E)", name_.Data()),
               200, 1, 180);
    hist_ediff_->GetXaxis()->SetTitle("sct(E) - nai_ext(E) [MeV]");
    hist_ediff_->GetYaxis()->SetTitle("Counts");

    hist_nais_bang_dist_ =
      new TH1F("nais_ydist",
               TString::Format("%s NaIs b angle dist.", name_.Data()),
               200, -250, 250);
    hist_nais_bang_dist_->GetYaxis()->SetTitle("Counts");
    hist_nais_bang_dist_->GetXaxis()->SetTitle("Vertex b angle [mrad]");

    hist_rdc_bang_dist_ =
      new TH1F("rdc_ydist",
               TString::Format("%s RDC b angle dist.", name_.Data()),
               200, -250, 250);
    hist_rdc_bang_dist_->GetYaxis()->SetTitle("Counts");
    hist_rdc_bang_dist_->GetXaxis()->SetTitle("Vertex a angle [mrad]");

    hist_rdc_nais_aang_ =
      new TH2F("rdc_nais_bang",
               TString::Format("%s RDC vs. Nai aang", name_.Data()),
               400, -250, 250, 200, -250, 250);
    hist_rdc_nais_aang_->GetYaxis()->SetTitle("RDC vertex b angle [mrad]");
    hist_rdc_nais_aang_->GetXaxis()->SetTitle("NaI vertex b angle [mrad]");

    hist_type2_theta_ =
      new TH1F("rdc_type2_theta",
               TString::Format("%s: type-2 events vs. #theta_{p}", name_.Data()),
               200, 50, 75);
    hist_type2_theta_->GetYaxis()->SetTitle("Counts");
    hist_type2_theta_->GetXaxis()->SetTitle("#theta_{p} [lab. deg.]");
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

  TVector3 find_extrp_nai_hit_pos(NaiHitEvt& nai_evt) {
    TVector3 position(nai_evt.xpos, nai_evt.ypos, 0);
    TVector3 direction(nai_evt.aang, nai_evt.bang, 1);
    TVector3 hit_position = nais_array_.find_hit_pos(position, direction);
    return hit_position;
  }

  double find_nearest_nai_dist(NaiHitEvt& nai_evt) {
    TVector3 position(nai_evt.xpos, nai_evt.ypos, 0);
    TVector3 direction(nai_evt.aang, nai_evt.bang, 1);
    double dist = nais_array_.find_nearest_nai_dist(position, direction);
    return dist;
  }

  double find_nearest_nai_id(NaiHitEvt& nai_evt) {
    TVector3 position(nai_evt.xpos, nai_evt.ypos, 0);
    TVector3 direction(nai_evt.aang, nai_evt.bang, 1);
    int id = nais_array_.find_nearest_nai_id(position, direction);
    return id;
  }

  double compute_nai_hit_aang(const TVector3& nai_hit_pos) {
    return nai_hit_pos.X() / (s13::gk_dist_rdc_target + nai_hit_pos.Z()) * 1000;
  }

  double compute_nai_hit_bang(const TVector3& nai_hit_pos) {
    return nai_hit_pos.Y() / (s13::gk_dist_rdc_target + nai_hit_pos.Z()) * 1000;
  }

  double compute_rdc_hit_aang(double rdc_xpos) {
    return rdc_xpos / s13::gk_dist_rdc_target * 1000;
  }

  double compute_rdc_hit_bang(double rdc_ypos) {
    return rdc_ypos / s13::gk_dist_rdc_target * 1000;
  }

  // misc stuff
  TString name_;
  ElasticScatteringCuts cuts_;
  NaiArrayBox nais_array_;

  // settings
  double thresh_e_;

  // stats
  double total_evts_;
  double good_evts_;
  double bad_type1_evts_;
  double bad_type1_pm1_evts_;
  double bad_type2_evts_;
  double bad_type3_evts_;
  double bad_type4_evts_;
  double extrp_evts_;
  double welldef_e_evts_;
  double bad_ext_evts_;
  double bad_pm1_ext_evts_;
  double no_extrp_naiid_evts_;

  // output histograms
  TH1* hist_ediff_;
  TH2* hist_nais_dist_;
  TH2* hist_rdc_dist_;
  TH1* hist_nais_bang_dist_;
  TH1* hist_rdc_bang_dist_;
  TH1* hist_rdc_nais_aang_;
  TH1* hist_type2_theta_;

  // utilities
  s13::misc::MessageLogger log_;
};


void espri_nais_ext() {
  RootScript script("espri_nais_ext", g_local_chain);
  gStyle->SetOptFit();
  ElasticScatteringCuts cuts;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_phi_cut = true;
  cuts.apply_theta_cut = true;
  cuts.apply_espri_aang_cut = true;
  cuts.apply_espri_min_e_de_cut = false;
  cuts.apply_espri_e_cut = false;

  cuts.vertexXY_radius = 10;
  cuts.theta_width = 3.2;
  cuts.phi_width = 37.8;
  cuts.espri_aang_width = 45;

  cuts.espri_selector = Espri::right;
  NaisHitPattern hit_pattern("ESR");
  hit_pattern.compute(*g_local_chain, cuts);
  hit_pattern.draw_all(script);
}

int main() {
  espri_nais_ext();
}
