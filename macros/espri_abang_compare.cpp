
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


class AbangCompare {

public:
  AbangCompare(TString name = "") : name_{name} {}

  AbangCompare(const AbangCompare&) = delete;
  AbangCompare(const AbangCompare&&) = delete;
  AbangCompare& operator=(const AbangCompare&) = delete;
  AbangCompare& operator=(const AbangCompare&&) = delete;

  void compute(TChain& chain, const ElasticScatteringCuts& cuts) {
    cuts_ = cuts;
    assert(cuts_.espri_selector != Espri::both);

    reset();

    ScatteringTreeWalker tree_walker;
    tree_walker.Walk(chain, *this);

    write_summary();
  }

  void operator() (ScatteringEvent& evt) {
    if (!cuts_.IsElasticEvent(evt)) {
      return;
    }
    total_evts_ += 1;

    double hit_xpos = 0, hit_ypos = 0;
    double rdc_aang = 0, rdc_bang = 0;
    if (cuts_.espri_selector == Espri::left) {
      hit_xpos = evt.esl_xpos;
      hit_ypos = evt.esl_ypos;
      rdc_aang = evt.esl_aang;
      rdc_bang = evt.esl_bang;
    } else {
      hit_xpos = evt.esr_xpos;
      hit_ypos = evt.esr_ypos;
      rdc_aang = evt.esr_aang;
      rdc_bang = evt.esr_bang;
    }

    if (std::abs(rdc_aang) > 1570) {
      bad_rdc_aang_evts_ += 1;
    }
    if (std::abs(rdc_bang) > 1570) {
      bad_rdc_bang_evts_ += 1;
    }

    // in mrads
    double vertex_aang = hit_xpos / 995 * 1000;
    double vertex_bang = hit_ypos / 995 * 1000;

    hist_rdc_abang_->Fill(rdc_aang, rdc_bang);
    hist_vertex_abang_->Fill(vertex_aang, vertex_bang);

    hist_rdc_vert_aang_->Fill(vertex_aang, rdc_aang);
    hist_rdc_vert_aang_1d_->Fill((vertex_aang - rdc_aang));
    if (std::abs(vertex_aang - rdc_aang) > 30) {
      hist_rdc_vert_aang_1d_2_->Fill((vertex_aang - rdc_aang));
    }

    hist_rdc_vert_bang_->Fill(vertex_bang, rdc_bang);
    hist_rdc_vert_bang_1d_->Fill((vertex_bang - rdc_bang));
    if (std::abs(vertex_bang - rdc_bang) > 30) {
      hist_rdc_vert_bang_1d_2_->Fill((vertex_bang - rdc_bang));
    }
  }

  void draw_all(RootScript& script) {
    script.NewPage(1, 1).cd();
    hist_rdc_abang_->Draw("colz");
    gPad->SetLogz();
    script.NewPage(1, 1).cd();
    hist_vertex_abang_->Draw("colz");
    gPad->SetLogz();

    script.NewPage(1, 1).cd();
    hist_rdc_vert_aang_->Draw("colz");
    gPad->SetLogz();

    script.NewPage(1, 1).cd();
    hist_rdc_vert_aang_1d_->Draw();
    gPad->SetLogy();

    script.NewPage(1, 1).cd();
    hist_rdc_vert_aang_1d_2_->Draw();
    gPad->SetLogy();

    script.NewPage(1, 1).cd();
    hist_rdc_vert_bang_->Draw("colz");
    gPad->SetLogz();

    script.NewPage(1, 1).cd();
    hist_rdc_vert_bang_1d_->Draw();
    gPad->SetLogy();

    script.NewPage(1, 1).cd();
    hist_rdc_vert_bang_1d_2_->Draw();
    gPad->SetLogy();
  }

private:

  void reset() {
    total_evts_ = 0;
    bad_rdc_aang_evts_ = 0;
    bad_rdc_bang_evts_ = 0;

    hist_rdc_abang_ = new TH2F("rdc_abang",
                           TString::Format("%s - RDC a&b angles", name_.Data()),
                           400, -350, 350, 400, -350, 350);
    hist_rdc_abang_->GetYaxis()->SetTitle("b angle [mrad]");
    hist_rdc_abang_->GetXaxis()->SetTitle("a angle [mrad]");

    hist_vertex_abang_ = new TH2F("vert_abang",
                               TString::Format("%s - Vertex a&b angles", name_.Data()),
                               400, -350, 350, 400, -350, 350);
    hist_vertex_abang_->GetYaxis()->SetTitle("b angle [mrad]");
    hist_vertex_abang_->GetXaxis()->SetTitle("a angle [mrad]");

    hist_rdc_vert_aang_ = new TH2F("rdcv_aang",
                                 TString::Format("%s - Vertex & RDC a angles", name_.Data()),
                                 400, -350, 350, 400, -350, 350);
    hist_rdc_vert_aang_->GetYaxis()->SetTitle("RDC a angle [mrad]");
    hist_rdc_vert_aang_->GetXaxis()->SetTitle("Vertex a angle [mrad]");

    hist_rdc_vert_bang_ = new TH2F("rdcv_bang",
                                   TString::Format("%s - Vertex & RDC b angles", name_.Data()),
                                   400, -350, 350, 400, -350, 350);
    hist_rdc_vert_bang_->GetYaxis()->SetTitle("RDC a angle [mrad]");
    hist_rdc_vert_bang_->GetXaxis()->SetTitle("Vertex a angle [mrad]");

    hist_rdc_vert_aang_1d_ = new TH1F("rdcv_aang_1d",
                                   TString::Format("%s - Vertex - RDC a angles", name_.Data()),
                                   400, -650, 650);
    hist_rdc_vert_aang_1d_->GetYaxis()->SetTitle("Counts");
    hist_rdc_vert_aang_1d_->GetXaxis()->SetTitle("Vert. - RDC [mrad]");

    hist_rdc_vert_aang_1d_2_ = new TH1F("rdcv_aang_1d_2",
                                      TString::Format("%s - Vertex - RDC a angles", name_.Data()),
                                      400, -650, 650);
    hist_rdc_vert_aang_1d_2_->GetYaxis()->SetTitle("Counts");
    hist_rdc_vert_aang_1d_2_->GetXaxis()->SetTitle("Vert. - RDC [mrad]");

     hist_rdc_vert_bang_1d_ = new TH1F("rdcv_bang_1d",
                                   TString::Format("%s - Vertex - RDC b angles", name_.Data()),
                                   400, -650, 650);
    hist_rdc_vert_bang_1d_->GetYaxis()->SetTitle("Counts");
    hist_rdc_vert_bang_1d_->GetXaxis()->SetTitle("Vert. - RDC [mrad]");

    hist_rdc_vert_bang_1d_2_ = new TH1F("rdcv_bang_1d_2",
                                      TString::Format("%s - Vertex - RDC b angles", name_.Data()),
                                      400, -650, 650);
    hist_rdc_vert_bang_1d_2_->GetYaxis()->SetTitle("Counts");
    hist_rdc_vert_bang_1d_2_->GetXaxis()->SetTitle("Vert. - RDC [mrad]");
  }

  void write_summary() {
    log_.info("Summary for %s...", name_.Data());
    log_.info("Total events: %d", total_evts_);
    log_.info("Bad aang events: %d", bad_rdc_aang_evts_);
    log_.info("Bad bang events: %d", bad_rdc_bang_evts_);
  }


  TString name_;
  ElasticScatteringCuts cuts_;

  // stats
  int total_evts_;
  int bad_rdc_aang_evts_;
  int bad_rdc_bang_evts_;

  // hists
  TH2* hist_rdc_abang_;
  TH2* hist_vertex_abang_;
  TH2* hist_rdc_vert_aang_;
  TH2* hist_rdc_vert_bang_;
  TH1* hist_rdc_vert_aang_1d_;
  TH1* hist_rdc_vert_bang_1d_;
  TH1* hist_rdc_vert_aang_1d_2_;
  TH1* hist_rdc_vert_bang_1d_2_;

  // utilities
  s13::misc::MessageLogger log_;
};

void espri_abang_compare() {
  RootScript script("espri_abang_compare", g_local_chain);
  gStyle->SetOptFit();
  ElasticScatteringCuts cuts;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_phi_cut = true;
  cuts.apply_theta_cut = true;
  cuts.apply_espri_aang_cut = false;
  cuts.apply_espri_min_e_de_cut = false;
  cuts.apply_espri_e_cut = false;

  cuts.espri_selector = Espri::left;
  AbangCompare abangs("ESL");
  abangs.compute(*g_local_chain, cuts);
  abangs.draw_all(script);
}

int main() {
  espri_abang_compare();
}
