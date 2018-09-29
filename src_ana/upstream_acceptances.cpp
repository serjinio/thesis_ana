
#include "TROOT.h"
#include "TThread.h"
#include "THStack.h"

#include "cli.hpp"
#include "consts.hpp"
#include "drawing.hpp"
#include "rootscript.hpp"

#include "acceptance.hpp"


static s13::misc::CliOptions cli_opts;


static double constexpr k_wnd_disp_z = 10;


THStack* make_hist_stack(std::vector<TH2*> hists, TString title = "") {
  THStack *hs = new THStack("hs", title);
  for (auto el : hists) {
    hs->Add(el);
  }
  return hs;
}

s13::ana::AcceptanceSet
make_symmetric(s13::ana::AcceptanceSet acc_set) {
  s13::ana::AcceptanceSet sym_set;
  for (auto& el : acc_set) {
    sym_set.insert(el);
    sym_set.insert(s13::ana::AcceptanceElement(el.x1, -el.x2));
  }
  return sym_set;
}

s13::ana::AcceptanceSet
acceptance_union(const s13::ana::AcceptanceSet& set1, const s13::ana::AcceptanceSet& set2) {
  s13::ana::AcceptanceSet union_acc;
  std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(),
                        std::inserter(union_acc, union_acc.end()));
  return union_acc;
}

void upstream_acceptances() {
  auto type = s13::ana::AcceptanceTypes::cm_theta_phi_recoil;
  double bin_size = 0.5;
  // auto hist_wnd_acc =
  //   s13::ana::compute_tgt_wnd_acceptance(s13::ana::AcceptanceTypes::cm_theta_phi_recoil, 8);
  auto rdc_acc = s13::ana::compute_rdc_acceptance(type, bin_size);
  rdc_acc = make_symmetric(rdc_acc);
  auto nai_acc =
    s13::ana::compute_nai_acceptance(type, bin_size);
  nai_acc = make_symmetric(nai_acc);
  auto fdc2_acc =
    s13::ana::compute_fdc2_acceptance(s13::ana::AcceptanceTypes::cm_theta_phi_fragment,
                                      bin_size);

  auto common_region = acceptance_union(nai_acc, fdc2_acc);

  std::vector<TH2*> hists;
  std::vector<TString> titles = {"Common region", "RDC", "FDC2"};
  std::vector<int> colors = {kRed, kBlue, kBlack};
  std::vector<s13::ana::AcceptanceSet> accs = {common_region, fdc2_acc, rdc_acc};
  for (size_t i = 0; i < accs.size(); i++) {
    auto h = s13::ana::make_upstream_acceptance_hist(type, accs[i],
                                                     bin_size, titles[i],
                                                     s13::ana::Espri::left);

    h->SetLineColor(colors[i]);
    hists.push_back(h);
  }

  auto hist_stack = make_hist_stack(hists, "Acceptances");
  s13::ana::RootScript script("upstream_acceptances");
  gStyle->SetOptStat(0000000);
  script.NewPage(1).cd();
  hist_stack->Draw("nostack box");
  hist_stack->GetXaxis()->SetTitle("x1 [deg.]");
  hist_stack->GetYaxis()->SetTitle("x2 [deg.]");
  gPad->BuildLegend(0.3, 0.75, 0.9, 0.9);
  gPad->Update();
}

int main(int argc, char** argv) {
  s13::misc::MessageLogger log("main()");
  auto styles = gROOT->GetListOfStyles();
  for (int i = 0; i <= styles->LastIndex(); i++) {
    log.debug("%s", styles->At(i)->GetName());
  }
  gROOT->SetStyle("Modern");
  TThread::Initialize();
  cli_opts.require_cuts_conf_opt(false);
  try {
    cli_opts.parse_options(argc, argv);
  } catch (std::exception& ex) {
    log.error("Invalid argument provided: %s", ex.what());
    return -1;
  }

  upstream_acceptances();
}
