#include <TChain.h>
#include <TGraph.h>
#include <TMultiGraph.h>

#include "init_ds.hpp"
#include "consts.hpp"
#include "fdc2algs.hpp"
#include "drawing.hpp"
#include "cli.hpp"
#include "cuts_conf.hpp"
#include "rootscript.hpp"


static s13::misc::CliOptions cli_opts;


void
draw_events_dist(s13::ana::RootScript& script,
                 s13::ana::Fdc2EvtsDistAlgPtr fdc2_events_dist_alg) {
  script.NewPage(1).cd();
  fdc2_events_dist_alg->hist_signal_xy_dist()->Draw("colz");
  script.NewPage(1).cd();
  fdc2_events_dist_alg->hist_bg_xy_dist()->Draw("colz");
  script.NewPage(1).cd();
  std::cout << "stuff" << std::endl;
  fdc2_events_dist_alg->hist_fdc2_ypos_s1dc_bang()->Draw("colz");
  script.NewPage(1).cd();
  fdc2_events_dist_alg->hist_fdc2_xpos_s1dc_aang()->Draw("colz");

  script.NewPage(3,2).cd(s13::ana::PadSeq::row);
  for (TH2* el : fdc2_events_dist_alg->histvec_fdc2_ypos_s1dc_bang()) {
    script.cd();
    el->Draw("colz");
  }

  script.NewPage(3,2).cd(s13::ana::PadSeq::row);
  for (TH2* el : fdc2_events_dist_alg->histvec_fdc2_ypos_s1dc_bang_raw()) {
    script.cd();
    el->Draw("colz");
  }

  script.NewPage(1).cd();
  fdc2_events_dist_alg->hist_s1dc_abang()->Draw("colz");
}

s13::ana::Fdc2EvtsDistAlgPtr
compute_events_dist(s13::ana::ElasticScatteringCuts cuts) {
  s13::misc::MessageLogger log("compute_events_dist()");
  log.debug("Computing events distribution...");
  auto fdc2_events_dist_alg = make_fdc2_events_dist_alg(cuts);

  std::shared_ptr<s13::ana::ScatteringTreeWalker>
    tree_walker(new s13::ana::ScatteringTreeWalker());
  tree_walker->add(fdc2_events_dist_alg);

  if (!cli_opts.use_mt()) {
    tree_walker->walk(*init_dataset_from_cli(cli_opts)[0]);
  } else {
    s13::ana::ParallelTreeWalker
      ptree_walker(tree_walker, init_dataset_from_cli(cli_opts));
    ptree_walker.walk();
  }
  std::cout << "Total elastic events: "
            << fdc2_events_dist_alg->elastic_events_num() << std::endl;

  return fdc2_events_dist_alg;
}

s13::ana::ElasticScatteringCuts
init_cuts(std::string cuts_config_fname) {
  std::fstream fs(cuts_config_fname);
  return s13::io::parse_cuts_config(fs);
}

std::string compute_out_filename(std::string basename,
                                 s13::ana::Espri espri_selector) {
  std::string out_filename = basename;
  if (espri_selector == s13::ana::Espri::both) {
    return out_filename;
  }

  if (espri_selector == s13::ana::Espri::left) {
    return out_filename + "_esl";
  } else if (espri_selector == s13::ana::Espri::right) {
    return out_filename + "_esr";
  } else {
    throw std::invalid_argument("Invalid ESPRI selector value! "
                                "Cannot compute filename.");
  }
}

void fdc2_events_dist() {
  s13::ana::ElasticScatteringCuts cuts = init_cuts(cli_opts.cuts_config());
  std::string out_filename =
    compute_out_filename("fdc2_events_dist", cuts.espri_selector);
  auto dist_alg = compute_events_dist(cuts);

  s13::ana::RootScript script(out_filename);
  draw_events_dist(script, dist_alg);
}

int main(int argc, char** argv) {
  s13::misc::MessageLogger logger("main()");
  TThread::Initialize();
  try {
    cli_opts.parse_options(argc, argv);
  } catch (std::exception& ex) {
    logger.error("Invalid argument provided: %s", ex.what());
    return -1;
  }

  fdc2_events_dist();
}
