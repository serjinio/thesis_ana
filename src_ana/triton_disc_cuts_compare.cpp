/**
 *   \file triton_disc_cuts.cpp
 *   \brief Comparison of several options for triton discrimination (cuts).
 *
 */

#include "cli.hpp"
#include "init_ds.hpp"
#include "consts.hpp"
#include "csalg.hpp"
#include "drawing.hpp"
#include "scattyield.hpp"
#include "cuts_conf.hpp"
#include "rootscript.hpp"


static s13::misc::CliOptions cli_opts;


using ElaCuts = s13::ana::ElasticScatteringCuts;
using Evt = s13::ana::ScatteringEvent;

class TritonDiscAlg : public s13::ana::TTreeAlgorithmBase {

public:

  TritonDiscAlg(ElaCuts cuts) :
    logger_{"TritonDiscAlg"}, cuts_{cuts} {
      init_hists();
    }

  virtual void setup() {

  }

  virtual void finalize() {

  }

  virtual void process_event(ScatteringEvent& evt) {

  }

  virtual void merge(const TTreeAlgorithmBase& other) {
    throw std::logic_error("Not implemented");
  }

  virtual TTreeAlgorithmBase* clone() const {
    throw std::logic_error("Not implemented");
  }

  virtual void serialize(std::ostream& os) {
    throw std::logic_error("Not implemented");
  }

  virtual void deserialize(std::istream& is) {
    throw std::logic_error("Not implemented");
  }

private:

  void init_hists() {

  }

  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;

};


void triton_disc_cuts() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  auto alg = check_gaps(cuts);

  s13::ana::RootScript script("hodf_acceptance");
  draw_hodf_acc(script, alg);
}

int main(int argc, char** argv) {
  s13::misc::MessageLogger log("main()");
  TThread::Initialize();
  try {
    cli_opts.parse_options(argc, argv);
  } catch (std::exception& ex) {
    log.error("Invalid argument provided: %s", ex.what());
    return -1;
  }

  triton_disc_cuts();
}
