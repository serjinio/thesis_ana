

#include <type_traits>

#include "beam.hpp"


std::unique_ptr<ROOT::Math::Interpolator>
s13::ana::load_beam_vs_tgtr(std::string filename) {
  auto csv_cols = io::load_csv(filename, 2, 2);
  auto* interp = new ROOT::Math::Interpolator(csv_cols[0], csv_cols[1]);
  return std::unique_ptr<std::remove_pointer<decltype(interp)>::type>(interp);
}

double
s13::ana::get_num_beam(DatasetType ds_type, PolarizationDirection pol_dir, double tgtR) {
  double num_incident = -1;

  if (ds_type == DatasetType::he4) {
    auto he4_beam_interp = load_beam_vs_tgtr("data/he4_beam_vs_tgtR.csv");
    if (pol_dir == PolarizationDirection::up) {
      he4_beam_interp = load_beam_vs_tgtr("data/he4_beam_pol_up_vs_tgtR.csv");
    } else if (pol_dir == PolarizationDirection::down) {
      he4_beam_interp = load_beam_vs_tgtr("data/he4_beam_pol_down_vs_tgtR.csv");
    }
    num_incident = he4_beam_interp->Eval(tgtR);
  } else if (ds_type == DatasetType::he6) {
    auto he6_beam_interp = load_beam_vs_tgtr("data/he6_beam_vs_tgtR.csv");
    if (pol_dir == PolarizationDirection::up) {
      he6_beam_interp = load_beam_vs_tgtr("data/he6_beam_pol_up_vs_tgtR.csv");
    } else if (pol_dir == PolarizationDirection::down) {
      he6_beam_interp = load_beam_vs_tgtr("data/he6_beam_pol_down_vs_tgtR.csv");
    }
    num_incident = he6_beam_interp->Eval(tgtR);
  } else if (ds_type == DatasetType::carbon) {
    auto carbon_beam_interp = load_beam_vs_tgtr("data/he6_carbon_runs_beam_vs_tgtR.csv");
    num_incident = carbon_beam_interp->Eval(tgtR);
  } else {
    throw std::invalid_argument("Unsupported beam type");
  }

  return num_incident;
}
