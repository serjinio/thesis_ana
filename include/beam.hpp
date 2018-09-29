/**
 *   \file beam.hpp
 *   \brief Beam-related computations
 *
 *
 */


#include "common.hpp"


namespace s13 {
  namespace ana {

    std::unique_ptr<ROOT::Math::Interpolator>
    load_beam_vs_tgtr(std::string filename);

    double
    get_num_beam(DatasetType ds_type, PolarizationDirection pol_dir, double tgtR);

  }
}
