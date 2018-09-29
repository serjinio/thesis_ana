/**
 *   \file protonede.hpp
 *   \brief Functions to compute simulated proton dE-E curve.
 *
 */

#pragma once


#include <string>

#include "TGraph.h"
#include "Math/Interpolator.h"

#include "common.hpp"


namespace s13 {
  namespace ana {

    std::unique_ptr<ROOT::Math::Interpolator>
    load_he_kin_e(std::string filename);

    std::unique_ptr<ROOT::Math::Interpolator>
    load_he_kin_e2(std::string filename);

    std::unique_ptr<ROOT::Math::Interpolator>
    load_stopping_powers(std::string filename);

    Double_t
    compute_proton_pla_de(Double_t proton_e,
                          std::unique_ptr<ROOT::Math::Interpolator>& stp_powers);

    /*
      Returns two interpolator objects to obtain dE & E as a function of proton angle.
    */
    std::pair<
      std::shared_ptr<ROOT::Math::Interpolator>,
      std::shared_ptr<ROOT::Math::Interpolator> >
    compute_he6_proton_ref_de_e();

    TGraph*
    make_de_e_graph(std::shared_ptr<ROOT::Math::Interpolator>& angle_de,
                    std::shared_ptr<ROOT::Math::Interpolator>& angle_e);

    TGraph*
    make_he6_proton_de_e_graph();

    TGraph*
    make_proton_angle_energy_graph(DatasetType ds_type,
                                   Espri espri);

    TGraph*
    make_he4_proton_angle_energy_graph();

    TGraph*
    make_he6_proton_angle_energy_graph();
  }
}
