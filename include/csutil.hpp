/**
 *   \file csutil.hpp
 *   \brief Utility class - helper computations for CS
 *
 */

#pragma once

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMath.h>

#include "consts.hpp"
#include "common.hpp"
#include "elacuts.hpp"
#include "commonalg.hpp"
#include "solidangle.hpp"
#include "scattyield.hpp"


namespace s13 {
  namespace ana {

    /*
      Makes interpolator for conversion of lab angle to CM.
    */
    TInterpPtr make_lab2cm_angle_interp(DatasetType ds_type);

    // DEPRECATED: those angle conversion functions are deprecated as
    // they use nonrelativistic conversion formulas
    // use make_lab2_cm_angle_interp() instead
    /*
      Convert theta2 angle lab to theta CM.
    */
    // Double_t lab_theta_angle_to_cm(Double_t theta_lab);
    /*
      Converts theta CM to theta2 lab.
    */
    // Double_t cm_theta_angle_to_lab(Double_t theta_cm);

    /*
      Convert solid angle from lab to CM.
    */
    double lab2cm_solid_angle_factor(double proton_cm_theta_rad,
                                     DatasetType ds_type);

    /*
      Converts yield from lab to center of mass frame.
    */
    ScatteringYield lab_yield_to_cm(const ScatteringYield& lab_yield,
                                    DatasetType ds_type);

    /*
      Converts CS from lab to center of mass frame.
     */
    ScatteringYield lab_cs_to_cm(const ScatteringYield& lab_cs,
                                 DatasetType ds_type);

    /*
      Computes cross section values from passed elastic scattering yields
      and experimental parameters. Essentialy it just converts
      counts into mbarns.

      Params:
        solid_angle: SolidAngle interface implementationwhich will return
          solid angle for the bin.
        num_beam_incident: Number of incident beam particles.
        tgt_areal_density: Target density [cm-2]
      Returns: ScatteringYield object containing yields in mbarn units.
    */
    ScatteringYield
    cs_from_yields(ScatteringYield& yields,
                   SolidAngle& solid_angle,
                   double num_beam_incident,
                   double tgt_areal_density);

  }
}
