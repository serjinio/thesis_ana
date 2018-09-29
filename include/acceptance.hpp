/**
 *
 *
 *
 */


#pragma once


#include <utility>
#include <algorithm>
#include <iostream>
#include <vector>

#include <TMath.h>
#include <TVector3.h>

#include "common.hpp"
#include "msgoutput.hpp"
#include "naiarraymodel.hpp"


namespace s13 {
  namespace ana {

    /*
      Abstract class represents geometrical acceptance
      of some detector.
    */
    class GeometricalAcceptance {

    public:

      // threshold in mm for comparisons, like detecting geometry crossings
      static constexpr double k_geom_thresh = 1;

      GeometricalAcceptance() {};

      virtual bool in_lab_acceptance(TVector3 vertex_pos,
                                     TVector3 direction) = 0;

    };

    /*
      Acceptance of target window.
    */
    class TargetWndAcceptance : public GeometricalAcceptance {

    public:

      // radius of the window
      static constexpr double k_window_radius = 112;
      // minimal & maximal Z which can get into window acceptance;
      static constexpr double k_min_z = 17.8;
      static constexpr double k_max_z = 87.8;
      static constexpr double k_max_abs_y = 33;


      TargetWndAcceptance(double dispZ) :
        logger_{"TargetWndAcceptance"},
        lcs_transform_vec_{0, 0, -dispZ}
      {}

      virtual bool in_lab_acceptance(TVector3 vertex_pos,
                                     TVector3 direction) {
        TVector3 particle_vertex = vertex_pos;
        TVector3 lcs_radius_vec = to_lcs(particle_vertex);

        // Step of intersection search - half of geometrical threshold, should
        // guarantee that we won't skip over.
        double pos_increment = k_geom_thresh / 2;
        // To make sure magnitude of direction vector does not affect the computation
        direction.SetMag(1);
        direction *= pos_increment;

        while (lcs_radius_vec.Mag() < k_window_radius * 1.2) {
          if (is_in_lab_acceptance_lcs(lcs_radius_vec)) {
            // logger_.debug("The following vertex_pos (X, Y, Z): (%.1f, %.1f, %.1f)"
            //               " and direction: (lab theta, lab phi):"
            //               " (%.1f, %.1f) is in acceptance.",
            //               vertex_pos.X(), vertex_pos.Y(), vertex_pos.Z(),
            //               direction.Theta() / s13::gk_d2r,
            //               (direction.Phi() / s13::gk_d2r) - 90);
            // logger_.debug("lcs Z: %.1f", lcs_radius_vec.Z());

            return true;
          }
          particle_vertex += direction;
          // logger_.debug("particle_vertex (X, Y, Z): (%.1f, %.1f, %.1f)",
          //               particle_vertex.X(), particle_vertex.Y(),
          //               particle_vertex.Z());
          lcs_radius_vec = to_lcs(particle_vertex);
          // logger_.debug("lcs_radius_vec (X, Y, Z): (%.1f, %.1f, %.1f)",
          //               lcs_radius_vec.X(), lcs_radius_vec.Y(),
          //               lcs_radius_vec.Z());
        }

        return false;
      }

    private:

      /*
        Converts vector from global to local, window coordinates.
      */
      TVector3 to_lcs(TVector3 vec) {
        return vec + lcs_transform_vec_;
      }

      bool is_at_window_plane_lcs(const TVector3& lcs_radius_vec) {
        // logger_.debug("is_at_window_plane_lcs(): lcs_radius_vec.Mag(): %.1f",
        //               lcs_radius_vec.Mag());
        if (std::abs(lcs_radius_vec.Mag() -  k_window_radius)
            < k_geom_thresh) {
          return true;
        }
        return false;
      }

      bool is_in_lab_acceptance_lcs(const TVector3& lcs_radius_vec) {
        if (!is_at_window_plane_lcs(lcs_radius_vec)) {
          return false;
        }

        if (std::abs(lcs_radius_vec.Y()) < k_max_abs_y &&
            lcs_radius_vec.Z() > k_min_z && lcs_radius_vec.Z() < k_max_z) {
          return true;
        }

        return false;
      }

      s13::misc::MessageLogger logger_;

      TVector3 lcs_transform_vec_;
    };

    class RdcAcceptance : public GeometricalAcceptance {

    public:
      RdcAcceptance(double distZ, double rotY) :
        logger_{"RdcAcceptance"}, distZ_{distZ}, rotY_{rotY} {}

      virtual bool in_lab_acceptance(TVector3 vertex_pos,
                                     TVector3 direction) {
        TVector3 particle_vertex = vertex_pos;
        direction.SetMag(1);
        // Advance particle vertex close to RDC plane
        particle_vertex += direction *
          ((distZ_ - k_geom_thresh) - particle_vertex.Mag());
        direction *= k_geom_thresh / 2;
        while (particle_vertex.Mag() < distZ_ * 1.2) {
          if (is_in_rdc_window(particle_vertex)) {
            return true;
          }
          particle_vertex += direction;
        }

        return false;
      }

    private:

      TVector3 to_lcs(TVector3 vec) {
        vec.RotateY(-rotY_);
        vec.SetZ(vec.Z() - distZ_);
        return vec;
      }

      bool is_in_rdc_window(TVector3& vec) {
        auto lcs_vec = to_lcs(vec);
        if (std::abs(lcs_vec.Z()) < k_geom_thresh) {
          if (std::abs(lcs_vec.Y()) < s13::gk_espri_rdc_width / 2
              && std::abs(lcs_vec.X()) < s13::gk_espri_rdc_width / 2) {
            return true;
          }
        }
        return false;
      }

      s13::misc::MessageLogger logger_;

      double distZ_;
      double rotY_;
    };

    class NaiAcceptance : public GeometricalAcceptance {

    public:
      NaiAcceptance(double distZ, double rotY) :
        logger_{"NaiAcceptance"}, distZ_{distZ}, rotY_{rotY} {}

      virtual bool in_lab_acceptance(TVector3 vertex_pos,
                                     TVector3 direction) {
        TVector3 particle_vertex = vertex_pos;
        direction.SetMag(1);
        // Advance particle vertex close to RDC plane
        particle_vertex += direction *
          ((distZ_ - k_geom_thresh) - particle_vertex.Mag());
        direction *= k_geom_thresh / 2;
        while (particle_vertex.Mag() < distZ_ * 1.5) {
          TVector3 lcs_radius_vec = to_lcs(particle_vertex);
          if (nai_box_.find_hit_nai_id(lcs_radius_vec) != -1) {
            return true;
          }
          particle_vertex += direction;
        }

        return false;
      }

    private:

      TVector3 to_lcs(TVector3 vec) {
        vec.RotateY(-rotY_);
        vec.SetZ(vec.Z() - distZ_);
        return vec;
      }

      s13::misc::MessageLogger logger_;

      double distZ_;
      double rotY_;

      sim::NaiArrayBox nai_box_;
    };

    /*
      Holds preferences to compute the acceptance.
    */
    struct LabAcceptancePrefs {
      double min_theta = 0;
      double max_theta = 90;
      double min_phi = -180;
      double max_phi = 180;
      double step_theta = 0.5;
      double step_phi = 0.5;
    };

    /*
      Preferences to compute acceptance in terms of a & b angles.
    */
    struct AbAngAcceptancePrefs {
      double min_aang = -90;
      double max_aang = 90;
      double min_bang = -90;
      double max_bang = 90;
      double step_aang = 0.5;
      double step_bang = 0.5;
    };


    /*
      Preferences to compute acceptance.
    */
    struct AcceptancePrefs {
      double min_x1 = -90;
      double max_x1 = 90;
      double min_x2 = -90;
      double max_x2 = 90;
      double step_x1 = 0.5;
      double step_x2 = 0.5;
    };

    /*
      Simple struct to hold acceptance points.
      Supports comparison and ordering - for storing it in sets.
    */
    struct PolarPoint {

      PolarPoint(double theta_, double phi_) :
        theta{theta_}, phi{phi_} {};

      double theta;
      double phi;
    };

    bool operator==(const PolarPoint& rhs, const PolarPoint& lhs);
    bool operator<(const PolarPoint& rhs, const PolarPoint& lhs);

    /*
      Struct to hold acceptane points.
      Supports comparison and ordering - for storing in sets.
    */
    struct AbAngPoint {

      AbAngPoint(double aang_, double bang_) :
        aang{aang_}, bang{bang_} {};

      double aang;
      double bang;
    };

    bool operator==(const AbAngPoint& rhs, const AbAngPoint& lhs);
    bool operator<(const AbAngPoint& rhs, const AbAngPoint& lhs);

    /*
      Types of angles in terms of which acceptance is given.
    */
    enum class AcceptanceTypes {
      lab_ab, lab_theta_phi, cm_theta_phi_fragment, cm_theta_phi_recoil
    };

    /*
      Holds cooridnates of the solid angle for detector acceptance.
      type member specifies type of angles held by x1 & x2.
    */
    struct AcceptanceElement {

      AcceptanceElement(double x1_, double x2_) :
        x1{x1_}, x2{x2_} {};

      double x1;
      double x2;
    };

    bool operator==(const AcceptanceElement& rhs, const AcceptanceElement& lhs);
    bool operator<(const AcceptanceElement& rhs, const AcceptanceElement& lhs);

    /*
      Type which holds results of acceptance computation.
    */
    using AcceptanceSet = std::set<AcceptanceElement>;

    /*
      Returns default settings for acceptance computation.
    */
    LabAcceptancePrefs default_lab_acceptance_prefs();
    AbAngAcceptancePrefs default_abang_acceptance_prefs();

    /*
      Computes acceptance of the target window by given reaction vertex cooridnates.

      step_theta and step_phi parameters define step size at which acceptance
      is calculated. If you perform logical operations on different
      acceptances, then you should check that they've been computed
      with the same step settings.
    */
    std::set<PolarPoint>
    lab_acceptance(TVector3 vertex_pos,
                   GeometricalAcceptance& detector_acceptance,
                   LabAcceptancePrefs* prefs = nullptr);

    std::set<AbAngPoint>
    lab_to_abang_acceptance(std::set<s13::ana::PolarPoint>& lab_acc);

    /*
      Converts given set of acceptance points from lab to CM angles.
    */
    std::set<PolarPoint>
    acceptance_lab_to_cm(std::set<PolarPoint> acceptance_lab);

    std::set<AbAngPoint>
    abang_acceptance(TVector3 vertex_pos,
                     GeometricalAcceptance& detector_acceptance,
                     AbAngAcceptancePrefs* prefs = nullptr);

    AcceptanceSet
    compute_acceptance(TVector3 vertex_pos,
                       GeometricalAcceptance& detector_acceptance,
                       AcceptanceTypes type,
                       AcceptancePrefs prefs);

    ///
    /// Functions to compute acceptance of specific detectors
    ///

    TH2*
    make_upstream_lababang_acceptance_hist(const AcceptanceSet& acc,
                                           double bin_size,
                                           TString title,
                                           Espri espri_selector = Espri::left);

    TH2*
    make_upstream_lab_acceptance_hist(const AcceptanceSet& acc,
                                      double bin_size,
                                      TString title,
                                      Espri espri_selector = Espri::left);

    TH2*
    make_upstream_cm_acceptance_hist(const AcceptanceSet& acc,
                                     double bin_size,
                                     TString title,
                                     Espri espri_selector = Espri::left);

    TH2*
    make_upstream_acceptance_hist(AcceptanceTypes type,
                                  const AcceptanceSet& acc_set,
                                  double bin_size,
                                  TString title,
                                  Espri espri_selector = Espri::left);

    TH2*
    make_downstream_lababang_acceptance_hist(const AcceptanceSet& acc,
                                           double bin_size,
                                           TString title);

    TH2*
    make_downstream_lab_acceptance_hist(const AcceptanceSet& acc,
                                      double bin_size,
                                      TString title);

    TH2*
    make_downstream_cm_acceptance_hist(const AcceptanceSet& acc,
                                     double bin_size,
                                     TString title);

    TH2*
    make_downstream_acceptance_hist(AcceptanceTypes type,
                                    const AcceptanceSet& acc_set,
                                    double bin_size,
                                    TString title);

    AcceptanceSet
    compute_tgt_wnd_acceptance(AcceptanceTypes type,
                               double vertexY,
                               double step = 0.5);

    TH2*
    compute_tgt_wnd_acceptance_hist(AcceptanceTypes type,
                                    double vertexY,
                                    double step = 0.5);

    AcceptanceSet
    compute_rdc_acceptance(AcceptanceTypes type,
                           double step = 0.5);

    AcceptanceSet
    compute_nai_acceptance(AcceptanceTypes type,
                           double step = 0.5);

    AcceptanceSet
    compute_fdc2_acceptance(AcceptanceTypes type,
                            double step = 0.5);


    TH2*
    compute_rdc_acceptance_hist(AcceptanceTypes type,
                                double step = 0.5);

    TH2*
    compute_nai_acceptance_hist(AcceptanceTypes type,
                                double step = 0.5);

    TH2*
    compute_fdc2_acceptance_hist(AcceptanceTypes type,
                                 double step = 0.5);

  }
}
