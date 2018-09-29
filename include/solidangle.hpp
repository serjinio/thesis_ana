/**
 *   \file solidangle.hpp
 *   \brief Contains code for computation of solid angles.
 *
 */

#pragma once

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMath.h>

#include "commonalg.hpp"
#include "common.hpp"
#include "msgoutput.hpp"
#include "consts.hpp"


namespace s13 {
  namespace ana {

    /*
      Interface for solid angle computation.
    */
    class SolidAngle {

    public:
      SolidAngle(CoordinateFrameType cf_type = CoordinateFrameType::lab) :
        cf_type_{cf_type} {};

      /*
        Computes solid angle for a bin at a given theta angle.
      */
      virtual double at(double theta) const = 0;

      CoordinateFrameType coordinate_frame() const {
        return cf_type_;
      }

    private:
      CoordinateFrameType cf_type_;
    };

    /*
      Abstract base class for solid angle computation for ESPRI
      defines convenience functions.
    */
    class EspriSolidAngle : public SolidAngle {

    public:
      EspriSolidAngle(CoordinateFrameType cf_type, s13::ana::Espri espri_selector,
                      int bin_num = s13::gk_cs_bin_number,
                      double theta_range_start = s13::gk_cs_range_start,
                      double theta_range_end = s13::gk_cs_range_end) :
        SolidAngle{cf_type}, bin_num_{bin_num},
        theta_range_start_{theta_range_start},
        theta_range_end_{theta_range_end}, espri_selector_{espri_selector} {};

      int bin_num() const {
        return bin_num_;
      }

      double min_theta() const {
        return theta_range_start_;
      }

      double max_theta() const {
        return theta_range_end_;
      }

      s13::ana::Espri espri_selector() const {
        return espri_selector_;
      }

      double bin_width() const {
        return std::abs(theta_range_start_ - theta_range_end_) / bin_num_;
      }

      double bin_min_theta(double center_theta) const {
        return center_theta - bin_width() / 2;
      }

      double bin_max_theta(double center_theta) const {
        return center_theta + bin_width() / 2;
      }

    private:
      int bin_num_;
      double theta_range_start_;
      double theta_range_end_;
      s13::ana::Espri espri_selector_;
    };

    /*
      Computes solid angle for RDC surface/plane.
    */
    class RdcSolidAngle : public EspriSolidAngle {

    public:
      RdcSolidAngle(CoordinateFrameType cf_type, s13::ana::Espri espri_selector,
                    double y_width,
                    int bin_num = s13::gk_cs_bin_number,
                    double theta_range_start = s13::gk_cs_range_start,
                    double theta_range_end = s13::gk_cs_range_end) :
        EspriSolidAngle{cf_type, espri_selector, bin_num,
          theta_range_start, theta_range_end},
        y_width_{y_width}, log_{"RdcSolidAngle"} {}

      virtual double at(double theta) const {
        if (coordinate_frame() == CoordinateFrameType::lab) {
          return compute_lab_sa(theta);
        } else {
          log_.error("Solid angle computation for CM frame not implemented yet!");
          throw "not implemented!";
        }
      }

    private:

      /*
        Computes area of RDC segment in mm2
      */
      Double_t compute_rdc_segment_area(Double_t theta) const {
        Double_t ksi_angle = std::abs(gk_rdc_center_angle - theta);
        Double_t bin_min_angle = ksi_angle - bin_width() / 2;
        Double_t bin_max_angle = ksi_angle + bin_width() / 2;
        Double_t rdc_min_x = gk_dist_rdc_target
          * TMath::Tan(bin_min_angle * gk_d2r);
        Double_t rdc_max_x = gk_dist_rdc_target
          * TMath::Tan(bin_max_angle * gk_d2r);
        log_.debug("bin_min_theta: %.1f; bin_max_theta: %.1f",
                   bin_min_theta(theta), bin_max_theta(theta));
        log_.debug("rdc_min_x: %.1f; rdc_max_x: %.1f", rdc_min_x, rdc_max_x);
        return (rdc_max_x - rdc_min_x) * y_width_;
      }

      double compute_lab_sa(double theta) const {
        Double_t ksi_angle = std::abs(gk_rdc_center_angle - theta);
        Double_t dist_to_bin_center = gk_dist_rdc_target / TMath::Cos(ksi_angle * gk_d2r);
        Double_t rdc_segment_area = compute_rdc_segment_area(theta);
        Double_t solid_angle = rdc_segment_area * TMath::Cos(ksi_angle * gk_d2r)
          / std::pow(dist_to_bin_center, 2);
        // account for two ESPRIs left & right
        if (espri_selector() == Espri::both) {
          solid_angle *= 2;
        }
        log_.debug("ksi_angle: %.2f; dist to bin center [mm]: %.2f; "
                   "RDC segment area [mm2]: %.1f; solid angle [sr]: %.8f",
                   ksi_angle, dist_to_bin_center, rdc_segment_area, solid_angle);
        return solid_angle;
      }


      double y_width_;

      s13::misc::MessageLogger log_;
    };

    /*
      Computes solid angle for RDC surface/plane.
    */
    class RdcSolidAngle2 : public EspriSolidAngle {

    public:
      RdcSolidAngle2(CoordinateFrameType cf_type, s13::ana::Espri espri_selector,
                     double y_width,
                     int bin_num = s13::gk_cs_bin_number,
                     double theta_range_start = s13::gk_cs_range_start,
                     double theta_range_end = s13::gk_cs_range_end) :
        EspriSolidAngle{cf_type, espri_selector, bin_num,
          theta_range_start, theta_range_end},
        y_width_{y_width}, log_{"RdcSolidAngle2 (polar coords)"} {}

      virtual double at(double theta) const {
        if (coordinate_frame() == CoordinateFrameType::lab) {
          return compute_lab_polar_sa(theta);
        } else {
          log_.error("Solid angle computation for CM frame not implemented yet!");
          throw "not implemented!";
        }
      }

    private:

      Double_t bin_width_deg() const {
        return std::abs(max_theta() - min_theta()) / bin_num();
      }

      Double_t bin_width_rad() const {
        return bin_width_deg() * D2R;
      }

      Double_t compute_lab_polar_sa(Double_t theta) const {
        TVector3 phi1_vec(0, -y_width_ / 2, gk_dist_rdc_target);
        phi1_vec.RotateY(theta * gk_d2r);
        TVector3 phi2_vec(0, y_width_ / 2, gk_dist_rdc_target);
        phi2_vec.RotateY(theta * gk_d2r);
        // to check curvature of theta=const dPhi segment
        TVector3 phi_curv_check(0, 0, gk_dist_rdc_target);
        phi_curv_check.RotateY(theta * gk_d2r);
        Double_t x0 = phi_curv_check.X();
        Double_t z0 = phi_curv_check.Z();
        // we use phi defined from Y axis, use -90 deg.
        // to convert to our definition
        Double_t delta_phi = std::abs((phi1_vec.Phi() - 90*gk_d2r) -
                                      (phi2_vec.Phi() - 90*gk_d2r));
        phi_curv_check.RotateZ((delta_phi / 2));
        Double_t xlim = phi_curv_check.X();
        Double_t zlim = phi_curv_check.Z();
        // Double_t delta_theta = bin_width_rad();
        Double_t min_theta_deg = theta - bin_width_deg() / 2;
        Double_t max_theta_deg = theta + bin_width_deg() / 2;
        // Double_t solid_angle = TMath::Sin(theta * gk_d2r)
        //   * delta_theta * delta_phi;
        Double_t solid_angle = delta_phi
          * (TMath::Cos(min_theta_deg * gk_d2r)
             - TMath::Cos(max_theta_deg * gk_d2r));

        // account for two ESPRIs left & right
        if (espri_selector() == Espri::both) {
          solid_angle *= 2;
        }

        log_.info("computing solid angle at theta [deg.]: %.2f...", theta);
        log_.info("bin theta width [deg.]: %.2f; phi width [deg]: %.2f; "
                   "solid angle [sr]: %.8f",
                   bin_width_deg(), delta_phi / gk_d2r, solid_angle);
        log_.debug("phi angles [deg.]: (%.2f; %.2f)",
                   phi2_vec.Phi() / gk_d2r - 90, phi1_vec.Phi() / gk_d2r - 90);
        log_.debug("X0;Xlim coordinates of check vector: %.2f; %.2f",
                   x0, xlim);
        log_.debug("X displacement for phi=0 to phi=phi_width/2 rotation: %.2f",
                   x0 - xlim);
        log_.debug("Z displacement for phi=0 to phi=phi_width/2 rotation: %.2f",
                   z0 - zlim);
        log_.debug("Y coordinate of X disp. check vector after phi rotation: %.2f",
                   phi_curv_check.Y());
        return solid_angle;
      }

      double y_width_;

      s13::misc::MessageLogger log_;
    };
  }
}
