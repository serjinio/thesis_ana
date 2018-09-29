/**
 *   \file naiarraymodel.hpp
 *   \brief Contains code to model NaI array.
 *
 *  The model is used to determine which crystal to use for E deposit extraction
 *  basing on proton direction vector.
 *
 */

#include <algorithm>
#include <iostream>
#include <utility>
#include <TVector3.h>

#include "consts.hpp"


namespace s13 {
  namespace sim {

   class NaiCrystalBox {

   public:

     static constexpr double k_widthX = 432;
     static constexpr double k_heightY = 46;
     static constexpr double k_depthZ = 46;

     explicit NaiCrystalBox(int id) : NaiCrystalBox(id, 0, 0, 0) {};

     /*
       Constructs NaI crystal box at distZ mm away from
       origin (along Z axis) and then rotated by rotationX
       degrees around X axis.
     */
     NaiCrystalBox(int id, double distY, double distZ, double rotationX) :
       id_{id}, distY_{distY}, distZ_{distZ}, rotationX_{rotationX} {
         std::cout << "[DEBUG:NaiCrystalBox()] Constructing NaI box with ID: "
                   << id_ << " and center coodinates (X, Y, Z): "
                   << "(0, " << distY_ << ", " << distZ_ << ")" << std::endl;
       }

     int id() const {
       return id_;
     }

     bool contains(const TVector3 &pos) const {
       TVector3 point = lcs_transform(pos);
       if (std::abs(point.X()) < k_widthX / 2 &&
           std::abs(point.Y()) < k_heightY /2 &&
           std::abs(point.Z()) < k_depthZ / 2) {
         return true;
       }

       return false;
     }

     /*
       Computes distance of a given point from the center of NaI.
     */
     double distance_from_center(const TVector3& pos) const {
       TVector3 point = lcs_transform(pos);
       if (std::abs(point.X()) > k_widthX / 2 * 1.1) {
         return 9999;
       }

       return std::sqrt(std::pow(point.Z(), 2) + std::pow(point.Y(), 2));
     }

   private:

     TVector3 lcs_transform(const TVector3& pos) const {
       TVector3 point = pos;
       // transform given point into local coordinate frame
       // of a crystal which center is placed at the origin
       // with its axis along X-axis or coordinate frame
       point.SetZ(point.Z() - distZ_);
       point.SetY(point.Y() - distY_);
       point.RotateX(-rotationX_ * gk_d2r);
       return point;
     }

     int id_;
     double distY_; /*!< displacement along Y axis */
     double distZ_; /*!< displacement along Z axis */
     double rotationX_;
   };

    /*
      Class to determine NaI crystal hit by a given proton
      position and direction vectors.
    */
    class NaiArrayBox {

    public:

      static double constexpr k_nai_hit_search_step = 1;

      explicit NaiArrayBox(double distZ = gk_dist_rdc_naiarray_center) :
        distZ_{distZ} {
        // array contains centers of NaI crystals along Y axis
        double naiY_disps[] = {215, 160, 69, -5, -80, -170, -220};
        double nai_tilt_angles[] = {10.5, 7, 3.5, 0, -3.5, -7, 10.5};
        // double naiY_disps[] = {215, 150, 69, -5, -80, -150, -220};
        std::cout << "[DEBUG:NaiArrayBox()] Constructing NaI array..." << std::endl;

        for (int i = 0; i < gk_espri_num_nais; i++) {
          double plane_disp = gk_espri_naiarray_planes_offset * std::pow(-1, i + 1);
          double naiY = distZ + plane_disp;
          auto nai_crystal = NaiCrystalBox(i, naiY_disps[i], naiY,
                                           nai_tilt_angles[i]);
          crystals_.push_back(nai_crystal);
        }
      };

      /*
        Returns point with coordinates at which a hit to NaI was detected.
        The position and direction vectors are assumed to be from RDC data.

        If no hit was detected then TVector3(0,0,0) will be returned.
      */
      TVector3 find_hit_pos(TVector3 position, TVector3 direction) {
        double dist_nai_plane1 = gk_dist_rdc_naiarray_center
          - gk_espri_naiarray_planes_offset;
        // double dist_nai_plane1 = 10; /*!< for debugging will search starting 10 mm away from RDC */
        double dist_nai_plane2 = distZ_ + gk_espri_naiarray_planes_offset;

        for (double dist = dist_nai_plane1 / 2;
             dist < 3 * dist_nai_plane2;
             dist += k_nai_hit_search_step) {
          TVector3 traj_pos = position + (direction * dist);
          int nai_id = find_hit_nai_id(traj_pos);
          if (nai_id != -1) {
            return traj_pos;
          }
        }

        return position + (direction * distZ_);
      }

      /*
        Returns which NaI crystal was hit by a given position
        of a particle.

        Call find_hit_pos() to get a position vector at first,
        then using it cal this method to get NaI ID.

        Returns: ID of NaI crystal [0, 6] or -1 if there is no NaI crystal
                 at this point.
      */
      int find_hit_nai_id(const TVector3& point) const {
        for (auto& el : crystals_) {
          if (el.contains(point)) {
            return el.id();
          }
        }
        return -1;
      }

      /*
        Finds nearest NaI ID by a given position and direction
        vectors (information from RDC).

        Used when find_hit_pos cannot find intersection with any
        of crystals in the array - for approximate search.
      */
      std::pair<int, double>
      find_nearest_nai_id_dist(const TVector3 position,
                               const TVector3 direction) const {
        double dist_nai_plane1 = gk_dist_rdc_naiarray_center
          - gk_espri_naiarray_planes_offset;
        double dist_nai_plane2 = gk_dist_rdc_naiarray_center
          + gk_espri_naiarray_planes_offset;
        dist_nai_plane1 /= direction.Z();
        dist_nai_plane2 /= direction.Z();
        TVector3 plane1_traj_vec(position + (direction * dist_nai_plane1));
        TVector3 plane2_traj_vec(position + (direction * dist_nai_plane2));

        double min_dist = 92;
        double nai_id = -1;
        for (auto& el : crystals_) {
          auto traj_vec = el.id() % 2 == 0 ? plane1_traj_vec : plane2_traj_vec;
          double dist = el.distance_from_center(traj_vec);
          if (min_dist > dist) {
            min_dist = dist;
            nai_id = el.id();
          }
        }

        return std::make_pair(nai_id, min_dist);
      }

      double find_nearest_nai_dist(const TVector3 position,
                                   const TVector3 direction) const {
        std::pair<int, double> retval =
          find_nearest_nai_id_dist(position, direction);
        return std::get<1>(retval);
      }

      int find_nearest_nai_id(const TVector3 position,
                              const TVector3 direction) const {
        std::pair<int, double> retval =
          find_nearest_nai_id_dist(position, direction);
        return std::get<0>(retval);
      }

    private:

      double distZ_;
      std::vector<NaiCrystalBox> crystals_;
    };

  }
}
