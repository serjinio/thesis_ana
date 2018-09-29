

#include <TMath.h>
#include <TVector3.h>

#include "consts.hpp"
#include "msgoutput.hpp"
#include "vertzalg.hpp"


double
s13::ana::compute_vertz(double rdc_xpos, double rdc_aang, Espri espri) {
  assert(espri != Espri::both);
  s13::misc::MessageLogger logger("s13::ana::compute_vertz");
  // define RDC hit point in RHS coordinate system
  double xpos = -rdc_xpos;
  double rdc_rot_angle = s13::gk_rdc_center_angle * s13::gk_d2r;
  if (espri == Espri::right) {
    rdc_rot_angle = -rdc_rot_angle;
  }
  double angle = -TMath::ATan(rdc_aang * 1e-3) + rdc_rot_angle;

  TVector3 hit_pos(xpos, 0, s13::gk_dist_rdc_target);
  hit_pos.RotateY(rdc_rot_angle);
  double vertz = hit_pos.Z() - (hit_pos.X() / TMath::Tan(angle));
  logger.trace("Computed vertZ: %.2f given RDC X: %.2f mm;"
               " RDC aang: %.2f mrad", vertz, rdc_xpos, rdc_aang);
  return vertz;
}
