
#include <TMath.h>

#include "bgsubtract.hpp"
#include "elacuts.hpp"
#include "drawing.hpp"


using ElaCuts = s13::ana::ElasticScatteringCuts;


s13::ana::OutPhiBgFactor
s13::ana::make_bg_norm_factor_interp2(DatasetType ds_type, Dsdcs dsdc_selector) {
  if (ds_type != DatasetType::he4) {
    throw std::logic_error("Only He4 supports out-phi BG extraction.");
  }
  if (dsdc_selector != Dsdcs::s1dc) {
    throw std::logic_error("Invalid DSDC chamber selected!");
  }

  // S1DC only for now
  // TODO: compute & implement FDC0 coefficients
  std::array<double, 2> a_consts = {{0.05800131907685552, -4.623628452986631e-4}};
  std::array<double, 2> b_consts = {{0.006294687071806979, 3.130789602420486e-4}};
  double c = -0.7791743451828663;
  OutPhiBgFactor bg_fact{a_consts, b_consts, c};
  return bg_fact;
}

s13::ana::ScatteringYield
s13::ana::convert_to_in_phi_bg(DatasetType ds_type,
                               const ElasticScatteringCuts& cuts,
                               ScatteringYield bg_yields) {
  misc::MessageLogger log{"s13::ana::convert_to_in_phi_bg()"};
  auto bg_factor = s13::ana::make_bg_norm_factor_interp2(ds_type, cuts.dsdc_selector);
  ScatteringYield bg_in_yields(bg_yields.GetBinsNumber(),
                               bg_yields.RangeStart(), bg_yields.RangeEnd());
  for (int i = 0; i < bg_yields.GetBinsNumber(); i++) {
    Double_t p_theta = bg_yields.GetBinArg(i);
    Double_t bin_yield_bg_out = bg_yields.GetBinValue(i);
    Double_t r_fact = bg_factor.eval(p_theta, cuts.GetPhiCutWidth(p_theta));
    Double_t bin_yield_bg_in = bin_yield_bg_out * (r_fact - 1);

    log.debug("\tProton theta: %.2f", p_theta);
    log.debug("\tNormalization factor: %.2f", r_fact);
    log.debug("\tbin_yield_bg_in: %.2f", bin_yield_bg_in);

    bg_in_yields.SetBinValue(p_theta, bin_yield_bg_in);
  }

  return bg_in_yields;
}

double
s13::ana::draw_kin_corr_hist_sn(double angle,
                                TH1* bin_hist,
                                Dsdcs dsdc_selector,
                                double theta_width,
                                UseAbsVal use_abs,
                                Degrees bg_min_theta_angle) {
  s13::misc::MessageLogger logger("draw_bin_sn()");
  double sigma = dsdc_selector == Dsdcs::s1dc ?
    s13::gk_theta_distribution_sigma_s1dc :
    s13::gk_theta_distribution_sigma_fdc0;
  const double bg_min_angle = bg_min_theta_angle.get();
  const double bg_max_angle = -3*sigma;

  double bg_int =
    bin_hist->Integral(bin_hist->GetXaxis()->FindBin(bg_min_angle),
                       bin_hist->GetXaxis()->FindBin(bg_max_angle));
  double signal_int =
    bin_hist->Integral(bin_hist->GetXaxis()->FindBin(-theta_width/2),
                       bin_hist->GetXaxis()->FindBin(theta_width/2));
  if (use_abs.get()) {
    bg_int = std::abs(bg_int);
    signal_int = std::abs(signal_int);
  }

  double sn = signal_int / bg_int;
  s13::ana::draw_text(s13::ana::tstrfmt("S/N: %.2f", sn));
  logger.debug("S/N for bin at %.2f deg.: %.2f (S: %.2f; N: %.2f)",
               angle, sn, signal_int, bg_int);

  double hist_ymin = bin_hist->GetMinimum();
  double hist_ymax = bin_hist->GetMaximum() * 0.7;
  s13::ana::draw_marker_line_ysym(theta_width / 2, hist_ymin, hist_ymax);
  s13::ana::draw_marker_line(bg_min_angle, hist_ymin, bg_min_angle, hist_ymax);
  s13::ana::draw_marker_line(bg_max_angle, hist_ymin, bg_max_angle, hist_ymax);
  return sn;
}
