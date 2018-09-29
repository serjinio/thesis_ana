


#include "TROOT.h"
#include "TPDF.h"
#include "TPaveStats.h"

#include "common.hpp"
#include "cli.hpp"
#include "init_ds.hpp"
#include "consts.hpp"
#include "csalg.hpp"
#include "drawing.hpp"
#include "scattyield.hpp"
#include "protonede.hpp"
#include "treewalk.hpp"
#include "cuts_conf.hpp"
#include "rootscript.hpp"
#include "vertzalg.hpp"


static s13::misc::CliOptions cli_opts;


using Script = s13::ana::RootScript;
using Yield = s13::ana::ScatteringYield;
using DsType = s13::ana::DatasetType;
using PolDir = s13::ana::PolarizationDirection;
using Dsdcs = s13::ana::Dsdcs;
using ElaCuts = s13::ana::ElasticScatteringCuts;
using Evt = s13::ana::ScatteringEvent;
using FitParams = s13::ana::FitParams;
using Espri = s13::ana::Espri;
using LinearInterp = s13::ana::LinearInterp;


class EspriECal : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  static const std::string k_he6_kin_fname;
  static constexpr double k_nai_cal_center_angle1 = 67.0;
  static constexpr double k_nai_cal_center_angle2 = 69.0;
  static constexpr double k_nai_cal_center_angle3 = 59.5;
  static constexpr double k_nai_cal_angle12_width = 1.0;
  static constexpr double k_nai_cal_angle3_width = 7.0;

  std::vector<TH1*> nai_adc_angle1_hists;
  std::vector<TH1*> nai_adc_angle2_hists;
  std::vector<TH1*> nai_adc_angle3_hists;
  std::vector<TH1*> nai_cal_consts_fit_hists;

private:

  static int hist_counter;

  using TInterp = std::unique_ptr<ROOT::Math::Interpolator>;

  TInterp espri_pla_stp_powers =
    s13::ana::load_stopping_powers("data/p_in_BC400.csv");
  double p_e_calib1 = 0;
  double p_e_calib2 = 0;
  double p_e_calib3 = 0;

  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;

public:

  EspriECal(const ElaCuts& cuts) :
    logger_{"EspriECal"}, cuts_{cuts} {
      assert(cuts_.espri_selector != Espri::both
             && "Can calibrate only L/R ESPRI separately!");
      compute_center_proton_e();
      init_hists();
    }

  EspriECal(const EspriECal& other) :
    logger_{"EspriECal"}, cuts_{other.cuts()} {
      assert(cuts_.espri_selector != Espri::both
             && "Can calibrate only L/R ESPRI separately!");
      clone_nai_hists(other.nai_adc_angle1_hists, nai_adc_angle1_hists);
      clone_nai_hists(other.nai_adc_angle2_hists, nai_adc_angle2_hists);
      clone_nai_hists(other.nai_adc_angle3_hists, nai_adc_angle3_hists);
    }

  EspriECal(EspriECal&&) = delete;
  EspriECal operator=(const EspriECal&) = delete;
  EspriECal operator=(EspriECal&&) = delete;

  virtual void setup() {
  }

  virtual void process_event(Evt& evt) {
    if (cuts_.IsSignalEvent(evt)) {
      fill_nai_adc_hists(evt);
    }
  }

  virtual void finalize() {
  }

  ElaCuts cuts() const {
    return cuts_;
  }

  virtual EspriECal* clone() const {
    auto this_clone = new EspriECal(*this);
    return this_clone;
  }

  virtual void merge(const TTreeAlgorithmBase& other) {
    const EspriECal& other_alg =
      static_cast<const EspriECal&>(other);
    assert(other_alg.nai_adc_angle1_hists.size() == nai_adc_angle1_hists.size());
    for (size_t i = 0; i < nai_adc_angle1_hists.size(); i++) {
      nai_adc_angle1_hists[i]->Add(other_alg.nai_adc_angle1_hists[i]);
      nai_adc_angle2_hists[i]->Add(other_alg.nai_adc_angle2_hists[i]);
      nai_adc_angle3_hists[i]->Add(other_alg.nai_adc_angle3_hists[i]);
    }
  }

  void write_nai_e_cal_consts() {
    fit_nai_adc_hists();

    std::string prefix = cuts_.espri_selector == Espri::left ? "esl" : "esr";
    std::string fname = "data/" + prefix + "_nai_e_gaus_means.csv";
    std::string fname2 = "data/" + prefix + "_nai_calib.csv";
    std::ofstream ofs(fname), ofs2(fname2);
    ofs << "#nai_id,ped,cal_const(" << k_nai_cal_center_angle1 << ")"
      "cal_const(" << k_nai_cal_center_angle2 << ")"
      "cal_const(" << k_nai_cal_center_angle3 << ")" << std::endl;
    ofs2 << "#nai_id,ped,cal_const" << std::endl;

    auto mean_consts1 =
      compute_nai_adc_mean_values(nai_adc_angle1_hists, k_nai_cal_center_angle1,
                                  "gaus_fn1");
    auto mean_consts2 =
      compute_nai_adc_mean_values(nai_adc_angle2_hists, k_nai_cal_center_angle2,
                                  "gaus_fn2");
    auto mean_consts3 =
      compute_nai_adc_mean_values(nai_adc_angle3_hists, k_nai_cal_center_angle3,
                                  "gaus_fn3");
    for (size_t i = 0; i < mean_consts1.size(); ++i) {
      ofs << i << "," << "1," << mean_consts1[i] << "," <<
        mean_consts2[i] << "," << mean_consts3[i] << std::endl;
    }

    std::vector<double> fit_consts =
      fit_nai_e_cal_consts({mean_consts1, mean_consts2, mean_consts3},
                           {p_e_calib1, p_e_calib2, p_e_calib3});
    for (size_t i = 0; i < fit_consts.size(); ++i) {
      ofs2 << i << "," << "1," << fit_consts[i] << std::endl;
    }
  }

private:

  std::vector<double>
  fit_nai_e_cal_consts(std::vector< std::vector<double> > mean_consts,
                       std::vector<double> cal_energies) {
    assert(mean_consts.size() == cal_energies.size() && !mean_consts.empty());
    static int fit_hist_no = 1;
    std::vector<double> fit_consts;
    auto fit_fn = [](Double_t* x, Double_t* param) -> Double_t {
      return param[0] * x[0];
    };

    for (size_t i = 0; i < mean_consts[0].size(); ++i) {
      std::vector<double> cs;
      for (size_t k = 0; k < mean_consts.size(); k++) {
        cs.push_back(mean_consts[k][i]);
      }

      TString hist_name = s13::ana::tstrfmt("nai_consts_fit_hist%i", fit_hist_no++);
      TString hist_title = s13::ana::tstrfmt("NaI#%i consts fit", i);
      TH1F *histo = new TH1F(hist_name, hist_title, 200, 0, 100);
      histo->SetMarkerStyle(21); histo->SetMarkerSize(0.8); histo->SetStats(0);

      for(size_t k=0; k < cal_energies.size();  k++) {
        int bin_no = histo->GetXaxis()->FindBin(cal_energies[k]);
        histo->SetBinContent(bin_no, cs[k]);
      }

      TF1* fit_tfn = new TF1("lin_fit", fit_fn, 0, 100, 1);
      fit_tfn->SetLineWidth(2); fit_tfn->SetLineColor(kMagenta);
      fit_tfn->SetParameters(&cs[0]);
      histo->Fit("lin_fit");
      nai_cal_consts_fit_hists.push_back(histo);
      fit_consts.push_back(1/fit_tfn->GetParameters()[0]);
    }

    return fit_consts;
  }

  void fit_nai_adc_hists() {
    for (size_t i = 0; i < nai_adc_angle1_hists.size(); i++) {
      TF1* fit_fn1, *fit_fn2, *fit_fn3;
      if (cuts_.espri_selector == Espri::left) {
        fit_fn1 = new TF1("gaus_fn1", "gaus", 1000, 1610);
        fit_fn1->SetParameters(3, 1300, 150);
        fit_fn2 = new TF1("gaus_fn2", "gaus", 600, 1320);
        fit_fn2->SetParameters(3, 1060, 170);
        fit_fn3 = new TF1("gaus_fn3", "gaus", 1100, 1900);
        fit_fn3->SetParameters(3, 1400, 170);

        nai_adc_angle1_hists[i]->Fit("gaus_fn1", "", "", 1000, 1610);
        nai_adc_angle2_hists[i]->Fit("gaus_fn2", "", "", 600, 1320);
        nai_adc_angle3_hists[i]->Fit("gaus_fn3", "", "", 1000, 1900);
      } else if (cuts_.espri_selector == Espri::right) {
        fit_fn1 = new TF1("gaus_fn1", "gaus", 850, 1400);
        fit_fn1->SetParameters(3, 1120, 150);
        fit_fn2 = new TF1("gaus_fn2", "gaus", 600, 1200);
        fit_fn2->SetParameters(3, 850, 170);
        fit_fn3 = new TF1("gaus_fn3", "gaus", 750, 1400);
        fit_fn3->SetParameters(3, 1100, 200);

        nai_adc_angle1_hists[i]->Fit("gaus_fn1", "", "", 850, 1400);
        nai_adc_angle2_hists[i]->Fit("gaus_fn2", "", "", 600, 1200);
        nai_adc_angle3_hists[i]->Fit("gaus_fn3", "", "", 750, 1400);
      } else {
        throw std::invalid_argument("ESPRI selector cannot be set to both!");
      }
    }
  }

  void clone_nai_hists(const std::vector<TH1*>& other_hists,
                       std::vector<TH1*>& my_hists) {
    for (size_t i = 0; i < other_hists.size(); i++) {
      TString new_name = TString::Format("naiadc%i", hist_counter++);
      auto hist_clone =
        static_cast<TH1*>(other_hists[i]->Clone(new_name));
      my_hists.push_back(hist_clone);
    }
  }

  void log_proton_e(double angle, double width,
                    double p_e_before_pla, double p_e_cal) {
    logger_.info("Will compute NaIs calibration constants at recoil proton "
                 "angle %.2f lab deg.", angle);
    logger_.info("\t - proton E before dE plastic: %.2f MeV", p_e_before_pla);
    logger_.info("\t - proton E after plastic (calibration E): %.2f MeV", p_e_cal);
    logger_.info("\t - proton angle range: +/-%.2f deg", width/2);
  }

  std::vector<double>
  compute_nai_adc_mean_values(const std::vector<TH1*>& nai_adc_hists,
                              double p_e_calib,
                              TString fit_fn_name) {
    std::vector<double> nai_means;
    for (size_t i = 0; i < nai_adc_hists.size(); ++i) {
      auto fit_fn = nai_adc_hists[i]->GetFunction(fit_fn_name);
      assert(fit_fn != nullptr);
      double* params = fit_fn->GetParameters();
      double gaus_mean = params[1];
      logger_.info("Fit params for '%s' (gaus: const, mean, sigma):"
                   " %.2f; %.2f; %.2f",
                   nai_adc_hists[i]->GetTitle(),
                   params[0], params[1], params[2]);

      // double cal_const = p_e_calib / gaus_mean;
      nai_means.push_back(gaus_mean);
    }
    return nai_means;
  }

  void compute_center_proton_e() {
    const std::string
      k_he6_esl_kin_fname = "data/p_he6_esl_p_e_with_e_losses.csv";
    const std::string
      k_he6_esr_kin_fname = "data/p_he6_esr_p_e_with_e_losses.csv";
    const std::string& kin_fname = cuts_.espri_selector == Espri::left ?
      k_he6_esl_kin_fname : k_he6_esr_kin_fname;
    TInterp e_interp = s13::ana::load_he_kin_e2(kin_fname);
    double p_e_angle1 = e_interp->Eval(k_nai_cal_center_angle1 * s13::gk_d2r);
    double p_e_angle2 = e_interp->Eval(k_nai_cal_center_angle2 * s13::gk_d2r);
    double p_e_angle3 = e_interp->Eval(k_nai_cal_center_angle3 * s13::gk_d2r);
    double pla_de_angle1 =
      s13::ana::compute_proton_pla_de(p_e_angle1, espri_pla_stp_powers);
    double pla_de_angle2 =
      s13::ana::compute_proton_pla_de(p_e_angle2, espri_pla_stp_powers);
    double pla_de_angle3 =
      s13::ana::compute_proton_pla_de(p_e_angle3, espri_pla_stp_powers);
    p_e_calib1 = p_e_angle1 - pla_de_angle1;
    p_e_calib2 = p_e_angle2 - pla_de_angle2;
    p_e_calib3 = p_e_angle3 - pla_de_angle3;
    log_proton_e(k_nai_cal_center_angle1, k_nai_cal_angle12_width,
                 p_e_angle1, p_e_calib1);
    log_proton_e(k_nai_cal_center_angle2, k_nai_cal_angle12_width,
                 p_e_angle2, p_e_calib2);
    log_proton_e(k_nai_cal_center_angle3, k_nai_cal_angle3_width,
                 p_e_angle3, p_e_calib3);
  }

  TH1* init_kin_corr_hist(TString title) {
    auto th1 = &s13::ana::make_th1;
    auto hist = th1(100, -8, 8, title, "#theta_{exp} - #theta_{theor}", "Counts");
    return hist;
  }

  TH2* init_2d_kin_corr_hist(TString title) {
    auto th2 = &s13::ana::make_th2;
    auto hist = th2(200, 0, 12, 200, 50, 75, title, "#theta_{He}", "#theta_{p}");
    return hist;
  }

  TH1* init_nai_adc_hist(TString user_title) {
    auto th1 = &s13::ana::make_th1;
    auto title = s13::ana::tstrfmt("%s", user_title.Data());
    auto hist = th1(70, 0, 2000, title, "ADC value [ch]", "Counts");
    hist_counter++;
    return hist;
  }

  void init_hists() {
    for (int i = 0; i < s13::gk_espri_num_nais; i++) {
      auto title = s13::ana::tstrfmt("NaI#%d (E: %.2f MeV @ %.2f deg)", i,
                                     p_e_calib1, k_nai_cal_center_angle1);
      auto hist = init_nai_adc_hist(title);
      nai_adc_angle1_hists.push_back(hist);

      title = s13::ana::tstrfmt("NaI#%d (E: %.2f MeV) @ %.2f deg", i,
                                p_e_calib2, k_nai_cal_center_angle2);
      hist = init_nai_adc_hist(title);
      nai_adc_angle2_hists.push_back(hist);

      title = s13::ana::tstrfmt("NaI#%d (E: %.2f MeV) @ %.2f deg", i,
                                p_e_calib3, k_nai_cal_center_angle3);
      hist = init_nai_adc_hist(title);
      nai_adc_angle3_hists.push_back(hist);
    }
  }

  virtual void fill_nai_adc_hists(Evt& evt) {
    if (!is_well_defined(evt.p_theta)
        || !is_well_defined(evt.he_theta)) {
      return;
    }
    double xpos = cuts_.IsLeftEspriEvent(evt) ? evt.esl_xpos : evt.esr_xpos;
    if (std::abs(xpos) > 9000) {
      return;
    }

    // basing on which angle to calibrate
    // evt.p_theta - measured RDC angle
    // evt.p_theta_eff - angle from kinematics
    // the choice depends on input energy value for calibration
    const double p_theta = evt.p_theta;
    if (p_theta < 53 || p_theta > 73) {
      return;
    }

    double* naie_raw = cuts_.IsLeftEspriEvent(evt) ? evt.esl_naie_raw : evt.esr_naie_raw;
    double min_theta1 = k_nai_cal_center_angle1 - k_nai_cal_angle12_width / 2;
    double max_theta1 = k_nai_cal_center_angle1 + k_nai_cal_angle12_width / 2;
    double min_theta2 = k_nai_cal_center_angle2 - k_nai_cal_angle12_width / 2;
    double max_theta2 = k_nai_cal_center_angle2 + k_nai_cal_angle12_width / 2;
    double min_theta3 = k_nai_cal_center_angle3 - k_nai_cal_angle3_width / 2;
    double max_theta3 = k_nai_cal_center_angle3 + k_nai_cal_angle3_width / 2;

    if (p_theta > min_theta1 && p_theta < max_theta1) {
      for (int i = 0; i < s13::gk_espri_num_nais; i++) {
        if (naie_raw[i] < 50) continue;
        auto hist = nai_adc_angle1_hists[i];
        hist->Fill(naie_raw[i]);
      }
    } else if (p_theta > min_theta2 && p_theta < max_theta2) {
      for (int i = 0; i < s13::gk_espri_num_nais; i++) {
        if (naie_raw[i] < 50) continue;
        auto hist = nai_adc_angle2_hists[i];
        hist->Fill(naie_raw[i]);
      }
    } else if (p_theta > min_theta3 && p_theta < max_theta3) {
      for (int i = 0; i < s13::gk_espri_num_nais; i++) {
        if (naie_raw[i] < 50) continue;
        auto hist = nai_adc_angle3_hists[i];
        hist->Fill(naie_raw[i]);
      }
    }
  }

};

constexpr double EspriECal::k_nai_cal_center_angle1;
constexpr double EspriECal::k_nai_cal_center_angle2;
constexpr double EspriECal::k_nai_cal_center_angle3;
constexpr double EspriECal::k_nai_cal_angle12_width;
constexpr double EspriECal::k_nai_cal_angle3_width;
int EspriECal::hist_counter = 0;


void draw_cut_marker_lines(double cut_width, TH1* hist) {
  s13::ana::draw_marker_line_ysym(cut_width / 2,
                                  hist->GetMinimum(),
                                  hist->GetMaximum());
}

void
draw_nais_adc_hists(s13::ana::RootScript& script, std::shared_ptr<EspriECal> alg) {
  alg->write_nai_e_cal_consts();

  script.NewPage(4, 2);
  for (size_t i = 0; i < alg->nai_adc_angle1_hists.size(); i++) {
    script.cd();
    alg->nai_adc_angle1_hists[i]->Draw();
    gPad->SetLogy();
  }

  script.NewPage(4, 2);
  for (size_t i = 0; i < alg->nai_adc_angle2_hists.size(); i++) {
    script.cd();
    alg->nai_adc_angle2_hists[i]->Draw();
    gPad->SetLogy();
  }

  script.NewPage(4, 2);
  for (size_t i = 0; i < alg->nai_adc_angle3_hists.size(); i++) {
    script.cd();
    alg->nai_adc_angle3_hists[i]->Draw();
    gPad->SetLogy();
  }

  script.NewPage(4, 2);
  for (size_t i = 0; i < alg->nai_cal_consts_fit_hists.size(); i++) {
    script.cd();
    alg->nai_cal_consts_fit_hists[i]->Draw();
    // gPad->SetLogy();
  }
}

void espri_e_cal() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  cuts.use_espri_aang_corr = true;
  auto alg =
    s13::ana::walk_alg<EspriECal>(cuts,
                                 init_dataset_from_cli(cli_opts),
                                 cli_opts.use_mt());

  auto script_fname =
    s13::ana::tstrfmt("%s_%s", cli_opts.output_file_basename().data(),
                      s13::ana::espri_selector_to_str(cuts.espri_selector).data());
  s13::ana::RootScript script(script_fname);
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit();
  draw_nais_adc_hists(script, alg);
}

int main(int argc, char** argv) {
  s13::misc::MessageLogger log("main()");
  // s13::ana::SetRootStyle(s13::ana::RootStyle::publication);

  TThread::Initialize();
  cli_opts.require_cuts_conf_opt(true);
  try {
    cli_opts.parse_options(argc, argv);
  } catch (std::exception& ex) {
    log.error("Invalid argument provided: %s", ex.what());
    return -1;
  }

  espri_e_cal();
}
