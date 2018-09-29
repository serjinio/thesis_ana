

#pragma once


#include <string>
#include <memory>
#include <array>

#include <TMath.h>

#include "consts.hpp"
#include "common.hpp"
#include "scattevt.hpp"


namespace s13 {
  namespace ana {

    class EDistParamsInterp {

      using TInterp = ROOT::Math::Interpolator;
      using TInterpPtr = std::unique_ptr<TInterp>;

      std::string ds_fname_;
      io::CsvColumnsVector dataset_;
      TInterpPtr p_interp_;

    public:

      EDistParamsInterp(std::string interp_ds_fname) :
        ds_fname_{interp_ds_fname} {
        init();
      }

      double operator()(double p_theta) const {
        return p_interp_->Eval(p_theta);
      }

    private:

      void init() {
        dataset_ = io::load_csv(ds_fname_, 2, 2);
        assert(dataset_.size() == 2 && dataset_[0].size() > 0 &&
               "Failed to load interpolation dataset");
        auto interp = new TInterp(dataset_[0].size(),
                                  ROOT::Math::Interpolation::kLINEAR);
        p_interp_.reset(interp);
        p_interp_->SetData(dataset_[0], dataset_[1]);
      }

    };


    std::shared_ptr<EDistParamsInterp>
    make_he6_e_dist_sigma_interp();

    std::shared_ptr<EDistParamsInterp>
    make_he6_e_dist_mean_interp(s13::ana::Espri espri);

    std::shared_ptr<EDistParamsInterp>
    make_he4_e_dist_sigma_interp();

    std::shared_ptr<EDistParamsInterp>
    make_he4_e_dist_mean_interp(s13::ana::Espri espri);


    class ERelMeanInterp {
    private:
      using Evt = s13::ana::ScatteringEvent;
      using NaisFitConsts = std::vector<std::vector<double> >;

      NaisFitConsts esl_nai_lin_fit_consts_;
      NaisFitConsts esr_nai_lin_fit_consts_;

    public:
      ERelMeanInterp(const NaisFitConsts& esl_nai_mean_lin_fit_consts,
                     const NaisFitConsts& esr_nai_mean_lin_fit_consts) :
        esl_nai_lin_fit_consts_{esl_nai_mean_lin_fit_consts},
        esr_nai_lin_fit_consts_{esr_nai_mean_lin_fit_consts} {
          if (esl_nai_lin_fit_consts_.size() != 3 ||
              esr_nai_lin_fit_consts_.size() != 3) {
            throw std::invalid_argument("Should be exactly 2 fit"
                                        " coeffs. for each NaI!");
          }
          if (esl_nai_lin_fit_consts_[0].size() != s13::gk_espri_num_nais ||
              esr_nai_lin_fit_consts_[0].size() != s13::gk_espri_num_nais) {
            throw std::invalid_argument("Number of coeffs. pairs must correspond"
                                        " to the number of NaIs in ESPRI!");
          }
        }


      double eval(const Evt& evt) const {
        int nai_idx = get_nai_idx(evt);
        if (!evt.IsEspriEvent() || nai_idx == -1) {
          return 0;
        }
        const NaisFitConsts& nais_consts = evt.IsLeftEspriEvent() ?
          esl_nai_lin_fit_consts_ : esr_nai_lin_fit_consts_;
        double c = nais_consts[1][nai_idx];
        double k = nais_consts[2][nai_idx];
        return c + k * evt.p_theta_eff;
      }

      double eval(double p_theta, Espri espri, int nai_idx) {
        if (nai_idx < 0 || nai_idx >= s13::gk_espri_num_nais) {
          throw std::invalid_argument("Invalid NaI index passed!");
        }
        if (espri == Espri::both) {
          throw std::invalid_argument("Invalid ESPRI selector!");
        }

        const NaisFitConsts& nais_consts = espri == Espri::left ?
          esl_nai_lin_fit_consts_ : esr_nai_lin_fit_consts_;
        double c = nais_consts[1][nai_idx];
        double k = nais_consts[2][nai_idx];
        return c + k * p_theta;
      }

    private:
      int get_nai_idx(const Evt& evt) const {
        if (!evt.IsEspriEvent()) {
          return -1;
        }

        int max_nai_idx = -1;
        double max_nai_adc = 0;
        const double* naie_raw = evt.IsLeftEspriEvent() ?
          evt.esl_naie_raw : evt.esr_naie_raw;
        for (int i = 0; i < s13::gk_espri_num_nais; i++) {
          if (naie_raw[i] > max_nai_adc) {
            max_nai_adc = naie_raw[i];
            max_nai_idx = i;
          }
        }
        return max_nai_idx;
      }

    };

    std::shared_ptr<ERelMeanInterp> make_he6_nai_erel_mean_interp();

    std::shared_ptr<ERelMeanInterp> make_he4_nai_erel_mean_interp();

  }
}
