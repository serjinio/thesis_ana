/**
 *   \file treewalker.hpp
 *   \brief Conatains TreeWalker - to traverse a TTree and related types.
 *
 */

#pragma once

#include <assert.h>

#include <TString.h>
#include <TThread.h>
#include <TChain.h>

#include <algorithm>
#include <iostream>
#include <vector>
#include <future>

#include "consts.hpp"
#include "scattevt.hpp"
#include "commonalg.hpp"
#include "msgoutput.hpp"


namespace s13 {
  namespace ana{

    class AbstractTreeWalker {

    public:

      AbstractTreeWalker() {};

      void add(std::shared_ptr<TTreeAlgorithmBase> alg) {
        algorithms.push_back(alg);
      }

      template <typename... Ts>
      void emplace(Ts&&... args) {
        algorithms.emplace_back(std::forward(args)...);
      }

      /*
        Returns algorithm reference at index i.
      */
      const std::shared_ptr<TTreeAlgorithmBase>& algorithm_at(int i) {
        return algorithms.at(i);
      }

      int algorithms_number() {
        return algorithms.size();
      }

    protected:
      std::vector<std::shared_ptr<TTreeAlgorithmBase> > algorithms;

    };


    /*
      Class to extract events from ROOT TTree and pass to provided algorithms.
     */
    class ScatteringTreeWalker : public AbstractTreeWalker {

    public:

      static constexpr double k_inf = -99999;

      ScatteringTreeWalker() : log_{"ScatteringTreeWalker"} {}

      ScatteringTreeWalker(const ScatteringTreeWalker& other) :
        log_{"ScatteringTreeWalker"} {
        for (auto el : other.algorithms) {
          std::shared_ptr<TTreeAlgorithmBase> alg_copy(el->clone());
          algorithms.push_back(alg_copy);
        }
      }

      ScatteringTreeWalker operator=(const ScatteringTreeWalker& other) {
        ScatteringTreeWalker tmp(other);
        swap(tmp);
        return *this;
      }

      void swap(ScatteringTreeWalker& other) noexcept(true) {
        std::swap(algorithms, other.algorithms);
      }

      /*
         Walks through a given tree, calls ProcessEvent of each algorithm.
      */
      void Walk(TTree& tree) {
        // reset state of all algorithms at first
        for (auto alg : algorithms) {
          alg->setup();
        }

        BindTreeVars(tree);
        CleanTreeVars();

        for (int i_entry = 0; i_entry < tree.GetEntries(); i_entry++) {
          CleanTreeVars();
          tree.GetEntry(i_entry);
          ConvertValues();

          auto evt = MakeScatteringEvent();
          for (auto alg : algorithms) {
            alg->process_event(evt);
          }
        }

        // reset state of all algorithms at first
        for (auto alg : algorithms) {
          alg->finalize();
        }

        tree.ResetBranchAddresses();
      }

      void walk(TTree& tree) {
        Walk(tree);
      }

      void operator() (TTree& tree) {
        walk(tree);
      }

    private:

      ScatteringEvent MakeScatteringEvent() {
        ScatteringEvent evt;

        evt.event_no = m_event_no;

        for (int i = 0; i < 8; i++) {
          evt.triggers[i] = m_triggers[i];
        }

        evt.p_phi = m_p_phi;
        evt.p_theta = m_p_theta;
        evt.p_aang = m_p_aang;
        evt.p_bang = m_p_bang;
        evt.p_theta_eff = m_p_theta_eff;
        evt.p_e_eff = m_p_e_eff;
        evt.p_e_sim = m_p_e_sim;
        evt.p_e = m_p_e;
        evt.p_de = m_p_de;
        evt.p_e_raw = m_p_e_raw;
        evt.p_de_raw = m_p_de_raw;
        evt.he_phi = m_s1dc_phi;
        evt.he_theta = m_s1dc_theta;
        evt.he_theta_theor = m_he_theta_theor;
        evt.vert_xpos = m_tgt_up_vert_xpos;
        evt.vert_ypos = m_tgt_up_vert_ypos;
        evt.tgt_up_aang = m_tgt_up_aang;
        evt.tgt_up_bang = m_tgt_up_bang;

        evt.bdc1_xpos = m_bdc1_xpos;
        evt.bdc1_ypos = m_bdc1_ypos;
        evt.bdc2_xpos = m_bdc2_xpos;
        evt.bdc2_ypos = m_bdc2_ypos;

        evt.esl_vert_zpos = m_esl_vertex_zpos;
        evt.esr_vert_zpos = m_esr_vertex_zpos;
        //evt.es_vert_zpos = m_es_vertex_zpos;

        evt.esl_aang = m_esl_aang;
        evt.esl_bang = m_esl_bang;
        evt.esr_aang = m_esr_aang;
        evt.esr_bang = m_esr_bang;

        evt.esl_xpos = m_esl_xpos;
        evt.esl_ypos = m_esl_ypos;
        evt.esr_xpos = m_esr_xpos;
        evt.esr_ypos = m_esr_ypos;

        evt.esl_e_cal = m_esl_e_cal;
        evt.esr_e_cal = m_esr_e_cal;
        evt.esl_e_raw = m_esl_e_raw;
        evt.esr_e_raw = m_esr_e_raw;
        evt.esl_e_eff = m_esl_e_eff;
        evt.esr_e_eff = m_esr_e_eff;

        evt.esl_de_raw = m_esl_de_raw;
        evt.esr_de_raw = m_esr_de_raw;

        evt.esl_nhitwires = m_esl_nhitwires;
        evt.esl_xnhitwires = m_esl_xnhitwires;
        evt.esl_ynhitwires = m_esl_ynhitwires;
        evt.esr_nhitwires = m_esr_nhitwires;
        evt.esr_xnhitwires = m_esr_xnhitwires;
        evt.esr_ynhitwires = m_esr_ynhitwires;

        evt.fdc0_phi = m_fdc0_phi;
        evt.fdc0_theta = m_fdc0_theta;
        evt.fdc0_xpos = m_fdc0_xpos;
        evt.fdc0_ypos = m_fdc0_ypos;
        evt.fdc0_aang = m_fdc0_aang;
        evt.fdc0_bang = m_fdc0_bang;

        evt.s1dc_phi = m_s1dc_phi;
        evt.s1dc_theta = m_s1dc_theta;
        evt.s1dc_xpos = m_s1dc_xpos;
        evt.s1dc_ypos = m_s1dc_ypos;
        evt.s1dc_aang = m_s1dc_aang;
        evt.s1dc_bang = m_s1dc_bang;

        evt.tgt_dwn_aang = m_tgt_dwn_aang;
        evt.tgt_dwn_bang = m_tgt_dwn_bang;
        evt.tgt_dwn_vert_xpos = m_tgt_dwn_vert_xpos;
        evt.tgt_dwn_vert_ypos = m_tgt_dwn_vert_ypos;

        evt.fdc2_xpos = m_fdc2_xpos;
        evt.fdc2_ypos = m_fdc2_ypos;
        evt.fdc2_aang = m_fdc2_aang;
        evt.fdc2_bang = m_fdc2_bang;

        for (int i = 0; i < 7; i++) {
          evt.esl_naie_cal[i] = m_esl_naie_cal[i];
          evt.esr_naie_cal[i] = m_esr_naie_cal[i];
          evt.esl_naie_raw[i] = m_esl_naie_raw[i];
          evt.esr_naie_raw[i] = m_esr_naie_raw[i];
        }

        for (int i = 0; i < 24; i++) {
          evt.hodf_q[i] = m_hodf_q[i];
          evt.hodf_t[i] = m_hodf_t[i];
        }

        return evt;
      }

      void BindTreeVars(TTree& tree) {
        tree.SetBranchAddress("event_no", &m_event_no);

        tree.SetBranchAddress("triggers", &m_triggers);
        tree.SetBranchAddress("p_phi", &m_p_phi);
        tree.SetBranchAddress("p_theta", &m_p_theta);
        tree.SetBranchAddress("p_theta_eff", &m_p_theta_eff);
        tree.SetBranchAddress("p_e_eff", &m_p_e_eff);
        tree.SetBranchAddress("p_e_sim", &m_p_e_sim);
        tree.SetBranchAddress("p_e", &m_p_e);
        tree.SetBranchAddress("p_de", &m_p_de);

        tree.SetBranchAddress("esr_e", &m_esr_e_raw);
        tree.SetBranchAddress("esl_e", &m_esl_e_raw);
        tree.SetBranchAddress("esr_de", &m_esr_de_raw);
        tree.SetBranchAddress("esl_de", &m_esl_de_raw);

        tree.SetBranchAddress("esl_aang", &m_esl_aang);
        tree.SetBranchAddress("esl_bang", &m_esl_bang);
        tree.SetBranchAddress("esr_aang", &m_esr_aang);
        tree.SetBranchAddress("esr_bang", &m_esr_bang);

        tree.SetBranchAddress("esl_xpos", &m_esl_xpos);
        tree.SetBranchAddress("esl_ypos", &m_esl_ypos);
        tree.SetBranchAddress("esr_xpos", &m_esr_xpos);
        tree.SetBranchAddress("esr_ypos", &m_esr_ypos);

        tree.SetBranchAddress("esl_e_cal", &m_esl_e_cal);
        tree.SetBranchAddress("esr_e_cal", &m_esr_e_cal);

        tree.SetBranchAddress("esl_e_eff", &m_esl_e_eff);
        tree.SetBranchAddress("esr_e_eff", &m_esr_e_eff);

        tree.SetBranchAddress("esl_rdc_nhitwires", &m_esl_nhitwires);
        tree.SetBranchAddress("esl_rdc_xnhitwires", &m_esl_xnhitwires);
        tree.SetBranchAddress("esl_rdc_ynhitwires", &m_esl_ynhitwires);
        tree.SetBranchAddress("esr_rdc_nhitwires", &m_esr_nhitwires);
        tree.SetBranchAddress("esr_rdc_xnhitwires", &m_esr_xnhitwires);
        tree.SetBranchAddress("esr_rdc_ynhitwires", &m_esr_ynhitwires);

        tree.SetBranchAddress("esl_naie_cal", &m_esl_naie_cal);
        tree.SetBranchAddress("esr_naie_cal", &m_esr_naie_cal);
        tree.SetBranchAddress("esl_naie", &m_esl_naie_raw);
        tree.SetBranchAddress("esr_naie", &m_esr_naie_raw);

        tree.SetBranchAddress("s1dc_phi", &m_s1dc_phi);
        tree.SetBranchAddress("s1dc_theta", &m_s1dc_theta);
        tree.SetBranchAddress("he_theta_theor", &m_he_theta_theor);
        tree.SetBranchAddress("tgt_up_xpos", &m_tgt_up_vert_xpos);
        tree.SetBranchAddress("tgt_up_ypos", &m_tgt_up_vert_ypos);
        tree.SetBranchAddress("tgt_up_aang", &m_tgt_up_aang);
        tree.SetBranchAddress("tgt_up_bang", &m_tgt_up_bang);
        tree.SetBranchAddress("bdc1_xpos", &m_bdc1_xpos);
        tree.SetBranchAddress("bdc1_ypos", &m_bdc1_ypos);
        tree.SetBranchAddress("bdc2_xpos", &m_bdc2_xpos);
        tree.SetBranchAddress("bdc2_ypos", &m_bdc2_ypos);

        tree.SetBranchAddress("es_vertex_zpos", &m_es_vertex_zpos);
        tree.SetBranchAddress("esl_vertex_zpos", &m_esl_vertex_zpos);
        tree.SetBranchAddress("esr_vertex_zpos", &m_esr_vertex_zpos);
        tree.SetBranchAddress("hodf_q", &m_hodf_q);
        tree.SetBranchAddress("hodf_t", &m_hodf_t);

        tree.SetBranchAddress("fdc0_phi", &m_fdc0_phi);
        tree.SetBranchAddress("fdc0_theta", &m_fdc0_theta);
        tree.SetBranchAddress("fdc0_xpos", &m_fdc0_xpos);
        tree.SetBranchAddress("fdc0_ypos", &m_fdc0_ypos);
        tree.SetBranchAddress("fdc0_aang", &m_fdc0_aang);
        tree.SetBranchAddress("fdc0_bang", &m_fdc0_bang);

        tree.SetBranchAddress("s1dc_xpos", &m_s1dc_xpos);
        tree.SetBranchAddress("s1dc_ypos", &m_s1dc_ypos);
        tree.SetBranchAddress("s1dc_aang", &m_s1dc_aang);
        tree.SetBranchAddress("s1dc_bang", &m_s1dc_bang);

        tree.SetBranchAddress("tgt_dwn_xpos", &m_tgt_dwn_vert_xpos);
        tree.SetBranchAddress("tgt_dwn_ypos", &m_tgt_dwn_vert_ypos);
        tree.SetBranchAddress("tgt_dwn_aang", &m_tgt_dwn_aang);
        tree.SetBranchAddress("tgt_dwn_bang", &m_tgt_dwn_bang);

        tree.SetBranchAddress("fdc2_xpos", &m_fdc2_xpos);
        tree.SetBranchAddress("fdc2_ypos", &m_fdc2_ypos);
        tree.SetBranchAddress("fdc2_aang", &m_fdc2_aang);
        tree.SetBranchAddress("fdc2_bang", &m_fdc2_bang);
      }

      void CleanTreeVars() {
        m_event_no = k_inf;

        for (int i = 0; i < 8; i++) {
          m_triggers[i] = false;
        }

        m_p_phi = k_inf;
        m_p_theta = k_inf;
        m_p_aang = k_inf;
        m_p_bang = k_inf;
        m_p_theta_eff = k_inf;
        m_p_e_eff = k_inf;
        m_p_e_sim = k_inf;
        m_p_e = k_inf;
        m_p_de = k_inf;
        m_he_theta_theor = k_inf;
        m_tgt_up_vert_xpos = k_inf;
        m_tgt_up_vert_ypos = k_inf;
        m_tgt_up_aang = k_inf;
        m_tgt_up_bang = k_inf;
        m_bdc1_xpos = k_inf;
        m_bdc1_ypos = k_inf;
        m_bdc2_xpos = k_inf;
        m_bdc2_ypos = k_inf;
        m_es_vertex_zpos = k_inf;

        m_esl_vertex_zpos = k_inf;
        m_esr_vertex_zpos = k_inf;

        for (int i = 0; i < 24; i++) {
          m_hodf_q[i] = k_inf;
          m_hodf_t[i] = k_inf;
        }

        m_p_e_raw = k_inf;
        m_p_de_raw = k_inf;
        m_esr_e_raw = k_inf;
        m_esl_e_raw = k_inf;
        m_esr_de_raw = k_inf;
        m_esl_de_raw = k_inf;

        m_esl_aang = k_inf;
        m_esl_bang = k_inf;
        m_esr_aang = k_inf;
        m_esr_bang = k_inf;

        m_esl_xpos = k_inf;
        m_esl_ypos = k_inf;
        m_esr_xpos = k_inf;
        m_esr_ypos = k_inf;

        m_esl_e_cal = k_inf;
        m_esr_e_cal = k_inf;

        m_esl_e_eff = k_inf;
        m_esr_e_eff = k_inf;

        m_esl_nhitwires = k_inf;
        m_esl_xnhitwires = k_inf;
        m_esl_ynhitwires = k_inf;
        m_esr_nhitwires = k_inf;
        m_esr_xnhitwires = k_inf;
        m_esr_ynhitwires = k_inf;

        for (int i = 0; i < 7; i++) {
          m_esl_naie_cal[i] = k_inf;
          m_esr_naie_cal[i] = k_inf;
          m_esl_naie_raw[i] = k_inf;
          m_esr_naie_raw[i] = k_inf;
        }

        m_fdc0_phi = k_inf;
        m_fdc0_theta = k_inf;
        m_fdc0_xpos = k_inf;
        m_fdc0_ypos = k_inf;
        m_fdc0_aang = k_inf;
        m_fdc0_bang = k_inf;

        m_s1dc_phi = k_inf;
        m_s1dc_theta = k_inf;
        m_s1dc_xpos = k_inf;
        m_s1dc_ypos = k_inf;
        m_s1dc_aang = k_inf;
        m_s1dc_bang = k_inf;

        m_tgt_dwn_aang = k_inf;
        m_tgt_dwn_bang = k_inf;
        m_tgt_dwn_vert_xpos = k_inf;
        m_tgt_dwn_vert_ypos = k_inf;

        m_fdc2_xpos = k_inf;
        m_fdc2_ypos = k_inf;
        m_fdc2_aang = k_inf;
        m_fdc2_bang = k_inf;
      }

      bool is_well_defined(double value) {
        if (value == k_inf) {
          return false;
        } else  {
          return true;
        }
      }

      Double_t rad_to_degrees(double value_rad) {
        if (!is_well_defined(value_rad)) {
          return value_rad;
        }
        return value_rad * (1 / s13::gk_d2r);
      }

      Double_t mrad_to_degrees(double value_mrad) {
        return rad_to_degrees(value_mrad * 1e-3);
      }

      void ConvertValues() {
        m_p_phi = rad_to_degrees(m_p_phi);
        m_p_theta = rad_to_degrees(m_p_theta);
        m_p_theta_eff = rad_to_degrees(m_p_theta_eff);
        m_s1dc_phi = rad_to_degrees(m_s1dc_phi);
        m_s1dc_theta = rad_to_degrees(m_s1dc_theta);
        m_fdc0_phi = rad_to_degrees(m_fdc0_phi);
        m_fdc0_theta = rad_to_degrees(m_fdc0_theta);
        m_he_theta_theor = rad_to_degrees(m_he_theta_theor);

        m_fdc0_aang = rad_to_degrees(m_fdc0_aang);
        m_fdc0_bang = rad_to_degrees(m_fdc0_bang);
        m_s1dc_aang = rad_to_degrees(m_s1dc_aang);
        m_s1dc_bang = rad_to_degrees(m_s1dc_bang);

        // TODO: remove this and get values of ESL&ESR directly from the scattree
        m_p_e_raw = m_p_phi < 0 ? m_esl_e_raw : m_esr_e_raw;
        m_p_de_raw = m_p_phi < 0 ? m_esl_de_raw : m_esr_de_raw;
        m_p_aang = m_p_phi < 0 ? m_esl_aang : m_esr_aang;
        m_p_bang = m_p_phi < 0 ? m_esl_bang : m_esr_bang;
      }

      // TTree data vars
      int m_event_no;

      Bool_t m_triggers[8];
      Double_t m_p_phi;
      Double_t m_p_theta;
      Double_t m_p_aang;
      Double_t m_p_bang;
      Double_t m_p_theta_eff;
      Double_t m_p_e_eff;
      Double_t m_p_e_sim;
      Double_t m_p_e;
      Double_t m_p_e_raw;
      Double_t m_p_de;
      Double_t m_p_de_raw;
      Double_t m_s1dc_phi;
      Double_t m_s1dc_theta;
      Double_t m_he_theta_theor;

      Double_t m_tgt_up_vert_xpos;
      Double_t m_tgt_up_vert_ypos;
      Double_t m_tgt_up_aang;
      Double_t m_tgt_up_bang;
      Double_t m_bdc1_xpos;
      Double_t m_bdc1_ypos;
      Double_t m_bdc2_xpos;
      Double_t m_bdc2_ypos;

      Double_t m_es_vertex_zpos;

      Double_t m_esl_vertex_zpos;
      Double_t m_esr_vertex_zpos;

      Double_t m_hodf_q[24];
      Double_t m_hodf_t[24];

      Double_t m_esr_e_raw;
      Double_t m_esl_e_raw;
      Double_t m_esr_de_raw;
      Double_t m_esl_de_raw;

      Double_t m_esl_aang;
      Double_t m_esl_bang;
      Double_t m_esr_aang;
      Double_t m_esr_bang;

      Double_t m_esl_xpos;
      Double_t m_esl_ypos;
      Double_t m_esr_xpos;
      Double_t m_esr_ypos;

      Double_t m_esl_e_cal;
      Double_t m_esr_e_cal;

      Double_t m_esl_e_eff;
      Double_t m_esr_e_eff;

      Double_t m_esl_naie_cal[7];
      Double_t m_esr_naie_cal[7];
      Double_t m_esl_naie_raw[7];
      Double_t m_esr_naie_raw[7];

      Int_t m_esl_nhitwires;
      Int_t m_esl_xnhitwires;
      Int_t m_esl_ynhitwires;
      Int_t m_esr_nhitwires;
      Int_t m_esr_xnhitwires;
      Int_t m_esr_ynhitwires;

      Double_t m_fdc0_phi;
      Double_t m_fdc0_theta;
      Double_t m_fdc0_xpos;
      Double_t m_fdc0_ypos;
      Double_t m_fdc0_aang;
      Double_t m_fdc0_bang;

      Double_t m_s1dc_xpos;
      Double_t m_s1dc_ypos;
      Double_t m_s1dc_aang;
      Double_t m_s1dc_bang;

      Double_t m_tgt_dwn_aang;
      Double_t m_tgt_dwn_bang;
      Double_t m_tgt_dwn_vert_xpos;
      Double_t m_tgt_dwn_vert_ypos;

      Double_t m_fdc2_xpos;
      Double_t m_fdc2_ypos;
      Double_t m_fdc2_aang;
      Double_t m_fdc2_bang;

      s13::misc::MessageLogger log_;
    };



    /*
      Helper structure for thread function.
    */
    struct ParallelTreeWalkerThreadArgs {
      ScatteringTreeWalker* p_tree_walker;
      TChain* p_chain;
    };

    /*
      Thread function
    */
    void* parallel_tree_walker_thread_fn(void *arg);

    /*
      Class to run multiple TreeWalkers in separate threads.
    */
    class ParallelTreeWalker {

    public:

      /*
        Constructs ParallelTreeWalker.

        Params:
          tree_walker: A reference instance of ScatteringTreeWalker to use for running
            on one of chains segments. It must support copy assignment.
          chains: Series of TChain* on which a dedicated instance of passed
            tree_walker will be run.
      */
      ParallelTreeWalker(std::shared_ptr<ScatteringTreeWalker> tree_walker,
                         std::vector<TChain*> chains) :
        reference_instance_{tree_walker}, chains_{chains},
        log_{"ParallelTreeWalker"}  {
          using tree_walker_ptr_t = std::shared_ptr<ScatteringTreeWalker>;
          log_.debug("Contructing instances of tree walkers...");
          for (size_t i = 0; i < chains_.size(); i++) {
            auto walker =
              tree_walker_ptr_t(new ScatteringTreeWalker(*tree_walker));
            walkers_.push_back(walker);
            log_.debug("Instance #%i constructed", i + 1);
          }
        }

      /*
        Walks through chain segments passed in constructor and
        merges the results into one ScatteringTreeWalker which can be
        obtained through a call to reference_instance().
      */
      void walk() {
        log_.debug("Starting walk()...");
        std::vector<std::future<void> > futs;
        std::vector<std::thread> threads;

        log_.debug("Creating async tasks for tree walkers...");
        long i = 0;
        for (std::shared_ptr<ScatteringTreeWalker> w : walkers_) {
          log_.debug("Starting thread #%d...", i + 1);
          TChain* chain_segment = chains_.at(i);
          std::thread t([w, chain_segment] () {
              w->walk(*chain_segment);
            });
          threads.push_back(std::move(t));

          ++i;
        }

        log_.debug("Waiting for results...");
        i = 0;
        for (auto& t : threads) {
          t.join();
          log_.debug("Joined thread #%d", i + 1);
          ++i;
        }

        merge();
      }

      /*
        Returns a reference instance of ScatteringTreeWalker. After calling walk()
        it will contain merged results of algorithms which were ran on a
        chains segments.
      */
      std::shared_ptr<ScatteringTreeWalker> reference_instance() {
        return reference_instance_;
      }

    private:

      void merge() {
        log_.debug("Merging results from async tasks walkers into"
                   " reference walker...");
        int counter = 1;
        for (auto w : walkers_) {
          for (int i = 0; i < w->algorithms_number(); i++) {
            auto& alg = w->algorithm_at(i);
            reference_instance_->algorithm_at(i)->merge(*alg);
          }
          log_.debug("Merged from tree walker #%d", counter);
          ++counter;
        }
      }

      std::shared_ptr<ScatteringTreeWalker> reference_instance_;
      std::vector<TChain*> chains_;

      std::vector<std::shared_ptr<ScatteringTreeWalker> > walkers_;

      s13::misc::MessageLogger log_;
    };

  }
}
