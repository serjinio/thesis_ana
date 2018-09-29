/**
 *   \file treewalk.hpp
 *   \brief Function for walking TTrees with algs
 *
 */

#pragma once

#include <assert.h>

#include <TString.h>
#include <TChain.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "consts.hpp"
#include "commonalg.hpp"
#include "elacuts.hpp"
#include "treewalker.hpp"
#include "msgoutput.hpp"


namespace s13 {
  namespace ana {

    template <typename AlgT>
    std::shared_ptr<AlgT> walk_alg(const ElasticScatteringCuts& cuts,
                                   const std::vector<TChain*>& dataset,
                                   bool use_mt = false) {
      using PtrTreeWalker = std::shared_ptr<s13::ana::ScatteringTreeWalker>;
      s13::misc::MessageLogger logger("s13::ana::walk_alg()");
      std::shared_ptr<AlgT> alg(new AlgT(cuts));
      auto tree_walker = PtrTreeWalker(new ScatteringTreeWalker());
      tree_walker->add(alg);

      if (!use_mt) {
        logger.debug("Starting algorithm walk (single thread)...");
        tree_walker->walk(*dataset[0]);
      } else {
        logger.debug("Starting algorithm walk (multiple threads)...");
        s13::ana::ParallelTreeWalker ptree_walker(tree_walker, dataset);
        ptree_walker.walk();
      }

      return alg;
    }

    template <typename AlgT, typename ...ArgTs>
    std::shared_ptr<AlgT> walk_alg2(const std::vector<TChain*>& dataset,
                                    bool use_mt,
                                    ArgTs&& ...args) {
      using PtrTreeWalker = std::shared_ptr<s13::ana::ScatteringTreeWalker>;
      s13::misc::MessageLogger logger("s13::ana::walk_alg(...)");
      std::shared_ptr<AlgT> alg(new AlgT(std::forward<ArgTs>(args)...));
      auto tree_walker = PtrTreeWalker(new ScatteringTreeWalker());
      tree_walker->add(alg);

      if (!use_mt) {
        logger.debug("Starting algorithm walk (single thread)...");
        tree_walker->walk(*dataset[0]);
      } else {
        logger.debug("Starting algorithm walk (multiple threads)...");
        s13::ana::ParallelTreeWalker ptree_walker(tree_walker, dataset);
        ptree_walker.walk();
      }

      return alg;
    }
  }
}
