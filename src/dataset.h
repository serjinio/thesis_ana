/*
  Dataset classes represent data selected for processing.

  Construction:
    A dataset could be constructed from a range of runs, like: 134-138,143,146.
    This should result in a dataset with all listed runs numbers included.
 */

#pragma once

#include <iostream>
#include <vector>

#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"


namespace s13 {
    namespace ana {

        enum class DetectorTypes {
            BripsPla, Bdcs, Dsdcs, Fdc2, Espri, Hodf
        };

        const std::vector<DetectorTypes> ALL_DETECTORS = {
            DetectorTypes::BripsPla, DetectorTypes::Bdcs, DetectorTypes::Dsdcs,
            DetectorTypes::Fdc2, DetectorTypes::Espri, DetectorTypes::Hodf };

        class SimpleDataset {
        public:
            SimpleDataset(std::vector<int> runNumbers, TTree* rootTree) :
                mp_RootTree(rootTree), m_RunNumbers(runNumbers) {}

            ~SimpleDataset() {
                delete mp_RootTree;
            }

            TTree* GetTree() {
                return mp_RootTree;
            }

            std::vector<int> GetRunNumbers() {
                return m_RunNumbers;
            }
        private:
            TTree *mp_RootTree;
            std::vector<int> m_RunNumbers;
        };

        SimpleDataset MakeDataset(std::vector<int> runNumbers,
                                  std::vector<DetectorTypes> detTypes =
                                   ALL_DETECTORS);

    }
}


// Local Variables:
// mode: c++
// End:
