

#include "dataset.h"
#include "srccommon.hpp"


using namespace s13::ana;
using namespace s13::common;

//common::basic_istream::char_traits::CsvRow

static TChain* MakeChain(std::string chainName, std::vector<int> runNumbers,
                 std::string runFilenameSuffix) {
    TChain *chain = new TChain(chainName.c_str());

    for (int runNumber : runNumbers) {
      std::string runFilename = "rootfiles/run";
      runFilename = runFilename + std::to_string(runNumber);
      runFilename = runFilename + "-" + runFilenameSuffix + ".root";
      std::cout << "Adding file to chain: " << runFilename << std::endl;
      chain->Add(runFilename.c_str());
    }

    return chain;
}

SimpleDataset s13::ana::MakeDataset(std::vector<int> runNumbers,
                                    std::vector<DetectorTypes> detTypes) {
  // main chain, contains trigger bit data, always present
  TChain *main_chain = MakeChain("BRIPSPlaTree", runNumbers, "BRIPSPla");
  std::vector<TChain*> chains;

  for (auto detType : detTypes) {
    TChain *chain = nullptr;
    switch (detType) {
    case DetectorTypes::BripsPla:
      // pass - this is out main chain - always present
      break;
    case DetectorTypes::Bdcs:
      chain = MakeChain("BDCsTree", runNumbers, "BDCs");
      main_chain->AddFriend(chain->GetName());
      break;
    case DetectorTypes::Dsdcs:
      chain = MakeChain("DSDCsTree", runNumbers, "DSDCs");
      main_chain->AddFriend(chain->GetName());
      break;
    case DetectorTypes::Fdc2:
      chain = MakeChain("FDC2Tree", runNumbers, "FDC2");
      main_chain->AddFriend(chain->GetName());
      break;
    case DetectorTypes::Espri:
      chain = MakeChain("ESPRITree", runNumbers, "ESPRI");
      main_chain->AddFriend(chain->GetName());
      break;
    case DetectorTypes::Hodf:
      chain = MakeChain("HODFTree", runNumbers, "HODF");
      main_chain->AddFriend(chain->GetName());
      break;
    default:
      std::cerr << "ERROR: wrong detector type passed!" << std::endl;
      exit(1);
      break;
    }
  }

  SimpleDataset ds{runNumbers, main_chain};
  return ds;
}

/*
int main(int argc, char **argv) {

    // if (argc != 2) {
    //     std::cerr << "ERROR: This program reqiures run number as first argument!" << std::endl;
    //     exit(1);
    // }
    // everything should be relative to our home directory

    auto ana_root_dir = GetAnaRootDir();
    std::cout << "Working directory will be set to ana root: " << ana_root_dir \
        << std::endl;
    if (chdir(ana_root_dir.c_str()) != 0) {
        std::cerr << "ERROR: Unable to set program working directory!" << \
            std::endl;
        exit(1);
    }

    TApplication theApp("tapp", &argc, argv);
    auto ds = MakeDataset(std::vector<int> {272,273});

    TCanvas *canvas = new TCanvas("c1");
    canvas->Divide(2,2);

    canvas->cd(1);
    canvas->SetLogz();
    ds.GetTree()->Draw("TargetUp_XPos:TargetUp_YPos >>h(400,-40,40,400,-40,40)",
                        "", "colz");
    canvas->cd(2);
    canvas->SetLogz();
    ds.GetTree()->Draw("TargetDwn_XPos:TargetDwn_YPos >>h(400,-40,40,400,-40,40)",
                        "", "colz");
    canvas->cd(3);
    canvas->SetLogz();
    ds.GetTree()->Draw("ESPRIPlaR_QCal:ESPRIPlaR_TCal >>h(400,0,12000,400,0,5000)",
                       "", "colz");

    canvas->Show();
    theApp.Run();

    return 0;
}
*/
