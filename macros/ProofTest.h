//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 25 20:10:17 2016 by ROOT version 6.09/01
// from TChain scattree/
//////////////////////////////////////////////////////////

#ifndef ProofTest_h
#define ProofTest_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class ProofTest : public TSelector {
public :
  TTreeReader fReader;  //!the tree reader
  TTree *fChain = 0;   //!pointer to the analyzed TTree or TChain

  // Readers to access the data (delete the ones you do not need).
  TTreeReaderValue<Int_t> event_no = {fReader, "event_no"};
  TTreeReaderArray<Bool_t> triggers = {fReader, "triggers"};
  TTreeReaderValue<Double_t> f7_q = {fReader, "f7_q"};
  TTreeReaderValue<Double_t> f7_t = {fReader, "f7_t"};
  TTreeReaderValue<Double_t> sbt1_q = {fReader, "sbt1_q"};
  TTreeReaderValue<Double_t> sbt1_t = {fReader, "sbt1_t"};
  TTreeReaderValue<Double_t> sbt2_q = {fReader, "sbt2_q"};
  TTreeReaderValue<Double_t> sbt2_t = {fReader, "sbt2_t"};
  TTreeReaderValue<Double_t> sbv_q = {fReader, "sbv_q"};
  TTreeReaderValue<Double_t> sbv_t = {fReader, "sbv_t"};
  TTreeReaderValue<Double_t> bdc1_xpos = {fReader, "bdc1_xpos"};
  TTreeReaderValue<Double_t> bdc1_ypos = {fReader, "bdc1_ypos"};
  TTreeReaderValue<Double_t> bdc2_xpos = {fReader, "bdc2_xpos"};
  TTreeReaderValue<Double_t> bdc2_ypos = {fReader, "bdc2_ypos"};
  TTreeReaderValue<Double_t> tgt_up_aang = {fReader, "tgt_up_aang"};
  TTreeReaderValue<Double_t> tgt_up_bang = {fReader, "tgt_up_bang"};
  TTreeReaderValue<Double_t> tgt_up_theta = {fReader, "tgt_up_theta"};
  TTreeReaderValue<Double_t> tgt_up_phi = {fReader, "tgt_up_phi"};
  TTreeReaderValue<Double_t> tgt_up_xpos = {fReader, "tgt_up_xpos"};
  TTreeReaderValue<Double_t> tgt_up_ypos = {fReader, "tgt_up_ypos"};
  TTreeReaderValue<Double_t> fdc0_xpos = {fReader, "fdc0_xpos"};
  TTreeReaderValue<Double_t> fdc0_ypos = {fReader, "fdc0_ypos"};
  TTreeReaderValue<Double_t> fdc0_aang = {fReader, "fdc0_aang"};
  TTreeReaderValue<Double_t> fdc0_bang = {fReader, "fdc0_bang"};
  TTreeReaderValue<Double_t> fdc0_theta = {fReader, "fdc0_theta"};
  TTreeReaderValue<Double_t> fdc0_phi = {fReader, "fdc0_phi"};
  TTreeReaderValue<Double_t> s1dc_xpos = {fReader, "s1dc_xpos"};
  TTreeReaderValue<Double_t> s1dc_ypos = {fReader, "s1dc_ypos"};
  TTreeReaderValue<Double_t> s1dc_aang = {fReader, "s1dc_aang"};
  TTreeReaderValue<Double_t> s1dc_bang = {fReader, "s1dc_bang"};
  TTreeReaderValue<Double_t> s1dc_theta = {fReader, "s1dc_theta"};
  TTreeReaderValue<Double_t> s1dc_phi = {fReader, "s1dc_phi"};
  TTreeReaderValue<Double_t> tgt_dwn_aang = {fReader, "tgt_dwn_aang"};
  TTreeReaderValue<Double_t> tgt_dwn_bang = {fReader, "tgt_dwn_bang"};
  TTreeReaderValue<Double_t> tgt_dwn_theta = {fReader, "tgt_dwn_theta"};
  TTreeReaderValue<Double_t> tgt_dwn_phi = {fReader, "tgt_dwn_phi"};
  TTreeReaderValue<Double_t> tgt_dwn_xpos = {fReader, "tgt_dwn_xpos"};
  TTreeReaderValue<Double_t> tgt_dwn_ypos = {fReader, "tgt_dwn_ypos"};
  TTreeReaderValue<Double_t> esl_xpos = {fReader, "esl_xpos"};
  TTreeReaderValue<Double_t> esl_ypos = {fReader, "esl_ypos"};
  TTreeReaderValue<Double_t> esl_aang = {fReader, "esl_aang"};
  TTreeReaderValue<Double_t> esl_bang = {fReader, "esl_bang"};
  TTreeReaderValue<Double_t> esl_aang_tgt = {fReader, "esl_aang_tgt"};
  TTreeReaderValue<Double_t> esl_bang_tgt = {fReader, "esl_bang_tgt"};
  TTreeReaderValue<Double_t> esl_de = {fReader, "esl_de"};
  TTreeReaderValue<Double_t> esl_de_cal = {fReader, "esl_de_cal"};
  TTreeReaderValue<Double_t> esl_tof = {fReader, "esl_tof"};
  TTreeReaderArray<Double_t> esl_naie = {fReader, "esl_naie"};
  TTreeReaderArray<Double_t> esl_naie_cal = {fReader, "esl_naie_cal"};
  TTreeReaderValue<Double_t> esl_e = {fReader, "esl_e"};
  TTreeReaderValue<Double_t> esl_e_cal = {fReader, "esl_e_cal"};
  TTreeReaderValue<Double_t> esl_vertex_zpos = {fReader, "esl_vertex_zpos"};
  TTreeReaderValue<Double_t> esl_p_theta = {fReader, "esl_p_theta"};
  TTreeReaderValue<Double_t> esl_p_phi = {fReader, "esl_p_phi"};
  TTreeReaderValue<Double_t> esr_xpos = {fReader, "esr_xpos"};
  TTreeReaderValue<Double_t> esr_ypos = {fReader, "esr_ypos"};
  TTreeReaderValue<Double_t> esr_aang = {fReader, "esr_aang"};
  TTreeReaderValue<Double_t> esr_bang = {fReader, "esr_bang"};
  TTreeReaderValue<Double_t> esr_aang_tgt = {fReader, "esr_aang_tgt"};
  TTreeReaderValue<Double_t> esr_bang_tgt = {fReader, "esr_bang_tgt"};
  TTreeReaderValue<Double_t> esr_de = {fReader, "esr_de"};
  TTreeReaderValue<Double_t> esr_de_cal = {fReader, "esr_de_cal"};
  TTreeReaderValue<Double_t> esr_tof = {fReader, "esr_tof"};
  TTreeReaderArray<Double_t> esr_naie = {fReader, "esr_naie"};
  TTreeReaderArray<Double_t> esr_naie_cal = {fReader, "esr_naie_cal"};
  TTreeReaderValue<Double_t> esr_e = {fReader, "esr_e"};
  TTreeReaderValue<Double_t> esr_e_cal = {fReader, "esr_e_cal"};
  TTreeReaderValue<Double_t> esr_vertex_zpos = {fReader, "esr_vertex_zpos"};
  TTreeReaderValue<Double_t> esr_p_theta = {fReader, "esr_p_theta"};
  TTreeReaderValue<Double_t> esr_p_phi = {fReader, "esr_p_phi"};
  TTreeReaderValue<Double_t> es_vertex_zpos = {fReader, "es_vertex_zpos"};
  TTreeReaderValue<Double_t> p_theta = {fReader, "p_theta"};
  TTreeReaderValue<Double_t> p_theta_eff = {fReader, "p_theta_eff"};
  TTreeReaderValue<Double_t> p_phi = {fReader, "p_phi"};
  TTreeReaderValue<Double_t> p_e = {fReader, "p_e"};
  //TTreeReaderValue<Double_t> p_e_eff = {fReader, "p_e_eff"};
  TTreeReaderValue<Double_t> p_e_sim = {fReader, "p_e_sim"};
  TTreeReaderValue<Double_t> p_de = {fReader, "p_de"};
  TTreeReaderValue<Double_t> p_de_sim = {fReader, "p_de_sim"};
  TTreeReaderValue<Double_t> p_tof = {fReader, "p_tof"};
  TTreeReaderValue<Double_t> fdc2_xpos = {fReader, "fdc2_xpos"};
  TTreeReaderValue<Double_t> fdc2_ypos = {fReader, "fdc2_ypos"};
  TTreeReaderValue<Double_t> fdc2_aang = {fReader, "fdc2_aang"};
  TTreeReaderValue<Double_t> fdc2_bang = {fReader, "fdc2_bang"};
  TTreeReaderArray<Double_t> hodf_q = {fReader, "hodf_q"};
  TTreeReaderArray<Double_t> hodf_t = {fReader, "hodf_t"};
  TTreeReaderValue<Double_t> he_theta_theor = {fReader, "he_theta_theor"};


  ProofTest(TTree * /*tree*/ = 0) {}

  virtual ~ProofTest() {}

  virtual Int_t Version() const { return 2; }

  virtual void Begin(TTree *tree);

  virtual void SlaveBegin(TTree *tree);

  virtual void Init(TTree *tree);

  virtual Bool_t Notify();

  virtual Bool_t Process(Long64_t entry);

  virtual Int_t GetEntry(Long64_t entry, Int_t getall = 0) {
    return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
  }

  virtual void SetOption(const char *option) { fOption = option; }

  virtual void SetObject(TObject *obj) { fObject = obj; }

  virtual void SetInputList(TList *input) { fInput = input; }

  virtual TList *GetOutputList() const { return fOutput; }

  virtual void SlaveTerminate();

  virtual void Terminate();

  TH1 *hist_p_theta;
  TH1 *hist_he_theta;

//ClassDef(ProofTest, 0);


};

#endif

#ifdef ProofTest_cxx
void ProofTest::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t ProofTest::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef ProofTest_cxx
