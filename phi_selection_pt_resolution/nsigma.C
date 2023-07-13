// A simple TTreeReader use: read data from hsimple.root (written by hsimple.C)
#include "FemtoPairFormat.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <TDirectory.h>
#include <TH2.h>
#include <TSystem.h>
#include <cmath>
#include <math.h>

// ClassImp( FemtoPair );

int ican = 0;
void makeCan() {
  TCanvas *can = new TCanvas(TString::Format("can%d", ican++), "", 900, 600);
  can->SetTopMargin(0.06);
  can->SetRightMargin(0.01);
}

void nsigma() {

  // histograms for all tracks
  TH1F *all_DCA = new TH1F("All Tracks", "Run 19 A+A DCA All Tracks;DCA;counts",
                           100, 0., 3.);

  TH1F *all_Npion = new TH1F(
      "All Tracks", "Run 19 A+A N#sigma#pi All Tracks;(N#sigmaPion);counts",
      100, -10, 10);

  TH2F *all_pt_NPion = new TH2F(
      "All Tracks", "Run 19 A+A N#sigmaPion All Tracks;P_{T};N#sigmaPion", 100,
      0.1, 1.2, 100, -4, 10);
  // Open the file containing the tree (INPUT data).
  TFile *myFile = TFile::Open("input.root");

  // This setup the reader, access the data
  TTreeReader myReader("PairDst", myFile);

  // We need this so ROOT understands the data format
  gSystem->Load("FemtoPairFormat_h.so");
  TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

  // Lorentz vectors (4-vectors) for kinematics in special relativity
  TLorentzVector lv1, lv2, lv;

  // Loop over all entries of the TTree or TChain.
  while (myReader.Next()) {

    all_Npion->Fill(pair->d1_mNSigmaPion);
    all_Npion->Fill(pair->d2_mNSigmaPion);

    all_pt_NPion->Fill(pair->d1_mPt, pair->d1_mNSigmaPion);
    all_pt_NPion->Fill(pair->d2_mPt, pair->d2_mNSigmaPion);
  } // loop on events

  makeCan();
  all_Npion->Draw();
  all_pt_NPion->Draw("colz");
  all_pt_NPion->FitSlicesY();
  TH1F *all_pt_NPion_1 = (TH1F *)gDirectory->Get("all_pt_NPion_1");
  if (all_pt_NPion_1) {
    std::cout << "all_pt_NPion_1 created successfully." << std::endl;
    makeCan();
    all_pt_NPion_1->Draw();
  } else {
    std::cout << "all_pt_NPion_1 is not created." << std::endl;
  }
}
