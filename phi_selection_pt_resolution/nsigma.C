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
      0., 1.2, 100, -8, 8);
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

    // kaon candidate filling with PId selection
    if (pair->d1_mPt > 0.06 && pair->d2_mPt > 0.06 &&
        abs(pair->d1_mNSigmaKaon) < 5 && abs(pair->d1_mNSigmaPion) >= 5 &&
        abs(pair->d2_mNSigmaKaon) < 5 && abs(pair->d2_mNSigmaPion) >= 5) {

      // reconstrauction of parent lorentz vector
      // assuming kaon mass
      TLorentzVector lv1, lv2, lv;
      lv1.SetPtEtaPhiM(pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.493);
      lv2.SetPtEtaPhiM(pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.493);
      lv = lv1 + lv2;

      // reco level filling
      if (pair->mChargeSum == 0) {
      } else if (pair->mChargeSum == -2) {

      } else if (pair->mChargeSum == 2) {
      }
    }

  } // loop on events

  makeCan();
  all_Npion->Draw();
  all_pt_NPion->Draw("colz");
  all_pt_NPion->FitSlicesY();
}
