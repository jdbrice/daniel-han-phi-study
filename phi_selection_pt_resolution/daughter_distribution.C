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

void daughter_distribution() {
  //    Create a file for output (saving) of histograms we compute
  TFile *fo = new TFile("inpuit.root", "RECREATE");

  // Create a ROOT 1D histogram
  TH1F *hdca = new TH1F("dca", "Run19 Au+Au data;dca (cm); counts", 500, 0, 5);
  // Create a ROOT 1D histogram for parent Mass
  TH1F *mass_rc = new TH1F(
      "parent_mass", "Run19 Au+Au data;parent_mass (GeV); counts", 100, 0.98, 1.1);

  // Create a ROOT 2D histogram for Pion/Kaon distributions
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
    // get their distance to closest approach ( info about track pair)
    double dca1 = pair->d1_mDCA;
    double dca2 = pair->d2_mDCA;

    hdca->Fill(dca1);
    hdca->Fill(dca2);

    //        if ( fabs(pair->d1_mNSigmaKaon) < 3. && fabs(pair->d2_mNSigmaKaon)
    //        < 3. && pair->d1_mPt > 0.01 && pair->d2_mPt > 0.01 ){ // event and
    //        track selection will go here

    //      only select += pair
    if (pair->mChargeSum == 0.) {
      // Setup the Lorentz Vectors of the daugher tracks
      lv1.SetPtEtaPhiM(pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi,
                       0.493); // we use Kaon mass here.
      lv2.SetPtEtaPhiM(pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.493);

      // compute parent particle lorentz vector from daughters
      // use conservation of momentum and energy to set the four-momentum
      // vector for parent particle.
      lv = lv1 + lv2;

      // selecting Kaons as daughter particles while rejecting Pions
      // selecting only low lv.Pt() as we investigated using the code blow
      // that lv.Pt() should be around 0.1 GeV
      // selecting only masses close to m_\phi
      if (fabs(pair->d1_mNSigmaKaon) < 5 && fabs(pair->d2_mNSigmaKaon) < 5 &&
          fabs(pair->d1_mNSigmaPion) > 5 && fabs(pair->d2_mNSigmaPion) > 5)
      {
        mass_rc->Fill(lv.M());
      }
    } // selection
  }   // loop on events

  fo->cd();

  makeCan();
  mass_rc->Draw();
  fo->Write();
}
