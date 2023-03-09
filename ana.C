// A simple TTreeReader use: read data from hsimple.root (written by hsimple.C)
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "TColor.h"

#include "FemtoPairFormat.h"

// ClassImp( FemtoPair );

int ican = 0;
void makeCan() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican++ ), "", 900, 600 );
    can->SetTopMargin(0.06);
    can->SetRightMargin(0.01);
}

void ana() {
    //    Create a file for output (saving) of histograms we compute
    TFile * fo = new TFile( "data.root", "RECREATE" );

    // Create a ROOT 1D histogram
    auto * hdca = new TH1F("dca", "Run19 Au+Au data;dca (cm); counts", 500, 0, 5);

    // Open the file containing the tree (INPUT data).
    TFile *myFile = TFile::Open("input.root");
    //This setsup the reader, access the data
    TTreeReader myReader("PairDst", myFile);

    // We need this so ROOT understands the data format
    gSystem->Load("FemtoPairFormat_h.so");
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

    // Lorentz vectors (4-vectors) for kinematics in special relativity
    TLorentzVector lv1, lv2, lv;

   // // Loop over all entries of the TTree or TChain.
    while (myReader.Next()) {
        // get their distance to closest approach ( info about track pair)
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;

        hdca->Fill( dca1 );
        hdca->Fill( dca2 );


        if ( true ){ // event and track selection will go here
            // Setup the Lorentz Vectors of the daugher tracks
            lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.493 );
            lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.493 );

            // compute parent particle lorentz vector from daughters

            
        } // selection
    } // loop on events

    fo->cd();

    makeCan();
    // draw histogram here
    
    hdca->Draw();
    gPad->SetLogy(1); // set semilog-y

    // save plot to a PDF
    gPad->Print( "plot0.pdf" );

    // write all histograms to output data file
    fo->Write();
  }             