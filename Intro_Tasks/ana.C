// A simple TTreeReader use: read data from hsimple.root (written by hsimple.C)
#include "TFile.h"
#include "TH1F.h"
#include <TSystem.h>
#include "TH2F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "TColor.h"
#include <cmath>
#include <math.h>
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
    TFile * fo = new TFile( "inpuit.root", "RECREATE" );

    // Create a ROOT 1D histogram
    TH1F * hdca = new TH1F("dca", "Run19 Au+Au data;dca (cm); counts", 500, 0, 5);
    // Create a ROOT 1D histogram for parent Mass
    TH1F * mass = new TH1F("parent_mass", "Run19 Au+Au data;parent_mass (GeV); counts", 500, 0, 5);
    TH2F * pt_mass = new TH2F("pt_mass",";mass(GeV);pt(GeV)", 5000,0,5, 500,0,5);
    // Open the file containing the tree (INPUT data).
    TFile *myFile = TFile::Open("input.root");
    //This setsup the reader, access the data
    TTreeReader myReader("PairDst", myFile);

    // We need this so ROOT understands the data format
    gSystem->Load("FemtoPairFormat_h.so");
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

    // Lorentz vectors (4-vectors) for kinematics in special relativity
    TLorentzVector lv1, lv2, lv;

    double parent_mass;
   // // Loop over all entries of the TTree or TChain.
    while (myReader.Next()) {
        // get their distance to closest approach ( info about track pair)
        double parent_mass = 0;
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;

        hdca->Fill( dca1 );
        hdca->Fill( dca2 );


//        if ( fabs(pair->d1_mNSigmaKaon) < 3. && fabs(pair->d2_mNSigmaKaon) < 3. && pair->d1_mPt > 0.01 && pair->d2_mPt > 0.01 ){ // event and track selection will go here
        if(pair->mChargeSum == 0.){
            // Setup the Lorentz Vectors of the daugher tracks
            lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.493 );
            lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.493 );

            // compute parent particle lorentz vector from daughters
            // chis is to be updated. I believe there are easier way to compute 
            //
            // compute the four-momentum vector for particle 1
            double px1 = cos(pair->d1_mPhi) * (pair->d1_mPt);
            double py1 = sin(pair->d1_mPhi) * (pair->d1_mPt);
            double pz1 = sinh(pair-> d1_mEta) * (pair-> d1_mPt);
            double e1 = sqrt(pow(px1, 2.)+pow(py1, 2.)+pow(pz1,2.) + pow(0.139,2.));
            // compute the four-momentum vector for particle 2
            double px2 = cos(pair->d2_mPhi) * (pair->d2_mPt);
            double py2 = sin(pair->d2_mPhi) * (pair->d2_mPt);
            double pz2 = sinh(pair-> d2_mEta) * (pair-> d2_mPt);
            double e2 = sqrt(pow(px2, 2.)+pow(py2, 2.)+pow(pz2,2.) + pow(0.139,2.));

            // use conservation of momentum and energy to set the four-momentum
            // vector for parent particle.
            lv.SetPxPyPzE(px1+px2, py1 + py2, pz1 + pz2, e1 + e2);
            lv = lv1 + lv2;// extract the parent mass from the parent four vector
            parent_mass = lv.E() - fabs(lv.Px()) - fabs(lv.Py()) - fabs(lv.Pz());
            parent_mass = sqrt(pow(lv.E(),2) - pow(lv.Px(), 2.)+pow(lv.Py(), 2.)+pow(lv.Pz(),2.));

            mass->Fill( parent_mass );

            pt_mass -> Fill(lv.M(), lv.Pt());
        } // selection
    } // loop on events

    fo->cd();

    makeCan();
    // draw histogram here
    // hdca->Draw();
    mass->Draw();
    gPad->Print( "parent_mass.pdf" );
    gPad->SetLogy(1); // set semilog-y

    // save plot to a PDF
    gPad->Print( "parent_mass_log.pdf" );
    pt_mass-> Draw("colz");
    gPad->SetLogy(0); // set semilog-y
    gPad->SetLogz(1); // set semilog-z
    // write all histograms to output data file
    fo->Write();
  }             
