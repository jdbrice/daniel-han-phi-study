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

void daughter_distribution() {
    //    Create a file for output (saving) of histograms we compute
    TFile * fo = new TFile( "inpuit.root", "RECREATE" );

    // Create a ROOT 1D histogram
    TH1F * hdca = new TH1F("dca", "Run19 Au+Au data;dca (cm); counts", 500, 0, 5);
    // Create a ROOT 1D histogram for parent Mass
    TH1F * mass = new TH1F("parent_mass", "Run19 Au+Au data;parent_mass (GeV); counts", 50, 1, 2);
    TH2F * pt_mass = new TH2F("pt_mass",";mass(GeV);pt(GeV)", 500,0,5, 500,0,1.5);

    // Create a ROOT 2D histogram for Pion/Kaon distributions
    TH2F * track1_sigma_KP= new TH2F("Track1",";d1_N_sigma_kaon;d1_N_sigma_pion", 500, -5 , 5, 500, -5, 5);
    TH2F * track2_sigma_KP= new TH2F("Track2",";d2_N_sigma_kaon;d2_N_sigma_pion", 500, -5 , 5, 500, -5, 5);
    TH1F * phi_pt = new TH1F("phi_pt", ";Phi_pt(GeV); Count", 50, 0, 1.5);
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
    // Loop over all entries of the TTree or TChain.
    while (myReader.Next()) {
        // get their distance to closest approach ( info about track pair)
        double parent_mass = 0;
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;

        hdca->Fill( dca1 );
        hdca->Fill( dca2 );


//        if ( fabs(pair->d1_mNSigmaKaon) < 3. && fabs(pair->d2_mNSigmaKaon) < 3. && pair->d1_mPt > 0.01 && pair->d2_mPt > 0.01 ){ // event and track selection will go here


//      only select += pair
        if(pair->mChargeSum == 0.){
            // Setup the Lorentz Vectors of the daugher tracks
            lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.493 ); // we use Kaon mass here.
            lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.493 );
            // compute parent particle lorentz vector from daughters

            // use conservation of momentum and energy to set the four-momentum
            // vector for parent particle.
            lv = lv1 + lv2;


            track1_sigma_KP -> Fill(pair->d1_mNSigmaKaon, pair->d1_mNSigmaPion);
            track2_sigma_KP -> Fill(pair->d2_mNSigmaKaon, pair->d2_mNSigmaPion);

            if(fabs(fabs(pair->d1_mNSigmaKaon)) < 5 && fabs(pair->d2_mNSigmaKaon) < 5 &&
                    lv.Pt()< 0.02 &&
                    fabs(pair->d1_mNSigmaPion) >5 && fabs(pair->d2_mNSigmaPion) >5 &&
                    fabs(pair->d1_mNSigmaProton) > 5 && fabs(pair->d2_mNSigmaProton) > 5
                    ) 
            {

                mass-> Fill(lv.M());
            }

            // We apply the selection rule to select Kaon, and fill the phi_pt values for phi mass += 15MeV 
            if(fabs(pair->d1_mNSigmaKaon) < 5 && fabs(pair->d2_mNSigmaKaon) < 5 &&
                    lv.M()>= 1. && lv.M()<= 1.035 &&
                    pair->d1_mPt > 0.01 && pair->d2_mPt > 0.01 && 
                    fabs(pair->d1_mNSigmaPion) > 5 && fabs(pair->d2_mNSigmaPion > 5)&&
                    fabs(pair->d1_mNSigmaProton) > 5 && fabs(pair->d2_mNSigmaProton)>5){

                phi_pt-> Fill(lv.Pt());
                pt_mass-> Fill(lv.Pt(),lv.M());
            }

        } // selection
    } // loop on events

    fo->cd();

    makeCan();

    // pt_mass-> Draw("colz");
    
    //track1_sigma_KP-> Draw("colz");
    //track2_sigma_KP-> Draw("colz");

    //phi_pt->Draw();


    mass->Draw();
    // pt_mass->Draw("colz");
    // write all histograms to output data file
    fo->Write();
  }             
