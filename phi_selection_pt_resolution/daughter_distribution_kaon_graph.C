// A simple TTreeReader use: read data from hsimple.root (written by hsimple.C)
#include "FemtoPairFormat.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPaveStats.h"
#include "TLorentzVector.h"
#include "TPaveStats.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <Rtypes.h>
#include <TH2.h>
#include <TSystem.h>
#include <TVirtualPad.h>
#include <cmath>
#include <complex>
#include <math.h>

const double NSigmaPionShift = 0.8;
const double NSigmaElectronShift = 0.424112;
const double NSigmaProtonShift = 2.36;
const double NSigmaKaonShift = 2.4;

// ClassImp( FemtoPair );

int ican = 0;
void makeCan(bool logz = false) {
  TCanvas *can = new TCanvas(TString::Format("can%d", ican++), "", 900, 600);
  can->SetTopMargin(0.10);
  can->SetRightMargin(0.03);

  if (logz) {
    can->SetLogz(1); // Set Z-axis to logarithmic scale
  }
}

void daughter_distribution_kaon_graph() {

  // histograms for all tracks
  TH1F *all_DCA = new TH1F(
      "All Tracks", "Run 19 A+A DCA All Tracks;DCA(cm);counts", 100, 0., 3.);

  TH1F *all_plus_DCA = new TH1F(
      "Positive Tracks", "Run 19 A+A DCA All Tracks;DCA(cm);counts", 100, 0., 3.);

  TH1F *all_minus_DCA = new TH1F(
      "Negative Tracks", "Run 19 A+A DCA All Tracks;DCA(cm);counts", 100, 0., 3.);

  TH1F *all_pt =
      new TH1F("All Tracks", "Run 19 A+A P_{T} All Tracks;P_{T}(GeV/c);counts",
               100, 0., 1.8);

  TH1F *all_plus_pt =
      new TH1F("Positive Tracks", "Run 19 A+A P_{T} All Tracks;P_{T}(GeV/c);counts",
               100, 0., 1.8);

  TH1F *all_minus_pt =
      new TH1F("Negative Tracks", "Run 19 A+A P_{T} All Tracks;P_{T}(GeV/c);counts",
               100, 0., 1.8);

  TH1F *all_eta = new TH1F(
      "All Tracks", "Run 19 A+A #eta All Tracks;#eta;counts", 100, -2.2, 2.2);

  TH1F *all_plus_eta = new TH1F(
      "Positive Tracks", "Run 19 A+A #eta All Tracks;#eta;counts", 100, -2.2, 2.2);

  TH1F *all_minus_eta = new TH1F(
      "Negative Tracks", "Run 19 A+A #eta All Tracks;#eta;counts", 100, -2.2, 2.2);

  TH1F *all_phi =
      new TH1F("All Tracks", "Run 19 A+A #phi All Tracks;#phi(rad);counts", 200,
               -6.3, 6.3);
  TH1F *all_plus_phi =
      new TH1F("Positive Tracks", "Run 19 A+A #phi All Tracks;#phi(rad);counts",
               200, -6.3, 6.3);
  TH1F *all_minus_phi =
      new TH1F("Negative Tracks", "Run 19 A+A #phi All Tracks;#phi(rad);counts",
               200, -6.3, 6.3);

  TH2F *plus_phi_eta =
      new TH2F("Positive Track",
               "Run 19 A+A #phi vs. #eta Positive Tracks;#eta;#phi(rad)", 100,
               -2., 2., 100, -3.15, 3.15);

  TH2F *minus_phi_eta =
      new TH2F("Negative Track",
               "Run 19 A+A #phi vs. #eta Negative Tracks;#eta;#phi(rad)", 100,
               -2., 2., 100, -3.15, 3.15);
  // Open the file containing the tree (INPUT data).
  TFile *myFile = TFile::Open("input.root");

  // This setup the reader, access the data
  TTreeReader myReader("PairDst", myFile);

  // We need this so ROOT understands the data format
  gSystem->Load("FemtoPairFormat_h.so");
  TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

  // Lorentz vectors (4-vectors) for kinematics in special relativity
  TLorentzVector lv1, lv2, lv;

  double num_all_track = 0;
  double num_chargesum = 0;
  double num_kaon_accept = 0;
  double num_pion_reject = 0;
  double num_electron_reject = 0;
  double num_proton_reject = 0;
  double num_mass_cut = 0;
  double num_pt_cut = 0;
  // Loop over all entries of the TTree or TChain.
  while (myReader.Next()) {
    num_all_track++;
    const double NSigmaPionShift = 0.8;
    const double NSigmaElectronShift = 0.424112;
    const double NSigmaProtonShift = 2.36;
    const double NSigmaKaonShift = 2.4;

    bool kaon_accept = abs(pair->d1_mNSigmaKaon - NSigmaKaonShift) < 5 &&
                       abs(pair->d2_mNSigmaKaon - NSigmaKaonShift) < 5;

    bool pion_reject = abs(pair->d1_mNSigmaPion - NSigmaPionShift) > 5 &&
                       abs(pair->d2_mNSigmaPion - NSigmaPionShift) > 5;

    bool electron_reject =
        abs(pair->d1_mNSigmaElectron - NSigmaElectronShift) > 5 &&
        abs(pair->d2_mNSigmaElectron - NSigmaElectronShift) > 5;

    bool proton_reject = abs(pair->d1_mNSigmaProton - NSigmaProtonShift) > 5 &&
                         abs(pair->d2_mNSigmaProton - NSigmaProtonShift) > 5;

    bool mass_cut = lv.M() >= 1.0 && lv.M() <= 1.04;

    bool pt_cut = lv.Pt() <= 0.1;

    if (pair->mChargeSum == 0) {
      num_chargesum++;
      if (kaon_accept) {
        num_kaon_accept++;
        if (pion_reject) {
          num_pion_reject++;
          if (electron_reject) {
            num_electron_reject++;
            if (proton_reject) {
              num_proton_reject ++;
              if (mass_cut){
                num_mass_cut ++;
                if (pt_cut){
                  num_pt_cut ++;
                }
              }
            }
          }
        }
      }
    }

      lv1.SetPtEtaPhiM(pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.493);
      lv2.SetPtEtaPhiM(pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.493);
      lv = lv1 + lv2;

      // all track histogram filling
      all_pt->Fill(pair->d1_mPt);
      all_pt->Fill(pair->d2_mPt);

      all_DCA->Fill(pair->d1_mDCA);
      all_DCA->Fill(pair->d2_mDCA);

      all_eta->Fill(pair->d1_mEta);
      all_eta->Fill(pair->d2_mEta);

      all_phi->Fill(pair->d1_mPhi);
      all_phi->Fill(pair->d2_mPhi);

      if (pair->mChargeSum == 0) {
        all_plus_pt->Fill(pair->d1_mPt);
        all_minus_pt->Fill(pair->d2_mPt);

        all_plus_eta->Fill(pair->d1_mEta);
        all_minus_eta->Fill(pair->d2_mEta);

        all_plus_phi->Fill(pair->d1_mPhi);
        all_minus_phi->Fill(pair->d2_mPhi);

        all_plus_DCA->Fill(pair->d1_mDCA);
        all_minus_DCA->Fill(pair->d2_mDCA);

        plus_phi_eta->Fill(pair->d1_mEta, pair->d1_mPhi);
        minus_phi_eta->Fill(pair->d2_mEta, pair->d2_mPhi);
      } else if (pair->mChargeSum == 2) {
        all_plus_pt->Fill(pair->d1_mPt);
        all_plus_pt->Fill(pair->d2_mPt);
        all_plus_eta->Fill(pair->d1_mEta);
        all_plus_eta->Fill(pair->d2_mEta);
        all_plus_phi->Fill(pair->d1_mPhi);
        all_plus_phi->Fill(pair->d2_mPhi);
        all_plus_DCA->Fill(pair->d1_mDCA);
        all_plus_DCA->Fill(pair->d2_mDCA);
        plus_phi_eta->Fill(pair->d1_mEta, pair->d1_mPhi);
        plus_phi_eta->Fill(pair->d2_mEta, pair->d2_mPhi);
      } else if (pair->mChargeSum == -2) {
        all_minus_pt->Fill(pair->d1_mPt);
        all_minus_pt->Fill(pair->d2_mPt);
        all_minus_eta->Fill(pair->d1_mEta);
        all_minus_eta->Fill(pair->d2_mEta);
        all_minus_phi->Fill(pair->d1_mPhi);
        all_minus_phi->Fill(pair->d2_mPhi);
        all_minus_DCA->Fill(pair->d1_mDCA);
        all_minus_DCA->Fill(pair->d2_mDCA);
        minus_phi_eta->Fill(pair->d1_mEta, pair->d1_mPhi);
        minus_phi_eta->Fill(pair->d2_mEta, pair->d2_mPhi);
      }

    } // loop on events

    std::cout << num_all_track << "    " << num_chargesum << "    " <<num_kaon_accept << "    "
              << num_pion_reject << "    " << num_electron_reject << "    "
              << num_proton_reject << "    " << num_mass_cut << "    " << num_pt_cut << std::endl;
    makeCan();
    all_plus_DCA->GetXaxis()->SetTitleSize(0.05);
    all_plus_DCA->GetXaxis()->CenterTitle();
    all_plus_DCA->GetXaxis()->SetTitleOffset(0.8);
    // all_plus_DCA->SetMarkerStyle(20);
    // all_plus_DCA->SetMarkerSize(0.5);
    // all_minus_DCA->SetMarkerStyle(20);
    // all_minus_DCA->SetMarkerSize(0.5);
    all_plus_DCA->SetFillColorAlpha(kRed, 0.3);
    all_minus_DCA->SetFillColorAlpha(kBlue, 0.3);
    all_DCA->SetFillColorAlpha(kGreen, 0.8);
    all_plus_DCA->Draw("hist");
    all_plus_DCA->Draw("e1;same");
    all_minus_DCA->Draw("sames;hist");
    all_minus_DCA->Draw("e1;same");

    TLegend *legend_DCA = new TLegend(0.78, 0.55, 0.98, 0.75);
    legend_DCA->AddEntry(all_plus_DCA, "Positive Track DCA", "f");
    legend_DCA->AddEntry(all_minus_DCA, "Negative Track DCA", "f");
    legend_DCA->Draw("same");
    // gPad->Print("./Plots_sigma/All/LFS/all_DCA.png");

    makeCan();
    all_plus_pt->GetXaxis()->SetTitleSize(0.05);
    all_plus_pt->GetXaxis()->CenterTitle();
    all_plus_pt->GetXaxis()->SetTitleOffset(0.8);
    // all_plus_pt->SetMarkerStyle(20);
    // all_plus_pt->SetMarkerSize(0.5);
    // all_minus_pt->SetMarkerStyle(20);
    // all_minus_pt->SetMarkerSize(0.5);
    all_plus_pt->SetFillColorAlpha(kRed, 0.3);
    all_minus_pt->SetFillColorAlpha(kBlue, 0.3);
    all_pt->SetFillColorAlpha(kGreen, 0.8);
    all_plus_pt->Draw("hist");
    all_plus_pt->Draw("e1;same");
    all_minus_pt->Draw("sames;hist");
    all_minus_pt->Draw("e1;same");

    TLegend *legend_pt = new TLegend(0.78, 0.55, 0.98, 0.75);
    legend_pt->AddEntry(all_plus_pt, "Positive Track P_{T}", "f");
    legend_pt->AddEntry(all_minus_pt, "Negative Track P_{T}", "f");
    legend_pt->Draw("same");
    // gPad->Print("./Plots_sigma/All/LFS/all_pt.png");

    makeCan();
    all_plus_eta->GetXaxis()->SetTitleSize(0.05);
    all_plus_eta->GetXaxis()->CenterTitle();
    all_plus_eta->GetXaxis()->SetTitleOffset(0.8);
    // all_plus_eta->SetMarkerStyle(20);
    // all_plus_eta->SetMarkerSize(0.5);
    // all_minus_eta->SetMarkerStyle(20);
    // all_minus_eta->SetMarkerSize(0.5);
    all_plus_eta->SetFillColorAlpha(kRed, 0.3);
    all_minus_eta->SetFillColorAlpha(kBlue, 0.3);
    all_eta->SetFillColorAlpha(kGreen, 0.8);
    all_plus_eta->Draw("hist");
    all_plus_eta->Draw("e1;same");
    all_minus_eta->Draw("sames;hist");
    all_minus_eta->Draw("e1;same");

    TLegend *legend_eta = new TLegend(0.78, 0.55, 0.98, 0.75);
    legend_eta->AddEntry(all_plus_eta, "Positive Track #eta", "f");
    legend_eta->AddEntry(all_minus_eta, "Negative Track #eta", "f");
    legend_eta->Draw("same");
    // gPad->Print("./Plots_sigma/All/LFS/all_eta.png");

    makeCan();
    all_plus_phi->GetXaxis()->SetTitleSize(0.05);
    all_plus_phi->GetXaxis()->CenterTitle();
    all_plus_phi->GetXaxis()->SetTitleOffset(0.8);
    // all_plus_phi->SetMarkerStyle(20);
    // all_plus_phi->SetMarkerSize(0.5);
    // all_minus_phi->SetMarkerStyle(20);
    // all_minus_phi->SetMarkerSize(0.5);
    all_plus_phi->SetFillColorAlpha(kRed, 0.3);
    all_minus_phi->SetFillColorAlpha(kBlue, 0.3);
    all_phi->SetFillColorAlpha(kGreen, 0.8);
    all_plus_phi->Draw("hist");
    all_plus_phi->Draw("e1;same");
    all_minus_phi->Draw("sames;hist");
    all_minus_phi->Draw("e1;same");

    TLegend *legend_phi = new TLegend(0.78, 0.55, 0.98, 0.75);
    legend_phi->AddEntry(all_plus_phi, "Positive Track #phi", "f");
    legend_phi->AddEntry(all_minus_phi, "Negative Track #phi", "f");
    legend_phi->Draw("same");
    // gPad->Print("./Plots_sigma/All/LFS/all_phi.png");
    //
    makeCan();
    plus_phi_eta->SetContour(100);
    plus_phi_eta->GetXaxis()->SetTitleSize(0.05);
    plus_phi_eta->GetXaxis()->CenterTitle();
    plus_phi_eta->GetXaxis()->SetTitleOffset(0.8);
    plus_phi_eta->GetYaxis()->SetTitleSize(0.05);
    plus_phi_eta->GetYaxis()->CenterTitle();
    plus_phi_eta->GetYaxis()->SetTitleOffset(0.8);
    plus_phi_eta->Draw("colz");
    gPad->Print("./Plots_sigma/All/LFS/plus_phi_eta.png");

    makeCan();
    minus_phi_eta->SetContour(100);
    minus_phi_eta->GetXaxis()->SetTitleSize(0.05);
    minus_phi_eta->GetXaxis()->CenterTitle();
    minus_phi_eta->GetXaxis()->SetTitleOffset(0.8);
    minus_phi_eta->GetYaxis()->SetTitleSize(0.05);
    minus_phi_eta->GetYaxis()->CenterTitle();
    minus_phi_eta->GetYaxis()->SetTitleOffset(0.8);
    minus_phi_eta->Draw("colz");
    gPad->Print("./Plots_sigma/All/LFS/minus_phi_eta.png");
  }
