// include header files
#include "FemtoPairFormat.h"
#include "root/Math/Vector4D.h"
#include "root/TCanvas.h"
#include "root/TEfficiency.h"
#include "root/TFile.h"
#include "root/TGraphAsymmErrors.h"
#include "root/TH1F.h"
#include "root/TLorentzVector.h"
#include "root/TSystem.h"
#include "root/TTreeReader.h"
#include "root/TTreeReaderValue.h"
#include <root/Math/Vector4Dfwd.h>
#include <root/TH1.h>
#include <root/TLegend.h>
#include <root/TVirtualPad.h>
//
// fitted PID shift for run 19 Au-Au data
const double NSigmaPionShift = 0.8;
const double NSigmaElectronShift = 0.424112;
const double NSigmaProtonShift = 2.36;
const double NSigmaKaonShift = 2.4;

// kaon_mass
const double KAON_MASS = 0.493;
const double PION_MASS = 0.139;

// create a 720p canvas
int ican = 0;
void makeCan(bool logz = false) {
  TCanvas *can = new TCanvas(TString::Format("can%d", ican++), "", 1280, 720);
  can->SetTopMargin(0.10);
  can->SetRightMargin(0.03);
}


double total_eff;

TH1F *effPt = new TH1F("eff_pt", "All Probes Reco P_T", 50, 0., 2.);
TH1F *effM = new TH1F("eff_mass", "All Probes Reco M", 50, 1., 1.04);
TH1F *effY = new TH1F("eff_y", "All Probes Reco y", 50, -1, 1);

TH1F *recoPt = new TH1F("eff_pt", "All Probes Reco P_T", 50, 0., 2.);
TH1F *recoM = new TH1F("eff_mass", "All Probes Reco M", 50, 0.67, 0.87);
TH1F *recoY = new TH1F("eff_y", "All Probes Reco Y", 50, -1, 1);
TH1F *dEdxProbe = new TH1F("All Probes", "All Probes dE/dx", 10, -50, 30);
TH1F *recoPtProbeCut =
    new TH1F("Selected Probes", "Selected Probes Reco P_T", 50, 0., 2.);
TH1F *recoMProbeCut =
    new TH1F("Selected Probes", "Selected Probes Reco M", 50, 0.67, 0.87);
TH1F *recoYProbeCut =
    new TH1F("Selected Probes", "Selected Probes Reco Y", 50, -1, 1);
// routine used to calcualte the pid effiency
void PidEff() {
  // loading run 19 pairdst file
  TFile *myFile = TFile::Open("~/Documents/Research Data/Run_19_Au_Au.root");

  TFile *output_file = new TFile("PidEff.root", "Recreate");

  // This setup the reader, access the data
  TTreeReader myReader("PairDst", myFile);

  // We need this so ROOT understands the data format
  gSystem->Load("./FemtoPairFormat_h.so");
  TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

  // Lorentz vectors (4-vectors) for kinematics in special relativity
  TLorentzVector lv1, lv2, lv;
  // Loop over all entries of the TTree or TChain.
  while (myReader.Next()) {
    // load in pair data to LorentzVector
    lv1.SetPtEtaPhiM(pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, PION_MASS);
    lv2.SetPtEtaPhiM(pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, PION_MASS);
    lv = lv1 + lv2;

    // selection pid conditions
    bool inMassRange = (0.67 < lv.M()) && (lv.M() < 0.87);
    // tag condition on d1
    bool tagPionD1 =
        ((abs(pair->d1_mNSigmaKaon - NSigmaKaonShift) > 5.) &&
         (abs(pair->d1_mNSigmaPion - NSigmaPionShift) < 1.) &&
         (abs(pair->d1_mNSigmaProton - NSigmaProtonShift) > 5.) &&
         (abs(pair->d1_mNSigmaElectron - NSigmaElectronShift) > 5.));
    // probe selection on d2
    bool probePionD2 = abs(pair->d2_mNSigmaPion - NSigmaPionShift) < 5.;

    if (inMassRange && tagPionD1) {
      dEdxProbe->Fill(pair->d2_mNSigmaKaon);
      recoM->Fill(lv.M());
      recoPt->Fill(lv.Pt());
      recoY->Fill(lv.Rapidity());
    }
    if (inMassRange && tagPionD1 && probePionD2) {
      dEdxProbe->Fill(pair->d2_mNSigmaKaon);
      recoMProbeCut->Fill(lv.M());
      recoPtProbeCut->Fill(lv.Pt());
      recoYProbeCut->Fill(lv.Rapidity());
    }
  }
  makeCan();
  dEdxProbe->Draw();
  makeCan();
  TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend1->AddEntry(recoM);
  legend1->AddEntry(recoMProbeCut);
  recoM->SetFillColor(kBlue);
  recoMProbeCut->SetFillColor(kRed);
  recoM->Draw();
  recoMProbeCut->Draw("same");
  legend1->Draw("same");
  //
  TEfficiency *ptEff = new TEfficiency(*recoPtProbeCut, *recoPt);
  ptEff->SetStatisticOption(TEfficiency::kFCP);
  TEfficiency *yEff = new TEfficiency(*recoYProbeCut, *recoY);
  yEff->SetStatisticOption(TEfficiency::kFCP);
  TEfficiency *mEff = new TEfficiency(*recoMProbeCut, *recoM);
  mEff->SetStatisticOption(TEfficiency::kFCP);
  //
  makeCan();
  ptEff->SetTitle("PID Effiency per Momentum Bin");
  ptEff->Draw("ap");
  makeCan();
  yEff->SetTitle("PID Effiency per y Bin");
  yEff->Draw("ap");
  gPad->Update();
  makeCan();
  TGraphAsymmErrors *yEffGraph = yEff->GetPaintedGraph();
  yEffGraph->SetMinimum(-0.2);
  yEffGraph->Draw("ap");
  makeCan();
  mEff->SetTitle("PID Effiency per inv. Mass Bin");
  mEff->Draw("ap");
  //
  total_eff = recoYProbeCut->GetEntries() / recoY->GetEntries();
  //
  ptEff->Write();
  mEff->Write();
  yEff->Write();
  //
  output_file->Close();
  //
  std::cout << "total effiency is " << total_eff << std::endl;
}
