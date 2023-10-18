// include header files
#include "FemtoPairFormat.h"
#include "Math/Vector4D.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TLegend.h>
#include <TVirtualPad.h>

// fitted PID shift for run 19 Au-Au data
const double NSigmaPionShift = 0.8;
const double NSigmaElectronShift = 0.424112;
const double NSigmaProtonShift = 2.36;
const double NSigmaKaonShift = 2.4;

// kaon_mass
const double KAON_MASS = 0.493;

// create a 720p canvas
int ican = 0;
void makeCan(bool logz = false) {
  TCanvas *can = new TCanvas(TString::Format("can%d", ican++), "", 1280, 720);
  can->SetTopMargin(0.10);
  can->SetRightMargin(0.03);
}

double total_eff;

TH1F *effPt =
    new TH1F("Selected/All probes", "All Probes Reco P_T", 10, 0., 2.);
TH1F *effM = new TH1F("Selected/All Probes", "All Probes Reco M", 10, 1., 1.04);
TH1F *effEta =
    new TH1F("Selected/All Probes", "All Probes Reco Eta", 10, -2, 2);

TH1F *recoPt = new TH1F("All probes", "All Probes Reco P_T", 10, 0., 2.);
TH1F *recoM = new TH1F("All Probes", "All Probes Reco M", 10, 1., 1.04);
TH1F *recoEta = new TH1F("All Probes", "All Probes Reco Eta", 10, -2, 2);
TH1F *dEdxProbe = new TH1F("All Probes", "All Probes dE/dx", 10, -50, 30);
TH1F *recoPtProbeCut =
    new TH1F("Selected Probes", "Selected Probes Reco P_T", 10, 0., 2.);
TH1F *recoMProbeCut =
    new TH1F("Selected Probes", "Selected Probes Reco M", 10, 1., 1.04);
TH1F *recoEtaProbeCut =
    new TH1F("Selected Probes", "Selected Probes Reco Eta", 10, -2, 2);
// routine used to calcualte the pid effiency
void PidEff() {
  // loading run 19 pairdst file
  TFile *myFile = TFile::Open("~/Documents/Research Data/Run_19_Au_Au.root");

  // This setup the reader, access the data
  TTreeReader myReader("PairDst", myFile);

  // We need this so ROOT understands the data format
  gSystem->Load("FemtoPairFormat_h.so");
  TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

  // Lorentz vectors (4-vectors) for kinematics in special relativity
  ROOT::Math::PtEtaPhiMVector lv1, lv2, lv;
  // Loop over all entries of the TTree or TChain.
  while (myReader.Next()) {
    // load in pair data to LorentzVector
    lv1.SetCoordinates(pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, KAON_MASS);
    lv2.SetCoordinates(pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, KAON_MASS);
    lv = lv1 + lv2;

    // selection pid conditions
    bool inMassRange = (1.0 < lv.M()) && (lv.M() < 1.04);
    bool inPtrange = (lv.Pt() < 1.0);
    // tag condition on d1
    bool tagKaonD1 =
        ((abs(pair->d1_mNSigmaKaon - NSigmaKaonShift) < 1) &&
         (abs(pair->d1_mNSigmaPion - NSigmaPionShift) > 5.) &&
         (abs(pair->d1_mNSigmaProton - NSigmaProtonShift) > 5.) &&
         (abs(pair->d1_mNSigmaElectron - NSigmaElectronShift) > 5.));
    // probe selection on d2
    bool probeKaonD2 = abs(pair->d2_mNSigmaKaon - NSigmaKaonShift) < 5.;

    if (inMassRange && tagKaonD1) {
      dEdxProbe->Fill(pair->d2_mNSigmaKaon);
      recoM->Fill(lv.M());
      recoPt->Fill(lv.Pt());
      recoEta->Fill(lv.Eta());
    }
    if (inMassRange && tagKaonD1 && probeKaonD2) {
      dEdxProbe->Fill(pair->d2_mNSigmaKaon);
      recoMProbeCut->Fill(lv.M());
      recoPtProbeCut->Fill(lv.Pt());
      recoEtaProbeCut->Fill(lv.Eta());
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
  makeCan();
  TLegend *legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend2->AddEntry(recoPt);
  legend2->AddEntry(recoPtProbeCut);
  recoPt->SetFillColor(kBlue);
  recoPtProbeCut->SetFillColor(kRed);
  recoPt->Draw();
  recoPtProbeCut->Draw("same");
  legend2->Draw("same");
  makeCan();
  TLegend *legend3 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend3->AddEntry(recoEta);
  legend3->AddEntry(recoEtaProbeCut);
  recoEta->SetFillColor(kBlue);
  recoEtaProbeCut->SetFillColor(kRed);
  recoEta->Draw();
  recoEtaProbeCut->Draw("same");
  legend3->Draw("same");

  TEfficiency *ptEff = new TEfficiency(*recoPtProbeCut, *recoPt);
  ptEff->SetStatisticOption(TEfficiency::kFCP);
  TEfficiency *etaEff = new TEfficiency(*recoEtaProbeCut, *recoEta);
  etaEff->SetStatisticOption(TEfficiency::kFCP);
  TEfficiency *mEff = new TEfficiency(*recoMProbeCut, *recoM);
  mEff->SetStatisticOption(TEfficiency::kFCP);

  makeCan();
  ptEff->SetTitle("PID Effiency per Momentum Bin");
  ptEff->Draw("ap");
  gPad->Update();
  makeCan();
  TGraphAsymmErrors* ptEffGraph = ptEff->GetPaintedGraph();
  ptEffGraph->SetMinimum(-0.2);
  ptEffGraph->Draw("ap");
  makeCan();
  etaEff->SetTitle("PID Effiency per #eta Bin");
  etaEff->Draw("ap");
  makeCan();
  mEff->SetTitle("PID Effiency per inv. Mass Bin");
  mEff->Draw("ap");

  total_eff = recoEtaProbeCut->GetEntries() / recoEta->GetEntries();

  std::cout << "total effiency is " << total_eff << std::endl;
}
