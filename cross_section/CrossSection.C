#include "TFile.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include <RtypesCore.h>
#include <TH1.h>
#include <iostream>
#include <vector>

using namespace std;

TFile *ReadFile(string fileName) {
  TFile *histogramFile = new TFile(fileName.c_str());
  return histogramFile;
}

TH1F *ReadHistogram(TFile *inFile, string histName) {
  TH1F *outputHist = (TH1F *)inFile->Get(histName.c_str());
  return outputHist;
}

TH1F *DivideHistograms(TH1F *h1, TH1F *h2) {
  if (!h1 || !h2) {
    std::cerr << "One or both histograms are null!" << std::endl;
    return nullptr;
  }

  // Check bin compatibility
  if (h1->GetNbinsX() != h2->GetNbinsX()) {
    std::cerr << "Histograms have different bin counts!" << std::endl;
    return nullptr;
  }

  // Clone the first histogram to keep it unchanged and to store the result
  TH1F *hResult = (TH1F *)h1->Clone("hResult");
  hResult->SetTitle("h1/h2");

  // Divide
  hResult->Divide(h2);

  return hResult;
}

Double_t CalculateTotalLuminosityFromCrossSection(Double_t numEvent,
                                                  Double_t totalCrossSection) {
  return numEvent / totalCrossSection;
}

TH1F *GetDifferentialCrossSectionFromDistribution(
    TH1F *distributionHist, Double_t totalLuminosity,
    vector<TH1F *> *efficiencyHist = nullptr,
    vector<TF1 *> *efficiencyFunc = nullptr) {
  TH1F *outputHist = (TH1F *)distributionHist->Clone();
  if (!efficiencyHist && !efficiencyFunc) {
    for (int i = 1; i <= outputHist->GetNbinsX(); i++) {
      Double_t differentialCrossSection =
          (outputHist->GetBinContent(i)) /
          (outputHist->GetBinWidth(i) * totalLuminosity);
      outputHist->SetBinContent(i, differentialCrossSection);
    }
  } else if (efficiencyHist) {
    for (int i = 1; i <= outputHist->GetNbinsX(); i++) {
      Double_t totalEffiency = 1;
      Double_t binCenter = distributionHist->GetBinCenter(i);
      for (TH1F *effHist : *efficiencyHist) {
        int correspondingBin = effHist->FindBin(binCenter);
        totalEffiency *= effHist->GetBinContent(correspondingBin);
      }
      Double_t differentialCrossSection =
          (outputHist->GetBinContent(i)) /
          (outputHist->GetBinWidth(i) * totalLuminosity * totalEffiency);
      outputHist->SetBinContent(i, differentialCrossSection);
    }
  }
  return outputHist;
}

void CrossSection() {
  TFile *slightFile = ReadFile("./histograms.root");
  TFile *iTPC_effFile = ReadFile("./itpc_eff.root");
  TH1F *ptDistr = ReadHistogram(slightFile, "PtKa");

  vector<TH1F *> efficiencyHist;
  TH1F *iTPC_effHist = ReadHistogram(iTPC_effFile, "Reco Effiency");
  efficiencyHist.push_back(iTPC_effHist);

  Double_t numEvent = ptDistr->GetEntries();
  Double_t totalCrossSection = 1.532;
  Double_t totalLuminosity =
      CalculateTotalLuminosityFromCrossSection(numEvent, totalCrossSection);

  TH1F *dSigmaDPt = GetDifferentialCrossSectionFromDistribution(
      ptDistr, totalLuminosity, &efficiencyHist);
  dSigmaDPt->SetTitle("d#sigma/dP_{T}");
  dSigmaDPt->GetXaxis()->SetTitle("P_{T} (GeV)");
  dSigmaDPt->GetYaxis()->SetTitle("d#sigma/dP_{T} (GeV^{-1} mb)");

  dSigmaDPt->Draw();
}
