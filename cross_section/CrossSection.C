#include "Rtypes.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TLegend.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
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
    vector<TF1 *> *efficiencyFunc = nullptr,
    vector<TEfficiency *> *tEfficiencies = nullptr) {

  TH1F *outputHist = (TH1F *)distributionHist->Clone();

  for (int i = 1; i < outputHist->GetNbinsX(); i++) {
    Double_t binCenter = distributionHist->GetBinCenter(i);
    Double_t totalEfficiency = 1.0;

    // Multiply with efficiencies from TEfficiency objects
    if (tEfficiencies) {
      for (TEfficiency *Efficiencies : *tEfficiencies) {
        int correspondingBin =
            static_cast<TH1F *>(
                const_cast<TH1 *>(Efficiencies->GetTotalHistogram()))
                ->FindBin(binCenter);
        totalEfficiency *= Efficiencies->GetEfficiency(correspondingBin);
      }
    }

    // Multiply with efficiencies from TH1F histograms
    if (efficiencyHist) {
      for (TH1F *effHist : *efficiencyHist) {
        int correspondingBin = effHist->FindBin(binCenter);
        totalEfficiency *= effHist->GetBinContent(correspondingBin);
      }
    }

    // Multiply with efficiencies from TF1 functions
    if (efficiencyFunc) {
      for (TF1 *effFunc : *efficiencyFunc) {
        totalEfficiency *= effFunc->Eval(binCenter);
      }
    }
    Double_t differentialCrossSection =
        (outputHist->GetBinContent(i)) /
        (outputHist->GetBinWidth(i) * totalLuminosity * totalEfficiency);

    outputHist->SetBinContent(i, differentialCrossSection);
  }

  return outputHist;
}

void CrossSection() {
  TFile *inputFile = ReadFile("./kaonDistribution.root");
  TFile *iTPC_effFile = ReadFile("./itpc_eff.root");
  TFile *PidEffFile = ReadFile("./PidEff.root");

  TH1F *ptDistr = ReadHistogram(inputFile, "reco_phi_pt");
  TH1F *massDistr = ReadHistogram(inputFile, "reco_phi_mass");
  TH1F *etaDistr = ReadHistogram(inputFile, "reco_phi_eta");

  vector<TEfficiency *> tEfficiencies;
  TEfficiency* massPidTEfficiency = (TEfficiency*) PidEffFile->Get("eff_pt_clone");
  tEfficiencies.push_back(massPidTEfficiency);

  vector<TH1F *> efficiencyHist;
  TH1F *iTPC_effHist = ReadHistogram(iTPC_effFile, "Reco Effiency");
  efficiencyHist.push_back(iTPC_effHist);

  vector<TF1 *> efficiencyFunctions;
  TF1 *fL = new TF1("constFunc", "[0]", -10, 10);
  fL->SetParameter(0,0.3);

  TF1 *mutEff = new TF1("constFunc", "[0]", -10, 10);
  mutEff->SetParameter(0,0.9);

  TF1 *dcaEff = new TF1("constFunc", "[0]", -10, 10);
  dcaEff->SetParameter(0,0.95);

  efficiencyFunctions.push_back(fL);
  efficiencyFunctions.push_back(mutEff);
  efficiencyFunctions.push_back(dcaEff);

  Double_t numEvent = ptDistr->GetEntries();
  Double_t totalLuminosity = 13.267;


  TH1F *dSigmaDPt = GetDifferentialCrossSectionFromDistribution(
      ptDistr, totalLuminosity, &efficiencyHist,&efficiencyFunctions,&tEfficiencies);
  dSigmaDPt->SetTitle("d#sigma/dP_{T}");
  dSigmaDPt->GetXaxis()->SetTitle("P_{T} (GeV)");
  dSigmaDPt->GetYaxis()->SetTitle("d#sigma/dP_{T} (GeV^{-1} #mu b)");
  dSigmaDPt->SetFillColor(kBlue);
  dSigmaDPt->Draw();
  std::cout<<dSigmaDPt->Integral()<<std::endl;

  TFile * slightCoherentFile = ReadFile("./slight/kaonCoherent.root");
  TH1F * slightSigmaPtCoherent =(TH1F*) slightCoherentFile->Get("PtKa");
  double slightCoherentTotalSigma = 1524;
  slightSigmaPtCoherent->Scale(slightCoherentTotalSigma/slightSigmaPtCoherent->Integral());
  std::cout<<slightSigmaPtCoherent->Integral() << std::endl;

  TFile * slightIncoherentFile = ReadFile("./slight/kaonIncoherent.root");
  TH1F * slightSigmaPtIncoherent =(TH1F*) slightIncoherentFile->Get("PtKa");
  double slightIncoherentTotalSigma = 931;
  slightSigmaPtIncoherent->Scale(slightIncoherentTotalSigma/slightSigmaPtIncoherent->Integral());
  std::cout<<slightSigmaPtIncoherent->Integral() << std::endl;

  TH1F * slightSigmaPtAll  = (TH1F*) slightSigmaPtCoherent->Clone();
  slightSigmaPtAll->Add(slightSigmaPtIncoherent);
  slightSigmaPtAll->SetTitle("d#sigma/dP_{T} Slight Coherent + Incoherent");
  slightSigmaPtAll->GetXaxis()->SetTitle("P_{T} (GeV)");
  slightSigmaPtAll->GetYaxis()->SetTitle("d#sigma/dP_{T} (GeV^{-1} #mu b)");
  slightSigmaPtAll->SetFillColor(kRed);
  slightSigmaPtAll->Draw("hist;same");

  TLegend * legend = new TLegend();
  legend->AddEntry(slightSigmaPtAll,"Slight Coherent + Incoherent");
  legend->AddEntry(dSigmaDPt,"Run 19 Au Au");
  legend->Draw("same");

}
