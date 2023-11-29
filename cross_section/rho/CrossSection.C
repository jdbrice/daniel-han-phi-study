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
#include <fstream>
#include <iomanip>

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

  std::ofstream outFile("out.dat");
  outFile << std::left;
  outFile << std::setw(10) << "Count"
    << std::setw(20) << "Luminosity(mb-1)"
    << std::setw(20) << "bin Width"
    << std::setw(20) << "Pid Eff"
    << std::setw(20) << "iTPC cutoff Eff"
    << std::setw(20) << "UPC fraction"
    << std::setw(20) << "DCA"
    << std::setw(20) << "Track Reco" 
    << std::setw(20) << "Differential Cross Section (mb)" << std::endl;

  TH1F *outputHist = (TH1F *)distributionHist->Clone();

  for (int i = 1; i < outputHist->GetNbinsX()+1; i++) {
    Double_t binCenter = distributionHist->GetBinCenter(i);
    Double_t totalEfficiency = 1.0;

    outFile << std::setw(20) << outputHist->GetBinContent(i);
    outFile << std::setw(20) << 13.267 * 1000.;
    outFile << std::setw(20) << outputHist->GetBinWidth(i);

    // Multiply with efficiencies from TEfficiency objects
    if (tEfficiencies) {
      for (TEfficiency *Efficiencies : *tEfficiencies) {
        int correspondingBin =
            static_cast<TH1F *>(
                const_cast<TH1 *>(Efficiencies->GetTotalHistogram()))
                ->FindBin(binCenter);
        totalEfficiency *= Efficiencies->GetEfficiency(correspondingBin);
        outFile << std::setw(10)<<Efficiencies->GetEfficiency(correspondingBin);
      }
    }

    // Multiply with efficiencies from TH1F histograms
    if (efficiencyHist) {
      for (TH1F *effHist : *efficiencyHist) {
        int correspondingBin = effHist->FindBin(binCenter);
        totalEfficiency *= effHist->GetBinContent(correspondingBin);
        outFile << std::setw(20)<<effHist->GetBinContent(correspondingBin);
      }
    }

    // Multiply with efficiencies from TF1 functions
    if (efficiencyFunc) {
      for (TF1 *effFunc : *efficiencyFunc) {
        totalEfficiency *= effFunc->Eval(binCenter);
        outFile << std::setw(20) << effFunc->Eval(binCenter);
      }
    }
    Double_t differentialCrossSection =
        (outputHist->GetBinContent(i)) /
        (outputHist->GetBinWidth(i) * totalLuminosity * totalEfficiency);
    outFile << setw(20) << differentialCrossSection << std::endl;

    outputHist->SetBinContent(i, differentialCrossSection);
  }


  return outputHist;
}

void CrossSection() {
  TFile *inputFile = ReadFile("./pionDistribution.root");
  TFile *PidEffFile = ReadFile("./PidEff.root");
  TFile *iTPC_effFile = ReadFile("./itpc_eff.root");

  TH1F *yDistr = ReadHistogram(inputFile, "reco_rho_y");

  vector<TEfficiency *> tEfficiencies;
  TEfficiency* yPidTEfficiency = (TEfficiency*) PidEffFile->Get("eff_y_clone");
  tEfficiencies.push_back(yPidTEfficiency);

  vector<TH1F *> efficiencyHist;
  TH1F *iTPC_effHist = ReadHistogram(iTPC_effFile, "Reco Effiency");
  efficiencyHist.push_back(iTPC_effHist);

  vector<TF1 *> efficiencyFunctions;
  TF1 *fL = new TF1("constUPCFraction", "[0]", -10, 10);
  fL->SetParameter(0,0.3);

  TF1 *mutEff = new TF1("constTrackRecoEff", "[0]", -10, 10);
  mutEff->SetParameter(0,0.9);

  TF1 *dcaEff = new TF1("constDCAEff", "[0]", -10, 10);
  dcaEff->SetParameter(0,0.95);

  efficiencyFunctions.push_back(fL);
  efficiencyFunctions.push_back(mutEff);
  efficiencyFunctions.push_back(dcaEff);

  Double_t numEvent = yDistr->GetEntries();
  Double_t totalLuminosity = 13.267*1000;

  TH1F *dSigmaDy = GetDifferentialCrossSectionFromDistribution(
      yDistr, totalLuminosity, &efficiencyHist,&efficiencyFunctions,&tEfficiencies);
  dSigmaDy->SetTitle("d#sigma/dy");
  dSigmaDy->GetXaxis()->SetTitle("y");
  dSigmaDy->GetYaxis()->SetTitle("d#sigma/dy ((#Delta y)^{-1} mb)");
  dSigmaDy->SetFillColor(kBlue);
  dSigmaDy->Draw();
  std::cout<<dSigmaDy->Integral("width")<<std::endl;

  TFile * slightCoherentFile = ReadFile("./rho_coherent.root");
  TH1F * slightSigmaRapCoherent =(TH1F*) slightCoherentFile->Get("RapPi");
  double slightCoherentTotalSigma = 39.964;
  slightSigmaRapCoherent->Scale(slightCoherentTotalSigma/slightSigmaRapCoherent->Integral(),"width");
  std::cout<<slightSigmaRapCoherent->Integral("width") << std::endl;
  slightSigmaRapCoherent->SetFillColor(kRed);
  slightSigmaRapCoherent->Draw("same;hist");

  TLegend * legend = new TLegend();
  legend->AddEntry(slightSigmaRapCoherent,"Slight Coherent");
  legend->AddEntry(dSigmaDy,"Run 19 Au Au");
  legend->Draw("same");
}
