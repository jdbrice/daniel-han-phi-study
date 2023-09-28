#ifndef ANALYZE_H
#define ANALYZE_H

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TH1F.h"
#include <string>

class Analyze 
{
 public:
  Analyze(); //Constructor
  Analyze(TString infile, Int_t nEvents); //Special constructor
  ~Analyze(); //Destructor
  Int_t Init();
  Int_t NextEvent();
  TParticle* NextParticle();
  void doAnalysis(); //Function doing the analysis
  
 private:
 
  TH1F *fPtEl;       //Histogram of pt of electrons
  TH1F *fRapEl;      //Histogram of rapidity of electrons
  TH1F *fInvMassEl;  //Histogram of ivariant mass of electrons
  TH1F *fAzimuthEl;  //Histogram of pt azimuthal angle of electrons

  TH1F *fPtMu;       //Histogram of pt of muons
  TH1F *fRapMu;      //Histogram of rapidity of muons
  TH1F *fInvMassMu;  //Histogram of ivariant mass of muons
  TH1F *fAzimuthMu;  //Histogram of pt azimuthal angle of muons

  TH1F *fPtPi;       //Histogram of pt of pions
  TH1F *fRapPi;      //Histogram of rapidity of pions
  TH1F *fInvMassPi;  //Histogram of ivariant mass of pions
  TH1F *fAzimuthPi;  //Histogram of pt azimuthal angle of pions

  TH1F *fPtKa;       //Histogram of pt of kaons
  TH1F *fPtKa_tpc;       //Histogram of pt of kaons
  TH1F *fPtKa_tpc_2;
  TH1F *fRapKa;      //Histogram of rapidity of kaons
  TH1F *fEtaKa;      //Histogram of rapidity of kaons
  TH1F *fEtaPi;      //Histogram of rapidity of kaons
  TH1F *fEtaEl;      //Histogram of rapidity of kaons
  TH1F *fInvMassKa;  //Histogram of ivariant mass of kaons
  TH1F *fAzimuthKa;  //Histogram of pt azimuthal angle of kaons

  TH1F *fPt1;
  TH1F *fPt2;

  TH1F *fRap1;
  TH1F *fRap2;

  TH1F *fEta1;
  TH1F *fEta2;

  TH1F *fPhi1;
  TH1F *fPhi2;
  
  FILE *filelist;
  TString fInfile;
  Int_t fNParticles;
  Int_t fNEvents;
  
};

#endif
