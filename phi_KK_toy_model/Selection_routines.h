#include "Random_routines.h"
#include <TFile.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <string>

#ifndef SELECTOR
#define SELECTOR

// This is an rudimentery attempt in coding a particle selector
// that tells the user, with a given particle four momentum, the number of
// sigma away of that four momentum to either Pion or Kaon
// using dE/dx vs Pt with a known blur(sigma). The "true" dE/dx file
// must be provided as a root file for this code to function.
class Selector {

public:
  // constructor for the particle selector. This requires the absolute file path
  // for the "true" dE/dx root file and a fixed blur amount
  Selector(double sigma);

  // overloaded constructor for using default values.
  Selector();
  // destructor
  ~Selector(){};

  // draw the blurred version of dEdx for both particles
  void draw_dEdx();
  // path to true dE/dx root file
  // It is assumed that the true dEdx file is named "dEdx.root"
  // in the directory where the analysis code is compiled.

  // Root object for the dEdx file
  TFile *dEdx_file = new TFile("dEdx.root");

  // Root histogram for dEdxKaon using most probable value
  TH1D *dEdxKaon = (TH1D *)dEdx_file->Get("mpmK");
  // Root histogram for dEdxPion using most probable value
  TH1D *dEdxPion = (TH1D *)dEdx_file->Get("mpmPi");
  TH1D *dEdxElectron = (TH1D *)dEdx_file->Get("mpmE");

  // calculate the nsigma values
  double compute_NSigmaPion(double pt, TH1 *particle_dedx_distr);
  double compute_NSigmaKaon(double pt, TH1 *particle_dedx_distr);
  double compute_NSigmaElectron(double pt, TH1 *particle_dedx_distr);

  // draw dEdx for both particles
  void draw_dEdx_blurred();

private:
  // return the smeared dedx value for a given particle type;
  double sample_dedx(double pt, TH1 *particle_dedx_distr);

private:
  // fixed blur. This is obtained by adding a gaussian blur of dEdx using
  // the precise method of "eyeballing" to match the dEdx measurement graph
  double sigma = 0.00025;

};

#endif
