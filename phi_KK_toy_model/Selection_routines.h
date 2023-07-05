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
class Selector
{

public:
  // constructor for the particle selector. This requires the absolute file path
  // for the "true" dE/dx root file and a fixed blur amount
  Selector(double sigma);

  // overloaded constructor for using default values.
  Selector();
  // destructor
  ~Selector(){};

  // calculate NSigmaKaon with a parent vector and its topological
  // reconstruction. The target is the exact dEdx of Kaon
  // , while the "measurement" is the blurred dEdx of the assumed particle
  double get_NSigmaKaon(TLorentzVector *mc_ptr);

  // calculate NSigmaPion with a parent vector and its topological
  // reconstruction. The target is the exact dEdx of Pion
  // , while the "measurement" is the blurred dEdx of the assumed particle
  double get_NSigmaPion(TLorentzVector *mc_ptr);

  double get_NSigmaElectron(TLorentzVector *mc_ptr);

  // draw the blurred version of dEdx for both particles
  void draw_dEdx();

  // draw dEdx for both particles
  void draw_dEdx_blurred();

private:
  // fixed blur. This is obtained by adding a gaussian blur of dEdx using
  // the precise method of "eyeballing" to match the dEdx measurement graph
  double sigma = 0.00025;

  // path to true dE/dx root file
  // It is assumed that the true dEdx file is named "dEdx.root"
  // in the directory where the analysis code is compiled.
  std::string file_path = "dEdx.root";

  // Root object for the dEdx file
  TFile *dEdx_file = new TFile(file_path.c_str());

  // Root histogram for dEdxKaon using most probable value
  TH1D *dEdxKaon = (TH1D *)dEdx_file->Get("mpmK");
  // Root histogram for dEdxPion using most probable value
  TH1D *dEdxPion = (TH1D *)dEdx_file->Get("mpmPi");
  TH1D *dEdxElectron = (TH1D *)dEdx_file->Get("mpmE");
};

#endif
