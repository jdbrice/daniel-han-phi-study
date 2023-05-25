#include "Random_routines.h"
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
  Selector(std::string file_path, double sigma);

  // destructor
  ~Selector(){};

  // calculate NSigmaKaon with a parent vector and its topological
  // reconstruction.
  double get_NSigmaKaon(TLorentzVector *mc_ptr);

  // calculate NSigmaPion with a parent vector and its topological
  // reconstruction.
  double get_NSigmaPion(TLorentzVector *mc_ptr);

  // blur the dE/dx file with a fixed amount of percent blur
  void blur_dEdx_file();

private:
  // fixed percent blur
  double sigma;

  // path to true dE/dx root file
  std::string file_path;


};

#endif
