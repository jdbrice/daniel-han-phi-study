#include "particle_routines.h"

namespace ParticleRoutines {

std::vector<std::shared_ptr<Particle>>
twoBodyDecay(std::shared_ptr<Particle> parent, double d1Mass, double d2Mass) {
  // TODO: Implement the function
  return {}; // Placeholder return
}

std::shared_ptr<Particle> generateParticlePtEtaPhiM(int pdgCode, double pt,
                                                    double eta, double phi,
                                                    double M) {
  // TODO: Implement the function
  return nullptr; // Placeholder return
}

std::shared_ptr<Particle> generateParticleSlight(int pdgCode,
                                                 TFile *slightFile) {
  // TODO: Implement the function
  return nullptr; // Placeholder return
                  //
                  //
}

double calculateNsigmaKaon(std::shared_ptr<Particle> target, TFile *dEdxFile) {
  // TODO: Implement the function
  return 0.0; // Placeholder return
}

double calculateNsigmaPion(std::shared_ptr<Particle> target, TFile *dEdxFile) {
  // TODO: Implement the function
  return 0.0; // Placeholder return
}

double calculateNsigmaElectron(std::shared_ptr<Particle> target,
                               TFile *dEdxFile) {
  // TODO: Implement the function
  return 0.0; // Placeholder return
}

double calculateNsigmaProton(std::shared_ptr<Particle> target,
                             TFile *dEdxFile) {
  // TODO: Implement the function
  return 0.0; // Placeholder return
}

void addTransverseMomentumLoss(std::shared_ptr<Particle> target) {
  // TODO: Implement the function
}

void addPtGaussianError(std::shared_ptr<Particle> target) {
  // TODO: Implement the function
}

void setPtEtaPhiM(std::shared_ptr<Particle> particle, double pt, double eta,
                  double phi, double M) {
  // TODO: Implement the function
}

void setPxPyPzE(std::shared_ptr<Particle> particle, double px, double py,
                double pz, double E) {
  // TODO: Implement the function
}

} // namespace ParticleRoutines
