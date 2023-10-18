#pragma once

#include "particle.h" // Assuming your Particle class is defined in Particle.h
#include <TFile.h>
#include <TH1F.h>
#include <memory>
#include <vector>

namespace ParticleRoutines {
// Function declarations for particle routines
std::vector<std::shared_ptr<Particle>> twoBodyDecay(std::shared_ptr<Particle> parent, double d1Mass,
                  double d2Mass);

std::shared_ptr<Particle> generateParticlePtEtaPhiM(int pdgCode, double pt,
                                                    double eta, double phi,
                                                    double M);

std::shared_ptr<Particle> generateParticleSlight(int pdgCode,
                                                 TFile *slightFile);

double calculateNsigmaKaon(std::shared_ptr<Particle> target, TFile *dEdxFile);

double calculateNsigmaPion(std::shared_ptr<Particle> target, TFile *dEdxFile);

double calculateNsigmaElectron(std::shared_ptr<Particle> target,
                               TFile *dEdxFile);

double calculateNsigmaProton(std::shared_ptr<Particle> target, TFile *dEdxFile);


void addTransverseMomentumLoss(std::shared_ptr<Particle> target);

void addPtGaussianError(std::shared_ptr<Particle> target);

void setPtEtaPhiM(std::shared_ptr<Particle> particle, double pt, double eta,
                  double phi, double M);

void setPxPyPzE(std::shared_ptr<Particle> particle, double px, double py,
                double pz, double E);

}; // namespace ParticleRoutines
