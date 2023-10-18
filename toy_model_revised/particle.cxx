#include "particle.h"

// Initialize the static counter for unique IDs
int Particle::counter = 0;

Particle::Particle(
    int pdgCode,
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiMVector> fourVector) {
  this->pdgCode = pdgCode;
  this->fourVector = fourVector;
  this->UID = counter++;
}
void Particle::addDaughter(std::shared_ptr<Particle> &daughter) {
  daughters.push_back(daughter);
}

std::vector<std::shared_ptr<Particle>> &Particle::getDaughters() {
  return daughters;
}

int Particle::getUID() { return this->UID; }

int Particle::getPdgCode() { return this->pdgCode; }

bool Particle::getRecoStatus() { return this->reconstructed; }

void Particle::setReconstructed(bool recoStatus) {
  this->reconstructed = recoStatus;
}
