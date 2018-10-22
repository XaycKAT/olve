#pragma once
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class DetectorConstruction;
class G4ParticleGun;
class G4Event;


class PrimaryGeneratorGun : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorGun(DetectorConstruction*);
   ~PrimaryGeneratorGun();

  public:
    void GeneratePrimaries(G4Event*);

  private:
    G4ParticleGun* particleGun;
    DetectorConstruction* myDetector;
};


