#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4GeneralParticleSource.hh"

#include "G4VUserPrimaryGeneratorAction.hh"

class DetectorConstruction;
class G4ParticleGun;
class G4Event;

 
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction(DetectorConstruction*);    
   ~PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);

  private:
    G4GeneralParticleSource* particleGun;
    //G4ParticleGun particleGun;
    DetectorConstruction* myDetector;
};


#endif


