#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include<G4Event.hh>
#include<G4ParticleGun.hh>
#include<G4ParticleTable.hh>
#include<G4ParticleDefinition.hh>
#include<globals.hh>
#include "G4GeneralParticleSource.hh"

#include "G4Proton.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* myDC):myDetector(myDC)
{
    particleGun = new G4GeneralParticleSource();

}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
	particleGun->GeneratePrimaryVertex(anEvent);
}
