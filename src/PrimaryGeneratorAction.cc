#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include<G4Event.hh>
#include<G4ParticleGun.hh>
#include<G4ParticleTable.hh>
#include<G4ParticleDefinition.hh>
#include<globals.hh>
#include<G4GeneralParticleSource.hh>

#include "G4Proton.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* myDC):myDetector(myDC)
{
	//particleGun = new G4GeneralParticleSource();
	G4int n_particle = 1;
	particleGun = new G4ParticleGun(n_particle);

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle("p");
  
    particleGun->SetParticleDefinition(G4Proton::ProtonDefinition());
	particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
    particleGun->SetParticlePosition(G4ThreeVector(0*cm, 0*cm, 35*cm));
	particleGun->SetParticleEnergy(130*GeV);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
	particleGun->GeneratePrimaryVertex(anEvent);
}
