#include "PrimaryGeneratorGun.hh"
#include "DetectorConstruction.hh"

#include<G4Event.hh>
#include<G4ParticleGun.hh>
#include<G4ParticleTable.hh>
#include<G4ParticleDefinition.hh>
#include<globals.hh>

#include "G4Proton.hh"

PrimaryGeneratorGun::PrimaryGeneratorGun(DetectorConstruction* myDC):myDetector(myDC)
{
    G4int n_particle = 1;
    particleGun = new G4ParticleGun(n_particle);

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    //G4ParticleDefinition* particle = particleTable->FindParticle("mu+");

    particleGun->SetParticleDefinition(G4Proton::ProtonDefinition());
    particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,-1.,0.));
    particleGun->SetParticlePosition(G4ThreeVector(2*cm, 200*cm, 0*cm));
    particleGun->SetParticleEnergy(130*GeV);
}

PrimaryGeneratorGun::~PrimaryGeneratorGun()
{
    delete  particleGun;
}

void PrimaryGeneratorGun::GeneratePrimaries(G4Event* anEvent)
{
    particleGun->GeneratePrimaryVertex(anEvent);
}
