﻿#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include<globals.hh>
#include<G4VUserDetectorConstruction.hh>
#include<G4VSolid.hh>
#include<G4LogicalVolume.hh>
#include<G4VPhysicalVolume.hh>
#include<G4Material.hh>
#include<G4VisAttributes.hh>
#include<G4RotationMatrix.hh>
#include "G4ThreeVector.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include<G4AssemblyVolume.hh>
#include<iostream>
#include<fstream>
using namespace std;
class G4Box;
class G4VPhysicalVolume;
class G4UniformMagField;
class G4GenericMessenger;

//extern int sizeDet;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  private:

    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
    //ofstream fileModelInfo{"modelInfo.dat"};
    ofstream filePos{"position.dat",ofstream::binary};

    G4VPhysicalVolume* fPlasPV;      // the plastic physical volume
    G4VPhysicalVolume* fWolfPV; // the wolfram physical volume

    G4VPhysicalVolume* fSilicPV;
    G4GenericMessenger*  fMessenger; // messenger 
    G4UniformMagField*   fMagField;  // magnetic field

    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps

    int WriteFile(G4AssemblyVolume*, int count );
    void SetAssemblyNames(G4AssemblyVolume*, int startnum,int aNum, vector<int> &sizeLayers);
  public:
    DetectorConstruction();
    
    virtual ~DetectorConstruction();
    int sizeDet;
    DetectorConstruction(int );

  public:
    virtual G4VPhysicalVolume* Construct();

    // set methods
    //
    void SetMagField(G4double fieldValue);
    void ConstructSDandField(G4double fieldValue);
    void setSize(int par){sizeDet=par;}
    
    // get methods
    //
    const G4VPhysicalVolume* GetWolfPV() const;
    const G4VPhysicalVolume* GetPlasPV() const;
    const G4VPhysicalVolume* GetSilicPV() const;

};

// inline functions

inline const G4VPhysicalVolume* DetectorConstruction::GetWolfPV() const { 
  return fWolfPV; 
}
inline const G4VPhysicalVolume* DetectorConstruction::GetPlasPV() const { 
  return fPlasPV; 
}
inline const G4VPhysicalVolume* DetectorConstruction::GetSilicPV() const {
  return fSilicPV;
}
  
#endif
