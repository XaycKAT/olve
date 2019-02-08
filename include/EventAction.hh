#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "DetectorConstruction.hh"
#include <iostream>
#include <vector>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

using namespace std;
class RunAction;

class G4GenericMessenger;

class EventAction : public G4UserEventAction, public DetectorConstruction
{
private:
    G4GenericMessenger*  fMessenger;
    RunAction*  fRunAction;
    ofstream filespec{"/opt/kurbanov_out/spectrumTest.dat",ofstream::binary};
    

    G4double  fEnergyWolf;
    G4double  fEnergyPlas;
    G4double  fTrackLWolf;
    G4double  fTrackLPlas;
    G4int     fPrintModulo;

    
public:
    bool checkEmptyEvent=false;
    EventAction();
    void SetSize ( int sizeDet );
    vector<G4double> plasEnergy;
    vector<G4ThreeVector> silicPos;
    vector<G4ThreeVector> parMomentum;
    vector<G4ThreeVector> parPosition;
    vector<G4ThreeVector> parPositionEnd;

    virtual ~EventAction();

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);
    
    void AddPlas(G4double de, int num);
    void AddSilicPos(G4ThreeVector pos,int num);
    void AddMomentum(G4ThreeVector vec, G4ThreeVector pos, G4ThreeVector endpos);
    void SetPrintModulo(G4int value);
    void WriteFileVec(ofstream &,G4ThreeVector &vec);

};


inline void EventAction::AddPlas(G4double de, int num) {

    plasEnergy[num]+=de;
    if(de > 0.001)
        checkEmptyEvent=true;

}
inline void EventAction::AddSilicPos(G4ThreeVector pos, int num){

    silicPos[num]=pos;

}

inline void EventAction::AddMomentum(G4ThreeVector vec, G4ThreeVector pos, G4ThreeVector endpos){
    parMomentum.push_back(vec);
    parPosition.push_back(pos);
    parPositionEnd.push_back(endpos);
}


inline void EventAction::SetPrintModulo(G4int value) {
    fPrintModulo = value;
}


#endif


