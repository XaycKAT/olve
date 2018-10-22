#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "DetectorConstruction.hh"
#include <iostream>
#include <vector>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "map"

using namespace std;
class RunAction;

class G4GenericMessenger;

class EventAction : public G4UserEventAction, public DetectorConstruction
{
private:
    G4GenericMessenger*  fMessenger;
    RunAction*  fRunAction;
    ofstream filespec{"spectrum.dat"};
    

    G4double  fEnergyWolf;
    G4double  fEnergyPlas;
    G4double  fTrackLWolf;
    G4double  fTrackLPlas;
    G4int     fPrintModulo;

    int PlasEnergySize = 0;
    
public:
    int checkEmptyEvent=false;
    EventAction();
    void SetSize ( int sizeDet );
    G4double *plasEnergy = nullptr;
    G4ThreeVector *silicPos = nullptr;
    vector<G4ThreeVector> parMomentum;
    vector<G4ThreeVector> parPosition;
    vector<G4ThreeVector> parPositionEnd;

    //map <G4String,G4double> mapSiPads;
    virtual ~EventAction();

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);
    
    void AddPlas(G4double de, int num);
    void AddSilicPos(G4ThreeVector pos,int num);
    void AddMomentum(G4ThreeVector vec, G4ThreeVector pos, G4ThreeVector endpos);
    void SetPrintModulo(G4int value);

};


//inline void EventAction::AddSilic(G4String name, G4double de) {
//    auto it = mapSiPads.find(name);
//    if(it == mapSiPads.end())
//    {
//        mapSiPads.insert(pair<G4String,G4double>(name,de));
//    }
//    else
//    {
//        it->second += de;
//    }

//}

inline void EventAction::AddPlas(G4double de, int num) {

    plasEnergy[num]+=de;
    if(de>0.0001)
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


