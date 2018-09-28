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
    ofstream filespec{"spectrum.dat"};
    
   
    G4double  fEnergyWolf;
    G4double  fEnergyPlas;
    G4double  fTrackLWolf; 
    G4double  fTrackLPlas;
    G4int     fPrintModulo;

    int PlasEnergySize = 0;
    
  public:
	
    EventAction();
    void SetSize ( int sizeDet );
    G4double *plasEnergy = nullptr;
    G4double *silicEnergy = nullptr;
    virtual ~EventAction();

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);
    
  //  void AddWolf(G4double de, G4double dl);
    void AddPlas(G4double de, int dl);
                     
    void SetPrintModulo(G4int value);
 
};

// inline functions

/*inline void EventAction::AddWolf(G4double de, G4double dl) {
  fEnergyWolf += de; 
  fTrackLWolf += dl;
}*/

inline void EventAction::AddPlas(G4double de, int dl) {

  plasEnergy[dl]+=de;

}

inline void EventAction::SetPrintModulo(G4int value) {
  fPrintModulo = value;
}
                     

#endif

    
