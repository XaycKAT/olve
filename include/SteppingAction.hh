#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include <G4Types.hh>

class DetectorConstruction;
class EventAction;

/// Stepping action class.
///
/// In UserSteppingAction() there are collected the energy deposit and track 
/// lengths of charged particles in Absober and Gap layers and
/// updated in B4aEventAction.

class SteppingAction : public G4UserSteppingAction
{
public:
    SteppingAction(const DetectorConstruction* detectorConstruction,
                   EventAction* eventAction);
    virtual ~SteppingAction();
    G4int pastEvt=-1;
    bool count = true;

    virtual void UserSteppingAction(const G4Step* step);
    
private:
    const DetectorConstruction* fDetConstruction;
    EventAction*  fEventAction;
};


#endif
