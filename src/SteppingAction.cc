
#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
using namespace std;

SteppingAction::SteppingAction(
        const DetectorConstruction* detectorConstruction,
        EventAction* eventAction)
    : G4UserSteppingAction(),
      fDetConstruction(detectorConstruction),
      fEventAction(eventAction)
{
}


SteppingAction::~SteppingAction()
{ 
}


void SteppingAction::UserSteppingAction(const G4Step* step)
{
    // Collect energy and track length step by step

    // get volume of the current step
    G4VPhysicalVolume* volume
            = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

    G4StepPoint* preStepPoint = step->GetPreStepPoint();
    G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
    //G4int copyNo = theTouchable->GetCopyNumber();
    G4int depth = theTouchable->GetHistoryDepth();
    if ( depth == 2 ) depth = 1;
    G4int cellCopyNo = theTouchable->GetCopyNumber(depth);
    G4String name=theTouchable->GetVolume()->GetName();
    G4ThreeVector worldPosition = preStepPoint->GetPosition();
    G4ThreeVector localPosition = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPosition);

    // energy deposit

    G4double edep = step->GetTotalEnergyDeposit();
    G4String VolName= step->GetTrack()->GetVolume()->GetName();
    //double charge=step->GetTrack()->GetDefinition()->GetPDGCharge();

    const G4Track* track = step->GetTrack();
    G4ThreeVector pMomentum;
    G4ThreeVector pPosition;
    G4ThreeVector pPositionEnd;

    //cout<<cellCopyNo<<'\t'<<fixed<<setprecision(2)<<VolName<<'\t'<<track->GetDefinition()->GetParticleName()
    //   <<'\t' <<track->GetTrackID()<<'\t' <<worldPosition<<'\t'<<"E:"<<edep<<endl;
    if (track->GetTrackID() == 1) {
    pMomentum =track->GetMomentumDirection();
    pPosition = track->GetVertexPosition();
    pPositionEnd = track->GetPosition();
    }

    if (volume == fDetConstruction->GetPlasPV() || volume == fDetConstruction->GetSilicPV()) {
        fEventAction->AddPlas(edep,cellCopyNo);
    }
    if (volume == fDetConstruction->GetSilicPV() && track->GetTrackID() == 1) {
        fEventAction->AddSilicPos(worldPosition,cellCopyNo);
    }
    G4int preEvt = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
    if(preEvt - aftEvt == 1){
        fEventAction->AddMomentum(pMomentum,pPosition,pPositionEnd);
    }
    aftEvt=G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
}

