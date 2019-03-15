
#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4AssemblyVolume.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
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

    G4VPhysicalVolume* volume
            = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

    G4StepPoint* preStepPoint = step->GetPreStepPoint();
    G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
    G4int depth = theTouchable->GetHistoryDepth();
    if ( depth == 2 ) depth = 1;
    G4int copyNo = theTouchable->GetCopyNumber(depth);
    G4ThreeVector worldPosition = preStepPoint->GetPosition();
    G4ThreeVector localPosition = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPosition);
    G4double edep = step->GetTotalEnergyDeposit();
    G4String VolName= step->GetTrack()->GetVolume()->GetName();
    const G4Track* track = step->GetTrack();
 //   G4LogicalVolume* lv=theTouchable->GetVolume()->GetLogicalVolume();
 //   lv->GetInstanceID();

//    auto logicVolname = lv->GetName();
//    auto PhysVolname = volume->GetName();
//    G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
//    lvs->GetVolume("CellLV");
//    lvs->GetInstance();
//    vector<G4LogicalVolume*>::const_iterator lvit;
//    string s;
//    for( lvit = lvs->begin(); lvit != lvs->end(); lvit++ )
//    {
//         s = (*lvit)->GetName();
//    }



    G4ThreeVector pMomentum;
    G4ThreeVector pPosition;
    G4ThreeVector pPositionEnd;

//   cout<<copyNo<<'\t'<<fixed<<setprecision(2)<<VolName<<'\t'<<track->GetDefinition()->GetParticleName()
//       <<'\t' <<track->GetTrackID()<<'\t' <<worldPosition<<'\t'<<"E:"<<edep<<endl;
    if (track->GetTrackID() == 1) {
        pMomentum =track->GetMomentumDirection();
        pPosition = track->GetVertexPosition();
        pPositionEnd = track->GetPosition();
    }
    G4int presentEvt = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
    if (volume == fDetConstruction->GetPlasPV() || volume == fDetConstruction->GetSilicPV()) {
        fEventAction->AddPlas(edep,copyNo);
    }
    if(presentEvt - pastEvt == 1)
        count=true;
    if(volume == fDetConstruction->GetPlasPV() )
        count=false;
    if (volume == fDetConstruction->GetSilicPV() && track->GetTrackID() == 1 && track->GetParentID() == 0 && count)  {
        fEventAction->AddSilicPos(worldPosition,copyNo);
    }
    if(presentEvt - pastEvt == 1){
        fEventAction->AddMomentum(pMomentum,pPosition,pPositionEnd);
    }
    pastEvt=G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
}

