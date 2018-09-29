
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
    G4int copyNo = theTouchable->GetCopyNumber();
    G4int depth = theTouchable->GetHistoryDepth();
    if ( depth == 2 ) depth = 1;
    G4int cellCopyNo = theTouchable->GetCopyNumber(depth);
    G4String name=theTouchable->GetVolume()->GetName();
    G4ThreeVector worldPosition = preStepPoint->GetPosition();
    G4ThreeVector localPosition = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPosition);

    // energy deposit

    G4double edep = step->GetTotalEnergyDeposit();

    G4String VolName= step->GetTrack()->GetVolume()->GetName();
    double charge=step->GetTrack()->GetDefinition()->GetPDGCharge();

    std::cout<<cellCopyNo<<'\t'<<VolName<<'\t'<<setprecision(2)<<worldPosition<<'\t'<<"Energy:"<<edep<<'\t'<<"Charge:"<<charge<<endl;


    if (volume == fDetConstruction->GetPlasPV() || volume == fDetConstruction->GetSilicPV()) {
        fEventAction->AddPlas(edep,cellCopyNo);
    }

//    if ( volume == fDetConstruction->GetSilicPV() ) {
//        fEventAction->AddSilic(VolName,edep);
//    }

}

