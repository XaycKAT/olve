#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
//#include "SteppingVerbose.hh"
//#include "PhysicsList.hh"
#include<G4RunManager.hh>
#include<G4UImanager.hh>
#include<G4UIterminal.hh>
#include<G4VisExecutive.hh>
#include<G4Material.hh>
#include<G4UserRunAction.hh>
#include<G4Run.hh>

#include<iostream>
#include<string>
#include<CLHEP/Random/Random.h>
#include<unistd.h>
#include<time.h>
#include "FTFP_BERT.hh"


using namespace std;

const char macros[]="vis.mac";

/*class RunAction: public G4UserRunAction
{
public:
  void BeginOfRunAction(const G4Run* aRun)
  {
    G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  }
};*/

int main(int argc,char** argv)
{

    //G4VSteppingVerbose::SetInstance(new SteppingVerbose);
    CLHEP::HepRandom::setTheSeed(time(0)+getpid());

    G4RunManager * runManager = new G4RunManager;
    DetectorConstruction* detConstruction = new DetectorConstruction;
    runManager->SetUserInitialization(detConstruction);

    //G4VUserPhysicsList *p = new PhysicsList;
    //runManager->SetUserInitialization(p);
    //runManager->SetUserInitialization(new FTFP_BERT);
    G4VModularPhysicsList* physicsList = new FTFP_BERT;
    runManager->SetUserInitialization(physicsList);
    // Set user action classes
    runManager->SetUserAction(new PrimaryGeneratorAction(detConstruction));
    runManager->SetUserAction(new RunAction);
    EventAction* eventAction = new EventAction();
    runManager->SetUserAction(eventAction);
    SteppingAction* steppingAction  = new SteppingAction(detConstruction, eventAction);
    runManager->SetUserAction(steppingAction);
  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();


  // get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();


    runManager->Initialize();

    eventAction->SetSize ( detConstruction->sizeDet );

    //cout<<"===============================================================";
    //cout<<endl;
    //cout<< *(G4Material::GetMaterialTable()) << endl;
    //cout<<"===============================================================";
    //cout<<endl;

    G4UImanager * UI = G4UImanager::GetUIpointer();
    G4UIsession * session = new G4UIterminal();

    UI->ExecuteMacroFile(macros);
    session->SessionStart();
    UImanager->ApplyCommand("/control/execute gps1.mac");
    delete session;
    delete visManager;
    delete runManager;

    return 0;
}
