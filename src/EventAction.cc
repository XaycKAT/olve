

#include "EventAction.hh"
#include "RunAction.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4GenericMessenger.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

EventAction::EventAction()
    : G4UserEventAction(),
      fMessenger(0),
      fEnergyWolf(0.),
      fEnergyPlas(0.),
      fTrackLWolf(0.),
      fTrackLPlas(0.),
      fPrintModulo(1)
{
    fMessenger = new G4GenericMessenger(this, "/Olve/event/", "Event control");

    G4GenericMessenger::Command& setPrintModulo
            = fMessenger->DeclareProperty("setPrintModulo",
                                          fPrintModulo,
                                          "Print events modulo n");
    setPrintModulo.SetRange("value>0");
}

void EventAction::SetSize(int sizeDet)
{
    PlasEnergySize = sizeDet;
    plasEnergy = new G4double[sizeDet]();
    silicPos = new G4ThreeVector[sizeDet]();
}

void EventAction::BeginOfEventAction(const G4Event* evt)
{  


    G4int eventID = evt->GetEventID();
    if ( eventID % fPrintModulo == 0 )  {
        G4cout << "\n---> Begin of event: " << eventID << G4endl;
        //CLHEP::HepRandom::showEngineStatus();
    }

    // initialisation per event
    fEnergyWolf = 0.;
    fEnergyPlas = 0.;
    fTrackLWolf = 0.;
    fTrackLPlas= 0.;
}


void EventAction::EndOfEventAction(const G4Event* evt)
{
    G4int eventID = evt->GetEventID();
    //    if ( eventID % fPrintModulo == 0 )  {
    //        G4cout << "\n---> End of event: " << eventID << '\t'<< G4endl;}
    if (!filespec.is_open())
    {
        std::runtime_error("Can't open output file");
    }
    if(checkEmptyEvent)
    {
        filespec << "#" << eventID <<setprecision(2)<< fixed <<' '<<
                    parMomentum[eventID]<<' '<<parPosition[eventID]<<' '<<parPositionEnd[eventID]<< endl;
        for(int i = 0; i < PlasEnergySize ; i++)
        {
            if(plasEnergy[i]!=0)
                filespec  <<i<< ' '<<setprecision(3) << plasEnergy[i]<<'\t' << silicPos[i]<< endl;
            else
                continue;
        }
        for(int i = 0; i < PlasEnergySize ; i++)
        {
            plasEnergy[i]=0;
        }
    }
    checkEmptyEvent=false;


}  
EventAction::~EventAction()
{
    cout<<"num of copys: "<<PlasEnergySize<<endl;
    if (silicPos)
    {
        delete[] silicPos;
    }
    if (plasEnergy)
    {
        delete[] plasEnergy;
    }
}

