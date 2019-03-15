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
    plasEnergy.resize(sizeDet);
    silicPos.resize(sizeDet);
}
void EventAction::WriteFileVec(ofstream &file, G4ThreeVector &vec)
{
    float x = vec.x();
    float y = vec.y();
    float z = vec.z();
    file.write((char*)&x,sizeof (x));
    file.write((char*)&y,sizeof (y));
    file.write((char*)&z,sizeof (z));
}

void EventAction::BeginOfEventAction(const G4Event* evt)
{  


    G4int eventID = evt->GetEventID();
    if ( eventID % fPrintModulo == 0 )  {
        G4cout << "\n---> Begin of event: " << eventID << G4endl;
        //CLHEP::HepRandom::showEngineStatus();
    }
    for(unsigned int i=0; i < plasEnergy.size(); i++)
    {
        plasEnergy[i]=0;
        silicPos[i]={0.,0.,0.};
    }
    // initialisation per event
    fEnergyWolf = 0.;
    fEnergyPlas = 0.;
    fTrackLWolf = 0.;
    fTrackLPlas= 0.;
}


void EventAction::EndOfEventAction(const G4Event* evt)
{
    int32_t eventID = evt->GetEventID();
    G4PrimaryVertex* primaryVertex = evt->GetPrimaryVertex();
    G4PrimaryParticle* primaryParticle = primaryVertex->GetPrimary();
    float primEnergy = primaryParticle->GetKineticEnergy();
    G4cout << primEnergy << endl;
    //    if ( eventID % fPrintModulo == 0 )  {
    //        G4cout << "\n---> End of event: " << eventID << '\t'<< G4endl;}
    if (!filespec.is_open())
    {
        std::runtime_error("Can't open output file");
    }

    if(checkEmptyEvent)
    {
        //filespec.write((char*)&eventID,sizeof (eventID));
        filespec.write((char*)&primEnergy,sizeof (primEnergy));
        WriteFileVec(filespec,parMomentum[eventID]);
        WriteFileVec(filespec,parPosition[eventID]);
        int32_t count=0;
        for(auto edep : plasEnergy)
        {
            if(edep != 0.)
                count++;
        }
        filespec.write((char*)&count,sizeof(count));
        for(int32_t i=0; i < plasEnergy.size(); i++)
        {
            if(plasEnergy[i]!=0)
            {
                filespec.write((char*)&i,sizeof (i));
                float edep = plasEnergy[i];
                filespec.write((char*)&edep,sizeof (edep));
                WriteFileVec(filespec,silicPos[i]);
            }
            else
                continue;
        }
    }
    checkEmptyEvent=false;


}  
EventAction::~EventAction()
{
    cout<<"num of copies: "<<plasEnergy.size()<<endl;
    filespec.close();
}

