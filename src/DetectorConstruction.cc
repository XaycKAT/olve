#include"DetectorConstruction.hh"
#include<G4NistManager.hh>
#include<G4Box.hh>
#include<G4Tubs.hh>
#include<G4Trap.hh>
#include<G4Polyhedra.hh>
#include<G4LogicalVolume.hh>
#include<G4PVPlacement.hh>
#include<G4SDManager.hh>
#include<G4VisAttributes.hh>
#include "G4UImessenger.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4AssemblyVolume.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4GenericMessenger.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Polyline.hh"
#include "G4Point3D.hh"
#include "G4String.hh"
#include "G4Circle.hh"

#include <math.h>
#include <vector>
#include <array>
#define PI 3.14159265
using namespace std;



DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction(),
      fWolfPV(0),
      fPlasPV(0),
      fSilicPV(0),
      fMessenger(0),
      fMagField(0),
      fCheckOverlaps(true)
{
    fMessenger
            = new G4GenericMessenger(this, "/B4/det/", "Detector construction control");

    G4GenericMessenger::Command& setMagFieldCmd
            = fMessenger->DeclareMethod("setMagField",
                                        &DetectorConstruction::SetMagField,
                                        "Define magnetic field value (in X direction");
    setMagFieldCmd.SetUnitCategory("Magnetic flux density");
}

DetectorConstruction::~DetectorConstruction()
{
}

int DetectorConstruction::WriteFile(G4AssemblyVolume* av, ofstream &file, int count)
{
    auto it = av->GetVolumesIterator();
    //advance(it,count);
    auto size = av->TotalImprintedVolumes();
    for (unsigned int i = count; i < size; it++, ++i)
    {
        //auto &volume = *it;
        G4ThreeVector t=(*it)->GetObjectTranslation();
        file <<(*it)->GetCopyNo()<<'\t'<<setprecision(4)<< t<<endl;
    }
    return size;
}
G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // Define materials
    DefineMaterials();

    // Define volumes
    return DefineVolumes();
}
void DetectorConstruction::DefineMaterials()
{
    // Lead material defined using NIST Manager
    G4NistManager* nistManager = G4NistManager::Instance();
    nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

    // wolfram material


    G4double a = 183.84*g/mole;  // mass of a mole;
    G4double z=74.;  // z=mean number of protons;
    G4double density = 17.1*g/cm3;
    new G4Material("wolfram", z, a, density);


    // Vacuum
    new G4Material("Galactic", 1., 1.01*g/mole,density= universe_mean_density,
                   kStateGas, 2.73*kelvin, 3.e-18*pascal);
    a = 28.09*g/mole;
    z=14.;
    density=2.33*g/cm3;
    new G4Material("silicon", z, a,density);

    //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
    // Geometry parameters

    G4double worldSizeXY = 3000*cm;
    G4double worldSizeZ  = 3000*cm;


    //double Pi=3.14159;
    G4double phiStart = 0;
    G4double phiTotal= 2*M_PI;
    G4double numSide = 6;	//шестигран призма
    //wolfram
    int numZPlanesWolf=2;
    G4double wolfOut=25.01/2.*sqrt(3)/2.*mm;   //внешний радиус вольфр призмы (вписанная окружность)
    G4double wolfIn=0;   //внутр радиус
    G4double wolfLengthF=2*mm;	//выступающий слой
    G4double *rwolfIn = new G4double[numZPlanesWolf]{ wolfIn,wolfIn};
    G4double *rwolfOut = new G4double[numZPlanesWolf]{ wolfOut, wolfOut};
    G4double *zPlaneWolf = new G4double[numZPlanesWolf]{ 0,wolfLengthF};

    //plastic
    int numZPlanesPlas=numZPlanesWolf;
    G4double plasLength=20*mm;
    G4double *zPlanePlas = new G4double[numZPlanesPlas]{ wolfLengthF,plasLength+wolfLengthF};
    G4double plasOut=wolfOut;
    G4double plasIn=0;
    G4double *rplasIn = new G4double[numZPlanesPlas]{ plasIn, plasIn};
    G4double *rplasOut = new G4double[numZPlanesPlas]{ plasOut, plasOut};

    //cell
    G4int nofCell=19;//static_cast<int>(1600./wolfOut/2);// кол-во ячеек в главной линии
    int nofCellLayers=1; //количество слоев ячеек
    int sideL=(nofCell+1)/2 ;
    int numZPlanesCell=numZPlanesWolf;
    G4double cellOut=wolfOut;
    G4double cellIn=0;
    //G4double cellR=2*wolfOut/ sqrt(3); //радиус описанной окружности
    G4double fullCellLength=plasLength+2*wolfLengthF;
    G4double *rcellIn = new G4double[numZPlanesCell]{ cellIn, cellIn};
    G4double *rcellOut = new G4double[numZPlanesCell]{ cellOut, cellOut};
    G4double *zPlaneCell = new G4double[numZPlanesCell]{ 0,fullCellLength};
    G4RotationMatrix* zRot = new G4RotationMatrix; // Rotates X and Z axes only
    zRot->rotateZ(M_PI/6.*rad);
    G4double centerCells=(nofCell-1)*wolfOut;

    //Pad Up & Bottom
    G4double padSizeX=20*mm; //размер одного пада
    G4double padSizeZ=20*mm;
    G4double padThick=0.5*mm;
    G4double padStep=3*cm;  //отступ между плитами падов
    int nofPadY=4;      // количестов плит
    G4double padOutStep=5*cm;   //отступ от призмы
    double p=2*padStep/sqrt(3)/padSizeX;
    double ch;
    if(p-floor(p)>0.5)
        ch=0.01;
    else
        ch=-0.01;
    while(p-floor(p)>0.01)
    {
        padStep+=ch;
        p=2*padStep/sqrt(3)/padSizeX;
    }
    p=(nofCell*wolfOut+2*(padOutStep)/sqrt(3))/padSizeX;
    if(p-floor(p)>0.5)
        ch=0.01;
    else
        ch=-0.01;
    while(p-floor(p)>0.01)
    {
        padOutStep+=ch;
        p=(nofCell*wolfOut+2*(padOutStep)/sqrt(3))/padSizeX;
    }
    int nofPadX=static_cast<int>((nofCell*wolfOut+2*(padOutStep)/sqrt(3))/padSizeX); //кол-во падов по Х
    int nofPadZ=static_cast<int>((nofCellLayers*fullCellLength+2*padOutStep)/padSizeZ);
    G4double putPlatesX=centerCells-nofPadX*padSizeX/2.+padSizeX/2.;
    G4double putPlatesY=nofPadX*sin(PI/3.)*padSizeX;
    G4double putPlatesZ=-((nofPadZ*padSizeZ - fullCellLength*nofCellLayers)/2. - padSizeZ/2.);
    G4double putPlatesZ1=-putPlatesZ+fullCellLength*nofCellLayers;

    //Pad Diagonlex
    G4double delta=nofPadX*padSizeX-wolfOut-padOutStep*2/sqrt(3)-centerCells;
    G4double putPlatesDX1=wolfOut+padOutStep*2/sqrt(3)-padSizeX/2.*cos(PI/3.);
    G4double putPlatesDY1=padSizeX*sqrt(3)/4.;
    G4double putPlatesDX2=(nofCell-1/2.)*wolfOut*2+padOutStep*2/sqrt(3)-padSizeX/2.*cos(PI/3.);


    //Pad Border
    int nofPadBX=static_cast<int>((nofCell*2*wolfOut+padOutStep*sqrt(5))/padSizeX+3);
    int nofPadBZ=static_cast<int>((putPlatesY*2)/padSizeZ+3);
    G4double putPlatesBX=centerCells+(1-nofPadBX)*padSizeX/2;
    G4double putPlatesBY=0;
    G4double putPlatesBZ=-(padOutStep-padSizeX/2);
    G4double putPlatesBZ1=((nofPadY-1)*(padStep)+nofCellLayers*fullCellLength+padOutStep-padSizeX/2.);

    // Get materials
    G4Material* defaultMaterial = G4Material::GetMaterial("Galactic");
    G4Material* wolfMaterial = G4Material::GetMaterial("wolfram");
    G4Material* plasticMaterial = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    G4Material* siliconMat = G4Material::GetMaterial("silicon");

    if ( ! defaultMaterial || ! wolfMaterial || ! plasticMaterial ) {
        G4ExceptionDescription msg;
        msg << "Cannot retrieve materials already defined.";
        G4Exception("B4DetectorConstruction::DefineVolumes()",
                    "MyCode0001", FatalException, msg);
    }
    //
    // WORLD
    //

    G4VSolid* world
            = new G4Box("World",           // its name
                        worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size

    G4LogicalVolume* worldLV
            = new G4LogicalVolume(
                world,           // its solid
                defaultMaterial,  // its material
                "WorldLV");         // its name

    G4VPhysicalVolume* worldPV
            = new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(),  // at (0,0,0)
                worldLV,          // its logical volume
                "WorldPV",          // its name
                0,                // its mother  volume
                false,            // no boolean operation
                0,                // copy number
                fCheckOverlaps);  // posFileSizeing overlaps


    //
    // Cell
    //


    G4Polyhedra* cell
            = new G4Polyhedra("cell",
                              phiStart+M_PI/6.*rad,phiTotal, numSide, numZPlanesCell,  zPlaneCell, rcellIn,rcellOut );

    G4LogicalVolume* cellLV
            = new G4LogicalVolume(
                cell,           // its cell
                defaultMaterial,  // its material
                "CellLV");         // its name



    //
    // plastic
    //


    G4Polyhedra* plastic
            = new G4Polyhedra("plastic",
                              phiStart,phiTotal, numSide, numZPlanesPlas,  zPlanePlas, rplasIn,rplasOut );

    G4LogicalVolume *plasLV = new G4LogicalVolume(plastic,plasticMaterial,"plasLV");

    fPlasPV = new G4PVPlacement(
                zRot,               // no rotation
                G4ThreeVector(0, 0., 0),  // at (0,0,0)
                plasLV,             // its logical volume
                "plasPV",           // its name
                cellLV,             // its mother  volume
                false,              // no boolean operation
                0,                  // copy number
                fCheckOverlaps);    // posFileSizeing overlaps



    //
    // wolfram
    //
    G4Polyhedra *sWolf
            = new G4Polyhedra("sWolf",            // its name
                              phiStart,phiTotal, numSide, numZPlanesWolf,  zPlaneWolf, rwolfIn,rwolfOut ); // its size*/


    G4LogicalVolume* swolfLV
            = new G4LogicalVolume(
                sWolf,        // its solid
                wolfMaterial, // its materialM_PI/6.*rad
                "WolfLV");          // its name


    fWolfPV
            = new G4PVPlacement(
                zRot,                // no rotation
                G4ThreeVector(0, 0., 0), // its position
                swolfLV,       // its logical volume
                "sWolfPV",           // its name
                cellLV,          // its mother  volume
                false,            // no boolean operation
                0,                // copy number
                fCheckOverlaps);  // posFileSizeing overlaps
    fWolfPV
            = new G4PVPlacement(
                zRot,                // no rotation
                G4ThreeVector(0, 0., plasLength+wolfLengthF), // its position
                swolfLV,       // its logical volume
                "sWolfPV",           // its name
                cellLV,          // its mother  volume
                false,            // no boolean operation
                0,                // copy number
                fCheckOverlaps);  // posFileSizeing overlaps



    G4AssemblyVolume* assemblyCell = new G4AssemblyVolume();
    G4RotationMatrix Ra;
    G4ThreeVector Tm;
    G4Transform3D Tr;
    G4RotationMatrix Rm;
    assemblyCell->AddPlacedVolume( cellLV, Tr );

    // создание массива ячеек в форме правильной призмы
    ofstream file("position.dat");
    for(  int j = 0; j < nofCellLayers; j++ )
    {
        for(  int k = -sideL+1; k < sideL; k++ )
        {
            for(  int i = 0; i < nofCell-abs(k); i++ )
            {
                Tm={ i*2*wolfOut + abs(k)*wolfOut, k*3*wolfOut/sqrt(3),j*fullCellLength};
                Tr = G4Transform3D(Rm,Tm);
                assemblyCell->MakeImprint( worldLV, Tr );
                //file << posFileSize+1 <<'\t'<<setprecision(4)<< Tm <<endl;

            }
        }
    }
    //пады

    G4Box* siPlate=new G4Box("siPlate", padThick/2., padSizeX/2., padSizeZ/2.);
    G4LogicalVolume* siPlateLV=new G4LogicalVolume(siPlate, defaultMaterial, "siPlateLV");

    G4Box* siPad= new G4Box("siPad",           // its name
                            padThick/2., padSizeX/2., padSizeZ/2.); // its size

    G4LogicalVolume* siPadLV
            = new G4LogicalVolume(
                siPad,           // its solid
                siliconMat,  // its material
                "siPadLV");         // its name
    fSilicPV
            = new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0, 0., 0), // its position
                siPadLV,       // its logical volume
                "silicPV",           // its name
                siPlateLV,          // its mother  volume
                false,            // no boolean operation
                0,                // copy number
                fCheckOverlaps);  // posFileSizeing overlaps

    G4AssemblyVolume* assemblyPlate = new G4AssemblyVolume();

    for( int j=0; j < nofPadY; j++)
    {
        for(  int i = 0; i < nofPadX+j*(2*padStep/padSizeX/sqrt(3))-1*bool(j); i++ )
        {
            for( int k=0; k < nofPadZ; k++)
            {
                Tm={ i*padSizeX-j*(padStep/sqrt(3)),j*padStep,k*padSizeZ};
                zRot = new G4RotationMatrix;
                zRot->rotateZ(M_PI/2.*rad);
                Rm=zRot->invert();
                Tr = G4Transform3D(Rm,Tm);
                assemblyPlate->AddPlacedVolume( siPlateLV, Tr );
            }
        }
    }
    G4AssemblyVolume* assemblyPlateB = new G4AssemblyVolume();

    for( int j = 0; j < nofPadY; j++)
    {
        for(  int i = -nofPadBZ/2+1; i < nofPadBZ/2; i++ )
        {
            for( int k = 0; k < nofPadBX-abs(i); k++)
            {
                Tm={k*padSizeX+abs(i)*padSizeX/2,j*padStep,i*padSizeZ };
                zRot = new G4RotationMatrix;
                zRot->rotateZ(M_PI/2.*rad);
                Rm=zRot->invert();
                Tr = G4Transform3D(Rm,Tm);
                assemblyPlateB->AddPlacedVolume( siPlateLV, Tr );

            }
        }
    }

    Tm={putPlatesBX,putPlatesBY,putPlatesBZ}; // пады с заднего  торца
    G4RotationMatrix* xRot = new G4RotationMatrix;
    xRot->rotateX(M_PI/2.*rad);
    Rm=xRot->invert();
    Tr = G4Transform3D(Rm,Tm);
    assemblyPlateB->MakeImprint(worldLV, Tr);

    Tm={putPlatesBX,putPlatesBY,putPlatesBZ1};   // пады с лицевого торца
    Tr = G4Transform3D(Rm,Tm);
    //assemblyPlateB->MakeImprint(worldLV, Tr);

    Tm={putPlatesX,putPlatesY,putPlatesZ};    //верхние пады
    G4RotationMatrix Rm1;
    Tr=G4Transform3D(Rm1,Tm);
    assemblyPlate->MakeImprint(worldLV, Tr);


    Tm={putPlatesX,-putPlatesY,putPlatesZ1};   //нижние пады
    zRot = new G4RotationMatrix;
    zRot->rotateX(M_PI*rad);
    Rm1=zRot->invert();
    Tr=G4Transform3D(Rm1,Tm);
    assemblyPlate->MakeImprint(worldLV, Tr);


    Tm={-putPlatesDX1,putPlatesDY1,putPlatesZ};  // пады сверху слева
    zRot = new G4RotationMatrix;
    zRot->rotateZ(-M_PI/3.*rad);
    Rm1=zRot->invert();
    Tr=G4Transform3D(Rm1,Tm);
    assemblyPlate->MakeImprint(worldLV, Tr);


    Tm={putPlatesDX2,-putPlatesDY1,putPlatesZ};  //пады снизу справа
    zRot = new G4RotationMatrix;
    zRot->rotateZ(2*M_PI/3.*rad);
    Rm1=zRot->invert();
    Tr=G4Transform3D(Rm1,Tm);
    assemblyPlate->MakeImprint(worldLV, Tr);


    Tm={-putPlatesDX1,-putPlatesDY1,putPlatesZ1}; //пады снизу слева
    zRot = new G4RotationMatrix;
    zRot->rotateZ(M_PI/3.*rad);
    zRot->rotateX(M_PI*rad);
    Rm1=zRot->invert();
    Tr=G4Transform3D(Rm1,Tm);
    assemblyPlate->MakeImprint(worldLV, Tr);


    Tm={putPlatesDX2,putPlatesDY1,putPlatesZ1};   //пады сверху справа
    zRot = new G4RotationMatrix;
    zRot->rotateZ(-2*M_PI/3.*rad);
    zRot->rotateX(M_PI*rad);
    Rm1=zRot->invert();
    Tr=G4Transform3D(Rm1,Tm);
    assemblyPlate->MakeImprint(worldLV, Tr);

    WriteFile(assemblyCell,file,0);
    file<<"#plates"<<endl;
    WriteFile(assemblyPlateB,file,0);
    WriteFile(assemblyPlate,file,0);
    int posFileSize=assemblyPlate->TotalImprintedVolumes()+assemblyPlateB->TotalImprintedVolumes()+
            assemblyCell->TotalImprintedVolumes();

    file.close();
    sizeDet=posFileSize;

    worldLV->SetVisAttributes(G4VisAttributes::Invisible);
    G4VisAttributes* plasColour= new G4VisAttributes(G4Colour(1,1,1));
    G4VisAttributes* wolfColour= new G4VisAttributes(G4Colour(0.5,0.5,0.5));
    G4VisAttributes* silicColour= new G4VisAttributes(G4Colour(1,1,0));
    plasColour->SetVisibility(true);
    wolfColour->SetVisibility(true);
    silicColour->SetVisibility(true);
    plasLV->SetVisAttributes(plasColour);
    swolfLV->SetVisAttributes(wolfColour);
    siPadLV->SetVisAttributes(silicColour);
    cellLV->SetVisAttributes(G4VisAttributes::Invisible);
    siPlateLV->SetVisAttributes(G4VisAttributes::Invisible);
    return worldPV;
}
void DetectorConstruction::SetMagField(G4double fieldValue)
{
    // Apply a global uniform magnetic field along X axis
    G4FieldManager* fieldManager
            = G4TransportationManager::GetTransportationManager()->GetFieldManager();

    // Delete the existing magnetic field
    if ( fMagField )  delete fMagField;

    if ( fieldValue != 0. ) {
        // create a new one if not null
        fMagField
                = new G4UniformMagField(G4ThreeVector(fieldValue, 0., 0.));

        fieldManager->SetDetectorField(fMagField);
        fieldManager->CreateChordFinder(fMagField);
    }
    else {
        fMagField = 0;
        fieldManager->SetDetectorField(fMagField);
    }
}
