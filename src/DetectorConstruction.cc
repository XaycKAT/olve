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

#include <math.h>
#include <vector>
using namespace std;



DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction(),
      fWolfPV(0),
      fPlasPV(0),
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
    new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                   kStateGas, 2.73*kelvin, 3.e-18*pascal);
    a = 28.09*g/mole;
    G4Element* elSi = new G4Element("elSi", "Si", 14., a);
    G4Material* Sci = new G4Material("silicon", density, 2);
    Sci->AddElement(elSi, 9);
    //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
    // Geometry parameters

    G4double worldSizeXY = 300*cm;
    G4double worldSizeZ  = 300*cm;


    //double Pi=3.14159;
    G4double phiStart = 0;
    G4double phiTotal= 2*M_PI;
    G4double numSide = 6;	//шестигран призма
    //wolfram
    int numZPlanesWolf=6;
    G4double wolfOut=19.101/2.*mm;   //внешний радиус вольфр призмы (вписанная окружность)
    G4double wolfIn=19.1/2.*mm;    //внутр радиус
    G4double wolfLength=20*mm;	//длина призмы
    G4double wolfLengthF=2*mm;	//выступающий слой
    G4double fullWolfLength=wolfLength+wolfLengthF;
    G4double fullWolfLengthF=wolfLength+2*wolfLengthF;
    G4double *rwolfIn = new G4double[numZPlanesWolf]{ 0,0,wolfIn,wolfIn,0,0};
    G4double *rwolfOut = new G4double[numZPlanesWolf]{ wolfOut, wolfOut,wolfOut,wolfOut, wolfOut, wolfOut};
    G4double *zPlaneWolf = new G4double[numZPlanesWolf]{ 0,wolfLengthF,wolfLengthF,fullWolfLength,fullWolfLength,fullWolfLengthF};

    //plastic
    int numZPlanesPlas=2;
    G4double *zPlanePlas = new G4double[numZPlanesPlas]{ wolfLengthF,fullWolfLength};
    G4double plasOut=19.1/2.*mm;
    G4double plasIn=0;
    G4double *rplasIn = new G4double[numZPlanesPlas]{ plasIn, plasIn};
    G4double *rplasOut = new G4double[numZPlanesPlas]{ plasOut, plasOut};

    //cell
    G4int nofCell=5;// кол-во ячеек в главной линии
    int nofCellLayers=2; //количество слоев ячеек
    int sideL=(nofCell+1)/2 ;
    int numZPlanesCell=2;
    G4double cellOut=wolfOut; G4double cellIn=0;
    G4double cellR=2*wolfOut/ sqrt(3); //радиус описанной окружности

    G4double *rcellIn = new G4double[numZPlanesCell]{ cellIn, cellIn};
    G4double *rcellOut = new G4double[numZPlanesCell]{ cellOut, cellOut};
    G4double *zPlaneCell = new G4double[numZPlanesCell]{ 0,fullWolfLengthF};
    G4RotationMatrix* zRot = new G4RotationMatrix; // Rotates X and Z axes only
    zRot->rotateZ(M_PI/6.*rad);
    G4double centerCells=(nofCell-1)*plasOut;

    //Pad Up & Bottom
    G4double padSizeX=10*mm; //размер одного пада
    G4double padSizeZ=10*mm;
    G4double padThick=2*mm;
    G4double padStep=1*cm;  //отступ между плитами падов
    int nofPadY=3;      // количестов плит
    G4double padOutStep=3*cm;   //отступ от призмы
    int nofPadX=static_cast<int>((nofCell/2.+1)*2*wolfOut/(padSizeX)-nofCell%2+1); //кол-во падов по Х
    int nofPadZ=static_cast<int>(nofCellLayers*fullWolfLengthF/padSizeZ+1);
    G4double putPlatesX=centerCells+(1-nofPadX)*padSizeX/2;
    G4double putPlatesY=(1+(1.5 *((nofCell+1)/2 -1)))*cellR + padOutStep;
    G4double putPlatesZ=(nofPadZ*padSizeZ - fullWolfLengthF*nofCellLayers)/2.-padSizeZ/2.;

    //Pad Border
    int nofPadB=static_cast<int>(nofCell*2*wolfOut/padSizeX+1);
    G4double putPlatesBX=centerCells+(1-nofPadB)*padSizeX/2;
    G4double putPlatesBY=-(nofPadB-1)*padSizeX/2.;
    G4double putPlatesBZ=nofCellLayers*fullWolfLengthF+padOutStep+(nofPadY-1)*(padThick+padStep);

    //Pad Diagonlex
    G4double putPlatesDX1=wolfOut+cellR/sqrt(3)/2.+padOutStep*sqrt(3)/2.;
    G4double putPlatesDY1=padOutStep/2.;//sqrt(3)*wolfOut*nofCell/4.;

    G4double putPlatesDX2=(nofCell-1/2.)*wolfOut*2+cellR/sqrt(3)/2.+padOutStep*sqrt(3)/2.;
    G4double putPlatesDY2=(padOutStep+(padStep)*(nofPadY-1))/2.;
    G4double putPlatesDX3=wolfOut+cellR/sqrt(3)/2.+(padOutStep+(padStep)*(nofPadY-1))*sqrt(3)/2.;
    G4double putPlatesDX4=(nofCell-1/2.)*wolfOut*2+cellR/sqrt(3)/2.+(padOutStep+(padStep)*(nofPadY-1))*sqrt(3)/2.;


    // Get materials
    G4Material* defaultMaterial = G4Material::GetMaterial("Galactic");
    G4Material* wolfMaterial = G4Material::GetMaterial("wolfram");
    G4Material* plasticMaterial = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    G4Material* siliconMat = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

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
                "World");         // its name

    G4VPhysicalVolume* worldPV
            = new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(),  // at (0,0,0)
                worldLV,          // its logical volume
                "World",          // its name
                0,                // its mother  volume
                false,            // no boolean operation
                0,                // copy number
                fCheckOverlaps);  // checking overlaps


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
                "Cell");         // its name



    //
    // plastic
    //


    G4Polyhedra* plastic
            = new G4Polyhedra("plastic",
                              phiStart,phiTotal, numSide, numZPlanesPlas,  zPlanePlas, rplasIn,rplasOut );

    G4LogicalVolume *plasLV = new G4LogicalVolume(plastic,plasticMaterial,"plastic");

    fPlasPV = new G4PVPlacement(
                zRot,                // no rotation
                G4ThreeVector(0, 0., 0),  // at (0,0,0)
                plasLV,          // its logical volume
                "plastic",    // its name
                cellLV,          // its mother  volume
                false,            // no boolean operation
                0,                // copy number
                fCheckOverlaps);  // checking overlaps



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
                "sWolf");          // its name


    fWolfPV
            = new G4PVPlacement(
                zRot,                // no rotation
                G4ThreeVector(0, 0., 0), // its position
                swolfLV,       // its logical volume
                "sWolf",           // its name
                cellLV,          // its mother  volume
                false,            // no boolean operation
                0,                // copy number
                fCheckOverlaps);  // checking overlaps


    // Define one layer as one assembly volume
    G4AssemblyVolume* assemblyCell = new G4AssemblyVolume();
    // Rotation and translation of a plate inside the assembly
    G4RotationMatrix Ra;
    G4ThreeVector Ta;
    G4Transform3D Tr;
    // Rotation of the assembly inside the world
    G4RotationMatrix Rm;

    // Fill the assembly by the plates
    assemblyCell->AddPlacedVolume( cellLV, Tr );

    int check=0;
    // создание массива ячеек в форме правильной призмы
    ofstream file("position.dat");

    for(  int j = 0; j < nofCellLayers; j++ )
    {
        for(  int k = -sideL+1; k < sideL; k++ )
        {
            for(  int i = 0; i < nofCell-abs(k); i++ )
            {
                // Translation of the assembly inside the world
                G4ThreeVector Tm( i*(0 + 2*wolfOut) + abs(k)*wolfOut, k*3*wolfOut/sqrt(3),j*(0 + fullWolfLengthF));
                Tr = G4Transform3D(Rm,Tm);
                assemblyCell->MakeImprint( worldLV, Tr );
                file << check <<'\t'<<setprecision(4)<< Tm <<endl;
                check++;
            }
        }
    }


    cout<<check<<endl;
    sizeDet=check;



    G4Box* siPlate=new G4Box("siPlate", 300*cm,300*cm,300*cm);
    G4LogicalVolume* siPlateLV=new G4LogicalVolume(siPlate, siliconMat, "siPlateLV");

    G4Box* siPad= new G4Box("siPad",           // its name
                            padThick/2., padSizeX/2., padSizeX/2.); // its size

    G4LogicalVolume* siPadLV
            = new G4LogicalVolume(
                siPad,           // its solid
                defaultMaterial,  // its material
                "siPadLV");         // its name

    G4AssemblyVolume* assemblyPlate = new G4AssemblyVolume();


    for( int j=0; j < nofPadY; j++)
    {

        for(  int i = 0; i < nofPadX; i++ )
        {
            for( int k=0; k < nofPadZ; k++)
            {
                G4ThreeVector Tm( i*padSizeX,j*padStep,k*padSizeZ);
                G4RotationMatrix* zRot = new G4RotationMatrix;
                zRot->rotateZ(M_PI/2.*rad);
                Rm=zRot->invert();
                Tr = G4Transform3D(Rm,Tm);
                assemblyPlate->AddPlacedVolume( siPadLV, Tr );
            }
        }
    }
    G4ThreeVector Tm(putPlatesX,putPlatesY,-putPlatesZ);    //верхние пады
    G4RotationMatrix Rm1;
    Tr=G4Transform3D(Rm1,Tm);
    assemblyPlate->MakeImprint(worldLV, Tr);

    Tm={putPlatesX,-putPlatesY-((nofPadY-1)*(padThick+padStep)),-putPlatesZ};   //нижние пады
    Tr=G4Transform3D(Rm1,Tm);
    assemblyPlate->MakeImprint(worldLV, Tr);

    {
        Tm={-putPlatesDX1,putPlatesDY1,0};  // пады сверху слева
        G4RotationMatrix* zRot = new G4RotationMatrix;
        zRot->rotateZ(-M_PI/3.*rad);
        Rm1=zRot->invert();
        Tr=G4Transform3D(Rm1,Tm);
        assemblyPlate->MakeImprint(worldLV, Tr);
    }
    {
        Tm={putPlatesDX2,-putPlatesDY1,0};  //пады снизу справа
        G4RotationMatrix* zt = new G4RotationMatrix;
        zt->rotateZ(2*M_PI/3.*rad);
        Rm1=zt->invert();
        Tr=G4Transform3D(Rm1,Tm);
        assemblyPlate->MakeImprint(worldLV, Tr);
    }
    {
        Tm={-putPlatesDX3,-putPlatesDY2,0}; //пады снизу слева
        G4RotationMatrix* zt = new G4RotationMatrix;
        zt->rotateZ(M_PI/3.*rad);
        Rm1=zt->invert();
        Tr=G4Transform3D(Rm1,Tm);
        assemblyPlate->MakeImprint(worldLV, Tr);
    }
    {
        Tm={putPlatesDX4,putPlatesDY2,0};   //пады сверху справа
        G4RotationMatrix* zt = new G4RotationMatrix;
        zt->rotateZ(-2*M_PI/3.*rad);
        Rm1=zt->invert();
        Tr=G4Transform3D(Rm1,Tm);
        assemblyPlate->MakeImprint(worldLV, Tr);
    }


    G4AssemblyVolume* assemblyPlateB = new G4AssemblyVolume();

    for( int j = 0; j < nofPadY; j++)
    {

        for(  int i = 0; i < nofPadB; i++ )
        {
            for( int k = 0; k < nofPadB; k++)
            {
                G4ThreeVector Tm(i*padSizeX,j*padStep,k*padSizeZ );
                G4RotationMatrix* zRot = new G4RotationMatrix;
                zRot->rotateZ(M_PI/2.*rad);
                Rm=zRot->invert();
                Tr = G4Transform3D(Rm,Tm);
                assemblyPlateB->AddPlacedVolume( siPadLV, Tr );

            }
        }
    }
    Tm={putPlatesBX,putPlatesBY,-padOutStep-padThick/2.}; // пады с лицевого торца
    G4RotationMatrix* xRot = new G4RotationMatrix;
    xRot->rotateX(M_PI/2.*rad);
    Rm=xRot->invert();
    Tr = G4Transform3D(Rm,Tm);
    assemblyPlateB->MakeImprint(worldLV, Tr);

    Tm={putPlatesBX,putPlatesBY,putPlatesBZ};   // пады с заднего торца
    Tr = G4Transform3D(Rm,Tm);
    assemblyPlateB->MakeImprint(worldLV, Tr);


    worldLV->SetVisAttributes(G4VisAttributes::Invisible);
    G4VisAttributes* plasColour= new G4VisAttributes(G4Colour(0,1.0,0));
    G4VisAttributes* wolfColour= new G4VisAttributes(G4Colour(1, 0, 0));
    G4VisAttributes* silicColour= new G4VisAttributes(G4Colour(0, 0, 1));
    plasColour->SetVisibility(true);
    wolfColour->SetVisibility(true);
    silicColour->SetVisibility(true);
    plasLV->SetVisAttributes(plasColour);
    swolfLV->SetVisAttributes(wolfColour);
    //silicLV->SetVisAttributes(silicColour);


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
