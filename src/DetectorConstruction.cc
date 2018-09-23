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
      // Define /B4/det commands using generic messenger class
      fMessenger
        = new G4GenericMessenger(this, "/B4/det/", "Detector construction control");

      // Define /B4/det/setMagField command
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
      G4double wolfOut=19.1/2*mm;   //внешний радиус вольфр призмы
      G4double wolfIn=19.1/2*mm;    //внутр радиус
      G4double wolfLength=20*mm;	//длина призмы
      G4double wolfLengthF=2*mm;	//выступающий слой
      G4double fullWolfLength=wolfLength+wolfLengthF;
      G4double fullWolfLengthF=wolfLength+2*wolfLengthF;
      G4double rwolfIn[numZPlanesWolf] = { 0,0,wolfIn,wolfIn,0,0};
      G4double rwolfOut[numZPlanesWolf] = { wolfOut, wolfOut,wolfOut,wolfOut, wolfOut, wolfOut};
      G4double zPlaneWolf[numZPlanesWolf] = { 0,wolfLengthF,wolfLengthF,fullWolfLength,fullWolfLength,fullWolfLengthF};

      //plastic
      int numZPlanesPlas=2;
      G4double zPlanePlas[numZPlanesPlas] = { wolfLengthF,fullWolfLength};
      G4double plasOut=19.1/2*mm;
      G4double plasIn=0;
      G4double rplasIn[numZPlanesPlas] = { plasIn, plasIn};
      G4double rplasOut[numZPlanesPlas] = { plasOut, plasOut};

      //cell
      G4int nofCell=3;// кол-во ячеек в главной линии
      int nofLayers=3; //количество слоев ячеек
      int sideL=(nofCell+1)/2. ;
      int numZPlanesCell=2;
      G4double cellOut=wolfOut; G4double cellIn=0;
      G4double rcellIn[numZPlanesCell] = { cellIn, cellIn};
      G4double rcellOut[numZPlanesCell] = { cellOut, cellOut};
      G4double zPlaneCell[numZPlanesCell] = { 0,fullWolfLengthF};
      G4double cellR=23.1*mm;

      //silicon
      int numZPlanesSilic=6;
      G4double silicThick=0.1*cm;
      G4double silicRange=1*cm;
      G4double silicOut=plasOut*nofCell+silicRange;
      G4double silicIn=silicOut-silicThick;
      G4double silicLength=wolfLength;	//длина призмы
      G4double silicLengthF=2*mm;	//выступающий слой
      G4double fullSilicLength=silicLength*nofLayers+silicLengthF*nofLayers+(nofLayers-1)*silicLengthF;
      G4double fullSilicLengthF=silicLength*nofLayers+2*silicLengthF*nofLayers;
      G4double rsilicIn[numZPlanesSilic] = { 0,0,silicIn,silicIn,0,0};
      G4double rsilicOut[numZPlanesSilic] = { silicOut, silicOut,silicOut,silicOut, silicOut, silicOut};
      G4double zPlaneSilic[numZPlanesSilic] = { 0,silicLengthF,silicLengthF,
                                                fullSilicLength+2*silicRange,fullSilicLength+2*silicRange,fullSilicLengthF+2*silicRange};

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
      G4RotationMatrix* zRot = new G4RotationMatrix; // Rotates X and Z axes only
      zRot->rotateZ(M_PI/6.*rad);

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
        G4AssemblyVolume* assemblyDetector = new G4AssemblyVolume();
        // Rotation and translation of a plate inside the assembly
        G4RotationMatrix Ra;
        G4ThreeVector Ta, Ta1;
        G4Transform3D Tr, Tr1;
        // Rotation of the assembly inside the world
        G4RotationMatrix Rm, Rm1;

        // Fill the assembly by the plates
        assemblyDetector->AddPlacedVolume( cellLV, Tr );

        int check=0;
      // создание массива ячеек в форме правильной призмы
      ofstream file("position.dat");

      for( unsigned int j = 0; j < nofLayers; j++ )
      {
        for(  int k = -sideL+1; k < sideL; k++ )
        {
            for( unsigned int i = 0; i < nofCell-abs(k); i++ )
            {
            // Translation of the assembly inside the world
            G4ThreeVector Tm( i*(0 + 2*wolfOut) + abs(k)*wolfOut, k*3*wolfOut/sqrt(3),j*(0 + fullWolfLengthF));
            Tr = G4Transform3D(Rm,Tm);
            assemblyDetector->MakeImprint( worldLV, Tr );
            file << check <<'\t'<<setprecision(4)<< Tm <<endl;
            check++;
            }
        }
      }


      cout<<check<<endl;
      sizeDet=check;

      G4Polyhedra* Silic
        = new G4Polyhedra("Silic",            // its name
                     phiStart,phiTotal, numSide, numZPlanesSilic,  zPlaneSilic, rsilicIn, rsilicOut ); // its size


      G4LogicalVolume* silicLV
        = new G4LogicalVolume(
                     Silic,         // its solid
                     siliconMat,    // its material
                     "silicLV");          // its name


      G4VPhysicalVolume* silicPV
        = new G4PVPlacement(
                     zRot,                // no rotation
                     G4ThreeVector((nofCell/2.-1/2.)*plasOut*2,0, -silicRange), // its position
                     silicLV,       // its logical volume
                     "silicPV",           // its name
                     worldLV,          // its mother  volume
                     false,            // no boolean operation
                     0,fCheckOverlaps);             // copy number */





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
