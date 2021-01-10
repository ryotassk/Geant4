//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Geometry.cc
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "Geometry.hh"
#include "G4Sphere.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"

#include "SensitiveVolume.hh"
#include "G4SDManager.hh"

//------------------------------------------------------------------------------
  Geometry::Geometry(){}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
  Geometry::~Geometry() {}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
  G4VPhysicalVolume* Geometry::Construct()
//------------------------------------------------------------------------------
{
// Get pointer to 'Material Manager'
   G4NistManager* materi_Man = G4NistManager::Instance();

// Define 'World Volume'
   // Define the shape of solid
   G4double leng_X_World = 10.0*m;         // X-full-length of world 
   G4double leng_Y_World = 10.0*m;         // Y-full-length of world 
   G4double leng_Z_World = 10.0*m;         // Z-full-length of world 
   G4Box* solid_World = 
     new G4Box("Solid_World", leng_X_World/2.0, leng_Y_World/2.0, leng_Z_World/2.0);

   // Define logical volume
   G4Material* materi_World = materi_Man->FindOrBuildMaterial("G4_AIR");
   G4LogicalVolume* logVol_World = 
     new G4LogicalVolume(solid_World, materi_World, "LogVol_World");
   logVol_World->SetVisAttributes (G4VisAttributes::Invisible);
   
   // Placement of logical volume
   G4int copyNum_World = 0;               // Set ID number of world
   G4PVPlacement* physVol_World  = 
     new G4PVPlacement(G4Transform3D(), "PhysVol_World", logVol_World, 0, 
                       false, copyNum_World);
       
   G4double naihankei = 1.0*m;         // X-full-length of world
   G4double gaihankei = 3.0*m;         // Y-full-length of world
   G4double startPhi = 0*rad;         // Z-full-length of world
   G4double endPhi = 2*acos(-1)*rad;
   G4double startTheta = 0*rad;
   G4double endTheta = acos(-1)*rad;
   G4Sphere* solid_Rock = 
     new G4Sphere("solid_Rock", naihankei, gaihankei, startPhi, endPhi, startTheta, endTheta);

   G4String symbol;
   G4double a,z,density,fractionmass;
   G4int natoms,ncomponents;

   //   G4Isotope* he3 = new G4Isotope(name="he3", iz=2, in=3, a=3.0160293191*g/mole);
   //   G4Isotope* U5 = new G4Isotope("U235",iz=92,n=235,a=235.01*g/mole);
   //   G4Isotope* U8 = new G4Isotope("U238",iz=92,n=238,a=238.03*g/mole);
   
   //   G4Element* U = new G4Element("enriched Uranium",symbol="U",ncomponents=2);
   //   U->AddIsotope(U5,abundance = 90.*perCent);
   //   U->AddIsotope(U8,abundance = 10.*perCent);

   G4Element* O = new G4Element("Oxygen",symbol="O",z=8.,a=16.00*g/mole);
   G4Element* Si = new G4Element("Silicon",symbol="Si",z=14.,a=28.09*g/mole);
   G4Element* Al = new G4Element("Aluminium",symbol="Al",z=13.,a=26.98*g/mole);
   G4Element* Fe = new G4Element("Iron",symbol="Fe",z=26.,a=55.85*g/mole);
   G4Element* Mn = new G4Element("Manganese",symbol="Mn",z=25.,a=54.94*g/mole);
   G4Element* Mg = new G4Element("Magnesium",symbol="Mg",z=12.,a=24.31*g/mole);
   G4Element* Ca = new G4Element("Calcium",symbol="Ca",z=20.,a=40.08*g/mole);
   G4Element* Na = new G4Element("Sodium",symbol="Na",z=11.,a=22.99*g/mole);
   G4Element* P = new G4Element("Phosphorus",symbol="P",z=15.,a=30.97*g/mole);
   G4Element* S = new G4Element("Sulfur",symbol="S",z=16.,a=32.07*g/mole);
   G4Element* Zn = new G4Element("Zinc",symbol="Zn",z=30.,a=65.39*g/mole);
   G4Element* H = new G4Element("Hydrogen",symbol="H",z=1.,a=1.01*g/mole);
   //   G4Element* He3 = new G4Element("He3",symbol="He3",ncomponents=1);
   //   He3->AddIsotope(he3,abundance=100.*perCent);
/*
   G4Material* SiO2 =
     new G4Material("quartz",density=2.650*g/cm3,ncomponents=2);
   SiO2->AddElement(Si,natoms=1);
   SiO2->AddElement(O,natoms=2);

   G4Material* Al2O3 =
     new G4Material("alo",density=3.950*g/cm3,ncomponents=2);
   Al2O3->AddElement(Al,natoms=2);
   Al2O3->AddElement(O,natoms=3);

   G4Material* Fe2O3 =
     new G4Material("feo",density=5.240*g/cm3,ncomponents=2);
   Fe2O3->AddElement(Fe,natoms=2);
   Fe2O3->AddElement(O,natoms=3);

   G4Material* MnO =
     new G4Material("mno",density=5.370*g/cm3,ncomponents=2);
   MnO->AddElement(Mn,natoms=1);
   MnO->AddElement(O,natoms=1);

   G4Material* MgO =
     new G4Material("mgo",density=3.580*g/cm3,ncomponents=2);
   MgO->AddElement(Mg,natoms=1);
   MgO->AddElement(O,natoms=1);

   G4Material* CaO =
     new G4Material("cao",density=3.340*g/cm3,ncomponents=2);
   CaO->AddElement(Ca,natoms=1);
   CaO->AddElement(O,natoms=1);

   G4Material* Na2O =
     new G4Material("nao",density=2.270*g/cm3,ncomponents=2);
   Na2O->AddElement(Na,natoms=2);
   Na2O->AddElement(O,natoms=1);

   G4Material* P2O5 =
     new G4Material("po",density=2.390*g/cm3,ncomponents=2);
   P2O5->AddElement(P,natoms=2);
   P2O5->AddElement(O,natoms=5);

   G4Material* SO3 =
     new G4Material("so",density=1.920*g/cm3,ncomponents=2);
   SO3->AddElement(S,natoms=1);
   SO3->AddElement(O,natoms=3);

   G4Material* ZnO =
     new G4Material("zno",density=5.610*g/cm3,ncomponents=2);
   ZnO->AddElement(Zn,natoms=1);
   ZnO->AddElement(O,natoms=1);

   G4Material* H2 =
     new G4Material("h2",density=0.071*g/cm3,ncomponents=1);
   H2->AddElement(H,natoms=2); */

   //   G4Material* Herium3 = new G4Material("Herium3",density=1.34644166*mg/cm3,ncomponents=1,kStateGas,temperature=300.15*kelvin,pressure=10.*atmosphere);
   //   Herium3->AddElement(He3,1);

   G4Material* materi_Rock =
     new G4Material("materi_Rock",density = 3.0*g/cm3,ncomponents=11);
    
    materi_Rock->AddElement(Al, fractionmass=0.0601);
    materi_Rock->AddElement(Ca, fractionmass=0.2814);
    materi_Rock->AddElement(Fe, fractionmass=0.0766);
    materi_Rock->AddElement(Mg, fractionmass=0.0060);
    materi_Rock->AddElement(Mn, fractionmass=0.0085);
    materi_Rock->AddElement(Na, fractionmass=0.0001);
    materi_Rock->AddElement(O, fractionmass=0.3981);
    materi_Rock->AddElement(P, fractionmass=0.0015);
    materi_Rock->AddElement(S, fractionmass=0.0004);
    materi_Rock->AddElement(Si, fractionmass=0.1671);
    materi_Rock->AddElement(Zn, fractionmass=0.0002);

/*    materi_Rock->AddElement(Al, fractionmass=0.0582);
    materi_Rock->AddElement(Ca, fractionmass=0.2729);
    materi_Rock->AddElement(Fe, fractionmass=0.0744);
    materi_Rock->AddElement(H, fractionmass=0.0300);
    materi_Rock->AddElement(Mg, fractionmass=0.0058);
    materi_Rock->AddElement(Mn, fractionmass=0.0082);
    materi_Rock->AddElement(Na, fractionmass=0.0001);
    materi_Rock->AddElement(O, fractionmass=0.3862);
    materi_Rock->AddElement(P, fractionmass=0.0015);
    materi_Rock->AddElement(S, fractionmass=0.0004);
    materi_Rock->AddElement(Si, fractionmass=0.1621);
    materi_Rock->AddElement(Zn, fractionmass=0.0002); */
    
   /*   materi_Rock->AddMaterial(SiO2,fractionmass=0.3575);
   materi_Rock->AddMaterial(Al2O3,fractionmass=0.1135);
   materi_Rock->AddMaterial(Fe2O3,fractionmass=0.1095);
   materi_Rock->AddMaterial(MnO,fractionmass=0.0109);
   materi_Rock->AddMaterial(MgO,fractionmass=0.0099);
   materi_Rock->AddMaterial(CaO,fractionmass=0.3937);
   materi_Rock->AddMaterial(Na2O,fractionmass=0.0002);
   materi_Rock->AddMaterial(P2O5,fractionmass=0.0035);
   materi_Rock->AddMaterial(SO3,fractionmass=0.0010);
   materi_Rock->AddMaterial(ZnO,fractionmass=0.0003);

   materi_Rock->AddMaterial(SiO2,fractionmass=0.3468);
   materi_Rock->AddMaterial(Al2O3,fractionmass=0.1101);
   materi_Rock->AddMaterial(Fe2O3,fractionmass=0.1062);
   materi_Rock->AddMaterial(MnO,fractionmass=0.0105);
   materi_Rock->AddMaterial(MgO,fractionmass=0.0096);
   materi_Rock->AddMaterial(CaO,fractionmass=0.3819);
   materi_Rock->AddMaterial(Na2O,fractionmass=0.0002);
   materi_Rock->AddMaterial(P2O5,fractionmass=0.0034);
   materi_Rock->AddMaterial(SO3,fractionmass=0.0010);
   materi_Rock->AddMaterial(ZnO,fractionmass=0.0003);
   //  materi_Rock->AddElement(H,fractionmass=0.03);
   materi_Rock->AddMaterial(H2,fractionmass=0.03); */

   G4LogicalVolume* logVol_Rock = 
     new G4LogicalVolume( solid_Rock, materi_Rock, "LogVol_Rock", 0, 0, 0 );

   G4double pos_X_LogV = 0.0*cm;           // X-location LogV 
   G4double pos_Y_LogV = 0.0*cm;           // Y-location LogV
   G4double pos_Z_LogV = 0.0*cm;           // Z-location LogV
   G4ThreeVector threeVect_LogV = G4ThreeVector(pos_X_LogV, pos_Y_LogV, pos_Z_LogV);
   G4RotationMatrix rotMtrx_LogV = G4RotationMatrix();
   G4Transform3D trans3D_LogV = G4Transform3D(rotMtrx_LogV, threeVect_LogV);

   G4int copyNum_LogV = 1000;                // Set ID number of LogV
   new G4PVPlacement(trans3D_LogV, "PhysVol_Rock", logVol_Rock, physVol_World, 
                     false, copyNum_LogV);
   
// Define 'lab'
   // Define the shape of solid
   G4double naihankei2 = 0.0*m;         // X-full-length of world
   G4double gaihankei2 = 1.0*m;         // Y-full-length of world
   G4double startPhi2 = 0*rad;         // Z-full-length of world
   G4double endPhi2 = 2*acos(-1)*rad;
   G4double startTheta2 = 0*rad;
   G4double endTheta2 = acos(-1)*rad;
   G4Sphere* solid_lab = 
     new G4Sphere("solid_lab", naihankei2, gaihankei2, startPhi2, endPhi2, startTheta2, endTheta2);
   // Define logical volume
   G4Material* materi_lab = materi_Man->FindOrBuildMaterial("G4_AIR");
   G4LogicalVolume* logVol_lab = 
     new G4LogicalVolume( solid_lab, materi_lab, "LogVol_lab", 0, 0, 0 );
   
   // Placement of logical volume
   G4double labpos_X_LogV = 0.0*cm;           // X-location LogV 
   G4double labpos_Y_LogV = 0.0*cm;           // Y-location LogV
   G4double labpos_Z_LogV = 0.0*cm;           // Z-location LogV
   G4ThreeVector labthreeVect_LogV = G4ThreeVector(labpos_X_LogV, labpos_Y_LogV, labpos_Z_LogV);
   G4RotationMatrix labrotMtrx_LogV = G4RotationMatrix();
   G4Transform3D labtrans3D_LogV = G4Transform3D(labrotMtrx_LogV, labthreeVect_LogV);
   
   copyNum_LogV = 2000;                // Set ID number of LogV
   new G4PVPlacement(labtrans3D_LogV, "PhysVol_lab", logVol_lab, physVol_World, 
                     false,copyNum_LogV);

   /*
// Define 'PPC Detector'
   // Define the shape of solid
   G4double radius_PPC = 5.1*cm;
   G4double leng_Z_PPC = 43.2*cm;
   G4Tubs* solid_PPC = new G4Tubs("Solid_PPC", 0., radius_PPC, leng_Z_PPC, 0., 360.*deg); 

   // Define logical volume
   G4Material* materi_PPC = materi_Man->FindOrBuildMaterial("Herium3");
   G4LogicalVolume* logVol_PPC = 
     new G4LogicalVolume( solid_PPC, materi_PPC, "LogVol_PPC", 0, 0, 0 );
   
   // Placement of logical volume
   G4double PPCpos_X_LogV = 0.0*cm;           // X-location LogV 
   G4double PPCpos_Y_LogV = 0.0*cm;           // Y-location LogV
   G4double PPCpos_Z_LogV = 0.0*cm;           // Z-location LogV
   G4ThreeVector PPCthreeVect_LogV = G4ThreeVector(PPCpos_X_LogV, PPCpos_Y_LogV, PPCpos_Z_LogV);
   G4RotationMatrix PPCrotMtrx_LogV = G4RotationMatrix();
   G4Transform3D PPCtrans3D_LogV = G4Transform3D(PPCrotMtrx_LogV, PPCthreeVect_LogV);
   
//   copyNum_LogV = 3000;                // Set ID number of LogV
   new G4PVPlacement(PPCtrans3D_LogV, "PhysVol_PPC", logVol_PPC, physVol_World, 
                     false, copyNum_LogV);
   */
// Sensitive volume
    auto aSV = new SensitiveVolume("SensitiveVolume");
    logVol_lab->SetSensitiveDetector(aSV);         // Add sensitivity to the logical volume
    auto SDman = G4SDManager::GetSDMpointer();
    SDman->AddNewDetector(aSV);

// Return the physical world
   return physVol_World;
}
