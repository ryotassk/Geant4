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
  Geometry::Geometry()
//: G4VUserDetectorConstruction(),
//  fScoringVolume(0)
{ }
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
   G4double leng_X_World = 100.0*m;         // X-full-length of world
   G4double leng_Y_World = 100.0*m;         // Y-full-length of world
   G4double leng_Z_World = 100.0*m;         // Z-full-length of world 
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

   G4String symbol;
   G4double a,z,n,density,fractionmass,temperature,abundance,pressure;
   G4int natoms,ncomponents;

    G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
   G4Isotope* he3 = new G4Isotope("he3", z=2, n=3, a=3.0160293191*g/mole);
     G4Element* B  = new G4Element("Boron",symbol="B" , z= 5., a= 10.81*g/mole);
    G4Element* C = new G4Element("Carbon",symbol="C",z=6.,a=12.01*g/mole);
    G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
    G4Element* Si = new G4Element("Silicon",symbol="Si",z=14.,a=28.09*g/mole);
    G4Element* P = new G4Element("Phosphorus",symbol="P",z=15.,a=30.97*g/mole);
    G4Element* S = new G4Element("Sulfur",symbol="S",z=16.,a=32.07*g/mole);
    G4Element* Cr = new G4Element("Chromium",symbol="Cr",z=24.,a=52.00*g/mole);
    G4Element* Mn = new G4Element("Manganese",symbol="Mn",z=25.,a=54.94*g/mole);
    G4Element* Fe = new G4Element("Iron",symbol="Fe",z=26.,a=55.85*g/mole);
    G4Element* Ni = new G4Element("Nickel",symbol="Ni",z=28.,a=58.69*g/mole);
    
   G4Element* He3 = new G4Element("He3",symbol="He3",ncomponents=1);
   He3->AddIsotope(he3,abundance=100.*perCent);

   G4Material* Herium3 = new G4Material("Helium3",density=1.346*mg/cm3,ncomponents=1,kStateGas,temperature=300.15*kelvin,pressure=10.*atmosphere);
   Herium3->AddElement(He3,1);
    
   
   G4Material* SUS_304 =
     new G4Material("SUS_304",density = 8.03*g/cm3,ncomponents=8);
    SUS_304->AddElement(C,fractionmass=0.0008);
    SUS_304->AddElement(Si,fractionmass=0.0010);
    SUS_304->AddElement(Mn,fractionmass=0.0020);
    SUS_304->AddElement(P,fractionmass=0.0004);
    SUS_304->AddElement(S,fractionmass=0.0003);
    SUS_304->AddElement(Ni,fractionmass=0.1050);
    SUS_304->AddElement(Cr,fractionmass=0.2000);
    SUS_304->AddElement(Fe,fractionmass=0.6905);
    
    G4Material* polyethylene=new G4Material("polyethylene",density=0.92*g/cm3,ncomponents=2);
    polyethylene->AddElement(C,1);
    polyethylene->AddElement(H,2);
    
    const G4double R_B4C  = 40;    // target value, wt-%
    const G4double R0_B4C = 20;    // original value, wt-%
    const G4double W0_H   = 0.065; // fraction Mass
    const G4double W0_O   = 0.173;
    const G4double W0_Si  = 0.303;
    const G4double W0_B   = 0.157;
    const G4double W0_C   = 0.302;
    
    G4double totalMass = R_B4C / R0_B4C * (W0_B + W0_C) + (100.0 - R_B4C) / (100.0 - R0_B4C) * (W0_H + W0_O + W0_Si);
    G4double W_H       = (100.0 - R_B4C) / (100.0 - R0_B4C) * W0_H  / totalMass;
    G4double W_O       = (100.0 - R_B4C) / (100.0 - R0_B4C) * W0_O  / totalMass;
    G4double W_Si      = (100.0 - R_B4C) / (100.0 - R0_B4C) * W0_Si / totalMass;
    G4double W_B       = R_B4C / R0_B4C * W0_B / totalMass;
    G4double W_C       = R_B4C / R0_B4C * W0_C / totalMass;
    
    // ref. density (2015/7/30 from T.Iida)
    G4Material* RubberB4C = new G4Material("RubberB4C", density = 1.42 *g/cm3, ncomponents = 5);
    RubberB4C -> AddElement(H,  fractionmass = W_H);
    RubberB4C -> AddElement(O,  fractionmass = W_O);
    RubberB4C -> AddElement(Si, fractionmass = W_Si);
    RubberB4C -> AddElement(B,  fractionmass = W_B);
    RubberB4C -> AddElement(C,  fractionmass = W_C);
    
/*    // Define 'lab'
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
    
    G4int copyNuml_LogV = 0;                // Set ID number of LogV
    new G4PVPlacement(labtrans3D_LogV, "PhysVol_lab", logVol_lab, physVol_World,
                      false,copyNuml_LogV);

 */
    // Define 'PPC case'
    // Define the shape of solid
      G4double inner_radius_case = 2.49*cm;
      G4double outer_radius_case = 2.54*cm;
      G4double leng_Z_case = 21.555*cm;
      G4Tubs* solid_case = new G4Tubs("solid_case", inner_radius_case, outer_radius_case, leng_Z_case, 0., 360.*deg);

     G4LogicalVolume* logVol_case =
       new G4LogicalVolume( solid_case, SUS_304, "logVol_case", 0, 0, 0 );

     G4double pos_X_LogV = 0.0*cm;           // X-location LogV
     G4double pos_Y_LogV = 0.0*cm;           // Y-location LogV
     G4double pos_Z_LogV = 0.0*cm;           // Z-location LogV
     G4ThreeVector threeVect_LogV = G4ThreeVector(pos_X_LogV, pos_Y_LogV, pos_Z_LogV);
     G4RotationMatrix rotMtrx_LogV = G4RotationMatrix();
     G4Transform3D trans3D_LogV = G4Transform3D(rotMtrx_LogV, threeVect_LogV);

    G4int copyNum_LogV = 0;                // Set ID number of LogV
    new G4PVPlacement(trans3D_LogV, "PhysVol_case", logVol_case, physVol_World,
                      false, copyNum_LogV);
    
     G4double radius_casefuta = 2.54*cm;
     G4double leng_Z_casefuta = 0.025*cm;
     G4Tubs* solid_casefuta = new G4Tubs("solid_casefuta", 0, radius_casefuta, leng_Z_casefuta, 0., 360.*deg);

    G4LogicalVolume* logVol_casefuta =
      new G4LogicalVolume( solid_casefuta, SUS_304, "logVol_casefuta", 0, 0, 0 );

    G4double posf_Z_LogV = 21.58*cm;           // Z-location LogV
    G4ThreeVector threeVectf_LogV = G4ThreeVector(pos_X_LogV, pos_Y_LogV, posf_Z_LogV);
    G4Transform3D trans3Df_LogV = G4Transform3D(rotMtrx_LogV, threeVectf_LogV);

    G4int copyNumf_LogV = 1000;                // Set ID number of LogV
    new G4PVPlacement(trans3Df_LogV, "PhysVol_casefuta", logVol_casefuta, physVol_World,
                      false, copyNumf_LogV);
    
    posf_Z_LogV = -21.58*cm;           // Z-location LogV
    threeVectf_LogV = G4ThreeVector(pos_X_LogV, pos_Y_LogV, posf_Z_LogV);
    trans3Df_LogV = G4Transform3D(rotMtrx_LogV, threeVectf_LogV);

    copyNumf_LogV = 2000;                // Set ID number of LogV
    new G4PVPlacement(trans3Df_LogV, "PhysVol_casefuta", logVol_casefuta, physVol_World,
                      false, copyNumf_LogV);

   
   
// Define 'PPC gas'
   // Define the shape of solid
   G4double radius_gas = 2.49*cm;
   G4double leng_Z_gas = 19.05*cm;
   G4Tubs* solid_gas = new G4Tubs("solid_gas", 0., radius_gas, leng_Z_gas, 0., 360.*deg);

   // Define logical volume
   G4Material* materi_gas = materi_Man->FindOrBuildMaterial("Helium3");
   G4LogicalVolume* logVol_gas =
     new G4LogicalVolume( solid_gas, materi_gas, "logVol_gas", 0, 0, 0 );
   
    G4int copyNumg_LogV = 0;                // Set ID number of LogV
   new G4PVPlacement(trans3D_LogV, "PhysVol_gas", logVol_gas, physVol_World,
                     false, copyNumg_LogV);
/*
    // Define 'polyethylene shield'
    // Define the shape of solid
      G4double inner_radius_p = 2.565*cm;
      G4double outer_radius_p = 9.065*cm;
      G4double leng_Z_p = 23.75*cm;
      G4Tubs* solid_p = new G4Tubs("solid_p", inner_radius_p, outer_radius_p, leng_Z_p, 0., 360.*deg);

     G4LogicalVolume* logVol_p =
       new G4LogicalVolume( solid_p, polyethylene, "logVol_p", 0, 0, 0 );
    
    G4int copyNump_LogV = 0;                // Set ID number of LogV
    new G4PVPlacement(trans3D_LogV, "PhysVol_p", logVol_p, physVol_World,
                      false, copyNump_LogV);
    
     G4double radius_pfuta = 9.065*cm;
     G4double leng_Z_pfuta = 3.25*cm;
     G4Tubs* solid_pfuta = new G4Tubs("solid_pfuta", 0, radius_pfuta, leng_Z_pfuta, 0., 360.*deg);

    G4LogicalVolume* logVol_pfuta =
      new G4LogicalVolume( solid_pfuta, polyethylene, "logVol_pfuta", 0, 0, 0 );

    G4double pospf_Z_LogV = 27.00*cm;           // Z-location LogV
    G4ThreeVector threeVectpf_LogV = G4ThreeVector(pos_X_LogV, pos_Y_LogV, pospf_Z_LogV);
    G4Transform3D trans3Dpf_LogV = G4Transform3D(rotMtrx_LogV, threeVectpf_LogV);

    G4int copyNumpf_LogV = 1000;                // Set ID number of LogV
    new G4PVPlacement(trans3Dpf_LogV, "PhysVol_pfuta", logVol_pfuta, physVol_World,
                      false, copyNumpf_LogV);
    
    pospf_Z_LogV = -27.00*cm;           // Z-location LogV
    threeVectpf_LogV = G4ThreeVector(pos_X_LogV, pos_Y_LogV, pospf_Z_LogV);
    trans3Dpf_LogV = G4Transform3D(rotMtrx_LogV, threeVectpf_LogV);

    copyNumpf_LogV = 2000;                // Set ID number of LogV
    new G4PVPlacement(trans3Dpf_LogV, "PhysVol_pfuta", logVol_pfuta, physVol_World,
                      false, copyNumpf_LogV);

    // Define 'B4C'
    // Define the shape of solid
      G4double inner_radius_b = 9.065*cm;
      G4double outer_radius_b = 9.565*cm;
      G4double leng_Z_b = 30.25*cm;
      G4Tubs* solid_b = new G4Tubs("solid_b", inner_radius_b, outer_radius_b, leng_Z_b, 0., 360.*deg);

     G4LogicalVolume* logVol_b =
       new G4LogicalVolume( solid_b, RubberB4C, "logVol_b", 0, 0, 0 );
    
    G4int copyNumb_LogV = 0;                // Set ID number of LogV
    new G4PVPlacement(trans3D_LogV, "PhysVol_b", logVol_b, physVol_World,
                      false, copyNumb_LogV);
    
     G4double radius_bfuta = 9.565*cm;
     G4double leng_Z_bfuta = 0.25*cm;
     G4Tubs* solid_bfuta = new G4Tubs("solid_bfuta", 0, radius_bfuta, leng_Z_bfuta, 0., 360.*deg);

    G4LogicalVolume* logVol_bfuta =
      new G4LogicalVolume( solid_bfuta, RubberB4C, "logVol_bfuta", 0, 0, 0 );

    G4double posbf_Z_LogV = 30.50*cm;           // Z-location LogV
    G4ThreeVector threeVectbf_LogV = G4ThreeVector(pos_X_LogV, pos_Y_LogV, posbf_Z_LogV);
    G4Transform3D trans3Dbf_LogV = G4Transform3D(rotMtrx_LogV, threeVectbf_LogV);

    G4int copyNumbf_LogV = 1000;                // Set ID number of LogV
    new G4PVPlacement(trans3Dbf_LogV, "PhysVol_bfuta", logVol_bfuta, physVol_World,
                      false, copyNumbf_LogV);

    pospf_Z_LogV = -30.50*cm;           // Z-location LogV
    threeVectbf_LogV = G4ThreeVector(pos_X_LogV, pos_Y_LogV, posbf_Z_LogV);
    trans3Dbf_LogV = G4Transform3D(rotMtrx_LogV, threeVectbf_LogV);

    copyNumbf_LogV = 2000;                // Set ID number of LogV
    new G4PVPlacement(trans3Dbf_LogV, "PhysVol_bfuta", logVol_bfuta, physVol_World,
                      false, copyNumbf_LogV);
*/
// Sensitive volume
    auto aSV = new SensitiveVolume("SensitiveVolume");
    logVol_gas->SetSensitiveDetector(aSV);         // Add sensitivity to the logical volume
    auto SDman = G4SDManager::GetSDMpointer();
    SDman->AddNewDetector(aSV);
    
//    fScoringVolume = logVol_bfuta;
    
// Return the physical world
   return physVol_World;
}
