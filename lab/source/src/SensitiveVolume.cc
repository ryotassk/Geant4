//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// SensitiveVolume.cc
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "SensitiveVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4HCofThisEvent.hh"
#include "Analysis.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


//------------------------------------------------------------------------------
  SensitiveVolume::SensitiveVolume(G4String name)
  : G4VSensitiveDetector(name)
{}
//------------------------------------------------------------------------------
  SensitiveVolume::~SensitiveVolume()
{}
//------------------------------------------------------------------------------
  void SensitiveVolume::Initialize(G4HCofThisEvent*)
{
     sum_eDep = 0.;
     //     sum_stepLength =0.;
     //   time = 0.;
}
//------------------------------------------------------------------------------
  void SensitiveVolume::EndOfEvent(G4HCofThisEvent*)
{
   //  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
   //   analysisManager->FillNtupleDColumn(0, sum_eDep);
   //   analysisManager->FillNtupleDColumn(1, time);
   //  analysisManager->AddNtupleRow();
/***
   G4cout <<  " eDep = "<< G4BestUnit(sum_eDep, "Energy")
          << " stepLength = " << G4BestUnit(sum_stepLength, "Length") 
          <<G4endl;
***/
}

//------------------------------------------------------------------------------
G4bool SensitiveVolume::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{

    
if(aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName()=="neutron"){
  G4StepPoint* point1 = aStep->GetPreStepPoint();
  if(point1->GetStepStatus()==fGeomBoundary){
  G4double kinE = point1->GetKineticEnergy();

  //  if(aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName()=="neutron"){return false;}
  // if(preStepPoint->GetStepStatus()==PhysVol_lab);

  //      G4double edep = aStep->GetTotalEnergyDeposit();
      //      G4double stepLength = aStep->GetNumberOfEvent();
   sum_eDep = kinE;
      //      sum_stepLength = sum_stepLength + stepLength;

  //  if(sum_eDep<kinE){
  //   sum_eDep = kinE;
  // }
   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
     analysisManager->FillNtupleDColumn(0, sum_eDep);
     analysisManager->AddNtupleRow();
  }
  }
 return true;
}






