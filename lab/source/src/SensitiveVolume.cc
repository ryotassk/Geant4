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
     neutron_energy = 0.;
}
//------------------------------------------------------------------------------
  void SensitiveVolume::EndOfEvent(G4HCofThisEvent*)
{

}

//------------------------------------------------------------------------------
G4bool SensitiveVolume::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    if(aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName()=="neutron"){
        G4StepPoint* point1 = aStep->GetPreStepPoint();
        if(point1->GetStepStatus()==fGeomBoundary){
            G4double kinE = point1->GetKineticEnergy();
            neutron_energy = kinE;
      
            G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
            analysisManager->FillNtupleDColumn(0, neutron_energy);
            analysisManager->AddNtupleRow();
        }
    }
 return true;
}






