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
    if(neutron_energy >= 0.000000001){
    
        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
        analysisManager->FillNtupleDColumn(0, neutron_energy);
        analysisManager->AddNtupleRow();
    
   }
}

//------------------------------------------------------------------------------

G4bool SensitiveVolume::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    G4double pre_kinE=0,post_kinE=0;
    
    if(aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName()=="neutron"){
        G4StepPoint* point1 = aStep->GetPreStepPoint();
        G4StepPoint* point2 = aStep->GetPostStepPoint();
        if(point1->GetStepStatus()==fGeomBoundary){
            pre_kinE = point1->GetKineticEnergy();
        }
        if(point2->GetStepStatus()==fGeomBoundary){
            post_kinE = point2->GetKineticEnergy();
        }
        if(pre_kinE != post_kinE){
            neutron_energy = pre_kinE ;
        }
      
    }
  
    return true;
}






