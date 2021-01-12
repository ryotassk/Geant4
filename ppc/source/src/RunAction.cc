#include "RunAction.hh"
#include "Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

RunAction::RunAction()
 : G4UserRunAction()
{
    // Get analysis manager
    auto analysisManager = G4AnalysisManager::Instance();
    //G4cout << "Using " << analysisManager->GetType() << G4endl;
    analysisManager->SetActivation(true);

     analysisManager->SetVerboseLevel(1);
     analysisManager->SetFileName("");
    
	// creating a Ntuple for all

    analysisManager-> CreateNtuple( "Flux",  "Energy");
    analysisManager->CreateNtupleDColumn( "Energy");

    analysisManager->FinishNtuple();

  }

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();  
}

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  analysisManager->OpenFile();
    //G4cout << "File " << fileName << " Open" <<G4endl;
}

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // save the file and close
  analysisManager->Write();
  analysisManager->CloseFile();
}
