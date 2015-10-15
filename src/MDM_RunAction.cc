//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// ********************************************************************
//
// $Id: MDM_RunAction.cc 75216 2013-10-29 16:08:11Z gcosmo $
//
/// \file MDM_RunAction.cc
/// \brief Implementation of the MDM_RunAction class

#include "MDM_RunAction.hh"
#include "MDM_PrimaryGeneratorAction.hh"
#include "MDM_DetectorConstruction.hh"
#include "MDM_Run.hh"
#include "G4Timer.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "g4root.hh"
#include <string>
//#include <string.h>
#include <sstream>

// from sample code
#include "G4SipmUiMessenger.hh"
//#include "ParticleSourceMessenger.hh"
//#include "DetectorConstruction.hh"
#include "persistency/PersistencyHandler.hh"
#include "persistency/PersistVisitorFactory.hh"


//#include "boost/lexical_cast.hpp"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_RunAction::MDM_RunAction(std::string _filename)
: G4UserRunAction(), fTimer(0), filename(_filename)
{ 

	//from sample code
persistencyHandler = new PersistencyHandler(PersistVisitorFactory::getInstance()->create(filename));

  fTimer = new G4Timer;

  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in B4Analysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(0);
  analysisManager->SetFirstHistoId(1);

  // Book histograms, ntuple
  //
  
  // Creating histograms
  analysisManager->CreateH1("h1_1","Edep in scintillator", 1000, 0.0, 10.0*eV,"eV");
  analysisManager->CreateH1("h1_2","Photons hit in scintillator", 1000, 0., 1000);
  analysisManager->CreateH1("h1_3","trackL in scintillator", 1000, 0.0, 0.5*cm,"cm");
  analysisManager->CreateH1("h1_4","Photons hit in sipm", 1000, 0., 100);
  analysisManager->CreateH2("h2_1","Photons hit position in sipm",1000,-1*mm,1*mm,1000,-1*mm,1*mm,"mm","mm","x_pos","y_pos");


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_RunAction::~MDM_RunAction()
{
	delete fTimer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* MDM_RunAction::GenerateRun()
{
  return new MDM_Run; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  //inform the runManager to save random number seed
	G4RunManager::GetRunManager()->SetRandomNumberStore(false);

	  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //  
  // what a pain in ass to get the G4int to G4string done !! there must be a better way1!

  int run_ID = (int)(aRun->GetRunID());
  G4String fileName = "MDM__run";
  std::ostringstream temp;
  temp << run_ID;
  fileName += (G4String)temp.str();

   analysisManager->OpenFile(fileName);
 
	G4cout <<"### Run"<< run_ID << " start."<< G4endl;
	fTimer->Start();
	
	//from sample code
	persistencyHandler->open(filename);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void MDM_RunAction::EndOfRunAction(const G4Run* run)
{

	fTimer->Stop();

  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;
  

          
  // print histogram statistics
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4cout<<"number of even = "<<nofEvents<<" "<< *fTimer << G4endl;

// Save histograms

analysisManager->Write();
analysisManager->CloseFile();

//from sample code
// Persist run settings.
//persistencyHandler->persist(G4SipmUiMessenger::getInstance());
//persistencyHandler->persist(ParticleSourceMessenger::getInstance());
// Get detector.
const MDM_DetectorConstruction* detector =
		(const MDM_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
persistencyHandler->persist(detector->getSipmHousing()->getSipm()->getModel());
persistencyHandler->persist(detector->getSipmHousing()->getSipm()->getModel()->getVoltageTraceModel());
// Close output.
persistencyHandler->close();


}

std::string MDM_RunAction::getFilename() const {
	return filename;
}

PersistencyHandler* MDM_RunAction::getPersistencyHandler() const {
	return persistencyHandler;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
