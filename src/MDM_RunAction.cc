//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
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
#include <string.h>
#include <sstream>

//#include "boost/lexical_cast.hpp"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_RunAction::MDM_RunAction()
: G4UserRunAction(), fTimer(0)
{ 
  // add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);  

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
  analysisManager->CreateH1("h1_1","Edep in scintillator", 1000, 0., 1*MeV);
  analysisManager->CreateH1("h1_2","Photons hit in scintillator", 10000, 0., 10000);
  analysisManager->CreateH1("h1_3","trackL in scintillator", 1000, 0., 0.5*cm);
  analysisManager->CreateH1("h1_4","Photons hit in sipm", 10000, 0., 10000);
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
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void MDM_RunAction::EndOfRunAction(const G4Run* run)
{

	fTimer->Stop();

  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;
  

  /*
  const MDM_Run* MDM_Run = static_cast<const MDM_Run*>(run);

  // Compute dose
  //
  G4double edep  = MDM_Run->GetEdep();
  G4double edep2 = MDM_Run->GetEdep2();
  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;

  const MDM_DetectorConstruction* detectorConstruction
   = static_cast<const MDM_DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
  G4double dose = edep/mass;
  G4double rmsDose = rms/mass;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const MDM_PrimaryGeneratorAction* generatorAction
   = static_cast<const MDM_PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }
          
  // Print
  //  
  if (IsMaster()) {
    G4cout
     << "\n--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << "\n--------------------End of Local Run------------------------";
  }
  
  G4cout
     << "\n The run consists of " << nofEvents << " "<< runCondition
     << "\n Dose in scoring volume : " 
     << G4BestUnit(dose,"Dose") << " +- " << G4BestUnit(rmsDose,"Dose")
     << "\n------------------------------------------------------------\n"
     << G4endl;

	  */


    // print histogram statistics
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->GetH1(1) ) {
    G4cout << "\n ----> print histograms statistic ";
    if(isMaster) {
      G4cout << "for the entire run \n" << G4endl; 
    }
    else {
      G4cout << "for the local thread \n" << G4endl; 
    }
    
    G4cout << " Escint : mean = " 
       << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy") 
       << " rms = " 
       << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;
    
    G4cout << " Esipm : mean = " 
       << G4BestUnit(analysisManager->GetH1(2)->mean(), "Energy") 
       << " rms = " 
       << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Energy") << G4endl;
    
    G4cout << " Lscint : mean = " 
      << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length") 
      << " rms = " 
      << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Length") << G4endl;

    G4cout << " Lsipm : mean = " 
      << G4BestUnit(analysisManager->GetH1(4)->mean(), "Length") 
      << " rms = " 
      << G4BestUnit(analysisManager->GetH1(4)->rms(),  "Length") << G4endl;
  }

  G4cout<<"number of even = "<<run->GetNumberOfEvent()<<" "<< *fTimer << G4endl;

// Save histograms

analysisManager->Write();
analysisManager->CloseFile();


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
