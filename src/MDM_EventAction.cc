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
// $Id: MDM_EventAction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file MDM_EventAction.cc
/// \brief Implementation of the MDM_EventAction class

#include "MDM_EventAction.hh"
#include "MDM_CalorimeterSD.hh"
#include "MDM_CalorHit.hh"
#include "MDM_SipmHit.hh"
#include "MDM_SipmSD.hh"

#include "MDM_Run.hh"



#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "g4root.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_EventAction::MDM_EventAction()
: G4UserEventAction(),
  fEdep(0.),
  fEnergyAbs(0.),
  fEnergyGap(0.),
  fTrackLAbs(0.),
  fTrackLGap(0.),
  fScintHCID(-1),
  fSipmHCID(-1)

{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_EventAction::~MDM_EventAction()
{	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_EventAction::BeginOfEventAction(const G4Event* event)
{    
  fEdep = 0.;
    // initialisation per event
  fEnergyAbs = 0.;
  fEnergyGap = 0.;
  fTrackLAbs = 0.;
  fTrackLGap = 0.;

  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  //getcollectionID(() method is a heavy operation, it should not be invoked for every event
  if(fScintHCID<0)
  {
	  fScintHCID = SDman->GetCollectionID("scintHitsCollection");
	  fSipmHCID = SDman->GetCollectionID("sipmHitsCollection");
	  G4cout << "Collection ID prints - scintHCID: " << fScintHCID << G4endl;
	  G4cout << "Collection ID prints - sipmHCID: " << fSipmHCID << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_EventAction::EndOfEventAction(const G4Event* event)
{   
	/*
  // Get hits collection IDs (only once)
	if( fScintHCID == -1)
	{
		fScintHCID = G4SDManager::GetSDMpointer()->GetCollectionID("scintHitsCollection");
		fSipmHCID = G4SDManager::GetSDMpointer()->GetCollectionID("sipmHitsCollection");
	}
	*/

	if(fScintHCID<0 || fSipmHCID<0) return;
	

	// get hits collections
	MDM_CalorHitsCollection* scintHC = GetHitsCollection(fScintHCID,event);
	MDM_SipmHitsCollection* sipmHC = GetHitsCollection_2(fSipmHCID,event);

	// get hit with total value
	G4int n_scint_hit=0;
	G4int n_sipm_hit=0;

	if(scintHC!=NULL)
		n_scint_hit = scintHC->entries();

	if(sipmHC!=NULL)
		n_sipm_hit = sipmHC->entries();

	MDM_CalorHit* scintHit;
	MDM_SipmHit* sipmHit;

	if(n_scint_hit>0)
		scintHit = (*scintHC)[n_scint_hit-1];


	if(n_sipm_hit>0)
		sipmHit = (*sipmHC)[n_sipm_hit-1];

	// Print per event (module n)
	
	G4cout << " Total Hits in the scintillator: " << n_scint_hit << G4endl;

	/*
	for(int i=0;i<n_scint_hit;i++)
	{
		G4cout << "----Hit #" <<i << G4endl;
		(*scintHC)[i]->Print();
	}
	*/
	G4cout << " Total Hits in the SiPM: " << n_sipm_hit << G4endl;

	/*
	for(int i=0;i<n_sipm_hit;i++)
	{
		G4cout << "----Hit #" <<i << G4endl;
		(*sipmHC)[i]->Print();
	}
	*/



	/*
	G4int eventID = event->GetEventID();
	G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
	if ( (printModulo > 0) && ( eventID % printModulo == 0)) {
		G4cout << "----> End of event: " << eventID << G4endl;
		PrintEventStatistics(scintHit->GetEdep(),scintHit->GetTrackLength(), sipmHit->GetEdep(), sipmHit->GetTrackLength());
	}
	*/

 
  // accumulate statistics in MDM_Run
 // MDM_Run* run = static_cast<MDM_Run*>(
  //      G4RunManager::GetRunManager()->GetNonConstCurrentRun());
 // run->AddEdep(fEdep);


  // get analysis manager
	// 2015-5-31 temperaly remove the root file generator to speed up the simulation
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //G4cout <<"debug1" << G4endl;

  // fill histograms
  if(n_scint_hit>0)
  {
	  analysisManager->FillH1(1, scintHit->GetTotalEdep());
	 // G4cout <<"debug2" << G4endl;
	  analysisManager->FillH1(2, G4double(n_scint_hit));
	 // G4cout <<"debug3" << G4endl;
	  analysisManager->FillH1(3, scintHit->GetTotalTrackLength());
	 // G4cout <<"debug4" << G4endl;
	  analysisManager->FillH1(4, G4double(n_sipm_hit));  // total photons per event
	  
  }else
  {
	  analysisManager->FillH1(1, 0.0);
	  analysisManager->FillH1(2, G4double(n_scint_hit));
	  analysisManager->FillH1(3, 0.0);
	  analysisManager->FillH1(4, G4double(n_sipm_hit));  // total photons per event
  }

  
  if(n_sipm_hit>0)
  {
	for(int i=0;i<n_sipm_hit;i++)
	{
		//G4cout << "----Hit #" <<i << G4endl;
		//(*sipmHC)[i]->Print();
		sipmHit = (*sipmHC)[i];
		G4ThreeVector pos = sipmHit->GetPos();
		analysisManager->FillH2(1,pos.x(),pos.y(),1.0);
	}
	//G4cout <<"debug5" << G4endl;
  }
  
	  

  // fill ntuple
 // analysisManager->FillNtupleDColumn(0, fEnergyAbs);
 // analysisManager->FillNtupleDColumn(1, fEnergyGap);
 // analysisManager->FillNtupleDColumn(2, fTrackLAbs);
 // analysisManager->FillNtupleDColumn(3, fTrackLGap);
 // analysisManager->AddNtupleRow();  

  PrintEventStatistics(scintHit->GetTotalEdep(),scintHit->GetTotalTrackLength(),sipmHit->GetPhotons(),sipmHit->GetTrackID());


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_CalorHitsCollection* MDM_EventAction::GetHitsCollection(G4int hcID, const G4Event* event) const
{
	MDM_CalorHitsCollection* hitsCollection = static_cast<MDM_CalorHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));

	if(!hitsCollection){
		G4ExceptionDescription msg;
		msg<< "Canot access hitsCollection ID" << hcID;
		G4Exception("MDM_EventAction::GetHitsCollection()","MyCode0003", FatalException, msg);
	}

	return hitsCollection;
}

MDM_SipmHitsCollection* MDM_EventAction::GetHitsCollection_2(G4int hcID, const G4Event* event) const
{
	MDM_SipmHitsCollection* hitsCollection = static_cast<MDM_SipmHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));

	if(!hitsCollection){
		G4ExceptionDescription msg;
		msg<< "Canot access hitsCollection ID" << hcID;
		G4Exception("MDM_EventAction::GetHitsCollection()","MyCode0004", FatalException, msg);
	}

	return hitsCollection;
}

void MDM_EventAction::PrintEventStatistics(G4double scintEdep, G4double scintTrackLength, G4double sipmPhotons, G4double sipmTrackID) const
{
	// print event statistics
	G4cout
		<<"   Scintillator: total energy: "
		<<std::setw(7) << G4BestUnit(scintEdep, "Energy")
		<<"           total track length: "
		<<std::setw(7) << G4BestUnit(scintTrackLength, "Length")
		<<G4endl
		<<"           SiPM: total photons: "
		<<std::setw(7) << sipmPhotons
		<<"           Track ID: "
		<<std::setw(7) << sipmTrackID
		<<G4endl;
}