//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
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

#include "MDM_RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Timer.hh"

#include "g4root.hh"

#include "persistency/PersistencyHandler.hh"
#include <G4DigiManager.hh>

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
	
	G4cout<<"write to the root "<<G4endl;

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

	// commented out on 26 Aug 2015 for reducing the output message
	/*
	for(int i=0;i<n_scint_hit;i++)
	{
		G4cout << "----Hit #" <<i << G4endl;
		(*scintHC)[i]->Print();
	}
	
	G4cout << " Total Hits in the SiPM: " << n_sipm_hit << G4endl;

	for(int i=0;i<n_sipm_hit;i++)
	{
		G4cout << "----Hit #" <<i << G4endl;
		(*sipmHC)[i]->Print();
	}
	*/


  // get analysis manager

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //G4cout <<"debug1" << G4endl;

  // fill histograms
  if(n_scint_hit>0)
  {
	  analysisManager->FillH1(1, scintHit->GetTotalEdep());
	  analysisManager->FillH1(2, G4double(n_scint_hit));
	  analysisManager->FillH1(3, scintHit->GetTotalTrackLength());
	  analysisManager->FillH1(4, G4double(n_sipm_hit));  // total photons detected per event
	  
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
		sipmHit = (*sipmHC)[i];
		G4ThreeVector pos = sipmHit->GetPos();
		analysisManager->FillH2(1,pos.getX(),pos.getY(),1.0);
	}

  }
  
	  

  PrintEventStatistics(scintHit->GetTotalEdep(),scintHit->GetTotalTrackLength(),sipmHit->GetPhotons(),sipmHit->GetTrackID());

  //G4Sipm sample code
	PersistencyHandler* persistency =
			((MDM_RunAction*) G4RunManager::GetRunManager()->GetUserRunAction())->getPersistencyHandler();
	// Run all digitizer modules.
	G4DigiManager* digiManager = G4DigiManager::GetDMpointer();
	G4DCtable* dcTable = digiManager->GetDCtable();
	G4cout<< "dcTable numbers "<<dcTable->entries()<<G4endl;
	for (int i = 0; i < dcTable->entries(); i++) {
		G4String dmName = dcTable->GetDMname(i);
		G4VDigitizerModule* dm = digiManager->FindDigitizerModule(dmName);
		if (dm) {
			dm->Digitize();
		}
	}
	G4Timer timer;
	timer.Start();
	// Process hits collections.

	G4HCofThisEvent* hCof = event->GetHCofThisEvent();
	if (hCof != NULL) {
		for (int i = 0; i < hCof->GetCapacity(); ++i) {
			G4VHitsCollection* hc = hCof->GetHC(i);
			if (hc != NULL) {
				if (dynamic_cast<G4SipmHitsCollection*>(hc)) {
					persistency->persist((G4SipmHitsCollection*) hc);
				}
			}
		}
	}

	// Process digi collections.
	G4DCofThisEvent* dCof = event->GetDCofThisEvent();
	if (dCof != NULL) {
		for (int i = 0; i < dCof->GetCapacity(); ++i) {
			G4VDigiCollection* dc = dCof->GetDC(i);
			if (dc != NULL) {
				if (dynamic_cast<G4SipmDigiCollection*>(dc)) {
					persistency->persist((G4SipmDigiCollection*) dc);
				}
				if (dynamic_cast<G4SipmVoltageTraceDigiCollection*>(dc)) {
					persistency->persist((G4SipmVoltageTraceDigiCollection*) dc);
				}
			}
		}
	}

	timer.Stop();
	std::cout << "EventAction::EndOfEventAction(): persist time (" << timer << ")." << std::endl;


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
