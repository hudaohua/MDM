//
// ********************************************************************
// * License and Disclaimer                                           *
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *
//
// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// ********************************************************************
//
// $Id: MDM_SipmSD.cc 75215 2013-10-29 16:07:06Z gcosmo $
//
/// \file MDM_SipmSD.cc
/// \brief Implementation of the MDM_SipmSD class

#include "MDM_SipmSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include <G4OpticalPhoton.hh>




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_SipmSD::MDM_SipmSD(
                            const G4String& name, 
                            const G4String& hitsCollectionName,
                            G4int nofCells)
 : G4VSensitiveDetector(name),
   fSipmHitsCollection(0),
   fNofCells(nofCells)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_SipmSD::~MDM_SipmSD() 
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_SipmSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fSipmHitsCollection 
    = new MDM_SipmHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce
  static G4int hcID1 = -1;
  if(hcID1<0)
{ hcID1 = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  G4cout << "SiPm hcID1: " << hcID1 << "- size "<< collectionName.size() <<G4endl;

	// hcID1 = GetCollectionID(0);
  }
  hce->AddHitsCollection( hcID1, fSipmHitsCollection ); 

  // Create hits
  // fNofCells for cells + one more for total sums 
 // for (G4int i=0; i<fNofCells+1; i++ ) {
 //   fHitsCollection->insert(new MDM_SipmHit());
  //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * The data processing depends on the hitting particle:
 * - mark that the sensitive detector has been hit
 * - <b> if the particle is an optical photons AND has been created outside the sensitive detector: </b>
 *   - save the time and position of the sensitive detector's hit as well as the particle's momentum when hitting the sensitive detector and the name of the sensitive detector's volume that was hit (PhotonDetectorDataStorage)
 *   - write the data needed to resimulate a photon that hit the sensitive detector (start time, start position, start momentum, polarisation) into a text file
 * - stop the particle
 */

G4bool MDM_SipmSD::ProcessHits(G4Step* step, 
                                     G4TouchableHistory*)
{  
	G4Track* theTrack = step->GetTrack();
	//DataStorage->FillEmptyEntriesUpToCurrentTrack(theTrack);

	// if the particle is no optical photon, stop here
	G4cout << theTrack->GetParticleDefinition()->GetParticleName() <<G4endl;
		
	if((theTrack->GetParticleDefinition())!=(G4OpticalPhoton::OpticalPhotonDefinition()))
	{
		G4int trackID = theTrack->GetTrackID();
		//DataStorage->SetPhotonDetectorWasHit(trackID);
		//Stop track.		
		theTrack->SetTrackStatus(fStopAndKill);
		return false;
	}

	G4cout << G4OpticalPhoton::OpticalPhotonDefinition()->GetParticleName() << G4endl;
	// if optical photon was created inside the sensitive detector, stop here
	G4StepPoint * thePreStepPoint = step->GetPreStepPoint();
	G4StepPoint * thePostStepPoint = step->GetPostStepPoint();

	G4cout << theTrack->GetLogicalVolumeAtVertex()->GetName() << G4endl;
	G4cout << thePostStepPoint->GetPhysicalVolume()->GetLogicalVolume()->GetName() << G4endl;
	if(theTrack->GetLogicalVolumeAtVertex() == thePostStepPoint->GetPhysicalVolume()->GetLogicalVolume())
	{
		// Stop track.
		theTrack->SetTrackStatus(fStopAndKill);

		return false;
	}

	// This is now a true hit
	// create a hit instance and store all the data
	/*
	if(fSipmHitsCollection->entries()>5000)
	{
		//G4cout<<"Too many sipm hit"<<G4endl;
		theTrack->SetTrackStatus(fStopAndKill);
		return false;
	}
	*/
	MDM_SipmHit* hit = new MDM_SipmHit();
	hit->SetTrackID(theTrack->GetTrackID());
	G4ThreeVector hitPoint = (thePreStepPoint->GetPosition() + thePostStepPoint->GetPosition())/2;  // take hit position at middle of pre and post step
	hit->SetPos(hitPoint);

	G4int n_entry;
	n_entry = fSipmHitsCollection->entries();

	if(n_entry>0)
	{
		MDM_SipmHit* last_hit = (*fSipmHitsCollection)[n_entry-1];
		hit->SetPhotons(last_hit->GetPhotons());
		hit->incPhotons();

	}else
	{
		//first hit
		hit->SetPhotons(0);
	}

	fSipmHitsCollection->insert(hit);
	

	// Stop track.
	theTrack->SetTrackStatus(fStopAndKill);

	return true;

	  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_SipmSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
     G4int nofHits = fSipmHitsCollection->entries();
     G4cout << "\n-------->Hits Collection: in this event they are " << nofHits 
            << " hits in the tracker chambers: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fSipmHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
