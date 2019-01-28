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
#include <G4UnitsTable.hh>



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
  G4cout << "SiPm sensitiveDetectorName: " << SensitiveDetectorName << "- size "<< collectionName.size() <<G4endl;

	// hcID1 = GetCollectionID(0);
  }
  hce->AddHitsCollection( hcID1, fSipmHitsCollection ); 


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
	//theTrack->SetTrackStatus(fStopAndKill);
	// try a empty function
	//return true;

	G4String volume_name = theTrack->GetLogicalVolumeAtVertex()->GetName();

	// if the particle is no optical photon, stop here

	//debug
	/*
			  G4cout
		     << "TotalE_Dep:" << std::setw(2) << G4BestUnit(step->GetTotalEnergyDeposit(),"Energy")
			// << " TrackTotalE:" << std::setw(2) << G4BestUnit(theTrack->GetTotalEnergy(),"Energy")
		    // << " Delta_Ene:" << std::setw(2) << G4BestUnit(step->GetDeltaEnergy(),"Energy")  // this function is going to be removed soon
			// << " Track_len:" << std::setw(2) << G4BestUnit(step->GetStepLength(),"Length")
			<< " Pos:"<< std::setw(2) << step->GetPostStepPoint()->GetPosition()
			// << " Particle :" << std::setw(5) << theTrack->GetParticleDefinition()->GetParticleName()
			 << " Track sta:" << std::setw(5) << theTrack->GetTrackStatus()
			 << " Parent ID:" << std::setw(2) << theTrack->GetParentID()
			 << " Track  ID:" << std::setw(2) << theTrack->GetTrackID()
			 << " step   no:" << std::setw(2) << theTrack->GetCurrentStepNumber()
			// << " LocalTime:" << std::setw(2) << G4BestUnit(theTrack->GetLocalTime(),"Time")
			// << " GlobalTim:" << std::setw(2) << G4BestUnit(theTrack->GetGlobalTime(),"Time")
			 << " VertexVolumn:" << std::setw(2) << theTrack->GetLogicalVolumeAtVertex()->GetName()
		     << G4endl;
  */
	if((theTrack->GetParticleDefinition())!=(G4OpticalPhoton::OpticalPhotonDefinition()))
	{
		//Stop track.		
		theTrack->SetTrackStatus(fStopAndKill);
		return false;
	}

	// if optical photon was created inside the sensitive detector, stop here
	if(volume_name == "Sensor" )
	{
		// Stop track.
		theTrack->SetTrackStatus(fStopAndKill);
		return false;
	}

	// debug removed later
		// print the information related to each hit
	/*
		  G4cout
		     << "TotalE_Dep:" << std::setw(2) << G4BestUnit(step->GetTotalEnergyDeposit(),"Energy")
			// << " TrackTotalE:" << std::setw(2) << G4BestUnit(theTrack->GetTotalEnergy(),"Energy")
		    // << " Delta_Ene:" << std::setw(2) << G4BestUnit(step->GetDeltaEnergy(),"Energy")  // this function is going to be removed soon
			// << " Track_len:" << std::setw(2) << G4BestUnit(step->GetStepLength(),"Length")
			<< " Pos:"<< std::setw(2) << step->GetPostStepPoint()->GetPosition()
			// << " Particle :" << std::setw(5) << theTrack->GetParticleDefinition()->GetParticleName()
			 << " Track sta:" << std::setw(5) << theTrack->GetTrackStatus()
			 << " Parent ID:" << std::setw(2) << theTrack->GetParentID()
			 << " Track  ID:" << std::setw(2) << theTrack->GetTrackID()
			 << " step   no:" << std::setw(2) << theTrack->GetCurrentStepNumber()
			// << " LocalTime:" << std::setw(2) << G4BestUnit(theTrack->GetLocalTime(),"Time")
			 << " GlobalTim:" << std::setw(2) << G4BestUnit(theTrack->GetGlobalTime(),"Time")
			 << " VertexVolumn:" << std::setw(2) << theTrack->GetLogicalVolumeAtVertex()->GetName()
		     << G4endl;
	*/
	// in csv format
	/*
		G4cout
		 << std::setw(7) << G4BestUnit(step->GetTotalEnergyDeposit(),"Energy")  //"TotalE_Dep:" 
		 <<","<< std::setw(7) << G4BestUnit(theTrack->GetTotalEnergy(),"Energy")  // " TrackTotalE:" 
	    //<<","<< std::setw(7) << G4BestUnit(step->GetDeltaEnergy(),"Energy") //" Delta_Ene:"
		 <<","<< std::setw(7) << G4BestUnit(step->GetStepLength(),"Length") // " Track_len:" 
		 <<","<< std::setw(7) << theTrack->GetParticleDefinition()->GetParticleName() //" Particle :" 
		 <<","<< std::setw(7) << theTrack->GetTrackStatus() //" Track sta:" 
		 <<","<< std::setw(7) << theTrack->GetParentID() //" Parent ID:"
		 <<","<< std::setw(7) << theTrack->GetTrackID() // " Track  ID:" 
		 <<","<< std::setw(7) << theTrack->GetCurrentStepNumber() //" step   no:" 
		 <<","<< std::setw(7) << G4BestUnit(theTrack->GetLocalTime(),"Time") //" LocalTime:" 
		 <<","<< std::setw(7) << G4BestUnit(theTrack->GetGlobalTime(),"Time")  //" GlobalTim:"
		 <<",Hit in SiPM"
	 << G4endl;
	*/

	// This is now a true hit

	MDM_SipmHit* hit = new MDM_SipmHit();
	// acquire interested information from the step

	G4StepPoint * thePreStepPoint = step->GetPreStepPoint();
	G4StepPoint * thePostStepPoint = step->GetPostStepPoint();
	G4int trackID = theTrack->GetTrackID();
	G4ThreeVector hitPoint = (thePreStepPoint->GetPosition() + thePostStepPoint->GetPosition())/2;  // take hit position at middle of pre and post step
	G4double HT = theTrack->GetGlobalTime();  // time since the event was created.
	//G4double HT = theTrack->GetLocalTime();  // time since the track was created.
	G4double trackEnergy = theTrack->GetTotalEnergy();

	// save to hit
	hit->SetHitTime(HT);
	hit->SetTrackID(trackID);
	hit->SetPos(hitPoint);
	hit->SetTrackEnergy(trackEnergy);
	hit->SetVertexVolumnName(volume_name);
	
	//hit->SetPhotons(theTrack->GetParentID()); //temperaly

	// insert the hit into collection
	fSipmHitsCollection->insert(hit);	

	// Stop track.
	theTrack->SetTrackStatus(fStopAndKill);

	return true;

	  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_SipmSD::EndOfEvent(G4HCofThisEvent* hitCollection)
{
  if ( verboseLevel>1 ) { 
     G4int nofHits = fSipmHitsCollection->entries();
     G4cout << "\n-------->Hits Collection: in this event they are " << nofHits 
            << " hits in the SiPMs: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fSipmHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
