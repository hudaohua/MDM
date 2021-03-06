//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// ********************************************************************
// $Id: MDM_CalorimeterSD.cc 75215 2013-10-29 16:07:06Z gcosmo $
//
/// \file MDM_CalorimeterSD.cc
/// \brief Implementation of the MDM_CalorimeterSD class

#include "MDM_CalorimeterSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4OpticalPhoton.hh"
#include "G4UnitsTable.hh"




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_CalorimeterSD::MDM_CalorimeterSD(
                            const G4String& name, 
                            const G4String& hitsCollectionName,
                            G4int nofCells)
 : G4VSensitiveDetector(name),
   fCalorHitsCollection(0),
   fNofCells(nofCells)
{

	//collectionName is a predefined vector
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_CalorimeterSD::~MDM_CalorimeterSD() 
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_CalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fCalorHitsCollection 
    = new MDM_CalorHitsCollection(SensitiveDetectorName, collectionName[0]); 
    //debug
  //G4cout <<"In CalorSD,the sensitivedetector name is"<< SensitiveDetectorName <<" the collectionName[0] is" <<collectionName[0] <<G4endl;

  // Add this collection in hce
  static G4int hcID = -1;
  if(hcID<0)
  {
	  hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
	  //debug
	  G4cout << "Scintillator - SensitiveDetectorName: " << SensitiveDetectorName << "- size "<< collectionName.size() <<" hcID="<<hcID<<G4endl;
	  //hcID = GetCollectionID(0);
  }

  hce->AddHitsCollection( hcID, fCalorHitsCollection );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



G4bool MDM_CalorimeterSD::ProcessHits(G4Step* step, 
                                     G4TouchableHistory* history)
{  

	G4double edep = step->GetTotalEnergyDeposit();

	if(edep == 0. ) return false;   // ignore if no energy deposit

	//don't record any opticalphoton hit
	G4Track* theTrack = step->GetTrack();
	if((theTrack->GetParticleDefinition())==(G4OpticalPhoton::OpticalPhotonDefinition()))
	{
		return false;
	}

	G4double tLen = step->GetStepLength();
// debug removed later
	// print the information related to each hit
	// use a CSV format
	/*
	  G4cout
	     << "TotalE_Dep:" << std::setw(2) << G4BestUnit(edep,"Energy")
		 << " TrackTotalE:" << std::setw(2) << G4BestUnit(theTrack->GetTotalEnergy(),"Energy")
	    // << " Delta_Ene:" << std::setw(2) << G4BestUnit(step->GetDeltaEnergy(),"Energy")
		 << " Track_len:" << std::setw(2) << G4BestUnit(step->GetStepLength(),"Length")
		 << " Particle :" << std::setw(5) << theTrack->GetParticleDefinition()->GetParticleName()
		 << " Track sta:" << std::setw(5) << theTrack->GetTrackStatus()
		 << " Parent ID:" << std::setw(2) << theTrack->GetParentID()
		 << " Track  ID:" << std::setw(2) << theTrack->GetTrackID()
		 << " step   no:" << std::setw(2) << theTrack->GetCurrentStepNumber()
		 << " LocalTime:" << std::setw(2) << G4BestUnit(theTrack->GetLocalTime(),"Time")
		 << " GlobalTim:" << std::setw(2) << G4BestUnit(theTrack->GetGlobalTime(),"Time")
	     << G4endl;
		 */
	/*
	G4cout
		 << std::setw(7) << G4BestUnit(edep,"Energy")  //"TotalE_Dep:" 
		 <<","<< std::setw(7) << G4BestUnit(theTrack->GetTotalEnergy(),"Energy")  // " TrackTotalE:" 
	    // <<  << std::setw(7) << G4BestUnit(step->GetDeltaEnergy(),"Energy") //" Delta_Ene:"
		 <<","<< std::setw(7) << G4BestUnit(step->GetStepLength(),"Length") // " Track_len:" 
		 <<","<< std::setw(7) << theTrack->GetParticleDefinition()->GetParticleName() //" Particle :" 
		 <<","<< std::setw(7) << theTrack->GetTrackStatus() //" Track sta:" 
		 <<","<< std::setw(7) << theTrack->GetParentID() //" Parent ID:"
		 <<","<< std::setw(7) << theTrack->GetTrackID() // " Track  ID:" 
		 <<","<< std::setw(7) << theTrack->GetCurrentStepNumber() //" step   no:" 
		 <<","<< std::setw(7) << G4BestUnit(theTrack->GetLocalTime(),"Time") //" LocalTime:" 
		 <<","<< std::setw(7) << G4BestUnit(theTrack->GetGlobalTime(),"Time")  //" GlobalTim:"
		 <<",Hit in scintillator"
	 << G4endl;
	 */
	G4StepPoint* preStep = step->GetPreStepPoint();

	MDM_CalorHit* hit = new MDM_CalorHit();

	hit->SetPos(preStep->GetPosition());
	hit->SetMomentum(preStep->GetMomentum());
	hit->SetEdep(edep);
	hit->SetParticle(theTrack->GetDefinition());
	hit->SetTrackLength(tLen); 
	hit->SetTrackID(theTrack->GetTrackID());

	//acquire the last deposit energy
	//G4int n_entry = 0;
	//n_entry = fCalorHitsCollection->entries();
	//if(n_entry>0)
	//{
	//	MDM_CalorHit* last_hit = (*fCalorHitsCollection)[n_entry-1];
	//	hit->SetTotalEdep(last_hit->GetTotalEdep());
	//	hit->IncEdep(edep);
	//	hit->SetTotalTrackLength(last_hit->GetTotalTrackLength());
	//	hit->IncTotalTrackLength(tLen);
	//	// redundant code. All the photon particle has been filtered at the start of this function.
	//	/*
	//	if(hit->GetParticle()==G4OpticalPhoton::OpticalPhotonDefinition())
	//	{
	//		hit->SetTotalPhotons(last_hit->GetTotalPhotons());
	//		hit->InctotalPhotons(1);  // add total photons by 1
	//		hit->SetHitTime(theTrack->GetGlobalTime());  //time since the event was created
	//	}
	//	*/

	//}else
	//{
	//	// first hit
	//	hit->SetTotalEdep(edep);
	//	hit->SetTotalTrackLength(tLen);
	//	// redundant code. All the photon particle has been filtered at the start of this function.
	//	/*
	//	if(hit->GetParticle()==G4OpticalPhoton::OpticalPhotonDefinition())
	//	{
	//		hit->SetTotalPhotons(1); // first photon generated.
	//		hit->SetHitTime(theTrack->GetGlobalTime());  //time since the event was created
	//	}
	//	*/
	//}
	fCalorHitsCollection->insert(hit);


	return true;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_CalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
     G4int nofHits = fCalorHitsCollection->entries();
     G4cout << "\n-------->Hits Collection: in this event they are " << nofHits 
            << " hits in the gamma detector: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fCalorHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
