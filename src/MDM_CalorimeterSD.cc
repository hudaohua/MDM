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
#include <G4OpticalPhoton.hh>




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

  // Add this collection in hce
  static G4int hcID = -1;
  if(hcID<0)
  {
	  hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
	  G4cout << "Calorimeter hcID: " << hcID << "- size "<< collectionName.size() <<G4endl;
	  //hcID = GetCollectionID(0);
  }
	hce->AddHitsCollection( hcID, fCalorHitsCollection ); 

  // Create hits
  // fNofCells for cells + one more for total sums 
  //for (G4int i=0; i<fNofCells+1; i++ ) {
   // fHitsCollection->insert(new MDM_CalorHit());
  //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



G4bool MDM_CalorimeterSD::ProcessHits(G4Step* step, 
                                     G4TouchableHistory*)
{  

	G4double edep = step->GetTotalEnergyDeposit();
	G4double tLen = step->GetStepLength();

	if(edep ==0) return false;   // ignore if no energy deposit

	/*
	if(fCalorHitsCollection->entries()>5000)
	{
		//G4cout<<"Too many scintillator hit"<<G4endl;
		return false;
	}
	*/

	G4StepPoint* preStep = step->GetPreStepPoint();

	MDM_CalorHit* hit = new MDM_CalorHit();

	hit->SetPos(preStep->GetPosition());
	hit->SetMomentum(preStep->GetMomentum());
	hit->SetEdep(edep);
	hit->SetParticle(step->GetTrack()->GetDefinition());
	hit->SetTrackLength(tLen); 

	//acquire the last deposit energy
	G4int n_entry;
	n_entry = fCalorHitsCollection->entries();
	if(n_entry>0)
	{
		MDM_CalorHit* last_hit = (*fCalorHitsCollection)[n_entry-1];
		hit->SetTotalEdep(last_hit->GetTotalEdep());
		hit->IncEdep(edep);
		hit->SetTotalTrackLength(last_hit->GetTotalTrackLength());
		hit->IncTotalTrackLength(tLen);
	}else
	{
		// first hit
		hit->SetTotalEdep(edep);
		hit->SetTotalTrackLength(tLen);
	}
	fCalorHitsCollection->insert(hit);


	return true;

	/*
  // energy deposit
  G4double edep = step->GetTotalEnergyDeposit();
  
  // step length
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
  }

  if ( edep==0. && stepLength == 0. ) return false;      

  G4TouchableHistory* touchable
    = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());
    
  // Get calorimeter cell id 
  G4int layerNumber = touchable->GetReplicaNumber(1);
  
  // Get hit accounting data for this cell
  MDM_CalorHit* hit = (*fHitsCollection)[layerNumber];
  if ( ! hit ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hit " << layerNumber; 
    G4Exception("MDM_CalorimeterSD::ProcessHits()",
      "MyCode0004", FatalException, msg);
  }         

  // Get hit for total accounting
  MDM_CalorHit* hitTotal 
    = (*fHitsCollection)[fHitsCollection->entries()-1];
  
  // Add values
  hit->Add(edep, stepLength);
  hitTotal->Add(edep, stepLength); 
      
  return true;

  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_CalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
     G4int nofHits = fCalorHitsCollection->entries();
     G4cout << "\n-------->Hits Collection: in this event they are " << nofHits 
            << " hits in the tracker chambers: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fCalorHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
