//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// ********************************************************************
//
// $Id: MDM_SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file MDM_SteppingAction.cc
/// \brief Implementation of the MDM_SteppingAction class

#include "MDM_SteppingAction.hh"
#include "MDM_EventAction.hh"
#include "MDM_DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_SteppingAction::MDM_SteppingAction(
	const MDM_DetectorConstruction* detectorConstruction,
	MDM_EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fDetConstruction(detectorConstruction),
  fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_SteppingAction::~MDM_SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_SteppingAction::UserSteppingAction(const G4Step* step)
{

	// no action here
	// all the data processing is handled by hit and sensitive detector

	// only process photons
/*
  G4Track* theTrack = step->GetTrack();

	// if the particle is no optical photon, stop here
	if(theTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) return;
	
	G4ThreeVector PreviousPreStepPointMomentum;
	G4ThreeVector PreviousPreStepPointPosition;
	G4double PreviousPreStepPointTime;
	G4ThreeVector PreviousPreStepPointPolarisation;

		G4StepPoint * preStepPoint = step->GetPreStepPoint();
		PreviousPreStepPointMomentum = preStepPoint->GetMomentum();
		PreviousPreStepPointPosition = preStepPoint->GetPosition();
		PreviousPreStepPointTime = preStepPoint->GetGlobalTime();
		PreviousPreStepPointPolarisation = preStepPoint->GetPolarization();

*/

// no need to extract information at stepping , all handled by sensitive detector and hit
	/*

//		G4double stepLength = step->GetStepLength();


	// get volume of the current step
	G4VPhysicalVolume* pvolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

	
	// energy deposit
	G4double edep = step->GetTotalEnergyDeposit();
	 // step length
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
  }
      
  if ( pvolume == fDetConstruction->GetShape1PV() ) {
    fEventAction->AddAbs(edep,stepLength);
  }
  
  if ( pvolume == fDetConstruction->GetShape2PV() ) {
    fEventAction->AddGap(edep,stepLength);
  }

  if (!fScoringVolume) { 
    const MDM_DetectorConstruction* detectorConstruction
      = static_cast<const MDM_DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   
  }

  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
      
  // check if we are in scoring volume
  if (volume != fScoringVolume) return;

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep);  
 
 */
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

