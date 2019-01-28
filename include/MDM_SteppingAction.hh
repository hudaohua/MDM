//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// ********************************************************************
//
// $Id: MDM_SteppingAction.hh 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file MDM_SteppingAction.hh
/// \brief Definition of the MDM_SteppingAction class

#ifndef MDM_SteppingAction_h
#define MDM_SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include <G4OpBoundaryProcess.hh>
#include "globals.hh"


class MDM_DetectorConstruction;
class MDM_EventAction;
class G4LogicalVolume;

/// Stepping action class
/// 

class MDM_SteppingAction : public G4UserSteppingAction
{
  public:
    MDM_SteppingAction(const MDM_DetectorConstruction* detectorConstruction,MDM_EventAction* eventAction);
    virtual ~MDM_SteppingAction();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);

  private:
    MDM_EventAction*  fEventAction;
    G4LogicalVolume* fScoringVolume;
	const MDM_DetectorConstruction* fDetConstruction;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
