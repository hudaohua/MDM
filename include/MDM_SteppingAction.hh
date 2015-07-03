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
