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
// $Id: MDM_ActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file MDM_ActionInitialization.cc
/// \brief Implementation of the MDM_ActionInitialization class

#include "MDM_ActionInitialization.hh"
#include "MDM_PrimaryGeneratorAction.hh"
#include "MDM_RunAction.hh"
#include "MDM_EventAction.hh"
#include "MDM_SteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_ActionInitialization::MDM_ActionInitialization(MDM_DetectorConstruction* detConstruction)
 : G4VUserActionInitialization(),
 fDetConstruction(detConstruction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_ActionInitialization::~MDM_ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_ActionInitialization::BuildForMaster() const
{
  SetUserAction(new MDM_RunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_ActionInitialization::Build() const
{
  SetUserAction(new MDM_PrimaryGeneratorAction);
  SetUserAction(new MDM_RunAction);
  
  MDM_EventAction* eventAction = new MDM_EventAction;
  SetUserAction(eventAction);
  
  SetUserAction(new MDM_SteppingAction(fDetConstruction,eventAction));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
