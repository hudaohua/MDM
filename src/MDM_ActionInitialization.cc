//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *
//
// * version log history
// * v0.1	08/05/2015	convert from old test program to MDM project  *
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
