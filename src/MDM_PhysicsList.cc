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
// $Id: MDM_PhysicsList.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file MDM_PhysicsList.cc
/// \brief Implementation of the MDM_PhysicsList class

#include "MDM_PhysicsList.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_PhysicsList::MDM_PhysicsList() 
: G4VModularPhysicsList(){
  SetVerboseLevel(0);
  //default cut value (1.0mm)
  defaultCutValue = 0.1*mm;  //cut range should be defined less than 1mm

  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  RegisterPhysics(opticalPhysics);

  /* not valid any more for Geant4.11.0; so comment out
   * shall use G4OpticalParameters class *
  opticalPhysics->SetScintillationYieldFactor(1.0);
  opticalPhysics->SetScintillationExcitationRatio(0.0);

  opticalPhysics->SetMaxNumPhotonsPerStep(100);
  opticalPhysics->SetMaxBetaChangePerStep(10.0);

  opticalPhysics->SetTrackSecondariesFirst(kCerenkov,true);
  opticalPhysics->SetTrackSecondariesFirst(kScintillation,true);
  */
  // Default physics
  RegisterPhysics(new G4DecayPhysics());

  // Radioactive decay
  RegisterPhysics(new G4RadioactiveDecayPhysics());

  // EM physics option 4
  // most accurate
   RegisterPhysics(new G4EmStandardPhysics_option4()); 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_PhysicsList::~MDM_PhysicsList()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_PhysicsList::SetCuts()
{
  //G4VUserPhysicsList::SetCuts();
	// cuts means the "secondary production threshold distance"
	// if primary no longer has enough energy to produce secondaries which travel at least 1mm(cut value),
	// two things will happen
	// - discrete energy loss ceases(no more secondaries produced)
	// - the primary is tracked down to zero energy using continuous energy loss
	SetCutsWithDefault();

	// experiment code to set cut to different value for different particals
	// has to set cut value for gamma first then e+/e-
	// not sure about proton
	
	SetCutValue(0.01*mm,"gamma");
	SetCutValue(0.01*mm,"e-");
	SetCutValue(0.01*mm,"e+");
	SetCutValue(0.05*mm,"proton");
	

	if(verboseLevel>0) DumpCutValuesTable();
}  
