//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// ********************************************************************
//
// $Id: MDM_PhysicsList.hh 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file MDM_PhysicsList.hh
/// \brief Definition of the MDM_PhysicsList class

#ifndef MDM_PhysicsList_h
#define MDM_PhysicsList_h 1

#include "G4VModularPhysicsList.hh"

/// Modular physics list
///
/// It includes the folowing physics builders
/// - G4DecayPhysics
/// - G4RadioactiveDecayPhysics
/// - G4EmStandardPhysics

class MDM_PhysicsList: public G4VModularPhysicsList
{
public:
  MDM_PhysicsList();
  virtual ~MDM_PhysicsList();

  virtual void SetCuts();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

