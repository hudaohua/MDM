//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// ********************************************************************
//
// $Id: MDM_PrimaryGeneratorAction.cc 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file MDM_PrimaryGeneratorAction.cc
/// \brief Implementation of the MDM_PrimaryGeneratorAction class
// Use G4ParticleGun example
// Try to create a 10x10 deg block over an one pi sphere

#include "MDM_PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"  //define Pi
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

// using GPS
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_PrimaryGeneratorAction::MDM_PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0)
{
	//Using GPS
	  fParticleGun = new G4GeneralParticleSource();
	
/* Using ParticleGun
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);
  // direction will be adjusted at each event
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  fParticleGun->SetParticleEnergy(10.*keV);
*/

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_PrimaryGeneratorAction::~MDM_PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
	//using GPS
	// no code is required here
  //using particle gun
	/*
  G4double envSizeXY = 0;
  G4double envSizeZ = 0;

  // source uniformly distributed in a square area
  G4double r = 1.5*mm;
  G4double x0 = r * (G4UniformRand()-0.5);
  G4double y0 = r * (G4UniformRand()-0.5);
  G4double z0 = +2.8*cm;   // fix the z position
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
 */
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

