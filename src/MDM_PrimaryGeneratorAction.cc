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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_PrimaryGeneratorAction::MDM_PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0), 
  fEnvelopeBox(0)
{
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
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.
  
  G4double envSizeXY = 0;
  G4double envSizeZ = 0;

  // This pieces of code may be the source of memory running out
  // or at least it triggers a bad_cast warning
  //ops maybe not right, bad_cast warning is still there.
  /*
  if (!fEnvelopeBox)
  {
    G4LogicalVolume* envLV
      = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
    if ( envLV ) fEnvelopeBox = dynamic_cast<G4Box*>(envLV->GetSolid());
  }

  if ( fEnvelopeBox ) {
    envSizeXY = fEnvelopeBox->GetXHalfLength()*2.;
    envSizeZ = fEnvelopeBox->GetZHalfLength()*2.;
  }  
  else  {
    G4ExceptionDescription msg;
    msg << "Envelope volume of box shape not found.\n"; 
    msg << "Perhaps you have changed geometry.\n";
    msg << "The gun will be place at the center.";
    G4Exception("MDM_PrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002",JustWarning,msg);
  }
  */



  // radiation source is a circular plate
  // locate at (0,0,1.5cm) of the evelop
  // shooting torwards (0,0,-1) direction
  // circular radius is 0.5cm
  //G4double radius = 0.15*cm; //match the sintillator size
  //G4double r = radius * G4UniformRand();
  //G4double phi =twopi*G4UniformRand();
  //G4double x0 = r * std::cos(phi);
 // G4double y0 = r * std::sin(phi);
 // G4double z0 = +0.4 * envSizeZ; 

  // source uniformly distributed in a square area
  G4double r = 1.5*mm;
  G4double x0 = r * (G4UniformRand()-0.5);
  G4double y0 = r * (G4UniformRand()-0.5);
  G4double z0 = +2.8*cm;   // fix the z position
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  //generate a random direction
  //G4double d_x = (G4UniformRand()-0.5)*3*mm;
 // G4double d_y = (G4UniformRand()-0.5)*3*mm;
 // G4double d_z = -1*cm; 
 // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(d_x,d_y,d_z));
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

