//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// ********************************************************************
//
// $Id: MDM_PrimaryGeneratorAction.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file MDM_PrimaryGeneratorAction.hh
/// \brief Definition of the MDM_PrimaryGeneratorAction class

#ifndef MDM_PrimaryGeneratorAction_h
#define MDM_PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

//using gun
//class G4ParticleGun;
//using GPS
class G4GeneralParticleSource;
class G4Event;
/// The primary generator action class with particle gun.
///


class MDM_PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    MDM_PrimaryGeneratorAction();    
    virtual ~MDM_PrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);         
  
    // method to access particle gun
	  //using Particle Gun
    //const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
	//using GPS
  const G4GeneralParticleSource* GetParticleGun() const { return fParticleGun; }

  private:
	  //using Particle Gun
   // G4ParticleGun*  fParticleGun; // pointer a to G4 gun class
	  //using GPS
	G4GeneralParticleSource*  fParticleGun;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


