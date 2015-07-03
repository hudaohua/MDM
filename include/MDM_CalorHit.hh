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
// $Id: MDM_CalorHit.hh 69223 2013-04-23 12:36:10Z gcosmo $
//
/// \file MDM_CalorHit.hh
/// \brief Definition of the MDM_CalorHit class

#ifndef MDM_CalorHit_h
#define MDM_CalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "tls.hh"

//Define the scintilator hit class
// This class defines the data collection when hit occurs


class MDM_CalorHit : public G4VHit
{
  public:
    MDM_CalorHit();
    MDM_CalorHit(const MDM_CalorHit&);
    virtual ~MDM_CalorHit();

    // operators
    const MDM_CalorHit& operator=(const MDM_CalorHit&);
    G4int operator==(const MDM_CalorHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // methods to handle data
    void IncEdep(G4double de);
	inline void IncTotalTrackLength(G4double len) {fTotalTrackLength += len;}

    // get methods
	G4double GetEdep() const; 
	G4double GetTrackLength() const;
	G4ThreeVector GetPos() const;
	G4ThreeVector GetMomentum() const;
	G4ParticleDefinition* GetParticle() const;
	G4double GetTotalEdep() const;
	inline G4double GetTotalTrackLength() const { return fTotalTrackLength;}

	// set methods
	void SetEdep(G4double de);
	void SetTrackLength(G4double track);
	void SetPos(G4ThreeVector xyz);
	void SetMomentum(G4ThreeVector mom);
	void SetParticle(G4ParticleDefinition* pdef);
	void SetTotalEdep(G4double tEdep);
	inline void SetTotalTrackLength(G4double len) { fTotalTrackLength = len; }
	      
  private:

	//data variable
    G4double fEdep;        ///< Energy deposit in the sensitive volume
	G4ThreeVector fPos;    // position info when hit
    G4ThreeVector fMomentum;  //
	G4double fTrackLength; ///< Track length in the  sensitive volume
	G4ParticleDefinition* fParticle; 
	G4double fTotalEdep;
	G4double fTotalTrackLength;
	
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<MDM_CalorHit> MDM_CalorHitsCollection;

extern G4ThreadLocal G4Allocator<MDM_CalorHit>* MDM_CalorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* MDM_CalorHit::operator new(size_t)
{
  if(!MDM_CalorHitAllocator)
    MDM_CalorHitAllocator = new G4Allocator<MDM_CalorHit>;
  void *hit;
  hit = (void *) MDM_CalorHitAllocator->MallocSingle();
  return hit;
}

inline void MDM_CalorHit::operator delete(void *hit)
{
  if(!MDM_CalorHitAllocator)
    MDM_CalorHitAllocator = new G4Allocator<MDM_CalorHit>;
  MDM_CalorHitAllocator->FreeSingle((MDM_CalorHit*) hit);
}

inline void MDM_CalorHit::IncEdep(G4double de) {
  fTotalEdep += de; 
}


//inline functions for accessing the privat data variables

inline G4double MDM_CalorHit::GetEdep() const { 
  return fEdep; 
}

inline G4double MDM_CalorHit::GetTotalEdep() const { 
  return fTotalEdep; 
}


inline G4double MDM_CalorHit::GetTrackLength() const { 
  return fTrackLength; 
}

inline G4ThreeVector MDM_CalorHit::GetPos() const { 
  return fPos; 
}

inline G4ThreeVector MDM_CalorHit::GetMomentum() const {
  return fMomentum;
}

inline G4ParticleDefinition* MDM_CalorHit::GetParticle() const {
	return fParticle;
}

inline void MDM_CalorHit::SetEdep(G4double de) { 
  fEdep = de; 
}

inline void MDM_CalorHit::SetTotalEdep(G4double tEdep) { 
  fTotalEdep = tEdep; 
}

inline void MDM_CalorHit::SetTrackLength(G4double track){ 
  fTrackLength = track; 
}

inline void MDM_CalorHit::SetPos(G4ThreeVector xyz) { 
  fPos = xyz; 
}

inline void MDM_CalorHit::SetMomentum(G4ThreeVector mom){
	fMomentum = mom;
}

inline void MDM_CalorHit::SetParticle(G4ParticleDefinition* pdef){
	fParticle = pdef;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
