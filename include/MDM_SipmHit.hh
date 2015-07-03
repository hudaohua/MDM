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
// $Id: MDM_SipmHit.hh 69223 2013-04-23 12:36:10Z gcosmo $
//
/// \file MDM_SipmHit.hh
/// \brief Definition of the MDM_SipmHit class

#ifndef MDM_SipmHit_h
#define MDM_SipmHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

/// Sipm hit class
///
/// It defines data members to store the the energy deposit and track lengths
/// of charged particles in a selected volume:
/// - fEdep, fTrackLength

class MDM_SipmHit : public G4VHit
{
  public:
    MDM_SipmHit();
    MDM_SipmHit(const MDM_SipmHit&);
    virtual ~MDM_SipmHit();

    // operators
    const MDM_SipmHit& operator=(const MDM_SipmHit&);
    G4int operator==(const MDM_SipmHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // methods to handle data
    void incPhotons();

    // get methods
    G4double GetPhotons() const;
    G4ThreeVector GetPos() const;
	G4int GetTrackID() const;

	// set methods
	void SetPhotons(G4int photons);
	void SetPos(G4ThreeVector pos);
	void SetTrackID(G4int tID);
      
  private:
	G4int fPhotons;
	G4ThreeVector fPos;
	G4int fTrackID;
   
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<MDM_SipmHit> MDM_SipmHitsCollection;

extern G4ThreadLocal G4Allocator<MDM_SipmHit>* MDM_SipmHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* MDM_SipmHit::operator new(size_t)
{
  if(!MDM_SipmHitAllocator)
    MDM_SipmHitAllocator = new G4Allocator<MDM_SipmHit>;
  void *hit;
  hit = (void *) MDM_SipmHitAllocator->MallocSingle();
  return hit;
}

inline void MDM_SipmHit::operator delete(void *hit)
{
  if(!MDM_SipmHitAllocator)
    MDM_SipmHitAllocator = new G4Allocator<MDM_SipmHit>;
  MDM_SipmHitAllocator->FreeSingle((MDM_SipmHit*) hit);
}


inline void MDM_SipmHit::incPhotons() {
  fPhotons++;
}

inline G4double MDM_SipmHit::GetPhotons() const { 
  return fPhotons; 
}

inline G4ThreeVector MDM_SipmHit::GetPos() const { 
  return fPos; 
}

inline G4int MDM_SipmHit::GetTrackID() const { 
  return fTrackID; 
}

inline void MDM_SipmHit::SetPhotons(G4int photons)  { 
  fPhotons = photons; 
}

inline void MDM_SipmHit::SetPos(G4ThreeVector pos) { 
  fPos = pos; 
}

inline void MDM_SipmHit::SetTrackID(G4int tID) { 
  fTrackID = tID; 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
