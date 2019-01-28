//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// ********************************************************************
//
// $Id: MDM_BetaHit.hh 69223 2013-04-23 12:36:10Z gcosmo $
//
/// \file MDM_BetaHit.hh
/// \brief Definition of the MDM_BetaHit class

#ifndef MDM_BetaHit_h
#define MDM_BetaHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "tls.hh"

//Define the scintilator hit class
// This class defines the data collection when hit occurs


class MDM_BetaHit : public G4VHit
{
  public:
    MDM_BetaHit();
    MDM_BetaHit(const MDM_BetaHit&);
    virtual ~MDM_BetaHit();

    // operators
    const MDM_BetaHit& operator=(const MDM_BetaHit&);
    G4int operator==(const MDM_BetaHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // methods to handle data
    inline void IncEdep(G4double de) { fTotalEdep += de; }
	inline void IncTotalTrackLength(G4double len) {fTotalTrackLength += len;}
	inline void InctotalPhotons(G4int dn) { iTotalPhotons += dn;}

    // get methods
	G4double GetEdep() const; 
	G4double GetTrackLength() const;
	G4ThreeVector GetPos() const;
	G4ThreeVector GetMomentum() const;
	G4ParticleDefinition* GetParticle() const;
	G4double GetTotalEdep() const;
	inline G4double GetTotalTrackLength() const { return fTotalTrackLength;}
	inline G4int GetTotalPhotons() const { return iTotalPhotons;}
	inline G4double GetHitTime() const { return fHitTime;}
	inline G4int GetTrackID() const { return iTrackID;}

	// set methods
	void SetEdep(G4double de);
	void SetTrackLength(G4double track);
	void SetPos(G4ThreeVector xyz);
	void SetMomentum(G4ThreeVector mom);
	void SetParticle(G4ParticleDefinition* pdef);
	void SetTotalEdep(G4double tEdep);
	inline void SetTotalTrackLength(G4double len) { fTotalTrackLength = len; }
	inline void SetTotalPhotons(G4int totalPhotons) {iTotalPhotons = totalPhotons;}
	inline void SetHitTime(G4double hitTime){ fHitTime = hitTime;}
	inline void SetTrackID(G4int tID){ iTrackID = tID;}
	      
  private:

	//data variable
    G4double fEdep;        ///< Energy deposit in the sensitive volume
	G4ThreeVector fPos;    // position info when hit
    G4ThreeVector fMomentum;  //
	G4double fTrackLength; ///< Track length in the  sensitive volume
	G4ParticleDefinition* fParticle; 
	G4double fTotalEdep;
	G4double fTotalTrackLength;
	G4int iTotalPhotons;
	G4double fHitTime;
	G4int iTrackID;
	
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<MDM_BetaHit> MDM_BetaHitsCollection;

extern G4ThreadLocal G4Allocator<MDM_BetaHit>* MDM_BetaHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* MDM_BetaHit::operator new(size_t)
{
  if(!MDM_BetaHitAllocator)
    MDM_BetaHitAllocator = new G4Allocator<MDM_BetaHit>;
  void *hit;
  hit = (void *) MDM_BetaHitAllocator->MallocSingle();
  return hit;
}

inline void MDM_BetaHit::operator delete(void *hit)
{
  if(!MDM_BetaHitAllocator)
    MDM_BetaHitAllocator = new G4Allocator<MDM_BetaHit>;
  MDM_BetaHitAllocator->FreeSingle((MDM_BetaHit*) hit);
}



//inline functions for accessing the privat data variables

inline G4double MDM_BetaHit::GetEdep() const { 
  return fEdep; 
}

inline G4double MDM_BetaHit::GetTotalEdep() const { 
  return fTotalEdep; 
}


inline G4double MDM_BetaHit::GetTrackLength() const { 
  return fTrackLength; 
}

inline G4ThreeVector MDM_BetaHit::GetPos() const { 
  return fPos; 
}

inline G4ThreeVector MDM_BetaHit::GetMomentum() const {
  return fMomentum;
}

inline G4ParticleDefinition* MDM_BetaHit::GetParticle() const {
	return fParticle;
}

inline void MDM_BetaHit::SetEdep(G4double de) { 
  fEdep = de; 
}

inline void MDM_BetaHit::SetTotalEdep(G4double tEdep) { 
  fTotalEdep = tEdep; 
}

inline void MDM_BetaHit::SetTrackLength(G4double track){ 
  fTrackLength = track; 
}

inline void MDM_BetaHit::SetPos(G4ThreeVector xyz) { 
  fPos = xyz; 
}

inline void MDM_BetaHit::SetMomentum(G4ThreeVector mom){
	fMomentum = mom;
}

inline void MDM_BetaHit::SetParticle(G4ParticleDefinition* pdef){
	fParticle = pdef;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
