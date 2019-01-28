//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
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
	
    // get methods
    G4ThreeVector GetPos() const;
	G4int GetTrackID() const;
	G4double GetHitTime() const;
	G4double GetTrackEnergy() const;
	G4String GetVertexVolumnName() const;

	// set methods
	void SetPos(G4ThreeVector pos);
	void SetTrackID(G4int tID);
	void SetHitTime(G4double TH);
	void SetTrackEnergy(G4double trackEne);
	void SetVertexVolumnName(G4String name);
      
  private:
	G4ThreeVector fPos;
	G4int fTrackID;
	G4double fHitTime;
	G4double fTrackEnergy;
	G4String sVertexVolumnName;
   
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




inline G4ThreeVector MDM_SipmHit::GetPos() const { 
  return fPos; 
}

inline G4int MDM_SipmHit::GetTrackID() const { 
  return fTrackID; 
}

inline G4double MDM_SipmHit::GetHitTime() const {
	return fHitTime;
}

inline G4double MDM_SipmHit::GetTrackEnergy() const {
	return fTrackEnergy;
}

inline G4String MDM_SipmHit::GetVertexVolumnName() const {
	return sVertexVolumnName;
}




inline void MDM_SipmHit::SetPos(G4ThreeVector pos) { 
  fPos = pos; 
}

inline void MDM_SipmHit::SetTrackID(G4int tID) { 
  fTrackID = tID; 
}

inline void MDM_SipmHit::SetHitTime(G4double TH){
	fHitTime = TH;
}

inline void MDM_SipmHit::SetTrackEnergy(G4double trackEne){
	fTrackEnergy = trackEne;
}

inline void MDM_SipmHit::SetVertexVolumnName(G4String name){
	sVertexVolumnName = name;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
