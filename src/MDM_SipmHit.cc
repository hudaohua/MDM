//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// ********************************************************************
//
// $Id: MDM_SipmHit.cc 69586 2013-05-08 14:20:11Z gcosmo $
//
/// \file MDM_SipmHit.cc
/// \brief Implementation of the MDM_SipmHit class

#include "MDM_SipmHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<MDM_SipmHit>* MDM_SipmHitAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_SipmHit::MDM_SipmHit()
 : G4VHit(),
   fPos(0.0,0.0,0.0),
   fTrackID(-1),
   fHitTime(0.0),
   fTrackEnergy(0.0),
   sVertexVolumnName("")
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_SipmHit::~MDM_SipmHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_SipmHit::MDM_SipmHit(const MDM_SipmHit& right)
  : G4VHit()
{
  fPos			  = right.fPos;
  fTrackID        = right.fTrackID;
  fHitTime		  = right.fHitTime;
  fTrackEnergy	  = right.fTrackEnergy;
  sVertexVolumnName	  = right.sVertexVolumnName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const MDM_SipmHit& MDM_SipmHit::operator=(const MDM_SipmHit& right)
{
  fPos			  = right.fPos;
  fTrackID        = right.fTrackID;
  fHitTime        = right.fHitTime;
  fTrackEnergy    = right.fTrackEnergy;
  sVertexVolumnName   = right.sVertexVolumnName;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int MDM_SipmHit::operator==(const MDM_SipmHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_SipmHit::Print()
{
  G4cout
	 << "TrackID: "
	 << std::setw(7) << fTrackID
     << " Position: " 
     << std::setw(7) << G4BestUnit( fPos,"Length")
	 << " HitTime: "
	 << std::setw(7) << G4BestUnit( fHitTime, "Time")
	 << " TrackEnergy: "
	 << std::setw(7) << G4BestUnit( fTrackEnergy, "Energy")
	 << " VertexVolumnName: "
	 << std::setw(7) << sVertexVolumnName
	 << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
