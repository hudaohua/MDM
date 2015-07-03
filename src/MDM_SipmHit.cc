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
   fPhotons(0.),
   fPos(0.),
   fTrackID(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_SipmHit::~MDM_SipmHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_SipmHit::MDM_SipmHit(const MDM_SipmHit& right)
  : G4VHit()
{
  fPhotons        = right.fPhotons;
  fPos			  = right.fPos;
  fTrackID        = right.fTrackID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const MDM_SipmHit& MDM_SipmHit::operator=(const MDM_SipmHit& right)
{
  fPhotons        = right.fPhotons;
  fPos			  = right.fPos;
  fTrackID        = right.fTrackID;

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
     << " Photons: " 
     << std::setw(7) << fPhotons
     << " Position: " 
     << std::setw(7) << G4BestUnit( fPos,"Length")
	 << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
