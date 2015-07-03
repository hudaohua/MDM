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
// $Id: MDM_CalorHit.cc 69586 2013-05-08 14:20:11Z gcosmo $
//
/// \file MDM_CalorHit.cc
/// \brief Implementation of the MDM_CalorHit class

#include "MDM_CalorHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<MDM_CalorHit>* MDM_CalorHitAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_CalorHit::MDM_CalorHit()
 : G4VHit(),
   fEdep(0.),
   fTotalEdep(0.),
   fPos(0.),
   fTrackLength(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_CalorHit::~MDM_CalorHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_CalorHit::MDM_CalorHit(const MDM_CalorHit& right)
  : G4VHit()
{
  fEdep        = right.fEdep;
  fTotalEdep   = right.fTotalEdep;
  fPos         = right.fPos;
  fTrackLength = right.fTrackLength;
  fMomentum    = right.fMomentum;  //
  fParticle    = right.fParticle; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const MDM_CalorHit& MDM_CalorHit::operator=(const MDM_CalorHit& right)
{
  fEdep        = right.fEdep;
  fTotalEdep   = right.fTotalEdep;
  fPos         = right.fPos;
  fTrackLength = right.fTrackLength;
  fMomentum    = right.fMomentum;  //
  fParticle    = right.fParticle; 
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int MDM_CalorHit::operator==(const MDM_CalorHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_CalorHit::Print()
{
  G4cout
     << "Edep: " 
     << std::setw(7) << G4BestUnit(fEdep,"Energy")
     << "   Track length: " 
     << std::setw(7) << G4BestUnit(fTrackLength,"Length")
	 << "   Position: "
	 << std::setw(7) << G4BestUnit(fPos,"Length")
     << "   Momentum: "
     << std::setw(7) << G4BestUnit(fMomentum,"Energy")
	 << "   Particle: "
	 << std::setw(10) << fParticle->GetParticleName()
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
