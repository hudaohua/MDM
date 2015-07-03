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
// $Id: MDM_DetectorConstruction.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file MDM_DetectorConstruction.hh
/// \brief Definition of the MDM_DetectorConstruction class

#ifndef MDM_DetectorConstruction_h
#define MDM_DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class MDM_DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    MDM_DetectorConstruction();
    virtual ~MDM_DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

	virtual void ConstructSDandField();
    
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
	// get methods
    //
    const G4VPhysicalVolume* GetShape1PV() const;
    const G4VPhysicalVolume* GetShape2PV() const;

  protected:
    G4LogicalVolume*  fScoringVolume;


  private:
	G4VPhysicalVolume* physDetector;
	G4VPhysicalVolume* physSensor;
};

// inline functions

inline const G4VPhysicalVolume* MDM_DetectorConstruction::GetShape1PV() const { 
  return physDetector; 
}

inline const G4VPhysicalVolume* MDM_DetectorConstruction::GetShape2PV() const  { 
  return physSensor; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

