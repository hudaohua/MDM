//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
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

