//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// ********************************************************************
//
// $Id: MDM_BetaSD.hh 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file MDM_BetaSD.hh
/// \brief Definition of the MDM_BetaSD class

#ifndef MDM_BetaSD_h
#define MDM_BetaSD_h 1

#include "G4VSensitiveDetector.hh"

#include "MDM_BetaHit.hh"

//#include <vector>

class G4Step;
class G4HCofThisEvent;

/// Beta sensitive detector class
///
/// In Initialize(), it creates one hit for each Beta layer and one more
/// hit for accounting the total quantities in all layers.
///
/// The values are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step.

class MDM_BetaSD : public G4VSensitiveDetector
{
  public:
    MDM_BetaSD(const G4String& name, 
                     const G4String& hitsCollectionName, 
                     G4int nofCells);
    virtual ~MDM_BetaSD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

  private:
    MDM_BetaHitsCollection* fBetaHitsCollection;
    G4int     fNofCells;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

