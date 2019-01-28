//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// ********************************************************************

//
// $Id: MDM_EventAction.hh 75216 2013-10-29 16:08:11Z gcosmo $
//
/// \file MDM_EventAction.hh
/// \brief Definition of the MDM_EventAction class

#ifndef MDM_EventAction_h
#define MDM_EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "MDM_CalorHit.hh"
#include "MDM_SipmHit.hh"
#include "MDM_BetaHit.hh"
#include "G4Timer.hh"

/// Event action class
///

class MDM_EventAction : public G4UserEventAction
{
  public:
    MDM_EventAction();
    virtual ~MDM_EventAction();
    
    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

private:
	MDM_CalorHitsCollection* GetHitsCollection(G4int hcID, const G4Event* event) const;
	MDM_SipmHitsCollection* GetHitsCollection_2(G4int hcID, const G4Event* event) const;
	MDM_BetaHitsCollection* GetHitsCollection_3(G4int hcID, const G4Event* event) const;
	void PrintEventStatistics(G4double scintEdep, G4double betaEdep,G4double scintTrackLength, G4int sipmGammaPhotons, G4int sipmBetaPhotons) const;

  private:
	G4Timer* fTimer;

	//data member
	G4int fScintHCID;
	G4int fSipmHCID;
	G4int fBetaHCID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
