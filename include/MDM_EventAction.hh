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

/// Event action class
///

class MDM_EventAction : public G4UserEventAction
{
  public:
    MDM_EventAction();
    virtual ~MDM_EventAction();
    
    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    void AddEdep(G4double edep) { fEdep += edep; }
	void AddAbs(G4double de, G4double dl);
    void AddGap(G4double de, G4double dl);

private:
	MDM_CalorHitsCollection* GetHitsCollection(G4int hcID, const G4Event* event) const;
	MDM_SipmHitsCollection* GetHitsCollection_2(G4int hcID, const G4Event* event) const;
	void PrintEventStatistics(G4double scintEdep, G4double scintTrackLength, G4double sipmPhotons, G4double sipmTrackID) const;

  private:
    G4double  fEdep;
	G4double  fEnergyAbs;
    G4double  fEnergyGap;
    G4double  fTrackLAbs; 
    G4double  fTrackLGap;

	//data member
	G4int fScintHCID;
	G4int fSipmHCID;
};

// inline functions

inline void MDM_EventAction::AddAbs(G4double de, G4double dl) {
  fEnergyAbs += de; 
  fTrackLAbs += dl;
}

inline void MDM_EventAction::AddGap(G4double de, G4double dl) {
  fEnergyGap += de; 
  fTrackLGap += dl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
