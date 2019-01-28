//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// ********************************************************************
//
// $Id: MDM_RunAction.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file MDM_RunAction.hh
/// \brief Definition of the MDM_RunAction class

#ifndef MDM_RunAction_h
#define MDM_RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

/* 21/11/15	removed the SiPM module. //
//from sample code
#include <string>
*/

class G4Timer;
class G4Run;
class G4LogicalVolume;

/* 21/11/15	removed the SiPM module. //
//from sample code
class PersistencyHandler;
*/

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume 
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

class MDM_RunAction : public G4UserRunAction
{
  public:
    MDM_RunAction();
    virtual ~MDM_RunAction();

    virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    /* 21/11/15	removed the SiPM module. //
    std::string getFilename() const;

    PersistencyHandler* getPersistencyHandler() const;
*/

  private:
	G4Timer* fTimer;

	/* 21/11/15	removed the SiPM module. //
	//from sample code
	std::string filename;
	PersistencyHandler* persistencyHandler;
	*/

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

