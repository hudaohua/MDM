//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// ********************************************************************
//
// $Id: MDM_ActionInitialization.hh 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file MDM_ActionInitialization.hh
/// \brief Definition of the MDM_ActionInitialization class

#ifndef MDM_ActionInitialization_h
#define MDM_ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
/* 21/11/15	removed the SiPM module. //
#include <string>
 */

class MDM_DetectorConstruction;
/// Action initialization class.

class MDM_ActionInitialization : public G4VUserActionInitialization
{
  public:
    MDM_ActionInitialization(MDM_DetectorConstruction*);
    virtual ~MDM_ActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

private:
	MDM_DetectorConstruction* fDetConstruction;
	/* 21/11/15	removed the SiPM module. //
	std::string filename;
	 */


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
