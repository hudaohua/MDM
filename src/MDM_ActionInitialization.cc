//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *
//
// * version log history
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// ********************************************************************
//
// $Id: MDM_ActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file MDM_ActionInitialization.cc
/// \brief Implementation of the MDM_ActionInitialization class

#include "MDM_ActionInitialization.hh"
#include "MDM_PrimaryGeneratorAction.hh"
#include "MDM_RunAction.hh"
#include "MDM_EventAction.hh"
#include "MDM_SteppingAction.hh"

/* 21/11/15	removed the SiPM module. //
// copy from G4Sipm sample code
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
*/

#include <G4RunManager.hh>
#include <G4GeometryTolerance.hh>

/* 21/11/15	removed the SiPM module. //
#include "G4SipmUiMessenger.hh"
*/



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_ActionInitialization::MDM_ActionInitialization(MDM_DetectorConstruction* detConstruction)
 : G4VUserActionInitialization(),
 fDetConstruction(detConstruction)
{
	/* 21/11/15	removed the SiPM module. //
	// from SiPm package sample code
	// create the data file
	if (boost::filesystem::is_directory(path)) {
		std::string time = boost::posix_time::to_iso_string(boost::posix_time::microsec_clock::local_time());
		filename = boost::str(boost::format("%s/g4sipm_%s") % path % time);
	} else {
		filename = path;
	}
	*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_ActionInitialization::~MDM_ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_ActionInitialization::BuildForMaster() const
{
  SetUserAction(new MDM_RunAction());
  G4cout << "Debug info: I'm the master thread!" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_ActionInitialization::Build() const
{
	/* 21/11/15	removed the SiPM module. //
	//copy from sample code
	// Load messengers.
	G4SipmUiMessenger::getInstance();
	*/

//	ParticleSourceMessenger::getInstance();
//	PersistencyHandlerMessenger::getInstance();


  SetUserAction(new MDM_PrimaryGeneratorAction);
  SetUserAction(new MDM_RunAction());

  MDM_EventAction* eventAction = new MDM_EventAction();
  SetUserAction(eventAction);
  SetUserAction(new MDM_SteppingAction(fDetConstruction,eventAction));

  
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
