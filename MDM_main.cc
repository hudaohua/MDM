//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// * v0.2	02/06/2015	optical surface configuration code updated.	  *
// * v0.3	30/06/2015	optical property parameters update.			  *
// * v0.4	24/11/2015  removed the G4SiPm package.					  *
// * v0.5	17/02/2016  update the geometry to introduce scintillator *
// *					array, dual layer and SiPm array.			  *
// * v0.6a      09/01/2018   update the scintillator material to LXSR  *
//                           update the program to multiple thread
//                           recomplie the code to geant4.10.03 patch3
// ********************************************************************
//
// $Id: MDM_main.cc 75216 2013-10-29 16:08:11Z gcosmo $
// version log
/// v0.3 30/06/2015


#include "MDM_DetectorConstruction.hh"
#include "MDM_ActionInitialization.hh"
#include "MDM_PhysicsList.hh"  //user defined physics list


#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "QBBC.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "Randomize.hh"
#include "G4ScoringManager.hh"

//#include "G4GDMLParser.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  G4cout << "Title: MDM Simulation Software" << G4endl
	     << "Author: Hubert Hu."<< G4endl
		 << "Version: 0.6a" <<G4endl
		 << "Date: 09/01/2018" <<G4endl;

  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
 
  // Construct the default run manager
  //
  // no multithreded for now
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  runManager -> SetNumberOfThreads(12); // if not set, default to be 2.
  G4cout << "Debug info: I'm running in multiple thread mode!" << G4endl;
#else
  G4RunManager* runManager = new G4RunManager;
  G4cout << "Debug info: I'm running in single thread mode!" << G4endl;
#endif

  // call this scroing manager immediately after the g4runmanager
  // so that the scoring UI can be used.
  // need to explore how to hard code the scoring function

 // G4ScoringManager::GetScoringManager();

  // Set mandatory initialization classes
  //
  // Detector construction
  MDM_DetectorConstruction* detConstruction = new MDM_DetectorConstruction();
  runManager->SetUserInitialization(detConstruction);

  // Physics list
  // use Customized PhysicsList
 // G4VModularPhysicsList* physicsList = new QBBC;
 // physicsList->SetVerboseLevel(1);
  runManager->SetUserInitialization(new MDM_PhysicsList);
    
  // User action initialization
  runManager->SetUserInitialization(new MDM_ActionInitialization(detConstruction));

  // Initialize G4 kernel
  //
  runManager->Initialize();

  	//write to GDML file
	//G4VPhysicalVolume* W = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume();
	//G4GDMLParser parser;
	//parser.Write("MDM_output.gdml",W,true);
  
#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
    
  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  //enable tracking verbose to level 2
  //UImanager->ApplyCommand("/tracking/verbose 1");
 
  if (argc!=1) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode : define UI session
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute init_vis.mac");
#else
    UImanager->ApplyCommand("/control/execute init.mac"); 
#endif
    ui->SessionStart();
    delete ui;
#endif
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !
  
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
