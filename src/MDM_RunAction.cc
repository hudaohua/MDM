//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *

// ********************************************************************
//
// $Id: MDM_RunAction.cc 75216 2013-10-29 16:08:11Z gcosmo $
//
/// \file MDM_RunAction.cc
/// \brief Implementation of the MDM_RunAction class

#include "MDM_RunAction.hh"
//#include "MDM_PrimaryGeneratorAction.hh"
//#include "MDM_DetectorConstruction.hh"
#include "MDM_Run.hh"
#include "MDM_Analysis.hh"

#include "G4Timer.hh"
//#include "G4RunManager.hh"
//#include "G4Run.hh"
//#include "G4LogicalVolumeStore.hh"
//#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <string>
#include <sstream>





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


MDM_RunAction::MDM_RunAction()
: G4UserRunAction(), fTimer(0)
{ 

  fTimer = new G4Timer;



  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in MDM_Analysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  analysisManager->SetVerboseLevel(0);
 // analysisManager->SetFirstHistoId(1);

  // Book histograms, ntuple
  //
  
  // Creating histograms
  
  analysisManager->CreateH1("h1_1","Total photon energy in SiPm", 100000, 0.0, 100*keV,"eV");  // max energy 10KeV
  analysisManager->CreateH1("h1_2","Single Photon energy in sipm", 1000, 0.0, 20*eV,"eV");
  analysisManager->CreateH1("h1_3","Highest track ID", 1000000, 0., 1E8);  // max counts 10,000,000
  analysisManager->CreateH1("h1_4","Gamma Photon hit time in sipm", 10000, 0., 200*ns, "ns");
  analysisManager->CreateH1("h1_5","Total Photon hit time in sipm", 10000, 0., 200*ns, "ns");
  analysisManager->CreateH1("h1_6","GPS energy distribution",10000,0,100*MeV,"MeV"); //resolution 100keV
  analysisManager->CreateH2("h2_1","Photon hit position in sipm",1000,-2.0*mm,2.0*mm,1000,-2.0*mm,2.0*mm,"mm","mm");
  analysisManager->CreateH2("h2_2","Envery vs Photons for Beta detector", 1000, 0.01*MeV, 100*MeV,1000,0.,1E5,"MeV");
  analysisManager->CreateH2("h2_3","Envery vs Photons for Gamma detector", 1000, 0.01*MeV, 100*MeV,1000,0.,1E8,"MeV");
  // not supported at Geant4.10.00 sp2, supported after Geant4.10.01
  analysisManager->CreateH3("h3_1","Primary Position",100,-40*mm,+40*mm,100,-40*mm,+40*mm,100,-40*mm,+40*mm,"mm","mm","mm");
  analysisManager->CreateH3("h3_2","Instrument Angluar Response",18,0,pi,36,-pi,pi,1000,0,1E8,"deg","deg");
  

  // creating ntuple
  //
  analysisManager->SetFirstNtupleId(0);
  analysisManager->CreateNtuple("Name-MDM","Title-Data");
  analysisManager->CreateNtupleIColumn(0,"Photon_in_SiPM");   // column Id = 0
  analysisManager->CreateNtupleIColumn(0,"Photon_from_Beta_Det"); // column Id = 1
  analysisManager->CreateNtupleIColumn(0,"Photon_from_Gamma_Det"); // column Id = 2
  analysisManager->CreateNtupleDColumn(0,"Edep_in_Beta_Det"); // column Id = 3
  analysisManager->CreateNtupleDColumn(0,"Edep_in_Gamma_Det"); // column Id = 4
  analysisManager->CreateNtupleDColumn(0,"Total_energy_of_Photon_collected_by_SiPM"); // column Id = 5
  analysisManager->CreateNtupleDColumn(0,"Primary_Energy"); // column Id = 6
  analysisManager->CreateNtupleDColumn(0,"Primary_pos_x");
  analysisManager->CreateNtupleDColumn(0,"Primary_pos_y");
  analysisManager->CreateNtupleDColumn(0,"Primary_pos_z");
  analysisManager->CreateNtupleDColumn(0,"Primary_dir_theta");
  analysisManager->CreateNtupleDColumn(0,"Primary_dir_phi");
  analysisManager->CreateNtupleSColumn(0,"Primary_type");
  analysisManager->CreateNtupleIColumn(0,"No_of_Hits_in_Beta");
  analysisManager->CreateNtupleIColumn(0,"No_of_Hits_in_Gamma");
  analysisManager->CreateNtupleDColumn(0,"Det_Efficiency");
  analysisManager->FinishNtuple(0);





}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_RunAction::~MDM_RunAction()
{
	delete fTimer;
	delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* MDM_RunAction::GenerateRun()
{
  return new MDM_Run; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  //inform the runManager to save random number seed
//	G4RunManager::GetRunManager()->SetRandomNumberStore(true;

	  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  // 
  analysisManager->SetNtupleMerging(true);

  // Open an output file
  //  
  // what a pain in ass to get the G4int to G4string done !! there must be a better way1!
  
	int run_ID = (int)(aRun->GetRunID());
	G4String fileName = "MDM_run";
	std::ostringstream temp;
	temp << run_ID;
	fileName += (G4String)temp.str();
    
	analysisManager->OpenFile(fileName);
 
	G4cout <<"==========================Run"<< run_ID << " start======================"<< G4endl;
	fTimer->Start();
	

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void MDM_RunAction::EndOfRunAction(const G4Run* run)
{

	fTimer->Stop();

	G4int nofEvents = run->GetNumberOfEvent();
	  //if (nofEvents == 0) return;

	  // print histogram statistics
	  //
	G4cout<<"number of event = "<<nofEvents<<" "<< *fTimer << G4endl;
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->Write();
	analysisManager->CloseFile();

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
