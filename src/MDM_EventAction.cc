//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// ********************************************************************
//
// $Id: MDM_EventAction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file MDM_EventAction.cc
/// \brief Implementation of the MDM_EventAction class

#include "MDM_EventAction.hh"
//#include "MDM_CalorimeterSD.hh"
#include "MDM_CalorHit.hh"
#include "MDM_SipmHit.hh"
//#include "MDM_SipmSD.hh"
#include "MDM_BetaHit.hh"
//#include "MDM_BetaSD.hh"

//#include "MDM_Run.hh"

//#include "MDM_RunAction.hh"
#include "MDM_Analysis.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Timer.hh"
#include "G4ios.hh"
#include "MDM_Global.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_EventAction::MDM_EventAction()
: G4UserEventAction(),
  fScintHCID(-1),
  fSipmHCID(-1),
  fBetaHCID(-1),
  fTimer(0)

{
  // set printing event number per each event
   G4RunManager::GetRunManager()->SetPrintProgress(0);     
   fTimer = new G4Timer;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_EventAction::~MDM_EventAction()
{	
	delete fTimer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_EventAction::BeginOfEventAction(const G4Event* event)
{    
	fTimer->Start();

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  //getcollectionID(() method is a heavy operation, it should not be invoked for every event
  if(fSipmHCID<0)
  {
	  fSipmHCID = SDman->GetCollectionID("sipmHitsCollection");
	  G4cout << "Collection ID prints - sipmHCID: " << fSipmHCID << G4endl;
	  fScintHCID = SDman->GetCollectionID("scintHitsCollection");	
	  G4cout << "Collection ID prints - scintHCID: " << fScintHCID << G4endl;
#ifdef ENABLE_BETASD_GAMMA_SD
	    
	  fBetaHCID = SDman->GetCollectionID("betaHitsCollection");
	  
	  G4cout << "Collection ID prints - betaHCID: " << fBetaHCID << G4endl;
#endif
  
  }
#ifdef OUTPUT_DATA_FORMAT
  G4cout<<event->GetEventID()<<",";
#else
  G4cout<<"=============Event "<<event->GetEventID() <<" starts=================="<<G4endl;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_EventAction::EndOfEventAction(const G4Event* event)
{   
	// do nothing if the hit collection is not found
#ifdef ENABLE_BETASD_GAMMA_SD
	if(fScintHCID<0 || fSipmHCID<0 || fBetaHCID<0) return;
#else
	if(fScintHCID<0 ||fSipmHCID<0) return;
#endif
	// get hits collections

	MDM_SipmHitsCollection* sipmHC = GetHitsCollection_2(fSipmHCID,event);
	MDM_CalorHitsCollection* scintHC = GetHitsCollection(fScintHCID,event);
#ifdef ENABLE_BETASD_GAMMA_SD
	MDM_BetaHitsCollection* betaHC = GetHitsCollection_3(fBetaHCID,event);
#endif

#ifdef ENABLE_BETASD_GAMMA_SD
	if ( (!scintHC) || (!sipmHC) || (!betaHC) ) 
    {
        G4ExceptionDescription msg;
        msg << "Some of hits collections of this event not found." << G4endl; 
        G4Exception("B5EventAction::EndOfEventAction()",
                    "B5Code001", JustWarning, msg);
        return;
    }
#else
	if ( (!scintHC) || (!sipmHC) ) 
    {
        G4ExceptionDescription msg;
        msg << "Some of hits collections of this event not found." << G4endl; 
        G4Exception("B5EventAction::EndOfEventAction()",
                    "B5Code001", JustWarning, msg);
        return;
    }
#endif


	// get hit with total value
	// declare all the in function parameters here
	G4int n_scint_hit=0;
	G4int n_sipm_hit=0;
	G4int n_beta_hit =0;
	G4int n_gamma_detector_photons=0;
	G4int n_beta_detector_photons=0;

	  // data collection for SiPm_hit
	  int i;
	  G4double totalTrackEnergy = 0.0;
	  G4int largestTrackID = 0;
	  G4int local_track_id = 0;
	  G4String VertexvolumnName;

	n_scint_hit = scintHC->entries();
	MDM_CalorHit* scintHit = new MDM_CalorHit();
		if(n_scint_hit>0)
	{
		scintHit = (*scintHC)[n_scint_hit-1];
	}

#ifdef ENABLE_BETASD_GAMMA_SD
	//if(betaHC!=NULL)
		n_beta_hit = betaHC->entries();
	//if(scintHC!=NULL)
		
	//if(sipmHC!=NULL)
	MDM_BetaHit* betaHit = new MDM_BetaHit();
	
	if(n_beta_hit>0)
	{
		betaHit = (*betaHC)[n_beta_hit-1];
	}	



#endif
	n_sipm_hit = sipmHC->entries();

	MDM_SipmHit* sipmHit = new MDM_SipmHit();

	if(n_sipm_hit>0)
	{
		sipmHit = (*sipmHC)[n_sipm_hit-1];
	}
	#ifdef OUTPUT_DATA_FORMAT

		#ifdef ENABLE_BETASD_GAMMA_SD
			G4cout<<n_beta_hit<<",";
			
		#endif
			G4cout<<n_scint_hit<<",";
			G4cout<<n_sipm_hit<<",";
	#else
		#ifdef ENABLE_BETASD_GAMMA_SD
			G4cout << "size of betaHit: " << sizeof(*betaHit) << ",n_beta_hit="<<n_scint_hit << G4endl;
			
		#endif
			G4cout << "size of scintHit: " << sizeof(*scintHit) << ",n_scint_hit="<<n_scint_hit<<G4endl;
			G4cout << " size of sipmHit: " << sizeof(*sipmHit) << ",n_sipm_hit="<<n_sipm_hit << G4endl;
	#endif	
  // get analysis manager

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(0);


  // comment out for testing the highest processing speed
  if(n_sipm_hit>0)
  {
	  // statistics per hit
	for(i=0;i<n_sipm_hit;i++)
	{
		sipmHit = (*sipmHC)[i];
		analysisManager->FillH1(1, sipmHit->GetTrackEnergy()); //single photon energy
		G4ThreeVector pos = sipmHit->GetPos();
		analysisManager->FillH2(0,pos.getX()/mm,pos.getY()/mm,1.0); //photon hit position		


		totalTrackEnergy += sipmHit->GetTrackEnergy();
		local_track_id = sipmHit->GetTrackID();  
		if(local_track_id>largestTrackID)
		{
			largestTrackID = local_track_id;
		}

		VertexvolumnName = sipmHit->GetVertexVolumnName();

		if(VertexvolumnName=="Detector")                               //(VertexvolumnName=="Detector")
		{
			//photons generated in the gamma detector
			analysisManager->FillH1(3, sipmHit->GetHitTime());
			n_gamma_detector_photons ++;
		}else if(VertexvolumnName=="Beta_Detector")                         //(VertexvolumnName=="Beta_Detector")
		{
			//photons generated in the beta detector
			//comment out 05/06/2022
			//analysisManager->FillH1(4, sipmHit->GetHitTime());
			n_beta_detector_photons ++;
		}
		// store all the particle into H4
		analysisManager->FillH1(4, sipmHit->GetHitTime());
	}
	analysisManager->FillH1(0, totalTrackEnergy);  //total energy of all the photons detected by the SiPm
	analysisManager->FillH1(2, largestTrackID);//Largest track ID collected by the SiPM

 }
 
  //data collection for Primary
 G4ThreeVector PrimaryPos=event->GetPrimaryVertex(0)->GetPosition();
 G4double pos_x = PrimaryPos.getX();
 G4double pos_y = PrimaryPos.getY();
 G4double pos_z = PrimaryPos.getZ();
 //G4double theta = PrimaryPos.getTheta();
 //G4double phi = PrimaryPos.getPhi();

 G4PrimaryParticle * Primary = event->GetPrimaryVertex(0)->GetPrimary(0);
 G4double PrimaryEne = Primary->GetKineticEnergy();  // 
 G4String PrimaryName = Primary->GetParticleDefinition()->GetParticleName();
 G4double dir_theta = Primary->GetMomentumDirection().getTheta();   // momentum direction x
 G4double dir_phi = Primary->GetMomentumDirection().getPhi();
 //G4double dir_z2 = Primary->GetMomentumDirection().getZ();

 analysisManager->FillH3(0,pos_x/mm,pos_y/mm,pos_z/mm,1.0);
 analysisManager->FillH3(1,dir_theta,dir_phi,n_sipm_hit,1.0);

 analysisManager->FillH1(5,PrimaryEne);

 analysisManager->FillH2(1,PrimaryEne,n_beta_detector_photons,1.0);
 analysisManager->FillH2(2,PrimaryEne,n_gamma_detector_photons,1.0);

 // data collection for beta and gamma

#ifdef ENABLE_BETASD_GAMMA_SD  
 G4double betaTotalEdep=0;
 for(i=0;i<n_beta_hit;i++)
 {
	 betaHit = (*betaHC)[i];
	 betaTotalEdep += betaHit->GetEdep();
 }

 
#endif
 G4double gammaTotalEdep=0;
 for(i=0;i<n_scint_hit;i++)
 {
	 scintHit=(*scintHC)[i];
	 gammaTotalEdep += scintHit->GetEdep();
 }

 // calculate detector efficiency
 G4double Det_Efficiency = (n_beta_detector_photons/500.0+n_gamma_detector_photons/32000.0)/PrimaryEne;
 // Fill up the Ntuple
 analysisManager->FillNtupleIColumn(0,0, n_sipm_hit);
 analysisManager->FillNtupleIColumn(0,1, n_beta_detector_photons);
 analysisManager->FillNtupleIColumn(0,2, n_gamma_detector_photons);
#ifdef ENABLE_BETASD_GAMMA_SD  
 analysisManager->FillNtupleDColumn(0,3, betaTotalEdep);  // even 0 keV is acceptable
#endif
 analysisManager->FillNtupleDColumn(0,4, gammaTotalEdep);  // even 0 keV is acceptable
 analysisManager->FillNtupleDColumn(0,5, totalTrackEnergy);
 analysisManager->FillNtupleDColumn(0,6, PrimaryEne);
 analysisManager->FillNtupleDColumn(0,7, pos_x);
 analysisManager->FillNtupleDColumn(0,8, pos_y);
 analysisManager->FillNtupleDColumn(0,9, pos_z);
 analysisManager->FillNtupleDColumn(0,10, dir_theta);
 analysisManager->FillNtupleDColumn(0,11, dir_phi);
 analysisManager->FillNtupleSColumn(0,12, PrimaryName);
 analysisManager->FillNtupleIColumn(0,13, n_beta_hit );
 analysisManager->FillNtupleIColumn(0,14, n_scint_hit);
 analysisManager->FillNtupleDColumn(0,15, Det_Efficiency);
 analysisManager->AddNtupleRow(0);

 
 fTimer->Stop();

 #ifdef OUTPUT_DATA_FORMAT
 G4cout << PrimaryEne <<","

#ifdef ENABLE_BETASD_GAMMA_SD 
		<< betaTotalEdep <<","    //total energy deposit in beta
		<< gammaTotalEdep <<","     //total energy deposit in gamma
#endif
		<< totalTrackEnergy <<","  // total photon energy deposit in the SiPM
		<< n_beta_detector_photons <<","
		<< n_gamma_detector_photons <<","
		<< dir_theta <<","
		<< dir_phi <<","
		<< fTimer->GetUserElapsed()<<","
		<< PrimaryName << G4endl;
#else
 //G4cout<<"===================Event "<<event->GetEventID() <<" ends=================="<<G4endl;
 //event complete, print the statistics information
 
 G4cout << " Primaries generated in this event    : " << event->GetNumberOfPrimaryVertex() << G4endl;
 G4cout << " Primary 0 position                   : " << G4BestUnit(pos_x,"Length") <<","<<G4BestUnit(pos_y,"Length") <<","<< G4BestUnit(pos_z,"Length") << G4endl;
 G4cout << " Primary 0 direction                  : " << G4BestUnit(dir_theta,"Angle") <<","<<G4BestUnit(dir_phi,"Angle")  << G4endl;
 G4cout << " Primary 0 energy                     : " << PrimaryEne << G4endl;
 G4cout << " Primary 0 name                       : " << PrimaryName <<G4endl;
// Primary->Print();
 G4cout << " Total Hits in the gamma detector     : " << n_scint_hit << G4endl;
 G4cout << " Total Hits in the beta detector      : " << n_beta_hit << G4endl;
 G4cout << " Total Photons collected sipm sensor  : " << n_sipm_hit << G4endl;
 G4cout << " Total Photons from the beta detector : " << n_beta_detector_photons << G4endl;
 G4cout << " Total Photons from the gamma detector: " << n_gamma_detector_photons << G4endl;
 G4cout << " Total energy deposit in beta detector: " << betaTotalEdep << G4endl;
 G4cout << " Total energy deposit in gamma detector: " <<gammaTotalEdep << G4endl;
 G4cout << " Total photons energy in sipm         : " << totalTrackEnergy << G4endl;
 G4cout << " Largest track ID                     : " << largestTrackID << G4endl;
 G4cout << " Processing speed                     : " << n_sipm_hit/fTimer->GetUserElapsed() <<G4endl;
 G4cout << " Detector efficiency                  : " << Det_Efficiency*100.0<<" %"<<G4endl;
 // PrintEventStatistics(scintHit->GetTotalEdep(),betaHit->GetTotalEdep(),totalTrackEnergy,n_gamma_detector_photons,n_beta_detector_photons);
 G4cout << "This event runs for: " << *fTimer << G4endl;
 
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_CalorHitsCollection* MDM_EventAction::GetHitsCollection(G4int hcID, const G4Event* event) const
{
	MDM_CalorHitsCollection* hitsCollection = static_cast<MDM_CalorHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));

	if(!hitsCollection){
		G4ExceptionDescription msg;
		msg<< "Canot access hitsCollection ID" << hcID;
		G4Exception("MDM_EventAction::GetHitsCollection()","MyCode0003", FatalException, msg);
	}

	return hitsCollection;
}

MDM_SipmHitsCollection* MDM_EventAction::GetHitsCollection_2(G4int hcID, const G4Event* event) const
{
	MDM_SipmHitsCollection* hitsCollection = static_cast<MDM_SipmHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));

	if(!hitsCollection){
		G4ExceptionDescription msg;
		msg<< "Canot access hitsCollection ID" << hcID;
		G4Exception("MDM_EventAction::GetHitsCollection()","MyCode0004", FatalException, msg);
	}

	return hitsCollection;
}

MDM_BetaHitsCollection* MDM_EventAction::GetHitsCollection_3(G4int hcID, const G4Event* event) const
{
	MDM_BetaHitsCollection* hitsCollection = static_cast<MDM_BetaHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));

	if(!hitsCollection){
		G4ExceptionDescription msg;
		msg<< "Canot access hitsCollection ID" << hcID;
		G4Exception("MDM_EventAction::GetHitsCollection()","MyCode0005", FatalException, msg);
	}

	return hitsCollection;
}

void MDM_EventAction::PrintEventStatistics(G4double scintEdep,G4double betaEdep, G4double scintTrackLength, G4int sipmGammaPhotons, G4int sipmBetaPhotons) const
{
	// print event statistics
	G4cout
		<<"Scintillator: total energy: "
		<<std::setw(5) << G4BestUnit(scintEdep, "Energy")
		<<"	Beta: total energy: "
		<<std::setw(5) << G4BestUnit(betaEdep, "Energy")
		<<"	SiPM: total energy: "
		<<std::setw(5) << G4BestUnit(scintTrackLength, "Energy")
		<<G4endl
		<<"	Gamma total photons: "
		<<std::setw(5) << sipmGammaPhotons
		<<"	Beta total photons: "
		<<std::setw(5) << sipmBetaPhotons
		<<G4endl;

}
