//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: MDM_DetectorConstruction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file MDM_DetectorConstruction.cc
/// \brief Implementation of the MDM_DetectorConstruction class

#include "MDM_DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include <G4ProductionCuts.hh>
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4SDManager.hh"
#include "MDM_CalorimeterSD.hh"
#include "MDM_SipmSD.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4OpBoundaryProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_DetectorConstruction::MDM_DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MDM_DetectorConstruction::~MDM_DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* MDM_DetectorConstruction::Construct()
{  
	G4double density;
	G4double temperature, pressure;
	G4int ncomponents;
	G4String name;
	G4double a, z;
	
	// Get nist material manager
	G4NistManager* nist = G4NistManager::Instance();

	/*********************************************************************************************************************************************
	* Define material Mylar and its optical property 
	*********************************************************************************************************************************************/
	G4Material* Mylar = nist->FindOrBuildMaterial("G4_MYLAR");
	//optical property of mylar
	G4MaterialPropertiesTable* mptMylar = new G4MaterialPropertiesTable();
	G4double AbsorptionLengthMylar = 100.*cm;
	G4double refractiveIndexMylar = 1.640;
	mptMylar->AddConstProperty("RINDEX",refractiveIndexMylar);
	mptMylar->AddConstProperty("ABSLENGTH",AbsorptionLengthMylar);
	Mylar->SetMaterialPropertiesTable(mptMylar);

	
	/*********************************************************************************************************************************************
	* Define material Silicon
	*********************************************************************************************************************************************/
	G4Material* Silicon = nist->FindOrBuildMaterial("G4_Si");
	
	/*********************************************************************************************************************************************
	* Define material vacuum and its optical property 
	*********************************************************************************************************************************************/
	density     = universe_mean_density;                //from PhysicalConstants.h
	pressure    = 3.e-18*pascal;
	temperature = 2.73*kelvin;
	G4Material* vacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole, density,kStateGas,temperature,pressure);

	const G4int nEntries_vacc = 5;
	G4double photonEnergies_vacc[nEntries_vacc] = { 2.761*eV,2.485*eV,2.258*eV,2.071*eV,1.911*eV};
	//optical property of vacuum
	G4double refractiveIndexVacuum[nEntries_vacc] = {1.00, 1.00, 1.00, 1.00, 1.00};
	G4MaterialPropertiesTable* mptVacuum = new G4MaterialPropertiesTable();
	mptVacuum->AddProperty("RINDEX",photonEnergies_vacc,refractiveIndexVacuum,nEntries_vacc);
	vacuum->SetMaterialPropertiesTable(mptVacuum);

	/*********************************************************************************************************************************************
	* Define material CsI(Tl) and its optical property 
	*********************************************************************************************************************************************/
	G4Material* CsI = nist->FindOrBuildMaterial("G4_CESIUM_IODIDE");
	G4Element* Tl = nist->FindOrBuildElement("Tl");
	G4Material* CsI_Tl = new G4Material("CsI_Tl", density=4.51*g/cm3,ncomponents=2);
	CsI_Tl->AddMaterial(CsI,99.6*perCent);
	CsI_Tl->AddElement(Tl,0.4*perCent); 

	/************************************************************************************************************************************************* 
	defination of scintillator material refer to example 5.6
	CsI crystal property refers to http://www.crystals.saint-gobain.com/CsI%28Tl%29_scintillator.aspx
	the property term definition refers to http://wiki.opengatecollaboration.org/index.php/Users_Guide_V6.2:Generating_and_tracking_optical_photons
	another property sources in geant4 code: http://hypernews.slac.stanford.edu/HyperNews/geant4/get/AUX/2011/06/01/03.27-76314-2DetectorConstruction.txt
	define the energy spectrum
	can use as many points as user wants
	**************************************************************************************************************************************************/

	/************************************************************************************************************************************************
	The refracive index can be sourced from: http://refractiveindex.info/?shelf=main&book=CsI&page=Li
	As MPPC from Hamamatsu is only responsible to 300nm to 900nm. The data range between 300nm and 1000nm is taken here.
	**************************************************************************************************************************************************/
	const G4int nEntries = 25;
	//energy values; currently this value corresponds to the maximum emission energy of CsI (550 nm)
	G4double photonEnergies[nEntries] = { 
		4.443*eV,4.202*eV,3.973*eV,3.758*eV,3.553*eV,
		3.360*eV,3.177*eV,3.005*eV,2.841*eV,2.687*eV,
		2.541*eV,2.403*eV,2.271*eV,2.148*eV,2.031*eV,
		1.921*eV,1.817*eV,1.718*eV,1.624*eV,1.536*eV,
		1.452*eV,1.374*eV,1.299*eV,1.228*eV,1.161*eV
	};
	// optical properties of Cesium Iodide
	G4double refractiveIndexCesiumIodide[nEntries] = {
		2.041,1.990,1.951,1.920,1.895,
		1.874,1.857,1.842,1.830,1.819,
		1.810,1.802,1.795,1.789,1.783,
		1.779,1.775,1.772,1.768,1.766,
		1.763,1.761,1.759,1.757,1.756
	};
	/************************************************************************************************************************************************
	This parameters defined the absorption lenth of incident photon at the energy range in CsI(Tl) crystal.
	**************************************************************************************************************************************************/
	G4double absorptionLengthCesiumIodide[nEntries] = {
		39.3*cm,39.3*cm,39.3*cm,39.3*cm,39.3*cm,
		39.3*cm,39.3*cm,39.3*cm,39.3*cm,39.3*cm,
		39.3*cm,39.3*cm,39.3*cm,39.3*cm,39.3*cm,
		39.3*cm,39.3*cm,39.3*cm,39.3*cm,39.3*cm,
		39.3*cm,39.3*cm,39.3*cm,39.3*cm,39.3*cm
	};

	/************************************************************************************************************************************************
	Data extracted from datasheet http://www.crystals.saint-gobain.com/uploadedFiles/SG-Crystals/Documents/CsI%28Tl%29%20and%20%28Na%29%20data%20sheet.pdf
	The emission graph in the datasheet onlys shows effecitve ragne from 300nm to 700nm.
	**************************************************************************************************************************************************/
	G4double Scnt_SLOW[nEntries] = {
		0.00,0.00,0.00,0.00,0.01,
		0.04,0.10,0.14,0.19,0.32,
		0.60,0.86,0.99,0.84,0.65,
		0.39,0.20,0.15,0.01,0.00,
		0.00,0.00,0.00,0.00,0.00
	};

	// defining an object of the G4MaterialPropertiesTable specific to CsI
	G4MaterialPropertiesTable* mptCesiumIodide = new G4MaterialPropertiesTable();
	// add the optical properties of Cesium Iodide to this object
	// mptCesiumIodide->AddProperty("FASTCOMPONENT", photonEnergies, Scnt_FAST, nEntries);
	mptCesiumIodide->AddProperty("SLOWCOMPONENT", photonEnergies, Scnt_SLOW, nEntries);
	mptCesiumIodide->AddProperty("RINDEX",photonEnergies,refractiveIndexCesiumIodide,nEntries);
	mptCesiumIodide->AddProperty("ABSLENGTH",photonEnergies,absorptionLengthCesiumIodide,nEntries);
	mptCesiumIodide->AddConstProperty("SCINTILLATIONYIELD", 54000./MeV);
	// no need for yield ratio if only one compound(fast or slow) is defined.
	// mptCesiumIodide->AddConstProperty("YIELDRATIO", 0.0);
	// mptCesiumIodide->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
	mptCesiumIodide->AddConstProperty("SLOWTIMECONSTANT", 1000.*ns);
	// RESOLUTIONSCALE can be calculated from the energy resolution of the scintillator.
	// here the energy resolution should be the intrinsic energy resolution of the scintillator
	// resolution scale here will produce a statistical fluctuation around the average yield set with scintillationyeild.
	// while values>1 broaden the fluctuation. 0 produces no fluctuation.
	mptCesiumIodide->AddConstProperty("RESOLUTIONSCALE", 1.0);

	// set these to Cesium Iodide
	CsI_Tl->SetMaterialPropertiesTable(mptCesiumIodide);



	/*********************************************************************************************************************************************
	* Define World
	* World defines the boundary of simulation
	*********************************************************************************************************************************************/
	// Option to switch on/off checking of volumes overlaps
	G4bool checkOverlaps = true;

	G4double world_sizeXY = 2*cm;
	G4double world_sizeZ  = 2*cm;
	G4Material* world_mat = vacuum;
  
	G4Box* solidWorld =    
		new G4Box("World",                       //its name
		0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size, actual size is 2x the parameter value
      
	G4LogicalVolume* logicWorld =                         
		new G4LogicalVolume(solidWorld,          //its solid
							vacuum,           //its material
							"World");            //its name
                                   
	G4VPhysicalVolume* physWorld = 
		new G4PVPlacement(0,                     //no rotation
						  G4ThreeVector(),       //at (0,0,0)
						  logicWorld,            //its logical volume
						  "World",               //its name
						  0,                     //its mother  volume
						  false,                 //no boolean operation
						  0,                     //copy number
						  checkOverlaps);        //overlaps checking
  
	G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //white
	simpleBoxVisAtt->SetVisibility(true);
	logicWorld->SetVisAttributes(simpleBoxVisAtt);

	/*********************************************************************************************************************************************
	* Define the main detector - CsI Scintillator shape
	* Scintillator size is a hemisphere with radius 3mm
	*********************************************************************************************************************************************/

	CsI_Tl->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

	// place the detector at the center of the world
	G4ThreeVector detectorPos = G4ThreeVector(0*mm, 0*mm, 0*mm); 
   
	// scintillator shape  
	G4Sphere* sphereDetector=new G4Sphere("Detector",0.*mm,1.5*mm,0*degree,360*degree,0,90*degree);

    // Define the logic volumne                  
	G4LogicalVolume* logicDetector =                         
	new G4LogicalVolume(sphereDetector,         //its shape
                        CsI_Tl,					//its material
                        "Detector");			//its name

	// set the detector colour to Green
	logicDetector->SetVisAttributes(new G4VisAttributes(G4Colour::Green()));

    // Define the physical volumne           
	physDetector = new G4PVPlacement(0,           //no rotation
                    detectorPos,            //at position
                    logicDetector,           //its logical volume
                    "Detector",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

	/*********************************************************************************************************************************************
	* Define the main detector's wrapping material optical property
	*********************************************************************************************************************************************/

	// Define detector optical surface property
	G4OpticalSurface* detectorWrap = new G4OpticalSurface("detectorWrap");
	//G4LogicalBorderSurface* surface1 = new G4LogicalBorderSurface("detectorWrap", physDetector, physWorld,detectorWrap);

	G4double sigma_alpha = 0.1;
	detectorWrap ->SetType (dielectric_metal);
	detectorWrap ->SetFinish(polished); // polished
	detectorWrap ->SetModel(unified);    // glisur was Geant3 model and expired.
	detectorWrap ->SetSigmaAlpha(sigma_alpha);

	const G4int num = 5;
	G4double pp[num] = { 2.761*eV,2.485*eV,2.258*eV,2.071*eV,1.911*eV};
	G4double reflectivity[num] = {0.99,0.99,0.99,0.99,0.99};   // QT throws exception when this is set to 1.0
	G4double efficiency[num] = {0.0, 0.0, 0.0, 0.0, 0.0};
	G4double rindex[num] = {1.56,1.56,1.56,1.56,1.56};

	G4MaterialPropertiesTable* detectorWrapProperty = new G4MaterialPropertiesTable();
	detectorWrapProperty->AddProperty("REFLECTIVITY",pp,reflectivity,num);
	detectorWrapProperty->AddProperty("EFFICIENCY",pp,efficiency,num);	
	detectorWrapProperty->AddProperty("RINDEX",pp,rindex,num);

	detectorWrap->SetMaterialPropertiesTable(detectorWrapProperty);


	/*********************************************************************************************************************************************
	* Define the photon sensor - SiPM shape
	* Use Hamamatsu MPPC, S10362-11-100C
	* ceramic type, sensitive area 1x1mm, 
	* Assume sensitive depth is 100um
	* Assume the material is silicon,
	* Define the physic size as 1mm x 1mm x 1mm  
	* Sensor size is a cube.
	*********************************************************************************************************************************************/
	// Define the position next to the scintillator
	G4ThreeVector sensorPos = G4ThreeVector(0*mm, 0*mm, -0.5*mm); //1.5mm + 0.5mm

	// Define the shape
	G4Box* solidSensor = new G4Box("Sensor",0.5*mm, 0.5*mm, 0.5*mm);  // size is half of the target value

	// Define the logic volumne
	G4LogicalVolume* logicSensor =                         
    new G4LogicalVolume(solidSensor,         //its solid
                        Silicon,			//its material
                        "Sensor");           //its name

	//set shape2 colour to blue
	logicSensor->SetVisAttributes(new G4VisAttributes(G4Colour::Blue()));

    // Define the physics volumne           
	physSensor = new G4PVPlacement(0,                       //no rotation
                    sensorPos,                    //at position
                    logicSensor,             //its logical volume
                    "Sensor",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
	// Define the optical surface between detector and sensor
	G4OpticalSurface* PhotonSensor_op_surf = new G4OpticalSurface("PhotonSensor_op_surf");
//	G4LogicalBorderSurface* surface2 = new G4LogicalBorderSurface("PhotonSensor_op_surf", physDetector, physSensor,PhotonSensor_op_surf);



	sigma_alpha = 0.1;
	PhotonSensor_op_surf ->SetType (dielectric_metal);
	PhotonSensor_op_surf ->SetFinish(polished);
	PhotonSensor_op_surf ->SetModel(unified);  //was glisur
	//PhotonSensor_op_surf ->SetSigmaAlpha(sigma_alpha);

	const G4int num1 = 5;
	G4double pp1[num1] = { 2.761*eV,2.485*eV,2.258*eV,2.071*eV,1.911*eV};
	G4double reflectivity1[num1] = {1.0,1.0,1.0,1.0,1.0};   // 10% reflectivity
	G4double efficiency1[num1] = {0.0, 0.0, 0.0, 0.0, 0.0};
	G4double rindex1[num1] = {1.56,1.56,1.56,1.56,1.56};

	G4MaterialPropertiesTable* PhotonSensor_op_surf_Property = new G4MaterialPropertiesTable();
	PhotonSensor_op_surf_Property->AddProperty("REFLECTIVITY",pp1,reflectivity1,num1);
//	PhotonSensor_op_surf_Property->AddProperty("EFFICIENCY",pp1,efficiency1,num1);	
	PhotonSensor_op_surf_Property->AddProperty("RINDEX",pp1,rindex1,num1);

	PhotonSensor_op_surf_Property->AddProperty("SPECULARLOBECONSTANT", pp1,reflectivity1, num1);   //0
	PhotonSensor_op_surf_Property->AddProperty("SPECULARSPIKECONSTANT", pp1,efficiency1, num1);    //1
	PhotonSensor_op_surf_Property->AddProperty("BACKSCATTERCONSTANT", pp1,reflectivity1, num1);    //0

	PhotonSensor_op_surf->SetMaterialPropertiesTable(PhotonSensor_op_surf_Property);
	
	/* create logical skin surfaces */
	G4LogicalSkinSurface* surface1 = new G4LogicalSkinSurface("detectorWrap",logicDetector,detectorWrap);

	G4LogicalSkinSurface* surface2 = new G4LogicalSkinSurface("PhotonSensor_op_surf",logicSensor,PhotonSensor_op_surf);


  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MDM_DetectorConstruction::ConstructSDandField()
{
	// sensitive detector
	MDM_CalorimeterSD* scintSD = new MDM_CalorimeterSD("Scintillator","scintHitsCollection",1);
	// register the sensitive detector with manager
	G4SDManager::GetSDMpointer()->AddNewDetector(scintSD);
	// attach the sensitive detector to the real detector logical volume
	SetSensitiveDetector("Detector",scintSD);	

	MDM_SipmSD* sipmSD = new MDM_SipmSD("SiPM","sipmHitsCollection",1);
	// register the sensitive detector with manager
	G4SDManager::GetSDMpointer()->AddNewDetector(sipmSD);
	SetSensitiveDetector("Sensor",sipmSD);

}