//
// **********************************************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history																		*
// * v0.1	08/05/2015	convert from old test program to MDM project							*
// * v0.2	02/06/2015	Code corrected for proper optical surface performance.					*
// *					The silicon detector now can detect the photons properly.				*
// *					The previous mistake is due to use the G4LogicalBoarderSurface class	*
// *					instead of G4LogicalSkinSurface class to define the surface property.	*
// * v0.5   12/10/2016	Add a alumni shielding, surrounding the SiPM							 *
// * v0.6   09/04/2017  Add 100% reflectivity to Si       -- cancelled				             *
// * v0.7   05/12/2017  Change the main scintillator to LXSR,								     *
// *                    change the wrapping material to 50nm aluminium                           *
// *                    change the Bc404 scintillator geometry to cylinder, D3mm, H1mm           *
// * v0.8	27/01/2019	update the main scintillator and beta scintillator to a new size,        *
// *					total: D6mm,H3.5mm														 *
// *                    beta scintillator: D6mm, H1mm,											 *
// *					main scintillator: D6mm, H2.5mm,                                         *
// * v0.9	26/05/2022	beta scintillator: BC404, D6MM, H15 um,                                  *
// *                    main scintillator: LXSR, D6mm, H4mm                                      *
// ***********************************************************************************************
//

#include "MDM_DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
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
#include "MDM_BetaSD.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4OpBoundaryProcess.hh"

#include "MDM_Global.hh"


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

	G4double DET_R = 3*mm;
	G4double BETA_DET_HZ = 0.017*mm; // estimated thickness of BC404 coating with +- 5um error)
	G4double GAMMA_DET_HZ = 2*mm;
	
	// Get nist material manager
	G4NistManager* nist = G4NistManager::Instance();

	/*********************************************************************************************************************************************
	* Define material Bc404 and its optical property
	* data copied from https://arxiv.org/pdf/2106.04734.pdf
	* data also from http://www.crystals.saint-gobain.com/uploadedFiles/SG-Crystals/Documents/SGC%20BC400-404-408-412-416%20Data%20Sheet.pdf

	*********************************************************************************************************************************************/
	G4Material* Bc404 = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
	
	const G4int NUMENTRIES = 12;
	G4double PhotonEnergy_Bc404[NUMENTRIES] =	{ 
		3.44*eV, 3.26*eV, 3.1*eV, 3.02*eV, 2.95*eV,
		2.92*eV, 2.82*eV, 2.76*eV, 2.7*eV, 2.58*eV,
		2.38*eV, 2.08*eV };
	G4double RINDEX_Bc404[NUMENTRIES] =	{
		1.58, 1.58, 1.58, 1.58, 1.58,
		1.58, 1.58, 1.58, 1.58, 1.58,
		1.58, 1.58 };
	G4double ABSORPTION_Bc404[NUMENTRIES] =	{
		140*cm, 140*cm, 140*cm, 140*cm, 140*cm,
		140*cm, 140*cm, 140*cm, 140*cm, 140*cm,
		140*cm, 140*cm }; 
	G4double SCINTILLATION_Bc404[NUMENTRIES] =	{
		0.04, 0.07, 0.20, 0.49, 0.84,
		1.00, 0.83, 0.55, 0.40, 0.17,
		0.03, 0 };

	G4MaterialPropertiesTable *Bc404_mt = new G4MaterialPropertiesTable();
	Bc404_mt->AddProperty("RINDEX", PhotonEnergy_Bc404, RINDEX_Bc404,  NUMENTRIES);
	Bc404_mt->AddProperty("ABSLENGTH", PhotonEnergy_Bc404, ABSORPTION_Bc404,  NUMENTRIES);
	Bc404_mt->AddProperty("SCINTILLATIONCOMPONENT1", PhotonEnergy_Bc404, SCINTILLATION_Bc404, NUMENTRIES);
	Bc404_mt->AddConstProperty("SCINTILLATIONYIELD",500./MeV); 
	Bc404_mt->AddConstProperty("RESOLUTIONSCALE",1.0);
	Bc404_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 1.8*ns);
	Bc404->SetMaterialPropertiesTable(Bc404_mt);  

	/*********************************************************************************************************************************************
	* Define material LXSR and its optical property from SEMICRO.ORG
	* data also from https://semicro.org/products/luminix-scintillator-lxsr-new

	*********************************************************************************************************************************************/
	G4Element* Yttrium = nist->FindOrBuildElement("Y"); // G4Element should use the pure element symbol not the material starting with "G4", e.g. G4_Y
	G4Element* Oxygen = nist->FindOrBuildElement("O");
	G4Element* Si = nist->FindOrBuildElement("Si");
	G4Element* Cerium = nist->FindOrBuildElement("Ce");
	G4Material* LXSR = new G4Material("LXSR",4.50*g/cm3,4); //
	LXSR->AddElement(Yttrium,99.2*perCent);
	LXSR->AddElement(Oxygen,0.4*perCent);
	LXSR->AddElement(Si,0.2*perCent);
	LXSR->AddElement(Cerium,0.2*perCent);

	const G4int LXSR_NUMENTRIES = 21;
	G4double PhotonEnergy_LXSR[LXSR_NUMENTRIES] =	{ 
		4.97*eV, 4.14*eV, 3.55*eV, 3.45*eV, 3.36*eV,
		3.31*eV, 3.27*eV, 3.23*eV, 3.19*eV, 3.15*eV,
		3.11*eV, 2.96*eV, 2.76*eV, 2.70*eV, 2.64*eV, 
		2.59*eV, 2.54*eV, 2.44*eV, 2.30*eV, 2.07*eV, 
		1.91*eV};

	G4double SCINTILLATION_LXSR[LXSR_NUMENTRIES] =	{
		0.00132, 0.00263, 0.00526, 0.13158, 0.26316,
		0.39474, 0.52632, 0.65789, 0.78947, 0.92105,
		1.0,     1.0,     0.92105, 0.78947, 0.65789,
		0.52631, 0.394737, 0.263158, 0.131579, 0.052632,
		0.005263};

	G4double RINDEX_LXSR[LXSR_NUMENTRIES] =	{
		1.8, 1.8, 1.8, 1.8, 1.8,
		1.8, 1.8, 1.8, 1.8, 1.8,
		1.8, 1.8, 1.8, 1.8, 1.8,
		1.8, 1.8, 1.8, 1.8, 1.8,
		1.8};

	G4double ABSLENGTH_LXSR[LXSR_NUMENTRIES] =	{
		39.3*cm,39.3*cm,39.3*cm,39.3*cm,39.3*cm,
		39.3*cm,39.3*cm,39.3*cm,39.3*cm,39.3*cm,
		39.3*cm,39.3*cm,39.3*cm,39.3*cm,39.3*cm,
		39.3*cm,39.3*cm,39.3*cm,39.3*cm,39.3*cm,
		39.3*cm};

	G4MaterialPropertiesTable *LXSR_mt = new G4MaterialPropertiesTable();
	LXSR_mt->AddProperty("RINDEX", PhotonEnergy_LXSR, RINDEX_LXSR, LXSR_NUMENTRIES);
	LXSR_mt->AddProperty("ABSLENGTH",PhotonEnergy_LXSR, ABSLENGTH_LXSR, LXSR_NUMENTRIES );
	LXSR_mt->AddProperty("SCINTILLATIONCOMPONENT1", PhotonEnergy_LXSR, SCINTILLATION_LXSR, LXSR_NUMENTRIES);
	LXSR_mt->AddConstProperty("SCINTILLATIONYIELD",32000./MeV);   //80% of NaI
	LXSR_mt->AddConstProperty("RESOLUTIONSCALE",1.0);
	LXSR_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 42.*ns);
	LXSR->SetMaterialPropertiesTable(LXSR_mt);  

	/*********************************************************************************************************************************************
	* Define material TEFLON and its optical property
	*********************************************************************************************************************************************/

	// use TEFLON material because it has the highest reflection coeffcient 0.99
	G4Material* Wrap = nist->FindOrBuildMaterial("G4_TEFLON");
	//optical property of Wrap_material
	/*
	G4MaterialPropertiesTable* mptWrap = new G4MaterialPropertiesTable();
	G4double AbsorptionLengthWrap = 100.*cm;
	G4double refractiveIndexWrap = 1.35;
	mptWrap->AddConstProperty("RINDEX",refractiveIndexWrap);
	mptWrap->AddConstProperty("REFLECTIVITY",1); // 0.99
	mptWrap->AddConstProperty("ABSLENGTH",AbsorptionLengthWrap);
	Wrap->SetMaterialPropertiesTable(mptWrap);
	*/
	//Only if the materials'optical property rindex is defined, the photon can enter the material//
	
	/*********************************************************************************************************************************************
	* Define material Silicon and its optical property
	*********************************************************************************************************************************************/
	G4Material* Silicon = nist->FindOrBuildMaterial("G4_Si");
	
	const G4int nEntries_si = 5;
	G4double photonEnergies_si[nEntries_si] = { 2.761*eV,2.485*eV,2.258*eV,2.071*eV,1.911*eV};
	//optical property of vacuum
	G4double refractiveIndexSi[nEntries_si] = {1.50, 1.50, 1.50, 1.50, 1.50};
	G4MaterialPropertiesTable* mptSi = new G4MaterialPropertiesTable();
	mptSi->AddProperty("RINDEX",photonEnergies_si,refractiveIndexSi,nEntries_si);
	//mptSi->AddConstProperty("REFLECTIVITY",1); //100%
	Silicon->SetMaterialPropertiesTable(mptSi);

	/*********************************************************************************************************************************************
	* Define material Aluminum and its optical property
	*********************************************************************************************************************************************/
	G4Material* Aluminum = nist->FindOrBuildMaterial("G4_Al");
	G4MaterialPropertiesTable* mptWrap = new G4MaterialPropertiesTable();

	G4double refractiveIndexWrap[LXSR_NUMENTRIES] =	{
		0.1846, 0.2642, 0.3667, 0.3886, 0.4122,
		0.4248, 0.4374, 0.45, 0.4627, 0.4753,
		0.4879, 0.5411, 0.6332, 0.6684, 0.7036,
		0.7393, 0.7759, 0.8492, 0.9727, 1.262,
		1.558};
	G4double reflectivityIndexWrap[LXSR_NUMENTRIES] =	{
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1};
	
	G4double AbsorptionLengthWrap[LXSR_NUMENTRIES] = {
		80.*cm,	80.*cm,	80.*cm,	80.*cm,	80.*cm,
		80.*cm,	80.*cm,	80.*cm,	80.*cm,	80.*cm,
		80.*cm,	80.*cm,	80.*cm,	80.*cm,	80.*cm,
		80.*cm,	80.*cm,	80.*cm,	80.*cm,	80.*cm,
		80.*cm};

	mptWrap->AddProperty("RINDEX",PhotonEnergy_LXSR,refractiveIndexWrap,LXSR_NUMENTRIES);
	mptWrap->AddProperty("REFLECTIVITY",PhotonEnergy_LXSR,reflectivityIndexWrap, LXSR_NUMENTRIES); // 0.99
	mptWrap->AddProperty("ABSLENGTH",PhotonEnergy_LXSR, AbsorptionLengthWrap, LXSR_NUMENTRIES);
	Wrap->SetMaterialPropertiesTable(mptWrap);
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
		0.01,0.01,0.01,0.01,0.01,
		0.04,0.10,0.14,0.19,0.32,
		0.60,0.86,0.99,0.84,0.65,
		0.39,0.20,0.15,0.01,0.01,
		0.01,0.01,0.01,0.01,0.01
	};
	
	// defining an object of the G4MaterialPropertiesTable specific to CsI
	G4MaterialPropertiesTable* mptCesiumIodide = new G4MaterialPropertiesTable();
	// add the optical properties of Cesium Iodide to this object
	// mptCesiumIodide->AddProperty("FASTCOMPONENT", photonEnergies, Scnt_FAST, nEntries);
	mptCesiumIodide->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergies, Scnt_SLOW, nEntries);
	mptCesiumIodide->AddProperty("RINDEX",photonEnergies,refractiveIndexCesiumIodide,nEntries);
	mptCesiumIodide->AddProperty("ABSLENGTH",photonEnergies,absorptionLengthCesiumIodide,nEntries);
	mptCesiumIodide->AddConstProperty("SCINTILLATIONYIELD", 54000./MeV);
	// no need for yield ratio if only one compound(fast or slow) is defined.
	// mptCesiumIodide->AddConstProperty("YIELDRATIO", 0.0);
	// mptCesiumIodide->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
	mptCesiumIodide->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 1000.*ns);
	// RESOLUTIONSCALE can be calculated from the energy resolution of the scintillator.
	// here the energy resolution should be the intrinsic energy resolution of the scintillator
	// resolution scale here will produce a statistical fluctuation around the average yield set with scintillationyeild.
	// while values>1 broaden the fluctuation. 0 produces no fluctuation.
	mptCesiumIodide->AddConstProperty("RESOLUTIONSCALE", 1.0);  // no fluctuatoin.

	// set these to Cesium Iodide
	CsI_Tl->SetMaterialPropertiesTable(mptCesiumIodide);
	

	// Define all the shapes

	/*********************************************************************************************************************************************
	* Define World
	* World defines the boundary of simulation
	*********************************************************************************************************************************************/
	// Option to switch on/off checking of volumes overlaps
	G4bool checkOverlaps = true;

	// define the world size, 
	G4double world_sizeXY = 4*cm;
	G4double world_sizeZ  = 4*cm;
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
	* Define the top coating, 50nm thick
	*********************************************************************************************************************************************/
	G4double top_coat_hz = 0.000025*mm;
	G4ThreeVector coat_top_Pos = G4ThreeVector(0*mm, 0*mm, top_coat_hz); //

	// Define the shape
	//G4Tubs* coat_top = new G4Tubs("coat_top",0*mm,1.5*mm,0.0000025*mm,0*degree,360*degree);  // make the size much bigger than the SiPM
	G4Tubs* coat_top = new G4Tubs("coat_top",0*mm,DET_R,top_coat_hz,0*degree,360*degree);

	// Define the logic volumne
	G4LogicalVolume* logic_coat_top =
    new G4LogicalVolume(coat_top,         //its solid
                        Aluminum,			//its material
                        "coat_top");           //its name

	//set shield_bot colour to white
	logic_coat_top->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));

    // Define the physics volumne
	G4VPhysicalVolume* phys_coat_top;
	phys_coat_top = new G4PVPlacement(0,                       //no rotation
                    coat_top_Pos,                    //at position
                    logic_coat_top,             //its logical volume
                    "coat_top",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

	/*********************************************************************************************************************************************
	* Define the first detector - Bc404 Scintillator shape
	* Scintillator size is a hemisphere with thickness of 0.2mm
	*********************************************************************************************************************************************/
#ifndef REMOVE_BETA_DET
	
	// place the detector at the center of the world
	G4ThreeVector Beta_detectorPos = G4ThreeVector(0*mm, 0*mm, -BETA_DET_HZ); 
   
	// fast scintillator shape  definition
	G4double innerRadius = 0.*mm;
	G4double outerRadius = DET_R;
	G4double hz = BETA_DET_HZ; // half z
	G4double startAngle = 0.*deg;
	G4double spanningAngle = 360.*deg;

	G4Tubs* fastDetectorTube
    = new G4Tubs("Beta_Detector",
                 innerRadius, 
                 outerRadius,
                 hz,
                 startAngle, 
                 spanningAngle);
	

    // Define the logic volumne                  
	G4LogicalVolume* logic_Beta_Detector =                         
	new G4LogicalVolume(fastDetectorTube,   //its shape
                        Bc404,					//its material
                        "Beta_Detector");			//its name

	// set the detector colour to white and transparent
	G4Colour white = G4Color::Brown();
	logic_Beta_Detector->SetVisAttributes(new G4VisAttributes(G4Color(white.GetRed(), white.GetGreen(), white.GetBlue(), 1.0)));

    // Define the physical volume
	physBeta_Detector = new G4PVPlacement(0,           //no rotation
                    Beta_detectorPos,            //at position
                    logic_Beta_Detector,           //its logical volume
                    "Beta_Detector",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

#endif

	/*********************************************************************************************************************************************
	* Define the main detector - CsI Scintillator shape
	* Scintillator size is a hemisphere with radius 3mm
	*********************************************************************************************************************************************/

	//CsI_Tl->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

	// place the detector
	G4ThreeVector detectorPos = G4ThreeVector(0*mm, 0*mm, -(BETA_DET_HZ*2+GAMMA_DET_HZ)); 

	// slow scintillator shape  definition
	innerRadius = 0.*mm;
	outerRadius = DET_R;
	hz = GAMMA_DET_HZ; // half z
	startAngle = 0.*deg;
	spanningAngle = 360.*deg;

	G4Tubs* slowDetectorTube
    = new G4Tubs("Detector",
                 innerRadius, 
                 outerRadius,
                 hz,
                 startAngle, 
                 spanningAngle);

    // Define the logic volumne                  
	G4LogicalVolume* logicDetector =                         
	new G4LogicalVolume(slowDetectorTube,         //its shape
                        LXSR,					//its material
                        "Detector");			//its name

	// set the detector colour to Green
	logicDetector->SetVisAttributes(new G4VisAttributes(G4Colour::Green()));

    // Define the physical volume
	physDetector = new G4PVPlacement(0,           //no rotation
                    detectorPos,            //at position
                    logicDetector,           //its logical volume
                    "Detector",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


	/*********************************************************************************************************************************************
	* Define the detector wrapping/coating with Aluminium
	* wapping thickness is defined at 50nm
	*********************************************************************************************************************************************/
	// place the detector at the center of the world


	G4ThreeVector wrapPos = G4ThreeVector(0*mm, 0*mm, -(BETA_DET_HZ+GAMMA_DET_HZ));
	// wrap shape
	G4Tubs* bundleWrap=new G4Tubs("Wrap",DET_R,DET_R+0.00005*mm,(BETA_DET_HZ+GAMMA_DET_HZ),0*degree,360*degree);

    // Define the logic volumne
	G4LogicalVolume* logicWrap =
	new G4LogicalVolume(bundleWrap,         //its shape
                        Wrap,					//its material
                        "Wrap");			//its name

	// set the wrapping colour to Cyan but transparent
	G4Colour cyan = G4Color::White();
	logicWrap->SetVisAttributes(new G4VisAttributes(G4Color(cyan.GetRed(), cyan.GetGreen(), cyan.GetBlue(), 0.5)));

    // Define the physical volume
	physWrap = new G4PVPlacement(0,           //no rotation
                    wrapPos,            //at position
                    logicWrap,           //its logical volume
                    "Wrap",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


	/*********************************************************************************************************************************************
	* Define the photon sensor - SiPM shape
	* Use Hamamatsu MPPC, S13360-1325CS
	* ceramic type, sensitive area 1.3x1.3mm,
	* Assume sensitive depth is 100um
	* Assume the material is silicon,
	* Define the physic size as 1.3mm x 1.3mm x 1mm
	* Sensor size is a cube.
	*********************************************************************************************************************************************/
	// Define the position next to the scintillator

	G4ThreeVector sensorPos = G4ThreeVector(0*mm, 0*mm, -((BETA_DET_HZ+GAMMA_DET_HZ)*2+0.5*mm)); //

	// Define the shape
	G4Box* solidSensor = new G4Box("Sensor",0.65*mm, 0.65*mm, 0.5*mm);  // size is half of the target value

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

	/*********************************************************************************************************************************************
	* Define the photon sensor - SiPM shielding
	* material: Al
	* thickness: 4mm
	* 4pcs around and 1pcs at the bottom
	*********************************************************************************************************************************************/
	// Define the position next to the scintillator
	G4ThreeVector shield_bot_Pos = G4ThreeVector(0*mm, 0*mm, -((BETA_DET_HZ+GAMMA_DET_HZ)*2+1*mm+2*mm)); //

	// Define the shape
	G4Box* shield_bot = new G4Box("Shield_bot",5*mm, 5*mm, 2*mm);  // make the size much bigger than the SiPM

	// Define the logic volumne
	G4LogicalVolume* logic_shield_bot =
    new G4LogicalVolume(shield_bot,         //its solid
                        Aluminum,			//its material
                        "shield_bot");           //its name

	//set shield_bot colour to white
	logic_shield_bot->SetVisAttributes(new G4VisAttributes(G4Colour::Gray()));

    // Define the physics volumne
	G4VPhysicalVolume* phys_shield_bot;
	phys_shield_bot = new G4PVPlacement(0,                       //no rotation
                    shield_bot_Pos,                    //at position
                    logic_shield_bot,             //its logical volume
                    "shield_bot",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

	/*********************************************************************************************************************************************
	* Define the main detector's wrapping material optical property
	*********************************************************************************************************************************************/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	//NOTE: possible Properties:
//         "RINDEX":			(spectrum (in dependence of the photon energy))		(obligatory property!)
//		defines the refraction index of the material, used for boundary processes, Cerenkov radiation and Rayleigh scattering
//         "ABSLENGTH":			(spectrum (in dependence of the photon energy))
//		defines the absorption length (absorption spectrum) of the material, used for the "normal" absorption of optical photons (default is infinity, i.e. no absorption)
//		("ABSLENGTH" should not be specified for WLS materials, the absorption length for the WLS process is specified by "WLSABSLENGTH")
//         "RAYLEIGH":			(spectrum (in dependence of the photon energy))
//		defines the absorption length of the material, used for the rayleigh scattering of optical photons (default is infinity, i.e. no scattering)
//
//         "SCINTILLATIONYIELD":	(constant value (energy independent))			(obligatory property for scintillator materials!)
//		defines the mean number of photons, emitted per MeV energy deposition in the scintillator material (the real number is Poisson/Gauss distributed)
//		(can also be specified separately for different particles by putting "ELECTRON...", "PROTON...", "DEUTERON...", "TRITON...", "ALPHA...", "ION..." infront of "SCINTILLATIONYIELD")
//		(default is 0, i.e. no scintillation process)
//         "RESOLUTIONSCALE":		(constant value (energy independent))
//		defines the intrinsic resolution of the scintillator material, used for the statistical distribution of the number of generated photons in the scintillation process
//		(values > 1 result in a wider distribution, values < 1 result in a narrower distribution -> 1 should be chosen as default)
//		(default is 0)
//         "FASTCOMPONENT":		(spectrum (in dependence of the photon energy))		(at least one "...COMPONENT" is obligatory for scintillator materials!)
//		defines the emission spectrum of the material, used for the fast scintillation process
//         "SLOWCOMPONENT":		(spectrum (in dependence of the photon energy))		(at least one "...COMPONENT" is obligatory for scintillator materials!)
//		defines the emission spectrum of the material, used for the slow scintillation process
//         "FASTTIMECONSTANT":		(constant value (energy independent))
//		defines the decay time (time between energy deposition and photon emission), used for the fast scintillation process (default is 0)
//         "SLOWTIMECONSTANT":		(constant value (energy independent))
//		defines the decay time (time between energy deposition and photon emission), used for the slow scintillation process (default is 0)
//         "FASTSCINTILLATIONRISETIME":	(constant value (energy independent))
//		defines the rise time (time between the start of the emission and the emission peak), used for the fast scintillation process (default is 0)
//         "SLOWSCINTILLATIONRISETIME":	(constant value (energy independent))
//		defines the rise time (time between the start of the emission and the emission peak), used for the slow scintillation process (default is 0)
//         "YIELDRATIO":		(constant value (energy independent))			(obligatory property for scintillator materials, if both "...COMPONENT"s are specified!)
//		defines relative strength of the fast scintillation process as a fraction of total scintillation yield (default is 0)
//
//         "WLSABSLENGTH":		(spectrum (in dependence of the photon energy))		(obligatory property for WLS materials!)
//		defines the absorption length (absorption spectrum) of the material, used for the WLS process (default is infinity, i.e. no WLS process)
//         "WLSCOMPONENT":		(spectrum (in dependence of the photon energy))		(obligatory property for WLS materials!)
//		defines the emission spectrum of the material, used for the WLS process
//         "WLSTIMECONSTANT":		(constant value (energy independent))
//		defines the decay time (time between absorption and emission), used for the WLS process (default is 0)
//         "WLSMEANNUMBERPHOTONS":	(constant value (energy independent))
//		defines the mean number of photons, emitted for each photon that was absorbed by the WLS material
//		(if specified, the real number of emitted photons is Poisson distributed, else the real number of emitted photons is 1)

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	//NOTE: !!!!!! In order to understand simulation of optical surface properties (and esspecially the relations between the different options and parameters),
//             you should definitely start with reading the diagram to be found at:
//             http://hypernews.slac.stanford.edu/HyperNews/geant4/get/AUX/2012/05/23/23.20-70533-ces_in_geant4_revised.png
//             (in http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/445.html)
//             The following comments are only a excerpt and CANNOT compete with the precision of this diagram !!!!!!
//
//      possible reflection models:
//         glisur:  original and obsolete GEANT3 model
//         unified: unified model (provides a range of different reflection mechanisms, cf. "possible Properties" below)
//         LUT:     using the Look-Up-Table for the surface simulation
//
//      possible surface types:
//         dielectric_dielectric: if both materials are dielectric, i.e. the reflection and refraction probabilities are defined by the refractive indices via the Fresnel equations
//                                (the reflection probability CANNOT be defined by the material's reflectivity, cf. "REFLECTIVITY" below)
//         dielectric_metal:      if one material is dielectric and the other is metal-like (i.e. the photons can only be reflected or absorbed but not refracted)
//                                this surface type has to be used if the reflection probability should be defined by the material's reflectivity
//                                (this does not work for dielectric_dielectric, cf. "REFLECTIVITY" below)
//         dielectric_LUT:        if Look-Up-Tables should be used
//         firsov:                for Firsov Process (O.B. Firsov, �Reflection of fast ions from a dense medium at glancing angles�, Sov. Phys.-Docklady, vol. 11, No. 8, pp.732-733, 1967)
//         x_ray:                 for x-ray mirror process
//
//      possible surface finishs: (cf. http://hypernews.slac.stanford.edu/HyperNews/geant4/get/AUX/2012/05/23/23.20-70533-ces_in_geant4_revised.png !!!)
//         for dielectric_metal:
//            polished:             perfectly polished surface
//                                  -> specular spike reflection
//                                     (requiring the reflectivity of the metal)
//            ground:               rough surface
//                                  -> specular spike, specular lobe reflection, (diffuse) Lambertian reflection or back scattering
//                                     (requiring the reflectivity of the metal and the angle alpha between a micro-facet normal and the average surface normal)
//
//         for dielectric_dielectric:
//            polished:             perfectly polished surface
//                                  -> specular spike reflection
//                                     (requiring the refraction indices of both materials)
//            ground:               rough surface
//                                  -> specular spike, specular lobe reflection, (diffuse) Lambertian reflection or back scattering
//                                     (requiring the refraction indices of both materials and the angle alpha between a micro-facet normal and the average surface normal)
//            polishedfrontpainted: the volume has a perfectly polished surface and an absorbing paint without air gap
//                                  -> specular spike reflection
//                                     (requiring the reflectivity of the paint)
//            groundfrontpainted:   the volume has a rough surface and an absorbing paint without air gap
//                                  -> (diffuse) Lambertian reflection
//                                     (requiring the reflectivity of the paint)
//            polishedbackpainted:  the volume has a rough surface and a perfectly polished, absorbing paint with air gap
//                                  -> volume-air surface: specular spike, specular lobe reflection, (diffuse) Lambertian reflection or back scattering
//                                     (requiring the refraction indices of the volume material and the surface as well as the angle alpha between a micro-facet normal and the average surface normal)
//                                  -> air-paint surface: specular spike reflection
//                                     (requiring the reflectivity of the paint)
//            groundbackpainted:    the volume has a rough surface and a rough, absorbing paint with air gap
//                                  -> volume-air surface: specular spike, specular lobe reflection, (diffuse) Lambertian reflection or back scattering
//                                     (requiring the refraction indices of the volume material and the surface as well as the angle alpha between a micro-facet normal and the average surface normal)
//                                  -> air-paint surface: (diffuse) Lambertian reflection
//                                     (requiring the reflectivity of the paint)
//
//         for the Look-Up-Table model:
//            polishedlumirrorair:
//            polishedlumirrorglue:
//            polishedair:
//            polishedteflonair:
//            polishedtioair:
//            polishedtyvekair:
//            polishedvm2000air:
//            polishedvm2000glue:
//            etchedlumirrorair:
//            etchedlumirrorglue:
//            etchedair:
//            etchedteflonair:
//            etchedtioair:
//            etchedtyvekair:
//            etchedvm2000air:
//            etchedvm2000glue:
//            groundlumirrorair:
//            groundlumirrorglue:
//            groundair:
//            groundteflonair:
//            groundtioair:
//            groundtyvekair:
//            groundvm2000air:
//            groundvm2000glue:
//
//      sigma alpha: defines the roughness of the surface (sigma alpha specifies the standard deviation of the distribution of the micro-facet normals in [rad])

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//NOTE: possible Properties (for reflection model "unified"):
//         "RINDEX":				(spectrum (in dependence of the photon energy))			(obligatory property for surfaces with surface finish "...backpainted"!)
//		in case of surface finish "...backpainted", defines the refraction index of the thin material layer between the surface and the paint
//		(only used for boundary processes and only in case of surface finish "...backpainted")
//         "REFLECTIVITY":			(spectrum (in dependence of the photon energy))
//		does NOT IN ANY CASE define the reflectivity of the surface but defines (1 - absorption coefficient), which may make a difference:
//			for dielectric_metal:      it is really the reflectivity of the metal, i.e. the probability that the photon is reflected at all
//						   (every photon not absorbed by the metal will be reflected)
//			for dielectric_dielectric: it is NOT the reflectivity but defines the absorption coefficient (absorption coefficient = 1 - "reflectivity") of a dirty surface
//						   (every photon not absorbed by dirt will be reflected or refracted as normal for "dielectric_dielectric" surfaces and the reflectivity
//						   for this process can be specified via "TRANSMITTANCE" or is calculated via the Fresnel equations, cf. "TRANSMITTANCE")
//		(default is 1., i.e. nothing is absorbed)
//         "REALRINDEX":			(spectrum (in dependence of the photon energy))
//		defines the real part of the complex refraction index of the surface   //FIXME surface <-> material ???
//		(in case that "REFLECTIVITY" is not specified and if "REALRINDEX" and "IMAGINARYRINDEX" are both specified,
//		the refectivity is calculated from the complex refraction index via the Fresnel equations   //FIXME shouldn't it be two refraction indices???
//		-> therefore, the complex refraction index should not be specified for surface type "dielectric_dielectric", cf. "REFLECTIVITY"
//		   (if one wants to define the reflectivity of a "dielectric_dielectric" surface from the complex refraction index,
//		   one has to calculate the transmittance via the Fresnel equations and specify it with "TRANSMITTANCE"))
//         "IMAGINARYRINDEX":			(spectrum (in dependence of the photon energy))
//		defines the imaginary part of the complex refraction index of the surface   //FIXME surface <-> material ???
//		(in case that "REFLECTIVITY" is not specified and if "REALRINDEX" and "IMAGINARYRINDEX" are both specified,
//		the refectivity is calculated from the complex refraction index via the Fresnel equations   //FIXME shouldn't it be two refraction indices???
//		-> therefore, the complex refraction index should not be specified for surface type "dielectric_dielectric", cf. "REFLECTIVITY"
//		   (if one wants to define the reflectivity of a "dielectric_dielectric" surface from the complex refraction index,
//		   one has to calculate the transmittance via the Fresnel equations and specify it with "TRANSMITTANCE"))
//         "EFFICIENCY":			(spectrum (in dependence of the photon energy))
//		defines the detection efficiency of absorbed photons (default is 0)
//         "TRANSMITTANCE":			(spectrum (in dependence of the photon energy))
//		in case of "dielectric_dielectric" surfaces, defines the fraction of photons (reaching the surface despite of dirt (cf. "REFLECTIVITY"))
//		which are refracted by the surface instead of being reflected
//		(if "TRANSMITTANCE" is not specified, the transmittance is calculated from the (real) refraction indices of the two materials forming the surface via the Fresnel equations)
//		(only used for boundary processes and only in case of surface type "dielectric_dielectric")
//         "SPECULARSPIKECONSTANT":		(spectrum (in dependence of the photon energy))
//		defines the probability for reflection at the average surface
//		(only used in case of surface finish "ground" or "...backpainted", default is 0)
//         "SPECULARLOBECONSTANT":		(spectrum (in dependence of the photon energy))
//		defines the probability for reflection at a micro facet surface, i.e. the direction is smeared around the direction of the specular spike reflection
//		(only used in case of surface finish "ground" or "...backpainted", default is 0)
//         "BACKSCATTERCONSTANT":		(spectrum (in dependence of the photon energy))
//		defines the probability of back scattering, caused by several reflections within a deep grove
//		(only used in case of surface finish "ground" or "...backpainted", default is 0)
//	    !!! the probability of internal (diffuse) Lambertian reflection can not be specified directly but is defined via
//		    100% = "SPECULARSPIKECONSTANT" + "SPECULARLOBECONSTANT" + "BACKSCATTERCONSTANT" + Lambertian (-> default is 1) !!!
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Define detector optical surface property
	G4OpticalSurface* detector2Wrap = new G4OpticalSurface("detector2Wrap");

	//no refraction, only reflection or absorption
	detector2Wrap ->SetType (dielectric_dielectric); // dielectric-dielectric for all the reflectors
	detector2Wrap ->SetFinish(polishedbackpainted); // only specular spike reflection takes place. The wrapping is a perfectly smooth mirror.
	detector2Wrap ->SetModel(unified);    // glisur was Geant3 model.
	detector2Wrap ->SetSigmaAlpha(0.0227);   // 1.3deg

	//define the wrap surface as a perfect reflection mirror
	const G4int num = 5;
	G4double pp[num] = { 2.761*eV,2.485*eV,2.258*eV,2.071*eV,1.911*eV};
	G4double reflectivity[num] = {1.00,1.00,1.00,1.00,1.00};   //0.99  QT throws exception when this is set to 1.0
	G4double efficiency[num] = {1.00, 1.00, 1.00, 1.00, 1.00};  //0.99 define the detection efficiency of absorbed photons, default is 0.
	G4double rindex[num] = {1.35,1.35,1.35,1.35,1.35};  // refraction index of Teflon
	//G4double specularlobe[num] = {0.3, 0.3, 0.3, 0.3, 0.3};  // temparaly data
	//G4double specularspike[num] = {0.2, 0.2, 0.2, 0.2, 0.2}; // temparaly data
	//G4double backscatter[num] = {0.1, 0.1, 0.1, 0.1, 0.1};   // temparaly data

	G4MaterialPropertiesTable* detectorWrapProperty = new G4MaterialPropertiesTable();
	detectorWrapProperty->AddProperty("REFLECTIVITY",pp,reflectivity,num);
	detectorWrapProperty->AddProperty("EFFICIENCY",pp,efficiency,num);	
	detectorWrapProperty->AddProperty("RINDEX",pp,rindex,num);

	//detectorWrapProperty -> AddProperty("SPECULARLOBECONSTANT",pp,specularlobe,num);
	//detectorWrapProperty -> AddProperty("SPECULARSPIKECONSTANT",pp,specularspike,num);
	//detectorWrapProperty -> AddProperty("BACKSCATTERCONSTANT",pp,backscatter,num) ;
	
	detector2Wrap->SetMaterialPropertiesTable(detectorWrapProperty);

	G4LogicalBorderSurface* surface1 = new G4LogicalBorderSurface("betadetectorWrap", physBeta_Detector, physWrap,detector2Wrap);

	G4LogicalBorderSurface* surface2 = new G4LogicalBorderSurface("detectorWrap", physDetector, physWrap,detector2Wrap);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// define optical boundary between beta_detector and main_detector

	G4OpticalSurface* beta_det2main_det = new G4OpticalSurface("beta_det2main_det");
	
	//primary refraction and reflection
	beta_det2main_det ->SetType(dielectric_dielectric); // dielectric-dielectric for all the reflectors
	beta_det2main_det ->SetFinish(polished); // Snell's law is applied basedon refractive index of the two media
	beta_det2main_det ->SetModel(unified);    // glisur was Geant3 model.
	//beta_det2main_det ->SetSigmaAlpha(0.0227);   // 1.3deg


	const G4int num2 = 5;
	G4double pp2[num2] = { 2.761*eV,2.485*eV,2.258*eV,2.071*eV,1.911*eV};
	G4double reflectivity2[num2] = {1.00,1.00,1.00,1.00,1.00};   //0.99  QT throws exception when this is set to 1.0
	G4double efficiency2[num2] = {1.00, 1.00, 1.00, 1.00, 1.00};  //0.99 define the detection efficiency of absorbed photons, default is 0.
	//G4double rindex2[num2] = {1.35,1.35,1.35,1.35,1.35};  // refraction index of Teflon


	G4MaterialPropertiesTable* beta_det2main_det_Property = new G4MaterialPropertiesTable();
	beta_det2main_det_Property->AddProperty("REFLECTIVITY",pp2,reflectivity2,num2);
	beta_det2main_det_Property->AddProperty("EFFICIENCY",pp2,efficiency2,num2);	
	//beta_det2main_det_Property->AddProperty("RINDEX",pp2,rindex2,num2);

	beta_det2main_det->SetMaterialPropertiesTable(beta_det2main_det_Property);

	G4LogicalBorderSurface* surface3 = new G4LogicalBorderSurface("beta_det2main_det", physBeta_Detector, physDetector,beta_det2main_det);
 
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Define the optical surface between detector and sensor
	// the surface is the optical grease from GE bayer Silicones Albany NY
	G4OpticalSurface* PhotonSensor_op_surf = new G4OpticalSurface("PhotonSensor_op_surf");

	//sigma_alpha = 0.9; //surface roughness
	PhotonSensor_op_surf ->SetType (dielectric_dielectric);
	PhotonSensor_op_surf ->SetFinish(polished);  // no refraction, only reflection or absorption
	PhotonSensor_op_surf ->SetModel(unified);  //was glisur
	//PhotonSensor_op_surf ->SetSigmaAlpha(0.209);  // 12deg

	const G4int num1 = 5;
	G4double pp1[num1] = { 2.761*eV,2.485*eV,2.258*eV,2.071*eV,1.911*eV};
	G4double reflectivity1[num1] = {1.0,1.0,1.0,1.0,1.0};   // 0% reflectivity, all absorbed
	G4double efficiency1[num1] = {1.0, 1.0, 1.0, 1.0, 1.0};
	//G4double rindex1[num1] = {1.4042,1.4042,1.4042,1.4042,1.4042};
	//G4double specularlobe1[num1] = {0.3, 0.3, 0.3, 0.3, 0.3};  // temparaly data
	//G4double specularspike1[num1] = {0.5, 0.5, 0.5, 0.5, 0.5}; // temparaly data
	//G4double backscatter1[num1] = {0.1, 0.1, 0.1, 0.1, 0.1};   // temparaly data

	G4MaterialPropertiesTable* PhotonSensor_op_surf_Property = new G4MaterialPropertiesTable();
	PhotonSensor_op_surf_Property->AddProperty("REFLECTIVITY",pp1,reflectivity1,num1);
	PhotonSensor_op_surf_Property->AddProperty("EFFICIENCY",pp1,efficiency1,num1);	
	//PhotonSensor_op_surf_Property->AddProperty("RINDEX",pp1,rindex1,num1);

	//PhotonSensor_op_surf_Property->AddProperty("SPECULARLOBECONSTANT", pp1,specularlobe1, num1);   //0
	//PhotonSensor_op_surf_Property->AddProperty("SPECULARSPIKECONSTANT", pp1,specularspike1, num1);    //1
	//PhotonSensor_op_surf_Property->AddProperty("BACKSCATTERCONSTANT", pp1,backscatter1, num1);    //0


	PhotonSensor_op_surf->SetMaterialPropertiesTable(PhotonSensor_op_surf_Property);
	
	G4LogicalBorderSurface* surface4 = new G4LogicalBorderSurface("PhotonSensor_op_surf", physDetector, physSensor,PhotonSensor_op_surf);


	
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

#ifdef ENABLE_BETASD_GAMMA_SD
	#ifndef REMOVE_BETA_DET
	MDM_BetaSD* betaSD = new MDM_BetaSD("Beta","betaHitsCollection",1);
	// register the sensitive detector with manager
	G4SDManager::GetSDMpointer()->AddNewDetector(betaSD);
	SetSensitiveDetector("Beta_Detector",betaSD);
	#endif
#endif

}

