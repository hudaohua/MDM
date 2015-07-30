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
// ***********************************************************************************************
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
//         firsov:                for Firsov Process (O.B. Firsov, “Reflection of fast ions from a dense medium at glancing angles”, Sov. Phys.-Docklady, vol. 11, No. 8, pp.732-733, 1967)
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
	G4OpticalSurface* detectorWrap = new G4OpticalSurface("detectorWrap");
	//G4LogicalBorderSurface* surface1 = new G4LogicalBorderSurface("detectorWrap", physDetector, physWorld,detectorWrap);

	//no refraction, only reflection or absorption
	//G4double sigma_alpha = 0.1;
	detectorWrap ->SetType (dielectric_metal);
	detectorWrap ->SetFinish(polished); // polished
	detectorWrap ->SetModel(glisur);    // glisur was Geant3 model.
	//detectorWrap ->SetSigmaAlpha(sigma_alpha);

	const G4int num = 5;
	G4double pp[num] = { 2.761*eV,2.485*eV,2.258*eV,2.071*eV,1.911*eV};
	G4double reflectivity[num] = {0.75,0.75,0.75,0.75,0.75};   // QT throws exception when this is set to 1.0
	G4double efficiency[num] = {0.8, 0.8, 0.8, 0.8, 0.8};  //define the detection efficiency of absorbed photons, default is 0.
	G4double rindex[num] = {1.56,1.56,1.56,1.56,1.56};
	//G4double specularlobe[num] = {0.3, 0.3, 0.3, 0.3, 0.3};  // temparaly data
	//G4double specularspike[num] = {0.2, 0.2, 0.2, 0.2, 0.2}; // temparaly data
	//G4double backscatter[num] = {0.1, 0.1, 0.1, 0.1, 0.1};   // temparaly data

	G4MaterialPropertiesTable* detectorWrapProperty = new G4MaterialPropertiesTable();
	detectorWrapProperty->AddProperty("REFLECTIVITY",pp,reflectivity,num);
	detectorWrapProperty->AddProperty("EFFICIENCY",pp,efficiency,num);	
	detectorWrapProperty->AddProperty("RINDEX",pp,rindex,num);

	//detectorWrapProperty -> AddProperty("SPECULARLOBECONSTANT",pp,specularlobe,NUM);
	//detectorWrapProperty -> AddProperty("SPECULARSPIKECONSTANT",pp,specularspike,NUM);
	//detectorWrapProperty -> AddProperty("BACKSCATTERCONSTANT",pp,backscatter,NUM);
	
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

	//sigma_alpha = 0.0; //surface roughness
	PhotonSensor_op_surf ->SetType (dielectric_metal);
	PhotonSensor_op_surf ->SetFinish(polished);
	PhotonSensor_op_surf ->SetModel(glisur);  //was glisur
	//PhotonSensor_op_surf ->SetSigmaAlpha(sigma_alpha);

	const G4int num1 = 5;
	G4double pp1[num1] = { 2.761*eV,2.485*eV,2.258*eV,2.071*eV,1.911*eV};
	G4double reflectivity1[num1] = {0.9,0.9,0.9,0.9,0.9};   // 10% reflectivity
	G4double efficiency1[num1] = {0.8, 0.8, 0.8, 0.8, 0.8};
//	G4double rindex1[num1] = {1.56,1.56,1.56,1.56,1.56};

	G4MaterialPropertiesTable* PhotonSensor_op_surf_Property = new G4MaterialPropertiesTable();
	PhotonSensor_op_surf_Property->AddProperty("REFLECTIVITY",pp1,reflectivity1,num1);
	PhotonSensor_op_surf_Property->AddProperty("EFFICIENCY",pp1,efficiency1,num1);	
//	PhotonSensor_op_surf_Property->AddProperty("RINDEX",pp1,rindex1,num1);

//	PhotonSensor_op_surf_Property->AddProperty("SPECULARLOBECONSTANT", pp1,reflectivity1, num1);   //0
//	PhotonSensor_op_surf_Property->AddProperty("SPECULARSPIKECONSTANT", pp1,efficiency1, num1);    //1
//	PhotonSensor_op_surf_Property->AddProperty("BACKSCATTERCONSTANT", pp1,reflectivity1, num1);    //0


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