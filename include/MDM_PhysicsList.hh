//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// * v0.2   01/10/2015  update to the same as example OpNovice        *
// ********************************************************************
//
// $Id: MDM_PhysicsList.hh 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file MDM_PhysicsList.hh
/// \brief Definition of the MDM_PhysicsList class

#ifndef MDM_PhysicsList_h
#define MDM_PhysicsList_h 1


#include "globals.hh"
#include "G4VUserPhysicsList.hh"

//class OpNovicePhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MDM_PhysicsList : public G4VUserPhysicsList
{
  public:
	MDM_PhysicsList();
    virtual ~MDM_PhysicsList();

  public:
    virtual void ConstructParticle();
    virtual void ConstructProcess();

    virtual void SetCuts();

    //these methods Construct physics processes and register them
    void ConstructDecay();
    void ConstructEM();
    void ConstructOp();

    //for the Messenger
    void SetVerbose(G4int);
    void SetNbOfPhotonsCerenkov(G4int);

  private:
    G4int                fVerboseLebel;
    //OpNovicePhysicsListMessenger* fMessenger;
    G4int fMaxNumPhotonStep;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

