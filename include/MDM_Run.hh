//
// ********************************************************************
// * Author: Hubert Hu												  *
// * Email: d.hu@ucl.ac.uk											  *
// * Electronics engineer @ Mullard Space Science Lab				  *

// * version log history											  *
// * v0.1	08/05/2015	convert from old test program to MDM project  *
// ********************************************************************
//
// $Id: MDM_Run.hh 66536 2012-12-19 14:32:36Z ihrivnac $
//
/// \file MDM_Run.hh
/// \brief Definition of the MDM_Run class

#ifndef MDM_Run_h
#define MDM_Run_h 1

#include "G4Run.hh"
#include "globals.hh"

class G4Event;

/// Run class
///

class MDM_Run : public G4Run
{
  public:
    MDM_Run();
    virtual ~MDM_Run();

    // method from the base class
    virtual void Merge(const G4Run*);
    
    void AddEdep (G4double edep); 

    // get methods
    G4double GetEdep()  const { return fEdep; }
    G4double GetEdep2() const { return fEdep2; }

  private:
    G4double  fEdep;
    G4double  fEdep2;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

