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

