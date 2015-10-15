/*
 * PersistencyHandler.hh
 *
 *  Created on: Feb 11, 2014
 *      Author: Tim Niggemann, III Phys. Inst. A, RWTH Aachen University
 */

#ifndef PERSISTENCYHANDLER_HH_
#define PERSISTENCYHANDLER_HH_

#include <boost/variant.hpp>

#include "PersistVisitor.hh"

/**
 * Variant for all supported, persistable types.
 */
typedef boost::variant<G4SipmHitsCollection*, G4SipmDigiCollection*, G4SipmVoltageTraceDigiCollection*,
		G4SipmUiMessenger*,  G4SipmModel*, G4SipmVoltageTraceModel*> Persistable;

class PersistencyHandler {
private:
	PersistVisitor* visitor;

public:
	PersistencyHandler(PersistVisitor* visitor);
	virtual ~PersistencyHandler();

	void open(std::string filename);
	void persist(Persistable p);
	void close();

	PersistVisitor* getVisitor() const;
};

#endif /* PERSISTENCYHANDLER_HH_ */
