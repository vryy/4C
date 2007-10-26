/*!
\file physics.H

\brief contains information about physical fields

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef PHYSICS_H
#define PHYSICS_H

#include <string>
#include <sstream>

namespace Physics
{

// unique name of a physical field
// the order they appear here determines the order in the unknown vectors!!!
enum PhysVar
{
    Dispx,              // displacements
    Dispy,
    Dispz,
    Velx,               // velocities
    Vely,
    Velz,
    Temp,               // temerature
    Pres,               // pressure
    ConcentrationH2O,   // concentration of solvent
    ConcentrationNO2,
    LMPLambdax,         // Lagrange multiplier unknowns
    LMPLambday,
    LMPLambdaz,
    LMPMux,
    LMPMuy,
    LMPMuz
};

std::string physVarToString(const PhysVar var);

};


#endif /*PHYSICS_H_*/
#endif  // #ifdef CCADISCRET
