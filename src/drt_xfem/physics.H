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

namespace XFEM
{
namespace PHYSICS
{

///
/// unique name of a physical field / primary unknown,
///
/// Note: the order they appear here determines the order in the unknown vectors!!!
enum Field
{
     Dispx              /// x displacement
    ,Dispy              /// y displacement
    ,Dispz              /// z displacement
    ,Velx               /// x velocity
    ,Vely               /// y velocity
    ,Velz               /// z velocity
    ,Temp               /// temerature
    ,Pres               /// pressure
    ,ConcentrationH2O   /// concentration of H2O
    ,ConcentrationNO2   /// concentration of NO2
    ,LMPLambdax         /// Lagrange multiplier unknown in x-direction
    ,LMPLambday         /// Lagrange multiplier unknown in y-direction
    ,LMPLambdaz         /// Lagrange multiplier unknown in z-direction
    ,LMPMux             /// another Lagrange multiplier unknown in x-direction
    ,LMPMuy             ///< another Lagrange multiplier unknown in y-direction
    ,LMPMuz             ///< another Lagrange multiplier unknown in z-direction
    ,Tauxx              ///< fluid deviatoric stress
    ,Tauxy              ///< fluid deviatoric stress
    ,Tauxz              ///< fluid deviatoric stress
    ,Tauyx              ///< fluid deviatoric stress
    ,Tauyy              ///< fluid deviatoric stress
    ,Tauyz              ///< fluid deviatoric stress
    ,Tauzx              ///< fluid deviatoric stress
    ,Tauzy              ///< fluid deviatoric stress
    ,Tauzz              ///< fluid deviatoric stress
};

std::string physVarToString(const Field var);

};
};


#endif  // PHYSICS_H
#endif  // #ifdef CCADISCRET
