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

namespace XFEM
{
namespace PHYSICS
{

/*!
 * \brief unique name of a physical field / primary unknown
 * 
 * If one has, for instance pressure at nodes AND element level, one has to define to distinct names for that purpose
 * like Pres and PresEle
 *
 * Note: the order they appear here determines the order in the unknown vectors!!!
 */
enum Field
{
    undefinedField,     ///< default vaue if nothing particular is defined
    Dispx,              ///< x displacement
    Dispy,              ///< y displacement
    Dispz,              ///< z displacement
    Velx,               ///< x velocity
    Vely,               ///< y velocity
    Velz,               ///< z velocity
    Temp,               ///< temperature
    Pres,               ///< pressure
    PresEle,            ///< pressure on element level
    ConcentrationH2O,   ///< concentration of H2O
    ConcentrationNO2,   ///< concentration of NO2
    LMPLambdax,         ///< Lagrange multiplier unknown in x-direction
    LMPLambday,         ///< Lagrange multiplier unknown in y-direction
    LMPLambdaz,         ///< Lagrange multiplier unknown in z-direction
    LMPMux,             ///< another Lagrange multiplier unknown in x-direction
    LMPMuy,             ///< another Lagrange multiplier unknown in y-direction
    LMPMuz,             ///< another Lagrange multiplier unknown in z-direction
    Tauxx,              ///< fluid deviatoric stress
    Tauxy,              ///< fluid deviatoric stress
    Tauxz,              ///< fluid deviatoric stress
    Tauyx,              ///< fluid deviatoric stress
    Tauyy,              ///< fluid deviatoric stress
    Tauyz,              ///< fluid deviatoric stress
    Tauzx,              ///< fluid deviatoric stress
    Tauzy,              ///< fluid deviatoric stress
    Tauzz,              ///< fluid deviatoric stress
    Sigmaxx,            ///< stress
    Sigmaxy,            ///< stress
    Sigmaxz,            ///< stress
    Sigmayx,            ///< stress
    Sigmayy,            ///< stress
    Sigmayz,            ///< stress
    Sigmazx,            ///< stress
    Sigmazy,            ///< stress
    Sigmazz             ///< stress
};

std::string physVarToString(const Field var);

}
}


#endif  // PHYSICS_H
#endif  // #ifdef CCADISCRET
