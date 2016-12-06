/*!
\file physics.H

\brief contains information about physical fields

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
*/


#ifndef PHYSICS_H
#define PHYSICS_H

#include <string>

namespace XFEM
{
//! contains physical variables and related routines
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
    undefinedField,     ///< default value if nothing particular is defined
    Dispx,              ///< x displacement
    Dispy,              ///< y displacement
    Dispz,              ///< z displacement
    Velx,               ///< x velocity
    Vely,               ///< y velocity
    Velz,               ///< z velocity
    Veln,               ///< normal velocity
    Temp,               ///< temperature
    Pres,               ///< pressure
    DiscPres,           ///< discontinuous pressure on element level
    ConcentrationH2O,   ///< concentration of H2O
    ConcentrationNO2,   ///< concentration of NO2
    HeatFlux_x,         ///< heat flux in x direction
    HeatFlux_y,         ///< heat flux in y direction
    HeatFlux_z,         ///< heat flux in z direction
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
    Epsilonxx,          ///< fluid strain
    Epsilonxy,          ///< fluid strain
    Epsilonxz,          ///< fluid strain
    Epsilonyx,          ///< fluid strain
    Epsilonyy,          ///< fluid strain
    Epsilonyz,          ///< fluid strain
    Epsilonzx,          ///< fluid strain
    Epsilonzy,          ///< fluid strain
    Epsilonzz,          ///< fluid strain
    Sigmaxx,            ///< stress
    Sigmaxy,            ///< stress
    Sigmaxz,            ///< stress
    Sigmayx,            ///< stress
    Sigmayy,            ///< stress
    Sigmayz,            ///< stress
    Sigmazx,            ///< stress
    Sigmazy,            ///< stress
    Sigmazz,            ///< stress
    Velxiface,          ///< x velocity interface
    Velyiface,          ///< y velocity interface
    Velziface,          ///< z velocity interface
    Dispxiface,         ///< x displacement interface
    Dispyiface,         ///< y displacement interface
    Dispziface          ///< z displacement interface
};

//! return a std::string representation for each enum value in XFEM::PHYSICS::Field
std::string physVarToString(const Field var);

}
}


#endif  // PHYSICS_H
