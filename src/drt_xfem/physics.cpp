/*!
\file physics.cpp

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


#include <iostream>
#include "physics.H"
#include "../drt_lib/drt_dserror.H"

// return a std::string representation for each enum value in XFEM::PHYSICS::Field
std::string XFEM::PHYSICS::physVarToString(const XFEM::PHYSICS::Field var)
{
  std::string text;
  switch (var){
  case Dispx:            text = "Dispx"; break;
  case Dispy:            text = "Dispy"; break;
  case Dispz:            text = "Dispz"; break;
  case Velx:             text = "Velx "; break;
  case Vely:             text = "Vely "; break;
  case Velz:             text = "Velz "; break;
  case Temp:             text = "Temp "; break;
  case Pres:             text = "Pres "; break;
  case DiscPres:         text = "DiscPres "; break;
  case HeatFlux_x:       text = "HeatFlux_x"; break;
  case HeatFlux_y:       text = "HeatFlux_y"; break;
  case HeatFlux_z:       text = "HeatFlux_z"; break;
  case LMPLambdax:       text = "LMPLambdax"; break;
  case LMPLambday:       text = "LMPLambday"; break;
  case LMPLambdaz:       text = "LMPLambdaz"; break;
  case Tauxx:            text = "Tauxx"; break;
  case Tauxy:            text = "Tauxy"; break;
  case Tauxz:            text = "Tauxz"; break;
  case Tauyx:            text = "Tauyx"; break;
  case Tauyy:            text = "Tauyy"; break;
  case Tauyz:            text = "Tauyz"; break;
  case Tauzx:            text = "Tauzx"; break;
  case Tauzy:            text = "Tauzy"; break;
  case Tauzz:            text = "Tauzz"; break;
  case Epsilonxx:        text = "Epsilonxx"; break;
  case Epsilonxy:        text = "Epsilonxy"; break;
  case Epsilonxz:        text = "Epsilonxz"; break;
  case Epsilonyx:        text = "Epsilonyx"; break;
  case Epsilonyy:        text = "Epsilonyy"; break;
  case Epsilonyz:        text = "Epsilonyz"; break;
  case Epsilonzx:        text = "Epsilonzx"; break;
  case Epsilonzy:        text = "Epsilonzy"; break;
  case Epsilonzz:        text = "Epsilonzz"; break;
  case Sigmaxx:          text = "Sigmaxx"; break;
  case Sigmaxy:          text = "Sigmaxy"; break;
  case Sigmaxz:          text = "Sigmaxz"; break;
  case Sigmayx:          text = "Sigmayx"; break;
  case Sigmayy:          text = "Sigmayy"; break;
  case Sigmayz:          text = "Sigmayz"; break;
  case Sigmazx:          text = "Sigmazx"; break;
  case Sigmazy:          text = "Sigmazy"; break;
  case Sigmazz:          text = "Sigmazz"; break;
  case Velxiface:        text = "velxiface"; break;
  case Velyiface:        text = "velyiface"; break;
  case Velziface:        text = "velziface"; break;
  default:
    std::cout << var << std::endl;
    dserror("no std::string defined for Field");
  };
  return text;
}

