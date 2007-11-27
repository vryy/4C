/*!
\file physics.cpp

\brief 

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET

#include "physics.H"
#include "../drt_lib/drt_dserror.H"

/// return a string representation for each enum value in XFEM::Physics::Field
std::string XFEM::Physics::physVarToString(const XFEM::Physics::Field var)
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
        case LMPLambdax:       text = "LMPLambdax"; break;
        case LMPLambday:       text = "LMPLambday"; break;
        case LMPLambdaz:       text = "LMPLambdaz"; break;
        default: dserror("no string defined for Field");
    };
    return text;
};

#endif  // #ifdef CCADISCRET
