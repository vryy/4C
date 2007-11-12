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


std::string XFEM::Physics::physVarToString(const Physics::PhysVar var)
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
        default: dserror("no string defined for PhysVar");
    };
    return text;
};

#endif  // #ifdef CCADISCRET
