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
#ifdef D_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef PHYSICS_H
#define PHYSICS_H


namespace Physics
{

// unique name of a physical field
enum PhysVar
{
    Dispx,
    Dispy,
    Dispz,
    Velx,
    Vely,
    Velz,
    Temp,
    Pres,
    ConcentrationH2O,
    ConcentrationNO2,
    LagrMultLambdax,
    LagrMultLambday,
    LagrMultLambdaz,
    LagrMultMux,
    LagrMultMuy,
    LagrMultMuz
};

};


#endif /*PHYSICS_H_*/
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_XFEM
