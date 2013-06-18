/*!-----------------------------------------------------------------------------------------------*
\file fluid_turbulence_defines.H

  \brief preprocessor define flags for turbulence model development

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifndef FLUID_TURBULENCE_DEFINES_H
#define FLUID_TURBULENCE_DEFINES_H

namespace FLD {

//#define TESTING

//----------------------------------//
// scale similarity model           //
//----------------------------------//

//> subgrid-scale model as proposed by Liu, Meneveau and Katz (1994)
          /*
                        /                                \
                       |             ^     ^   ^          |
                       | nabla o ( (u*u) - u * u ) ,  v   |
                       |                                  |
                        \                                /
          */
//> only right hand side
//#define SGS_TERM_PREFILTERING


//> modified version of the subgrid-scale model proposed by Liu, Meneveau and Katz (1994)
          /*
                        /                                \
                       |                   ^   ^          |
                       | nabla o ( (u*u) - u * u ) ,  v   |
                       |                                  |
                        \                                /
          */
//> only right hand side
//#define SGS_TERM

//> reformulated version of the modified scale similarity model
//> allows consideration of only one of the following terms
//> only partially integrated
//> only on right hand side
              /*
                            /                               \
                           |        ^   ^                    |
                           | ( du * u - u * du ) ,  eps(v)   |
                           |                                 |
                            \                               /
              */
//#define SGS_TERM_CROSS
              /*
                            /                      \
                           |                        |
                           | ( du * du ) , eps(v)   |
                           |                        |
                            \                      /
              */
//#define SGS_TERM_REYNOLDS

//----------------------------------//
// multifractal subgrid-scales      //
//----------------------------------//

// !! CAUTION !!: these options are not available for coupling with scalar transport

//> multifractal subgrid-scale model as proposed by Burton and Dahm (2005)
//> cross-stress term
//> currently conservative form (may be switched off)
      /*
               /                                      \
              |                                        |
              | ( du o nabla u - u o nabla du ) ,  v   |
              |                                        |
               \                                      /
      */
//#define SGS_TERM_CROSS
//> Reynolds-stress term
//> currently conservative form (may be switched off)
      /*
               /                       \
              |                         |
              | ( du o nabla du) ,  v   |
              |                         |
               \                       /
      */
//#define SGS_TERM_REYNOLDS

//> turn direction-dependence on: can be used for different Ns or weights for N
//> caution: smoothing not included
//#define DIR_N
//> parameter depending on the ratio of resolved scales and viscous scales
//> same values for x, y and z if DIRNUM is turned of
//#define NUMX 2
//#define NUMY 2
//#define NUMZ 2
//> weights for N depending on direction
//#define WEIGHT_NX 1
//#define WEIGHT_NY 1
//#define WEIGHT_NZ 0.5

//> calculate N depending on the element Reynolds-number for each element
//> length of element edge (minimum or maximum)
//#define HMIN
//#define HMAX

//> use constant B instead of calculation
//> Warning: overwrites calculation of B
//#define CONST_B
//> high Reynolds number limit: B = 0.42
//#define B_CONST 0.21


}
#endif // #ifndef FLUID_TURBULENCE_DEFINES_H
