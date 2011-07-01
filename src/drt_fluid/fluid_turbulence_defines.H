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
#ifdef CCADISCRET
#ifndef FLUID_TURBULENCE_DEFINES_H
#define FLUID_TURBULENCE_DEFINES_H


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

//> subgrid-stress term is considered in residual of momentum equation for stabilization terms
//> only original version
//#define TAU_SUBGRID_IN_RES_MOM

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

//> partially integrated form
//#define PARTINT


//----------------------------------//
// multifractal subgrid-scales      //
//----------------------------------//

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
#define SGS_TERM_CROSS
//> Reynolds-stress term
//> currently conservative form (may be switched off)
      /*
               /                       \
              |                         |
              | ( du o nabla du) ,  v   |
              |                         |
               \                       /
      */
#define SGS_TERM_REYNOLDS

//> include lhs contribution (as far as possible)
//#define LHS

//> partially integrated form (only rhs side)
//#define PARTINT

//> filter width ratio
    /*
                    Test Filter
           alpha = -------------
                    Grid Filter
     */
//> algebraic_multigrid_operator: alpha = 3
//> box_filter: alpha = 2
#define ALPHA 3.0
//> proportionality constant
#define C 0.125

//> parameter depending on the ratio of resolved scales and viscous scales
//> same values for x, y and z if DIRNUM is turned of
#define NUMX 2
#define NUMY 2
#define NUMZ 2
//> turn direction-dependence on
//#define DIRNUM

//> calculate N depending on the element Reynolds-number for each element
#define RE_ELE
//> set reference length for Reynolds-number
//> stream length
//#define STREAM_LENGTH
//> volume equivalent diameter
//#define VOLUME_EQ_DIA
//> length of cube with same volume
#define CUBIC_LENGTH
//> set reference velocity
//> use resolved velocity at element center
//#define RESOLVED_VEL
//> use fine-scale velocity at element center
//#define FINESCALE_VEL
//> scale element Reynolds-number
//#define FINESCALE_VEL_SCALED
//> factor for scaling
//#define FAC 0.75
//> use strain rate
#define STRAINRATE
//> set N equal to zero if N < Tol
//#define CLIPPING
//#define TOL 1.0
//> set N equal to zero in all elements adjacent to the wall
//#define ZERO_N


//> additional eddy viscosity models

//> small-scale eddy viscosity (similar to avm3)
#define FINE_SCALE_EDDY_VISC
//> set Smagorinsky-constant
#define MYCS 0.06
//> set velocity for calculation of subgrid-viscosity
//> use resolved velocity
//#define ALL_SCALE
//> use fine-scale velocity
#define SMALL_SCALE
//> use subgrid-scale velocity of multifractal subgrid-scale model
//#define SUBGRID_SCALE

//> calculate subgrid-viscosity based on subgid-energy
//#define SUBGRID_ENR_EDDY_VISC

#endif // #ifndef FLUID_TURBULENCE_DEFINES_H
#endif // #ifdef CCADISCRET
