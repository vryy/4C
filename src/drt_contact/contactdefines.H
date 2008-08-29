/*!----------------------------------------------------------------------
\file contactdefines.H

\brief A set of preprocessor defines for contact methods

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef CONTACTDEFINES_H
#define CONTACTDEFINES_H

/************************************************************************/
/* Contact algorithm options                                            */
/************************************************************************/

//#define CONTACTWNORMAL         /* flag for weighted nodal normal calc. */
//#define CONTACTBOUNDMOD        /* flag for modification of boundary nodes */
//#define CONTACTONEMORTARLOOP   /* flag for computation of D,M matrices */

/************************************************************************/
/* Debugging options                                                    */
/************************************************************************/

//#define CONTACTDIMOUTPUT       /* flag for output of global matrix dim. */
//#define CONTACTCHECKHUEEBER    /* flag for S.Hüeber's quad. conv. check */
//#define CONTACTFDNORMAL        /* flag for FD check of normal derivative */
//#define CONTACTFDMORTARD       /* flag for FD check of mortar matrix D */
//#define CONTACTFDMORTARM       /* flag for FD check of mortar matrix M */
//#define CONTACTFDGAP           /* flag for FD check of weighted gap g */
//#define CONTACTFDTANGLM        /* flag for FD check of tangential LM */
//#define CONTACTGMSH1           /* flag for gmsh output in each time step */
//#define CONTACTGMSH2           /* flag for gmsh output in each iteration */

/************************************************************************/
/* Contact algorithm parameters                                         */
/************************************************************************/

#define CONTACTPROJTOL 0.05      /* projection tol. for overlap decision */
#define CONTACTPROJLIM 1.0e-8    /* no problems with exact node projections */
#define CONTACTMAXITER 10        /* max. no. iterations for local Newton */
#define CONTACTCONVTOL 1.0e-12   /* convergence tolerance for local Newton */
#define CONTACTCRITDIST 0.5      /* critical distance for potential contact */
#define CONTACTDEPTHOCTREE 5     /* max. depth for contact search octree */

namespace CONTACT {
  // static solution routine
  void contact_stru_static_drt();
}

#endif  // #ifndef CONTACTDEFINES_H
#endif  // #ifdef CCADISCRET
