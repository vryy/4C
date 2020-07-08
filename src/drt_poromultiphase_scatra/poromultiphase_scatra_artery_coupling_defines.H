/*----------------------------------------------------------------------------*/
/*! \file

\brief A set of defines for coupling between poromultiphase_scatra-framework
       and flow in artery networks including scalar transport (mainly for
       projection algorithms)

\level 3

*/
/*----------------------------------------------------------------------------*/
#ifndef SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_DEFINES_H_
#define SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_DEFINES_H_
/************************************************************************/
/* Projection                                                           */
/************************************************************************/
#define PROJMAXITER 10 /* max number of steps for local Newton */
// default: 10
#define PROJOUTPUT false /* output of projection algorithm */
// default: false
#define XIETATOL                                                                               \
  1.0e-9 /* used for check if two coordinates xi, eta in element parameter space are identical \
default: 1.0e-9  --> segments (in parameter space) cannot be smaller than this */
#define COLINEARTOL 1.0e-08 /* used for check if elements are colinear */
// default: 1.0e-08
#define CONVTOLNEWTONPROJ 1.0e-10 /* tolerance for Newton projection algorithm */
// default: 1.0e-10
/************************************************************************/
/* Random                                                               */
/************************************************************************/
#define KAPPAINVTOL 1.0e-10 /* used if kappa vector can be inverted */
// default: 1.0e-10
#define CONVTOLNEWTONMESHMOVEMENT                                                                                       \
  1.0e-09                      /* used for convergence check of Newton for                                              \
//default: 1.0e-9                             recomputation of eta and xi in current configuration \
                     */
#define MESHMOVEMENTMAXITER 10 /* max number of steps for local Newton (mesh movement) */
// default: 10
#endif /* SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_DEFINES_H_ */
