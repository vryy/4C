/*----------------------------------------------------------------------------*/
/*!
\file poromultiphase_scatra_artery_coupling_defines.H

\brief A set of defines for coupling between poromultiphase_scatra-framework
       and flow in artery networks including scalar transport (mainly for
       projection algorithms)

\level 3

\maintainer Johannes Kremheller
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
#define SMALLESTSEGMENT 1.0e-10 /* segments cannot be smaller than this */
// default: 1.0e-12
#define ETAABTOL 1.0e-10 /* used for check if segments are identical */
// default: 1.0e-10
#define VALIDPROJTOL                                                 \
  1.0e-9                    /* used for check if projection is valid \
//default: 1.0e-9                (lies inside parameter space of this element) */
#define COLINEARTOL 1.0e-12 /* used for check if elements are colinear */
// default: 1.0e-12
#define CONVTOLNEWTONPROJ 1.0e-14 /* tolerance for Newton projection algorithm */
// default: 1.0e-14
#define NUMPROJCHECKS 1001 /* number of checks performed in interval [-0.999, 0.999] */
// default: 51
#define THROW_ERR_AT_DEGENERATE_CASE false /* throw error if degenerate case is detected */
// default: false
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
#define MAXNUMSEGPERELE 7 /* maximum number of segments per artery element */
// default: 5
#endif /* SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_DEFINES_H_ */
