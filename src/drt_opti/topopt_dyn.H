/*---------------------------------------------------------------------*/
/*! \file

\brief control routine of topology optimization

\maintainer Martin Kronbichler

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef TOPOPT_DYN_H
#define TOPOPT_DYN_H


/*!

  \brief control routine for topology optimization on fluidic domains

  \detailed

  Basically two approaches exist to optimize the topology of fluidic domains:
  1) A more shape optimization-like approach where the fluidic boundary is moved due
  to an topology sensitivity analysis
  2) Density based approaches where a density function rho_opt(x) ranging in [0 1]
  is used to compute a pseudo-porosity alpha. An according porosity term is added
  to the NS-equations resulting in zero-velocities when rho_opt is 0 and
  accordingly alpha >> 1 (=solid) and standard NS-equations when rho_opt is 1 and
  alpha << 1 (=fluid)

  Going the second approach, the optimization is performed via a sensitivity
  analysis using the adjoint equations of the extended NS-equations


  \author winklmaier
  \date 12/11
*/

void fluid_topopt_dyn();


#endif /* #ifndef TOPOPT_DYN_H */
