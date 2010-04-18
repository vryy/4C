/*----------------------------------------------------------------------*/
/*!
\file inpar_tsi.H

\brief

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                              dano 11/09 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef INPAR_TSI_H
#define INPAR_TSI_H

/*----------------------------------------------------------------------*
 |                                                           dano 11/09 |
 *----------------------------------------------------------------------*/
namespace INPAR
{

/// Input parameters for TSI
namespace TSI
{

//  /// Type of partitioned coupling for TSI problems
//  enum PartitionedCouplingMethod
//  {
//    // 23.11.09 check which method to use!!!
//    DirichletNeumann,
//    RobinNeumann,
//    DirichletRobin,
//    RobinRobin
//  };

}  // namespace TSI

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_TSI_H
#endif  // #ifdef CCADISCRET
