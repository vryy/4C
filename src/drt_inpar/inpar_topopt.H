/*!------------------------------------------------------------------------------------------------*
\file inpar_topopt.H

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef INPAR_TOPOPT_H_
#define INPAR_TOPOPT_H_



namespace INPAR
{
  /*!
  \brief input parameter section for topology optimization

  The enumerations of input parameters required within the topology optimization
  are listed here

  */
  namespace TOPOPT
  {

    /// type of convergence check
    enum ConvCheck
    {
      inc,
      res,
      inc_and_res
    };

    /// reentry point of reinitialization
    enum Restart_Point
    {
      fluid,
      adjoint,
      gradient,
      opti_step
    };

    /// initial field for scalar transport problem
    enum InitialField
    {
      initfield_zero_field,
      initfield_field_by_function
    };
  }
}
#endif /* INPAR_TOPOPT_H_ */
#endif  // #ifdef CCADISCRET
