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

    /// type of adjoint equations
    enum AdjointType
    {
      discrete_adjoint,
      cont_adjoint
    };

    /// initial field for optimization variable
    enum InitialDensityField
    {
      initdensfield_zero_field,
      initdensfield_field_by_function
    };

    /// type of initial fluid field
    enum InitialAdjointField
    {
      initadjointfield_zero_field,
      initadjointfield_field_by_function
    };

    /// test configurations for adjoint equations
    enum AdjointCase
    {
      adjointtest_no,
      adjointtest_primal,
      adjointtest_stat_const_vel_lin_pres,
      adjointtest_stat_lin_vel_quad_pres,
      adjointtest_stat_quad_vel_lin_pres,
      adjointtest_stat_all_terms_all_constants,
      adjointtest_instat_varying_theta,
      adjointtest_instat_all_terms_all_constants,
      adjointtest_instat_primal_and_dual
    };

    /// test configuration for optimizer
    enum OptiCase
    {
      optitest_no,
      optitest_snake_one_constr,
      optitest_snake_multiple_constr,
      optitest_workflow_without_fluiddata
    };
  }
}
#endif /* INPAR_TOPOPT_H_ */
