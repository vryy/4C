/*--------------------------------------------------------------------------*/
/*!
\file lubrication_ele_action.H

\brief provides enum of actions. Motivation was to have a quick enum
comparison instead of a cumbersome std::string compare.

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_ACTION_H_
#define SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_ACTION_H_


namespace LUBRICATION
{
  /*--------------------------------------------------------------------------*/
  /*!
   * \brief enum that provides all possible lubrication actions
   *///                                                            wirtz 11/15
  /*--------------------------------------------------------------------------*/
  enum Action
  {
    // domain action
    set_general_lubrication_parameter,  // set general parameters for element evaluation
    set_time_parameter,                 // set time-integration parameters for element evaluation
    calc_mat_and_rhs,                   // calc_condif_systemmat_and_residual,
    calc_mean_pressures,                // calc_mean_pressures,
    calc_error,                         // calc_error
    calc_lubrication_coupltang          // calculate off-diagonal tangent matrix term
  };                                    // enum Action

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief enum that provides all possible lubrication actions on a boundary
   *///                                                            wirtz 11/15
  /*--------------------------------------------------------------------------*/
  enum BoundaryAction
  {
    // new action                       //PRE_BACI_CLEANING NAME:
    bd_calc_weak_Dirichlet,  // WeakDirichlet,
    bd_calc_Neumann,         // n/a
  };                         // enum LUBRICATION::BoundaryAction
}  // namespace LUBRICATION



#endif /* SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_ACTION_H_ */
