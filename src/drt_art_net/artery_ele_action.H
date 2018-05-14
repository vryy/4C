/*----------------------------------------------------------------------*/
/*!
 \file artery_ele_action.H

 \brief element actions for evaluation of artery element

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_ART_NET_ARTERY_ELE_ACTION_H_
#define SRC_DRT_ART_NET_ARTERY_ELE_ACTION_H_


namespace ARTERY
{

/*--------------------------------------------------------------------------*/
/*!
 * \brief enum that provides all possible artery actions
 *///                                                   kremheller 03/18
/*--------------------------------------------------------------------------*/
enum Action {
  none,
  calc_sys_matrix_rhs,
  get_initial_artery_state,
  solve_riemann_problem,
  set_term_bc,
  set_scatra_term_bc,
  set_scatra_bc,
  calc_postpro_vals,
  calc_scatra_sys_matrix_rhs,
  calc_scatra_from_scatra_fb,
  evaluate_wf_wb,
  evaluate_scatra_analytically
};

}

#endif /* SRC_DRT_ART_NET_ARTERY_ELE_ACTION_H_ */
