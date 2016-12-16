/*----------------------------------------------------------------------*/
/*!
 \file porofluidmultiphase_ele_action.H

 \brief element actions for evaluation of porofluidmultiphase element

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_ACTION_H_
#define SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_ACTION_H_


namespace POROFLUIDMULTIPHASE
{

/*--------------------------------------------------------------------------*/
/*!
 * \brief enum that provides all possible porofluidmultiphase actions
 *///                                                            vuong 08/16
/*--------------------------------------------------------------------------*/
enum Action {
  // domain action
  set_general_parameter,     // set general parameters for element evaluation
  set_timestep_parameter,    // set time-integration parameters for element evaluation
  calc_mat_and_rhs,          // calc system matrix and residual,
  calc_pres_and_sat,         // calculate pressure and saturation
  calc_solidpressure,        // calculate solid pressure
  calc_error,                // calc_error
  recon_flux_at_nodes        // reconstruct flux at nodes
}; // enum Action

/*--------------------------------------------------------------------------*/
/*!
 * \brief enum that provides all possible porofluidmultiphase actions on a boundary
 *///                                                            vuong 08/16
/*--------------------------------------------------------------------------*/
enum BoundaryAction {
  bd_calc_Neumann,                  // evaluate neumann loads
}; // enum POROFLUIDMULTIPHASE::BoundaryAction
} // namespace POROFLUIDMULTIPHASE



#endif /* SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_ACTION_H_ */
