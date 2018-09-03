/*!---------------------------------------------------------------------------
\file acou_ele_action.H

\brief Provides enum of actions for acoustic elements

<pre>
\level 2

\maintainer Luca Berardocco
            berardoccoo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
*---------------------------------------------------------------------------*/

#ifndef ACOU_ELE_ACTION_H
#define ACOU_ELE_ACTION_H

namespace ACOU
{
  /*--------------------------------------------------------------------------
   | enum that provides all possible acou actions
   *--------------------------------------------------------------------------*/
  enum Action
  {
    none,                                         // 0
    calc_systemmat_and_residual,                  // 1
    calc_abc,                                     // 2
    calc_pressuremon,                             // 3
    calc_acou_error,                              // 4
    project_field,                                // 5
    project_optical_field,                        // 6
    update_secondary_solution,                    // 7
    interpolate_hdg_to_node,                      // 8
    update_secondary_solution_and_calc_residual,  // 9
    interpolate_psi_to_node,                      // 10
    ele_init,                                     // 11
    fill_restart_vecs,                            // 12
    ele_init_from_restart,                        // 13
    project_dirich_field,                         // 14
    bd_integrate,                                 // 15
    calc_average_pressure,                        // 16
    calc_ader_sysmat_and_residual,                // 17
    update_ader_solution,                         // 18
    prepare_ader_postprocessing,                  // 19
    ader_postprocessing,                          // 20
    ader_postpro_gradanddiv,                      // 21
    get_gauss_points                              // 22
  };                                              // enum Action

}  // namespace ACOU


#endif /* ACOU_ELE_ACTION_H */
