/*!---------------------------------------------------------------------------
\file elemag_ele_action.H

\brief Provides enum of actions for electromagnetics elements

<pre>
\level 2

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            089 - 289-15244
</pre>
*---------------------------------------------------------------------------*/

#ifndef ELEMAG_ELE_ACTION_H
#define ELEMAG_ELE_ACTION_H

namespace ELEMAG
{
  /*--------------------------------------------------------------------------
   | enum that provides all possible elemag actions
   *--------------------------------------------------------------------------*/
  enum Action
  {
    none,                                         // 0
    calc_systemmat_and_residual,                  // 1
    calc_abc,                                     // 2
    calc_elemag_error,                            // 3
    project_field,                                // 4
    update_secondary_solution,                    // 5
    interpolate_hdg_to_node,                      // 6
    update_secondary_solution_and_calc_residual,  // 7
    ele_init,                                     // 8
    fill_restart_vecs,                            // 9
    ele_init_from_restart,                        // 10
    project_dirich_field,                         // 11
    bd_integrate,                                 // 12
    get_gauss_points,                             // 13
    project_field_test,                           // 14
    project_field_test_trace,                     // 15
    compute_error                                 // 16
  };                                              // enum Action

}  // namespace ELEMAG


#endif /* ELEMAG_ELE_ACTION_H */
