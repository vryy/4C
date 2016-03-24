/*!---------------------------------------------------------------------------
\file acou_ele_action.H

\brief Provides enum of actions.

<pre>
\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>

*---------------------------------------------------------------------------*/

#ifndef ACOU_ELE_ACTION_H
#define ACOU_ELE_ACTION_H

namespace ACOU
{

/*--------------------------------------------------------------------------
 | enum that provides all possible acou actions
 *--------------------------------------------------------------------------*/
enum Action {
  none,                               // 0
  calc_systemmat_and_residual,        // 1
  calc_abc,                           // 2
  calc_pressuremon,                   // 3
  calc_pmon_nodevals,                 // 4
  calc_acou_error,                    // 5
  project_field,                      // 6
  project_optical_field,              // 7
  update_secondary_solution,          // 8
  interpolate_hdg_to_node,            // 9
  update_secondary_solution_and_calc_residual, // 10
  interpolate_psi_to_node,            // 11
  ele_init,                           // 12
  fill_restart_vecs,                  // 13
  ele_init_from_restart,              // 14
  project_dirich_field,               // 15
  ele_init_from_acouoptrestart,       // 16
  bd_integrate,                       // 17
  calc_average_pressure               // 18
}; // enum Action

} // namespace ACOU


#endif /* ACOU_ELE_ACTION_H */
