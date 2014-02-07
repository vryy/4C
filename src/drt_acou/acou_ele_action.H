/*!---------------------------------------------------------------------------
\file acou_ele_action.H

\brief Provides enum of actions.

<pre>
Maintainer: Svenja Schoeder
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
  none,
  calc_systemmat_and_residual,
  calc_abc,
  project_field,
  project_optical_field,
  update_secondary_solution,
  interpolate_hdg_to_node,
  calc_integr_objf,
  update_secondary_solution_and_calc_residual,
  update_secondary_solution_and_calc_mat_residual_noli
}; // enum Action

} // namespace ACOU


#endif /* ACOU_ELE_ACTION_H */
