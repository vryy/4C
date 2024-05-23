/*----------------------------------------------------------------------*/
/*! \file

\brief manages a mesh coupling object with knowledge of a level set field

\level 3

*/
/*----------------------------------------------------------------------*/

#include "4C_xfem_coupling_mesh_coupled_levelset.hpp"

#include "4C_xfem_coupling_mesh.hpp"

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
//! constructor
XFEM::MeshCouplingNavierSlipTwoPhase::MeshCouplingNavierSlipTwoPhase(
    Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
    const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                   ///< discretization is identified
    Teuchos::RCP<DRT::Discretization>&
        cond_dis,           ///< discretization from which cutter discretization can be derived
    const int coupling_id,  ///< id of composite of coupling conditions
    const double time,      ///< time
    const int step,         ///< time step
    bool marked_geometry    ///< is this a marked geometry mesh boundary
    )
    : MeshCouplingNavierSlip(bg_dis, cond_name, cond_dis, coupling_id, time, step, marked_geometry)
{
}


void XFEM::MeshCouplingNavierSlipTwoPhase::set_condition_specific_parameters()
{
  XFEM::MeshCouplingNavierSlip::set_condition_specific_parameters();

  std::vector<CORE::Conditions::Condition*> conditions_NS;
  cutter_dis_->GetCondition(cond_name_, conditions_NS);
}

FOUR_C_NAMESPACE_CLOSE
