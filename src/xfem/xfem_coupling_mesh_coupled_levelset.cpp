/*----------------------------------------------------------------------*/
/*! \file

\brief manages a mesh coupling object with knowledge of a level set field

\level 3

*/
/*----------------------------------------------------------------------*/

#include "xfem_coupling_mesh.H"
#include "xfem_coupling_mesh_coupled_levelset.H"

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
    : MeshCouplingNavierSlip(bg_dis, cond_name, cond_dis, coupling_id, time, step, marked_geometry),
      twophase_coupling_object_(Teuchos::null)
{
}

void XFEM::MeshCouplingNavierSlipTwoPhase::SetTwoPhaseCouplingPointer(
    Teuchos::RCP<XFEM::LevelSetCouplingTwoPhase> twophase_coupling_object)
{
  twophase_coupling_object_ = twophase_coupling_object;
  if (twophase_coupling_object_ == Teuchos::null)
    dserror("twophase_coupling_object_, object is null-pointer.");
}

void XFEM::MeshCouplingNavierSlipTwoPhase::SetConditionSpecificParameters()
{
  XFEM::MeshCouplingNavierSlip::SetConditionSpecificParameters();

  std::vector<DRT::Condition*> conditions_NS;
  cutter_dis_->GetCondition(cond_name_, conditions_NS);

  // Create maps for easy extraction at gausspoint level
  for (std::vector<DRT::Condition*>::iterator i = conditions_NS.begin(); i != conditions_NS.end();
       ++i)
  {
    int cond_int = (*i)->Id();

    double sliplength_smear = (*i)->GetDouble("slipsmear");
    if (!slipsmear_map_.insert(std::make_pair(cond_int, sliplength_smear)).second)
      dserror("ID already existing! For slipsmear_map_.");


    double normalpenalty_scaling = (*i)->GetDouble("normalpen_scaling");
    if (!scalednormalpenalty_map_.insert(std::make_pair(cond_int, normalpenalty_scaling)).second)
      dserror("ID already existing! For scalednormalpenalty_map_.");
  }
}
