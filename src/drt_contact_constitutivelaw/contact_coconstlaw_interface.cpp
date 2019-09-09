/*---------------------------------------------------------------------*/
/*! \file

\brief Contact interface used for constitutive laws

\level 3

\maintainer Nora Hagmeyer

*/
/*---------------------------------------------------------------------*/
#include "contact_coconstlaw_interface.H"

#include <Epetra_CrsMatrix.h>
#include <Epetra_FEVector.h>
#include <Epetra_Time.h>
#include "../drt_contact/contact_node.H"
#include "../drt_contact/contact_element.H"
#include "../drt_contact/contact_integrator.H"
#include "../drt_contact/contact_interpolator.H"
#include "../drt_contact/contact_coupling2d.H"
#include "../drt_contact/contact_coupling3d.H"
#include "../drt_contact/contact_defines.H"
#include "../drt_contact/contact_line_coupling.H"
#include "../drt_contact/friction_node.H"
#include "../drt_mortar/mortar_binarytree.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_mortar/mortar_projector.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/drt_utils_parmetis.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_contact_constitutivelaw/coconstitutivelaw.H"
#include "../drt_contact_constitutivelaw/coconstlaw_parameter.H"

#include "../drt_mortar/mortar_coupling3d_classes.H"
#include "../drt_contact/contact_nitsche_utils.H"

#include <Teuchos_TimeMonitor.hpp>
#include "../drt_contact/selfcontact_binarytree_unbiased.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                                       |
 *----------------------------------------------------------------------*/
CONTACT::ConstitutivelawInterface::ConstitutivelawInterface(
    const Teuchos::RCP<MORTAR::IDataContainer>& idata_ptr, const int id, const Epetra_Comm& comm,
    const int dim, const Teuchos::ParameterList& icontact, bool selfcontact,
    INPAR::MORTAR::RedundantStorage redundant)
    : CoInterface(idata_ptr, id, comm, dim, icontact, selfcontact, redundant)
{
  Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> coconstlaw =
      CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::Factory(
          id + 1);  // todo This is temporary, until we can hand the interface its roughness law id
                    // via condition
  coconstlaw_ = coconstlaw;
  return;
}
/*----------------------------------------------------------------------*
 |  Evaluate regularized normal forces (nodes)                          |
 *----------------------------------------------------------------------*/
void CONTACT::ConstitutivelawInterface::AssembleRegNormalForces(
    bool& localisincontact, bool& localactivesetchange)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over all slave row nodes on the current interface
  for (int i = 0; i < SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    int dim = cnode->NumDof();
    double gap = cnode->CoData().Getg();

    double kappa = cnode->CoData().Kappa();

    double lmuzawan = 0.0;
    for (int k = 0; k < dim; ++k) lmuzawan += cnode->MoData().lmuzawa()[k] * cnode->MoData().n()[k];

#ifdef CONTACTFDPENALTYKC1
    // set lagrangian multipliers explicitly to constant
    // and corresponding derivatives to zero

    for (int j = 0; j < dim; ++j) cnode->MoData().lm()[j] = i * j;

    cnode->CoData().GetDerivZ().clear();

    continue;
#endif

    // Activate/Deactivate node and notice any change
    if ((cnode->Active() == false) && (coconstlaw_->Parameter()->GetOffset() - kappa * gap >= 0))
    {
      cnode->Active() = true;
      localactivesetchange = true;
    }

    else if ((cnode->Active() == true) && (coconstlaw_->Parameter()->GetOffset() - kappa * gap < 0))
    {
      cnode->Active() = false;
      localactivesetchange = true;

      // std::cout << "node #" << gid << " is now inactive, gap=" << gap << std::endl;
    }
    //********************************************************************

    // Compute derivZ-entries with the Macauley-Bracket
    // of course, this is only done for active constraints in order
    // for linearization and r.h.s to match!
    if (cnode->Active() == true)
    {
      localisincontact = true;

      double* normal = cnode->MoData().n();

      // compute lagrange multipliers and store into node
      for (int j = 0; j < dim; ++j)
        cnode->MoData().lm()[j] = (lmuzawan - coconstlaw_->Evaluate(kappa * gap)) * normal[j];

      // compute derivatives of lagrange multipliers and store into node
      // contribution of derivative of weighted gap
      std::map<int, double>& derivg = cnode->CoData().GetDerivG();
      std::map<int, double>::iterator gcurr;
      // printf("lm=%f\n", -coconstlaw_->Evaluate(kappa * gap));

      // contribution of derivative of normal
      std::vector<GEN::pairedvector<int, double>>& derivn = cnode->CoData().GetDerivN();
      GEN::pairedvector<int, double>::iterator ncurr;

      for (int j = 0; j < dim; ++j)
      {
        for (gcurr = derivg.begin(); gcurr != derivg.end(); ++gcurr)
          cnode->AddDerivZValue(j, gcurr->first,
              -kappa * coconstlaw_->EvaluateDeriv(kappa * gap) * (gcurr->second) * normal[j]);
        for (ncurr = (derivn[j]).begin(); ncurr != (derivn[j]).end(); ++ncurr)
          cnode->AddDerivZValue(
              j, ncurr->first, -coconstlaw_->Evaluate(kappa * gap) * ncurr->second);
        for (ncurr = (derivn[j]).begin(); ncurr != (derivn[j]).end(); ++ncurr)
          cnode->AddDerivZValue(j, ncurr->first, +lmuzawan * ncurr->second);
      }
    }

    // be sure to remove all LM-related stuff from inactive nodes
    else
    {
      // clear lagrange multipliers
      for (int j = 0; j < dim; ++j) cnode->MoData().lm()[j] = 0.0;

      // clear derivz
      cnode->CoData().GetDerivZ().clear();

    }  // Macauley-Bracket
  }    // loop over slave nodes

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate regularized tangential forces                              |
 *----------------------------------------------------------------------*/
void CONTACT::ConstitutivelawInterface::AssembleRegTangentForcesPenalty()
{
  dserror("Frictional contact not yet implemented for rough surfaces\n");
}
