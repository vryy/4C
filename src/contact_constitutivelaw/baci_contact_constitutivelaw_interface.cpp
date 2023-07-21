/*---------------------------------------------------------------------*/
/*! \file

\brief Contact interface used for constitutive laws

\level 3


*/
/*---------------------------------------------------------------------*/
#include "baci_contact_constitutivelaw_interface.H"
#include "baci_contact_constitutivelaw_contactconstitutivelaw.H"
#include "baci_contact_constitutivelaw_contactconstitutivelaw_parameter.H"

#include <Epetra_CrsMatrix.h>
#include <Epetra_FEVector.h>

#include "baci_contact_node.H"
#include "baci_contact_element.H"
#include "baci_contact_defines.H"

#include "baci_inpar_mortar.H"
#include "baci_inpar_contact.H"

#include "baci_lib_discret.H"
#include "baci_utils_exceptions.H"
#include "baci_lib_node.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                                       |
 *----------------------------------------------------------------------*/
CONTACT::ConstitutivelawInterface::ConstitutivelawInterface(
    const Teuchos::RCP<MORTAR::InterfaceDataContainer>& interfaceData, const int id,
    const Epetra_Comm& comm, const int dim, const Teuchos::ParameterList& icontact,
    bool selfcontact, const int contactconstitutivelawid)
    : CoInterface(interfaceData, id, comm, dim, icontact, selfcontact)
{
  Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> coconstlaw =
      CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::Factory(contactconstitutivelawid);
  coconstlaw_ = coconstlaw;
  return;
}
/*----------------------------------------------------------------------*
 |  Evaluate regularized normal forces (nodes)                          |
 *----------------------------------------------------------------------*/
void CONTACT::ConstitutivelawInterface::AssembleRegNormalForces(
    bool& localisincontact, bool& localactivesetchange)
{
  // loop over all slave row nodes on the current interface
  for (int i = 0; i < SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
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
    if ((cnode->Active() == false) && (-coconstlaw_->Parameter()->GetOffset() - kappa * gap >= 0))
    {
      cnode->Active() = true;
      localactivesetchange = true;
    }

    else if ((cnode->Active() == true) &&
             (-coconstlaw_->Parameter()->GetOffset() - kappa * gap < 0))
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
      std::vector<CORE::GEN::pairedvector<int, double>>& derivn = cnode->CoData().GetDerivN();
      CORE::GEN::pairedvector<int, double>::iterator ncurr;

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
