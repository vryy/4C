/*---------------------------------------------------------------------*/
/*! \file

\brief Contact interface used for constitutive laws

\level 3


*/
/*---------------------------------------------------------------------*/
#include "4C_contact_constitutivelaw_interface.hpp"

#include "4C_contact_constitutivelaw_contactconstitutivelaw.hpp"
#include "4C_contact_constitutivelaw_contactconstitutivelaw_parameter.hpp"
#include "4C_contact_defines.hpp"
#include "4C_contact_element.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_inpar_mortar.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_CrsMatrix.h>
#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                                       |
 *----------------------------------------------------------------------*/
CONTACT::ConstitutivelawInterface::ConstitutivelawInterface(
    const Teuchos::RCP<Mortar::InterfaceDataContainer>& interfaceData, const int id,
    const Epetra_Comm& comm, const int dim, const Teuchos::ParameterList& icontact,
    bool selfcontact, const int contactconstitutivelawid)
    : Interface(interfaceData, id, comm, dim, icontact, selfcontact)
{
  Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> coconstlaw =
      CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::Factory(contactconstitutivelawid);
  coconstlaw_ = coconstlaw;
  return;
}
/*----------------------------------------------------------------------*
 |  Evaluate regularized normal forces (nodes)                          |
 *----------------------------------------------------------------------*/
void CONTACT::ConstitutivelawInterface::assemble_reg_normal_forces(
    bool& localisincontact, bool& localactivesetchange)
{
  // loop over all slave row nodes on the current interface
  for (int i = 0; i < SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    Core::Nodes::Node* node = Discret().gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    int dim = cnode->NumDof();
    double gap = cnode->Data().Getg();

    double kappa = cnode->Data().Kappa();

    double lmuzawan = 0.0;
    for (int k = 0; k < dim; ++k) lmuzawan += cnode->MoData().lmuzawa()[k] * cnode->MoData().n()[k];

#ifdef CONTACTFDPENALTYKC1
    // set lagrangian multipliers explicitly to constant
    // and corresponding derivatives to zero

    for (int j = 0; j < dim; ++j) cnode->MoData().lm()[j] = i * j;

    cnode->Data().GetDerivZ().clear();

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
      // Evaluate pressure
      double pressure = coconstlaw_->evaluate(kappa * gap, cnode);
      // Evaluate pressure derivative
      double pressurederiv = coconstlaw_->EvaluateDeriv(kappa * gap, cnode);

      localisincontact = true;

      double* normal = cnode->MoData().n();

      // compute lagrange multipliers and store into node
      for (int j = 0; j < dim; ++j) cnode->MoData().lm()[j] = (lmuzawan - pressure) * normal[j];

      // compute derivatives of lagrange multipliers and store into node
      // contribution of derivative of weighted gap
      std::map<int, double>& derivg = cnode->Data().GetDerivG();
      std::map<int, double>::iterator gcurr;
      // printf("lm=%f\n", -coconstlaw_->evaluate(kappa * gap));

      // contribution of derivative of normal
      std::vector<Core::Gen::Pairedvector<int, double>>& derivn = cnode->Data().GetDerivN();
      Core::Gen::Pairedvector<int, double>::iterator ncurr;

      for (int j = 0; j < dim; ++j)
      {
        for (gcurr = derivg.begin(); gcurr != derivg.end(); ++gcurr)
          cnode->AddDerivZValue(
              j, gcurr->first, -kappa * pressurederiv * (gcurr->second) * normal[j]);
        for (ncurr = (derivn[j]).begin(); ncurr != (derivn[j]).end(); ++ncurr)
          cnode->AddDerivZValue(j, ncurr->first, -pressure * ncurr->second);
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
      cnode->Data().GetDerivZ().clear();

    }  // Macauley-Bracket
  }    // loop over slave nodes

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate regularized tangential forces                              |
 *----------------------------------------------------------------------*/
void CONTACT::ConstitutivelawInterface::assemble_reg_tangent_forces_penalty()
{
  FOUR_C_THROW("Frictional contact not yet implemented for rough surfaces\n");
}

FOUR_C_NAMESPACE_CLOSE
