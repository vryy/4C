/*---------------------------------------------------------------------*/
/*!
\file contact_coconstlaw_interface.cpp

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
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::ConstitutivelawInterface::ConstitutivelawInterface(
    const Teuchos::RCP<MORTAR::IDataContainer>& idata_ptr, const int id, const Epetra_Comm& comm,
    const int dim, const Teuchos::ParameterList& icontact, bool selfcontact,
    INPAR::MORTAR::RedundantStorage redundant)
    : CoInterface(idata_ptr, id, comm, dim, icontact, selfcontact, redundant)
{
  // For now, this should be the location of a CoLawDataContainer Factory (until input can be read
  // in) Only the container should then be passed to the factory.
  Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> coconstlawdata =
      CONTACT::CONSTITUTIVELAW::Container::ContainerFactory(
          INPAR::CONTACT::ConstitutiveLawType::colaw_penalty,
          IParams().get<double>("PENALTYPARAM"));
  Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> coconstlaw =
      CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::Factory(coconstlawdata);
  coconstlaw_ = coconstlaw;
  // printf("deriv=%f \n", coconstlaw_->Evaluate(0.0));
  return;
}
/*----------------------------------------------------------------------*
 |  Evaluate regularized normal forces (nodes)                popp 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::ConstitutivelawInterface::AssembleRegNormalForces(
    bool& localisincontact, bool& localactivesetchange)
{
  // if (coconstlaw_->GetCoConstLawType() != CONTACT::MAT::colaw_linear) printf("law is not
  // linear\n");
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

    // modified gap for zero initial gap
    // (if gap is below zero, it is explicitly set to zero)
    // double modgap = cnode->CoData().Getg();
    // if (abs(modgap) < 1.0e-10) modgap=0.0;

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

    //********************************************************************
    // Decision on active /  inactive nodes (regularization)
    //
    // CASE 1: Penalty approach
    // A node is activated if its weighted gap is negative or deactivated
    // if its gap is equal zero or positive.
    // -> the regularization reads: lambda_n = kappa * pp * < -gap >
    //
    // CASE 2: Uzawa augmented Lagrange approach
    // A node is activated if its Lagrange multiplier, stemming from the
    // last Uzawa Lagrange multiplier AND the current regularization is
    // negative or deactivated if its LM is equal zero or positive.
    // -> the regularization reads: lambda_n = < lmuzawa_n - kappa * pp * gap >
    //
    // As the Uzawa Lagrange multipliers are zero in the penalty approach,
    // the two cases can formally be treted identically, see below.
    // We do not need an explicit separation of cases!
    //
    //********************************************************************

    // Activate/Deactivate node and notice any change
    if ((cnode->Active() == false) && (coconstlaw_->Parameter()->GetOffset() - kappa * gap >= 0))
    {
      cnode->Active() = true;
      localactivesetchange = true;

      //      std::cout << "node #" << gid << " is now active (";
      //      for (int j = 0; j < dim; j++) std::cout << " " << cnode->Dofs()[j] << " ";
      //      std::cout << ") gap=" << gap * kappa << std::endl;
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
      //      std::cout << "GID " << gid << std::endl;
      //      std::cout << "LMUZAWAN " << lmuzawan << std::endl;
      // std::cout << "Gap " << kappa * gap << std::endl;

      localisincontact = true;

      double* normal = cnode->MoData().n();

      // compute lagrange multipliers and store into node
      for (int j = 0; j < dim; ++j)
        cnode->MoData().lm()[j] = (lmuzawan - coconstlaw_->Evaluate(kappa * gap)) * normal[j];
      // std::cout << "lm= " << cnode->MoData().lm()[1] << std::endl;
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
 |  Evaluate regularized tangential forces                gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::ConstitutivelawInterface::AssembleRegTangentForcesPenalty()
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // penalty parameter in tangential direction
  double pptan = IParams().get<double>("PENALTYPARAMTAN");
  double frcoeff = IParams().get<double>("FRCOEFF");

  INPAR::CONTACT::FrictionType ftype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(IParams(), "FRICTION");

  // loop over all slave row nodes on the current interface
  for (int i = 0; i < SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    // get some informatiom form the node
    double gap = cnode->CoData().Getg();
    int dim = cnode->NumDof();
    double kappa = cnode->CoData().Kappa();
    double* n = cnode->MoData().n();

    // Lagrange multiplier from Uzawa algorithm
    Epetra_SerialDenseMatrix lmuzawa(dim, 1);
    for (int k = 0; k < dim; ++k) lmuzawa(k, 0) = cnode->MoData().lmuzawa()[k];

    // Lagrange multiplier in normal direction
    double lmuzawan = 0.0;
    for (int k = 0; k < dim; ++k) lmuzawan += cnode->MoData().lmuzawa()[k] * cnode->MoData().n()[k];

    // tangential plane
    Epetra_SerialDenseMatrix tanplane(dim, dim);
    if (dim == 3)
    {
      tanplane(0, 0) = 1 - (n[0] * n[0]);
      tanplane(0, 1) = -(n[0] * n[1]);
      tanplane(0, 2) = -(n[0] * n[2]);
      tanplane(1, 0) = -(n[1] * n[0]);
      tanplane(1, 1) = 1 - (n[1] * n[1]);
      tanplane(1, 2) = -(n[1] * n[2]);

      tanplane(2, 0) = -(n[2] * n[0]);
      tanplane(2, 1) = -(n[2] * n[1]);
      tanplane(2, 2) = 1 - (n[2] * n[2]);
    }
    else if (dim == 2)
    {
      tanplane(0, 0) = 1 - (n[0] * n[0]);
      tanplane(0, 1) = -(n[0] * n[1]);

      tanplane(1, 0) = -(n[1] * n[0]);
      tanplane(1, 1) = 1 - (n[1] * n[1]);
    }
    else
      dserror("Error in AssembleTangentForces: Unknown dimension.");

    // evaluate traction
    Epetra_SerialDenseMatrix jumpvec(dim, 1);

    for (int i = 0; i < dim; i++) jumpvec(i, 0) = cnode->FriData().jump()[i];

    // evaluate kappa.pptan.jumptan
    Epetra_SerialDenseMatrix temptrac(dim, 1);
    temptrac.Multiply('N', 'N', kappa * pptan, tanplane, jumpvec, 0.0);

    // fill vector tractionold
    std::vector<double> tractionold(dim);
    for (int i = 0; i < dim; i++) tractionold[i] = cnode->FriData().tractionold()[i];

    // Evaluate trailtraction (tractionold+temptrac in penalty case)
    std::vector<double> trailtraction(dim);
    double magnitude = 0;
    for (int i = 0; i < dim; i++)
    {
      trailtraction[i] = tractionold[i] + temptrac(i, 0);
      magnitude += (trailtraction[i] * trailtraction[i]);
    }

    // evaluate magnitude of trailtraction
    magnitude = sqrt(magnitude);

    // evaluate maximal tangential traction
    double maxtantrac = frcoeff * (lmuzawan - coconstlaw_->Evaluate(kappa * gap));

    if (cnode->Active() == false)
    {
      // do nothing
      cnode->FriData().Slip() = false;
    }
    else if (cnode->Active() == true &&
             ((abs(maxtantrac) - magnitude >= 0) or ftype == INPAR::CONTACT::friction_stick))
    {
      // std::cout << "Node " << gid << " is stick" << std::endl;
      cnode->FriData().Slip() = false;

      // in the stick case, traction is trailtraction
      for (int i = 0; i < dim; i++) cnode->FriData().traction()[i] = trailtraction[i];

      // compute lagrange multipliers and store into node
      for (int j = 0; j < dim; ++j)
        cnode->MoData().lm()[j] = n[j] * (-coconstlaw_->Evaluate(kappa * gap)) + trailtraction[j];
    }
    else
    {
      // std::cout << "Node " << gid << " is slip" << std::endl;
      cnode->FriData().Slip() = true;

      // in the slip case, traction is evaluated with a return map algorithm
      for (int i = 0; i < dim; i++)
        cnode->FriData().traction()[i] = maxtantrac / magnitude * trailtraction[i];

      // compute lagrange multipliers and store into node
      for (int j = 0; j < dim; ++j)
        cnode->MoData().lm()[j] = n[j] * (-coconstlaw_->Evaluate(kappa * gap)) +
                                  maxtantrac / magnitude * trailtraction[j];
    }

    // linearization of contact forces (lagrange multipliers)
    // this consists the linearization of the tangential part,
    // the normal part was already done in AssembleRegNormalTraction

    // stick nodes
    if (cnode->Active() == true && cnode->FriData().Slip() == false)
    {
      /***************************************** tanplane.deriv(jump) ***/
      std::vector<std::map<int, double>>& derivjump = cnode->FriData().GetDerivJump();
      std::map<int, double>::iterator colcurr;
      GEN::pairedvector<int, double>::iterator _colcurr;

      // loop over dimensions
      for (int dimrow = 0; dimrow < cnode->NumDof(); ++dimrow)
      {
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr = derivjump[dim].begin(); colcurr != derivjump[dim].end(); ++colcurr)
          {
            int col = colcurr->first;
            double val = pptan * kappa * (colcurr->second) * tanplane(dimrow, dim);
            cnode->AddDerivZValue(dimrow, col, val);
          }
        }
      }

      /**************************************** deriv(tanplane).jump  ***/
      std::vector<GEN::pairedvector<int, double>>& derivn = cnode->CoData().GetDerivN();

      // loop over dimensions
      for (int dimrow = 0; dimrow < cnode->NumDof(); ++dimrow)
      {
        // loop over all entries of the current derivative map
        for (_colcurr = derivn[dimrow].begin(); _colcurr != derivn[dimrow].end(); ++_colcurr)
        {
          for (int dim = 0; dim < cnode->NumDof(); ++dim)
          {
            int col = _colcurr->first;
            double val =
                -pptan * kappa * (_colcurr->second) * n[dim] * cnode->FriData().jump()[dim];
            cnode->AddDerivZValue(dimrow, col, val);
          }
        }
      }

      // loop over dimensions
      for (int dim = 0; dim < cnode->NumDof(); ++dim)
      {
        // loop over all entries of the current derivative map
        for (_colcurr = derivn[dim].begin(); _colcurr != derivn[dim].end(); ++_colcurr)
        {
          for (int dimrow = 0; dimrow < cnode->NumDof(); ++dimrow)
          {
            int col = _colcurr->first;
            double val =
                -pptan * kappa * (_colcurr->second) * n[dimrow] * cnode->FriData().jump()[dim];
            cnode->AddDerivZValue(dimrow, col, val);
          }
        }
      }
    }
    // slip nodes
    else if (cnode->Active() == true && cnode->FriData().Slip() == true)
    {
      /******************** tanplane.deriv(jump).maxtantrac/magnidude ***/

      std::vector<std::map<int, double>>& derivjump = cnode->FriData().GetDerivJump();
      std::map<int, double>::iterator colcurr;
      GEN::pairedvector<int, double>::iterator _colcurr;

      // loop over dimensions
      for (int dimrow = 0; dimrow < cnode->NumDof(); ++dimrow)
      {
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr = derivjump[dim].begin(); colcurr != derivjump[dim].end(); ++colcurr)
          {
            int col = colcurr->first;
            double val =
                pptan * kappa * (colcurr->second) * tanplane(dimrow, dim) * maxtantrac / magnitude;
            cnode->AddDerivZValue(dimrow, col, val);
          }
        }
      }

      /******************** deriv(tanplane).jump.maxtantrac/magnitude ***/
      std::vector<GEN::pairedvector<int, double>>& derivn = cnode->CoData().GetDerivN();

      // loop over dimensions
      for (int dimrow = 0; dimrow < cnode->NumDof(); ++dimrow)
      {
        // loop over all entries of the current derivative map
        for (_colcurr = derivn[dimrow].begin(); _colcurr != derivn[dimrow].end(); ++_colcurr)
        {
          for (int dim = 0; dim < cnode->NumDof(); ++dim)
          {
            int col = _colcurr->first;
            double val = -pptan * kappa * (_colcurr->second) * n[dim] *
                         cnode->FriData().jump()[dim] * maxtantrac / magnitude;
            cnode->AddDerivZValue(dimrow, col, val);
          }
        }
      }
      // loop over dimensions
      for (int dim = 0; dim < cnode->NumDof(); ++dim)
      {
        // loop over all entries of the current derivative map
        for (_colcurr = derivn[dim].begin(); _colcurr != derivn[dim].end(); ++_colcurr)
        {
          for (int dimrow = 0; dimrow < cnode->NumDof(); ++dimrow)
          {
            int col = _colcurr->first;
            double val = -pptan * kappa * (_colcurr->second) * n[dimrow] *
                         cnode->FriData().jump()[dim] * maxtantrac / magnitude;
            cnode->AddDerivZValue(dimrow, col, val);
          }
        }
      }

      /******************** tanplane.jump.deriv(maxtantrac)/magnitude ***/
      std::map<int, double>& derivg = cnode->CoData().GetDerivG();
      std::map<int, double>::iterator gcurr;

      for (int j = 0; j < cnode->NumDof(); ++j)
      {
        for (gcurr = derivg.begin(); gcurr != derivg.end(); ++gcurr)
        {
          cnode->AddDerivZValue(j, gcurr->first,
              -frcoeff * coconstlaw_->EvaluateDeriv(kappa * gap) * kappa * (gcurr->second) *
                  trailtraction[j] / magnitude);
        }
      }

      /******************** tanplane.jump.maxtantrac/deriv(magnitude) ***/
      // vector double temp
      std::vector<double> temp(cnode->NumDof());
      for (int dim = 0; dim < cnode->NumDof(); ++dim)
        temp[dim] = -maxtantrac / (magnitude * magnitude) * trailtraction[dim];

      // loop over dimensions
      for (int dimout = 0; dimout < cnode->NumDof(); ++dimout)
      {
        double traction = 0;
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
          traction += tanplane(dimout, dim) * cnode->FriData().jump()[dim] * kappa * pptan;

        traction += tractionold[dimout];

        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr = derivjump[dim].begin(); colcurr != derivjump[dim].end(); ++colcurr)
          {
            int col = colcurr->first;
            double val =
                tanplane(dimout, dim) * pptan * kappa * (colcurr->second) * traction / magnitude;

            for (int dimrow = 0; dimrow < cnode->NumDof(); ++dimrow)
            {
              double val1 = val * temp[dimrow];
              cnode->AddDerivZValue(dimrow, col, val1);
            }
          }
        }
      }

      // loop over dimensions
      for (int dimout = 0; dimout < cnode->NumDof(); ++dimout)
      {
        double traction = 0;
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
          traction += tanplane(dimout, dim) * cnode->FriData().jump()[dim] * kappa * pptan;

        traction += tractionold[dimout];

        // loop over all entries of the current derivative map
        for (_colcurr = derivn[dimout].begin(); _colcurr != derivn[dimout].end(); ++_colcurr)
        {
          int col = _colcurr->first;

          for (int dim = 0; dim < cnode->NumDof(); ++dim)
          {
            double val = -_colcurr->second * n[dim] * cnode->FriData().jump()[dim] * traction /
                         magnitude * pptan * kappa;
            for (int dimrow = 0; dimrow < cnode->NumDof(); ++dimrow)
            {
              double val1 = val * temp[dimrow];
              cnode->AddDerivZValue(dimrow, col, val1);
            }
          }
        }
      }

      // loop over dimensions
      for (int dimout = 0; dimout < cnode->NumDof(); ++dimout)
      {
        double traction = 0;
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
          traction += tanplane(dimout, dim) * cnode->FriData().jump()[dim] * kappa * pptan;

        traction += tractionold[dimout];

        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          // loop over all entries of the current derivative map
          for (_colcurr = derivn[dim].begin(); _colcurr != derivn[dim].end(); ++_colcurr)
          {
            int col = _colcurr->first;

            double val = -_colcurr->second * n[dimout] * cnode->FriData().jump()[dim] * traction /
                         magnitude * pptan * kappa;

            for (int dimrow = 0; dimrow < cnode->NumDof(); ++dimrow)
            {
              double val1 = val * temp[dimrow];
              cnode->AddDerivZValue(dimrow, col, val1);
            }
          }
        }
      }
    }  // if Slip == true
    else
    {
      // clear tractions
      for (int j = 0; j < dim; ++j) cnode->MoData().lm()[j] = 0;
      // clear derivz
      cnode->CoData().GetDerivZ().clear();
    }
  }  // loop over active nodes
  return;
}
