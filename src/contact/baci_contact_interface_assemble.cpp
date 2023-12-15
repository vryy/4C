/*---------------------------------------------------------------------*/
/*! \file
\brief Contains assemble routines for multiple types of matrices and cases

\level 2


*/
/*---------------------------------------------------------------------*/


#include "baci_contact_element.H"
#include "baci_contact_friction_node.H"
#include "baci_contact_interface.H"
#include "baci_coupling_adapter.H"
#include "baci_linalg_utils_densematrix_communication.H"
#include "baci_linalg_utils_densematrix_multiply.H"
#include "baci_linalg_utils_sparse_algebra_assemble.H"
#include "baci_linalg_utils_sparse_algebra_manipulation.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------*
 |  Assemble slave coordinates (xs)                 gitterle 10/09|
 *----------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleSlaveCoord(Teuchos::RCP<Epetra_Vector>& xsmod)
{
  // loop over all slave nodes
  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    int dim = cnode->NumDof();

    CORE::LINALG::SerialDenseVector xspatial(dim);
    std::vector<int> dof(dim);
    std::vector<int> owner(dim);

    for (int k = 0; k < dim; ++k)
    {
      xspatial(k) = cnode->xspatial()[k];
      dof[k] = cnode->Dofs()[k];
      owner[k] = cnode->Owner();
    }

    // do assembly
    CORE::LINALG::Assemble(*xsmod, xspatial, dof, owner);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate regularized normal forces (nodes)                popp 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleRegNormalForces(
    bool& localisincontact, bool& localactivesetchange)
{
  // penalty parameter
  double pp = InterfaceParams().get<double>("PENALTYPARAM");

  // loop over all slave row nodes on the current interface
  for (int i = 0; i < SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
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
    if ((cnode->Active() == false) && (lmuzawan - kappa * pp * gap >= 0))
    {
      cnode->Active() = true;
      localactivesetchange = true;

      // std::cout << "node #" << gid << " is now active (";
      // for( int j=0; j<dim; j++)
      //  std::cout << " " << cnode->Dofs()[j] << " ";
      // std::cout << ") gap=" << gap << std::endl;
    }

    else if ((cnode->Active() == true) && (lmuzawan - kappa * pp * gap < 0))
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
      //      std::cout << "GAP " << gap << std::endl;

      localisincontact = true;

      double* normal = cnode->MoData().n();

      // compute lagrange multipliers and store into node
      for (int j = 0; j < dim; ++j)
        cnode->MoData().lm()[j] = (lmuzawan - kappa * pp * gap) * normal[j];

      // compute derivatives of lagrange multipliers and store into node

      // contribution of derivative of weighted gap
      std::map<int, double>& derivg = cnode->CoData().GetDerivG();
      std::map<int, double>::iterator gcurr;

      // contribution of derivative of normal
      std::vector<CORE::GEN::pairedvector<int, double>>& derivn = cnode->CoData().GetDerivN();
      CORE::GEN::pairedvector<int, double>::iterator ncurr;

      for (int j = 0; j < dim; ++j)
      {
        for (gcurr = derivg.begin(); gcurr != derivg.end(); ++gcurr)
          cnode->AddDerivZValue(j, gcurr->first, -kappa * pp * (gcurr->second) * normal[j]);
        for (ncurr = (derivn[j]).begin(); ncurr != (derivn[j]).end(); ++ncurr)
          cnode->AddDerivZValue(j, ncurr->first, -kappa * pp * gap * ncurr->second);
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
void CONTACT::CoInterface::AssembleRegTangentForcesPenalty()
{
  // penalty parameter in tangential direction
  double ppnor = InterfaceParams().get<double>("PENALTYPARAM");
  double pptan = InterfaceParams().get<double>("PENALTYPARAMTAN");
  double frcoeff = InterfaceParams().get<double>("FRCOEFF");

  INPAR::CONTACT::FrictionType ftype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(InterfaceParams(), "FRICTION");

  // loop over all slave row nodes on the current interface
  for (int i = 0; i < SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    // get some informatiom form the node
    double gap = cnode->CoData().Getg();
    int numdof = cnode->NumDof();
    double kappa = cnode->CoData().Kappa();
    double* n = cnode->MoData().n();

    // Lagrange multiplier from Uzawa algorithm
    CORE::LINALG::SerialDenseMatrix lmuzawa(numdof, 1);
    for (int k = 0; k < numdof; ++k) lmuzawa(k, 0) = cnode->MoData().lmuzawa()[k];

    // Lagrange multiplier in normal direction
    double lmuzawan = 0.0;
    for (int k = 0; k < numdof; ++k)
      lmuzawan += cnode->MoData().lmuzawa()[k] * cnode->MoData().n()[k];

    // tangential plane
    CORE::LINALG::SerialDenseMatrix tanplane(numdof, numdof);
    if (numdof == 3)
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
    else if (numdof == 2)
    {
      tanplane(0, 0) = 1 - (n[0] * n[0]);
      tanplane(0, 1) = -(n[0] * n[1]);

      tanplane(1, 0) = -(n[1] * n[0]);
      tanplane(1, 1) = 1 - (n[1] * n[1]);
    }
    else
      dserror("Error in AssembleTangentForces: Unknown dimension.");

    // evaluate traction
    CORE::LINALG::SerialDenseMatrix jumpvec(numdof, 1);

    for (int i = 0; i < numdof; i++) jumpvec(i, 0) = cnode->FriData().jump()[i];

    // evaluate kappa.pptan.jumptan
    CORE::LINALG::SerialDenseMatrix temptrac(numdof, 1);
    CORE::LINALG::multiply(0.0, temptrac, kappa * pptan, tanplane, jumpvec);

    // fill vector tractionold
    std::vector<double> tractionold(numdof);
    for (int i = 0; i < numdof; i++) tractionold[i] = cnode->FriData().tractionold()[i];

    // Evaluate trailtraction (tractionold+temptrac in penalty case)
    std::vector<double> trailtraction(numdof);
    double magnitude = 0;
    for (int i = 0; i < numdof; i++)
    {
      trailtraction[i] = tractionold[i] + temptrac(i, 0);
      magnitude += (trailtraction[i] * trailtraction[i]);
    }

    // evaluate magnitude of trailtraction
    magnitude = sqrt(magnitude);

    // evaluate maximal tangential traction
    double maxtantrac = frcoeff * (lmuzawan - kappa * ppnor * gap);

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
      for (int i = 0; i < numdof; i++) cnode->FriData().traction()[i] = trailtraction[i];

      // compute lagrange multipliers and store into node
      for (int j = 0; j < numdof; ++j)
        cnode->MoData().lm()[j] = n[j] * (-kappa * ppnor * gap) + trailtraction[j];
    }
    else
    {
      // std::cout << "Node " << gid << " is slip" << std::endl;
      cnode->FriData().Slip() = true;

      // in the slip case, traction is evaluated with a return map algorithm
      for (int i = 0; i < numdof; i++)
        cnode->FriData().traction()[i] = maxtantrac / magnitude * trailtraction[i];

      // compute lagrange multipliers and store into node
      for (int j = 0; j < numdof; ++j)
        cnode->MoData().lm()[j] =
            n[j] * (-kappa * ppnor * gap) + maxtantrac / magnitude * trailtraction[j];
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
      CORE::GEN::pairedvector<int, double>::iterator _colcurr;

      // loop over dimensions
      for (int dimrow = 0; dimrow < numdof; ++dimrow)
      {
        for (int dim = 0; dim < numdof; ++dim)
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
      std::vector<CORE::GEN::pairedvector<int, double>>& derivn = cnode->CoData().GetDerivN();

      // loop over dimensions
      for (int dimrow = 0; dimrow < numdof; ++dimrow)
      {
        // loop over all entries of the current derivative map
        for (_colcurr = derivn[dimrow].begin(); _colcurr != derivn[dimrow].end(); ++_colcurr)
        {
          for (int dim = 0; dim < numdof; ++dim)
          {
            int col = _colcurr->first;
            double val =
                -pptan * kappa * (_colcurr->second) * n[dim] * cnode->FriData().jump()[dim];
            cnode->AddDerivZValue(dimrow, col, val);
          }
        }
      }

      // loop over dimensions
      for (int dim = 0; dim < numdof; ++dim)
      {
        // loop over all entries of the current derivative map
        for (_colcurr = derivn[dim].begin(); _colcurr != derivn[dim].end(); ++_colcurr)
        {
          for (int dimrow = 0; dimrow < numdof; ++dimrow)
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
      CORE::GEN::pairedvector<int, double>::iterator _colcurr;

      // loop over dimensions
      for (int dimrow = 0; dimrow < numdof; ++dimrow)
      {
        for (int dim = 0; dim < numdof; ++dim)
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
      std::vector<CORE::GEN::pairedvector<int, double>>& derivn = cnode->CoData().GetDerivN();

      // loop over dimensions
      for (int dimrow = 0; dimrow < numdof; ++dimrow)
      {
        // loop over all entries of the current derivative map
        for (_colcurr = derivn[dimrow].begin(); _colcurr != derivn[dimrow].end(); ++_colcurr)
        {
          for (int dim = 0; dim < numdof; ++dim)
          {
            int col = _colcurr->first;
            double val = -pptan * kappa * (_colcurr->second) * n[dim] *
                         cnode->FriData().jump()[dim] * maxtantrac / magnitude;
            cnode->AddDerivZValue(dimrow, col, val);
          }
        }
      }
      // loop over dimensions
      for (int dim = 0; dim < numdof; ++dim)
      {
        // loop over all entries of the current derivative map
        for (_colcurr = derivn[dim].begin(); _colcurr != derivn[dim].end(); ++_colcurr)
        {
          for (int dimrow = 0; dimrow < numdof; ++dimrow)
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

      for (int j = 0; j < numdof; ++j)
      {
        for (gcurr = derivg.begin(); gcurr != derivg.end(); ++gcurr)
        {
          cnode->AddDerivZValue(j, gcurr->first,
              -frcoeff * kappa * ppnor * (gcurr->second) * trailtraction[j] / magnitude);
        }
      }

      /******************** tanplane.jump.maxtantrac/deriv(magnitude) ***/
      // vector double temp
      std::vector<double> temp(numdof);
      for (int dim = 0; dim < numdof; ++dim)
        temp[dim] = -maxtantrac / (magnitude * magnitude) * trailtraction[dim];

      // loop over dimensions
      for (int dimout = 0; dimout < numdof; ++dimout)
      {
        double traction = 0;
        for (int dim = 0; dim < numdof; ++dim)
          traction += tanplane(dimout, dim) * cnode->FriData().jump()[dim] * kappa * pptan;

        traction += tractionold[dimout];

        for (int dim = 0; dim < numdof; ++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr = derivjump[dim].begin(); colcurr != derivjump[dim].end(); ++colcurr)
          {
            int col = colcurr->first;
            double val =
                tanplane(dimout, dim) * pptan * kappa * (colcurr->second) * traction / magnitude;

            for (int dimrow = 0; dimrow < numdof; ++dimrow)
            {
              double val1 = val * temp[dimrow];
              cnode->AddDerivZValue(dimrow, col, val1);
            }
          }
        }
      }

      // loop over dimensions
      for (int dimout = 0; dimout < numdof; ++dimout)
      {
        double traction = 0;
        for (int dim = 0; dim < numdof; ++dim)
          traction += tanplane(dimout, dim) * cnode->FriData().jump()[dim] * kappa * pptan;

        traction += tractionold[dimout];

        // loop over all entries of the current derivative map
        for (_colcurr = derivn[dimout].begin(); _colcurr != derivn[dimout].end(); ++_colcurr)
        {
          int col = _colcurr->first;

          for (int dim = 0; dim < numdof; ++dim)
          {
            double val = -_colcurr->second * n[dim] * cnode->FriData().jump()[dim] * traction /
                         magnitude * pptan * kappa;
            for (int dimrow = 0; dimrow < numdof; ++dimrow)
            {
              double val1 = val * temp[dimrow];
              cnode->AddDerivZValue(dimrow, col, val1);
            }
          }
        }
      }

      // loop over dimensions
      for (int dimout = 0; dimout < numdof; ++dimout)
      {
        double traction = 0;
        for (int dim = 0; dim < numdof; ++dim)
          traction += tanplane(dimout, dim) * cnode->FriData().jump()[dim] * kappa * pptan;

        traction += tractionold[dimout];

        for (int dim = 0; dim < numdof; ++dim)
        {
          // loop over all entries of the current derivative map
          for (_colcurr = derivn[dim].begin(); _colcurr != derivn[dim].end(); ++_colcurr)
          {
            int col = _colcurr->first;

            double val = -_colcurr->second * n[dimout] * cnode->FriData().jump()[dim] * traction /
                         magnitude * pptan * kappa;

            for (int dimrow = 0; dimrow < numdof; ++dimrow)
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
      for (int j = 0; j < numdof; ++j) cnode->MoData().lm()[j] = 0;
      // clear derivz
      cnode->CoData().GetDerivZ().clear();
    }
  }  // loop over active nodes
}

/*----------------------------------------------------------------------*
 |  Evaluate regularized tangential forces                gitterle 10/09|
 |  (Uzawa Aug. Lagr.)                                                  |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleRegTangentForcesUzawa()
{
  // penalty parameter in tangential direction
  double ppnor = InterfaceParams().get<double>("PENALTYPARAM");
  double pptan = InterfaceParams().get<double>("PENALTYPARAMTAN");
  double frcoeff = InterfaceParams().get<double>("FRCOEFF");

  INPAR::CONTACT::FrictionType ftype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(InterfaceParams(), "FRICTION");

  // loop over all slave row nodes on the current interface
  for (int i = 0; i < SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    // get some informatiom form the node
    double gap = cnode->CoData().Getg();
    int dim = cnode->NumDof();
    double kappa = cnode->CoData().Kappa();
    double* n = cnode->MoData().n();

    // Lagrange multiplier from Uzawa algorithm
    CORE::LINALG::SerialDenseMatrix lmuzawa(dim, 1);
    for (int k = 0; k < dim; ++k) lmuzawa(k, 0) = cnode->MoData().lmuzawa()[k];

    // Lagrange multiplier in normal direction
    double lmuzawan = 0.0;
    for (int k = 0; k < dim; ++k) lmuzawan += cnode->MoData().lmuzawa()[k] * cnode->MoData().n()[k];

    // tangential plane
    CORE::LINALG::SerialDenseMatrix tanplane(dim, dim);
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

    // Lagrange multiplier in tangential direction
    CORE::LINALG::SerialDenseMatrix lmuzawatan(dim, 1);
    CORE::LINALG::multiply(lmuzawatan, tanplane, lmuzawa);

    // evaluate traction
    CORE::LINALG::SerialDenseMatrix jumpvec(dim, 1);

    for (int i = 0; i < dim; i++) jumpvec(i, 0) = cnode->FriData().jump()[i];

    // evaluate kappa.pptan.jumptan
    CORE::LINALG::SerialDenseMatrix temptrac(dim, 1);
    CORE::LINALG::multiply(0.0, temptrac, kappa * pptan, tanplane, jumpvec);

    // Evaluate trailtraction
    std::vector<double> trailtraction(dim);
    double magnitude = 0;
    for (int i = 0; i < dim; i++)
    {
      trailtraction[i] = lmuzawatan(i, 0) + temptrac(i, 0);
      magnitude += (trailtraction[i] * trailtraction[i]);
    }

    // evaluate magnitude of trailtraction
    magnitude = sqrt(magnitude);

    // evaluate maximal tangential traction
    double maxtantrac = frcoeff * (lmuzawan - kappa * ppnor * gap);

    if (cnode->Active() == false)
    {
    }
    else if (cnode->Active() == true &&
             ((abs(maxtantrac) - magnitude >= 0) or ftype == INPAR::CONTACT::friction_stick))
    {
      // std::cout << "Node " << gid << " is stick" << std::endl;
      cnode->FriData().Slip() = false;

      // compute lagrange multipliers and store into node
      for (int j = 0; j < dim; ++j)
        cnode->MoData().lm()[j] = n[j] * (lmuzawan - kappa * ppnor * gap) + trailtraction[j];
    }
    else
    {
      // std::cout << "Node " << gid << " is slip" << std::endl;
      cnode->FriData().Slip() = true;

      // compute lagrange multipliers and store into node
      for (int j = 0; j < dim; ++j)
        cnode->MoData().lm()[j] =
            n[j] * (lmuzawan - kappa * ppnor * gap) + trailtraction[j] * maxtantrac / magnitude;
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
      CORE::GEN::pairedvector<int, double>::iterator _colcurr;
      const int numdof = cnode->NumDof();

      // loop over dimensions
      for (int dimrow = 0; dimrow < numdof; ++dimrow)
      {
        for (int dim = 0; dim < numdof; ++dim)
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

      /******************************* deriv(tanplane).(lmuzawa+jump) ***/
      std::vector<CORE::GEN::pairedvector<int, double>>& derivn = cnode->CoData().GetDerivN();

      // loop over dimensions
      for (int dimrow = 0; dimrow < numdof; ++dimrow)
      {
        // loop over all entries of the current derivative map
        for (_colcurr = derivn[dimrow].begin(); _colcurr != derivn[dimrow].end(); ++_colcurr)
        {
          for (int dim = 0; dim < numdof; ++dim)
          {
            int col = _colcurr->first;
            double val =
                -pptan * kappa * (_colcurr->second) * n[dim] * (cnode->FriData().jump()[dim]);
            val = val - (_colcurr->second) * n[dim] * (cnode->MoData().lmuzawa()[dim]);
            cnode->AddDerivZValue(dimrow, col, val);
          }
        }
      }

      // loop over dimensions
      for (int dim = 0; dim < numdof; ++dim)
      {
        // loop over all entries of the current derivative map
        for (_colcurr = derivn[dim].begin(); _colcurr != derivn[dim].end(); ++_colcurr)
        {
          for (int dimrow = 0; dimrow < numdof; ++dimrow)
          {
            int col = _colcurr->first;
            double val =
                -pptan * kappa * (_colcurr->second) * n[dimrow] * (cnode->FriData().jump()[dim]);
            val = val - (_colcurr->second) * n[dimrow] * (cnode->MoData().lmuzawa()[dim]);
            cnode->AddDerivZValue(dimrow, col, val);
          }
        }
      }
    }

    // slip nodes
    else if (cnode->Active() == true && cnode->FriData().Slip() == true)
    {
      /***************************************** tanplane.deriv(jump) ***/
      std::vector<std::map<int, double>>& derivjump = cnode->FriData().GetDerivJump();
      std::map<int, double>::iterator colcurr;
      CORE::GEN::pairedvector<int, double>::iterator _colcurr;
      const int numdof = cnode->NumDof();

      // loop over dimensions
      for (int dimrow = 0; dimrow < numdof; ++dimrow)
      {
        for (int dim = 0; dim < numdof; ++dim)
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

      /******************************* deriv(tanplane).(lmuzawa+jump) ***/
      std::vector<CORE::GEN::pairedvector<int, double>>& derivn = cnode->CoData().GetDerivN();

      // loop over dimensions
      for (int dimrow = 0; dimrow < numdof; ++dimrow)
      {
        // loop over all entries of the current derivative map
        for (_colcurr = derivn[dimrow].begin(); _colcurr != derivn[dimrow].end(); ++_colcurr)
        {
          for (int dim = 0; dim < numdof; ++dim)
          {
            int col = _colcurr->first;
            double val =
                -pptan * kappa * (_colcurr->second) * n[dim] * cnode->FriData().jump()[dim];
            val = (val - (_colcurr->second) * n[dim] * (cnode->MoData().lmuzawa()[dim])) *
                  maxtantrac / magnitude;
            cnode->AddDerivZValue(dimrow, col, val);
          }
        }
      }

      // loop over dimensions
      for (int dim = 0; dim < numdof; ++dim)
      {
        // loop over all entries of the current derivative map
        for (_colcurr = derivn[dim].begin(); _colcurr != derivn[dim].end(); ++_colcurr)
        {
          for (int dimrow = 0; dimrow < numdof; ++dimrow)
          {
            int col = _colcurr->first;
            double val =
                -pptan * kappa * (_colcurr->second) * n[dimrow] * cnode->FriData().jump()[dim];
            val = (val - (_colcurr->second) * n[dimrow] * (cnode->MoData().lmuzawa()[dim])) *
                  maxtantrac / magnitude;
            cnode->AddDerivZValue(dimrow, col, val);
          }
        }
      }

      /******************** tanplane.jump.deriv(maxtantrac)/magnitude ***/
      std::map<int, double>& derivg = cnode->CoData().GetDerivG();
      std::map<int, double>::iterator gcurr;

      for (int j = 0; j < numdof; ++j)
      {
        for (gcurr = derivg.begin(); gcurr != derivg.end(); ++gcurr)
        {
          cnode->AddDerivZValue(j, gcurr->first,
              -frcoeff * kappa * ppnor * (gcurr->second) * trailtraction[j] / magnitude);
        }
      }

      for (int j = 0; j < numdof; ++j)
      {
        for (_colcurr = (derivn[j]).begin(); _colcurr != (derivn[j]).end(); ++_colcurr)
        {
          for (int k = 0; k < numdof; ++k)
          {
            double val =
                frcoeff * (_colcurr->second) * lmuzawa(j, 0) * trailtraction[k] / magnitude;
            cnode->AddDerivZValue(k, _colcurr->first, val);
          }
        }
      }

      /******************** tanplane.jump.maxtantrac/deriv(magnitude) ***/
      // vector double temp
      std::vector<double> temp(numdof);
      for (int dim = 0; dim < numdof; ++dim)
        temp[dim] = -maxtantrac / (magnitude * magnitude) * trailtraction[dim];

      // loop over dimensions
      for (int dimout = 0; dimout < numdof; ++dimout)
      {
        double traction = 0;
        for (int dim = 0; dim < numdof; ++dim)
          traction += tanplane(dimout, dim) *
                      (lmuzawa(dim, 0) + cnode->FriData().jump()[dim] * kappa * pptan);

        for (int dim = 0; dim < numdof; ++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr = derivjump[dim].begin(); colcurr != derivjump[dim].end(); ++colcurr)
          {
            int col = colcurr->first;
            double val =
                tanplane(dimout, dim) * pptan * kappa * (colcurr->second) * traction / magnitude;

            for (int dimrow = 0; dimrow < numdof; ++dimrow)
            {
              double val1 = val * temp[dimrow];
              cnode->AddDerivZValue(dimrow, col, val1);
            }
          }
        }
      }

      // loop over dimensions
      for (int dimout = 0; dimout < numdof; ++dimout)
      {
        double traction = 0;
        for (int dim = 0; dim < numdof; ++dim)
          traction += tanplane(dimout, dim) *
                      (lmuzawa(dim, 0) + cnode->FriData().jump()[dim] * kappa * pptan);

        // loop over all entries of the current derivative map
        for (_colcurr = derivn[dimout].begin(); _colcurr != derivn[dimout].end(); ++_colcurr)
        {
          int col = _colcurr->first;

          for (int dim = 0; dim < numdof; ++dim)
          {
            double val = -_colcurr->second * n[dim] *
                         (lmuzawa(dim, 0) + cnode->FriData().jump()[dim] * pptan * kappa) *
                         traction / magnitude;
            for (int dimrow = 0; dimrow < numdof; ++dimrow)
            {
              double val1 = val * temp[dimrow];
              cnode->AddDerivZValue(dimrow, col, val1);
            }
          }
        }
      }

      // loop over dimensions
      for (int dimout = 0; dimout < numdof; ++dimout)
      {
        double traction = 0;
        for (int dim = 0; dim < numdof; ++dim)
          traction += tanplane(dimout, dim) *
                      (lmuzawa(dim, 0) + cnode->FriData().jump()[dim] * kappa * pptan);

        for (int dim = 0; dim < numdof; ++dim)
        {
          // loop over all entries of the current derivative map
          for (_colcurr = derivn[dim].begin(); _colcurr != derivn[dim].end(); ++_colcurr)
          {
            int col = _colcurr->first;
            double val = -_colcurr->second * n[dimout] *
                         (lmuzawa(dim, 0) + cnode->FriData().jump()[dim] * pptan * kappa) *
                         traction / magnitude;

            for (int dimrow = 0; dimrow < numdof; ++dimrow)
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
}

/*----------------------------------------------------------------------*
 |  Assemble derivatives of lagrange multipliers              popp 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinZ(CORE::LINALG::SparseMatrix& linzglobal)
{
  // loop over all slave nodes (row map)
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID()) dserror("AssembleLinZ: Node ownership inconsistency!");

    // derivz is the std::vector<map> we want to assemble
    std::vector<std::map<int, double>>& derivz = cnode->CoData().GetDerivZ();

    if ((int)derivz.size() > 0)
    {
      int rowsize = cnode->NumDof();
      int colsize = (int)derivz[0].size();

      // consistency check
      for (int j = 0; j < rowsize - 1; ++j)
        if ((int)derivz[j].size() != (int)derivz[j + 1].size())
          dserror("AssembleLinZ: Column dim. of nodal derivz-map is inconsistent!");

      std::map<int, double>::iterator colcurr;

      // loop over dofs
      for (int k = 0; k < rowsize; ++k)
      {
        int row = cnode->Dofs()[k];  // row index equals global dof index of this #i node's dof k
        int l = 0;

        // loop over all directional derivative entries using the map iterator
        for (colcurr = derivz[k].begin(); colcurr != derivz[k].end(); ++colcurr)
        {
          int col =
              colcurr->first;  // col index equals global id of directional derivative component ,l
          double val = colcurr->second;
          linzglobal.Assemble(val, row, col);
          l++;
        }

        if (l != colsize) dserror("AssembleLinZ: l = %i but colsize = %i", k, colsize);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix with nodal tangents or/and normals         popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleTN(Teuchos::RCP<CORE::LINALG::SparseMatrix> tglobal,
    Teuchos::RCP<CORE::LINALG::SparseMatrix> nglobal)
{
  // nothing to do if no active nodes
  if (activenodes_ == Teuchos::null) return;

  if (Dim() != 2 && Dim() != 3) dserror("Dim() must be either 2 or 3!");

  // loop over all active slave nodes of the interface
  for (int i = 0; i < activenodes_->NumMyElements(); ++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    auto* cnode = dynamic_cast<CoNode*>(node);
    const int numdof = cnode->NumDof();

    if (cnode->Owner() != Comm().MyPID()) dserror("AssembleTN: Node ownership inconsistency!");

    if (tglobal != Teuchos::null)
    {
      if (constr_direction_ == INPAR::CONTACT::constr_xyz)
      {
        if (Dim() == 2)
        {
          // prepare assembly
          std::vector<int> lmrowT(cnode->NumDof());
          std::vector<int> lmrowownerT(cnode->NumDof());
          std::vector<int> lmcol(cnode->NumDof());

          CORE::LINALG::SerialDenseMatrix Tnode(numdof, numdof);
          for (int i = 0; i < numdof; ++i)
          {
            lmrowT[i] = cnode->Dofs()[i];
            lmrowownerT[i] = cnode->Owner();
            for (int j = 0; j < numdof; ++j)
            {
              lmcol[j] = cnode->Dofs()[j];
              Tnode(i, j) = cnode->CoData().txi()[i] * cnode->CoData().txi()[j];
            }
          }
          tglobal->Assemble(-1, Tnode, lmrowT, lmrowownerT, lmcol);
        }
        else if (Dim() == 3)
        {
          std::vector<int> lmrowT(cnode->NumDof());
          std::vector<int> lmrowownerT(cnode->NumDof());
          std::vector<int> lmcol(cnode->NumDof());

          CORE::LINALG::SerialDenseMatrix Tnode(numdof, numdof);

          for (int i = 0; i < numdof; ++i)
          {
            lmrowT[i] = cnode->Dofs()[i];
            lmrowownerT[i] = cnode->Owner();
            for (int j = 0; j < numdof; ++j)
            {
              lmcol[j] = cnode->Dofs()[j];
              Tnode(i, j) = cnode->CoData().txi()[i] * cnode->CoData().txi()[j] +
                            cnode->CoData().teta()[i] * cnode->CoData().teta()[j];
            }
          }
          tglobal->Assemble(-1, Tnode, lmrowT, lmrowownerT, lmcol);
        }
        else
          dserror("Dim() must be either 2D or 3D");
      }
      else
      {
        if (Dim() == 2)
        {
          // prepare assembly
          int colsize = numdof;
          std::vector<int> lmrowT(1);
          std::vector<int> lmrowownerT(1);
          std::vector<int> lmcol(colsize);

          lmrowT[0] = activet_->GID(i);
          lmrowownerT[0] = cnode->Owner();

          /**************************************************** T-matrix ******/
          CORE::LINALG::SerialDenseMatrix Tnode(1, colsize);

          for (int j = 0; j < colsize; ++j)
          {
            lmcol[j] = cnode->Dofs()[j];
            Tnode(0, j) = cnode->CoData().txi()[j];
          }

          // assemble into matrix of normal vectors T
          tglobal->Assemble(-1, Tnode, lmrowT, lmrowownerT, lmcol);
        }

        else if (Dim() == 3)
        {
          // prepare assembly
          int colsize = numdof;
          std::vector<int> lmrowT(2);
          std::vector<int> lmrowownerT(2);
          std::vector<int> lmcol(colsize);

          lmrowT[0] = activet_->GID(2 * i);
          lmrowT[1] = activet_->GID(2 * i + 1);
          lmrowownerT[0] = cnode->Owner();
          lmrowownerT[1] = cnode->Owner();

          /**************************************************** T-matrix ******/
          CORE::LINALG::SerialDenseMatrix Tnode(2, colsize);

          for (int j = 0; j < colsize; ++j)
          {
            lmcol[j] = cnode->Dofs()[j];
            Tnode(0, j) = cnode->CoData().txi()[j];
            Tnode(1, j) = cnode->CoData().teta()[j];
          }

          // assemble into matrix of normal vectors T
          tglobal->Assemble(-1, Tnode, lmrowT, lmrowownerT, lmcol);
        }
        else
          dserror("Dim() must be either 2D or 3D");
      }
    }

    if (nglobal != Teuchos::null)
    {
      // nodal normal
      double* nodalnormal = cnode->MoData().n();

      int row = activen_->GID(i);

      // add normal to corresponding row in global matrix
      for (int k = 0; k < numdof; ++k)
        nglobal->Assemble(
            nodalnormal[k], row, cnode->Dofs()[k]);  // use the first dof for normal direction!!!
    }
  }
}

/*----------------------------------------------------------------------*
 |  Assemble matrix S containing gap g~ derivatives           popp 02/09|
 |  PS: "AssembleS" is an outdated name which could make                |
 |  you confused.                                                       |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleS(CORE::LINALG::SparseMatrix& sglobal)
{
  // nothing to do if no active nodes
  if (activenodes_ == Teuchos::null) return;

  // loop over all active slave nodes of the interface
  for (int i = 0; i < activenodes_->NumMyElements(); ++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    auto* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID()) dserror("AssembleS: Node ownership inconsistency!");

    // prepare assembly
    std::map<int, double>& dgmap = cnode->CoData().GetDerivG();
    std::map<int, double>::iterator colcurr;
    int row = activen_->GID(i);

    for (colcurr = dgmap.begin(); colcurr != dgmap.end(); ++colcurr)
    {
      int col = colcurr->first;
      double val = colcurr->second;

      // do not assemble zeros into s matrix
      if (constr_direction_ == INPAR::CONTACT::constr_xyz)
      {
        for (int j = 0; j < cnode->NumDof(); j++)
          if (abs(val * cnode->MoData().n()[j]) > 1.0e-12)
            sglobal.Assemble(val * cnode->MoData().n()[j], cnode->Dofs()[j], col);
      }
      else if (abs(val) > 1.0e-12)
        sglobal.Assemble(val, row, col);
    }
  }  // for (int i=0;i<activenodes_->NumMyElements();++i)
}

/*----------------------------------------------------------------------*
 |  Assemble tangent deriv or/and normal deriv matrix         popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleTNderiv(Teuchos::RCP<CORE::LINALG::SparseMatrix> tderivglobal,
    Teuchos::RCP<CORE::LINALG::SparseMatrix> nderivglobal, bool usePoroLM)
{
  // nothing to do if no active nodes
  if (activenodes_ == Teuchos::null) return;

  if (Dim() != 2 && Dim() != 3) dserror("Dim() must be either 2 or 3!");

  // loop over all active slave nodes of the interface
  for (int i = 0; i < activenodes_->NumMyElements(); ++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())  // move this check into debug?
      dserror("AssembleTNderiv: Node ownership inconsistency!");

    if (tderivglobal != Teuchos::null)  // assemble tangential derivs?
    {
      std::vector<CORE::GEN::pairedvector<int, double>>& dtmap = cnode->CoData().GetDerivTxi();
      CORE::GEN::pairedvector<int, double>::iterator colcurr;

      for (int dim = 0; dim < Dim() - 1; ++dim)  // for both tangents
      {
        if (dim == 1)  // just for 3d case, 2nd tangent
          dtmap = cnode->CoData().GetDerivTeta();
        int colsize = (int)dtmap[0].size();
        int mapsize = (int)dtmap.size();

        int row = activet_->GID((Dim() - 1) * i + dim);

        if (Dim() == 2 && mapsize == 3) mapsize = 2;  //??

        for (int j = 0; j < mapsize - 1; ++j)  // move this check into debug?
          if ((int)dtmap[j].size() != (int)dtmap[j + 1].size())
          {
            std::cout << "size j = " << dtmap[j].size() << "   size j+1 = " << dtmap[j + 1].size()
                      << std::endl;
            dserror("AssembleTNderiv: Column dim. of nodal DerivT-map is inconsistent!");
          }

        // begin assembly of Tderiv-matrix
        // std::cout << std::endl << "->Assemble P for Node ID: " << cnode->Id() << std::endl;

        // loop over all derivative maps (=dimensions)
        for (int j = 0; j < mapsize; ++j)
        {
          int k = 0;

          // loop over all entries of the current derivative map
          for (colcurr = dtmap[j].begin(); colcurr != dtmap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
            double val;
            if (!usePoroLM)
              val = cnode->MoData().lm()[j] * (colcurr->second);
            else
              val = cnode->CoPoroData().poroLM()[j] * (colcurr->second);
            // std::cout << "lm[" << j << "]=" << cnode->MoData().lm()[j] << " deriv=" <<
            // colcurr->second << std::endl; std::cout << "Assemble P: " << row << " " << col << " "
            // << val << std::endl;
            // do not assemble zeros into P matrix

            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int i = 0; i < cnode->NumDof(); ++i)
              {
                if (abs(val) > 1.0e-12)
                {
                  double t;
                  if (dim == 0)
                    t = cnode->CoData().txi()[i];
                  else
                    t = cnode->CoData().teta()[i];
                  tderivglobal->Assemble(val * t, cnode->Dofs()[i], col);
                }
              }
            }
            else if (abs(val) > 1.0e-12)
              tderivglobal->Assemble(val, row, col);
            ++k;
          }

          if (k != colsize) dserror("AssembleTNderiv: k = %i but colsize = %i", k, colsize);
        }
      }
    }

    if (nderivglobal != Teuchos::null)  // assemble normal derivs?
    {
      std::vector<CORE::GEN::pairedvector<int, double>>& dnmap = cnode->CoData().GetDerivN();
      CORE::GEN::pairedvector<int, double>::iterator colcurr;

      int colsize = (int)dnmap[0].size();
      int mapsize = (int)dnmap.size();

      int row = activen_->GID(i);

      if (Dim() == 2 && mapsize == 3) mapsize = 2;  //??

      for (int j = 0; j < mapsize - 1; ++j)  // move this check into debug?
        if ((int)dnmap[j].size() != (int)dnmap[j + 1].size())
          dserror("AssembleTNderiv: Column dim. of nodal DerivN-map is inconsistent!");

      // loop over all derivative maps (=dimensions)
      for (int j = 0; j < mapsize; ++j)
      {
        int k = 0;

        // loop over all entries of the current derivative map
        for (colcurr = dnmap[j].begin(); colcurr != dnmap[j].end(); ++colcurr)
        {
          int col = colcurr->first;
          double val;
          if (!usePoroLM)
            val = cnode->MoData().lm()[j] * (colcurr->second);
          else
            val = cnode->CoPoroData().poroLM()[j] * (colcurr->second);

          // do not assemble zeros into P matrix

          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int i = 0; i < cnode->NumDof(); ++i)
              if (abs(val) > 1.0e-12)
                nderivglobal->Assemble(val * cnode->MoData().n()[i], cnode->Dofs()[i], col);
          }
          else if (abs(val) > 1.0e-12)
            nderivglobal->Assemble(val, row, col);
          ++k;
        }

        if (k != colsize) dserror("AssembleTNderiv: k = %i but colsize = %i", k, colsize);
      }
    }

  }  // for (int i=0;i<activenodes_->NumMyElements();++i)

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrices LinD containing fc derivatives          popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinD(CORE::LINALG::SparseMatrix& lindglobal, bool usePoroLM)
{
  /**********************************************************************/
  // NEW VERSION (09/2010): No more communication, thanks to FE_MATRIX!
  /**********************************************************************/
  // we have: D_jk,c with j = Lagrange multiplier slave dof
  //                 with k = Displacement slave dof
  //                 with c = Displacement slave or master dof
  // we compute (LinD)_kc = D_jk,c * z_j
  /**********************************************************************/

  // loop over all LM slave nodes (row map)
  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);
    int dim = cnode->NumDof();

    // Mortar matrix D and M derivatives
    std::map<int, std::map<int, double>>& dderiv = cnode->CoData().GetDerivD();

    // current Lagrange multipliers
    double* lm;
    if (!usePoroLM)
      lm = cnode->MoData().lm();
    else
      lm = cnode->CoPoroData().poroLM();

    // get sizes and iterator start
    int slavesize = (int)dderiv.size();
    std::map<int, std::map<int, double>>::iterator scurr = dderiv.begin();

    /********************************************** LinDMatrix **********/
    // loop over all DISP slave nodes in the DerivD-map of the current LM slave node
    for (int k = 0; k < slavesize; ++k)
    {
      int sgid = scurr->first;
      ++scurr;

      DRT::Node* snode = idiscret_->gNode(sgid);
      if (!snode) dserror("Cannot find node with gid %", sgid);
      CoNode* csnode = dynamic_cast<CoNode*>(snode);

      // Mortar matrix D derivatives
      std::map<int, double>& thisdderiv = cnode->CoData().GetDerivD()[sgid];
      int mapsize = (int)(thisdderiv.size());

      // inner product D_{jk,c} * z_j for index j
      for (int prodj = 0; prodj < dim; ++prodj)
      {
        int row = csnode->Dofs()[prodj];
        std::map<int, double>::iterator scolcurr = thisdderiv.begin();

        // loop over all directional derivative entries
        for (int c = 0; c < mapsize; ++c)
        {
          int col = scolcurr->first;
          double val = lm[prodj] * (scolcurr->second);
          ++scolcurr;

          // owner of LM slave node can do the assembly, although it actually
          // might not own the corresponding rows in lindglobal (DISP slave node)
          // (FE_MATRIX automatically takes care of non-local assembly inside!!!)
          // std::cout << "Assemble LinD: " << row << " " << col << " " << val << std::endl;
          if (abs(val) > 1.0e-12) lindglobal.FEAssemble(val, row, col);
        }

        // check for completeness of DerivD-Derivatives-iteration
        if (scolcurr != thisdderiv.end())
          dserror("AssembleLinDM: Not all derivative entries of DerivD considered!");
      }
    }

    // check for completeness of DerivD-Slave-iteration
    if (scurr != dderiv.end())
      dserror("AssembleLinDM: Not all DISP slave entries of DerivD considered!");
    /******************************** Finished with LinDMatrix **********/
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Assemble matrices LinM containing fc derivatives          popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinM(CORE::LINALG::SparseMatrix& linmglobal, bool usePoroLM)
{
  /**********************************************************************/
  // NEW VERSION (09/2010): No more communication, thanks to FE_MATRIX!
  /**********************************************************************/
  // we have: M_jl,c with j = Lagrange multiplier slave dof
  //                 with l = Displacement master dof
  //                 with c = Displacement slave or master dof
  // we compute (LinM)_lc = M_jl,c * z_j
  /**********************************************************************/

  // loop over all LM slave nodes (row map)
  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);
    int dim = cnode->NumDof();

    // Mortar matrix D and M derivatives
    std::map<int, std::map<int, double>>& mderiv = cnode->CoData().GetDerivM();

    // current Lagrange multipliers
    double* lm;
    if (!usePoroLM)
      lm = cnode->MoData().lm();
    else
      lm = cnode->CoPoroData().poroLM();

    // get sizes and iterator start
    int mastersize = (int)mderiv.size();
    std::map<int, std::map<int, double>>::iterator mcurr = mderiv.begin();

    /********************************************** LinMMatrix **********/
    // loop over all master nodes in the DerivM-map of the current LM slave node
    for (int l = 0; l < mastersize; ++l)
    {
      int mgid = mcurr->first;
      ++mcurr;

      DRT::Node* mnode = idiscret_->gNode(mgid);
      if (!mnode) dserror("Cannot find node with gid %", mgid);
      CoNode* cmnode = dynamic_cast<CoNode*>(mnode);

      // Mortar matrix M derivatives
      std::map<int, double>& thismderiv = cnode->CoData().GetDerivM()[mgid];
      int mapsize = (int)(thismderiv.size());

      // inner product M_{jl,c} * z_j for index j
      for (int prodj = 0; prodj < dim; ++prodj)
      {
        int row = cmnode->Dofs()[prodj];
        std::map<int, double>::iterator mcolcurr = thismderiv.begin();

        // loop over all directional derivative entries
        for (int c = 0; c < mapsize; ++c)
        {
          int col = mcolcurr->first;
          double val = lm[prodj] * (mcolcurr->second);
          ++mcolcurr;

          // owner of LM slave node can do the assembly, although it actually
          // might not own the corresponding rows in lindglobal (DISP slave node)
          // (FE_MATRIX automatically takes care of non-local assembly inside!!!)
          // std::cout << "Assemble LinM: " << row << " " << col << " " << val << std::endl;
          if (abs(val) > 1.0e-12) linmglobal.FEAssemble(-val, row, col);
        }

        // check for completeness of DerivM-Derivatives-iteration
        if (mcolcurr != thismderiv.end())
          dserror("AssembleLinDM: Not all derivative entries of DerivM considered!");
      }
    }

    // check for completeness of DerivM-Master-iteration
    if (mcurr != mderiv.end())
      dserror("AssembleLinDM: Not all master entries of DerivM considered!");
    /******************************** Finished with LinMMatrix **********/
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Assemble matrices LinDM containing fc derivatives        farah 02/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinDM(
    CORE::LINALG::SparseMatrix& lindglobal, CORE::LINALG::SparseMatrix& linmglobal, bool usePoroLM)
{
  // call both sub functions
  AssembleLinD(lindglobal, usePoroLM);
  AssembleLinM(linmglobal, usePoroLM);

  return;
}


/*----------------------------------------------------------------------*
 |  Assemble normal weighted gap                              popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleG(Epetra_Vector& gglobal)
{
  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID()) dserror("AssembleG: Node ownership inconsistency!");

    /**************************************************** g-vector ******/
    if (cnode->CoData().Getg() != 0.0)
    {
      double gap = cnode->CoData().Getg();

      // std::cout << "Node ID: " << cnode->Id() << " HasProj: " << cnode->HasProj()
      //      << " IsActive: " << cnode->Active() << " Gap: " << gap << std::endl;

      // check if this inactive node has a feasible projection
      // else, it cannot be in contact and weighted gap should be positive
      // (otherwise wrong results possible for g~ because of non-positivity
      // of dual shape functions!!!)
      //******************************************************************
      // TODO: This is only necessary for quadratic LM shape functions!
      // By the way, it makes the method slightly inconsistent
      // (e.g. patch tests with slave side being wider than master side).
      // However, we are able to solve many problems with this little trick.
      // But not all problems, e.g. dropping edge problems would still fail!!!
      // To solve this dilemma, we need a clever modification of the LM shape
      // functions such that they have positive integral values on the
      // "projecting" element part. Once we have this, the following trick
      // can (and should) also be removed in order to make the method
      // consistent again! (08/2013)
      //******************************************************************

      if (!nurbs_)  // only for Lagrange elements
      {
        bool node_has_quad_element = false;
        for (int i = 0; i < cnode->NumElement(); i++)
        {
          if (dynamic_cast<MORTAR::MortarElement*>(cnode->Elements()[i])->IsQuad() == true)
          {
            node_has_quad_element = true;
            break;
          }
        }
        if (!cnode->HasProj() && !cnode->Active() && node_has_quad_element)
        {
          gap = 1.0e12;
          cnode->CoData().Getg() = gap;
        }
      }

      if (constr_direction_ == INPAR::CONTACT::constr_xyz)
      {
        CORE::LINALG::SerialDenseVector gnode(Dim());
        std::vector<int> lm(Dim());
        std::vector<int> lmowner(Dim());
        for (int i = 0; i < Dim(); i++)
        {
          gnode(i) = gap * cnode->MoData().n()[i];
          lm[i] = cnode->Dofs()[i];
          lmowner[i] = cnode->Owner();
        }
        CORE::LINALG::Assemble(gglobal, gnode, lm, lmowner);
      }
      else
      {
        static CORE::LINALG::SerialDenseVector gnode(1);
        static std::vector<int> lm(1);
        static std::vector<int> lmowner(1);

        gnode(0) = gap;
        lm[0] = cnode->Id();
        lmowner[0] = cnode->Owner();
        CORE::LINALG::Assemble(gglobal, gnode, lm, lmowner);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble inactive right hand side                    hiermeier 08/13|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleInactiverhs(Epetra_Vector& inactiverhs)
{
  // FIXME It's possible to improve the performance, if only recently active nodes of the inactive
  // node set, i.e. nodes, which were active in the last iteration, are considered. Since you know,
  // that the lagrange multipliers of former inactive nodes are still equal zero.

  Teuchos::RCP<Epetra_Map> inactivenodes = CORE::LINALG::SplitMap(*snoderowmap_, *activenodes_);
  Teuchos::RCP<Epetra_Map> inactivedofs = CORE::LINALG::SplitMap(*sdofrowmap_, *activedofs_);

  static std::vector<int> lm_gid(Dim());
  static std::vector<int> lm_owner(Dim());
  static CORE::LINALG::SerialDenseVector lm_i(Dim());

  for (int i = 0; i < inactivenodes->NumMyElements(); ++i)
  {
    int gid = inactivenodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("AssembleInactiverhs: Node ownership inconsistency!");
    if (Dim() == 2)
    {
      // calculate the tangential rhs
      for (int j = 0; j < Dim(); ++j)
      {
        lm_owner[j] = cnode->Owner();
        lm_i[j] = -cnode->MoData().lm()[j];  // already negative rhs!!!
      }
      lm_gid[0] = inactivedofs->GID(2 * i);
      lm_gid[1] = inactivedofs->GID(2 * i + 1);

      CORE::LINALG::Assemble(inactiverhs, lm_i, lm_gid, lm_owner);
    }
    else if (Dim() == 3)
    {
      // calculate the tangential rhs
      for (int j = 0; j < Dim(); ++j)
      {
        lm_owner[j] = cnode->Owner();
        lm_i[j] = -cnode->MoData().lm()[j];  // already negative rhs!!!
      }
      lm_gid[0] = inactivedofs->GID(3 * i);
      lm_gid[1] = inactivedofs->GID(3 * i + 1);
      lm_gid[2] = inactivedofs->GID(3 * i + 2);

      CORE::LINALG::Assemble(inactiverhs, lm_i, lm_gid, lm_owner);
    }
  }
}

/*----------------------------------------------------------------------*
 |  Assemble tangential right-hand side                  hiermeier 08/13|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleTangrhs(Epetra_Vector& tangrhs)
{
  static std::vector<int> lm_gid(Dim() - 1);
  static std::vector<int> lm_owner(Dim() - 1);
  static CORE::LINALG::SerialDenseVector lm_t(Dim() - 1);

  for (int i = 0; i < activenodes_->NumMyElements(); ++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID()) dserror("AssembleTangrhs: Node ownership inconsistency!");
    if (constr_direction_ == INPAR::CONTACT::constr_xyz)
    {
      if (Dim() == 2)
      {
        std::vector<int> lm_gid(2);
        std::vector<int> lm_owner(2);
        CORE::LINALG::SerialDenseVector lm_t(2);
        for (int i = 0; i < Dim(); ++i)
        {
          lm_gid[i] = cnode->Dofs()[i];
          lm_owner[i] = cnode->Owner();
          lm_t[i] = 0.;
          for (int j = 0; j < Dim(); ++j)
          {
            lm_t[i] -=
                cnode->CoData().txi()[i] * cnode->CoData().txi()[j] * cnode->MoData().lm()[j];
          }
        }
        CORE::LINALG::Assemble(tangrhs, lm_t, lm_gid, lm_owner);
      }
      else if (Dim() == 3)
      {
        std::vector<int> lm_gid(3);
        std::vector<int> lm_owner(3);
        CORE::LINALG::SerialDenseVector lm_t(3);

        for (int i = 0; i < Dim(); ++i)
        {
          lm_gid[i] = cnode->Dofs()[i];
          lm_owner[i] = cnode->Owner();
          lm_t[i] = 0.;
          for (int j = 0; j < Dim(); ++j)
          {
            lm_t[i] -=
                cnode->CoData().txi()[i] * cnode->CoData().txi()[j] * cnode->MoData().lm()[j];
            lm_t[i] -=
                cnode->CoData().teta()[i] * cnode->CoData().teta()[j] * cnode->MoData().lm()[j];
          }
        }
        CORE::LINALG::Assemble(tangrhs, lm_t, lm_gid, lm_owner);
      }
    }
    else
    {
      if (Dim() == 2)
      {
        lm_gid[0] = activet_->GID(i);
        lm_owner[0] = cnode->Owner();

        lm_t[0] = 0.0;
        for (int j = 0; j < Dim(); ++j)
          lm_t[0] -= cnode->CoData().txi()[j] * cnode->MoData().lm()[j];  // already negative rhs!!!

        CORE::LINALG::Assemble(tangrhs, lm_t, lm_gid, lm_owner);
      }
      else if (Dim() == 3)
      {
        lm_gid[0] = activet_->GID(2 * i);      // even
        lm_gid[1] = activet_->GID(2 * i + 1);  // odd
        lm_owner[0] = cnode->Owner();
        lm_owner[1] = cnode->Owner();

        // calculate the tangential rhs
        lm_t[0] = 0.0;
        lm_t[1] = 0.0;
        for (int j = 0; j < Dim(); ++j)
        {
          lm_t[0] -= cnode->CoData().txi()[j] * cnode->MoData().lm()[j];  // already negative rhs!!!
          lm_t[1] -=
              cnode->CoData().teta()[j] * cnode->MoData().lm()[j];  // already negative rhs!!!
        }
        CORE::LINALG::Assemble(tangrhs, lm_t, lm_gid, lm_owner);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 |  Assemble matrix LinStick with tangential+D+M derivatives  mgit 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinStick(CORE::LINALG::SparseMatrix& linstickLMglobal,
    CORE::LINALG::SparseMatrix& linstickDISglobal, Epetra_Vector& linstickRHSglobal)
{
  // FIXGIT: Assemble LinStick is containing a matrix for the de-
  // rivatives of the Lagrange multipliers. This is according to Hueeber.
  // Because of worse convergence, this is not implemented, but the
  // code is commented after the algorithm.

  // create map of stick nodes
  Teuchos::RCP<Epetra_Map> sticknodes = CORE::LINALG::SplitMap(*activenodes_, *slipnodes_);
  Teuchos::RCP<Epetra_Map> stickt = CORE::LINALG::SplitMap(*activet_, *slipt_);

  // nothing to do if no stick nodes
  if (sticknodes->NumMyElements() == 0) return;

  // information from interface contact parameter list
  INPAR::CONTACT::FrictionType ftype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(InterfaceParams(), "FRICTION");
  double frcoeff_in =
      InterfaceParams().get<double>("FRCOEFF");  // the friction coefficient from the input
  bool gp_slip = DRT::INPUT::IntegralValue<int>(InterfaceParams(), "GP_SLIP_INCR");
  bool frilessfirst = DRT::INPUT::IntegralValue<int>(InterfaceParams(), "FRLESS_FIRST");

  double frcoeff = 0.;  // the friction coefficient actually used
  bool consistent = false;

#if defined(CONSISTENTSTICK) && defined(CONSISTENTSLIP)
  dserror(
      "It's not reasonable to activate both, the consistent stick and slip branch, "
      "because both together will lead again to an inconsistent formulation!");
#endif
#ifdef CONSISTENTSTICK
  consistent = true;
#endif

  if (consistent && DRT::INPUT::IntegralValue<int>(InterfaceParams(), "REGULARIZED_NORMAL_CONTACT"))
    dserror("no consistent stick for regularized contact");

  // loop over all stick nodes of the interface
  for (int i = 0; i < sticknodes->NumMyElements(); ++i)
  {
    int gid = sticknodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("AssembleLinStick: Node ownership inconsistency!");

    // get friction coefficient for this node
    // in case of TSI, the nodal temperature influences the local friction coefficient
    frcoeff = cnode->FrCoeff(frcoeff_in);

    double cn_input = GetCnRef()[GetCnRef().Map().LID(cnode->Id())];
    double ct_input = GetCtRef()[GetCtRef().Map().LID(cnode->Id())];

    // more information from node
    double* n = cnode->MoData().n();
    double* txi = cnode->CoData().txi();
    double* teta = cnode->CoData().teta();
    double* z = cnode->MoData().lm();

    // evaluation of specific components of entries to assemble
    double znor = 0.0;
    double ztxi = 0.0;
    double zteta = 0.0;

    for (int j = 0; j < Dim(); j++)
    {
      znor += n[j] * z[j];
      ztxi += txi[j] * z[j];
      zteta += teta[j] * z[j];
    }

    // prepare assembly, get information from node
    std::vector<CORE::GEN::pairedvector<int, double>> dnmap = cnode->CoData().GetDerivN();
    std::vector<CORE::GEN::pairedvector<int, double>> dtximap = cnode->CoData().GetDerivTxi();
    std::vector<CORE::GEN::pairedvector<int, double>> dtetamap = cnode->CoData().GetDerivTeta();
    std::map<int, double> dscmap;

    // iterator for maps
    std::map<int, double>::iterator colcurr;
    CORE::GEN::pairedvector<int, double>::iterator _colcurr;

    // row number of entries
    std::vector<int> row(Dim() - 1);
    if (Dim() == 2)
    {
      row[0] = stickt->GID(i);
    }
    else if (Dim() == 3)
    {
      row[0] = stickt->GID(2 * i);
      row[1] = stickt->GID(2 * i) + 1;
    }
    else
      dserror("AssemblelinStick: Dimension not correct");

    //****************************************************************
    // CONSISTENT TREATMENT OF CASE FRCOEFF=0 (FRICTIONLESS)
    //****************************************************************
    // popp 08/2012
    //
    // There is a problem with our frictional nonlinear complementarity
    // function when applied to the limit case frcoeff=0 (frictionless).
    // In this case, the simple frictionless sliding condition should
    // be consistently recovered, which unfortunately is not the case.
    // This fact is well-known (see PhD thesis S. Hueeber) and now
    // taken care of by a special treatment as can be seen below
    //
    // seitz 11/2013
    // In case of a consistent treatment of the stick condition (i.e. using
    // the full linearization as in Diss. Gitterle eq.4.21) the case of
    // a vanishing frictional bound (frcoeff*znor==0) needs to be treated
    // as frictionless contact as well.
    // The FRLESS_FIRST option might end up here, depending on if nodes newly
    // in contact are initialized as stick or slip
    //
    //****************************************************************
    if (frcoeff * znor == 0.0 || (cnode->CoData().ActiveOld() == false && frilessfirst))
    {
      //**************************************************************
      // calculation of matrix entries of linearized slip condition
      //**************************************************************

      // 1) Entries from differentiation with respect to LM
      /******************************************************************/

      // loop over the dimension
      for (int dim = 0; dim < cnode->NumDof(); ++dim)
      {
        int col = cnode->Dofs()[dim];
        double valtxi = txi[dim];

        double valteta = 0.0;
        if (Dim() == 3) valteta = teta[dim];

        if (constr_direction_ == INPAR::CONTACT::constr_xyz)
        {
          for (int j = 0; j < Dim(); j++)
          {
            if (abs(valtxi * txi[j]) > 1.e-12)
              linstickLMglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
            if (Dim() == 3)
              if (abs(valteta * teta[j]) > 1.e-12)
                linstickLMglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
          }
        }
        else
        {
          if (abs(valtxi) > 1.0e-12) linstickLMglobal.Assemble(valtxi, row[0], col);
          if (Dim() == 3)
            if (abs(valteta) > 1.0e-12) linstickLMglobal.Assemble(valteta, row[1], col);
        }
      }

      // 2) Entries on right hand side
      /******************************************************************/
      if (constr_direction_ == INPAR::CONTACT::constr_xyz)
      {
        CORE::LINALG::SerialDenseVector rhsnode(Dim());
        std::vector<int> lm(Dim());
        std::vector<int> lmowner(Dim());
        for (int j = 0; j < Dim(); j++)
        {
          lm[j] = cnode->Dofs()[j];
          lmowner[j] = cnode->Owner();
          rhsnode[j] -= ztxi * txi[j];
          if (Dim() == 3) rhsnode[j] -= zteta * teta[j];
        }
        CORE::LINALG::Assemble(linstickRHSglobal, rhsnode, lm, lmowner);
      }
      else
      {
        CORE::LINALG::SerialDenseVector rhsnode(Dim() - 1);
        std::vector<int> lm(Dim() - 1);
        std::vector<int> lmowner(Dim() - 1);

        lm[0] = cnode->Dofs()[1];
        lmowner[0] = cnode->Owner();

        rhsnode[0] = -ztxi;  // already negative rhs!!!

        if (Dim() == 3)
        {
          rhsnode[1] = -zteta;  // already negative rhs!!!

          lm[1] = cnode->Dofs()[2];
          lmowner[1] = cnode->Owner();
        }
        CORE::LINALG::Assemble(linstickRHSglobal, rhsnode, lm, lmowner);
      }

      // 3) Entries from differentiation with respect to displacements
      /******************************************************************/

      // loop over dimensions
      for (int j = 0; j < Dim(); ++j)
      {
        // loop over all entries of the current derivative map (txi)
        for (_colcurr = dtximap[j].begin(); _colcurr != dtximap[j].end(); ++_colcurr)
        {
          int col = _colcurr->first;
          double val = (_colcurr->second) * z[j];

          // do not assemble zeros into matrix
          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int j = 0; j < Dim(); j++)
              if (abs(val * txi[j]) > 1.e-12)
                linstickDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
          }
          else if (abs(val) > 1.0e-12)
            linstickDISglobal.Assemble(val, row[0], col);
        }

        if (Dim() == 3)
        {
          // loop over all entries of the current derivative map (teta)
          for (_colcurr = dtetamap[j].begin(); _colcurr != dtetamap[j].end(); ++_colcurr)
          {
            int col = _colcurr->first;
            double val = (_colcurr->second) * z[j];

            // do not assemble zeros into s matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(val * teta[j]) > 1.e-12)
                  linstickDISglobal.Assemble(val * teta[j], cnode->Dofs()[j], col);
            }
            else if (abs(val) > 1.0e-12)
              linstickDISglobal.Assemble(val, row[1], col);
          }
        }
      }
    }
    else if (consistent && ftype == INPAR::CONTACT::friction_coulomb)
    {
      std::map<int, double>& dgmap = cnode->CoData().GetDerivG();

      // check for Dimension of derivative maps
      for (int j = 0; j < Dim() - 1; ++j)
        if ((int)dnmap[j].size() != (int)dnmap[j + 1].size())
          dserror("AssembleLinStick: Column dim. of nodal DerivN-map is inconsistent!");

      for (int j = 0; j < Dim() - 1; ++j)
        if ((int)dtximap[j].size() != (int)dtximap[j + 1].size())
          dserror("AssembleLinStick: Column dim. of nodal DerivTxi-map is inconsistent!");

      if (Dim() == 3)
      {
        for (int j = 0; j < Dim() - 1; ++j)
          if ((int)dtximap[j].size() != (int)dtximap[j + 1].size())
            dserror("AssembleLinStick: Column dim. of nodal DerivTeta-map is inconsistent!");
      }

      double cn = cn_input;
      double ct = ct_input;

      double& wgap = cnode->CoData().Getg();

      // evaluation of specific components of entries to assemble
      double jumptxi = 0;
      double jumpteta = 0;

      double* jump = cnode->FriData().jump();

      // choose slip increment
      if (gp_slip)
      {
        jumptxi = cnode->FriData().jump_var()[0];

        if (Dim() == 3) jumpteta = cnode->FriData().jump_var()[1];
      }
      else
      {
        // more information from node
        for (int i = 0; i < Dim(); i++)
        {
          jumptxi += txi[i] * jump[i];
          jumpteta += teta[i] * jump[i];
        }
      }

      // check for dimensions
      if (Dim() == 2 and (jumpteta != 0.0))
        dserror("AssembleLinStick: jumpteta must be zero in 2D");

      //**************************************************************
      // calculation of matrix entries of linearized stick condition
      //**************************************************************

      // 1) Entries from differentiation with respect to LM
      /******************************************************************/

      // loop over the dimension
      for (int dim = 0; dim < cnode->NumDof(); ++dim)
      {
        double valtxi = 0.0;
        double valteta = 0.0;
        int col = cnode->Dofs()[dim];
        valtxi = -frcoeff * ct * jumptxi * n[dim];

        if (Dim() == 3)
        {
          valteta = -frcoeff * ct * jumpteta * n[dim];
        }

        // do not assemble zeros into matrix
        if (constr_direction_ == INPAR::CONTACT::constr_xyz)
        {
          for (int j = 0; j < Dim(); j++)
          {
            if (abs(valtxi * txi[j]) > 1.e-12)
              linstickLMglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
            if (Dim() == 3)
              if (abs(valteta * teta[j]) > 1.e-12)
                linstickLMglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
          }
        }
        else
        {
          if (abs(valtxi) > 1.0e-12) linstickLMglobal.Assemble(valtxi, row[0], col);
          if (Dim() == 3)
            if (abs(valteta) > 1.0e-12) linstickLMglobal.Assemble(valteta, row[1], col);
        }
      }

      // Entries on right hand side ****************************
      if (constr_direction_ == INPAR::CONTACT::constr_xyz)
      {
        CORE::LINALG::SerialDenseVector rhsnode(Dim());
        std::vector<int> lm(Dim());
        std::vector<int> lmowner(Dim());

        for (int j = 0; j < Dim(); j++)
        {
          lm[j] = cnode->Dofs()[j];
          lmowner[j] = cnode->Owner();
          rhsnode(j) += frcoeff * (znor - cn * wgap) * ct * jumptxi * txi[j];
          if (Dim() == 3) rhsnode(j) += frcoeff * (znor - cn * wgap) * ct * jumpteta * teta[j];
        }
        CORE::LINALG::Assemble(linstickRHSglobal, rhsnode, lm, lmowner);
      }
      else
      {
        CORE::LINALG::SerialDenseVector rhsnode(Dim() - 1);
        std::vector<int> lm(Dim() - 1);
        std::vector<int> lmowner(Dim() - 1);
        rhsnode(0) = frcoeff * (znor - cn * wgap) * ct * jumptxi;

        lm[0] = cnode->Dofs()[1];
        lmowner[0] = cnode->Owner();
        if (Dim() == 3)
        {
          rhsnode(1) = frcoeff * (znor - cn * wgap) * ct * jumpteta;
          lm[1] = cnode->Dofs()[2];
          lmowner[1] = cnode->Owner();
        }
        CORE::LINALG::Assemble(linstickRHSglobal, rhsnode, lm, lmowner);
      }


      // 3) Entries from differentiation with respect to displacements
      /******************************************************************/

      if (gp_slip)
      {
        std::vector<std::map<int, double>> derivjump_ = cnode->FriData().GetDerivVarJump();

        // txi
        for (colcurr = derivjump_[0].begin(); colcurr != derivjump_[0].end(); ++colcurr)
        {
          int col = colcurr->first;
          double valtxi = -frcoeff * (znor - cn * wgap) * ct * (colcurr->second);

          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int j = 0; j < Dim(); j++)
              if (abs(valtxi * txi[j]) > 1.e-12)
                linstickDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
          }
          else if (abs(valtxi) > 1.0e-12)
            linstickDISglobal.Assemble(valtxi, row[0], col);
        }
        // teta
        for (colcurr = derivjump_[1].begin(); colcurr != derivjump_[1].end(); ++colcurr)
        {
          int col = colcurr->first;
          double valteta = -frcoeff * (znor - cn * wgap) * ct * (colcurr->second);

          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int j = 0; j < Dim(); j++)
              if (abs(valteta * teta[j]) > 1.e-12)
                linstickDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
          }
          else if (abs(valteta) > 1.0e-12)
            linstickDISglobal.Assemble(valteta, row[1], col);
        }

        // ... old slip
        // get linearization of jump vector
        std::vector<std::map<int, double>> derivjump = cnode->FriData().GetDerivJump();

        // loop over dimensions
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          // linearization of normal direction *****************************************
          // loop over all entries of the current derivative map
          for (_colcurr = dnmap[dim].begin(); _colcurr != dnmap[dim].end(); ++_colcurr)
          {
            int col = _colcurr->first;
            double valtxi = 0.0;
            valtxi = -frcoeff * z[dim] * _colcurr->second * ct * jumptxi;

            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(valtxi * txi[j]) > 1.e-12)
                  linstickDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(valtxi) > 1.0e-12)
              linstickDISglobal.Assemble(valtxi, row[0], col);

            if (Dim() == 3)
            {
              double valteta = 0.0;
              valteta = -frcoeff * z[dim] * _colcurr->second * ct * jumpteta;

              // do not assemble zeros into matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(valteta * teta[j]) > 1.e-12)
                    linstickDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(valteta) > 1.0e-12)
                linstickDISglobal.Assemble(valteta, row[1], col);
            }
          }
        }  // loop over all dimensions
      }
      else  // std slip
      {
        // ... old slip
        // get linearization of jump vector
        std::vector<std::map<int, double>> derivjump = cnode->FriData().GetDerivJump();

        // loop over dimensions
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          // loop over all entries of the current derivative map (jump)
          for (colcurr = derivjump[dim].begin(); colcurr != derivjump[dim].end(); ++colcurr)
          {
            int col = colcurr->first;

            double valtxi = 0.0;
            valtxi = -frcoeff * (znor - cn * wgap) * ct * txi[dim] * (colcurr->second);

            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(valtxi * txi[j]) > 1.e-12)
                  linstickDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(valtxi) > 1.0e-12)
              linstickDISglobal.Assemble(valtxi, row[0], col);

            if (Dim() == 3)
            {
              double valteta = 0.0;
              valteta = -frcoeff * (znor - cn * wgap) * ct * teta[dim] * (colcurr->second);

              // do not assemble zeros into matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(valteta * teta[j]) > 1.e-12)
                    linstickDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(valteta) > 1.0e-12)
                linstickDISglobal.Assemble(valteta, row[1], col);
            }
          }

          // linearization first tangential direction *********************************
          // loop over all entries of the current derivative map (txi)
          for (_colcurr = dtximap[dim].begin(); _colcurr != dtximap[dim].end(); ++_colcurr)
          {
            int col = _colcurr->first;
            double valtxi = 0.0;
            valtxi = -frcoeff * (znor - cn * wgap) * ct * jump[dim] * _colcurr->second;

            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(valtxi * txi[j]) > 1.e-12)
                  linstickDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(valtxi) > 1.0e-12)
              linstickDISglobal.Assemble(valtxi, row[0], col);
          }
          // linearization second tangential direction *********************************
          if (Dim() == 3)
          {
            // loop over all entries of the current derivative map (teta)
            for (_colcurr = dtetamap[dim].begin(); _colcurr != dtetamap[dim].end(); ++_colcurr)
            {
              int col = _colcurr->first;
              double valteta = 0.0;
              valteta = -frcoeff * (znor - cn * wgap) * ct * jump[dim] * _colcurr->second;

              // do not assemble zeros into matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(valteta * teta[j]) > 1.e-12)
                    linstickDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(valteta) > 1.0e-12)
                linstickDISglobal.Assemble(valteta, row[1], col);
            }
          }

          // linearization of normal direction *****************************************
          // loop over all entries of the current derivative map
          for (_colcurr = dnmap[dim].begin(); _colcurr != dnmap[dim].end(); ++_colcurr)
          {
            int col = _colcurr->first;
            double valtxi = 0.0;
            valtxi = -frcoeff * z[dim] * _colcurr->second * ct * jumptxi;

            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(valtxi * txi[j]) > 1.e-12)
                  linstickDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(valtxi) > 1.0e-12)
              linstickDISglobal.Assemble(valtxi, row[0], col);

            if (Dim() == 3)
            {
              double valteta = 0.0;
              valteta = -frcoeff * z[dim] * _colcurr->second * ct * jumpteta;

              // do not assemble zeros into matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(valteta * teta[j]) > 1.e-12)
                    linstickDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(valteta) > 1.0e-12)
                linstickDISglobal.Assemble(valteta, row[1], col);
            }
          }
        }  // loop over all dimensions
      }

      // linearization of weighted gap**********************************************
      // loop over all entries of the current derivative map fixme
      for (colcurr = dgmap.begin(); colcurr != dgmap.end(); ++colcurr)
      {
        int col = colcurr->first;
        double valtxi = 0.0;
        valtxi = frcoeff * colcurr->second * ct * cn * jumptxi;

        // do not assemble zeros into matrix
        if (constr_direction_ == INPAR::CONTACT::constr_xyz)
        {
          for (int j = 0; j < Dim(); j++)
            if (abs(valtxi * txi[j]) > 1.e-12)
              linstickDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
        }
        else if (abs(valtxi) > 1.0e-12)
          linstickDISglobal.Assemble(valtxi, row[0], col);

        if (Dim() == 3)
        {
          double valteta = 0.0;
          valteta = frcoeff * colcurr->second * ct * cn * jumpteta;

          // do not assemble zeros into matrix
          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int j = 0; j < Dim(); j++)
              if (abs(valteta * teta[j]) > 1.e-12)
                linstickDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
          }
          else if (abs(valteta) > 1.0e-12)
            linstickDISglobal.Assemble(valteta, row[1], col);
        }
      }
    }
    else  // not consistent stick
    {
      for (int j = 0; j < Dim() - 1; ++j)
        if ((int)dtximap[j].size() != (int)dtximap[j + 1].size())
          dserror("AssembleLinStick: Column dim. of nodal DerivTxi-map is inconsistent!");

      if (Dim() == 3)
      {
        for (int j = 0; j < Dim() - 1; ++j)
          if ((int)dtximap[j].size() != (int)dtximap[j + 1].size())
            dserror("AssembleLinStick: Column dim. of nodal DerivTeta-map is inconsistent!");
      }

      // evaluation of specific components of entries to assemble
      double jumptxi = 0;
      double jumpteta = 0;

      // more information from node
      double* jump = cnode->FriData().jump();

      if (gp_slip)
      {
        jumptxi = cnode->FriData().jump_var()[0];

        if (Dim() == 3) jumpteta = cnode->FriData().jump_var()[1];
      }
      else
      {
        for (int i = 0; i < Dim(); i++)
        {
          jumptxi += txi[i] * jump[i];
          jumpteta += teta[i] * jump[i];
        }
      }

      // check for dimensions
      if (Dim() == 2 and (jumpteta != 0.0))
        dserror("AssembleLinStick: jumpteta must be zero in 2D");

      // Entries on right hand side
      /************************************************ (-utxi, -uteta) ***/
      if (constr_direction_ == INPAR::CONTACT::constr_xyz)
      {
        CORE::LINALG::SerialDenseVector rhsnode(Dim());
        std::vector<int> lm(Dim());
        std::vector<int> lmowner(Dim());
        for (int j = 0; j < Dim(); j++)
        {
          lm[j] = cnode->Dofs()[j];
          lmowner[j] = cnode->Owner();
          rhsnode(j) -= jumptxi * txi[j];
          if (Dim() == 3) rhsnode(j) -= jumpteta * teta[j];
        }
        CORE::LINALG::Assemble(linstickRHSglobal, rhsnode, lm, lmowner);
      }
      else
      {
        CORE::LINALG::SerialDenseVector rhsnode(Dim() - 1);
        std::vector<int> lm(Dim() - 1);
        std::vector<int> lmowner(Dim() - 1);

        // modification to stabilize the convergence of the lagrange multiplier incr (hiermeier
        // 08/13)
        if (abs(jumptxi) < 1e-15)
          rhsnode(0) = 0.0;
        else
          rhsnode(0) = -jumptxi;

        lm[0] = cnode->Dofs()[1];
        lmowner[0] = cnode->Owner();

        if (Dim() == 3)
        {
          // modification to stabilize the convergence of the lagrange multiplier incr (hiermeier
          // 08/13)
          if (abs(jumpteta) < 1e-15)
            rhsnode(1) = 0.0;
          else
            rhsnode(1) = -jumpteta;

          lm[1] = cnode->Dofs()[2];
          lmowner[1] = cnode->Owner();
        }
        CORE::LINALG::Assemble(linstickRHSglobal, rhsnode, lm, lmowner);
      }


      // The routine "ApplyDirichlet" in SaddlepointSolve can only set ones on the diagonal
      // if there has already been a diagonal entry in the sparse matrix
      for (int j = 0; j < Dim(); j++)
      {
        if (cnode->DbcDofs()[j] == true)
        {
          linstickLMglobal.Assemble(1.e-12, cnode->Dofs()[j], cnode->Dofs()[j]);
        }
      }

      // Entries from differentiation with respect to displacements
      /*** 1 ************************************** tangent.deriv(jump) ***/
      if (gp_slip)
      {
        std::map<int, double> derivjump1 = cnode->FriData().GetDerivVarJump()[0];
        std::map<int, double> derivjump2 = cnode->FriData().GetDerivVarJump()[1];

        for (colcurr = derivjump1.begin(); colcurr != derivjump1.end(); ++colcurr)
        {
          int col = colcurr->first;
          double valtxi = colcurr->second;

          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int j = 0; j < Dim(); j++)
              if (abs(valtxi * txi[j]) > 1.e-12)
                linstickDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
          }
          else if (abs(valtxi) > 1.0e-12)
            linstickDISglobal.Assemble(valtxi, row[0], col);
        }

        if (Dim() == 3)
        {
          for (colcurr = derivjump2.begin(); colcurr != derivjump2.end(); ++colcurr)
          {
            int col = colcurr->first;
            double valteta = colcurr->second;

            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(valteta * teta[j]) > 1.e-12)
                  linstickDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
            }
            else if (abs(valteta) > 1.0e-12)
              linstickDISglobal.Assemble(valteta, row[1], col);
          }
        }
      }
      else  // std. slip
      {
        // get linearization of jump vector
        std::vector<std::map<int, double>> derivjump = cnode->FriData().GetDerivJump();

        if (derivjump.size() < 1) dserror("AssembleLinStick: Derivative of jump is not exiting!");

        // loop over dimensions
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          // loop over all entries of the current derivative map (jump)
          for (colcurr = derivjump[dim].begin(); colcurr != derivjump[dim].end(); ++colcurr)
          {
            int col = colcurr->first;
            double valtxi = txi[dim] * colcurr->second;

            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(valtxi * txi[j]) > 1.e-12)
                  linstickDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(valtxi) > 1.0e-12)
              linstickDISglobal.Assemble(valtxi, row[0], col);

            if (Dim() == 3)
            {
              double valteta = teta[dim] * colcurr->second;

              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(valteta * teta[j]) > 1.e-12)
                    linstickDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(valteta) > 1.0e-12)
                linstickDISglobal.Assemble(valteta, row[1], col);
            }
          }
        }

        /*** 2 ************************************** deriv(tangent).jump ***/
        // loop over dimensions
        for (int j = 0; j < Dim(); ++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (_colcurr = dtximap[j].begin(); _colcurr != dtximap[j].end(); ++_colcurr)
          {
            int col = _colcurr->first;
            double val = jump[j] * _colcurr->second;

            // do not assemble zeros into s matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(val * txi[j]) > 1.e-12)
                  linstickDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(val) > 1.0e-12)
              linstickDISglobal.Assemble(val, row[0], col);
          }

          if (Dim() == 3)
          {
            // loop over all entries of the current derivative map (teta)
            for (_colcurr = dtetamap[j].begin(); _colcurr != dtetamap[j].end(); ++_colcurr)
            {
              int col = _colcurr->first;
              double val = jump[j] * _colcurr->second;

              // do not assemble zeros into matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(val * teta[j]) > 1.e-12)
                    linstickDISglobal.Assemble(val * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(val) > 1.0e-12)
                linstickDISglobal.Assemble(val, row[1], col);
            }
          }
        }
      }
    }
  }
  return;
}

/*---------------------------------------------------------------------*
 | Assemble matrix LinSlip with tangential+D+M derivatives  mgit 02/09 |
 *---------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinSlip(CORE::LINALG::SparseMatrix& linslipLMglobal,
    CORE::LINALG::SparseMatrix& linslipDISglobal, Epetra_Vector& linslipRHSglobal)
{
  // nothing to do if no slip nodes
  if (slipnodes_->NumMyElements() == 0) return;

  // information from interface contact parameter list
  INPAR::CONTACT::FrictionType ftype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(InterfaceParams(), "FRICTION");
  double frbound = InterfaceParams().get<double>("FRBOUND");
  double frcoeff_in =
      InterfaceParams().get<double>("FRCOEFF");  // the friction coefficient from the input
  bool gp_slip = DRT::INPUT::IntegralValue<int>(InterfaceParams(), "GP_SLIP_INCR");
  bool frilessfirst = DRT::INPUT::IntegralValue<int>(InterfaceParams(), "FRLESS_FIRST");

  // the friction coefficient adapted by every node (eg depending on the local temperature)
  double frcoeff = 0.;

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
  // Coulomb Friction
  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
  if (ftype == INPAR::CONTACT::friction_coulomb)
  {
    // loop over all slip nodes of the interface
    for (int i = 0; i < slipnodes_->NumMyElements(); ++i)
    {
      int gid = slipnodes_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("Cannot find node with gid %", gid);
      FriNode* cnode = dynamic_cast<FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("AssembleLinSlip: Node ownership inconsistency!");

      // get friction coefficient for this node
      // in case of TSI, the nodal temperature influences the local friction coefficient
      frcoeff = cnode->FrCoeff(frcoeff_in);

      double cn_input = GetCnRef()[GetCnRef().Map().LID(cnode->Id())];
      double ct_input = GetCtRef()[GetCtRef().Map().LID(cnode->Id())];

      // prepare assembly, get information from node
      std::vector<CORE::GEN::pairedvector<int, double>> dnmap = cnode->CoData().GetDerivN();
      std::vector<CORE::GEN::pairedvector<int, double>> dtximap = cnode->CoData().GetDerivTxi();
      std::vector<CORE::GEN::pairedvector<int, double>> dtetamap = cnode->CoData().GetDerivTeta();

      double cn = cn_input;
      double ct = ct_input;

      // check for Dimension of derivative maps
      for (int j = 0; j < Dim() - 1; ++j)
        if ((int)dnmap[j].size() != (int)dnmap[j + 1].size())
          dserror("AssembleLinSlip: Column dim. of nodal DerivTxi-map is inconsistent!");

      for (int j = 0; j < Dim() - 1; ++j)
        if ((int)dtximap[j].size() != (int)dtximap[j + 1].size())
          dserror("AssembleLinSlip: Column dim. of nodal DerivTxi-map is inconsistent!");

      if (Dim() == 3)
      {
        for (int j = 0; j < Dim() - 1; ++j)
          if ((int)dtximap[j].size() != (int)dtximap[j + 1].size())
            dserror("AssembleLinSlip: Column dim. of nodal DerivTeta-map is inconsistent!");
      }

      // more information from node
      double* n = cnode->MoData().n();
      double* txi = cnode->CoData().txi();
      double* teta = cnode->CoData().teta();
      double* z = cnode->MoData().lm();
      double wgap = cnode->CoData().Getg();

      // iterator for maps
      CORE::GEN::pairedvector<int, double>::iterator colcurr;
      std::map<int, double>::iterator _colcurr;

      // row number of entries
      std::vector<int> row(Dim() - 1);
      if (Dim() == 2)
      {
        row[0] = slipt_->GID(i);
      }
      else if (Dim() == 3)
      {
        row[0] = slipt_->GID(2 * i);
        row[1] = slipt_->GID(2 * i) + 1;
      }
      else
        dserror("AssemblelinSlip: Dimension not correct");

      // evaluation of specific components of entries to assemble
      double znor = 0.0;
      double ztxi = 0.0;
      double zteta = 0.0;
      double jumptxi = 0.0;
      double jumpteta = 0.0;
      double euclidean = 0.0;
      double* jump = cnode->FriData().jump();

      if (gp_slip)
      {
        jumptxi = cnode->FriData().jump_var()[0];

        if (Dim() == 3) jumpteta = cnode->FriData().jump_var()[1];

        for (int i = 0; i < Dim(); i++)
        {
          znor += n[i] * z[i];
          ztxi += txi[i] * z[i];
          zteta += teta[i] * z[i];
        }
      }
      else  // std. slip
      {
        for (int i = 0; i < Dim(); i++)
        {
          znor += n[i] * z[i];
          ztxi += txi[i] * z[i];
          zteta += teta[i] * z[i];
          jumptxi += txi[i] * jump[i];
          jumpteta += teta[i] * jump[i];
        }
      }

      // evaluate euclidean norm ||vec(zt)+ct*vec(jumpt)||
      std::vector<double> sum1(Dim() - 1, 0);
      sum1[0] = ztxi + ct * jumptxi;
      if (Dim() == 3) sum1[1] = zteta + ct * jumpteta;
      if (Dim() == 2) euclidean = abs(sum1[0]);
      if (Dim() == 3) euclidean = sqrt(sum1[0] * sum1[0] + sum1[1] * sum1[1]);

      // check of dimensions
      if (Dim() == 2 and (zteta != 0.0 or jumpteta != 0.0))
        dserror("AssemblelinSlip: zteta and jumpteta must be zero in 2D");

      //****************************************************************
      // CONSISTENT TREATMENT OF CASE FRCOEFF=0 (FRICTIONLESS)
      //****************************************************************
      // popp 08/2012
      //
      // There is a problem with our frictional nonlinear complementarity
      // function when applied to the limit case frcoeff=0 (frictionless).
      // In this case, the simple frictionless sliding condition should
      // be consistently recovered, which unfortunately is not the case.
      // This fact is well-known (see PhD thesis S. Hueeber) and now
      // taken care of by a special treatment as can be seen below
      //
      // seitz11/13
      // Like that FRLESS_FIRST results in a frictionless constraint.
      // Nodes with a vanishing slip and tangential Lagrange multiplier
      // are treated as frictionless as well.
      //
      //****************************************************************
      if (frcoeff == 0.0 || (cnode->CoData().ActiveOld() == false && frilessfirst) ||
          euclidean == 0.0)
      {
        //**************************************************************
        // calculation of matrix entries of linearized slip condition
        //**************************************************************

        // 1) Entries from differentiation with respect to LM
        /******************************************************************/

        // loop over the dimension
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          int col = cnode->Dofs()[dim];
          double valtxi = txi[dim];

          double valteta = 0.0;
          if (Dim() == 3) valteta = teta[dim];

          // do not assemble zeros into matrix
          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int j = 0; j < Dim(); j++)
            {
              if (abs(valtxi * txi[j]) > 1.e-12)
                linslipLMglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
              if (Dim() == 3)
                if (abs(valteta * teta[j]) > 1.e-12)
                  linslipLMglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
            }
          }
          else
          {
            if (abs(valtxi) > 1.0e-12) linslipLMglobal.Assemble(valtxi, row[0], col);
            if (Dim() == 3)
              if (abs(valteta) > 1.0e-12) linslipLMglobal.Assemble(valteta, row[1], col);
          }
        }

        // 2) Entries on right hand side
        /******************************************************************/
        if (constr_direction_ == INPAR::CONTACT::constr_xyz)
        {
          CORE::LINALG::SerialDenseVector rhsnode(Dim());
          std::vector<int> lm(Dim());
          std::vector<int> lmowner(Dim());
          for (int j = 0; j < Dim(); j++)
          {
            lm[j] = cnode->Dofs()[j];
            lmowner[j] = cnode->Owner();
            rhsnode[j] -= ztxi * txi[j];
            if (Dim() == 3) rhsnode[j] -= zteta * teta[j];
          }
          CORE::LINALG::Assemble(linslipRHSglobal, rhsnode, lm, lmowner);
        }
        else
        {
          CORE::LINALG::SerialDenseVector rhsnode(Dim() - 1);
          std::vector<int> lm(Dim() - 1);
          std::vector<int> lmowner(Dim() - 1);

          lm[0] = cnode->Dofs()[1];
          lmowner[0] = cnode->Owner();

          rhsnode[0] = -ztxi;  // already negative rhs!!!

          if (Dim() == 3)
          {
            rhsnode[1] = -zteta;  // already negative rhs!!!

            lm[1] = cnode->Dofs()[2];
            lmowner[1] = cnode->Owner();
          }
          CORE::LINALG::Assemble(linslipRHSglobal, rhsnode, lm, lmowner);
        }

        // 3) Entries from differentiation with respect to displacements
        /******************************************************************/

        // loop over dimensions
        for (int j = 0; j < Dim(); ++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (colcurr = dtximap[j].begin(); colcurr != dtximap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
            double val = (colcurr->second) * z[j];

            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(val * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(val) > 1.0e-12)
              linslipDISglobal.Assemble(val, row[0], col);
          }

          if (Dim() == 3)
          {
            // loop over all entries of the current derivative map (teta)
            for (colcurr = dtetamap[j].begin(); colcurr != dtetamap[j].end(); ++colcurr)
            {
              int col = colcurr->first;
              double val = (colcurr->second) * z[j];

              // do not assemble zeros into s matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(val * teta[j]) > 1.e-12)
                    linslipDISglobal.Assemble(val * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(val) > 1.0e-12)
                linslipDISglobal.Assemble(val, row[1], col);
            }
          }
        }
      }

      //****************************************************************
      // STANDARD TREATMENT OF CASE FRCOEFF!=0 (FRICTIONAL)
      //****************************************************************
      else
      {
        //**************************************************************
        // calculation of matrix entries of linearized slip condition
        //**************************************************************

        // 1) Entries from differentiation with respect to LM
        /******************************************************************/

        // loop over the dimension
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          double valtxi = 0.0;
          int col = cnode->Dofs()[dim];

          double valtxi0 = euclidean * txi[dim];
          double valtxi1 = ((ztxi + ct * jumptxi) / euclidean * ztxi) * txi[dim];
          double valtxi3 = (zteta + ct * jumpteta) / euclidean * ztxi * teta[dim];

#ifdef CONSISTENTSLIP
          valtxi0 = valtxi0 / (znor - cn * wgap);
          valtxi1 = valtxi1 / (znor - cn * wgap);
          valtxi3 = valtxi3 / (znor - cn * wgap);

          // Additional term
          valtxi0 -= euclidean * ztxi / pow(znor - cn * wgap, 2.0) * n[dim];

          double valtxi2 = -frcoeff * txi[dim];
#else
          double valtxi2 =
              -frcoeff * (znor - cn * wgap) * txi[dim] - frcoeff * (ztxi + ct * jumptxi) * n[dim];
#endif
          valtxi = valtxi0 + valtxi1 + valtxi2 + valtxi3;

          double valteta = 0.0;
          if (Dim() == 3)
          {
            double valteta0 = euclidean * teta[dim];
            double valteta1 = ((ztxi + ct * jumptxi) / euclidean * zteta) * txi[dim];
            double valteta3 = (zteta + ct * jumpteta) / euclidean * zteta * teta[dim];
#ifdef CONSISTENTSLIP
            valteta0 = valteta0 / (znor - cn * wgap);
            valteta1 = valteta1 / (znor - cn * wgap);
            valteta3 = valteta3 / (znor - cn * wgap);

            // Additional term
            valteta0 -= euclidean * zteta / pow(znor - cn * wgap, 2.0) * n[dim];

            double valteta2 = -frcoeff * teta[dim];
#else
            double valteta2 = -frcoeff * (znor - cn * wgap) * teta[dim] -
                              frcoeff * (zteta + ct * jumpteta) * n[dim];
#endif
            valteta = valteta0 + valteta1 + valteta2 + valteta3;
          }

          // do not assemble zeros into matrix
          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int j = 0; j < Dim(); j++)
            {
              if (abs(valtxi * txi[j]) > 1.e-12)
                linslipLMglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
              if (Dim() == 3)
                if (abs(valteta * teta[j]) > 1.e-12)
                  linslipLMglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
            }
          }
          else
          {
            if (abs(valtxi) > 1.0e-12) linslipLMglobal.Assemble(valtxi, row[0], col);
            if (Dim() == 3)
              if (abs(valteta) > 1.0e-12) linslipLMglobal.Assemble(valteta, row[1], col);
          }
        }

        // 2) Entries on right hand side
        /******************************************************************/
        if (constr_direction_ == INPAR::CONTACT::constr_xyz)
        {
          CORE::LINALG::SerialDenseVector rhsnode(Dim());
          std::vector<int> lm(Dim());
          std::vector<int> lmowner(Dim());

#ifdef CONSISTENTSLIP
          double valuetxi1 =
              -(euclidean)*ztxi / (znor - cn * wgap) + frcoeff * (ztxi + ct * jumptxi);
#else
          double valuetxi1 =
              -(euclidean)*ztxi + (frcoeff * (znor - cn * wgap)) * (ztxi + ct * jumptxi);
#endif

          for (int j = 0; j < Dim(); j++)
          {
            lm[j] = cnode->Dofs()[j];
            lmowner[j] = cnode->Owner();
            rhsnode(j) += valuetxi1 * txi[j];
          }

          if (Dim() == 3)
          {
#ifdef CONSISTENTSLIP
            double valueteta1 =
                -(euclidean)*zteta / (znor - cn * wgap) + frcoeff * (zteta + ct * jumpteta);
#else
            double valueteta1 =
                -(euclidean)*zteta + (frcoeff * (znor - cn * wgap)) * (zteta + ct * jumpteta);
#endif

            for (int j = 0; j < Dim(); j++) rhsnode(j) += valueteta1 * teta[j];
          }
          CORE::LINALG::Assemble(linslipRHSglobal, rhsnode, lm, lmowner);
        }
        else
        {
          CORE::LINALG::SerialDenseVector rhsnode(Dim() - 1);
          std::vector<int> lm(Dim() - 1);
          std::vector<int> lmowner(Dim() - 1);
#ifdef CONSISTENTSLIP
          double valuetxi1 =
              -(euclidean)*ztxi / (znor - cn * wgap) + frcoeff * (ztxi + ct * jumptxi);
#else
          double valuetxi1 =
              -(euclidean)*ztxi + (frcoeff * (znor - cn * wgap)) * (ztxi + ct * jumptxi);
#endif
          rhsnode(0) = valuetxi1;
          lm[0] = cnode->Dofs()[1];
          lmowner[0] = cnode->Owner();

          if (Dim() == 3)
          {
#ifdef CONSISTENTSLIP
            double valueteta1 =
                -(euclidean)*zteta / (znor - cn * wgap) + frcoeff * (zteta + ct * jumpteta);
#else
            double valueteta1 =
                -(euclidean)*zteta + (frcoeff * (znor - cn * wgap)) * (zteta + ct * jumpteta);
#endif
            rhsnode(1) = valueteta1;

            lm[1] = cnode->Dofs()[2];
            lmowner[1] = cnode->Owner();
          }
          CORE::LINALG::Assemble(linslipRHSglobal, rhsnode, lm, lmowner);
        }

        // 3) Entries from differentiation with respect to displacements
        /******************************************************************/
        std::map<int, double> derivjump1, derivjump2;  // for gp slip
        std::vector<std::map<int, double>> derivjump;  // for dm slip

        /*** 01  ********* -Deriv(euclidean).ct.tangent.deriv(u)*ztan ***/
        if (gp_slip)
        {
          derivjump1 = cnode->FriData().GetDerivVarJump()[0];
          derivjump2 = cnode->FriData().GetDerivVarJump()[1];

          for (_colcurr = derivjump1.begin(); _colcurr != derivjump1.end(); ++_colcurr)
          {
            int col = _colcurr->first;
            double valtxi1 = (ztxi + ct * jumptxi) / euclidean * ct * _colcurr->second * ztxi;
            double valteta1 = (ztxi + ct * jumptxi) / euclidean * ct * _colcurr->second * zteta;

            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
              {
                if (abs(valtxi1 * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(valtxi1 * txi[j], cnode->Dofs()[j], col);
                if (abs(valteta1 * teta[j]) > 1.e-12)
                  linslipDISglobal.Assemble(valteta1 * teta[j], cnode->Dofs()[j], col);
              }
            }
            else
            {
              if (abs(valtxi1) > 1.0e-12) linslipDISglobal.Assemble(valtxi1, row[0], col);
              if (abs(valteta1) > 1.0e-12) linslipDISglobal.Assemble(valteta1, row[1], col);
            }
          }

          if (Dim() == 3)
          {
            for (_colcurr = derivjump2.begin(); _colcurr != derivjump2.end(); ++_colcurr)
            {
              int col = _colcurr->first;
              double valtxi2 = (zteta + ct * jumpteta) / euclidean * ct * _colcurr->second * ztxi;
              double valteta2 = (zteta + ct * jumpteta) / euclidean * ct * _colcurr->second * zteta;

              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                {
                  if (abs(valtxi2 * txi[j]) > 1.e-12)
                    linslipDISglobal.Assemble(valtxi2 * txi[j], cnode->Dofs()[j], col);
                  if (abs(valteta2 * teta[j]) > 1.e-12)
                    linslipDISglobal.Assemble(valteta2 * teta[j], cnode->Dofs()[j], col);
                }
              }
              else
              {
                if (abs(valtxi2) > 1.0e-12) linslipDISglobal.Assemble(valtxi2, row[0], col);
                if (abs(valteta2) > 1.0e-12) linslipDISglobal.Assemble(valteta2, row[1], col);
              }
            }
          }
        }
        else  // std. slip increment
        {
          // get linearization of jump vector
          derivjump = cnode->FriData().GetDerivJump();

          // loop over dimensions
          for (int dim = 0; dim < cnode->NumDof(); ++dim)
          {
            // loop over all entries of the current derivative map (jump)
            for (_colcurr = derivjump[dim].begin(); _colcurr != derivjump[dim].end(); ++_colcurr)
            {
              int col = _colcurr->first;

              double valtxi1 =
                  (ztxi + ct * jumptxi) / euclidean * ct * txi[dim] * _colcurr->second * ztxi;
              double valteta1 =
                  (ztxi + ct * jumptxi) / euclidean * ct * txi[dim] * _colcurr->second * zteta;
              double valtxi2 =
                  (zteta + ct * jumpteta) / euclidean * ct * teta[dim] * _colcurr->second * ztxi;
              double valteta2 =
                  (zteta + ct * jumpteta) / euclidean * ct * teta[dim] * _colcurr->second * zteta;

#ifdef CONSISTENTSLIP
              valtxi1 = valtxi1 / (znor - cn * wgap);
              valteta1 = valteta1 / (znor - cn * wgap);
              valtxi2 = valtxi2 / (znor - cn * wgap);
              valteta2 = valteta2 / (znor - cn * wgap);
#endif

              // do not assemble zeros into matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                {
                  if (abs((valtxi1 + valtxi2) * txi[j]) > 1.e-12)
                    linslipDISglobal.Assemble((valtxi1 + valtxi2) * txi[j], cnode->Dofs()[j], col);
                  if (Dim() == 3)
                    if (abs((valteta1 + valteta2) * teta[j]) > 1.e-12)
                      linslipDISglobal.Assemble(
                          (valteta1 + valteta2) * teta[j], cnode->Dofs()[j], col);
                }
              }
              else
              {
                if (abs(valtxi1) > 1.0e-12) linslipDISglobal.Assemble(valtxi1, row[0], col);
                if (abs(valteta1) > 1.0e-12) linslipDISglobal.Assemble(valteta1, row[1], col);
                if (abs(valtxi2) > 1.0e-12) linslipDISglobal.Assemble(valtxi2, row[0], col);
                if (abs(valteta2) > 1.0e-12) linslipDISglobal.Assemble(valteta2, row[1], col);
              }
            }

#ifdef CONSISTENTSLIP
            /*** Additional Terms ***/
            // normal derivative
            for (colcurr = dnmap[dim].begin(); colcurr != dnmap[dim].end(); ++colcurr)
            {
              int col = colcurr->first;
              double valtxi =
                  -euclidean * ztxi * z[dim] / pow(znor - cn * wgap, 2.0) * (colcurr->second);
              double valteta =
                  -euclidean * zteta * z[dim] / pow(znor - cn * wgap, 2.0) * (colcurr->second);

              // do not assemble zeros into s matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                {
                  if (abs(valtxi * txi[j]) > 1.e-12)
                    linslipDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
                  if (Dim() == 3)
                    if (abs(valteta * teta[j]) > 1.e-12)
                      linslipDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
                }
              }
              else
              {
                if (abs(valtxi) > 1.0e-12) linslipDISglobal.Assemble(valtxi, row[0], col);
                if (abs(valteta) > 1.0e-12) linslipDISglobal.Assemble(valteta, row[1], col);
              }
            }
#endif
          }
        }

#ifdef CONSISTENTSLIP
        /*** Additional Terms ***/
        // wgap derivative
        std::map<int, double>& dgmap = cnode->CoData().GetDerivG();

        for (colcurr = dgmap.begin(); colcurr != dgmap.end(); ++colcurr)
        {
          int col = colcurr->first;
          double valtxi = +euclidean * ztxi / pow(znor - cn * wgap, 2.0) * cn * (colcurr->second);
          double valteta = +euclidean * zteta / pow(znor - cn * wgap, 2.0) * cn * (colcurr->second);

          // do not assemble zeros into matrix
          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int j = 0; j < Dim(); j++)
            {
              if (abs(valtxi * txi[j]) > 1.e-12)
                linslipDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
              if (Dim() == 3)
                if (abs(valteta * teta[j]) > 1.e-12)
                  linslipDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
            }
          }
          else
          {
            if (abs(valtxi) > 1.0e-12) linslipDISglobal.Assemble(valtxi, row[0], col);
            if (abs(valteta) > 1.0e-12) linslipDISglobal.Assemble(valteta, row[1], col);
          }
        }
#endif


        /*** 02 ***************** frcoeff*znor*ct*tangent.deriv(jump) ***/
        if (gp_slip)
        {
          for (_colcurr = derivjump1.begin(); _colcurr != derivjump1.end(); ++_colcurr)
          {
            int col = _colcurr->first;
            double valtxi = -frcoeff * (znor - cn * wgap) * ct * _colcurr->second;
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(valtxi * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(valtxi) > 1.0e-12)
              linslipDISglobal.Assemble(valtxi, row[0], col);
          }

          if (Dim() == 3)
          {
            for (_colcurr = derivjump2.begin(); _colcurr != derivjump2.end(); ++_colcurr)
            {
              int col = _colcurr->first;
              double valteta = -frcoeff * (znor - cn * wgap) * ct * _colcurr->second;

              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(valteta * teta[j]) > 1.e-12)
                    linslipDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(valteta) > 1.0e-12)
                linslipDISglobal.Assemble(valteta, row[1], col);
            }
          }
        }
        else
        {
          // loop over dimensions
          for (int dim = 0; dim < cnode->NumDof(); ++dim)
          {
            // loop over all entries of the current derivative map (jump)
            for (_colcurr = derivjump[dim].begin(); _colcurr != derivjump[dim].end(); ++_colcurr)
            {
              int col = _colcurr->first;

              // std::cout << "val " << colcurr->second << std::endl;
#ifdef CONSISTENTSLIP
              double valtxi = -frcoeff * ct * txi[dim] * colcurr->second;
              double valteta = -frcoeff * ct * teta[dim] * colcurr->second;
#else
              double valtxi =
                  (-1) * (frcoeff * (znor - cn * wgap)) * ct * txi[dim] * _colcurr->second;
              double valteta =
                  (-1) * (frcoeff * (znor - cn * wgap)) * ct * teta[dim] * _colcurr->second;
#endif
              // do not assemble zeros into matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(valtxi * txi[j]) > 1.e-12)
                    linslipDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
              }
              else if (abs(valtxi) > 1.0e-12)
                linslipDISglobal.Assemble(valtxi, row[0], col);

              if (Dim() == 3)
              {
                if (constr_direction_ == INPAR::CONTACT::constr_xyz)
                {
                  for (int j = 0; j < Dim(); j++)
                    if (abs(valteta * teta[j]) > 1.e-12)
                      linslipDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
                }
                else if (abs(valteta) > 1.0e-12)
                  linslipDISglobal.Assemble(valteta, row[1], col);
              }
            }
          }
        }
        /*** 1 ********************************* euclidean.deriv(T).z ***/
        // loop over dimensions
        for (int j = 0; j < Dim(); ++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (colcurr = dtximap[j].begin(); colcurr != dtximap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
            double val = euclidean * (colcurr->second) * z[j];

#ifdef CONSISTENTSLIP
            val = val / (znor - cn * wgap);
#endif

            // do not assemble zeros into s matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(val * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(val) > 1.0e-12)
              linslipDISglobal.Assemble(val, row[0], col);
          }

          if (Dim() == 3)
          {
            // loop over all entries of the current derivative map (teta)
            for (colcurr = dtetamap[j].begin(); colcurr != dtetamap[j].end(); ++colcurr)
            {
              int col = colcurr->first;
              double val = euclidean * (colcurr->second) * z[j];

#ifdef CONSISTENTSLIP
              val = val / (znor - cn * wgap);
#endif

              // do not assemble zeros into s matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(val * teta[j]) > 1.e-12)
                    linslipDISglobal.Assemble(val * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(val) > 1.0e-12)
                linslipDISglobal.Assemble(val, row[1], col);
            }
          }
        }

        /*** 2 ********************* deriv(euclidean).deriv(T).z.ztan ***/
        // loop over dimensions
        for (int j = 0; j < Dim(); ++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (colcurr = dtximap[j].begin(); colcurr != dtximap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
            double valtxi = (ztxi + ct * jumptxi) / euclidean * (colcurr->second) * z[j] * ztxi;
            double valteta = (ztxi + ct * jumptxi) / euclidean * (colcurr->second) * z[j] * zteta;

#ifdef CONSISTENTSLIP
            valtxi = valtxi / (znor - cn * wgap);
            valteta = valteta / (znor - cn * wgap);
#endif

            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(valtxi * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(valtxi) > 1.0e-12)
              linslipDISglobal.Assemble(valtxi, row[0], col);

            if (Dim() == 3)
            {
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(valteta * teta[j]) > 1.e-12)
                    linslipDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(valteta) > 1.0e-12)
                linslipDISglobal.Assemble(valteta, row[1], col);
            }
          }

          if (Dim() == 3)
          {
            // 3D loop over all entries of the current derivative map (teta)
            for (colcurr = dtetamap[j].begin(); colcurr != dtetamap[j].end(); ++colcurr)
            {
              int col = colcurr->first;
              double valtxi = (zteta + ct * jumpteta) / euclidean * (colcurr->second) * z[j] * ztxi;
              double valteta =
                  (zteta + ct * jumpteta) / euclidean * (colcurr->second) * z[j] * zteta;

#ifdef CONSISTENTSLIP
              valtxi = valtxi / (znor - cn * wgap);
              valteta = valteta / (znor - cn * wgap);
#endif

              // do not assemble zeros into matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(valtxi * txi[j]) > 1.e-12)
                    linslipDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
              }
              else if (abs(valtxi) > 1.0e-12)
                linslipDISglobal.Assemble(valtxi, row[0], col);

              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(valteta * teta[j]) > 1.e-12)
                    linslipDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(valteta) > 1.0e-12)
                linslipDISglobal.Assemble(valteta, row[1], col);
            }
          }
        }

        /*** 3 ****************** deriv(euclidean).deriv(T).jump.ztan ***/
        if (gp_slip)
        {
          //!!!!!!!!!!!!!!! DO NOTHING !!!!!!!
        }
        else
        {
          // loop over dimensions
          for (int j = 0; j < Dim(); ++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr = dtximap[j].begin(); colcurr != dtximap[j].end(); ++colcurr)
            {
              int col = colcurr->first;
              double valtxi =
                  (ztxi + ct * jumptxi) / euclidean * ct * (colcurr->second) * jump[j] * ztxi;
              double valteta =
                  (ztxi + ct * jumptxi) / euclidean * ct * (colcurr->second) * jump[j] * zteta;

#ifdef CONSISTENTSLIP
              valtxi = valtxi / (znor - cn * wgap);
              valteta = valteta / (znor - cn * wgap);
#endif

              // do not assemble zeros into s matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(valtxi * txi[j]) > 1.e-12)
                    linslipDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
                for (int j = 0; j < Dim(); j++)
                  if (abs(valteta * teta[j]) > 1.e-12)
                    linslipDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
              }
              else
              {
                if (abs(valtxi) > 1.0e-12) linslipDISglobal.Assemble(valtxi, row[0], col);
                if (abs(valteta) > 1.0e-12) linslipDISglobal.Assemble(valteta, row[1], col);
              }
            }

            if (Dim() == 3)
            {
              // loop over all entries of the current derivative map (teta)
              for (colcurr = dtetamap[j].begin(); colcurr != dtetamap[j].end(); ++colcurr)
              {
                int col = colcurr->first;
                double valtxi =
                    (zteta + ct * jumpteta) / euclidean * ct * (colcurr->second) * jump[j] * ztxi;
                double valteta =
                    (zteta + ct * jumpteta) / euclidean * ct * (colcurr->second) * jump[j] * zteta;

#ifdef CONSISTENTSLIP
                valtxi = valtxi / (znor - cn * wgap);
                valteta = valteta / (znor - cn * wgap);
#endif

                // do not assemble zeros into matrix
                if (constr_direction_ == INPAR::CONTACT::constr_xyz)
                {
                  for (int j = 0; j < Dim(); j++)
                    if (abs(valtxi * txi[j]) > 1.e-12)
                      linslipDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
                  for (int j = 0; j < Dim(); j++)
                    if (abs(valteta * teta[j]) > 1.e-12)
                      linslipDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
                }
                else
                {
                  if (abs(valtxi) > 1.0e-12) linslipDISglobal.Assemble(valtxi, row[0], col);
                  if (abs(valteta) > 1.0e-12) linslipDISglobal.Assemble(valteta, row[1], col);
                }
              }
            }
          }
        }
        /*** 4 ************************** (frcoeff*znor).deriv(T).z ***/
        // loop over all dimensions
        for (int j = 0; j < Dim(); ++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (colcurr = dtximap[j].begin(); colcurr != dtximap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
#ifdef CONSISTENTSLIP
            double val = -frcoeff * (colcurr->second) * z[j];
#else
            double val = (-1) * (frcoeff * (znor - cn * wgap)) * (colcurr->second) * z[j];
#endif
            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(val * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(val) > 1.0e-12)
              linslipDISglobal.Assemble(val, row[0], col);
          }

          if (Dim() == 3)
          {
            // loop over all entries of the current derivative map (teta)
            for (colcurr = dtetamap[j].begin(); colcurr != dtetamap[j].end(); ++colcurr)
            {
              int col = colcurr->first;
#ifdef CONSISTENTSLIP
              double val = -frcoeff * (colcurr->second) * z[j];
#else
              double val = (-1.0) * (frcoeff * (znor - cn * wgap)) * (colcurr->second) * z[j];
#endif
              // do not assemble zeros into matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(val * teta[j]) > 1.e-12)
                    linslipDISglobal.Assemble(val * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(val) > 1.0e-12)
                linslipDISglobal.Assemble(val, row[1], col);
            }
          }
        }

        /*** 5 *********************** (frcoeff*znor).deriv(T).jump ***/
        if (gp_slip)
        {
          //!!!!!!!!!!!!!!! DO NOTHING !!!!!!!!!!
        }
        else
        {
          // loop over all dimensions
          for (int j = 0; j < Dim(); ++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr = dtximap[j].begin(); colcurr != dtximap[j].end(); ++colcurr)
            {
              int col = colcurr->first;
#ifdef CONSISTENTSLIP
              double val = -frcoeff * ct * (colcurr->second) * jump[j];
#else
              double val = (-1) * (frcoeff * (znor - cn * wgap)) * ct * (colcurr->second) * jump[j];
#endif
              // do not assemble zeros into matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(val * txi[j]) > 1.e-12)
                    linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
              }
              else if (abs(val) > 1.0e-12)
                linslipDISglobal.Assemble(val, row[0], col);
            }

            if (Dim() == 3)
            {
              // loop over all entries of the current derivative map (teta)
              for (colcurr = dtetamap[j].begin(); colcurr != dtetamap[j].end(); ++colcurr)
              {
                int col = colcurr->first;
#ifdef CONSISTENTSLIP
                double val = -frcoeff * ct * (colcurr->second) * jump[j];
#else
                double val =
                    (-1) * (frcoeff * (znor - cn * wgap)) * ct * (colcurr->second) * jump[j];
#endif
                // do not assemble zeros into s matrix
                if (constr_direction_ == INPAR::CONTACT::constr_xyz)
                {
                  for (int j = 0; j < Dim(); j++)
                    if (abs(val * teta[j]) > 1.e-12)
                      linslipDISglobal.Assemble(val * teta[j], cnode->Dofs()[j], col);
                }
                else if (abs(val) > 1.0e-12)
                  linslipDISglobal.Assemble(val, row[1], col);
              }
            }
          }
        }

#ifndef CONSISTENTSLIP
        /*** 6 ******************* -frcoeff.Deriv(n).z(ztan+ct*utan) ***/
        // loop over all dimensions
        for (int j = 0; j < Dim(); ++j)
        {
          // loop over all entries of the current derivative map
          for (colcurr = dnmap[j].begin(); colcurr != dnmap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
            double valtxi = (-1) * (ztxi + ct * jumptxi) * frcoeff * (colcurr->second) * z[j];
            double valteta = (-1) * (zteta + ct * jumpteta) * frcoeff * (colcurr->second) * z[j];

            // do not assemble zeros into s matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(valtxi * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
              for (int j = 0; j < Dim(); j++)
                if (abs(valteta * teta[j]) > 1.e-12)
                  linslipDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
            }
            else
            {
              if (abs(valtxi) > 1.0e-12) linslipDISglobal.Assemble(valtxi, row[0], col);
              if (abs(valteta) > 1.0e-12) linslipDISglobal.Assemble(valteta, row[1], col);
            }
          }
        }

        /*** 7 ****************** frcoeff*cn*deriv (g).(ztan+ct*utan) ***/
        // prepare assembly
        std::map<int, double>& dgmap = cnode->CoData().GetDerivG();

        // loop over all entries of the current derivative map
        for (_colcurr = dgmap.begin(); _colcurr != dgmap.end(); ++_colcurr)
        {
          int col = _colcurr->first;
          double valtxi = frcoeff * cn * (_colcurr->second) * (ztxi + ct * jumptxi);
          double valteta = frcoeff * cn * (_colcurr->second) * (zteta + ct * jumpteta);

          // do not assemble zeros into matrix
          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int j = 0; j < Dim(); j++)
              if (abs(valtxi * txi[j]) > 1.e-12)
                linslipDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
            for (int j = 0; j < Dim(); j++)
              if (abs(valteta * teta[j]) > 1.e-12)
                linslipDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
          }
          else
          {
            if (abs(valtxi) > 1.0e-12) linslipDISglobal.Assemble(valtxi, row[0], col);
            if (abs(valteta) > 1.0e-12) linslipDISglobal.Assemble(valteta, row[1], col);
          }
        }
#endif
      }  // if (frcoeff==0.0)
    }    // loop over all slip nodes of the interface
  }      // Coulomb friction

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
  // Tresca Friction
  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
  if (ftype == INPAR::CONTACT::friction_tresca)
  {
    // loop over all slip nodes of the interface
    for (int i = 0; i < slipnodes_->NumMyElements(); ++i)
    {
      int gid = slipnodes_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("Cannot find node with gid %", gid);
      FriNode* cnode = dynamic_cast<FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("AssembleLinSlip: Node ownership inconsistency!");

      double ct = GetCtRef()[GetCtRef().Map().LID(cnode->Id())];

      // preparation of assembly
      // get Deriv N and calculate DerivD form DerivN

      // only for 2D so far, in this case calculation is very easy
      // dty =  dnx
      // dtx = -dny
      // FIXGIT: in the future DerivD will be called directly form node

      std::vector<CORE::GEN::pairedvector<int, double>> dnmap = cnode->CoData().GetDerivN();

      // iterator
      CORE::GEN::pairedvector<int, double>::iterator _colcurr;
      std::map<int, double>::iterator colcurr;

      std::vector<std::map<int, double>> dtmap(Dim());

      for (_colcurr = dnmap[0].begin(); _colcurr != dnmap[0].end(); _colcurr++)
        dtmap[1].insert(std::pair<int, double>(_colcurr->first, _colcurr->second));

      for (_colcurr = dnmap[1].begin(); _colcurr != dnmap[1].end(); _colcurr++)
        dtmap[0].insert(std::pair<int, double>(_colcurr->first, (-1) * _colcurr->second));

      // get more information from node
      double* jump = cnode->FriData().jump();
      double* txi = cnode->CoData().txi();
      double* xi = cnode->xspatial();
      double* z = cnode->MoData().lm();
      int row = slipt_->GID(i);

      int colsize = (int)dtmap[0].size();
      int mapsize = (int)dtmap.size();

      for (int j = 0; j < mapsize - 1; ++j)
        if ((int)dtmap[j].size() != (int)dtmap[j + 1].size())
          dserror("AssembleLinSlip: Column dim. of nodal DerivT-map is inconsistent!");

      // calculation of parts of the complementary function
      double ztan = txi[0] * z[0] + txi[1] * z[1];
      double jumptan = txi[0] * jump[0] + txi[1] * jump[1];
      // double temp = ztan + ct*jumptan;
      // double epk = frbound/abs(temp);
      // double Fpk = ztan*temp/(frbound*abs(temp));
      // double Mpk = epk*(1-Fpk);
      // double fac = 1/(abs(ztan+ct*jumptan))*1/(1-Mpk)*(-1);

      // calculation of |ztan+ct*utan|
      double sum = 0;
      int prefactor = 1;
      for (int dim = 0; dim < Dim(); dim++) sum += txi[dim] * z[dim] + ct * txi[dim] * jump[dim];

      // calculate |sum| and prefactor
      if (sum < 0)
      {
        sum = -sum;
        prefactor = (-1);
      }

      //****************************************************************
      // CONSISTENT TREATMENT OF CASE FRBOUND=0 (FRICTIONLESS)
      //****************************************************************
      // popp 08/2012
      //
      // There is a problem with our frictional nonlinear complementarity
      // function when applied to the limit case frbound=0 (frictionless).
      // In this case, the simple frictionless sliding condition should
      // be consistently recovered, which unfortunately is not the case.
      // This fact is well-known (see PhD thesis S. Hueeber) and now
      // taken care of by a special treatment as can be seen below
      //
      //****************************************************************
      if (frbound == 0.0 || (cnode->CoData().ActiveOld() == false && frilessfirst))
      {
        //**************************************************************
        // calculation of matrix entries of linearized slip condition
        //**************************************************************

        // 1) Entries from differentiation with respect to LM
        /******************************************************************/

        // loop over the dimension
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          int col = cnode->Dofs()[dim];
          double valtxi = txi[dim];

          // do not assemble zeros into matrix
          if (abs(valtxi) > 1.0e-12) linslipLMglobal.Assemble(valtxi, row, col);
        }

        // 2) Entries on right hand side
        /******************************************************************/

        CORE::LINALG::SerialDenseVector rhsnode(1);
        std::vector<int> lm(1);
        std::vector<int> lmowner(1);

        rhsnode(0) = -ztan;
        lm[0] = cnode->Dofs()[1];
        lmowner[0] = cnode->Owner();

        CORE::LINALG::Assemble(linslipRHSglobal, rhsnode, lm, lmowner);

        // 3) Entries from differentiation with respect to displacements
        /******************************************************************/

        // loop over dimensions
        for (int j = 0; j < Dim(); ++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (colcurr = dtmap[j].begin(); colcurr != dtmap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
            double val = (colcurr->second) * z[j];

            // do not assemble zeros into matrix
            if (abs(val) > 1.0e-12) linslipDISglobal.Assemble(val, row, col);
          }
        }
      }

      //****************************************************************
      // STANDARD TREATMENT OF CASE FRBOUND!=0 (FRICTIONAL)
      //****************************************************************
      else
      {
        //****************************************************************
        // calculation of matrix entries of the linearized slip condition
        //****************************************************************

        // 1) Entries from differentiation with respect to LM
        /**************** (Deriv(abs)*ztan+|ztan+ct*jumptan|-frbound).tan ***/

        // loop over the dimension
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          int col = cnode->Dofs()[dim];
          double val = (prefactor * ztan + sum - frbound) * txi[dim];

          // do not assemble zeros into matrix
          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int j = 0; j < Dim(); j++)
              if (abs(val * txi[j]) > 1.e-12)
                linslipLMglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
          }
          else if (abs(val) > 1.0e-12)
            linslipLMglobal.Assemble(val, row, col);
        }

        // 2) Entries on right hand side
        /************ -C + entries from writing Delta(z) as z(k+1)-z(k) ***/

        // -C and remaining terms
        double value1 = -(abs(ztan + ct * jumptan)) * ztan + frbound * (ztan + ct * jumptan);

        if (constr_direction_ == INPAR::CONTACT::constr_xyz)
        {
          CORE::LINALG::SerialDenseVector rhsnode(Dim());
          std::vector<int> lm(Dim());
          std::vector<int> lmowner(Dim());

          for (int j = 0; j < Dim(); j++)
          {
            lm[j] = cnode->Dofs()[j];
            lmowner[j] = cnode->Owner();
            rhsnode(j) = value1 * txi[j];
          }

          CORE::LINALG::Assemble(linslipRHSglobal, rhsnode, lm, lmowner);
        }
        else
        {
          CORE::LINALG::SerialDenseVector rhsnode(1);
          rhsnode(0) = value1;

          std::vector<int> lm(1);
          std::vector<int> lmowner(1);

          lm[0] = cnode->Dofs()[1];
          lmowner[0] = cnode->Owner();

          CORE::LINALG::Assemble(linslipRHSglobal, rhsnode, lm, lmowner);
        }

        // 3) Entries from differentiation with respect to displacements
        /***************************** -Deriv(abs)*ct*tan.(D-Dn-1)*ztan ***/

        // we need the nodal entries of the D-matrix and the old one
        double D = cnode->MoData().GetD()[cnode->Id()];
        double Dold = cnode->FriData().GetDOld()[cnode->Id()];

        if (abs(Dold) < 0.0001) dserror("Error:No entry for Dold");

        // loop over all derivative maps (=dimensions)
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          int col = cnode->Dofs()[dim];
          double val = prefactor * (-1) * ct * txi[dim] * (D - Dold) * ztan;
          // std::cout << "01 GID " << gid << " row " << row << " col " << col << " val " << val <<
          // std::endl;

          // do not assemble zeros into matrix
          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int j = 0; j < Dim(); j++)
              if (abs(val * txi[j]) > 1.e-12)
                linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
          }
          else if (abs(val) > 1.0e-12)
            linslipDISglobal.Assemble(val, row, col);
        }

        /***************************** -Deriv(abs)*ct*tan.(M-Mn-1)*ztan ***/

        // we need the nodal entries of the M-matrix and the old one
        std::map<int, double>& mmap = cnode->MoData().GetM();
        std::map<int, double>& mmapold = cnode->FriData().GetMOld();

        // create a set of nodes including nodes according to M entries
        // from current and previous time step
        std::set<int> mnodes;

        // iterator
        std::set<int>::iterator mcurr;

        std::set<int> mnodescurrent = cnode->FriData().GetMNodes();
        std::set<int> mnodesold = cnode->FriData().GetMNodesOld();

        for (mcurr = mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
          mnodes.insert(*mcurr);

        for (mcurr = mnodesold.begin(); mcurr != mnodesold.end(); mcurr++) mnodes.insert(*mcurr);

        // loop over all master nodes (find adjacent ones to this stick node)
        for (mcurr = mnodes.begin(); mcurr != mnodes.end(); mcurr++)
        {
          int gid = *mcurr;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("Cannot find node with gid %", gid);
          FriNode* cmnode = dynamic_cast<FriNode*>(mnode);

          double mik = mmap[cmnode->Id()];
          double mikold = mmapold[cmnode->Id()];

          // compute linstick-matrix entry of the current active node / master node pair
          // loop over all derivative maps (=dimensions)
          for (int dim = 0; dim < cnode->NumDof(); ++dim)
          {
            int col = cmnode->Dofs()[dim];
            double val = prefactor * (+1) * ct * txi[dim] * (mik - mikold) * ztan;
            // std::cout << "02 GID " << gid << " row " << row << " col " << col << " val " << val
            // << std::endl;

            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(val * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(val) > 1.0e-12)
              linslipDISglobal.Assemble(val, row, col);
          }
        }

        /************************************** frbound*ct*tan.(D-Dn-1) ***/

        // loop over all derivative maps (=dimensions)
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          int col = cnode->Dofs()[dim];
          double val = frbound * ct * txi[dim] * (D - Dold);
          // std::cout << "03 GID " << gid << " row " << row << " col " << col << " val " << val <<
          // std::endl;

          // do not assemble zeros into matrix
          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int j = 0; j < Dim(); j++)
              if (abs(val * txi[j]) > 1.e-12)
                linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
          }
          else if (abs(val) > 1.0e-12)
            linslipDISglobal.Assemble(val, row, col);
        }

        /********************************** -frbound*ct*tan.(M-Mn-1).xm ***/

        // loop over all master nodes
        for (mcurr = mnodes.begin(); mcurr != mnodes.end(); mcurr++)
        {
          int gid = *mcurr;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("Cannot find node with gid %", gid);
          FriNode* cmnode = dynamic_cast<FriNode*>(mnode);

          double mik = mmap[cmnode->Id()];
          double mikold = mmapold[cmnode->Id()];

          // loop over all derivative maps (=dimensions)
          for (int dim = 0; dim < cnode->NumDof(); ++dim)
          {
            int col = cmnode->Dofs()[dim];
            double val = frbound * (-1) * ct * txi[dim] * (mik - mikold);
            // std::cout << "04 GID " << gid << " row " << row << " col " << col << " val " << val
            // << std::endl;
            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(val * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(val) > 1.0e-12)
              linslipDISglobal.Assemble(val, row, col);
          }
        }

        /************************************ |ztan+ct*utan|.DerivT.z ***/

        // loop over all derivative maps (=dimensions)
        for (int j = 0; j < mapsize; ++j)
        {
          int k = 0;

          // loop over all entries of the current derivative map
          for (colcurr = dtmap[j].begin(); colcurr != dtmap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
            double val = sum * (colcurr->second) * z[j];
            // std::cout << "1 GID " << gid << " row " << row << " col " << col << " val " << val <<
            // std::endl;

            // do not assemble zeros into s matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(val * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(val) > 1.0e-12)
              linslipDISglobal.Assemble(val, row, col);

            ++k;
          }

          if (k != colsize) dserror("AssembleLinSlip: k = %i but colsize = %i", k, colsize);
        }

        /*********************************** Deriv(abs)*DerivT.z*ztan ***/

        // loop over all derivative maps (=dimensions)
        for (int j = 0; j < mapsize; ++j)
        {
          int k = 0;

          // loop over all entries of the current derivative map
          for (colcurr = dtmap[j].begin(); colcurr != dtmap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor * (colcurr->second) * z[j] * ztan;
            // std::cout << "2 GID " << gid << " row " << row << " col " << col << " val " << val <<
            // std::endl;

            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(val * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(val) > 1.0e-12)
              linslipDISglobal.Assemble(val, row, col);

            ++k;
          }

          if (k != colsize) dserror("AssembleLinSlip: k = %i but colsize = %i", k, colsize);
        }

        /******************************* Deriv(abs)*DerivT.jump+*ztan ***/

        // loop over all derivative maps (=dimensions)
        for (int j = 0; j < mapsize; ++j)
        {
          int k = 0;

          // loop over all entries of the current derivative map
          for (colcurr = dtmap[j].begin(); colcurr != dtmap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor * ct * (colcurr->second) * jump[j] * ztan;
            // std::cout << "3 GID " << gid << " row " << row << " col " << col << " val " << val <<
            // std::endl;

            // do not assemble zeros into s matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(val * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(val) > 1.0e-12)
              linslipDISglobal.Assemble(val, row, col);

            ++k;
          }

          if (k != colsize) dserror("AssembleLinSlip: k = %i but colsize = %i", k, colsize);
        }

        /*************************** -Deriv(abs).ct.tan.DerivD.x*ztan ***/

        // we need the dot product t*x of this node
        double tdotx = 0.0;
        for (int dim = 0; dim < cnode->NumDof(); ++dim) tdotx += txi[dim] * xi[dim];

        // prepare assembly
        std::map<int, double>& ddmap = cnode->CoData().GetDerivD()[gid];

        // loop over all entries of the current derivative map
        for (colcurr = ddmap.begin(); colcurr != ddmap.end(); ++colcurr)
        {
          int col = colcurr->first;
          double val = (-1) * prefactor * ct * tdotx * colcurr->second * ztan;
          // std::cout << "4 GID " << gid << " row " << row << " col " << col << " val " << val <<
          // std::endl;

          // do not assemble zeros into matrix
          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int j = 0; j < Dim(); j++)
              if (abs(val * txi[j]) > 1.e-12)
                linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
          }
          else if (abs(val) > 1.0e-12)
            linslipDISglobal.Assemble(val, row, col);
        }

        /**************************** Deriv(abs).ct.tan.DerivM.x*ztan ***/

        // we need the Lin(M-matrix) entries of this node
        std::map<int, std::map<int, double>>& dmmap = cnode->CoData().GetDerivM();
        std::map<int, std::map<int, double>>::iterator dmcurr;

        // loop over all master nodes in the DerivM-map of the active slave node
        for (dmcurr = dmmap.begin(); dmcurr != dmmap.end(); ++dmcurr)
        {
          int gid = dmcurr->first;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("Cannot find node with gid %", gid);
          FriNode* cmnode = dynamic_cast<FriNode*>(mnode);
          double* mxi = cmnode->xspatial();

          // we need the dot product ns*xm of this node pair
          double tdotx = 0.0;
          for (int dim = 0; dim < cnode->NumDof(); ++dim) tdotx += txi[dim] * mxi[dim];

          // compute entry of the current active node / master node pair
          std::map<int, double>& thisdmmap = cnode->CoData().GetDerivM(gid);

          // loop over all entries of the current derivative map
          for (colcurr = thisdmmap.begin(); colcurr != thisdmmap.end(); ++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor * ct * tdotx * colcurr->second * ztan;
            // std::cout << "5 GID " << gid << " row " << row << " col " << col << " val " << val <<
            // std::endl;

            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(val * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(val) > 1.0e-12)
              linslipDISglobal.Assemble(val, row, col);
          }
        }

        /****************************************** -frbound.DerivT.z ***/

        // loop over all derivative maps (=dimensions)
        for (int j = 0; j < mapsize; ++j)
        {
          int k = 0;

          // loop over all entries of the current derivative map
          for (colcurr = dtmap[j].begin(); colcurr != dtmap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
            double val = (-1) * frbound * (colcurr->second) * z[j];
            // std::cout << "6 GID " << gid << " row " << row << " col " << col << " val " << val <<
            // std::endl;

            // do not assemble zeros into s matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(val * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(val) > 1.0e-12)
              linslipDISglobal.Assemble(val, row, col);

            ++k;
          }

          if (k != colsize) dserror("AssembleLinSlip: k = %i but colsize = %i", k, colsize);
        }

        /************************************ -frbound.ct.DerivT.jump ***/

        // loop over all derivative maps (=dimensions)
        for (int j = 0; j < mapsize; ++j)
        {
          int k = 0;

          // loop over all entries of the current derivative map
          for (colcurr = dtmap[j].begin(); colcurr != dtmap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
            double val = (-1) * frbound * ct * (colcurr->second) * jump[j];
            // std::cout << "7 GID " << gid << " row " << row << " col " << col << " val " << val <<
            // std::endl;

            // do not assemble zeros into s matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(val * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(val) > 1.0e-12)
              linslipDISglobal.Assemble(val, row, col);

            ++k;
          }

          if (k != colsize) dserror("AssembleLinSlip: k = %i but colsize = %i", k, colsize);
        }

        /************************************* +frbound.ct.T.DerivD.x ***/

        // we need the dot product t*x of this node
        tdotx = 0.0;
        for (int dim = 0; dim < cnode->NumDof(); ++dim) tdotx += txi[dim] * xi[dim];

        // loop over all entries of the current derivative map
        for (colcurr = ddmap.begin(); colcurr != ddmap.end(); ++colcurr)
        {
          int col = colcurr->first;
          double val = (-1) * (-1) * frbound * ct * tdotx * colcurr->second;
          // std::cout << "8 GID " << gid << " row " << row << " col " << col << " val " << val <<
          // std::endl;

          // do not assemble zeros into matrix
          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int j = 0; j < Dim(); j++)
              if (abs(val * txi[j]) > 1.e-12)
                linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
          }
          else if (abs(val) > 1.0e-12)
            linslipDISglobal.Assemble(val, row, col);
        }

        /********************************  -frbound.ct.T.DerivM.x ******/

        // loop over all master nodes in the DerivM-map of the active slave node
        for (dmcurr = dmmap.begin(); dmcurr != dmmap.end(); ++dmcurr)
        {
          int gid = dmcurr->first;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("Cannot find node with gid %", gid);
          FriNode* cmnode = dynamic_cast<FriNode*>(mnode);
          double* mxi = cmnode->xspatial();

          // we need the dot product ns*xm of this node pair
          double tdotx = 0.0;
          for (int dim = 0; dim < cnode->NumDof(); ++dim) tdotx += txi[dim] * mxi[dim];

          // compute entry of the current active node / master node pair
          std::map<int, double>& thisdmmap = cnode->CoData().GetDerivM(gid);

          // loop over all entries of the current derivative map
          for (colcurr = thisdmmap.begin(); colcurr != thisdmmap.end(); ++colcurr)
          {
            int col = colcurr->first;
            double val = (-1) * frbound * ct * tdotx * colcurr->second;
            // std::cout << "9 GID " << gid << " row " << row << " col " << col << " val " << val <<
            // std::endl;


            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(val * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(val) > 1.0e-12)
              linslipDISglobal.Assemble(val, row, col);
          }
        }
      }
    }
  }  // Tresca friction

  return;
}

void CONTACT::CoInterface::AssembleNormalContactRegularization(
    CORE::LINALG::SparseMatrix& d_disp, CORE::LINALG::SparseMatrix& d_lm, Epetra_Vector& f)
{
  const bool regularization =
      DRT::INPUT::IntegralValue<int>(InterfaceParams(), "REGULARIZED_NORMAL_CONTACT");
  if (!regularization) dserror("you should not be here");
  const double k = 1. / InterfaceParams().get<double>("REGULARIZATION_STIFFNESS");
  const double gmax = InterfaceParams().get<double>("REGULARIZATION_THICKNESS");
  const int dim = Dim();
  static const INPAR::CONTACT::ConstraintDirection constr_direction =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::ConstraintDirection>(
          InterfaceParams(), "CONSTRAINT_DIRECTIONS");

  for (int i = 0; i < ActiveNodes()->NumMyElements(); ++i)
  {
    DRT::Node* node = Discret().gNode(ActiveNodes()->GID(i));
    if (!node) dserror("node not found");
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
    if (!cnode) dserror("not a contact node");

    CORE::LINALG::Matrix<3, 1> n(cnode->MoData().n(), true);
    CORE::LINALG::Matrix<3, 1> lm(cnode->MoData().lm(), true);
    const double lm_n = lm.Dot(n);

    const double gLM = gmax * (1. - exp(-k / gmax * lm_n));
    const double d_gLM = k * exp(-k / gmax * lm_n);

    if (cnode->MoData().GetD().size() != 1)
      dserror(
          "we need to have a D-value for active contact nodes\nAnd exactly one due to "
          "biorthogonality");
    double dval = cnode->MoData().GetD().at(cnode->Id());

    if (constr_direction == INPAR::CONTACT::constr_xyz)
      for (int d = 0; d < dim; ++d) f[f.Map().LID(cnode->Dofs()[d])] += n(d) * dval * gLM;
    else
      f[f.Map().LID(cnode->Dofs()[0])] += dval * gLM;

    for (int l = 0; l < dim; ++l)
    {
      if (constr_direction == INPAR::CONTACT::constr_xyz)
        for (int d = 0; d < dim; ++d)
          d_lm.Assemble(n(d) * dval * n(l) * d_gLM, cnode->Dofs()[d], cnode->Dofs()[l]);
      else
        d_lm.Assemble(dval * d_gLM * n(l), cnode->Dofs()[0], cnode->Dofs()[l]);
    }

    for (auto p = cnode->CoData().GetDerivD().at(cnode->Id()).begin();
         p != cnode->CoData().GetDerivD().at(cnode->Id()).end(); ++p)
    {
      const int col = p->first;
      const double val = gLM * p->second;
      if (constr_direction == INPAR::CONTACT::constr_xyz)
        for (int d = 0; d < dim; ++d) d_disp.Assemble(val * n(d), cnode->Dofs()[d], col);
      else
        d_disp.Assemble(val, cnode->Dofs()[0], col);
    }

    for (int l = 0; l < dim; ++l)
      for (auto p = cnode->CoData().GetDerivN().at(l).begin();
           p != cnode->CoData().GetDerivN().at(l).end(); ++p)
      {
        const int col = p->first;
        const double val = dval * d_gLM * lm(l) * p->second;
        if (constr_direction == INPAR::CONTACT::constr_xyz)
          for (int d = 0; d < dim; ++d) d_disp.Assemble(val * n(d), cnode->Dofs()[d], col);
        else
          d_disp.Assemble(val, cnode->Dofs()[0], col);
      }
  }

  return;
}

void CONTACT::CoInterface::AssembleLinSlipNormalRegularization(
    CORE::LINALG::SparseMatrix& linslipLMglobal, CORE::LINALG::SparseMatrix& linslipDISglobal,
    Epetra_Vector& linslipRHSglobal)
{
  // nothing to do if no slip nodes
  if (slipnodes_->NumMyElements() == 0) return;

  // information from interface contact parameter list
  INPAR::CONTACT::FrictionType ftype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(InterfaceParams(), "FRICTION");
  double frcoeff_in =
      InterfaceParams().get<double>("FRCOEFF");  // the friction coefficient from the input
  bool gp_slip = DRT::INPUT::IntegralValue<int>(InterfaceParams(), "GP_SLIP_INCR");
  if (gp_slip) dserror("not implemented");
  bool frilessfirst = DRT::INPUT::IntegralValue<int>(InterfaceParams(), "FRLESS_FIRST");

  // the friction coefficient adapted by every node (eg depending on the local temperature)
  double frcoeff = 0.;

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
  // Coulomb Friction
  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
  if (ftype == INPAR::CONTACT::friction_coulomb)
  {
    // loop over all slip nodes of the interface
    for (int i = 0; i < slipnodes_->NumMyElements(); ++i)
    {
      int gid = slipnodes_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("Cannot find node with gid %", gid);
      FriNode* cnode = dynamic_cast<FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("AssembleLinSlip: Node ownership inconsistency!");

      // get friction coefficient for this node
      // in case of TSI, the nodal temperature influences the local friction coefficient
      frcoeff = cnode->FrCoeff(frcoeff_in);

      double cn = GetCnRef()[GetCnRef().Map().LID(cnode->Id())];
      double ct = GetCtRef()[GetCtRef().Map().LID(cnode->Id())];

      // prepare assembly, get information from node
      std::vector<CORE::GEN::pairedvector<int, double>> dnmap = cnode->CoData().GetDerivN();
      std::vector<CORE::GEN::pairedvector<int, double>> dtximap = cnode->CoData().GetDerivTxi();
      std::vector<CORE::GEN::pairedvector<int, double>> dtetamap = cnode->CoData().GetDerivTeta();

      // check for Dimension of derivative maps
      for (int j = 0; j < Dim() - 1; ++j)
        if ((int)dnmap[j].size() != (int)dnmap[j + 1].size())
          dserror("AssembleLinSlip: Column dim. of nodal DerivTxi-map is inconsistent!");

      for (int j = 0; j < Dim() - 1; ++j)
        if ((int)dtximap[j].size() != (int)dtximap[j + 1].size())
          dserror("AssembleLinSlip: Column dim. of nodal DerivTxi-map is inconsistent!");

      if (Dim() == 3)
      {
        for (int j = 0; j < Dim() - 1; ++j)
          if ((int)dtximap[j].size() != (int)dtximap[j + 1].size())
            dserror("AssembleLinSlip: Column dim. of nodal DerivTeta-map is inconsistent!");
      }

      // more information from node
      double* n = cnode->MoData().n();
      double* txi = cnode->CoData().txi();
      double* teta = cnode->CoData().teta();
      double* z = cnode->MoData().lm();
      double wgap = cnode->CoData().Getg();

      // iterator for maps
      CORE::GEN::pairedvector<int, double>::iterator colcurr;
      std::map<int, double>::iterator _colcurr;

      // row number of entries
      std::vector<int> row(Dim() - 1);
      if (Dim() == 2)
      {
        row[0] = slipt_->GID(i);
      }
      else if (Dim() == 3)
      {
        row[0] = slipt_->GID(2 * i);
        row[1] = slipt_->GID(2 * i) + 1;
      }
      else
        dserror("AssemblelinSlip: Dimension not correct");

      // evaluation of specific components of entries to assemble
      double znor = 0.0;
      double ztxi = 0.0;
      double zteta = 0.0;
      double jumptxi = 0.0;
      double jumpteta = 0.0;
      double euclidean = 0.0;
      double* jump = cnode->FriData().jump();

      for (int i = 0; i < Dim(); i++)
      {
        znor += n[i] * z[i];
        ztxi += txi[i] * z[i];
        zteta += teta[i] * z[i];
        jumptxi += txi[i] * jump[i];
        jumpteta += teta[i] * jump[i];
      }

      // setup regularization
      static const bool regularization =
          DRT::INPUT::IntegralValue<int>(InterfaceParams(), "REGULARIZED_NORMAL_CONTACT");
      if (!regularization) dserror("you should not be here");
      static const double k = 1. / InterfaceParams().get<double>("REGULARIZATION_STIFFNESS");
      static const double gmax = InterfaceParams().get<double>("REGULARIZATION_THICKNESS");
      if (cnode->MoData().GetD().size() != 1)
        dserror(
            "we need to have a D-value for active contact nodes\nAnd exactly one due to "
            "biorthogonality");
      double dval = cnode->MoData().GetD().at(cnode->Id());
      const double gLM = gmax * (1. - exp(-k / gmax * znor));
      const double d_gLM = k * exp(-k / gmax * znor);
      const double frbound = frcoeff * std::max(0., znor - cn * (wgap + dval * gLM));

      //      std::cout << "node: " << cnode->Id() << " wgap: " << wgap << " dval*gLM: " << dval*gLM
      //      << std::endl;

      // evaluate euclidean norm ||vec(zt)+ct*vec(jumpt)||
      std::vector<double> sum1(Dim() - 1, 0);
      sum1[0] = ztxi + ct * jumptxi;
      if (Dim() == 3) sum1[1] = zteta + ct * jumpteta;
      if (Dim() == 2) euclidean = abs(sum1[0]);
      if (Dim() == 3) euclidean = sqrt(sum1[0] * sum1[0] + sum1[1] * sum1[1]);

      // check of dimensions
      if (Dim() == 2 and (zteta != 0.0 or jumpteta != 0.0))
        dserror("AssemblelinSlip: zteta and jumpteta must be zero in 2D");

      //****************************************************************
      // CONSISTENT TREATMENT OF CASE FRCOEFF=0 (FRICTIONLESS)
      //****************************************************************
      // popp 08/2012
      //
      // There is a problem with our frictional nonlinear complementarity
      // function when applied to the limit case frcoeff=0 (frictionless).
      // In this case, the simple frictionless sliding condition should
      // be consistently recovered, which unfortunately is not the case.
      // This fact is well-known (see PhD thesis S. Hueeber) and now
      // taken care of by a special treatment as can be seen below
      //
      // seitz11/13
      // Like that FRLESS_FIRST results in a frictionless constraint.
      // Nodes with a vanishing slip and tangential Lagrange multiplier
      // are treated as frictionless as well.
      //
      //****************************************************************
      if (frbound == 0.0 || (cnode->CoData().ActiveOld() == false && frilessfirst) ||
          euclidean == 0.0)
      {
        //**************************************************************
        // calculation of matrix entries of linearized slip condition
        //**************************************************************

        // 1) Entries from differentiation with respect to LM
        /******************************************************************/

        // loop over the dimension
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          int col = cnode->Dofs()[dim];
          double valtxi = txi[dim];

          double valteta = 0.0;
          if (Dim() == 3) valteta = teta[dim];

          // do not assemble zeros into matrix
          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int j = 0; j < Dim(); j++)
            {
              if (abs(valtxi * txi[j]) > 1.e-12)
                linslipLMglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
              if (Dim() == 3)
                if (abs(valteta * teta[j]) > 1.e-12)
                  linslipLMglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
            }
          }
          else
          {
            if (abs(valtxi) > 1.0e-12) linslipLMglobal.Assemble(valtxi, row[0], col);
            if (Dim() == 3)
              if (abs(valteta) > 1.0e-12) linslipLMglobal.Assemble(valteta, row[1], col);
          }
        }

        // 2) Entries on right hand side
        /******************************************************************/
        if (constr_direction_ == INPAR::CONTACT::constr_xyz)
        {
          CORE::LINALG::SerialDenseVector rhsnode(Dim());
          std::vector<int> lm(Dim());
          std::vector<int> lmowner(Dim());
          for (int j = 0; j < Dim(); j++)
          {
            lm[j] = cnode->Dofs()[j];
            lmowner[j] = cnode->Owner();
            rhsnode[j] -= ztxi * txi[j];
            if (Dim() == 3) rhsnode[j] -= zteta * teta[j];
          }
          CORE::LINALG::Assemble(linslipRHSglobal, rhsnode, lm, lmowner);
        }
        else
        {
          CORE::LINALG::SerialDenseVector rhsnode(Dim() - 1);
          std::vector<int> lm(Dim() - 1);
          std::vector<int> lmowner(Dim() - 1);

          lm[0] = cnode->Dofs()[1];
          lmowner[0] = cnode->Owner();

          rhsnode[0] = -ztxi;  // already negative rhs!!!

          if (Dim() == 3)
          {
            rhsnode[1] = -zteta;  // already negative rhs!!!

            lm[1] = cnode->Dofs()[2];
            lmowner[1] = cnode->Owner();
          }
          CORE::LINALG::Assemble(linslipRHSglobal, rhsnode, lm, lmowner);
        }

        // 3) Entries from differentiation with respect to displacements
        /******************************************************************/

        // loop over dimensions
        for (int j = 0; j < Dim(); ++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (colcurr = dtximap[j].begin(); colcurr != dtximap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
            double val = (colcurr->second) * z[j];

            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(val * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(val) > 1.0e-12)
              linslipDISglobal.Assemble(val, row[0], col);
          }

          if (Dim() == 3)
          {
            // loop over all entries of the current derivative map (teta)
            for (colcurr = dtetamap[j].begin(); colcurr != dtetamap[j].end(); ++colcurr)
            {
              int col = colcurr->first;
              double val = (colcurr->second) * z[j];

              // do not assemble zeros into s matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(val * teta[j]) > 1.e-12)
                    linslipDISglobal.Assemble(val * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(val) > 1.0e-12)
                linslipDISglobal.Assemble(val, row[1], col);
            }
          }
        }
      }

      //****************************************************************
      // STANDARD TREATMENT OF CASE FRCOEFF!=0 (FRICTIONAL)
      //****************************************************************
      else
      {
        //**************************************************************
        // calculation of matrix entries of linearized slip condition
        //**************************************************************


        // 1) Entries from differentiation with respect to LM
        /******************************************************************/

        // loop over the dimension
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          double valtxi = 0.0;
          int col = cnode->Dofs()[dim];

          double valtxi0 = euclidean * txi[dim];
          double valtxi1 = ((ztxi + ct * jumptxi) / euclidean * ztxi) * txi[dim];
          double valtxi3 = (zteta + ct * jumpteta) / euclidean * ztxi * teta[dim];

          double valtxi2 = -frbound * txi[dim] -
                           frcoeff * (1. - cn * dval * d_gLM) * (ztxi + ct * jumptxi) * n[dim];
          valtxi = valtxi0 + valtxi1 + valtxi2 + valtxi3;

          double valteta = 0.0;
          if (Dim() == 3)
          {
            double valteta0 = euclidean * teta[dim];
            double valteta1 = ((ztxi + ct * jumptxi) / euclidean * zteta) * txi[dim];
            double valteta3 = (zteta + ct * jumpteta) / euclidean * zteta * teta[dim];
            double valteta2 = -frbound * teta[dim] -
                              frcoeff * (1. - cn * dval * d_gLM) * (zteta + ct * jumpteta) * n[dim];
            valteta = valteta0 + valteta1 + valteta2 + valteta3;
          }

          // do not assemble zeros into matrix
          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int j = 0; j < Dim(); j++)
            {
              if (abs(valtxi * txi[j]) > 1.e-12)
                linslipLMglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
              if (Dim() == 3)
                if (abs(valteta * teta[j]) > 1.e-12)
                  linslipLMglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
            }
          }
          else
          {
            if (abs(valtxi) > 1.0e-12) linslipLMglobal.Assemble(valtxi, row[0], col);
            if (Dim() == 3)
              if (abs(valteta) > 1.0e-12) linslipLMglobal.Assemble(valteta, row[1], col);
          }
        }

        // 2) Entries on right hand side
        /******************************************************************/
        if (constr_direction_ == INPAR::CONTACT::constr_xyz)
        {
          CORE::LINALG::SerialDenseVector rhsnode(Dim());
          std::vector<int> lm(Dim());
          std::vector<int> lmowner(Dim());

          double valuetxi1 = -(euclidean)*ztxi + (frbound) * (ztxi + ct * jumptxi);

          for (int j = 0; j < Dim(); j++)
          {
            lm[j] = cnode->Dofs()[j];
            lmowner[j] = cnode->Owner();
            rhsnode(j) += valuetxi1 * txi[j];
          }

          if (Dim() == 3)
          {
            double valueteta1 = -(euclidean)*zteta + (frbound) * (zteta + ct * jumpteta);

            for (int j = 0; j < Dim(); j++) rhsnode(j) += valueteta1 * teta[j];
          }
          CORE::LINALG::Assemble(linslipRHSglobal, rhsnode, lm, lmowner);
        }
        else
        {
          CORE::LINALG::SerialDenseVector rhsnode(Dim() - 1);
          std::vector<int> lm(Dim() - 1);
          std::vector<int> lmowner(Dim() - 1);
          double valuetxi1 = -(euclidean)*ztxi + (frbound) * (ztxi + ct * jumptxi);
          rhsnode(0) = valuetxi1;
          lm[0] = cnode->Dofs()[1];
          lmowner[0] = cnode->Owner();

          if (Dim() == 3)
          {
            double valueteta1 = -(euclidean)*zteta + (frbound) * (zteta + ct * jumpteta);
            rhsnode(1) = valueteta1;

            lm[1] = cnode->Dofs()[2];
            lmowner[1] = cnode->Owner();
          }
          CORE::LINALG::Assemble(linslipRHSglobal, rhsnode, lm, lmowner);
        }

        // 3) Entries from differentiation with respect to displacements
        /******************************************************************/
        std::map<int, double> derivjump1, derivjump2;  // for gp slip
        std::vector<std::map<int, double>> derivjump;  // for dm slip

        // get linearization of jump vector
        derivjump = cnode->FriData().GetDerivJump();

        // loop over dimensions
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          // loop over all entries of the current derivative map (jump)
          for (_colcurr = derivjump[dim].begin(); _colcurr != derivjump[dim].end(); ++_colcurr)
          {
            int col = _colcurr->first;

            double valtxi1 =
                (ztxi + ct * jumptxi) / euclidean * ct * txi[dim] * _colcurr->second * ztxi;
            double valteta1 =
                (ztxi + ct * jumptxi) / euclidean * ct * txi[dim] * _colcurr->second * zteta;
            double valtxi2 =
                (zteta + ct * jumpteta) / euclidean * ct * teta[dim] * _colcurr->second * ztxi;
            double valteta2 =
                (zteta + ct * jumpteta) / euclidean * ct * teta[dim] * _colcurr->second * zteta;

            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
              {
                if (abs((valtxi1 + valtxi2) * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble((valtxi1 + valtxi2) * txi[j], cnode->Dofs()[j], col);
                if (Dim() == 3)
                  if (abs((valteta1 + valteta2) * teta[j]) > 1.e-12)
                    linslipDISglobal.Assemble(
                        (valteta1 + valteta2) * teta[j], cnode->Dofs()[j], col);
              }
            }
            else
            {
              if (abs(valtxi1) > 1.0e-12) linslipDISglobal.Assemble(valtxi1, row[0], col);
              if (abs(valteta1) > 1.0e-12) linslipDISglobal.Assemble(valteta1, row[1], col);
              if (abs(valtxi2) > 1.0e-12) linslipDISglobal.Assemble(valtxi2, row[0], col);
              if (abs(valteta2) > 1.0e-12) linslipDISglobal.Assemble(valteta2, row[1], col);
            }
          }
        }

        // loop over dimensions
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          // loop over all entries of the current derivative map (jump)
          for (_colcurr = derivjump[dim].begin(); _colcurr != derivjump[dim].end(); ++_colcurr)
          {
            int col = _colcurr->first;

            // std::cout << "val " << colcurr->second << std::endl;
            double valtxi = (-1) * (frbound)*ct * txi[dim] * _colcurr->second;
            double valteta = (-1) * (frbound)*ct * teta[dim] * _colcurr->second;

            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(valtxi * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(valtxi) > 1.0e-12)
              linslipDISglobal.Assemble(valtxi, row[0], col);

            if (Dim() == 3)
            {
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(valteta * teta[j]) > 1.e-12)
                    linslipDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(valteta) > 1.0e-12)
                linslipDISglobal.Assemble(valteta, row[1], col);
            }
          }
        }
        /*** 1 ********************************* euclidean.deriv(T).z ***/
        // loop over dimensions
        for (int j = 0; j < Dim(); ++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (colcurr = dtximap[j].begin(); colcurr != dtximap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
            double val = euclidean * (colcurr->second) * z[j];

            // do not assemble zeros into s matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(val * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(val) > 1.0e-12)
              linslipDISglobal.Assemble(val, row[0], col);
          }

          if (Dim() == 3)
          {
            // loop over all entries of the current derivative map (teta)
            for (colcurr = dtetamap[j].begin(); colcurr != dtetamap[j].end(); ++colcurr)
            {
              int col = colcurr->first;
              double val = euclidean * (colcurr->second) * z[j];

              // do not assemble zeros into s matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(val * teta[j]) > 1.e-12)
                    linslipDISglobal.Assemble(val * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(val) > 1.0e-12)
                linslipDISglobal.Assemble(val, row[1], col);
            }
          }
        }

        /*** 2 ********************* deriv(euclidean).deriv(T).z.ztan ***/
        // loop over dimensions
        for (int j = 0; j < Dim(); ++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (colcurr = dtximap[j].begin(); colcurr != dtximap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
            double valtxi = (ztxi + ct * jumptxi) / euclidean * (colcurr->second) * z[j] * ztxi;
            double valteta = (ztxi + ct * jumptxi) / euclidean * (colcurr->second) * z[j] * zteta;

            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(valtxi * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(valtxi) > 1.0e-12)
              linslipDISglobal.Assemble(valtxi, row[0], col);

            if (Dim() == 3)
            {
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(valteta * teta[j]) > 1.e-12)
                    linslipDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(valteta) > 1.0e-12)
                linslipDISglobal.Assemble(valteta, row[1], col);
            }
          }

          if (Dim() == 3)
          {
            // 3D loop over all entries of the current derivative map (teta)
            for (colcurr = dtetamap[j].begin(); colcurr != dtetamap[j].end(); ++colcurr)
            {
              int col = colcurr->first;
              double valtxi = (zteta + ct * jumpteta) / euclidean * (colcurr->second) * z[j] * ztxi;
              double valteta =
                  (zteta + ct * jumpteta) / euclidean * (colcurr->second) * z[j] * zteta;

              // do not assemble zeros into matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(valtxi * txi[j]) > 1.e-12)
                    linslipDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
              }
              else if (abs(valtxi) > 1.0e-12)
                linslipDISglobal.Assemble(valtxi, row[0], col);

              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(valteta * teta[j]) > 1.e-12)
                    linslipDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(valteta) > 1.0e-12)
                linslipDISglobal.Assemble(valteta, row[1], col);
            }
          }
        }

        // loop over dimensions
        for (int j = 0; j < Dim(); ++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (colcurr = dtximap[j].begin(); colcurr != dtximap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
            double valtxi =
                (ztxi + ct * jumptxi) / euclidean * ct * (colcurr->second) * jump[j] * ztxi;
            double valteta =
                (ztxi + ct * jumptxi) / euclidean * ct * (colcurr->second) * jump[j] * zteta;

            // do not assemble zeros into s matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(valtxi * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
              for (int j = 0; j < Dim(); j++)
                if (abs(valteta * teta[j]) > 1.e-12)
                  linslipDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
            }
            else
            {
              if (abs(valtxi) > 1.0e-12) linslipDISglobal.Assemble(valtxi, row[0], col);
              if (abs(valteta) > 1.0e-12) linslipDISglobal.Assemble(valteta, row[1], col);
            }
          }

          if (Dim() == 3)
          {
            // loop over all entries of the current derivative map (teta)
            for (colcurr = dtetamap[j].begin(); colcurr != dtetamap[j].end(); ++colcurr)
            {
              int col = colcurr->first;
              double valtxi =
                  (zteta + ct * jumpteta) / euclidean * ct * (colcurr->second) * jump[j] * ztxi;
              double valteta =
                  (zteta + ct * jumpteta) / euclidean * ct * (colcurr->second) * jump[j] * zteta;

              // do not assemble zeros into matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(valtxi * txi[j]) > 1.e-12)
                    linslipDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
                for (int j = 0; j < Dim(); j++)
                  if (abs(valteta * teta[j]) > 1.e-12)
                    linslipDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
              }
              else
              {
                if (abs(valtxi) > 1.0e-12) linslipDISglobal.Assemble(valtxi, row[0], col);
                if (abs(valteta) > 1.0e-12) linslipDISglobal.Assemble(valteta, row[1], col);
              }
            }
          }
        }
        /*** 4 ************************** (frcoeff*znor).deriv(T).z ***/
        // loop over all dimensions
        for (int j = 0; j < Dim(); ++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (colcurr = dtximap[j].begin(); colcurr != dtximap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
            double val = (-1) * (frbound) * (colcurr->second) * z[j];

            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(val * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(val) > 1.0e-12)
              linslipDISglobal.Assemble(val, row[0], col);
          }

          if (Dim() == 3)
          {
            // loop over all entries of the current derivative map (teta)
            for (colcurr = dtetamap[j].begin(); colcurr != dtetamap[j].end(); ++colcurr)
            {
              int col = colcurr->first;
              double val = (-1.0) * (frbound) * (colcurr->second) * z[j];

              // do not assemble zeros into matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(val * teta[j]) > 1.e-12)
                    linslipDISglobal.Assemble(val * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(val) > 1.0e-12)
                linslipDISglobal.Assemble(val, row[1], col);
            }
          }
        }

        /*** 5 *********************** (frcoeff*znor).deriv(T).jump ***/
        // loop over all dimensions
        for (int j = 0; j < Dim(); ++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (colcurr = dtximap[j].begin(); colcurr != dtximap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
            double val = (-1) * (frbound)*ct * (colcurr->second) * jump[j];

            // do not assemble zeros into matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(val * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(val * txi[j], cnode->Dofs()[j], col);
            }
            else if (abs(val) > 1.0e-12)
              linslipDISglobal.Assemble(val, row[0], col);
          }

          if (Dim() == 3)
          {
            // loop over all entries of the current derivative map (teta)
            for (colcurr = dtetamap[j].begin(); colcurr != dtetamap[j].end(); ++colcurr)
            {
              int col = colcurr->first;
              double val = (-1) * (frbound)*ct * (colcurr->second) * jump[j];

              // do not assemble zeros into s matrix
              if (constr_direction_ == INPAR::CONTACT::constr_xyz)
              {
                for (int j = 0; j < Dim(); j++)
                  if (abs(val * teta[j]) > 1.e-12)
                    linslipDISglobal.Assemble(val * teta[j], cnode->Dofs()[j], col);
              }
              else if (abs(val) > 1.0e-12)
                linslipDISglobal.Assemble(val, row[1], col);
            }
          }
        }

        /*** 6 ******************* -frcoeff.Deriv(n).z(ztan+ct*utan) ***/
        // loop over all dimensions
        for (int j = 0; j < Dim(); ++j)
        {
          // loop over all entries of the current derivative map
          for (colcurr = dnmap[j].begin(); colcurr != dnmap[j].end(); ++colcurr)
          {
            int col = colcurr->first;
            double valtxi = (-1) * (ztxi + ct * jumptxi) * frcoeff * (1. - cn * dval * d_gLM) *
                            (colcurr->second) * z[j];
            double valteta = (-1) * (zteta + ct * jumpteta) * frcoeff * (1. - cn * dval * d_gLM) *
                             (colcurr->second) * z[j];

            // do not assemble zeros into s matrix
            if (constr_direction_ == INPAR::CONTACT::constr_xyz)
            {
              for (int j = 0; j < Dim(); j++)
                if (abs(valtxi * txi[j]) > 1.e-12)
                  linslipDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
              for (int j = 0; j < Dim(); j++)
                if (abs(valteta * teta[j]) > 1.e-12)
                  linslipDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
            }
            else
            {
              if (abs(valtxi) > 1.0e-12) linslipDISglobal.Assemble(valtxi, row[0], col);
              if (abs(valteta) > 1.0e-12) linslipDISglobal.Assemble(valteta, row[1], col);
            }
          }
        }

        /*** 7 ****************** deriv (g) ***/
        // prepare assembly
        std::map<int, double>& dgmap = cnode->CoData().GetDerivG();

        // loop over all entries of the current derivative map
        for (_colcurr = dgmap.begin(); _colcurr != dgmap.end(); ++_colcurr)
        {
          int col = _colcurr->first;
          double valtxi = frcoeff * cn * (_colcurr->second) * (ztxi + ct * jumptxi);
          double valteta = frcoeff * cn * (_colcurr->second) * (zteta + ct * jumpteta);

          // do not assemble zeros into matrix
          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int j = 0; j < Dim(); j++)
              if (abs(valtxi * txi[j]) > 1.e-12)
                linslipDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
            for (int j = 0; j < Dim(); j++)
              if (abs(valteta * teta[j]) > 1.e-12)
                linslipDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
          }
          else
          {
            if (abs(valtxi) > 1.0e-12) linslipDISglobal.Assemble(valtxi, row[0], col);
            if (abs(valteta) > 1.0e-12) linslipDISglobal.Assemble(valteta, row[1], col);
          }
        }

        std::map<int, double>& dd = cnode->CoData().GetDerivD()[cnode->Id()];
        for (auto p = dd.begin(); p != dd.end(); ++p)
        {
          int col = p->first;
          double valtxi = frcoeff * cn * gLM * p->second * (ztxi + ct * jumptxi);
          double valteta = frcoeff * cn * gLM * p->second * (zteta + ct * jumpteta);

          // do not assemble zeros into matrix
          if (constr_direction_ == INPAR::CONTACT::constr_xyz)
          {
            for (int j = 0; j < Dim(); j++)
              if (abs(valtxi * txi[j]) > 1.e-12)
                linslipDISglobal.Assemble(valtxi * txi[j], cnode->Dofs()[j], col);
            for (int j = 0; j < Dim(); j++)
              if (abs(valteta * teta[j]) > 1.e-12)
                linslipDISglobal.Assemble(valteta * teta[j], cnode->Dofs()[j], col);
          }
          else
          {
            if (abs(valtxi) > 1.0e-12) linslipDISglobal.Assemble(valtxi, row[0], col);
            if (abs(valteta) > 1.0e-12) linslipDISglobal.Assemble(valteta, row[1], col);
          }
        }
      }  // if (frcoeff==0.0)
    }    // loop over all slip nodes of the interface
  }      // Coulomb friction

  else
  {
    dserror("no tresca friction with regularized normal constraint");
  }  // Tresca friction

  return;
}

/*---------------------------------------------------------------------------*
 |  Assemble normal coupling weighted condition for poro contact   ager 07/14|
 *--------------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleNCoup(Epetra_Vector& gglobal)
{
  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i = 0; i < activenodes_->NumMyElements(); ++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    CoNode* mrtnode = dynamic_cast<CoNode*>(node);

    if (mrtnode->Owner() != Comm().MyPID()) dserror("AssembleDMG: Node ownership inconsistency!");

    /**************************************************** nCoup-vector ******/
    if (mrtnode->CoPoroData().GetnCoup() != 0.0)
    {
      double nCoup = mrtnode->CoPoroData().GetnCoup();

      CORE::LINALG::SerialDenseVector gnode(1);
      std::vector<int> lm(1);
      std::vector<int> lmowner(1);

      gnode(0) = nCoup;
      lm[0] = activen_->GID(i);

      lmowner[0] = mrtnode->Owner();

      CORE::LINALG::Assemble(gglobal, gnode, lm, lmowner);
    }
  }

  return;
}

/*---------------------------------------------------------------------*
 |  Assemble linearisation of normal coupling                          |
 |          weighted condition for poro contact              ager 07/14|
 *--------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleNCoupLin(
    CORE::LINALG::SparseMatrix& sglobal, CORE::ADAPTER::Coupling& coupfs, bool AssembleVelocityLin)
{
  // nothing to do if no active nodes
  if (activenodes_ == Teuchos::null) return;

  Teuchos::RCP<const Epetra_Map> MasterDofMap_full;
  Teuchos::RCP<const Epetra_Map> PermSlaveDofMap_full;

  if (AssembleVelocityLin)
  {
    // store map on all processors, simple but expensive
    MasterDofMap_full = CORE::LINALG::AllreduceEMap(*coupfs.MasterDofMap());
    PermSlaveDofMap_full = CORE::LINALG::AllreduceEMap(*coupfs.PermSlaveDofMap());
  }

  for (int i = 0; i < activenodes_->NumMyElements(); ++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID()) dserror("AssembleS: Node ownership inconsistency!");

    std::map<int, double>::iterator colcurr;
    int row = activen_->GID(i);

    std::map<int, double>& dgmap = cnode->CoPoroData().GetDerivnCoup();
    if (AssembleVelocityLin)
    {  // Assign fluid velocity linearization to matrix

      dgmap = cnode->CoPoroData().GetVelDerivnCoup();
      for (colcurr = dgmap.begin(); colcurr != dgmap.end(); ++colcurr)
      {
        int col = PermSlaveDofMap_full->GID(MasterDofMap_full->LID(colcurr->first));
        double val = colcurr->second;

        // do not assemble zeros into s matrix
        if (abs(val) > 1.0e-12) sglobal.Assemble(val, row, col);
      }
      // Assign pressure linearization to matrix
      dgmap = cnode->CoPoroData().GetPresDerivnCoup();
      for (colcurr = dgmap.begin(); colcurr != dgmap.end(); ++colcurr)
      {
        int col = PermSlaveDofMap_full->GID(MasterDofMap_full->LID(colcurr->first)) + Dim();
        double val = colcurr->second;

        // do not assemble zeros into s matrix
        if (abs(val) > 1.0e-12) sglobal.Assemble(val, row, col);
      }
    }
    else
    {  // Assign skeleton displacement linearization to matrix
      for (colcurr = dgmap.begin(); colcurr != dgmap.end(); ++colcurr)
      {
        int col = colcurr->first;
        double val = colcurr->second;

        // do not assemble zeros into s matrix
        if (abs(val) > 1.0e-12) sglobal.Assemble(val, row, col);
      }
    }
  }
  return;
}

/*------------------------------------------------------------------------*
 | Derivative of D-matrix multiplied with a slave dof vector              |
 *------------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleCoupLinD(
    CORE::LINALG::SparseMatrix& CoupLin, const Teuchos::RCP<Epetra_Vector> x)
{
  // we have: D_jk,c with j = Slave dof
  //                 with k = Displacement slave dof
  //                 with c = Displacement slave or master dof
  // we compute (LinD)_kc = D_jk,c * x_j

  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    // Mortar matrix D derivatives
    std::map<int, std::map<int, double>>& dderiv = cnode->CoData().GetDerivD();

    // get sizes and iterator start
    int slavesize = (int)dderiv.size();
    std::map<int, std::map<int, double>>::iterator scurr = dderiv.begin();

    /********************************************** LinDMatrix **********/
    // loop over all DISP slave nodes in the DerivD-map of the current slave node
    for (int k = 0; k < slavesize; ++k)
    {
      int sgid = scurr->first;
      ++scurr;

      DRT::Node* snode = idiscret_->gNode(sgid);
      if (!snode) dserror("Cannot find node with gid %", sgid);
      CoNode* csnode = dynamic_cast<CoNode*>(snode);  // current slave node

      // Mortar matrix D derivatives
      std::map<int, double>& thisdderiv = cnode->CoData().GetDerivD()[sgid];
      int mapsize = (int)(thisdderiv.size());

      if (cnode->NumDof() != csnode->NumDof())
        dserror("Mortar Nodes on interface must have same number of dofs!");

      // inner product D_{jk,c} * z_j for index j
      for (int prodj = 0; prodj < cnode->NumDof(); ++prodj)
      {
        int row = csnode->Dofs()[prodj];
        std::map<int, double>::iterator scolcurr = thisdderiv.begin();

        // loop over all directional derivative entries
        for (int c = 0; c < mapsize; ++c)
        {
          int col = scolcurr->first;
          int slavedofgid = cnode->Dofs()[prodj];

          int slavedoflid = x->Map().LID(slavedofgid);
          if (slavedoflid < 0) dserror("invalid slave dof lid");
          double val = (*x)[slavedoflid] * (scolcurr->second);

          ++scolcurr;


          //************   ASSEMBLY INTO THE MATRIX    **********************
          if (abs(val) > 1.0e-12) CoupLin.FEAssemble(val, row, col);
        }
        // check for completeness of DerivD-Derivatives-iteration
        if (scolcurr != thisdderiv.end())
          dserror("AssembleCoupLin: Not all derivative entries of Lin(D*z_s) considered!");
      }
    }

    if (scurr != dderiv.end())
      dserror("AssembleCoupLin: Not all connected slave nodes to Lin(D*z_s) considered!");
  }

  return;
}

/*-----------------------------------------------------------------------------------*
 | Derivative of transposed M-matrix multiplied with a slave dof vector  seitz 01/18 |
 *-----------------------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleCoupLinM(
    CORE::LINALG::SparseMatrix& CoupLin, const Teuchos::RCP<Epetra_Vector> x)
{
  // we have: M_jl,c with j = Slave dof
  //                 with l = Displacement master dof
  //                 with c = Displacement slave or master dof
  // we compute (CoupLin)_lc = M_jl,c * x_j

  // loop over all slave nodes (row map)
  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    // Mortar matrix M derivatives
    std::map<int, std::map<int, double>>& mderiv = cnode->CoData().GetDerivM();

    // get sizes and iterator start
    int mastersize = (int)mderiv.size();
    std::map<int, std::map<int, double>>::iterator mcurr = mderiv.begin();

    /********************************************** LinMMatrix **********/
    // loop over all master nodes in the DerivM-map of the current LM slave node
    for (int l = 0; l < mastersize; ++l)
    {
      int mgid = mcurr->first;
      ++mcurr;

      DRT::Node* mnode = idiscret_->gNode(mgid);
      if (!mnode) dserror("Cannot find node with gid %", mgid);
      CoNode* cmnode = dynamic_cast<CoNode*>(mnode);

      // Mortar matrix M derivatives
      std::map<int, double>& thismderiv = cnode->CoData().GetDerivM()[mgid];
      int mapsize = (int)(thismderiv.size());

      if (cnode->NumDof() != cmnode->NumDof())
        dserror("Mortar Nodes on interface must have same number of dofs!");

      // inner product M_{jl,c} * z_j for index j
      for (int prodj = 0; prodj < cmnode->NumDof(); ++prodj)
      {
        int row = cmnode->Dofs()[prodj];
        std::map<int, double>::iterator mcolcurr = thismderiv.begin();

        // loop over all directional derivative entries
        for (int c = 0; c < mapsize; ++c)
        {
          int col = mcolcurr->first;
          int slavedofgid = cnode->Dofs()[prodj];

          int slavedoflid = x->Map().LID(slavedofgid);
          double val = (*x)[slavedoflid] * (mcolcurr->second);

          ++mcolcurr;

          ///************   ASSEMBLY INTO THE MATRIX    **********************
          if (abs(val) > 1.0e-12) CoupLin.FEAssemble(val, row, col);
        }

        // check for completeness of DerivM-Derivatives-iteration
        if (mcolcurr != thismderiv.end())
          dserror("AssembleCoupLin: Not all derivative entries of DerivM considered!");
      }
    }  // loop over all master nodes, connected to this slave node => Lin(M*z_s) finished

    // check for completeness of DerivM-Master-iteration
    if (mcurr != mderiv.end())
      dserror("AssembleCoupLin: Not all master entries of Lin(M*z_s) considered!");

    /*********************************************************************/
    /*******************       Finished with Lin(M*z_s)         **********/
    /*********************************************************************/
  }
}
BACI_NAMESPACE_CLOSE
