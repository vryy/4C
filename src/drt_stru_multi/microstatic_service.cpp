/*---------------------------------------------------------------------*/
/*!

\brief Service for microstructural problems

\maintainer Martin Kronbichler

\level 2

*/
/*---------------------------------------------------------------------*/


#include "microstatic.H"
#include "../drt_structure/stru_aux.H"
#include "../linalg/linalg_utils.H"


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
void STRUMULTI::MicroStatic::DetermineToggle()
{
  int np = 0;  // local number of prescribed (=boundary) dofs needed for the
               // creation of vectors and matrices for homogenization
               // procedure

  std::vector<DRT::Condition*> conds;
  discret_->GetCondition("MicroBoundary", conds);
  for (unsigned i = 0; i < conds.size(); ++i)
  {
    const std::vector<int>* nodeids = conds[i]->Get<std::vector<int>>("Node Ids");
    if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
    const int nnode = (*nodeids).size();

    for (int i = 0; i < nnode; ++i)
    {
      // do only nodes in my row map
      if (!discret_->NodeRowMap()->MyGID((*nodeids)[i])) continue;
      DRT::Node* actnode = discret_->gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node %d", (*nodeids)[i]);
      std::vector<int> dofs = discret_->Dof(actnode);
      const unsigned numdf = dofs.size();

      for (unsigned j = 0; j < numdf; ++j)
      {
        const int gid = dofs[j];

        const int lid = disn_->Map().LID(gid);
        if (lid < 0) dserror("Global id %d not on this proc in system vector", gid);

        if ((*dirichtoggle_)[lid] != 1.0)  // be careful not to count dofs more
                                           // than once since nodes belong to
                                           // several surfaces simultaneously
          ++np;

        (*dirichtoggle_)[lid] = 1.0;
      }
    }
  }

  np_ = np;
}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
void STRUMULTI::MicroStatic::SetUpHomogenization()
{
  int indp = 0;
  int indf = 0;

  ndof_ = discret_->DofRowMap()->NumMyElements();

  std::vector<int> pdof(np_);
  std::vector<int> fdof(ndof_ - np_);  // changed this, previously this
                                       // has been just fdof(np_),
                                       // but how should that
                                       // work for ndof_-np_>np_???

  for (int it = 0; it < ndof_; ++it)
  {
    if ((*dirichtoggle_)[it] == 1.0)
    {
      pdof[indp] = discret_->DofRowMap()->GID(it);
      ++indp;
    }
    else
    {
      fdof[indf] = discret_->DofRowMap()->GID(it);
      ++indf;
    }
  }

  // create map based on the determined dofs of prescribed and free nodes
  pdof_ = Teuchos::rcp(new Epetra_Map(-1, np_, &pdof[0], 0, discret_->Comm()));
  fdof_ = Teuchos::rcp(new Epetra_Map(-1, ndof_ - np_, &fdof[0], 0, discret_->Comm()));

  // create importer
  importp_ = Teuchos::rcp(new Epetra_Import(*pdof_, *(discret_->DofRowMap())));
  importf_ = Teuchos::rcp(new Epetra_Import(*fdof_, *(discret_->DofRowMap())));

  // create vector containing material coordinates of prescribed nodes
  Epetra_Vector Xp_temp(*pdof_);

  std::vector<DRT::Condition*> conds;
  discret_->GetCondition("MicroBoundary", conds);
  for (unsigned i = 0; i < conds.size(); ++i)
  {
    const std::vector<int>* nodeids = conds[i]->Get<std::vector<int>>("Node Ids");
    if (!nodeids) dserror("MicroBoundary condition does not have nodal cloud");
    const int nnode = (*nodeids).size();

    for (int i = 0; i < nnode; ++i)
    {
      // do only nodes in my row map
      if (!discret_->NodeRowMap()->MyGID((*nodeids)[i])) continue;
      DRT::Node* actnode = discret_->gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node %d", (*nodeids)[i]);

      // nodal coordinates
      const double* x = actnode->X();

      std::vector<int> dofs = discret_->Dof(actnode);

      for (int k = 0; k < 3; ++k)
      {
        const int gid = dofs[k];

        const int lid = disn_->Map().LID(gid);
        if (lid < 0) dserror("Global id %d not on this proc in system vector", gid);

        for (int l = 0; l < np_; ++l)
        {
          if (pdof[l] == gid) Xp_temp[l] = x[k];
        }
      }
    }
  }

  Xp_ = LINALG::CreateVector(*pdof_, true);
  *Xp_ = Xp_temp;

  // now create D and its transpose DT (following Miehe et al., 2002)
  // NOTE: D_ has the same row GIDs (0-8), but different col IDs on different procs (corresponding
  // to pdof_).
  D_ = Teuchos::rcp(new Epetra_MultiVector(*pdof_, 9));

  for (int n = 0; n < np_ / 3; ++n)
  {
    (*(*D_)(0))[3 * n] = (*Xp_)[3 * n];
    (*(*D_)(3))[3 * n] = (*Xp_)[3 * n + 1];
    (*(*D_)(6))[3 * n] = (*Xp_)[3 * n + 2];

    (*(*D_)(1))[3 * n + 1] = (*Xp_)[3 * n + 1];
    (*(*D_)(4))[3 * n + 1] = (*Xp_)[3 * n + 2];
    (*(*D_)(7))[3 * n + 1] = (*Xp_)[3 * n];

    (*(*D_)(2))[3 * n + 2] = (*Xp_)[3 * n + 2];
    (*(*D_)(5))[3 * n + 2] = (*Xp_)[3 * n];
    (*(*D_)(8))[3 * n + 2] = (*Xp_)[3 * n + 1];
  }

  Epetra_MultiVector DT(*pdof_, 9);

  for (int n = 0; n < np_ / 3; ++n)
  {
    (*(DT(0)))[3 * n] = (*Xp_)[3 * n];
    (*(DT(1)))[3 * n + 1] = (*Xp_)[3 * n + 1];
    (*(DT(2)))[3 * n + 2] = (*Xp_)[3 * n + 2];
    (*(DT(3)))[3 * n] = (*Xp_)[3 * n + 1];
    (*(DT(4)))[3 * n + 1] = (*Xp_)[3 * n + 2];
    (*(DT(5)))[3 * n + 2] = (*Xp_)[3 * n];
    (*(DT(6)))[3 * n] = (*Xp_)[3 * n + 2];
    (*(DT(7)))[3 * n + 1] = (*Xp_)[3 * n];
    (*(DT(8)))[3 * n + 2] = (*Xp_)[3 * n + 1];
  }

  rhs_ = Teuchos::rcp(new Epetra_MultiVector(*(discret_->DofRowMap()), 9));

  for (int i = 0; i < 9; ++i)
  {
    ((*rhs_)(i))->Export(*(DT(i)), *importp_, Insert);
  }
}


/*----------------------------------------------------------------------*
 |  check convergence of Newton iteration (public)              lw 12/07|
 *----------------------------------------------------------------------*/
bool STRUMULTI::MicroStatic::Converged()
{
  // check for single norms
  bool convdis = false;
  bool convfres = false;

  // residual displacement
  switch (normtypedisi_)
  {
    case INPAR::STR::convnorm_abs:
      convdis = normdisi_ < toldisi_;
      break;
    case INPAR::STR::convnorm_rel:
      convdis = normdisi_ / normchardis_ < toldisi_;
      break;
    case INPAR::STR::convnorm_mix:
      convdis = ((normdisi_ < toldisi_) or (normdisi_ / normchardis_ < toldisi_));
      break;
    default:
      dserror("Cannot check for convergence of residual displacements!");
      break;
  }

  // residual forces
  switch (normtypefres_)
  {
    case INPAR::STR::convnorm_abs:
      convfres = normfres_ < tolfres_;
      break;
    case INPAR::STR::convnorm_rel:
      convfres = normfres_ / normcharforce_ < tolfres_;
      break;
    case INPAR::STR::convnorm_mix:
      convfres = ((normfres_ < tolfres_) or (normfres_ / normcharforce_ < tolfres_));
      break;
    default:
      dserror("Cannot check for convergence of residual forces!");
      break;
  }

  // combine displacement-like and force-like residuals
  bool conv = false;
  if (combdisifres_ == INPAR::STR::bop_and)
    conv = convdis and convfres;
  else if (combdisifres_ == INPAR::STR::bop_or)
    conv = convdis or convfres;
  else
    dserror("Something went terribly wrong with binary operator!");

  return conv;
}


/*----------------------------------------------------------------------*
 |  calculate reference norms for relative convergence checks   lw 12/07|
 *----------------------------------------------------------------------*/
void STRUMULTI::MicroStatic::CalcRefNorms()
{
  // The reference norms are used to scale the calculated iterative
  // displacement norm and/or the residual force norm. For this
  // purpose we only need the right order of magnitude, so we don't
  // mind evaluating the corresponding norms at possibly different
  // points within the timestep (end point, generalized midpoint).
  // In the beginning (construction of macroscale time integrator and
  // first macroscale predictor), macro displacements are generally
  // 0 leading to no load on the micro problem. Consequently, the
  // microscale reference norms are 0 in case of displacements, and
  // near 0 in case of the residual (sum of ndof numerical near zero
  // values). To enable convergence in these cases, the reference norm
  // is automatically set to 1 if the calculated values are below
  // the chosen tolerances. Simply testing against 0 only works for
  // the displacements, but not for the residual!

  normchardis_ = STR::AUX::CalculateVectorNorm(iternorm_, dis_);
  if (normchardis_ < toldisi_) normchardis_ = 1.0;

  double fintnorm = STR::AUX::CalculateVectorNorm(iternorm_, fintn_);
  double freactnorm = STR::AUX::CalculateVectorNorm(iternorm_, freactn_);
  normcharforce_ = std::max(fintnorm, freactnorm);
  if (normcharforce_ < tolfres_) normcharforce_ = 1.0;
}


/*----------------------------------------------------------------------*
 |  print to screen and/or error file                           lw 12/07|
 *----------------------------------------------------------------------*/
void STRUMULTI::MicroStatic::PrintNewton(bool print_unconv, Epetra_Time timer)
{
  bool relres = (normtypefres_ == INPAR::STR::convnorm_rel);

  bool relres_reldis = ((normtypedisi_ == INPAR::STR::convnorm_rel) && relres);

  if (relres)
  {
    normfres_ /= normcharforce_;
  }

  if (relres_reldis)
  {
    normfres_ /= normcharforce_;
    normdisi_ /= normchardis_;
  }

  if (print_unconv)
  {
    if (printscreen_ and (stepn_ % printscreen_ == 0))
    {
      if (relres)
      {
        printf("      MICROSCALE numiter %2d scaled res-norm %10.5e absolute dis-norm %20.15E\n",
            numiter_ + 1, normfres_, normdisi_);
        fflush(stdout);
      }
      else if (relres_reldis)
      {
        printf("      MICROSCALE numiter %2d scaled res-norm %10.5e scaled dis-norm %20.15E\n",
            numiter_ + 1, normfres_, normdisi_);
        fflush(stdout);
      }
      else
      {
        printf("      MICROSCALE numiter %2d absolute res-norm %10.5e absolute dis-norm %20.15E\n",
            numiter_ + 1, normfres_, normdisi_);
        fflush(stdout);
      }
    }
  }
  else
  {
    double timepernlnsolve = timer.ElapsedTime();

    if (relres)
    {
      printf(
          "      MICROSCALE Newton iteration converged: numiter %d scaled res-norm %e absolute "
          "dis-norm %e time %10.5f\n\n",
          numiter_, normfres_, normdisi_, timepernlnsolve);
      fflush(stdout);
    }
    else if (relres_reldis)
    {
      printf(
          "      MICROSCALE Newton iteration converged: numiter %d scaled res-norm %e scaled "
          "dis-norm %e time %10.5f\n\n",
          numiter_, normfres_, normdisi_, timepernlnsolve);
      fflush(stdout);
    }
    else
    {
      printf(
          "      MICROSCALE Newton iteration converged: numiter %d absolute res-norm %e absolute "
          "dis-norm %e time %10.5f\n\n",
          numiter_, normfres_, normdisi_, timepernlnsolve);
      fflush(stdout);
    }
  }
}


/*----------------------------------------------------------------------*
 |  print to screen                                             lw 12/07|
 *----------------------------------------------------------------------*/
void STRUMULTI::MicroStatic::PrintPredictor()
{
  if (normtypefres_ == INPAR::STR::convnorm_rel)
  {
    normfres_ /= normcharforce_;
    std::cout << "      MICROSCALE Predictor scaled res-norm " << normfres_ << std::endl;
  }
  else
  {
    std::cout << "      MICROSCALE Predictor absolute res-norm " << normfres_ << std::endl;
  }
  fflush(stdout);
}
