/*!----------------------------------------------------------------------
\file contact_augmented_strategy_tools.cpp

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Created on: Jun 3, 2014

Maintainer: Michael Hiermeier
            hiermeier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15268
</pre>

*----------------------------------------------------------------------*/

#include "contact_augmented_strategy.H"
#include "contact_augmented_interface.H"
#include "../drt_contact/contact_defines.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*
 | Finite Difference check at GP level                   hiermeier 05/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::AugFDCheckGP()
{
#ifdef CONTACTFD_KAPPA
  // FD check of the linearization of the scaling parameter kappa
  for (int i=0;i<(int)interface_.size();++i)
    interface_[i]->FDCheckKappaLin();
#endif

#ifdef CONTACTFD_AWGAP
  // FD check of the linearization
  for (int i=0;i<(int)interface_.size();++i)
    interface_[i]->FDCheckAWGapLin();
#endif

#ifdef CONTACTFD_VARWGAPSL
  // FD check of the linearization
  // of the variation of WGap (slave)
  for (int i=0;i<(int)interface_.size();++i)
    interface_[i]->FDCheckVarWGapLinSl();
#endif

#ifdef CONTACTFD_VARWGAPMA
  // FD check of the linearization
  // of the variation of WGap (master)
  for (int i=0;i<(int)interface_.size();++i)
    interface_[i]->FDCheckVarWGapLinMa();
#endif

#ifdef CONTACTFD_AUGA
  // FD check of the linearization of the scaling parameter kappa
  for (int i=0;i<(int)interface_.size();++i)
    interface_[i]->FDCheckAugALin();
#endif

  return;
}

/*----------------------------------------------------------------------*
 | Finite Difference check at global level               hiermeier 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::AugFDCheckGlobal()
{
  // *** linearization w.r.t. displ. *********************************
  /*---------------------------------*
   | Active normal constraint        |
   *---------------------------------*/
#ifdef CONTACTFD_DLMNWGAPLINMATRIX
  if (!gFdCheck_) EvalFDCheckGlobalDispl(*dLmNWGapLinMatrix_,dLmNWGapRhs_);
#endif
  /*---------------------------------*
   | Active tangential constraint    |
   *---------------------------------*/
#ifdef CONTACTFD_DLMTLMTLINMATRIX
  if (!gFdCheck_) EvalFDCheckGlobalDispl(*dLmTLmTLinMatrix_,dLmTLmTRhs_);
#endif
  /*---------------------------------*
   | Active tangential constraint    |
   *---------------------------------*/
#ifdef CONTACTFD_AUGINACTIVELINMATRIX
  if (!gFdCheck_) EvalFDCheckGlobalDispl(*augInactiveLinMatrix_,augInactiveRhs_);
#endif
  /*---------------------------------*
   | Std force terms                 |
   *---------------------------------*/
#if defined(CONTACTFD_DGLMSLLINMATRIX) || defined(CONTACTFD_DGLMMALINMATRIX)
  Teuchos::RCP<Epetra_Vector>  z_exp = Teuchos::rcp(new Epetra_Vector(*gAugActiveSlaveNDofs_));
  LINALG::Export(*z_,*z_exp);
  // *** slave side ***
#ifdef CONTACTFD_DGLMSLLINMATRIX
  augDnMatrix_->Multiply(true,*z_exp,*augfs_lm_);
  if (!gFdCheck_) EvalFDCheckGlobalDispl(*dGLmSlLinMatrix_,augfs_lm_);
#endif
  // *** master side ***
#ifdef CONTACTFD_DGLMMALINMATRIX
  augMnMatrix_->Multiply(true,*z_exp,*augfm_lm_);
  if (!gFdCheck_) EvalFDCheckGlobalDispl(*dGLmMaLinMatrix_,augfm_lm_);
#endif
#endif
  /*---------------------------------*
   | Regularization term             |
   *---------------------------------*/
#if defined(CONTACTFD_DGGSLLINMATRIX) || defined(CONTACTFD_DGGMALINMATRIX)
  double cn= Params().get<double>("SEMI_SMOOTH_CN");
  // *** slave side ***
#ifdef CONTACTFD_DGGSLLINMATRIX
  augDnMatrix_->Multiply(true,*aWGapRhs_,*augfs_g_);
  augfs_g_->Scale(-cn);
  if (!gFdCheck_) EvalFDCheckGlobalDispl(*dGGSlLinMatrix_,augfs_g_);
#endif
  // *** master side ***
#ifdef CONTACTFD_DGGMALINMATRIX
  augMnMatrix_->Multiply(true,*aWGapRhs_,*augfm_g_);
  augfm_g_->Scale(-cn);
  if (!gFdCheck_) EvalFDCheckGlobalDispl(*dGGMaLinMatrix_,augfm_g_);
#endif
#endif

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for global matrices           hiermeier 06/14|
 | w.r.t. to the displacements                                          |
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::EvalFDCheckGlobalDispl(LINALG::SparseMatrix derivMatrix,
                                                            Teuchos::RCP<Epetra_Vector>& rhsVector)
{
  // activate global finite difference indicator
  gFdCheck_ = true;

  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");
//  if ((int) interface_.size()>1) dserror("ERROR: FD check only for ONE interface");

  const Epetra_Map rowMap = derivMatrix.RowMap();
  const Epetra_Map colMap = derivMatrix.ColMap();

  const Epetra_BlockMap rhsMap = rhsVector->Map();

  LINALG::SparseMatrix fdMatrixRef = LINALG::SparseMatrix(rowMap,100);
  LINALG::SparseMatrix fdMatrixNew = LINALG::SparseMatrix(rowMap,100);
  int dim = Dim();

  // create reference Matrix:
  // loop over all columns
  for (int c=0;c<colMap.NumMyElements();++c)
  {
    int colId = colMap.GID(c);
    // loop over all rows of the right hand side vector
    for (int r=0;r<rhsMap.NumMyElements();++r)
    {
      int rowId = rhsMap.GID(r);
      if (rowMap.LID(rowId)==-1)
        dserror("ERROR: Row gids of the corresponding rhs-vector and the derivative matrix do not match!");

      double val = (*rhsVector)[r];

      fdMatrixRef.Assemble(val,rowId,colId);
    }
  }
  fdMatrixRef.Complete(colMap,rowMap);

  double delta = 1.0e-8;
  // loop over all columns of the reference Matrix
  for (int fd=0;fd<colMap.NumMyElements();++fd)
  {
    int colId = colMap.GID(fd);
    int gid = colId/dim;
    int dof = colId%dim;
    bool iCheck = false;

    // do the finite difference step
    for (int i=0;i<(int) interface_.size();++i)
    {
      iCheck = interface_[i]->UpdateInterfaces(gid,dof,delta,true);
      if (iCheck) break;
    }
    if (!iCheck) dserror("ERROR: Node % not found!",gid);

    // Update matrix and rhs
    Initialize();
    Teuchos::RCP<LINALG::SparseOperator> kteff = Teuchos::null;
    Teuchos::RCP<Epetra_Vector> feff = Teuchos::null;
    EvaluateContact(kteff,feff);

    // loop over all rows of the updated right hand side vector
    // and save the values in a new matrix
    for (int r=0; r<rhsMap.NumMyElements();++r)
    {
      int rowId = rhsMap.GID(r);
      double val = (*rhsVector)[r];
      fdMatrixNew.Assemble(val,rowId,colId);
    }

    // Undo finite difference step
    iCheck = false;
    for (int i=0;i<(int) interface_.size();++i)
    {
      iCheck = interface_[i]->UpdateInterfaces(gid,dof,delta,false);
      if (iCheck) break;
    }
    if (!iCheck) dserror("ERROR: Node % not found!",gid);
    Initialize();
    EvaluateContact(kteff,feff);
  }

  // calculate the finite difference
  fdMatrixNew.Add(fdMatrixRef,false,-1.0/delta,1.0/delta);
  fdMatrixNew.Complete(colMap,rowMap);

  // loop over all rows
  for (int r=0;r<rowMap.NumMyElements();++r)
  {
    int rowId = rowMap.GID(r);
    // check if the row belongs to the slave or master side
    std::string rSlMa = "(S)";
    if (gsdofrowmap_->LID(rowId)==-1)
      rSlMa = "(M)";

    int w = 0;

    // *** finite differences ***
    // get all non-zero values and the corresponding ids of the current row
    int rLengthFD = fdMatrixNew.EpetraMatrix()->NumGlobalEntries(rowId);
    int numEntriesFD = 0;
    std::vector<double> rValFD(rLengthFD);
    std::vector<int>    cIdsFD(rLengthFD);
    fdMatrixNew.EpetraMatrix()->ExtractGlobalRowCopy(rowId,rLengthFD,numEntriesFD,&rValFD[0],&cIdsFD[0]);

    // *** analytical solution ***
    // get all non-zero values and the corresponding ids of the current row
    int rLengthAna = derivMatrix.EpetraMatrix()->NumGlobalEntries(rowId);
    int numEntriesAna = 0;
    std::vector<double> rValAna(rLengthAna);
    std::vector<int>    cIdsAna(rLengthAna);
    derivMatrix.EpetraMatrix()->ExtractGlobalRowCopy(rowId,rLengthAna,numEntriesAna,&rValAna[0],&cIdsAna[0]);

    /*-------------------------------------------------------------*
     |   Compare analytical and finite difference solution         |
     *-------------------------------------------------------------*/
    std::cout << "\nFINITE DIFFERENCE CHECK FOR ROW-ID # " << rowId << rSlMa << std::endl;
    for (int c=0; c<rLengthFD;++c)
    {
      // Do the finite difference check only for values which are greater than some threshold
      if (rValFD[c]>1e-12)
      {
        // check if the column belongs to the slave or master side
        std::string cSlMa = "(S)";
        if (gsdofrowmap_->LID(cIdsFD[c])==-1)
          cSlMa = "(M)";

        // search for entry in the analytical solution
        int anaId = -1;
        for (int cc=0;cc<rLengthAna;++cc)
          if (cIdsFD[c]==cIdsAna[cc])
          {
            anaId = cc;
            break;
          }

        if(anaId==-1)
          std::cout << "*** WARNING: Global column #" << cIdsFD[c] << " in global row #" << rowId <<
            " could not be found in the analytical derivative matrix! (fd= " << rValFD[c] << ") ***" << std::endl;
        else
        {
          double dev = rValFD[c]-rValAna[anaId];

          std::cout << cIdsFD[c] << cSlMa << ":"
              "   fd="  << std::setw(14) << std::setprecision(5) << std::scientific << rValFD[c] <<
              "   ana=" << std::setw(14) << std::setprecision(5) << std::scientific << rValAna[anaId] <<
              "   DEVIATION=" << std::setw(14) << std::setprecision(5) << std::scientific << dev <<
              "   REL-ERROR [%]=" << std::setw(14) << std::setprecision(5) << std::scientific << abs(dev/rValFD[c])*100;

          if( abs(dev) > 1.0e-4 )
          {
            std::cout << " ***** WARNING ***** ";
            w++;
          }
          else if( abs(dev) > 1.0e-5 )
          {
            std::cout << " ***** warning ***** ";
            w++;
          }

          std::cout << std::endl;
        }
      }
    }
    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** " << std::endl;
  }

  // deactivate global finite difference indicator
  gFdCheck_=false;

  return;
}

