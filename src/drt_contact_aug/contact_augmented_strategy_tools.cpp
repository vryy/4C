/*---------------------------------------------------------------------*/
/*!
\file contact_augmented_strategy_tools.cpp

\brief Tools for the augmented contact solving strategy with
       standard Lagrangian multipliers.

\level 3

\maintainer Michael Hiermeier

\date Jun 3, 2014

*/
/*---------------------------------------------------------------------*/

#include "contact_augmented_strategy.H"
#include "contact_augmented_interface.H"
#include "../drt_contact/contact_defines.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::AugFDCheckGP(
    CONTACT::ParamsInterface& cparams)
{
#ifdef CONTACTFD_KAPPA
  // FD check of the linearization of the scaling parameter kappa
  for (int i=0;i<(int)interface_.size();++i)
    interface_[i]->FDCheckKappaLin(cparams);
#endif

#ifdef CONTACTFD_AWGAP
  // FD check of the linearization
  for (int i=0;i<(int)interface_.size();++i)
    interface_[i]->FDCheckAWGapLin(cparams);
#endif

#ifdef CONTACTFD_VARWGAPSL
  // FD check of the linearization
  // of the variation of WGap (slave)
  for (int i=0;i<(int)interface_.size();++i)
    interface_[i]->FDCheckVarWGapLinSl(cparams);
#endif

#ifdef CONTACTFD_VARWGAPMA
  // FD check of the linearization
  // of the variation of WGap (master)
  for (int i=0;i<(int)interface_.size();++i)
    interface_[i]->FDCheckVarWGapLinMa(cparams);
#endif

#ifdef CONTACTFD_AUGA
  // FD check of the linearization of the scaling parameter kappa
  for (int i=0;i<(int)interface_.size();++i)
    interface_[i]->FDCheckAugALin(cparams);
#endif

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::AugFDCheckGlobal(
    CONTACT::ParamsInterface& cparams)
{
  // *** linearization w.r.t. displ. *********************************
  /*---------------------------------*
   | Active normal constraint        |
   *---------------------------------*/
#ifdef CONTACTFD_DLMNWGAPLINMATRIX
  if (!Data().FiniteDifferenceIndicator())
    EvalFDCheckGlobalDispl(*dLmNWGapLinMatrix_,dLmNWGapRhs_,cparams);
#endif
  /*---------------------------------*
   | Active tangential constraint    |
   *---------------------------------*/
#ifdef CONTACTFD_DLMTLMTLINMATRIX
  if (!Data().FiniteDifferenceIndicator())
    EvalFDCheckGlobalDispl(*dLmTLmTLinMatrix_,dLmTLmTRhs_,cparams);
#endif
  /*---------------------------------*
   | Active tangential constraint    |
   *---------------------------------*/
#ifdef CONTACTFD_AUGINACTIVELINMATRIX
  if (!Data().FiniteDifferenceIndicator())
    EvalFDCheckGlobalDispl(*augInactiveLinMatrix_,augInactiveRhs_,cparams);
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
  if (!Data().FiniteDifferenceIndicator())
    EvalFDCheckGlobalDispl(*dGLmSlLinMatrix_,augfs_lm_,cparams);
#endif
  // *** master side ***
#ifdef CONTACTFD_DGLMMALINMATRIX
  augMnMatrix_->Multiply(true,*z_exp,*augfm_lm_);
  if (!Data().FiniteDifferenceIndicator())
    EvalFDCheckGlobalDispl(*dGLmMaLinMatrix_,augfm_lm_,cparams);
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
  if (!Data().FiniteDifferenceIndicator())
    EvalFDCheckGlobalDispl(*dGGSlLinMatrix_,augfs_g_,cparams);
#endif
  // *** master side ***
#ifdef CONTACTFD_DGGMALINMATRIX
  augMnMatrix_->Multiply(true,*aWGapRhs_,*augfm_g_);
  augfm_g_->Scale(-cn);
  if (!Data().FiniteDifferenceIndicator())
    EvalFDCheckGlobalDispl(*dGGMaLinMatrix_,augfm_g_,cparams);
#endif
#endif

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::EvalFDCheckGlobalDispl(
    LINALG::SparseMatrix derivMatrix,
    Teuchos::RCP<Epetra_Vector>& rhsVector,
    CONTACT::ParamsInterface& cparams)
{
  // activate global finite difference indicator
 Data().FiniteDifferenceIndicator() = true;

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
      iCheck = interface_[i]->UpdateInterfaces(gid,dof,delta,true,cparams);
      if (iCheck) break;
    }
    if (!iCheck) dserror("ERROR: Node % not found!",gid);

    // Update matrix and rhs
    Initialize();
    AssembleContactRHS();
    AssembleContactStiff();

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
      iCheck = interface_[i]->UpdateInterfaces(gid,dof,delta,false,cparams);
      if (iCheck) break;
    }
    if (!iCheck) dserror("ERROR: Node % not found!",gid);
    Initialize();
    AssembleContactRHS();
    AssembleContactStiff();
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
  Data().FiniteDifferenceIndicator() = false;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::MultiplyElementwise(
    const Epetra_Vector& source,
    const Epetra_Map& source2targetMap,
    Epetra_Vector& target,
    const bool inverse) const
{
  // consistency check
  if (source2targetMap.NumMyElements()!=target.Map().NumMyElements())
  {
    dserror("The number of local elements of the source2targetMap and the "
        "target.Map() have to be equal! \n"
        ".........................................................\n"
        "source2targetMap.NumMyElements() = %i on proc %i \n"
        "target.Map().NumMyElements() = %i on proc %i \n"
        ".........................................................\n",
        source2targetMap.NumMyElements(),Comm().MyPID(),
        target.Map().NumMyElements(),Comm().MyPID());
  }

  // nothing to do, if the target map size is equal zero
  if (source2targetMap.NumGlobalElements()==0)
    return;

  Epetra_Vector source_exp = Epetra_Vector(source2targetMap);
  LINALG::Export(source,source_exp);
  int error = 0;
  error = source_exp.ReplaceMap(target.Map());
  if (error) dserror("The source map couldn't be replaced by the target map! (error=%i)", error);

  if (inverse)
    error = target.ReciprocalMultiply(1.0,source_exp,target,0.0);
  else
    error = target.Multiply(1.0,source_exp,target,0.0);

  if (error) dserror("The element-wise multiplication failed! (error=%i)",error);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::RedistributeRowMap(
    const Epetra_Map& refMap,
    Epetra_Map& modMap)
{
  // only for the parallel redistribution case
  if (!ParRedist()) return;

  int count = 0;
  std::vector<int> myGids(refMap.NumMyElements());

  const Teuchos::RCP<Epetra_Map> allreducedMap = LINALG::AllreduceEMap(modMap);

  for (int i=0; i<refMap.NumMyElements(); ++i)
  {
    int gid = refMap.GID(i);
    if (allreducedMap->LID(gid)>=0)
    {
      myGids[count] = gid;
      ++count;
    }
  }

  myGids.resize(count);
  int gCount=0;
  Comm().SumAll(&count,&gCount,1);
  modMap = Epetra_Map(gCount,count,&myGids[0],0,Comm());

  return;
}
