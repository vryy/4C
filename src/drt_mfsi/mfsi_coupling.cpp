#ifdef CCADISCRET

#include "mfsi_coupling.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::Coupling::SetupCouplingMatrices(const Epetra_Map& masterdomainmap, const Epetra_Map& slavedomainmap)
{
  // we always use the masterdofmap for the domain
  matmm_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*MasterDofMap(),1,true));
  matsm_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*MasterDofMap(),1,true));

  matmm_trans_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,masterdomainmap,1,true));
  matsm_trans_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,slavedomainmap,1,true));

  int length = MasterDofMap()->NumMyElements();
  double one = 1.;
  for (int i=0; i<length; ++i)
  {
    int sgid = SlaveDofMap()->GID(i);
    int mgid = MasterDofMap()->GID(i);

    int err = matmm_->InsertGlobalValues(mgid, 1, &one, &mgid);
    if (err!=0)
      dserror("InsertGlobalValues for entry (%d,%d) failed with err=%d",mgid,mgid,err);

    err = matsm_->InsertGlobalValues(mgid, 1, &one, &sgid);
    if (err!=0)
      dserror("InsertGlobalValues for entry (%d,%d) failed with err=%d",mgid,sgid,err);

    err = matmm_trans_->InsertGlobalValues(mgid, 1, &one, &mgid);
    if (err!=0)
      dserror("InsertGlobalValues for entry (%d,%d) failed with err=%d",mgid,mgid,err);

    err = matsm_trans_->InsertGlobalValues(sgid, 1, &one, &mgid);
    if (err!=0)
      dserror("InsertGlobalValues for entry (%d,%d) failed with err=%d",sgid,mgid,err);
  }

  matmm_->FillComplete(masterdomainmap,*MasterDofMap());
  matsm_->FillComplete(slavedomainmap,*MasterDofMap());

  matmm_trans_->FillComplete(*MasterDofMap(),masterdomainmap);
  matsm_trans_->FillComplete(*MasterDofMap(),slavedomainmap);

}


#endif
