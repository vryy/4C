#ifdef CCADISCRET

#include <vector>
#include <algorithm>

#include "mfsi_coupling.H"
#include "../drt_fsi/fsi_utils.H"


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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::Coupling::SetupInnerDofMaps(const Epetra_Map& masterdomainmap, const Epetra_Map& slavedomainmap)
{
  // do the master inner map
  int* mastergids = masterdomainmap.MyGlobalElements();
  int masterelements = masterdomainmap.NumMyElements();

  std::vector<int> master;
  master.reserve(masterelements);

  std::remove_copy_if(mastergids,
                      mastergids+masterelements,
                      back_inserter(master),
                      FSI::Utils::MyGID(&*MasterDofMap()));

  masterinnerdofmap_ = Teuchos::rcp(new Epetra_Map(-1,master.size(),&master[0],0,masterdomainmap.Comm()));

  // do the slave inner map
  int* slavegids = slavedomainmap.MyGlobalElements();
  int slaveelements = slavedomainmap.NumMyElements();

  std::vector<int> slave;
  slave.reserve(slaveelements);

  std::remove_copy_if(slavegids,
                      slavegids+slaveelements,
                      back_inserter(slave),
                      FSI::Utils::MyGID(&*SlaveDofMap()));

  slaveinnerdofmap_ = Teuchos::rcp(new Epetra_Map(-1,slave.size(),&slave[0],0,slavedomainmap.Comm()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix> MFSI::Coupling::RCSlaveToMaster(const Teuchos::RCP<Epetra_CrsMatrix> s) const
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix> MFSI::Coupling::RSlaveToMaster(const Teuchos::RCP<Epetra_CrsMatrix> s) const
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix> MFSI::Coupling::CSlaveToMaster(const Teuchos::RCP<Epetra_CrsMatrix> s) const
{
  return Teuchos::null;
}

#endif
