/*----------------------------------------------------------------------*/
/*!
\file matpar_manager_uniform.cpp
\brief Creating patches from inputfile
<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_manager_uniform.H"

#include "invana_utils.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_element.H"

#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"

#include <queue>

typedef std::map<int, std::vector<int> > PATCHES;

/*----------------------------------------------------------------------*/
INVANA::MatParManagerUniform::MatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret)
   : MatParManagerPerElement(discret)
{}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerUniform::Setup()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "-----------------------------" << std::endl;
    std::cout << "MatParManagerUniform Setup:" << std::endl;
  }

  // call setup of the Base class to have all the
  // layout of the elementwise distribution
  MatParManagerPerElement::Setup();
  optparams_elewise_ = Teuchos::rcp(new Epetra_MultiVector(*optparams_));
  elewise_map_ = Teuchos::rcp(new Epetra_Map(*paramlayoutmap_));

  // create projection from elements to patches
  CreateProjection();

  // initialize parameters
  InitParameters();

  // Some user information
  if (Comm().MyPID() == 0)
    std::cout << std::endl;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerUniform::FillParameters(Teuchos::RCP<Epetra_MultiVector> params)
{
  params->PutScalar(0.0);

  // Inject into the elementwise solution space
  int err = projector_->Multiply(true,*optparams_,*optparams_elewise_);
  if (err!=0)
    dserror("Application of prolongator failed.");

  // loop the parameter blocks
  for (int k=0; k<paramapextractor_->NumMaps(); k++)
  {
    Teuchos::RCP<Epetra_Vector> tmp = paramapextractor_->ExtractVector(*(*optparams_elewise_)(0),k);
    for (int i=0; i< tmp->MyLength(); i++)
    {
      int pgid = tmp->Map().GID(i); // !! the local id of the partial map is not the local parameter id!!
      int plid = paramapextractor_->FullMap()->LID(pgid);
      params->ReplaceGlobalValue(ParamsLIDtoeleGID()[plid],k,(*tmp)[i]);
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerUniform::ApplyParametrization(
    DcsMatrix& matrix, Teuchos::RCP<Epetra_MultiVector> diagonals)
{
  // this is ok here since we have a sparse approximation
  Teuchos::RCP<Epetra_CrsMatrix> fullmatrix = matrix.FillMatrix();

  // todo: this is not ok! loop over the single columns of matrix
  // and extract only the diagonal component.
  // matrix * projector_
  Teuchos::RCP<Epetra_CrsMatrix> mr = LINALG::Multiply(fullmatrix,false,projector_,false);
  // projector'*matrix*projector
  Teuchos::RCP<Epetra_CrsMatrix> pmr = LINALG::Multiply(projector_,true,mr,false);

  Epetra_Vector diagonal(pmr->RowMap(),true);
  pmr->ExtractDiagonalCopy(diagonal);

  // loop the parameter blocks
  for (int k=0; k<paramapextractor_->NumMaps(); k++)
  {
    Teuchos::RCP<Epetra_Vector> tmp = paramapextractor_->ExtractVector(diagonal,k);
    for (int i=0; i< tmp->MyLength(); i++)
    {
      int pgid = tmp->Map().GID(i); // !! the local id of the partial map is not the local parameter id!!
      int plid = paramapextractor_->FullMap()->LID(pgid);
      diagonals->ReplaceGlobalValue(ParamsLIDtoeleGID()[plid],k,(*tmp)[i]);
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerUniform::InitParameters()
{
  // sanity checks
  if ( not projector_->DomainMap().PointSameAs(optparams_elewise_->Map()))
    dserror("Restrictor->DomainMap error.");

  if ( not projector_->RangeMap().PointSameAs(optparams_->Map()))
    dserror("Restrictor->RangeMap error");

  // parameters are not initialized from input but
  // from the elementwise layout
  int err = projector_->Multiply(false,*optparams_elewise_,*optparams_);
  if (err!=0)
    dserror("Application of restrictor failed.");

  // set initial values
  optparams_initial_->Scale(1.0,*optparams_);

}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerUniform::ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
                                                            double val,
                                                            int elepos,
                                                            int paraposglobal,
                                                            int paraposlocal)
{
  if (EleGIDtoparamsLID().find(elepos) == EleGIDtoparamsLID().end())
    dserror("proc %d, ele %d not in this map", Discret()->Comm().MyPID(), elepos);

  // parameter in the 'elementwise' (there can be multiple
  // per element) parameter layout
  int plid = EleGIDtoparamsLID()[elepos].at(paraposlocal);

  // check to which patch it belongs
  int patchid = pidtopatch_[plid];

  int success = dfint->SumIntoMyValue(patchid,0,val);
  if (success!=0) dserror("gid %d is not on this processor", plid);
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerUniform::Finalize(Teuchos::RCP<Epetra_MultiVector> source,
    Teuchos::RCP<Epetra_MultiVector> target)
{
  // sum across processor
  std::vector<double> val(source->MyLength(),0.0);
  Discret()->Comm().SumAll((*source)(0)->Values(),&val[0],source->MyLength());

  for (int i=0; i<target->MyLength(); i++)
    target->SumIntoGlobalValue(i,0,val[i]);

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerUniform::CreateProjection()
{
  // create optimization parameter maps anew
  int numeleperproc=0;
  if (Comm().MyPID()==0) numeleperproc = NumParams();
  paramlayoutmapunique_ = Teuchos::rcp(new Epetra_Map(-1,numeleperproc,0,Comm()));
  paramlayoutmap_ = LINALG::AllreduceEMap(*paramlayoutmapunique_);

  // build correspondence elewiseoptparam <-> patch
  for (int i=0; i<paramapextractor_->NumMaps(); i++)
  {
    const Teuchos::RCP<const Epetra_Map>& map = paramapextractor_->Map(i);
    int mylength = map->NumMyElements();
    for (int j=0; j<mylength; j++)
    {
      int lid = paramapextractor_->FullMap()->LID(map->GID(j));
      pidtopatch_[lid]=i;
    }
  }

  //check bandwidth of restrictor and prolongator
  int maxbw = 0;
  for (int i=0; i<paramapextractor_->NumMaps(); i++)
  {
    int current = paramapextractor_->Map(i)->NumMyElements();
    if (current > maxbw)
      maxbw = current;
  }

  Teuchos::RCP<Epetra_Map> colmap = LINALG::AllreduceEMap(*paramapextractor_->FullMap(),0);
  projector_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*paramlayoutmapunique_,*colmap,maxbw,false));

  // the projections will have a range only on proc 0
  // but the partial maps might be distributed
  // -> reduce them before inserting into projection
  std::vector<Teuchos::RCP<Epetra_Map> > maps;
  for (int i=0; i<paramapextractor_->NumMaps(); i++)
    maps.push_back(LINALG::AllreduceEMap(*paramapextractor_->Map(i),0));

  for (int i=0; i<projector_->NumMyRows(); i++)
  {
    // scale length of the row to 1
    int numentries = maps[i]->NumMyElements();
    std::vector<double> values(numentries, 1.0/sqrt(numentries));

    int err = projector_->InsertGlobalValues(i,numentries,&values[0],maps[i]->MyGlobalElements());
    if (err < 0 )
      dserror("Restrictor/Prolongator insertion failed.");
  }
  int err = projector_->FillComplete(*paramapextractor_->FullMap(), *paramlayoutmapunique_,true);
  if (err != 0 )
    dserror("Restrictor/Prolongator FillComplete failed.");


  optparams_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmapunique_,1,true));
  optparams_initial_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmapunique_,1,true));

  return;
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix> INVANA::MatParManagerUniform::InitialCovariance()
{
  // the best inital guess to get from datfile input is a
  // unit diagonal.
  Teuchos::RCP<Epetra_CrsMatrix> cov = Teuchos::rcp(new
      Epetra_CrsMatrix(Copy,*paramlayoutmapunique_,*paramlayoutmapunique_,1,false));

  // insert ones onto the diagonal
  double values=1.0;
  for (int i=0; i<cov->NumMyRows(); i++)
  {
    int gid = paramlayoutmapunique_->GID(i);
    cov->InsertGlobalValues(gid,1,&values,&gid);
  }
  cov->FillComplete();

  return cov;
}
