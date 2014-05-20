/*----------------------------------------------------------------------*/
/*!

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>

!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_manager.H"
#include "matpar_manager_uniform.H"

#include "invana_utils.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_element.H"

#include "../drt_mat/material.H"
#include "../drt_mat/aaaneohooke_stopro.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_inpar/inpar_material.H"
#include "../drt_comm/comm_utils.H"

STR::INVANA::MatParManagerUniform::MatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret)
   :MatParManager(discret)
{
  paramlayoutmap_ = Teuchos::rcp(new Epetra_Map(1,1,0,*(DRT::Problem::Instance()->GetNPGroup()->LocalComm())));
  paramlayoutmapunique_ = LINALG::AllreduceEMap(*paramlayoutmap_,0);

  optparams_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,NumParams(),true));
  optparams_o_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,NumParams(),true));

  //initialize parameter vector from material parameters given in the input file
  InitParams();

  // set the parameters to be available for the elements
  SetParams();

}

void STR::INVANA::MatParManagerUniform::FillParameters(Teuchos::RCP<Epetra_MultiVector> params)
{
  // zero out to make sure
  params->PutScalar(0.0);

  for (int i=0; i<NumParams(); i++)
    (*params)(i)->PutScalar((*(*optparams_)(i))[0]);
}

void STR::INVANA::MatParManagerUniform::InitParameters(int parapos, double val)
{
  (*optparams_)(parapos)->PutScalar(val);
}

void STR::INVANA::MatParManagerUniform::ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
                                                         double val,
                                                         int elepos,
                                                         int paraposglobal,
                                                         int paraposlocal)
{
  // only row elements contribute
  if (not discret_->ElementRowMap()->MyGID(elepos)) return;

  // every proc can do the 'product rule' on his own since the uniformly ditributed optparams are kept redudantly
  int success = dfint->SumIntoGlobalValue(0,paraposglobal,val);
  if (success!=0) dserror("gid %d is not on this processor", elepos);
}

void STR::INVANA::MatParManagerUniform::Consolidate(Teuchos::RCP<Epetra_MultiVector> dfint)
{

  double val=0.0;
  for (int i=0; i<dfint->NumVectors(); i++)
  {
    val = 0.0;
    discret_->Comm().SumAll((*dfint)(i)->Values(),&val,1);
    dfint->ReplaceGlobalValue(0,i,val);
  }
}

STR::INVANA::MatParManagerPerElement::MatParManagerPerElement(Teuchos::RCP<DRT::Discretization> discret)
   :MatParManager(discret)
{
  // temp map to keep correspondence of parameter block position and eleids
  // used to build the mapextractor and the various maps to keep track of parameters and elements
  std::map< int,std::vector<int> > elemap;
  int nummyparams=0;

  // fill it
  for (int i=0; i<discret_->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = discret_->lRowElement(i);
    int elematid = actele->Material()->Parameter()->Id();

    if (paramap_.find(elematid) == paramap_.end() )
      continue;

    std::vector<int> actparapos = parapos_.at( elematid );
    std::vector<int>::const_iterator it;
    for ( it=actparapos.begin(); it!=actparapos.end(); it++)
    {
      elemap[*it].push_back(actele->Id());
      nummyparams++;
    }
  }

  // generate global ids plus build map paramsLIDtoeleGID_
  std::map<int, std::vector<int> > gids;
  int count = 0;
  for (int i=0; i<discret_->Comm().NumProc(); i++)
  {
    if (discret_->Comm().MyPID() == i)
    {
      for (int j=0; j<NumParams(); j++)
      {
        for (int k=0; k<(int)elemap[j].size(); k++)
        {
          gids[j].push_back(count);
          paramsLIDtoeleGID_.push_back(elemap[j].at(k));
          count++;
        }
      }
    }
    discret_->Comm().Broadcast(&count,1,i);
  }

  //build map eleGIDtoparamsLID_
  for (int i=0; i<(int)paramsLIDtoeleGID_.size(); i++)
  {
    // the blocks are ordered so this just
    eleGIDtoparamsLID_[paramsLIDtoeleGID_[i]].push_back(i);
  }

  // the full map of the vector layout
  paramlayoutmap_ = Teuchos::rcp(new Epetra_Map(-1,nummyparams,0,*(DRT::Problem::Instance()->GetNPGroup()->LocalComm())));
  paramlayoutmapunique_ = Teuchos::rcp(new Epetra_Map(*paramlayoutmap_));

  // the partial maps:
  std::vector< Teuchos::RCP<const Epetra_Map> > partials;
  for (int i=0; i<NumParams(); i++)
  {
    partials.push_back(Teuchos::rcp(new Epetra_Map(-1,gids[i].size(),gids[i].data(),0,discret_->Comm())));
  }

  // finally build the MapExtractor
  paramapextractor_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*paramlayoutmap_,partials));

  optparams_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,NumVectors(),true));
  optparams_o_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,NumVectors(),true));

  //initialize parameter vector from material parameters given in the input file
  InitParams();

  // set the parameters to be available for the elements
  SetParams();

}

void STR::INVANA::MatParManagerPerElement::FillParameters(Teuchos::RCP<Epetra_MultiVector> params)
{
  params->PutScalar(0.0);

  // loop the parameter blocks
  for (int k=0; k<paramapextractor_->NumMaps(); k++)
  {
    Teuchos::RCP<Epetra_Vector> tmp = paramapextractor_->ExtractVector(*(*optparams_)(0),k);
    for (int i=0; i< tmp->MyLength(); i++)
    {
      int pgid = tmp->Map().GID(i);
      int plid = paramlayoutmap_->LID(pgid);
      params->ReplaceGlobalValue(paramsLIDtoeleGID_[plid],k,(*tmp)[i]);
    }
  }
}

void STR::INVANA::MatParManagerPerElement::InitParameters(int parapos, double val)
{
  Teuchos::RCP<Epetra_Vector> tmp = Teuchos::rcp(new Epetra_Vector(*paramapextractor_->Map(parapos), false));
  tmp->PutScalar(val);

  paramapextractor_->InsertVector(tmp,parapos,Teuchos::rcp((*optparams_)(0),false));

}

void STR::INVANA::MatParManagerPerElement::ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
                                                            double val,
                                                            int elepos,
                                                            int paraposglobal,
                                                            int paraposlocal)
{
  if (eleGIDtoparamsLID_.find(elepos) == eleGIDtoparamsLID_.end())
    dserror("proc %d, ele %d not in this map", discret_->Comm().MyPID(), elepos);

  int plid = eleGIDtoparamsLID_[elepos].at(paraposlocal);
  int success = dfint->SumIntoMyValue(plid,0,val);
  if (success!=0) dserror("gid %d is not on this processor", plid);
}
