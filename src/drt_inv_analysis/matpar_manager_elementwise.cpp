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
#include "matpar_manager_elementwise.H"

#include "invana_utils.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_element.H"

#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"

STR::INVANA::MatParManagerPerElement::MatParManagerPerElement(Teuchos::RCP<DRT::Discretization> discret)
   :MatParManager(discret)
{
  // temp map to keep correspondence of parameter block position and eleids
  // used to build the mapextractor and the various maps to keep track of parameters and elements
  std::map< int,std::vector<int> > elemap;
  int nummyparams=0;

  // fill it
  for (int i=0; i<Discret()->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = Discret()->lRowElement(i);
    int elematid = actele->Material()->Parameter()->Id();

    if (ParaMap().find(elematid) == ParaMap().end() )
      continue;

    std::vector<int> actparapos = ParaPos().at( elematid );
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
  for (int i=0; i<Discret()->Comm().NumProc(); i++)
  {
    if (Discret()->Comm().MyPID() == i)
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
    Discret()->Comm().Broadcast(&count,1,i);
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
    partials.push_back(Teuchos::rcp(new Epetra_Map(-1,gids[i].size(),gids[i].data(),0,Discret()->Comm())));
  }

  // finally build the MapExtractor
  paramapextractor_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*paramlayoutmap_,partials));

  optparams_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,NumVectors(),true));
  optparams_o_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,NumVectors(),true));
  optparams_initial_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,NumVectors(),true));

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
    dserror("proc %d, ele %d not in this map", Discret()->Comm().MyPID(), elepos);

  int plid = eleGIDtoparamsLID_[elepos].at(paraposlocal);
  int success = dfint->SumIntoMyValue(plid,0,val);
  if (success!=0) dserror("gid %d is not on this processor", plid);
}
