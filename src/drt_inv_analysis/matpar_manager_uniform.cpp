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
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"

STR::INVANA::MatParManagerUniform::MatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret)
   :MatParManager(discret)
{
  paramlayoutmap_ = Teuchos::rcp(new Epetra_Map(1,1,0,*(DRT::Problem::Instance()->GetNPGroup()->LocalComm())));
  int numeleperproc=0;
  if (Discret()->Comm().MyPID()==0) numeleperproc = 1;
  paramlayoutmapunique_ = Teuchos::rcp(new Epetra_Map(1,numeleperproc,0,*(DRT::Problem::Instance()->GetNPGroup()->LocalComm())));

  // build the mapextractor
  // the partial maps (only one map in case of uniform material parameters)
  std::vector< Teuchos::RCP<const Epetra_Map> > partials;
  partials.push_back(paramlayoutmapunique_);
  paramapextractor_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*paramlayoutmapunique_,partials));

  optparams_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,NumParams(),true));
  optparams_o_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,NumParams(),true));
  optparams_initial_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,NumParams(),true));

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
  // only this proc's row elements contribute
  if (not Discret()->ElementRowMap()->MyGID(elepos)) return;

  // every proc can do the 'product rule' on his own since the uniformly ditributed optparams are kept redundantly
  int success = dfint->SumIntoGlobalValue(0,paraposglobal,val);
  if (success!=0) dserror("error code %d", success);
}

void STR::INVANA::MatParManagerUniform::Consolidate(Teuchos::RCP<Epetra_MultiVector> dfint)
{

  double val=0.0;
  for (int i=0; i<dfint->NumVectors(); i++)
  {
    val = 0.0;
    Discret()->Comm().SumAll((*dfint)(i)->Values(),&val,1);
    dfint->ReplaceGlobalValue(0,i,val);
  }
}
