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

INVANA::MatParManagerUniform::MatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret)
   :MatParManager(discret)
{}

/*----------------------------------------------------------------------*/
/* Setup                                                    keh 10/14   */
/*----------------------------------------------------------------------*/
void INVANA::MatParManagerUniform::Setup()
{
  paramlayoutmap_ = Teuchos::rcp(new Epetra_Map(NumParams(),NumParams(),0,*(DRT::Problem::Instance()->GetNPGroup()->LocalComm())));
  int numeleperproc=0;
  if (Discret()->Comm().MyPID()==0) numeleperproc = NumParams();
  paramlayoutmapunique_ = Teuchos::rcp(new Epetra_Map(-1,numeleperproc,0,*(DRT::Problem::Instance()->GetNPGroup()->LocalComm())));

  // build the mapextractor
  // the partial maps (only one map in case of uniform material parameters)
  std::vector< Teuchos::RCP<const Epetra_Map> > partials;
  partials.push_back(paramlayoutmapunique_);
  paramapextractor_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*paramlayoutmapunique_,partials));

  optparams_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,NumVectors(),true));
  optparams_initial_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,NumVectors(),true));

  //initialize parameter vector from material parameters given in the input file
  InitParams();
}

void INVANA::MatParManagerUniform::FillParameters(Teuchos::RCP<Epetra_MultiVector> params)
{
  for (int i=0; i<NumParams(); i++)
    (*params)(i)->PutScalar((*(*optparams_)(0))[i]);
}

void INVANA::MatParManagerUniform::InitParameters(int parapos, double val)
{
  optparams_->ReplaceGlobalValue(parapos,0,val);
}

void INVANA::MatParManagerUniform::ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
                                                         double val,
                                                         int elepos,
                                                         int paraposglobal,
                                                         int paraposlocal)
{
  // only this proc's row elements contribute
  if (not Discret()->ElementRowMap()->MyGID(elepos)) return;

  // every proc can do the 'product rule' on his own since the uniformly distributed optparams are kept redundantly
  int success = dfint->SumIntoGlobalValue(paraposglobal,0,val);
  if (success!=0) dserror("error code %d", success);
}

void INVANA::MatParManagerUniform::Finalize(Teuchos::RCP<Epetra_MultiVector> dfint)
{
  std::vector<double> val(dfint->MyLength(),0.0);
  Discret()->Comm().SumAll((*dfint)(0)->Values(),&val[0],dfint->MyLength());

  for (int i=0; i<NumParams(); i++)
    dfint->ReplaceGlobalValue(i,0,val[i]);

}
