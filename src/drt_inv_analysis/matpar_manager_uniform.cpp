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


void STR::INVANA::MatParManagerUniform::ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
                                                         double val,
                                                         int elepos,
                                                         int parapos)
{
  // only row elements contribute
  if (not discret_->ElementRowMap()->MyGID(elepos)) return;

  // no every proc should get this
  int success = dfint->SumIntoGlobalValue(0,parapos,val);
  if (success!=0) dserror("gid %d is not on this processor", elepos);
}

void STR::INVANA::MatParManagerUniform::Consolidate(Teuchos::RCP<Epetra_MultiVector> dfint)
{

  double val=0.0;
  for (int i=0; i<dfint->NumVectors(); i++)
  {
    discret_->Comm().SumAll((*dfint)(i)->Values(),&val,1);
    dfint->ReplaceGlobalValue(0,i,val);
  }
}

STR::INVANA::MatParManagerPerElement::MatParManagerPerElement(Teuchos::RCP<DRT::Discretization> discret)
   :MatParManager(discret)
{
  paramlayoutmap_ = Teuchos::rcp(new Epetra_Map(*(discret->ElementColMap())));
  paramlayoutmapunique_ = Teuchos::rcp(new Epetra_Map(*(discret->ElementRowMap())));

  optparams_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,NumParams(),true));
  optparams_o_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,NumParams(),true));

  //initialize parameter vector from material parameters given in the input file
  InitParams();

  // set the parameters to be available for the elements
  SetParams();

}

void STR::INVANA::MatParManagerPerElement::FillParameters(Teuchos::RCP<Epetra_MultiVector> params)
{
  params->Scale(1.0,*optparams_);
}

void STR::INVANA::MatParManagerPerElement::ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
                                                         double val,
                                                         int elepos,
                                                         int parapos)
{
  int success = dfint->SumIntoGlobalValue(elepos,parapos,val);
  if (success!=0) dserror("gid %d is not on this processor", elepos);
}
