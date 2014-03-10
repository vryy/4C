/*----------------------------------------------------------------------*/
/*!

<pre>
Maintainer: Jonas Biehler
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_manager.H"
#include "matpar_manager_uniform.H"

#include "invana_utils.H"
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

  optparams_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,numparams_,true));
  optparams_o_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,numparams_,true));

  //initialize parameter vector from material parameters given in the input file
  InitParams();

  // set the parameters to be available for the elements
  SetParams();

}

void STR::INVANA::MatParManagerUniform::FillParameters(Teuchos::RCP<Epetra_MultiVector> params)
{
  params->PutScalar(0.0);

  for (int i=0; i<numparams_; i++)
    (*params)(i)->PutScalar((*(*optparams_)(i))[0]);
}


void STR::INVANA::MatParManagerUniform::ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
                                                         Teuchos::RCP<Epetra_Vector> dfintp,
                                                         int parapos)
{
  double val;
  Epetra_Vector ones(*(discret_->ElementRowMap()),true);
  ones.PutScalar(1.0);

  STR::INVANA::MVDotProduct(dfintp,Teuchos::rcp(&ones),&val,discret_->ElementRowMap());
  dfint->SumIntoGlobalValue(0,parapos,val);
}

STR::INVANA::MatParManagerPerElement::MatParManagerPerElement(Teuchos::RCP<DRT::Discretization> discret)
   :MatParManager(discret)
{
  paramlayoutmap_ = Teuchos::rcp(new Epetra_Map(*(discret->ElementColMap())));

  optparams_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,numparams_,true));
  optparams_o_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,numparams_,true));

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
                                                            Teuchos::RCP<Epetra_Vector> dfintp,
                                                            int parapos)
{
  (*dfint)(parapos)->Update(1.0,*dfintp,1.0);
}
