/*----------------------------------------------------------------------*/
/*!
\file pat_matpar_manager.H
\brief manage material parameters during optimization

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            089 - 289-15271
</pre>

!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "pat_matpar_manager.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_scatra_ele/scatra_ele_action.H"

/*----------------------------------------------------------------------*/
ACOU::PatMatParManagerUniform::PatMatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret)
:MatParManagerUniform(discret)
{
// TODO think about different inheritance
}

/*----------------------------------------------------------------------*/
ACOU::PatMatParManagerPerElement::PatMatParManagerPerElement(Teuchos::RCP<DRT::Discretization> discret, bool scaleele)
:MatParManagerPerElement(discret)
{
  scalegradele_ = scaleele;
}

/*----------------------------------------------------------------------*/
void ACOU::PatMatParManagerUniform::Evaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint)
{
  Teuchos::RCP<const Epetra_Vector> phi = discret_->GetState("phi");

  for (int i=0; i<discret_->NumMyColElements(); i++)
  {
    DRT::Element* actele;
    actele = discret_->lColElement(i);
    int elematid = actele->Material()->Parameter()->Id();

    if (paramap_.find(elematid) == paramap_.end() )
      continue;

    // list to define routines at elementlevel
    Teuchos::ParameterList p;
    p.set<int>("action",SCATRA::calc_integr_grad_reac);
    p.set<int>("scatratype",INPAR::SCATRA::scatratype_condif);

    // this works if we optimize only with respect to reac
    p.set<bool>("signum_mu",actele->Material()->Parameter()->GetParameter(1,discret_->ElementColMap()->LID(actele->Id()))<0.0);
    p.set<bool>("scaleele",false);

    std::vector<int> actparams = paramap_.at( elematid );
    std::vector<int>::const_iterator it;
    for ( it=actparams.begin(); it!=actparams.end(); it++)
    {
      //initialize element vectors
      int ndof = actele->NumNode();
      Epetra_SerialDenseMatrix elematrix1(ndof,ndof,false);
      Epetra_SerialDenseMatrix elematrix2(ndof,ndof,false);
      Epetra_SerialDenseVector elevector1(ndof);
      Epetra_SerialDenseVector elevector2(ndof);
      Epetra_SerialDenseVector elevector3(ndof);

      DRT::Element::LocationArray la(discret_->NumDofSets());
      actele->LocationVector(*discret_,la,false);
      actele->Evaluate(p,*discret_,la,elematrix1,elematrix2,elevector1,elevector2,elevector3);

      //reuse elevector2
      for (int l=0; l<(int)la[0].lm_.size(); l++)
      {
        int lid=phi->Map().LID(la[0].lm_.at(l));
        if (lid==-1) dserror("not found on this processor");
        elevector2[l] = (*phi)[lid];
      }
      double val2 = elevector2.Dot(elevector1);

      // Assemble the final gradient; this is parametrization class business
      // (i.e contraction to (optimization)-parameter space:
      ContractGradient(dfint,val2,actele->Id(),parapos_.at(elematid).at(it-actparams.begin()));

    }//loop this elements material parameters (only the ones to be optimized)

  }//loop elements

  Consolidate(dfint);
  return;
}
/*----------------------------------------------------------------------*/
void ACOU::PatMatParManagerPerElement::Evaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint)
{
  Teuchos::RCP<const Epetra_Vector> phi = discret_->GetState("phi");

  for (int i=0; i<discret_->NumMyColElements(); i++)
  {
    DRT::Element* actele;
    actele = discret_->lColElement(i);
    int elematid = actele->Material()->Parameter()->Id();

    if (paramap_.find(elematid) == paramap_.end() )
      continue;

    // list to define routines at elementlevel
    Teuchos::ParameterList p;
    p.set<int>("action",SCATRA::calc_integr_grad_reac);
    p.set<int>("scatratype",INPAR::SCATRA::scatratype_condif);

    // this works if we optimize only with respect to reac
    p.set<bool>("signum_mu",actele->Material()->Parameter()->GetParameter(1,discret_->ElementColMap()->LID(actele->Id()))<0.0);
    p.set<bool>("scaleele",scalegradele_);

    std::vector<int> actparams = paramap_.at( elematid );
    std::vector<int>::const_iterator it;
    for ( it=actparams.begin(); it!=actparams.end(); it++)
    {
      //initialize element vectors
      int ndof = actele->NumNode();
      Epetra_SerialDenseMatrix elematrix1(ndof,ndof,false);
      Epetra_SerialDenseMatrix elematrix2(ndof,ndof,false);
      Epetra_SerialDenseVector elevector1(ndof);
      Epetra_SerialDenseVector elevector2(ndof);
      Epetra_SerialDenseVector elevector3(ndof);

      DRT::Element::LocationArray la(discret_->NumDofSets());
      actele->LocationVector(*discret_,la,false);
      actele->Evaluate(p,*discret_,la,elematrix1,elematrix2,elevector1,elevector2,elevector3);

      //reuse elevector2
      for (int l=0; l<(int)la[0].lm_.size(); l++)
      {
        int lid=phi->Map().LID(la[0].lm_.at(l));
        if (lid==-1) dserror("not found on this processor");
        elevector2[l] = (*phi)[lid];
      }
      double val2 = elevector2.Dot(elevector1);

      // Assemble the final gradient; this is parametrization class business
      // (i.e contraction to (optimization)-parameter space:
      ContractGradient(dfint,val2,actele->Id(),parapos_.at(elematid).at(it-actparams.begin()));

    }//loop this elements material parameters (only the ones to be optimized)

  }//loop elements

  Consolidate(dfint);

  return;
}
