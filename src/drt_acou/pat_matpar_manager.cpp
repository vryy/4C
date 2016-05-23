/*----------------------------------------------------------------------*/
/*!
\file pat_matpar_manager.cpp
\brief manage material parameters during optimization

<pre>
\level 3
\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            089 - 289-15271
</pre>

!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "pat_matpar_manager.H"
#include "acou_ele.H"
#include "acou_ele_action.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*/
ACOU::OptMatParManagerUniform::OptMatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret)
:MatParManagerUniform(discret)
{}

/*----------------------------------------------------------------------*/
ACOU::AcouMatParManagerUniform::AcouMatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret)
:MatParManagerUniform(discret)
{}

/*----------------------------------------------------------------------*/
ACOU::OptMatParManagerPerElement::OptMatParManagerPerElement(Teuchos::RCP<DRT::Discretization> discret)
:MatParManagerPerElement(discret)
{}

/*----------------------------------------------------------------------*/
ACOU::AcouMatParManagerPerElement::AcouMatParManagerPerElement(Teuchos::RCP<DRT::Discretization> discret)
:MatParManagerPerElement(discret)
{}

/*----------------------------------------------------------------------*/
void ACOU::OptMatParManagerUniform::AddEvaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint)
{
  // get the actual set of elementwise material parameters from the derived classes
  Teuchos::RCP<Epetra_MultiVector> getparams = Teuchos::rcp(new Epetra_MultiVector(*(Discret()->ElementRowMap()),NumParams(),false));
  FillParameters(getparams);

  // export to column layout to be able to run column elements
  Discret()->Comm().Barrier();
  LINALG::Export(*getparams,*WriteParamsVec());

  // get solution
  Teuchos::RCP<const Epetra_Vector> phi = Discret()->GetState("phi");
  for (int i=0; i<Discret()->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = Discret()->lRowElement(i);
    int elematid = actele->Material()->Parameter()->Id();

    if (ParaMap().find(elematid) == ParaMap().end() )
      continue;

    // list to define routines at elementlevel
    Teuchos::ParameterList p;
    p.set<int>("action",SCATRA::calc_integr_grad_reac);

    std::vector<int> actparams = ParaMap().at(elematid);
    std::vector<int>::const_iterator it;
    for ( it=actparams.begin(); it!=actparams.end(); it++)
    {
      if((*it) == 1)
        p.set<int>("action",SCATRA::calc_integr_grad_reac);
      else if((*it) == 0)
        p.set<int>("action",SCATRA::calc_integr_grad_diff);
      else
        dserror("unkown elematid provided");

      //initialize element vectors
      int ndof = actele->NumNode();
      Epetra_SerialDenseMatrix elematrix1(ndof,ndof,false);
      Epetra_SerialDenseMatrix elematrix2(ndof,ndof,false);
      Epetra_SerialDenseVector elevector1(ndof);
      Epetra_SerialDenseVector elevector2(ndof);
      Epetra_SerialDenseVector elevector3(ndof);

      DRT::Element::LocationArray la(Discret()->NumDofSets());
      actele->LocationVector(*Discret(),la,false);
      actele->Evaluate(p,*Discret(),la,elematrix1,elematrix2,elevector1,elevector2,elevector3);

      double metaval = (*(*GetMatParamsVec())( ParaPos().at(elematid).at(it-actparams.begin()) ))[actele->LID()];
      double val1 = metaparams_.DMaterialDMeta(metaval);
      elevector1.Scale(val1);

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
      ContractGradient(dfint,val2,actele->Id(),ParaPos().at(elematid).at(it-actparams.begin()), it-actparams.begin());

    }//loop this elements material parameters (only the ones to be optimized)

  }//loop elements

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::AcouMatParManagerUniform::AddEvaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint)
{
  // loop the row elements
  for (int i=0; i<Discret()->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = Discret()->lRowElement(i);
    int elematid = actele->Material()->Parameter()->Id();

    if (ParaMap().find(elematid) == ParaMap().end() )
    {
      //std::cout<<"Warning, skipping elematid "<<elematid<<" in ele "<<actele->Id()<<std::endl;
      continue;
    }
    const DRT::ELEMENTS::Acou * hdgele = dynamic_cast<const DRT::ELEMENTS::Acou*>(actele);
    double val = 0.0;
    std::vector<int> actparams = ParaMap().at( elematid );
    std::vector<int>::const_iterator it;
    for ( it=actparams.begin(); it!=actparams.end(); it++)
    {
      // get the correct value
      if((*it) == 0)
        val = hdgele->GetDensityGradient();
      else if((*it) == 1)
        val = hdgele->GetSoSGradient();
      // write it to the gradient
      ContractGradient(dfint,val,actele->Id(),ParaPos().at(elematid).at(it-actparams.begin()), it-actparams.begin());
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::OptMatParManagerPerElement::AddEvaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint)
{
  // get the actual set of elementwise material parameters from the derived classes
  Teuchos::RCP<Epetra_MultiVector> getparams = Teuchos::rcp(new Epetra_MultiVector(*(Discret()->ElementRowMap()),NumParams(),false));
  FillParameters(getparams);

  // export to column layout to be able to run column elements
  Discret()->Comm().Barrier();
  LINALG::Export(*getparams,*WriteParamsVec());

  // get solution
  Teuchos::RCP<const Epetra_Vector> phi = Discret()->GetState("phi");
  for (int i=0; i<Discret()->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = Discret()->lRowElement(i);
    int elematid = actele->Material()->Parameter()->Id();

    if (ParaMap().find(elematid) == ParaMap().end() )
    {
      //std::cout<<"Warning, skipping elematid "<<elematid<<" in ele "<<actele->Id()<<std::endl;
      continue;
    }
    // list to define routines at elementlevel
    Teuchos::ParameterList p;
    // set parameters
    p.set<bool>("signum_mu",actele->Material()->Parameter()->GetParameter(1,Discret()->ElementColMap()->LID(actele->Id()))<0.0);
    p.set<bool>("signum_D", actele->Material()->Parameter()->GetParameter(0,Discret()->ElementColMap()->LID(actele->Id()))<0.0);

    std::vector<int> actparams = ParaMap().at( elematid );
    std::vector<int>::const_iterator it;
    for ( it=actparams.begin(); it!=actparams.end(); it++)
    {
      if((*it) == 1)
        p.set<int>("action",SCATRA::calc_integr_grad_reac);
      else if((*it) == 0)
        p.set<int>("action",SCATRA::calc_integr_grad_diff);
      else
        dserror("unkown elematid provided");

      //initialize element vectors
      int ndof = actele->NumNode();
      Epetra_SerialDenseMatrix elematrix1(ndof,ndof,false);
      Epetra_SerialDenseMatrix elematrix2(ndof,ndof,false);
      Epetra_SerialDenseVector elevector1(ndof);
      Epetra_SerialDenseVector elevector2(ndof);
      Epetra_SerialDenseVector elevector3(ndof);

      DRT::Element::LocationArray la(Discret()->NumDofSets());
      actele->LocationVector(*Discret(),la,false);
      actele->Evaluate(p,*Discret(),la,elematrix1,elematrix2,elevector1,elevector2,elevector3);

      double metaval = (*(*GetMatParamsVec())( ParaPos().at(elematid).at(it-actparams.begin()) ))[actele->LID()];
      double val1 = metaparams_.DMaterialDMeta(metaval);
      elevector1.Scale(val1);

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
      ContractGradient(dfint,val2,actele->Id(),ParaPos().at(elematid).at(it-actparams.begin()), it-actparams.begin());

    }//loop this elements material parameters (only the ones to be optimized)

  }//loop elements

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::AcouMatParManagerPerElement::AddEvaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint)
{
  // loop the row elements
  for (int i=0; i<Discret()->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = Discret()->lRowElement(i);
    int elematid = actele->Material()->Parameter()->Id();

    if (ParaMap().find(elematid) == ParaMap().end() )
    {
      //std::cout<<"Warning, skipping elematid "<<elematid<<" in ele "<<actele->Id()<<std::endl;
      continue;
    }
    const DRT::ELEMENTS::Acou * hdgele = dynamic_cast<const DRT::ELEMENTS::Acou*>(actele);
    double val = 0.0;
    std::vector<int> actparams = ParaMap().at( elematid );
    std::vector<int>::const_iterator it;
    for ( it=actparams.begin(); it!=actparams.end(); it++)
    {
      // get the correct value
      if((*it) == 0)
        val = hdgele->GetDensityGradient();
      else if((*it) == 1)
        val = hdgele->GetSoSGradient();

      double metaval = (*(*GetMatParamsVec())( ParaPos().at(elematid).at(it-actparams.begin()) ))[actele->LID()];
      double val1 = metaparams_.DMaterialDMeta(metaval);
      val *= val1;

      // write it to the gradient
      ContractGradient(dfint,val,actele->Id(),ParaPos().at(elematid).at(it-actparams.begin()), it-actparams.begin());
    }
  }

  return;
}

//! write parameter values to the correct positions in the given vector
void ACOU::OptMatParManagerPerElement::WriteValuesToVector(int elematid, int geleid, double valit0, double valit1, Teuchos::RCP<Epetra_MultiVector> params)
{
  std::vector<int> actparams = ParaMap().at( elematid );
  std::vector<int>::const_iterator it;
  double val = 0.0;
  for ( it=actparams.begin(); it!=actparams.end(); it++)
  {
    // get the correct value
    if((*it) == 0)
      val = valit0;
    else if((*it) == 1)
      val = valit1;
    // write it to the gradient
    ContractGradient(params,val,geleid,ParaPos().at(elematid).at(it-actparams.begin()), it-actparams.begin());
  }
  return;
}

//! write parameter values to the correct positions in the given vector
void ACOU::AcouMatParManagerPerElement::WriteValuesToVector(int elematid, int geleid, double valit0, double valit1, Teuchos::RCP<Epetra_MultiVector> params)
{
  std::vector<int> actparams = ParaMap().at( elematid );
  std::vector<int>::const_iterator it;
  double val = 0.0;
  for ( it=actparams.begin(); it!=actparams.end(); it++)
  {
    // get the correct value
    if((*it) == 0)
      val = valit0;
    else if((*it) == 1)
      val = valit1;
    // write it to the gradient
    ContractGradient(params,val,geleid,ParaPos().at(elematid).at(it-actparams.begin()), it-actparams.begin());
  }
  return;
}

/*----------------------------------------------------------------------*/
void ACOU::OptMatParManagerPerElement::SetAction(Teuchos::ParameterList& p)
{
  p.set<int>("action",SCATRA::bd_integrate_shape_functions);
  return;
}


/*----------------------------------------------------------------------*/
void ACOU::AcouMatParManagerPerElement::SetAction(Teuchos::ParameterList& p)
{
  p.set<int>("action", ACOU::bd_integrate);
  return;
}
