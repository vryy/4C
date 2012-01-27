/*!------------------------------------------------------------------------------------------------*
\file topopt_optimizer.cpp

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "topopt_optimizer.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"



TOPOPT::Optimizer::Optimizer(
    RCP<DRT::Discretization> discret,
    const ParameterList& params
) :
discret_(discret),
params_(params.sublist("TOPOLOGY OPTIMIZER"))
{
  // topology density fields
  dens_i_ = rcp(new Epetra_Vector(*discret_->NodeRowMap(),false));
  dens_ip_ = rcp(new Epetra_Vector(*discret_->NodeRowMap(),false));

  // fluid fields for optimization
  vel_ = rcp(new map<int,Epetra_Vector>);
  adjointvel_ = rcp(new map<int,Epetra_Vector>);

  SetInitialDensityField(DRT::INPUT::IntegralValue<INPAR::TOPOPT::InitialField>(params_,"INITIALFIELD"),
      params_.get<int>("INITFUNCNO"));

  return;
}



void TOPOPT::Optimizer::SetInitialDensityField(
    const INPAR::TOPOPT::InitialField initfield,
    const int startfuncno
)
{
  switch (initfield)
  {
  case INPAR::TOPOPT::initfield_zero_field:
  {
    dens_ip_->PutScalar(0.0);
    break;
  }
  case INPAR::TOPOPT::initfield_field_by_function:
  {
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);

      // evaluate component k of spatial function
      double initialval = DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(0,lnode->X(),0.0,NULL); // scalar
      dens_ip_->ReplaceMyValues(1,&initialval,&lnodeid); // lnodeid = ldofid
    }
    break;
  }
  default:
    dserror("unknown initial field");
  }

  dens_i_->Update(1.0,*dens_ip_,0.0);
}



void TOPOPT::Optimizer::ImportFluidData(
    RCP<Epetra_Vector> vel,
    int step)
{
  vel_->insert(pair<int,Epetra_Vector>(step,*vel));
}



void TOPOPT::Optimizer::ImportAdjointFluidData(
    RCP<Epetra_Vector> vel,
    int step)
{
  adjointvel_->insert(pair<int,Epetra_Vector>(step,*vel));
}



const Epetra_Map* TOPOPT::Optimizer::RowMap()
{
  return discret_->NodeRowMap();
}



const Epetra_Map* TOPOPT::Optimizer::ColMap()
{
  return discret_->NodeColMap();
}



#endif  // #ifdef CCADISCRET
