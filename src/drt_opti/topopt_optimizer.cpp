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
  dens_i_ = rcp(new Epetra_Vector(*discret_->DofRowMap(),false));
  dens_ip_ = rcp(new Epetra_Vector(*discret_->DofRowMap(),false));

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
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);

      int numdofs = nodedofset.size();
      for (int k=0;k< numdofs;++k)
      {
        const int dofgid = nodedofset[k];
        int doflid = dofrowmap->LID(dofgid);
        // evaluate component k of spatial function
        double initialval = DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(k,lnode->X(),0.0,NULL);
        dens_ip_->ReplaceMyValues(1,&initialval,&doflid);
      }
    }
    break;
  }
  default:
    dserror("unknown initial field");
  }

  dens_i_->Update(1.0,*dens_ip_,0.0);
}



#endif  // #ifdef CCADISCRET
