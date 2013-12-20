/*----------------------------------------------------------------------*/
/*!
\file fluid_timint_bdf2.cpp
\brief BDF2 time-integration scheme

<pre>
Maintainers: Ursula Rasthofer & Martin Kronbichler
             {rasthofer,kronbichler}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-235
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_timint_bdf2.H"
#include "../drt_fluid_turbulence/scale_sep_gmo.H"
#include "../drt_io/io.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_turbulence/dyn_smag.H"
#include "../drt_fluid_turbulence/dyn_vreman.H"
#include "../drt_fluid_turbulence/boxfilter.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntBDF2::TimIntBDF2(
    const Teuchos::RCP<DRT::Discretization>&      actdis,
        const Teuchos::RCP<LINALG::Solver>&           solver,
        const Teuchos::RCP<Teuchos::ParameterList>&   params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output,
        bool                                          alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis,solver,params,output,alefluid),
  theta_(1.0)
{

  //check, if starting algorithm is desired
  if (numstasteps_ > 0)
    dserror("no starting algorithm supported for schemes other than af-gen-alpha");

  SetElementTimeParameter();

  Initialize();

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                    bk 11/13 |
*----------------------------------------------------------------------*/
FLD::TimIntBDF2::~TimIntBDF2()
{
  return;
}

/*----------------------------------------------------------------------*
| Print information about current time step to screen          bk 11/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntBDF2::PrintTimeStepInfo()
{
  if (myrank_==0)
  {
    printf("TIME: %11.4E/%11.4E  DT = %11.4E       BDF2          STEP = %4d/%4d \n",
           time_,maxtime_,dta_,step_,stepmax_);
  }
  return;
}

/*----------------------------------------------------------------------*
| calculate pseudo-theta for startalgo_                        bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntBDF2::SetTheta()
{
  // for BDF2, theta is set by the time-step sizes, 2/3 for const. dt

  if (step_ > 1)
    theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);
  else
  {
    // use backward Euler for the first time step
    velnm_->Update(1.0,*veln_,0.0); // results in hist_ = veln_
    theta_ = 1.0;
  }

  return;
}

/*----------------------------------------------------------------------*
| set old part of right hand side                              bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntBDF2::SetOldPartOfRighthandside()
{
  /*
     BDF2: for constant time step:

                   mom: hist_ = 4/3 veln_  - 1/3 velnm_
                  (con: hist_ = 4/3 densn_ - 1/3 densnm_)

  */

      hist_->Update(4./3., *veln_, -1./3., *velnm_, 0.0);

  return;
}

/*----------------------------------------------------------------------*
| set integration-scheme-specific state                        bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntBDF2::SetStateTimInt()
{

  discret_->SetState("velaf",velnp_);

  return;
}

/*----------------------------------------------------------------------*
| calculate acceleration                                       bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntBDF2::CalculateAcceleration(
    const Teuchos::RCP<const Epetra_Vector>    velnp,
    const Teuchos::RCP<const Epetra_Vector>    veln,
    const Teuchos::RCP<const Epetra_Vector>    velnm,
    const Teuchos::RCP<const Epetra_Vector>    accn,
    const Teuchos::RCP<Epetra_Vector>          accnp
)
{

  /*

  BDF2:

                 2*dt(n)+dt(n-1)                  dt(n)+dt(n-1)
   acc(n+1) = --------------------- vel(n+1) - --------------- vel(n)
               dt(n)*[dt(n)+dt(n-1)]              dt(n)*dt(n-1)

                       dt(n)
             + ----------------------- vel(n-1)
               dt(n-1)*[dt(n)+dt(n-1)]

  */

  if (dta_*dtp_ < EPS15) dserror("Zero time step size!!!!!");
  const double sum = dta_ + dtp_;

  accnp->Update((2.0*dta_+dtp_)/(dta_*sum),*velnp, -sum/(dta_*dtp_),*veln ,0.0);
  accnp->Update(dta_/(dtp_*sum),*velnm,1.0);

  return;
}

/*----------------------------------------------------------------------*
| set gamma                                                    bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntBDF2::SetGamma(Teuchos::ParameterList& eleparams)
{

  eleparams.set("gamma"  ,1.0);
  return;
}

/*----------------------------------------------------------------------*
| scale separation                                             bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntBDF2::Sep_Multiply()
{
  Sep_->Multiply(false,*velnp_,*fsvelaf_);
  return;
}

/*----------------------------------------------------------------------*
 | paraview output of filtered velocity                  rasthofer 02/11|
 *----------------------------------------------------------------------*/
void FLD::TimIntBDF2::OutputofFilteredVel(
     Teuchos::RCP<Epetra_Vector> outvec,
     Teuchos::RCP<Epetra_Vector> fsoutvec)
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  RCP<Epetra_Vector> row_finescaleveltmp;
  row_finescaleveltmp = Teuchos::rcp(new Epetra_Vector(*dofrowmap,true));

  // get fine scale velocity
  if (scale_sep_ == INPAR::FLUID::algebraic_multigrid_operator)
    Sep_->Multiply(false,*velnp_,*row_finescaleveltmp);
  else
    ScaleSepGMO_->ApplyScaleSeparation(velnp_,row_finescaleveltmp);
  // get filtered or coarse scale velocity
  outvec->Update(1.0,*velnp_,-1.0,*row_finescaleveltmp,0.0);

  fsoutvec->Update(1.0,*row_finescaleveltmp,0.0);

  return;
}

// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------
void FLD::TimIntBDF2::SetElementTimeParameter()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action",FLD::set_time_parameter);

  //set time integration scheme
  eleparams.set<int>("TimeIntegrationScheme", timealgo_);

  // set general element parameters
  eleparams.set("dt",dta_);
  eleparams.set("theta",theta_);
  eleparams.set("omtheta",0.0);

  // set scheme-specific element parameters and vector values
  eleparams.set("total time",time_);


  // call standard loop over elements
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  return;
}

/*----------------------------------------------------------------------*
| return time integration factor                               bk 12/13 |
*-----------------------------------------------------------------------*/
const double FLD::TimIntBDF2::TimIntParam() const
{
  double retval = 0.0;
  return retval;
}

/*----------------------------------------------------------------------*
 | filtered quantities for classical LES models          rasthofer 02/11|
 *----------------------------------------------------------------------*/
void FLD::TimIntBDF2::ApplyScaleSeparationForLES()
{
  if (turbmodel_==INPAR::FLUID::dynamic_smagorinsky)
  {
    // perform filtering and computation of Cs
    // compute averaged values for LijMij and MijMij
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = Dirichlet();
    DynSmag_->ApplyFilterForDynamicComputationOfCs(velnp_,scaaf_,0.0,dirichtoggle);
  }
  else if (turbmodel_==INPAR::FLUID::scale_similarity or turbmodel_==INPAR::FLUID::scale_similarity_basic)
  {
    switch (scale_sep_)
    {
    case INPAR::FLUID::box_filter:
    {
      // perform filtering
      const Teuchos::RCP<const Epetra_Vector> dirichtoggle = Dirichlet();
      // call only filtering
      Boxf_->ApplyFilter(velnp_,scaaf_,0.0,dirichtoggle);

      // get filtered fields
      filteredvel_->PutScalar(0.0);
      filteredreystr_->PutScalar(0.0);

      Boxf_->GetFilteredVelocity(filteredvel_);
      Boxf_->GetFilteredReynoldsStress(filteredreystr_);
      Boxf_->GetFineScaleVelocity(finescalevel_);

      // store fine-scale velocity
      Boxf_->OutputofFineScaleVel(fsvelaf_);
      break;
    }
    case INPAR::FLUID::algebraic_multigrid_operator:
    {
      const Epetra_Map* dofrowmap = discret_->DofRowMap();
      const Epetra_Map* dofcolmap = discret_->DofColMap();

      RCP<Epetra_Vector> row_filteredveltmp;
      row_filteredveltmp = Teuchos::rcp(new Epetra_Vector(*dofrowmap,true));
      RCP<Epetra_Vector> col_filteredveltmp;
      col_filteredveltmp = Teuchos::rcp(new Epetra_Vector(*dofcolmap,true));

      RCP<Epetra_Vector> row_finescaleveltmp;
      row_finescaleveltmp = Teuchos::rcp(new Epetra_Vector(*dofrowmap,true));
      RCP<Epetra_Vector> col_finescaleveltmp;
      col_finescaleveltmp = Teuchos::rcp(new Epetra_Vector(*dofcolmap,true));

      RCP<Epetra_MultiVector> row_filteredreystretmp;
      row_filteredreystretmp = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap,3,true));
      RCP<Epetra_MultiVector> col_filteredreystretmp;
      col_filteredreystretmp = Teuchos::rcp(new Epetra_MultiVector(*dofcolmap,3,true));
      RCP<Epetra_MultiVector> row_reystretmp;
      row_reystretmp = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap,3,true));
      RCP<Epetra_MultiVector> row_finescalereystretmp;
      row_finescalereystretmp = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap,3,true));

      /*-------------------------------------------------------------------
       * remark:
       * - first, the fine scale velocity is computed
       * - second, we get the coarse scale velocity by subtracting
       *   the fine scale velocity from the resolved velocity
       * - the same procedure is applied to the reynolds stresses
       * - One might think of directly computing the coarse scale quantities
       *   by applying the respective scale-separation operator. However, in
       *   doing so, coarse scale quantities are set to zero on dirichlet
       *   boundaries due to the way the scale-separation operators are
       *   constructed. Setting the fine scale part of the velocity equal zero
       *   is a reasonable choice. However, this strategy is not reasonable for
       *   coarse scale quantities. The coarse scale quantities should rather
       *   contain the exact values. This is ensured by the proposed approach.
       *-------------------------------------------------------------------*/
        //std::cout << "one-step theta" << std::endl;
        // get fine scale velocity
        Sep_->Multiply(false,*velnp_,*row_finescaleveltmp);
        // get filtered or coarse scale velocity
        row_filteredveltmp->Update(1.0,*velnp_,-1.0,*row_finescaleveltmp,0.0);


        // calculate reynoldsstress
        // loop all nodes on this proc
        for (int nid=0;nid<discret_->NumMyRowNodes();++nid)
        {
          // get the node
          DRT::Node* node = discret_->lRowNode(nid);
          // get the dofs of the node
          std::vector<int> dofs= discret_->Dof(node);
          //we only loop over all velocity dofs
          for(int di=0;di<discret_->NumDof(node)-1;++di)
          {
            // get global id of the dof
            int gidi = dofs[di];
            // get local id of the dof
            int lidi = discret_->DofRowMap()->LID(gidi);
            // get the velocity
            double veli=(*velnp_)[lidi];
            for(int dj=0;dj<discret_->NumDof(node)-1;++dj)
            {
              // get the second velocity in the same way
              int gidj =dofs[dj];
              int lidj = discret_->DofRowMap()->LID(gidj);
              double velj=(*velnp_)[lidj];
              // multiply the velocity to get the final component of the reynoldsstress tensor
              double velivelj = veli*velj;
              // store it
              ((*row_reystretmp)(dj))->ReplaceGlobalValues(1,&velivelj,&gidi);
            }
          }
        }

        //get filtered reynoldsstress
        Sep_->Multiply(false,*row_reystretmp,*row_finescalereystretmp);
        row_filteredreystretmp->Update(1.0,*row_reystretmp,-1.0,*row_finescalereystretmp,0.0);

      // export quantities in dofrowmap format to dofcolumnmap format
      LINALG::Export(*row_filteredveltmp,*col_filteredveltmp);
      LINALG::Export(*row_finescaleveltmp,*col_finescaleveltmp);
      LINALG::Export(*row_filteredreystretmp,*col_filteredreystretmp);

      // transfer quantities from dofcolumnmap to nodecolmap
      // filtered velocity and fine scale subgrid velocity
      // loop all nodes on this proc (including ghosted ones)
      for (int nid=0;nid<discret_->NumMyColNodes();++nid)
      {
        // get the node
        DRT::Node* node = discret_->lColNode(nid);
        // get global ids of all dofs of the node
        std::vector<int> dofs= discret_->Dof(node);

        //we only loop over all velocity dofs
        for(int di=0;di<discret_->NumDof(node)-1;++di)
        {
          // get global id of the dof
          int gidi = dofs[di];
          // get local dof id corresponding to the global id
          int lidi = discret_->DofColMap()->LID(gidi);
          // get the values of the dof
          double valvel = (*col_filteredveltmp)[lidi];
          double valfsvel = (*col_finescaleveltmp)[lidi];
          // store them
          int err = 0;
          err += ((*filteredvel_)(di))->ReplaceMyValues(1,&valvel,&nid);
          err += ((*finescalevel_)(di))->ReplaceMyValues(1,&valfsvel,&nid);
          if (err!=0) dserror("dof not on proc");
        }
      }

      // transfer filtered reynoldsstress
      // loop all nodes on this proc (including ghosted ones)
      for (int nid=0;nid<discret_->NumMyColNodes();++nid)
      {
        // get the node
        DRT::Node* node = discret_->lColNode(nid);
        // get global ids of all dofs of the node
        std::vector<int> dofs= discret_->Dof(node);

        //we only loop over all velocity dofs
        for(int di=0;di<discret_->NumDof(node)-1;++di)
        {
          // get global id of the dof
          int gidi = dofs[di];
          // get local dof id corresponding to the global id
          int lidi = discret_->DofColMap()->LID(gidi);
          //loop over all components
          for(int dj=0;dj<discret_->NumDof(node)-1;++dj)
          {
            // get the values
            double val = (*((*col_filteredreystretmp)(dj)))[lidi];
            // and store it
            const int ij = di*3+dj;
            int err = ((*filteredreystr_)(ij))->ReplaceMyValues(1,&val,&nid);
            if (err!=0) dserror("dof not on proc");
          }
        }
      }

      // store fine-scale velocity
      fsvelaf_->Update(1.0,*row_finescaleveltmp,0.0);

      break;
    }
    case INPAR::FLUID::geometric_multigrid_operator:
    {
      dserror("Not available for scale-similarity type models!");
      break;
    }
    default:
    {
      dserror("Unknown filter type!");
      break;
    }
    }
  }
  else if (turbmodel_==INPAR::FLUID::multifractal_subgrid_scales)
  {
    switch (scale_sep_)
    {
    case INPAR::FLUID::box_filter:
    {
      // perform filtering
      const Teuchos::RCP<const Epetra_Vector> dirichtoggle = Dirichlet();
      // call only filtering
      Boxf_->ApplyFilter(velnp_,scaaf_,0.0,dirichtoggle);
      // get fine-scale velocity
      Boxf_->OutputofFineScaleVel(fsvelaf_);

      break;
    }
    case INPAR::FLUID::algebraic_multigrid_operator:
    {
      // get fine-scale part of velocity at time n+alpha_F or n+1
      Sep_->Multiply(false,*velnp_,*fsvelaf_);

      // set fine-scale velocity for parallel nigthly tests
      // separation matrix depends on the number of proc here
      if (turbmodel_==INPAR::FLUID::multifractal_subgrid_scales and
          (DRT::INPUT::IntegralValue<int>(params_->sublist("MULTIFRACTAL SUBGRID SCALES"),"SET_FINE_SCALE_VEL")))
        fsvelaf_->PutScalar(0.01);

      break;
    }
    case INPAR::FLUID::geometric_multigrid_operator:
    {

      ScaleSepGMO_->ApplyScaleSeparation(velnp_,fsvelaf_);

      break;
    }
    default:
    {
      dserror("Unknown filter type!");
      break;
    }
    }

    // set fine-scale vector
    discret_->SetState("fsvelaf",fsvelaf_);
  }
  else if (turbmodel_==INPAR::FLUID::dynamic_vreman)
  {
    // perform filtering
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = Dirichlet();

    Vrem_->ApplyFilterForDynamicComputationOfCv(velnp_,scaaf_,0.0,dirichtoggle);

  }
  else
    dserror("Unknown turbulence model!");

  return;
}
