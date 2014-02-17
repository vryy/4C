/*----------------------------------------------------------------------*/
/*!
\file fluid_timint_genalpha.cpp
\brief Generalized-alpha time-integration scheme

<pre>
Maintainers: Ursula Rasthofer & Martin Kronbichler
             {rasthofer,kronbichler}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-235
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_timint_genalpha.H"
#include "../drt_fluid_turbulence/scale_sep_gmo.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_turbulence/dyn_smag.H"
#include "../drt_fluid_turbulence/dyn_vreman.H"
#include "../drt_fluid_turbulence/boxfilter.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntGenAlpha::TimIntGenAlpha(
    const Teuchos::RCP<DRT::Discretization>&      actdis,
    const Teuchos::RCP<LINALG::Solver>&           solver,
    const Teuchos::RCP<Teuchos::ParameterList>&   params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output,
    bool                                          alefluid /*= false*/)
: FluidImplicitTimeInt(actdis,solver,params,output,alefluid),
  alphaM_(params_->get<double>("alpha_M")),
  alphaF_(params_->get<double>("alpha_F")),
  gamma_ (params_->get<double>("gamma")),
  startalgo_(false)
  // af-generalized-alpha parameters: gamma_ = 0.5 + alphaM_ - alphaF_
  // (may be reset below when starting algorithm is used)
{

  // starting algorithm only for af-generalized-alpha so far
  // -> check for time-integration scheme and reasonability of number of steps
  if (numstasteps_ > 0)
  {
    if (timealgo_ != INPAR::FLUID::timeint_afgenalpha)
      dserror("no starting algorithm supported for schemes other than af-gen-alpha");
    else startalgo_= true;
    if (numstasteps_>stepmax_)
      dserror("more steps for starting algorithm than steps overall");
  }

  SetElementTimeParameter();

  Initialize();

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                     bk 11/13 |
*-----------------------------------------------------------------------*/
FLD::TimIntGenAlpha::~TimIntGenAlpha()
{
  return;
}


/*----------------------------------------------------------------------*
| Print information about current time step to screen          bk 11/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::PrintTimeStepInfo()
{
  if (myrank_==0)
  {
    switch (timealgo_)
    {
    case INPAR::FLUID::timeint_afgenalpha:
      printf("TIME: %11.4E/%11.4E  DT = %11.4E  Af-Generalized-Alpha  STEP = %4d/%4d \n",
             time_,maxtime_,dta_,step_,stepmax_);
      break;
    case INPAR::FLUID::timeint_npgenalpha:
      printf("TIME: %11.4E/%11.4E  DT = %11.4E  Np-Generalized-Alpha  STEP = %4d/%4d \n",
             time_,maxtime_,dta_,step_,stepmax_);
      break;
    default:
      dserror("parameter out of range: IOP\n");
      break;
    } /* end of switch(timealgo) */
  }
  return;
}

/*----------------------------------------------------------------------*
| calculate pseudo-theta for startalgo_                        bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::SetTheta()
{
  // -------------------------------------------------------------------
  //  For af-generalized-alpha time-integration scheme:
  //  set "pseudo-theta", calculate initial accelerations according to
  //  prescribed Dirichlet values for generalized-alpha time
  //  integration and values at intermediate time steps
  // -------------------------------------------------------------------
  // starting algorithm
  if (startalgo_)
  {
    // use backward-Euler-type parameter combination
    if (step_<=numstasteps_)
    {
      if (myrank_==0)
      {
        std::cout<<"Starting algorithm for Af_GenAlpha active."
            <<"Performing step "<<step_ <<" of "<<numstasteps_
            <<" Backward Euler starting steps"<<std::endl;
      }
      alphaM_ = 1.0;
      alphaF_ = 1.0;
      gamma_  = 1.0;
    }
    else
    {
      // recall original user wish
      alphaM_ = params_->get<double>("alpha_M");
      alphaF_ = params_->get<double>("alpha_F");
      gamma_  = params_->get<double>("gamma");
      // do not enter starting algorithm section in the future
      startalgo_ = false;
    }
  }

  // compute "pseudo-theta" for af-generalized-alpha scheme
  theta_ = alphaF_*gamma_/alphaM_;

  return;
}

/*----------------------------------------------------------------------*
| set old part of right hand side                              bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::SetOldPartOfRighthandside()
{
  /*
     for low-Mach-number flow: distinguish momentum and continuity part
     (continuity part only meaningful for low-Mach-number flow)

      af-generalized-alpha:

                   mom: hist_ = 0.0
                  (con: hist_ = 0.0)

  */

  hist_->PutScalar(0.0);

  return;
}

/*----------------------------------------------------------------------*
 | update acceleration for generalized-alpha time integration  vg 02/09 |
 *----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::GenAlphaUpdateAcceleration()
{

  //                                  n+1     n
  //                               vel   - vel
  //       n+1      n  gamma-1.0      (i)
  //    acc    = acc * --------- + ------------
  //       (i)           gamma      gamma * dt
  //

  // compute factors
  const double fact1 = 1.0/(gamma_*dta_);
  const double fact2 = 1.0 - (1.0/gamma_);

  // consider both velocity and pressure degrees of freedom in case of
  // artificial compressibility or
  // extract and update only velocity degrees of freedom, since in
  // low-Mach-number flow, 'pressure' components are used to store
  // temporal derivatives of scalar/temperature values
  if (physicaltype_ == INPAR::FLUID::artcomp)
  {
    accnp_->Update(fact2,*accn_,0.0);
    accnp_->Update(fact1,*velnp_,-fact1,*veln_,1.0);
  }
  else
  {
    Teuchos::RCP<Epetra_Vector> onlyaccn  = velpressplitter_.ExtractOtherVector(accn_ );
    Teuchos::RCP<Epetra_Vector> onlyveln  = velpressplitter_.ExtractOtherVector(veln_ );
    Teuchos::RCP<Epetra_Vector> onlyvelnp = velpressplitter_.ExtractOtherVector(velnp_);

    Teuchos::RCP<Epetra_Vector> onlyaccnp = Teuchos::rcp(new Epetra_Vector(onlyaccn->Map()));

    onlyaccnp->Update(fact2,*onlyaccn,0.0);
    onlyaccnp->Update(fact1,*onlyvelnp,-fact1,*onlyveln,1.0);

    // copy back into global vector
    LINALG::Export(*onlyaccnp,*accnp_);
  }

} // TimIntGenAlpha::GenAlphaUpdateAcceleration

/*----------------------------------------------------------------------*
 | compute values at intermediate time steps for gen.-alpha    vg 02/09 |
 *----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::GenAlphaIntermediateValues()
{
  // set intermediate values for acceleration and potential temporal
  // derivatives
  //
  //       n+alphaM                n+1                      n
  //    acc         = alpha_M * acc     + (1-alpha_M) *  acc
  //       (i)                     (i)

  // consider both velocity and pressure degrees of freedom in case of
  // artificial compressibility or
  // extract and update only velocity degrees of freedom, since in
  // low-Mach-number flow, 'pressure' components are used to store
  // temporal derivatives of scalar/temperature values
  if (physicaltype_ == INPAR::FLUID::artcomp)
  {
    accam_->Update((alphaM_),*accnp_,(1.0-alphaM_),*accn_,0.0);
  }
  else
  {
    Teuchos::RCP<Epetra_Vector> onlyaccn  = velpressplitter_.ExtractOtherVector(accn_ );
    Teuchos::RCP<Epetra_Vector> onlyaccnp = velpressplitter_.ExtractOtherVector(accnp_);

    Teuchos::RCP<Epetra_Vector> onlyaccam = Teuchos::rcp(new Epetra_Vector(onlyaccnp->Map()));

    onlyaccam->Update((alphaM_),*onlyaccnp,(1.0-alphaM_),*onlyaccn,0.0);

    // copy back into global vector
    LINALG::Export(*onlyaccam,*accam_);
  }

  // set intermediate values for velocity
  //
  //       n+alphaF              n+1                   n
  //      u         = alpha_F * u     + (1-alpha_F) * u
  //       (i)                   (i)
  //
  // and pressure
  //
  //       n+alphaF              n+1                   n
  //      p         = alpha_F * p     + (1-alpha_F) * p
  //       (i)                   (i)
  //
  // note that its af-genalpha with mid-point treatment of the pressure,
  // not implicit treatment as for the genalpha according to Whiting
  velaf_->Update((alphaF_),*velnp_,(1.0-alphaF_),*veln_,0.0);

} // TimIntGenAlpha::GenAlphaIntermediateValues

/*----------------------------------------------------------------------*
| set integration-scheme-specific state                        bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::SetStateTimInt()
{

  discret_->SetState("velaf",velaf_);
  if (timealgo_==INPAR::FLUID::timeint_npgenalpha)
    discret_->SetState("velnp",velnp_);


  return;
}

/*----------------------------------------------------------------------*
| return alphaF_                                               bk 12/13 |
*-----------------------------------------------------------------------*/
double FLD::TimIntGenAlpha::SetTimeFac()
{
  return alphaF_;
}

/*----------------------------------------------------------------------*
| calculate acceleration                                       bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::CalculateAcceleration(
    const Teuchos::RCP<const Epetra_Vector>    velnp,
    const Teuchos::RCP<const Epetra_Vector>    veln,
    const Teuchos::RCP<const Epetra_Vector>    velnm,
    const Teuchos::RCP<const Epetra_Vector>    accn,
    const Teuchos::RCP<Epetra_Vector>          accnp
)
{
  // do nothing: new acceleration is calculated at beginning of next time step

  return;
}

/*----------------------------------------------------------------------*
| set gamma                                                    bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::SetGamma(Teuchos::ParameterList& eleparams)
{
  eleparams.set("gamma"  ,gamma_);
  return;
}

/*----------------------------------------------------------------------*
| scale separation                                             bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::Sep_Multiply()
{
  Sep_->Multiply(false,*velaf_,*fsvelaf_);
  return;
}

/*----------------------------------------------------------------------*
| update velaf_                                                bk 12/13 |
*-----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::UpdateVelafGenAlpha()
{
  velaf_->Update((alphaF_),*velnp_,(1.0-alphaF_),*veln_,0.0);
  return;
}

/*----------------------------------------------------------------------*
 | paraview output of filtered velocity                  rasthofer 02/11|
 *----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::OutputofFilteredVel(
     Teuchos::RCP<Epetra_Vector> outvec,
     Teuchos::RCP<Epetra_Vector> fsoutvec)
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  Teuchos::RCP<Epetra_Vector> row_finescaleveltmp;
  row_finescaleveltmp = Teuchos::rcp(new Epetra_Vector(*dofrowmap,true));

  // get fine scale velocity
  if (scale_sep_ == INPAR::FLUID::algebraic_multigrid_operator)
    Sep_->Multiply(false,*velaf_,*row_finescaleveltmp);
  else
    ScaleSepGMO_->ApplyScaleSeparation(velaf_,row_finescaleveltmp);
  // get filtered or coarse scale velocity
  outvec->Update(1.0,*velaf_,-1.0,*row_finescaleveltmp,0.0);

  fsoutvec->Update(1.0,*row_finescaleveltmp,0.0);

  return;
}

// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------
void FLD::TimIntGenAlpha::SetElementTimeParameter()
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
  if (time_ >0.0)
  {
    eleparams.set("total time",time_-(1-alphaF_)*dta_);
  }
  else
  {
    eleparams.set("total time",time_);
  }
  eleparams.set("alphaF",alphaF_);
  eleparams.set("alphaM",alphaM_);
  eleparams.set("gamma",gamma_);

  // call standard loop over elements
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  return;
}

/*----------------------------------------------------------------------*
| return time integration factor                               bk 12/13 |
*-----------------------------------------------------------------------*/
const double FLD::TimIntGenAlpha::TimIntParam() const
{
  double retval = 0.0;
    // this is the interpolation weight for quantities from last time step
    retval = 1.0 - alphaF_;

  return retval;
}

/*----------------------------------------------------------------------*
 | filtered quantities for classical LES models          rasthofer 02/11|
 *----------------------------------------------------------------------*/
void FLD::TimIntGenAlpha::ApplyScaleSeparationForLES()
{
  if (turbmodel_==INPAR::FLUID::dynamic_smagorinsky)
  {
    // perform filtering and computation of Cs
    // compute averaged values for LijMij and MijMij
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = Dirichlet();

    DynSmag_->ApplyFilterForDynamicComputationOfCs(velaf_,scaaf_,ReturnThermpressaf(),dirichtoggle);
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
        Boxf_->ApplyFilter(velaf_,scaaf_,ReturnThermpressaf(),dirichtoggle);

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

      Teuchos::RCP<Epetra_Vector> row_filteredveltmp;
      row_filteredveltmp = Teuchos::rcp(new Epetra_Vector(*dofrowmap,true));
      Teuchos::RCP<Epetra_Vector> col_filteredveltmp;
      col_filteredveltmp = Teuchos::rcp(new Epetra_Vector(*dofcolmap,true));

      Teuchos::RCP<Epetra_Vector> row_finescaleveltmp;
      row_finescaleveltmp = Teuchos::rcp(new Epetra_Vector(*dofrowmap,true));
      Teuchos::RCP<Epetra_Vector> col_finescaleveltmp;
      col_finescaleveltmp = Teuchos::rcp(new Epetra_Vector(*dofcolmap,true));

      Teuchos::RCP<Epetra_MultiVector> row_filteredreystretmp;
      row_filteredreystretmp = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap,3,true));
      Teuchos::RCP<Epetra_MultiVector> col_filteredreystretmp;
      col_filteredreystretmp = Teuchos::rcp(new Epetra_MultiVector(*dofcolmap,3,true));
      Teuchos::RCP<Epetra_MultiVector> row_reystretmp;
      row_reystretmp = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap,3,true));
      Teuchos::RCP<Epetra_MultiVector> row_finescalereystretmp;
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
      // get fine scale velocity
      Sep_->Multiply(false,*velaf_,*row_finescaleveltmp);
      // get filtered or coarse scale velocity
      row_filteredveltmp->Update(1.0,*velaf_,-1.0,*row_finescaleveltmp,0.0);

      /*-------------------------------------------------------------------
       * idea:
       * - calculate reynolds stress tensor at each node
       * - filter the reynolds stress tensor by multiplication with
       *   large scale separation operator
       *-------------------------------------------------------------------*/
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
          double veli=(*velaf_)[lidi];
          for(int dj=0;dj<discret_->NumDof(node)-1;++dj)
          {
            // get the second velocity in the same way
            int gidj =dofs[dj];
            int lidj = discret_->DofRowMap()->LID(gidj);
            double velj=(*velaf_)[lidj];
            // multiply the velocity to get the final component of the reynoldsstress tensor
            double velivelj = veli*velj;
            // store it
            ((*row_reystretmp)(dj))->ReplaceGlobalValues(1,&velivelj,&gidi);
          }
        }
      }

      //get the filtered reynoldsstress
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
      Boxf_->ApplyFilter(velaf_,scaaf_,ReturnThermpressaf(),dirichtoggle);

      // get fine-scale velocity
      Boxf_->OutputofFineScaleVel(fsvelaf_);

      break;
    }
    case INPAR::FLUID::algebraic_multigrid_operator:
    {
      // get fine-scale part of velocity at time n+alpha_F or n+1

      Sep_->Multiply(false,*velaf_,*fsvelaf_);


      // set fine-scale velocity for parallel nigthly tests
      // separation matrix depends on the number of proc here
      if (turbmodel_==INPAR::FLUID::multifractal_subgrid_scales and
          (DRT::INPUT::IntegralValue<int>(params_->sublist("MULTIFRACTAL SUBGRID SCALES"),"SET_FINE_SCALE_VEL")))
        fsvelaf_->PutScalar(0.01);

      break;
    }
    case INPAR::FLUID::geometric_multigrid_operator:
    {

    ScaleSepGMO_->ApplyScaleSeparation(velaf_,fsvelaf_);


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

    Vrem_->ApplyFilterForDynamicComputationOfCv(velaf_,scaaf_,ReturnThermpressaf(),dirichtoggle);

  }
  else
    dserror("Unknown turbulence model!");

  return;
}

/*----------------------------------------------------------------------*
| extrapolate end point                                        bk 12/13 |
*-----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::TimIntGenAlpha::ExtrapolateEndPoint
(
  Teuchos::RCP<Epetra_Vector> vecn,
  Teuchos::RCP<Epetra_Vector> vecm
)
{
  Teuchos::RCP<Epetra_Vector> vecnp=
  FluidImplicitTimeInt::ExtrapolateEndPoint(vecn,vecm);

  // For gen-alpha extrapolate mid-point quantities to end-point.
  // Otherwise, equilibrium time level is already end-point.

  vecnp->Update((alphaF_-1.0)/alphaF_,*vecn,1.0/alphaF_);

  return vecnp;
}
