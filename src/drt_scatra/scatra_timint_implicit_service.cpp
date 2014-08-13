/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_implicit_service.cpp
\brief Service routines of the scalar transport time integration class

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_timint_implicit.H"
#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../linalg/linalg_solver.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_lib/drt_timecurve.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_fluid/fluid_rotsym_periodicbc_utils.H"
#include "../drt_fluid_turbulence/dyn_smag.H"
#include "../drt_fluid_turbulence/dyn_vreman.H"
#include "turbulence_hit_scalar_forcing.H"

// for AVM3 solver:
#include <MLAPI_Workspace.h>
#include <MLAPI_Aggregation.h>

// for printing electrode status to file
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_pstream.H"


/*==========================================================================*
 |                                                                          |
 | public:                                                                  |
 |                                                                          |
 *==========================================================================*/

/*==========================================================================*/
//! @name general framework
/*==========================================================================*/

/*----------------------------------------------------------------------*
 |  calculate mass / heat flux vector                        gjb   04/08|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> SCATRA::ScaTraTimIntImpl::CalcFlux(const bool writetofile, const int num)
{
  switch(writeflux_)
  {
  case INPAR::SCATRA::flux_total_domain:
  case INPAR::SCATRA::flux_diffusive_domain:
  {
    return CalcFluxInDomain(writeflux_);
    break;
  }
  case INPAR::SCATRA::flux_total_boundary:
  case INPAR::SCATRA::flux_diffusive_boundary:
  case INPAR::SCATRA::flux_convective_boundary:
  {
    // calculate normal flux vector field only for the user-defined boundary conditions:
    std::vector<std::string> condnames;
    condnames.push_back("ScaTraFluxCalc");

    return CalcFluxAtBoundary(condnames, writetofile, num);
    break;
  }
  default:
    break;
  }
  // else: we just return a zero vector field (needed for result testing)
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  return Teuchos::rcp(new Epetra_MultiVector(*dofrowmap,3,true));
}

/*----------------------------------------------------------------------*
 |  calculate mass / heat flux vector field in comp. domain    gjb 06/09|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> SCATRA::ScaTraTimIntImpl::CalcFluxInDomain
(const INPAR::SCATRA::FluxType fluxtype)
{
  // get a vector layout from the discretization to construct matching
  // vectors and matrices    local <-> global dof numbering
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // empty vector for (normal) mass or heat flux vectors (always 3D)
  Teuchos::RCP<Epetra_MultiVector> flux = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap,3,true));

  // We have to treat each spatial direction separately
  Teuchos::RCP<Epetra_Vector> fluxx = LINALG::CreateVector(*dofrowmap,true);
  Teuchos::RCP<Epetra_Vector> fluxy = LINALG::CreateVector(*dofrowmap,true);
  Teuchos::RCP<Epetra_Vector> fluxz = LINALG::CreateVector(*dofrowmap,true);

  // we need a vector for the integrated shape functions
  Teuchos::RCP<Epetra_Vector> integratedshapefcts = LINALG::CreateVector(*dofrowmap,true);

  {
    Teuchos::ParameterList eleparams;
    eleparams.set<int>("action",SCATRA::integrate_shape_functions);
    eleparams.set<int>("scatratype",scatratype_);
    // we integrate shape functions for the first numscal_ dofs per node!!
    Epetra_IntSerialDenseVector dofids(7); // make it big enough!
    for(int rr=0;rr<7;rr++)
    {
      dofids(rr) = -1; // do not integrate shape functions for these dofs
    }
    for (std::vector<int>::iterator it = writefluxids_->begin(); it!=writefluxids_->end(); ++it)
    {
      dofids((*it)-1) = (*it); // do not integrate shape functions for these dofs
    }
    eleparams.set("dofids",dofids);

    if (isale_)
      discret_->AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);
    // evaluate fluxes in the whole computational domain
    // (e.g., for visualization of particle path-lines) or L2 projection for better consistency
    discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,integratedshapefcts,Teuchos::null,Teuchos::null);
  }

  // set action for elements
  Teuchos::ParameterList params;
  params.set<int>("action",SCATRA::calc_flux_domain);
  params.set<int>("scatratype",scatratype_);

  // provide velocity field and potentially acceleration/pressure field
  // (export to column map necessary for parallel evaluation)
  discret_->AddMultiVectorToParameterList(params,"convective velocity field",convel_);
  discret_->AddMultiVectorToParameterList(params,"velocity field",vel_);
  discret_->AddMultiVectorToParameterList(params,"acceleration/pressure field",accpre_);

  //provide displacement field in case of ALE
  params.set("isale",isale_);
  if (isale_)
    discret_->AddMultiVectorToParameterList(params,"dispnp",dispnp_);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);

  // evaluate fluxes in the whole computational domain (e.g., for visualization of particle path-lines)
  discret_->Evaluate(params,Teuchos::null,Teuchos::null,fluxx,fluxy,fluxz);

  // insert values into final flux vector for visualization
  // we do not solve a global equation system for the flux values here
  // but perform a lumped mass matrix approach, i.e., dividing by the values of
  // integrated shape functions
  for (int i = 0;i<flux->MyLength();++i)
  {
    const double intshapefct = (*integratedshapefcts)[i];
    // is zero at electric potential dofs
    if (abs(intshapefct) > EPS13)
    {
      flux->ReplaceMyValue(i,0,((*fluxx)[i])/intshapefct);
      flux->ReplaceMyValue(i,1,((*fluxy)[i])/intshapefct);
      flux->ReplaceMyValue(i,2,((*fluxz)[i])/intshapefct);
    }
  }

  // clean up
  discret_->ClearState();

  return flux;
} // SCATRA::ScaTraTimIntImpl::CalcFluxInDomain


/*----------------------------------------------------------------------*
 |  calculate mass / heat normal flux at specified boundaries  gjb 06/09|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> SCATRA::ScaTraTimIntImpl::CalcFluxAtBoundary(
    std::vector<std::string>& condnames,
    const bool writetofile,
    const int num,
    bool fixtotflux)
{
  // The normal flux calculation is based on the idea proposed in
  // GRESHO ET AL.,
  // "THE CONSISTENT GALERKIN FEM FOR COMPUTING DERIVED BOUNDARY
  // QUANTITIES IN THERMAL AND/OR FLUIDS PROBLEMS",
  // INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN FLUIDS, VOL. 7, 371-394 (1987)
  // For the moment, we are lumping the 'boundary mass matrix' instead of solving
  // a small linear system!

  if (fixtotflux ==1) writeflux_= INPAR::SCATRA::flux_total_boundary;

  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // empty vector for (normal) mass or heat flux vectors (always 3D)
  Teuchos::RCP<Epetra_MultiVector> flux = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap,3,true));

  // determine the averaged normal vector field for indicated boundaries
  // used for the output of the normal flux as a vector field
  // is computed only once; for ALE formulation recalculation is necessary
  if ((normals_ == Teuchos::null) or (isale_== true))
    normals_ = ComputeNormalVectors(condnames);

  if (writeflux_==INPAR::SCATRA::flux_convective_boundary)
  {
    // zero out trueresidual vector -> we do not need this info
    trueresidual_->PutScalar(0.0);
  }
  else
  {
    // was the residual already prepared?
    if ((solvtype_!=INPAR::SCATRA::solvertype_nonlinear) and (lastfluxoutputstep_ != step_))
    {
      lastfluxoutputstep_ = step_;

      // For nonlinear problems we already have the actual residual vector
      // from the last convergence test!
      // For linear problems we have to compute this information first, since
      // the residual (w.o. Neumann boundary) has not been computed after the last solve!

      // zero out matrix entries
      sysmat_->Zero();

      // zero out residual vector
      residual_->PutScalar(0.0);

      Teuchos::ParameterList eleparams;
      // action for elements
      eleparams.set<int>("action",SCATRA::calc_mat_and_rhs);
      eleparams.set<int>("scatratype",scatratype_);

      // other parameters that might be needed by the elements
      // the resulting system has to be solved incrementally
      SetElementTimeParameter(true);

      // provide velocity field and potentially acceleration/pressure field
      // (export to column map necessary for parallel evaluation)
      discret_->AddMultiVectorToParameterList(eleparams,"convective velocity field",convel_);
      discret_->AddMultiVectorToParameterList(eleparams,"velocity field",vel_);
      discret_->AddMultiVectorToParameterList(eleparams,"acceleration/pressure field",accpre_);
      // and provide fine-scale velocity for multifractal subgrid-scale modeling only
      if (turbmodel_==INPAR::FLUID::multifractal_subgrid_scales or fssgd_ == INPAR::SCATRA::fssugrdiff_smagorinsky_small)
        discret_->AddMultiVectorToParameterList(eleparams,"fine-scale velocity field",fsvel_);

      //provide displacement field in case of ALE
      if (isale_)
        discret_->AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

      // clear state
      discret_->ClearState();

      // AVM3 separation for incremental solver: get fine-scale part of scalar
      if (incremental_ and
          (fssgd_ != INPAR::SCATRA::fssugrdiff_no or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales))
       AVM3Separation();

      // add element parameters according to time-integration scheme
      // we want to have an incremental solver!!
      AddTimeIntegrationSpecificVectors(true);

      // add element parameters according to specific problem
      AddProblemSpecificParametersAndVectors(eleparams);

      {
        // call standard loop over elements
        discret_->Evaluate(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null);
        discret_->ClearState();
      }

      // scaling to get true residual vector for all time integration schemes
      trueresidual_->Update(ResidualScaling(),*residual_,0.0);

      // undo potential changes
      SetElementTimeParameter();

    } // if ((solvtype_!=INPAR::SCATRA::solvertype_nonlinear) && (lastfluxoutputstep_ != step_))
  } // if (writeflux_==INPAR::SCATRA::flux_convective_boundary)

  // if desired add the convective flux contribution
  // to the trueresidual_ now.
  if((writeflux_==INPAR::SCATRA::flux_total_boundary)
      or (writeflux_==INPAR::SCATRA::flux_convective_boundary))
  {
    if (myrank_==0)
      std::cout<<"Convective flux contribution is added to trueresidual_ vector.\n"
      "Be sure not to address the same boundary part twice!\n Two flux calculation boundaries "
      "should also not share a common node!"<<std::endl;

    // now we evaluate the conditions and separate via ConditionID
    for (unsigned int i=0; i < condnames.size(); i++)
    {
      std::vector<DRT::Condition*> cond;
      discret_->GetCondition(condnames[i],cond);

      discret_->ClearState();
      Teuchos::ParameterList params;

      params.set<int>("action",SCATRA::bd_add_convective_mass_flux);
      params.set<int>("scatratype",scatratype_);

      // add element parameters according to time-integration scheme
      AddTimeIntegrationSpecificVectors();

      // provide velocity field
      // (export to column map necessary for parallel evaluation)
      discret_->AddMultiVectorToParameterList(params,"velocity field",convel_);

      //provide displacement field in case of ALE
      params.set("isale",isale_);
      if (isale_)
        discret_->AddMultiVectorToParameterList(params,"dispnp",dispnp_);

      // call loop over boundary elements and add integrated fluxes to trueresidual_
      discret_->EvaluateCondition(params,trueresidual_,condnames[i]);
      discret_->ClearState();
    }
  }

  // vector for effective flux over all defined boundary conditions
  // maximal number of fluxes  (numscal+1 -> ionic species + potential) is generated
  // for OUTPUT standard -> last entry is not used
  std::vector<double> normfluxsum(numdofpernode_);

  for (unsigned int i=0; i < condnames.size(); i++)
  {
    std::vector<DRT::Condition*> cond;
    discret_->GetCondition(condnames[i],cond);

    // go to the next condition type, if there's nothing to do!
    if (!cond.size()) continue;

    if (myrank_ == 0)
    {
      std::cout<<"Normal fluxes at boundary '"<<condnames[i]<<"':\n"
      <<"+----+-----+-------------------------+------------------+--------------------------+"<<std::endl;
      printf("| ID | DOF | Integral of normal flux | Area of boundary | Mean normal flux density |\n");
    }

    // first, add to all conditions of interest a ConditionID
    for (int condid = 0; condid < (int) cond.size(); condid++)
    {
      // is there already a ConditionID?
      const std::vector<int>*    CondIDVec  = cond[condid]->Get<std::vector<int> >("ConditionID");
      if (CondIDVec)
      {
        if ((*CondIDVec)[0] != condid)
          dserror("Condition %s has non-matching ConditionID",condnames[i].c_str());
      }
      else
      {
        // let's add a ConditionID
        cond[condid]->Add("ConditionID",condid);
      }
    }

    // now we evaluate the conditions and separate via ConditionID
    for (int condid = 0; condid < (int) cond.size(); condid++)
    {
      Teuchos::ParameterList params;

      int addflux=0;
      if(condnames[i]=="ScaTraFluxCalc")
      {
        if((cond[condid])->GetInt("output")==INPAR::SCATRA::fluxeval_alldof)
        {
          // flux for additional dof
          addflux=1;
          params.set<int>("alldof",INPAR::SCATRA::fluxeval_alldof);
        }
      }
      // maximal number of fluxes  (numscal+1 -> ionic species + potential) is used if it is
      // specified in the BC
      const int numfluxeval=numscal_+addflux;

      // calculate integral of shape functions over indicated boundary and it's area
      params.set("boundaryint",0.0);
      params.set<int>("action",SCATRA::bd_integrate_shape_functions);
      params.set<int>("scatratype",scatratype_);

      //provide displacement field in case of ALE
      params.set("isale",isale_);
      if (isale_)
        discret_->AddMultiVectorToParameterList(params,"dispnp",dispnp_);

      // create vector (+ initialization with zeros)
      Teuchos::RCP<Epetra_Vector> integratedshapefunc = LINALG::CreateVector(*dofrowmap,true);

      // call loop over elements
      discret_->ClearState();
      discret_->EvaluateCondition(params,integratedshapefunc,condnames[i],condid);
      discret_->ClearState();

      // maximal number of fluxes  (numscal+1 -> ionic species + potential)
      std::vector<double> normfluxintegral(numfluxeval);

      // insert values into final flux vector for visualization
      int numrownodes = discret_->NumMyRowNodes();
      for (int lnodid = 0; lnodid < numrownodes; ++lnodid )
      {
        DRT::Node* actnode = discret_->lRowNode(lnodid);
        for (int idof = 0; idof < numfluxeval; ++idof)
        {
          int dofgid = discret_->Dof(0,actnode,idof);
          int doflid = dofrowmap->LID(dofgid);

          if ((*integratedshapefunc)[doflid] != 0.0)
          {
            // this is the value of the normal flux density
            double normflux = ((*trueresidual_)[doflid])/(*integratedshapefunc)[doflid];
            // compute integral value for every degree of freedom
            normfluxintegral[idof] += (*trueresidual_)[doflid];

            // care for the slave nodes of rotationally symm. periodic boundary conditions
            double rotangle(0.0);
            bool havetorotate = FLD::IsSlaveNodeOfRotSymPBC(actnode,rotangle);

            // do not insert slave node values here, since they would overwrite the
            // master node values owning the same dof
            // (rotation of slave node vectors is performed later during output)
            if (not havetorotate)
            {
              // for visualization, we plot the normal flux with
              // outward pointing normal vector
              for (int idim = 0; idim < 3; idim++)
              {
                Epetra_Vector* normalcomp = (*normals_)(idim);
                double normalveccomp =(*normalcomp)[lnodid];
                int err=flux->ReplaceMyValue(doflid,idim,normflux*normalveccomp);
                if (err!=0) dserror("Detected error in ReplaceMyValue");
              }
            }
          }
        }
      }

      // get area of the boundary on this proc
      double boundaryint = params.get<double>("boundaryint");

      // care for the parallel case
      std::vector<double> parnormfluxintegral(numfluxeval);
      discret_->Comm().SumAll(&normfluxintegral[0],&parnormfluxintegral[0],numfluxeval);
      double parboundaryint = 0.0;
      discret_->Comm().SumAll(&boundaryint,&parboundaryint,1);

      for (int idof = 0; idof < numfluxeval; ++idof)
      {
        // print out results
        if (myrank_ == 0)
        {
          printf("| %2d | %2d  |       %+10.4E       |    %10.4E    |        %+10.4E       |\n",
              condid,idof,parnormfluxintegral[idof],parboundaryint,parnormfluxintegral[idof]/parboundaryint);
        }
        normfluxsum[idof]+=parnormfluxintegral[idof];
      }

      // statistics section for normfluxintegral
      if (DoBoundaryFluxStatistics())
      {
        // add current flux value to the sum!
        (*sumnormfluxintegral_)[condid] += parnormfluxintegral[0]; // only first scalar!
        int samstep = step_-samstart_+1;

        // dump every dumperiod steps (i.e., write to screen)
        bool dumpit(false);
        if (dumperiod_==0)
          {dumpit=true;}
        else
          {if(samstep%dumperiod_==0) dumpit=true;}

        if(dumpit)
        {
          double meannormfluxintegral = (*sumnormfluxintegral_)[condid]/samstep;
          // dump statistical results
          if (myrank_ == 0)
          {
            printf("| %2d | Mean normal-flux integral (step %5d -- step %5d) :   %12.5E |\n", condid,samstart_,step_,meannormfluxintegral);
          }
        }
      }

      // print out results to file as well (only if really desired)
      if ((myrank_ == 0) and (writetofile==true))
      {
        std::ostringstream temp;
        temp << condid;
        temp << num;
        const std::string fname
        = DRT::Problem::Instance()->OutputControlFile()->FileName()+".boundaryflux_"+condnames[i]+"_"+temp.str()+".txt";

        std::ofstream f;
        if (Step() <= 1)
        {
          f.open(fname.c_str(),std::fstream::trunc);
          f << "#| ID | Step | Time | Area of boundary |";
          for(int idof = 0; idof < numfluxeval; ++idof)
          {
            f<<" Integral of normal flux "<<idof<<" | Mean normal flux density "<<idof<<" |";
          }
          f<<"\n";
        }
        else
          f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

        f << condid << " " << Step() << " " << Time() << " "<< parboundaryint<< " ";
        for (int idof = 0; idof < numfluxeval; ++idof)
        {
          f << parnormfluxintegral[idof] << " "<< parnormfluxintegral[idof]/parboundaryint<< " ";
        }
        f << "\n";
        f.flush();
        f.close();
      } // write to file

    } // loop over condid

    if (myrank_==0)
      std::cout<<"+----+-----+-------------------------+------------------+--------------------------+"<<std::endl;
  }

  // print out the accumulated normal flux over all indicated boundaries
  // the accumulated normal flux over all indicated boundaries is not printed for numscal+1
  if (myrank_ == 0)
  {
    for (int idof = 0; idof < numscal_; ++idof)
    {
      printf("| Sum of all normal flux boundary integrals for scalar %d: %+10.5E             |\n"
          ,idof,normfluxsum[idof]);
    }
    std::cout<<"+----------------------------------------------------------------------------------+"<<std::endl;
  }
  // clean up
  discret_->ClearState();

  return flux;
} // SCATRA::ScaTraTimIntImpl::CalcFluxAtBoundary


/*----------------------------------------------------------------------*
 | calculate initial time derivative of phi at t=t_0           gjb 08/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::CalcInitialPhidt()
{
  // assemble system: M phidt^0 = f^n - K\phi^n - C(u_n)\phi^n
  CalcInitialPhidtAssemble();
  // solve for phidt_0
  CalcInitialPhidtSolve();
} // SCATRA::ScaTraTimIntImpl::CalcInitialPhidt

/*----------------------------------------------------------------------*
 | calculate initial time derivative of phi at t=t_0 (assembly)gjb 08/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::CalcInitialPhidtAssemble()
{
  // time measurement:
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + calc initial phidt");

  double initialtime = Time();
  if (myrank_ == 0)
  {
    IO::cout<<"SCATRA: calculating initial time derivative of phi (step "
    << Step() <<","<<" time "<<initialtime<<")"<< IO::endl;
  }

  // zero out matrix and rhs entries !
  sysmat_->Zero();
  residual_->PutScalar(0.0);

  // evaluate Dirichlet boundary conditions at time t=0
  // the values should match your initial field at the boundary!
  //ApplyDirichletBC(time_,phin_,phidtn_);
  ApplyDirichletBC(time_,phin_,Teuchos::null);

  {
    // evaluate Neumann boundary conditions at time t=0
    neumann_loads_->PutScalar(0.0);
    Teuchos::ParameterList p;
    p.set<int>("scatratype",scatratype_);
    p.set("isale",isale_);
    // provide displacement field in case of ALE
    if (isale_) discret_->AddMultiVectorToParameterList(p,"dispnp",dispnp_);
    discret_->ClearState();
    discret_->EvaluateNeumann(p,*neumann_loads_);
    discret_->ClearState();

    // add potential Neumann boundary condition at time t=0
    // and zero out the residual_ vector!
    residual_->Update(1.0,*neumann_loads_,0.0);
  }

  // call elements to calculate matrix and right-hand-side
  {
    // create the parameters for the discretization
    Teuchos::ParameterList eleparams;

    // action for elements
    eleparams.set<int>("action",SCATRA::calc_initial_time_deriv);
    // set type of scalar transport problem (after preevaluate evaluate, which need scatratype is called)
    eleparams.set<int>("scatratype",scatratype_);

    // add additional parameters
    AddTimeIntegrationSpecificVectors();

    // provide velocity field and potentially acceleration/pressure field
    // (export to column map necessary for parallel evaluation)
    discret_->AddMultiVectorToParameterList(eleparams,"convective velocity field",convel_);
    discret_->AddMultiVectorToParameterList(eleparams,"velocity field",vel_);
    discret_->AddMultiVectorToParameterList(eleparams,"acceleration/pressure field",accpre_);

    //provide displacement field in case of ALE
    if (isale_)
      discret_->AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("hist",zeros_); // we need a zero vector here!!!!
    discret_->SetState("phinp",phin_); // that's the initial field phi0

    // set external volume force (required, e.g., for forced homogeneous isotropic turbulence)
    if (forcing_!=Teuchos::null)
    {
      eleparams.set("forcing",true);
      discret_->SetState("forcing",forcing_);
    }

    // add problem specific time-integration parameters
    AddProblemSpecificParametersAndVectorsForCalcInitialPhiDt(eleparams);

    // call loop over elements
    discret_->Evaluate(eleparams,sysmat_,residual_);
    discret_->ClearState();

    // finalize the complete matrix
    sysmat_->Complete();
  }

  // what's next?
  return;
} // SCATRA::ScaTraTimIntImpl::CalcInitialPhidtAssemble

/*----------------------------------------------------------------------*
 | calculate initial time derivative of phi at t=t_0 (solver) gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::CalcInitialPhidtSolve()
{
  // are we really at step 0?
  //dsassert(step_==0,"Step counter is not 0");

  // We determine phidtn at every node (including boundaries)
  // To be consistent with time integration scheme we do not prescribe
  // any values at Dirichlet boundaries for phidtn !!!
  // apply Dirichlet boundary conditions to system matrix
  // LINALG::ApplyDirichlettoSystem(sysmat_,phidtn_,residual_,phidtn_,*(dbcmaps_->CondMap()));

  // solve for phidtn
  solver_->Solve(sysmat_->EpetraOperator(),phidtn_,residual_,true,true);

  // copy solution also to phidtnp
  phidtnp_->Update(1.0,*phidtn_,0.0);

  // reset the matrix (and its graph!) since we solved
  // a very special problem here that has a different sparsity pattern
  if (DRT::INPUT::IntegralValue<int>(*params_,"BLOCKPRECOND"))
    BlockSystemMatrix()->Reset();
  else
    SystemMatrix()->Reset();
  // reset the solver as well
  solver_->Reset();

  // that's it
  return;
} // SCATRA::ScaTraTimIntImpl::CalcInitialPhidtSolve


/*----------------------------------------------------------------------*
 | compute density from concentration(s)                      gjb 07/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ComputeDensity()
{
  double newdensity(0.0);
  int err(0);

  // loop over all local nodes
  for(int lnodeid=0; lnodeid<discret_->NumMyRowNodes(); lnodeid++)
  {
    // get the processor's local node
    DRT::Node* lnode = discret_->lRowNode(lnodeid);

    // get the degrees of freedom associated with this node
    std::vector<int> nodedofs;
    nodedofs = discret_->Dof(lnode);
    int numdof = nodedofs.size();

    newdensity= 1.0;
    // loop over all species
    for(int k=0; k<numscal_; k++)
    {
      /*
        //                  k=numscal_-1
        //          /       ----                         \
        //         |        \                            |
        // rho_0 * | 1 +    /       alfa_k * (c_k - c_0) |
        //         |        ----                         |
        //          \       k=0                          /
        //
        // For use of molar mass M_k:  alfa_k = M_k/rho_0  !!
       */

      // global and processor's local DOF ID
      const int globaldofid = nodedofs[k];
      const int localdofid = phinp_->Map().LID(globaldofid);
      if (localdofid < 0)
        dserror("localdofid not found in map for given globaldofid");

      // compute contribution to density due to ionic species k
      newdensity += densific_[k]*((*phinp_)[localdofid]-c0_[k]);
    }

    // insert the current density value for this node
    // (has to be at position of el potential/ the position of the last dof!
    const int globaldofid = nodedofs[numdof-1];
    const int localdofid = phinp_->Map().LID(globaldofid);
    if (localdofid < 0)
      dserror("localdofid not found in map for given globaldofid");

    err = densnp_->ReplaceMyValue(localdofid,0,newdensity);

    if (err != 0) dserror("error while inserting a value into densnp_");
  } // loop over all local nodes
  return;
} // SCATRA::ScaTraTimIntImpl::ComputeDensity


/*----------------------------------------------------------------------*
 |  output of some mean values                               gjb   01/09|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::OutputMeanScalars(const int num)
{
  if (outmean_)
  {
    // set scalar values needed by elements
    discret_->ClearState();
    discret_->SetState("phinp",phinp_);
    // set action for elements
    Teuchos::ParameterList eleparams;
    eleparams.set<int>("action",SCATRA::calc_mean_scalars);
    eleparams.set("inverting",false);
    eleparams.set<int>("scatratype",scatratype_);

    //provide displacement field in case of ALE
    if (isale_)
      discret_->AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

    // evaluate integrals of scalar(s) and domain
    Teuchos::RCP<Epetra_SerialDenseVector> scalars
    = Teuchos::rcp(new Epetra_SerialDenseVector(numscal_+1));
    discret_->EvaluateScalars(eleparams, scalars);
    discret_->ClearState();   // clean up

    const double domint = (*scalars)[numscal_];

    // print out values
    if (myrank_ == 0)
    {
      if (scatratype_==INPAR::SCATRA::scatratype_loma)
        std::cout << "Mean scalar: " << std::setprecision (9) << (*scalars)[0]/domint << std::endl;
      else
      {
        std::cout << "Domain integral:          " << std::setprecision (9) << domint << std::endl;
        for (int k = 0; k < numscal_; k++)
        {
          std::cout << "Total concentration (c_"<<k+1<<"): "<< std::setprecision (9) << (*scalars)[k] << std::endl;
          std::cout << "Mean concentration (c_"<<k+1<<"): "<< std::setprecision (9) << (*scalars)[k]/domint << std::endl;
        }
      }
    }

    // print out results to file as well
    if (myrank_ == 0)
    {
      std::stringstream number;
      number << num;
      const std::string fname
      = DRT::Problem::Instance()->OutputControlFile()->FileName()+number.str()+".meanvalues.txt";

      std::ofstream f;
      if (Step() <= 1)
      {
        f.open(fname.c_str(),std::fstream::trunc);
        if (scatratype_==INPAR::SCATRA::scatratype_loma)
          f << "#| Step | Time | Mean scalar |\n";
        else
        {
          f << "#| Step | Time | Domain integral ";
          for (int k = 0; k < numscal_; k++)
          {
            f << "Total concentration (c_"<<k+1<<")| Mean concentration (c_"<<k+1<<") ";
          }
          f << "\n";
        }
      }
      else
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

      f << Step() << " " << Time() << " ";
      if (scatratype_==INPAR::SCATRA::scatratype_loma)
        f << (*scalars)[0]/domint << "\n";
      else
      {
        f << std::setprecision (9) << domint << " ";
        for (int k = 0; k < numscal_; k++)
        {
          f << std::setprecision (9) << (*scalars)[k] << " ";
          f << std::setprecision (9) << (*scalars)[k]/domint << " ";
        }
        f << "\n";
      }
      f.flush();
      f.close();
    }

  } // if(outmean_)

  return;
} // SCATRA::ScaTraTimIntImpl::OutputMeanScalars



// TODO: SCATRA_ELE_CLEANING: BIOFILM
/*----------------------------------------------------------------------*
 | Evaluate surface/interface permeability                     BIOFILM  |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SurfacePermeability(
    Teuchos::RCP<LINALG::SparseOperator> matrix,
    Teuchos::RCP<Epetra_Vector>          rhs)
{
  // time measurement: evaluate condition 'SurfacePermeability'
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'ScaTraCoupling'");

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action",SCATRA::bd_calc_surface_permeability);
  condparams.set<int>("scatratype",scatratype_);

  // provide displacement field in case of ALE
  condparams.set("isale",isale_);
  if (isale_) discret_->AddMultiVectorToParameterList(condparams,"dispnp",dispnp_);

  // set vector values needed by elements
  discret_->ClearState();

  // add element parameters according to time-integration scheme
  AddTimeIntegrationSpecificVectors();

  std::string condstring("ScaTraCoupling");
  discret_->EvaluateCondition(condparams,matrix,Teuchos::null,rhs,Teuchos::null,Teuchos::null,condstring);
  discret_->ClearState();

  matrix->Complete();

  double scaling = 1.0/ResidualScaling();

  rhs->Scale(scaling);
  matrix->Scale(scaling);

  return;
} // SCATRA::ScaTraTimIntImpl::SurfacePermeability


/*----------------------------------------------------------------------*
 |                                                      rasthofer 03/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::PrintScatraType()
{
  std::cout << "# YOUR SCATRA TYPE: ";

  switch (scatratype_)
  {
  case INPAR::SCATRA::scatratype_undefined:
    std::cout << "Undefined" << std::endl;
    break;
  case INPAR::SCATRA::scatratype_condif:
    std::cout << "ConvectionDiffusion" << std::endl;
    break;
  case INPAR::SCATRA::scatratype_loma:
    std::cout << "LowMachNumberFlow" << std::endl;
    break;
  case INPAR::SCATRA::scatratype_elch:
    std::cout << "Elch" << std::endl;
    break;
  case INPAR::SCATRA::scatratype_levelset:
    std::cout << "LevelSet" << std::endl;
    break;
  case INPAR::SCATRA::scatratype_poro:
    std::cout << "Poroscatra" << std::endl;
    break;
  case INPAR::SCATRA::scatratype_advreac:
    std::cout << "Advanced_Reaction" << std::endl;
    break;
  case INPAR::SCATRA::scatratype_pororeac:
    std::cout << "Poro_Scatra_Reaction" << std::endl;
    break;
  case INPAR::SCATRA::scatratype_anisotrop:
    std::cout << "Anisotropic Diffusion" << std::endl;
    break;
  case INPAR::SCATRA::scatratype_cardiac_monodomain:
    std::cout << "Cardiac_Monodomain" << std::endl;
    break;
  default:
    dserror("Fix your scatratype!");
    break;
  }

  return;
}


/*==========================================================================*
 |                                                                          |
 | protected:                                                               |
 |                                                                          |
 *==========================================================================*/

/*==========================================================================*/
// general framework
/*==========================================================================*/

/*----------------------------------------------------------------------*
 | add approximation to flux vectors to a parameter list      gjb 05/10 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AddFluxApproxToParameterList(
    Teuchos::ParameterList& p,
    const enum INPAR::SCATRA::FluxType fluxtype
)
{
  Teuchos::RCP<Epetra_MultiVector> flux
    = CalcFluxInDomain(fluxtype);

  // post_drt_ensight does not support multivectors based on the dofmap
  // for now, I create single vectors that can be handled by the filters

  // get the noderowmap
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  Teuchos::RCP<Epetra_MultiVector> fluxk = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3,true));
  for(int k=0;k<numscal_;++k)
  {
    std::ostringstream temp;
    temp << k;
    std::string name = "flux_phi_"+temp.str();
    for (int i = 0;i<fluxk->MyLength();++i)
    {
      DRT::Node* actnode = discret_->lRowNode(i);
      int dofgid = discret_->Dof(actnode,k);
      fluxk->ReplaceMyValue(i,0,((*flux)[0])[(flux->Map()).LID(dofgid)]);
      fluxk->ReplaceMyValue(i,1,((*flux)[1])[(flux->Map()).LID(dofgid)]);
      fluxk->ReplaceMyValue(i,2,((*flux)[2])[(flux->Map()).LID(dofgid)]);
    }
    discret_->AddMultiVectorToParameterList(p,name,fluxk);
  }
  return;
} // SCATRA::ScaTraTimIntImpl::AddFluxApproxToParameterList

/*----------------------------------------------------------------------*
 | compute outward pointing unit normal vectors at given b.c.  gjb 01/09|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> SCATRA::ScaTraTimIntImpl::ComputeNormalVectors(
    const std::vector<std::string>& condnames
)
{
  // create vectors for x,y and z component of average normal vector field
  // get noderowmap of discretization
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  Teuchos::RCP<Epetra_MultiVector> normal = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3,true));

  discret_->ClearState();

  // set action for elements
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",SCATRA::bd_calc_normal_vectors);
  eleparams.set<int>("scatratype",scatratype_);
  eleparams.set<Teuchos::RCP<Epetra_MultiVector> >("normal vectors",normal);

  //provide displacement field in case of ALE
  eleparams.set("isale",isale_);
  if (isale_)
    discret_->AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

  // loop over all intended types of conditions
  for (unsigned int i=0; i < condnames.size(); i++)
  {
    discret_->EvaluateCondition(eleparams,condnames[i]);
  }

  // clean up
  discret_->ClearState();

  // the normal vector field is not properly scaled up to now. We do this here
  int numrownodes = discret_->NumMyRowNodes();
  Epetra_Vector* xcomp = (*normal)(0);
  Epetra_Vector* ycomp = (*normal)(1);
  Epetra_Vector* zcomp = (*normal)(2);
  for (int i=0; i<numrownodes; ++i)
  {
    double x = (*xcomp)[i];
    double y = (*ycomp)[i];
    double z = (*zcomp)[i];
    double norm = sqrt(x*x + y*y + z*z);
    // form the unit normal vector
    if (norm > EPS15)
    {
      normal->ReplaceMyValue(i,0,x/norm);
      normal->ReplaceMyValue(i,1,y/norm);
      normal->ReplaceMyValue(i,2,z/norm);
    }
  }

  return normal;
} // SCATRA:: ScaTraTimIntImpl::ComputeNormalVectors

/*----------------------------------------------------------------------*
 | evaluate Neumann inflow boundary condition                  vg 03/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ComputeNeumannInflow(
    Teuchos::RCP<LINALG::SparseOperator> matrix,
    Teuchos::RCP<Epetra_Vector>          rhs)
{
  // time measurement: evaluate condition 'Neumann inflow'
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'TransportNeumannInflow'");

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action",SCATRA::bd_calc_Neumann_inflow);
  condparams.set<int>("scatratype",scatratype_);
  condparams.set("isale",isale_);

  // provide velocity field and potentially acceleration/pressure field
  // as well as displacement field in case of ALE
  // (export to column map necessary for parallel evaluation)
  discret_->AddMultiVectorToParameterList(condparams,"convective velocity field",convel_);
  discret_->AddMultiVectorToParameterList(condparams,"velocity field",vel_);
  discret_->AddMultiVectorToParameterList(condparams,"acceleration/pressure field",accpre_);
  if (isale_) discret_->AddMultiVectorToParameterList(condparams,"dispnp",dispnp_);

  // clear state
  discret_->ClearState();

  // add element parameters and vectors according to time-integration scheme
  AddTimeIntegrationSpecificVectors();

  std::string condstring("TransportNeumannInflow");
  discret_->EvaluateCondition(condparams,matrix,Teuchos::null,rhs,Teuchos::null,Teuchos::null,condstring);
  discret_->ClearState();

  return;
} // SCATRA::ScaTraTimIntImpl::ComputeNeumannInflow

/*----------------------------------------------------------------------*
 | evaluate boundary cond. due to convective heat transfer     vg 10/11 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::EvaluateConvectiveHeatTransfer(
    Teuchos::RCP<LINALG::SparseOperator> matrix,
    Teuchos::RCP<Epetra_Vector>          rhs)
{
  // time measurement: evaluate condition 'ThermoConvections'
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'ThermoConvections'");

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action",SCATRA::bd_calc_convective_heat_transfer);
  condparams.set<int>("scatratype",scatratype_);
  condparams.set("isale",isale_);

  // clear state
  discret_->ClearState();

  // add element parameters and vectors according to time-integration scheme
  AddTimeIntegrationSpecificVectors();

  std::string condstring("ThermoConvections");
  discret_->EvaluateCondition(condparams,matrix,Teuchos::null,rhs,Teuchos::null,Teuchos::null,condstring);
  discret_->ClearState();

  return;
} // SCATRA::ScaTraTimIntImpl::EvaluateConvectiveHeatTransfer

/*----------------------------------------------------------------------*
 | write state vectors to Gmsh postprocessing files        henke   12/09|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::OutputToGmsh(
    const int step,
    const double time
    ) const
{
  // turn on/off screen output for writing process of Gmsh postprocessing file
  const bool screen_out = true;

  // create Gmsh postprocessing file
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_scalar", step, 500, screen_out, discret_->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
//  {
//    // add 'View' to Gmsh postprocessing file
//    gmshfilecontent << "View \" " << "Phin \" {" << std::endl;
//    // draw scalar field 'Phindtp' for every element
//    IO::GMSH::ScalarFieldToGmsh(discret_,phin_,gmshfilecontent);
//    gmshfilecontent << "};" << std::endl;
//  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Phinp \" {" << std::endl;
    // draw scalar field 'Phinp' for every element
    IO::GMSH::ScalarFieldToGmsh(discret_,phinp_,gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }
//  {
//    // add 'View' to Gmsh postprocessing file
//    gmshfilecontent << "View \" " << "Phidtn \" {" << std::endl;
//    // draw scalar field 'Phindtn' for every element
//    IO::GMSH::ScalarFieldToGmsh(discret_,phidtn_,gmshfilecontent);
//    gmshfilecontent << "};" << std::endl;
//  }
//  {
//    // add 'View' to Gmsh postprocessing file
//    gmshfilecontent << "View \" " << "Phidtnp \" {" << std::endl;
//    // draw scalar field 'Phindtp' for every element
//    IO::GMSH::ScalarFieldToGmsh(discret_,phidtnp_,gmshfilecontent);
//    gmshfilecontent << "};" << std::endl;
//  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Convective Velocity \" {" << std::endl;
    // draw vector field 'Convective Velocity' for every element
    IO::GMSH::VectorFieldNodeBasedToGmsh(discret_,convel_,gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }
  gmshfilecontent.close();
  if (screen_out) std::cout << " done" << std::endl;
} // ScaTraTimIntImpl::OutputToGmsh

/*----------------------------------------------------------------------*
 |  write mass / heat flux vector to BINIO                   gjb   08/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::OutputFlux(Teuchos::RCP<Epetra_MultiVector> flux)
{
  //safety check
  if (flux == Teuchos::null)
    dserror("Null pointer for flux vector output. Output() called before Update() ??");

  // WORK-AROUND FOR NURBS DISCRETIZATIONS
  // using noderowmap is problematic. Thus, we do not add normal vectors
  // to the scalar result field (scalar information is enough anyway)
  DRT::NURBS::NurbsDiscretization* nurbsdis
  = dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*discret_));
  if(nurbsdis!=NULL)
  {
    Teuchos::RCP<Epetra_Vector> normalflux = Teuchos::rcp(((*flux)(0)),false);
    output_->WriteVector("normalflux", normalflux, IO::DiscretizationWriter::dofvector);
    return; // leave here
  }

  // post_drt_ensight does not support multivectors based on the dofmap
  // for now, I create single vectors that can be handled by the filter

  // get the noderowmap
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  Teuchos::RCP<Epetra_MultiVector> fluxk = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3,true));
  for (std::vector<int>::iterator it = writefluxids_->begin(); it!=writefluxids_->end(); ++it)
  {
    int k=(*it);

    std::ostringstream temp;
    temp << k;
    std::string name = "flux_phi_"+temp.str();
    for (int i = 0;i<fluxk->MyLength();++i)
    {
      DRT::Node* actnode = discret_->lRowNode(i);
      int dofgid = discret_->Dof(0,actnode,k-1);
      // get value for each component of flux vector
      double xvalue = ((*flux)[0])[(flux->Map()).LID(dofgid)];
      double yvalue = ((*flux)[1])[(flux->Map()).LID(dofgid)];
      double zvalue = ((*flux)[2])[(flux->Map()).LID(dofgid)];
      // care for the slave nodes of rotationally symm. periodic boundary conditions
      double rotangle(0.0); //already converted to radians
      bool havetorotate = FLD::IsSlaveNodeOfRotSymPBC(actnode,rotangle);
      if (havetorotate)
      {
        double xvalue_rot = (xvalue*cos(rotangle)) - (yvalue*sin(rotangle));
        double yvalue_rot = (xvalue*sin(rotangle)) + (yvalue*(cos(rotangle)));
        xvalue = xvalue_rot;
        yvalue = yvalue_rot;
      }
      // insert values
      int err = fluxk->ReplaceMyValue(i,0,xvalue);
      err += fluxk->ReplaceMyValue(i,1,yvalue);
      err += fluxk->ReplaceMyValue(i,2,zvalue);
      if (err!=0) dserror("Detected error in ReplaceMyValue");
    }
    if (numdofpernode_==1)
      output_->WriteVector("flux", fluxk, IO::DiscretizationWriter::nodevector);
    else
      output_->WriteVector(name, fluxk, IO::DiscretizationWriter::nodevector);
  }
  // that's it
  return;
} // ScaTraTimIntImpl::OutputFlux


//TODO: SCATRA_ELE_CLEANING: BIOFILM
// This action is not called at element level!!
// Can it be integrated in CalcFlux_domain?
/*----------------------------------------------------------------------*
 |  output of integral reaction                               mc   03/13|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::OutputIntegrReac(const int num)
{
  if (outintegrreac_)
  {
    if (!discret_->Filled()) dserror("FillComplete() was not called");
    if (!discret_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

    // set scalar values needed by elements
    discret_->ClearState();
    discret_->SetState("phinp",phinp_);
    // set action for elements
    Teuchos::ParameterList eleparams;
    eleparams.set<int>("action",SCATRA::calc_integr_reaction);
    eleparams.set<int>("scatratype",scatratype_);
    Teuchos::RCP<std::vector<double> > myreacnp = Teuchos::rcp(new std::vector<double>(numscal_,0.0));
    eleparams.set<Teuchos::RCP<std::vector<double> > >("local reaction integral",myreacnp);

    discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

    myreacnp = eleparams.get<Teuchos::RCP<std::vector<double> > >("local reaction integral");
    // global integral of reaction terms
    std::vector<double> intreacterm(numscal_,0.0);
    for (int k=0; k<numscal_; ++k)
      phinp_->Map().Comm().SumAll(&((*myreacnp)[k]),&intreacterm[k],1);
    discret_->ClearState();   // clean up

    // print out values
    if (myrank_ == 0)
    {
      for (int k = 0; k < numscal_; k++)
      {
        std::cout << "Total reaction (r_"<<k<<"): "<< std::setprecision (9) << intreacterm[k] << std::endl;
      }
    }

    // print out results to file as well
    if (myrank_ == 0)
    {
      std::stringstream number;
      number << num;
      const std::string fname
      = DRT::Problem::Instance()->OutputControlFile()->FileName()+number.str()+".integrreacvalues.txt";

      std::ofstream f;
      if (Step() <= 1)
      {
        f.open(fname.c_str(),std::fstream::trunc);
        f << "#| Step | Time ";
        for (int k = 0; k < numscal_; k++)
        {
          f << "| Total reaction (r_"<<k<<") ";
        }
        f << "\n";

      }
      else
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

      f << Step() << " " << Time() << " ";
      for (int k = 0; k < numscal_; k++)
      {
        f << std::setprecision (9) << intreacterm[k] << " ";
      }
      f << "\n";
      f.flush();
      f.close();
    }

  } // if(outintegrreac_)

  return;
} // SCATRA::ScaTraTimIntImpl::OutputIntegrReac


/*==========================================================================*/
// ELCH
/*==========================================================================*/



/*----------------------------------------------------------------------*
 |  prepare AVM3-based scale separation                        vg 10/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AVM3Preparation()
{
  // time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:            + avm3");

  // create normalized all-scale subgrid-diffusivity matrix
  sysmat_sd_->Zero();

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // action for elements, time factor and stationary flag
  eleparams.set<int>("action",SCATRA::calc_subgrid_diffusivity_matrix);

  // set type of scalar transport problem
  eleparams.set<int>("scatratype",scatratype_);

  //provide displacement field in case of ALE
  eleparams.set("isale",isale_);
  if (isale_)
    discret_->AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

  // add element parameters according to time-integration scheme
  AddTimeIntegrationSpecificVectors();

  // call loop over elements
  discret_->Evaluate(eleparams,sysmat_sd_,residual_);
  discret_->ClearState();

  // finalize the normalized all-scale subgrid-diffusivity matrix
  sysmat_sd_->Complete();

  // apply DBC to normalized all-scale subgrid-diffusivity matrix
  LINALG::ApplyDirichlettoSystem(sysmat_sd_,phinp_,residual_,phinp_,*(dbcmaps_->CondMap()));

  // get normalized fine-scale subgrid-diffusivity matrix
  {
    // this is important to have!!!
    // MLAPI::Init() without arguments uses internally MPI_COMM_WOLRD
    MLAPI::Init();

    // extract the ML parameters
    Teuchos::ParameterList&  mlparams = solver_->Params().sublist("ML Parameters");
    // remark: we create a new solver with ML preconditioner here, since this allows for also using other solver setups
    // to solve the system of equations
    // get the solver number used form the multifractal subgrid-scale model parameter list
    const int scale_sep_solvernumber = extraparams_->sublist("MULTIFRACTAL SUBGRID SCALES").get<int>("ML_SOLVER");
    if (scale_sep_solvernumber != (-1))     // create a dummy solver
    {
      Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(scale_sep_solvernumber),
                                            discret_->Comm(),
                                            DRT::Problem::Instance()->ErrorFile()->Handle()));
      // compute the null space,
      discret_->ComputeNullSpaceIfNecessary(solver->Params(),true);
      // and, finally, extract the ML parameters
      mlparams = solver->Params().sublist("ML Parameters");
    }

    // get toggle vector for Dirchlet boundary conditions
    const Teuchos::RCP<const Epetra_Vector> dbct = DirichletToggle();

    // get nullspace parameters
    double* nullspace = mlparams.get("null space: vectors",(double*)NULL);
    if (!nullspace) dserror("No nullspace supplied in parameter list");
    int nsdim = mlparams.get("null space: dimension",1);

    // modify nullspace to ensure that DBC are fully taken into account
    if (nullspace)
    {
      const int length = sysmat_sd_->OperatorRangeMap().NumMyElements();
      for (int i=0; i<nsdim; ++i)
        for (int j=0; j<length; ++j)
          if ((*dbct)[j]!=0.0) nullspace[i*length+j] = 0.0;
    }

    // get plain aggregation Ptent
    Teuchos::RCP<Epetra_CrsMatrix> crsPtent;
    MLAPI::GetPtent(*sysmat_sd_->EpetraMatrix(),mlparams,nullspace,crsPtent);
    LINALG::SparseMatrix Ptent(crsPtent);

    // compute scale-separation matrix: S = I - Ptent*Ptent^T
    Sep_ = LINALG::Multiply(Ptent,false,Ptent,true);
    Sep_->Scale(-1.0);
    Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(Sep_->RowMap(),false);
    tmp->PutScalar(1.0);
    Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(Sep_->RowMap(),false);
    Sep_->ExtractDiagonalCopy(*diag);
    diag->Update(1.0,*tmp,1.0);
    Sep_->ReplaceDiagonalValues(*diag);

    //complete scale-separation matrix and check maps
    Sep_->Complete(Sep_->DomainMap(),Sep_->RangeMap());
    if (!Sep_->RowMap().SameAs(sysmat_sd_->RowMap())) dserror("rowmap not equal");
    if (!Sep_->RangeMap().SameAs(sysmat_sd_->RangeMap())) dserror("rangemap not equal");
    if (!Sep_->DomainMap().SameAs(sysmat_sd_->DomainMap())) dserror("domainmap not equal");

    // precomputation of unscaled diffusivity matrix:
    // either two-sided S^T*M*S: multiply M by S from left- and right-hand side
    // or one-sided M*S: only multiply M by S from left-hand side
    if (not incremental_)
    {
      Mnsv_ = LINALG::Multiply(*sysmat_sd_,false,*Sep_,false);
      //Mnsv_ = LINALG::Multiply(*Sep_,true,*Mnsv_,false);
    }
  }

  return;
} // ScaTraTimIntImpl::AVM3Preparation

/*----------------------------------------------------------------------*
 |  scaling of AVM3-based subgrid-diffusivity matrix           vg 10/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AVM3Scaling(Teuchos::ParameterList& eleparams)
{
  // time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:            + avm3");

  // some necessary definitions
  int ierr;
  double* sgvsqrt = 0;
  int length = subgrdiff_->MyLength();

  // square-root of subgrid-viscosity-scaling vector for left and right scaling
  sgvsqrt = (double*)subgrdiff_->Values();
  for (int i = 0; i < length; ++i)
  {
    sgvsqrt[i] = sqrt(sgvsqrt[i]);
    int err = subgrdiff_->ReplaceMyValues(1,&sgvsqrt[i],&i);
    if (err != 0) dserror("index not found");
  }

  // get unscaled S^T*M*S from Sep
  sysmat_sd_ = Teuchos::rcp(new LINALG::SparseMatrix(*Mnsv_));

  // left and right scaling of normalized fine-scale subgrid-viscosity matrix
  ierr = sysmat_sd_->LeftScale(*subgrdiff_);
  if (ierr) dserror("Epetra_CrsMatrix::LeftScale returned err=%d",ierr);
  ierr = sysmat_sd_->RightScale(*subgrdiff_);
  if (ierr) dserror("Epetra_CrsMatrix::RightScale returned err=%d",ierr);

  // add the subgrid-viscosity-scaled fine-scale matrix to obtain complete matrix
  Teuchos::RCP<LINALG::SparseMatrix> sysmat = SystemMatrix();
  sysmat->Add(*sysmat_sd_,false,1.0,1.0);

  // set subgrid-diffusivity vector to zero after scaling procedure
  subgrdiff_->PutScalar(0.0);

  return;
} // ScaTraTimIntImpl::AVM3Scaling

/*----------------------------------------------------------------------*
 | construct toggle vector for Dirichlet dofs                  gjb 11/08|
 | assures backward compatibility for avm3 solver; should go away once  |
 *----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> SCATRA::ScaTraTimIntImpl::DirichletToggle()
{
  if (dbcmaps_ == Teuchos::null)
    dserror("Dirichlet map has not been allocated");
  Teuchos::RCP<Epetra_Vector> dirichones = LINALG::CreateVector(*(dbcmaps_->CondMap()),false);
  dirichones->PutScalar(1.0);
  Teuchos::RCP<Epetra_Vector> dirichtoggle = LINALG::CreateVector(*(discret_->DofRowMap()),true);
  dbcmaps_->InsertCondVector(dirichones, dirichtoggle);
  return dirichtoggle;
} // ScaTraTimIntImpl::DirichletToggle


/*========================================================================*/
// turbulence and related
/*========================================================================*/


/*----------------------------------------------------------------------*
 | provide access to the dynamic Smagorinsky filter     rasthofer 08/12 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AccessDynSmagFilter(
  Teuchos::RCP<FLD::DynSmagFilter> dynSmag
)
{
  DynSmag_ = Teuchos::rcp(dynSmag.get(), false);

  // access to the dynamic Smagorinsky class is provided
  // by the fluid scatra coupling algorithm
  // therefore, we only have to add some scatra specific parameters here
  if (turbmodel_ == INPAR::FLUID::dynamic_smagorinsky)
  {
    if (myrank_ == 0)
    {
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
      std::cout << "Dynamic Smagorinsky model: provided access for ScaTra       " << std::endl;
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    }
    DynSmag_->AddScatra(discret_,scatratype_);
  }

  return;
}

/*----------------------------------------------------------------------*
 | provide access to the dynamic Vreman class               krank 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AccessVreman(
  Teuchos::RCP<FLD::Vreman> vrem
)
{
  Vrem_ = Teuchos::rcp(vrem.get(), false);

  // access to the dynamic Vreman class is provided
  // by the fluid scatra coupling algorithm
  // therefore, we only have to add some scatra specific parameters here
  if (turbmodel_ == INPAR::FLUID::dynamic_vreman)
  {
    if (myrank_ == 0)
    {
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
      std::cout << "Dynamic Vreman model: provided access for ScaTra            " << std::endl;
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    }
    Vrem_->AddScatra(discret_,scatratype_);
  }

  return;
}


/*-------------------------------------------------------------------------*
 | calculate mean CsgsB to estimate CsgsD                                  |
 | for multifractal subgrid-scale model                    rasthofer 08/12 |
 *-------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::RecomputeMeanCsgsB()
{
  if (DRT::INPUT::IntegralValue<int>(extraparams_->sublist("MULTIFRACTAL SUBGRID SCALES"),"ADAPT_CSGS_PHI"))
  {
    // mean Cai
    double meanCai = 0.0;

    // variables required for calculation
    // local sums
    double local_sumCai = 0.0;
    double local_sumVol = 0.0;
    // global sums
    double global_sumCai = 0.0;
    double global_sumVol = 0.0;

    // define element matrices and vectors --- dummies
    Epetra_SerialDenseMatrix emat1;
    Epetra_SerialDenseMatrix emat2;
    Epetra_SerialDenseVector evec1;
    Epetra_SerialDenseVector evec2;
    Epetra_SerialDenseVector evec3;

    // generate a parameterlist for communication and control
    Teuchos::ParameterList myparams;
    // action for elements
    myparams.set<int>("action",SCATRA::calc_mean_Cai);
    myparams.set<int>("scatratype",scatratype_);
    // add convective velocity
    discret_->AddMultiVectorToParameterList(myparams,"convective velocity field",convel_);

    // add element parameters according to time-integration scheme
    AddTimeIntegrationSpecificVectors();
    // set thermodynamic pressure in case of loma
    AddProblemSpecificParametersAndVectors(myparams);

    // loop all elements on this proc (excluding ghosted ones)
    for (int nele=0;nele<discret_->NumMyRowElements();++nele)
    {
      // get the element
      DRT::Element* ele = discret_->lRowElement(nele);

      // get element location vector, dirichlet flags and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      ele->LocationVector(*discret_,lm,lmowner,lmstride);

      // call the element evaluate method to integrate functions
      int err = ele->Evaluate(myparams,*discret_,lm,
                            emat1,emat2,
                            evec1,evec2,evec2);
      if (err) dserror("Proc %d: Element %d returned err=%d",myrank_,ele->Id(),err);

      // get contributions of this element and add it up
      local_sumCai += myparams.get<double>("Cai_int");
      local_sumVol += myparams.get<double>("ele_vol");
    }

    // gather contibutions of all procs
    discret_->Comm().SumAll(&local_sumCai,&global_sumCai,1);
    discret_->Comm().SumAll(&local_sumVol,&global_sumVol,1);

    // calculate mean Cai
    meanCai = global_sumCai/global_sumVol;

    //std::cout << "Proc:  " << myrank_ << "  local vol and Cai   "
    //<< local_sumVol << "   " << local_sumCai << "  global vol and Cai   "
    //<< global_sumVol << "   " << global_sumCai << "  mean   " << meanCai << std::endl;

    if (myrank_ == 0)
    {
      std::cout << "\n+--------------------------------------------------------------------------------------------+" << std::endl;
      std::cout << "Multifractal subgrid scales: adaption of CsgsD from near-wall limit of CsgsB:  " << std::setprecision (8) << meanCai << std::endl;
      std::cout << "+--------------------------------------------------------------------------------------------+\n" << std::endl;
    }

    // set meanCai via pre-evaluate call
    myparams.set<int>("action",SCATRA::set_mean_Cai);
    myparams.set<int>("INPAR::SCATRA::ScaTraType",scatratype_);
    myparams.set<double>("meanCai",meanCai);
    // call standard loop over elements
    discret_->Evaluate(myparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  }

  return;
}


/*----------------------------------------------------------------------*
 | calculate intermediate solution                       rasthofer 05/13|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::CalcIntermediateSolution()
{
  if (special_flow_ == "scatra_forced_homogeneous_isotropic_turbulence" and
      extraparams_->sublist("TURBULENCE MODEL").get<std::string>("SCALAR_FORCING")=="isotropic" and
      DRT::INPUT::IntegralValue<INPAR::FLUID::ForcingType>(extraparams_->sublist("TURBULENCE MODEL"),"FORCING_TYPE")
      == INPAR::FLUID::linear_compensation_from_intermediate_spectrum)
  {
    bool activate = true;

    if (activate)
    {
      if (homisoturb_forcing_ == Teuchos::null)
        dserror("Forcing expected!");

      if (myrank_ == 0)
      {
        std::cout << "+------------------------------------------------------------------------+\n";
        std::cout << "|     calculate intermediate solution\n";
        std::cout << "|"<< std::endl;
      }

      // turn off forcing in NonlinearSolve()
      homisoturb_forcing_->ActivateForcing(false);

      // temporary store velnp_ since it will be modified in NonlinearSolve()
      const Epetra_Map* dofrowmap = discret_->DofRowMap();
      Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*dofrowmap,true);
      tmp->Update(1.0,*phinp_,0.0);

      // compute intermediate solution without forcing
      forcing_->PutScalar(0.0); // just to be sure
      NonlinearSolve();

      // calculate required forcing
      homisoturb_forcing_->CalculateForcing(step_);

      // reset velnp_
      phinp_->Update(1.0,*tmp,0.0);

      // recompute intermediate values, since they have been likewise overwritten
      // only for gen.-alpha
      ComputeIntermediateValues();

      homisoturb_forcing_->ActivateForcing(true);

      if (myrank_ == 0)
      {
        std::cout << "|\n";
        std::cout << "+------------------------------------------------------------------------+\n";
        std::cout << "+------------------------------------------------------------------------+\n";
        std::cout << "|" << std::endl;
      }
    }
    else
      // set force to zero
      forcing_->PutScalar(0.0);
  }

  return;
}


/*==========================================================================*/
// Biofilm related
/*==========================================================================*/

/*----------------------------------------------------------------------*/
/* set scatra-fluid displacement vector due to biofilm growth          */
void SCATRA::ScaTraTimIntImpl::SetScFldGrDisp(Teuchos::RCP<Epetra_MultiVector> scatra_fluid_growth_disp)
{
  scfldgrdisp_= scatra_fluid_growth_disp;

  return;
}

/*----------------------------------------------------------------------*/
/* set scatra-structure displacement vector due to biofilm growth          */
void SCATRA::ScaTraTimIntImpl::SetScStrGrDisp(Teuchos::RCP<Epetra_MultiVector> scatra_struct_growth_disp)
{
  scstrgrdisp_= scatra_struct_growth_disp;

  return;
}

/*==========================================================================*/
// Reconstruct gradients
/*==========================================================================*/
//! Calculate the reconstructed nodal gradient of phi
Teuchos::RCP<Epetra_MultiVector> SCATRA::ScaTraTimIntImpl::GetGradientAtNodes(){
  Teuchos::RCP<Epetra_MultiVector> gradPhi = Teuchos::rcp(new Epetra_MultiVector(*(discret_->DofRowMap()), 3, true));

  ReconstructGradientAtNodes(gradPhi);
  return gradPhi;
}


/*----------------------------------------------------------------------*
 | Calculate the reconstructed nodal gradient of phi        winter 04/14|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ReconstructGradientAtNodes(
    Teuchos::RCP<Epetra_MultiVector> gradPhi)
{
  // zero out matrix entries
  sysmat_->Zero();

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action",SCATRA::recon_gradients_at_nodes);

  // set type of scalar transport problem
  eleparams.set<int>("scatratype",scatratype_);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",EvaluationPhi());

  Teuchos::RCP<Epetra_MultiVector> rhs_vectors = Teuchos::rcp(new Epetra_MultiVector(*(discret_->DofRowMap()), 3, true));

  discret_->Evaluate(eleparams,sysmat_,Teuchos::null,Teuchos::rcp((*rhs_vectors)(0),false),Teuchos::rcp((*rhs_vectors)(1),false),Teuchos::rcp((*rhs_vectors)(2),false));

  discret_->ClearState();

  // finalize the complete matrix
  sysmat_->Complete();

  //TODO Remove for-loop and add Belos solver instead
  for (int idim=0; idim<3 ; idim++)
    solver_->Solve(sysmat_->EpetraOperator(),Teuchos::rcp((*gradPhi)(idim),false),Teuchos::rcp((*rhs_vectors)(idim),false),true,true);

}

/*--------------------------------------------------------------------------------------------------------------------*
 | convergence check (only for two-way coupled problems e.g. low-Mach-number flow, two phase flow)           vg 09/11 |
 *--------------------------------------------------------------------------------------------------------------------*/
bool SCATRA::ScaTraTimIntImpl::ConvergenceCheck(int          itnum,
                                                int          itmax,
                                                const double ittol)
{
  bool stopnonliniter = false;

  // define L2-norm of residual, incremental scalar and scalar
  double resnorm_L2(0.0);
  double phiincnorm_L2(0.0);
  double phinorm_L2(0.0);

  if (numscal_>1)
    dserror("Convergence check for two way coupled ScaTra problems, does not support multiple scalar fields. Yet!");

  // for the time being, only one scalar considered for low-Mach-number flow
  /*if (numscal_>1)
  {
    Teuchos::RCP<Epetra_Vector> onlyphi = splitter_->ExtractCondVector(increment_);
    onlyphi->Norm2(&phiincnorm_L2);

    splitter_->ExtractCondVector(phinp_,onlyphi);
    onlyphi->Norm2(&phinorm_L2);
  }
  else*/

  residual_ ->Norm2(&resnorm_L2);
  increment_->Norm2(&phiincnorm_L2);
  phinp_    ->Norm2(&phinorm_L2);

  // check for any INF's and NaN's
  if (std::isnan(resnorm_L2) or
      std::isnan(phiincnorm_L2) or
      std::isnan(phinorm_L2))
    dserror("At least one of the calculated vector norms is NaN.");

  if (abs(std::isinf(resnorm_L2)) or
      abs(std::isinf(phiincnorm_L2)) or
      abs(std::isinf(phinorm_L2)))
    dserror("At least one of the calculated vector norms is INF.");

  // for scalar norm being (close to) zero, set to one
  if (phinorm_L2 < 1e-5) phinorm_L2 = 1.0;

  if (myrank_==0)
  {
    printf("+------------+-------------------+--------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|- residual   -|- scalar-inc -|\n");
    printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |",
         itnum,itmax,ittol,resnorm_L2,phiincnorm_L2/phinorm_L2);
    printf("\n");
    printf("+------------+-------------------+--------------+--------------+\n");
  }

  if ((resnorm_L2 <= ittol) and
      (phiincnorm_L2/phinorm_L2 <= ittol)) stopnonliniter=true;

  // warn if itemax is reached without convergence, but proceed to next timestep
  if ((itnum == itmax) and
      ((resnorm_L2 > ittol) or (phiincnorm_L2/phinorm_L2 > ittol)))
  {
    stopnonliniter=true;
    if (myrank_==0)
    {
      printf("|            >>>>>> not converged in itemax steps!             |\n");
      printf("+--------------------------------------------------------------+\n");
    }
  }

  return stopnonliniter;
} // SCATRA::ScaTraTimIntImplicit::ConvergenceCheck
