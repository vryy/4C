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
#include "scatra_ele_action.H"
#include "scatra_utils.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_fluid/fluid_rotsym_periodicbc_utils.H"
#include "../drt_fluid_turbulence/dyn_smag.H"
#include "turbulence_hit_scalar_forcing.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
// for AVM3 solver:
#include <MLAPI_Workspace.h>
#include <MLAPI_Aggregation.h>
// for printing electrode status to file
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_pstream.H"
//access to the material data (ELCH)
#include "../drt_mat/material.H"
#include "../drt_mat/ion.H"
#include "../drt_mat/matlist.H"
#include "../drt_inpar/inpar_elch.H"


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
    for(int rr=0;rr < numscal_;rr++)
    {
      dofids(rr) = rr;
    }
    for(int rr=numscal_;rr<7;rr++)
    {
      dofids(rr) = -1; // do not integrate shape functions for these dofs
    }
    eleparams.set("dofids",dofids);
    eleparams.set("isale",isale_);
    if (isale_)
      AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);
    // evaluate fluxes in the whole computational domain
    // (e.g., for visualization of particle path-lines) or L2 projection for better consistency
    discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,integratedshapefcts,Teuchos::null,Teuchos::null);
  }

  // set action for elements
  Teuchos::ParameterList params;
  params.set<int>("action",SCATRA::calc_flux_domain);
  params.set<int>("scatratype",scatratype_);
  params.set("frt",frt_);
  params.set<int>("fluxtype",fluxtype);
  params.set<double>("time-step length",dta_);

  // provide velocity field and potentially acceleration/pressure field
  // (export to column map necessary for parallel evaluation)
  AddMultiVectorToParameterList(params,"convective velocity field",convel_);
  AddMultiVectorToParameterList(params,"velocity field",vel_);
  AddMultiVectorToParameterList(params,"acceleration/pressure field",accpre_);

  //provide displacement field in case of ALE
  params.set("isale",isale_);
  if (isale_)
    AddMultiVectorToParameterList(params,"dispnp",dispnp_);

  // parameters for stabilization
  params.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

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
    bool biogrowth)
{
  // The normal flux calculation is based on the idea proposed in
  // GRESHO ET AL.,
  // "THE CONSISTENT GALERKIN FEM FOR COMPUTING DERIVED BOUNDARY
  // QUANTITIES IN THERMAL AND/OR FLUIDS PROBLEMS",
  // INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN FLUIDS, VOL. 7, 371-394 (1987)
  // For the moment, we are lumping the 'boundary mass matrix' instead of solving
  // a small linear system!

  if (biogrowth ==1) writeflux_= INPAR::SCATRA::flux_total_boundary;

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

      // other parameters that might be needed by the elements
      eleparams.set("time-step length",dta_);
      eleparams.set<int>("scatratype",scatratype_);
      eleparams.set("incremental solver",true); // say yes and you get the residual!!
      eleparams.set<int>("form of convective term",convform_);
      eleparams.set<int>("fs subgrid diffusivity",fssgd_);
      //eleparams.set("turbulence model",turbmodel_);
      // set general parameters for turbulent flow
      eleparams.sublist("TURBULENCE MODEL") = extraparams_->sublist("TURBULENCE MODEL");
      // set model-dependent parameters
      eleparams.sublist("SUBGRID VISCOSITY") = extraparams_->sublist("SUBGRID VISCOSITY");
      // and set parameters for multifractal subgrid-scale modeling
      if (turbmodel_==INPAR::FLUID::multifractal_subgrid_scales)
        eleparams.sublist("MULTIFRACTAL SUBGRID SCALES") = extraparams_->sublist("MULTIFRACTAL SUBGRID SCALES");
      eleparams.set("frt",frt_);
      if (scatratype_ == INPAR::SCATRA::scatratype_loma)
        eleparams.set<bool>("update material",(&(extraparams_->sublist("LOMA")))->get<bool>("update material",false));

      // provide velocity field and potentially acceleration/pressure field
      // (export to column map necessary for parallel evaluation)
      AddMultiVectorToParameterList(eleparams,"convective velocity field",convel_);
      AddMultiVectorToParameterList(eleparams,"velocity field",vel_);
      AddMultiVectorToParameterList(eleparams,"acceleration/pressure field",accpre_);
      // and provide fine-scale velocity for multifractal subgrid-scale modeling only
      if (turbmodel_==INPAR::FLUID::multifractal_subgrid_scales or fssgd_ == INPAR::SCATRA::fssugrdiff_smagorinsky_small)
        AddMultiVectorToParameterList(eleparams,"fine-scale velocity field",fsvel_);

      //provide displacement field in case of ALE
      eleparams.set("isale",isale_);
      if (isale_)
        AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

      // parameters for stabilization
      eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

      // clear state
      discret_->ClearState();

      // AVM3 separation for incremental solver: get fine-scale part of scalar
      if (incremental_ and
          (fssgd_ != INPAR::SCATRA::fssugrdiff_no or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales))
       AVM3Separation();

      // we have to perform some dirty action here...
      bool incremental_old = incremental_;
      incremental_ = true;
      // add element parameters according to time-integration scheme
      AddSpecificTimeIntegrationParameters(eleparams);
      // undo
      incremental_ = incremental_old;

      {
        // call standard loop over elements
        discret_->Evaluate(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null);
        discret_->ClearState();
      }

      // scaling to get true residual vector for all time integration schemes
      trueresidual_->Update(ResidualScaling(),*residual_,0.0);

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
      AddSpecificTimeIntegrationParameters(params);

      // provide velocity field
      // (export to column map necessary for parallel evaluation)
      AddMultiVectorToParameterList(params,"velocity field",convel_);

      //provide displacement field in case of ALE
      params.set("isale",isale_);
      if (isale_)
        AddMultiVectorToParameterList(params,"dispnp",dispnp_);

      // call loop over boundary elements and add integrated fluxes to trueresidual_
      discret_->EvaluateCondition(params,trueresidual_,condnames[i]);
      discret_->ClearState();
    }
  }

  // vector for effective flux over all defined boundary conditions
  // maximal number of fluxes  (numscal+1 -> ionic species + potential) is generated
  // for OUTPUT standard -> last entry is not used
  std::vector<double> normfluxsum(numscal_+1);

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
      if(condnames[i]=="ScaTraFluxCalc" and IsElch(scatratype_))
      {
        if((cond[condid])->GetInt("output")==INPAR::SCATRA::fluxeval_alldof)
        {
          // flux for additional dof
          addflux+=1;
          params.set<int>("alldof",INPAR::SCATRA::fluxeval_alldof);
        }
      }
      // maximal number of fluxes  (numscal+1 -> ionic species + potential) is used if it is
      // specified in the BC
      const int numfluxeval = numscal_ + addflux;

      // calculate integral of shape functions over indicated boundary and it's area
      params.set("boundaryint",0.0);
      params.set<int>("action",SCATRA::bd_integrate_shape_functions);
      params.set<int>("scatratype",scatratype_);

      //provide displacement field in case of ALE
      params.set("isale",isale_);
      if (isale_)
        AddMultiVectorToParameterList(params,"dispnp",dispnp_);

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
                flux->ReplaceMyValue(doflid,idim,normflux*normalveccomp);
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
 |  calculate error compared to analytical solution            gjb 10/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::EvaluateErrorComparedToAnalyticalSol()
{
  const INPAR::SCATRA::CalcError calcerr
    = DRT::INPUT::IntegralValue<INPAR::SCATRA::CalcError>(*params_,"CALCERROR");

  switch (calcerr)
  {
  case INPAR::SCATRA::calcerror_no: // do nothing (the usual case)
    break;
  case INPAR::SCATRA::calcerror_Kwok_Wu:
  {
    //   References:

    //   Kwok, Yue-Kuen and Wu, Charles C. K.
    //   "Fractional step algorithm for solving a multi-dimensional
    //   diffusion-migration equation"
    //   Numerical Methods for Partial Differential Equations
    //   1995, Vol 11, 389-397

    //   G. Bauer, V. Gravemeier, W.A. Wall, A 3D finite element approach for the coupled
    //   numerical simulation of electrochemical systems and fluid flow,
    //   International Journal for Numerical Methods in Engineering, 86
    //   (2011) 1339–1359. DOI: 10.1002/nme.3107

    // create the parameters for the error calculation
    Teuchos::ParameterList p;
    p.set<int>("action",SCATRA::calc_error);
    p.set<int>("scatratype",scatratype_);
    p.set("total time",time_);
    p.set("frt",frt_);
    p.set<int>("calcerrorflag",calcerr);
    //provide displacement field in case of ALE
    p.set("isale",isale_);
    if (isale_)
      AddMultiVectorToParameterList(p,"dispnp",dispnp_);
    // parameters for Elch/DiffCond formulation
    if(IsElch(scatratype_))
      p.sublist("DIFFCOND") = extraparams_->sublist("ELCH CONTROL").sublist("DIFFCOND");

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("phinp",phinp_);

    // get (squared) error values
    Teuchos::RCP<Epetra_SerialDenseVector> errors
      = Teuchos::rcp(new Epetra_SerialDenseVector(3));
    discret_->EvaluateScalars(p, errors);
    discret_->ClearState();

    double conerr1 = 0.0;
    double conerr2 = 0.0;
    // for the L2 norm, we need the square root
    if(numscal_==2)
    {
      conerr1 = sqrt((*errors)[0]);
      conerr2 = sqrt((*errors)[1]);
    }
    else if(numscal_==1)
    {
      conerr1 = sqrt((*errors)[0]);
      conerr2 = 0.0;
    }
    else
      dserror("The analytical solution of Kwok and Wu is only defined for two species");

    double poterr  = sqrt((*errors)[2]);

    if (myrank_ == 0)
    {
      printf("\nL2_err for Kwok and Wu (time = %f):\n", time_);
      printf(" concentration1 %15.8e\n concentration2 %15.8e\n potential      %15.8e\n\n",
             conerr1,conerr2,poterr);
    }
#if 0
    if (myrank_ == 0)
    {
      // append error of the last time step to the error file
      if ((step_==stepmax_) or (time_==maxtime_))// write results to file
      {
        ostd::stringstream temp;
        const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
        //const std::string fname = simulation+".relerror";
        const std::string fname = "XXX_kwok_xele.relerror";

        double elelength=0.0;
        if(simulation.find("5x5")!=std::string::npos)
          elelength=0.2;
        else if(simulation.find("10x10")!=std::string::npos)
          elelength=0.1;
        else if(simulation.find("20x20")!=std::string::npos)
          elelength=0.05;
        else if(simulation.find("40x40")!=std::string::npos)
          elelength=0.025;
        else if(simulation.find("50x50")!=std::string::npos)
          elelength=0.02;
        else if(simulation.find("80x80")!=std::string::npos)
          elelength=0.0125;
        else std::cout << "Warning: file name did not allow a evaluation of the element size!!!" << std::endl;

        std::ofstream f;
        f.precision(12);
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
        f << "#| " << simulation << "\n";
        //f << "#| Step | Time | rel. L2-error velocity mag |  rel. L2-error pressure  |\n";
        f << "#| Step | Time | c1 abs. error L2 | c2 abs. error L2 | phi abs. error L2 |  element length  |\n";
        //f << step_ << " " << time_ << " " << velerr/velint << " " << preerr/pint << " "<<"\n";
        f << step_ << " " << time_ << " " << conerr1 << " " << conerr2 << " " <<  poterr << " "<< elelength << "" << "\n";
        f.flush();
        f.close();
      }
    }
#endif
  }
  break;
  case INPAR::SCATRA::calcerror_cylinder:
  {
    //   Reference:
    //   G. Bauer, V. Gravemeier, W.A. Wall, A 3D finite element approach for the coupled
    //   numerical simulation of electrochemical systems and fluid flow,
    //   International Journal for Numerical Methods in Engineering, 2011

    // create the parameters for the error calculation
    Teuchos::ParameterList p;
    p.set<int>("action",SCATRA::calc_error);
    p.set<int>("scatratype",scatratype_);
    p.set("total time",time_);
    p.set("frt",frt_);
    p.set<int>("calcerrorflag",calcerr);
    //provide displacement field in case of ALE
    p.set("isale",isale_);
    if (isale_)
      AddMultiVectorToParameterList(p,"dispnp",dispnp_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("phinp",phinp_);

    // get (squared) error values
    Teuchos::RCP<Epetra_SerialDenseVector> errors
      = Teuchos::rcp(new Epetra_SerialDenseVector(3));
    discret_->EvaluateScalars(p, errors);
    discret_->ClearState();

    // for the L2 norm, we need the square root
    double conerr1 = sqrt((*errors)[0]);
    double conerr2 = sqrt((*errors)[1]);
    double poterr  = sqrt((*errors)[2]);

    if (myrank_ == 0)
    {
      printf("\nL2_err for concentric cylinders (time = %f):\n", time_);
      printf(" concentration1 %15.8e\n concentration2 %15.8e\n potential      %15.8e\n\n",
             conerr1,conerr2,poterr);
    }
  }
  break;
  case INPAR::SCATRA::calcerror_electroneutrality:
  {
    // compute L2 norm of electroneutrality condition

    // create the parameters for the error calculation
    Teuchos::ParameterList p;
    p.set<int>("action",SCATRA::calc_error);
    p.set<int>("scatratype",scatratype_);
    p.set("total time",time_);
    p.set("frt",frt_);
    p.set<int>("calcerrorflag",calcerr);
    //provide displacement field in case of ALE
    p.set("isale",isale_);
    if (isale_)
      AddMultiVectorToParameterList(p,"dispnp",dispnp_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("phinp",phinp_);

    // get (squared) error values
    Teuchos::RCP<Epetra_SerialDenseVector> errors
      = Teuchos::rcp(new Epetra_SerialDenseVector(1));
    discret_->EvaluateScalars(p, errors);
    discret_->ClearState();

    // for the L2 norm, we need the square root
    double err = sqrt((*errors)[0]);

    if (myrank_ == 0)
    {
      printf("\nL2_err for electroneutrality (time = %f):\n", time_);
      printf(" Deviation from ENC: %15.8e\n\n",err);
    }
  }
  break;
  default:
    dserror("Cannot calculate error. Unknown type of analytical test problem"); break;
  }
  return;
} // SCATRA::ScaTraTimIntImpl::EvaluateErrorComparedToAnalyticalSol

/*----------------------------------------------------------------------*
 | calculate initial time derivative of phi at t=t_0           gjb 08/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::CalcInitialPhidt()
{
  // TODO
  // Calculation of initial derivative yields in different results for the uncharged particle and
  // the binary electrolyte solution
  // -> Check calculation procedure of the method

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

//  // are we really at step 0?
//  if(Step() !=0 ) //and scatratype_!= INPAR::SCATRA::scatratype_levelset)
//   ;// dserror("Step counter is not 0");

  // zero out matrix and rhs entries !
  sysmat_->Zero();
  residual_->PutScalar(0.0);

  // evaluate Dirichlet boundary conditions at time t=0
  // the values should match your initial field at the boundary!
  //ApplyDirichletBC(time_,phin_,phidtn_);
  ApplyDirichletBC(time_,phin_,Teuchos::null);

  {
    // TODO: Add coupling condition to potential
    // evaluate Neumann boundary conditions at time t=0
    neumann_loads_->PutScalar(0.0);
    Teuchos::ParameterList p;
    p.set("total time",time_);
    p.set<int>("scatratype",scatratype_);
    p.set("isale",isale_);
    // provide displacement field in case of ALE
    if (isale_) AddMultiVectorToParameterList(p,"dispnp",dispnp_);
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

    // set flag whether to consider stabilization terms or not
    // for backward compatibility we do it like this:
    if (Step() == 0)
      eleparams.set<bool>("onlySGFEM",true); // initial calculation for step 0
    else
      eleparams.set<bool>("onlySGFEM",false);// for reinitialization of level-set fields

    // set type of scalar transport problem
    eleparams.set<int>("scatratype",scatratype_);

    // add additional parameters
    AddSpecificTimeIntegrationParameters(eleparams);

    // other parameters that are needed by the elements
    eleparams.set("incremental solver",true); // we need an incremental formulation for this
    eleparams.set<int>("form of convective term",convform_);
    if (IsElch(scatratype_))
      eleparams.set("frt",frt_); // factor F/RT
    else if (scatratype_==INPAR::SCATRA::scatratype_loma)
    {
      eleparams.set("thermodynamic pressure",thermpressn_);
      eleparams.set("time derivative of thermodynamic pressure",thermpressdtn_);
    }

    // provide velocity field and potentially acceleration/pressure field
    // (export to column map necessary for parallel evaluation)
    AddMultiVectorToParameterList(eleparams,"convective velocity field",convel_);
    AddMultiVectorToParameterList(eleparams,"velocity field",vel_);
    AddMultiVectorToParameterList(eleparams,"acceleration/pressure field",accpre_);

    // set time step length
    eleparams.set("time-step length",dta_);
    // give correct time (override things set by AddSpecificTimeIntegrationParameters)
    eleparams.set("total time",initialtime);

    // ensure backward compatability (override genalpha settings):
    eleparams.set("using generalized-alpha time integration",false);
    // genalpha should only be supported for alpha_F = alpha_M = 1
    // otherwise it is a two-step scheme that needs a start-up algorithm !!

    // other parameters that might be needed by the elements
    eleparams.set<int>("fs subgrid diffusivity",INPAR::SCATRA::fssugrdiff_no); // no fssgd for this Evaluate() call!!!!
    // set general parameters for turbulent flow
    eleparams.sublist("TURBULENCE MODEL") = extraparams_->sublist("TURBULENCE MODEL");
    // set model-dependent parameters
    eleparams.sublist("SUBGRID VISCOSITY") = extraparams_->sublist("SUBGRID VISCOSITY");
    // and set parameters for multifractal subgrid-scale modeling
    eleparams.set("turbulent inflow",turbinflow_);

    if (scatratype_ == INPAR::SCATRA::scatratype_loma)
      eleparams.set<bool>("update material",(&(extraparams_->sublist("LOMA")))->get<bool>("update material",false));

    // set switch for reinitialization
    eleparams.set("reinitswitch",reinitswitch_);

    // parameters for stabilization (here required for material evaluation location)
    eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

    // parameters for Elch/DiffCond formulation
    if(IsElch(scatratype_))
      eleparams.sublist("DIFFCOND") = extraparams_->sublist("ELCH CONTROL").sublist("DIFFCOND");

    //provide displacement field in case of ALE
    eleparams.set("isale",isale_);
    if (isale_)
      AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

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
    eleparams.set("isale",isale_);
    if (isale_)
      AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

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


/*==========================================================================*/
// ELCH
/*==========================================================================*/

/*----------------------------------------------------------------------*
 | compute density from ion concentrations                    gjb 07/09 |
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
    // loop over all ionic species
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
    err = elchdensnp_->ReplaceMyValue(localdofid,0,newdensity);
    if (err != 0) dserror("error while inserting a value into elchdensnp_");

  } // loop over all local nodes
  return;
} // SCATRA::ScaTraTimIntImpl::ComputeDensity

/*----------------------------------------------------------------------*
 |  output of electrode status information to screen         gjb  01/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::OutputElectrodeInfo(
    bool printtoscreen,
    bool printtofile)
{
  // evaluate the following type of boundary conditions:
  std::string condname("ElectrodeKinetics");
  std::vector<DRT::Condition*> cond;
  discret_->GetCondition(condname,cond);

  // leave method, if there's nothing to do!
  if (!cond.size()) return;

  double sum(0.0);

  if ((myrank_ == 0) and printtoscreen)
  {
    std::cout<<"Status of '"<<condname<<"':\n"
    <<"++----+---------------------+------------------+----------------------+--------------------+----------------+----------------+"<<std::endl;
    printf("|| ID |    Total current    | Area of boundary | Mean current density | Mean overpotential | Electrode pot. | Mean Concentr. |\n");
  }

  // first, add to all conditions of interest a ConditionID
  for (int condid = 0; condid < (int) cond.size(); condid++)
  {
    // is there already a ConditionID?
    const std::vector<int>*    CondIDVec  = cond[condid]->Get<std::vector<int> >("ConditionID");
    if (CondIDVec)
    {
      if ((*CondIDVec)[0] != condid)
        dserror("Condition %s has non-matching ConditionID",condname.c_str());
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
    double currtangent(0.0); // this value remains unused here!
    double currresidual(0.0); // this value remains unused here!
    double electrodesurface(0.0); // this value remains unused here!
    double electrodepot(0.0); // this value remains unused here!
    double meanoverpot(0.0); // this value remains unused here!

    OutputSingleElectrodeInfo(
        cond[condid],
        condid,
        printtoscreen,
        printtofile,
        sum,
        currtangent,
        currresidual,
        electrodesurface,
        electrodepot,
        meanoverpot);
  } // loop over condid

  if ((myrank_==0) and printtoscreen)
  {
    std::cout<<"++----+---------------------+------------------+----------------------+--------------------+----------------+----------------+"<<std::endl;
    // print out the net total current for all indicated boundaries
    printf("Net total current over boundary: %10.3E\n\n",sum);
  }

  // clean up
  discret_->ClearState();

  return;
} // ScaTraImplicitTimeInt::OutputElectrodeInfo

/*----------------------------------------------------------------------*
 |  get electrode status for single boundary condition       gjb  11/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::OutputSingleElectrodeInfo(
    DRT::Condition* condition,
    const int condid,
    const bool printtoscreen,
    const bool printtofile,
    double& currentsum,
    double& currtangent,
    double& currresidual,
    double& electrodesurface,
    double& electrodepot,
    double& meanoverpot)
{
  // safety check: is there already a ConditionID?
  const std::vector<int>* CondIDVec  = condition->Get<std::vector<int> >("ConditionID");
  if (not CondIDVec) dserror("Condition has not yet a ConditionID");

  // set vector values needed by elements
  discret_->ClearState();
  // needed for double-layer capacity!
  discret_->SetState("timederivative",phidtnp_);

  // set action for elements
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",SCATRA::bd_calc_elch_electrode_kinetics);
  eleparams.set<int>("scatratype",scatratype_);
  eleparams.set("calc_status",true); // just want to have a status ouput!
  eleparams.set("frt",frt_);

  // parameters for Elch/DiffCond formulation
  if(IsElch(scatratype_))
    eleparams.sublist("DIFFCOND") = extraparams_->sublist("ELCH CONTROL").sublist("DIFFCOND");

  //provide displacement field in case of ALE
  eleparams.set("isale",isale_);
  if (isale_)
    AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

  // Since we just want to have the status ouput for t_{n+1},
  // we have to take care for Gen.Alpha!
  // AddSpecificTimeIntegrationParameters cannot be used since we do not want
  // an evaluation for t_{n+\alpha_f} !!!

  // TODO
  // Warning:
  // Specific time integration parameter are set in the following function.
  // In the case of a genalpha-time integration scheme the solution vector phiaf_ at time n+af
  // is passed to the element evaluation routine. Therefore, the electrode status is evaluate at a
  // different time (n+af) than our output routine (n+1), resulting in slightly different values at the electrode.
  // A different approach is not possible (without major hacks) since the time-integration scheme is
  // necessary to perform galvanostatic simulations, for instance.
  // Think about: double layer effects for genalpha time-integratio scheme

  // add element parameters according to time-integration scheme
  AddSpecificTimeIntegrationParameters(eleparams);

  // values to be computed
  eleparams.set("currentintegral",0.0);
  eleparams.set("boundaryintegral",0.0);
  eleparams.set("overpotentialintegral",0.0);
  eleparams.set("electrodedifferencepotentialintegral",0.0);
  eleparams.set("opencircuitpotentialintegral",0.0);
  eleparams.set("concentrationintegral",0.0);
  eleparams.set("currentderiv",0.0);
  eleparams.set("currentresidual",0.0);

  // would be nice to have a EvaluateScalar for conditions!!!
  discret_->EvaluateCondition(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,"ElectrodeKinetics",condid);

  // get integral of current on this proc
  double currentintegral = eleparams.get<double>("currentintegral");
  // get area of the boundary on this proc
  double boundaryint = eleparams.get<double>("boundaryintegral");
  // get integral of overpotential on this proc
  double overpotentialint = eleparams.get<double>("overpotentialintegral");
  // get integral of overpotential on this proc
  double epdint = eleparams.get<double>("electrodedifferencepotentialintegral");
  // get integral of overpotential on this proc
  double ocpint = eleparams.get<double>("opencircuitpotentialintegral");
  // get integral of reactant concentration on this proc
  double cint = eleparams.get<double>("concentrationintegral");
  // tangent of current w.r.t. electrode potential on this proc
  double currderiv = eleparams.get<double>("currentderiv");
  // get negative current residual (rhs of galvanostatic balance equation)
  double currentresidual = eleparams.get<double>("currentresidual");

  // care for the parallel case
  double parcurrentintegral = 0.0;
  discret_->Comm().SumAll(&currentintegral,&parcurrentintegral,1);
  double parboundaryint = 0.0;
  discret_->Comm().SumAll(&boundaryint,&parboundaryint,1);
  double paroverpotentialint = 0.0;
  discret_->Comm().SumAll(&overpotentialint,&paroverpotentialint,1);
  double parepdint = 0.0;
  discret_->Comm().SumAll(&epdint,&parepdint,1);
  double parocpint = 0.0;
  discret_->Comm().SumAll(&ocpint,&parocpint,1);
  double parcint = 0.0;
  discret_->Comm().SumAll(&cint,&parcint,1);
  double parcurrderiv = 0.0;
  discret_->Comm().SumAll(&currderiv,&parcurrderiv ,1);
  double parcurrentresidual = 0.0;
  discret_->Comm().SumAll(&currentresidual,&parcurrentresidual ,1);

  // access some parameters of the actual condition
  double pot = condition->GetDouble("pot");
  const int curvenum = condition->GetInt("curve");
  if (curvenum>=0)
  {
    const double curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time_);
    // adjust potential at metal side accordingly
    pot *= curvefac;
  }

  // specify some return values
  currentsum += parcurrentintegral; // sum of currents
  currtangent  = parcurrderiv;      // tangent w.r.t. electrode potential on metal side
  currresidual = parcurrentresidual;
  electrodesurface = parboundaryint;
  electrodepot = pot;
  meanoverpot = paroverpotentialint/parboundaryint;

  // clean up
  discret_->ClearState();

  // print out results to screen/file if desired
  if (myrank_ == 0)
  {
    if (printtoscreen) // print out results to screen
    {
      printf("|| %2d |     %10.3E      |    %10.3E    |      %10.3E      |     %10.3E     |   %10.3E   |   %10.3E   |\n",
          condid,parcurrentintegral,parboundaryint,parcurrentintegral/parboundaryint,paroverpotentialint/parboundaryint, pot, parcint/parboundaryint);
    }

    if (printtofile)// write results to file
    {
      std::ostringstream temp;
      temp << condid;
      const std::string fname
      = DRT::Problem::Instance()->OutputControlFile()->FileName()+".electrode_status_"+temp.str()+".txt";

      std::ofstream f;
      if (Step() <= 1)
      {
        f.open(fname.c_str(),std::fstream::trunc);
        f << "#| ID | Step | Time | Total current | Area of boundary | Mean current density | Mean overpotential | Mean electrode pot. diff. | Mean opencircuit pot. | Electrode pot. | Mean Concentr. |\n";
      }
      else
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

      f << condid << "  " << Step() << "  " << Time() << "  " << parcurrentintegral << "  " << parboundaryint
      << "  " << parcurrentintegral/parboundaryint << "  " << paroverpotentialint/parboundaryint << "  " <<
      parepdint/parboundaryint << "  " << parocpint/parboundaryint << "  " << pot << "  "
      << parcint/parboundaryint << "  " <<"\n";
      f.flush();
      f.close();
    }
  } // if (myrank_ == 0)

  return;
} // SCATRA::ScaTraTimIntImpl::OutputSingleElectrodeInfo

/*-------------------------------------------------------------------------*
 | parameters valid for the diffusion-conduction formulation    ehrl 12/12 |
 *-------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ValidParameterDiffCond()
{
  if(myrank_==0)
  {
    // Parameters defined in "ELCH CONTROL"
    Teuchos::ParameterList& elchparams = extraparams_->sublist("ELCH CONTROL");

    if(DRT::INPUT::IntegralValue<INPAR::ELCH::ElchMovingBoundary>(elchparams,"MOVINGBOUNDARY")!=INPAR::ELCH::elch_mov_bndry_no)
      dserror("Moving boundaries are not supported in the ELCH diffusion-conduction framework!!");

    if(DRT::INPUT::IntegralValue<int>(elchparams,"NATURAL_CONVECTION"))
      dserror("Natural convection is not supported in the ELCH diffusion-conduction framework!!");

    if((elchparams.get<int>("MAGNETICFIELD_FUNCNO"))>0)
      dserror("Simulations including a magnetic field are not supported in the ELCH diffusion-conduction framework!!");

    // Parameters defined in "SCALAR TRANSPORT DYNAMIC"
    Teuchos::ParameterList& scatraparams = *params_;

    if((DRT::INPUT::IntegralValue<INPAR::SCATRA::SolverType>(scatraparams,"SOLVERTYPE"))!= INPAR::SCATRA::solvertype_nonlinear)
      dserror("The only solvertype supported by the ELCH diffusion-conduction framework is the non-linar solver!!");

    if((DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatraparams,"VELOCITYFIELD"))!= INPAR::SCATRA::velocity_zero)
      dserror("Convective ion transport is neglected so far!!");

    if((DRT::INPUT::IntegralValue<INPAR::SCATRA::FluxType>(scatraparams,"WRITEFLUX"))== INPAR::SCATRA::flux_diffusive_domain or
       (DRT::INPUT::IntegralValue<INPAR::SCATRA::FluxType>(scatraparams,"WRITEFLUX"))== INPAR::SCATRA::flux_total_domain)
      dserror("This feature is needed -> Think about!!");

    if((DRT::INPUT::IntegralValue<int>(scatraparams,"SKIPINITDER"))!= true)
       std::cout << "Initial phidt - What is it doing in case of diffcond??" << std::endl;

    if((DRT::INPUT::IntegralValue<INPAR::SCATRA::ConvForm>(scatraparams,"CONVFORM"))!= INPAR::SCATRA::convform_convective)
      dserror("Only the convective formulation is supported so far!!");

    if((DRT::INPUT::IntegralValue<int>(scatraparams,"NEUMANNINFLOW"))== true)
      dserror("Neuman inflow BC's are not supported by the ELCH diffusion-conduction framework!!");

    if((DRT::INPUT::IntegralValue<int>(scatraparams,"CONV_HEAT_TRANS"))== true)
      dserror("Convective heat transfer BC's are not supported by the ELCH diffusion-conduction framework!!");

    if((DRT::INPUT::IntegralValue<INPAR::SCATRA::FSSUGRDIFF>(scatraparams,"FSSUGRDIFF"))!= INPAR::SCATRA::fssugrdiff_no)
      dserror("Subgrid diffusivity is not supported by the ELCH diffusion-conduction framework!!");

    if((DRT::INPUT::IntegralValue<int>(scatraparams,"BLOCKPRECOND"))== true)
          dserror("Block preconditioner is not supported so far!!");

    if((DRT::INPUT::IntegralValue<INPAR::FLUID::MeshTying>(scatraparams,"MESHTYING"))!= INPAR::FLUID::no_meshtying)
      dserror("Subgrid diffusivity is not supported by the ELCH diffusion-conduction framework!!");

    // Parameters defined in "SCALAR TRANSPORT DYNAMIC"
    Teuchos::ParameterList& scatrastabparams = params_->sublist("STABILIZATION");

    if((DRT::INPUT::IntegralValue<INPAR::SCATRA::StabType>(scatrastabparams,"STABTYPE"))!= INPAR::SCATRA::stabtype_no_stabilization)
      dserror("No stabilization is necessary for solving the ELCH diffusion-conduction framework!!");

    if((DRT::INPUT::IntegralValue<INPAR::SCATRA::TauType>(scatrastabparams,"DEFINITION_TAU"))!= INPAR::SCATRA::tau_zero)
      dserror("No stabilization is necessary for solving the ELCH diffusion-conduction framework!!");

    if((DRT::INPUT::IntegralValue<INPAR::SCATRA::EvalTau>(scatrastabparams,"EVALUATION_TAU"))!= INPAR::SCATRA::evaltau_integration_point)
      dserror("Evaluation of stabilization parameter only at Gauss points!!");

    if((DRT::INPUT::IntegralValue<INPAR::SCATRA::EvalMat>(scatrastabparams,"EVALUATION_MAT"))!= INPAR::SCATRA::evalmat_integration_point)
      dserror("Evaluation of material only at Gauss points!!");

    if((DRT::INPUT::IntegralValue<INPAR::SCATRA::Consistency>(scatrastabparams,"CONSISTENCY"))!= INPAR::SCATRA::consistency_no)
          dserror("Consistence formulation is not in the ELCH diffusion-conduction framework!!");

    if((DRT::INPUT::IntegralValue<int>(scatrastabparams,"SUGRVEL"))== true)
          dserror("Subgrid velocity is not incoperated in the ELCH diffusion-conduction framework!!");

    if((DRT::INPUT::IntegralValue<int>(scatrastabparams,"ASSUGRDIFF"))== true)
          dserror("Subgrid diffusivity is not incoperated in the ELCH diffusion-conduction framework!!");
  }

  return;
}


/*==========================================================================*/
// low-Mach-number flow
/*==========================================================================*/

/*----------------------------------------------------------------------*
 | set initial thermodynamic pressure                          vg 07/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetInitialThermPressure()
{
  // get thermodynamic pressure and gas constant from material parameters
  // (if no temperature equation, zero values are returned)
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",SCATRA::get_material_parameters);
  eleparams.set<int>("scatratype",scatratype_);
  eleparams.set("isale",isale_);
  // provide displacement field in case of ALE
  if (isale_) AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  thermpressn_ = eleparams.get("thermodynamic pressure", 98100.0);

  // initialize also value at n+1
  // (computed if not constant, otherwise prescribed value remaining)
  thermpressnp_ = thermpressn_;

  // initialize time derivative of thermodynamic pressure at n+1 and n
  // (computed if not constant, otherwise remaining zero)
  thermpressdtnp_ = 0.0;
  thermpressdtn_  = 0.0;

  // compute values at intermediate time steps
  // (only for generalized-alpha time-integration scheme)
  // -> For constant thermodynamic pressure, this is done here once and
  // for all simulation time.
  ComputeThermPressureIntermediateValues();

  return;
} // SCATRA::ScaTraTimIntImpl::SetInitialThermPressure

/*----------------------------------------------------------------------*
 | compute initial time derivative of thermodynamic pressure   vg 07/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ComputeInitialThermPressureDeriv()
{
  // define element parameter list
  Teuchos::ParameterList eleparams;

  // DO THIS BEFORE PHINP IS SET (ClearState() is called internally!!!!)
  // compute flux approximation and add it to the parameter list
  AddFluxApproxToParameterList(eleparams,INPAR::SCATRA::flux_diffusive_domain);

  // set scalar vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phin_);

  // provide velocity field and potentially acceleration/pressure field
  // (export to column map necessary for parallel evaluation)
  AddMultiVectorToParameterList(eleparams,"convective velocity field",convel_);
  AddMultiVectorToParameterList(eleparams,"velocity field",vel_);
  AddMultiVectorToParameterList(eleparams,"acceleration/pressure field",accpre_);

  // provide displacement field in case of ALE
  eleparams.set("isale",isale_);
  if (isale_) AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

  // set parameters for element evaluation
  eleparams.set<int>("action",SCATRA::calc_domain_and_bodyforce);
  eleparams.set<int>("scatratype",scatratype_);
  eleparams.set("total time",0.0);

  // variables for integrals of domain and bodyforce
  Teuchos::RCP<Epetra_SerialDenseVector> scalars
    = Teuchos::rcp(new Epetra_SerialDenseVector(2));

  // evaluate domain and bodyforce integral
  discret_->EvaluateScalars(eleparams, scalars);

  // get global integral values
  double pardomint  = (*scalars)[0];
  double parbofint  = (*scalars)[1];

  // set action for elements
  eleparams.set<int>("action",SCATRA::bd_calc_loma_therm_press);

  // variables for integrals of normal velocity and diffusive flux
  double normvelint      = 0.0;
  double normdifffluxint = 0.0;
  eleparams.set("normal velocity integral",normvelint);
  eleparams.set("normal diffusive flux integral",normdifffluxint);

  // evaluate velocity-divergence and diffusive (minus sign!) flux on boundaries
  // We may use the flux-calculation condition for calculation of fluxes for
  // thermodynamic pressure, since it is usually at the same boundary.
  std::vector<std::string> condnames;
  condnames.push_back("ScaTraFluxCalc");
  for (unsigned int i=0; i < condnames.size(); i++)
  {
    discret_->EvaluateCondition(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,condnames[i]);
  }

  // get integral values on this proc
  normvelint      = eleparams.get<double>("normal velocity integral");
  normdifffluxint = eleparams.get<double>("normal diffusive flux integral");

  // get integral values in parallel case
  double parnormvelint      = 0.0;
  double parnormdifffluxint = 0.0;
  discret_->Comm().SumAll(&normvelint,&parnormvelint,1);
  discret_->Comm().SumAll(&normdifffluxint,&parnormdifffluxint,1);

  // clean up
  discret_->ClearState();

  // compute initial time derivative of thermodynamic pressure
  // (with specific heat ratio fixed to be 1.4)
  const double shr = 1.4;
  thermpressdtn_ = (-shr*thermpressn_*parnormvelint
                    + (shr-1.0)*(-parnormdifffluxint+parbofint))/pardomint;

  // set time derivative of thermodynamic pressure at n+1 equal to the one at n
  // for following evaluation of intermediate values
  thermpressdtnp_ = thermpressdtn_;

  // compute values at intermediate time steps
  // (only for generalized-alpha time-integration scheme)
  ComputeThermPressureIntermediateValues();

  return;
} // SCATRA::ScaTraTimIntImpl::ComputeInitialThermPressureDeriv

/*----------------------------------------------------------------------*
 | compute initial total mass in domain                        vg 01/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ComputeInitialMass()
{
  // set scalar values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phin_);
  // set action for elements
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",SCATRA::calc_mean_scalars);
  eleparams.set<int>("scatratype",scatratype_);
  // inverted scalar values are required here
  eleparams.set("inverting",true);

  //provide displacement field in case of ALE
  eleparams.set("isale",isale_);
  if (isale_) AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

  // evaluate integral of inverse temperature
  Teuchos::RCP<Epetra_SerialDenseVector> scalars
    = Teuchos::rcp(new Epetra_SerialDenseVector(numscal_+1));
  discret_->EvaluateScalars(eleparams, scalars);
  discret_->ClearState();   // clean up

  // compute initial mass times gas constant: R*M_0 = int(1/T_0)*tp
  initialmass_ = (*scalars)[0]*thermpressn_;

  // print out initial total mass
  if (myrank_ == 0)
  {
    std::cout << std::endl;
    std::cout << "+--------------------------------------------------------------------------------------------+" << std::endl;
    std::cout << "Initial total mass in domain (times gas constant): " << initialmass_ << std::endl;
    std::cout << "+--------------------------------------------------------------------------------------------+" << std::endl;
  }

  return;
} // SCATRA::ScaTraTimIntImpl::ComputeInitialMass

/*----------------------------------------------------------------------*
 | compute thermodynamic pressure from mass conservation       vg 01/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ComputeThermPressureFromMassCons()
{
  // set scalar values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);
  // set action for elements
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",SCATRA::calc_mean_scalars);
  eleparams.set<int>("scatratype",scatratype_);
  // inverted scalar values are required here
  eleparams.set("inverting",true);

  //provide displacement field in case of ALE
  eleparams.set("isale",isale_);
  if (isale_) AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

  // evaluate integral of inverse temperature
  Teuchos::RCP<Epetra_SerialDenseVector> scalars
    = Teuchos::rcp(new Epetra_SerialDenseVector(numscal_+1));
  discret_->EvaluateScalars(eleparams, scalars);
  discret_->ClearState();   // clean up

  // compute thermodynamic pressure: tp = R*M_0/int(1/T)
  thermpressnp_ = initialmass_/(*scalars)[0];

  // print out thermodynamic pressure
  if (myrank_ == 0)
  {
    std::cout << std::endl;
    std::cout << "+--------------------------------------------------------------------------------------------+" << std::endl;
    std::cout << "Thermodynamic pressure from mass conservation: " << thermpressnp_ << std::endl;
    std::cout << "+--------------------------------------------------------------------------------------------+" << std::endl;
  }

  // compute time derivative of thermodynamic pressure at time step n+1
  ComputeThermPressureTimeDerivative();

  // compute values at intermediate time steps
  // (only for generalized-alpha time-integration scheme)
  ComputeThermPressureIntermediateValues();

  return;
} // SCATRA::ScaTraTimIntImpl::ComputeThermPressureFromMassCons

/*----------------------------------------------------------------------*
 | convergence check (only for low-Mach-number flow)           vg 09/11 |
 *----------------------------------------------------------------------*/
bool SCATRA::ScaTraTimIntImpl::ConvergenceCheck(int          itnum,
                                                int          itmax,
                                                const double ittol)
{
  bool stopnonliniter = false;

  // define L2-norm of residual, incremental scalar and scalar
  double resnorm_L2(0.0);
  double phiincnorm_L2(0.0);
  double phinorm_L2(0.0);

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
} // SCATRA::ScaTraTimIntImpl::ConvergenceCheck

/*----------------------------------------------------------------------*
 | Evaluate surface/interface permeability                              |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SurfacePermeability(
    RCP<LINALG::SparseOperator> matrix,
    RCP<Epetra_Vector>          rhs)
{
  // time measurement: evaluate condition 'SurfacePermeability'
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'ScaTraCoupling'");

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action",SCATRA::bd_calc_surface_permeability);
  condparams.set<int>("scatratype",scatratype_);
  condparams.set("incremental solver",incremental_);

  // provide displacement field in case of ALE
  condparams.set("isale",isale_);
  if (isale_) AddMultiVectorToParameterList(condparams,"dispnp",dispnp_);

  // set vector values needed by elements
  discret_->ClearState();

  // add element parameters according to time-integration scheme
  AddSpecificTimeIntegrationParameters(condparams);

  std::string condstring("ScaTraCoupling");
  discret_->EvaluateCondition(condparams,matrix,Teuchos::null,rhs,Teuchos::null,Teuchos::null,condstring);
  discret_->ClearState();

  matrix->Complete();

  double scaling = 1.0/ResidualScaling();

  rhs->Scale(scaling);
  matrix->Scale(scaling);

  return;
} // SCATRA::ScaTraTimIntImpl::SurfacePermeability


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
    AddMultiVectorToParameterList(p,name,fluxk);
  }
  return;
} // SCATRA::ScaTraTimIntImpl::AddFluxApproxToParameterList

/*----------------------------------------------------------------------*
 | compute outward pointing unit normal vectors at given b.c.  gjb 01/09|
 *----------------------------------------------------------------------*/
RCP<Epetra_MultiVector> SCATRA::ScaTraTimIntImpl::ComputeNormalVectors(
    const std::vector<std::string>& condnames
)
{
  // create vectors for x,y and z component of average normal vector field
  // get noderowmap of discretization
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  RCP<Epetra_MultiVector> normal = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3,true));

  discret_->ClearState();

  // set action for elements
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",SCATRA::bd_calc_normal_vectors);
  eleparams.set<int>("scatratype",scatratype_);
  eleparams.set<RCP<Epetra_MultiVector> >("normal vectors",normal);

  //provide displacement field in case of ALE
  eleparams.set("isale",isale_);
  if (isale_)
    AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

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
    RCP<LINALG::SparseOperator> matrix,
    RCP<Epetra_Vector>          rhs)
{
  // time measurement: evaluate condition 'Neumann inflow'
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'TransportNeumannInflow'");

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action",SCATRA::bd_calc_Neumann_inflow);
  condparams.set<int>("scatratype",scatratype_);
  condparams.set("incremental solver",incremental_);
  condparams.set("isale",isale_);

  // provide velocity field and potentially acceleration/pressure field
  // as well as displacement field in case of ALE
  // (export to column map necessary for parallel evaluation)
  AddMultiVectorToParameterList(condparams,"convective velocity field",convel_);
  AddMultiVectorToParameterList(condparams,"velocity field",vel_);
  AddMultiVectorToParameterList(condparams,"acceleration/pressure field",accpre_);
  if (isale_) AddMultiVectorToParameterList(condparams,"dispnp",dispnp_);

  // clear state
  discret_->ClearState();

  // add element parameters and vectors according to time-integration scheme
  AddSpecificTimeIntegrationParameters(condparams);

  std::string condstring("TransportNeumannInflow");
  discret_->EvaluateCondition(condparams,matrix,Teuchos::null,rhs,Teuchos::null,Teuchos::null,condstring);
  discret_->ClearState();

  return;
} // SCATRA::ScaTraTimIntImpl::ComputeNeumannInflow

/*----------------------------------------------------------------------*
 | evaluate boundary cond. due to convective heat transfer     vg 10/11 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::EvaluateConvectiveHeatTransfer(
    RCP<LINALG::SparseOperator> matrix,
    RCP<Epetra_Vector>          rhs)
{
  // time measurement: evaluate condition 'ThermoConvections'
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'ThermoConvections'");

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action",SCATRA::bd_calc_convective_heat_transfer);
  condparams.set<int>("scatratype",scatratype_);
  condparams.set("incremental solver",incremental_);
  condparams.set("isale",isale_);

  // clear state
  discret_->ClearState();

  // add element parameters and vectors according to time-integration scheme
  AddSpecificTimeIntegrationParameters(condparams);

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
void SCATRA::ScaTraTimIntImpl::OutputFlux(RCP<Epetra_MultiVector> flux)
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
    RCP<Epetra_Vector> normalflux = Teuchos::rcp(((*flux)(0)),false);
    output_->WriteVector("normalflux", normalflux, IO::DiscretizationWriter::dofvector);
    return; // leave here
  }

  // post_drt_ensight does not support multivectors based on the dofmap
  // for now, I create single vectors that can be handled by the filter

  // get the noderowmap
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  Teuchos::RCP<Epetra_MultiVector> fluxk = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3,true));
  for (std::vector<int>::iterator it = writefluxids_.begin(); it!=writefluxids_.end(); ++it)
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
    if (numscal_==1)
      output_->WriteVector("flux", fluxk, IO::DiscretizationWriter::nodevector);
    else
      output_->WriteVector(name, fluxk, IO::DiscretizationWriter::nodevector);
  }
  // that's it
  return;
} // ScaTraTimIntImpl::OutputFlux

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
    eleparams.set("time-step length",dta_);
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
 | perform setup of natural convection applications (ELCH)    gjb 07/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetupElchNatConv()
{
  // loads densification coefficients and the initial mean concentration

  // only required for ELCH with natural convection
  if (extraparams_->isSublist("ELCH CONTROL"))
  {
    if (DRT::INPUT::IntegralValue<int>(extraparams_->sublist("ELCH CONTROL"),"NATURAL_CONVECTION") == true)
    {
      // allocate denselch_ with *dofrowmap and initialize it
      const Epetra_Map* dofrowmap = discret_->DofRowMap();
      elchdensnp_ = LINALG::CreateVector(*dofrowmap,true);
      elchdensnp_->PutScalar(1.0);

      // Calculate the initial mean concentration value
      if (numscal_ < 1) dserror("Error since numscal = %d. Not allowed since < 1",numscal_);
      c0_.resize(numscal_);

      discret_->ClearState();
      discret_->SetState("phinp",phinp_);
      // set action for elements
      Teuchos::ParameterList eleparams;
      eleparams.set<int>("action",SCATRA::calc_mean_scalars);
      eleparams.set<int>("scatratype",scatratype_);
      eleparams.set("inverting",false);

      //provide displacement field in case of ALE
      eleparams.set("isale",isale_);
      if (isale_)
        AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

      // evaluate integrals of concentrations and domain
      Teuchos::RCP<Epetra_SerialDenseVector> scalars
      = Teuchos::rcp(new Epetra_SerialDenseVector(numscal_+1));
      discret_->EvaluateScalars(eleparams, scalars);
      discret_->ClearState();   // clean up

      // calculate mean_concentration
      const double domint  = (*scalars)[numscal_];
      for(int k=0;k<numscal_;k++)
      {
        c0_[k] = (*scalars)[k]/domint;
      }

      //initialization of the densification coefficient vector
      densific_.resize(numscal_);
      DRT::Element*   element = discret_->lRowElement(0);
      RCP<MAT::Material>  mat = element->Material();

      if (mat->MaterialType() == INPAR::MAT::m_matlist)
      {
        const MAT::MatList* actmat = static_cast<const MAT::MatList*>(mat.get());

        for (int k = 0;k<numscal_;++k)
        {
          const int matid = actmat->MatID(k);
          Teuchos::RCP<const MAT::Material> singlemat = actmat->MaterialById(matid);

          if (singlemat->MaterialType() == INPAR::MAT::m_ion)
          {
            const MAT::Ion* actsinglemat = static_cast<const MAT::Ion*>(singlemat.get());
            densific_[k] = actsinglemat->Densification();
            if (densific_[k] < 0.0) dserror("received negative densification value");
          }
          else
            dserror("material type is not allowed");
        }
      }
      if (mat->MaterialType() == INPAR::MAT::m_ion) // for a single species calculation
      {
        const MAT::Ion* actmat = static_cast<const MAT::Ion*>(mat.get());
        densific_[0] = actmat->Densification();
        if (densific_[0] < 0.0) dserror("received negative densification value");
        if (numscal_ > 1) dserror("Single species calculation but numscal = %d > 1",numscal_);
      }
    }
  }

  return;
} // ScaTraTimIntImpl::SetupElchNatConv

/*----------------------------------------------------------------------*
 | define the magnetic field                                  gjb 05/10 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetMagneticField(const int funcno)
{
  if (funcno > 0)
  {
    int err(0);
    const int numdim = 3; // the magnetic field is always 3D

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode = discret_->lRowNode(lnodeid);
      for(int index=0;index<numdim;++index)
      {
        double value = DRT::Problem::Instance()->Funct(funcno-1).Evaluate(index,lnode->X(),time_,NULL);
        // no time-dependency included, yet!
        err = magneticfield_->ReplaceMyValue(lnodeid, index, value);
        if (err!=0) dserror("error while inserting a value into magneticfield_");
      }
    }
  }
  return;

} // ScaTraImplicitTimeInt::SetMagneticField

/*----------------------------------------------------------------------*
 | calculate initial electric potential field at t=t_0         gjb 04/10|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::CalcInitialPotentialField()
{
  if (IsElch(scatratype_))
  {
    if (DRT::INPUT::IntegralValue<int>(*params_,"INITPOTCALC"))
    {
      // time measurement:
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + calc initial potential field");
      if (myrank_ == 0)
        std::cout<<"SCATRA: calculating initial field for electric potential"<<std::endl;

      // are we really at step 0?
      dsassert(step_==0,"Step counter is not 0");

      // construct intermediate vectors
      const Epetra_Map* dofrowmap = discret_->DofRowMap();
      RCP<Epetra_Vector> rhs = LINALG::CreateVector(*dofrowmap,true);
      RCP<Epetra_Vector> phi0 = LINALG::CreateVector(*dofrowmap,true);

      // zero out matrix entries
      sysmat_->Zero();

      // evaluate Dirichlet boundary conditions at time t=0
      // the values should match your initial field at the boundary!
      //ApplyDirichletBC(time_,phin_,phidtn_);
      ApplyDirichletBC(time_,phin_,Teuchos::null);

      // ToDo:
      // contributions due to Neumann b.c. or ElectrodeKinetics b.c.
      // have to be summed up here, and applied
      // as a current flux condition at the potential field!

      // evaluate Neumann boundary conditions at time t=0
      /*{
        neumann_loads_->PutScalar(0.0);
        Teuchos::ParameterList p;
        p.set("total time",time_);
        p.set<int>("scatratype",scatratype_);
        p.set("isale",isale_);
        discret_->ClearState();
        discret_->EvaluateNeumann(p,*neumann_loads_);
        discret_->ClearState();

        // add potential Neumann boundary condition at time t=0
        // and zero out the residual_ vector!
        // residual_->Update(1.0,*neumann_loads_,0.0);
      }*/

      // call elements to calculate matrix and right-hand-side
      {

        // create the parameters for the discretization
        Teuchos::ParameterList eleparams;

        // action for elements
        eleparams.set<int>("action",SCATRA::calc_elch_initial_potential);

        // set type of scalar transport problem
        eleparams.set<int>("scatratype",scatratype_);

        // factor F/RT
        eleparams.set("frt",frt_);

        // parameters for stabilization
        eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

        // parameters for Elch/DiffCond formulation
        if(IsElch(scatratype_))
          eleparams.sublist("DIFFCOND") = extraparams_->sublist("ELCH CONTROL").sublist("DIFFCOND");

        //provide displacement field in case of ALE
        eleparams.set("isale",isale_);
        if (isale_)
          AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

        // set vector values needed by elements
        discret_->ClearState();
        discret_->SetState("phi0",phin_);

        // call loop over elements
        discret_->Evaluate(eleparams,sysmat_,rhs);
        discret_->ClearState();

        // finalize the complete matrix
        sysmat_->Complete();
      }

      // apply Dirichlet boundary conditions to system matrix
      LINALG::ApplyDirichlettoSystem(sysmat_,phi0,rhs,phi0,*(dbcmaps_->CondMap()));

      // solve
      solver_->Solve(sysmat_->EpetraOperator(),phi0,rhs,true,true);

      // copy solution of initial potential field to the solution vectors
      Teuchos::RCP<Epetra_Vector> onlypot = splitter_->ExtractCondVector(phi0);
      // insert values into the whole solution vectors
      splitter_->InsertCondVector(onlypot,phinp_);
      splitter_->InsertCondVector(onlypot,phin_);

      // reset the matrix (and its graph!) since we solved
      // a very special problem here that has a different sparsity pattern
      if (DRT::INPUT::IntegralValue<int>(*params_,"BLOCKPRECOND"))
        BlockSystemMatrix()->Reset();
      else
        SystemMatrix()->Reset();
    }
  }
  // go on!
  return;
} // ScaTraTimIntImpl::CalcInitialPotentialField

/*----------------------------------------------------------------------*
 |  calculate conductivity of electrolyte solution             gjb 07/09|
 *----------------------------------------------------------------------*/
Epetra_SerialDenseVector SCATRA::ScaTraTimIntImpl::ComputeConductivity()
{
  // we perform the calculation on element level hiding the material access!
  // the initial concentration distribution has to be uniform to do so!!

  // create the parameters for the elements
  Teuchos::ParameterList p;
  p.set<int>("action",SCATRA::calc_elch_conductivity);
  p.set<int>("scatratype",scatratype_);
  p.set("frt",frt_);

  // parameters for Elch/DiffCond formulation
  if(IsElch(scatratype_))
    p.sublist("DIFFCOND") = extraparams_->sublist("ELCH CONTROL").sublist("DIFFCOND");

  //provide displacement field in case of ALE
  p.set("isale",isale_);
  if (isale_)
    AddMultiVectorToParameterList(p,"dispnp",dispnp_);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);

  // pointer to current element
  DRT::Element* actele = discret_->lRowElement(0);

  // get element location vector, dirichlet flags and ownerships
  std::vector<int> lm;  // location vector
  std::vector<int> lmowner;  // processor which owns DOFs
  std::vector<int> lmstride;  // nodal block sizes in element matrices

  actele->LocationVector(*discret_,lm,lmowner,lmstride);

  // define element matrices and vectors
  // -- which are empty and unused, just to satisfy element Evaluate()
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // define element vector
  Epetra_SerialDenseVector sigma(numscal_+1);

  // call the element evaluate method of the first row element
  int err = actele->Evaluate(p,*discret_,lm,elematrix1,elematrix2,sigma,elevector2,elevector3);
  if (err) dserror("error while computing conductivity");
  discret_->ClearState();

  return sigma;
} // ScaTraTimIntImpl::ComputeConductivity

/*----------------------------------------------------------------------*
 | apply galvanostatic control                                gjb 11/09 |
 *----------------------------------------------------------------------*/
bool SCATRA::ScaTraTimIntImpl::ApplyGalvanostaticControl()
{
  // for galvanostatic ELCH applications we have to adjust the
  // applied cell voltage and continue Newton-Raphson iterations until
  // we reach the desired value for the electric current.

  // leave method, if there's nothing to do!
  if (extraparams_->isSublist("ELCH CONTROL") == false) return true;

  if (DRT::INPUT::IntegralValue<int>(extraparams_->sublist("ELCH CONTROL"),"GALVANOSTATIC"))
  {
    std::vector<DRT::Condition*> cond;
    discret_->GetCondition("ElectrodeKinetics",cond);
    if (!cond.empty())
    {
      const unsigned condid_cathode = extraparams_->sublist("ELCH CONTROL").get<int>("GSTATCONDID_CATHODE");
      const unsigned condid_anode = extraparams_->sublist("ELCH CONTROL").get<int>("GSTATCONDID_ANODE");
      int gstatitemax = (extraparams_->sublist("ELCH CONTROL").get<int>("GSTATITEMAX"));
      double gstatcurrenttol = (extraparams_->sublist("ELCH CONTROL").get<double>("GSTATCURTOL"));
      const int curvenum = extraparams_->sublist("ELCH CONTROL").get<int>("GSTATCURVENO");
      const double tol = extraparams_->sublist("ELCH CONTROL").get<double>("GSTATCONVTOL");
      const double effective_length = extraparams_->sublist("ELCH CONTROL").get<double>("GSTAT_LENGTH_CURRENTPATH");

      const double potold = cond[condid_cathode]->GetDouble("pot");
      double potnew = potold;
      double actualcurrent(0.0);
      double currtangent(0.0);
      double currresidual(0.0);
      double electrodesurface(0.0);
      //Assumption: Residual at BV1 is the negative of the value at BV2, therefore only the first residual is calculated
      double newtonrhs(0.0);

      // for all time integration schemes, compute the current value for phidtnp
      // this is needed for evaluating charging currents due to double-layer capacity
      // This may only be called here and not inside OutputSingleElectrodeInfo!!!!
      // Otherwise you modify your output to file called during Output()
      ComputeTimeDerivative();

      double targetcurrent = DRT::Problem::Instance()->Curve(curvenum-1).f(time_);
      double timefac = 1.0/ResidualScaling();

      double currtangent_anode(0.0);
      double currtangent_cathode(0.0);
      double potinc_ohm(0.0);
      double electrodepot(0.0);
      double meanoverpot(0.0);

      double potdiffbulk(0.0);
      double potdiffcell(0.0);

      // loop over all BV
      // degenerated to a loop over 2 (user-specified) BV conditions
      for (unsigned int icond = 0; icond < cond.size(); icond++)
      {
        // consider only the specified electrode kinetics boundaries!
        if ((icond != condid_cathode)and((icond != condid_anode)))
          continue;

        actualcurrent = 0.0;
        currtangent = 0.0;
        currresidual = 0.0;
        electrodepot = 0.0;
        meanoverpot = 0.0;

        // note: only the potential at the boundary with id condid_cathode will be adjusted!
        OutputSingleElectrodeInfo(
            cond[icond],
            icond,
            false,
            false,
            actualcurrent,
            currtangent,
            currresidual,
            electrodesurface,
            electrodepot,
            meanoverpot
            );

        // bulk voltage loss = V_A - eta_A - V_C + eta_C
        // cell voltage loss = V_A - V_C
        if (icond==condid_anode)
        {
          potdiffcell += electrodepot;
          potdiffbulk += (electrodepot - (meanoverpot));
        }
        if (icond==condid_cathode)
        {
          potdiffcell -= electrodepot;
          potdiffbulk -= (electrodepot - (meanoverpot));
        }

        // store the tangent for later usage
        if (icond==condid_cathode)
          currtangent_cathode=currtangent;
        if (icond==condid_anode)
          currtangent_anode=currtangent;

        if (icond==condid_cathode)
        {
          //Assumption: Residual at BV1 is the negative of the value at BV2, therefore only the first residual is calculated
          // newtonrhs = -residual, with the definition:  residual := timefac*(-I + I_target)
          newtonrhs = + currresidual - (timefac*targetcurrent); // newtonrhs is stored only from cathode!
          if (myrank_==0)
          {
            std::cout<<"\nGALVANOSTATIC MODE:\n";
            std::cout<<"iteration "<<gstatnumite_<<" / "<<gstatitemax<<std::endl;
            std::cout<<"  actual reaction current = "<<std::scientific<<actualcurrent<<std::endl;
            std::cout<<"  required total current  = "<<targetcurrent<<std::endl;
            std::cout<<"  negative residual (rhs) = "<<newtonrhs<<std::endl<<std::endl;
          }

          if (gstatnumite_ > gstatitemax)
          {
            if (myrank_==0) std::cout<< std::endl <<"  --> maximum number iterations reached. Not yet converged!"<<std::endl<<std::endl;
            return true; // we proceed to next time step
          }
          else if (abs(newtonrhs)< gstatcurrenttol)
          {
            if (myrank_==0) std::cout<< std::endl <<"  --> Newton-RHS-Residual is smaller than " << gstatcurrenttol<< "!" <<std::endl<<std::endl;
            return true; // we proceed to next time step
          }
          // electric potential increment of the last iteration
          else if ((gstatnumite_ > 1) and (abs(gstatincrement_)< (1+abs(potold))*tol)) // < ATOL + |pot|* RTOL
          {
            if (myrank_==0) std::cout<< std::endl <<"  --> converged: |"<<gstatincrement_<<"| < "<<(1+abs(potold))*tol<<std::endl<<std::endl;
            return true; // galvanostatic control has converged
          }

          // update applied electric potential
          // potential drop ButlerVolmer conditions (surface ovepotential) and in the electrolyte (ohmic overpotential) are conected in parallel:
          //
          // I_0 = I_BV1 = I_ohmic = I_BV2
          // R(I_target, I) = R_BV1(I_target, I) = R_ohmic(I_target, I) = -R_BV2(I_target, I)
          // delta E_0 = delta U_BV1 + delta U_ohmic - (delta U_BV2)
          // => delta E_0 = (R_BV1(I_target, I)/J) + (R_ohmic(I_target, I)/J) - (-R_BV2(I_target, I)/J)

          potinc_ohm=(-1.0*effective_length*newtonrhs)/(sigma_(numscal_)*timefac*electrodesurface);

          // print additional information
          if (myrank_==0)
          {
            std::cout<< "  area (condid " << icond <<")               = " << electrodesurface << std::endl;
            std::cout<< "  actualcurrent - targetcurrent = " << (actualcurrent-targetcurrent) << std::endl;
            std::cout<< "  conductivity                  = " << sigma_(numscal_) << std::endl;
          }
        }

        // safety check
        if (abs(currtangent)<EPS13)
          dserror("Tangent in galvanostatic control is near zero: %lf",currtangent);

      }
      // end loop over electrode kinetics

      // actual potential difference is used to calculate the current path length
      // -> it is possible to compute the new ohmic potential step
      //    without the input parameter GSTAT_LENGTH_CURRENTPATH
      //
      // for only one BV boundary condition, keep the original approach using GSTAT_LENGTH_CURRENTPATH
      if (cond.size()>=2)
      {
        if(myrank_==0)
        {
          std::cout << std::endl <<"  cell potential difference = "<<potdiffcell<<std::endl;
          std::cout<<"  bulk potential difference = "<<potdiffbulk<<std::endl;
        }
        if (abs(actualcurrent) > EPS10)
        {
          if(myrank_==0)
          {
            std::cout<<"  Defined GSTAT_LENGTH_CURRENTPATH:   "<< effective_length << std::endl;
            std::cout<<"  Suggested GSTAT_LENGTH_CURRENTPATH: "<< (potdiffbulk/(actualcurrent))*(sigma_(numscal_)*electrodesurface)<<std::endl;
            std::cout<<"  dV/dI =  "<< (-1.0)*potdiffbulk/actualcurrent << std::endl;
            std::cout<<"  potinc_ohm (CURRENTPATH based): "<< potinc_ohm <<std::endl;
            std::cout<<"  New guess for potinc_ohm:       "<< (-1.0)*(potdiffbulk*newtonrhs)/(timefac*actualcurrent)<<std::endl;
          }
          potinc_ohm = (-1.0)*(potdiffbulk*newtonrhs)/(timefac*(actualcurrent));
        }
        else
        {
          potinc_ohm = 0.0;
        }
      }

      // Newton step:  Jacobian * \Delta pot = - Residual
      std::cout << "currtangent_cathode:  " << currtangent_cathode << std::endl;
      const double potinc_cathode = newtonrhs/((-1)*currtangent_cathode);
      double potinc_anode = 0.0;
      if (abs(currtangent_anode)>EPS13) // anode surface overpotential is optional
        potinc_anode = newtonrhs/((-1)*currtangent_anode);
      gstatincrement_ = (potinc_cathode+potinc_anode+potinc_ohm);
      // update electric potential
      potnew += gstatincrement_;

      // print info to screen
      if (myrank_==0)
      {
        std::cout<<std::endl<< "  ohmic overpotential                        = " << potinc_ohm << std::endl;
        std::cout<< "  overpotential increment cathode (condid " << condid_cathode <<") = " << potinc_cathode << std::endl;
        if (abs(potinc_anode)>EPS12) // prevents output if an anode is not considered
          std::cout<< "  overpotential increment anode   (condid " << condid_anode <<") = " << potinc_anode << std::endl;

        std::cout<< "  total increment for potential              = " << potinc_cathode+potinc_anode+potinc_ohm << std::endl;
        std::cout<< std::endl;
        std::cout<< "  old electrode potential (condid "<<condid_cathode <<") = "<<potold<<std::endl;
        std::cout<< "  new electrode potential (condid "<<condid_cathode <<") = "<<potnew<<std::endl<<std::endl;
      }
      // replace potential value of the boundary condition (on all processors)
      cond[condid_cathode]->Add("pot",potnew);
      gstatnumite_++;
      return false; // not yet converged -> continue Newton iteration with updated potential
    }
  }
  return true; //default

} // end ApplyGalvanostaticControl()

/*----------------------------------------------------------------------*
 | evaluate contribution of electrode kinetics to eq. system  gjb 02/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::EvaluateElectrodeKinetics(
    RCP<LINALG::SparseOperator> matrix,
    RCP<Epetra_Vector>          rhs
)
{
  // time measurement: evaluate condition 'ElectrodeKinetics'
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'ElectrodeKinetics'");

  discret_->ClearState();

  // create an parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action",SCATRA::bd_calc_elch_electrode_kinetics);
  condparams.set<int>("scatratype",scatratype_);
  condparams.set("frt",frt_); // factor F/RT
  condparams.set("isale",isale_);

  // parameters for Elch/DiffCond formulation
  if(IsElch(scatratype_))
    condparams.sublist("DIFFCOND") = extraparams_->sublist("ELCH CONTROL").sublist("DIFFCOND");

  if (isale_)   //provide displacement field in case of ALE
    AddMultiVectorToParameterList(condparams,"dispnp",dispnp_);

  // add element parameters and set state vectors according to time-integration scheme
  AddSpecificTimeIntegrationParameters(condparams);

  std::string condstring("ElectrodeKinetics");
  // evaluate ElectrodeKinetics conditions at time t_{n+1} or t_{n+alpha_F}
  discret_->EvaluateCondition(condparams,matrix,Teuchos::null,rhs,Teuchos::null,Teuchos::null,condstring);
  discret_->ClearState();

  return;
} // ScaTraTimIntImpl::EvaluateElectrodeKinetics

/*----------------------------------------------------------------------*
 | Nernst BC                                                 ehrl 03/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AdaptDC(
    Teuchos::RCP<LINALG::SparseOperator> matrix,
    RCP<Epetra_Vector>          residual,
    RCP<Epetra_Vector>          phinp
)
{
  discret_->ClearState();

  // create an parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action",SCATRA::bd_calc_elch_adapt_DC);
  condparams.set<int>("scatratype",scatratype_);
  condparams.set("frt",frt_); // factor F/RT
  condparams.set("isale",isale_);

  // parameters for Elch/DiffCond formulation
  if(IsElch(scatratype_))
    condparams.sublist("DIFFCOND") = extraparams_->sublist("ELCH CONTROL").sublist("DIFFCOND");

  // add element parameters and set state vectors according to time-integration scheme
  // we need here concentration at t+np
  discret_->SetState("phinp",phinp_);

  std::string condstring("ElectrodeKinetics");
  // evaluate ElectrodeKinetics conditions at time t_{n+1} or t_{n+alpha_F}
  // phinp (view to phinp)
  discret_->EvaluateCondition(condparams,matrix,Teuchos::null,residual,phinp,Teuchos::null,condstring);
  discret_->ClearState();

  return;
} // ScaTraTimIntImpl::EvaluateElectrodeKinetics

/*----------------------------------------------------------------------*
 | check for zero/negative concentration values               gjb 01/10 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::CheckConcentrationValues(RCP<Epetra_Vector> vec)
{
  // action only for ELCH applications
  if (IsElch(scatratype_))
  {
    // for NURBS discretizations we skip the following check.
    // Control points (i.e., the "nodes" and its associated dofs can be located
    // outside the domain of interest. Thus, they can have negative
    // concentration values although the concentration solution is positive
    // in the whole computational domain!
    if(dynamic_cast<DRT::NURBS::NurbsDiscretization*>(discret_.get())!=NULL)
      return;

    // this option can be helpful in some rare situations
    bool makepositive(false);

    std::vector<int> numfound(numscal_,0);
#if 0
    std::stringstream myerrormessage;
#endif
    for (int i = 0; i < discret_->NumMyRowNodes(); i++)
    {
      DRT::Node* lnode = discret_->lRowNode(i);
      std::vector<int> dofs;
      dofs = discret_->Dof(lnode);

      for (int k = 0; k < numscal_; k++)
      {
        const int lid = discret_->DofRowMap()->LID(dofs[k]);
        if (((*vec)[lid]) < EPS13 )
        {
          numfound[k]++;
          if (makepositive)
            ((*vec)[lid]) = EPS13;
#if 0
          myerrormessage<<"PROC "<<myrank_<<" dof index: "<<k<<setprecision(7)<<scientific<<
              " val: "<<((*vec)[lid])<<" node gid: "<<lnode->Id()<<
              " coord: [x] "<< lnode->X()[0]<<" [y] "<< lnode->X()[1]<<" [z] "<< lnode->X()[2]<<std::endl;
#endif
        }
      }
    }

    // print warning to screen
    for (int k = 0; k < numscal_; k++)
    {
      if (numfound[k] > 0)
      {
        std::cout<<"WARNING: PROC "<<myrank_<<" has "<<numfound[k]<<
        " nodes with zero/neg. concentration values for species "<<k;
        if (makepositive)
          std::cout<<"-> were made positive (set to 1.0e-13)"<<std::endl;
        else
          std::cout<<std::endl;
      }
    }

#if 0
    // print detailed info to error file
    for(int p=0; p < discret_->Comm().NumProc(); p++)
    {
      if (p==myrank_) // is it my turn?
      {
        // finish error message
        myerrormessage.flush();

        // write info to error file
        if ((errfile_!=NULL) and (myerrormessage.str()!=""))
        {
          fprintf(errfile_,myerrormessage.str().c_str());
          // std::cout<<myerrormessage.str()<<std::endl;
        }
      }
      // give time to finish writing to file before going to next proc ?
      discret_->Comm().Barrier();
    }
#endif

  }
  // so much code for a simple check!
  return;
} // ScaTraTimIntImpl::CheckConcentrationValues

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
    AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

  // add element parameters according to time-integration scheme
  AddSpecificTimeIntegrationParameters(eleparams);

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
    RCP<Epetra_CrsMatrix> crsPtent;
    MLAPI::GetPtent(*sysmat_sd_->EpetraMatrix(),mlparams,nullspace,crsPtent);
    LINALG::SparseMatrix Ptent(crsPtent);

    // compute scale-separation matrix: S = I - Ptent*Ptent^T
    Sep_ = LINALG::Multiply(Ptent,false,Ptent,true);
    Sep_->Scale(-1.0);
    RCP<Epetra_Vector> tmp = LINALG::CreateVector(Sep_->RowMap(),false);
    tmp->PutScalar(1.0);
    RCP<Epetra_Vector> diag = LINALG::CreateVector(Sep_->RowMap(),false);
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
    DynSmag_->AddScatra(discret_,scatratype_,pbcmapmastertoslave_);
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
    AddMultiVectorToParameterList(myparams,"convective velocity field",convel_);
    // add turbulence and mutifracatl subgrid-scales parameters
    myparams.sublist("TURBULENCE MODEL") = extraparams_->sublist("TURBULENCE MODEL");
    myparams.sublist("MULTIFRACTAL SUBGRID SCALES") = extraparams_->sublist("MULTIFRACTAL SUBGRID SCALES");
    // parameters for stabilization (required for evaluation of material)
    myparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");
    // add element parameters according to time-integration scheme
    AddSpecificTimeIntegrationParameters(myparams);

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

    extraparams_->sublist("MULTIFRACTAL SUBGRID SCALES").set<double>("meanCai",meanCai);
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
