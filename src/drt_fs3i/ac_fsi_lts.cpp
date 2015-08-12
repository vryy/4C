/*!----------------------------------------------------------------------
\file ac_fsi_lts.cpp

\brief cpp-file associated with algorithmic routines for two-way coupled partitioned
       solution approaches to fluid-structure-scalar-scalar interaction
       (FS3I). Specifically related version for multiscale approches. This file thereby holds
       all functions related with the large time scale simulation and
       the small to large to small time scale 'communication'.

\date 2015-07-29

\maintainer Moritz Thon
            thon@mhpc.mw.tum.de
            089/289-10364

\level 3
----------------------------------------------------------------------*/


#include "ac_fsi.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fsi/fsi_monolithic.H"
#include "../drt_adapter/ad_fld_fluid_ac_fsi.H"
#include "../drt_adapter/ad_str_fpsiwrapper.H"
#include "../drt_scatra/scatra_algorithm.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_inpar/inpar_material.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/growth_scd.H"
#include "../drt_mat/growth_law.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_mapextractor.H"

/*----------------------------------------------------------------------*
 | timeloop for small time scales                            Thon 07/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::LargeTimeScaleLoop()
{
  PrepareLargeTimeScaleLoop();

  while (LargeTimeScaleLoopNotFinished())
  {
    LargeTimeScalePrepareTimeStep();

    LargeTimeScaleOuterLoop();

    LargeTimeScaleUpdateAndOutput();
  }

  FinishLargeTimeScaleLoop();
}

/*----------------------------------------------------------------------*
 | Prepare the large time scale loop                         Thon 08/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::PrepareLargeTimeScaleLoop()
{
  // print info
  if (Comm().MyPID()==0)
  {
    std::cout<<"\n************************************************************************"
               "\n                         LARGE TIME SCALE LOOP"
               "\n************************************************************************"<<std::endl;
  }

  //Set large time scale time step in the struct scatra field
  scatravec_[1]->ScaTraField()->SetDt( dt_large_ );

  //Set zeros velocities since we assume that the large time scale can not see the deformation of the small time scale
  Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(fsi_->StructureField()->Velnp()->Map(),true));
  scatravec_[1]->ScaTraField()->SetVelocityField(zeros,
                                         Teuchos::null,
                                         zeros,
                                         Teuchos::null,
                                         Teuchos::null,
                                         fsi_->StructureField()->Discretization() );

  fsineedsupdate_=false;

  // Save the phinp vector at the beginning of the large time scale loop in
  // in order to estimate the so far induced growth
  *structurephinp_bltsl_ = *scatravec_[1]->ScaTraField()->Phinp();
  *structurephinp_blgu_ = *scatravec_[1]->ScaTraField()->Phinp();
}

/*----------------------------------------------------------------------*
 | Finish the large time scale loop                          Thon 08/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::FinishLargeTimeScaleLoop()
{
  //Set small time scale time step size
  scatravec_[1]->ScaTraField()->SetDt( dt_ );

  // Fix time and step in fsi and fluid scatra field

  // We start the small time scale with a new cycle. But since dt_large is a
  // multiple of fsiperiod_ we are already at the this time_
  // We do not modify the step_ counter; we just keep counting..

  SetTimeStepInFSI(time_,step_);
  // we now have to fix the time_ and step_ of the structure field, since this is not shifted
  // in PrepareTimeStep(), but in Update() which we here will not call. So..
  fsi_->StructureField()->SetTime(time_);
  fsi_->StructureField()->SetTimen(time_+fsi_->StructureField()->Dt());
  fsi_->StructureField()->SetStep(step_);
  fsi_->StructureField()->SetStepn(step_+1);

  scatravec_[0]->ScaTraField()->SetTimeStep(time_,step_);

  //clean up large time scale loop
  structurephinp_bltsl_->PutScalar(0.0);
  structurephinp_blgu_->PutScalar(0.0);

  //we start with a clean small time scale loop
  fsiisperiodic_ = false;
  scatraisperiodic_ = false;
}

/*----------------------------------------------------------------------*
 | timeloop for large time scales                            Thon 07/15 |
 *----------------------------------------------------------------------*/
bool FS3I::ACFSI::LargeTimeScaleLoopNotFinished()
{
  return NotFinished() and not fsineedsupdate_;
}

/*----------------------------------------------------------------------*
 | Prepare small time scale time step                        Thon 07/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::LargeTimeScalePrepareTimeStep()
{
  //Increment time and step
  step_ += 1;
  time_ += dt_large_;

  //Print to screen
  if (Comm().MyPID()==0)
  {
    std::cout << "\n\n"<< "TIME:  "    << std::scientific <<std::setprecision(8)<< time_ << "/" << std::scientific << timemax_
             << "     DT = " << std::scientific << dt_large_
             << "     STEP = " << std::setw(4) << step_ << "/" << std::setw(4) << numstep_
             << "\n\n";
  }

  //prepare structure scatra field
  scatravec_[1]->ScaTraField()->PrepareTimeStep();
}

/*----------------------------------------------------------------------*
 | OuterLoop for sequentially staggered FS3I scheme          Thon 08/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::LargeTimeScaleOuterLoop()
{
  DoStructScatraStep();

  if ( DoesGrowthNeedsUpdate() ) //includes the check for fsineedsupdate_
  {
    LargeTimeScaleDoGrowthUpdate();
  }

}

/*----------------------------------------------------------------------*
 | Do a large time scale structe scatra step                 Thon 08/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::DoStructScatraStep()
{
  if (Comm().MyPID()==0)
  {
    std::cout<<"\n************************************************************************"
               "\n                       AC STRUCTURE SCATRA SOLVER"
               "\n************************************************************************\n"<<std::endl;

    std::cout<<"+- step/max -+- tol ---- [norm] -+-- scal-res --+-- scal-inc --+"<<std::endl;
  }

  bool stopnonliniter=false;
  int itnum = 0;

  while (stopnonliniter==false)
  {
    StructScatraEvaluateSolveIterUpdate();
    itnum++;
    if (StructScatraConvergenceCheck(itnum))
      break;
  }
}

/*--------------------------------------------------------------------------------*
 | evaluate, solver and iteratively update structure scalar problem    Thon 08/15 |
 *--------------------------------------------------------------------------------*/
void FS3I::ACFSI::StructScatraEvaluateSolveIterUpdate()
{
  if (infperm_)
      dserror("This not a valid option!"); //just for safety

  const Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra = scatravec_[1]->ScaTraField(); //structure scatra

  //----------------------------------------------------------------------
  //evaluate the structure scatra field
  //----------------------------------------------------------------------
  scatra->PrepareLinearSolve();

  //----------------------------------------------------------------------
  // calculate contributions due to finite interface permeability
  //----------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> rhs_struct_scal = scatracoupforce_[1];
  Teuchos::RCP<LINALG::SparseMatrix> mat_struct_scal = scatracoupmat_[1];

  rhs_struct_scal->PutScalar(0.0);
  mat_struct_scal->Zero();

  scatra->SurfacePermeability(mat_struct_scal,rhs_struct_scal);

  // apply Dirichlet boundary conditions to coupling matrix and vector
  const Epetra_Map* dofrowmap = scatra->Discretization()->DofRowMap();
  Teuchos::RCP<Epetra_Vector> zeros = LINALG::CreateVector(*dofrowmap,true);

  const Teuchos::RCP< const Epetra_Map > dbcmap = scatra->DirichMaps()->CondMap();
  mat_struct_scal->ApplyDirichlet(*dbcmap,false);
  LINALG::ApplyDirichlettoSystem(rhs_struct_scal,zeros,*dbcmap);

  //----------------------------------------------------------------------
  // add coupling to the resiudal
  //----------------------------------------------------------------------
  const Teuchos::RCP<Epetra_Vector> residual = scatra->Residual();
  residual->Update(1.0,*rhs_struct_scal,1.0);

  // add contribution of the fluid field
  Teuchos::RCP<Epetra_Vector> rhs_fluid_scal_boundary = scatrafieldexvec_[0]->ExtractVector(scatracoupforce_[0],1);
  Teuchos::RCP<Epetra_Vector> rhs_fluid_scal = scatrafieldexvec_[1]->InsertVector(Scatra1ToScatra2(rhs_fluid_scal_boundary),1);
  //Note: we have to care for the fact that scatracoupforce_[0] has been evaluatet with the small time step:
  rhs_fluid_scal->Scale(dt_large_/dt_);

  residual->Update(-1.0,*rhs_fluid_scal,1.0);

  //----------------------------------------------------------------------
  // add coupling to the sysmat
  //----------------------------------------------------------------------
  const Teuchos::RCP<LINALG::SparseMatrix> sysmat = scatra->SystemMatrix();
  sysmat->Add(*mat_struct_scal,false,1.0,1.0);

  //----------------------------------------------------------------------
  // solve the scatra problem
  //----------------------------------------------------------------------
  const Teuchos::RCP<Epetra_Vector> structurescatraincrement = LINALG::CreateVector(*scatra->DofRowMap(),true);
  scatra->Solver()->Solve( sysmat->EpetraOperator(),
                           structurescatraincrement,
                           residual,
                           true,
                           true);

  //----------------------------------------------------------------------
  // update the strucutre scatra increment
  //----------------------------------------------------------------------
  scatra->UpdateIter(structurescatraincrement);
}

/*----------------------------------------------------------------------*
 | check convergence of structure scatra field               Thon 08/15 |
 *----------------------------------------------------------------------*/
bool FS3I::ACFSI::StructScatraConvergenceCheck(const int itnum)
{

  const Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra = scatravec_[1]->ScaTraField(); //structure scatra

  // some input parameters for the scatra fields
  const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();
  const int scatraitemax = scatradyn.sublist("NONLINEAR").get<int>("ITEMAX");
  const double scatraittol = scatradyn.sublist("NONLINEAR").get<double>("CONVTOL");
  const double scatraabstolres = scatradyn.sublist("NONLINEAR").get<double>("ABSTOLRES");


  double conresnorm(0.0);
  scatra->Residual()->Norm2(&conresnorm);
  double incconnorm(0.0);
  scatra->Increment()->Norm2(&incconnorm);
  double phinpnorm(0.0);
  scatra->Phinp()->Norm2(&phinpnorm);

  // care for the case that nothing really happens in the concentration field
  if (phinpnorm < 1e-5)
    phinpnorm = 1.0;

  // print the screen info
  if (Comm().MyPID()==0)
  {
    printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |\n",
           itnum,scatraitemax,scatraittol,conresnorm,incconnorm/phinpnorm);
  }

  // this is the convergence check
  // We always require at least one solve. We test the L_2-norm of the
  // current residual. Norm of residual is just printed for information
  if (conresnorm <= scatraittol and incconnorm/phinpnorm <= scatraittol)
  {
    if (Comm().MyPID()==0)
    {
      // print 'finish line'
      printf("+------------+-------------------+--------------+--------------+\n\n");
    }
    return true;
  }

  // abort iteration, when there's nothing more to do! -> more robustness
  else if (conresnorm < scatraabstolres)
  {
    // print 'finish line'
    if (Comm().MyPID()==0)
    {
      printf("+------------+-------------------+--------------+--------------+\n\n");
    }
    return true;
  }

  // if itemax is reached without convergence stop the simulation
  else if (itnum == scatraitemax)
  {
    if (Comm().MyPID()==0)
    {
      printf("+---------------------------------------------------------------+\n");
      printf("|    scalar-scalar field did not converge in itemax steps!     |\n");
      printf("+---------------------------------------------------------------+\n");
    }
    // yes, we stop!
    dserror("Structure scatra not converged in itemax steps!");
    return true;
  }
  else
    return false;
}

/*----------------------------------------------------------------------*
 | Do we need to update the structure scatra displacments               |
 | due to growth                                             Thon 08/15 |
 *----------------------------------------------------------------------*/
bool FS3I::ACFSI::DoesGrowthNeedsUpdate()
{
  bool growthneedsupdate = false;

  // check if the structure material is a growth material. We assume here
  // that the structure has the same material for the whole discretiazation.
  // Hence we check only the first element:
  Teuchos::RCP<DRT::Discretization> structuredis = fsi_->StructureField()->Discretization();
  const int GID = structuredis->ElementColMap()->GID(0); //global element ID

  Teuchos::RCP<MAT::Material> structurematerial = structuredis->gElement(GID)->Material();;

  if ( structurematerial->MaterialType() == INPAR::MAT::m_growth_volumetric_scd )
  {
    //get alpha and growth inducing scalar.
    double alpha = 0.0;
    int sc1 = 1;

    Teuchos::RCP<MAT::GrowthMandel> growthscdmaterial = Teuchos::rcp_dynamic_cast<MAT::GrowthMandel>(structurematerial);

    if (growthscdmaterial==Teuchos::null)
      dserror("Dynamic cast to MAT::GrowthMandel failed!");

    Teuchos::RCP<MAT::GrowthLaw> growthlaw = growthscdmaterial->Parameter()->growthlaw_;

    switch (growthlaw->MaterialType())
    {
      case INPAR::MAT::m_growth_ac:
      {
        Teuchos::RCP<MAT::GrowthLawAC> growthlawac = Teuchos::rcp_dynamic_cast<MAT::GrowthLawAC>(growthlaw);
        if (growthscdmaterial==Teuchos::null)
          dserror("Dynamic cast to MAT::GrowthLawAC failed!");
        alpha = growthlawac->Parameter()->alpha_;
        sc1 = growthlawac->Parameter()->Sc1_;
        break;
      }
      case INPAR::MAT::m_growth_ac_radial:
      {
        Teuchos::RCP<MAT::GrowthLawACRadial> growthlawacradial = Teuchos::rcp_dynamic_cast<MAT::GrowthLawACRadial>(growthlaw);
        if (growthlawacradial==Teuchos::null)
          dserror("Dynamic cast to MAT::GrowthLawACRadial failed!");
        alpha = growthlawacradial->Parameter()->alpha_;
        sc1 = growthlawacradial->Parameter()->Sc1_;
        break;
      }
      default:
      {
        dserror("Growth law not supported in AC-FS3I!");
        break;
      }
    }
    //Puh! That was exhausting. But we have to keep going.

    const Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra = scatravec_[1]->ScaTraField(); //structure scatra
    const Teuchos::RCP<const Epetra_Vector> phinp = scatra->Phinp(); //fluidscatra

    //first we check the error, iff we need a growth update
    {
      const Teuchos::RCP<Epetra_Vector> phidiff_blgu = LINALG::CreateVector(*scatra->DofRowMap(),true);

      //build difference vector with the reference
      phidiff_blgu->Update(1.0,*phinp,-1.0,*structurephinp_blgu_,0.0);

      //Extract the dof of interest
      Teuchos::RCP<Epetra_Vector> phidiff_blgu_j = extractjthscalar_[sc1-1]->ExtractCondVector(phidiff_blgu);

      //get the maximum
      double max_phidiff_blgu = 0.0;
      phidiff_blgu_j->MaxValue(&max_phidiff_blgu);

      if (Comm().MyPID()==0)
        std::cout<<std::scientific<<std::setprecision(2)<<" The maximal local growth since the last growth update is "<<alpha*max_phidiff_blgu<<std::endl;

      //now the actual comparison:
      const double growth_tol = DRT::Problem::Instance()->FS3IDynamicParams().sublist("AC").get<double>("GROWTH_TOL");

      if (max_phidiff_blgu * alpha >= growth_tol)
      {
        growthneedsupdate = true;
      }
    }

    //now check the error, iff we need a small time scale update
    {
      const Teuchos::RCP<Epetra_Vector> phidiff_bltsl_ = LINALG::CreateVector(*scatra->DofRowMap(),true);

      //build difference vector with the reference
      phidiff_bltsl_->Update(1.0,*phinp,-1.0,*structurephinp_bltsl_,0.0);

      //Extract the dof of interest
      Teuchos::RCP<Epetra_Vector> phidiff_bltsl_j = extractjthscalar_[sc1-1]->ExtractCondVector(phidiff_bltsl_);

      //get the maximum
      double max_phidiff_bltsl = 0.0;
      phidiff_bltsl_j->MaxValue(&max_phidiff_bltsl);

      if (Comm().MyPID()==0)
        std::cout<<std::scientific<<std::setprecision(2)<<" The maximal local growth since the beginning of the large time scale loop is "<<alpha*max_phidiff_bltsl<<std::endl;

      //now the actual comparison:
      const double fsi_update_tol = DRT::Problem::Instance()->FS3IDynamicParams().sublist("AC").get<double>("FSI_UPDATE_TOL");

      if (max_phidiff_bltsl * alpha >= fsi_update_tol)
      {
        growthneedsupdate = true; //we want to update the displacements before we go back to the small time scale
        fsineedsupdate_ = true;
      }
    }
  }
  else
    dserror("In AC-FS3I we want growth, so use a growth material like MAT_GrowthVolumetricScd!");

  return growthneedsupdate;
}

/*-------------------------------------------------------------------------*
 | update the structure scatra displacments due to growth       Thon 08/15 |
 *-------------------------------------------------------------------------*/
void FS3I::ACFSI::LargeTimeScaleDoGrowthUpdate()
{
  const Teuchos::RCP<SCATRA::ScaTraTimIntImpl> fluidscatra = scatravec_[0]->ScaTraField();
  const Teuchos::RCP<SCATRA::ScaTraTimIntImpl> structurescatra = scatravec_[1]->ScaTraField();

  // Note: we never do never proceed with time_ and step_, so this really just about updating the growth, i.e.
  // the displacments of the structure scatra fiels

  //----------------------------------------------------------------------
  // print to screen
  //----------------------------------------------------------------------
  if (Comm().MyPID()==0)
  {
    std::cout<<"\n************************************************************************"
               "\n                            AC GROWTH UPDATE"
               "\n************************************************************************"<<std::endl;
  }


  //----------------------------------------------------------------------
  // finish present structure scatra time step (no output)
  //----------------------------------------------------------------------
  structurescatra->Update(1);

  //----------------------------------------------------------------------
  // Switch time step of structure scatra
  //----------------------------------------------------------------------
  //Switch back the time step to do the update with the same (small) timestep as the fsi (subcycling time step possible!)
  const double dt_fluid = fsi_->FluidField()->Dt();
  fluidscatra->SetDt( dt_fluid );
  structurescatra->SetDt( dt_fluid );

  //----------------------------------------------------------------------
  // Prepare time steps
  //----------------------------------------------------------------------
  // fsi problem
  fsi_->PrepareTimeStep();
  //scatra fields
  fluidscatra->PrepareTimeStep();
  structurescatra->PrepareTimeStep();

  //----------------------------------------------------------------------
  // Fix time_ and step_ counters
  //----------------------------------------------------------------------
  // fsi problem
  SetTimeStepInFSI(time_,step_);
  //scatra fields
  fluidscatra->SetTimeStep(time_,step_);
  structurescatra->SetTimeStep(time_,step_);

  //----------------------------------------------------------------------
  // do the growth update
  //----------------------------------------------------------------------
  //NOTE: in there SetFSISolution() does not mess up our WSS since fsiisperiodic_ == true! But let's check for safety
  if ( not fsiisperiodic_)
    dserror("Here fsiisperiodic_ must be true to not mess up the WSS in ExtractWSS()");

  //the actuall calculations
  SmallTimeScaleOuterLoopIterStagg();

  //----------------------------------------------------------------------
  // write the output
  //----------------------------------------------------------------------
  //write fsi output. Scatra outputs are done later
  // fsi output
  fsi_->PrepareOutput();
  fsi_->Update();
  FsiOutput();
  //fluid scatra update. Structure scatra is done later
  fluidscatra->Update(0);

  //----------------------------------------------------------------------
  // Switch back time step and velocity field
  //----------------------------------------------------------------------
  //Now set the time step back:
  fluidscatra->SetDt( dt_ );
  structurescatra->SetDt( dt_large_ );

  //Set zeros velocities since we assume that the large time scale can not see the deformation of the small time scale
  Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(fsi_->StructureField()->Velnp()->Map(),true));
  structurescatra->SetVelocityField(zeros,
                                    Teuchos::null,
                                    zeros,
                                    Teuchos::null,
                                    Teuchos::null,
                                    fsi_->StructureField()->Discretization() );

  //----------------------------------------------------------------------
  // Update growth check vector
  //----------------------------------------------------------------------
  *structurephinp_blgu_ = *structurescatra->Phinp();
}

/*----------------------------------------------------------------------*
 | Update and output the large time scale                    Thon 08/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::LargeTimeScaleUpdateAndOutput()
{
  //keep fsi time and fluid scatra field up to date
  SetTimeStepInFSI(time_,step_);
  scatravec_[0]->ScaTraField()->SetTimeStep(time_,step_);

  //write fsi and fluid scatra output. Structure scatra output is done later anyways
  //FSI update and output
//  fsi_->PrepareOutput();
//  fsi_->Update();
//  FsiOutput();

  //fluid scatra update and output. structure scatra is done later
//  scatravec_[0]->ScaTraField()->Update(0);
  scatravec_[0]->ScaTraField()->Output(0);

  //now update and output the structure scatra field
  scatravec_[1]->ScaTraField()->Update(1);
  scatravec_[1]->ScaTraField()->Output(1);
}

/*----------------------------------------------------------------------*
 | Build map extractor which extracts the j-th dof           Thon 08/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<LINALG::MapExtractor> > FS3I::ACFSI::BuildMapExtractor()
{
  std::vector<Teuchos::RCP<LINALG::MapExtractor> > extractjthscalar;

  const Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra = scatravec_[1]->ScaTraField(); //structure scatra
  const int numscal = scatra->NumScal();
  const Teuchos::RCP<const DRT::Discretization> dis = scatra->Discretization();

  for (int k=0; k<numscal; k++)
  {
    std::set<int> conddofset;
    std::set<int> otherdofset;

    int numrownodes = dis->NumMyRowNodes();
    for (int i=0; i<numrownodes; ++i)
    {
      DRT::Node* node = dis->lRowNode(i);

      std::vector<int> dof = dis->Dof(0,node);
      if ( dof.size()!=(unsigned)scatravec_[1]->ScaTraField()->NumScal() )
        dserror("There was some error building the Map Extractor!");
      for (unsigned j=0; j<dof.size(); ++j)
      {
        // test for dof position
        if (j != static_cast<unsigned>(k))
        {
          otherdofset.insert(dof[j]);
        }
        else
        {
          conddofset.insert(dof[j]);
        }
      }
    }
    std::vector<int> conddofmapvec;
    conddofmapvec.reserve(conddofset.size());
    conddofmapvec.assign(conddofset.begin(), conddofset.end());
    conddofset.clear();
    Teuchos::RCP<Epetra_Map> conddofmap = Teuchos::rcp(new Epetra_Map(-1,
                                  conddofmapvec.size(),
                                  &conddofmapvec[0],
                                  0,
                                  dis->Comm()));
    conddofmapvec.clear();

    std::vector<int> otherdofmapvec;
    otherdofmapvec.reserve(otherdofset.size());
    otherdofmapvec.assign(otherdofset.begin(), otherdofset.end());
    otherdofset.clear();
    Teuchos::RCP<Epetra_Map> otherdofmap =
      Teuchos::rcp(new Epetra_Map(-1,
                                  otherdofmapvec.size(),
                                  &otherdofmapvec[0],
                                  0,
                                  dis->Comm()));
    otherdofmapvec.clear();

    Teuchos::RCP<LINALG::MapExtractor> getjdof= Teuchos::rcp(new LINALG::MapExtractor);
    getjdof->Setup(*dis->DofRowMap(),conddofmap,otherdofmap);
    extractjthscalar.push_back(getjdof);
  }

  return extractjthscalar;
}

/*----------------------------------------------------------------------*
 | Compare if two doubles are relatively equal               Thon 08/15 |
 *----------------------------------------------------------------------*/
bool FS3I::ACFSI::IsRealtiveEqualTo(const double A, const double B, const double Ref)
{
  return ( (std::abs(A-B)/Ref) < 1e-12 );
}

/*----------------------------------------------------------------------*
 | Compare if A mod B is relatively equal to zero            Thon 08/15 |
 *----------------------------------------------------------------------*/
bool FS3I::ACFSI::ModuloIsRealtiveZero(const double value, const double modulo, const double Ref)
{
  return IsRealtiveEqualTo(fmod(value+modulo/2,modulo)-modulo/2,0.0,Ref);
}
