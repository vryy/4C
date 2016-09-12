/*!----------------------------------------------------------------------
\file strtimint_statmech.cpp

\brief time integration for structural problems with statistical mechanics

\level 3

\maintainer Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276

*----------------------------------------------------------------------*/

#include <Teuchos_Time.hpp>
#include "Teuchos_RCP.hpp"

#include "strtimint_statmech.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_locsys.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_krylov_projector.H"
#include "../drt_statmech/statmech_manager.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_constraint/constraint_manager.H"
#include "../drt_constraint/constraintsolver.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_beamcontact.H"
#include "../drt_plastic_ssn/plastic_ssn_manager.H"
#include "../drt_structure/strtimint.H"
#include "../drt_beamcontact/beam3contact_manager.H"
#include "stru_aux.H"

#include "../drt_beam3/beam3.H"
#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/beam3eb.H"
#include "../drt_beam3/beam3cl.H"
#include "../drt_truss3/truss3.H"
#include "../drt_truss3cl/truss3cl.H"
#include "../drt_truss2/truss2.H"
#include "../drt_discsh3/discsh3.H"

#include <fenv.h>


//#define GMSHPTCSTEPS
//#define MEASURETIME

/*----------------------------------------------------------------------*
 |  ctor (public)                                          mueller 03/12|
 *----------------------------------------------------------------------*/
STR::TimIntStatMech::TimIntStatMech(const Teuchos::ParameterList& params,
                                    const Teuchos::ParameterList& sdynparams,
                                    const Teuchos::ParameterList& xparams,
                                    Teuchos::RCP<DRT::Discretization> actdis,
                                    Teuchos::RCP<LINALG::Solver> solver,
                                    Teuchos::RCP<LINALG::Solver> contactsolver,
                                    Teuchos::RCP<IO::DiscretizationWriter> output) :
TimIntOneStepTheta(sdynparams,params,sdynparams, xparams, actdis,solver,contactsolver,output),
isconverged_(false),
iterges_(0)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g. redistribution of elements.
  // Only then call the setup to this class. This will call the setup to all classes in the inheritance hierarchy.
  // This way, this class may also override a method that is called during Setup() in a base class.
  return;
} // STR::TimIntStatMech::TimIntStatMech()

/*----------------------------------------------------------------------------------------------*
 * Initialize this class                                                            rauch 09/16 |
 *----------------------------------------------------------------------------------------------*/
void STR::TimIntStatMech::Init
(
    const Teuchos::ParameterList& timeparams,
    const Teuchos::ParameterList& sdynparams,
    const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<LINALG::Solver> solver
)
{
  // call Init() in base class
  STR::TimIntOneStepTheta::Init(timeparams,sdynparams,xparams,actdis,solver);

  //getting number of dimensions for diffusion coefficient calculation
  ndim_= DRT::Problem::Instance()->NDim();

  // print dbc type for this simulation to screen
  if(HaveStatMech())
    StatMechPrintBCType();

  if(!discret_->Comm().MyPID())
  {
    if(HaveStatMech())
    {
      std::cout<<"StatMech output path: "<<statmechman_->StatMechRootPath()<<"/StatMechOutput/"<<std::endl;
      std::cout<<"================================================================"<<std::endl;
    }
    else if (HaveStatMechBilayer())
    {
      std::cout<<"StatMech output path: "<<statmechmanBilayer_->StatMechRootPath()<<"/StatMechOutput/"<<std::endl;
      std::cout<<"================================================================"<<std::endl;
    }
  }

  //suppress all output printed to screen in case of single filament studies in order not to generate too much output on the cluster
  if(HaveStatMech())
    SuppressOutput();

  //in case that beam contact is activated by respective input parameter, a Beam3cmanager object is created
  if(HaveStatMech())
    InitializeBeamContact();

  //temporary safety check -> will be removed when time loop is reorganized
  if(DRT::Problem::Instance()->ProblemType() != prb_statmech)
    dserror("Due to current reorganizing of the time loop: Set your PROBLEMTYP to StatMech");

  return;
}

/*----------------------------------------------------------------------------------------------*
 * Setup this class                                                                 rauch 09/16 |
 *----------------------------------------------------------------------------------------------*/
void STR::TimIntStatMech::Setup()
{
  // call Setup() in base class
  STR::TimIntOneStepTheta::Setup();

  // retrieve number of random numbers per element and store them in randomnumbersperelement_
  RandomNumbersPerElement();

  // retrieve number of random numbers per node and store them in randomnumberspernode_
  RandomNumbersPerNode();

  return;
}

/*----------------------------------------------------------------------*
 |  print Boundary condition type to screen      (public) mueller 03/12 |
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::StatMechPrintBCType()
{
  if(!discret_->Comm().MyPID())
  {
    INPAR::STATMECH::DBCType dbctype = DRT::INPUT::IntegralValue<INPAR::STATMECH::DBCType>(statmechman_->GetStatMechParams(),"DBCTYPE");
    switch(dbctype)
    {
      // standard DBC application
      case INPAR::STATMECH::dbctype_std:
        std::cout<<"- Conventional Input file based application of DBCs"<<std::endl;
      break;
      // shear with a fixed Dirichlet node set
      case INPAR::STATMECH::dbctype_shearfixed:
        std::cout<<"- DBCs for rheological measurements applied: fixed node set"<<std::endl;
      break;
      case INPAR::STATMECH::dbctype_shearfixeddel:
        std::cout<<"- DBCs for rheological measurements applied: fixed node set"<<std::endl;
      break;
      // shear with an updated Dirichlet node set (only DOF in direction of oscillation is subject to BC, others free)
      case INPAR::STATMECH::dbctype_sheartrans:
        std::cout<<"- DBCs for rheological measurements applied: transient node set"<<std::endl;
      break;
      // pin down and release individual nodes
      case INPAR::STATMECH::dbctype_pinnodes:
        std::cout<<"- Special DBCs pinning down selected nodes"<<std::endl;
      break;
      // apply affine shear deformation
      case INPAR::STATMECH::dbctype_affineshear:
        std::cout<<"- DBCs for affine shear deformation"<<std::endl;
      break;
      case INPAR::STATMECH::dbctype_affinesheardel:
        std::cout<<"- DBCs for affine shear deformation"<<std::endl;
      break;
      case INPAR::STATMECH::dbctype_movablesupport1d:
        std::cout<<"- DBCs 1D movable support"<<std::endl;
      break;
      // no DBCs at all
      case INPAR::STATMECH::dbctype_none:
        std::cout<<"- No application of DBCs (i.e. also no DBCs by Input file)"<<std::endl;
      break;
      // default: everything involving periodic boundary conditions
      default:
        dserror("Check your DBC type! %i", dbctype);
      break;
    }
    INPAR::STATMECH::NBCType nbctype = DRT::INPUT::IntegralValue<INPAR::STATMECH::NBCType>(statmechman_->GetStatMechParams(), "NBCTYPE");
    switch(nbctype)
    {
      case INPAR::STATMECH::nbctype_constcreep:
        std::cout<<"- constant shear Neumann boundary condition on a nodal subset" << std::endl;
        break;
      case INPAR::STATMECH::nbctype_randompointforce:
        std::cout<<"- random time-dependent Neumann point force on single nodes" << std::endl;
        break;
      default:
        std::cout<<"- standard input file based definition of Neumann boundary conditions / no NBCs"<< std::endl;
        break;
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  random number per element                    (public) mueller 03/12 |
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::RandomNumbersPerElement()
{
  //maximal number of random numbers to be generated per time step for any column map element of this processor
  int randomnumbersperlocalelement = 0;

  /*check maximal number of nodes of an element with stochastic forces on this processor*/
  for (int i=0; i<  discret_->NumMyColElements(); ++i)
  {
    const DRT::ElementType & eot = discret_->lColElement(i)->ElementType();
    /*stochastic forces implemented so far only for the following elements:*/
    if ( eot == DRT::ELEMENTS::Beam3Type::Instance() )
    {
      //see whether current element needs more random numbers per time step than any other before
      randomnumbersperlocalelement = std::max(randomnumbersperlocalelement,dynamic_cast<DRT::ELEMENTS::Beam3*>(discret_->lColElement(i))->HowManyRandomNumbersINeed());

      //in case of periodic boundary conditions beam3 elements require a special initialization if they are broken by the periodic boundaries in the initial configuration
      if((statmechman_->GetPeriodLength())->at(0) > 0.0)
        statmechman_->PeriodicBoundaryBeam3Init(discret_->lColElement(i));
    }
    else if ( eot == DRT::ELEMENTS::Beam3rType::Instance() )
    {
      //see whether current element needs more random numbers per time step than any other before
      randomnumbersperlocalelement = std::max(randomnumbersperlocalelement,dynamic_cast<DRT::ELEMENTS::Beam3r*>(discret_->lColElement(i))->HowManyRandomNumbersINeed());

      //in case of periodic boundary conditions beam3 elements require a special initialization if they are broken by the periodic boundaries in the initial configuration
      if((statmechman_->GetPeriodLength())->at(0) > 0.0)
        statmechman_->PeriodicBoundaryBeam3rInit(discret_->lColElement(i));
    }
    else if ( eot == DRT::ELEMENTS::BeamCLType::Instance() )
    {
      //see whether current element needs more random numbers per time step than any other before
      randomnumbersperlocalelement = std::max(randomnumbersperlocalelement,dynamic_cast<DRT::ELEMENTS::BeamCL*>(discret_->lColElement(i))->HowManyRandomNumbersINeed());

      //in case of periodic boundary conditions beam3 elements require a special initialization if they are broken by the periodic boundaries in the initial configuration
      if((statmechman_->GetPeriodLength())->at(0) > 0.0)
        statmechman_->PeriodicBoundaryBeamCLInit(discret_->lColElement(i));
    }
    else if ( eot == DRT::ELEMENTS::Beam3ebType::Instance() )
    {
//      //see whether current element needs more random numbers per time step than any other before
      randomnumbersperlocalelement = std::max(randomnumbersperlocalelement,dynamic_cast<DRT::ELEMENTS::Beam3eb*>(discret_->lColElement(i))->HowManyRandomNumbersINeed());
//
//      //in case of periodic boundary conditions beam3 elements require a special initialization if they are broken by the periodic boundaries in the initial configuration
      if((statmechman_->GetPeriodLength())->at(0) > 0.0)
        statmechman_->PeriodicBoundaryBeam3ebInit(discret_->lColElement(i));
    }
    else if ( eot == DRT::ELEMENTS::Truss3Type::Instance() )
    {
      //see whether current element needs more random numbers per time step than any other before
      randomnumbersperlocalelement = std::max(randomnumbersperlocalelement,dynamic_cast<DRT::ELEMENTS::Truss3*>(discret_->lColElement(i))->HowManyRandomNumbersINeed());

      //in case of periodic boundary conditions truss3 elements require a special initialization if they are broken by the periodic boundaries in the initial configuration
      if((statmechman_->GetPeriodLength())->at(0) > 0.0)
        statmechman_->PeriodicBoundaryTruss3Init(discret_->lColElement(i));
    }
    else if ( eot == DRT::ELEMENTS::Truss3CLType::Instance() )
    {
      //see whether current element needs more random numbers per time step than any other before
      randomnumbersperlocalelement = std::max(randomnumbersperlocalelement,dynamic_cast<DRT::ELEMENTS::Truss3CL*>(discret_->lColElement(i))->HowManyRandomNumbersINeed());

      //in case of periodic boundary conditions Truss3 elements require a special initialization if they are broken by the periodic boundaries in the initial configuration
      if((statmechman_->GetPeriodLength())->at(0) > 0.0)
        statmechman_->PeriodicBoundaryTruss3CLInit(discret_->lColElement(i));
    }
    else if ( eot == DRT::ELEMENTS::DiscSh3Type::Instance() )
    {
      //see whether current element needs more random numbers per time step than any other before
      randomnumbersperlocalelement=std::max(randomnumbersperlocalelement,dynamic_cast<DRT::ELEMENTS::DiscSh3*>(discret_->lColElement(i))->HowManyRandomNumbersINeed());
    }
    else
      continue;
  } //for (int i=0; i<dis_.NumMyColElements(); ++i)

  /*so far the maximal number of random numbers required per element has been checked only locally on this processor;
   *now we compare the results of each processor and store the maximal one in maxrandomnumbersperglobalelement_*/
  discret_->Comm().MaxAll(&randomnumbersperlocalelement,&maxrandomnumbersperglobalelement_ ,1);

  return;
}

/*----------------------------------------------------------------------*
 |  random number per node                    (public) mueller 03/12 |
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::RandomNumbersPerNode()
{
  //maximal number of random numbers to be generated per time step for any column map element of this processor
  int randomnumbersperlocalnode = 0;

  /*check maximal number of nodes of an element with stochastic forces on this processor*/
  for (int i=0; i<  discret_->NumMyColElements(); ++i)
  {
    const DRT::ElementType & eot = discret_->lColElement(i)->ElementType();
    /*stochastic forces implemented so far only for the following elements:*/
     if ( eot == DRT::ELEMENTS::DiscSh3Type::Instance() )
    {
      //see whether current element needs more random numbers per time step than any other before
      randomnumbersperlocalnode=std::max(randomnumbersperlocalnode,dynamic_cast<DRT::ELEMENTS::DiscSh3*>(discret_->lColElement(i))->HowManyRandomNumbersPerNode());
    }
    else
      continue;
  } //for (int i=0; i<dis_.NumMyColElements(); ++i)

  /*so far the maximal number of random numbers required per element has been checked only locally on this processor;
   *now we compare the results of each processor and store the maximal one in maxrandomnumbersperglobalelement_*/
  discret_->Comm().MaxAll(&randomnumbersperlocalnode,&maxrandomnumbersperglobalnode_ ,1);

  return;
}

/*----------------------------------------------------------------------*
 |  suppress output in some cases                (public) mueller 03/12 |
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::SuppressOutput()
{
  if( DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechman_->GetStatMechParams(), "SPECIAL_OUTPUT") == INPAR::STATMECH::statout_endtoendlog ||
      DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechman_->GetStatMechParams(), "SPECIAL_OUTPUT") == INPAR::STATMECH::statout_endtoendconst ||
      DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechman_->GetStatMechParams(), "SPECIAL_OUTPUT") == INPAR::STATMECH::statout_orientationcorrelation ||
      DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechman_->GetStatMechParams(), "SPECIAL_OUTPUT") == INPAR::STATMECH::statout_anisotropic)
  {
    printscreen_ = 0;
    std::cout<<"\n\nPay Attention: from now on regular output to screen suppressed !!!\n\n";
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Initialize beam contact                      (public) mueller 03/12 |
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::InitializeBeamContact()
{
  if(HaveBeamContact())
  {
    //check wheter appropriate parameters are set in the parameter list "CONTACT DYNAMIC"
    //initialize beam contact detection strategy (for network simulations, octree is the choice)
    const Teuchos::ParameterList& beamcontact = DRT::Problem::Instance()->BeamContactParams();

    // conditions for potential-based beam interaction
    std::vector<DRT::Condition*> beampotconditions(0);
    discret_->GetCondition("BeamPotentialLineCharge",beampotconditions);

    if (DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Strategy>(beamcontact,"BEAMS_STRATEGY") != INPAR::BEAMCONTACT::bstr_none or (int)beampotconditions.size()!=0)
      buildoctree_ = true;
    else
      dserror("Check your input parameters in CONTACT DYNAMIC!");
  }
  return;
}

///*----------------------------------------------------------------------*
// |  integrate in time          (static/public)             mueller 03/12|
// *----------------------------------------------------------------------*/
//int STR::TimIntStatMech::Integrate()
//{
//  dserror("time loop of statmech problems moved to ADAPTER::StructureStatMech::Integrate()"
//      "this method will be removed soon");
//  // set statmech internal time and time step size
//  statmechman_->UpdateTimeAndStepSize((*dt_)[0],(*time_)[0],true);
//  // this is necessary in case the time step size changed with the initial step (originally timen_ is set in strtimint.cpp)
//  timen_ = (*time_)[0] + (*dt_)[0];
//
//  double eps = 1.0e-12;
//  while( (timen_ <= timemax_+eps) and (stepn_ <= stepmax_) )
//  {
//#ifdef MEASURETIME
//    const double t0 = Teuchos::Time::wallTime();
//#endif
//
//    // preparations for statistical mechanics in this time step
//    StatMechPrepareStep();
//#ifdef MEASURETIME
//    const double t1 = Teuchos::Time::wallTime();
//#endif
//    //redo time step in case of bad random configuration
//    do
//    {
//      // Update of statmech specific quantities as well as new set of random numbers
//      StatMechUpdate();
//
//      // Solve system of equations according to parameters and methods chosen by input file
//      if(HaveBeamContact())
//        BeamContactNonlinearSolve();
//      else // standard procedure without beam contact
//      {
//        //pay attention: for a constant predictor an incremental velocity update is necessary, which has
//        //been deleted out of the code in order to simplify it
////        if(!discret_->Comm().MyPID())
////          std::cout<<"target time = "<<timen_<<", time step = "<<(*dt_)[0]<<std::endl;
//        Predict();
//
//        Solve();
//      }
//
//      /*if iterations have not converged a new trial requires setting all intern element variables, statmechmanager class variables
//       *and the state of the discretization to status at the beginning of this time step*/
//      StatMechRestoreConvState();
//    }
//    while(!isconverged_);
//
//#ifdef MEASURETIME
//    const double t2 = Teuchos::Time::wallTime();
//#endif
//
//    //periodic shift of configuration at the end of the time step in order to avoid improper output
//    statmechman_->PeriodicBoundaryShift(*disn_, ndim_, timen_, (*dt_)[0]);
//
//#ifdef MEASURETIME
//    const double t3 = Teuchos::Time::wallTime();
//#endif
//
//    // update all that is relevant
//    UpdateAndOutput();
//
//#ifdef MEASURETIME
//    const double t4 = Teuchos::Time::wallTime();
//#endif
//
//    //special output for statistical mechanics
//    StatMechOutput();
//
//#ifdef MEASURETIME
//    if(!discret_->Comm().MyPID())
//    {
//      std::cout<<"\n=================Time  Measurement================"<<std::endl;
//      std::cout<<"TimIntStatMech::Integrate"<<std::endl;
//      std::cout<<"StatMechPrepareStep :\t"<<std::setprecision(4)<<t1-t0<<"\ts"<<std::endl;
//      std::cout<<"Newton              :\t"<<std::setprecision(4)<<t2-t1<<"\ts"<<std::endl;
//      std::cout<<"PeriodicShift       :\t"<<std::setprecision(4)<<t3-t2<<"\ts"<<std::endl;
//      std::cout<<"UpdateAndOutput     :\t"<<std::setprecision(4)<<t4-t3<<"\ts"<<std::endl;
//      std::cout<<"StatMechOutput      :\t"<<std::setprecision(4)<<Teuchos::Time::wallTime()-t4<<"\ts"<<std::endl;
//      std::cout<<"=================================================="<<std::endl;
//      std::cout<<"total time         :\t"<<std::setprecision(4)<<Teuchos::Time::wallTime()-t0<<"\ts"<<std::endl;
//    }
//#endif
//  }
//  return 0;
//} // void STR::TimIntStatMech::Integrate()


/*----------------------------------------------------------------------*
 |does what it says it does                       (public) mueller 03/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::UpdateAndOutput()
{
  // calculate stresses, strains and energies
  // note: this has to be done before the update since otherwise a potential
  // material history is overwritten
  PrepareOutput();

  // update displacements, velocities, accelerations (UpdateStepBeamContact() happens within)
  // after this call we will have disn_==dis_, etc
  UpdateStepState();

  // hand over time step size and time of the latest converged time step
  if(HaveStatMech())
    statmechman_->UpdateTimeAndStepSize((*dt_)[0], timen_);
  else if(HaveStatMechBilayer())
  {
    statmechmanBilayer_->UpdateTimeAndStepSize((*dt_)[0], timen_);
  }

  // update time and step
//  if(!discret_->Comm().MyPID())
//    std::cout<<"pre UpdateStepTime           : time = "<<(*time_)[0]<<", timen_ = "<<timen_<<", dt = "<<(*dt_)[0]<<std::endl;
  UpdateStepTime();

  // update time and step
//  if(!discret_->Comm().MyPID())
//    std::cout<<"pre UpdateStepElement        : time = "<<(*time_)[0]<<", timen_ = "<<timen_<<", dt = "<<(*dt_)[0]<<"\n\n"<<std::endl;
  // update everything on the element level
  UpdateStepElement();

  // write output
  Output();

  // print info about finished time step
  PrintStep();

  return;
}//UpdateAndOutput()

/*----------------------------------------------------------------------*
 |do consistent predictor step for Brownian dynamics (public)cyron 10/09|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::Predict()
{
  //Set iter_ to zero for the predictor step
  iter_=0;
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  //consistent predictor for backward Euler time integration scheme (theta==1.0)
  if(pred_==INPAR::STR::pred_constdis)
  {
    PredictConstDisConsistVel();
    // note: currently just c&p, have to look into it...
    normdisi_ = 1.0e6;
    normpres_ = 1.0e6;
  }
  else
    dserror("Trouble in determining predictor %i", pred_);

  // Apply Dirichlet Boundary Conditions
  ApplyDirichletBC(timen_, disn_, veln_);

  // calculate internal force and stiffness matrix
  Teuchos::ParameterList params;
  params.set<bool>("predict",true);
  EvaluateForceStiffResidual(params);

  // reactions are negative to balance residual on DBC
  freact_->Update(-1.0, *fres_, 0.0);
// not necessary, since fres_ DBC DOFs were already blanked
//  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);

  //------------------------------------------------ build residual norm
  CalcRefNorms();

  PrintPredictor();

  return;
} //STR::TimIntStatMech::Predict()

/*----------------------------------------------------------------------*
 |  predictor                                     (public) mueller 03/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::PredictConstDisConsistVel()
{
  // time step size
  const double dt = (*dt_)[0];

  /*special part for STATMECH: initialize disn_ and veln_ with zero; this is necessary only for the following case:
   * assume that an iteration step did not converge and is repeated with new random numbers; if the failure of conver
   * gence lead to disn_ = NaN and veln_ = NaN this would affect also the next trial as e.g. disn_->Update(1.0,*((*dis_)(0)),0.0);
   * would set disn_ to NaN as even 0*NaN = NaN!; this would defeat the purpose of the repeated iterations with new
   * random numbers and has thus to be avoided; therefore we initialized disn_ and veln_ with zero which has no effect
   * in any other case*/
  disn_->PutScalar(0.0);
  veln_->PutScalar(0.0);

  // constant predictor : displacement in domain
  disn_->Update(1.0, *(*dis_)(0), 0.0);

  // new end-point velocities
  veln_->Update(1.0/(theta_*dt), *disn_,
                -1.0/(theta_*dt), *(*dis_)(0),
                0.0);
  veln_->Update(-(1.0-theta_)/theta_, *(*vel_)(0),
                1.0);

  // no accelerations (1st order differential equation)
//  accn_->Update(1.0/(theta_*theta_*dt*dt), *disn_,
//                -1.0/(theta_*theta_*dt*dt), *(*dis_)(0),
//                0.0);
//  accn_->Update(-1.0/(theta_*theta_*dt), *(*vel_)(0),
//                -(1.0-theta_)/theta_, *(*acc_)(0),
//                1.0);
  return;
}//STR::TimIntStatMech::PredictConstDisConsistVel()

/*----------------------------------------------------------------------*
 | solve equilibrium                              (public) mueller 03/12|
 *----------------------------------------------------------------------*/
INPAR::STR::ConvergenceStatus STR::TimIntStatMech::Solve()
{
  if(HaveBeamContact())
  {
    BeamContactNonlinearSolve();
  }
  else
  {
    if(itertype_==INPAR::STR::soltech_ptc)
      PTC();
    else if(itertype_==INPAR::STR::soltech_newtonfull)
      NewtonFull();
    else if(itertype_==INPAR::STR::soltech_newtonls)
      NewtonLS();
    else
      dserror("itertype %d not implemented for StatMech applications! Choose ptc/fullnewton/lsnewton !", itertype_);
  }

  INPAR::STR::ConvergenceStatus status = INPAR::STR::conv_success;

  if(isconverged_)
    status = INPAR::STR::conv_success;
  else
    status = INPAR::STR::conv_nonlin_fail;

  // Only relevant, if the input parameter DIVERCONT is used and set to divcontype_ == adapt_step:
  // In this case, the time step size is halved as consequence of a non-converging nonlinear solver.
  // After a prescribed number of converged time steps, the time step is doubled again. The following
  // methods checks, if the time step size can be increased again.
  CheckForTimeStepIncrease(status);

  return status;
}

/*-----------------------------------------------------------------------------------*
 | error action in case of non-convergence of nonlinear solver  (public)   meier 01/15|
 *-----------------------------------------------------------------------------------*/
bool STR::TimIntStatMech::PerformErrorAction()
{
  // what to do when nonlinear solver does not converge

  /*if iterations have not converged a new trial requires setting all intern element variables, statmechmanager class variables
   *and the state of the discretization to status at the beginning of this time step*/
  StatMechRestoreConvState();

  switch (divcontype_)
  {
    case INPAR::STR::divcont_stop:
    {
      // write restart output of last converged step before stopping
      OutputStep(true);

      // we should not get here, dserror for safety
      dserror("Nonlinear solver did not converge! ");
      return false;
    }
    case INPAR::STR::divcont_continue:
    {
      dserror("DIVERCONT-Type divcont_continue not supported by statmech so far! ");
      return false;
    }
    break;
    case INPAR::STR::divcont_repeat_step:
    {
      IO::cout << "Nonlinear solver failed to converge. Repeat time step with new random numbers!"
               << IO::endl;
      // do nothing since we didn't update yet
      return true;
    }
    break;
    case INPAR::STR::divcont_halve_step:
    {
      IO::cout << "Nonlinear solver failed to converge. Divide timestep in half!"
               << IO::endl;

      // undo last time step increment
      timen_-=(*dt_)[0];
      // halve the time step size
      (*dt_)[0]=(*dt_)[0]*0.5;
      // update the number of max time steps
      stepmax_= stepmax_ + (stepmax_-stepn_)+1;
      // reset timen_ according to new time step increment
      timen_+= (*dt_)[0];

      return false;
    }
    break;
    case INPAR::STR::divcont_adapt_step:
    {
      IO::cout << "Nonlinear solver failed to converge. Divide timestep in half!"
               << IO::endl;

      // undo last time step increment
      timen_-=(*dt_)[0];
      // halve the time step size
      (*dt_)[0]=(*dt_)[0]*0.5;
      // update the number of max time steps
      stepmax_= stepmax_ + (stepmax_-stepn_)+1;
      // reset timen_ according to new time step increment
      timen_+= (*dt_)[0];

      divconrefinementlevel_++;
      divconnumfinestep_=0;

      return false;
    }
    break;
    case INPAR::STR::divcont_rand_adapt_step:
    {
      // generate random number between 0.51 and 1.99 (as mean value of random
      // numbers generated on all processors), alternating values larger
      // and smaller than 1.0
      double proc_randnum_get = ((double) rand()/(double) RAND_MAX);
      double proc_randnum = proc_randnum_get;
      double randnum = 1.0;
      discret_->Comm().SumAll(&proc_randnum,&randnum,1);
      randnum /= discret_->Comm().NumProc();
      if (rand_tsfac_ > 1.0)      rand_tsfac_ = randnum*0.49+0.51;
      else if (rand_tsfac_ < 1.0) rand_tsfac_ = randnum*0.99+1.0;
      else                        rand_tsfac_ = randnum*1.48+0.51;
      if (myrank_ == 0) IO::cout << "Nonlinear solver failed to converge: modifying time-step size by random number between 0.51 and 1.99 -> here: " << rand_tsfac_ << " !" << IO::endl;
      // multiply time-step size by random number
      (*dt_)[0]=(*dt_)[0]*rand_tsfac_;
      // update maximum number of time steps
      stepmax_= (1.0/rand_tsfac_)*stepmax_ + (1.0-(1.0/rand_tsfac_))*stepn_ + 1;
      // reset timen_ because it is set in the constructor
      timen_ = (*time_)[0] + (*dt_)[0];;
      return INPAR::STR::conv_success;
    }
    break;
    case INPAR::STR::divcont_repeat_simulation:
    {
      dserror("DIVERCONT-Type divcont_repeat_simulation not supported by statmech so far! ");
      return false;
    }
    break;
    default:
      dserror("Unknown DIVER_CONT case");
    return false;
    break;
  }

  return false; // make compiler happy
}

/*----------------------------------------------------------------------*
 | apply Dirichlet Boundary Conditions            (public) mueller 03/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::ApplyDirichletBC(const double                time,
                                           Teuchos::RCP<Epetra_Vector> dis,
                                           Teuchos::RCP<Epetra_Vector> vel)
{
  // needed parameters
  Teuchos::ParameterList p;
  p.set("total time", time);  // target time (i.e. timen_)
  p.set("delta time", (*dt_)[0]);

  // set vector values needed by elements
  discret_->ClearState();

  //discret_->SetState("displacement",disn_);
  //discret_->SetState("velocity",veln_);
  // predicted dirichlet values
  // disn then also holds prescribed new dirichlet displacements

  // determine DBC evaluation mode (new vs. old)
  if (HaveStatMech())
    statmechman_->EvaluateDirichletStatMech(p, dis, vel, dbcmaps_);
  if (HaveStatMechBilayer())
    statmechmanBilayer_->EvaluateDirichletStatMech(p, dis, vel, dbcmaps_);

  discret_->ClearState();

  return;
} //STR::TimIntStatMech::ApplyDirichletBC()

/*----------------------------------------------------------------------*
 |  evaluate residual                             (public) mueller 03/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::EvaluateForceStiffResidual(Teuchos::ParameterList& params)
{
  // get info about prediction step from parameter list
  bool predict = false;
  if(params.isParameter("predict"))
    predict = params.get<bool>("predict");

  // initialise stiffness matrix to zero
  stiff_->Zero();

  // theta-interpolate state vectors
  EvaluateMidState();

  // build new external forces
  fextn_->PutScalar(0.0);
  ApplyForceExternal(timen_, (*dis_)(0), disn_, (*vel_)(0), fextn_, stiff_);

  // additional external forces are added (e.g. interface forces)
  fextn_->Update(1.0, *fifc_, 1.0);

  // initialise internal forces
  fintn_->PutScalar(0.0);
  // internal forces and stiffness matrix
  ApplyForceStiffInternal(timen_, (*dt_)[0], disn_, disi_, veln_, fintn_, stiff_);

  // note: neglected in statmech...
  // inertial forces #finertt_
  //mass_->Multiply(false, *acct_, *finertt_);

  // viscous forces due Rayleigh damping
  if (damping_ == INPAR::STR::damp_rayleigh)
    damp_->Multiply(false, *velt_, *fvisct_);

  // Build residual
  BuildResidual();

  // evaluate beam contact
  if(HaveBeamContact())
    ApplyForceStiffBeamContact(stiff_, fres_, disn_, predict);

  // blank residual at DOFs on Dirichlet BC already here (compare with strtimint_impl: there, this is first called on freact_, then again on fres_ which seems unnecessary)
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);

  // build tangent matrix : effective dynamic stiffness matrix
  //    K_{Teffdyn} = 1/(theta*dt^2) M
  //                + 1/dt C
  //                + theta K_{T}
  // note : no mass terms in statmech case -> zeros in matrix -> Comment it(?)
  //stiff_->Add(*mass_, false, 1.0/(theta_*(*dt_)[0]*(*dt_)[0]), theta_);
  if (damping_ == INPAR::STR::damp_rayleigh)
    stiff_->Add(*damp_, false, 1.0/(*dt_)[0], 1.0);

  // close stiffness matrix
  stiff_->Complete();

//  // blank residual at DOFs on Dirichlet BC
//  {
//    Epetra_Vector frescopy(*fres_);
//    fres_->Multiply(1.0,*invtoggle_,frescopy,0.0);
//  }

  return;
}//STR::TimIntStatMech::EvaluateForceStiffResidual()

/*----------------------------------------------------------------------*
 | evaluate theta-state vectors by averaging end-point vectors          |
 |                                                (public) mueller 03/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::EvaluateMidState()
{
  // mid-displacements D_{n+1-alpha_f} (dism)
  //    D_{n+theta} := theta * D_{n+1} + (1-theta) * D_{n}
  dist_->Update(theta_, *disn_, 1.0-theta_, *(*dis_)(0), 0.0);

  // mid-velocities V_{n+1-alpha_f} (velm)
  //    V_{n+theta} := theta * V_{n+1} + (1-theta) * V_{n}
  velt_->Update(theta_, *veln_, 1.0-theta_, *(*vel_)(0), 0.0);

  // note: no accelerations in statmech...
  // mid-accelerations A_{n+1-alpha_m} (accm)
  //    A_{n+theta} := theta * A_{n+1} + (1-theta) * A_{n}
  // acct_->Update(theta_, *accn_, 1.0-theta_, *(*acc_)(0), 0.0);

  return;
}


/*----------------------------------------------------------------------*
 | evaluate external forces at t_{n+1} (public)            mueller 09/13|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::ApplyForceExternal(const double                          time,
                                             const Teuchos::RCP<Epetra_Vector>     dis,
                                             const Teuchos::RCP<Epetra_Vector>     disn,
                                             const Teuchos::RCP<Epetra_Vector>     vel,
                                             Teuchos::RCP<Epetra_Vector>&          fext,
                                             Teuchos::RCP<LINALG::SparseOperator>& fextlin)
{
  Teuchos::ParameterList p;
  // other parameters needed by the elements
  p.set("total time", time);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState(0,"displacement", dis);

  if (damping_ == INPAR::STR::damp_material)
    discret_->SetState(0,"velocity", vel);

  if (HaveStatMech())
    statmechman_->EvaluateNeumannStatMech(p,disn, fext, fextlin);
  else if (HaveStatMechBilayer())
    statmechmanBilayer_->EvaluateNeumannStatMech(p,disn, fext, fextlin);

  // go away
  return;
}

/*----------------------------------------------------------------------*
 | internal forces and stiffness (public)                  mueller 03/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::ApplyForceStiffInternal(const double                         time,
                                                  const double                         dt,
                                                  const Teuchos::RCP<Epetra_Vector>    dis,  // displacement state
                                                  const Teuchos::RCP<Epetra_Vector>    disi,  // residual displacements
                                                  const Teuchos::RCP<Epetra_Vector>    vel,  // velocity state
                                                  Teuchos::RCP<Epetra_Vector>          fint,  // internal force
                                                  Teuchos::RCP<LINALG::SparseOperator> stiff, // stiffness matrix
                                                  double                               t_eval) // time
{
  double t_evaluate = Teuchos::Time::wallTime();

  // create the parameters for the discretization
  Teuchos::ParameterList p;

  //passing statistical mechanics parameters to elements
  if(HaveStatMech())
    statmechman_->AddStatMechParamsTo(p, randomnumbers_);
  else if(HaveStatMechBilayer())
  {
    statmechmanBilayer_->AddStatMechParamsTo(p, randomnumbers_);

    // get reference volume
    p.set("action", "calc_struct_refvol");
    Teuchos::RCP<Epetra_SerialDenseVector> vol_ref
    = Teuchos::rcp(new Epetra_SerialDenseVector(1));
    discret_->EvaluateScalars(p, vol_ref);

    // get reference CG
    p.set("action", "calc_struct_refCG");
    Teuchos::RCP<Epetra_SerialDenseVector> CG_ref
    = Teuchos::rcp(new Epetra_SerialDenseVector(3));
    discret_->EvaluateScalars(p, CG_ref);

    // get reference area
    p.set("action", "calc_struct_refarea");
    Teuchos::RCP<Epetra_SerialDenseVector> area_ref
    = Teuchos::rcp(new Epetra_SerialDenseVector(1));
    discret_->EvaluateScalars(p, area_ref);

    // get current volume
    discret_->SetState("displacement", disn_);
    p.set("action", "calc_struct_currvol");
    Teuchos::RCP<Epetra_SerialDenseVector> vol_curr
    = Teuchos::rcp(new Epetra_SerialDenseVector(1));
    discret_->EvaluateScalars(p, vol_curr);

    // get current CG
    p.set("action", "calc_struct_currCG");
    Teuchos::RCP<Epetra_SerialDenseVector> CG_curr
    = Teuchos::rcp(new Epetra_SerialDenseVector(3));
    discret_->EvaluateScalars(p, CG_curr);


    // get current area
    p.set("action", "calc_struct_currarea");
    Teuchos::RCP<Epetra_SerialDenseVector> area_curr
    = Teuchos::rcp(new Epetra_SerialDenseVector(1));
    discret_->EvaluateScalars(p, area_curr);

    discret_->ClearState();

    p.set("reference volume",double((*vol_ref)(0)));
    p.set("current volume",double((*vol_curr)(0)));
    p.set<Teuchos::RCP<Epetra_SerialDenseVector> >("reference CG", CG_ref);
    p.set<Teuchos::RCP<Epetra_SerialDenseVector> >("current CG", CG_curr);
    p.set("reference area",double((*area_ref)(0)));
    p.set("current area",double((*area_curr)(0)));
  }

  // action for elements
  const std::string action = "calc_struct_nlnstiff";
  p.set("action", action);
  p.set("total time", time);
  p.set("delta time", dt);

  // set vector values needed by elements
  discret_->ClearState();
  // extended SetState(0,...) in case of multiple dofsets (e.g. TSI)
  discret_->SetState(0,"residual displacement", disi);
  discret_->SetState(0,"displacement", dis);
  discret_->SetState(0,"velocity", vel);
  discret_->Evaluate(p, stiff, Teuchos::null, fint, Teuchos::null, Teuchos::null);
  discret_->ClearState();
  if (HaveFaceDiscret())
  {
    AssembleEdgeBasedMatandRHS(p,fint,dis,vel);
  }

  t_eval += timer_->WallTime() - t_evaluate;

  // that's it
  return;
}

/*----------------------------------------------------------------------*
 |  calculate reference norms for relative convergence checks           |
 |                                                (public) mueller 03/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::BuildResidual()
{
  // build residual  Res = M . A_{n+theta}
  //                     + C . V_{n+theta}
  //                     + F_{int;n+theta}
  //                     - F_{ext;n+theta}
  fres_->Update(-theta_, *fextn_, -(1.0-theta_), *fext_, 0.0);
  fres_->Update(theta_, *fintn_, (1.0-theta_), *fint_, 1.0);

  if (damping_ == INPAR::STR::damp_rayleigh)
    fres_->Update(1.0, *fvisct_, 1.0);
  // note: finertt_ is zero vector in statmech case
  // fres_->Update(1.0, *finertt_, 1.0);
  return;
}

/*----------------------------------------------------------------------*
 |  calculate reference norms for relative convergence checks           |
 |                                                (public) mueller 03/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::CalcRefNorms()
{
  normfres_ = STR::AUX::CalculateVectorNorm(iternorm_, fres_);
  // determine characteristic norms
  // we set the minumum of CalcRefNormForce() and #tolfres_, because
  // we want to prevent the case of a zero characteristic fnorm
  normcharforce_ = CalcRefNormForce();
  if (normcharforce_ == 0.0) normcharforce_ = tolfres_;
  normchardis_ = CalcRefNormDisplacement();
  if (normchardis_ == 0.0) normchardis_ = toldisi_;

  return;
}

/*----------------------------------------------------------------------*
 |  Newton-Raphson iteration                      (public) mueller 04/13|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::NewtonFull()
{
  // we do a Newton-Raphson iteration here.
  // the specific time integration has set the following
  // --> On #fres_ is the positive force residuum
  // --> On #stiff_ is the effective dynamic stiffness matrix

  // check whether we have a sanely filled stiffness matrix
  if (not stiff_->Filled()){ dserror("Effective stiffness matrix must be filled here");}

  if (outputeveryiter_)
  {
    int restart = DRT::Problem::Instance()->Restart();
    if (stepn_ == (restart + 1)) outputcounter_ = 0;
    OutputEveryIter(true);
  }

  // initialise equilibrium loop
  iter_ = 1;
  normfres_ = CalcRefNormForce();
  // normdisi_ was already set in predictor; this is strictly >0
  timer_->ResetStartTime();

  // create out-of-balance force for 2nd, 3rd, ... Uzawa iteration
  if(HaveBeamContact())
    InitializeNewtonUzawa();

  while ( ( (not Converged()) and (iter_ <= itermax_) ) or (iter_ <= itermin_) )
  {
    // increment total number of iterations
    iterges_++;
    // make negative residual
    fres_->Scale(-1.0);

    // apply Dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_,
                                   GetLocSysTrafo(), zeros_, *(dbcmaps_->CondMap()));

    // *********** time measurement ***********
    double dtcpu = timer_->WallTime();
    // *********** time measurement ***********

//    // STC preconditioning (no effect on StatMech)
//    STCPreconditioning();

    // solve for disi_
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (solveradapttol_ and (iter_ > 1))
    {
      double worst = normfres_;
      double wanted = tolfres_;
      solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
    }

    // linear solver call
    solver_->Solve(stiff_->EpetraOperator(), disi_, fres_, true, iter_==1, projector_);
    solver_->ResetTolerance();

    //In beam contact applications it can be necessary to limit the Newton step size (scaled residual displacements)
    if(HaveStatMech())
      LimitStepsizeBeamContact(disi_);

    //Biopolymer network applications (for filaments with beam3eb elements and spring crosslinkers)
    if (HaveStatMech())
      LimitStepsizeBeam(disi_);

//    // recover standard displacements (no effect on StatMech)
//    RecoverSTCSolution();

    // *********** time measurement ***********
    dtsolve_ = timer_->WallTime() - dtcpu;
    // *********** time measurement ***********

    // update end-point displacements etc
    UpdateIter(iter_);

    // compute residual forces #fres_ and stiffness #stiff_
    // whose components are globally oriented
    Teuchos::ParameterList params;
    EvaluateForceStiffResidual(params);

    // extract reaction forces
    // reactions are negative to balance residual on DBC
    freact_->Update(-1.0, *fres_, 0.0);
    dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);

    // blank residual at DOFs on Dirichlet BC
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);

//    // cancel in residual those forces that would excite rigid body modes and
//    // that thus vanish in the Krylov space projection (no effect on StatMech)
//    if (projector_!=Teuchos::null)
//      projector_->ApplyPT(*fres_);

    // (trivial)
    if (pressure_ != Teuchos::null)
    {
      Teuchos::RCP<Epetra_Vector> pres = pressure_->ExtractCondVector(fres_);
      Teuchos::RCP<Epetra_Vector> disp = pressure_->ExtractOtherVector(fres_);
      normpfres_ = STR::AUX::CalculateVectorNorm(iternorm_, pres);
      normfres_ = STR::AUX::CalculateVectorNorm(iternorm_, disp);

      pres = pressure_->ExtractCondVector(disi_);
      disp = pressure_->ExtractOtherVector(disi_);
      normpres_ = STR::AUX::CalculateVectorNorm(iternorm_, pres);
      normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disp);
    }
    else
    {
      // build residual force norm
      normfres_ = STR::AUX::CalculateVectorNorm(iternorm_, fres_);
      // build residual displacement norm
      normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);
    }

    // print stuff
    if(printscreen_)
      PrintNewtonIter();

    // increment equilibrium loop index
    iter_ += 1;

    // DEBUG: Gmsh Output
//    GmshOutputEveryIter();
    //    if(HaveStatMech())
//      if()
    // Print GMSHOUTPUT
    //********************Begin: GMSH Output*******************************
        // STEP 1: OUTPUT OF TIME STEP INDEX
//        std::ostringstream filename;
//        filename << "/home/mukherjee/workspace/baci/release/GmshOutput/";
//
//        if (iter_<1000000)
//          filename << "network_timestep" << std::setw(4) << std::setfill('0')<<step_+1<<"Iter"<<std::setw(2) << std::setfill('0') << iter_-1;
//        else /*(timestep>=1000000)*/
//          dserror("ERROR: Gmsh output implemented for max 999.999 time steps");
//
//        // finish filename
//        filename<<".pos";
//
////        statmechman_->GmshOutput((*disn_),filename,step_,(*time_)[0]);
//        statmechmanBilayer_->GmshOutput((*disn_),filename,step_,(*time_)[0]);
   //********************End: GMSH Output:*******************************

        // Print total surface area
//        PrintTotalSurfAreaperIter();

        // Print total enclosed Volume
//        PrintTotalVolperIter();

        // Print total EdgeLength
//        PrintTotalEdgeLenperIter();

    // leave the loop without going to maxiter iteration because most probably, the process will not converge anyway from here on
    if(normfres_>1.0e4 && iter_>4)
    {
      std::cout << "Attention: Left nonlinear equilibrium loop since normfres_>1.0e4!" << std::endl;
      break;
    }
  }  // end equilibrium loop

  // correct iteration counter
  iter_ -= 1;

  // call monitor (no effect on StatMech)
  if (conman_->HaveMonitor())
    conman_->ComputeMonitorValues(disn_);

  ConvergenceStatusUpdate(Converged());

  // test whether max iterations was hit
  if ( isconverged_ && !discret_->Comm().MyPID() && printscreen_)
  {
    std::cout<<"Newton-Raphson-iteration converged with..."<<std::endl;
    PrintNewtonIter();
  }

 // std::cout << "iterges_: " << iterges_ << std::endl;

  return;

} // STR::TimIntStatMech::FullNewton()

/*----------------------------------------------------------------------*
 |  initialize Newton for 2nd, 3rd, ... Uzawa iteration      cyron 12/10|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::InitializeNewtonUzawa()
{
//  bool  loadlin    = params_.get<bool>("LOADLIN",false);

  // create out-of-balance force for 2nd, 3rd, ... Uzawa iteration
  if (beamcman_->GetUzawaIter() > 1)
  {
    //--------------------------- recompute external forces if nonlinear
    // at state n, the external forces and linearization are interpolated at
    // time 1-alphaf in a TR fashion
//    if (loadlin)
//    {
//      Teuchos::ParameterList p;
//      // action for elements
//      p.set("action","calc_struct_eleload");
//      // other parameters needed by the elements
//      p.set("total time",timen_);
//      p.set("delta time",(*dt_)[0]);
//      p.set("alpha f",1-theta_);
//      // set vector values needed by elements
//      discret_->ClearState();
//      discret_->SetState("displacement",disn_);
//      discret_->SetState("velocity",veln_);
////      discret_->SetState("displacement",dism_); // mid point
////      discret_->SetState("velocity",velm_);
//      fextn_->PutScalar(0.0); // TR
////      fextm_->PutScalar(0.0);
//      fextlin_->Zero();
////      discret_->EvaluateNeumann(p,fextm_,fextlin_);
//      discret_->EvaluateNeumann(p,fextn_,fextlin_);
//      fextlin_->Complete();
//      discret_->ClearState();
//      fextm_->Update(1.0,*fextn_,0.0,*fext_,0.0);
//    }

    //---------------------------- compute internal forces and stiffness
    {
      // zero out stiffness
      stiff_->Zero();
      // create the parameters for the discretization
      Teuchos::ParameterList p;
      // action for elements
      p.set("action","calc_struct_nlnstiff");
      // other parameters that might be needed by the elements
      p.set("total time",timen_);
      p.set("delta time",(*dt_)[0]);
      p.set("alpha f",1-theta_);



      if(HaveStatMech())
        statmechman_->AddStatMechParamsTo(p, randomnumbers_);
      else if(HaveStatMechBilayer())
        statmechmanBilayer_->AddStatMechParamsTo(p, randomnumbers_);
      else
        dserror("No Statistical Mechanics module found!");

      // set vector values needed by elements
      discret_->ClearState();

      // scale IncD_{n+1} by (1-alphaf) to obtain mid residual displacements IncD_{n+1-alphaf}
//      disi_->Scale(1.-alphaf);

      discret_->SetState("residual displacement",disi_);
      discret_->SetState("displacement",disn_);
      discret_->SetState("velocity",veln_);

      fint_->PutScalar(0.0);  // initialise internal force vector
      discret_->Evaluate(p,stiff_,Teuchos::null,fint_,Teuchos::null,Teuchos::null);

      discret_->ClearState();
    }

    //------------------------------------------ compute residual forces
    fres_->Update(-1.0,*fint_,1.0,*fextn_,0.0); // fext_ oder fextn_ ???
    //**********************************************************************
    //**********************************************************************
    // evaluate beam contact
    if(HaveBeamContact())
    {
      Teuchos::ParameterList beamcontactparams;
      beamcontactparams.set("iter", iter_);
      beamcontactparams.set("dt", (*dt_)[0]);
      beamcontactparams.set("numstep", step_);

      beamcman_->Evaluate(*SystemMatrix(),*fres_,*disn_,beamcontactparams);
    }
    //**********************************************************************
    //**********************************************************************

    // blank residual DOFs that are on Dirichlet BC
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);
  }

  return;
}//STR::TimIntStatMech::InitializeNewtonUzawa()

/*----------------------------------------------------------------------*
 |  read restart (public)                                    cyron 12/08|
 *----------------------------------------------------------------------*/
//void STR::TimIntStatMech::ReadRestart(int step)
//{
//  Teuchos::RCP<DRT::Discretization> Teuchos::rcpdiscret = Teuchos::rcp(&discret_,false);
//  IO::DiscretizationReader reader(rcpdiscret,step);
//  double time  = reader.ReadDouble("time");
//  int    rstep = reader.ReadInt("step");
//  if (rstep != step) dserror("Time step on file not equal to given step");
//
//  reader.ReadVector(dis_, "displacement");
//  reader.ReadVector(vel_, "velocity");
//  reader.ReadVector(acc_, "acceleration");
//  reader.ReadVector(fext_,"fexternal");
//  reader.ReadMesh(step);
//
//  // read restart information for contact
//  statmechmanager_->ReadRestart(reader);
//
//#ifdef INVERSEDESIGNUSE
//  int idrestart = -1;
//  idrestart = reader.ReadInt("InverseDesignRestartFlag");
//  if (idrestart==-1) dserror("expected inverse design restart flag not on file");
//  // if idrestart==0 then the file is from a INVERSEDESIGCREATE phase
//  // and we have to zero out the inverse design displacements.
//  // The stored reference configuration is on record at the element level
//  if (!idrestart)
//  {
//    dis_->PutScalar(0.0);
//    vel_->PutScalar(0.0);
//    acc_->PutScalar(0.0);
//  }
//#endif
//
//  // override current time and step with values from file
//  params_.set<double>("total time",time);
//  params_.set<int>   ("step",rstep);
//
//  if (surf_stress_man_->HaveSurfStress())
//    surf_stress_man_->ReadRestart(rstep, DRT::Problem::Instance()->InputControlFile()->FileName());
//
//  if (constrMan_->HaveConstraint())
//  {
//    double uzawatemp = reader.ReadDouble("uzawaparameter");
//    constrSolv_->SetUzawaParameter(uzawatemp);
//    Teuchos::RCP<Epetra_Map> constrmap=constrMan_->GetConstraintMap();
//    Teuchos::RCP<Epetra_Vector> tempvec = LINALG::CreateVector(*constrmap,true);
//    reader.ReadVector(tempvec, "lagrmultiplier");
//    constrMan_->SetLagrMultVector(tempvec);
//    reader.ReadVector(tempvec, "refconval");
//    constrMan_->SetRefBaseValues(tempvec,time);
//  }
//
//  return;
//}//STR::TimIntStatMech::ReadRestart()

/*----------------------------------------------------------------------*
 |  do output including statistical mechanics data(public)    cyron 12/08|
 *----------------------------------------------------------------------*/
//void STR::TimIntStatMech::Output()
//{
//  // -------------------------------------------------------------------
//  // get some parameters from parameter list
//  // -------------------------------------------------------------------
//  double timen         = params_.get<double>("total time"             ,0.0);
//  double dt            = params_.get<double>("delta time"             ,0.01);
//  double alphaf        = 1.0-theta_;
//  int    istep         = params_.get<int>   ("step"                   ,0);
//  int    nstep         = params_.get<int>   ("nstep"                  ,5);
//  int    numiter       = params_.get<int>   ("num iterations"         ,-1);
//
//  bool   iodisp        = params_.get<bool>  ("io structural disp"     ,true);
//  int    updevrydisp   = params_.get<int>   ("io disp every nstep"    ,10);
//  INPAR::STR::StressType iostress = DRT::INPUT::get<INPAR::STR::StressType>(params_, "io structural stress",INPAR::STR::stress_none);
//  int    updevrystress = params_.get<int>   ("io stress every nstep"  ,10);
//  INPAR::STR::StrainType iostrain      = DRT::INPUT::get<INPAR::STR::StrainType>(params_, "io structural strain",INPAR::STR::strain_none);
//  bool   iosurfactant  = params_.get<bool>  ("io surfactant"          ,false);
//
//  int    writeresevry  = params_.get<int>   ("write restart every"    ,0);
//
//  bool   printscreen   = params_.get<bool>  ("print to screen"        ,true);
//  bool   printerr      = params_.get<bool>  ("print to err"           ,true);
//  FILE*  errfile       = params_.get<FILE*> ("err file"               ,NULL);
//  if (!errfile) printerr = false;
//
//  bool isdatawritten = false;
//
//  //------------------------------------------------- write restart step
//  if ((writeresevry and istep%writeresevry==0) or istep==nstep)
//  {
//    output_.WriteMesh(istep,timen);
//    output_.NewStep(istep, timen);
//    output_.WriteVector("displacement",dis_);
//    output_.WriteVector("velocity",vel_);
//    output_.WriteVector("acceleration",acc_);
//    output_.WriteVector("fexternal",fext_);
//
//#ifdef INVERSEDESIGNCREATE // indicate that this restart is from INVERSEDESIGCREATE phase
//    output_.WriteInt("InverseDesignRestartFlag",0);
//#endif
//#ifdef INVERSEDESIGNUSE // indicate that this restart is from INVERSEDESIGNUSE phase
//    output_.WriteInt("InverseDesignRestartFlag",1);
//#endif
//
//    isdatawritten = true;
//
////____________________________________________________________________________________________________________
////note:the following block is the only difference to Output() in strugenalpha.cpp-----------------------------
///* write restart information for statistical mechanics problems; all the information is saved as class variables
// * of StatMechManager*/
//    statmechman_->WriteRestart(output_);
////------------------------------------------------------------------------------------------------------------
////____________________________________________________________________________________________________________
//
//    if (surf_stress_man_->HaveSurfStress())
//      surf_stress_man_->WriteRestart(istep, timen);
//
//    if (constrMan_->HaveConstraint())
//    {
//      output_.WriteDouble("uzawaparameter",constrSolv_->GetUzawaParameter());
//      output_.WriteVector("lagrmultiplier",constrMan_->GetLagrMultVector());
//      output_.WriteVector("refconval",constrMan_->GetRefBaseValues());
//    }
//
//    if (discret_->Comm().MyPID()==0 and printscreen)
//    {
//      std::cout << "====== Restart written in step " << istep << std::endl;
//      fflush(stdout);
//    }
//    if (errfile and printerr)
//    {
//      fprintf(errfile,"====== Restart written in step %d\n",istep);
//      fflush(errfile);
//    }
//  }
//
//  //----------------------------------------------------- output results
//  if (iodisp and updevrydisp and istep%updevrydisp==0 and !isdatawritten)
//  {
//    output_.NewStep(istep, timen);
//    output_.WriteVector("displacement",dis_);
//    output_.WriteVector("velocity",vel_);
//    output_.WriteVector("acceleration",acc_);
//    output_.WriteVector("fexternal",fext_);
//    output_.WriteElementData();
//
//    if (surf_stress_man_->HaveSurfStress() and iosurfactant)
//      surf_stress_man_->WriteResults(istep,timen);
//
//    isdatawritten = true;
//  }
//
//  //------------------------------------- do stress calculation and output
//  if (updevrystress and !(istep%updevrystress) and iostress!=INPAR::STR::stress_none)
//  {
//    // create the parameters for the discretization
//    Teuchos::ParameterList p;
//    // action for elements
//    p.set("action","calc_struct_stress");
//    // other parameters that might be needed by the elements
//    p.set("total time",timen);
//    p.set("delta time",dt);
//    p.set("alpha f",1-theta_);
//    Teuchos::RCP<std::vector<char> > stress = Teuchos::rcp(new std::vector<char>());
//    Teuchos::RCP<std::vector<char> > strain = Teuchos::rcp(new std::vector<char>());
//    p.set("stress", stress);
//    p.set<int>("iostress", iostress);
//    p.set("strain", strain);
//    p.set<int>("iostrain", iostrain);
//    // set vector values needed by elements
//    discret_->ClearState();
//    discret_->SetState("residual displacement",zeros_);
//    discret_->SetState("displacement",dis_);
//    discret_->SetState("velocity",vel_);
//    discret_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
//    discret_->ClearState();
//    if (!isdatawritten) output_.NewStep(istep, timen);
//    isdatawritten = true;
//
//    switch (iostress)
//    {
//    case INPAR::STR::stress_cauchy:
//      output_.WriteVector("gauss_cauchy_stresses_xyz",*stress,*discret_->ElementRowMap());
//      break;
//    case INPAR::STR::stress_2pk:
//      output_.WriteVector("gauss_2PK_stresses_xyz",*stress,*discret_->ElementRowMap());
//      break;
//    case INPAR::STR::stress_none:
//      break;
//    default:
//      dserror ("requested stress type not supported");
//    }
//
//    switch (iostrain)
//    {
//    case INPAR::STR::strain_ea:
//      output_.WriteVector("gauss_EA_strains_xyz",*strain,*discret_->ElementRowMap());
//      break;
//    case INPAR::STR::strain_gl:
//      output_.WriteVector("gauss_GL_strains_xyz",*strain,*discret_->ElementRowMap());
//      break;
//    case INPAR::STR::strain_none:
//      break;
//    default:
//      dserror("requested strain type not supported");
//    }
//  }
//
//  //---------------------------------------------------------- print out
//  if (!myrank_)
//  {
//    if (printscreen)
//    {
//      printf("step %6d | nstep %6d | time %-14.8E | dt %-14.8E | numiter %3d\n",
//             istep,nstep,timen,dt,numiter);
//      printf("----------------------------------------------------------------------------------\n");
//      fflush(stdout);
//    }
//    if (printerr)
//    {
//      fprintf(errfile,"step %6d | nstep %6d | time %-14.8E | dt %-14.8E | numiter %3d\n",
//              istep,nstep,timen,dt,numiter);
//      fprintf(errfile,"----------------------------------------------------------------------------------\n");
//      fflush(errfile);
//    }
//  }
//  return;
//}//STR::TimIntStatMech::Output()

/*----------------------------------------------------------------------*
 |  Pseudo Transient Continuation                 (public)   cyron 12/10|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::PTC()
{
  //---------------------------------------------------------------some sanity checks
  if (not stiff_->Filled())
    dserror("Effective stiffness matrix must be filled here");
  //------------------------------------------------------------ for time measurement
  double sumsolver     = 0;
  double sumevaluation = 0;
  double sumptc = 0;

  const double tbegin = Teuchos::Time::wallTime();

  //--------------------create out-of-balance force for 2nd, 3rd, ... Uzawa iteration
  if(HaveBeamContact())
    InitializeNewtonUzawa();

  //=================================================== equilibrium loop
  // initialise equilibrium loop
  iter_ = 1;
  normfres_ = CalcRefNormForce();
  // normdisi_ was already set in predictor; this is strictly >0
  timer_->ResetStartTime();

  //-----------------------------parameters from statistical mechanics parameter list
  // hard wired ptc parameters
  Teuchos::ParameterList statmechparams;
  if(HaveStatMech())
    statmechparams = statmechman_->GetStatMechParams();
  if(HaveStatMechBilayer())
    statmechparams = statmechmanBilayer_->GetStatMechBilayerParams();
  double ctransptc = statmechparams.get<double>("CTRANSPTC0",0.0);
  // crotptc is used here as equivalent to dti of the PTC scheme
  double crotptc   = statmechparams.get<double>("CROTPTC0",0.145);
  double alphaptc  = statmechparams.get<double>("ALPHAPTC",6.0);

  // PTC parameters
  double nc;
  fres_->NormInf(&nc);
  double resinit = nc;

  if(nc==0.0)
    dserror("nc == 0.0! PTC scheme not applicable! Choose fullnewton");

  //printf("fresnorm %10.5e disinorm %10.5e nc %10.5e\n",normfres_,normdisi_,nc);

  //parameters to make sure that in last iteration step botch PTC parameters have reached zero
  double ctransptcold = ctransptc;
  double crotptcold   = crotptc;

  while (((!Converged() || ctransptcold > 0.0 || crotptcold > 0.0) and iter_<=itermax_) or (iter_ <= itermin_))
  {
    // increment total number of iterations
    iterges_++;

    //save PTC parameters of the so far last iteration step
    ctransptcold = ctransptc;
    crotptcold   = crotptc;

    // make negative residual
    fres_->Scale(-1.0);

    //backward Euler
    stiff_->Complete();

    //the following part was especially introduced for Brownian dynamics. Set to zero for std Newton-Raphson
    if(statmechparams.get<int>("MAXITERPTC",5)>0)
      PTCBrownianForcesAndDamping((*dt_)[0],crotptc,ctransptc, sumptc);

    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more

    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fres_,zeros_,*(dbcmaps_->CondMap()));

    //--------------------------------------------------- solve for disi
    const double t_solver = Teuchos::Time::wallTime();
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (solveradapttol_ and (iter_ > 1))
    {
      double worst = normfres_;
      double wanted = tolfres_;
      solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
    }
    solver_->Solve(stiff_->EpetraOperator(),disi_,fres_,true,iter_==1);
    solver_->ResetTolerance();

    //In beam contact applications it can be necessary to limit the Newton step size (scaled residual displacements)
    LimitStepsizeBeamContact(disi_);

    sumsolver += Teuchos::Time::wallTime() - t_solver;

    // update displacements and velocities for this iteration step
    UpdateIter(iter_);

    // empty parameter list
    Teuchos::ParameterList params;

    //---------------- compute internal forces, stiffness and residual
    EvaluateForceStiffResidual(params);

    // reactions are negative to balance residual on DBC
    // note: due to the use of the old "dirichtoggle_" vector, fres_ dofs with DBCs have already been blanked
    freact_->Update(-1.0, *fres_, 0.0);

    //---------------------------------------------- build residual norm
    // build residual displacement norm
    normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);
    // build residual force norm
    normfres_ = STR::AUX::CalculateVectorNorm(iternorm_, fres_);

    // update dti_ of the PTC scheme
    dti_ = crotptc;
    PrintNewtonIter();

    //------------------------------------ PTC update of artificial time
    // compute inf norm of residual
    PTCStatMechUpdate(ctransptc,crotptc,nc,resinit,alphaptc);
#ifdef GMSHPTCSTEPS
    // GmshOutput
    std::ostringstream filename;
    if(DRT::INPUT::IntegralValue<int>(statmechman_->GetStatMechParams(),"GMSHOUTPUT") && HaveBeamContact())
    {
      filename << statmechman_->StatMechRootPath() <<"/GmshOutput/network"<< (*time_)[0] <<"_u"<<std::setw(2) << std::setfill('0')<<beamcman_->GetUzawaIter()<<"_n"<<std::setw(2) << std::setfill('0')<<iter_<<".pos";
      statmechman_->GmshOutput(*disn_,filename,stepn_,timen_,beamcman_);
    }
    else
    {
      filename << statmechman_->StatMechRootPath()<< "/GmshOutput/network"<< (*time_)[0] <<"_n"<<std::setw(2) << std::setfill('0')<<iter_<<".pos";
      statmechman_->GmshOutput(*disn_,filename,stepn_,timen_);
    }
#endif
    //--------------------------------- increment equilibrium loop index
    ++iter_;

    // leave the loop without going to maxiter iteration because most probably, the process will not converge anyway from here on
    if(normfres_>1.0e4 && iter_>4)
    {
      std::cout << "Attention: Left nonlinear equilibrium loop since normfres_>1.0e4!" << std::endl;
      break;
    }
  }
  //============================================= end equilibrium loop
  // iter_ started at 1, so "--"
  iter_--;
//  print_unconv = false;

  ConvergenceStatusUpdate(Converged());

  // test whether max iterations was hit
  if ( isconverged_ && !discret_->Comm().MyPID() && printscreen_)
  {
    std::cout<<"PTC-iteration converged with..."<<std::endl;
    PrintNewtonIter();

    std::cout << "\n***\nevaluation time: " << sumevaluation<< " seconds\nptc time: "<< sumptc <<" seconds\nsolver time: "<< sumsolver <<" seconds\ntotal solution time: "<<Teuchos::Time::wallTime() - tbegin<<" seconds\n***\n";
  }

  return;
} // STR::TimIntStatMech::PTC()

/*----------------------------------------------------------------------*
 |  evaluate outcome of PTC and chose action accordingly   mueller 02/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::PTCBrownianForcesAndDamping(double& dt, double& crotptc, double& ctransptc, double& sumptc)
{
  const double t_ptc = Teuchos::Time::wallTime();
  // create the parameters for the discretization
  Teuchos::ParameterList p;

  p.set("action","calc_struct_ptcstiff");
  p.set("delta time",dt);
  p.set("crotptc",crotptc);
  p.set("ctransptc",ctransptc);

  //add statistical vector to parameter list for statistical forces and damping matrix computation
  if(HaveStatMech())
    statmechman_->AddStatMechParamsTo(p);
  else if(HaveStatMechBilayer())
    statmechmanBilayer_->AddStatMechParamsTo(p);

  //evaluate ptc stiffness contribution in all the elements
  discret_->Evaluate(p,stiff_,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  sumptc += Teuchos::Time::wallTime() - t_ptc;

  return;
}

/*----------------------------------------------------------------------*
 |  update of rot. and transl. ptc damping for statmech    mueller 03/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::PTCStatMechUpdate(double& ctransptc, double& crotptc, double& nc, double& resinit, double& alphaptc)
{
  double np;
  fres_->NormInf(&np);

  if(nc==0.0)
    dserror("nc == 0.0! PTC scheme not applicable! Choose fullnewton");

  // SER step size control
  crotptc *= std::pow((np/nc),alphaptc);
  ctransptc *= std::pow((np/nc),alphaptc);
  nc = np;

  Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
  int maxptciter = statmechparams.get<int>("MAXITERPTC",5);
  double resfrac = statmechparams.get<double>("RESLOWPTC",0.001);

  // modification: turn off ptc once residual is small enough
  if(np < resfrac*resinit || iter_ > maxptciter)
  {
    ctransptc = 0.0;
    crotptc = 0.0;
  }
  return;
}//STR::TimIntStatMech::PTCStatMechUpdate()



/*---------------------------------------------------------------------
 solution with line search algorithm                  hiermeier 08/13
---------------------------------------------------------------------*/
int STR::TimIntStatMech::NewtonLS()
{
  // The specific time integration has set the following
  // --> On #fres_ is the positive force residuum
  // --> On #stiff_ is the effective dynamic stiffness matrix

  int linsolve_error= 0;
  int fscontrol = 0;             // integer for a first step control (equal 1: deactivation) // fixme check if this is necessary for structural mechanics
  bool eval_error = false;       // an error occurred in the structure evaluation

  // check whether we have a sanely filled stiffness matrix
  if (not stiff_->Filled())
    dserror("Effective stiffness matrix must be filled here");

  if (outputeveryiter_)
  {
    int restart = DRT::Problem::Instance()->Restart();
    if (stepn_ == (restart + 1)) outputcounter_ = 0;
    OutputEveryIter(true);
  }

  // initialize equilibrium loop (outer Full Newton loop)
  iter_ = 1;
  normfres_ = CalcRefNormForce();
  // normdisi_ was already set in predictor; this is strictly >0
  timer_->ResetStartTime();

  // Merit function at current stage and for ls step
  std::vector<double> merit_fct (2);

  // Temporal copies of different vectors. Necessary for the sufficient decrease check.
  Teuchos::RCP<Epetra_Vector> tdisn = Teuchos::rcp(new Epetra_Vector(*disn_));
  Teuchos::RCP<Epetra_Vector> tveln = Teuchos::rcp(new Epetra_Vector(*veln_));
  Teuchos::RCP<Epetra_Vector> taccn = Teuchos::rcp(new Epetra_Vector(*accn_));

  // equilibrium iteration loop (outer full Newton loop)
  while ( ( (not Converged() and (not linsolve_error)) and (iter_ <= itermax_) ) or (iter_ <= itermin_) )
  {
    // initialize the Newton line search iteration counter
    int iter_ls  = 0;
    double step_red = 1.0;

    /*************************************************************
    ***           Save successful iteration state               ***
    **************************************************************/

    // It's necessary to save a temporal copy of the end-point displacements,
    // before any update is performed (because of the pseudo energy norm):
    tdisn->Update(1.0, *disn_, 0.0);                  // copy of the displ vector
    tveln->Update(1.0, *veln_, 0.0);                  // copy of the velocity vector
    taccn->Update(1.0, *accn_, 0.0);                  // copy of the acceleration vector

    /*************************************************************
    ***                       Solver Call                       ***
    **************************************************************/
    linsolve_error = LsSolveNewtonStep();

    // Evaluate merit function
    if (iter_==1)
      LsEvalMeritFct(merit_fct[0]);
    else
      merit_fct[0]=merit_fct[1];

    // Check if pred_constdis is used. If yes, the first step is not controlled.
    if (pred_ == INPAR::STR::pred_constdis or pred_ == INPAR::STR::pred_constdisvelacc)
      fscontrol = 1;
    else if ((pred_==INPAR::STR::pred_tangdis || pred_==INPAR::STR::pred_constacc ||
             pred_==INPAR::STR::pred_constvel)|| (iter_ > 1))
      fscontrol=0;
    else
      dserror("The behavior of the chosen predictor is not yet tested in the line search framework.");

    /*************************************************************
    ***      Update right-hand side and stiffness matrix        ***
    **************************************************************/
    Teuchos::ParameterList params;
    params.set<bool>("tolerate_errors",true);
    params.set<bool>("eval_error",false);
    if (fresn_str_!=Teuchos::null)
    {
      // attention: though it is called rhs_norm it actually contains sum x_i^2, i.e. the square of the L2-norm
      params.set<double>("cond_rhs_norm",0.);
      // need to know the processor id
      params.set<int>("MyPID",myrank_);
    }
    {
      int exceptcount = 0;
      fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
      EvaluateForceStiffResidual(params);
      if (fetestexcept(FE_INVALID) || fetestexcept(FE_OVERFLOW)
          || fetestexcept(FE_DIVBYZERO) || params.get<bool>("eval_error")==true)
        exceptcount  = 1;
      int tmp=0;
      discret_->Comm().SumAll(&exceptcount,&tmp,1);
      if (tmp)
        eval_error=true;
      feclearexcept(FE_ALL_EXCEPT);
      feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    }

    // get residual of condensed variables (e.g. EAS) for NewtonLS
    if (fresn_str_!=Teuchos::null)
    {
      double loc=params.get<double>("cond_rhs_norm");
      discret_->Comm().SumAll(&loc,&cond_res_,1);
    }

    // blank residual at (locally oriented) Dirichlet DOFs
    // rotate to local co-ordinate systems
    if (locsysman_ != Teuchos::null)
      locsysman_->RotateGlobalToLocal(fres_);

    // extract reaction forces
    // reactions are negative to balance residual on DBC
    freact_->Update(-1.0, *fres_, 0.0);
    dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);
    // rotate reaction forces back to global co-ordinate system
    if (locsysman_ != Teuchos::null)
      locsysman_->RotateLocalToGlobal(freact_);

    // blank residual at DOFs on Dirichlet BC
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);
    // rotate back to global co-ordinate system
    if (locsysman_ != Teuchos::null)
      locsysman_->RotateLocalToGlobal(fres_);

    // cancel in residual those forces that would excite rigid body modes and
    // that thus vanish in the Krylov space projection
    if (projector_!=Teuchos::null)
      projector_->ApplyPT(*fres_);

    /*************************************************************
    ***           merit function (current iteration)            ***
    **************************************************************/
    int err=LsEvalMeritFct(merit_fct[1]);
    eval_error = (eval_error || err);

    if (outputeveryiter_) OutputEveryIter(true);

    /*************************************************************
    ***          1st inner LINE SEARCH loop                     ***
    **************************************************************/

    while ((iter_-fscontrol > 0) && ((!LsConverged(& merit_fct[0], step_red) || eval_error) && (iter_ls < ls_maxiter_)))
    {
      /*************************************************************
      ***           Display line search information               ***
      **************************************************************/
      if (iter_ls==0)
        LsPrintLineSearchIter(&merit_fct[0],iter_ls,step_red);

      // increase inner loop count
      ++iter_ls;

      /*************************************************************
      ***                   Step size control                     ***
      **************************************************************/
      step_red *= alpha_ls_;
      // >>>> displacement, velocity, acceleration <<<<<<<<<<<<<<<
      // scale displ. increment
      disi_->Scale(alpha_ls_);
      // load old displ. vector
      disn_->Update(1.0, *tdisn, 0.0);
      // load old vel. vector
      veln_->Update(1.0, *tveln, 0.0);
      // load old acc. vector
      accn_->Update(1.0, *taccn, 0.0);

      // Update nodal displ., vel., acc., etc.
      UpdateIter(iter_);
      /*************************************************************
      ***   Update right-hand side (and part. stiffness matrix)   ***
      **************************************************************/
      LsUpdateStructuralRHSandStiff(eval_error,merit_fct[1]);

      /*************************************************************
      ***           Display line search information               ***
      **************************************************************/
      LsPrintLineSearchIter(&merit_fct[0],iter_ls,step_red);

      if (!(eval_error) && (outputeveryiter_)) OutputEveryIter(true, true);
    }

    if (iter_ls!=0)
    {
      if ( (myrank_ == 0) and printscreen_ and (StepOld()%printscreen_==0) and  printiter_ )
      {
        std::ostringstream oss;
        std::string dashline;
        dashline.assign(64,'-');
        oss << dashline ;
        // print to screen (could be done differently...)
        if (printerrfile_)
        {
          fprintf(errfile_, "%s\n", oss.str().c_str());
          fflush(errfile_);
        }

        fprintf(stdout, "%s\n", oss.str().c_str());
        fflush(stdout);
      }
    }

    /*************************************************************
    ***      Print Newton Step information                      ***
    **************************************************************/

    // build residual force norm
    normfres_ = STR::AUX::CalculateVectorNorm(iternorm_, fres_);
    // build residual displacement norm
    normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);

    PrintNewtonIter();

    // increment equilibrium loop index
    iter_ += 1;

    // DEBUG
   // GmshOutputEveryIter();
  } // end equilibrium loop

  // correct iteration counter
  iter_ -= 1;

  // call monitor
  if (conman_->HaveMonitor())
    conman_->ComputeMonitorValues(disn_);

  //do nonlinear solver error check
  return NewtonFullErrorCheck(linsolve_error,0);
}


/*----------------------------------------------------------------------
   Solver Call (line search)                          hiermeier 09/13
----------------------------------------------------------------------*/
int STR::TimIntStatMech::LsSolveNewtonStep()
{
  int linsolve_error = 0;
  /*************************************************************
  ***           Prepare the solution procedure                ***
  **************************************************************/
  // make negative residual
  fres_->Scale(-1.0);

  // transform to local co-ordinate systems
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateGlobalToLocal(SystemMatrix(), fres_);

  // STC preconditioning
  STCPreconditioning();

  // apply Dirichlet BCs to system of equations
  disi_->PutScalar(0.0);  // Useful? depends on solver and more
  LINALG::ApplyDirichlettoSystem(stiff_, disi_, fres_,
                                 GetLocSysTrafo(), zeros_, *(dbcmaps_->CondMap()));

  /*************************************************************
  ***                     Solver Call                         ***
  **************************************************************/
  // *********** time measurement ***********
  double dtcpu = timer_->WallTime();
  // *********** time measurement ***********

  // solve for disi_
  // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
  if (solveradapttol_ and (iter_ > 1))
  {
    double worst = normfres_;
    double wanted = tolfres_;
    solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
  }

  linsolve_error = solver_->Solve(stiff_->EpetraOperator(), disi_, fres_, true, iter_==1, projector_);
  // check for problems in linear solver
  // however we only care about this if we have a fancy divcont action (meaning function will return 0 )
  linsolve_error = LinSolveErrorCheck(linsolve_error);

  //In beam contact applications it can be necessary to limit the Newton step size (scaled residual displacements)
  LimitStepsizeBeamContact(disi_);

  solver_->ResetTolerance();

  // recover standard displacements
  RecoverSTCSolution();

  // *********** time measurement ***********
  dtsolve_ = timer_->WallTime() - dtcpu;
  // *********** time measurement ***********

  // update end-point displacements etc
  UpdateIter(iter_);

  return(linsolve_error);
}

/*----------------------------------------------------------------------
   Update structural RHS and stiff (line search)      hiermeier 09/13
----------------------------------------------------------------------*/
void STR::TimIntStatMech::LsUpdateStructuralRHSandStiff(bool& isexcept,double& merit_fct)
{
  // --- Checking for floating point exceptions
  fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

  // compute residual forces #fres_ and stiffness #stiff_
  // whose components are globally oriented
  int exceptcount = 0;
  Teuchos::ParameterList params;
  // elements may tolerate errors usually leading to dserrors
  // in such cases the elements force the line search to reduce
  // the step size by setting "eval_error" to true
  params.set<bool>("tolerate_errors",true);
  params.set<bool>("eval_error",false);
  // condensed degrees of freedom need to know the step reduction
  params.set<double>("alpha_ls",alpha_ls_);
  // line search needs to know the residuals of additional condensed dofs
  if (fresn_str_!=Teuchos::null)
  {
    params.set<double>("cond_rhs_norm",0.);
    // need to know the processor id
    params.set<int>("MyPID",myrank_);
  }
  EvaluateForceStiffResidual(params);

  // get residual of condensed variables (e.g. EAS) for NewtonLS
  if (fresn_str_!=Teuchos::null)
  {
    double loc=params.get<double>("cond_rhs_norm");
    discret_->Comm().SumAll(&loc,&cond_res_,1);
  }

  if (fetestexcept(FE_INVALID) || fetestexcept(FE_OVERFLOW)
      || fetestexcept(FE_DIVBYZERO) || params.get<bool>("eval_error")==true)
    exceptcount  = 1;

  // synchronize the exception flag isexcept on all processors
  int exceptsum = 0;
  discret_->Comm().SumAll(& exceptcount, & exceptsum, 1);
  if (exceptsum > 0)
    isexcept = true;
  else
    isexcept=false;

  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
  feclearexcept(FE_ALL_EXCEPT);
  // blank residual at (locally oriented) Dirichlet DOFs
  // rotate to local co-ordinate systems
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateGlobalToLocal(fres_);

  // extract reaction forces
  // reactions are negative to balance residual on DBC
  freact_->Update(-1.0, *fres_, 0.0);
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);
  // rotate reaction forces back to global co-ordinate system
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateLocalToGlobal(freact_);

  // blank residual at DOFs on Dirichlet BC
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);
  // rotate back to global co-ordinate system
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateLocalToGlobal(fres_);

  // cancel in residual those forces that would excite rigid body modes and
  // that thus vanish in the Krylov space projection
  if (projector_!=Teuchos::null)
    projector_->ApplyPT(*fres_);

  /*************************************************************
  ***          merit function (current iteration)             ***
  **************************************************************/
  int err=LsEvalMeritFct(merit_fct);
  isexcept = (isexcept || err);

  return;
}


/*----------------------------------------------------------------------
   Evaluate the merit function (line search)          hiermeier 08/13
----------------------------------------------------------------------*/
int STR::TimIntStatMech::LsEvalMeritFct(double& merit_fct)
{
  fedisableexcept(FE_OVERFLOW);
  int err=0;
  // Calculate the quadratic norm of the right-hand side as merit function
  // Calculate the merit function value: (1/2) * <RHS,RHS>
  if (fresn_str_==Teuchos::null)
  {
    err=fres_->Dot(*fres_,& merit_fct);
  }
  else
  {
    merit_fct=0.;
    err=fresn_str_->Dot(*fresn_str_,&merit_fct);
    merit_fct+=cond_res_;
  }
  merit_fct *= 0.5 ;

  int exceptcount=0;
  if (fetestexcept(FE_OVERFLOW))
      exceptcount  = 1;
  int exceptsum = 0;
  discret_->Comm().SumAll(& exceptcount, & exceptsum, 1);
  if (exceptsum!=0)
    return err;
  feclearexcept(FE_ALL_EXCEPT);
  feenableexcept(FE_OVERFLOW);

  return 0;
}

/*----------------------------------------------------------------------
   Print information about the last line search step  hiermeier 09/13
----------------------------------------------------------------------*/
void STR::TimIntStatMech::LsPrintLineSearchIter(double* mf_value, int iter_ls, double step_red)
{
    normdisi_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);
  // print to standard out
  if ( (myrank_ == 0) and printscreen_ and (StepOld()%printscreen_==0) and  printiter_ )
  {
    std::ostringstream oss;
    if (iter_ls== 0)
    {
      std::string dashline;
      dashline.assign(64,'-');
      oss << dashline << std::endl;
      oss << std::setw(6) << "ls_iter";
      oss << std::setw(16) << "step_scale";
      oss << std::setw(16) << "abs-dis-norm";
      oss << std::setw(16) << "merit-fct";
      oss << std::setw(10)<< "te";
      if (HaveSemiSmoothPlasticity())
        oss << std::setw(10) << "#active";
      oss << std::endl;
    }

    oss << std::setw(7) << iter_ls;
    oss << std::setw(16) << std::setprecision(5) << std::scientific << step_red;
    // build residual displacement norm
    oss << std::setw(16) << std::setprecision(5) << std::scientific << normdisi_;
    if (iter_ls==0)
      oss << std::setw(16) << std::setprecision(5) << std::scientific << mf_value[0];
    else
      oss << std::setw(16) << std::setprecision(5) << std::scientific << mf_value[1];
    oss << std::setw(10) << std::setprecision(2) << std::scientific << dtele_;
    if (HaveSemiSmoothPlasticity())
      oss << std::setw(10) << std::scientific << plastman_->NumActivePlasticGP();

    // finish oss
    oss << std::ends;

      // print to screen (could be done differently...)
    if (printerrfile_)
    {
      fprintf(errfile_, "%s\n", oss.str().c_str());
      fflush(errfile_);
    }

    fprintf(stdout, "%s\n", oss.str().c_str());
    fflush(stdout);
  }

  // see you
  return;
}

/*----------------------------------------------------------------------
   Inner convergence check (line search)              hiermeier 08/13
----------------------------------------------------------------------*/
bool STR::TimIntStatMech::LsConverged(double* mf_value,double step_red)
{
  bool check_ls_mf = false;

  /*************************************************************
  ***           Check for sufficient descent                  ***
  **************************************************************/
  // mf_value[1]: NEW merit function value
  //            --> f(x + alpha_ls * dx)
  // mf_value[0]: OLD merit function value (initial value at the beginning of the time step
  //              or function value of the last converged iteration step. Converged means that
  //              the last step fulfilled the LsConverged test.)
  //            --> f(x)
  // The check follows to
  //            f(x + alpha_ls * dx) - f(x) <= - 2 * sigma_ls * step_red_ * f(x).
  check_ls_mf = ((mf_value[1] - mf_value[0]) <= - 2.0 * sigma_ls_ * step_red * mf_value[0]);


//  std::cout<<"step_red="<<step_red<<std::endl;
//  std::cout<<"check_ls_mf="<<check_ls_mf<<std::endl;
  return (check_ls_mf);
}


/*----------------------------------------------------------------------*
 |  incremental iteration update of state                  mueller 03/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::UpdateIter(const int iter)  //!< iteration counter
{
  if (iter <= 1)
  {
    UpdateIterIncrementally();
  }
  else
  {
    UpdateIterIteratively();
  }
  // morning is broken
  return;
}

/*----------------------------------------------------------------------*
 |  incremental iteration update of state                  mueller 03/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::UpdateIterIncrementally()
{
  // Auxiliary vector holding new velocities and accelerations
  // by extrapolation/scheme on __all__ DOFs. This includes
  // the Dirichlet DOFs as well. Thus we need to protect those
  // DOFs of overwriting; they already hold the
  // correctly 'predicted', final values.

  // this version leads to a segmentation fault if the time step has to be repeated...
  //Teuchos::RCP<Epetra_Vector> aux = LINALG::CreateVector(*DofRowMapView(), false);

  Teuchos::RCP<Epetra_Vector> aux = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()), false));


  // new end-point displacements
  // D_{n+1}^{<k+1>} := D_{n+1}^{<k>} + IncD_{n+1}^{<k>}
  disn_->Update(1.0, *disi_, 1.0);

  // new end-point velocities
  aux->Update(1.0/(theta_*(*dt_)[0]), *disn_,
               -1.0/(theta_*(*dt_)[0]), *(*dis_)(0),
               0.0);
  aux->Update(-(1.0-theta_)/theta_, *(*vel_)(0), 1.0);
  // put only to free/non-DBC DOFs
  // new version of updating velocity vector
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(aux), veln_);

  // note: no accelerations in statmech...
//  // new end-point accelerations
//  aux->Update(1.0/(theta_*theta_*(*dt_)[0]*(*dt_)[0]), *disn_,
//              -1.0/(theta_*theta_*(*dt_)[0]*(*dt_)[0]), *(*dis_)(0),
//              0.0);
//  aux->Update(-1.0/(theta_*theta_*(*dt_)[0]), *(*vel_)(0),
//              -(1.0-theta_)/theta_, *(*acc_)(0),
//              1.0);
//  // put only to free/non-DBC DOFs
//  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(aux), accn_);

  // bye
  return;
}

/*----------------------------------------------------------------------*
 |  iterative iteration update of state                    mueller 03/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::UpdateIterIteratively()
{
  // new end-point displacements
  // D_{n+1}^{<k+1>} := D_{n+1}^{<k>} + IncD_{n+1}^{<k>}
  disn_->Update(1.0, *disi_, 1.0);

  // new end-point velocities
  veln_->Update(1.0/(theta_*(*dt_)[0]), *disi_, 1.0);

  // note: no accelerations in statmech...
//  // new end-point accelerations
//  accn_->Update(1.0/((*dt_)[0]*(*dt_)[0]*theta_*theta_), *disi_, 1.0);

  return;
}

/*----------------------------------------------------------------------*
 | set relevant variables signaling divergence of PTC      mueller 02/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::ConvergenceStatusUpdate(bool converged)
{
  if(converged)
  {
    isconverged_ = true;
  }
  else
  {
    isconverged_ = false;
    if(HaveStatMechBilayer())
    {
      statmechmanBilayer_->UpdateNumberOfUnconvergedSteps();
    }
    else
      statmechman_->UpdateNumberOfUnconvergedSteps();
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate beam contact according to solution strategy                |
 |                                            (private)    mueller 02/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::BeamContactNonlinearSolve()
{
  INPAR::BEAMCONTACT::Strategy strategy = DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Strategy>(beamcman_->BeamContactParameters(),"BEAMS_STRATEGY");
  switch (strategy)
  {
    //solving strategy using regularization with penalty method (nonlinear solution approach: ordinary NEWTON (PTC))
    case INPAR::BEAMCONTACT::bstr_penalty:
      BeamContactPenalty();
    break;
    //solving strategy using regularization with augmented Lagrange method (nonlinear solution approach: nested UZAWA NEWTON (PTC))
    case INPAR::BEAMCONTACT::bstr_uzawa:
    {
      BeamContactAugLag();
    }
    break;
    default:
      dserror("Only penalty and Uzawa augmented Lagrange implemented in statmech_time.cpp for beam contact");
      break;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate beam contact using Penalty appr. (private)    mueller 02/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::BeamContactPenalty()
{

  if(itertype_==INPAR::STR::soltech_ptc)
    PTC();
  else if(itertype_==INPAR::STR::soltech_newtonfull)
    NewtonFull();
  else
    dserror("itertype %d not implemented for StatMech applications! Choose either ptc or fullnewton!", itertype_);

  beamcman_->UpdateConstrNorm();
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate beam contact using Augmented Lagrange                      |
 |                                            (private)    mueller 02/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::BeamContactAugLag()
{
  // get tolerance and maximum number of Uzawa steps from input file
  double eps = beamcman_->BeamContactParameters().get<double>("BEAMS_BTBUZAWACONSTRTOL");
  int maxuzawaiter = beamcman_->BeamContactParameters().get<int>("BEAMS_BTBUZAWAMAXSTEPS");

  Teuchos::ParameterList ioparams = DRT::Problem::Instance()->IOParams();

  // LOOP2: augmented Lagrangian (Uzawa)
  do
  {
    // increase iteration index
    beamcman_->UpdateUzawaIter();

    // if unconverged
    if(BeamContactExitUzawaAt(maxuzawaiter))
    {
      dserror("Maximal number of uzawa iterations reached. Adapt your tolerance or allow for larger number of uzawa iterations!");
    }

    if (discret_->Comm().MyPID() == 0 && ioparams.get<int>("STDOUTEVRY",0))
      std::cout << std::endl << "Starting Uzawa step No. " << beamcman_->GetUzawaIter() << std::endl;

    if(itertype_==INPAR::STR::soltech_ptc)
      PTC();
    else if(itertype_==INPAR::STR::soltech_newtonfull)
      NewtonFull();
    else
      dserror("itertype %d not implemented for StatMech applications! Choose either ptc or fullnewton!", itertype_);

    // in case uzawa step did not converge
    if(!isconverged_)
    {
      if(!discret_->Comm().MyPID() && ioparams.get<int>("STDOUTEVRY",0))
        std::cout<<"\n\nNewton iteration in Uzawa Step "<<beamcman_->GetUzawaIter()<<" unconverged - leaving Uzawa loop and restarting time step...!\n\n";
      break;
    }

    // update constraint norm and penalty parameter
    beamcman_->UpdateConstrNormUzawa();
    // update Uzawa Lagrange multipliers
    beamcman_->UpdateAlllmuzawa();
  } while (abs(beamcman_->GetConstrNorm()) >= eps);

  // Reset Penalty parameter, Lagrange Multipliers and Uzawa iteration counter for the Augmented Lagrangian loop
  beamcman_->ResetCurrentpp();
  beamcman_->ResetUzawaIter();
  beamcman_->ResetAlllmuzawa();

  beamcman_->UpdateConstrNorm();
  return;
}

/*----------------------------------------------------------------------*
 |  Check Uzawa convergence                      (private) mueller 02/12|
 *----------------------------------------------------------------------*/
bool STR::TimIntStatMech::BeamContactExitUzawaAt(int& maxuzawaiter)
{
  if (beamcman_->GetUzawaIter() > maxuzawaiter)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------*
 | Precautions for Contact during one time step (private) mueller 03/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::StatMechPrepareStep()
{
  if(HaveStatMech() || HaveStatMechBilayer())
  {
    Teuchos::ParameterList statmechparams;
    if(HaveStatMech())
      statmechparams = statmechman_->GetStatMechParams();
    else if(HaveStatMechBilayer())
      statmechparams = statmechmanBilayer_->GetStatMechBilayerParams();

    switch (DRT::INPUT::IntegralValue<INPAR::STATMECH::SimulationType>(statmechparams,"SIMULATION_TYPE"))
    {
    case INPAR::STATMECH::simulation_type_biopolymer_network:
    case INPAR::STATMECH::simulation_type_none:
    {
//
      // special preparations for the very first step
      if(step_ == 0)
      {
        /* In case we add an initial amount of already linked crosslinkers, we have to build the octree
         * even before the first statmechman_->Update() call because the octree is needed to decide
         * whether links can be set...*/
        if(statmechparams.get<int>("INITOCCUPIEDBSPOTS",0)>0)
          statmechman_->SetInitialCrosslinkers(beamcman_);

        // Initialize Statistical Mechanics Output
        statmechman_->InitOutput(DRT::Problem::Instance()->NDim(),*((*dis_)(0)),step_,(*dt_)[0]);
      }

      //save relevant class variables at the beginning of this time step
      statmechman_->WriteConv(beamcman_);

      double randnumtimeinc = statmechparams.get<double>("RANDNUMTIMEINT",-1.0);
      int randnumupdatestep = static_cast<int>( (timen_-(*dt_)[0]) / randnumtimeinc + 1.0e-8);

      //seed random generators of statmechman_ to generate the same random numbers even if the simulation was interrupted by a restart
      statmechman_->SeedRandomGenerators(step_, randnumupdatestep);

      if(!discret_->Comm().MyPID() && printscreen_)
      {
        std::cout<<"\nbegin time step "<<stepn_<<":";
      }
    }
    break;
    case INPAR::STATMECH::simulation_type_lipid_bilayer:
    {
      // special preparations for the very first step
      if(step_ == 0)
      {
        // Initialize Statistical Mechanics Output
        // TODO: create output for lipid bilayer analysis
        statmechmanBilayer_->InitOutput(DRT::Problem::Instance()->NDim(),*((*dis_)(0)),step_,(*dt_)[0]);
      }


      double randnumtimeinc = statmechparams.get<double>("RANDNUMTIMEINT",-1.0);
      int randnumupdatestep = static_cast<int>( (timen_-(*dt_)[0]) / randnumtimeinc + 1.0e-8);

      //seed random generators of statmechman_ to generate the same random numbers even if the simulation was interrupted by a restart
      statmechmanBilayer_->SeedRandomGenerators(step_, randnumupdatestep);

      if(!discret_->Comm().MyPID() && printscreen_)
      {
        std::cout<<"\nbegin time step "<<stepn_<<":";
      }
    }
    break;
    case INPAR::STATMECH::simulation_type_network_bilayer:
    {
      // special preparations for the very first step
      if(step_ == 0)
      {
        // Initialize Statistical Mechanics Output
        // TODO: create output for network & lipid bilayer analysis
        statmechmanBilayer_->InitOutput(DRT::Problem::Instance()->NDim(),*((*dis_)(0)),step_,(*dt_)[0]);
      }


      double randnumtimeinc = statmechparams.get<double>("RANDNUMTIMEINT",-1.0);
      int randnumupdatestep = static_cast<int>( (timen_-(*dt_)[0]) / randnumtimeinc + 1.0e-8);

      //seed random generators of statmechman_ to generate the same random numbers even if the simulation was interrupted by a restart
      statmechmanBilayer_->SeedRandomGenerators(step_, randnumupdatestep);
      statmechman_->SeedRandomGenerators(step_, randnumupdatestep);
      if(!discret_->Comm().MyPID() && printscreen_)
      {
        std::cout<<"\nbegin time step "<<stepn_<<":";
      }
    }
    break;
    default:
      dserror("Simulation type not found! Please check input parameters.");
      break;
    }
  }

  return;
}//StatMechPrepareStep()

/*----------------------------------------------------------------------*
 |  Call statmechmanager Update according to options chosen             |
 |                                            (private)    mueller 02/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::StatMechUpdate(bool newrandomnumbers)
{
  if(!HaveStatMech() && !HaveStatMechBilayer())
    dserror("You should not be here!");
//

//  Teuchos::ParameterList statmechparams = statmechman_->GetStatMechParams();
  Teuchos::ParameterList statmechparams;
  if(HaveStatMech())
    statmechparams = statmechman_->GetStatMechParams();
  if(HaveStatMechBilayer()) // Get parameter for bilayer
    statmechparams = statmechmanBilayer_->GetStatMechBilayerParams();
  //assuming that iterations will converge
//

  isconverged_ = true;
  const double t_admin = Teuchos::Time::wallTime();
  if(HaveStatMech())
  {
    if(HaveBeamContact())
      statmechman_->Update(step_, timen_, (*dt_)[0], *((*dis_)(0)), stiff_,ndim_,beamcman_,buildoctree_, printscreen_);
    else
      statmechman_->Update(step_, timen_, (*dt_)[0], *((*dis_)(0)), stiff_,ndim_, Teuchos::null,false,printscreen_);
  }
  else if(HaveStatMechBilayer())
    statmechmanBilayer_->Update(step_, timen_, (*dt_)[0], *((*dis_)(0)), stiff_, false, printscreen_);
  // print to screen
  StatMechPrintUpdate(t_admin);
//

  //Only generate new random numbers if necessary (for repeated time steps this depends on the type of divercont_ action)
  if(newrandomnumbers)
  {
    /*multivector for stochastic forces evaluated by each element; the numbers of vectors in the multivector equals the maximal
     *number of random numbers required by any element in the discretization per time step; therefore this multivector is suitable
     *for synchrinisation of these random numbers in parallel computing*/
    if(HaveStatMech())
      randomnumbers_ = Teuchos::rcp( new Epetra_MultiVector(*(discret_->ElementColMap()),maxrandomnumbersperglobalelement_,true) );
    else if(HaveStatMechBilayer())
      randomnumbers_ = Teuchos::rcp( new Epetra_MultiVector(*(discret_->NodeColMap()),maxrandomnumbersperglobalnode_,true) );

    /*pay attention: for a constant predictor an incremental velocity update is necessary, which has been deleted out of the code in oder to simplify it*/
    //generate gaussian random numbers for parallel use with mean value 0 and standard deviation (2KT / dt)^0.5
    double randnumtimeinc = statmechparams.get<double>("RANDNUMTIMEINT",-1.0);
    if(randnumtimeinc==-1.0)
    {
      if(HaveStatMech())
        statmechman_->GenerateGaussianRandomNumbers(randomnumbers_,0,pow(2.0 * statmechparams.get<double>("KT",0.0) / (*dt_)[0],0.5));
      else if(HaveStatMechBilayer())
        statmechmanBilayer_->GenerateGaussianRandomNumbers(randomnumbers_,0,pow(2.0 * statmechparams.get<double>("KT",0.0) / (*dt_)[0],0.5));
    }
    else
    {
      if(HaveStatMech())
        statmechman_->GenerateGaussianRandomNumbers(randomnumbers_,0,pow(2.0 * statmechparams.get<double>("KT",0.0) / randnumtimeinc ,0.5));
      else if(HaveStatMechBilayer())
        statmechmanBilayer_->GenerateGaussianRandomNumbers(randomnumbers_,0,pow(2.0 * statmechparams.get<double>("KT",0.0) / randnumtimeinc ,0.5));
    }
  }

  return;
} //StatMechUpdate()

/*----------------------------------------------------------------------*
 | Print Statistical Mechanics Update To Screen (private)  mueller 03/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::StatMechPrintUpdate(const double& t_admin)
{
  if(!myrank_ && printscreen_)
  {
    if(HaveStatMech())
    {
      INPAR::STATMECH::LinkerModel linkermodel = DRT::INPUT::IntegralValue<INPAR::STATMECH::LinkerModel>(DRT::Problem::Instance()->StatisticalMechanicsParams(),"LINKERMODEL");
      if(linkermodel != INPAR::STATMECH::linkermodel_none)
      {
        std::cout<<"\nTime for update of crosslinkers                   : " << Teuchos::Time::wallTime() - t_admin<< " seconds";
        std::cout<<"\nTotal number of elements after crosslinker update : "<<discret_->NumGlobalElements();
      }
      std::cout<<"\nNumber of unconverged steps since simulation start: "<<statmechman_->NumberOfUnconvergedSteps()<<"\n"<<std::endl;
    }
    else if(HaveStatMechBilayer())
    {
      std::cout<<"\nNumber of unconverged steps since simulation start: "<<statmechmanBilayer_->NumberOfUnconvergedSteps()<<"\n"<<std::endl;
    }
  }
  return;
}//StatMechPrintUpdate()

/*----------------------------------------------------------------------*
 |  Call statmechmanager Output according to options chosen             |
 |                                            (private)    mueller 02/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::StatMechOutput()
{
  // note: "step_ - 1" in order to make the modulo operations within Output() work properly.
  if(HaveStatMech())
  {
    if(HaveBeamContact())
      statmechman_->Output(ndim_,(*time_)[0],step_,(*dt_)[0],*((*dis_)(0)),*fint_,beamcman_, printscreen_);
    else
      statmechman_->Output(ndim_,(*time_)[0],step_,(*dt_)[0],*((*dis_)(0)),*fint_, Teuchos::null, printscreen_);
  }
  if(HaveStatMechBilayer())
  {
    statmechmanBilayer_->Output(ndim_,(*time_)[0],step_,(*dt_)[0],*((*dis_)(0)),*fint_, printscreen_);
  }
  return;
}// StatMechOutput()

/*----------------------------------------------------------------------*
 |  Reset relevant values, vectors, discretization before repeating     |
 |  the time step                             (private)    mueller 02/12|
 *----------------------------------------------------------------------*/
void STR::TimIntStatMech::StatMechRestoreConvState()
{
  if(HaveStatMech())
  {
    if(!isconverged_)
    {
      Teuchos::ParameterList p;
      p.set("action","calc_struct_reset_istep");
      discret_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
      statmechman_->RestoreConv((*dis_)(0), stiff_, beamcman_,printscreen_);
      buildoctree_ = true;
    }
  }
  return;
} // StatMechRestoreConvState()

void STR::TimIntStatMech::PrintTotalSurfAreaperIter()
{
  //we need displacements also of ghost nodes and hence export displacment vector to column map format
  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(*disn_, discol);

  double area_ele=0;
  double area_total=0;
  for (int i=0; i < discret_->NumMyRowElements() - 1; i++)
  {
    //getting pointer to current element
    DRT::Element* element = discret_->lRowElement(i);
    const DRT::ElementType & eot = element->ElementType();

    if(eot == DRT::ELEMENTS::DiscSh3Type::Instance())    // discrete shell element
    {
      DRT::ELEMENTS::DiscSh3* ele = dynamic_cast<DRT::ELEMENTS::DiscSh3*>(element);
      area_ele=ele->CalcSurfArea(*discret_,discol,false);
    }
    discret_->Comm().SumAll(&area_ele, &area_total, 1);
  }

  //proc 0 write complete output into file, all other proc inactive
  if(!discret_->Comm().MyPID())
  {

    FILE* fp = NULL; //file pointer for statistical output file

    //name of output file
    std::ostringstream outputfilename;
    outputfilename.str("");
    outputfilename << "/home/mukherjee/workspace/baci/release/StatMechOutput/SurfaceAreaIterations" << ".dat";

    fp = fopen(outputfilename.str().c_str(), "a");
    std::stringstream filecontent;
    //          filecontent << iter_-1;
    filecontent << std::scientific << std::setprecision(10);
    filecontent << " " << area_total;

    filecontent << std::endl;
    fputs(filecontent.str().c_str(), fp);
    fclose(fp);
  }

  return;

}

void STR::TimIntStatMech::PrintTotalVolperIter()
{
  //we need displacements also of ghost nodes and hence export displacment vector to column map format
  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(*disn_, discol);

  double vol_ele=0;
  double vol_total=0;
  for (int i=0; i < discret_->NumMyRowElements() - 1; i++)
  {
    //getting pointer to current element
    DRT::Element* element = discret_->lRowElement(i);
    const DRT::ElementType & eot = element->ElementType();

    if(eot == DRT::ELEMENTS::DiscSh3Type::Instance())    // discrete shell element
    {
      DRT::ELEMENTS::DiscSh3* ele = dynamic_cast<DRT::ELEMENTS::DiscSh3*>(element);
      vol_ele=ele->CalcVolume(*discret_,discol,false);
    }
    discret_->Comm().SumAll(&vol_ele, &vol_total, 1);
  }

  //proc 0 write complete output into file, all other proc inactive
  if(!discret_->Comm().MyPID())
  {

    FILE* fp = NULL; //file pointer for statistical output file

    //name of output file
    std::ostringstream outputfilename;
    outputfilename.str("");
    outputfilename << "/home/mukherjee/workspace/baci/release/StatMechOutput/VolumeperIterations" << ".dat";

    fp = fopen(outputfilename.str().c_str(), "a");
    std::stringstream filecontent;
    //          filecontent << iter_-1;
    filecontent << std::scientific << std::setprecision(10);
    filecontent << " " << vol_total;

    filecontent << std::endl;
    fputs(filecontent.str().c_str(), fp);
    fclose(fp);
  }

  return;

}

void STR::TimIntStatMech::PrintTotalEdgeLenperIter()
{
  //we need displacements also of ghost nodes and hence export displacment vector to column map format
  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(*disn_, discol);

  double Edge_ele=0;
  double Edge_total=0;
  for (int i=0; i < facediscret_->NumMyRowFaces() - 1; i++)
  {
    //getting pointer to current element
    DRT::Element* element = facediscret_->lRowFace(i);
    const DRT::ElementType & eot = element->ElementType();

    if(eot == DRT::ELEMENTS::DiscSh3LineType::Instance())    // discrete shell Edge element
    {
      DRT::ELEMENTS::DiscSh3Line* ele = dynamic_cast<DRT::ELEMENTS::DiscSh3Line*>(element);
      Edge_ele=ele->GetCurrEdgeLength(*discret_,discol);
    }
    discret_->Comm().SumAll(&Edge_ele, &Edge_total, 1);
  }

  //proc 0 write complete output into file, all other proc inactive
  if(!discret_->Comm().MyPID())
  {

    FILE* fp = NULL; //file pointer for statistical output file

    //name of output file
    std::ostringstream outputfilename;
    outputfilename.str("");
    outputfilename << "/home/mukherjee/workspace/baci/release/StatMechOutput/EdgeperIterations" << ".dat";

    fp = fopen(outputfilename.str().c_str(), "a");
    std::stringstream filecontent;
    //          filecontent << iter_-1;
    filecontent << std::scientific << std::setprecision(10);
    filecontent << " " << Edge_total;

    filecontent << std::endl;
    fputs(filecontent.str().c_str(), fp);
    fclose(fp);
  }

  return;
}

void STR::TimIntStatMech::GmshOutputEveryIter()
{
      // STEP 1: OUTPUT OF TIME STEP INDEX
      std::ostringstream filename;
      filename << "/home/mukherjee/workspace/baci/release/GmshOutput/";

      if (iter_<1000000)
        filename << "Timestep" << std::setw(4) << std::setfill('0')<<step_+1<<"Iter"<<std::setw(2) << std::setfill('0') << iter_-1;
      else /*(timestep>=1000000)*/
        dserror("ERROR: Gmsh output implemented for max 999.999 time steps");

      // finish filename
      filename<<".pos";

//        statmechman_->GmshOutput((*disn_),filename,step_,(*time_)[0]);
      statmechmanBilayer_->GmshOutput((*disn_),filename,step_,(*time_)[0]);

  return;

}

