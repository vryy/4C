/*----------------------------------------------------------------------*/
/*!
\file scatra_cardiac_monodomain_dyn.cpp

\brief entry point for cardiac monodomain scalar transport problems

\level 2

<pre>
\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*/
/*----------------------------------------------------------------------*/
#include "scatra_cardiac_monodomain_dyn.H"

#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_particle/binning_strategy.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_scatra_ele/scatra_ele.H"

#include <iostream>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "scatra_algorithm.H"
#include "scatra_resulttest.H"
#include "scatra_utils_clonestrategy.H"


/*----------------------------------------------------------------------*
 * Main control routine for scalar transport problems, incl. various solvers
 *
 *        o Laplace-/ Poisson equation (zero velocity field)
 *          (with linear and nonlinear boundary conditions)
 *        o transport of passive scalar in velocity field given by spatial function
 *        o transport of passive scalar in velocity field given by Navier-Stokes
 *          (one-way coupling)
 *        o scalar transport in velocity field given by Navier-Stokes with natural convection
 *          (two-way coupling)
 *
 *----------------------------------------------------------------------*/
void scatra_cardiac_monodomain_dyn(int restart)
{
  // pointer to problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // access the communicator
  const Epetra_Comm& comm = problem->GetDis("fluid")->Comm();

//  // print problem type
  if (comm.MyPID() == 0)
  {
    std::cout << "###################################################"<< std::endl;
    std::cout << "# YOUR PROBLEM TYPE: " << DRT::Problem::Instance()->ProblemName() << std::endl;
    std::cout << "###################################################" << std::endl;
  }


  // print CardiacMonodomain-Logo to screen
  if (comm.MyPID()==0) printheartlogo();

  // access the problem-specific parameter list
  const Teuchos::ParameterList& scatradyn = problem->ScalarTransportDynamicParams();

  // access the fluid discretization
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  // access the scatra discretization
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");

  // ensure that all dofs are assigned in the right order; this creates dof numbers with
  //       fluid dof < scatra dof
  fluiddis->FillComplete();
  scatradis->FillComplete();

  // set velocity field
  const INPAR::SCATRA::VelocityField veltype
    = DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatradyn,"VELOCITYFIELD");
  switch (veltype)
  {
    case INPAR::SCATRA::velocity_zero:  // zero  (see case 1)
    case INPAR::SCATRA::velocity_function:  // function
    case INPAR::SCATRA::velocity_function_and_curve: // function and time curve
    {
      // we directly use the elements from the scalar transport elements section
      if (scatradis->NumGlobalNodes()==0)
        dserror("No elements in the ---TRANSPORT ELEMENTS section");

      // add proxy of velocity related degrees of freedom to scatra discretization
      if (scatradis->BuildDofSetAuxProxy(DRT::Problem::Instance()->NDim()+1, 0, 0, true ) != 1)
        dserror("Scatra discretization has illegal number of dofsets!");

      // allow TRANSPORT conditions, too
      // NOTE: we can not use the conditions given by 'conditions_to_copy = clonestrategy.ConditionsToCopy()'
      // since we than may have some scatra condition twice. So we only copy the Dirichlet and Neumann conditions:
      std::map<std::string,std::string> conditions_to_copy;
      conditions_to_copy.insert(std::pair<std::string,std::string>("TransportDirichlet","Dirichlet"));
      conditions_to_copy.insert(std::pair<std::string,std::string>("TransportPointNeumann","PointNeumann"));
      conditions_to_copy.insert(std::pair<std::string,std::string>("TransportLineNeumann","LineNeumann"));
      conditions_to_copy.insert(std::pair<std::string,std::string>("TransportSurfaceNeumann","SurfaceNeumann"));
      conditions_to_copy.insert(std::pair<std::string,std::string>("TransportVolumeNeumann","VolumeNeumann"));

      DRT::UTILS::DiscretizationCreatorBase creator;
      creator.CopyConditions(*scatradis,*scatradis,conditions_to_copy);

      // finalize discretization
      scatradis->FillComplete();



      // We have do use the binningstrategy and extended ghosting when we use p-adaptivity. This guarantees that the elements
      // at the boarder between the processor calculate correctly since one face is shared with the neighboring element
      // (which is owned by an other processor = ghosted element) which again is sharing other faces with elements on other
      // processors (extended ghosted element)
      if(DRT::INPUT::IntegralValue<bool>(scatradyn,"PADAPTIVITY"))
      {
        // redistribute discr. with help of binning strategy
        if(scatradis->Comm().NumProc()>1)
        {
          // create vector of discr.
          std::vector<Teuchos::RCP<DRT::Discretization> > dis;
          dis.push_back(scatradis);

          //binning strategy for parallel redistribution
          Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy = Teuchos::null;

          std::vector<Teuchos::RCP<Epetra_Map> > stdelecolmap;
          std::vector<Teuchos::RCP<Epetra_Map> > stdnodecolmap;

          /// binning strategy is created and parallel redistribution is performed
          binningstrategy = Teuchos::rcp(new BINSTRATEGY::BinningStrategy(dis,stdelecolmap,stdnodecolmap));
        }
      }

      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == (-1))
        dserror("no linear solver defined for SCALAR_TRANSPORT problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

      // create instance of scalar transport basis algorithm (empty fluid discretization)
      Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatraonly = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());

      // now we can call Init() on the base algo.
      // time integrator is constructed and initialized inside
      scatraonly->Init(
          scatradyn,
          scatradyn,
          DRT::Problem::Instance()->SolverParams(linsolvernumber));

      // NOTE : At this point we may redistribute and/or
      //        ghost our discretizations at will.

      // now we must call Setup()
      scatraonly->Setup();

      // read the restart information, set vectors and variables
      if (restart) scatraonly->ScaTraField()->ReadRestart(restart);

      // set initial velocity field
      // note: The order ReadRestart() before SetVelocityField() is important here!!
      // for time-dependent velocity fields, SetVelocityField() is additionally called in each PrepareTimeStep()-call
      (scatraonly->ScaTraField())->SetVelocityField(1);

      // enter time loop to solve problem with given convective velocity
      (scatraonly->ScaTraField())->TimeLoop();

      // perform the result test if required
      DRT::Problem::Instance()->AddFieldTest(scatraonly->CreateScaTraFieldTest());
      DRT::Problem::Instance()->TestAll(comm);

      break;
    }
    case INPAR::SCATRA::velocity_Navier_Stokes:  // Navier_Stokes
    {
      dserror("Navier Stokes case not implemented for cardiac monodomain scalar transport problem");

      // we use the fluid discretization as layout for the scalar transport discretization
      if (fluiddis->NumGlobalNodes()==0) dserror("Fluid discretization is empty!");

      const INPAR::SCATRA::FieldCoupling fieldcoupling = DRT::INPUT::IntegralValue<INPAR::SCATRA::FieldCoupling>(DRT::Problem::Instance()->ScalarTransportDynamicParams(),"FIELDCOUPLING");

      // create scatra elements if the scatra discretization is empty
      if (scatradis->NumGlobalNodes()==0)
      {
        if (fieldcoupling != INPAR::SCATRA::coupling_match)
          dserror("If you want matching fluid and scatra meshes, do clone you fluid mesh and use FIELDCOUPLING match!");

        fluiddis->FillComplete();
        scatradis->FillComplete();

        // fill scatra discretization by cloning fluid discretization
        DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(fluiddis,scatradis);

        // set implementation type of cloned scatra elements
        for(int i=0; i<scatradis->NumMyColElements(); ++i)
        {
          DRT::ELEMENTS::Transport* element = dynamic_cast<DRT::ELEMENTS::Transport*>(scatradis->lColElement(i));
          if(element == NULL)
            dserror("Invalid element type!");
          else
            element->SetImplType(INPAR::SCATRA::impltype_std);
        }

        // add proxy of fluid transport degrees of freedom to scatra discretization
        if(scatradis->AddDofSet(fluiddis->GetDofSetProxy()) != 1)
          dserror("Scatra discretization has illegal number of dofsets!");
      }
      else
      {
        if (fieldcoupling != INPAR::SCATRA::coupling_volmortar)
          dserror("If you want non-matching fluid and scatra meshes, you need to use FIELDCOUPLING volmortar!");

        // allow TRANSPORT conditions, too
        std::map<std::string,std::string> conditions_to_copy;
        SCATRA::ScatraFluidCloneStrategy clonestrategy;
        conditions_to_copy = clonestrategy.ConditionsToCopy();
        DRT::UTILS::DiscretizationCreatorBase creator;
        creator.CopyConditions(*scatradis,*scatradis,conditions_to_copy);

        //first call FillComplete for single discretizations.
        //This way the physical dofs are numbered successively
        fluiddis->FillComplete();
        scatradis->FillComplete();

        //build auxiliary dofsets, i.e. pseudo dofs on each discretization
        const int ndofpernode_scatra = scatradis->NumDof(0,scatradis->lRowNode(0));
        const int ndofperelement_scatra  = 0;
        const int ndofpernode_fluid = fluiddis->NumDof(0,fluiddis->lRowNode(0));
        const int ndofperelement_fluid = 0;
        if (fluiddis->BuildDofSetAuxProxy(ndofpernode_scatra, ndofperelement_scatra, 0, true ) != 1)
          dserror("unexpected dof sets in fluid field");
        if (scatradis->BuildDofSetAuxProxy(ndofpernode_fluid, ndofperelement_fluid, 0, true) != 1)
          dserror("unexpected dof sets in scatra field");

        //call AssignDegreesOfFreedom also for auxiliary dofsets
        //note: the order of FillComplete() calls determines the gid numbering!
        // 1. fluid dofs
        // 2. scatra dofs
        // 3. fluid auxiliary dofs
        // 4. scatra auxiliary dofs
        fluiddis->FillComplete(true, false,false);
        scatradis->FillComplete(true, false,false);

        //NOTE: we have do use the binningstrategy here since we build our fluid and scatra problems by inheritance,
        //i.e. by calling the constructor of the corresponding class. But since we have to use the binning-strategy before
        //creating the single field we have to do it here :-( We would prefer to to it like the SSI since than we could
        //extended ghosting
        //TODO (thon): make this if-case obsolete and allow for redistribution within volmortar->Setup() by removing inheitance-building of fields
        {
          // redistribute discr. with help of binning strategy
          if(fluiddis->Comm().NumProc()>1)
          {
            // create vector of discr.
            std::vector<Teuchos::RCP<DRT::Discretization> > dis;
            dis.push_back(fluiddis);
            dis.push_back(scatradis);

            //binning strategy for parallel redistribution
            Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy = Teuchos::null;

            std::vector<Teuchos::RCP<Epetra_Map> > stdelecolmap;
            std::vector<Teuchos::RCP<Epetra_Map> > stdnodecolmap;

            /// binning strategy is created and parallel redistribution is performed
            binningstrategy = Teuchos::rcp(new BINSTRATEGY::BinningStrategy(dis,stdelecolmap,stdnodecolmap));
          }
        }
      }

      // support for turbulent flow statistics
      const Teuchos::ParameterList& fdyn = (DRT::Problem::Instance()->FluidDynamicParams());

      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == (-1))
        dserror("no linear solver defined for SCALAR_TRANSPORT problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

      // create a scalar transport algorithm instance
      Teuchos::RCP<SCATRA::ScaTraAlgorithm> algo =
          Teuchos::rcp(new SCATRA::ScaTraAlgorithm(
              comm,
              scatradyn,
              fdyn,
              "scatra",
              DRT::Problem::Instance()->SolverParams(linsolvernumber)));

      // init algo (init fluid time integrator and scatra time integrator inside)
      algo->Init(
          scatradyn,
          scatradyn,
          DRT::Problem::Instance()->SolverParams(linsolvernumber) );

      // setup algo (setup fluid time integrator and scatra time integrator inside)
      algo->Setup();

      // read restart information
      // in case a inflow generation in the inflow section has been performed, there are not any
      // scatra results available and the initial field is used
      if (restart)
      {
        if ((DRT::INPUT::IntegralValue<int>(fdyn.sublist("TURBULENT INFLOW"),"TURBULENTINFLOW")==true) and
           (restart==fdyn.sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP")))
          algo->ReadInflowRestart(restart);
        else
          algo->ReadRestart(restart);
      }
      else
        if(DRT::INPUT::IntegralValue<int>(fdyn.sublist("TURBULENT INFLOW"),"TURBULENTINFLOW")==true)
          dserror("Turbulent inflow generation for passive scalar transport should be performed as fluid problem!");

      // solve the whole scalar transport problem
      algo->TimeLoop();

      // summarize the performance measurements
      Teuchos::TimeMonitor::summarize();

      // perform the result test
      DRT::Problem::Instance()->AddFieldTest(algo->FluidField()->CreateFieldTest());
      DRT::Problem::Instance()->AddFieldTest(algo->CreateScaTraFieldTest());
      DRT::Problem::Instance()->TestAll(comm);

      break;
    } // case 2
    default:
    {
      dserror("unknown velocity field type for transport of passive scalar");
      break;
    }
  }

  return;

} // end of scatra_cardiac_monodomain_dyn()


/*----------------------------------------------------------------------*/
// print ELCH-Module logo
/*----------------------------------------------------------------------*/
void printheartlogo()
{
  // more at http://www.ascii-art.de
  std::cout<<"                                                         "<<std::endl;
  std::cout<<"               |  \\ \\ | |/ /                           "<<std::endl;
  std::cout<<"               |  |\\ `' ' /                             "<<std::endl;
  std::cout<<"               |  ;'aorta \\      / , pulmonary          "<<std::endl;
  std::cout<<"               | ;    _,   |    / / ,  arteries          "<<std::endl;
  std::cout<<"      superior | |   (  `-.;_,-' '-' ,                   "<<std::endl;
  std::cout<<"     vena cava | `,   `-._       _,-'_                   "<<std::endl;
  std::cout<<"               |,-`.    `.)    ,<_,-'_, pulmonary                     ______ _____    "<<std::endl;
  std::cout<<"              ,'    `.   /   ,'  `;-' _,  veins                      |  ____|  __ \\   "<<std::endl;
  std::cout<<"             ;        `./   /`,    \\-'                               | |__  | |__) |  "<<std::endl;
  std::cout<<"             | right   /   |  ;\\   |\\                                |  __| |  ___/   "<<std::endl;
  std::cout<<"             | atrium ;_,._|_,  `, ' \\                               | |____| |       "<<std::endl;
  std::cout<<"             |        \\    \\ `       `,                              |______|_|       "<<std::endl;
  std::cout<<"             `      __ `    \\   left  ;,                                              "<<std::endl;
  std::cout<<"              \\   ,'  `      \\,  ventricle                                            "<<std::endl;
  std::cout<<"               \\_(            ;,      ;;                "<<std::endl;
  std::cout<<"               |  \\           `;,     ;;                "<<std::endl;
  std::cout<<"      inferior |  |`.          `;;,   ;'                 "<<std::endl;
  std::cout<<"     vena cava |  |  `-.        ;;;;,;'                  "<<std::endl;
  std::cout<<"               |  |    |`-.._  ,;;;;;'                   "<<std::endl;
  std::cout<<"               |  |    |   | ``';;;'                     "<<std::endl;
  std::cout<<"                       aorta                             "<<std::endl;
  std::cout<<"                                                         "<<std::endl;

}
