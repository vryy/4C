/*----------------------------------------------------------------------*/
/*! \file
  \brief control routine of fluid-porous-structure-interaction problems
 \level 2
 *------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  headers                                                             |
 *----------------------------------------------------------------------*/
#include "4C_fpsi_dyn.hpp"

#include "4C_fpsi.hpp"
#include "4C_fpsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fpsi.hpp"
#include "4C_lib_discret.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*------------------------------------------------------------------------------------------------*
 | main control routine for fluid-porous-structure-interaction problems                rauch 11/12 |
 *------------------------------------------------------------------------------------------------*/
void fpsi_drt()
{
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  // 1.- Get Communicator
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  // print the chuck
  if (comm.MyPID() == 0)
  {
    std::cout << "                            ``;@@@@@@@@@#@`               " << std::endl;
    std::cout << "                            .@#@#@@@@@@@@@@,+.`           " << std::endl;
    std::cout << "                          '##@@@@@@@@@@@@@@@+;`           " << std::endl;
    std::cout << "                         :@#@@@@@@@@@@@@@@@@'`#           " << std::endl;
    std::cout << "                         @@@#@@@@@@@@@@@@@@@#@#:          " << std::endl;
    std::cout << "                        #@@@@@@@@@@@@@@@@@@@@@@@          " << std::endl;
    std::cout << "                       ,##@@@@@@@@@@@@@@@@@@@@@@:         " << std::endl;
    std::cout << "                       ;###@@@@@@@@#@@@@@@@@@@@@#         " << std::endl;
    std::cout << "                        @@@@@@@@@@@@@#@#@@@@@@@@@         " << std::endl;
    std::cout << "                       `@#@@@@@@@@@@@@@@@@@@@@@@@         " << std::endl;
    std::cout << "                       ,@#@@@@@@@@@#':#@@##@@@@@@`        " << std::endl;
    std::cout << "                       +##@@@@@@@@@+: '#@@#@@@@@@@        " << std::endl;
    std::cout << "                       :@#@@@@@@@@@;    .`#@@@@@@@        " << std::endl;
    std::cout << "                       ,@@@@@@@@##;.     `@@#@@@##.       " << std::endl;
    std::cout << "                       .@##@@;.@@@        @@#@@@@#+       " << std::endl;
    std::cout << "                        @@@@#@#@#`         @##@@#@#       " << std::endl;
    std::cout << "                        @@@@@  `           #@#@@@##       " << std::endl;
    std::cout << "                        ,#@@````           +@@#@#@#       " << std::endl;
    std::cout << "                         @,@.   `..        +@#@,@##       " << std::endl;
    std::cout << "                         #@@,@# #@@@@#      @@ # @#       " << std::endl;
    std::cout << "                         '##@@@ ;@@@#@'     ` + :##       " << std::endl;
    std::cout << "                         `@@@@@  @ ``      ';   #@`       " << std::endl;
    std::cout << "                          @@@##            .@   @@        " << std::endl;
    std::cout << "                          #  @+            @:  :##        " << std::endl;
    std::cout << "                          '  @.            @,  @@,        " << std::endl;
    std::cout << "                          #  @            `@` '@@         " << std::endl;
    std::cout << "                          , `@            .#: ;#.         " << std::endl;
    std::cout << "                           +'@     #       #   ,          " << std::endl;
    std::cout << "                           @@@   +       :.`   .          " << std::endl;
    std::cout << "                           @##@@        ,@'@   @;         " << std::endl;
    std::cout << "                           @@@@@@#': #. @#@+   ;@         " << std::endl;
    std::cout << "                     @;+#@@@@####@@@@@@@@@@`   +#'        " << std::endl;
    std::cout << "                    @@#@@@#@@@@+;   ;@@#:#'    #@@        " << std::endl;
    std::cout << "                    ;#@@@#@@@@@'.   :#@:@@    ;#@@        " << std::endl;
    std::cout << "                     @#@@@@@@@@@@@@@@###@@   ,@@@@;       " << std::endl;
    std::cout << "                      @#@@@@@@@@#@@#@@@@@+   @@@#@#+      " << std::endl;
    std::cout << "                     :.#@@@@@@@#@#@##@###  `@#@@@@@@#     " << std::endl;
    std::cout << "                      _______ ______ _______  ______      " << std::endl;
    std::cout << "                        ||____ ||___|||_____    ||        " << std::endl;
    std::cout << "                        ||     ||     _____|| __||__      " << std::endl;
    std::cout << std::endl << std::endl;
  }

  // 2.- Parameter reading
  const Teuchos::ParameterList& fpsidynparams = problem->FPSIDynamicParams();
  const Teuchos::ParameterList& poroelastdynparams = problem->poroelast_dynamic_params();

  Teuchos::RCP<FPSI::Utils> FPSI_UTILS = FPSI::Utils::Instance();

  // 3.- Creation of Poroelastic + Fluid problem. (Discretization called inside)
  Teuchos::RCP<FPSI::FpsiBase> fpsi = Teuchos::null;
  fpsi = FPSI_UTILS->setup_discretizations(comm, fpsidynparams, poroelastdynparams);

  // 3.1- Read restart if needed.
  const int restartstep = problem->Restart();
  if (restartstep)
  {
    fpsi->read_restart(restartstep);
  }

  // 3.2.- redistribute the FPSI interface
  fpsi->redistribute_interface();

  //////////////////////////////////
  // 4.- Run of the actual problem.//
  //////////////////////////////////

  // 4.1.- Coupling and creation of combined dofmap
  fpsi->SetupSystem();
  // setup the linear solver, if necessary
  fpsi->SetupSolver();
  // 4.2.- Solve the whole problem
  fpsi->Timeloop();
  Teuchos::TimeMonitor::summarize();

  // 5. - perform the result test
  fpsi->TestResults(comm);


  return;
}  // fpsi_drt()

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
