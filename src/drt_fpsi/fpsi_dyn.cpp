/*!------------------------------------------------------------------------------------------------*
 \file fpsi.cpp

 \brief control routine of fluid-porous-structure-interaction problems

 <pre>
   Maintainer: Andreas Rauch
               rauch@lnm.mw.tum.de
 </pre>
 *------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  headers                                                             |
 *----------------------------------------------------------------------*/
#include "fpsi.H"
#include "fpsi_dyn.H"
#include "fpsi_utils.H"
#include "fpsi_monolithic.H"
#include "fpsi_monolithic_plain.H"
#include "fpsi_partitioned.H"
#include "../drt_inpar/inpar_fpsi.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include <Teuchos_TimeMonitor.hpp>

/*------------------------------------------------------------------------------------------------*
 | main control routine for fluid-porous-structure-interaction problems                rauch 11/12 |
 *------------------------------------------------------------------------------------------------*/
void fpsi_drt()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  //1.- Get Communicator
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
    std::cout <<  std::endl << std::endl;
  }

  //2.- Parameter reading
  const Teuchos::ParameterList& fpsidynparams       = problem->FPSIDynamicParams();
  const Teuchos::ParameterList& poroelastdynparams  = problem->PoroelastDynamicParams();

  Teuchos::RCP<FPSI::UTILS> FPSI_UTILS = FPSI::UTILS::Instance();

  //3.- Creation of Poroelastic + Fluid problem. (Discretization called inside)
  Teuchos::RCP<FPSI::FPSI_Base> fpsi_ = Teuchos::null;
  fpsi_ = FPSI_UTILS->SetupDiscretizations(comm, fpsidynparams,poroelastdynparams);

  //3.1- Read restart if needed.
  const int restartstep = DRT::Problem::Instance()->Restart();
  if (restartstep)
    {
      fpsi_->ReadRestart(restartstep);
    }

  //////////////////////////////////
 //4.- Run of the actual problem.//
//////////////////////////////////

  // 4.1.- Coupling and creation of combined dofmap
  fpsi_->SetupSystem();
  // 4.2.- Solve the whole problem
  fpsi_->Timeloop();
  Teuchos::TimeMonitor::summarize();

  // 5. - perform the result test
  fpsi_->TestResults(comm);


  return;
}//fpsi_drt()

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
