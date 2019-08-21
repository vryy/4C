/*---------------------------------------------------------------------*/
/*! \file

\brief student's c++/baci tutorial control algorithm

\maintainer  Martin Kronbichler

\level 2

*/
/*---------------------------------------------------------------------*/

#include "tutorial_dyn.H"
#include "inpar_tutorial.H"

#include "tutorial_nln_truss.H"
#include "tutorial_fixedpoint_scheme.H"
#include "tutorial_struct_ale_coupling.H"

#include "../drt_lib/drt_globalproblem.H"


void tutorial_drt()
{
  /// print welcome message to screen
  printwelcomemessage();

  /// Get pointer to global problem
  /// The class 'Problem' is a global singleton class.
  /// The instance is callable everywhere, as soon as
  ///  he header 'drt_globalproblem.H' is included.
  /// 'Instance()' is a method common to such singleton classes.
  /// It does no more than either creating a new unique 'Problem'-object,
  /// or - if it had already been built - it only returns a pointer to the
  /// global 'Problem'.
  DRT::Problem* problem = DRT::Problem::Instance();

  /// get ParameterList for the tutorial.
  /// this list contains all parameters specified in your .dat- (input-) file.
  const Teuchos::ParameterList& tutorialparams = problem->TutorialParams();

  /// choose algorithm : this must be specified in your .dat- (input-) file next to the parameter
  /// 'TYPE'
  int tutorial_type = DRT::INPUT::IntegralValue<int>(tutorialparams, "TYPE");


  /// detect type of tutorial and build corresponding objects
  switch (tutorial_type)
  {
    case INPAR::TUTORIAL::nonlinear_truss:
    {
      /// initialize rcp pointer to nonlinear truss tutorial object
      Teuchos::RCP<TUTORIAL::NonlinearTruss> algorithm_object = Teuchos::null;

      /// build object containing algorithm for nonlinear truss tutorial
      algorithm_object = Teuchos::rcp(new TUTORIAL::NonlinearTruss());

      /// define geometry and boundary conditions, and print to screen
      algorithm_object->SetupProblem();

      /// execute the time loop
      algorithm_object->TimeLoop();

      /// show the results
      algorithm_object->PrintResults();

      /// clean up the tutorial
      algorithm_object->TutorialDone();

      break;
    }
    case INPAR::TUTORIAL::partitioned_fixed_point_scheme:
    {
      /// initialize rcp pointer to struct-ale coupling tutorial object
      Teuchos::RCP<TUTORIAL::FixedPointScheme> algorithm_object = Teuchos::null;

      /// build object containing algorithm for struct-ale coupling tutorial
      algorithm_object = Teuchos::rcp(new TUTORIAL::FixedPointScheme());

      /// execute the time loop
      algorithm_object->TimeLoop();

      break;
    }
    case INPAR::TUTORIAL::struct_ale_coupling:
    {
      /// initialize rcp pointer to struct-ale coupling tutorial object
      Teuchos::RCP<TUTORIAL::StructAleCoupling> algorithm_object = Teuchos::null;

      /// build object containing algorithm for struct-ale coupling tutorial
      algorithm_object = Teuchos::rcp(new TUTORIAL::StructAleCoupling());

      /// define geometry and boundary conditions, and print to screen
      algorithm_object->SetupProblem();

      /// execute the time loop
      algorithm_object->TimeLoop();

      /// show the results
      algorithm_object->PrintResults();

      /// clean up the tutorial
      algorithm_object->TutorialDone();

      break;
    }
    default:
    {
      dserror("no valid tutorial type specified");
      break;
    }  // default

  }  // switch tutorial_type

  return;
}  // tutorial_drt


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void printwelcomemessage()
{
  std::cout << "\n*****************************************************************" << std::endl;
  std::cout << "*                                                               *" << std::endl;
  std::cout << "*          WELCOME TO THE HANDS-ON C++/BACI TUTORIAL            *" << std::endl;
  std::cout << "*                                                               *" << std::endl;
  std::cout << "*****************************************************************\n" << std::endl;
}
