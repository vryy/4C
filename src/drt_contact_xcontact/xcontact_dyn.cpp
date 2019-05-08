/*----------------------------------------------------------------------------*/
/**
\file xcontact_dyn.cpp

\brief inequality level-set approach a.k.a. xcontact or extended contact

\maintainer Matthias Mayr

\date Jun 14, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "../drt_contact_xcontact/xcontact_dyn.H"

#include "../drt_contact_xcontact/xcontact_algorithm_partitioned.H"
#include "../drt_inpar/inpar_contact_xcontact.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret_xfem.H"
#include "../drt_lib/drt_dofset_fixed_size.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void xcontact_dyn(const int& restart)
{
  DRT::Problem* problem_ptr = DRT::Problem::Instance();

  // create a communicator
  const Epetra_Comm& comm = problem_ptr->GetDis("structure")->Comm();

  if (comm.MyPID() == 0)
  {
    std::cout << "                    XXX     XXX                    \n";
    std::cout << "                     XXX   XXX                     \n";
    std::cout << "                      XXX XXX                      \n";
    std::cout << "                       XXXXX                       \n";
    std::cout << "                       XXXXX                       \n";
    std::cout << "                      XXX XXX                      \n";
    std::cout << "                     XXX   XXX                     \n";
    std::cout << "                    XXX     XXX                    \n";
    std::cout << "                                                   \n";
    std::cout << "       +++++ eXtended Contact Dynamics +++++       \n";
    std::cout << "       +++ INEQUALITY LEVEL-SET APPROACH +++       \n";
    std::cout << "                                                   \n";
  }

  // access parameters for eXtended contact dynamics
  const Teuchos::ParameterList& p_xcontact_dyn = problem_ptr->XContactDynamicParams();
  const Teuchos::ParameterList& p_structure_dyn = problem_ptr->StructuralDynamicParams();
  const Teuchos::ParameterList& p_scatra_dyn = problem_ptr->ScalarTransportDynamicParams();
  const Teuchos::ParameterList& p_xfem_general = problem_ptr->XFEMGeneralParams();

  // ---------------------------------------------------------------------------
  // Setup the X-FEM discretization
  // ---------------------------------------------------------------------------
  // access the structure discretization
  Teuchos::RCP<DRT::Discretization> struct_dis_ptr = problem_ptr->GetDis("structure");
  if (not struct_dis_ptr->Filled()) struct_dis_ptr->FillComplete();

  // ---------------------------------------------------------------------------
  // Build the coupling algorithm
  // ---------------------------------------------------------------------------
  // create an empty XCONTACT::Algorithm instance
  Teuchos::RCP<XCONTACT::ALGORITHM::Base> xcontact_alg_ptr = Teuchos::null;

  // coupling between scatra and structure field
  const INPAR::XCONTACT::FieldCouplingAlgorithm coupling =
      DRT::INPUT::IntegralValue<INPAR::XCONTACT::FieldCouplingAlgorithm>(
          p_xcontact_dyn, "COUPL_ALGO");

  // build the desired coupling algorithm
  switch (coupling)
  {
    case INPAR::XCONTACT::field_coupl_algo_partitioned:
    {
      xcontact_alg_ptr = Teuchos::rcp(new XCONTACT::ALGORITHM::Partitioned());
      break;
    }
    case INPAR::XCONTACT::field_coupl_algo_monolithic:
    {
      dserror(
          "The monolithic coupling scheme is currently unsupported "
          "for the eXtended contact approach!");
      break;
    }
    case INPAR::XCONTACT::field_coupl_algo_undefined:
    default:
    {
      dserror(
          "The coupling scheme for the eXtended contact approach "
          "is not defined. Fix your input file! (COUPL_ALGO = %d | %s)",
          coupling, INPAR::XCONTACT::FieldCouplingAlgorithm2String(coupling).c_str());
      break;
    }
  }

  // initialize and setup the chosen algorithm
  xcontact_alg_ptr->Init(
      p_xcontact_dyn, p_structure_dyn, p_scatra_dyn, p_xfem_general, struct_dis_ptr);
  xcontact_alg_ptr->Setup();

  if (restart) xcontact_alg_ptr->ReadRestart(restart);

  // start simulating
  xcontact_alg_ptr->Timeloop();

  return;
}
