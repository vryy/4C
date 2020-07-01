/*----------------------------------------------------------------------------*/
/** \file
\brief partitioned algorithm for the inequality level-set approach for contact
       problems a.k.a. xcontact or extended contact


\level 3
*/
/*----------------------------------------------------------------------------*/

#include "xcontact_algorithm_partitioned.H"
#include "xcontact_multi_discretization_wrapper.H"
#include "xcontact_levelset_algorithm.H"

#include "../drt_adapter/ad_str_xcontact.H"

#include "../drt_io/io_pstream.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XCONTACT::ALGORITHM::Partitioned::Partitioned()
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XCONTACT::ALGORITHM::Partitioned::Setup()
{
  // call the Setup() routine of the base algorithm first
  XCONTACT::ALGORITHM::Base::Setup();

  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XCONTACT::ALGORITHM::Partitioned::OuterLoop()
{
  int itnum = 0;
  bool isconverged = false;

  if (Comm().MyPID() == 0)
  {
    std::cout << "*=======================================================*\n";
    std::cout << "||     PARTITIONED OUTER ITERATION LOOP                ||\n";
    std::cout << "*=======================================================*\n";
    std::cout << std::flush;

    printf("TIME: %11.4E/%11.4E  DT = %11.4E  %s  STEP = %4d/%4d\n", Time(), MaxTime(), Dt(),
        ScaTraField().MethodTitle().c_str(), Step(), NumStep());
  }

  /* initially solve the structural problem and detect a possible overlap,
   * the current boundary conditions were applied in the PrepareTimeStep()
   * call. */
  DoStructureField();

  // Prepare variables for convergence check.
  PrepareOuterIteration();

  while (not isconverged)
  {
    ++itnum;

    // solve the level-set problem
    DoScaTraField();

    // solve the structural problem
    DoStructureField();

    // check convergence and stop iteration loop if convergence is achieved
    isconverged = ConvergenceCheck(itnum);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XCONTACT::ALGORITHM::Partitioned::PrepareOuterIteration()
{
  if (not StructureField().IsComingIntoContact()) return;

  // access the weighted gap vector with slave normal dof row map layout
  const Epetra_Vector& wgap = StructureField().GetWeightedGap();

  /* ToDo If more interfaces are involved, we need a map extractor to
   * extract the wgap vector of each contact interface and transfer
   * these partial wgap vectors to the corresponding ScaTra
   * discretizations. */
  MultiDiscret().Contact2ScaTra(wgap, *ScaTraField().Phinp(), XFEM::xstructure, true);

  ScaTraField().Reinitialization();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool XCONTACT::ALGORITHM::Partitioned::ConvergenceCheck(const int& itnum)
{
  IO::cout << __LINE__ << " -- " << __PRETTY_FUNCTION__ << ": Not yet implemented!" << IO::endl;
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XCONTACT::ALGORITHM::Partitioned::DoScaTraField()
{
  // Set relevant structure values in ScaTra field.
  if (not SetStructureValuesInScaTra()) return;

  if (Comm().MyPID() == 0)
  {
    std::cout << "*-------------------------------------------------------*\n";
    std::cout << "|              LEVEL-SET SOLVER                         |\n";
    std::cout << "*-------------------------------------------------------*\n";
    std::cout << std::flush;
  }

  // Solve the ScaTra field.
  ScaTraField().Solve();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XCONTACT::ALGORITHM::Partitioned::DoStructureField()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "*-------------------------------------------------------*\n";
    std::cout << "|        STRUCTURE / XCONTACT SOLVER                    |\n";
    std::cout << "*-------------------------------------------------------*\n";
    std::cout << std::flush;
  }

  // Set relevant ScaTra values in structure field.
  SetScaTraValuesInStructure();

  // Solve the structure field.
  StructureField().Solve();

  return;
}
