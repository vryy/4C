/*!----------------------------------------------------------------------
\file pasi_str_model_evaluator_partitioned.H

\brief model evaluator for structure part of partitioned pasi

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                           sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
#ifndef PASI_STR_MODEL_EVALUATOR_PARTITIONED_H_
#define PASI_STR_MODEL_EVALUATOR_PARTITIONED_H_

/*----------------------------------------------------------------------*
 | headers                                               sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
#include "../drt_structure_new/str_model_evaluator_generic.H"

/*----------------------------------------------------------------------*
 | forward declarations                                  sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
namespace ADAPTER
{
class Structure;
} // namespace ADAPTER

/*----------------------------------------------------------------------*
 | partitioned pasi structure model evaluator            sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
namespace STR {

namespace MODELEVALUATOR {

class PartitionedPASI: public Generic
{

public:

  //! constructor
  PartitionedPASI();

  //! destructor
  virtual ~PartitionedPASI(){};

  //! setup class variables [derived]
  virtual void Setup();

  //! @name Functions which are derived from the base generic class
  //! @{
  //! [derived]
  INPAR::STR::ModelType Type() const
  { return INPAR::STR::model_partitioned_coupling; }

  //! reset class variables (without jacobian) [derived]
  virtual void Reset(const Epetra_Vector& x)
  {return;};

  //! [derived]
  virtual bool EvaluateForce()
  {return true;};

  //! [derived]
  virtual bool EvaluateStiff()
  {return true;};

  //! [derived] not needed in partitioned scheme
  virtual bool EvaluateForceStiff()
  {return true;};

  //! derived
  virtual bool AssembleForce(Epetra_Vector& f,
      const double & timefac_np) const;

  //! Assemble the jacobian at \f$t_{n+1}\f$ not needed in partitioned scheme
  virtual bool AssembleJacobian(
      LINALG::SparseOperator& jac,
      const double & timefac_np) const
  {return true;};

  //! [derived]
  virtual void WriteRestart(
      IO::DiscretizationWriter& iowriter,
      const bool& forced_writerestart) const
  {return;};

  //! [derived]
  virtual void ReadRestart(IO::DiscretizationReader& ioreader)
  {return;};

  //! recover condensed Lagrange multipliers
  virtual void RecoverState(
      const Epetra_Vector& xold,
      const Epetra_Vector& dir,
      const Epetra_Vector& xnew)
  {return;};

  //! [derived]
  virtual void UpdateStepState(const double& timefac_n);

  //! [derived]
  virtual void UpdateStepElement()
  {return;};

  //! [derived]
  virtual void DetermineStressStrain()
  {return;};

  //! [derived]
  virtual void DetermineEnergy()
  {return;};

  //! [derived]
  virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const
  {return;};

  //! derived
  virtual void ResetStepState()
  {return;};

  //! [derived]
  virtual void PostOutput()
  {return;};

  //! [derived]
  virtual Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const;

  //! [derived]
  virtual Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const;

  //! [derived]
  virtual Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const;
  //! @}

}; // class PartitionedPASI: public Generic

} // namespace MODELEVALUATOR

} // namespace STR

/*----------------------------------------------------------------------*/
#endif // PASI_STR_MODEL_EVALUATOR_PARTITIONED_H_
