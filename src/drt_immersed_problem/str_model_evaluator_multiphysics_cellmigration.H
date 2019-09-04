/*----------------------------------------------------------------------*/
/*! \file

\brief Generic class for all model evaluators.

\level 3

\maintainer Jonas Eichinger

*----------------------------------------------------------------------*/


#ifndef STR_MODEL_EVALUATOR_MULTIPHYSICS_CELLMIGRATION_H_
#define STR_MODEL_EVALUATOR_MULTIPHYSICS_CELLMIGRATION_H_


#include "../drt_structure_new/str_model_evaluator_multiphysics.H"

// forward declaration
class Epetra_Vector;
class Epetra_Map;

namespace LINALG
{
  class SparseOperator;
}  // namespace LINALG

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

namespace IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace IO

namespace DRT
{
  class Discretization;
}  // namespace DRT

namespace STR
{
  class Integrator;
  namespace TIMINT
  {
    class BaseDataGlobalState;
    class BaseDataIO;
    class Base;
  }  // namespace TIMINT
  namespace MODELEVALUATOR
  {
    class Data;
    /*! \brief This is the model evaluator for cell migration.
     *
     *  This object wraps the model evaluators needed for cell migration.
     *
     *  \date 11/16
     *  \author rauch */
    class CellMigration : public Multiphysics
    {
     public:
      //! constructor
      CellMigration();

      //! destructor
      virtual ~CellMigration(){};

      //! initialize the class variables
      virtual void Init(const Teuchos::RCP<STR::MODELEVALUATOR::Data>& eval_data_ptr,
          const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
          const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
          const Teuchos::RCP<STR::Integrator>& int_ptr,
          const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr, const int& dof_offset);

      //! setup class variables
      virtual void Setup();

    };  // class CellMigration

  }  // namespace MODELEVALUATOR
}  // namespace STR

#endif /* STR_MODEL_EVALUATOR_MULTIPHYSICS_CELLMIGRATION_H_ */
