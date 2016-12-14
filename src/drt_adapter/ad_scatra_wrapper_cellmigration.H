/*----------------------------------------------------------------------*/
/*!
\file ad_scatra_wrapper_cellmigration.H
\brief Cell Migration specific wrapper for the scatra time integrator.
\level 1
\maintainer Andreas Rauch
*/
/*----------------------------------------------------------------------*/


#ifndef AD_SCATRA_WRAPPER_CELLMIGRATION_
#define AD_SCATRA_WRAPPER_CELLMIGRATION_


#include "ad_scatra_wrapper.H"

// forward declarations
namespace LINALG {
  class SparseOperator;
}


namespace ADAPTER {

/*! \brief Specialized adapter for cell migration simulations.
 *
 * \date 12/2016
 * \author rauch (rauch@lnm.mw.tum.de)    */
class AdapterScatraWrapperCellMigration : public AdapterScatraWrapper
{

public:

  /// constructor
  explicit AdapterScatraWrapperCellMigration(Teuchos::RCP<ScatraInterface> scatra);

  /// virtual to get polymorph destruction
  virtual ~AdapterScatraWrapperCellMigration() {};

  /// compute contribution of mechanical state to eq. system
  virtual void EvaluateAdditionalSolutionDependingModels(
      Teuchos::RCP<LINALG::SparseOperator> systemmatrix,      //!< system matrix
      Teuchos::RCP<Epetra_Vector>          rhs     //!< rhs vector
  );

  /// add vector to rhs
  void AddContributionToRHS(
      const Teuchos::RCP<const Epetra_Vector>& contributing_vector);

private:

  Teuchos::RCP<Epetra_Vector> rates_;

}; // class AdapterScatraWrapperCellMigration
} // namespace ADAPTER


#endif /* AD_SCATRA_WRAPPER_CELLMIGRATION_ */
