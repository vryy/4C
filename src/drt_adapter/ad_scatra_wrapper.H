/*----------------------------------------------------------------------*/
/*!
\file ad_scatra_wrapper.H
\brief Wrapper for the scatra time integrator.
\level 1
\maintainer Andreas Rauch
 */
/*----------------------------------------------------------------------*/


#ifndef AD_SCATRA_WRAPPER_
#define AD_SCATRA_WRAPPER_


#include "ad_scatra_interface.H"

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>


// forward declarations

namespace LINALG {
  class SparseOperator;
}

/*! \brief Just wrap the time integrator. Do nothing new.
 *
 * \date 12/2016
 * \author rauch (rauch@lnm.mw.tum.de)    */
namespace ADAPTER {

class AdapterScatraWrapper : public ScatraInterface
{

public:

  /// constructor
  explicit AdapterScatraWrapper(Teuchos::RCP<ScatraInterface> scatra);

  /// virtual to get polymorph destruction
  virtual ~AdapterScatraWrapper() {};

  /// compute contribution of mechanical state to eq. system
  virtual void EvaluateAdditionalSolutionDependingModels(
      Teuchos::RCP<LINALG::SparseOperator> systemmatrix,      //!< system matrix
      Teuchos::RCP<Epetra_Vector>          rhs                //!< rhs vector
  );


protected:

  //! return discretization
  virtual Teuchos::RCP<DRT::Discretization> Discretization() const
  {return scatra_timint_ -> Discretization(); };

  //! add parameters specific for time-integration scheme
  virtual void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver=false)
  {scatra_timint_ -> AddTimeIntegrationSpecificVectors(forcedincrementalsolver); };

  //! return number of dofset associated with displacement dofs
  virtual int NdsDisp()
  {return scatra_timint_ -> NdsDisp(); };

  /// return rcp ptr to neumann loads vector
  virtual Teuchos::RCP<Epetra_Vector> GetNeumannLoadsPtr()
  {return scatra_timint_->GetNeumannLoadsPtr();};

private:

  Teuchos::RCP<ScatraInterface> scatra_timint_; ///< underlying structural time integration

}; // class AdapterScatraWrapper
} // namespace ADAPTER


#endif /* AD_SCATRA_WRAPPER_ */
