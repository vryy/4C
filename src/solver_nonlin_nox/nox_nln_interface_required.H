/*-----------------------------------------------------------*/
/*!
\file nox_nln_interface_required.H

\brief %NOX::NLN extension of the %NOX::Epetra required
       interface.

\maintainer Michael Hiermeier

\date Jun 29, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_INTERFACE_REQUIRED_H_
#define NOX_NLN_INTERFACE_REQUIRED_H_

#include "nox_nln_enum_lists.H"
#include "../drt_lib/drt_dserror.H"

#include <NOX_Epetra_Interface_Required.H>    // base class
#include <NOX_Epetra_Vector.H>

namespace NOX {
namespace Abstract {
  class Group;
} // namespace Abstract
namespace NLN {
namespace Interface {
  class Required : public virtual NOX::Epetra::Interface::Required
  {
  public:
    //! Constructor
    Required() {};

    //! Destructor
    virtual ~Required() {};

    //! returns the right-hand-side norms of the primary DoF fields
    virtual double GetPrimaryRHSNorms(
        const Epetra_Vector& F,
        const NOX::NLN::StatusTest::QuantityType& checkQuantity,
        const NOX::Abstract::Vector::NormType& type = NOX::Abstract::Vector::TwoNorm,
        const bool& isScaled = false) const = 0;

    //! Returns the Root Mean Squares (abbr.: RMS) of the primary solution updates
    virtual double GetPrimarySolutionUpdateRMS(
        const Epetra_Vector& xNew, const Epetra_Vector& xOld,
        const double& aTol, const double& rTol,
        const NOX::NLN::StatusTest::QuantityType& checkQuantity,
        const bool& disable_implicit_weighting = false) const = 0;

    //! Returns the increment norm of the primary DoF fields
    virtual double GetPrimarySolutionUpdateNorms(
        const Epetra_Vector& xNew, const Epetra_Vector& xOld,
        const NOX::NLN::StatusTest::QuantityType& checkQuantity,
        const NOX::Abstract::Vector::NormType& type = NOX::Abstract::Vector::TwoNorm,
        const bool& isScaled = false) const = 0;

    //! Returns the previous solution norm of primary DoF fields
    virtual double GetPreviousPrimarySolutionNorms(
        const Epetra_Vector& xOld,
        const NOX::NLN::StatusTest::QuantityType& checkQuantity,
        const NOX::Abstract::Vector::NormType& type = NOX::Abstract::Vector::TwoNorm,
        const bool& isScaled = false) const = 0;

    //! compute and return some energy representative
    virtual double GetModelValue(
        const Epetra_Vector& x,
        const Epetra_Vector& F,
        const enum MeritFunction::MeritFctName merit_func_type ) const = 0;

    virtual double GetLinearizedModelTerms(
        const NOX::Abstract::Group* group,
        const Epetra_Vector& dir,
        const enum NOX::NLN::MeritFunction::MeritFctName mf_type,
        const enum NOX::NLN::MeritFunction::LinOrder linorder,
        const enum NOX::NLN::MeritFunction::LinType lintype ) const
    {
      dserror( "Not implemented!" );
      exit( EXIT_FAILURE );
    }

    //! calculate characteristic/reference norms for forces
    virtual double CalcRefNormForce() = 0;

    virtual Teuchos::RCP<const Epetra_Vector> GetLumpedMassMatrixPtr() const
    {
      dserror("The evaluation of the lumped mass matrix is not implemented!");
      return Teuchos::null;
    }

    virtual void CreateBackupState( const Epetra_Vector& dir )
    {
      dserror( "There is no meaningful implementation for this method!" );
    }

    virtual void RecoverFromBackupState()
    {
      dserror( "There is no meaningful implementation for this method!" );
    }
  };
} // namespace Interface
} // namespace STR
} // namespace NOX


#endif /* NOX_NLN_INTERFACE_REQUIRED_H_ */
