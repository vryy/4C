/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN extension of the %NOX::Epetra required
       interface.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_SOLVER_NONLIN_NOX_INTERFACE_REQUIRED_H
#define BACI_SOLVER_NONLIN_NOX_INTERFACE_REQUIRED_H

#include "baci_config.H"

#include "baci_solver_nonlin_nox_enum_lists.H"
#include "baci_utils_exceptions.H"

#include <NOX_Epetra_Interface_Required.H>  // base class
#include <NOX_Epetra_Vector.H>

#include <set>
#include <vector>

namespace NOX
{
  namespace Abstract
  {
    class Group;
  }  // namespace Abstract
}  // namespace NOX

BACI_NAMESPACE_OPEN

namespace NOX
{
  namespace NLN
  {
    namespace Interface
    {
      class Required : public virtual ::NOX::Epetra::Interface::Required
      {
       public:
        //! Constructor
        Required(){};

        //! returns the right-hand-side norms of the primary DoF fields
        virtual double GetPrimaryRHSNorms(const Epetra_Vector& F,
            const NOX::NLN::StatusTest::QuantityType& checkQuantity,
            const ::NOX::Abstract::Vector::NormType& type = ::NOX::Abstract::Vector::TwoNorm,
            const bool& isScaled = false) const = 0;

        //! Returns the Root Mean Squares (abbr.: RMS) of the primary solution updates
        virtual double GetPrimarySolutionUpdateRMS(const Epetra_Vector& xNew,
            const Epetra_Vector& xOld, const double& aTol, const double& rTol,
            const NOX::NLN::StatusTest::QuantityType& checkQuantity,
            const bool& disable_implicit_weighting = false) const = 0;

        //! Returns the increment norm of the primary DoF fields
        virtual double GetPrimarySolutionUpdateNorms(const Epetra_Vector& xNew,
            const Epetra_Vector& xOld, const NOX::NLN::StatusTest::QuantityType& checkQuantity,
            const ::NOX::Abstract::Vector::NormType& type = ::NOX::Abstract::Vector::TwoNorm,
            const bool& isScaled = false) const = 0;

        //! Returns the previous solution norm of primary DoF fields
        virtual double GetPreviousPrimarySolutionNorms(const Epetra_Vector& xOld,
            const NOX::NLN::StatusTest::QuantityType& checkQuantity,
            const ::NOX::Abstract::Vector::NormType& type = ::NOX::Abstract::Vector::TwoNorm,
            const bool& isScaled = false) const = 0;

        //! compute and return some energy representative
        virtual double GetModelValue(const Epetra_Vector& x, const Epetra_Vector& F,
            const enum MeritFunction::MeritFctName merit_func_type) const = 0;

        //! return model terms of a linear model (optional)
        virtual double GetLinearizedModelTerms(const ::NOX::Abstract::Group* group,
            const Epetra_Vector& dir, const enum NOX::NLN::MeritFunction::MeritFctName mf_type,
            const enum NOX::NLN::MeritFunction::LinOrder linorder,
            const enum NOX::NLN::MeritFunction::LinType lintype) const
        {
          dserror("Not implemented!");
          exit(EXIT_FAILURE);
        }

        //! calculate characteristic/reference norms for forces
        virtual double CalcRefNormForce() = 0;

        //! access the lumped mass matrix
        virtual Teuchos::RCP<const Epetra_Vector> GetLumpedMassMatrixPtr() const
        {
          dserror("The evaluation of the lumped mass matrix is not implemented!");
          return Teuchos::null;
        }

        //! create a backup state (optional)
        virtual void CreateBackupState(const Epetra_Vector& dir)
        {
          dserror("There is no meaningful implementation for this method!");
        }

        //! recover from a backup state (optional)
        virtual void RecoverFromBackupState()
        {
          dserror("There is no meaningful implementation for this method!");
        }

        //! compute element volumes (optional)
        virtual bool computeElementVolumes(
            const Epetra_Vector& x, Teuchos::RCP<Epetra_Vector>& ele_vols) const
        {
          dserror("There is no meaningful implementation for this method!");
          exit(EXIT_FAILURE);
        }

        //! access dofs of specific elements (optional)
        virtual void getDofsFromElements(
            const std::vector<int>& my_ele_gids, std::set<int>& my_ele_dofs) const
        {
          dserror("There is no meaningful implementation for this method!");
        };
      };
    }  // namespace Interface
  }    // namespace NLN
}  // namespace NOX


BACI_NAMESPACE_CLOSE

#endif  // SOLVER_NONLIN_NOX_INTERFACE_REQUIRED_H
