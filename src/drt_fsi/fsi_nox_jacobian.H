/*----------------------------------------------------------------------*/
/*!

\brief Matrix-free Newton-Krylov for FSI

\level 1

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/

#ifndef FSI_NOX_JACOBIAN_H
#define FSI_NOX_JACOBIAN_H

#include <Epetra_Comm.h>
#include <Epetra_Operator.h>
#include <Epetra_Import.h>

#include <NOX_Abstract_Group.H>
#include <NOX_Common.H>
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Utils.H>

#include <Teuchos_RCP.hpp>

// Forward Declarations
class Epetra_Comm;
class Epetra_Map;

namespace NOX
{
  namespace FSI
  {
    /// Matrix Free Newton Krylov based on an approximation of the residuum derivatives
    class FSIMatrixFree : public Epetra_Operator, public virtual NOX::Epetra::Interface::Jacobian
    {
     public:
      /*! \brief Constructor

        The vector \c x is used to clone the solution vector.
      */
      FSIMatrixFree(Teuchos::ParameterList& printParams,
          const Teuchos::RCP<NOX::Epetra::Interface::Required>& i,
          const NOX::Epetra::Vector& cloneVector);

      //! Pure virtual destructor
      virtual ~FSIMatrixFree();

      //! If set true, transpose of this operator will be applied.
      /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this
        flag affects only the Apply() and ApplyInverse() methods.  If the implementation of this
        interface does not support transpose use, this method should return a value of -1. \param
        UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

        \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does
        not support transpose.
      */
      virtual int SetUseTranspose(bool UseTranspose);

      //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
      /*!
        \param     X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
        \param     Y - A Epetra_MultiVector of dimension NumVectors containing result.

        \return Integer error code, set to 0 if successful.
      */
      virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

      //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
      /*!
        \param     X - A Epetra_MultiVector of dimension NumVectors to solve for.
        \param     Y -A Epetra_MultiVector of dimension NumVectors containing result.

        \return Integer error code, set to 0 if successful.

        \warning In order to work with AztecOO, any implementation of this method must
        support the case where X and Y are the same object.
      */
      virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

      //! Returns the infinity norm of the global matrix.
      /* Returns the quantity \f$ \| A \|_\infty\f$ such that
         \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].

         \warning This method must not be called unless HasNormInf() returns true.    */
      virtual double NormInf() const;

      //! Returns a character string describing the operator
      virtual const char* Label() const;

      //! Returns the current UseTranspose setting.
      virtual bool UseTranspose() const;

      //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
      virtual bool HasNormInf() const;

      //! Returns a reference to the Epetra_Comm communicator associated with this operator.
      virtual const Epetra_Comm& Comm() const;

      //! Returns the Epetra_BlockMap object associated with the domain of this matrix operator.
      virtual const Epetra_Map& OperatorDomainMap() const;

      //! Returns the Epetra_BlockMap object associated with the range of this matrix operator.
      virtual const Epetra_Map& OperatorRangeMap() const;

      //! Compute Jacobian given the specified input vector, x.  Returns true if computation was
      //! successful.
      virtual bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

      //! Clone a NOX::Abstract::Group derived object and use the computeF() method of that group
      //! for the perturbation instead of the NOX::Epetra::Interface::Required::computeF() method.
      //! This is required for LOCA to get the operators correct during homotopy.
      void setGroupForComputeF(const NOX::Abstract::Group& group);

     protected:
      //! Label for matrix
      std::string label;

      //! User provided interface function
      Teuchos::RCP<NOX::Epetra::Interface::Required> interface;

      //! The current solution vector
      NOX::Epetra::Vector currentX;

      //! Perturbed solution vector
      mutable NOX::Epetra::Vector perturbX;

      //! Perturbed solution vector
      mutable NOX::Epetra::Vector perturbY;

      //! Epetra_Map object used in the returns of the Epetra_Operator derived methods.
      /*! If the user is using Epetra_BlockMaps, then NOX::Epetra::MatrixFree must create an
       * equivalent Epetra_Map from the Epetra_BlockMap that can be used as the return object of the
       * OperatorDomainMap() and OperatorRangeMap() methods.
       */
      Teuchos::RCP<const Epetra_Map> epetraMap;

      //! Flag to enables the use of a group instead of the interface for the computeF() calls in
      //! the directional difference calculation.
      bool useGroupForComputeF;

      //! Pointer to the group for possible use in computeF() calls.
      Teuchos::RCP<NOX::Abstract::Group> groupPtr;

      //! Printing utilities.
      NOX::Utils utils;
    };

  }  // namespace FSI
}  // namespace NOX

#endif
