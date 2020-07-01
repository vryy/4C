/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_DIRECTION_DEFAULTSTEPTEST_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_DIRECTION_DEFAULTSTEPTEST_H_

#include <set>

#include <Teuchos_RCP.hpp>
#include <NOX_Abstract_Group.H>

class Epetra_Vector;
namespace NOX
{
  class Utils;
  namespace Abstract
  {
    class Vector;
  }
  namespace NLN
  {
    class Group;
    namespace Direction
    {
      namespace Test
      {
        /// Generic base class of all direction tests
        class Generic
        {
         public:
          /// constructor
          Generic(const Teuchos::RCP<NOX::Utils>& utils) : utils_(utils){};

          /// destructor
          virtual ~Generic(){};

          /// this routine is called in the first attempt
          virtual bool initAndCheckTest(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp) = 0;

          /// this routine is called in all follow-up checking attempts
          virtual bool checkTest(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp) = 0;

          /** \brief return diagonal vector
           *
           *  The method offers the possibility to return a diagonal vector which
           *  is specifically created based on the information collect throughout
           *  the test. */
          virtual Teuchos::RCP<Epetra_Vector> getCurrentDiagonal(
              const NOX::Abstract::Group& grp) const = 0;

          /// fill diagonal vector
          virtual void fillDiagonal(Epetra_Vector& diagonal) const = 0;

         protected:
          /// output object
          Teuchos::RCP<NOX::Utils> utils_;
        };

        /** \brief Default step test strategy called Volume Change (VC)
         *
         *  \author hiermeier */
        class VolumeChange : public Generic
        {
         public:
          /// constructor
          VolumeChange(const Teuchos::RCP<NOX::Utils>& utils) : Generic(utils){};

          /// destructor
          ~VolumeChange() override{};

         protected:
          /** \brief return true if the internally used test combination indicates acceptance
           *
           *  The test consists of three parts:
           *
           *  If either the step length increases AND bad/invalid elements occur,
           *  OR
           *  if the system matrix is negative/semi definite in search direction,
           *  the current direction will be rejected. Otherwise, if there are bad
           *  elements and the step length decreases or there are now bad elements
           *  and step length increases, while the system matrix stays always
           *  positive definite , the direction will be accepted. This state
           *  of acceptance is always reached after a certain number of
           *  modifications as long all entries on the primal diagonal are modified.
           *
           *  \author hiermeier */
          bool isAccepted();

          /// reset class members and perform the default check
          bool initAndCheckTest(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp) override;

          /** \brief perform the check
           *
           *  1) identify bad elements
           *  2) compute different scalar quality quantities based on the direction vector
           *     and the system matrix.
           *  3) Test for acceptance */
          bool checkTest(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp) override;

          /// return diagonal with the value one at the location of bad element DOFs
          Teuchos::RCP<Epetra_Vector> getCurrentDiagonal(
              const NOX::Abstract::Group& grp) const override;
          Teuchos::RCP<Epetra_Vector> getCurrentDiagonal(const NOX::NLN::Group& grp) const;

          /// fill the diagonal with the value one at the location of bad element DOFs
          inline void fillDiagonal(Epetra_Vector& diagonal) const override
          {
            fillDiagonalAtBadDofs(diagonal);
          }

          /// return true if all volumes are valid under the specified criteria
          bool isValidElementVolumes() const;

          /// return true if step length decreases
          bool isValidDirectionLength() const;

          /** \brief return true if the direction indicates a positive definite system matrix
           *
           *  \note The here used simple test does not tell you that the matrix is truly
           *  positive definite. To say that, the criterion must hold for ALL vectors
           *  beside the zero vector. However, since we are only interested in the current
           *  direction this test seems to be sufficient. */
          bool isPositiveDefinite() const;

          /// compute scalar quality quantities based on the primal search direction
          void computePrimalDirectionMeasures(
              NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp);

          /// compute element volumes
          NOX::Abstract::Group::ReturnType computeElementVolumes(
              NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp);

          /// detected bad elements based on the user-defined criteria
          void identifyBadElements(NOX::Abstract::Group& grp, int& gnew_num_bad_eles);

          /// fill my bad dofs set and return the number of newly bad elements
          int fillMyBadDofs(NOX::NLN::Group& grp);

          /// return an empty diagonal vector
          Teuchos::RCP<Epetra_Vector> getEmptyDiagonal(const NOX::NLN::Group& grp) const;

          /// fill the diagonal vector at all dofs of bad elements
          void fillDiagonalAtBadDofs(Epetra_Vector& diagonal) const;

         private:
          /** reference element volumes (volumes corresponding to previously
           * accepted Newton iterate) */
          Teuchos::RCP<Epetra_Vector> ref_ele_vols_;

          /** current element volumes (volumes corresponding to current
           * trial Newton default step) */
          Teuchos::RCP<Epetra_Vector> curr_ele_vols_;

          /// set containing all bad dofs on this proc
          std::set<int> my_bad_dofs_;

          /// quadratic l2-norm of the default primal step
          double dirdir_ = 100.0;

          /** inner product of primal direction, primal-primal system matrix
           *  block and primal direction */
          double dirres_ = 0.0;

          /// quadratic l2-norm of the default primal step (last accepted Newton iterate)
          double dirdir_last_ = 100.0;

          /// global (over all procs) number of bad elements
          int gnum_bad_eles_ = 0;
        };
      }  // namespace Test
    }    // namespace Direction
  }      // namespace NLN
}  // namespace NOX


#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_DIRECTION_DEFAULTSTEPTEST_H_ */
