/*!----------------------------------------------------------------------
\file constraintpenalty.H

\brief Basic constraint class, dealing with constraints living on boundaries

\maintainer Marc Hirschvogel

\level 2

*----------------------------------------------------------------------*/

#ifndef CONSTRAINTPENALTY_H
#define CONSTRAINTPENALTY_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Teuchos_ParameterList.hpp>

#include "constraint.H"



/*!
\brief UTILS: namespace containing general utils for BACI, as constraints, potentials, surface
stress, ...
*/

namespace UTILS
{
  /*!
  \brief Basic constraint class, dealing with constraint and monitor boundary conditions.
  This class cannot handle multi point constraints, they will be dealt with by a derived class.
  */
  class ConstraintPenalty : public Constraint

  {
   public:
    /*!
    \brief Constructor of a constraint based on a conditions with a given name.
    */

    ConstraintPenalty(
        Teuchos::RCP<DRT::Discretization> discr,  ///< Discretization constraint lives on
        const std::string& CondName               ///< Name of condition to create constraints from
    );

    /*!
        \brief Destructor

     */
    ~ConstraintPenalty() { return; };

    /// unused
    void Initialize(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<Epetra_Vector> systemvector3  ///< distributed vector that may be filled by
                                                   ///< assembly of element contributions
    );

    /// initialization routine called by the manager ctor
    void Initialize(Teuchos::ParameterList&
            params  ///< parameter list to communicate between elements and discretization
    );

    /// initialization routine called at restart to activate the right conditions
    void Initialize(const double& time  ///< current time
    );

    //! Evaluate routine to call from outside. In here the right action is determined and the
    //! #EvaluateConstraint routine is called
    void Evaluate(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<LINALG::SparseOperator>
            systemmatrix1,  ///< sparse matrix that may be filled by assembly of element
                            ///< contributions
        Teuchos::RCP<LINALG::SparseOperator>
            systemmatrix2,  ///< sparse (rectangular) matrix that may be filled by assembly of
                            ///< element contributions
        Teuchos::RCP<Epetra_Vector> systemvector1,  ///< distributed vector that may be filled by
                                                    ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> systemvector2,  ///< distributed vector that may be filled by
                                                    ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> systemvector3   ///< distributed vector that may be filled by
                                                    ///< assembly of element contributions
    );

   protected:
    std::map<int, double> penalties_;          ///< map containing penalty values
    std::map<int, double> rho_;                ///< map containing rhos for augmented Lagrange
    Teuchos::RCP<Epetra_Export> errorexport_;  ///< exporter for redundant and non-overlapping maps
    Teuchos::RCP<Epetra_Import> errorimport_;  ///< importer for redundant and non-overlapping maps
    Teuchos::RCP<Epetra_Map> rederrormap_;     ///< redundant map of errors
    Teuchos::RCP<Epetra_Map> errormap_;        ///< non-overlapping map of errors
    Teuchos::RCP<Epetra_Vector> initerror_;    ///< initial value of bc
    Teuchos::RCP<Epetra_Vector> acterror_;     ///< current value of bc
    Teuchos::RCP<Epetra_Vector>
        lagrvalues_;  ///< value of Lagrange multiplier in augmented Lagrange
    Teuchos::RCP<Epetra_Vector>
        lagrvalues_force_;  ///< value of Lagrange multiplier in augmented Lagrange


   private:
    // don't want = operator, cctor and destructor

    ConstraintPenalty operator=(const ConstraintPenalty& old);
    ConstraintPenalty(const ConstraintPenalty& old);


    //! Evaluate constraint conditions and assemble the results
    void EvaluateConstraint(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<LINALG::SparseOperator>
            systemmatrix1,  ///< sparse matrix that may be filled by assembly of element
                            ///< contributions
        Teuchos::RCP<LINALG::SparseOperator>
            systemmatrix2,  ///< sparse (rectangular) matrix that may be filled by assembly of
                            ///< element contributions
        Teuchos::RCP<Epetra_Vector> systemvector1,  ///< distributed vector that may be filled by
                                                    ///< aasembly of element contributions
        Teuchos::RCP<Epetra_Vector> systemvector2,  ///< distributed vector that may be filled by
                                                    ///< aasembly of element contributions
        Teuchos::RCP<Epetra_Vector> systemvector3   ///< distributed vector that may be filled by
                                                    ///< aasembly of element contributions
    );

    //! Compute and assemble initial constraint values (depending on user specific activation times)
    void EvaluateError(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<Epetra_Vector> systemvector  ///< distributed vector that may be filled by
                                                  ///< aasembly of element contributions
    );
  };  // class
}  // namespace UTILS

#endif /*CONSTRAINT_H*/
