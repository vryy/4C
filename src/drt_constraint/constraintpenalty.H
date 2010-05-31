/*!----------------------------------------------------------------------
\file constraintpenalty.H

\brief Basic constraint class, dealing with constraints living on boundaries

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifndef CONSTRAINTPENALTY_H
#define CONSTRAINTPENALTY_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Teuchos_ParameterList.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_sparsematrix.H"
#include "constraint.H"


using namespace std;
using namespace Teuchos;

/*!
\brief UTILS: namespace containing general utils for BACI, as constraints, potentials, surface stress, ...
*/

namespace UTILS
{
  /*!
  \brief Basic constraint class, dealing with constraint and monitor boundary conditions.
  This class cannot handle multi point constraints, they will be dealt with by a derived class.

  \author Kloeppel (kloeppel@lnm.mw.tum.de)
  */
  class ConstraintPenalty : public Constraint

  {
  public:


    /*!
    \brief Constructor of a constraint based on a conditions with a given name.
    */

    ConstraintPenalty
    (
      RCP<DRT::Discretization> discr,  ///< Discretization constraint lives on
      const string& CondName ///< Name of condition to create constraints from
    );

    /*!
        \brief Destructor

     */
    ~ConstraintPenalty()
    {
        return;
    };

    /// initialization routine called by the manager ctor to get correct reference base values and activating the right conditions at the beginning
    void Initialize
    (
      ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      RCP<Epetra_Vector>    systemvector3  ///< distributed vector that may be filled by assembly of element contributions
    );

    /// initialization routine called at restart to activate the right conditions
    void Initialize
    (
      const double& time  ///< current time
    );

    //! Evaluate routine to call from outside. In here the right action is determined and the #EvaluateConstraint routine is called
    void Evaluate
    (
      ParameterList&  params,  ///< parameter list to communicate between elements and discretization
      RCP<LINALG::SparseOperator> systemmatrix1, ///< sparse matrix that may be filled by assembly of element contributions
      RCP<LINALG::SparseOperator> systemmatrix2, ///< sparse (rectangular) matrix that may be filled by assembly of element contributions
      RCP<Epetra_Vector> systemvector1, ///< distributed vector that may be filled by assembly of element contributions
      RCP<Epetra_Vector> systemvector2, ///< distributed vector that may be filled by assembly of element contributions
      RCP<Epetra_Vector> systemvector3  ///< distributed vector that may be filled by assembly of element contributions
    );

  protected:

    map<int,double>  penalties_; ///< map with times at which constraint is supposed to become active
    RCP<Epetra_Export> errorexport_;
    RCP<Epetra_Import> errorimport_;
    RCP<Epetra_Map> rederrormap_;
    RCP<Epetra_Map> errormap_;
    RCP<Epetra_Vector> initerror_;
    RCP<Epetra_Vector> acterror_;

  private:

    // don't want = operator, cctor and destructor

    ConstraintPenalty operator = (const ConstraintPenalty& old);
    ConstraintPenalty(const ConstraintPenalty& old);


    //! Evaluate constraint conditions and assemble the results
    void EvaluateConstraint
    (
      ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      RCP<LINALG::SparseOperator> systemmatrix1, ///< sparse matrix that may be filled by assembly of element contributions
      RCP<LINALG::SparseOperator> systemmatrix2, ///< sparse (rectangular) matrix that may be filled by assembly of element contributions
      RCP<Epetra_Vector>    systemvector1, ///< distributed vector that may be filled by aasembly of element contributions
      RCP<Epetra_Vector>    systemvector2, ///< distributed vector that may be filled by aasembly of element contributions
      RCP<Epetra_Vector>    systemvector3  ///< distributed vector that may be filled by aasembly of element contributions
      );

    //! Compute and assemble initial constraint values (depending on user specific activation times)
    void EvaluateError
    (
      ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      RCP<Epetra_Vector>    systemvector ///< distributed vector that may be filled by aasembly of element contributions
    );
  }; //class
} //namespace

#endif /*CONSTRAINT_H*/
