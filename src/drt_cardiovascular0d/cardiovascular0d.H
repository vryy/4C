/*!----------------------------------------------------------------------
\file cardiovascular0d.H

\brief Monolithic coupling of 3D structure Cardiovascular0D models

\level 2

<pre>
\maintainer Amadeus Gebauer
</pre>

*----------------------------------------------------------------------*/

#ifndef CARDIOVASCULAR0D_H
#define CARDIOVASCULAR0D_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Teuchos_ParameterList.hpp>
#include "Epetra_FECrsMatrix.h"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_inpar/inpar_cardiovascular0d.H"


// forward declarations
namespace DRT
{
  class Condition;
  class Discretization;
}  // namespace DRT

namespace LINALG
{
  class SparseMatrix;
  class SparseOperator;
}  // namespace LINALG

namespace UTILS
{
  class Cardiovascular0D

  {
   public:
    //! 0D cardiovascular types
    enum Cardiovascular0DType
    {
      none,
      cardvasc0d_4elementwindkessel,
      cardvasc0d_arterialproxdist,
      cardvasc0d_syspulcirculation,
      cardvascrespir0d_syspulperiphcirculation
    };

    /*!
    \brief Constructor of a Cardiovascular0D based on conditions with a given name. It also
    takes care of the Cardiovascular0D IDs.
    */

    Cardiovascular0D(Teuchos::RCP<DRT::Discretization>
                         discr,       ///< Discretization where Cardiovascular0D lives on
        const std::string& CondName,  ///< Name of condition to create Cardiovascular0D from
        std::vector<int>& curID       ///< current ID
    );

    /*!
    \brief Constructor of a Cardiovascular0D based on a conditions with a given name.
    */

    Cardiovascular0D(Teuchos::RCP<DRT::Discretization>
                         discr,  ///< Discretization where Cardiovascular0D funtion lives on
        const std::string&
            CondName  ///< Name of condition to create Cardiovascular0D functions from
    );

    /*!
        \brief Destructor

     */
    virtual ~Cardiovascular0D() { return; };

    /*!
     \brief Return if there are Cardiovascular0D functions
    */
    bool HaveCardiovascular0D() { return cardiovascular0dtype_ != none; };

    /// Set state of the underlying discretization
    void SetState(const std::string& state,  ///< name of state to set
        Teuchos::RCP<Epetra_Vector> V        ///< values to set
    );


    /// initialization routine called by the manager ctor to get correct reference base values and
    /// activating the right conditions at the beginning
    virtual void Initialize(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<Epetra_Vector> sysvec1,  ///< distributed vector that may be filled by assembly
                                              ///< of element contributions
        Teuchos::RCP<Epetra_Vector>
            sysvec2  ///< distributed vector that may be filled by assembly of element contributions
    );

    //! Evaluate routine to call from outside. In here the right action is determined and the
    //! #EvaluateCardiovascular0D routine is called
    virtual void Evaluate(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<LINALG::SparseMatrix> sysmat1,  ///< Cardiovascular0D stiffness matrix
        Teuchos::RCP<LINALG::SparseOperator>
            sysmat2,  ///< Cardiovascular0D offdiagonal matrix dV/dd
        Teuchos::RCP<LINALG::SparseOperator>
            sysmat3,                          ///< Cardiovascular0D offdiagonal matrix dfext/dp
        Teuchos::RCP<Epetra_Vector> sysvec1,  ///< distributed vectors that may be filled by
                                              ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> sysvec2, Teuchos::RCP<Epetra_Vector> sysvec3,
        const Teuchos::RCP<Epetra_Vector> sysvec4, Teuchos::RCP<Epetra_Vector> sysvec5);

    /// Return type of Cardiovascular0D function
    Cardiovascular0DType Type() { return cardiovascular0dtype_; }

    std::vector<DRT::Condition*> GetCardiovascular0DCondition() { return cardiovascular0dcond_; }
    std::vector<DRT::Condition*> GetCardiovascular0DStructureCouplingCondition()
    {
      return cardiovascular0dstructcoupcond_;
    }

    INPAR::CARDIOVASCULAR0D::Cardvasc0DRespiratoryModel GetRespiratoryModel()
    {
      return respiratory_model_;
    }

    //! Evaluate Cardiovascular0D conditions and assemble the results
    void EvaluateDStructDp(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<LINALG::SparseOperator>
            sysmat  ///< Cardiovascular0D offdiagonal matrix dfext/dp
    );

   protected:
    Teuchos::RCP<DRT::Discretization> actdisc_;          ///< standard discretization
    std::vector<DRT::Condition*> cardiovascular0dcond_;  ///< 0D cardiovascular conditions
    std::vector<DRT::Condition*>
        cardiovascular0dstructcoupcond_;  ///< 0D cardiovascular structure coupling conditions
    Cardiovascular0DType cardiovascular0dtype_;  ///< Cardiovascular0D type
    const INPAR::CARDIOVASCULAR0D::Cardvasc0DAtriumModel atrium_model_;
    const INPAR::CARDIOVASCULAR0D::Cardvasc0DVentricleModel ventricle_model_;
    const INPAR::CARDIOVASCULAR0D::Cardvasc0DRespiratoryModel respiratory_model_;
    //! gaussian integration to be used
    DRT::UTILS::GaussRule2D gaussrule_;

   private:
    // don't want = operator, cctor and destructor

    Cardiovascular0D operator=(const Cardiovascular0D& old);
    Cardiovascular0D(const Cardiovascular0D& old);

    //! Return the Cardiovascular0DType based on the condition name
    Cardiovascular0DType GetCardiovascular0DType(const std::string& Name  ///< condition name
    );

  };  // class
}  // namespace UTILS

#endif /*CARDIOVASCULAR0D_H*/
