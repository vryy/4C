/*----------------------------------------------------------------------*/
/*! \file

\brief Monolithic coupling of 3D structure Cardiovascular0D models

\level 2


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CARDIOVASCULAR0D_HPP
#define FOUR_C_CARDIOVASCULAR0D_HPP

#include "4C_config.hpp"

#include "4C_discretization_condition.hpp"
#include "4C_discretization_fem_general_utils_integration.hpp"
#include "4C_inpar_cardiovascular0d.hpp"

#include <Epetra_FECrsMatrix.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Discret
{
  class Discretization;
}  // namespace Discret

namespace Core::LinAlg
{
  class SparseMatrix;
  class SparseOperator;
}  // namespace Core::LinAlg

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

    Cardiovascular0D(Teuchos::RCP<Discret::Discretization>
                         discr,            ///< discretization where Cardiovascular0D lives on
        const std::string& conditionname,  ///< Name of condition to create Cardiovascular0D from
        std::vector<int>& curID            ///< current ID
    );

    /*!
    \brief Constructor of a Cardiovascular0D based on a conditions with a given name.
    */

    Cardiovascular0D(Teuchos::RCP<Discret::Discretization>
                         discr,  ///< discretization where Cardiovascular0D funtion lives on
        const std::string&
            CondName  ///< Name of condition to create Cardiovascular0D functions from
    );

    /*!
        \brief Destructor

     */
    virtual ~Cardiovascular0D() = default;

    /*!
     \brief Return if there are Cardiovascular0D functions
    */
    bool have_cardiovascular0_d() { return cardiovascular0dtype_ != none; };

    /// Set state of the underlying discretization
    void set_state(const std::string& state,  ///< name of state to set
        Teuchos::RCP<Epetra_Vector> V         ///< values to set
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
        Teuchos::RCP<Core::LinAlg::SparseMatrix> sysmat1,  ///< Cardiovascular0D stiffness matrix
        Teuchos::RCP<Core::LinAlg::SparseOperator>
            sysmat2,  ///< Cardiovascular0D offdiagonal matrix dV/dd
        Teuchos::RCP<Core::LinAlg::SparseOperator>
            sysmat3,                          ///< Cardiovascular0D offdiagonal matrix dfext/dp
        Teuchos::RCP<Epetra_Vector> sysvec1,  ///< distributed vectors that may be filled by
                                              ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> sysvec2, Teuchos::RCP<Epetra_Vector> sysvec3,
        const Teuchos::RCP<Epetra_Vector> sysvec4, Teuchos::RCP<Epetra_Vector> sysvec5);

    /// Return type of Cardiovascular0D function
    Cardiovascular0DType Type() { return cardiovascular0dtype_; }

    std::vector<Core::Conditions::Condition*> get_cardiovascular0_d_condition()
    {
      return cardiovascular0dcond_;
    }
    std::vector<Core::Conditions::Condition*> get_cardiovascular0_d_structure_coupling_condition()
    {
      return cardiovascular0dstructcoupcond_;
    }

    Inpar::CARDIOVASCULAR0D::Cardvasc0DRespiratoryModel GetRespiratoryModel()
    {
      return respiratory_model_;
    }

    //! Evaluate Cardiovascular0D conditions and assemble the results
    void EvaluateDStructDp(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<Core::LinAlg::SparseOperator>
            sysmat  ///< Cardiovascular0D offdiagonal matrix dfext/dp
    );

   protected:
    Teuchos::RCP<Discret::Discretization> actdisc_;  ///< standard discretization
    std::vector<Core::Conditions::Condition*>
        cardiovascular0dcond_;  ///< 0D cardiovascular conditions
    std::vector<Core::Conditions::Condition*>
        cardiovascular0dstructcoupcond_;  ///< 0D cardiovascular structure coupling conditions
    Cardiovascular0DType cardiovascular0dtype_;  ///< Cardiovascular0D type
    const Inpar::CARDIOVASCULAR0D::Cardvasc0DAtriumModel atrium_model_;
    const Inpar::CARDIOVASCULAR0D::Cardvasc0DVentricleModel ventricle_model_;
    const Inpar::CARDIOVASCULAR0D::Cardvasc0DRespiratoryModel respiratory_model_;
    //! gaussian integration to be used
    Core::FE::GaussRule2D gaussrule_;

   private:
    // don't want = operator, cctor and destructor

    Cardiovascular0D operator=(const Cardiovascular0D& old);
    Cardiovascular0D(const Cardiovascular0D& old);

    //! Return the Cardiovascular0DType based on the condition name
    Cardiovascular0DType get_cardiovascular0_d_type(const std::string& Name  ///< condition name
    );

  };  // class
}  // namespace UTILS

FOUR_C_NAMESPACE_CLOSE

#endif
