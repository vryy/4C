/*----------------------------------------------------------------------*/
/*! \file

\brief Monolithic coupling of 3D structure Cardiovascular0D models

\level 2


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CARDIOVASCULAR0D_4ELEMENTWINDKESSEL_HPP
#define FOUR_C_CARDIOVASCULAR0D_4ELEMENTWINDKESSEL_HPP

#include "baci_config.hpp"

#include "baci_cardiovascular0d.hpp"
#include "baci_discretization_fem_general_utils_integration.hpp"
#include "baci_inpar_cardiovascular0d.hpp"

#include <Epetra_FECrsMatrix.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace DRT
{
  class Condition;
  class Discretization;
}  // namespace DRT

namespace CORE::LINALG
{
  class SparseMatrix;
  class SparseOperator;
}  // namespace CORE::LINALG

namespace UTILS
{
  /*! \brief Four-element Windkessel model without valves
   *
   * Boundary condition: DESIGN SURF CARDIOVASCULAR 0D WINDKESSEL ONLY CONDITIONS
   *
   * Original equation is:
   * C * dp/dt + (p - p_ref)/R_p - (1 + Z_c/R_p) * Q - (C Z_c  + * L/R_p) * dQ/dt - L C * d2Q/dt2 =
   * 0 where Q := -dV/dt
   *
   * -> we reformulate the ODE to three first-order ODEs with variables p, q, s:
   *
   *       [C * dp/dt + (p - p_ref)/R_p + (1 + Z_c/R_p) * q + (C Z_c  + L/R_p) * s + L C * ds/dt] [
   * 0 ] Res = [dV/dt - q ] = [ 0 ] [dq/dt - s ]   [ 0 ]
   *
   * The classical 3- or 2-element windkessel models are reproduced by setting L, or L and Z_c to
   * zero, respectively
   *
   */
  class Cardiovascular0D4ElementWindkessel : public Cardiovascular0D

  {
   public:
    /*!
    \brief Constructor of a Cardiovascular0D based on conditions with a given name. It also
    takes care of the Cardiovascular0D IDs.
    */

    Cardiovascular0D4ElementWindkessel(
        Teuchos::RCP<DRT::Discretization>
            discr,                         ///< Discretization where Cardiovascular0D lives on
        const std::string& conditionname,  ///< Name of condition to create Cardiovascular0D from
        std::vector<int>& curID            ///< current ID
    );



    /// initialization routine called by the manager ctor to get correct reference base values and
    /// activating the right conditions at the beginning
    void Initialize(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<Epetra_Vector> sysvec1,  ///< distributed vector that may be filled by assembly
                                              ///< of element contributions
        Teuchos::RCP<Epetra_Vector>
            sysvec2  ///< distributed vector that may be filled by assembly of element contributions
        ) override;

    //! Evaluate routine to call from outside. In here the right action is determined and the
    //! #EvaluateCardiovascular0D routine is called
    void Evaluate(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat1,  ///< Cardiovascular0D stiffness matrix
        Teuchos::RCP<CORE::LINALG::SparseOperator>
            sysmat2,  ///< Cardiovascular0D offdiagonal matrix dV/dd
        Teuchos::RCP<CORE::LINALG::SparseOperator>
            sysmat3,                          ///< Cardiovascular0D offdiagonal matrix dfext/dp
        Teuchos::RCP<Epetra_Vector> sysvec1,  ///< distributed vectors that may be filled by
                                              ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> sysvec2, Teuchos::RCP<Epetra_Vector> sysvec3,
        const Teuchos::RCP<Epetra_Vector> sysvec4, Teuchos::RCP<Epetra_Vector> sysvec5) override;

   private:
    // don't want = operator, cctor and destructor

    Cardiovascular0D4ElementWindkessel operator=(const Cardiovascular0D4ElementWindkessel& old);
    Cardiovascular0D4ElementWindkessel(const Cardiovascular0D4ElementWindkessel& old);



  };  // class
}  // namespace UTILS

FOUR_C_NAMESPACE_CLOSE

#endif
