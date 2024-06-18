/*----------------------------------------------------------------------*/
/*! \file

\brief Monolithic coupling of 3D structure Cardiovascular0D models

\level 2


An arterial 0D flow model derived from physical considerations of mass and momentum balance in the
proximal and distal arterial part, incl. valves (formulation proposed by Cristobal Bertoglio)
(DESIGN SURF CARDIOVASCULAR 0D ARTERIAL PROX DIST CONDITIONS):

-> variables are p_v, p_arp, q_arp, p_ard

      [dV_v/dt + (p_v - p_at)/R_atv(p_v,p_at) + (p_v - p_arp)/R_arv(p_v,p_arp)]   [ 0 ]
Res = [C_arp * d(p_arp)/dt + q_arp - (p_v - p_arp)/R_arv                      ]   [ 0 ]
      [(L_arp/R_arp) * d(q_arp)/dt + q_arp + (p_ard - p_arp)/R_arp            ] = [ 0 ]
      [C_ard * d(p_ard)/dt + (p_ard - p_ref)/R_ard - q_arp                    ]   [ 0 ]

with nonlinear valve resistances R_atv(p_v,p_at), R_arv(p_v,p_arp) - caution when using this since
its physical correctness is doubted by the code author! - reproduce classical piecewise linear
valves with k_p -> 0

*----------------------------------------------------------------------*/

#ifndef FOUR_C_CARDIOVASCULAR0D_ARTERIALPROXDIST_HPP
#define FOUR_C_CARDIOVASCULAR0D_ARTERIALPROXDIST_HPP

#include "4C_config.hpp"

#include "4C_cardiovascular0d.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_inpar_cardiovascular0d.hpp"

#include <Epetra_FECrsMatrix.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::LinAlg
{
  class SparseMatrix;
  class SparseOperator;
}  // namespace Core::LinAlg

namespace UTILS
{
  class Cardiovascular0DArterialProxDist : public Cardiovascular0D

  {
   public:
    /*!
    \brief Constructor of a Cardiovascular0D based on conditions with a given name. It also
    takes care of the Cardiovascular0D IDs.
    */

    Cardiovascular0DArterialProxDist(Teuchos::RCP<Core::FE::Discretization>
                                         discr,  ///< discretization where Cardiovascular0D lives on
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
    void evaluate(
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
        const Teuchos::RCP<Epetra_Vector> sysvec4, Teuchos::RCP<Epetra_Vector> sysvec5) override;

   private:
    // don't want = operator, cctor and destructor

    Cardiovascular0DArterialProxDist operator=(const Cardiovascular0DArterialProxDist& old);
    Cardiovascular0DArterialProxDist(const Cardiovascular0DArterialProxDist& old);



  };  // class
}  // namespace UTILS

FOUR_C_NAMESPACE_CLOSE

#endif
