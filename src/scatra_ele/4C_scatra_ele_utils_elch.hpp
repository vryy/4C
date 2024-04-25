/*----------------------------------------------------------------------*/
/*! \file

\brief utility class supporting element evaluation for electrochemistry problems


\level 2
 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_UTILS_ELCH_HPP
#define FOUR_C_SCATRA_ELE_UTILS_ELCH_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_inpar_elch.hpp"
#include "4C_lib_element.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    class ScaTraEleDiffManagerElch;

    // class implementation
    template <CORE::FE::CellType distype>
    class ScaTraEleUtilsElch
    {
     protected:
      //! number of element nodes
      static constexpr int nen_ = CORE::FE::num_nodes<distype>;

      //! number of space dimensions
      static constexpr int nsd_ = CORE::FE::dim<distype>;

     public:
      //! singleton access method
      static ScaTraEleUtilsElch<distype>* Instance(
          const int numdofpernode,    ///< number of degrees of freedom per node
          const int numscal,          ///< number of transported scalars per node
          const std::string& disname  ///< name of discretization
      );

      //! destructor
      virtual ~ScaTraEleUtilsElch() = default;

      //! evaluation of electrochemistry kinetics at integration point on domain or boundary element
      void EvaluateElchKineticsAtIntegrationPoint(const DRT::Element* ele,  ///< current element
          CORE::LINALG::SerialDenseMatrix& emat,                            ///< element matrix
          CORE::LINALG::SerialDenseVector& erhs,  ///< element right-hand side vector
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
              ephinp,  ///< state variables at element nodes
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
              ehist,                                   ///< history variables at element nodes
          const double timefac,                        ///< time factor
          const double fac,                            ///< Gauss integration factor
          const CORE::LINALG::Matrix<nen_, 1>& funct,  ///< shape functions at int. point
          const Teuchos::RCP<DRT::Condition>& cond,    ///< condition
          const int nume,                              ///< number of transferred electrons
          const std::vector<int>& stoich,              ///< stoichiometry of the reaction
          const double valence_k,                      ///< valence of the single reactant
          const int kinetics,                          ///< desired electrode kinetics model
          const double pot0,                           ///< actual electrode potential on metal side
          const double frt,                            ///< factor F/RT
          const double fns,     ///< factor fns = s_k / (nume * faraday * (-1))
          const double scalar,  ///< scaling factor for element matrix and right-hand side vector
                                ///< contributions
          const int k           ///< index of evaluated scalar
      ) const;

      //! evaluate electrode kinetics status information at integration point on domain or boundary
      //! element
      void EvaluateElectrodeStatusAtIntegrationPoint(const DRT::Element* ele,  ///< current element
          CORE::LINALG::SerialDenseVector& scalars,  ///< scalars to be computed
          const Teuchos::ParameterList& params,      ///< parameter list
          const Teuchos::RCP<DRT::Condition>& cond,  ///< condition
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
              ephinp,  ///< nodal values of concentration and electric potential
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
              ephidtnp,                                ///< nodal time derivative vector
          const CORE::LINALG::Matrix<nen_, 1>& funct,  ///< shape functions at integration point
          const int zerocur,                           ///< flag for zero current
          const int kinetics,                          ///< desired electrode kinetics model
          const std::vector<int>& stoich,              ///< stoichiometry of the reaction
          const int nume,                              ///< number of transferred electrons
          const double pot0,     ///< actual electrode potential on metal side at t_{n+1}
          const double frt,      ///< factor F/RT
          const double timefac,  ///< factor due to time discretization
          const double fac,      ///< integration factor
          const double scalar,   ///< scaling factor for current related quantities
          const int k            ///< index of evaluated scalar
      ) const;

      //! evaluate ion material
      void MatIon(const Teuchos::RCP<const MAT::Material> material,  //!< ion material
          const int k,                                               //!< ID of ion material
          const INPAR::ELCH::EquPot equpot,  //!< type of closing equation for electric potential
          const Teuchos::RCP<ScaTraEleDiffManagerElch>& diffmanager  //!< diffusion manager
      );

     protected:
      //! protected constructor for singletons
      ScaTraEleUtilsElch(const int numdofpernode,  ///< number of degrees of freedom per node
          const int numscal,                       ///< number of transported scalars per node
          const std::string& disname               ///< name of discretization
      );

     private:
      //! number of degrees of freedom per node
      const int numdofpernode_;

      //! number of transported scalars
      const int numscal_;
    };  // class ScaTraEleUtilsElch
  }     // namespace ELEMENTS
}  // namespace DRT
FOUR_C_NAMESPACE_CLOSE

#endif
