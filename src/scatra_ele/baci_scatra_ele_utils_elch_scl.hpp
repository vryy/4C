/*----------------------------------------------------------------------*/
/*! \file

\brief utility class supporting element evaluation with scls

\level 2

 */
/*----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_ELE_UTILS_ELCH_SCL_HPP
#define BACI_SCATRA_ELE_UTILS_ELCH_SCL_HPP
#include "baci_config.hpp"

#include "baci_inpar_elch.hpp"
#include "baci_scatra_ele_parameter_elch.hpp"
#include "baci_scatra_ele_utils_elch_diffcond.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declarations
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerScl;
    class ScaTraEleDiffManagerElchScl;

    // class implementation
    template <CORE::FE::CellType distype>
    class ScaTraEleUtilsElchScl : public ScaTraEleUtilsElchDiffCond<distype>
    {
      //! abbreviations
      using myelch = ScaTraEleUtilsElch<distype>;
      using mydiffcond = ScaTraEleUtilsElchDiffCond<distype>;

     public:
      //! singleton access method
      static ScaTraEleUtilsElchScl<distype>* Instance(
          const int numdofpernode,    ///< number of degrees of freedom per node
          const int numscal,          ///< number of transported scalars per node
          const std::string& disname  ///< name of discretization
      );

      //! evaluate electrolyte material
      void MatElchMat(Teuchos::RCP<const MAT::Material> material,  //!< electrolyte material
          const std::vector<double>& concentrations,               //!< local concentration values
          double temperature,                                      //!< temperature
          Teuchos::RCP<ScaTraEleDiffManagerElchScl> diffmanager,   //!< diffusion manager
          INPAR::ELCH::DiffCondMat& diffcondmat                    //!< ion type
      );
      //! evaluate electrolyte phase
      void MatElchPhase(Teuchos::RCP<const MAT::Material> material,  //!< electrolyte phase
          const std::vector<double>& concentrations,                 //!< local concentration values
          double temperature,                                        //!< temperature
          Teuchos::RCP<ScaTraEleDiffManagerElchScl> diffmanager,     //!< diffusion manager
          INPAR::ELCH::DiffCondMat& diffcondmat                      //!< ion type
      );

      //! evaluate Space Charge Layer Material
      void MatScl(Teuchos::RCP<const MAT::Material> material, const double concentration,
          const double temperature, Teuchos::RCP<ScaTraEleDiffManagerElchScl> diffmanager);

     protected:
      //! private constructor for singletons
      ScaTraEleUtilsElchScl(const int numdofpernode,  ///< number of degrees of freedom per node
          const int numscal,                          ///< number of transported scalars per node
          const std::string& disname                  ///< name of discretization
      );
    };
  }  // namespace ELEMENTS
}  // namespace DRT
BACI_NAMESPACE_CLOSE

#endif
