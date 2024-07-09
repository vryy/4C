/*----------------------------------------------------------------------*/
/*! \file

\brief utility class supporting element evaluation for concentrated electrolytes

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_UTILS_ELCH_DIFFCOND_HPP
#define FOUR_C_SCATRA_ELE_UTILS_ELCH_DIFFCOND_HPP

#include "4C_config.hpp"

#include "4C_inpar_elch.hpp"
#include "4C_scatra_ele_utils_elch_electrode.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    // forward declarations
    template <int nsd, int nen>
    class ScaTraEleInternalVariableManagerElchDiffCond;
    class ScaTraEleDiffManagerElchDiffCond;

    // class implementation
    template <Core::FE::CellType distype>
    class ScaTraEleUtilsElchDiffCond : public ScaTraEleUtilsElchElectrode<distype>
    {
      //! abbreviations
      typedef ScaTraEleUtilsElch<distype> myelch;
      typedef ScaTraEleUtilsElchElectrode<distype> myelectrode;

     public:
      //! singleton access method
      static ScaTraEleUtilsElchDiffCond<distype>* instance(
          const int numdofpernode,    ///< number of degrees of freedom per node
          const int numscal,          ///< number of transported scalars per node
          const std::string& disname  ///< name of discretization
      );

      //! evaluate electrolyte material
      void mat_elch_mat(Teuchos::RCP<const Core::Mat::Material> material,  //!< electrolyte material
          const std::vector<double>& concentrations,  //!< local concentration values
          double temperature,                         //!< temperature
          Inpar::ElCh::EquPot equpot,  //!< type of closing equation for electric potential
          double ffrt,                 //!< factor F^2/RT
          Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> diffmanager,  //!< diffusion manager
          Inpar::ElCh::DiffCondMat& diffcondmat                        //!< ion type
      );

      //! evaluate electrolyte phase
      void mat_elch_phase(Teuchos::RCP<const Core::Mat::Material> material,  //!< electrolyte phase
          const std::vector<double>& concentrations,  //!< local concentration values
          double temperature,                         //!< temperature
          const Inpar::ElCh::EquPot& equpot,  //!< type of closing equation for electric potential
          const double& ffrt,                 //!< factor F^2/RT
          Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> diffmanager,  //!< diffusion manager
          Inpar::ElCh::DiffCondMat& diffcondmat                        //!< ion type
      );

      //! evaluate standard Newman material
      void mat_newman(Teuchos::RCP<const Core::Mat::Material> material,  //!< Newman material
          double concentration,                                       //!< local concentration value
          double temperature,                                         //!< temperature
          Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> diffmanager  //!< diffusion manager
      );

     protected:
      //! private constructor for singletons
      ScaTraEleUtilsElchDiffCond(
          const int numdofpernode,    ///< number of degrees of freedom per node
          const int numscal,          ///< number of transported scalars per node
          const std::string& disname  ///< name of discretization
      );

    };  // class ScaTraEleUtilsElchDiffCond
  }     // namespace ELEMENTS
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
