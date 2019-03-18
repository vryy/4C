/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_utils.H

 \brief helper functions/classes for scalar transport within multiphase porous medium

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_UTILS_H_
#define SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_UTILS_H_

#include "../drt_inpar/inpar_poromultiphase_scatra.H"

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"

// forward declarations
namespace LINALG
{
  class MultiMapExtractor;
}
namespace DRT
{
  class Discretization;
}
namespace ADAPTER
{
  class Coupling;
}

namespace POROMULTIPHASESCATRA
{
  class PoroMultiPhaseScaTraBase;
  class PoroMultiPhaseScaTraArtCouplBase;

  namespace UTILS
  {
    //! setup discretizations and dofsets
    void SetupDiscretizationsAndFieldCoupling(const Epetra_Comm& comm,
        const std::string& struct_disname, const std::string& fluid_disname,
        const std::string& scatra_disname, int& ndsporo_disp, int& ndsporo_vel,
        int& ndsporo_solidpressure, int& ndsporofluid_scatra, const bool artery_coupl);

    //! exchange material pointers of discretizations
    void AssignMaterialPointers(const std::string& struct_disname, const std::string& fluid_disname,
        const std::string& scatra_disname, const bool artery_coupl);

    //! create solution algorithm depending on input file
    Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase>
    CreatePoroMultiPhaseScatraAlgorithm(INPAR::POROMULTIPHASESCATRA::SolutionSchemeOverFields
                                            solscheme,  //!< solution scheme to build (i)
        const Teuchos::ParameterList& timeparams,       //!< problem parameters (i)
        const Epetra_Comm& comm                         //!< communicator(i)
    );

    //! create coupling strategy for coupling with 1D network depending on input file
    Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplBase>
    CreateAndInitArteryCouplingStrategy(Teuchos::RCP<DRT::Discretization> arterydis,
        Teuchos::RCP<DRT::Discretization> contdis, const Teuchos::ParameterList& meshtyingparams,
        const std::string& condname, const std::string& artcoupleddofname,
        const std::string& contcoupleddofname);

    /**
     * Determine norm of vector
     * @param norm [in]: norm to use
     * @param vect [in]: the vector of interest
     * @return: the norm
     */
    double CalculateVectorNorm(const enum INPAR::POROMULTIPHASESCATRA::VectorNorm norm,
        const Teuchos::RCP<const Epetra_Vector> vect);

    /**
     * Get oxygen partial pressure from oxygen concentration via numerical inversion
     * templated to FAD and double
     *
     * @param Pb [out]: oxygen partial pressure
     * @param CaO2 [in]: oxygen concentration
     * @param CaO2_max [in]: maximum oxygen concentration
     * @param Pb50 [in]: partial pressure at 50% maximum oxygen concentration
     * @param n [in]: exponent in Hill equation
     * @param alpha_eff [in]: effective solubility of oxygen in blood
     */
    template <typename T>
    void GetOxyPartialPressureFromConcentration(T& Pb, const T& CaO2, const double& CaO2_max,
        const double& Pb50, const double& n, const double& alpha_eff)
    {
      // start value
      Pb = Pb50 * 2.0 * CaO2 / CaO2_max;

      bool converged = false;
      // Newton loop
      for (int i = 0; i < 20; i++)
      {
        // function
        T f = (std::pow(Pb50, n) + std::pow(Pb, n)) * CaO2 - CaO2_max * std::pow(Pb, n) -
              Pb * (std::pow(Pb, n) + std::pow(Pb50, n)) * alpha_eff;
        if (fabs(f) < 1.0e-10)
        {
          converged = true;
          break;
        }
        // deriv
        T dfdPb = n * std::pow(Pb, (n - 1)) * CaO2 - CaO2_max * n * std::pow(Pb, (n - 1)) -
                  ((n + 1) * std::pow(Pb, n) + std::pow(Pb50, n)) * alpha_eff;
        // update
        Pb = Pb - f / dfdPb;
      }
      if (!converged)
        dserror("local Newton for computation of oxygen partial pressure unconverged");
    }

  }  // namespace UTILS
  //! Print the logo
  void PrintLogo();
}  // namespace POROMULTIPHASESCATRA



#endif /* SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_UTILS_H_ */
