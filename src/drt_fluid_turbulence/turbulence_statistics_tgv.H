/*----------------------------------------------------------------------*/
/*! \file

\brief Compute dissipation rates for Taylor-Green vortex
       and write them to files.


\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef TGV_STATISTICS_H
#define TGV_STATISTICS_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_scatra.H"


namespace FLD
{
  class TurbulenceStatisticsTgv
  {
   public:
    /*!
    \brief Standard Constructor (public)
    */
    TurbulenceStatisticsTgv(Teuchos::RCP<DRT::Discretization> actdis,
        Teuchos::ParameterList& params, const std::string& statistics_outfilename);

    /*!
    \brief Destructor

    */
    virtual ~TurbulenceStatisticsTgv();

    //! @name functions for (spatial) averaging

    /*!
    \brief evaluate dissipation
    */
    void EvaluateResiduals(std::map<std::string, Teuchos::RCP<Epetra_Vector>> statevecs,
        std::map<std::string, Teuchos::RCP<Epetra_MultiVector>> statetenss,
        const double thermpressaf, const double thermpressam, const double thermpressdtaf,
        const double thermpressdtam);

    //@}

    //! @name Miscellaneous

    /*!
    \brief Dump the result to file.
    */

    void DumpStatistics(const int step);

    /*!
    \brief reset values
    */

    void ClearStatistics();


   private:
    //! time step size
    double dt_;

    //! number of elements in sample plane
    int numele_;

    //! The discretization (required for nodes, dofs etc;)
    Teuchos::RCP<DRT::Discretization> discret_;

    //! contains plane normal direction etc --- this is the original
    //! fluid dynamic parameterlist
    Teuchos::ParameterList& params_;
    //! name of statistics output file, despite the ending
    const std::string statistics_outfilename_;

    //! parameterlist for the element call when averages of residuals
    //! are calculated --- used for communication between element
    //! and averaging methods --- for fluid field
    Teuchos::ParameterList eleparams_;

    //! the dim_-coordinates of the homogeneous planes containing nodes
    Teuchos::RCP<std::vector<double>> nodeplanes_;

    //!--------------------------------------------------
    //!  averaged resdiuals and subscale quantities
    //!--------------------------------------------------

    //! sum over all in plane element sizes
    Teuchos::RCP<std::vector<double>> sumhk_;
    //! sum over all in plane element sizes for the Bazilevs
    //! parameter in the viscous regime
    Teuchos::RCP<std::vector<double>> sumhbazilevs_;
    //! sum over all in plane stream lengths
    Teuchos::RCP<std::vector<double>> sumstrle_;
    //! sum over all in plane stream lengths
    Teuchos::RCP<std::vector<double>> sumgradle_;

    //! sum over all in plane residuals
    Teuchos::RCP<std::vector<double>> sumtauM_;
    //! sum over all in plane squared residuals
    Teuchos::RCP<std::vector<double>> sumtauC_;

    //! sum over all in plane mk (parameter for stabilization parameter, 1/3 for lin ele)
    Teuchos::RCP<std::vector<double>> summk_;

    //! sum over all in plane residuals
    Teuchos::RCP<std::vector<double>> sumres_;
    //! sum over all in plane squared residuals
    Teuchos::RCP<std::vector<double>> sumres_sq_;
    //! sum over all in plane residuals norms
    Teuchos::RCP<std::vector<double>> sumabsres_;
    //! sum over all in plane values of svel/tau
    Teuchos::RCP<std::vector<double>> sumtauinvsvel_;
    //! sum over all in plane subscale velocities
    Teuchos::RCP<std::vector<double>> sumsvelaf_;
    //! sum over all in plane squared subscale velocities
    Teuchos::RCP<std::vector<double>> sumsvelaf_sq_;
    //! sum over all in plane subscale velocities norms
    Teuchos::RCP<std::vector<double>> sumabssvelaf_;

    //! sum over all in plane residuals of the continuity equation
    Teuchos::RCP<std::vector<double>> sumresC_;
    //! sum over all in plane squared residuals of the continuity equation
    Teuchos::RCP<std::vector<double>> sumresC_sq_;
    //! sum over all in plane subscale pressure values at current timestep
    Teuchos::RCP<std::vector<double>> sumspressnp_;
    //! sum over all in plane squared subscale pressure values at current timestep
    Teuchos::RCP<std::vector<double>> sumspressnp_sq_;

    //! sum over all in plane averaged dissipation rates from pspg stabilisation
    Teuchos::RCP<std::vector<double>> sum_eps_pspg_;
    //! sum over all in plane averaged dissipation rates from supg stabilisation
    Teuchos::RCP<std::vector<double>> sum_eps_supg_;
    //! sum over all in plane averaged dissipation rates from cross term
    Teuchos::RCP<std::vector<double>> sum_eps_cross_;
    //! sum over all in plane averaged dissipation rates from reynolds term
    Teuchos::RCP<std::vector<double>> sum_eps_rey_;
    //! sum over all in plane averaged dissipation rates from least squares continuity term
    Teuchos::RCP<std::vector<double>> sum_eps_graddiv_;
    //! sum over all in plane averaged dissipation rates from eddy viscosity model (Smagorinsky)
    Teuchos::RCP<std::vector<double>> sum_eps_eddyvisc_;
    //! sum over all in plane averaged dissipation rates from eddy viscosity model (AVM3)
    Teuchos::RCP<std::vector<double>> sum_eps_avm3_;
    //! sum over all in plane averaged dissipation rates from mfs model
    Teuchos::RCP<std::vector<double>> sum_eps_mfs_;
    //! sum over all in plane averaged dissipation rates from ,fs model (forwardscatter)
    Teuchos::RCP<std::vector<double>> sum_eps_mfscross_;
    //! sum over all in plane averaged dissipation rates from mfs model (backscatter)
    Teuchos::RCP<std::vector<double>> sum_eps_mfsrey_;
    //! sum over all in plane averaged dissipation rates from Galerkin viscous term
    Teuchos::RCP<std::vector<double>> sum_eps_visc_;
    //! sum over all in plane averaged dissipation rates from Galerkin convective term
    Teuchos::RCP<std::vector<double>> sum_eps_conv_;

    //! sum over all in plane averaged supg+cross stress
    Teuchos::RCP<std::vector<double>> sum_crossstress_;
    //! sum over all in plane averaged reynolds stress
    Teuchos::RCP<std::vector<double>> sum_reystress_;
  };

}  // namespace FLD

#endif  // not TGV_STATISTICS_H
