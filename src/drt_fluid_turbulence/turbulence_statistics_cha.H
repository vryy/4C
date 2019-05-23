/*----------------------------------------------------------------------*/
/*!

\brief Compute (time and space) averaged values for turbulent channel
       flows and write them to files.


o Create set of all available homogeneous planes
  (Construction based on a round robin communication pattern)

o loop planes (e.g. plane coordinates)

  - pointwise in-plane average of first- and second order moments
  - integral in-plane average of first- and second order moments
  - in-plane average of wall force

o in plane mean values are averaged in time over all steps between two
  outputs (by computation of the arithmetic mean)

  - time average pointwise values
  - time average integral values
  - time average forces
  - time average Smagorinsky stuff
  - time average residuals, subscale quantities etc.

o Write pointwise and integral statistics for first and second
  order moments
  ->   .flow_statistic

o Write statistics for the Smagorinsky "constant" Cs if a dynamic
  procedure to determine it is applied

  ->  .Cs_statistic

o Write statistics for subscale quantitites and residuals
  (subscale quantitites and residuals are averaged over element
   layers)
  ->  .res_statistic

Required parameters are the number of velocity degrees of freedom (3),
the normal direction to the plane, in which the average values in space
should be computed, and the basename of the statistics outfile. These
parameters are expected to be contained in the fluid time integration
parameter list given on input.

This method is intended to be called every upres_ steps during fluid
output.

\maintainer Martin Kronbichler

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef TURBULENCE_STATISTICS_H
#define TURBULENCE_STATISTICS_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_control_point.H"
#include "../linalg/linalg_utils.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_scatra.H"


namespace FLD
{
  class XWall;

  class TurbulenceStatisticsCha
  {
   public:
    /*!
    \brief Standard Constructor (public)

        o Create vector of homogeneous plane coordinates

    o Allocate 4 distributed toggle vectors and one distributed vector
      for squares, both for integral and pointwise means

    o allocate all sum_something vectors

    o initialise the output (open/clear files, print header)


    */
    TurbulenceStatisticsCha(Teuchos::RCP<DRT::Discretization> actdis, bool alefluid,
        Teuchos::RCP<Epetra_Vector> dispnp, Teuchos::ParameterList& params,
        const std::string& statistics_outfilename, bool subgrid_dissipation,
        Teuchos::RCP<FLD::XWall> xwallobj);

    /*!
    \brief Destructor

    */
    virtual ~TurbulenceStatisticsCha();

    //! @name functions for (spatial) averaging


    /*!
    \brief Store scatra discretization for passive scalar transport
    to caculate mean values of multifractal subgrid-scales modeling
    parameters
    */
    void StoreScatraDiscretAndParams(Teuchos::RCP<DRT::Discretization> scatradis,
        Teuchos::RCP<Teuchos::ParameterList> scatraparams,
        Teuchos::RCP<Teuchos::ParameterList> scatraextraparams,
        Teuchos::RCP<Teuchos::ParameterList> scatratimeparams);


    /*!
    \brief Compute the in-plane mean values of first and second order
    moments for velocities, pressure and Cs are added to global
    'sum' vectors.
    */
    void DoTimeSample(const Teuchos::RCP<const Epetra_Vector> velnp,
        const Teuchos::RCP<const Epetra_Vector> force);


    /*!
    \brief The values of velocity, pressure, temperature and its squared
    values are added to global vectors. This method allows to do the time
    average of the nodal values after a certain amount of timesteps.
    */
    void DoLomaTimeSample(const Teuchos::RCP<const Epetra_Vector> velnp,
        const Teuchos::RCP<const Epetra_Vector> scanp,
        const Teuchos::RCP<const Epetra_Vector> force, const double eosfac);


    /*!
    \brief The values of velocity, pressure, scalar and its squared
    values are added to global vectors. This method allows to do the time
    average of the nodal values after a certain amount of timesteps.
    */
    void DoScatraTimeSample(const Teuchos::RCP<const Epetra_Vector> velnp,
        const Teuchos::RCP<const Epetra_Vector> scanp,
        const Teuchos::RCP<const Epetra_Vector> force);


    /*!
    \brief Compute in plane means of u,u^2 etc. (integral version)

    The averages here are calculated by integration.

    The calculated value is added to the sum**,sumsq** variables in the
    component corresponding to the plane.

    Further documentation is provided in the element subroutines
    */
    void EvaluateIntegralMeanValuesInPlanes();

    /*!
    \brief Compute in-plane means (integral version) for low-Mach-number flow
    */
    void EvaluateLomaIntegralMeanValuesInPlanes(const double eosfac);

    /*!
    \brief Compute in-plane means (integral version) for turbulent passive scalar transport
    */
    void EvaluateScatraIntegralMeanValuesInPlanes();

    /*!
    \brief Compute in plane means of u,u^2 etc. (nodal quantities)

    The averages here are calculated as the arithmetic mean of
    point values:

    - generate 4 toggle vectors (u,v,w,p), for example

                              /  1  u dof in homogeneous plane
                   toggleu_  |
                              \  0  elsewhere

    - 2 * 4 scalarproducts for in plane mean values

    - apply toggle vectors to pointwise multiplied velnp.*velnp
      for second order moments


    The calculated values are added to the pointsum**,pointsumsq** variables
    in the component corresponding to the plane.

    velnp is the solution vector provided by the time integration
    algorithm
    */
    void EvaluatePointwiseMeanValuesInPlanes();

    /*!
      \brief Add computed dynamic Smagorinsky quantities (Smagorinsky
             constant, effective viscosity and (Cs_delta)^2 used
             during the computation)

      The increment is computed during the computation of Cs
      (in the filtering part of time integration, i.e. during
      an element call in the nonlinear iteration)

      We just store it here and add it to the sum as soon as we do
      the time sample.
    */

    void AddDynamicSmagorinskyQuantities();

    /*!
      \brief Add parameters of multifractal
             subgrid-scales model
    */

    void AddModelParamsMultifractal(const Teuchos::RCP<const Epetra_Vector> velnp,
        const Teuchos::RCP<const Epetra_Vector> fsvelnp, const bool withscatra);

    /*!
      \brief do averaging of residuals, dissipation rates etc
             (all gausspoint-quantities)

    */

    void EvaluateResiduals(std::map<std::string, Teuchos::RCP<Epetra_Vector>> statevecs,
        std::map<std::string, Teuchos::RCP<Epetra_MultiVector>> statetenss,
        const double thermpressaf, const double thermpressam, const double thermpressdtaf,
        const double thermpressdtam,
        std::map<std::string, Teuchos::RCP<Epetra_Vector>> scatrastatevecs,
        std::map<std::string, Teuchos::RCP<Epetra_MultiVector>> scatrafieldvecs,
        const int scatrandsvel);

    //@}

    //! @name Miscellaneous

    /*!
    \brief Compute a time average of the mean values over all steps
    since the last output. Dump the result to file.

    step on input is used to print the timesteps which belong to the
    statistic to the file

    */

    void TimeAverageMeansAndOutputOfStatistics(const int step);

    /*!
    \brief Compute a time average of the mean values over all steps
    of the sampling period so far. Dump the result to file.

    */

    void DumpStatistics(const int step);

    /*!
    \brief Compute a time average of the mean values for low-Mach-number
    flow over all steps of the sampling period so far. Dump the result to
    file.

    */

    void DumpLomaStatistics(const int step);

    /*!
    \brief Compute a time average of the mean values for turbulent
    passive scalar transprt over all steps of the sampling period so far.
    Dump the result to file.

    */

    void DumpScatraStatistics(const int step);

    /*!
    \brief Reset sums and number of samples to 0

    */

    void ClearStatistics();

    /*!
    \brief Provide the coordinates of the homogeneous planes for a
    turbulent channel flow

    */
    std::vector<double> ReturnNodePlaneCoords() { return (*nodeplanes_); };

    //@}


   protected:
    /*!
    \brief sort criterium for double values up to a tolerance of 10-9

    This is used to create sets of doubles (e.g. coordinates)

    */
    class PlaneSortCriterion
    {
     public:
      bool operator()(const double& p1, const double& p2) const { return (p1 < p2 - 1E-9); }

     protected:
     private:
    };

   private:
    //! direction normal to homogenous plane
    int dim_;

    //! number of elements in sample plane
    int numele_;

    //! number of samples taken
    int numsamp_;

    //! number of records written
    int countrecord_;

    //! flag for physical type of fluid flow (standard: incompressible)
    INPAR::FLUID::PhysicalType physicaltype_;

    //! The discretization (required for nodes, dofs etc;)
    Teuchos::RCP<DRT::Discretization> discret_;

    //! scatra discretization (required for additional multifractal subgrid-scale output in case on
    //! passive scalar transport)
    Teuchos::RCP<DRT::Discretization> scatradiscret_;

    //! flag for ale discretization
    bool alefluid_;

    //! node displacements due to mesh motion
    Teuchos::RCP<Epetra_Vector> dispnp_;

    //! contains plane normal direction etc --- this is the original
    //! fluid dynamic parameterlist
    Teuchos::ParameterList& params_;
    //! name of statistics output file, despite the ending
    const std::string statistics_outfilename_;
    //! contains plane normal direction etc --- this is the original
    //! scatra dynamic parameterlist
    Teuchos::RCP<Teuchos::ParameterList> scatraparams_;
    //! additional parameter list in scatra that contains
    //! required parameters form other sublists such as turbulence models, etc.
    Teuchos::RCP<Teuchos::ParameterList> scatraextraparams_;
    //! additional parameterlist in scatra containing the parameters of the
    //! resepctive time-integration scheme
    Teuchos::RCP<Teuchos::ParameterList> scatratimeparams_;

    //! toggle evaluation of dynamic Smagorinsky/Smagorinsky with
    //! wall damping quantities
    bool smagorinsky_;

    //! toggle evaluation of scale multifractal quantities
    bool multifractal_;

    //! toggle whether to evaluate residuals, taus, dissipation rates etc
    bool subgrid_dissipation_;

    //! boolean indicating turbulent inflow channel discretization
    const bool inflowchannel_;
    //! x-coordinate of outflow of inflow channel
    const double inflowmax_;

    //! parameterlist for the element call when averages of residuals
    //! are calculated --- used for communication between element
    //! and averaging methods --- for fluid field
    Teuchos::ParameterList eleparams_;

    //! parameterlist for the element call when averages of residuals
    //! are calculated --- used for communication between element
    //! and averaging methods --- for scalar field
    Teuchos::ParameterList scatraeleparams_;

    //! pointer to mean vel/pres and scalar field
    Teuchos::RCP<Epetra_Vector> meanvelnp_;
    Teuchos::RCP<Epetra_Vector> meanscanp_;

    //! pointer to vel/pres^2 field (space allocated in constructor)
    Teuchos::RCP<Epetra_Vector> squaredvelnp_;

    //! toogle vectors --- sums are computed by scalarproducts
    //  with these toggle vectors
    Teuchos::RCP<Epetra_Vector> toggleu_;
    Teuchos::RCP<Epetra_Vector> togglev_;
    Teuchos::RCP<Epetra_Vector> togglew_;
    Teuchos::RCP<Epetra_Vector> togglep_;

    //! the dim_-coordinates of the homogeneous planes containing nodes
    Teuchos::RCP<std::vector<double>> nodeplanes_;

    //! the dim_-coordinates of the homogeneous planes --- including
    // additional sampling planes
    Teuchos::RCP<std::vector<double>> planecoordinates_;

    //! a bounding box for the channel
    Teuchos::RCP<Epetra_SerialDenseMatrix> boundingbox_;

    //! viscosity to calculate l_tau, y+ etc.
    double dens_;
    double visc_;
    //! specific heat capacity to calculate Temp_tau (loma)
    double shc_;
    //! Schmidt number (passive scalar)
    double scnum_;

    Teuchos::RCP<FLD::XWall> myxwall_;

    //!--------------------------------------------------
    //!         integration based averaging
    //!--------------------------------------------------
    //
    //! sum over u (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumu_;
    //! sum over v (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumv_;
    //! sum over w (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumw_;
    //! sum over p (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sump_;
    //! sum over density (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumrho_;
    //! sum over T (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumT_;
    //! sum over density*u (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumrhou_;
    //! sum over density*u*T (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumrhouT_;

    //! sum over u^2 (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumsqu_;
    //! sum over v^2 (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumsqv_;
    //! sum over w^2 (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumsqw_;
    //! sum over p^2 (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumsqp_;
    //! sum over density^2 (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumsqrho_;
    //! sum over T^2 (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumsqT_;

    //! sum over uv (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumuv_;
    //! sum over uw (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumuw_;
    //! sum over vw (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumvw_;
    //! sum over uv (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumuT_;
    //! sum over uw (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumvT_;
    //! sum over vw (over one plane in each component)
    Teuchos::RCP<std::vector<double>> sumwT_;

    //!--------------------------------------------------
    //!       the pointwise averaged stuff
    //!--------------------------------------------------
    //
    //! vector of squared velocities and pressures
    Teuchos::RCP<Epetra_Vector> pointsquaredvelnp_;

    //! sum over u (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumu_;
    //! sum over v (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumv_;
    //! sum over w (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumw_;
    //! sum over p (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsump_;
    //! sum over T (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumT_;

    //! sum over u^2 (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumsqu_;
    //! sum over v^2 (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumsqv_;
    //! sum over w^2 (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumsqw_;
    //! sum over p^2 (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumsqp_;
    //! sum over T^2 (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumsqT_;

    //!--------------------------------------------------
    //!   averaged forces (mean, bottom and top)
    //!--------------------------------------------------

    //! sum over nodal forces on boundary in u direction
    double sumforceu_;
    //! sum over nodal forces on boundary in v direction
    double sumforcev_;
    //! sum over nodal forces on boundary in w direction
    double sumforcew_;

    //! sum over nodal forces on boundary in u direction
    double sumforcebu_;
    //! sum over nodal forces on boundary in v direction
    double sumforcebv_;
    //! sum over nodal forces on boundary in w direction
    double sumforcebw_;

    //! sum over nodal forces on boundary in u direction
    double sumforcetu_;
    //! sum over nodal forces on boundary in v direction
    double sumforcetv_;
    //! sum over nodal forces on boundary in w direction
    double sumforcetw_;

    //! heat flux on bottom and top boundary
    double sumqwb_;
    double sumqwt_;

    //!--------------------------------------------------
    //!  averaged quantities from dynamic Smagorinsky
    //!--------------------------------------------------

    //! sum over Cs --- used only for dynamic Smagorinsky model
    Teuchos::RCP<std::vector<double>> sumCs_;
    //! sum over (Cs*delta)^2 --- used only for dynamic Smagorinsky model
    Teuchos::RCP<std::vector<double>> sumCs_delta_sq_;
    //! sum over effective viscosity --- used only for dynamic Smagorinsky model
    Teuchos::RCP<std::vector<double>> sumvisceff_;
    //! increment of sumCs over in one timestep
    Teuchos::RCP<std::vector<double>> incrsumCs_;
    //! increment of sumCs_delta_sq over in one timestep
    Teuchos::RCP<std::vector<double>> incrsumCs_delta_sq_;
    //! increment of sumvisceff over in one timestep
    Teuchos::RCP<std::vector<double>> incrsumvisceff_;
    //! sum over Prt --- used only for dynamic Smagorinsky model
    Teuchos::RCP<std::vector<double>> sumPrt_;
    //! sum over (Cs*delta)^2/Prt --- used only for dynamic Smagorinsky model
    Teuchos::RCP<std::vector<double>> sumCs_delta_sq_Prt_;
    //! sum over effective diffusivity --- used only for dynamic Smagorinsky model
    Teuchos::RCP<std::vector<double>> sumdiffeff_;
    //! increment of sumCs_delta_sq_Prt over in one timestep
    Teuchos::RCP<std::vector<double>> incrsumCs_delta_sq_Prt_;
    //! increment of sumPrt over in one timestep
    Teuchos::RCP<std::vector<double>> incrsumPrt_;
    //! increment of sumdiffeff over in one timestep
    Teuchos::RCP<std::vector<double>> incrsumdiffeff_;
    //! sum over Ci --- used only for dynamic Smagorinsky model for loma
    Teuchos::RCP<std::vector<double>> sumCi_;
    //! sum over (Ci*delta)^2 --- used only for dynamic Smagorinsky model for loma
    Teuchos::RCP<std::vector<double>> sumCi_delta_sq_;
    //! increment of sumCi over in one timestep
    Teuchos::RCP<std::vector<double>> incrsumCi_;
    //! increment of sumCi_delta_sq over in one timestep
    Teuchos::RCP<std::vector<double>> incrsumCi_delta_sq_;

    //!--------------------------------------------------
    //!  averaged quantities from multifractal subgid-scales
    //!--------------------------------------------------

    //! sum over parameter N --- used only for multifractal subgid-scale model
    Teuchos::RCP<std::vector<double>> sumN_stream_;
    Teuchos::RCP<std::vector<double>> sumN_normal_;
    Teuchos::RCP<std::vector<double>> sumN_span_;
    //! increment of parameter N over in one time step
    Teuchos::RCP<std::vector<double>> incrsumN_stream_;
    Teuchos::RCP<std::vector<double>> incrsumN_normal_;
    Teuchos::RCP<std::vector<double>> incrsumN_span_;
    //! sum over parameter B --- used only for multifractal subgid-scale model
    Teuchos::RCP<std::vector<double>> sumB_stream_;
    Teuchos::RCP<std::vector<double>> sumB_normal_;
    Teuchos::RCP<std::vector<double>> sumB_span_;
    //! increment of parameter B over in one time step
    Teuchos::RCP<std::vector<double>> incrsumB_stream_;
    Teuchos::RCP<std::vector<double>> incrsumB_normal_;
    Teuchos::RCP<std::vector<double>> incrsumB_span_;
    //! sum over parameter Csgs --- used only for multifractal subgid-scale model
    Teuchos::RCP<std::vector<double>> sumCsgs_;
    //! increment of Csgs over in one time step
    Teuchos::RCP<std::vector<double>> incrsumCsgs_;
    //! sum over subgrid viscosity --- used only for multifractal subgid-scale model in combination
    //! with eddy viscosity model
    Teuchos::RCP<std::vector<double>> sumsgvisc_;
    //! increment of subgrid viscosity over in one time step
    Teuchos::RCP<std::vector<double>> incrsumsgvisc_;
    //! sum over parameter Nphi --- used only for multifractal subgid-scale model
    Teuchos::RCP<std::vector<double>> sumNphi_;
    //! increment of parameter Nphi over in one time step
    Teuchos::RCP<std::vector<double>> incrsumNphi_;
    //! sum over parameter Dphi --- used only for multifractal subgid-scale model
    Teuchos::RCP<std::vector<double>> sumDphi_;
    //! increment of parameter Dphi over in one time step
    Teuchos::RCP<std::vector<double>> incrsumDphi_;
    //! sum over parameter Csgs_phi --- used only for multifractal subgid-scale model
    Teuchos::RCP<std::vector<double>> sumCsgs_phi_;
    //! increment of Csgs_phi over in one time step
    Teuchos::RCP<std::vector<double>> incrsumCsgs_phi_;

    //!--------------------------------------------------
    //!  averaged resudiuals and subscale quantities
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
    //! sum over all in plane averaged dissipation rates from mfs model (forwardscatter)
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

    //!--------------------------------------------------
    //!  averaged resudiuals and subscale quantities of scalar field
    //!--------------------------------------------------

    //! sum over all in plane stabilization parameters
    Teuchos::RCP<std::vector<double>> sumtauS_;

    //! sum over all in plane residuals of the convection-diffusion equation
    Teuchos::RCP<std::vector<double>> sumresS_;
    //! sum over all in plane squared residuals of the convection-diffusion equation
    Teuchos::RCP<std::vector<double>> sumresS_sq_;

    //! sum over all in plane averaged dissipation rates from supg stabilisation
    Teuchos::RCP<std::vector<double>> sum_scatra_eps_supg_;
    //! sum over all in plane averaged dissipation rates from cross term
    Teuchos::RCP<std::vector<double>> sum_scatra_eps_cross_;
    //! sum over all in plane averaged dissipation rates from reynolds term
    Teuchos::RCP<std::vector<double>> sum_scatra_eps_rey_;
    //! sum over all in plane averaged dissipation rates from eddy viscosity model (Smagorinsky)
    Teuchos::RCP<std::vector<double>> sum_scatra_eps_eddyvisc_;
    //! sum over all in plane averaged dissipation rates from eddy viscosity model (AVM3)
    Teuchos::RCP<std::vector<double>> sum_scatra_eps_avm3_;
    //! sum over all in plane averaged dissipation rates from mfs model
    Teuchos::RCP<std::vector<double>> sum_scatra_eps_mfs_;
    //! sum over all in plane averaged dissipation rates from mfs model (forwardscatter)
    Teuchos::RCP<std::vector<double>> sum_scatra_eps_mfscross_;
    //! sum over all in plane averaged dissipation rates from mfs model (backscatter)
    Teuchos::RCP<std::vector<double>> sum_scatra_eps_mfsrey_;
    //! sum over all in plane averaged dissipation rates from Galerkin viscous term
    Teuchos::RCP<std::vector<double>> sum_scatra_eps_visc_;
    //! sum over all in plane averaged dissipation rates from Galerkin convective term
    Teuchos::RCP<std::vector<double>> sum_scatra_eps_conv_;

    //! number of sampling planes inside the element
    const int numsubdivisions_;
  };

}  // namespace FLD

#endif  // not TURBULENCE_STATISTICS_H
