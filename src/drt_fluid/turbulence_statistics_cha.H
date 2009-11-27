/*!----------------------------------------------------------------------
\file turbulence_statistics_cha.H

\brief Compute (time and space) averaged values for turbulent channel
       flows and write them to files.

<pre>

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

Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef TURBULENCE_STATISTICS_H
#define TURBULENCE_STATISTICS_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_control_point.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_inpar/inpar_fluid.H"

#include "../drt_mat/newtonianfluid.H"

using namespace std;
using namespace Teuchos;

namespace FLD
{

class TurbulenceStatisticsCha
{

public:
  /*!
  \brief Standard Constructor (public)

  <pre>
  o Create vector of homogeneous plane coordinates

  o Allocate 4 distributed toggle vectors and one distributed vector
    for squares, both for integral and pointwise means

  o allocate all sum_something vectors

  o initialise the output (open/clear files, print header)

  </pre>

  */
  TurbulenceStatisticsCha(RefCountPtr<DRT::Discretization> actdis             ,
                          bool                             alefluid           ,
                          RefCountPtr<Epetra_Vector>       dispnp             ,
                          ParameterList&                   params             ,
                          bool                             smagorinsky        ,
                          bool                             subgrid_dissipation);

  /*!
  \brief Destructor

  */
  virtual ~TurbulenceStatisticsCha();


  //! @name functions for (spatial) averaging

  /*!
  \brief Compute the in-plane mean values of first and second order
  moments for velocities, pressure and Cs are added to global
  'sum' vectors.
  */
  void DoTimeSample(Teuchos::RefCountPtr<Epetra_Vector> velnp,
                    Epetra_Vector & force);


  /*!
  \brief The values of velocity, pressure, temperature and its squared
  values are added to global vectors. This method allows to do the time
  average of the nodal values after a certain amount of timesteps.
  */
  void DoLomaTimeSample(Teuchos::RefCountPtr<Epetra_Vector> velnp,
                        Teuchos::RefCountPtr<Epetra_Vector> scanp,
                        Epetra_Vector &                     force,
                        const double                        eosfac);


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
    \brief do averaging of residuals, dissipation rates etc
           (all gausspoint-quantities)

  */
  void EvaluateResiduals(map<string,RCP<Epetra_Vector> > statevecs,
                         double                          time     );


  //@}

  //! @name Miscellaneous

  /*!
  \brief Compute a time average of the mean values over all steps
  since the last output. Dump the result to file.

  step on input is used to print the timesteps which belong to the
  statistic to the file

  */

  void TimeAverageMeansAndOutputOfStatistics(int step);

  /*!
  \brief Compute a time average of the mean values over all steps
  of the sampling period so far. Dump the result to file.

  */

  void DumpStatistics(int step);

  /*!
  \brief Compute a time average of the mean values for low-Mach-number
  flow over all steps of the sampling period so far. Dump the result to
  file.

  */

  void DumpLomaStatistics(int step);

  /*!
  \brief Reset sums and number of samples to 0

  */

  void ClearStatistics();

  /*!
  \brief Provide the coordinates of the homogeneous planes for a
  turbulent channel flow

  */
  vector<double> ReturnNodePlaneCoords()
    {
      return(*nodeplanes_);
    };

  //@}


protected:

  /*!
  \brief sort criterium for double values up to a tolerance of 10-9

  This is used to create sets of doubles (e.g. coordinates)

  */
  class PlaneSortCriterion
  {
  public:
    bool operator() (
      const double& p1,
      const double& p2) const {
      return (p1 < p2 - 1E-9);
    }

  protected:
  private:

  };

private:
  //! direction normal to homogenous plane
  int                              dim_;

  //! number of elements in sample plane
  int                              numele_;

  //! number of samples taken
  int                              numsamp_;

  //! number of records written
  int                              countrecord_;

  /// flag for physical type of fluid flow (standard: incompressible)
  INPAR::FLUID::PhysicalType 		physicaltype_;

  //! The discretisation (required for nodes, dofs etc;)
  RefCountPtr<DRT::Discretization> discret_;

  //! flag for ale discretisation
  bool                             alefluid_;

  //! node displacements due to mesh motion
  RefCountPtr<Epetra_Vector>       dispnp_;

  //! contains plane normal direction etc --- this is the original
  //! fluid dynamic parameterlist
  ParameterList&                   params_;

  //! toggle evaluation of dynamic Smagorinsky/Smagorinsky with
  //! wall damping quantities
  bool                             smagorinsky_;

  //! toggle whether to evaluate residuals, taus, dissipation rates etc
  bool                             subgrid_dissipation_;

  //! parameterlist for the element call when averages of residuals
  //! are calculated --- used for communication between element
  //! and averaging methods
  ParameterList                    eleparams_;

  //! pointer to mean vel/pres and scalar field
  RefCountPtr<Epetra_Vector>       meanvelnp_;
  RefCountPtr<Epetra_Vector>       meanscanp_;

   //! pointer to vel/pres^2 field (space allocated in constructor)
  RefCountPtr<Epetra_Vector>       squaredvelnp_;

  //! toogle vectors --- sums are computed by scalarproducts
  //  with these toggle vectors
  RefCountPtr<Epetra_Vector>       toggleu_;
  RefCountPtr<Epetra_Vector>       togglev_;
  RefCountPtr<Epetra_Vector>       togglew_;
  RefCountPtr<Epetra_Vector>       togglep_;

  //! the dim_-coordinates of the homogeneous planes containing nodes
  RefCountPtr<vector<double> >     nodeplanes_;

  //! the dim_-coordinates of the homogeneous planes --- including
  //additional sampling planes
  RefCountPtr<vector<double> >     planecoordinates_;

  //! a bounding box for the channel
  RefCountPtr<Epetra_SerialDenseMatrix> boundingbox_;

  //! viscosity to calculate l_tau, y+ etc.
  double                           visc_;

  //!--------------------------------------------------
  //!         integration based averaging
  //!--------------------------------------------------
  //
  //! sum over u (over one plane in each component)
  RefCountPtr<vector<double> >     sumu_;
  //! sum over v (over one plane in each component)
  RefCountPtr<vector<double> >     sumv_;
  //! sum over w (over one plane in each component)
  RefCountPtr<vector<double> >     sumw_;
  //! sum over p (over one plane in each component)
  RefCountPtr<vector<double> >     sump_;
  //! sum over density (over one plane in each component)
  RefCountPtr<vector<double> >     sumrho_;
  //! sum over T (over one plane in each component)
  RefCountPtr<vector<double> >     sumT_;
  //! sum over density*u (over one plane in each component)
  RefCountPtr<vector<double> >     sumrhou_;
  //! sum over density*u*T (over one plane in each component)
  RefCountPtr<vector<double> >     sumrhouT_;

  //! sum over u^2 (over one plane in each component)
  RefCountPtr<vector<double> >     sumsqu_;
  //! sum over v^2 (over one plane in each component)
  RefCountPtr<vector<double> >     sumsqv_;
  //! sum over w^2 (over one plane in each component)
  RefCountPtr<vector<double> >     sumsqw_;
  //! sum over p^2 (over one plane in each component)
  RefCountPtr<vector<double> >     sumsqp_;
  //! sum over density^2 (over one plane in each component)
  RefCountPtr<vector<double> >     sumsqrho_;
  //! sum over T^2 (over one plane in each component)
  RefCountPtr<vector<double> >     sumsqT_;

  //! sum over uv (over one plane in each component)
  RefCountPtr<vector<double> >     sumuv_;
  //! sum over uw (over one plane in each component)
  RefCountPtr<vector<double> >     sumuw_;
  //! sum over vw (over one plane in each component)
  RefCountPtr<vector<double> >     sumvw_;
  //! sum over uv (over one plane in each component)
  RefCountPtr<vector<double> >     sumuT_;
  //! sum over uw (over one plane in each component)
  RefCountPtr<vector<double> >     sumvT_;
  //! sum over vw (over one plane in each component)
  RefCountPtr<vector<double> >     sumwT_;

  //!--------------------------------------------------
  //!       the pointwise averaged stuff
  //!--------------------------------------------------
  //
  //! vector of squared velocities and pressures
  RefCountPtr<Epetra_Vector>       pointsquaredvelnp_;

  //! sum over u (over one plane in each component)
  RefCountPtr<vector<double> >     pointsumu_;
  //! sum over v (over one plane in each component)
  RefCountPtr<vector<double> >     pointsumv_;
  //! sum over w (over one plane in each component)
  RefCountPtr<vector<double> >     pointsumw_;
  //! sum over p (over one plane in each component)
  RefCountPtr<vector<double> >     pointsump_;
  //! sum over T (over one plane in each component)
  RefCountPtr<vector<double> >     pointsumT_;

  //! sum over u^2 (over one plane in each component)
  RefCountPtr<vector<double> >     pointsumsqu_;
  //! sum over v^2 (over one plane in each component)
  RefCountPtr<vector<double> >     pointsumsqv_;
  //! sum over w^2 (over one plane in each component)
  RefCountPtr<vector<double> >     pointsumsqw_;
  //! sum over p^2 (over one plane in each component)
  RefCountPtr<vector<double> >     pointsumsqp_;
  //! sum over T^2 (over one plane in each component)
  RefCountPtr<vector<double> >     pointsumsqT_;

  //!--------------------------------------------------
  //!   averaged forces (mean, bottom and top)
  //!--------------------------------------------------

  //! sum over nodal forces on boundary in u direction
  double                           sumforceu_;
  //! sum over nodal forces on boundary in v direction
  double                           sumforcev_;
  //! sum over nodal forces on boundary in w direction
  double                           sumforcew_;

  //! sum over nodal forces on boundary in u direction
  double                           sumforcebu_;
  //! sum over nodal forces on boundary in v direction
  double                           sumforcebv_;
  //! sum over nodal forces on boundary in w direction
  double                           sumforcebw_;

  //! sum over nodal forces on boundary in u direction
  double                           sumforcetu_;
  //! sum over nodal forces on boundary in v direction
  double                           sumforcetv_;
  //! sum over nodal forces on boundary in w direction
  double                           sumforcetw_;

  //! heat flux on bottom and top boundary
  double                           sumqwb_;
  double                           sumqwt_;

  //!--------------------------------------------------
  //!  averaged quantities from dynamic Smagorinsky
  //!--------------------------------------------------

  //! sum over Cs --- used only for dynamic Smagorinsky model
  RefCountPtr<vector<double> >     sumCs_;
  //! sum over (Cs*delta)^2 --- used only for dynamic Smagorinsky model
  RefCountPtr<vector<double> >     sumCs_delta_sq_;
  //! sum over effective viscosity --- used only for dynamic Smagorinsky model
  RefCountPtr<vector<double> >     sumvisceff_;
  //! increment of sumCs over in one timestep
  RefCountPtr<vector<double> >     incrsumCs_;
  //! increment of sumCs_delta_sq over in one timestep
  RefCountPtr<vector<double> >     incrsumCs_delta_sq_;
  //! increment of sumvisceff over in one timestep
  RefCountPtr<vector<double> >     incrsumvisceff_;

  //!--------------------------------------------------
  //!  averaged resudiuals and subscale quantities
  //!--------------------------------------------------

  //! sum over all in plane element sizes
  RefCountPtr<vector<double> >     sumhk_;
  //! sum over all in plane element sizes for the Bazilevs
  //! parameter in the viscous regime
  RefCountPtr<vector<double> >     sumhbazilevs_;
  //! sum over all in plane stream lengths
  RefCountPtr<vector<double> >     sumstrle_;
  //! sum over all in plane stream lengths
  RefCountPtr<vector<double> >     sumgradle_;

  //! sum over all in plane residuals
  RefCountPtr<vector<double> >     sumtauM_;
  //! sum over all in plane squared residuals
  RefCountPtr<vector<double> >     sumtauC_;

  //! sum over all in plane residuals
  RefCountPtr<vector<double> >     sumres_;
  //! sum over all in plane squared residuals
  RefCountPtr<vector<double> >     sumres_sq_;
  //! sum over all in plane residuals norms
  RefCountPtr<vector<double> >     sumabsres_;
  //! sum over all in plane values of svel/tau
  RefCountPtr<vector<double> >     sumtauinvsvel_;
  //! sum over all in plane subscale accelerations
  RefCountPtr<vector<double> >     sumsacc_;
  //! sum over all in plane squared subscale accelerations
  RefCountPtr<vector<double> >     sumsacc_sq_;
  //! sum over all in plane subscale accelerations norms
  RefCountPtr<vector<double> >     sumabssacc_;
  //! sum over all in plane subscale velocities
  RefCountPtr<vector<double> >     sumsvelaf_;
  //! sum over all in plane squared subscale velocities
  RefCountPtr<vector<double> >     sumsvelaf_sq_;
  //! sum over all in plane subscale velocities norms
  RefCountPtr<vector<double> >     sumabssvelaf_;

  //! sum over all in plane residuals of the continuity equation
  RefCountPtr<vector<double> >     sumresC_;
  //! sum over all in plane squared residuals of the continuity equation
  RefCountPtr<vector<double> >     sumresC_sq_;
  //! sum over all in plane subscale pressure values at current timestep
  RefCountPtr<vector<double> >     sumspressnp_   ;
  //! sum over all in plane squared subscale pressure values at current timestep
  RefCountPtr<vector<double> >     sumspressnp_sq_;

  //! sum over all in plane averaged dissipation rates from subscale acceleration
  RefCountPtr<vector<double> >     sum_eps_sacc_;
  //! sum over all in plane averaged dissipation rates from pspg stabilisation
  RefCountPtr<vector<double> >     sum_eps_pspg_;
  //! sum over all in plane averaged dissipation rates from supg stabilisation
  RefCountPtr<vector<double> >     sum_eps_supg_;
  //! sum over all in plane averaged dissipation rates from cross term
  RefCountPtr<vector<double> >     sum_eps_cross_;
  //! sum over all in plane averaged dissipation rates from reynolds term
  RefCountPtr<vector<double> >     sum_eps_rey_;
  //! sum over all in plane averaged dissipation rates from least squares continuity term
  RefCountPtr<vector<double> >     sum_eps_cstab_;
  //! sum over all in plane averaged dissipation rates from viscous stabilisation
  RefCountPtr<vector<double> >     sum_eps_vstab_;
  //! sum over all in plane averaged dissipation rates from eddy viscosity model (Smagorinsky)
  RefCountPtr<vector<double> >     sum_eps_eddyvisc_;
  //! sum over all in plane averaged dissipation rates from Galerkin viscous term
  RefCountPtr<vector<double> >     sum_eps_visc_;
  //! sum over all in plane averaged dissipation rates from Galerkin convective term
  RefCountPtr<vector<double> >     sum_eps_conv_;

  //! sum over all in plane averaged supg+cross stress
  RefCountPtr<vector<double> >     sum_crossstress_;
  //! sum over all in plane averaged reynolds stress
  RefCountPtr<vector<double> >     sum_reystress_;

};

}

#endif  // not TURBULENCE_STATISTICS_H

#endif /* CCADISCRET       */
