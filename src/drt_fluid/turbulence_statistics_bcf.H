/*!----------------------------------------------------------------------
\file turbulence_statistics_cha_mult_phase.H

\brief Compute (time and space) averaged values for turbulent channel
       flows with multiple phases and write them to files.

<pre>

This class has the essential features of the regular channel
flow. All averages are calculated separately for each phase.
In addition the fraction distribution is calculated.

o Create set of all available homogeneous planes
  (Construction based on a round robin communication pattern)

o loop planes (e.g. plane coordinates)
  - pointwise in-plane average of first- and second order moments
  - in-plane average of wall force

o in plane mean values are averaged in time over all steps between two
  outputs (by computation of the arithmetic mean)
  - time average pointwise values
  - time average forces

o Write pointwise and integral statistics for first and second
  order moments
  ->   .flow_statistic

Required parameters are the number of velocity degrees of freedom (3),
the normal direction to the plane, in which the average values in space
should be computed, and the basename of the statistics outfile. These
parameters are expected to be contained in the fluid time integration
parameter list given on input.

In addition the multiple phase version needs information about which
phase the velocity dofs belong to.

Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#ifndef TURBULENCE_STATISTICS_BCF_H
#define TURBULENCE_STATISTICS_BCF_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
//#include "turbulence_statistics_cha.H"
#include "../linalg/linalg_utils.H"
#include "../drt_inpar/inpar_fluid.H"

#include "../drt_mat/newtonianfluid.H"

using namespace std;
using namespace Teuchos;

//forward declaration
namespace DRT
{
  class DofSet;
}

namespace COMBUST
{

class TurbulenceStatisticsBcf
{

public:
  /*!
  \brief Standard Constructor (public)

  <pre>
  o Create vector of homogeneous plane coordinates

  o allocate all vectors of sum_something vectors

  o initialize the output (open/clear files, print header)

  </pre>

  */
  TurbulenceStatisticsBcf(
                          RCP<DRT::Discretization>      actdis,
                          Teuchos::ParameterList&                        params);


  /*!
  \brief Destructor

  */
  virtual ~TurbulenceStatisticsBcf();


  //! @name functions for (spatial) averaging

  /*!
  \brief Compute the in-plane mean values of first and second order
  moments for velocities and pressure are added to global
  'sum' vectors.

  \param stdvelnp             (in) : The current stdvelnp vector created with a stddofset,
                                     including only velocity and pressure dofs.
  \param stdforce             (in) : The current force or true residuum vector created
                                      with the same dofset as stdvelnp.
  \param stdforce             (in) : The dofset used to create the above mentioned vectors.
  \param discretmatchingvelnp (in) : The current velnp vector created with the dofset
                                     as returned by the discretization discret_.
  \param phinp                (in) : The phi vector based on the fluids node row map. This
                                     is a hack!
  */
  void DoTimeSample(
      Teuchos::RCP<const Epetra_Vector> stdvelnp                            ,
      Teuchos::RCP<const Epetra_Vector> stdforce                            ,
      Teuchos::RCP<const DRT::DofSet>   stddofset            = Teuchos::null,
      Teuchos::RCP<const Epetra_Vector> discretmatchingvelnp = Teuchos::null,
      Teuchos::RCP<const Epetra_Vector> phinp                = Teuchos::null);


  /*!
  \brief Compute in plane means of u,u^2 etc. (integral version)

  The averages here are calculated by integration.

  The calculated value is added to the sum**,sumsq** variables in the
  component corresponding to the plane.

  Further documentation is provided in the element subroutines
  */
  void EvaluateIntegralMeanValuesInPlanes();

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
  */
  void EvaluatePointwiseMeanValuesInPlanes();

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
  //====================================================
  // Members which are used only during the DoTimeSample
  //====================================================

  //! DofSet for the result Epetra vectors
  RCP<const DRT::DofSet>             stddofset_;

  //! pointer to mean vel/pres and scalar field
  RCP<const Epetra_Vector>           stdvelnp_;

  //! vector of squared velocities and pressures
  RCP<Epetra_Vector>                 pointsquaredvelnp_;

  //! pointer to the full velnp. Must match the discretization
  RCP<const Epetra_Vector>           fullvelnp_;

  //! pointer to the phinp vector. This is a dirty hack!!!
  RCP<const Epetra_Vector>           phinp_;

  //! toogle vectors --- sums are computed by scalarproducts
  //  with these toggle vectors
  std::vector<RCP<Epetra_Vector> >   toggleu_;
  std::vector<RCP<Epetra_Vector> >   togglev_;
  std::vector<RCP<Epetra_Vector> >   togglew_;
  std::vector<RCP<Epetra_Vector> >   togglep_;

  //! parameterlist for the element call when averages of residuals
  //! are calculated --- used for communication between element
  //! and averaging methods
  ParameterList                              eleparams_;

  //====================================================
  // Permanent members
  //====================================================

  //! direction normal to homogenous plane
  int                                        dim_;

  //! number of phases
  size_t                                     numphase_;

  //! number of elements in sample plane
  int                                        numele_;

  //! number of nodes in sample plane
  int                                        numnode_;

  //! number of samples taken
  int                                        numsamp_;

  //! number of records written
  int                                        countrecord_;

  /// flag for physical type of fluid flow (standard: incompressible) // only one required? can you combine different types in multiphase?
  INPAR::FLUID::PhysicalType                 physicaltype_;

  //! The discretisation (required for nodes, etc; dofs are in a separate dofset)
  RCP<DRT::Discretization>           discret_;

  //! contains plane normal direction etc --- this is the original
  //! fluid dynamic parameterlist
  Teuchos::ParameterList&                             params_;

  //! the dim_-coordinates of the homogeneous planes containing nodes
  RCP<std::vector<double> >               nodeplanes_;

  //! a bounding box for the channel
  RCP<Epetra_SerialDenseMatrix>      boundingbox_;

  //! viscosity to calculate l_tau, y+ etc.
  std::vector<double>                        visc_;

  //! maps the global material ids to the position in
  //  the vectors which belongs to that material
  RCP<map<int,int> >                 matidtoindex_;

  //!--------------------------------------------------
  //!         integration based averaging
  //!--------------------------------------------------

  //! sum over the volume of each phase (in one plane)
  std::vector<RCP<std::vector<double> > > sumvol_;


  //!--------------------------------------------------
  //!       the pointwise averaged stuff
  //!--------------------------------------------------
  //
  //! pointer to the sum of nodes per plane
  std::vector<RCP<std::vector<double> > > pointsumnode_;

  //! sum over u (over one plane in each component)
  std::vector<RCP<std::vector<double> > > pointsumu_;
  //! sum over v (over one plane in each component)
  std::vector<RCP<std::vector<double> > > pointsumv_;
  //! sum over w (over one plane in each component)
  std::vector<RCP<std::vector<double> > > pointsumw_;
  //! sum over p (over one plane in each component)
  std::vector<RCP<std::vector<double> > > pointsump_;

  //! sum over u (over one plane in each component)
  std::vector<RCP<std::vector<double> > > pointsumuv_;
  //! sum over v (over one plane in each component)
  std::vector<RCP<std::vector<double> > > pointsumuw_;
  //! sum over w (over one plane in each component)
  std::vector<RCP<std::vector<double> > > pointsumvw_;

  //! sum over u^2 (over one plane in each component)
  std::vector<RCP<std::vector<double> > > pointsumsqu_;
  //! sum over v^2 (over one plane in each component)
  std::vector<RCP<std::vector<double> > > pointsumsqv_;
  //! sum over w^2 (over one plane in each component)
  std::vector<RCP<std::vector<double> > > pointsumsqw_;
  //! sum over p^2 (over one plane in each component)
  std::vector<RCP<std::vector<double> > > pointsumsqp_;

  //!--------------------------------------------------
  //!   averaged forces (mean, bottom and top)
  //!--------------------------------------------------

  //! sum over nodal forces on boundary in u direction
  double                                     sumforceu_;
  //! sum over nodal forces on boundary in v direction
  double                                     sumforcev_;
  //! sum over nodal forces on boundary in w direction
  double                                     sumforcew_;

};

}

#endif  // not TURBULENCE_STATISTICS_BCF_H

