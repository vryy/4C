/*----------------------------------------------------------------------*/
/*! \file

\brief Compute (time and space) averaged values for turbulent flows
       around a rotating cylinder and write them to files.

o Create set of all available homogeneous shells
  (Construction based on a round robin communication pattern)

o loop shells (e.g. radial shell coordinates)

  - pointwise in-plane average of first- and second order moments

o Write pointwise statistics for first and second order moments
  ->   .flow_statistic

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

#ifndef TURBULENCE_STATISTICS_CCY_H
#define TURBULENCE_STATISTICS_CCY_H


#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseMatrix.h"



namespace DRT
{
  class Discretization;
}

namespace FLD
{
  class TurbulenceStatisticsCcy
  {
   public:
    /*!
    \brief Standard Constructor (public)

        o Create vector of radial coordinates of homogeneous shells

    o allocate all sum_something vectors

    o initialise the output (open/clear files, print header)


    */
    TurbulenceStatisticsCcy(Teuchos::RCP<DRT::Discretization> actdis, bool alefluid,
        Teuchos::RCP<Epetra_Vector> dispnp, Teuchos::ParameterList& params,
        const std::string& statistics_outfilename, const bool withscatra);

    /*!
    \brief Destructor

    */
    virtual ~TurbulenceStatisticsCcy();


    //! @name functions for (spatial) averaging

    /*!
    \brief Compute the in-shell mean values of first and second order
    moments for velocities, pressure (and transported scalar fields).
    */
    void DoTimeSample(Teuchos::RCP<Epetra_Vector> velnp, Teuchos::RCP<Epetra_Vector> scanp,
        Teuchos::RCP<Epetra_Vector> fullphinp);


    /*!
    \brief Compute in plane means of u,u^2 etc. (nodal quantities)

    The averages here are calculated as the arithmetic mean of
    point values (computed by interpolation)

    The calculated values are added to the pointsum**,pointsumsq** variables
    in the component corresponding to the plane.

    velnp is the solution vector provided by the time integration
    algorithm
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
    \brief Reset sums and number of samples to 0

    */

    void ClearStatistics();

    /*!
    \brief Provide the radius of the homogeneous shell for a
    flow around a rotating circular cylinder

    */
    std::vector<double> ReturnShellPlaneRadius() { return (*nodeshells_); };

    //@}

    // Add results from scalar transport field solver to statistics
    void AddScaTraResults(
        Teuchos::RCP<DRT::Discretization> scatradis, Teuchos::RCP<Epetra_Vector> phinp);

   protected:
    /*!
    \brief sort criterium for double values up to a tolerance of 10-6

    This is used to create sets of doubles (e.g. coordinates)

    */
    class PlaneSortCriterion
    {
     public:
      bool operator()(const double& p1, const double& p2) const { return (p1 < p2 - 1E-6); }

     protected:
     private:
    };

   private:
    //! direction normal to homogenous plane
    int dim_;

    //! number of samples taken
    int numsamp_;

    //! number of records written
    int countrecord_;

    //! The discretisation (required for nodes, dofs etc;)
    Teuchos::RCP<DRT::Discretization> discret_;

    //! the scatra discretization
    Teuchos::RCP<DRT::Discretization> scatradis_;

    //! node displacements due to mesh motion
    Teuchos::RCP<Epetra_Vector> dispnp_;

    //! contains plane normal direction etc --- this is the original
    //! fluid dynamic parameterlist
    Teuchos::ParameterList& params_;

    //! name of statistics output file, despite the ending
    const std::string statistics_outfilename_;


    //! parameterlist for the element call when averages of residuals
    //! are calculated --- used for communication between element
    //! and averaging methods
    Teuchos::ParameterList eleparams_;

    //! pointer to mean vel/pres field
    Teuchos::RCP<Epetra_Vector> meanvelnp_;

    //! the dim_-coordinates of the homogeneous planes containing nodes
    Teuchos::RCP<std::vector<double>> nodeshells_;

    //! the dim_-coordinates of the homogeneous planes --- including
    // additional sampling planes
    Teuchos::RCP<std::vector<double>> shellcoordinates_;

    //!--------------------------------------------------
    //!       the pointwise averaged stuff
    //!--------------------------------------------------
    //

    //! sum over u (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumu_;
    //! sum over v (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumv_;
    //! sum over w (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumw_;
    //! sum over p (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsump_;

    //! sum over u^2 (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumuu_;
    //! sum over v^2 (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumvv_;
    //! sum over w^2 (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumww_;
    //! sum over p^2 (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumpp_;

    //! sum over uv (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumuv_;
    //! sum over uw (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumuw_;
    //! sum over vw (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumvw_;

    // for additional scalar transport

    //! flag for mass transport statistics
    const bool withscatra_;

    //! number of scatra dofs per node
    int numscatradofpernode_;

    //! pointer to mean scalar field
    Teuchos::RCP<Epetra_Vector> meanscanp_;

    //! pointer to mean field of all scatra results
    Teuchos::RCP<Epetra_Vector> meanfullphinp_;

    //! sum over c (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumc_;
    //! sum over c^2 (over one plane in each component)
    Teuchos::RCP<std::vector<double>> pointsumcc_;

    //! sum over c (over one plane in each component)
    Teuchos::RCP<Epetra_SerialDenseMatrix> pointsumphi_;
    //! sum over c^2 (over one plane in each component)
    Teuchos::RCP<Epetra_SerialDenseMatrix> pointsumphiphi_;
  };

}  // namespace FLD

#endif  // not TURBULENCE_STATISTICS_CCY_H
