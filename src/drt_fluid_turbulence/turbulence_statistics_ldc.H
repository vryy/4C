/*----------------------------------------------------------------------*/
/*!

\brief Write (time and space) averaged values to file.

o Create sets for centerlines in x1-, x2- and x3-direction
  (Construction based on a round robin communication pattern)

o loop nodes closest to centerlines

  - generate 4 toggle vectors (u,v,w,p), for example

                            /  1  u dof in homogeneous plane
                 toggleu_  |
                            \  0  elsewhere

  - pointwise multiplication velnp.*velnp for second order
    moments

o values on centerlines are averaged in time over all steps between two
  outputs

Required parameters are the number of velocity degrees of freedom (3)
and the basename of the statistics outfile. These parameters are
expected to be contained in the fluid time integration parameter list
given on input.

This method is intended to be called every upres_ steps during fluid
output.

\maintainer Martin Kronbichler

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef TURBULENCE_STATISTICS_LDC_H
#define TURBULENCE_STATISTICS_LDC_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_inpar/inpar_fluid.H"

namespace IO
{
  class DiscretizationReader;
  class DiscretizationWriter;
}  // namespace IO

namespace FLD
{
  class TurbulenceStatisticsLdc
  {
   public:
    /*!
    \brief Standard Constructor (public)

        o Create sets for centerlines in x1-, x2- and x3-direction

    o Allocate distributed vector for squares

    */
    TurbulenceStatisticsLdc(Teuchos::RCP<DRT::Discretization> actdis,
        Teuchos::ParameterList& params, const std::string& statistics_outfilename);

    /*!
    \brief Destructor

    */
    virtual ~TurbulenceStatisticsLdc();


    //! @name functions for averaging

    /*!
    \brief The values of velocity, pressure and its squared values are
    added to global vectors. This method allows to do the time average
    of the nodal values after a certain amount of timesteps.
    */
    void DoTimeSample(Teuchos::RCP<Epetra_Vector> velnp);

    /*!
    \brief The values of velocity, pressure, temperature and its squared
    values are added to global vectors. This method allows to do the time
    average of the nodal values after a certain amount of timesteps.
    */
    void DoLomaTimeSample(Teuchos::RCP<Epetra_Vector> velnp, Teuchos::RCP<Epetra_Vector> scanp,
        Epetra_Vector& force, const double eosfac);

    /*!
    \brief Dump the result to file for incompressible flow.

    step on input is used to print the timesteps which belong to the
    statistic to the file
    */

    void DumpStatistics(int step);

    /*!
    \brief Dump the result to file for low-Mach-number flow.

    step on input is used to print the timesteps which belong to the
    statistic to the file
    */

    void DumpLomaStatistics(int step);

    /*!
    \brief Reset sums and number of samples to zero
    */

    void ClearStatistics();


    /*!
    \brief Input of statistics data after restart
    */

    void ReadRestart(IO::DiscretizationReader& reader);

    /*!
    \brief Write output file of statistics data to allow restart
    */
    void WriteRestart(IO::DiscretizationWriter& writer);

   protected:
    /*!
    \brief sort criterium for double values up to a tolerance of 10-9

    This is used to create sets of doubles (e.g. coordinates)

    */
    class LineSortCriterion
    {
     public:
      bool operator()(const double& p1, const double& p2) const { return (p1 < p2 - 1E-9); }

     protected:
     private:
    };

   private:
    //! number of samples taken
    int numsamp_;

    //! bounds for extension of cavity in x1-direction
    double x1min_;
    double x1max_;
    //! bounds for extension of cavity in x2-direction
    double x2min_;
    double x2max_;
    //! bounds for extension of cavity in x3-direction
    double x3min_;
    double x3max_;

    //! The discretisation (required for nodes, dofs etc;)
    Teuchos::RCP<DRT::Discretization> discret_;

    //! parameter list
    Teuchos::ParameterList& params_;

    //! name of statistics output file, despite the ending
    const std::string statistics_outfilename_;

    //! toogle vectors: sums are computed by scalarproducts
    Teuchos::RCP<Epetra_Vector> toggleu_;
    Teuchos::RCP<Epetra_Vector> togglev_;
    Teuchos::RCP<Epetra_Vector> togglew_;
    Teuchos::RCP<Epetra_Vector> togglep_;

    //! the coordinates of the centerlines in x1-, x2- and x3-direction
    Teuchos::RCP<std::vector<double>> x1coordinates_;
    Teuchos::RCP<std::vector<double>> x2coordinates_;
    Teuchos::RCP<std::vector<double>> x3coordinates_;

    //! sum over u (over the centerlines in x1-, x2- and x3-direction)
    Teuchos::RCP<std::vector<double>> x1sumu_;
    Teuchos::RCP<std::vector<double>> x2sumu_;
    Teuchos::RCP<std::vector<double>> x3sumu_;
    //! sum over v (over the centerlines in x1-, x2- and x3-direction)
    Teuchos::RCP<std::vector<double>> x1sumv_;
    Teuchos::RCP<std::vector<double>> x2sumv_;
    Teuchos::RCP<std::vector<double>> x3sumv_;
    //! sum over w (over the centerlines in x1-, x2- and x3-direction)
    Teuchos::RCP<std::vector<double>> x1sumw_;
    Teuchos::RCP<std::vector<double>> x2sumw_;
    Teuchos::RCP<std::vector<double>> x3sumw_;
    //! sum over p (over the centerlines in x1-, x2- and x3-direction)
    Teuchos::RCP<std::vector<double>> x1sump_;
    Teuchos::RCP<std::vector<double>> x2sump_;
    Teuchos::RCP<std::vector<double>> x3sump_;
    //! sum over density (over the centerlines in x1-, x2- and x3-direction)
    Teuchos::RCP<std::vector<double>> x1sumrho_;
    Teuchos::RCP<std::vector<double>> x2sumrho_;
    Teuchos::RCP<std::vector<double>> x3sumrho_;
    //! sum over T (over the centerlines in x1-, x2- and x3-direction)
    Teuchos::RCP<std::vector<double>> x1sumT_;
    Teuchos::RCP<std::vector<double>> x2sumT_;
    Teuchos::RCP<std::vector<double>> x3sumT_;

    //! sum over u^2 (over the centerlines in x1-, x2- and x3-direction)
    Teuchos::RCP<std::vector<double>> x1sumsqu_;
    Teuchos::RCP<std::vector<double>> x2sumsqu_;
    Teuchos::RCP<std::vector<double>> x3sumsqu_;
    //! sum over v^2 (over the centerlines in x1-, x2- and x3-direction)
    Teuchos::RCP<std::vector<double>> x1sumsqv_;
    Teuchos::RCP<std::vector<double>> x2sumsqv_;
    Teuchos::RCP<std::vector<double>> x3sumsqv_;
    //! sum over w^2 (over the centerlines in x1-, x2- and x3-direction)
    Teuchos::RCP<std::vector<double>> x1sumsqw_;
    Teuchos::RCP<std::vector<double>> x2sumsqw_;
    Teuchos::RCP<std::vector<double>> x3sumsqw_;
    //! sum over p^2 (over the centerlines in x1-, x2- and x3-direction)
    Teuchos::RCP<std::vector<double>> x1sumsqp_;
    Teuchos::RCP<std::vector<double>> x2sumsqp_;
    Teuchos::RCP<std::vector<double>> x3sumsqp_;
    //! sum over density^2 (over the centerlines in x1-, x2- and x3-direction)
    Teuchos::RCP<std::vector<double>> x1sumsqrho_;
    Teuchos::RCP<std::vector<double>> x2sumsqrho_;
    Teuchos::RCP<std::vector<double>> x3sumsqrho_;
    //! sum over T^2 (over the centerlines in x1-, x2- and x3-direction)
    Teuchos::RCP<std::vector<double>> x1sumsqT_;
    Teuchos::RCP<std::vector<double>> x2sumsqT_;
    Teuchos::RCP<std::vector<double>> x3sumsqT_;

    //! sum over uv (over the centerlines in x1-, x2- and x3-direction)
    Teuchos::RCP<std::vector<double>> x1sumuv_;
    Teuchos::RCP<std::vector<double>> x2sumuv_;
    Teuchos::RCP<std::vector<double>> x3sumuv_;
    //! sum over uw (over the centerlines in x1-, x2- and x3-direction)
    Teuchos::RCP<std::vector<double>> x1sumuw_;
    Teuchos::RCP<std::vector<double>> x2sumuw_;
    Teuchos::RCP<std::vector<double>> x3sumuw_;
    //! sum over vw (over the centerlines in x1-, x2- and x3-direction)
    Teuchos::RCP<std::vector<double>> x1sumvw_;
    Teuchos::RCP<std::vector<double>> x2sumvw_;
    Teuchos::RCP<std::vector<double>> x3sumvw_;
    //! sum over uT (over the centerlines in x1-, x2- and x3-direction)
    Teuchos::RCP<std::vector<double>> x1sumuT_;
    Teuchos::RCP<std::vector<double>> x2sumuT_;
    Teuchos::RCP<std::vector<double>> x3sumuT_;
    //! sum over vT (over the centerlines in x1-, x2- and x3-direction)
    Teuchos::RCP<std::vector<double>> x1sumvT_;
    Teuchos::RCP<std::vector<double>> x2sumvT_;
    Teuchos::RCP<std::vector<double>> x3sumvT_;
    //! sum over wT (over the centerlines in x1-, x2- and x3-direction)
    Teuchos::RCP<std::vector<double>> x1sumwT_;
    Teuchos::RCP<std::vector<double>> x2sumwT_;
    Teuchos::RCP<std::vector<double>> x3sumwT_;
  };

}  // namespace FLD

#endif  // not TURBULENCE_STATISTICS_LDC_H
