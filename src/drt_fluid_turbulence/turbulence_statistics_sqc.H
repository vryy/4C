/*----------------------------------------------------------------------*/
/*!

\brief Write (time and space) averaged values to file for
turbulent flow past a square cylinder

o Create sets for various evaluation lines in domain
  (Construction based on a round robin communication pattern):
  - centerline in x1-direction
  - centerline (with respect to cylinder center) in x2-direction
  - lines in wake at x1=7.5 and x1=11.5 in x2-direction
  - lines around cylinder

o loop nodes closest to centerlines

  - generate 4 toggle vectors (u,v,w,p), for example

                            /  1  u dof in homogeneous plane
                 toggleu_  |
                            \  0  elsewhere

  - pointwise multiplication velnp.*velnp for second order
    moments

o values on lines are averaged in time over all steps between two
  outputs

Required parameters are the number of velocity degrees of freedom (3)
and the basename of the statistics outfile. These parameters are
expected to be contained in the fluid time integration parameter list
given on input.

This method is intended to be called every upres_ steps during fluid
output

\maintainer Martin Kronbichler

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef TURBULENCE_STATISTICS_SQC_H
#define TURBULENCE_STATISTICS_SQC_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"


namespace FLD
{
  class TurbulenceStatisticsSqc
  {
   public:
    /*!
    \brief Standard Constructor (public)

        o Create sets for lines in x1- and x2-direction

    o Allocate distributed vector for squares

    */
    TurbulenceStatisticsSqc(Teuchos::RCP<DRT::Discretization> actdis,
        Teuchos::ParameterList& params, const std::string& statistics_outfilename);

    /*!
    \brief Destructor

    */
    virtual ~TurbulenceStatisticsSqc();


    //! @name functions for averaging

    /*!
    \brief The values of lift and drag and its squared values are added.
     This method allows to do the time average after a certain amount of
     timesteps.
    */
    void DoLiftDragTimeSample(double dragforce, double liftforce);

    /*!
    \brief The values of velocity and its squared values are added to
    global vectors. This method allows to do the time average of the
    nodal values after a certain amount of timesteps.
    */
    void DoTimeSample(Teuchos::RCP<Epetra_Vector> velnp);

    /*!
    \brief Dump the result to file.

    step on input is used to print the timesteps which belong to the
    statistic to the file
    */

    void DumpStatistics(int step);

    /*!
    \brief Reset sums and number of samples to 0
    */

    void ClearStatistics();


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
    //! homogeneous direction for sampling
    std::string homdir_;

    //! bounds for extension of cavity in x3-direction
    double x3min_;
    double x3max_;

    //! sums over lift and drag values
    double lift_;
    double drag_;
    double liftsq_;
    double dragsq_;

    //! The discretisation (required for nodes, dofs etc;)
    Teuchos::RCP<DRT::Discretization> discret_;

    //! parameter list
    Teuchos::ParameterList& params_;

    //! name of statistics output file, despite the ending
    const std::string statistics_outfilename_;


    //! pointer to vel/pres^2 field (space allocated in constructor)
    Teuchos::RCP<Epetra_Vector> squaredvelnp_;

    //! toogle vectors: sums are computed by scalarproducts
    Teuchos::RCP<Epetra_Vector> toggleu_;
    Teuchos::RCP<Epetra_Vector> togglev_;
    Teuchos::RCP<Epetra_Vector> togglew_;
    Teuchos::RCP<Epetra_Vector> togglep_;

    //! the coordinates of various lines
    Teuchos::RCP<std::vector<double>> x1ccoordinates_;
    Teuchos::RCP<std::vector<double>> x2ccoordinates_;
    Teuchos::RCP<std::vector<double>> x2wcoordinates_;
    Teuchos::RCP<std::vector<double>> clrcoordinates_;
    Teuchos::RCP<std::vector<double>> ctbcoordinates_;
    //! all coordinates in x1- and x2-direction (required for averaging of Smagorinsky constant)
    Teuchos::RCP<std::vector<double>> x1coordinates_;
    Teuchos::RCP<std::vector<double>> x2coordinates_;

    //! sum over u
    Teuchos::RCP<std::vector<double>> x1csumu_;
    Teuchos::RCP<std::vector<double>> x2csumu_;
    Teuchos::RCP<std::vector<double>> x2w1sumu_;
    Teuchos::RCP<std::vector<double>> x2w2sumu_;
    Teuchos::RCP<std::vector<double>> cyllsumu_;
    Teuchos::RCP<std::vector<double>> cyltsumu_;
    Teuchos::RCP<std::vector<double>> cylrsumu_;
    Teuchos::RCP<std::vector<double>> cylbsumu_;
    //! sum over v
    Teuchos::RCP<std::vector<double>> x1csumv_;
    Teuchos::RCP<std::vector<double>> x2csumv_;
    Teuchos::RCP<std::vector<double>> x2w1sumv_;
    Teuchos::RCP<std::vector<double>> x2w2sumv_;
    Teuchos::RCP<std::vector<double>> cyllsumv_;
    Teuchos::RCP<std::vector<double>> cyltsumv_;
    Teuchos::RCP<std::vector<double>> cylrsumv_;
    Teuchos::RCP<std::vector<double>> cylbsumv_;
    //! sum over w
    Teuchos::RCP<std::vector<double>> x1csumw_;
    Teuchos::RCP<std::vector<double>> x2csumw_;
    Teuchos::RCP<std::vector<double>> x2w1sumw_;
    Teuchos::RCP<std::vector<double>> x2w2sumw_;
    Teuchos::RCP<std::vector<double>> cyllsumw_;
    Teuchos::RCP<std::vector<double>> cyltsumw_;
    Teuchos::RCP<std::vector<double>> cylrsumw_;
    Teuchos::RCP<std::vector<double>> cylbsumw_;
    //! sum over p
    Teuchos::RCP<std::vector<double>> x1csump_;
    Teuchos::RCP<std::vector<double>> x2csump_;
    Teuchos::RCP<std::vector<double>> x2w1sump_;
    Teuchos::RCP<std::vector<double>> x2w2sump_;
    Teuchos::RCP<std::vector<double>> cyllsump_;
    Teuchos::RCP<std::vector<double>> cyltsump_;
    Teuchos::RCP<std::vector<double>> cylrsump_;
    Teuchos::RCP<std::vector<double>> cylbsump_;

    //! sum over u^2
    Teuchos::RCP<std::vector<double>> x1csumsqu_;
    Teuchos::RCP<std::vector<double>> x2csumsqu_;
    Teuchos::RCP<std::vector<double>> x2w1sumsqu_;
    Teuchos::RCP<std::vector<double>> x2w2sumsqu_;
    Teuchos::RCP<std::vector<double>> cyllsumsqu_;
    Teuchos::RCP<std::vector<double>> cyltsumsqu_;
    Teuchos::RCP<std::vector<double>> cylrsumsqu_;
    Teuchos::RCP<std::vector<double>> cylbsumsqu_;
    //! sum over v^2
    Teuchos::RCP<std::vector<double>> x1csumsqv_;
    Teuchos::RCP<std::vector<double>> x2csumsqv_;
    Teuchos::RCP<std::vector<double>> x2w1sumsqv_;
    Teuchos::RCP<std::vector<double>> x2w2sumsqv_;
    Teuchos::RCP<std::vector<double>> cyllsumsqv_;
    Teuchos::RCP<std::vector<double>> cyltsumsqv_;
    Teuchos::RCP<std::vector<double>> cylrsumsqv_;
    Teuchos::RCP<std::vector<double>> cylbsumsqv_;
    //! sum over w^2
    Teuchos::RCP<std::vector<double>> x1csumsqw_;
    Teuchos::RCP<std::vector<double>> x2csumsqw_;
    Teuchos::RCP<std::vector<double>> x2w1sumsqw_;
    Teuchos::RCP<std::vector<double>> x2w2sumsqw_;
    Teuchos::RCP<std::vector<double>> cyllsumsqw_;
    Teuchos::RCP<std::vector<double>> cyltsumsqw_;
    Teuchos::RCP<std::vector<double>> cylrsumsqw_;
    Teuchos::RCP<std::vector<double>> cylbsumsqw_;
    //! sum over uv
    Teuchos::RCP<std::vector<double>> x1csumuv_;
    Teuchos::RCP<std::vector<double>> x2csumuv_;
    Teuchos::RCP<std::vector<double>> x2w1sumuv_;
    Teuchos::RCP<std::vector<double>> x2w2sumuv_;
    Teuchos::RCP<std::vector<double>> cyllsumuv_;
    Teuchos::RCP<std::vector<double>> cyltsumuv_;
    Teuchos::RCP<std::vector<double>> cylrsumuv_;
    Teuchos::RCP<std::vector<double>> cylbsumuv_;
    //! sum over uw
    Teuchos::RCP<std::vector<double>> x1csumuw_;
    Teuchos::RCP<std::vector<double>> x2csumuw_;
    Teuchos::RCP<std::vector<double>> x2w1sumuw_;
    Teuchos::RCP<std::vector<double>> x2w2sumuw_;
    Teuchos::RCP<std::vector<double>> cyllsumuw_;
    Teuchos::RCP<std::vector<double>> cyltsumuw_;
    Teuchos::RCP<std::vector<double>> cylrsumuw_;
    Teuchos::RCP<std::vector<double>> cylbsumuw_;
    //! sum over vw
    Teuchos::RCP<std::vector<double>> x1csumvw_;
    Teuchos::RCP<std::vector<double>> x2csumvw_;
    Teuchos::RCP<std::vector<double>> x2w1sumvw_;
    Teuchos::RCP<std::vector<double>> x2w2sumvw_;
    Teuchos::RCP<std::vector<double>> cyllsumvw_;
    Teuchos::RCP<std::vector<double>> cyltsumvw_;
    Teuchos::RCP<std::vector<double>> cylrsumvw_;
    Teuchos::RCP<std::vector<double>> cylbsumvw_;
    //! sum over p^2
    Teuchos::RCP<std::vector<double>> x1csumsqp_;
    Teuchos::RCP<std::vector<double>> x2csumsqp_;
    Teuchos::RCP<std::vector<double>> x2w1sumsqp_;
    Teuchos::RCP<std::vector<double>> x2w2sumsqp_;
    Teuchos::RCP<std::vector<double>> cyllsumsqp_;
    Teuchos::RCP<std::vector<double>> cyltsumsqp_;
    Teuchos::RCP<std::vector<double>> cylrsumsqp_;
    Teuchos::RCP<std::vector<double>> cylbsumsqp_;
  };

}  // namespace FLD

#endif  // not TURBULENCE_STATISTICS_SQC_H
