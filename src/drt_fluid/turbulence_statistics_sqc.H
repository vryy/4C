/*!----------------------------------------------------------------------
\file turbulence_statistics_sqc.H

\brief Write (time and space) averaged values to file for
turbulent flow past a square cylinder


<pre>
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
output.

Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/

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


using namespace std;
using namespace Teuchos;

namespace FLD
{

class TurbulenceStatisticsSqc
{

public:
  /*!
  \brief Standard Constructor (public)

  <pre>
  o Create sets for lines in x1- and x2-direction

  o Allocate distributed vector for squares
  </pre>

  */
  TurbulenceStatisticsSqc(RCP<DRT::Discretization> actdis,
                          ParameterList&                   params);

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
    bool operator() (
      const double& p1,
      const double& p2) const {
      return (p1 < p2 - 1E-9);
    }

  protected:
  private:

  };

private:

  //! number of samples taken
  int                              numsamp_;
  //! homogeneous direction for sampling
  string                           homdir_;

  //! bounds for extension of cavity in x3-direction
  double                           x3min_;
  double                           x3max_;

  //! sums over lift and drag values
  double                           lift_;
  double                           drag_;
  double                           liftsq_;
  double                           dragsq_;

  //! The discretisation (required for nodes, dofs etc;)
  RCP<DRT::Discretization> discret_;

  //! parameter list
  ParameterList&                   params_;

   //! pointer to vel/pres^2 field (space allocated in constructor)
  RCP<Epetra_Vector>       squaredvelnp_;

  //! toogle vectors: sums are computed by scalarproducts
  RCP<Epetra_Vector>       toggleu_;
  RCP<Epetra_Vector>       togglev_;
  RCP<Epetra_Vector>       togglew_;
  RCP<Epetra_Vector>       togglep_;

  //! the coordinates of various lines
  RCP<vector<double> >     x1ccoordinates_;
  RCP<vector<double> >     x2ccoordinates_;
  RCP<vector<double> >     x2wcoordinates_;
  RCP<vector<double> >     clrcoordinates_;
  RCP<vector<double> >     ctbcoordinates_;
  //! all coordinates in x1- and x2-direction (required for averaging of Smagorinsky constant)
  RCP<vector<double> >     x1coordinates_;
  RCP<vector<double> >     x2coordinates_;

  //! sum over u
  RCP<vector<double> >     x1csumu_;
  RCP<vector<double> >     x2csumu_;
  RCP<vector<double> >     x2w1sumu_;
  RCP<vector<double> >     x2w2sumu_;
  RCP<vector<double> >     cyllsumu_;
  RCP<vector<double> >     cyltsumu_;
  RCP<vector<double> >     cylrsumu_;
  RCP<vector<double> >     cylbsumu_;
  //! sum over v
  RCP<vector<double> >     x1csumv_;
  RCP<vector<double> >     x2csumv_;
  RCP<vector<double> >     x2w1sumv_;
  RCP<vector<double> >     x2w2sumv_;
  RCP<vector<double> >     cyllsumv_;
  RCP<vector<double> >     cyltsumv_;
  RCP<vector<double> >     cylrsumv_;
  RCP<vector<double> >     cylbsumv_;
  //! sum over w
  RCP<vector<double> >     x1csumw_;
  RCP<vector<double> >     x2csumw_;
  RCP<vector<double> >     x2w1sumw_;
  RCP<vector<double> >     x2w2sumw_;
  RCP<vector<double> >     cyllsumw_;
  RCP<vector<double> >     cyltsumw_;
  RCP<vector<double> >     cylrsumw_;
  RCP<vector<double> >     cylbsumw_;
  //! sum over p
  RCP<vector<double> >     x1csump_;
  RCP<vector<double> >     x2csump_;
  RCP<vector<double> >     x2w1sump_;
  RCP<vector<double> >     x2w2sump_;
  RCP<vector<double> >     cyllsump_;
  RCP<vector<double> >     cyltsump_;
  RCP<vector<double> >     cylrsump_;
  RCP<vector<double> >     cylbsump_;

  //! sum over u^2
  RCP<vector<double> >     x1csumsqu_;
  RCP<vector<double> >     x2csumsqu_;
  RCP<vector<double> >     x2w1sumsqu_;
  RCP<vector<double> >     x2w2sumsqu_;
  RCP<vector<double> >     cyllsumsqu_;
  RCP<vector<double> >     cyltsumsqu_;
  RCP<vector<double> >     cylrsumsqu_;
  RCP<vector<double> >     cylbsumsqu_;
  //! sum over v^2
  RCP<vector<double> >     x1csumsqv_;
  RCP<vector<double> >     x2csumsqv_;
  RCP<vector<double> >     x2w1sumsqv_;
  RCP<vector<double> >     x2w2sumsqv_;
  RCP<vector<double> >     cyllsumsqv_;
  RCP<vector<double> >     cyltsumsqv_;
  RCP<vector<double> >     cylrsumsqv_;
  RCP<vector<double> >     cylbsumsqv_;
  //! sum over w^2
  RCP<vector<double> >     x1csumsqw_;
  RCP<vector<double> >     x2csumsqw_;
  RCP<vector<double> >     x2w1sumsqw_;
  RCP<vector<double> >     x2w2sumsqw_;
  RCP<vector<double> >     cyllsumsqw_;
  RCP<vector<double> >     cyltsumsqw_;
  RCP<vector<double> >     cylrsumsqw_;
  RCP<vector<double> >     cylbsumsqw_;
  //! sum over uv
  RCP<vector<double> >     x1csumuv_;
  RCP<vector<double> >     x2csumuv_;
  RCP<vector<double> >     x2w1sumuv_;
  RCP<vector<double> >     x2w2sumuv_;
  RCP<vector<double> >     cyllsumuv_;
  RCP<vector<double> >     cyltsumuv_;
  RCP<vector<double> >     cylrsumuv_;
  RCP<vector<double> >     cylbsumuv_;
  //! sum over uw
  RCP<vector<double> >     x1csumuw_;
  RCP<vector<double> >     x2csumuw_;
  RCP<vector<double> >     x2w1sumuw_;
  RCP<vector<double> >     x2w2sumuw_;
  RCP<vector<double> >     cyllsumuw_;
  RCP<vector<double> >     cyltsumuw_;
  RCP<vector<double> >     cylrsumuw_;
  RCP<vector<double> >     cylbsumuw_;
  //! sum over vw
  RCP<vector<double> >     x1csumvw_;
  RCP<vector<double> >     x2csumvw_;
  RCP<vector<double> >     x2w1sumvw_;
  RCP<vector<double> >     x2w2sumvw_;
  RCP<vector<double> >     cyllsumvw_;
  RCP<vector<double> >     cyltsumvw_;
  RCP<vector<double> >     cylrsumvw_;
  RCP<vector<double> >     cylbsumvw_;
  //! sum over p^2
  RCP<vector<double> >     x1csumsqp_;
  RCP<vector<double> >     x2csumsqp_;
  RCP<vector<double> >     x2w1sumsqp_;
  RCP<vector<double> >     x2w2sumsqp_;
  RCP<vector<double> >     cyllsumsqp_;
  RCP<vector<double> >     cyltsumsqp_;
  RCP<vector<double> >     cylrsumsqp_;
  RCP<vector<double> >     cylbsumsqp_;

};

}

#endif  // not TURBULENCE_STATISTICS_SQC_H

