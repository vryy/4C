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
#ifdef CCADISCRET

#ifndef TURBULENCE_STATISTICS_SQC_H
#define TURBULENCE_STATISTICS_SQC_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_utils.H"


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
  TurbulenceStatisticsSqc(RefCountPtr<DRT::Discretization> actdis,
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
  void DoTimeSample(Teuchos::RefCountPtr<Epetra_Vector> velnp);

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

  //! bounds for extension of cavity in x3-direction
  double                           x3min_;
  double                           x3max_;

  //! sums over lift and drag values
  double                           lift_;
  double                           drag_;
  double                           liftsq_;
  double                           dragsq_;

  //! The discretisation (required for nodes, dofs etc;)
  RefCountPtr<DRT::Discretization> discret_;

  //! parameter list
  ParameterList&                   params_;

   //! pointer to vel/pres^2 field (space allocated in constructor)
  RefCountPtr<Epetra_Vector>       squaredvelnp_;

  //! toogle vectors: sums are computed by scalarproducts
  RefCountPtr<Epetra_Vector>       toggleu_;
  RefCountPtr<Epetra_Vector>       togglev_;
  RefCountPtr<Epetra_Vector>       togglew_;
  RefCountPtr<Epetra_Vector>       togglep_;

  //! the coordinates of various lines
  RefCountPtr<vector<double> >     x1ccoordinates_;
  RefCountPtr<vector<double> >     x2ccoordinates_;
  RefCountPtr<vector<double> >     x2wcoordinates_;
  RefCountPtr<vector<double> >     clrcoordinates_;
  RefCountPtr<vector<double> >     ctbcoordinates_;

  //! sum over u
  RefCountPtr<vector<double> >     x1csumu_;
  RefCountPtr<vector<double> >     x2csumu_;
  RefCountPtr<vector<double> >     x2w1sumu_;
  RefCountPtr<vector<double> >     x2w2sumu_;
  RefCountPtr<vector<double> >     cyllsumu_;
  RefCountPtr<vector<double> >     cyltsumu_;
  RefCountPtr<vector<double> >     cylrsumu_;
  RefCountPtr<vector<double> >     cylbsumu_;
  //! sum over v
  RefCountPtr<vector<double> >     x1csumv_;
  RefCountPtr<vector<double> >     x2csumv_;
  RefCountPtr<vector<double> >     x2w1sumv_;
  RefCountPtr<vector<double> >     x2w2sumv_;
  RefCountPtr<vector<double> >     cyllsumv_;
  RefCountPtr<vector<double> >     cyltsumv_;
  RefCountPtr<vector<double> >     cylrsumv_;
  RefCountPtr<vector<double> >     cylbsumv_;
  //! sum over w
  RefCountPtr<vector<double> >     x1csumw_;
  RefCountPtr<vector<double> >     x2csumw_;
  RefCountPtr<vector<double> >     x2w1sumw_;
  RefCountPtr<vector<double> >     x2w2sumw_;
  RefCountPtr<vector<double> >     cyllsumw_;
  RefCountPtr<vector<double> >     cyltsumw_;
  RefCountPtr<vector<double> >     cylrsumw_;
  RefCountPtr<vector<double> >     cylbsumw_;
  //! sum over p
  RefCountPtr<vector<double> >     x1csump_;
  RefCountPtr<vector<double> >     x2csump_;
  RefCountPtr<vector<double> >     x2w1sump_;
  RefCountPtr<vector<double> >     x2w2sump_;
  RefCountPtr<vector<double> >     cyllsump_;
  RefCountPtr<vector<double> >     cyltsump_;
  RefCountPtr<vector<double> >     cylrsump_;
  RefCountPtr<vector<double> >     cylbsump_;

  //! sum over u^2
  RefCountPtr<vector<double> >     x1csumsqu_;
  RefCountPtr<vector<double> >     x2csumsqu_;
  RefCountPtr<vector<double> >     x2w1sumsqu_;
  RefCountPtr<vector<double> >     x2w2sumsqu_;
  RefCountPtr<vector<double> >     cyllsumsqu_;
  RefCountPtr<vector<double> >     cyltsumsqu_;
  RefCountPtr<vector<double> >     cylrsumsqu_;
  RefCountPtr<vector<double> >     cylbsumsqu_;
  //! sum over v^2
  RefCountPtr<vector<double> >     x1csumsqv_;
  RefCountPtr<vector<double> >     x2csumsqv_;
  RefCountPtr<vector<double> >     x2w1sumsqv_;
  RefCountPtr<vector<double> >     x2w2sumsqv_;
  RefCountPtr<vector<double> >     cyllsumsqv_;
  RefCountPtr<vector<double> >     cyltsumsqv_;
  RefCountPtr<vector<double> >     cylrsumsqv_;
  RefCountPtr<vector<double> >     cylbsumsqv_;
  //! sum over w^2
  RefCountPtr<vector<double> >     x1csumsqw_;
  RefCountPtr<vector<double> >     x2csumsqw_;
  RefCountPtr<vector<double> >     x2w1sumsqw_;
  RefCountPtr<vector<double> >     x2w2sumsqw_;
  RefCountPtr<vector<double> >     cyllsumsqw_;
  RefCountPtr<vector<double> >     cyltsumsqw_;
  RefCountPtr<vector<double> >     cylrsumsqw_;
  RefCountPtr<vector<double> >     cylbsumsqw_;
  //! sum over uv
  RefCountPtr<vector<double> >     x1csumuv_;
  RefCountPtr<vector<double> >     x2csumuv_;
  RefCountPtr<vector<double> >     x2w1sumuv_;
  RefCountPtr<vector<double> >     x2w2sumuv_;
  RefCountPtr<vector<double> >     cyllsumuv_;
  RefCountPtr<vector<double> >     cyltsumuv_;
  RefCountPtr<vector<double> >     cylrsumuv_;
  RefCountPtr<vector<double> >     cylbsumuv_;
  //! sum over uw
  RefCountPtr<vector<double> >     x1csumuw_;
  RefCountPtr<vector<double> >     x2csumuw_;
  RefCountPtr<vector<double> >     x2w1sumuw_;
  RefCountPtr<vector<double> >     x2w2sumuw_;
  RefCountPtr<vector<double> >     cyllsumuw_;
  RefCountPtr<vector<double> >     cyltsumuw_;
  RefCountPtr<vector<double> >     cylrsumuw_;
  RefCountPtr<vector<double> >     cylbsumuw_;
  //! sum over vw
  RefCountPtr<vector<double> >     x1csumvw_;
  RefCountPtr<vector<double> >     x2csumvw_;
  RefCountPtr<vector<double> >     x2w1sumvw_;
  RefCountPtr<vector<double> >     x2w2sumvw_;
  RefCountPtr<vector<double> >     cyllsumvw_;
  RefCountPtr<vector<double> >     cyltsumvw_;
  RefCountPtr<vector<double> >     cylrsumvw_;
  RefCountPtr<vector<double> >     cylbsumvw_;
  //! sum over p^2
  RefCountPtr<vector<double> >     x1csumsqp_;
  RefCountPtr<vector<double> >     x2csumsqp_;
  RefCountPtr<vector<double> >     x2w1sumsqp_;
  RefCountPtr<vector<double> >     x2w2sumsqp_;
  RefCountPtr<vector<double> >     cyllsumsqp_;
  RefCountPtr<vector<double> >     cyltsumsqp_;
  RefCountPtr<vector<double> >     cylrsumsqp_;
  RefCountPtr<vector<double> >     cylbsumsqp_;

};

}

#endif  // not TURBULENCE_STATISTICS_SQC_H

#endif /* CCADISCRET       */
