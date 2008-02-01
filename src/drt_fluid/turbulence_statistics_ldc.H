/*!----------------------------------------------------------------------
\file turbulence_statistics_ldc.H

\brief Write (time and space) averaged values to file.


<pre>
o Create sets for centerlines in x1- and x2-direction
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

Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef TURBULENCE_STATISTICS_LDC_H
#define TURBULENCE_STATISTICS_LDC_H

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

class TurbulenceStatisticsLdc
{
  
public:
  /*!
  \brief Standard Constructor (public)

  <pre>
  o Create sets for centerlines in x1- and x2-direction

  o Allocate distributed vector for squares
  </pre>

  */
  TurbulenceStatisticsLdc(RefCountPtr<DRT::Discretization> actdis,
                          ParameterList&                   params);

  /*!
  \brief Destructor

  */
  virtual ~TurbulenceStatisticsLdc();


  //! @name functions for averaging

  /*!
  \brief The values of velocity and its squared velues are added to
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

  //! the coordinates of the centerlines in x1- and x2-direction
  RefCountPtr<vector<double> >     x1coordinates_;
  RefCountPtr<vector<double> >     x2coordinates_;

  //! sum over u (over the centerlines in x1- and x2-direction)
  RefCountPtr<vector<double> >     x1sumu_;
  RefCountPtr<vector<double> >     x2sumu_;
  //! sum over v (over the centerlines in x1- and x2-direction)
  RefCountPtr<vector<double> >     x1sumv_;
  RefCountPtr<vector<double> >     x2sumv_;
  //! sum over w (over the centerlines in x1- and x2-direction)
  RefCountPtr<vector<double> >     x1sumw_; 
  RefCountPtr<vector<double> >     x2sumw_;
  //! sum over p (over the centerlines in x1- and x2-direction)
  RefCountPtr<vector<double> >     x1sump_;
  RefCountPtr<vector<double> >     x2sump_;

  //! sum over u^2 (over the centerlines in x1- and x2-direction)
  RefCountPtr<vector<double> >     x1sumsqu_;
  RefCountPtr<vector<double> >     x2sumsqu_;
  //! sum over v^2 (over the centerlines in x1- and x2-direction)
  RefCountPtr<vector<double> >     x1sumsqv_;
  RefCountPtr<vector<double> >     x2sumsqv_;
  //! sum over w^2 (over the centerlines in x1- and x2-direction)
  RefCountPtr<vector<double> >     x1sumsqw_;
  RefCountPtr<vector<double> >     x2sumsqw_;
  //! sum over uv (over the centerlines in x1- and x2-direction)
  RefCountPtr<vector<double> >     x1sumuv_;
  RefCountPtr<vector<double> >     x2sumuv_;
  //! sum over uw (over the centerlines in x1- and x2-direction)
  RefCountPtr<vector<double> >     x1sumuw_;
  RefCountPtr<vector<double> >     x2sumuw_;
  //! sum over vw (over the centerlines in x1- and x2-direction)
  RefCountPtr<vector<double> >     x1sumvw_;
  RefCountPtr<vector<double> >     x2sumvw_;
  //! sum over p^2 (over the centerlines in x1- and x2-direction)
  RefCountPtr<vector<double> >     x1sumsqp_;
  RefCountPtr<vector<double> >     x2sumsqp_;

};

#endif  // not TURBULENCE_STATISTICS_LDC_H

#endif /* CCADISCRET       */
