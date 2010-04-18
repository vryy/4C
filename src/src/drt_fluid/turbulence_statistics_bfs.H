/*!----------------------------------------------------------------------
\file turbulence_statistics_sqc.H

\brief Write (time and space) averaged values to file for
turbulent flow over a backward-facing step


<pre>
o Create sets for various evaluation lines in domain
  (Construction based on a round robin communication pattern):
  - 21 lines in x2-direction
  - lines along upper and lower wall

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

#ifndef TURBULENCE_STATISTICS_BFS_H
#define TURBULENCE_STATISTICS_BFS_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_inpar/inpar_fluid.H"


using namespace std;
using namespace Teuchos;

namespace FLD
{

class TurbulenceStatisticsBfs
{

public:
  /*!
  \brief Standard Constructor (public)

  <pre>
  o Create sets for lines

  o Allocate distributed vector for squares
  </pre>

  */
  TurbulenceStatisticsBfs(RCP<DRT::Discretization> actdis,
                          ParameterList&                   params);

  /*!
  \brief Destructor

  */
  virtual ~TurbulenceStatisticsBfs();


  //! @name functions for averaging

  /*!
  \brief The values of velocity and its squared values are added to
  global vectors. This method allows to do the time average of the
  nodal values after a certain amount of timesteps.
  */
  void DoTimeSample(Teuchos::RCP<Epetra_Vector> velnp);

  /*!
  \brief The values of velocity, pressure, temperature and its squared
  values are added to global vectors. This method allows to do the time
  average of the nodal values after a certain amount of timesteps.
  */
  void DoLomaTimeSample(Teuchos::RCP<Epetra_Vector> velnp,
                        Teuchos::RCP<Epetra_Vector> scanp,
                        const double                eosfac);

  /*!
  \brief Dump the result to file.

  step on input is used to print the timesteps which belong to the
  statistic to the file
  */

  void DumpStatistics(int step);

  /*!
  \brief Dump the result to file for low-Mach-number flow.

  step on input is used to print the timesteps which belong to the
  statistic to the file
  */

  void DumpLomaStatistics(int step, const double eosfac);


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
  int  numsamp_;

  //! number of coordinates in x1- and x2-direction
  int  numx1coor_;
  int  numx2coor_;

  //! number of locations in x1- and x2-direction for statistical evaluation
  int  numx1statlocations_;
  int  numx2statlocations_;

  //! bounds for extension of backward-facing step in x3-direction
  double  x3min_;
  double  x3max_;

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

  //! available x1- and x2-coordinates
  RCP<vector<double> >     x1coordinates_;
  RCP<vector<double> >     x2coordinates_;

  //! coordinates of locations in x1- and x2-direction for statistical evaluation
  LINALG::Matrix<21,1> x1statlocations_;
  LINALG::Matrix<2,1>  x2statlocations_;

  //! coordinates of supplementary locations in x2-direction for velocity derivative
  LINALG::Matrix<2,1>  x2supplocations_;

  //! matrices containing values
  RCP<Epetra_SerialDenseMatrix>     x1sumu_;
  RCP<Epetra_SerialDenseMatrix>     x1sump_;

  RCP<Epetra_SerialDenseMatrix>     x1sumrho_;
  RCP<Epetra_SerialDenseMatrix>     x1sumT_;

  RCP<Epetra_SerialDenseMatrix>     x2sumu_;
  RCP<Epetra_SerialDenseMatrix>     x2sumv_;
  RCP<Epetra_SerialDenseMatrix>     x2sumw_;
  RCP<Epetra_SerialDenseMatrix>     x2sump_;

  RCP<Epetra_SerialDenseMatrix>     x2sumrho_;
  RCP<Epetra_SerialDenseMatrix>     x2sumT_;

  RCP<Epetra_SerialDenseMatrix>     x2sumsqu_;
  RCP<Epetra_SerialDenseMatrix>     x2sumsqv_;
  RCP<Epetra_SerialDenseMatrix>     x2sumsqw_;
  RCP<Epetra_SerialDenseMatrix>     x2sumsqp_;

  RCP<Epetra_SerialDenseMatrix>     x2sumsqrho_;
  RCP<Epetra_SerialDenseMatrix>     x2sumsqT_;

  RCP<Epetra_SerialDenseMatrix>     x2sumuv_;
  RCP<Epetra_SerialDenseMatrix>     x2sumuw_;
  RCP<Epetra_SerialDenseMatrix>     x2sumvw_;

  RCP<Epetra_SerialDenseMatrix>     x2sumrhou_;
  RCP<Epetra_SerialDenseMatrix>     x2sumrhouT_;
  RCP<Epetra_SerialDenseMatrix>     x2sumrhov_;
  RCP<Epetra_SerialDenseMatrix>     x2sumrhovT_;
};

}

#endif  // not TURBULENCE_STATISTICS_SQC_H

#endif /* CCADISCRET       */
