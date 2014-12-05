/*!----------------------------------------------------------------------
\file turbulence_statistics_bfs.H

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

#ifndef TURBULENCE_STATISTICS_PH_H
#define TURBULENCE_STATISTICS_PH_H

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


namespace FLD
{

class TurbulenceStatisticsPh
{

public:
  /*!
  \brief Standard Constructor (public)

  <pre>
  o Create sets for lines

  o Allocate distributed vector for squares
  </pre>

  */
  TurbulenceStatisticsPh(Teuchos::RCP<DRT::Discretization> actdis,
                          Teuchos::ParameterList&           params);

  /*!
  \brief Destructor

  */
  virtual ~TurbulenceStatisticsPh();


  //! @name functions for averaging

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

  //! geometry of DNS of incompressible flow over Ph by Le, Moin and Kim or geometry of Avancha and Pletcher of LES of flow over Ph with heating
  enum GeoType{
   none,
   geometry_DNS_incomp_flow,
   geometry_LES_flow_with_heating
  };

  //! number of samples taken
  int  numsamp_;

  //! number of coordinates in x1- and x2-direction
  int  numx1coor_;
  int  numx2coor_;

  //! number of locations in x1- and x2-direction for statistical evaluation
  int  numx1statlocations_;
  int  numx2statlocations_;
//  int  numx1supplocations_;


  //! bounds for extension of backward-facing step in x1-direction
  double  x1min_;
  double  x1max_;

  //! bounds for extension of backward-facing step in x2-direction
  double  x2min_;
  double  x2max_;

  //! bounds for extension of backward-facing step in x3-direction
  double  x3min_;
  double  x3max_;

  //! length of domain in x1-direction
//  double x1dist_;

  //! half of the element length in x1-direction
  double x1elemlengthhalf;
//  double * x1elemlengthhalf_;
//  x1elemlengthhalf_ = &x1elemlengthhalf;

  //! number of elements in x1-direction
//  int numx1elem_;

  //! distance between x1-coordinates and x1setstatlocations
  double dist;

  //! allowed minimum distance between x1-coordinates and x1setstatlocations = x1elemlengthhalf_
  double mindist;


  //! The discretisation (required for nodes, dofs etc;)
  Teuchos::RCP<DRT::Discretization> discret_;

  //! parameter list
  Teuchos::ParameterList&                   params_;


//  //! boolean indicating turbulent inflow channel discretization
//  const bool inflowchannel_;
//  //! x-coordinate of outflow of inflow channel
//  const double inflowmax_;

   //! pointer to vel/pres^2 field (space allocated in constructor)
  Teuchos::RCP<Epetra_Vector>       squaredvelnp_;
  //! pointer to T^2 field (space allocated in constructor)
  Teuchos::RCP<Epetra_Vector>       squaredscanp_;
  //! pointer to 1/T field (space allocated in constructor)
  Teuchos::RCP<Epetra_Vector>       invscanp_;
  //! pointer to (1/T)^2 field (space allocated in constructor)
  Teuchos::RCP<Epetra_Vector>       squaredinvscanp_;

  //! toogle vectors: sums are computed by scalarproducts
  Teuchos::RCP<Epetra_Vector>       toggleu_;
  Teuchos::RCP<Epetra_Vector>       togglev_;
  Teuchos::RCP<Epetra_Vector>       togglew_;
  Teuchos::RCP<Epetra_Vector>       togglep_;

  //! available x1- and x2-coordinates
  Teuchos::RCP<std::vector<double> >     x1coordinates_;
  Teuchos::RCP<std::vector<double> >     x2coordinates_;

  //! coordinates of locations in x1- and x2-direction for statistical evaluation
  LINALG::Matrix<4,1> x1statlocations_;
//  LINALG::Matrix<17,1> x2statlocations_;

  LINALG::Matrix<4,17> x2statlocations;


  //! set coordinates of locations in x1-direction for statistical evaluation
  Teuchos::RCP<std::vector<double> > x1setstatlocations_;

  //! coordinates of supplementary locations in x2-direction for velocity derivative
//  LINALG::Matrix<2,1>  x2supplocations_;

  //! coordinates of supplementary locations in x1-direction for statistical evaluation
  //! (check of inflow profile)
  LINALG::Matrix<10,1>  x1supplocations_;

  //! matrices containing values
  Teuchos::RCP<Epetra_SerialDenseMatrix>     x1sumu_;
  Teuchos::RCP<Epetra_SerialDenseMatrix>     x1sump_;

  Teuchos::RCP<Epetra_SerialDenseMatrix>     x1sumrho_;
  Teuchos::RCP<Epetra_SerialDenseMatrix>     x1sumT_;

  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sumu_;
  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sumv_;
  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sumw_;
  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sump_;

  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sumrho_;
  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sumT_;

  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sumsqu_;
  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sumsqv_;
  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sumsqw_;
  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sumsqp_;

  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sumsqrho_;
  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sumsqT_;

  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sumuv_;
  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sumuw_;
  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sumvw_;

  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sumrhou_;
  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sumuT_;
  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sumrhov_;
  Teuchos::RCP<Epetra_SerialDenseMatrix>     x2sumvT_;

};

}

#endif  // not TURBULENCE_STATISTICS_PH_H

