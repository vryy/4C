/*!----------------------------------------------------------------------
\file turbulence_statistics.H

\brief Write (time and space) averaged values to file.


<pre>
o Create set of all available homogeneous planes
  (Construction based on a round robin communication pattern)

o loop planes (e.g. plane coordinates)

  - generate 4 toggle vectors (u,v,w,p), for example

                            /  1  u dof in homogeneous plane
                 toggleu_  |
                            \  0  elsewhere

  - pointwise multiplication velnp.*velnp for second order
    moments

  - 2 * 4 scalarproducts for in plane mean values

o in plane mean values are averaged in time over all steps between two
  outputs

o Write statistics for the Smagorinsky "constant" Cs if a dynamic
  procedure to determine it is applied
  
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
#include "../drt_lib/linalg_utils.H"

#include "../drt_mat/newtonianfluid.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


using namespace std;
using namespace Teuchos;

class TurbulenceStatistics
{
  
public:
  /*!
  \brief Standard Constructor (public)

  <pre>
  o Create vector of homogeneous plane coordinates

  o Allocate 4 distributed toggle vectors and one distributed vector
    for squares
  </pre>

  */
  TurbulenceStatistics(RefCountPtr<DRT::Discretization> actdis,
                       ParameterList&                   params);

  /*!
  \brief Destructor

  */
  virtual ~TurbulenceStatistics();


  //! @name functions for averaging

  /*!
  \brief The values of velocity and its squared velues are added to
  global vectors. This method allows to do the time average of the
  nodal values after a certain amount of timesteps.
  
  */
  void DoTimeSample(
    Teuchos::RefCountPtr<Epetra_Vector> velnp,
    Epetra_Vector & force
    );

  
  /*!
  \brief Compute in plane means of u,u^2 etc. (nodal quantities)

  The calculated value is added to the sum**,sumsq** variables in the
  component corresponding to the plane.

  velnp is the solution vector provided by the time integration
  algorithm 
  */
  void EvaluateMeanValuesInPlanes();

  /*!
  \brief 

  */
  void TimeAverageMeans();

  
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
  \brief Provide the coordinates of the homogeneous planes for a
  turbulent channel flow
  
  */
  vector<double> ReturnNodePlaneCoords()
    {
      return(*nodeplanes_);
    };

#if 0
  /*!
  \brief Return pointer to vector of values of the Smagorinsky constant.
  It will be passed to the element in order to compute average values.
  
  */
  RefCountPtr<vector<double> > ReturnCsAverage()
    {
      return(incrsumCs_);
    };
#endif
  /*!
  \brief Replace increment to compute average of Smagorinsky Constant, effective
  viscosity and (Cs_delta)^2
  
  */
   void ReplaceCsIncrement(
     RefCountPtr<vector<double> > incrsumCs,
     RefCountPtr<vector<double> > incr_Cs_delta_sq_sum,
     RefCountPtr<vector<double> > incr_visceff_sum
    )
    {
      for (unsigned rr=0;rr<incrsumCs->size();++rr)
      {
        (*incrsumCs_         )[rr] =(*incrsumCs           )[rr];
        (*incrsumCs_delta_sq_)[rr] =(*incr_Cs_delta_sq_sum)[rr];
        (*incrsumvisceff_    )[rr] =(*incr_visceff_sum    )[rr];
      }
      return;
    };
  
  /*!
  \brief Add computed sums to compute average of residuals
  
  */
   void AddToResAverage(
     RefCountPtr<vector<double> > incrres,
     RefCountPtr<vector<double> > incrres_sq,
     RefCountPtr<vector<double> > incrsacc,
     RefCountPtr<vector<double> > incrsacc_sq
    )
    {
      for (unsigned rr=0;rr<incrres->size();++rr)
      {
          (*sumres_    )[rr]+=(*incrres)    [rr];
          (*sumres_sq_ )[rr]+=(*incrres_sq) [rr];
          (*sumsacc_   )[rr]+=(*incrsacc)   [rr];
          (*sumsacc_sq_)[rr]+=(*incrsacc_sq)[rr];
      }
      return;
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

  //! The discretisation (required for nodes, dofs etc;)
  RefCountPtr<DRT::Discretization> discret_;

  //! contains plane normal direction etc
  ParameterList&                   params_;

   //! pointer to vel/pres field
  RefCountPtr<Epetra_Vector>       meanvelnp_;

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
  
  //! sum over u (over one plane in each component)
  RefCountPtr<vector<double> >     sumu_;
  //! sum over v (over one plane in each component)
  RefCountPtr<vector<double> >     sumv_;
  //! sum over w (over one plane in each component)
  RefCountPtr<vector<double> >     sumw_; 
  //! sum over p (over one plane in each component)
  RefCountPtr<vector<double> >     sump_;

  //! sum over u^2 (over one plane in each component)
  RefCountPtr<vector<double> >     sumsqu_;
  //! sum over v^2 (over one plane in each component)
  RefCountPtr<vector<double> >     sumsqv_;
  //! sum over w^2 (over one plane in each component)
  RefCountPtr<vector<double> >     sumsqw_;
  //! sum over uv (over one plane in each component)
  RefCountPtr<vector<double> >     sumuv_;
  //! sum over uw (over one plane in each component)
  RefCountPtr<vector<double> >     sumuw_;
  //! sum over vw (over one plane in each component)
  RefCountPtr<vector<double> >     sumvw_;
  //! sum over p^2 (over one plane in each component)
  RefCountPtr<vector<double> >     sumsqp_;


#if 1
  //! vector of squared velocities (pointwise)
  RefCountPtr<Epetra_Vector>       pointsquaredvelnp_;


//! sum over u (over one plane in each component)
  RefCountPtr<vector<double> >     pointsumu_;
  //! sum over v (over one plane in each component)
  RefCountPtr<vector<double> >     pointsumv_;
  //! sum over w (over one plane in each component)
  RefCountPtr<vector<double> >     pointsumw_; 
  //! sum over p (over one plane in each component)
  RefCountPtr<vector<double> >     pointsump_;

  //! sum over u^2 (over one plane in each component)
  RefCountPtr<vector<double> >     pointsumsqu_;
  //! sum over v^2 (over one plane in each component)
  RefCountPtr<vector<double> >     pointsumsqv_;
  //! sum over w^2 (over one plane in each component)
  RefCountPtr<vector<double> >     pointsumsqw_;
  //! sum over p^2 (over one plane in each component)
  RefCountPtr<vector<double> >     pointsumsqp_;
#endif
  //! sum over nodal forces on boundary in u direction
  double                           sumforceu_;
  //! sum over nodal forces on boundary in v direction
  double                           sumforcev_;
  //! sum over nodal forces on boundary in w direction
  double                           sumforcew_;

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

  //! sum over all in plane residuals
  RefCountPtr<vector<double> >     sumres_;
  //! sum over all in plane squared residuals
  RefCountPtr<vector<double> >     sumres_sq_;
  //! sum over all in plane subscale accelerations
  RefCountPtr<vector<double> >     sumsacc_;
  //! sum over all in plane squared subscale accelerations
  RefCountPtr<vector<double> >     sumsacc_sq_;
  
  
  //! viscosity to calculate l_tau, y+ etc.
  double                           visc_;

  
};

#endif  // not TURBULENCE_STATISTICS_H

#endif /* CCADISCRET       */
