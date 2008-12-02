/*!
\file turbulence_statistics_mean_general.H

\brief Computation of mean values of nodal/cp quantities. The
means are computed as time averages. 

An additional method to do an averaging in a homogeneous, axis 
aligned, direction is provided. A requirement for this averaging
procedure is that all nodes are repeated in identical planes
along the axis along we do the averaging (like it is the case
for a flow around a square cylinder, a backward facing step,
a diffuser or a plane channel).  
--- all space averaged information is distributed back to all 
nodes (for easy visualisation in paraview). Note that this 
involves heavy communication --- including two round robin 
loops --- and produces additional output.

<pre>

Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235

</pre>
*/
#ifdef CCADISCRET
#ifndef TURBULENCE_STATISTICS_MEAN_GENERAL_H
#define TURBULENCE_STATISTICS_MEAN_GENERAL_H

#include "Epetra_Vector.h"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_io/io.H"

using namespace std;
using namespace Teuchos;


namespace FLD
{

class TurbulenceStatisticsGeneralMean
{
public:

  /*!
  \brief Constructor (public)

  \param (in) the discretisation (containing nodes, dofs etc.)

  */
  TurbulenceStatisticsGeneralMean(RCP<DRT::Discretization> discret,string homdir);

  /*!
  \brief Destructor (public)

  */
  virtual ~TurbulenceStatisticsGeneralMean();

  //! @name Averaging

  /*!
  \brief Add vector to current time average

  \param dt  (in) time contribution corresponding to this sample
  \param vec (in) vector to add to time average

  */
  void AddToCurrentTimeAverage(const double             dt , 
                               const RCP<Epetra_Vector> vec);

  /*!
  \brief Perform a averaging of the current, already time averaged
  vector, in space in a homogeneous direction.

  \param dim (in) dimension to average in

  */
  void SpaceAverageInOneDirection(const int dim);

  /*!
  \brief Add vector to time average from previous steps

  */
  void AddToTotalTimeAverage();

  //! @name IO

  /*!
  \brief Read previous statistics from a file (for restart)

  \param (in) input reader to allow restart

  */
  void ReadOldStatistics(IO::DiscretizationReader&  input);


  /*!
  \brief Write the statistics to a file

  \param (in) output context

  */
  void WriteOldAverageVec(IO::DiscretizationWriter&  output);

  //! @name Misc

  /*!
  \brief Clear all statistics collected in the current period

  */
  void TimeReset();


  /*!
  \brief Clear all statistics collected up to now

  */
  void ResetComplete();

  //@}

private:

  //! the discretization
  RCP<DRT::Discretization>     discret_      ;

  //! vector containing homogeneous directions
  vector<int>                  homdir_       ;

  //! previous averages, done in time and space
  RCP<Epetra_Vector>           prev_avg_     ;
  //! number of time steps included in the previous average
  int                          prev_n_       ;
  //! time covered by previous average
  double                       prev_avg_time_;


  //! current averages, done in time and space
  RCP<Epetra_Vector>           curr_avg_     ;
  //! number of time steps included in the current average
  int                          curr_n_       ;
  //! time covered by current average
  double                       curr_avg_time_;

  //! compare operator for doubles up to a precision of 1e-8
  struct doublecomp {
    bool operator() (const double& lhs, const double& rhs) const
      {return lhs<rhs-1e-8;}
  };
}; // end class TurbulenceStatisticsGeneralMean

} // end namespace FLD

#endif  // #ifndef TURBULENCE_STATISTICS_MEAN_GENERAL_H
#endif  // #ifdef CCADISCRET
