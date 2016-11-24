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

\level 2

\maintainer Benjamin Krank

*/
#ifndef TURBULENCE_STATISTICS_MEAN_GENERAL_H
#define TURBULENCE_STATISTICS_MEAN_GENERAL_H

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

namespace LINALG
{
  class MapExtractor;
}


// forward declarations
namespace DRT
{
  class Discretization;
  class DofSet;
}
namespace IO
{
  class DiscretizationReader;
  class DiscretizationWriter;
}


namespace FLD
{

class TurbulenceStatisticsGeneralMean
{
public:

  /*!
  \brief Constructor (public)

  \param (in) the discretisation (containing nodes, dofs etc.)

  */
  TurbulenceStatisticsGeneralMean(
    Teuchos::RCP<DRT::Discretization> discret        ,
    std::string                   homdir         ,
    LINALG::MapExtractor&    velpressplitter,
    const bool               withscatra
    );

  TurbulenceStatisticsGeneralMean(
    Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<const DRT::DofSet>   standarddofset ,
    std::string                            homdir,
    LINALG::MapExtractor&             velpressplitter,
    const bool                        withscatra
    );

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
                               const Teuchos::RCP<Epetra_Vector> vec,
                               const Teuchos::RCP<Epetra_Vector> scavec = Teuchos::null,
                               const Teuchos::RCP<Epetra_Vector> scatravec = Teuchos::null);

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
  \brief Read previous scatra statistics from a file (for restart)

  \param (in) input reader to allow restart

  */
  void ReadOldStatisticsScaTra(IO::DiscretizationReader&  input);


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
  \brief Clear all statistics vectors based on fluid maps collected in the current period

  */
  void TimeResetFluidAvgVectors(const Epetra_Map& dofrowmap);

  /*!
  \brief Clear all statistics collected up to now

  */
  void ResetComplete();

  /*!
  \brief Clear statistics vectors based on fluid maps

  */
  void ResetFluidAvgVectors(const Epetra_Map& dofrowmap);

  /*!
  \brief Redistribute all statistics vectors

  */
  void Redistribute(
      Teuchos::RCP<const DRT::DofSet>   standarddofset,
      Teuchos::RCP<DRT::Discretization> discret);

  /*!
  \brief Add results from scalar transport field solver to statistics

  */
  void AddScaTraResults(
      Teuchos::RCP<DRT::Discretization> scatradis,
      Teuchos::RCP<Epetra_Vector> myphinp);

  /*!
  \brief Do output of ScaTra mean field for visualization/restart
         (statistics was already written during call of DoOutput())

  */
  void
  DoOutputForScaTra(
      IO::DiscretizationWriter& output,
      int                       step);

  //@}

private:

  //! the fluid discretization
  Teuchos::RCP<DRT::Discretization> discret_;

  //! dofset containing fluid standard dofs (no XFEM dofs)
  Teuchos::RCP<const DRT::DofSet> standarddofset_;

  //! the scatra discretization
  Teuchos::RCP<DRT::Discretization> scatradis_;

  //! a splitter between velocities and pressure dofs
  LINALG::MapExtractor&        velpressplitter_;

  //! vector containing homogeneous directions
  std::vector<int>             homdir_         ;

  //! previous averages, done in time and space
  Teuchos::RCP<Epetra_Vector>  prev_avg_       ;
  //! number of time steps included in the previous average
  int                          prev_n_         ;
  //! time covered by previous average
  double                       prev_avg_time_  ;


  //! current averages, done in time and space
  Teuchos::RCP<Epetra_Vector>  curr_avg_       ;
  //! number of time steps included in the current average
  int                          curr_n_         ;
  //! time covered by current average
  double                       curr_avg_time_  ;

  //! flag for additional averaging of scalar field
  bool                         withscatra_     ;
  //! previous scalar field averages, done in time and space
  Teuchos::RCP<Epetra_Vector>  prev_avg_sca_   ;
  //! current scalar field averages, done in time and space
  Teuchos::RCP<Epetra_Vector>  curr_avg_sca_   ;
  //! previous scalar field averages, done in time and space
  Teuchos::RCP<Epetra_Vector>  prev_avg_scatra_;
  //! current scalar field averages, done in time and space
  Teuchos::RCP<Epetra_Vector>  curr_avg_scatra_;

  //! compare operator for doubles up to a precision of 1e-8
  struct doublecomp {
    bool operator() (const double lhs, const double rhs) const
      {return lhs<(rhs-1e-8);}
  };
}; // end class TurbulenceStatisticsGeneralMean

} // end namespace FLD

#endif  // #ifndef TURBULENCE_STATISTICS_MEAN_GENERAL_H
