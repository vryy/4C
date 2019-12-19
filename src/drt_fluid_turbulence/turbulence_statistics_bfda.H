/*----------------------------------------------------------------------*/
/*! \file

\brief routines for bfda test case (Semesterarbeit Sebastian Buechner)
  ___________          ______________
 |           \ ______|              |
 |            _______O              |     ----------> z direction
 |___________/       |______________|

 (Origin of coordinate system is marked by O)


\maintainer Martin Kronbichler

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef TURBULENCE_STATISTICS_BFDA_H
#define TURBULENCE_STATISTICS_BFDA_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../drt_inpar/inpar_fluid.H"


#include "Teuchos_SerialDenseMatrix.hpp"

namespace FLD
{
  class TurbulenceStatisticsBfda
  {
   public:
    /*!
    \brief Standard Constructor (public)

        o Create vector "zcoordinates" with coordinates along z-axis
    o Create matrix "rcoordinates_" with coordinates of radial coordinates of evaluation planes
      columns are evaluation planes corresponding to the positions in "posEvaluation_"

    */
    TurbulenceStatisticsBfda(Teuchos::RCP<DRT::Discretization> actdis,
        Teuchos::ParameterList& params, const std::string& statistics_outfilename);

    /*!
    \brief Destructor

    */
    virtual ~TurbulenceStatisticsBfda();


    //! @name functions for averaging


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
      bool operator()(const double& p1, const double& p2) const { return (p1 < p2 - 1E-9); }

     protected:
     private:
    };

   private:
    //! number of samples taken
    int numsamp_;


    //! number of coordinates in z-direction
    int numzcoor_;

    //! number in radial direction for statistical evaluation
    int numrstatlocations_;

    //! The discretisation (required for nodes, dofs etc;)
    Teuchos::RCP<DRT::Discretization> discret_;

    //! parameter list
    Teuchos::ParameterList& params_;

    //! name of statistics output file, despite the ending
    const std::string statistics_outfilename_;

    //! positions for evaluation
    std::vector<double> posEvaluation_;

    //! positions for evaluation of actual mesh
    std::vector<double> actPosEvaluation_;

    //! Toggle vectors
    Teuchos::RCP<Epetra_Vector> togglew_;
    Teuchos::RCP<Epetra_Vector> togglep_;

    //! vector for z-coordinates
    Teuchos::RCP<std::vector<double>> zcoordinates_;
    //! matrix for r-coordinates (columns are evaluation planes
    Teuchos::SerialDenseMatrix<int, double> rcoordinates_;

    //! matrices for mean values along z-axis
    Teuchos::RCP<Epetra_SerialDenseMatrix> zsumw_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> zsump_;

    //! matrices for mean values of evaluation planes
    Teuchos::RCP<Epetra_SerialDenseMatrix> rsumw_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> rsump_;
  };

}  // namespace FLD

#endif  // not TURBULENCE_STATISTICS_BFS_H
