/*----------------------------------------------------------------------*/
/*! \file
\brief Interface class for random fields using FFT algorithms

\maintainer Jonas Nitzler

\level 2
*/
/*----------------------------------------------------------------------*/

#ifdef HAVE_FFTW
#ifndef RANDOMFIELD_H_
#define RANDOMFIELD_H_

// headers
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_RCP.hpp>
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "../drt_lib/drt_discret.H"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_LAPACK.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_RCP.hpp>
#include <Epetra_SerialDenseVector.h>
#include <math.h>


namespace UQ
{
  class RandomField
  {
   public:
    RandomField(Teuchos::RCP<DRT::Discretization> discret, const Teuchos::ParameterList& rfp);
    //! \brief Standard constructor

    virtual ~RandomField(){};

    //! \brief generates new realization of random field
    virtual void CreateNewSample(unsigned int seed) = 0;

    //! \brief get value of random field at a specific location
    virtual double EvalFieldAtLocation(const std::vector<double> location,
        const double paracont_parameter, const bool writetofile, const bool output) = 0;

    //! \brief get value of random field at a specific location
    virtual double EvalFieldAtLocation(const std::vector<double> location,
        const double local_median, const double paracont_parameter, const bool writetofile,
        const bool output) = 0;

    //! \brief write random field To file
    // virtual  void WriteRandomFieldToFile()=0;

    //! \brief write random field to file
    virtual void WriteRandomVariablesToFile(std::string filename) = 0;

    //! \brief return dimension of random field (either 2 or 3)
    virtual int Dimension() = 0;

    //! \brief compute bounding box of discretization
    virtual void ComputeBoundingBox(Teuchos::RCP<DRT::Discretization> discret);

    bool IsGaussian() const
    {
      if (marginal_pdf_ == normal)
        return true;
      else
        return false;
    }

    double LowerBound() const { return rf_lower_bound_; }

    double UpperBound() const { return rf_upper_bound_; }

    bool IsBounded() const { return is_bounded_; }

    bool HasSpatialMedian() const { return has_spatially_variable_median_; }

    //! brief get min corner of bounding box
    std::vector<double> GetMinBB() { return bb_min_; }

    //! brief get min corner of bounding box
    std::vector<double> GetMaxBB() { return bb_max_; }

   protected:
    //! \brief which proc am I on
    int myrank_;

    //! \brief parameters for evaluation i.e corners of the bounding box of the discretization
    //! \brief min bounding box
    std::vector<double> bb_min_;

    //! \brief max bounding box
    std::vector<double> bb_max_;

    double largestlength_;

    //! \brief do we truncate the field above and below a threshold
    bool is_bounded_;

    //! \brief lower cutoff threshold
    double rf_lower_bound_;

    //! \brief upper cutoff threshold
    double rf_upper_bound_;

    enum PROBABILITY_DISTRIBUTIONS
    {
      normal,
      lognormal,
      beta,
      none,
    };
    //! \brief target marginal pdf of random field
    PROBABILITY_DISTRIBUTIONS marginal_pdf_;

    //! brief field has spatially varying median
    bool has_spatially_variable_median_;

   private:
  };  // class RandomField
}  // namespace UQ
#endif  /* RANDOMFIELD_H_ */
#endif  // HAVE_FFTw
