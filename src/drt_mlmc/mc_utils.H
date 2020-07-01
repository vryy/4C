/*----------------------------------------------------------------------*/
/*! \file
\brief Utility functions for UQ analyis


\level 3
*/
/*----------------------------------------------------------------------*/

#ifdef HAVE_FFTW
#ifndef MC_UTILS_H_
#define MC_UTILS_H_

#include "Teuchos_RCP.hpp"


namespace UQ
{
  //! brief compares the first entry of a pair<double,pair<int,double> >. Used in sorting algorithm.
  //! Min for minimum search
  inline static bool MyComparePairsMin(const std::pair<double, std::pair<int, double>>& first,
      const std::pair<double, std::pair<int, double>>& second)
  {
    if (first.first < second.first)
      return true;
    else
      return false;
  }

  //! brief compares the first entry of a pair<double,pair<int,double> >. Used in sorting algorithm.
  //! Max for maximum search
  inline static bool MyComparePairsMax(const std::pair<double, std::pair<int, double>>& first,
      const std::pair<double, std::pair<int, double>>& second)
  {
    if (first.first > second.first)
      return true;
    else
      return false;
  }

  //! brief compute peak and quantile over all processors (for OutputUQ of red_airway)
  void ComputePeakAndQuantile(Teuchos::RCP<std::vector<double>> values1,
      Teuchos::RCP<std::vector<double>> values2, Teuchos::RCP<std::vector<int>> ele_ids,
      Teuchos::RCP<std::pair<double, std::pair<int, double>>> peak,
      Teuchos::RCP<std::pair<double, std::pair<int, double>>> quantile99, int tot_num_elements,
      std::string mode, const Epetra_Comm& mycomm);
}  // namespace UQ


#endif  /* MC_UTILS_H_ */
#endif  // HAVE_FFTw
