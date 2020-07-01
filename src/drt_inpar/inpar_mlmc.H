/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for multi-level monte carlo


\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef INPAR_MLMC_H
#define INPAR_MLMC_H

#include "inpar_parameterlist_utils.H"
/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace MLMC
  {
    enum CorrStruct
    {
      corr_gaussian,  ///< gaussian correlation structure
      // more to follow
    };
    enum MarginalPdf
    {
      pdf_gaussian,   // normal distribution
      pdf_beta,       // beta distribution
      pdf_lognormal,  // lognormal distribution
      // more to follow
    };
    enum CalcMethod
    {
      calc_m_fft,      // using fftw
      calc_m_cos,      // as cosinus series
      calc_m_fourier,  // as fourier series
    };
    enum UQStrategy
    {
      mc_plain,        // do plain and simple MC
      mc_paracont,     //  use parameter continuation scheme in MC
      mc_scaledthick,  // compute mechanical quantities based on scaled thickness
    };

    enum FWDProblem
    {
      structure,    // structure problem
      red_airways,  // reduced airwarys problem
    };
    /// To be filled with input params


    /// set the mlmc parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace MLMC

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_MLMC_H
