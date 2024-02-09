/*----------------------------------------------------------------------*/
/*! \file
 \brief main file containing routines for calculation of scatra element with chemotactic terms

\level 2

 *----------------------------------------------------------------------*/


#ifndef BACI_SCATRA_ELE_CALC_CHEMO_HPP
#define BACI_SCATRA_ELE_CALC_CHEMO_HPP

#include "baci_config.hpp"

#include "baci_mat_scatra_chemotaxis_mat.hpp"
#include "baci_scatra_ele_calc.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype>>
    class ScaTraEleCalcChemo : public virtual ScaTraEleCalc<distype, probdim>
    {
     protected:
      //! (private) protected constructor, since we are a Singleton.
      ScaTraEleCalcChemo(const int numdofpernode, const int numscal, const std::string& disname);

     private:
      typedef ScaTraEleCalc<distype, probdim> my;
      using my::nen_;
      using my::nsd_;
      using varmanager = ScaTraEleInternalVariableManager<nsd_, nen_>;

     public:
      //! Singleton access method
      static ScaTraEleCalcChemo<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

     protected:
      //! calculation of chemotactic element matrix
      void CalcMatChemo(CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac,
          const double timetaufac, const double densnp, const double scatrares,
          const CORE::LINALG::Matrix<nen_, 1>& sgconv,
          const CORE::LINALG::Matrix<nen_, 1>& diff) override;

      //! calculation of chemotactic RHS
      void CalcRHSChemo(CORE::LINALG::SerialDenseVector& erhs, const int k, const double rhsfac,
          const double rhstaufac, const double scatrares, const double densnp) override;

      //! get the material parameters
      void GetMaterialParams(const DRT::Element* ele,  //!< the element we are dealing with
          std::vector<double>& densn,                  //!< density at t_(n)
          std::vector<double>& densnp,                 //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,                 //!< density at t_(n+alpha_M)
          double& visc,                                //!< fluid viscosity
          const int iquad                              //!< id of current gauss point
          ) override;

      //! Clear all chemotaxtis related class variable (i.e. set them to zero)
      void ClearChemotaxisTerms();

      //! Get and save all chemotaxtis related class variable
      virtual void GetChemotaxisCoefficients(
          const Teuchos::RCP<const MAT::Material> material  //!< pointer to current material
      );

      //! Get ID of attractant (i.e. scalar which gradient the current scalar shall follow)
      virtual int GetPartner(const std::vector<int> pair);

      //! calculation of strong residual for stabilization
      void CalcStrongResidual(const int k,  //!< index of current scalar
          double& scatrares,                //!< residual of convection-diffusion-reaction eq
          const double densam,              //!< density at t_(n+am)
          const double densnp,              //!< density at t_(n+1)
          const double rea_phi,             //!< reactive contribution
          const double rhsint,              //!< rhs at integration point
          const double tau                  //!< the stabilisation parameter
          ) override;


      int numcondchemo_;                    //!< number of chemotactic conditions
      std::vector<std::vector<int>> pair_;  //!< vector containing the pairings
      std::vector<double> chemocoeff_;  //!< constants by which the chemotactic terms are multiplied

    };  // end Chemotaxis

  }  // namespace ELEMENTS

}  // namespace DRT


BACI_NAMESPACE_CLOSE

#endif  // SCATRA_ELE_CALC_CHEMO_H
