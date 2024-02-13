/*----------------------------------------------------------------------*/
/*! \file
\brief main file containing routines for calculation of scatra element formulated in reference
concentrations and with advanced reaction terms

\level 3

 *----------------------------------------------------------------------*/


#ifndef BACI_SCATRA_ELE_CALC_REFCONC_REAC_HPP
#define BACI_SCATRA_ELE_CALC_REFCONC_REAC_HPP

#include "baci_config.hpp"

#include "baci_scatra_ele_calc_advanced_reaction.hpp"

BACI_NAMESPACE_OPEN


namespace DRT
{
  namespace ELEMENTS
  {
    template <CORE::FE::CellType distype>
    class ScaTraEleCalcRefConcReac : public ScaTraEleCalcAdvReac<distype>
    {
     private:
      /// private constructor, since we are a Singleton.
      ScaTraEleCalcRefConcReac(
          const int numdofpernode, const int numscal, const std::string& disname);

      typedef ScaTraEleCalc<distype> my;
      typedef ScaTraEleCalcAdvReac<distype> advreac;
      using my::nen_;
      using my::nsd_;

     public:
      /// Singleton access method
      static ScaTraEleCalcRefConcReac<distype>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

     protected:
      //! Set reac. body force, reaction coefficient and derivatives
      void SetAdvancedReactionTerms(const int k,                  //!< index of current scalar
          const Teuchos::RCP<MAT::MatListReactions> matreaclist,  //!< index of current scalar
          const double* gpcoord  //!< current Gauss-point coordinates
          ) override;

      //! calculation of convective element matrix: add conservative contributions
      void CalcMatConvAddCons(
          CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                            //!< index of current scalar
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double vdiv,        //!< velocity divergence
          const double densnp       //!< density at time_(n+1)
          ) override;

      //! set internal variables
      void SetInternalVariablesForMatAndRHS() override;

      //! calculation of diffusive element matrix
      void CalcMatDiff(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                         //!< index of current scalar
          const double timefacfac  //!< domain-integration factor times time-integration factor
          ) override;

      //! calculate the Laplacian (weak form)
      void GetLaplacianWeakForm(double& val,                   //!< ?
          const CORE::LINALG::Matrix<nsd_, nsd_>& difftensor,  //!< ?
          const int vi,                                        //!< ?
          const int ui                                         //!< ?
      )
      {
        val = 0.0;
        for (unsigned j = 0; j < nsd_; j++)
        {
          for (unsigned i = 0; i < nsd_; i++)
          {
            val += my::derxy_(j, vi) * difftensor(j, i) * my::derxy_(i, ui);
          }
        }
        return;
      };

      //! standard Galerkin diffusive term on right hand side
      void CalcRHSDiff(CORE::LINALG::SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                         //!< index of current scalar
          const double rhsfac  //!< time-integration factor for rhs times domain-integration factor
          ) override;

      //! calculate the Laplacian (weak form)
      void GetLaplacianWeakFormRHS(double& val,                //!< ?
          const CORE::LINALG::Matrix<nsd_, nsd_>& difftensor,  //!< ?
          const CORE::LINALG::Matrix<nsd_, 1>& gradphi,        //!< ?
          const int vi                                         //!< ?
      )
      {
        val = 0.0;
        for (unsigned j = 0; j < nsd_; j++)
        {
          for (unsigned i = 0; i < nsd_; i++)
          {
            val += my::derxy_(j, vi) * difftensor(j, i) * gradphi(i);
          }
        }
        return;
      };

      // add nodal displacements to point coordinates
      void UpdateNodeCoordinates() override
      { /*nothing to to since we want reference coordinates*/
        return;
      };

     private:
      /// determinante of deformation gradient
      double J_;

      /// inverse of cauchy-green deformation gradient
      CORE::LINALG::Matrix<nsd_, nsd_> C_inv_;

      /// derivative dJ/dX by finite differences
      CORE::LINALG::Matrix<nsd_, 1> dJdX_;

    };  // end ScaTraEleCalcRefConcReac


  }  // namespace ELEMENTS

}  // namespace DRT


BACI_NAMESPACE_CLOSE

#endif  // SCATRA_ELE_CALC_REFCONC_REAC_H
