/*----------------------------------------------------------------------*/
/*! \file

\brief scatra_ele_calc_aniso.H

\level 3

 *----------------------------------------------------------------------*/


#ifndef FOUR_C_SCATRA_ELE_CALC_ANISO_HPP
#define FOUR_C_SCATRA_ELE_CALC_ANISO_HPP

#include "baci_config.hpp"

#include "baci_scatra_ele_calc.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    template <int NSD>
    class ScaTraEleDiffManagerAniso;

    template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype>>
    class ScaTraEleCalcAniso : public virtual ScaTraEleCalc<distype, probdim>
    {
     protected:
      /// (private) protected constructor, since we are a Singleton.
      /// this constructor is called from a derived class
      /// -> therefore, it has to be protected instead of private
      ScaTraEleCalcAniso(const int numdofpernode, const int numscal, const std::string& disname);

     private:
      typedef ScaTraEleCalc<distype, probdim> my;

     protected:
      using my::nen_;
      using my::nsd_;

     private:
      using varmanager = ScaTraEleInternalVariableManager<nsd_, nen_>;

     public:
      /// Singleton access method
      static ScaTraEleCalcAniso<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

     protected:
      /*========================================================================*/
      //! @name material and related and related functions
      /*========================================================================*/

      //! evaluate material
      void Materials(
          const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
          const int k,                                       //!< id of current scalar
          double& densn,                                     //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
          ) override;

      //! material ScaTra
      void MatScaTraAniso(
          const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
          const int k,                                       //!< id of current scalar
          double& densn,                                     //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
      );

      //! get diffusion manager for anisotropic diffusivity
      Teuchos::RCP<ScaTraEleDiffManagerAniso<nsd_>> DiffManager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleDiffManagerAniso<nsd_>>(my::diffmanager_);
      };

      /*========================================================================*/
      //! @name methods for evaluation of individual terms
      /*========================================================================*/

      //! calculate the Laplacian for all shape functions(strong form)
      void GetLaplacianStrongForm(
          CORE::LINALG::Matrix<nen_, 1>& diff  //!< laplace term to be computed
      );

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

      //! calculation of diffusive element matrix
      void CalcMatDiff(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                         //!< index of current scalar
          const double timefacfac  //!< domain-integration factor times time-integration factor
          ) override;

      //! standard Galerkin diffusive term on right hand side
      void CalcRHSDiff(CORE::LINALG::SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                         //!< index of current scalar
          const double rhsfac  //!< time-integration factor for rhs times domain-integration factor
          ) override;
    };


    /// Scatra anisotropic diffusion manager
    /*!
    This is a class to handle anisotropic diffusion. It is derived from the
    basic isotropic diffusion manager class in scatra_ele_calc.H. It doesn't
    contain subgrid diffusion.
    */
    template <int NSD>
    class ScaTraEleDiffManagerAniso : public ScaTraEleDiffManager
    {
     public:
      ScaTraEleDiffManagerAniso(int numscal) : ScaTraEleDiffManager(numscal), difftensor_(numscal)
      {
        return;
      }

      //! Set the anisotropic diffusion coefficient
      virtual void SetAnisotropicDiff(const CORE::LINALG::Matrix<NSD, NSD> difftensor, const int k)
      {
        difftensor_[k] = difftensor;
        return;
      }

      //! Return the stored anisotropic diffusion coefficient
      virtual CORE::LINALG::Matrix<NSD, NSD> GetAnisotropicDiff(const int k)
      {
        return difftensor_[k];
      }

     protected:
      //! tensor valued diffusion coefficient
      std::vector<CORE::LINALG::Matrix<NSD, NSD>> difftensor_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT


BACI_NAMESPACE_CLOSE

#endif  // SCATRA_ELE_CALC_ANISO_H
