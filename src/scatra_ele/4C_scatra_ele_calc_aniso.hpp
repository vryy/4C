/*----------------------------------------------------------------------*/
/*! \file

\brief scatra_ele_calc_aniso.H

\level 3

 *----------------------------------------------------------------------*/


#ifndef FOUR_C_SCATRA_ELE_CALC_ANISO_HPP
#define FOUR_C_SCATRA_ELE_CALC_ANISO_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    template <int nsd>
    class ScaTraEleDiffManagerAniso;

    template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype>>
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
      static ScaTraEleCalcAniso<distype, probdim>* instance(
          const int numdofpernode, const int numscal, const std::string& disname);

     protected:
      /*========================================================================*/
      //! @name material and related and related functions
      /*========================================================================*/

      //! get diffusion manager for anisotropic diffusivity
      Teuchos::RCP<ScaTraEleDiffManagerAniso<nsd_>> diff_manager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleDiffManagerAniso<nsd_>>(my::diffmanager_);
      };

      /*========================================================================*/
      //! @name methods for evaluation of individual terms
      /*========================================================================*/

      //! calculate the Laplacian for all shape functions(strong form)
      void get_laplacian_strong_form(
          Core::LinAlg::Matrix<nen_, 1>& diff  //!< laplace term to be computed
      );

      //! calculate the Laplacian (weak form)
      void get_laplacian_weak_form(double& val,                //!< ?
          const Core::LinAlg::Matrix<nsd_, nsd_>& difftensor,  //!< ?
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
      void get_laplacian_weak_form_rhs(double& val,            //!< ?
          const Core::LinAlg::Matrix<nsd_, nsd_>& difftensor,  //!< ?
          const Core::LinAlg::Matrix<nsd_, 1>& gradphi,        //!< ?
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
      void calc_mat_diff(Core::LinAlg::SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                           //!< index of current scalar
          const double timefacfac  //!< domain-integration factor times time-integration factor
          ) override;

      //! standard Galerkin diffusive term on right hand side
      void calc_rhs_diff(Core::LinAlg::SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                           //!< index of current scalar
          const double rhsfac  //!< time-integration factor for rhs times domain-integration factor
          ) override;
    };


    /// Scatra anisotropic diffusion manager
    /*!
    This is a class to handle anisotropic diffusion. It is derived from the
    basic isotropic diffusion manager class in scatra_ele_calc.H. It doesn't
    contain subgrid diffusion.
    */
    template <int nsd>
    class ScaTraEleDiffManagerAniso : public ScaTraEleDiffManager
    {
     public:
      ScaTraEleDiffManagerAniso(int numscal) : ScaTraEleDiffManager(numscal), difftensor_(numscal)
      {
        return;
      }

      //! Set the anisotropic diffusion coefficient
      virtual void set_anisotropic_diff(
          const Core::LinAlg::Matrix<nsd, nsd> difftensor, const int k)
      {
        difftensor_[k] = difftensor;
        return;
      }

      //! Return the stored anisotropic diffusion coefficient
      virtual Core::LinAlg::Matrix<nsd, nsd> get_anisotropic_diff(const int k)
      {
        return difftensor_[k];
      }

     protected:
      //! tensor valued diffusion coefficient
      std::vector<Core::LinAlg::Matrix<nsd, nsd>> difftensor_;
    };

  }  // namespace ELEMENTS
}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
