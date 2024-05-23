/*--------------------------------------------------------------------------*/
/*! \file

\brief Element evaluations for loma problems

       REMARK:
       Note, parameter "SGS_MATERIAL_UPDATE" is not anymore supported. The idea of
       this parameter is to update all material parameters (density, viscosity, conductivity, ...)
       using the temperature obtained from the solution, i.e., T^h, plus the subgrid-scale
       temperature (obtained from residual-based or multifractal subgrid-scales). This procedure has
originally been suggested by Avila et al. 2012. Elaborate testing of this option has not shown any
notable differences or improvements compared to the usual way, i.e., merely using T^h. For reuse of
this option see file scatra_ele_impl.cpp of revision 18836. If this option should be re-implemented,
it should not be combined with subgrid diffusivity approaches or artificial diffusion methods, since
a potential update of the material parameters may overwrite the subgrid diffusivity added to the
physical one. For avoiding this potential combination, it should be excluded in a loma parameter
list, where a FOUR_C_THROW() should be provided.

\level 2

*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_ELE_CALC_LOMA_HPP
#define FOUR_C_SCATRA_ELE_CALC_LOMA_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN


namespace DRT
{
  namespace ELEMENTS
  {
    template <CORE::FE::CellType distype>
    class ScaTraEleCalcLoma : public ScaTraEleCalc<distype>
    {
     private:
      //! private constructor for singletons
      ScaTraEleCalcLoma(const int numdofpernode, const int numscal, const std::string& disname);

      typedef ScaTraEleCalc<distype> my;
      using my::nen_;
      using my::nsd_;

     public:
      /// Singleton access method
      static ScaTraEleCalcLoma<distype>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

      int EvaluateAction(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const SCATRA::Action& action,
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) override;

     private:
      /*========================================================================*/
      //! @name material and related and related functions
      /*========================================================================*/

      //! evaluate material
      void Materials(
          const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          double& densn,                                           //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
          ) override;

      //! material mixfrac
      void MatMixFrac(
          const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          double& densn,                                           //!< density at t_(n)
          double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,  //!< density at t_(n+alpha_M)
          double& visc     //!< fluid viscosity
      );

      //! material Sutherland
      void MatSutherland(
          const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          double& densn,                                           //!< density at t_(n)
          double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,  //!< density at t_(n+alpha_M)
          double& visc     //!< fluid viscosity
      );

      //! material temperature-dependent water
      void MatTempDepWater(
          const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          double& densn,                                           //!< density at t_(n)
          double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,  //!< density at t_(n+alpha_M)
          double& visc     //!< fluid viscosity
      );

      //! material Arrhenius PV
      void MatArrheniusPV(
          const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          double& densn,                                           //!< density at t_(n)
          double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,  //!< density at t_(n+alpha_M)
          double& visc     //!< fluid viscosity
      );

      //! material Arrhenius Spec
      void MatArrheniusSpec(
          const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          double& densn,                                           //!< density at t_(n)
          double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,  //!< density at t_(n+alpha_M)
          double& visc     //!< fluid viscosity
      );

      //! material Arrhenius Temp
      void MatArrheniusTemp(
          const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          double& densn,                                           //!< density at t_(n)
          double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,  //!< density at t_(n+alpha_M)
          double& visc     //!< fluid viscosity
      );

      //! material Ferech PV
      void MatFerechPV(
          const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          double& densn,                                           //!< density at t_(n)
          double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,  //!< density at t_(n+alpha_M)
          double& visc     //!< fluid viscosity
      );

      //! material thermo St. Venant Kirchhoff
      void mat_thermo_st_venant_kirchhoff(
          const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          double& densn,                                           //!< density at t_(n)
          double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,  //!< density at t_(n+alpha_M)
          double& visc     //!< fluid viscosity
      );

      //! material Yoghurt
      void MatYoghurt(
          const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          double& densn,                                           //!< density at t_(n)
          double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,  //!< density at t_(n+alpha_M)
          double& visc     //!< fluid viscosity
      );

      /*========================================================================*/
      //! @name methods for evaluation of individual terms
      /*========================================================================*/

      //! compute rhs containing bodyforce
      void GetRhsInt(double& rhsint,  //!< rhs containing bodyforce at Gauss point
          const double densnp,        //!< density at t_(n+1)
          const int k                 //!< index of current scalar
          ) override;

      //! calculation of convective element matrix: add conservative contributions
      void CalcMatConvAddCons(
          CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                            //!< index of current scalar
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const CORE::LINALG::Matrix<nsd_, 1>& convelint,  //!< convective velocity at Gauss point
          const CORE::LINALG::Matrix<nsd_, 1>& gradphi,    //!< scalar gradient at Gauss point
          const double vdiv,                               //!< velocity divergence
          const double densnp,                             //!< density at time_(n+1)
          const double visc                                //!< viscosity
      );

      //! adaption of convective term for rhs
      void recompute_conv_phi_for_rhs(double& conv_phi,   //!< convective contribution
          const int k,                                    //!< index of current scalar
          const CORE::LINALG::Matrix<nsd_, 1>& sgvelint,  //!< subgrid-scale velocity at Gauss point
          const CORE::LINALG::Matrix<nsd_, 1>& gradphi,   //!< scalar gradient at Gauss point
          const double densnp,                            //!< density at time_(n+1)
          const double densn,                             //!< density at time_(n)
          const double phinp,                             //!< scalar at time_(n+1)
          const double phin,                              //!< scalar at time_(n)
          const double vdiv                               //!< velocity divergence
      );

      /*========================================================================*/
      //! @name additional service routines
      /*========================================================================*/

      //! calculate domain integral
      void calculate_domain_and_bodyforce(
          CORE::LINALG::SerialDenseVector& scalars, const DRT::Element* ele);

      //! extract element based or nodal values and return extracted values of phinp
      void extract_element_and_node_values(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la) override;

      //! get density at integration point
      double GetDensity(const DRT::Element* ele, Teuchos::RCP<const CORE::MAT::Material> material,
          Teuchos::ParameterList& params, const double tempnp) override;

      //! calculate viscous part of subgrid-scale velocity
      void calc_subgr_velocity_visc(CORE::LINALG::Matrix<nsd_, 1>& epsilonvel) override;

      /*========================================================================*/
      //! @name scalar degrees of freedom and related
      /*========================================================================*/

      //! scalar at t_(n+alpha_M)
      std::vector<CORE::LINALG::Matrix<nen_, 1>> ephiam_;
      //!
      std::vector<double> densgradfac_;

      //! thermodynamic pressure at t_(n+1) or t_(n+alpha_F) (LOMA specific)
      double thermpressnp_;
      //! thermodynamic pressure at t_(n+alpha_M) (LOMA specific)
      double thermpressam_;
      //! time derivative of thermodynamic pressure (LOMA specific)
      double thermpressdt_;

      //! specific heat capacity (either at constant pressure or at constant volume)
      double shc_;
    };

    /// Scatra reaction manager
    /*!
        special reaction manager for loma
    */
    class ScaTraEleReaManagerLoma : public ScaTraEleReaManager
    {
     public:
      ScaTraEleReaManagerLoma(int numscal) : ScaTraEleReaManager(numscal), reatemprhs_(numscal, 0.0)
      {
        return;
      }

      //! Set reaction / temperature term for rhs
      void SetReaTempRhs(const double reatemprhs, const int k)
      {
        reatemprhs_[k] = reatemprhs;
        return;
      }

      //! Return reaction / temperature term for rhs
      double GetReaTempRhs(const int k) { return reatemprhs_[k]; }

     private:
      //! reaction / temperature term for rhs
      std::vector<double> reatemprhs_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
