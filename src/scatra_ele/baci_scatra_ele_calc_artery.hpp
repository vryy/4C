/*----------------------------------------------------------------------*/
/*! \file
 \brief evaluation class containing routines for calculation of scalar transport
        within 1D-arteries (blood vessels)
        only pressure-based formulation supports this

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_ELE_CALC_ARTERY_HPP
#define FOUR_C_SCATRA_ELE_CALC_ARTERY_HPP

#include "baci_config.hpp"

#include "baci_lib_discret.hpp"
#include "baci_mat_cnst_1d_art.hpp"
#include "baci_scatra_ele_calc.hpp"

BACI_NAMESPACE_OPEN


namespace DRT
{
  namespace ELEMENTS
  {
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerArtery;

    template <CORE::FE::CellType distype, int probdim>
    class ScaTraEleCalcArtery : public ScaTraEleCalc<distype, probdim>
    {
     private:
     protected:
      //! (private) protected constructor, since we are a Singleton.
      ScaTraEleCalcArtery(const int numdofpernode, const int numscal, const std::string& disname);

     private:
     public:
      typedef ScaTraEleCalc<distype, probdim> my;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;


      //! Singleton access method
      static ScaTraEleCalcArtery<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

      /// Setup element evaluation
      int SetupCalc(DRT::Element* ele, DRT::Discretization& discretization) override;

     protected:
      //! extract element based or nodal values
      //  return extracted values of phinp
      void ExtractElementAndNodeValues(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la) override;

      //! evaluate shape functions and their derivatives at current integration point
      double EvalShapeFuncAndDerivsAtIntPoint(
          const CORE::FE::IntPointsAndWeights<nsd_ele_>& intpoints,  //!< integration points
          const int iquad                                            //!< id of current Gauss point
          ) override
      {
        const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

        // scale fac with the area of the artery pi*D^2/4
        return fac * M_PI * VarManager()->Diam() * VarManager()->Diam() / 4.0;
      }

      //! evaluate shape functions and their derivatives at element center
      double EvalShapeFuncAndDerivsAtEleCenter() override
      {
        // use one-point Gauss rule to do calculations at the element center
        const CORE::FE::IntPointsAndWeights<nsd_ele_> intpoints_tau(
            SCATRA::DisTypeToStabGaussRule<distype>::rule);

        // volume of the element (2D: element surface area; 1D: element length)
        // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
        // base class has to be called since we do not want the scaling with artery area here
        const double vol = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints_tau, 0);

        return vol;
      }

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

      //! calculation of convective element matrix in convective form (off diagonal term fluid)
      void CalcMatConvODFluid(
          CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                            //!< index of current scalar
          const int ndofpernodefluid,             //!< number of dofs per node of fluid element
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double densnp,      //!< density at time_(n+1)
          const CORE::LINALG::Matrix<nsd_, 1>& gradphi  //!< scalar gradient
          ) override;

      //! calculation of convective element matrix: add conservative contributions (off diagonal
      //! term fluid)
      void CalcMatConvAddConsODFluid(
          CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                            //!< index of current scalar
          const int ndofpernodefluid,             //!< number of dofs per node of fluid element
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double densnp,      //!< density at time_(n+1)
          const double phinp        //!< scalar at time_(n+1)
          ) override
      {
        dserror("not possible");
        return;
      }

      //! get internal variable manager for multiporo formulation
      Teuchos::RCP<ScaTraEleInternalVariableManagerArtery<nsd_, nen_>> VarManager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleInternalVariableManagerArtery<nsd_, nen_>>(
            my::scatravarmanager_);
      };

      //! set internal variables
      void SetInternalVariablesForMatAndRHS() override;

      //! artery pressure values at t_(n+1)
      CORE::LINALG::Matrix<nen_, 1> earterypressurenp_;
    };

    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerArtery : public ScaTraEleInternalVariableManager<NSD, NEN>
    {
      typedef ScaTraEleInternalVariableManager<NSD, NEN> my;

     public:
      ScaTraEleInternalVariableManagerArtery(int numscal)
          : ScaTraEleInternalVariableManager<NSD, NEN>(numscal), materialset_(false)
      {
        return;
      }

      // compute and set internal variables -- no L2-projection but evaluation at GP
      void SetInternalVariablesArtery(
          const CORE::LINALG::Matrix<NEN, 1>& funct,  //! array for shape functions
          const CORE::LINALG::Matrix<NSD, NEN>&
              derxy,  //! global derivatives of shape functions w.r.t x,y,z
          const CORE::LINALG::Matrix<NSD, NEN>&
              deriv,  //! global derivatives of shape functions w.r.t r,s,t
          const CORE::LINALG::Matrix<NSD, NSD>& xjm,
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>&
              ephinp,  //! scalar at t_(n+1) or t_(n+alpha_F)
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>& ephin,  //! scalar at t_(n)
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>&
              ehist,  //! history vector of transported scalars
          const CORE::LINALG::Matrix<NEN, 1>& earterypressure)
      {
        // call base class (scatra) with dummy variable econvelnp
        const CORE::LINALG::Matrix<NSD, NEN> econvelnp(true);
        const CORE::LINALG::Matrix<NSD, NEN> eforcevelocity(true);
        my::SetInternalVariables(funct, derxy, ephinp, ephin, econvelnp, ehist, eforcevelocity);

        static CORE::LINALG::Matrix<NSD, 1> pressuregrad(true);
        pressuregrad.Multiply(derxy, earterypressure);

        for (int k = 0; k < my::numscal_; ++k)
        {
          // convective velocity
          my::convelint_[k].Update(-Diam() * Diam() / 32.0 / Visc(), pressuregrad, 0.0);
          // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
          my::conv_[k].MultiplyTN(derxy, my::convelint_[k]);
          // overwrite convective term
          // - k/\mu*grad p * grad phi
          my::conv_phi_[k] = my::convelint_[k].Dot(my::gradphi_[k]);
        }
      };

      // Set the artery material in the scatra-Varmanager
      void SetArteryMaterial(DRT::Element* ele)
      {
        // check if we actually have two materials
        if (ele->NumMaterial() < 2) dserror("no second material available");
        // check for artery material
        if (ele->Material(1)->MaterialType() != INPAR::MAT::MaterialType::m_cnst_art)
          dserror("Secondary material is not of type m_cnst_art, but %d",
              ele->Material(1)->MaterialType());

        // here we rely that the Artery material has been added as second material
        arterymat_ = Teuchos::rcp_dynamic_cast<MAT::Cnst_1d_art>(ele->Material(1));

        materialset_ = true;
      }

      // return density
      double Dens() { return ArteryMat()->Density(); }

      // return diameter
      double Diam() { return ArteryMat()->Diam(); }

      // return viscosity
      double Visc() { return ArteryMat()->Viscosity(); }

      //! return artery material
      Teuchos::RCP<MAT::Cnst_1d_art> ArteryMat()
      {
        if (!materialset_) dserror("Artery Material has not yet been set in Variablemanager");

        return arterymat_;
      }

     private:
      //! artery material
      Teuchos::RCP<MAT::Cnst_1d_art> arterymat_;

      //! check if artery material has been set
      bool materialset_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT



BACI_NAMESPACE_CLOSE

#endif
