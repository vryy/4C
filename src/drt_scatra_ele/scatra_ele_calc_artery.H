/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_artery.H

 \brief evaluation class containing routines for calculation of scalar transport
        within 1D-arteries (blood vessels)
        only pressure-based formulation supports this

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_SCATRA_ELE_SCATRA_ELE_CALC_ARTERY_H_
#define SRC_DRT_SCATRA_ELE_SCATRA_ELE_CALC_ARTERY_H_

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "scatra_ele_calc.H"
#include "../drt_mat/cnst_1d_art.H"


namespace DRT
{
  namespace ELEMENTS
  {
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerArtery;

    template <DRT::Element::DiscretizationType distype, int probdim>
    class ScaTraEleCalcArtery : public ScaTraEleCalc<distype, probdim>
    {
     private:
     protected:
      //! (private) protected constructor, since we are a Singleton.
      ScaTraEleCalcArtery(const int numdofpernode, const int numscal, const std::string& disname);

     private:
     public:
      typedef ScaTraEleCalc<distype, probdim> my;

      virtual ~ScaTraEleCalcArtery() {}

      //! Singleton access method
      static ScaTraEleCalcArtery<distype, probdim>* Instance(const int numdofpernode,
          const int numscal, const std::string& disname,
          const ScaTraEleCalcArtery* delete_me = NULL);

      //! called upon destruction
      virtual void Done();

      /// Setup element evaluation
      virtual int SetupCalc(DRT::Element* ele, DRT::Discretization& discretization);

     protected:
      //! extract element based or nodal values
      //  return extracted values of phinp
      virtual void ExtractElementAndNodeValues(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la);

      //! evaluate shape functions and their derivatives at current integration point
      virtual double EvalShapeFuncAndDerivsAtIntPoint(
          const DRT::UTILS::IntPointsAndWeights<my::nsd_ele_>& intpoints,  //!< integration points
          const int iquad  //!< id of current Gauss point
      )
      {
        const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

        // scale fac with the area of the artery pi*D^2/4
        return fac * M_PI * VarManager()->Diam() * VarManager()->Diam() / 4.0;
      }

      //! evaluate shape functions and their derivatives at element center
      virtual double EvalShapeFuncAndDerivsAtEleCenter()
      {
        // use one-point Gauss rule to do calculations at the element center
        const DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints_tau(
            SCATRA::DisTypeToStabGaussRule<distype>::rule);

        // volume of the element (2D: element surface area; 1D: element length)
        // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
        // base class has to be called since we do not want the scaling with artery area here
        const double vol = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints_tau, 0);

        return vol;
      }

      //! evaluate material
      virtual void Materials(
          const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
          const int k,                                       //!< id of current scalar
          double& densn,                                     //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
      );

      //! calculation of convective element matrix in convective form (off diagonal term fluid)
      virtual void CalcMatConvODFluid(
          Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                     //!< index of current scalar
          const int ndofpernodefluid,      //!< number of dofs per node of fluid element
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double densnp,      //!< density at time_(n+1)
          const LINALG::Matrix<my::nsd_, 1>& gradphi  //!< scalar gradient
      );

      //! calculation of convective element matrix: add conservative contributions (off diagonal
      //! term fluid)
      virtual void CalcMatConvAddConsODFluid(
          Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                     //!< index of current scalar
          const int ndofpernodefluid,      //!< number of dofs per node of fluid element
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double densnp,      //!< density at time_(n+1)
          const double phinp        //!< scalar at time_(n+1)
      )
      {
        dserror("not possible");
        return;
      }

      //! get internal variable manager for multiporo formulation
      Teuchos::RCP<ScaTraEleInternalVariableManagerArtery<my::nsd_, my::nen_>> VarManager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleInternalVariableManagerArtery<my::nsd_, my::nen_>>(
            my::scatravarmanager_);
      };

      //! set internal variables
      virtual void SetInternalVariablesForMatAndRHS();

      //! artery pressure values at t_(n+1)
      LINALG::Matrix<my::nen_, 1> earterypressurenp_;
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

      virtual ~ScaTraEleInternalVariableManagerArtery(){};

      // compute and set internal variables -- no L2-projection but evaluation at GP
      void SetInternalVariablesArtery(
          const LINALG::Matrix<NEN, 1>& funct,  //! array for shape functions
          const LINALG::Matrix<NSD, NEN>&
              derxy,  //! global derivatives of shape functions w.r.t x,y,z
          const LINALG::Matrix<NSD, NEN>&
              deriv,  //! global derivatives of shape functions w.r.t r,s,t
          const LINALG::Matrix<NSD, NSD>& xjm,
          const std::vector<LINALG::Matrix<NEN, 1>>& ephinp,  //! scalar at t_(n+1) or t_(n+alpha_F)
          const std::vector<LINALG::Matrix<NEN, 1>>& ephin,   //! scalar at t_(n)
          const std::vector<LINALG::Matrix<NEN, 1>>&
              ehist,  //! history vector of transported scalars
          const LINALG::Matrix<NEN, 1>& earterypressure)
      {
        // call base class (scatra) with dummy variable econvelnp
        const LINALG::Matrix<NSD, NEN> econvelnp(true);
        my::SetInternalVariables(funct, derxy, ephinp, ephin, econvelnp, ehist);

        static LINALG::Matrix<NSD, 1> pressuregrad(true);
        pressuregrad.Multiply(derxy, earterypressure);

        // convective velocity
        my::convelint_.Update(-Diam() * Diam() / 32.0 / Visc(), pressuregrad, 0.0);
        // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
        my::conv_.MultiplyTN(derxy, my::convelint_);

        for (int k = 0; k < my::numscal_; ++k)
        {
          // overwrite convective term
          // - k/\mu*grad p * grad phi
          my::conv_phi_[k] = my::convelint_.Dot(my::gradphi_[k]);
        }

        return;
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



#endif /* SRC_DRT_SCATRA_ELE_SCATRA_ELE_CALC_ARTERY_H_ */
