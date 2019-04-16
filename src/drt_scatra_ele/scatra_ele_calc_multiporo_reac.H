/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_multiporo_reac.H

 \brief evaluation class containing routines for calculation of scalar transport
        within multiphase porous medium

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_SCATRA_ELE_SCATRA_ELE_CALC_MULTIPORO_REAC_H_
#define SRC_DRT_SCATRA_ELE_SCATRA_ELE_CALC_MULTIPORO_REAC_H_


#include "scatra_ele_calc_poro_reac.H"
#include "../drt_mat/fluidporo_multiphase.H"
#include "../drt_mat/fluidporo_singlephase.H"
#include "../drt_mat/fluidporo_multiphase_reactions.H"
#include "../drt_mat/scatra_mat_multiporo.H"
#include "../drt_porofluidmultiphase_ele/porofluidmultiphase_ele_calc_utils.H"

#include "../drt_porofluidmultiphase_ele/porofluid_phasemanager.H"
#include "../drt_porofluidmultiphase_ele/porofluidmultiphase_ele_action.H"
#include "../drt_porofluidmultiphase_ele/porofluidmultiphase_ele_calc.H"
#include "../drt_porofluidmultiphase_ele/porofluid_variablemanager.H"
#include "../drt_porofluidmultiphase_ele/porofluid_evaluator.H"
#include "../drt_porofluidmultiphase_ele/porofluidmultiphase_ele_parameter.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

// forward declarations
namespace MAT
{
  class ScatraMat;
}

namespace DRT
{
  namespace ELEMENTS
  {
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerMultiPoro;

    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleCalcMultiPoroReac : public ScaTraEleCalcPoroReac<distype>
    {
     private:
      /// private constructor, since we are a Singleton.
      ScaTraEleCalcMultiPoroReac(
          const int numdofpernode, const int numscal, const std::string& disname);

      typedef ScaTraEleCalc<distype> my;
      typedef ScaTraEleCalcPoroReac<distype> pororeac;
      typedef ScaTraEleCalcPoro<distype> poro;
      typedef ScaTraEleCalcAdvReac<distype> advreac;

     public:
      virtual ~ScaTraEleCalcMultiPoroReac() {}

      /// Singleton access method
      static ScaTraEleCalcMultiPoroReac<distype>* Instance(const int numdofpernode,
          const int numscal, const std::string& disname,
          const ScaTraEleCalcMultiPoroReac* delete_me = NULL);

      /// called upon destruction
      virtual void Done();

      /// Setup element evaluation
      virtual int SetupCalc(DRT::Element* ele, DRT::Discretization& discretization);

     protected:
      //! extract element based or nodal values
      //  return extracted values of phinp
      virtual void ExtractElementAndNodeValues(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la);

      //! extract element based or nodal values --> L2-projection case: called within
      //! ExtractElementAndNodeValues
      //  return extracted values of phinp
      virtual void ExtractNodalFlux(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          const int numfluidphases);

      //! set internal variables
      virtual void SetInternalVariablesForMatAndRHS();

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

      //! material MatMultiPoroFluid
      virtual void MatMultiPoroFluid(
          const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
          const int k,                                       //!< id of current scalar
          double& densn,                                     //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
      );

      //! material MatMultiPoroVolFrac
      virtual void MatMultiPoroVolFrac(
          const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
          const int k,                                       //!< id of current scalar
          double& densn,                                     //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
      );

      //! material MatMultiPoroSolid
      virtual void MatMultiPoroSolid(
          const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
          const int k,                                       //!< id of current scalar
          double& densn,                                     //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
      );

      //! material MatMultiPoroTemperature
      virtual void MatMultiPoroTemperature(
          const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
          const int k,                                       //!< id of current scalar
          double& densn,                                     //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
      );

      //! Set advanced reaction terms and derivatives
      virtual void SetAdvancedReactionTerms(const int k,          //!< index of current scalar
          const Teuchos::RCP<MAT::MatListReactions> matreaclist,  //!< index of current scalar
          const double* gpcoord  //!< current Gauss-point coordinates
      );

      //! compute pore pressure
      virtual double ComputePorePressure();

      //! get internal variable manager for multiporo formulation
      Teuchos::RCP<ScaTraEleInternalVariableManagerMultiPoro<my::nsd_, my::nen_>> VarManager()
      {
        return Teuchos::rcp_static_cast<
            ScaTraEleInternalVariableManagerMultiPoro<my::nsd_, my::nen_>>(my::scatravarmanager_);
      };

      //! calculation of convective element matrix in convective form
      //! the only difference to the base class version is, that there is no scaling with the
      //! density
      virtual void CalcMatConv(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                          //!< index of current scalar
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double densnp,      //!< density at time_(n+1)
          const LINALG::Matrix<my::nen_, 1>& sgconv  //!< subgrid-scale convective operator
      );

      //! adaption of convective term for rhs
      //! the only difference to the base class version is, that there is no scaling with the
      //! density
      virtual void RecomputeConvPhiForRhs(const int k,  //!< index of current scalar
          const LINALG::Matrix<my::nsd_, 1>& sgvelint,  //!< subgrid-scale velocity at Gauss point
          const double densnp,                          //!< density at time_(n+1)
          const double densn,                           //!< density at time_(n)
          const double vdiv                             //!< velocity divergence
      );

      //! calculation of convective element matrix: add conservative contributions
      //! the only difference to the base class version is, that there is no scaling with the
      //! density
      virtual void CalcMatConvAddCons(
          Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                     //!< index of current scalar
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double vdiv,        //!< velocity divergence
          const double densnp       //!< density at time_(n+1)
      );

      //! calculation of mass element matrix (standard shape functions)
      virtual void CalcMatMass(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int& k,                                         //!< index of current scalar
          const double& fac,                                    //!< domain-integration factor
          const double& densam                                  //!< density at time_(n+am)
      );

      //! calculation of convective element matrix (OD term structure coupling)
      //! difference to base class: linearization of mesh motion + shapederivatives pressure
      //! gradient have to be included
      virtual void CalcConvODMesh(Epetra_SerialDenseMatrix& emat, const int k,
          const int ndofpernodemesh, const double fac, const double rhsfac, const double densnp,
          const double J, const LINALG::Matrix<my::nsd_, 1>& gradphi,
          const LINALG::Matrix<my::nsd_, 1>& convelint);

      //! calculation of linearized mass (off diagonal/shapederivative term mesh)
      //! difference to base class: linearization of porosity is included
      virtual void CalcLinMassODMesh(
          Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                     //!< index of current scalar
          const int ndofpernodemesh,       //!< number of dofs per node of ale element
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double fac,     //!< domain-integration factor
          const double densam,  //!< density at time_(n+am)
          const double densnp,  //!< density at time_(n+1)
          const double phinp,   //!< scalar at time_(n+1)
          const double hist,    //!< history of time integartion
          const double J,       //!< determinant of Jacobian det(dx/ds)
          const LINALG::Matrix<1, my::nsd_ * my::nen_>&
              dJ_dmesh  //!< derivative of det(dx/ds) w.r.t. mesh displacement
      );

      //! standard Galerkin transient, old part of rhs and source term (off diagonal/shapederivative
      //! term mesh) difference to base class: linearization of porosity and advanced reaction terms
      //! are included
      virtual void CalcHistAndSourceODMesh(
          Epetra_SerialDenseMatrix& emat,  //!< element current to be filled
          const int k,                     //!< index of current scalar
          const int ndofpernodemesh,       //!< number of dofs per node of ale element
          const double fac,                //!< domain-integration factor
          const double rhsint,             //!< rhs at Gauss point
          const double J,                  //!< determinant of Jacobian det(dx/ds)
          const LINALG::Matrix<1, my::nsd_ * my::nen_>&
              dJ_dmesh,        //!< derivative of det(dx/ds) w.r.t. mesh displacement
          const double densnp  //!< density
      );

      //! standard Galerkin diffusive term (off diagonal/shapederivative term mesh)
      //! difference to base class: linearization of porosity and effective diffusivity
      virtual void CalcDiffODMesh(Epetra_SerialDenseMatrix& emat,  //!< element current to be filled
          const int k,                                             //!< index of current scalar
          const int ndofpernodemesh,  //!< number of dofs per node of ale element
          const double diffcoeff,     //!< diffusion coefficient
          const double fac,           //!< domain-integration factor
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double J,       //!< determinant of Jacobian det(dx/ds)
          const LINALG::Matrix<my::nsd_, 1>& gradphi,    //!< scalar gradient at Gauss point
          const LINALG::Matrix<my::nsd_, 1>& convelint,  //!< convective velocity
          const LINALG::Matrix<1, my::nsd_ * my::nen_>&
              dJ_dmesh  //!< derivative of det(dx/ds) w.r.t. mesh displacement
      );

      //! reactive terms (standard Galerkin) (off diagonal/shapederivative term mesh)
      //! difference to base class: linearization of porosity
      virtual void CalcReactODMesh(
          Epetra_SerialDenseMatrix& emat,  //!< element current to be filled
          const int k,                     //!< index of current scalar
          const int ndofpernodemesh,       //!< number of dofs per node of ale element
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double rea_phi,  //!< reactive term
          const double J,        //!< determinant of Jacobian det(dx/ds)
          const LINALG::Matrix<1, my::nsd_ * my::nen_>&
              dJ_dmesh  //!< derivative of det(dx/ds) w.r.t. mesh displacement
      );

      //! calculation of convective element matrix in convective form (off diagonal term fluid)
      virtual void CalcMatConvODFluid(
          Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                     //!< index of current scalar
          const int ndofpernodefluid,  //!< number of dofs per node of fluid element // only a dummy
                                       //!< variable
          const double timefacfac,     //!< domain-integration factor times time-integration factor
          const double densnp,         //!< density at time_(n+1)
          const LINALG::Matrix<my::nsd_, 1>& gradphi  //!< scalar gradient
      );

      //! calculation of convective element matrix in convective form -- additional conservative
      //! contributions (off diagonal term fluid) not yet implemented --> dserror
      virtual void CalcMatConvAddConsODFluid(
          Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                     //!< index of current scalar
          const int ndofpernodefluid,      //!< number of dofs per node of fluid element
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double densnp,      //!< density at time_(n+1)
          const double phinp        //!< scalar at time_(n+1)
      );

      //! calculation of linearized mass (off diagonal terms fluid)
      //! linearization of porosity*saturation*vtrans
      virtual void CalcLinMassODFluid(
          Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                     //!< index of current scalar
          const int ndofpernodemesh,  //!< number of dofs per node of fluid element // only a dummy
                                      //!< variable
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double fac,     //!< domain-integration factor
          const double densam,  //!< density at time_(n+am)
          const double densnp,  //!< density at time_(n+1)
          const double phinp,   //!< scalar at time_(n+1)
          const double hist     //!< history of time integartion
      );

      // calculate linearization of a mass matrix type matrix (OD-fluid terms)
      virtual void CalcLinMassMatrixTypeODFluid(Epetra_SerialDenseMatrix& emat, const int k,
          const std::vector<double>* prefaclinmassodfluid, const int totalnummultiphasedofpernode,
          double prefac);

      //! standard Galerkin transient, old part of rhs and source term (off diagonal terms fluid)
      //! linearization of porosity*saturation*vtrans + advanced reaction terms
      virtual void CalcHistAndSourceODFluid(
          Epetra_SerialDenseMatrix& emat,  //!< element current to be filled
          const int k,                     //!< index of current scalar
          const int ndofpernodemesh,       //!< number of dofs per node of ale element
          const double fac,                //!< domain-integration factor
          const double rhsint,             //!< rhs at Gauss point
          const double densnp              //!< density
      );

      //! standard Galerkin reactive term (off diagonal terms fluid)
      //! linearization of porosity*saturation*vreact
      virtual void CalcReactODFluid(
          Epetra_SerialDenseMatrix& emat,  //!< element current to be filled
          const int k,                     //!< index of current scalar
          const int ndofpernodemesh,       //!< number of dofs per node of ale element
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double rea_phi  //!< rhs at Gauss point
      );

      //! standard Galerkin diffusive term (off diagonal terms fluid)
      //! linearization of porosity*saturation*d_eff
      virtual void CalcDiffODFluid(
          Epetra_SerialDenseMatrix& emat,  //!< element current to be filled
          const int k,                     //!< index of current scalar
          const int ndofpernodemesh,       //!< number of dofs per node of ale element
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const LINALG::Matrix<my::nsd_, 1>& gradphi  //!< scalar gradient at Gauss point
      );

      //! fill coupling vector and add variables to reaction in order to compute reaction values and
      //! derivatives
      void FillCouplingVectorAndAddVariables(const int k,
          const Teuchos::RCP<MAT::MatListReactions> matreaclist,
          const Teuchos::RCP<ScaTraEleReaManagerAdvReac> remanager);

      //! Get right hand side including reaction bodyforce term
      virtual void GetRhsInt(double& rhsint,  //!< rhs containing bodyforce at Gauss point
          const double densnp,                //!< density at t_(n+1)
          const int k                         //!< index of current scalar
      );

      //! calculation of reactive element matrix
      virtual void CalcMatReact(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                           //!< index of current scalar
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double
              timetaufac,  //!< domain-integration factor times time-integration factor times tau
          const double taufac,                        //!< domain-integration factor times tau
          const double densnp,                        //!< density at time_(n+1)
          const LINALG::Matrix<my::nen_, 1>& sgconv,  //!< subgrid-scale convective operator
          const LINALG::Matrix<my::nen_, 1>& diff     //!< laplace term
      );

      //! "shapederivatives" pressure gradient
      virtual void ApplyShapeDerivsPressureGrad(Epetra_SerialDenseMatrix& emat, const int k,
          const double vrhs, const LINALG::Matrix<my::nsd_, 1>& gradphi,
          const LINALG::Matrix<my::nsd_, 1> refgradpres);

     protected:
      //! nodal flux values at t_(n+1)
      std::vector<LINALG::Matrix<my::nsd_, my::nen_>> efluxnp_;

      //! a vector containing all quantities, the equation is coupled with
      //! (i.e. pressures, saturations and porosity)
      std::vector<std::pair<std::string, double>> couplingvalues_;

      //! penalty factor to avoid very small "densities"
      const double penalty_ = 1.0;

      //! do we use L2-projection or evaluation at GP
      bool L2_projection_ = false;
    };


    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerMultiPoro
        : public ScaTraEleInternalVariableManager<NSD, NEN>
    {
      typedef ScaTraEleInternalVariableManager<NSD, NEN> my;

     public:
      ScaTraEleInternalVariableManagerMultiPoro(int numscal)
          : ScaTraEleInternalVariableManager<NSD, NEN>(numscal),
            multiconvelint_(0),
            multiconv_(0),
            zeroconvelint_(true),
            zeroconv_(true),
            temperatureconvelint_(true),
            temperatureconv_(true),
            pressure_(0),
            saturation_(0),
            density_(0),
            solidpressure_(0.0),
            delta_(numscal, 0.0),
            heatcapacity_(0),
            thermaldiffusivity_(0),
            heatcapacityeff_(0.0),
            min_val_of_phase_(numscal, 0.0),
            evaluate_scalar_(numscal, true),
            scalartophasemap_(
                numscal, {-1, MAT::ScatraMatMultiPoro::SpeciesType::species_undefined}),
            materialset_(false),
            myaction_(SCATRA::Action::calc_mat_and_rhs)
      {
        return;
      }

      virtual ~ScaTraEleInternalVariableManagerMultiPoro(){};

      // compute and set internal variables -- no L2-projection but evaluation at GP
      void SetInternalVariablesMultiPoro(
          const LINALG::Matrix<NEN, 1>& funct,  //! array for shape functions
          const LINALG::Matrix<NSD, NEN>&
              derxy,  //! global derivatives of shape functions w.r.t x,y,z
          const LINALG::Matrix<NSD, NEN>&
              deriv,  //! global derivatives of shape functions w.r.t r,s,t
          const LINALG::Matrix<NSD, NSD>& xjm, const LINALG::Matrix<NSD, NEN>& xyze0,
          const std::vector<LINALG::Matrix<NEN, 1>>& ephinp,  //! scalar at t_(n+1) or t_(n+alpha_F)
          const std::vector<LINALG::Matrix<NEN, 1>>& ephin,   //! scalar at t_(n)
          const std::vector<LINALG::Matrix<NEN, 1>>&
              ehist  //! history vector of transported scalars
      )
      {
        // call base class (scatra) with dummy variable econvelnp
        const LINALG::Matrix<NSD, NEN> econvelnp(true);
        my::SetInternalVariables(funct, derxy, ephinp, ephin, econvelnp, ehist);

        //------------------------get determinant of Jacobian dX / ds
        // transposed jacobian "dX/ds"
        LINALG::Matrix<NSD, NSD> xjm0;
        xjm0.MultiplyNT(deriv, xyze0);

        // inverse of transposed jacobian "ds/dX"
        const double det0 = xjm0.Determinant();

        const double det = xjm.Determinant();

        // determinant of deformationgradient det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds)
        // )^-1
        const double JacobianDefGradient = det / det0;

        // clear current gauss point data for safety
        phasemanager_->ClearGPState();
        variablemanager_->EvaluateGPVariables(funct, derxy);
        // access from outside to the phasemanager: scatra-discretization has fluid-dis on dofset 2
        phasemanager_->EvaluateGPState(
            JacobianDefGradient, *variablemanager_, ndsscatra_porofluid_);

        const int numfluidphases = phasemanager_->NumFluidPhases();
        const int numvolfrac = phasemanager_->NumVolFrac();

        // resize all phase related vectors
        multiconvelint_.resize(numfluidphases + numvolfrac);
        multiconv_.resize(numfluidphases + numvolfrac);
        pressure_.resize(numfluidphases);
        saturation_.resize(numfluidphases);
        density_.resize(numfluidphases + numvolfrac);
        heatcapacity_.resize(numfluidphases + numvolfrac + 1);
        thermaldiffusivity_.resize(numfluidphases + numvolfrac + 1);
        pressuregrad_.resize(numfluidphases + numvolfrac);
        difftensorsfluid_.resize(numfluidphases + numvolfrac);
        abspressuregrad_.resize(numfluidphases);
        volfrac_.resize(numvolfrac);
        volfracpressure_.resize(numvolfrac);

        const std::vector<LINALG::Matrix<NSD, 1>>& fluidgradphi = *(variablemanager_->GradPhinp());

        volfrac_ = phasemanager_->VolFrac();
        volfracpressure_ = phasemanager_->VolFracPressure();

        temperatureconv_ = 0.0;
        temperatureconvelint_ = 0.0;

        for (int k = 0; k < numfluidphases; ++k)
        {
          // current pressure gradient
          pressuregrad_[k].Clear();

          // phase density
          density_[k] = phasemanager_->Density(k);

          // compute the pressure gradient from the phi gradients
          for (int idof = 0; idof < numfluidphases; ++idof)
            pressuregrad_[k].Update(phasemanager_->PressureDeriv(k, idof), fluidgradphi[idof], 1.0);

          // compute the absolute value of the pressure gradient from the phi gradients
          abspressuregrad_[k] = 0.0;
          for (int i = 0; i < NSD; i++)
            abspressuregrad_[k] += pressuregrad_[k](i) * pressuregrad_[k](i);
          abspressuregrad_[k] = sqrt(abspressuregrad_[k]);

          // diffusion tensor
          difftensorsfluid_[k].Clear();
          phasemanager_->PermeabilityTensor(k, difftensorsfluid_[k]);
          difftensorsfluid_[k].Scale(
              phasemanager_->RelPermeability(k) /
              phasemanager_->DynViscosity(k, abspressuregrad_[k], ndsscatra_porofluid_));

          // Insert Darcy's law: porosity*S^\pi*(v^\pi - v_s) = - k/\mu * grad p
          multiconvelint_[k].Multiply(-1.0, difftensorsfluid_[k], pressuregrad_[k]);

          temperatureconvelint_.Update(heatcapacity_[k] * density_[k], multiconvelint_[k], 1.0);

          // phase pressure
          pressure_[k] = phasemanager_->Pressure(k);
          // phase saturation
          saturation_[k] = phasemanager_->Saturation(k);
          // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
          multiconv_[k].MultiplyTN(derxy, multiconvelint_[k]);

          temperatureconv_.Update(heatcapacity_[k] * density_[k], multiconv_[k], 1.0);
        }

        for (int k = numfluidphases; k < numfluidphases + numvolfrac; ++k)
        {
          // current pressure gradient
          pressuregrad_[k].Update(1.0, fluidgradphi[k + numvolfrac], 0.0);

          // vol frac density
          density_[k] = phasemanager_->VolFracDensity(k - numfluidphases);

          // diffusion tensor
          difftensorsfluid_[k].Clear();
          phasemanager_->PermeabilityTensorVolFracPressure(
              k - numfluidphases, difftensorsfluid_[k]);
          difftensorsfluid_[k].Scale(
              1.0 / phasemanager_->DynViscosityVolFracPressure(k - numfluidphases, -1.0,
                        ndsscatra_porofluid_));  // -1.0 --> don't need abspressgrad

          // Insert Darcy's law: porosity*(v^\pi - v_s) = - k/\mu * grad p
          multiconvelint_[k].Multiply(-1.0, difftensorsfluid_[k], pressuregrad_[k]);

          temperatureconvelint_.Update(heatcapacity_[k] * density_[k], multiconvelint_[k], 1.0);

          // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
          multiconv_[k].MultiplyTN(derxy, multiconvelint_[k]);

          temperatureconv_.Update(heatcapacity_[k] * density_[k], multiconv_[k], 1.0);
        }

        // solid pressure
        solidpressure_ = phasemanager_->SolidPressure();

        for (int k = 0; k < my::numscal_; ++k)
        {
          // overwrite convective term
          // - rho * k/\mu*grad p * grad phi
          switch (scalartophasemap_[k].species_type)
          {
            case MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid:
            case MAT::ScatraMatMultiPoro::SpeciesType::species_in_volfrac:
              my::conv_phi_[k] = multiconvelint_[scalartophasemap_[k].phaseID].Dot(my::gradphi_[k]);
              break;

            case MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid:
              my::conv_phi_[k] = 0;
              break;

            case MAT::ScatraMatMultiPoro::SpeciesType::species_temperature:
              my::conv_phi_[k] = temperatureconvelint_.Dot(my::gradphi_[k]);
              break;

            default:
              dserror(
                  "unknown species type %i for species %i!", scalartophasemap_[k].species_type, k);
              break;
          }
          // set flag if we actually have to evaluate the species
          SetEvaluateScalarFlag(k);
        }

        return;
      };

      // adapt convective term in case of L2-projection
      void AdaptConvectiveTermForL2(
          const LINALG::Matrix<NEN, 1>& funct,  //! array for shape functions
          const LINALG::Matrix<NSD, NEN>&
              derxy,  //! global derivatives of shape functions w.r.t x,y,z
          const std::vector<LINALG::Matrix<NSD, NEN>>&
              econvelnp  //! nodal convective velocity values at t_(n+1) or t_(n+alpha_F)
      )
      {
        const int numfluidphases = econvelnp.size();

        // resize all phase related vectors
        multiconvelint_.resize(numfluidphases);
        multiconv_.resize(numfluidphases);

        for (int k = 0; k < numfluidphases; ++k)
        {
          // convective velocity
          multiconvelint_[k].Multiply(1.0, econvelnp[k], funct);
          // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
          multiconv_[k].MultiplyTN(derxy, multiconvelint_[k]);
        }

        // check if the IDs make sense
        for (unsigned int k = 0; k < scalartophasemap_.size(); ++k)
        {
          if (scalartophasemap_[k].phaseID < 0 or scalartophasemap_[k].phaseID >= numfluidphases)
            dserror("Invalid phase ID %i", scalartophasemap_[k].phaseID);
        }

        // overwrite convective term
        for (int k = 0; k < my::numscal_; ++k)
        {
          switch (scalartophasemap_[k].species_type)
          {
            case MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid:
            case MAT::ScatraMatMultiPoro::SpeciesType::species_in_volfrac:
            case MAT::ScatraMatMultiPoro::SpeciesType::species_temperature:
              my::conv_phi_[k] = multiconvelint_[scalartophasemap_[k].phaseID].Dot(my::gradphi_[k]);
              break;

            case MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid:
              my::conv_phi_[k] = 0;
              break;

            default:
              dserror(
                  "unknown species type %i for species %i!", scalartophasemap_[k].species_type, k);
              break;
          }
        }
        return;
      };

      // Set the fluid-material in the scatra-Varmanager
      void SetFluidPoromultiphaseMaterial(DRT::Element* ele)
      {
        // check if we actually have three materials
        if (ele->NumMaterial() < 3) dserror("no third material available");

        // here we rely that the PoroMultiPhase material has been added as third material
        multiphasemat_ = Teuchos::rcp_dynamic_cast<MAT::FluidPoroMultiPhase>(
            ele->Material(ndsscatra_porofluid_));
        if (multiphasemat_ == Teuchos::null) dserror("cast to MAT::FluidPoroMultiPhase failed!");

        materialset_ = true;
      }

      /*========================================================================*/
      //! @name return methods for internal variables
      /*========================================================================*/

      //! return convective velocity
      virtual const LINALG::Matrix<NSD, 1>& ConVel(const int k) const
      {
        switch (scalartophasemap_[k].species_type)
        {
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid:
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_volfrac:
            return multiconvelint_[scalartophasemap_[k].phaseID];

          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid:
            return zeroconvelint_;

          case MAT::ScatraMatMultiPoro::SpeciesType::species_temperature:
            return temperatureconvelint_;

          case MAT::ScatraMatMultiPoro::SpeciesType::species_undefined:
          default:
            dserror(
                "unknown species type %i for species %i!", scalartophasemap_[k].species_type, k);
        }
        return zeroconvelint_;
      };

      //! return convective part in convective form
      virtual const LINALG::Matrix<NEN, 1>& Conv(const int k) const
      {
        switch (scalartophasemap_[k].species_type)
        {
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid:
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_volfrac:
            return multiconv_[scalartophasemap_[k].phaseID];

          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid:
            return zeroconv_;

          case MAT::ScatraMatMultiPoro::SpeciesType::species_temperature:
            return temperatureconv_;

          default:
            dserror(
                "unknown species type %i for species %i!", scalartophasemap_[k].species_type, k);
            return zeroconv_;
        }
      };

      //! return pressure associated with scalar k
      double Pressure(const int k) const
      {
        switch (scalartophasemap_[k].species_type)
        {
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid:
            return pressure_[scalartophasemap_[k].phaseID];

          default:
            dserror("ScalartophaseID = %i for species %i", scalartophasemap_[k].phaseID, k);
            return 0;
        }
      };

      //! return pressures
      const std::vector<double>& Pressure() const { return pressure_; };

      //! return saturation associated with scalar k
      double Saturation(const int k) const
      {
        switch (scalartophasemap_[k].species_type)
        {
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid:
            return saturation_[scalartophasemap_[k].phaseID];

          default:
            dserror("ScalartophaseID = %i for species %i", scalartophasemap_[k].phaseID, k);
            return 0;
        }
      };

      //! return saturation associated with scalar k
      const std::vector<double>& Saturation() const { return saturation_; };

      //! return density associated with scalar k
      double Density(const int k) const
      {
        switch (scalartophasemap_[k].species_type)
        {
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid:
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_volfrac:
            return density_[scalartophasemap_[k].phaseID];

          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid:
            return phasemanager_->SolidDensity();

          case MAT::ScatraMatMultiPoro::SpeciesType::species_temperature:
            // set to 1.0 because densities are included in the effective heat capacity
            return 1.0;

          default:
            dserror(
                "unknown species type %i for species %i!", scalartophasemap_[k].species_type, k);
            return 0;
        }
      };

      //! return density vector
      const std::vector<double>& Density() const { return density_; };

      //! return volfrac vector
      const std::vector<double>& VolFrac() const { return volfrac_; };

      //! return volume fraction associated with scalar k
      double VolFrac(const int k) const
      {
        switch (scalartophasemap_[k].species_type)
        {
          // case MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid:
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_volfrac:
            return volfrac_[GetPhaseID(k) - phasemanager_->NumFluidPhases()];

          default:
            dserror("ScalartophaseID = %i for species %i", scalartophasemap_[k].phaseID, k);
            return 0;
        }
      };

      //! return volfrac pressure vector
      const std::vector<double>& VolFracPressure() const { return volfracpressure_; };

      //! return solid pressure
      double SolidPressure() const { return solidpressure_; };

      //! set scalar ID to phase ID mapping and species type
      void SetPhaseIDAndSpeciesType(const int scalarID, const int phaseID,
          const MAT::ScatraMatMultiPoro::SpeciesType& spectype)
      {
        scalartophasemap_[scalarID].phaseID = phaseID;
        scalartophasemap_[scalarID].species_type = spectype;
      };

      //! get phase ID from scalar ID
      int GetPhaseID(const int scalarID) const
      {
        if (scalartophasemap_[scalarID].phaseID < 0)
          dserror(
              "ScalartophaseID = %i for species %i", scalartophasemap_[scalarID].phaseID, scalarID);

        return scalartophasemap_[scalarID].phaseID;
      };

      //! get species type of scalar 'k'
      MAT::ScatraMatMultiPoro::SpeciesType GetSpeciesType(const int k)
      {
        return scalartophasemap_[k].species_type;
      };

      //! set delta for evaluation of effective diffusivity
      void SetDelta(const double delta, const int k) { delta_[k] = delta; };

      //! set heat capacity
      void SetHeatCapacity(std::vector<double> cp) { heatcapacity_ = cp; };

      //! set effective heat capacity
      void SetEffectiveHeatCapacity(double cp_eff) { heatcapacityeff_ = cp_eff; };

      //! set thermal diffusivity
      void SetThermalDiffusivity(std::vector<double> kappa) { thermaldiffusivity_ = kappa; };

      //! set minimum value of corresponding phase under which we assume that mass fraction is
      //! also zero
      void SetMinValOfPhase(const double minvalofphase, const int k)
      {
        min_val_of_phase_[k] = minvalofphase;
      };

      //! set action
      void SetAction(const SCATRA::Action action) { myaction_ = action; };

      //! get delta
      double GetDelta(const int k) { return delta_[k]; };

      //! get heat capacity
      //! heat capacity order [ <fluid>  <volfrac>  <solid> ]
      double GetHeatCapacity(const int j) { return heatcapacity_[j]; }

      //! get effective heat capacity
      double GetEffectiveHeatCapacity() { return heatcapacityeff_; }

      //! get thermal diffusivity
      //! thermal diffusivity order [ <fluid>  <volfrac>  <solid> ]
      double GetThermalDiffusivity(const int j) { return thermaldiffusivity_[j]; }

      //! get evaluate scalar flag
      bool EvaluateScalar(const int k) { return evaluate_scalar_[k]; };

      //! get minimum value of corresponding phase under which we assume that mass fraction is
      //! also zero
      double GetMinValOfPhase(const int k) { return min_val_of_phase_[k]; };

      //! get pre-factor needed for OD-mesh-linearization of mass matrix
      double GetPreFactorForMassMatrixODMesh(const int k, const double fac)
      {
        double prefac = fac;

        switch (scalartophasemap_[k].species_type)
        {
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid:
          {
            // linearization of porosity
            //------------------------------------------------dporosity/dd = dporosity/dJ * dJ/dd =
            // dporosity/dJ * J * N_x
            // J denotes the determinant of the deformation gradient, i.e. det F = det ( d x / d X )
            // = det (dx/ds) * ( det(dX/ds) )^-1 in our case: prefac is scaled with density, i.e.
            // porosity --> scale it with 1.0/porosity here

            // we will pass fac/poro*J*dporosity/dJ into CalcMatMassODMesh, where it will internally
            // be scaled correctly

            if (phasemanager_->PorosityDependsOnStruct())
              prefac += fac / (phasemanager_->Porosity()) * phasemanager_->JacobianDefGrad() *
                        phasemanager_->PorosityDerivWrtJacobianDefGrad();
            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_volfrac:
          {
            // do nothing: correct prefac = fac
            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid:
          {
            if (phasemanager_->PorosityDependsOnStruct())
              prefac -= fac / (1 - phasemanager_->Porosity() - phasemanager_->SumAddVolFrac()) *
                        phasemanager_->JacobianDefGrad() *
                        phasemanager_->PorosityDerivWrtJacobianDefGrad();

            // else
            // do nothing: correct prefac = fac

            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_temperature:
          {
            if (phasemanager_->PorosityDependsOnStruct())
            {
              // evaluate the effective heat capacity factor
              const int numdofs = phasemanager_->NumFluidPhases() + phasemanager_->NumVolFrac();
              const int numfluidphases = phasemanager_->NumFluidPhases();

              double cp_eff = GetEffectiveHeatCapacity();

              prefac += fac / cp_eff * phasemanager_->SolidDensity() * GetHeatCapacity(numdofs) *
                        (-1.0) * phasemanager_->JacobianDefGrad() *
                        phasemanager_->PorosityDerivWrtJacobianDefGrad();

              for (int phase = 0; phase < numfluidphases; ++phase)
              {
                prefac += fac / cp_eff * phasemanager_->Density(phase) * GetHeatCapacity(phase) *
                          phasemanager_->JacobianDefGrad() * phasemanager_->Saturation(phase) *
                          phasemanager_->PorosityDerivWrtJacobianDefGrad();
              }
            }
            break;
          }
          default:
            dserror(
                "unknown species type %i for species %i!", scalartophasemap_[k].species_type, k);
            break;
        }  // switch(species type)

        return prefac;
      };

      //! get pre-factor needed for OD-fluid-linearization of mass matrix
      void GetPreFacLinMassODFluid(const int k, std::vector<double>* prefaclinmassodfluid)
      {
        // reset to zero
        std::fill(prefaclinmassodfluid->begin(), prefaclinmassodfluid->end(), 0.0);

        const int numfluidphases = phasemanager_->NumFluidPhases();

        switch (scalartophasemap_[k].species_type)
        {
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid:
          {
            const int curphase = GetPhaseID(k);
            const int totalnummultiphasedofpernode = MultiphaseMat()->NumMat();

            // d poro / d S_j = SaturationDeriv; scaling with 1.0/S since term is re-scaled with
            // densnp = poro*S*rho
            for (int idof = 0; idof < numfluidphases; ++idof)
              (*prefaclinmassodfluid)[idof] += phasemanager_->SaturationDeriv(curphase, idof) /
                                               phasemanager_->Saturation(curphase);

            // linearization of porosity only if porosity is pressure-dependent
            // d poro / d psi_j = porosityDeriv; scaling with 1.0/poro since term is re-scaled with
            // densnp = poro*S*rho
            if (phasemanager_->PorosityDependsOnFluid())
            {
              for (int idof = 0; idof < totalnummultiphasedofpernode; ++idof)
                (*prefaclinmassodfluid)[idof] +=
                    phasemanager_->PorosityDeriv(idof) / phasemanager_->Porosity();
            }
            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_volfrac:
          {
            const int curphase = GetPhaseID(k);
            // d volfrac_j / d volfrac_j = 1.0; scaling with 1.0/volfrac since term is re-scaled
            // with densnp = rho*volfrac
            (*prefaclinmassodfluid)[curphase] += 1.0 / VolFrac(k);

            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid:
          {
            // porosity of solid epsilon_s does not depend on fluid primary variables

            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_temperature:
          {
            // evaluate the effective heat capacity factor
            const int numdofs = phasemanager_->NumFluidPhases() + phasemanager_->NumVolFrac();
            const int totalnummultiphasedofpernode = MultiphaseMat()->NumMat();

            double cp_eff = GetEffectiveHeatCapacity();

            for (int phase = 0; phase < numfluidphases; ++phase)
            {
              for (int idof = 0; idof < numfluidphases; ++idof)
                (*prefaclinmassodfluid)[idof] += phasemanager_->SaturationDeriv(phase, idof) /
                                                 cp_eff * phasemanager_->Density(phase) *
                                                 GetHeatCapacity(phase) * phasemanager_->Porosity();

              if (phasemanager_->PorosityDependsOnFluid())
              {
                for (int idof = 0; idof < totalnummultiphasedofpernode; ++idof)
                  (*prefaclinmassodfluid)[idof] +=
                      phasemanager_->PorosityDeriv(idof) / cp_eff * phasemanager_->Density(phase) *
                      GetHeatCapacity(phase) * phasemanager_->Saturation(phase);
              }
            }

            for (int phase = numfluidphases; phase < numdofs; ++phase)
            {
              (*prefaclinmassodfluid)[phase] +=
                  GetHeatCapacity(phase) * phasemanager_->VolFracDensity(phase - numfluidphases) /
                  cp_eff;
            }
            break;
          }
          default:
          {
            dserror(
                "unknown species type %i for species %i!", scalartophasemap_[k].species_type, k);
            break;
          }
        }  // switch(species type)

        return;
      };

      //! get pre-factor needed for OD-fluid-linearization of convective term
      void GetPreFacLinConvODFluid(const int k, const unsigned ui,
          std::vector<double>* prefaclinconvodfluid, const LINALG::Matrix<NSD, 1>& gradphi,
          const LINALG::Matrix<1, NSD>& gradphiTdifftensor, const LINALG::Matrix<NEN, 1>& funct,
          const LINALG::Matrix<NSD, NEN>& derxy, const int phase)
      {
        // reset to zero
        std::fill(prefaclinconvodfluid->begin(), prefaclinconvodfluid->end(), 0.0);

        // get correct factor
        double laplawf(0.0);
        for (unsigned j = 0; j < NSD; j++) laplawf += derxy(j, ui) * gradphiTdifftensor(0, j);

        const int numfluidphases = phasemanager_->NumFluidPhases();
        const int numvolfrac = phasemanager_->NumVolFrac();

        // fluid phase
        if (phase >= 0 && phase < numfluidphases)
        {
          // derivative after fluid pressures
          for (int idof = 0; idof < numfluidphases; ++idof)
            (*prefaclinconvodfluid)[idof] += laplawf * phasemanager_->PressureDeriv(phase, idof);

          // derivative after relative permeabilty
          if (not phasemanager_->HasConstantRelPermeability(phase))
          {
            const LINALG::Matrix<NSD, 1> gradpres = PressureGradient(phase);
            const double abspressgrad = AbsPressureGradient(phase);

            static LINALG::Matrix<NSD, NSD> difftensor(true);
            phasemanager_->PermeabilityTensor(phase, difftensor);
            difftensor.Scale(phasemanager_->RelPermeabilityDeriv(phase) /
                             phasemanager_->DynViscosity(phase, abspressgrad, 2));

            static LINALG::Matrix<1, NSD> gradphiTdifftensor(true);
            gradphiTdifftensor.MultiplyTN(gradphi, difftensor);

            double laplawf(0.0);
            for (unsigned j = 0; j < NSD; j++) laplawf += gradpres(j) * gradphiTdifftensor(0, j);

            for (int idof = 0; idof < numfluidphases; ++idof)
              (*prefaclinconvodfluid)[idof] +=
                  laplawf * funct(ui) * phasemanager_->SaturationDeriv(phase, idof);
          }
        }
        // volfrac
        else if (phase < numfluidphases + numvolfrac)
        {
          // derivative after volfrac-pressure
          (*prefaclinconvodfluid)[phase + numvolfrac] += laplawf;
        }
        else
          dserror("GetPreFacLinConvODFluid has been called with phase %i", phase);

        return;
      };

      //! get pre-factor needed for OD-mesh-linearization of diffusive term
      double GetPreFactorForDiffMatrixODMesh(
          const int k, const double rhsfac, const double diffcoeff)
      {
        double prefac = 0.0;

        switch (scalartophasemap_[k].species_type)
        {
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid:
          {
            if (phasemanager_->PorosityDependsOnStruct())
            {
              const double delta = GetDelta(k);

              // linearization of porosity
              //------------------------------------------------dporosity/dd = dporosity/dJ * dJ/dd
              //= dporosity/dJ * J * N_x
              // J denotes the determinant of the deformation gradient, i.e. det F = det ( d x / d X
              // ) = det (dx/ds) * ( det(dX/ds) )^-1 in our case: diffusivity is scaled with
              // porosity^(delta + 1) --> scale it with 1.0/porosity^(delta + 1) here
              //              and build derivative d diff/d porosity = (delta + 1) * porosity^delta
              prefac = rhsfac * diffcoeff / std::pow(phasemanager_->Porosity(), delta + 1.0) *
                       (delta + 1.0) * std::pow(phasemanager_->Porosity(), delta) *
                       phasemanager_->JacobianDefGrad() *
                       phasemanager_->PorosityDerivWrtJacobianDefGrad();
            }
            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_volfrac:
          {
            // do nothing: linearization performed in base class
            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid:
          {
            if (phasemanager_->PorosityDependsOnStruct())
            {
              const double delta = GetDelta(k);

              // linearization of porosity
              //------------------------------------------------dporosity/dd =
              // dporosity/dJ * dJ/dd = dporosity/dJ * J * N_x
              // J denotes the determinant of the deformation gradient, i.e. det F =
              // det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1 in our case:
              // diffusivity is scaled with (1-porosity-sumaddvolfrac)^(delta + 1) -->
              // scale it with 1.0/(1-porosity-sumaddvolfrac)^(delta + 1) here
              //              and build derivative d diff/d porosity = (delta + 1) *
              //              (1-porosity-sumaddvolfrac)^delta
              prefac = -1.0 * rhsfac * diffcoeff /
                       std::pow(1.0 - phasemanager_->Porosity() - phasemanager_->SumAddVolFrac(),
                           delta + 1.0) *
                       (delta + 1.0) *
                       std::pow(1.0 - phasemanager_->Porosity() - phasemanager_->SumAddVolFrac(),
                           delta) *
                       phasemanager_->JacobianDefGrad() *
                       phasemanager_->PorosityDerivWrtJacobianDefGrad();
            }
            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_temperature:
          {
            if (phasemanager_->PorosityDependsOnStruct())
            {
              const int numdofs = phasemanager_->NumFluidPhases() + phasemanager_->NumVolFrac();
              const int numfluidphases = phasemanager_->NumFluidPhases();

              prefac = rhsfac * GetThermalDiffusivity(numdofs) * (-1.0) *
                       phasemanager_->JacobianDefGrad() *
                       phasemanager_->PorosityDerivWrtJacobianDefGrad();

              for (int phase = 0; phase < numfluidphases; ++phase)
                prefac += rhsfac * GetThermalDiffusivity(phase) * phasemanager_->JacobianDefGrad() *
                          phasemanager_->Saturation(phase) *
                          phasemanager_->PorosityDerivWrtJacobianDefGrad();
            }
            break;
          }
          default:
            dserror(
                "unknown species type %i for species %i!", scalartophasemap_[k].species_type, k);
            break;
        }  // switch(species type)

        return prefac;
      };

      //! get pre-factor needed for OD-fluid-linearization of diffusive term
      void GetPreFacDiffODFluid(const int k, const double rhsfac, const double diffcoeff,
          std::vector<double>* prefacdiffodfluid)
      {
        // reset to zero
        std::fill(prefacdiffodfluid->begin(), prefacdiffodfluid->end(), 0.0);

        switch (scalartophasemap_[k].species_type)
        {
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid:
          {
            const int curphase = GetPhaseID(k);
            const int numfluidphases = phasemanager_->NumFluidPhases();
            const int totalnummultiphasedofpernode = MultiphaseMat()->NumMat();

            const double delta = GetDelta(k);

            // linearization of saturation*porosity*d_eff = (saturation * porosity)^(delta+1) w.r.t
            // saturation
            //-------------------------------------------------------------------------------------------------------------
            // in our case: diffusivity is scaled with saturation^(delta + 1) --> scale it
            // with 1.0/saturation^(delta + 1) here
            //              and build derivative d diff/d saturation = (delta + 1) *
            //              saturation^delta
            const double vrhs_sat =
                rhsfac * diffcoeff / std::pow(phasemanager_->Saturation(curphase), delta + 1.0) *
                (delta + 1.0) * std::pow(phasemanager_->Saturation(curphase), delta);

            // linearization of saturation*porosity*d_eff = (saturation * porosity)^(delta+1) w.r.t
            // porosity
            //-------------------------------------------------------------------------------------------------------------
            // in our case: diffusivity is scaled with porosity^(delta + 1) --> scale it
            // with 1.0/porosity^(delta + 1) here
            //              and build derivative d diff/d porosity = (delta + 1) * porosity^delta
            const double vrhs_poro = rhsfac * diffcoeff /
                                     std::pow(phasemanager_->Porosity(), delta + 1.0) *
                                     (delta + 1.0) * std::pow(phasemanager_->Porosity(), delta);

            for (int idof = 0; idof < numfluidphases; ++idof)
              (*prefacdiffodfluid)[idof] +=
                  vrhs_sat * phasemanager_->SaturationDeriv(curphase, idof);

            // linearization of porosity only if porosity is pressure-dependent
            if (phasemanager_->PorosityDependsOnFluid())
            {
              for (int idof = 0; idof < totalnummultiphasedofpernode; ++idof)
                (*prefacdiffodfluid)[idof] += vrhs_poro * phasemanager_->PorosityDeriv(idof);
            }
            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_volfrac:
          {
            const int curphase = GetPhaseID(k);
            const double delta = GetDelta(k);
            // diffusivity is scaled with volfrac^(delta + 1) --> scale it
            // with 1.0/volfrac^(delta + 1) here
            // and build derivative d diff/d volfrac = (delta + 1) * volfrac^delta
            (*prefacdiffodfluid)[curphase] += rhsfac * diffcoeff /
                                              std::pow(VolFrac(k), delta + 1.0) * (delta + 1.0) *
                                              std::pow(VolFrac(k), delta);
            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid:
          {
            // porosity of solid epsilon_s does not depend on fluid primary variables
            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_temperature:
          {
            const int numfluidphases = phasemanager_->NumFluidPhases();
            const int numdofs = numfluidphases + phasemanager_->NumVolFrac();
            const int totalnummultiphasedofpernode = MultiphaseMat()->NumMat();

            for (int phase = 0; phase < numfluidphases; ++phase)
            {
              for (int idof = 0; idof < numfluidphases; ++idof)
                (*prefacdiffodfluid)[idof] += rhsfac * phasemanager_->SaturationDeriv(phase, idof) *
                                              GetThermalDiffusivity(phase) *
                                              phasemanager_->Porosity();

              if (phasemanager_->PorosityDependsOnFluid())
                for (int idof = 0; idof < totalnummultiphasedofpernode; ++idof)
                  (*prefacdiffodfluid)[idof] += rhsfac * phasemanager_->PorosityDeriv(idof) *
                                                GetThermalDiffusivity(phase) *
                                                phasemanager_->Saturation(phase);
            }

            for (int phase = numfluidphases; phase < numdofs; ++phase)
              (*prefacdiffodfluid)[phase] += rhsfac * GetThermalDiffusivity(phase);

            break;
          }
          default:
          {
            dserror(
                "unknown species type %i for species %i!", scalartophasemap_[k].species_type, k);
            break;
          }
        }  // switch(species type)

        return;
      };

      //! get difftensor of fluid phases (permeability*relpermeability/mu)
      void GetDiffTensorFluid(const int k, LINALG::Matrix<NSD, NSD>& difftensor, const int phase)
      {
        switch (scalartophasemap_[k].species_type)
        {
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid:
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_volfrac:
          case MAT::ScatraMatMultiPoro::SpeciesType::species_temperature:
          {
            difftensor = difftensorsfluid_[phase];
            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid:
          default:
            dserror(
                "unknown species type %i for species %i!", scalartophasemap_[k].species_type, k);
            break;
        }  // switch(species type)

        return;
      }

      //! compute gradient of pressure in reference configuration
      void GetRefGradPres(const int k, const LINALG::Matrix<NSD, NSD>& xjm,
          LINALG::Matrix<NSD, 1>& refgradpres, const int phase)
      {
        refgradpres.Clear();

        const int numfluidphases = phasemanager_->NumFluidPhases();
        const int numvolfrac = phasemanager_->NumVolFrac();

        // fluid phase
        if (phase >= 0 && phase < numfluidphases)
        {
          const std::vector<LINALG::Matrix<NSD, 1>>& fluidgradphi =
              *(variablemanager_->GradPhinp());

          // gradient of phi w.r.t. reference coordinates
          std::vector<LINALG::Matrix<NSD, 1>> reffluidgradphi(
              numfluidphases, LINALG::Matrix<NSD, 1>(true));
          for (int idof = 0; idof < numfluidphases; ++idof)
            reffluidgradphi[idof].Multiply(xjm, fluidgradphi[idof]);

          // compute the pressure gradient from the phi gradients
          for (int idof = 0; idof < numfluidphases; ++idof)
            refgradpres.Update(
                phasemanager_->PressureDeriv(phase, idof), reffluidgradphi[idof], 1.0);
        }
        // volfrac
        else if (phase < numfluidphases + numvolfrac)
        {
          refgradpres.Multiply(xjm, pressuregrad_[phase]);
        }
        else
          dserror("GetRefGradPres has been called with phase %i", phase);

        return;
      }

      // set flag if we actually evaluate the scalar or set it to zero
      void SetEvaluateScalarFlag(const int k)
      {
        const int numfluidphases = phasemanager_->NumFluidPhases();

        switch (scalartophasemap_[k].species_type)
        {
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid:
          {
            // if smaller than threshold (at GP) we do not evaluate
            evaluate_scalar_[k] =
                (fabs(saturation_[scalartophasemap_[k].phaseID]) > min_val_of_phase_[k]);
            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_volfrac:
          {
            // if smaller than minimum nodal volume fraction in element we do not evaluate
            evaluate_scalar_[k] = variablemanager_->ElementHasValidVolFracSpecies(
                scalartophasemap_[k].phaseID - numfluidphases);
            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid:
          {
            evaluate_scalar_[k] = true;
            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_temperature:
          {
            evaluate_scalar_[k] = true;
            break;
          }
          default:
          {
            dserror(
                "unknown species type %i for species %i!", scalartophasemap_[k].species_type, k);
            break;
          }
        }  // switch(species type)

        return;
      }

      //! get pre-factor needed for OD-mesh-linearization of hist and source term
      double GetPreFactorForHistAndSourceODMesh(const int k, const double fac, const double densnp,
          const double hist, const double rhsint)
      {
        // 1) linearization of mesh motion:
        //    call base class with fac*rhsint
        double prefac = fac * rhsint;

        switch (scalartophasemap_[k].species_type)
        {
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid:
          {
            // 2) linearization of history:
            //    prefactor is densnp = porosity * rho * S
            //    --> porosity has to be linearized
            // with
            // linearization of porosity
            //------------------------------------------------dporosity/dd = dporosity/dJ * dJ/dd =
            // dporosity/dJ * J * N_x
            // J denotes the determinant of the deformation gradient, i.e. det F = det ( d x / d X )
            // = det (dx/ds) * ( det(dX/ds) )^-1 in our case: prefac is scaled with density, i.e.
            // porosity --> scale it with 1.0/porosity here

            if (phasemanager_->PorosityDependsOnStruct())
              prefac += 1.0 * fac * hist * densnp / (phasemanager_->Porosity()) *
                        phasemanager_->JacobianDefGrad() *
                        phasemanager_->PorosityDerivWrtJacobianDefGrad();
            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_volfrac:
          {
            // do nothing: correct prefac = fac*rhsint
            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid:
          {
            if (phasemanager_->PorosityDependsOnStruct())
              prefac -= fac * hist * densnp /
                        (1 - phasemanager_->Porosity() - phasemanager_->SumAddVolFrac()) *
                        phasemanager_->JacobianDefGrad() *
                        phasemanager_->PorosityDerivWrtJacobianDefGrad();
            // else
            // do nothing: prefac = fac
            break;
          }
          case MAT::ScatraMatMultiPoro::SpeciesType::species_temperature:
          {
            // evaluate the effective heat capacity factor
            const int numdofs = phasemanager_->NumFluidPhases() + phasemanager_->NumVolFrac();
            const int numfluidphases = phasemanager_->NumFluidPhases();

            double cp_eff = GetEffectiveHeatCapacity();

            prefac += fac * hist * densnp / cp_eff * phasemanager_->SolidDensity() *
                      GetHeatCapacity(numdofs) * (-1.0) * phasemanager_->JacobianDefGrad() *
                      phasemanager_->PorosityDerivWrtJacobianDefGrad();

            for (int phase = 0; phase < numfluidphases; ++phase)
              prefac += fac * hist * densnp / cp_eff * phasemanager_->Density(phase) *
                        GetHeatCapacity(phase) * phasemanager_->Saturation(phase) *
                        phasemanager_->PorosityDerivWrtJacobianDefGrad() *
                        phasemanager_->JacobianDefGrad();

            break;
          }

          default:
            dserror(
                "unknown species type %i for species %i!", scalartophasemap_[k].species_type, k);
            break;
        }  // switch(species type)

        return prefac;
      };

      //! get action
      SCATRA::Action GetAction()
      {
        // return
        return myaction_;
      };

      //! return pressure gradient of current phase
      const LINALG::Matrix<NSD, 1>& PressureGradient(const int curphase) const
      {
        return pressuregrad_[curphase];
      };

      //! return absolute value of pressure gradient of current phase
      double AbsPressureGradient(const int curphase) const { return abspressuregrad_[curphase]; };

      //! return fluid material
      Teuchos::RCP<MAT::FluidPoroMultiPhase> MultiphaseMat()
      {
        if (!materialset_)
          dserror("Fluid-Multiphase Material has not yet been set in Variablemanager");

        return multiphasemat_;
      }

      //! Setup phasemanager and variablemanager of fluid
      void SetupPoroFluidManagers(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          const int numfluidphases, const int totalnummultiphasedofpernode)
      {
        // dummy parameter list
        DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter* para =
            DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter::Instance(discretization.Name());

        // create phase-manager
        phasemanager_ = DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface::CreatePhaseManager(
            *para, NSD, MultiphaseMat()->MaterialType(),
            POROFLUIDMULTIPHASE::Action::get_access_from_scatra, totalnummultiphasedofpernode,
            numfluidphases);

        // access from outside to the phasemanager: scatra-discretization has fluid-dis on dofset 2
        phasemanager_->Setup(ele, ndsscatra_porofluid_);

        // create variablemanager
        variablemanager_ = DRT::ELEMENTS::POROFLUIDMANAGER::VariableManagerInterface<NSD,
            NEN>::CreateVariableManager(*para, POROFLUIDMULTIPHASE::Action::get_access_from_scatra,
            MultiphaseMat(), totalnummultiphasedofpernode, numfluidphases);

        return;
      }

      // extract the element and node values of the poro-fluid --> extract them from its
      // variablemanager
      void ExtractElementAndNodeValuesOfPoroFluid(DRT::Element* ele,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          LINALG::Matrix<NSD, NEN>& xyze)
      {
        // access from outside to the variablemananger: scatra-discretization has fluid-dis on
        // dofset 2
        variablemanager_->ExtractElementAndNodeValues(
            *ele, discretization, la, xyze, ndsscatra_porofluid_);

        return;
      }

      // get the phasemanager of the fluid
      Teuchos::RCP<DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface> FluidPhaseManager()
      {
        return phasemanager_;
      }

      // get the variablemanager of the fluid
      Teuchos::RCP<DRT::ELEMENTS::POROFLUIDMANAGER::VariableManagerInterface<NSD, NEN>>
      FluidVarManager()
      {
        return variablemanager_;
      }

     private:
      //! convective velocity
      std::vector<LINALG::Matrix<NSD, 1>> multiconvelint_;
      //! convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
      std::vector<LINALG::Matrix<NEN, 1>> multiconv_;
      //! zero convective velocity
      LINALG::Matrix<NSD, 1> zeroconvelint_;
      //! zero convective part in convective form
      LINALG::Matrix<NEN, 1> zeroconv_;
      //! temperature convective velocity
      LINALG::Matrix<NSD, 1> temperatureconvelint_;
      //! temperature convective part in convective form
      LINALG::Matrix<NEN, 1> temperatureconv_;

      //! phase pressure
      std::vector<double> pressure_;
      //! phase saturation
      std::vector<double> saturation_;
      //! phase density
      std::vector<double> density_;
      //! solid pressure
      double solidpressure_;
      //! pressure gradient
      std::vector<LINALG::Matrix<NSD, 1>> pressuregrad_;
      //! pressure gradient
      std::vector<LINALG::Matrix<NSD, NSD>> difftensorsfluid_;
      //! norm of pressure-gradient
      std::vector<double> abspressuregrad_;
      //! volume fraction
      std::vector<double> volfrac_;
      //! volume fraction pressure
      std::vector<double> volfracpressure_;

      //! fluid-poro multiphase material
      Teuchos::RCP<MAT::FluidPoroMultiPhase> multiphasemat_;

      //! delta for effective diffusivity
      std::vector<double> delta_;
      //! heat capacity
      //! order [ <fluid>  <volfrac>  <solid> ]
      std::vector<double> heatcapacity_;
      //! thermal diffusivity
      //! order [ <fluid>  <volfrac>  <solid> ]
      std::vector<double> thermaldiffusivity_;
      //! effective heat capacity
      double heatcapacityeff_;

      //! minimum saturation under which the corresponding mass fraction is assumed to be zero
      std::vector<double> min_val_of_phase_;
      //! flag to check if we have to evaluate the species equation in this element
      std::vector<bool> evaluate_scalar_;
      //! scalar to phase map
      std::vector<MAT::ScatraMatMultiPoro::ScalarToPhaseMap> scalartophasemap_;
      //! check if multiphase material has been set
      bool materialset_;

      SCATRA::Action myaction_;

      //! phase manager of the fluid
      Teuchos::RCP<DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface> phasemanager_;

      //! variable manager of the fluid
      Teuchos::RCP<DRT::ELEMENTS::POROFLUIDMANAGER::VariableManagerInterface<NSD, NEN>>
          variablemanager_;

      //! dofset of fluid field on scatra dis
      // TODO: find a better way to do this
      const int ndsscatra_porofluid_ = 2;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif /* SRC_DRT_SCATRA_ELE_SCATRA_ELE_CALC_MULTIPORO_REAC_H_ */
