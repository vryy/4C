/*----------------------------------------------------------------------*/
/*! \file
 \brief evaluation class containing routines for calculation of scalar transport
        within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_ELE_CALC_MULTIPORO_REAC_HPP
#define FOUR_C_SCATRA_ELE_CALC_MULTIPORO_REAC_HPP


#include "baci_config.hpp"

#include "baci_global_data.hpp"
#include "baci_lib_discret.hpp"
#include "baci_mat_fluidporo_multiphase.hpp"
#include "baci_mat_fluidporo_multiphase_reactions.hpp"
#include "baci_mat_fluidporo_singlephase.hpp"
#include "baci_mat_scatra_multiporo.hpp"
#include "baci_porofluidmultiphase_ele_action.hpp"
#include "baci_porofluidmultiphase_ele_calc.hpp"
#include "baci_porofluidmultiphase_ele_calc_utils.hpp"
#include "baci_porofluidmultiphase_ele_evaluator.hpp"
#include "baci_porofluidmultiphase_ele_parameter.hpp"
#include "baci_porofluidmultiphase_ele_phasemanager.hpp"
#include "baci_porofluidmultiphase_ele_variablemanager.hpp"
#include "baci_scatra_ele_calc_poro_reac.hpp"
#include "baci_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN

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

    template <CORE::FE::CellType distype>
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
      using my::nen_;
      using my::nsd_;

     public:
      /// Singleton access method
      static ScaTraEleCalcMultiPoroReac<distype>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

      /// Setup element evaluation
      int SetupCalc(DRT::Element* ele, DRT::Discretization& discretization) override;

     protected:
      //! extract element based or nodal values
      //  return extracted values of phinp
      void ExtractElementAndNodeValues(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la) override;

      //! extract element based or nodal values --> L2-projection case: called within
      //! ExtractElementAndNodeValues
      //  return extracted values of phinp
      virtual void ExtractNodalFlux(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          const int numfluidphases);

      //! set internal variables
      void SetInternalVariablesForMatAndRHS() override;

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
      void SetAdvancedReactionTerms(const int k,                  //!< index of current scalar
          const Teuchos::RCP<MAT::MatListReactions> matreaclist,  //!< index of current scalar
          const double* gpcoord  //!< current Gauss-point coordinates
          ) override;

      //! compute pore pressure
      double ComputePorePressure() override;

      //! get internal variable manager for multiporo formulation
      Teuchos::RCP<ScaTraEleInternalVariableManagerMultiPoro<nsd_, nen_>> VarManager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleInternalVariableManagerMultiPoro<nsd_, nen_>>(
            my::scatravarmanager_);
      };

      //! calculation of convective element matrix in convective form
      //! the only difference to the base class version is, that there is no scaling with the
      //! density
      void CalcMatConv(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                         //!< index of current scalar
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double densnp,      //!< density at time_(n+1)
          const CORE::LINALG::Matrix<nen_, 1>& sgconv  //!< subgrid-scale convective operator
          ) override;

      //! adaption of convective term for rhs
      //! the only difference to the base class version is, that there is no scaling with the
      //! density
      void RecomputeConvPhiForRhs(const int k,            //!< index of current scalar
          const CORE::LINALG::Matrix<nsd_, 1>& sgvelint,  //!< subgrid-scale velocity at Gauss point
          const double densnp,                            //!< density at time_(n+1)
          const double densn,                             //!< density at time_(n)
          const double vdiv                               //!< velocity divergence
          ) override;

      //! calculation of convective element matrix: add conservative contributions
      //! the only difference to the base class version is, that there is no scaling with the
      //! density
      void CalcMatConvAddCons(
          CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                            //!< index of current scalar
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double vdiv,        //!< velocity divergence
          const double densnp       //!< density at time_(n+1)
          ) override;

      //! calculation of mass element matrix (standard shape functions)
      void CalcMatMass(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int& k,                                        //!< index of current scalar
          const double& fac,                                   //!< domain-integration factor
          const double& densam                                 //!< density at time_(n+am)
          ) override;

      //! calculation of convective element matrix (OD term structure coupling)
      //! difference to base class: linearization of mesh motion + shapederivatives pressure
      //! gradient have to be included
      void CalcConvODMesh(CORE::LINALG::SerialDenseMatrix& emat, const int k,
          const int ndofpernodemesh, const double fac, const double rhsfac, const double densnp,
          const double J, const CORE::LINALG::Matrix<nsd_, 1>& gradphi,
          const CORE::LINALG::Matrix<nsd_, 1>& convelint) override;

      //! calculation of linearized mass (off diagonal/shapederivative term mesh)
      //! difference to base class: linearization of porosity is included
      void CalcLinMassODMesh(
          CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                            //!< index of current scalar
          const int ndofpernodemesh,              //!< number of dofs per node of ale element
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double fac,     //!< domain-integration factor
          const double densam,  //!< density at time_(n+am)
          const double densnp,  //!< density at time_(n+1)
          const double phinp,   //!< scalar at time_(n+1)
          const double hist,    //!< history of time integartion
          const double J,       //!< determinant of Jacobian det(dx/ds)
          const CORE::LINALG::Matrix<1, nsd_ * nen_>&
              dJ_dmesh  //!< derivative of det(dx/ds) w.r.t. mesh displacement
          ) override;

      //! standard Galerkin transient, old part of rhs and source term (off diagonal/shapederivative
      //! term mesh) difference to base class: linearization of porosity and advanced reaction terms
      //! are included
      void CalcHistAndSourceODMesh(
          CORE::LINALG::SerialDenseMatrix& emat,  //!< element current to be filled
          const int k,                            //!< index of current scalar
          const int ndofpernodemesh,              //!< number of dofs per node of ale element
          const double fac,                       //!< domain-integration factor
          const double rhsint,                    //!< rhs at Gauss point
          const double J,                         //!< determinant of Jacobian det(dx/ds)
          const CORE::LINALG::Matrix<1, nsd_ * nen_>&
              dJ_dmesh,        //!< derivative of det(dx/ds) w.r.t. mesh displacement
          const double densnp  //!< density
          ) override;

      //! standard Galerkin diffusive term (off diagonal/shapederivative term mesh)
      //! difference to base class: linearization of porosity and effective diffusivity
      void CalcDiffODMesh(CORE::LINALG::SerialDenseMatrix& emat,  //!< element current to be filled
          const int k,                                            //!< index of current scalar
          const int ndofpernodemesh,  //!< number of dofs per node of ale element
          const double diffcoeff,     //!< diffusion coefficient
          const double fac,           //!< domain-integration factor
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double J,       //!< determinant of Jacobian det(dx/ds)
          const CORE::LINALG::Matrix<nsd_, 1>& gradphi,    //!< scalar gradient at Gauss point
          const CORE::LINALG::Matrix<nsd_, 1>& convelint,  //!< convective velocity
          const CORE::LINALG::Matrix<1, nsd_ * nen_>&
              dJ_dmesh  //!< derivative of det(dx/ds) w.r.t. mesh displacement
          ) override;

      //! reactive terms (standard Galerkin) (off diagonal/shapederivative term mesh)
      //! difference to base class: linearization of porosity
      void CalcReactODMesh(CORE::LINALG::SerialDenseMatrix& emat,  //!< element current to be filled
          const int k,                                             //!< index of current scalar
          const int ndofpernodemesh,  //!< number of dofs per node of ale element
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double rea_phi,  //!< reactive term
          const double J,        //!< determinant of Jacobian det(dx/ds)
          const CORE::LINALG::Matrix<1, nsd_ * nen_>&
              dJ_dmesh  //!< derivative of det(dx/ds) w.r.t. mesh displacement
          ) override;

      //! calculation of convective element matrix in convective form (off diagonal term fluid)
      void CalcMatConvODFluid(
          CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                            //!< index of current scalar
          const int ndofpernodefluid,  //!< number of dofs per node of fluid element // only a dummy
                                       //!< variable
          const double rhsfac,         //!< domain-integration factor times time-integration factor
          const double densnp,         //!< density at time_(n+1)
          const CORE::LINALG::Matrix<nsd_, 1>& gradphi  //!< scalar gradient
          ) override;

      //! calculation of convective element matrix in convective form -- additional conservative
      //! contributions (off diagonal term fluid) not yet implemented --> dserror
      void CalcMatConvAddConsODFluid(
          CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                            //!< index of current scalar
          const int ndofpernodefluid,             //!< number of dofs per node of fluid element
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double densnp,      //!< density at time_(n+1)
          const double phinp        //!< scalar at time_(n+1)
          ) override;

      //! calculation of linearized mass (off diagonal terms fluid)
      //! linearization of porosity*saturation*vtrans
      void CalcLinMassODFluid(
          CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                            //!< index of current scalar
          const int ndofpernodemesh,  //!< number of dofs per node of fluid element // only a dummy
                                      //!< variable
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double fac,     //!< domain-integration factor
          const double densam,  //!< density at time_(n+am)
          const double densnp,  //!< density at time_(n+1)
          const double phinp,   //!< scalar at time_(n+1)
          const double hist     //!< history of time integartion
          ) override;

      // calculate linearization of a mass matrix type matrix (OD-fluid terms)
      virtual void CalcLinMassMatrixTypeODFluid(CORE::LINALG::SerialDenseMatrix& emat, const int k,
          const std::vector<double>* prefaclinmassodfluid, const int totalnummultiphasedofpernode,
          double prefac);

      //! standard Galerkin transient, old part of rhs and source term (off diagonal terms fluid)
      //! linearization of porosity*saturation*vtrans + advanced reaction terms
      void CalcHistAndSourceODFluid(
          CORE::LINALG::SerialDenseMatrix& emat,  //!< element current to be filled
          const int k,                            //!< index of current scalar
          const int ndofpernodemesh,              //!< number of dofs per node of ale element
          const double fac,                       //!< domain-integration factor
          const double rhsint,                    //!< rhs at Gauss point
          const double densnp                     //!< density
          ) override;

      //! standard Galerkin reactive term (off diagonal terms fluid)
      //! linearization of porosity*saturation*vreact
      void CalcReactODFluid(
          CORE::LINALG::SerialDenseMatrix& emat,  //!< element current to be filled
          const int k,                            //!< index of current scalar
          const int ndofpernodemesh,              //!< number of dofs per node of ale element
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double rea_phi  //!< rhs at Gauss point
          ) override;

      //! standard Galerkin diffusive term (off diagonal terms fluid)
      //! linearization of porosity*saturation*d_eff
      void CalcDiffODFluid(CORE::LINALG::SerialDenseMatrix& emat,  //!< element current to be filled
          const int k,                                             //!< index of current scalar
          const int ndofpernodemesh,  //!< number of dofs per node of ale element
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const CORE::LINALG::Matrix<nsd_, 1>& gradphi  //!< scalar gradient at Gauss point
          ) override;

      //! fill coupling vector and add variables to reaction in order to compute reaction values and
      //! derivatives
      void FillCouplingVectorAndAddVariables(const int k,
          const Teuchos::RCP<MAT::MatListReactions> matreaclist,
          const Teuchos::RCP<ScaTraEleReaManagerAdvReac> remanager);

      //! Get right hand side including reaction bodyforce term
      void GetRhsInt(double& rhsint,  //!< rhs containing bodyforce at Gauss point
          const double densnp,        //!< density at t_(n+1)
          const int k                 //!< index of current scalar
          ) override;

      //! calculation of reactive element matrix
      void CalcMatReact(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                          //!< index of current scalar
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double
              timetaufac,  //!< domain-integration factor times time-integration factor times tau
          const double taufac,                          //!< domain-integration factor times tau
          const double densnp,                          //!< density at time_(n+1)
          const CORE::LINALG::Matrix<nen_, 1>& sgconv,  //!< subgrid-scale convective operator
          const CORE::LINALG::Matrix<nen_, 1>& diff     //!< laplace term
          ) override;

      //! "shapederivatives" pressure gradient
      virtual void ApplyShapeDerivsPressureGrad(CORE::LINALG::SerialDenseMatrix& emat, const int k,
          const double vrhs, const CORE::LINALG::Matrix<nsd_, 1>& gradphi,
          const CORE::LINALG::Matrix<nsd_, 1> refgradpres);

     protected:
      //! nodal flux values at t_(n+1)
      std::vector<CORE::LINALG::Matrix<nsd_, nen_>> efluxnp_;

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
            pressure_(0),
            saturation_(0),
            density_(0),
            solidpressure_(0.0),
            delta_(numscal, 0.0),
            relative_mobility_funct_id_(numscal),
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

      // compute and set internal variables -- no L2-projection but evaluation at GP
      void SetInternalVariablesMultiPoro(
          const CORE::LINALG::Matrix<NEN, 1>& funct,  //! array for shape functions
          const CORE::LINALG::Matrix<NSD, NEN>&
              derxy,  //! global derivatives of shape functions w.r.t x,y,z
          const CORE::LINALG::Matrix<NSD, NEN>&
              deriv,  //! global derivatives of shape functions w.r.t r,s,t
          const CORE::LINALG::Matrix<NSD, NSD>& xjm, const CORE::LINALG::Matrix<NSD, NEN>& xyze0,
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>&
              ephinp,  //! scalar at t_(n+1) or t_(n+alpha_F)
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>& ephin,  //! scalar at t_(n)
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>&
              ehist,  //! history vector of transported scalars
          const CORE::LINALG::Matrix<NSD, NEN>&
              eforcevelocity  //! nodal velocity due to external force
      )
      {
        // call base class (scatra) with dummy variable
        const CORE::LINALG::Matrix<NSD, NEN> dummy_econv(true);
        my::SetInternalVariables(funct, derxy, ephinp, ephin, dummy_econv, ehist, dummy_econv);

        // velocity due to the external force
        CORE::LINALG::Matrix<NSD, 1> force_velocity;
        force_velocity.Multiply(eforcevelocity, funct);

        //------------------------get determinant of Jacobian dX / ds
        // transposed jacobian "dX/ds"
        CORE::LINALG::Matrix<NSD, NSD> xjm0;
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
        relative_mobility_funct_id_.resize(my::numscal_);

        const std::vector<CORE::LINALG::Matrix<NSD, 1>>& fluidgradphi =
            *(variablemanager_->GradPhinp());

        volfrac_ = phasemanager_->VolFrac();
        volfracpressure_ = phasemanager_->VolFracPressure();

        //! convective velocity
        std::vector<CORE::LINALG::Matrix<NSD, 1>> phase_fluid_velocity(0.0);
        phase_fluid_velocity.resize(numfluidphases + numvolfrac);
        //! convective part in convective form: u_x*N,x + u_y*N,y
        std::vector<CORE::LINALG::Matrix<NEN, 1>> phase_fluid_velocity_conv(0.0);
        phase_fluid_velocity_conv.resize(numfluidphases + numvolfrac);

        //! temperature convective velocity
        CORE::LINALG::Matrix<NSD, 1> temperatureconvelint(true);
        //! temperature convective part in convective form
        CORE::LINALG::Matrix<NEN, 1> temperatureconv(true);

        for (int i_phase = 0; i_phase < numfluidphases; ++i_phase)
        {
          // current pressure gradient
          pressuregrad_[i_phase].Clear();

          // phase density
          density_[i_phase] = phasemanager_->Density(i_phase);

          // compute the pressure gradient from the phi gradients
          for (int idof = 0; idof < numfluidphases; ++idof)
            pressuregrad_[i_phase].Update(
                phasemanager_->PressureDeriv(i_phase, idof), fluidgradphi[idof], 1.0);

          // compute the absolute value of the pressure gradient from the phi gradients
          abspressuregrad_[i_phase] = 0.0;
          for (int i = 0; i < NSD; i++)
            abspressuregrad_[i_phase] += pressuregrad_[i_phase](i) * pressuregrad_[i_phase](i);
          abspressuregrad_[i_phase] = sqrt(abspressuregrad_[i_phase]);

          // diffusion tensor
          difftensorsfluid_[i_phase].Clear();
          phasemanager_->PermeabilityTensor(i_phase, difftensorsfluid_[i_phase]);
          difftensorsfluid_[i_phase].Scale(phasemanager_->RelPermeability(i_phase) /
                                           phasemanager_->DynViscosity(i_phase,
                                               abspressuregrad_[i_phase], ndsscatra_porofluid_));

          // Insert Darcy's law: porosity*S^\pi*(v^\pi - v_s) = - phase/\mu * grad p
          phase_fluid_velocity[i_phase].Multiply(
              -1.0, difftensorsfluid_[i_phase], pressuregrad_[i_phase]);
          temperatureconvelint.Update(
              heatcapacity_[i_phase] * density_[i_phase], phase_fluid_velocity[i_phase], 1.0);
          // in convective form: u_x*N,x + u_y*N,y
          phase_fluid_velocity_conv[i_phase].MultiplyTN(derxy, phase_fluid_velocity[i_phase]);
          temperatureconv.Update(
              heatcapacity_[i_phase] * density_[i_phase], phase_fluid_velocity_conv[i_phase], 1.0);

          // phase pressure
          pressure_[i_phase] = phasemanager_->Pressure(i_phase);
          // phase saturation
          saturation_[i_phase] = phasemanager_->Saturation(i_phase);
        }

        for (int i_volfrac = numfluidphases; i_volfrac < numfluidphases + numvolfrac; ++i_volfrac)
        {
          // current pressure gradient
          pressuregrad_[i_volfrac].Update(1.0, fluidgradphi[i_volfrac + numvolfrac], 0.0);

          // vol frac density
          density_[i_volfrac] = phasemanager_->VolFracDensity(i_volfrac - numfluidphases);

          // diffusion tensor
          difftensorsfluid_[i_volfrac].Clear();
          phasemanager_->PermeabilityTensorVolFracPressure(
              i_volfrac - numfluidphases, difftensorsfluid_[i_volfrac]);
          difftensorsfluid_[i_volfrac].Scale(
              1.0 / phasemanager_->DynViscosityVolFracPressure(i_volfrac - numfluidphases, -1.0,
                        ndsscatra_porofluid_));  // -1.0 --> don't need abspressgrad

          // Insert Darcy's law: porosity*(v^\pi - v_s) = - k/\mu * grad p
          phase_fluid_velocity[i_volfrac].Multiply(
              -1.0, difftensorsfluid_[i_volfrac], pressuregrad_[i_volfrac]);
          temperatureconvelint.Update(
              heatcapacity_[i_volfrac] * density_[i_volfrac], phase_fluid_velocity[i_volfrac], 1.0);
          // in convective form: u_x*N,x + u_y*N,y
          phase_fluid_velocity_conv[i_volfrac].MultiplyTN(derxy, phase_fluid_velocity[i_volfrac]);
          temperatureconv.Update(heatcapacity_[i_volfrac] * density_[i_volfrac],
              phase_fluid_velocity_conv[i_volfrac], 1.0);
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
              my::convelint_[k] = phase_fluid_velocity[scalartophasemap_[k].phaseID];
              // if the scalar reacts to the external force, add the velocity due to the external
              // force scaled with the relative mobility and the porosity * saturation
              if (my::reacts_to_force_[k])
              {
                const auto prefactor = EvaluateRelativeMobility(k) * phasemanager_->Porosity() *
                                       phasemanager_->Saturation(scalartophasemap_[k].phaseID);
                my::convelint_[k].Update(prefactor, force_velocity, 1.0);
              }
              my::conv_[k].MultiplyTN(derxy, my::convelint_[k]);
              my::conv_phi_[k] = my::convelint_[k].Dot(my::gradphi_[k]);
              break;

            case MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid:
              my::convelint_[k] = CORE::LINALG::Matrix<NSD, 1>(0.0);
              my::conv_[k] = CORE::LINALG::Matrix<NEN, 1>(0.0);
              my::conv_phi_[k] = 0;
              break;

            case MAT::ScatraMatMultiPoro::SpeciesType::species_temperature:
              my::convelint_[k] = temperatureconvelint;
              my::conv_[k] = temperatureconv;
              my::conv_phi_[k] = temperatureconvelint.Dot(my::gradphi_[k]);
              break;

            default:
              dserror(
                  "unknown species type %i for species %i!", scalartophasemap_[k].species_type, k);
              break;
          }
          // set flag if we actually have to evaluate the species
          SetEvaluateScalarFlag(k);
        }
      };

      // adapt convective term in case of L2-projection
      void AdaptConvectiveTermForL2(
          const CORE::LINALG::Matrix<NEN, 1>& funct,  //! array for shape functions
          const CORE::LINALG::Matrix<NSD, NEN>&
              derxy,  //! global derivatives of shape functions w.r.t x,y,z
          const std::vector<CORE::LINALG::Matrix<NSD, NEN>>&
              efluxnp  //! nodal flux values at t_(n+1) or t_(n+alpha_F)
      )
      {
        const int numfluidphases = efluxnp.size();

        std::vector<CORE::LINALG::Matrix<NSD, 1>> flux(0.0);
        flux.resize(numfluidphases);

        // in convective form: q_x*N,x + q_y*N,y
        std::vector<CORE::LINALG::Matrix<NEN, 1>> flux_conv(0.0);
        flux_conv.resize(numfluidphases);

        for (int i_phase = 0; i_phase < numfluidphases; ++i_phase)
        {
          flux[i_phase].Multiply(1.0, efluxnp[i_phase], funct);
          flux_conv[i_phase].MultiplyTN(derxy, flux[i_phase]);
        }

        for (auto& k : scalartophasemap_)
        {
          if (k.phaseID < 0 or k.phaseID >= numfluidphases)
            dserror("Invalid phase ID %i", k.phaseID);
        }

        // set convective term
        for (int k = 0; k < my::numscal_; ++k)
        {
          switch (scalartophasemap_[k].species_type)
          {
            case MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid:
            case MAT::ScatraMatMultiPoro::SpeciesType::species_in_volfrac:
              my::conv_phi_[k] = flux[scalartophasemap_[k].phaseID].Dot(my::gradphi_[k]);
              break;

            case MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid:
            case MAT::ScatraMatMultiPoro::SpeciesType::species_temperature:
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

      //! set relative mobility function ID
      void SetRelativeMobilityFunctionId(
          const int relative_mobility_funct_id, const int current_scalar)
      {
        relative_mobility_funct_id_[current_scalar] = relative_mobility_funct_id;
      };

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

      //! evaluate relative mobility
      auto EvaluateRelativeMobility(const int current_scalar)
      {
        std::vector<std::pair<std::string, double>> varfunction_variables;
        std::vector<std::pair<std::string, double>> varfunction_constants;
        const auto number_fluidphases = phasemanager_->NumFluidPhases();

        for (int i = 0; i < number_fluidphases; i++)
        {
          std::ostringstream temp;
          temp << i + 1;

          varfunction_variables.emplace_back("S" + temp.str(), phasemanager_->Saturation(i));
          varfunction_variables.emplace_back("p" + temp.str(), phasemanager_->Pressure(i));
        }
        varfunction_variables.emplace_back("porosity", phasemanager_->Porosity());

        const auto relative_mobility =
            GLOBAL::Problem::Instance()
                ->FunctionById<CORE::UTILS::FunctionOfAnything>(
                    relative_mobility_funct_id_[current_scalar] - 1)
                .Evaluate(varfunction_variables, varfunction_constants, 0);

        return relative_mobility;
      }

      //! get pre-factor needed for OD-fluid-linearization of convective term
      void GetPreFacLinConvODFluid(const int k, const unsigned ui,
          std::vector<double>* prefaclinconvodfluid, const CORE::LINALG::Matrix<NSD, 1>& gradphi,
          const CORE::LINALG::Matrix<1, NSD>& gradphiTdifftensor,
          const CORE::LINALG::Matrix<NEN, 1>& funct, const CORE::LINALG::Matrix<NSD, NEN>& derxy,
          const int phase)
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
            const CORE::LINALG::Matrix<NSD, 1> gradpres = PressureGradient(phase);
            const double abspressgrad = AbsPressureGradient(phase);

            static CORE::LINALG::Matrix<NSD, NSD> difftensor(true);
            phasemanager_->PermeabilityTensor(phase, difftensor);
            difftensor.Scale(phasemanager_->RelPermeabilityDeriv(phase) /
                             phasemanager_->DynViscosity(phase, abspressgrad, 2));

            static CORE::LINALG::Matrix<1, NSD> gradphiTdifftensor(true);
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
      void GetDiffTensorFluid(
          const int k, CORE::LINALG::Matrix<NSD, NSD>& difftensor, const int phase)
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
      void GetRefGradPres(const int k, const CORE::LINALG::Matrix<NSD, NSD>& xjm,
          CORE::LINALG::Matrix<NSD, 1>& refgradpres, const int phase)
      {
        refgradpres.Clear();

        const int numfluidphases = phasemanager_->NumFluidPhases();
        const int numvolfrac = phasemanager_->NumVolFrac();

        // fluid phase
        if (phase >= 0 && phase < numfluidphases)
        {
          const std::vector<CORE::LINALG::Matrix<NSD, 1>>& fluidgradphi =
              *(variablemanager_->GradPhinp());

          // gradient of phi w.r.t. reference coordinates
          std::vector<CORE::LINALG::Matrix<NSD, 1>> reffluidgradphi(
              numfluidphases, CORE::LINALG::Matrix<NSD, 1>(true));
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
      const CORE::LINALG::Matrix<NSD, 1>& PressureGradient(const int curphase) const
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
          CORE::LINALG::Matrix<NSD, NEN>& xyze)
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
      //! phase pressure
      std::vector<double> pressure_;
      //! phase saturation
      std::vector<double> saturation_;
      //! phase density
      std::vector<double> density_;
      //! solid pressure
      double solidpressure_;
      //! pressure gradient
      std::vector<CORE::LINALG::Matrix<NSD, 1>> pressuregrad_;
      //! pressure gradient
      std::vector<CORE::LINALG::Matrix<NSD, NSD>> difftensorsfluid_;
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

      //! function IDs of relative mobility functions
      std::vector<int> relative_mobility_funct_id_;

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

FOUR_C_NAMESPACE_CLOSE

#endif
