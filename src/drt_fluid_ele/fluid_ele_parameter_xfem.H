/*-----------------------------------------------------------*/
/*! \file

\brief Setting of specific XFEM based fluid parameter for element evaluation

\maintainer Martin Kronbichler

\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FLUID_ELE_PARAMETER_XFEM_H
#define FLUID_ELE_PARAMETER_XFEM_H

#include "fluid_ele_parameter_std.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_cut.H"


namespace DRT
{
  namespace ELEMENTS
  {
    class FluidEleParameterXFEM : public FluidEleParameterStd
    {
     public:
      virtual ~FluidEleParameterXFEM() {}

      /// Singleton access method
      static FluidEleParameterXFEM* Instance(bool create = true);

      /// called upon destruction
      virtual void Done();

      /// set all the XFEM specific parameters
      void SetElementXFEMParameter(Teuchos::ParameterList& params,  ///< parameter list
          int myrank                                                ///< pid (for output purpose)
      );

      /*----------------------------------------------------*/
      //! @name access routines for integration on cut elements
      /*----------------------------------------------------*/
      //@{

      //! get the volumecell integration method used for integrating cut elements?
      INPAR::CUT::VCellGaussPts VolumeCellGaussPoints() const { return vcellgausspts_; };

      //! get the boundarycell integration method used for integrating the surface in cut elements?
      INPAR::CUT::BCellGaussPts BoundaryCellGaussPoints() const { return bcellgausspts_; };

      //@}

      /*----------------------------------------------------*/
      //! @name interface coupling approach
      /*----------------------------------------------------*/
      //@{

      //! coupling approach (Nitsche or hybrid stress-based LM)
      INPAR::XFEM::CouplingMethod GetCouplingMethod() const { return coupling_method_; }

      //! get information, whether L2-projection between stress fields is accomplished on whole cut
      //! element or on physical volume
      INPAR::XFEM::Hybrid_LM_L2_Proj HybridLM_L2Proj() const { return hybrid_lm_l2_proj_; }

      //@}

      /*----------------------------------------------------*/
      //! @name access routines for the viscous stabilization in Nitsche's method and MixedHybrid_LM
      //! methods
      /*----------------------------------------------------*/
      //@{

      //! get the type of how to estimate the scaling of the trace inequality used for the viscous
      //! part of Nitsche's method?
      INPAR::XFEM::ViscStab_TraceEstimate ViscStabTracEstimate() const
      {
        return visc_stab_trace_estimate_;
      };

      //! get the element length definition used for viscous part of the penalty term in Nitsche's
      //! method
      INPAR::XFEM::ViscStab_hk ViscStabHK() const { return visc_stab_hk_; };

      //! get the dimensionless user defined scaling for the penalty term in Nitsche's method and
      //! scaling factor for the MHVS method (then gamma = 1/n, see publications)
      double NITStabScaling() const { return nit_stab_gamma_; };

      //! get Nitsche's penalty scaling for tangential terms
      double NITStabScalingTang() const { return nit_stab_gamma_tang_; };

      //! get information, whether the formulation should be symmetric/skew-symmetric in the adjoint
      //! viscous terms
      //! @return true, in case of a symmetric adjoint term
      bool IsViscousAdjointSymmetric() const
      {
        if (visc_adjoint_scaling_ == INPAR::XFEM::adj_none)
          dserror("Do not call IsViscousAdjointSymmetric with adj_none");
        return visc_adjoint_scaling_ == INPAR::XFEM::adj_sym;
      }

      //! get information, whether the formulation should be symmetric/skew-symmetric/none in the
      //! adjoint viscous terms
      //! @return double with scaling
      double GetViscousAdjointScaling() const
      {
        switch (visc_adjoint_scaling_)
        {
          case INPAR::XFEM::adj_none:
            return 0.0;
          case INPAR::XFEM::adj_sym:
            return 1.0;
          case INPAR::XFEM::adj_skew:
            return -1.0;
          default:
            dserror("Unkown type of AdjointScaling!");
        }
        return 0.0;  // make compiler happy
      }

      //! get the flag if the simulation is run as pseudo 2D simulation with only one element in the
      //! third dimension and strong Dirichlet condition to fix u_z = 0
      bool IsPseudo2D() const { return is_pseudo_2D_; };

      //@}


      /*----------------------------------------------------*/
      //! @name access routines for convective inflow stabilization
      //! parameters in Nitsche's method and MixedHybrid_LM methods
      /*----------------------------------------------------*/
      //@{
      //! get the type of scaling for convective/inflow stabilization term for xfluid-fluid problems
      INPAR::XFEM::XFF_ConvStabScaling XffConvStabScaling() const { return xff_conv_stab_scaling_; }

      //! get the type of scaling for convective/inflow stabilization term for classical xfluid
      //! problem
      INPAR::XFEM::ConvStabScaling ConvStabScaling() const { return conv_stab_scaling_; }

      //@}

      /*----------------------------------------------------*/
      //! @name access routines for definitions to compute
      //! Nitsche's penalty parameter for different flow
      //! regimes
      /*----------------------------------------------------*/
      //@{
      //! get information, whether we take the maximum from the viscous and convective penalty
      //! scaling or the sum of them
      INPAR::XFEM::MassConservationCombination MassConservationCombination() const
      {
        return mass_conservation_combo_;
      }

      //! get the type of scaling for convective/inflow stabilization term for classical xfluid
      //! problem
      INPAR::XFEM::MassConservationScaling MassConservationScaling() const
      {
        return mass_conservation_scaling_;
      }

      //@}

      //! for new OST-implementation: which interface terms to be evaluated for previous time step
      INPAR::XFEM::InterfaceTermsPreviousState InterfaceTermsPreviousState() const
      {
        return intterms_prev_state_;
      }

      //! assure a valid combination of input parameters for certain weighting
      void CheckParameterConsistencyForAveragingStrategy(
          int myrank, INPAR::XFEM::AveragingStrategy averaging_strategy) const;

     private:
      //! assure a valid combination of input parameters
      void CheckParameterConsistency(int myrank) const;

      /*----------------------------------------------------*/
      //! @name parameters for integration on cut elements
      /*----------------------------------------------------*/
      //@{

      //! which volumecell integration is used for integrating cut elements?
      INPAR::CUT::VCellGaussPts vcellgausspts_;

      //! which boundarycell integration is used for integrating the surface in cut elements?
      INPAR::CUT::BCellGaussPts bcellgausspts_;

      //@}


      /*----------------------------------------------------*/
      //! @name parameters describing the coupling approach
      /*----------------------------------------------------*/
      //@{

      //! coupling approach (Nitsche or hybrid stress-based LM)
      INPAR::XFEM::CouplingMethod coupling_method_;

      //! for coupling using stress-based LM: L2-projection between stress fields on whole cut
      //! element or on physical volume
      INPAR::XFEM::Hybrid_LM_L2_Proj hybrid_lm_l2_proj_;

      //@}


      /*----------------------------------------------------*/
      //! @name parameters for the viscous stabilization in Nitsche's method and MixedHybrid_LM
      //! methods
      /*----------------------------------------------------*/
      //@{

      //! how to estimate the scaling of the trace inequality used for the viscous part of Nitsche's
      //! method?
      INPAR::XFEM::ViscStab_TraceEstimate visc_stab_trace_estimate_;

      //! element length definition used for viscous part of the penalty term in Nitsche's method
      INPAR::XFEM::ViscStab_hk visc_stab_hk_;

      //! dimensionless user defined scaling for the penalty term in Nitsche's method and
      //! scaling factor for the MHVS method (then gamma = 1/n, see publications),
      //! this parameter shall be independent of element type, shape and polynomial degree
      double nit_stab_gamma_;

      //! scaling for penalty term in Nitsche's method in case normal and tangential direction are
      //! split
      double nit_stab_gamma_tang_;

      //! enum, that indicates matching signs between the viscous and adjoint viscous Nitsche terms
      //! or between the equivalent MHVS-terms
      INPAR::XFEM::AdjointScaling visc_adjoint_scaling_;

      //@}


      /*----------------------------------------------------*/
      //! @name pseudo 2D flag
      /*----------------------------------------------------*/
      //@{

      //! pseudo 2D flag for 2D simulation with one element in z-direction
      double is_pseudo_2D_;

      //@}

      /*----------------------------------------------------*/
      //! @name parameters for the convective interface stabilizations
      /*----------------------------------------------------*/
      //@{

      //! type of convective scaling for xfluid-fluid problem
      INPAR::XFEM::XFF_ConvStabScaling xff_conv_stab_scaling_;

      //! type of convective scaling for xfluid problem
      INPAR::XFEM::ConvStabScaling conv_stab_scaling_;

      //@}

      /*----------------------------------------------------*/
      //! @name definitions to compute Nitsche's penalty
      //! parameter for different flow regimes
      /*----------------------------------------------------*/
      //@{

      //! get information, whether we take the maximum from the viscous and convective penalty
      //! scaling or the sum of them
      INPAR::XFEM::MassConservationCombination mass_conservation_combo_;

      //! get the type of scaling for convective/inflow stabilization term for classical xfluid
      //! problem
      INPAR::XFEM::MassConservationScaling mass_conservation_scaling_;

      //@}

      //! which interface terms to be evaluated for previous time step (new OST)
      INPAR::XFEM::InterfaceTermsPreviousState intterms_prev_state_;

     protected:
      /// protected Constructor since we are a Singleton.
      FluidEleParameterXFEM();
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
