/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_parameter_intface.H

\brief Setting of general fluid parameter for internal faces evaluation

This file has to contain all parameters called in fluid_ele_intfaces.cpp.
Additional parameters required in derived classes of FluidIntFaceImpl have to
be set in problem specific parameter lists derived from this class.

\level 2

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249

*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_PARAMETER_INTFACE_H
#define FLUID_ELE_PARAMETER_INTFACE_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_singletondestruction.H"
//#include "../drt_fluid_ele/fluid_ele_parameter_std.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_xfem.H"
//#include "../drt_inpar/inpar_turbulence.H"
#include "fluid_ele_parameter_timint.H"


namespace DRT
{
  namespace ELEMENTS
  {
    class FluidEleParameterIntFace : public SingletonDestruction
    {
     public:
      //! Singleton access method
      static FluidEleParameterIntFace* Instance(bool create = true);

      //! called upon destruction
      void Done();


      ~FluidEleParameterIntFace() {}



      /*========================================================================*/
      //! @name set-routines, individually for each face
      /*========================================================================*/

      //! set the EOS pressure stabilization flag
      void Set_Face_EOS_Pres(const bool is_face_EOS_Pres) { is_face_EOS_Pres_ = is_face_EOS_Pres; };

      //! set the EOS convective streamline stabilization flag
      void Set_Face_EOS_Conv_Stream(const bool is_face_EOS_Conv_Stream)
      {
        is_face_EOS_Conv_Stream_ = is_face_EOS_Conv_Stream;
      };

      //! set the EOS convective cross wind stabilization flag
      void Set_Face_EOS_Conv_Cross(const bool is_face_EOS_Conv_Cross)
      {
        is_face_EOS_Conv_Cross_ = is_face_EOS_Conv_Cross;
      };

      //! set the EOS divergence vel-jump  flag
      void Set_Face_EOS_Div_vel_jump(const bool is_face_EOS_Div_vel_jump)
      {
        is_face_EOS_Div_vel_jump_ = is_face_EOS_Div_vel_jump;
      };

      //! set the EOS divergence div-jump flag
      void Set_Face_EOS_Div_div_jump(const bool is_face_EOS_Div_div_jump)
      {
        is_face_EOS_Div_div_jump_ = is_face_EOS_Div_div_jump;
      };

      //! set the viscous ghost-penalty stabilization flag
      void Set_Face_GP_visc(const bool is_face_GP_visc) { is_face_GP_visc_ = is_face_GP_visc; };

      //! set the transient ghost-penalty stabilization flag
      void Set_Face_GP_trans(const bool is_face_GP_trans) { is_face_GP_trans_ = is_face_GP_trans; };

      //! set the 2nd order derivatives ghost-penalty stabilization flag
      void Set_Face_GP_u_p_2nd(const bool is_face_GP_u_p_2nd)
      {
        is_face_GP_u_p_2nd_ = is_face_GP_u_p_2nd;
      };

      //! set the EOS pattern for the assembly of the current face
      void Set_Face_EOS_GP_Pattern(const INPAR::FLUID::EOS_GP_Pattern face_eos_gp_pattern)
      {
        face_eos_gp_pattern_ = face_eos_gp_pattern;
      };

      //! general fluid parameters are set
      void SetFaceGeneralFluidParameter(Teuchos::ParameterList& params,  //> parameter list
          int myrank);                                                   //> proc id

      //! general xfem (ghost-penalty) parameters are set
      void SetFaceGeneralXFEMParameter(Teuchos::ParameterList& params,  //> parameter list
          int myrank);                                                  //> proc id


      //! specific fluid xfem (ghost-penalty) parameters are set for the face and return if
      //! stabilization for current face is required
      bool SetFaceSpecificFluidXFEMParameter(
          const INPAR::XFEM::FaceType& face_type,  ///< which type of face std, ghost, ghost-penalty
          Teuchos::ParameterList& params           ///< parameter list
      );

      //! set flag if step is a ghost-penalty reconstruction step for xfluid time integration
      void Set_GhostPenaltyReconstruction(const bool is_ghost_penalty_reconstruction_step)
      {
        is_ghost_penalty_reconstruction_step_ = is_ghost_penalty_reconstruction_step;
      };



      /*========================================================================*/
      //! @name access-routines
      /*========================================================================*/

      /*----------------------------------------------------*/
      //! @name general parameters
      /*----------------------------------------------------*/

      //! Flag for physical type of the fluid flow (incompressible, loma, varying_density,
      //! Boussinesq, poro)
      INPAR::FLUID::PhysicalType PhysicalType() const { return physicaltype_; };

      //! Return function number of Oseen advective field
      int OseenFieldFuncNo() const { return oseenfieldfuncno_; };


      /*----------------------------------------------------*/
      //! @name general stabilization parameters (non-individually for each face, see below)
      /*----------------------------------------------------*/

      /// parameter for edge-based (EOS,CIP) stabilizations
      //! Flag to (de)activate pressure stabilization
      INPAR::FLUID::EOS_Pres EOS_Pres() const { return EOS_pres_; };
      //! Flag to (de)activate convective streamline stabilization
      INPAR::FLUID::EOS_Conv_Stream EOS_Conv_Stream() const { return EOS_conv_stream_; };
      //! Flag to (de)activate convective crosswind stabilization
      INPAR::FLUID::EOS_Conv_Cross EOS_Conv_Cross() const { return EOS_conv_cross_; };
      //! Flag to (de)activate divergence stabilization
      INPAR::FLUID::EOS_Div EOS_Div() const { return EOS_div_; };
      //! Flag to define element length
      INPAR::FLUID::EOS_ElementLength EOS_element_length() const { return EOS_element_length_; };
      //! Flag to define tau for edge-based stabilization
      INPAR::FLUID::EOS_TauType EOS_WhichTau() const { return EOS_whichtau_actual_; };

      bool Is_EOS_Pres() const { return EOS_pres_ != INPAR::FLUID::EOS_PRES_none; };
      //! Flag to (de)activate convective streamline stabilization
      bool Is_Conv_Stream() const { return EOS_conv_stream_; };
      //! Flag to (de)activate convective crosswind stabilization
      bool Is_Conv_Cross() const { return EOS_conv_cross_; };
      //! Flag to (de)activate divergence stabilization
      bool Is_Div() const { return EOS_div_; };
      //! Flag to define element length
      bool Is_element_length() const { return EOS_element_length_; };
      //! Flag to define tau for edge-based stabilization
      bool Is_WhichTau() const { return EOS_whichtau_actual_; };

      //! Flag to activate special least-squares condition for pseudo 2D examples where pressure
      //! level is determined via Krylov-projection
      bool presKrylov2Dz() const { return presKrylov2Dz_; };

      //! get the viscous ghost-penalty stabilization flag
      bool Is_General_Ghost_Penalty_visc() { return ghost_penalty_visc_; };

      //! get the transient ghost-penalty stabilization flag
      bool Is_General_Ghost_Penalty_trans() { return ghost_penalty_trans_; };

      //! get the 2nd order derivatives ghost-penalty stabilization flag
      bool Is_General_Ghost_Penalty_u_p_2nd() { return ghost_penalty_u_p_2nd_; };

      //! get the 2nd order derivatives ghost-penalty stabilization flag - normal derivatives or
      //! full derivatives
      bool Is_General_Ghost_Penalty_u_p_2nd_Normal() { return ghost_penalty_u_p_2nd_normal_; };

      //! get the viscous ghost-penalty stabilization factor
      double Ghost_Penalty_visc_fac() { return ghost_penalty_visc_fac; };

      //! get the transient ghost-penalty stabilization factor
      double Ghost_Penalty_trans_fac() { return ghost_penalty_trans_fac; };

      //! get the viscous ghost-penalty stabilization factor
      double Ghost_Penalty_visc_2nd_fac() { return ghost_penalty_visc_2nd_fac; };

      //! get the viscous ghost-penalty stabilization factor
      double Ghost_Penalty_press_2nd_fac() { return ghost_penalty_press_2nd_fac; };
      //
      /*----------------------------------------------------*/
      //! @name stabilization parameters, individually set for each face
      /*----------------------------------------------------*/

      //! get the EOS pressure stabilization flag
      bool Face_EOS_Pres() { return is_face_EOS_Pres_; };

      //! get the EOS convective streamline stabilization flag
      bool Face_EOS_Conv_Stream() { return is_face_EOS_Conv_Stream_; };

      //! get the EOS convective cross wind stabilization flag
      bool Face_EOS_Conv_Cross() { return is_face_EOS_Conv_Cross_; };

      //! get the EOS divergence vel-jump  flag
      bool Face_EOS_Div_vel_jump() { return is_face_EOS_Div_vel_jump_; };

      //! get the EOS divergence div-jump flag
      bool Face_EOS_Div_div_jump() { return is_face_EOS_Div_div_jump_; };

      //! get the viscous ghost-penalty stabilization flag
      bool Face_GP_visc() { return is_face_GP_visc_; };

      //! get the transient ghost-penalty stabilization flag
      bool Face_GP_trans() { return is_face_GP_trans_; };

      //! get the 2nd order derivatives ghost-penalty stabilization flag
      bool Face_GP_u_p_2nd() { return is_face_GP_u_p_2nd_; };

      //! get the EOS pattern for the assembly of the current face
      INPAR::FLUID::EOS_GP_Pattern Face_EOS_GP_Pattern() { return face_eos_gp_pattern_; };

      //! get flag if step is a ghost-penalty reconstruction step for xfluid time integration
      bool Is_GhostPenaltyReconstruction() { return is_ghost_penalty_reconstruction_step_; };


     protected:
      /*----------------------------------------------------*/
      //! @name general parameters
      /*----------------------------------------------------*/

      //! Flag SetFaceGeneralFluidParameter was called
      bool set_face_general_fluid_parameter_;

      //! Flag SetFaceGeneralXFEMParameter was called
      bool set_face_general_XFEM_parameter_;

      //! Flag for physical type of the fluid flow (incompressible, loma, varying_density,
      //! Boussinesq, Poro)
      INPAR::FLUID::PhysicalType physicaltype_;

      //! Which stabilization type
      INPAR::FLUID::StabType stabtype_;

      //! function number for advective velocity for Oseen problems
      int oseenfieldfuncno_;


      /*----------------------------------------------------*/
      //! @name general edge-based/face-oriented fluid stabilization parameters (not indiviually set
      //! for each face, see below)
      /*----------------------------------------------------*/

      //! which EOS pressure stabilization
      INPAR::FLUID::EOS_Pres EOS_pres_;

      //! which EOS convective streamline stabilization
      INPAR::FLUID::EOS_Conv_Stream EOS_conv_stream_;

      //! which EOS convective crosswind stabilization
      INPAR::FLUID::EOS_Conv_Cross EOS_conv_cross_;

      //! which EOS divergence stabilization
      INPAR::FLUID::EOS_Div EOS_div_;


      //! which EOS stabilization parameter definition
      INPAR::FLUID::EOS_TauType EOS_whichtau_;
      INPAR::FLUID::EOS_TauType EOS_whichtau_actual_;

      //! which EOS characteristic element length definition
      INPAR::FLUID::EOS_ElementLength EOS_element_length_;

      //! flag to active special least-squares condition for pseudo 2D examples where pressure level
      //! is determined via Krylov-projection
      bool presKrylov2Dz_;


      /*----------------------------------------------------*/
      //! @name general XFEM ghost penalty stabilization parameters (not indiviually set for each
      //! face, see below)
      /*----------------------------------------------------*/

      //! factor for viscous ghost penalty XFEM stabilization terms
      double ghost_penalty_visc_fac;

      //! factor for transient ghost penalty XFEM stabilization terms
      double ghost_penalty_trans_fac;

      //! general flag for viscous ghost-penalty stabilization terms
      bool ghost_penalty_visc_;

      //! general flag for transient ghost-penalty stabilization terms
      bool ghost_penalty_trans_;

      //! general flag for 2nd order derivatives in u- and p- ghost-penalty stabilization terms
      bool ghost_penalty_u_p_2nd_;

      //! general flag for 2nd order derivatives in u- and p- ghost-penalty stabilization terms
      bool ghost_penalty_u_p_2nd_normal_;

      //! factor for viscous ghost penalty XFEM stabilization terms
      double ghost_penalty_visc_2nd_fac;

      //! factor for pressure ghost penalty XFEM stabilization terms
      double ghost_penalty_press_2nd_fac;

      /*----------------------------------------------------*/
      //! @name EOS stabilization parameter, individually set for each face
      /*----------------------------------------------------*/

      bool is_face_EOS_Pres_;
      bool is_face_EOS_Conv_Stream_;
      bool is_face_EOS_Conv_Cross_;
      bool is_face_EOS_Div_vel_jump_;
      bool is_face_EOS_Div_div_jump_;

      /*----------------------------------------------------*/
      //! @name Ghost-Penalty parameter, individually set for each face
      /*----------------------------------------------------*/
      bool is_face_GP_visc_;
      bool is_face_GP_trans_;
      bool is_face_GP_u_p_2nd_;


      /*----------------------------------------------------*/
      //! @name Combined EOS/Ghost-Penalty parameters, individually set for each face
      /*----------------------------------------------------*/

      INPAR::FLUID::EOS_GP_Pattern face_eos_gp_pattern_;

      /*----------------------------------------------------*/
      //! @name Flag for ghost-penalty reconstruction steps
      /*----------------------------------------------------*/
      bool is_ghost_penalty_reconstruction_step_;

      /// private Constructor since we are a Singleton.
      FluidEleParameterIntFace();

     private:
      //! access time-integration parameters
      DRT::ELEMENTS::FluidEleParameterTimInt* fldparatimint_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
