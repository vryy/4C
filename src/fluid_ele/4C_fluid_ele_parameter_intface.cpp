/*----------------------------------------------------------------------*/
/*! \file

\brief Setting of general fluid parameter for internal faces evaluation

This file has to contain all parameters called in fluid_ele_intfaces.cpp.
Additional parameters required in derived classes of FluidIntFaceImpl have to
be set in problem specific parameter lists derived from this class.

\level 2


*/
/*----------------------------------------------------------------------*/

#include "4C_fluid_ele_parameter_intface.hpp"

#include "4C_io_pstream.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_TimeMonitor.hpp>

#include <iostream>
#include <string>

FOUR_C_NAMESPACE_OPEN


Discret::ELEMENTS::FluidEleParameterIntFace* Discret::ELEMENTS::FluidEleParameterIntFace::Instance(
    Core::UTILS::SingletonAction action)
{
  static auto singleton_owner = Core::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<Discret::ELEMENTS::FluidEleParameterIntFace>(
            new Discret::ELEMENTS::FluidEleParameterIntFace());
      });

  return singleton_owner.Instance(action);
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
Discret::ELEMENTS::FluidEleParameterIntFace::FluidEleParameterIntFace()
    : set_face_general_fluid_parameter_(false),
      set_face_general_XFEM_parameter_(false),
      physicaltype_(Inpar::FLUID::physicaltype_undefined),
      stabtype_(Inpar::FLUID::stabtype_edgebased),
      oseenfieldfuncno_(-1),
      EOS_pres_(Inpar::FLUID::EOS_PRES_none),
      EOS_conv_stream_(Inpar::FLUID::EOS_CONV_STREAM_none),
      EOS_conv_cross_(Inpar::FLUID::EOS_CONV_CROSS_none),
      EOS_div_(Inpar::FLUID::EOS_DIV_none),
      EOS_whichtau_(Inpar::FLUID::EOS_tau_burman_fernandez),
      EOS_element_length_(Inpar::FLUID::EOS_he_max_dist_to_opp_surf),
      presKrylov2Dz_(false),
      ghost_penalty_visc_fac_(0.0),
      ghost_penalty_trans_fac_(0.0),
      ghost_penalty_visc_(false),
      ghost_penalty_trans_(false),
      ghost_penalty_u_p_2nd_(false),
      ghost_penalty_u_p_2nd_normal_(false),
      ghost_penalty_visc_2nd_fac_(0.0),
      ghost_penalty_press_2nd_fac_(0.0),
      is_face_EOS_Pres_(false),
      is_face_EOS_Conv_Stream_(false),
      is_face_EOS_Conv_Cross_(false),
      is_face_EOS_Div_vel_jump_(false),
      is_face_EOS_Div_div_jump_(false),
      is_face_GP_visc_(false),
      is_face_GP_trans_(false),
      is_face_GP_u_p_2nd_(false),
      face_eos_gp_pattern_(Inpar::FLUID::EOS_GP_Pattern_up),
      is_ghost_penalty_reconstruction_step_(false)
{
  // we have to know the time parameters here to check for illegal combinations
  fldparatimint_ = Discret::ELEMENTS::FluidEleParameterTimInt::Instance();
}

//----------------------------------------------------------------------*
//  set general parameters                                 schott Jun14 |
//---------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidEleParameterIntFace::set_face_general_fluid_parameter(
    Teuchos::ParameterList& params, int myrank)
{
  if (set_face_general_fluid_parameter_ == false) set_face_general_fluid_parameter_ = true;
  // For turbulent inflow generation,
  // this function is indeed two times called.
  // In this sepcial case, calling this function twice
  // is ok!
  else
  {
    if (myrank == 0)
      std::cout
          << std::endl
          << (" Warning: general face fluid XFEM parameters should be set only once!!\n "
              " If you run a turbulent inflow generation, calling this function twice is ok!\n ")
          << std::endl
          << std::endl;
    //    FOUR_C_THROW(" general face fluid XFEM parameters should be set only once!! -> Check
    //    this?!");
  }


  // set flag for physical type of fluid flow
  physicaltype_ = Core::UTILS::GetAsEnum<Inpar::FLUID::PhysicalType>(params, "Physical Type");
  if ((physicaltype_ != Inpar::FLUID::incompressible) and
      (physicaltype_ != Inpar::FLUID::stokes) and (physicaltype_ != Inpar::FLUID::oseen) and
      (physicaltype_ != Inpar::FLUID::poro))
    FOUR_C_THROW("physical type is not supported for face stabilizations.");

  // get function number of given Oseen advective field if necessary
  if (physicaltype_ == Inpar::FLUID::oseen) oseenfieldfuncno_ = params.get<int>("OSEENFIELDFUNCNO");

  //---------------------------------
  // which basic stabilization type?
  // residual-based: for residualbased standard fluid or residual-based XFEM fluid in combination
  // with edge-based ghost penalty stabilization edge-based:     for pure edge-based/ghost-penalty
  // stabilization
  stabtype_ = Core::UTILS::GetAsEnum<Inpar::FLUID::StabType>(params, "STABTYPE");

  // --------------------------------
  // edge-based fluid stabilization can be used as standard fluid stabilization or
  // as ghost-penalty stabilization in addition to residual-based stabilizations in the XFEM

  // set parameters if single stabilization terms are switched on/off or which type of stabilization
  // is chosen

  Teuchos::ParameterList& stablist_edgebased = params.sublist("EDGE-BASED STABILIZATION");

  EOS_pres_ = Core::UTILS::IntegralValue<Inpar::FLUID::EosPres>(stablist_edgebased, "EOS_PRES");
  EOS_conv_stream_ = Core::UTILS::IntegralValue<Inpar::FLUID::EosConvStream>(
      stablist_edgebased, "EOS_CONV_STREAM");
  EOS_conv_cross_ =
      Core::UTILS::IntegralValue<Inpar::FLUID::EosConvCross>(stablist_edgebased, "EOS_CONV_CROSS");
  EOS_div_ = Core::UTILS::IntegralValue<Inpar::FLUID::EosDiv>(stablist_edgebased, "EOS_DIV");

  if (physicaltype_ == Inpar::FLUID::stokes and EOS_conv_stream_)
    FOUR_C_THROW("no EOS_CONV_STREAM stabilization required for Stokes problems");
  if (physicaltype_ == Inpar::FLUID::stokes and EOS_conv_cross_)
    FOUR_C_THROW("no EOS_CONV_CROSS stabilization required for Stokes problems");

  // activate special least-squares condition for pseudo 2D examples where pressure level is
  // determined via Krylov-projection
  presKrylov2Dz_ = (bool)Core::UTILS::IntegralValue<int>(stablist_edgebased, "PRES_KRYLOV_2Dz");

  // check for reasonable combinations of non-edgebased fluid stabilizations with edge-based
  // stabilizations
  if (stabtype_ != Inpar::FLUID::stabtype_edgebased)
  {
    // case of xfem check whether additional xfem-stabilization terms in the form of
    // edge-based terms are activated (i.e., ghost penalties)

    if (EOS_pres_ != Inpar::FLUID::EOS_PRES_xfem_gp and EOS_pres_ != Inpar::FLUID::EOS_PRES_none)
      FOUR_C_THROW(
          "check the combination of edge-based pressure-EOS and non-edge-based stabtype! Do you "
          "really want to use this combination?");
    if (EOS_conv_stream_ != Inpar::FLUID::EOS_CONV_STREAM_xfem_gp and
        EOS_conv_stream_ != Inpar::FLUID::EOS_CONV_STREAM_none)
      FOUR_C_THROW(
          "check the combination of edge-based pressure-EOS and non-edge-based stabtype! Do you "
          "really want to use this combination?");
    if (EOS_conv_cross_ != Inpar::FLUID::EOS_CONV_CROSS_xfem_gp and
        EOS_conv_cross_ != Inpar::FLUID::EOS_CONV_CROSS_none)
      FOUR_C_THROW(
          "check the combination of edge-based pressure-EOS and non-edge-based stabtype! Do you "
          "really want to use this combination?");
    if (EOS_div_ != Inpar::FLUID::EOS_DIV_div_jump_xfem_gp and
        EOS_div_ != Inpar::FLUID::EOS_DIV_vel_jump_xfem_gp and
        EOS_div_ != Inpar::FLUID::EOS_DIV_none)
      FOUR_C_THROW(
          "check the combination of edge-based pressure-EOS and non-edge-based stabtype! Do you "
          "really want to use this combination?");

    if (presKrylov2Dz_)
      FOUR_C_THROW(
          "pressure Krylov 2Dz condition not reasonable for non-pure edge-based stabilizations");
  }

  if (presKrylov2Dz_ and EOS_pres_ != Inpar::FLUID::EOS_PRES_std_eos)
    FOUR_C_THROW(
        "pressure Krylov 2Dz condition only reasonable for full p-EOS: EOS_PRES = std_eos");

  EOS_element_length_ = Core::UTILS::IntegralValue<Inpar::FLUID::EosElementLength>(
      stablist_edgebased, "EOS_H_DEFINITION");
  EOS_whichtau_ = Core::UTILS::IntegralValue<Inpar::FLUID::EosTauType>(
      stablist_edgebased, "EOS_DEFINITION_TAU");


  // set correct stationary definition of stabilization parameter automatically
  if (fldparatimint_->IsStationary())
  {
    if (EOS_whichtau_ == Inpar::FLUID::EOS_tau_burman_fernandez_hansbo)
      EOS_whichtau_ = Inpar::FLUID::EOS_tau_burman_fernandez_hansbo_wo_dt;
    if (EOS_whichtau_ == Inpar::FLUID::EOS_tau_burman_hansbo_dangelo_zunino)
      EOS_whichtau_ = Inpar::FLUID::EOS_tau_burman_hansbo_dangelo_zunino_wo_dt;
    if (EOS_whichtau_ == Inpar::FLUID::EOS_tau_schott_massing_burman_dangelo_zunino)
      EOS_whichtau_ = Inpar::FLUID::EOS_tau_schott_massing_burman_dangelo_zunino_wo_dt;
  }

  return;
}


//----------------------------------------------------------------------*
//  set general parameters                                 schott Jun14 |
//---------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidEleParameterIntFace::set_face_general_xfem_parameter(
    Teuchos::ParameterList& params, int myrank)
{
  if (set_face_general_XFEM_parameter_ == false)
    set_face_general_XFEM_parameter_ = true;
  else
  {
    FOUR_C_THROW(" general face fluid XFEM parameters should be set only once!! -> Check this?!");
  }

  // --------------------------------
  // XFEM-specific ghost-penalty stabilizations

  Teuchos::ParameterList& stablist_xfem = params.sublist("XFLUID DYNAMIC/STABILIZATION");

  ghost_penalty_visc_fac_ = stablist_xfem.get<double>("GHOST_PENALTY_FAC", 0.0);
  ghost_penalty_trans_fac_ = stablist_xfem.get<double>("GHOST_PENALTY_TRANSIENT_FAC", 0.0);

  ghost_penalty_visc_ = (bool)Core::UTILS::IntegralValue<int>(stablist_xfem, "GHOST_PENALTY_STAB");
  ghost_penalty_trans_ =
      (bool)Core::UTILS::IntegralValue<int>(stablist_xfem, "GHOST_PENALTY_TRANSIENT_STAB");
  ghost_penalty_visc_2nd_fac_ = stablist_xfem.get<double>("GHOST_PENALTY_2nd_FAC", 0.0);
  ghost_penalty_press_2nd_fac_ = stablist_xfem.get<double>("GHOST_PENALTY_PRESSURE_2nd_FAC", 0.0);

  // safety check
  if (fldparatimint_->IsStationary() and ghost_penalty_trans_)
    FOUR_C_THROW("Do not use transient ghost penalties for stationary problems");

  ghost_penalty_u_p_2nd_ =
      (bool)Core::UTILS::IntegralValue<int>(stablist_xfem, "GHOST_PENALTY_2nd_STAB");
  ghost_penalty_u_p_2nd_normal_ =
      (bool)Core::UTILS::IntegralValue<int>(stablist_xfem, "GHOST_PENALTY_2nd_STAB_NORMAL");

  return;
}

//----------------------------------------------------------------------*
// specific fluid xfem (ghost-penalty) parameters are set for the face  |
// and return if stabilization for current face is required             |
//                                                         schott Jun14 |
//---------------------------------------------------------------------*/
bool Discret::ELEMENTS::FluidEleParameterIntFace::set_face_specific_fluid_xfem_parameter(
    const Inpar::XFEM::FaceType& face_type,  ///< which type of face std, ghost, ghost-penalty
    Teuchos::ParameterList& params           ///< parameter list
)
{
  TEUCHOS_FUNC_TIME_MONITOR("XFEM::Edgestab EOS: set_face_specific_fluid_xfem_parameter");

  set_ghost_penalty_reconstruction(params.get<bool>("ghost_penalty_reconstruct", false));

  // do not assemble non-ghost-penalty faces if not necessary!
  // remark: flags for ghost penalty have to be set
  if (is_ghost_penalty_reconstruction())
  {
    if (face_type != Inpar::XFEM::face_type_ghost_penalty)
      return false;  // do not assembly (stabilizing just the ghost-penalty faces is sufficient)

    // otherwise continue and set the required flags for the face
  }

  //----------------------------------------------------------------------------
  // decide which terms have to be assembled and decide the assembly pattern

  // final decision which terms are assembled for the current face
  // set the values into the fldpara_intface_-list

  EOS_whichtau_actual_ = EOS_whichtau_;
  if (face_type == Inpar::XFEM::face_type_std)
  {
    Set_Face_EOS_Pres((EOS_Pres() == Inpar::FLUID::EOS_PRES_std_eos));
    set_face_eos_conv_stream((EOS_Conv_Stream() == Inpar::FLUID::EOS_CONV_STREAM_std_eos));
    set_face_eos_conv_cross((EOS_Conv_Cross() == Inpar::FLUID::EOS_CONV_CROSS_std_eos));
    set_face_eos_div_vel_jump((EOS_Div() == Inpar::FLUID::EOS_DIV_vel_jump_std_eos));
    set_face_eos_div_div_jump((EOS_Div() == Inpar::FLUID::EOS_DIV_div_jump_std_eos));

    Set_Face_GP_visc(false);
    Set_Face_GP_trans(false);
    Set_Face_GP_u_p_2nd(false);
  }
  else if (face_type == Inpar::XFEM::face_type_ghost_penalty)
  {
    Set_Face_EOS_Pres((EOS_Pres() != Inpar::FLUID::EOS_PRES_none));
    set_face_eos_conv_stream((EOS_Conv_Stream() != Inpar::FLUID::EOS_CONV_STREAM_none));
    set_face_eos_conv_cross((EOS_Conv_Cross() != Inpar::FLUID::EOS_CONV_CROSS_none));
    set_face_eos_div_vel_jump((EOS_Div() == Inpar::FLUID::EOS_DIV_vel_jump_std_eos or
                               EOS_Div() == Inpar::FLUID::EOS_DIV_vel_jump_xfem_gp));
    set_face_eos_div_div_jump((EOS_Div() == Inpar::FLUID::EOS_DIV_div_jump_std_eos or
                               EOS_Div() == Inpar::FLUID::EOS_DIV_div_jump_xfem_gp));

    Set_Face_GP_visc(is_general_ghost_penalty_visc());
    Set_Face_GP_trans(is_general_ghost_penalty_trans());
    Set_Face_GP_u_p_2nd(is_general_ghost_penalty_u_p_2nd());
  }
  else if (face_type == Inpar::XFEM::face_type_boundary_ghost_penalty)
  {
    // TODO: this can be improved if only pressure is assembled later on

    Set_Face_EOS_Pres((EOS_Pres() == Inpar::FLUID::EOS_PRES_xfem_gp));
    set_face_eos_conv_stream(false);
    set_face_eos_conv_cross(false);
    set_face_eos_div_vel_jump(false);
    set_face_eos_div_div_jump(false);

    Set_Face_GP_visc(false);
    Set_Face_GP_trans(false);
    Set_Face_GP_u_p_2nd(false);
  }
  else if (face_type == Inpar::XFEM::face_type_ghost)
  {
    Set_Face_EOS_Pres(false);
    set_face_eos_conv_stream(false);
    set_face_eos_conv_cross(false);
    set_face_eos_div_vel_jump(false);
    set_face_eos_div_div_jump(false);
    Set_Face_GP_visc(false);
    Set_Face_GP_trans(false);
    Set_Face_GP_u_p_2nd(false);
  }
  else if (face_type == Inpar::XFEM::face_type_porof)
  {
    EOS_whichtau_actual_ = Inpar::FLUID::EOS_tau_poroelast_fluid;
    Set_Face_EOS_Pres(true);
    set_face_eos_conv_stream(false);
    set_face_eos_conv_cross(false);
    set_face_eos_div_vel_jump(false);
    set_face_eos_div_div_jump(true);
    Set_Face_GP_visc(false);
    Set_Face_GP_trans(false);
    Set_Face_GP_u_p_2nd(false);
  }
  else
    FOUR_C_THROW("unknown face_type!!!");

  // which pattern has to be activated?
  // TODO: this can be improved if only pressure is assembled and so on!

  if (face_eos_div_div_jump())
  {
    set_face_eos_gp_pattern(Inpar::FLUID::EOS_GP_Pattern_up);
  }
  else
  {
    set_face_eos_gp_pattern(Inpar::FLUID::EOS_GP_Pattern_uvwp);
  }

  // return false if no stabilization is required
  if (!Face_EOS_Pres() and !face_eos_conv_stream() and !Face_EOS_Conv_Cross() and
      !face_eos_div_vel_jump() and !face_eos_div_div_jump() and !Face_GP_visc() and
      !Face_GP_trans() and !Face_GP_u_p_2nd())
  {
    return false;
  }

  return true;
}

FOUR_C_NAMESPACE_CLOSE
