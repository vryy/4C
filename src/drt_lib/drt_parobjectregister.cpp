
#include <sstream>
#include <string>
#include <iostream>

#include "drt_parobjectregister.H"

#include "../drt_nurbs_discret/drt_control_point.H"
#include "../drt_beam2/beam2.H"
#include "../drt_beam2r/beam2r.H"
#include "../drt_beam3/beam3.H"
#include "../drt_beam3ii/beam3ii.H"
#include "../drt_smoothrod/smoothrod.H"
#include "../drt_truss3/truss3.H"
#include "../drt_trusslm/trusslm.H"
#include "../drt_truss2/truss2.H"
#include "../drt_torsion3/torsion3.H"
#include "../drt_torsion2/torsion2.H"
#include "../drt_s8/shell8.H"
#include "../drt_scatra/scatra_element.H"
#include "../drt_f3/fluid3.H"
#include "../drt_f3/xfluid3.H"
#include "../drt_xdiff3/xdiff3.H"
#include "../drt_combust/combust3.H"
#include "../drt_ale2/ale2.H"
#include "../drt_ale2/ale2_nurbs.H"
#include "../drt_ale3/ale3.H"
#include "../drt_ale3/ale3_nurbs.H"
#include "../drt_bele3/bele3_4.H"
#include "../drt_bele3/bele3.H"
#include "../drt_bele3/vele3.H"
#include "../drt_bele3/bele2.H"
#include "../drt_constraint/constraint_element2.H"
#include "../drt_constraint/constraint_element3.H"
#include "../drt_w1/wall1.H"
#include "../drt_w1/wall1_nurbs.H"
#include "../drt_so3/so_hex8.H"
#include "../drt_so3/so_hex20.H"
#include "../drt_so3/so_hex27.H"
#include "../drt_so3/so_nurbs27.H"
#include "../drt_so3/so_sh8.H"
#include "../drt_so3/so_sh8p8.H"
#include "../drt_so3/so_tet4.H"
#include "../drt_so3/so_ptet.H"
#include "../drt_so3/so_nstet.H"
#include "../drt_so3/so_tet10.H"
#include "../drt_so3/so_weg6.H"
#include "../drt_so3/so_shw6.H"
#include "../drt_so3/so_disp.H"
#include "../drt_so3/so_hex8p1j1.H"
#include "../drt_so3/so_hex8fbar.H"
#include "../drt_thermo/thermo_element.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/thermostvenantkirchhoff.H"
#include "../drt_mat/thermoplasticlinelast.H"
#include "../drt_mat/micromaterial.H"
#include "../drt_mat/neohooke.H"
#include "../drt_mat/logneohooke.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/aaaneohooke_stopro.H"
#include "../drt_mat/aaaraghavanvorp_damage.H"
#include "../drt_mat/aaa_mixedeffects.H"
#include "../drt_mat/aaagasser.H"
#include "../drt_mat/lung_penalty.H"
#include "../drt_mat/lung_ogden.H"
#include "../drt_mat/anisotropic_balzani.H"
#include "../drt_mat/mooneyrivlin.H"
#include "../drt_mat/yeoh.H"
#include "../drt_mat/visconeohooke.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/contchainnetw.H"
#include "../drt_mat/artwallremod.H"
#include "../drt_mat/scatra_mat.H"
#include "../drt_mat/ion.H"
#include "../drt_mat/mixfrac.H"
#include "../drt_mat/sutherland.H"
#include "../drt_mat/arrhenius_spec.H"
#include "../drt_mat/arrhenius_temp.H"
#include "../drt_mat/arrhenius_pv.H"
#include "../drt_mat/ferech_pv.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"
#include "../drt_mat/yoghurt.H"
#include "../drt_mat/biocell.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/charmm.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_mat/cnst_1d_art.H"
#include "../drt_mat/fourieriso.H"
#include "../drt_mat/holzapfelcardiovascular.H"
#include "../drt_mat/humphreycardiovascular.H"
#include "../drt_mat/growth_ip.H"
#include "../drt_mat/constraintmixture.H"
#include "../drt_mat/itskov.H"
#include "../drt_mat/plasticneohooke.H"
#include "../drt_mat/plasticlinelast.H"
#include "../drt_mat/robinson.H"
#include "../drt_mat/biofilm.H"
#include "../drt_mortar/mortar_node.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_contact/contact_node.H"
#include "../drt_contact/friction_node.H"
#include "../drt_contact/contact_element.H"
#include "../drt_art_net/artery.H"
#include "../drt_red_airways/red_airway.H"


std::string DRT::ParObjectList()
{
  std::stringstream s;

  s << DRT::ContainerType::Instance().Name() << " "
    << DRT::ConditionObjectType::Instance().Name() << " "
    << DRT::NodeType::Instance().Name() << " "
    << DRT::NURBS::ControlPointType::Instance().Name() << " "
#ifdef D_BEAM2
    << DRT::ELEMENTS::Beam2Type::Instance().Name() << " "
#endif
#ifdef D_BEAM2R
    << DRT::ELEMENTS::Beam2rType::Instance().Name() << " "
#endif
#ifdef D_BEAM3
    << DRT::ELEMENTS::Beam3Type::Instance().Name() << " "
#endif
#ifdef D_BEAM3II
    << DRT::ELEMENTS::Beam3iiType::Instance().Name() << " "
#endif
#ifdef D_SMOOTHROD
    << DRT::ELEMENTS::SmoothrodType::Instance().Name() << " "
#endif
#ifdef D_TRUSS3
    << DRT::ELEMENTS::Truss3Type::Instance().Name() << " "
    << DRT::ELEMENTS::TrussLmType::Instance().Name() << " "
#endif
#ifdef D_TRUSS2
    << DRT::ELEMENTS::Truss2Type::Instance().Name() << " "
#endif
#ifdef D_TORSION3
    << DRT::ELEMENTS::Torsion3Type::Instance().Name() << " "
#endif
#ifdef D_TORSION2
    << DRT::ELEMENTS::Torsion2Type::Instance().Name() << " "
#endif
#ifdef D_SHELL8
    << DRT::ELEMENTS::Shell8Type::Instance().Name() << " "
#endif
#ifdef D_WALL1
    << DRT::ELEMENTS::Wall1Type::Instance().Name() << " "
    << DRT::ELEMENTS::NURBS::Wall1NurbsType::Instance().Name() << " "
#endif
#ifdef D_FLUID3
    << DRT::ELEMENTS::Combust3Type::Instance().Name() << " "
    << DRT::ELEMENTS::Fluid3Type::Instance().Name() << " "
    << DRT::ELEMENTS::XDiff3Type::Instance().Name() << " "
    << DRT::ELEMENTS::XFluid3Type::Instance().Name() << " "
#endif
#ifdef D_ALE
    << DRT::ELEMENTS::Ale3Type::Instance().Name() << " "
    << DRT::ELEMENTS::NURBS::Ale3_NurbsType::Instance().Name() << " "
#endif
#ifdef D_ALE
    << DRT::ELEMENTS::Ale2Type::Instance().Name() << " "
    << DRT::ELEMENTS::NURBS::Ale2_NurbsType::Instance().Name() << " "
#endif
    << DRT::ELEMENTS::Bele2Type::Instance().Name() << " "
    << DRT::ELEMENTS::Bele3Type::Instance().Name() << " "
    << DRT::ELEMENTS::Bele3_4Type::Instance().Name() << " "
    << DRT::ELEMENTS::Vele3Type::Instance().Name() << " "
#ifdef D_SOLID3
    << DRT::ELEMENTS::NStetType::Instance().Name() << " "
    << DRT::ELEMENTS::NURBS::So_nurbs27Type::Instance().Name() << " "
    << DRT::ELEMENTS::PtetType::Instance().Name() << " "
    << DRT::ELEMENTS::SoDispType::Instance().Name() << " "
    << DRT::ELEMENTS::So_Hex8P1J1Type::Instance().Name() << " "
    << DRT::ELEMENTS::So_hex8fbarType::Instance().Name() << " "
    << DRT::ELEMENTS::So_hex20Type::Instance().Name() << " "
    << DRT::ELEMENTS::So_hex27Type::Instance().Name() << " "
    << DRT::ELEMENTS::So_hex8Type::Instance().Name() << " "
    << DRT::ELEMENTS::So_sh8Type::Instance().Name() << " "
    << DRT::ELEMENTS::So_sh8p8Type::Instance().Name() << " "
    << DRT::ELEMENTS::So_shw6Type::Instance().Name() << " "
    << DRT::ELEMENTS::So_tet10Type::Instance().Name() << " "
    << DRT::ELEMENTS::So_tet4Type::Instance().Name() << " "
    << DRT::ELEMENTS::So_weg6Type::Instance().Name() << " "
#endif
#ifdef D_ARTNET //_1D_ARTERY_
    << DRT::ELEMENTS::ArteryType::Instance().Name() << " "
    << MAT::Cnst_1d_artType::Instance().Name() << " "
#endif
    << MAT::AAAgasserType::Instance().Name() << " "
    << MAT::AAAneohookeType::Instance().Name() << " "
    << MAT::AAAneohooke_stoproType::Instance().Name() << " "
    << MAT::AAAraghavanvorp_damageType::Instance().Name() << " "
    << MAT::AAA_mixedeffectsType::Instance().Name() << " "
    << MAT::AnisotropicBalzaniType::Instance().Name() << " "
    << MAT::ArrheniusPVType::Instance().Name() << " "
    << MAT::ArrheniusSpecType::Instance().Name() << " "
    << MAT::ArrheniusTempType::Instance().Name() << " "
    << MAT::ArtWallRemodType::Instance().Name() << " "
    << MAT::BioCellType::Instance().Name() << " "
    << MAT::BiofilmType::Instance().Name() << " "
    << MAT::CHARMMType::Instance().Name() << " "
    << MAT::CarreauYasudaType::Instance().Name() << " "
    << MAT::ConstraintMixtureType::Instance().Name() << " "
    << MAT::ConstraintMixtureHistoryType::Instance().Name() << " "
    << MAT::ContChainNetwType::Instance().Name() << " "
    << MAT::ElastHyperType::Instance().Name() << " "
    << MAT::FerEchPVType::Instance().Name() << " "
    << MAT::FourierIsoType::Instance().Name() << " "
    << MAT::GrowthType::Instance().Name() << " "
    << MAT::HolzapfelCardioType::Instance().Name() << " "
    << MAT::HumphreyCardioType::Instance().Name() << " "
    << MAT::IonType::Instance().Name() << " "
    << MAT::ItskovType::Instance().Name() << " "
    << MAT::LogNeoHookeType::Instance().Name() << " "
    << MAT::LungOgdenType::Instance().Name() << " "
    << MAT::LungPenaltyType::Instance().Name() << " "
    << MAT::MatListType::Instance().Name() << " "
    << MAT::MicroMaterialType::Instance().Name() << " "
    << MAT::MixFracType::Instance().Name() << " "
    << MAT::ModPowerLawType::Instance().Name() << " "
    << MAT::MooneyRivlinType::Instance().Name() << " "
    << MAT::NeoHookeType::Instance().Name() << " "
    << MAT::NewtonianFluidType::Instance().Name() << " "
    << MAT::ScatraMatType::Instance().Name() << " "
    << MAT::StVenantKirchhoffType::Instance().Name() << " "
    << MAT::SutherlandType::Instance().Name() << " "
    << MAT::ThermoStVenantKirchhoffType::Instance().Name() << " "
    << MAT::ThermoPlasticLinElastType::Instance().Name() << " "
    << MAT::ViscoAnisotropicType::Instance().Name() << " "
    << MAT::ViscoNeoHookeType::Instance().Name() << " "
    << MAT::YeohType::Instance().Name() << " "
    << MAT::YoghurtType::Instance().Name() << " "
    << MORTAR::MortarNodeType::Instance().Name() << " "
    << MORTAR::MortarElementType::Instance().Name() << " "
    << CONTACT::CoNodeType::Instance().Name() << " "
    << CONTACT::FriNodeType::Instance().Name() << " "
    << CONTACT::CoElementType::Instance().Name() << " "
    << DRT::ELEMENTS::ConstraintElement2Type::Instance().Name() << " "
    << DRT::ELEMENTS::ConstraintElement3Type::Instance().Name() << " "
#if defined(D_FLUID3)
    << DRT::ELEMENTS::TransportType::Instance().Name() << " "
#endif
#ifdef D_THERMO
    << DRT::ELEMENTS::ThermoType::Instance().Name() << " "
#endif
    << MAT::PlasticNeoHookeType::Instance().Name() << " "
    << MAT::PlasticLinElastType::Instance().Name() << " "
    << MAT::RobinsonType::Instance().Name() << " "
#ifdef D_RED_AIRWAYS
    << DRT::ELEMENTS::RedAirwayType::Instance().Name() << " "
#endif
    ;
  return s.str();
}


void PrintParObjectList()
{
  std::cout << "defined parobject types: " << DRT::ParObjectList() << "\n";
}
