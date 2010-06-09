/*!----------------------------------------------------------------------
\file drt_utils.cpp
\brief A collection of helper methods for namespace DRT

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>
<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifdef PARALLEL
#ifdef LINUX_MUENCH
typedef int idxtype;
extern "C"
{
  void METIS_PartGraphRecursive(int *, idxtype *, idxtype *,
				idxtype *, idxtype *, int *, int *, int *,
				int *, int *, idxtype *);
  void METIS_PartGraphKway(int *, idxtype *, idxtype *, idxtype *,
			   idxtype *, int *, int *, int *, int *, int *,
			   idxtype *);
}
#endif
#include <mpi.h>
#endif // PARALLEL

#include <algorithm>
#include <numeric>
#include <vector>
#include <set>
#include <map>
#include <string>

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_mapextractor.H"

#include "drt_utils.H"
#include "drt_node.H"
#include "../drt_nurbs_discret/drt_control_point.H"
#include "drt_dofset.H"
#include "drt_discret.H"

#if 1
#include "../drt_beam2/beam2.H"
#include "../drt_beam2r/beam2r.H"
#include "../drt_beam3/beam3.H"
#include "../drt_beam3ii/beam3ii.H"
#include "../drt_smoothrod/smoothrod.H"
#include "../drt_truss3/truss3.H"
#include "../drt_truss2/truss2.H"
#include "../drt_torsion3/torsion3.H"
#include "../drt_torsion2/torsion2.H"
#include "../drt_s8/shell8.H"
#include "../drt_f2/fluid2.H"
#include "../drt_f2/fluid2_nurbs.H"
#include "../drt_scatra/scatra_element.H"
#include "../drt_f3/fluid3.H"
#include "../drt_f3/fluid3_nurbs.H"
#include "../drt_f3/xfluid3.H"
#include "../drt_xdiff3/xdiff3.H"
#include "../drt_combust/combust3.H"
#include "../drt_ale2/ale2.H"
#include "../drt_ale2/ale2_nurbs.H"
#include "../drt_ale3/ale3.H"
#include "../drt_ale3/ale3_nurbs.H"
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
#include "../drt_thermo/thermo_element.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/thermostvenantkirchhoff.H"
#include "../drt_mat/micromaterial.H"
#include "../drt_mat/neohooke.H"
#include "../drt_mat/logneohooke.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/aaaraghavanvorp_damage.H"
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
#include "../drt_mat/biocell.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/charmm.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_mat/cnst_1d_art.H"
#include "../drt_mat/fourieriso.H"
#include "../drt_mat/holzapfelcardiovascular.H"
#include "../drt_mat/humphreycardiovascular.H"
#include "../drt_mat/itskov.H"
#include "../drt_mat/plasticneohooke.H"
#include "../drt_mortar/mortar_node.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_contact/contact_node.H"
#include "../drt_contact/friction_node.H"
#include "../drt_contact/contact_element.H"
#include "../drt_art_net/artery.H"
#include "../drt_red_airways/red_airway.H"
#endif

#include "drt_dserror.H"
#include "standardtypes_cpp.H"

#include "drt_parobjectfactory.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

void ReferenceParObjectTypes()
{
  static int linker_hack;
  if ( linker_hack==0 )
  {
    linker_hack = 1;

    DRT::ContainerType::Instance().Name();
    DRT::ConditionObjectType::Instance().Name();
    DRT::NodeType::Instance().Name();
    DRT::NURBS::ControlPointType::Instance().Name();
#ifdef D_BEAM2
    DRT::ELEMENTS::Beam2Type::Instance().Name();
    DRT::ELEMENTS::Beam2RegisterType::Instance().Name();
#endif
#ifdef D_BEAM2R
    DRT::ELEMENTS::Beam2rType::Instance().Name();
    DRT::ELEMENTS::Beam2rRegisterType::Instance().Name();
#endif
#ifdef D_BEAM3
    DRT::ELEMENTS::Beam3Type::Instance().Name();
    DRT::ELEMENTS::Beam3RegisterType::Instance().Name();
#endif
#ifdef D_BEAM3II
    DRT::ELEMENTS::Beam3iiType::Instance().Name();
    DRT::ELEMENTS::Beam3iiRegisterType::Instance().Name();
#endif
#ifdef D_SMOOTHROD
    DRT::ELEMENTS::SmoothrodType::Instance().Name();
    DRT::ELEMENTS::SmoothrodRegisterType::Instance().Name();
#endif
#ifdef D_TRUSS3
    DRT::ELEMENTS::Truss3Type::Instance().Name();
    DRT::ELEMENTS::Truss3RegisterType::Instance().Name();
#endif
#ifdef D_TRUSS2
    DRT::ELEMENTS::Truss2Type::Instance().Name();
    DRT::ELEMENTS::Truss2RegisterType::Instance().Name();
#endif
#ifdef D_TORSION3
    DRT::ELEMENTS::Torsion3Type::Instance().Name();
    DRT::ELEMENTS::Torsion3RegisterType::Instance().Name();
#endif
#ifdef D_TORSION2
    DRT::ELEMENTS::Torsion2Type::Instance().Name();
    DRT::ELEMENTS::Torsion2RegisterType::Instance().Name();
#endif
#ifdef D_SHELL8
    DRT::ELEMENTS::Shell8Type::Instance().Name();
    DRT::ELEMENTS::Shell8RegisterType::Instance().Name();
#endif
#ifdef D_WALL1
    DRT::ELEMENTS::Wall1Type::Instance().Name();
    DRT::ELEMENTS::Wall1RegisterType::Instance().Name();
    DRT::ELEMENTS::NURBS::Wall1NurbsType::Instance().Name();
#endif
#ifdef D_FLUID2
    DRT::ELEMENTS::Fluid2Type::Instance().Name();
    DRT::ELEMENTS::NURBS::Fluid2NurbsType::Instance().Name();
    DRT::ELEMENTS::Fluid2RegisterType::Instance().Name();
#endif
#ifdef D_FLUID3
    DRT::ELEMENTS::Fluid3Type::Instance().Name();
    DRT::ELEMENTS::NURBS::Fluid3NurbsType::Instance().Name();
    DRT::ELEMENTS::Fluid3RegisterType::Instance().Name();
    DRT::ELEMENTS::XFluid3Type::Instance().Name();
    DRT::ELEMENTS::XFluid3RegisterType::Instance().Name();
    DRT::ELEMENTS::XDiff3Type::Instance().Name();
    DRT::ELEMENTS::XDiff3RegisterType::Instance().Name();
    DRT::ELEMENTS::Combust3Type::Instance().Name();
    DRT::ELEMENTS::Combust3RegisterType::Instance().Name();
#endif
#ifdef D_ALE
    DRT::ELEMENTS::Ale3Type::Instance().Name();
    DRT::ELEMENTS::NURBS::Ale3_NurbsType::Instance().Name();
#endif
#ifdef D_ALE
    DRT::ELEMENTS::Ale2Type::Instance().Name();
    DRT::ELEMENTS::NURBS::Ale2_NurbsType::Instance().Name();
#endif
    DRT::ELEMENTS::Bele3Type::Instance().Name();
    DRT::ELEMENTS::Bele3RegisterType::Instance().Name();
    DRT::ELEMENTS::Vele3Type::Instance().Name();
    DRT::ELEMENTS::Vele3RegisterType::Instance().Name();
    DRT::ELEMENTS::Bele2Type::Instance().Name();
    DRT::ELEMENTS::Bele2RegisterType::Instance().Name();
#ifdef D_SOLID3
    //DRT::ELEMENTS::So_hex8Type::Instance().Name();
    DRT::ELEMENTS::Soh8RegisterType::Instance().Name();
    DRT::ELEMENTS::So_sh8Type::Instance().Name();
    DRT::ELEMENTS::Sosh8RegisterType::Instance().Name();
    DRT::ELEMENTS::So_sh8p8Type::Instance().Name();
    DRT::ELEMENTS::Sosh8p8RegisterType::Instance().Name();
    DRT::ELEMENTS::So_hex27Type::Instance().Name();
    DRT::ELEMENTS::Soh27RegisterType::Instance().Name();
    DRT::ELEMENTS::NURBS::So_nurbs27Type::Instance().Name();
    DRT::ELEMENTS::NURBS::So_Nurbs27_RegisterType::Instance().Name();
    DRT::ELEMENTS::So_hex20Type::Instance().Name();
    DRT::ELEMENTS::Soh20RegisterType::Instance().Name();
    DRT::ELEMENTS::So_weg6Type::Instance().Name();
    DRT::ELEMENTS::Sow6RegisterType::Instance().Name();
    DRT::ELEMENTS::So_shw6Type::Instance().Name();
    DRT::ELEMENTS::Soshw6RegisterType::Instance().Name();
    DRT::ELEMENTS::SoDispType::Instance().Name();
    DRT::ELEMENTS::SoDispRegisterType::Instance().Name();
    DRT::ELEMENTS::So_tet10Type::Instance().Name();
    DRT::ELEMENTS::Sotet10RegisterType::Instance().Name();
    DRT::ELEMENTS::So_tet4Type::Instance().Name();
    DRT::ELEMENTS::Sotet4RegisterType::Instance().Name();
    DRT::ELEMENTS::PtetType::Instance().Name();
    DRT::ELEMENTS::PtetRegisterType::Instance().Name();
    DRT::ELEMENTS::NStetType::Instance().Name();
    DRT::ELEMENTS::NStetRegisterType::Instance().Name();
    DRT::ELEMENTS::So_Hex8P1J1Type::Instance().Name();
    DRT::ELEMENTS::SoHex8P1J1RegisterType::Instance().Name();
#endif
#ifdef D_ARTNET //_1D_ARTERY_
    DRT::ELEMENTS::ArteryType::Instance().Name();
    MAT::Cnst_1d_artType::Instance().Name();
#endif
    DRT::ElementRegisterType::Instance().Name();
    MAT::NewtonianFluidType::Instance().Name();
    MAT::StVenantKirchhoffType::Instance().Name();
    MAT::ThermoStVenantKirchhoffType::Instance().Name();
    MAT::LungPenaltyType::Instance().Name();
    MAT::LungOgdenType::Instance().Name();
    MAT::AnisotropicBalzaniType::Instance().Name();
    MAT::MooneyRivlinType::Instance().Name();
    MAT::YeohType::Instance().Name();
    MAT::ViscoNeoHookeType::Instance().Name();
    MAT::ViscoAnisotropicType::Instance().Name();
    MAT::ContChainNetwType::Instance().Name();
    MAT::ArtWallRemodType::Instance().Name();
    MAT::MicroMaterialType::Instance().Name();
    MAT::NeoHookeType::Instance().Name();
    MAT::LogNeoHookeType::Instance().Name();
    MAT::AAAneohookeType::Instance().Name();
    MAT::AAAraghavanvorp_damageType::Instance().Name();
    MAT::AAAgasserType::Instance().Name();
    MAT::ScatraMatType::Instance().Name();
    MAT::IonType::Instance().Name();
    MAT::MixFracType::Instance().Name();
    MAT::SutherlandType::Instance().Name();
    MAT::ArrheniusSpecType::Instance().Name();
    MAT::ArrheniusTempType::Instance().Name();
    MAT::ArrheniusPVType::Instance().Name();
    MAT::FerEchPVType::Instance().Name();
    MAT::CarreauYasudaType::Instance().Name();
    MAT::ModPowerLawType::Instance().Name();
    MAT::BioCellType::Instance().Name();
    MAT::CHARMMType::Instance().Name();
    MAT::ItskovType::Instance().Name();
    MAT::MatListType::Instance().Name();
    MAT::ElastHyperType::Instance().Name();
    MAT::FourierIsoType::Instance().Name();
    MAT::HolzapfelCardioType::Instance().Name();
    MAT::HumphreyCardioType::Instance().Name();
    MORTAR::MortarNodeType::Instance().Name();
    MORTAR::MortarElementType::Instance().Name();
    CONTACT::CoNodeType::Instance().Name();
    CONTACT::FriNodeType::Instance().Name();
    CONTACT::CoElementType::Instance().Name();
    DRT::ELEMENTS::ConstraintElement2Type::Instance().Name();
    DRT::ELEMENTS::ConstraintElement2RegisterType::Instance().Name();
    DRT::ELEMENTS::ConstraintElement3Type::Instance().Name();
    DRT::ELEMENTS::ConstraintElement3RegisterType::Instance().Name();
#if defined(D_FLUID2) || defined(D_FLUID3)
    DRT::ELEMENTS::TransportType::Instance().Name();
#endif
#ifdef D_THERMO
    DRT::ELEMENTS::ThermoType::Instance().Name();
#endif
    MAT::PlasticNeoHookeType::Instance().Name();
#ifdef D_RED_AIRWAYS
    DRT::ELEMENTS::RedAirwayType::Instance().Name();
    DRT::ELEMENTS::RedAirwayRegisterType::Instance().Name();
#endif
  }
}

/*----------------------------------------------------------------------*
 |  allocate an instance of a specific impl. of ParObject (public) mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::UTILS::Factory(const vector<char>& data)
{
#if 1

  ReferenceParObjectTypes();
  return ParObjectFactory::Instance().Create( data );

#else
  // mv ptr behind the size record
  const int* ptr = (const int*)(&data[0]);
  // get the type
  const int type = *ptr;
  switch(type)
  {
    case ParObject_Container:
    {
      DRT::Container* object = new DRT::Container();
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Condition:
    {
      DRT::Condition* object = new DRT::Condition();
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Node:
    {
      double dummycoord[3] = {999.,999.,999.};
      DRT::Node* object = new DRT::Node(-1,dummycoord,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_ControlPoint:
    {
      double dummycoord[3] = {999.,999.,999.};
      double dummyweight   =  999.;
      DRT::NURBS::ControlPoint* object = new DRT::NURBS::ControlPoint(-1,dummycoord,dummyweight,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Element:
    {
      dserror("DRT::Element is pure virtual, cannot create instance");
    }
    break;
#ifdef D_BEAM2
    case ParObject_Beam2:
    {
      DRT::ELEMENTS::Beam2* object = new DRT::ELEMENTS::Beam2(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Beam2Register:
    {
      DRT::ELEMENTS::Beam2Register* object =
                      new DRT::ELEMENTS::Beam2Register(DRT::Element::element_beam2);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_BEAM2R
    case ParObject_Beam2r:
    {
      DRT::ELEMENTS::Beam2r* object = new DRT::ELEMENTS::Beam2r(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Beam2rRegister:
    {
      DRT::ELEMENTS::Beam2rRegister* object =
                      new DRT::ELEMENTS::Beam2rRegister(DRT::Element::element_beam2r);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_BEAM3
    case ParObject_Beam3:
    {
      DRT::ELEMENTS::Beam3* object = new DRT::ELEMENTS::Beam3(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Beam3Register:
    {
      DRT::ELEMENTS::Beam3Register* object =
                      new DRT::ELEMENTS::Beam3Register(DRT::Element::element_beam3);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_BEAM3II
    case ParObject_Beam3ii:
    {
      DRT::ELEMENTS::Beam3ii* object = new DRT::ELEMENTS::Beam3ii(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Beam3iiRegister:
    {
      DRT::ELEMENTS::Beam3iiRegister* object =
                      new DRT::ELEMENTS::Beam3iiRegister(DRT::Element::element_beam3ii);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_SMOOTHROD
    case ParObject_Smoothrod:
    {
      DRT::ELEMENTS::Smoothrod* object = new DRT::ELEMENTS::Smoothrod(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_SmoothrodRegister:
    {
      DRT::ELEMENTS::SmoothrodRegister* object =
                      new DRT::ELEMENTS::SmoothrodRegister(DRT::Element::element_smoothrod);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_TRUSS3
    case ParObject_Truss3:
    {
      DRT::ELEMENTS::Truss3* object = new DRT::ELEMENTS::Truss3(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Truss3Register:
    {
      DRT::ELEMENTS::Truss3Register* object =
                      new DRT::ELEMENTS::Truss3Register(DRT::Element::element_truss3);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_TRUSS2
    case ParObject_Truss2:
    {
      DRT::ELEMENTS::Truss2* object = new DRT::ELEMENTS::Truss2(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Truss2Register:
    {
      DRT::ELEMENTS::Truss2Register* object =
                      new DRT::ELEMENTS::Truss2Register(DRT::Element::element_truss2);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_TORSION3
    case ParObject_Torsion3:
    {
      DRT::ELEMENTS::Torsion3* object = new DRT::ELEMENTS::Torsion3(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Torsion3Register:
    {
      DRT::ELEMENTS::Torsion3Register* object =
                      new DRT::ELEMENTS::Torsion3Register(DRT::Element::element_torsion3);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_TORSION2
    case ParObject_Torsion2:
    {
      DRT::ELEMENTS::Torsion2* object = new DRT::ELEMENTS::Torsion2(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Torsion2Register:
    {
      DRT::ELEMENTS::Torsion2Register* object =
                      new DRT::ELEMENTS::Torsion2Register(DRT::Element::element_torsion2);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_SHELL8
    case ParObject_Shell8:
    {
      DRT::ELEMENTS::Shell8* object = new DRT::ELEMENTS::Shell8(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Shell8Register:
    {
      DRT::ELEMENTS::Shell8Register* object =
                      new DRT::ELEMENTS::Shell8Register(DRT::Element::element_shell8);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_WALL1
    case ParObject_Wall1:
    {
      DRT::ELEMENTS::Wall1* object = new DRT::ELEMENTS::Wall1(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Wall1Register:
    {
      DRT::ELEMENTS::Wall1Register* object =
                      new DRT::ELEMENTS::Wall1Register(DRT::Element::element_wall1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Wall1Nurbs:
    {
      DRT::ELEMENTS::NURBS::Wall1Nurbs* object = new DRT::ELEMENTS::NURBS::Wall1Nurbs(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_FLUID2
    case ParObject_Fluid2:
    {
      DRT::ELEMENTS::Fluid2* object = new DRT::ELEMENTS::Fluid2(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Fluid2Nurbs:
    {
      DRT::ELEMENTS::NURBS::Fluid2Nurbs* object = new DRT::ELEMENTS::NURBS::Fluid2Nurbs(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Fluid2Register:
    {
      DRT::ELEMENTS::Fluid2Register* object =
                      new DRT::ELEMENTS::Fluid2Register(DRT::Element::element_fluid2);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_FLUID3
    case ParObject_Fluid3:
    {
      DRT::ELEMENTS::Fluid3* object = new DRT::ELEMENTS::Fluid3(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Fluid3Nurbs:
    {
      DRT::ELEMENTS::NURBS::Fluid3Nurbs* object = new DRT::ELEMENTS::NURBS::Fluid3Nurbs(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Fluid3Register:
    {
      DRT::ELEMENTS::Fluid3Register* object =
                      new DRT::ELEMENTS::Fluid3Register(DRT::Element::element_fluid3);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_XFluid3:
    {
      DRT::ELEMENTS::XFluid3* object = new DRT::ELEMENTS::XFluid3(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_XFluid3Register:
    {
      DRT::ELEMENTS::XFluid3Register* object =
                      new DRT::ELEMENTS::XFluid3Register(DRT::Element::element_xfluid3);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_XDiff3:
    {
      DRT::ELEMENTS::XDiff3* object = new DRT::ELEMENTS::XDiff3(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_XDiff3Register:
    {
      DRT::ELEMENTS::XDiff3Register* object =
                      new DRT::ELEMENTS::XDiff3Register(DRT::Element::element_xdiff3);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Combust3:
    {
      DRT::ELEMENTS::Combust3* object = new DRT::ELEMENTS::Combust3(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Combust3Register:
    {
      DRT::ELEMENTS::Combust3Register* object =
    	              new DRT::ELEMENTS::Combust3Register(DRT::Element::element_combust3);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_ALE
    case ParObject_Ale3:
    {
      DRT::ELEMENTS::Ale3* object = new DRT::ELEMENTS::Ale3(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Ale3Register:
    {
      DRT::ELEMENTS::Ale3Register* object =
                      new DRT::ELEMENTS::Ale3Register(DRT::Element::element_ale3);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Ale3_Nurbs:
    {
      DRT::ELEMENTS::NURBS::Ale3Nurbs* object = new DRT::ELEMENTS::NURBS::Ale3Nurbs(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_ALE
    case ParObject_Ale2:
    {
      DRT::ELEMENTS::Ale2* object = new DRT::ELEMENTS::Ale2(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Ale2Register:
    {
      DRT::ELEMENTS::Ale2Register* object =
                      new DRT::ELEMENTS::Ale2Register(DRT::Element::element_ale2);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Ale2_Nurbs:
    {
      DRT::ELEMENTS::NURBS::Ale2Nurbs* object = new DRT::ELEMENTS::NURBS::Ale2Nurbs(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
#endif
    case ParObject_Bele3:
    {
      DRT::ELEMENTS::Bele3* object = new DRT::ELEMENTS::Bele3(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Bele3Register:
    {
      DRT::ELEMENTS::Bele3Register* object =
                  new DRT::ELEMENTS::Bele3Register(DRT::Element::element_bele3);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Vele3:
    {
      DRT::ELEMENTS::Vele3* object = new DRT::ELEMENTS::Vele3(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Vele3Register:
    {
      DRT::ELEMENTS::Vele3Register* object =
                  new DRT::ELEMENTS::Vele3Register(DRT::Element::element_vele3);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Bele2:
    {
      DRT::ELEMENTS::Bele2* object = new DRT::ELEMENTS::Bele2(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Bele2Register:
    {
      DRT::ELEMENTS::Bele2Register* object =
                  new DRT::ELEMENTS::Bele2Register(DRT::Element::element_bele2);
      object->Unpack(data);
      return object;
    }
    break;
#ifdef D_SOLID3
    case ParObject_So_hex8:
    {
      DRT::ELEMENTS::So_hex8* object = new DRT::ELEMENTS::So_hex8(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Soh8Register:
    {
      DRT::ELEMENTS::Soh8Register* object =
                new DRT::ELEMENTS::Soh8Register(DRT::Element::element_so_hex8);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_So_sh8:
    {
      DRT::ELEMENTS::So_sh8* object = new DRT::ELEMENTS::So_sh8(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Sosh8Register:
    {
      DRT::ELEMENTS::Sosh8Register* object =
                new DRT::ELEMENTS::Sosh8Register(DRT::Element::element_sosh8);
      object->Unpack(data);
      return object;
    }
    case ParObject_So_sh8p8:
    {
      DRT::ELEMENTS::So_sh8p8* object = new DRT::ELEMENTS::So_sh8p8(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Sosh8p8Register:
    {
      DRT::ELEMENTS::Sosh8p8Register* object =
                new DRT::ELEMENTS::Sosh8p8Register(DRT::Element::element_sosh8p8);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_So_hex27:
    {
      DRT::ELEMENTS::So_hex27* object = new DRT::ELEMENTS::So_hex27(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Soh27Register:
    {
      DRT::ELEMENTS::Soh27Register* object =
                new DRT::ELEMENTS::Soh27Register(DRT::Element::element_so_hex27);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_So_Nurbs27:
    {
      DRT::ELEMENTS::NURBS::So_nurbs27* object = new DRT::ELEMENTS::NURBS::So_nurbs27(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_So_Nurbs27_Register:
    {
      DRT::ELEMENTS::NURBS::Sonurbs27Register* object =
                new DRT::ELEMENTS::NURBS::Sonurbs27Register(DRT::Element::element_so_nurbs27);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_So_hex20:
    {
      DRT::ELEMENTS::So_hex20* object = new DRT::ELEMENTS::So_hex20(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Soh20Register:
    {
      DRT::ELEMENTS::Soh20Register* object =
                new DRT::ELEMENTS::Soh20Register(DRT::Element::element_so_hex20);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_So_weg6:
    {
      DRT::ELEMENTS::So_weg6* object =
                new DRT::ELEMENTS::So_weg6(-1,-1);
      object->Unpack(data);
      return object;
    }
    case ParObject_Sow6Register:
    {
      DRT::ELEMENTS::Sow6Register* object =
                new DRT::ELEMENTS::Sow6Register(DRT::Element::element_so_weg6);
      object->Unpack(data);
      return object;
    }
    case ParObject_So_shw6:
    {
      DRT::ELEMENTS::So_shw6* object =
                new DRT::ELEMENTS::So_shw6(-1,-1);
      object->Unpack(data);
      return object;
    }
    case ParObject_Soshw6Register:
    {
      DRT::ELEMENTS::Soshw6Register* object =
                new DRT::ELEMENTS::Soshw6Register(DRT::Element::element_so_shw6);
      object->Unpack(data);
      return object;
    }
    case ParObject_SoDisp:
    {
      DRT::ELEMENTS::SoDisp* object =
                new DRT::ELEMENTS::SoDisp(-1,-1);
      object->Unpack(data);
      return object;
    }
    case ParObject_SoDispRegister:
    {
      DRT::ELEMENTS::SoDispRegister* object =
                new DRT::ELEMENTS::SoDispRegister(DRT::Element::element_sodisp);
      object->Unpack(data);
      return object;
    }
    case ParObject_So_tet10:
    {
      DRT::ELEMENTS::So_tet10* object = new DRT::ELEMENTS::So_tet10(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Sotet10Register:
    {
      DRT::ELEMENTS::Sotet10Register* object =
                new DRT::ELEMENTS::Sotet10Register(DRT::Element::element_so_tet10);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_So_tet4:
    {
      DRT::ELEMENTS::So_tet4* object = new DRT::ELEMENTS::So_tet4(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Sotet4Register:
    {
      DRT::ELEMENTS::Sotet4Register* object =
                new DRT::ELEMENTS::Sotet4Register(DRT::Element::element_so_tet4);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Ptet:
    {
      DRT::ELEMENTS::Ptet* object = new DRT::ELEMENTS::Ptet(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_PtetRegister:
    {
      DRT::ELEMENTS::PtetRegister* object =
                new DRT::ELEMENTS::PtetRegister(DRT::Element::element_ptet);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_NStet:
    {
      DRT::ELEMENTS::NStet* object = new DRT::ELEMENTS::NStet(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_NStetRegister:
    {
      DRT::ELEMENTS::NStetRegister* object =
                new DRT::ELEMENTS::NStetRegister(DRT::Element::element_nstet);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_So_Hex8P1J1:
    {
      DRT::ELEMENTS::So_Hex8P1J1* object = new DRT::ELEMENTS::So_Hex8P1J1(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_SoHex8P1J1Register:
    {
      DRT::ELEMENTS::SoHex8P1J1Register* object =
                new DRT::ELEMENTS::SoHex8P1J1Register(DRT::Element::element_so_hex8p1j1);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_ARTNET //_1D_ARTERY_
    case ParObject_Artery:
    {
      DRT::ELEMENTS::Artery* object = new DRT::ELEMENTS::Artery(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Cnst_1d_art:
    {
      MAT::Cnst_1d_art* cnst_art = new MAT::Cnst_1d_art();
      cnst_art->Unpack(data);
      return cnst_art;
    }
    break;
    case ParObject_ArteryRegister:
    {
      DRT::ELEMENTS::ArteryRegister* object = new DRT::ELEMENTS::ArteryRegister(DRT::Element::element_artery);
      object->Unpack(data);
      return object;
    }
    break;
#endif
    case ParObject_ElementRegister:
    {
      DRT::ElementRegister* object =
                new DRT::ElementRegister(DRT::Element::element_none);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_NewtonianFluid:
    {
      MAT::NewtonianFluid* fluid = new MAT::NewtonianFluid();
      fluid->Unpack(data);
      return fluid;
    }
    case ParObject_StVenantKirchhoff:
    {
      MAT::StVenantKirchhoff* stvenantk = new MAT::StVenantKirchhoff();
      stvenantk->Unpack(data);
      return stvenantk;
    }
    case ParObject_ThermoStVenantKirchhoff:
    {
      MAT::ThermoStVenantKirchhoff* thrstvenantk = new MAT::ThermoStVenantKirchhoff();
      thrstvenantk->Unpack(data);
      return thrstvenantk;
    }
    case ParObject_LungPenalty:
    {
      MAT::LungPenalty* lungpen = new MAT::LungPenalty();
      lungpen->Unpack(data);
      return lungpen;
    }
    case ParObject_LungOgden:
    {
      MAT::LungOgden* lungog = new MAT::LungOgden();
      lungog->Unpack(data);
      return lungog;
    }
    case ParObject_AnisotropicBalzani:
    {
      MAT::AnisotropicBalzani* anba = new MAT::AnisotropicBalzani();
      anba->Unpack(data);
      return anba;
    }
    case ParObject_MooneyRivlin:
    {
      MAT::MooneyRivlin* moon = new MAT::MooneyRivlin();
      moon->Unpack(data);
      return moon;
    }
    case ParObject_Yeoh:
    {
      MAT::Yeoh* yeoh = new MAT::Yeoh();
      yeoh->Unpack(data);
      return yeoh;
    }
    case ParObject_ViscoNeoHooke:
    {
      MAT::ViscoNeoHooke* visco = new MAT::ViscoNeoHooke();
      visco->Unpack(data);
      return visco;
    }
    case ParObject_ViscoAnisotropic:
    {
      MAT::ViscoAnisotropic* visco = new MAT::ViscoAnisotropic();
      visco->Unpack(data);
      return visco;
    }
    case ParObject_ContChainNetw:
    {
      MAT::ContChainNetw* chain = new MAT::ContChainNetw();
      chain->Unpack(data);
      return chain;
    }
    case ParObject_ArtWallRemod:
    {
      MAT::ArtWallRemod* remod = new MAT::ArtWallRemod();
      remod->Unpack(data);
      return remod;
    }
    case ParObject_MicroMaterial:
    {
      MAT::MicroMaterial* micro = new MAT::MicroMaterial();
      micro->Unpack(data);
      return micro;
    }
    case ParObject_NeoHooke:
    {
      MAT::NeoHooke* neo = new MAT::NeoHooke();
      neo->Unpack(data);
      return neo;
    }
    case ParObject_LogNeoHooke:
    {
      MAT::LogNeoHooke* logneo = new MAT::LogNeoHooke();
      logneo->Unpack(data);
      return logneo;
    }
    case ParObject_AAAneohooke:
    {
      MAT::AAAneohooke* aaa = new MAT::AAAneohooke();
      aaa->Unpack(data);
      return aaa;
    }
    case ParObject_AAAraghavanvorp_damage:
    {
      MAT::AAAraghavanvorp_damage* aaadamage = new MAT::AAAraghavanvorp_damage();
      aaadamage->Unpack(data);
      return aaadamage; //aaadam;
    }
    case ParObject_AAAgasser:
    {
      MAT::AAAgasser* aaa = new MAT::AAAgasser();
      aaa->Unpack(data);
      return aaa;
    }
    case ParObject_ScatraMat:
    {
      MAT::ScatraMat* scatra_mat = new MAT::ScatraMat();
      scatra_mat->Unpack(data);
      return scatra_mat;
    }
    case ParObject_Ion:
    {
      MAT::Ion* ion = new MAT::Ion();
      ion->Unpack(data);
      return ion;
    }
    case ParObject_MixFrac:
    {
      MAT::MixFrac* mixfrac = new MAT::MixFrac();
      mixfrac->Unpack(data);
      return mixfrac;
    }
    case ParObject_Sutherland:
    {
      MAT::Sutherland* sutherland = new MAT::Sutherland();
      sutherland->Unpack(data);
      return sutherland;
    }
    case ParObject_ArrheniusSpec:
    {
      MAT::ArrheniusSpec* arrhenius_spec = new MAT::ArrheniusSpec();
      arrhenius_spec->Unpack(data);
      return arrhenius_spec;
    }
    case ParObject_ArrheniusTemp:
    {
      MAT::ArrheniusTemp* arrhenius_temp = new MAT::ArrheniusTemp();
      arrhenius_temp->Unpack(data);
      return arrhenius_temp;
    }
    case ParObject_ArrheniusPV:
    {
      MAT::ArrheniusPV* arrhenius_pv = new MAT::ArrheniusPV();
      arrhenius_pv->Unpack(data);
      return arrhenius_pv;
    }
    case ParObject_FerEchPV:
    {
      MAT::FerEchPV* ferech_pv = new MAT::FerEchPV();
      ferech_pv->Unpack(data);
      return ferech_pv;
    }
    case ParObject_CarreauYasuda:
    {
      MAT::CarreauYasuda* carYas = new MAT::CarreauYasuda();
      carYas->Unpack(data);
      return carYas;
    }
    case ParObject_ModPowerLaw:
    {
      MAT::ModPowerLaw* powLaw = new MAT::ModPowerLaw();
      powLaw->Unpack(data);
      return powLaw;
    }
    case ParObject_BioCell:
    {
      MAT::BioCell* biocell = new MAT::BioCell();
      biocell->Unpack(data);
      return biocell;
    }
    case ParObject_CHARMM:
    {
      MAT::CHARMM* charmm = new MAT::CHARMM();
      charmm->Unpack(data);
      return charmm;
    }
    case ParObject_Itskov:
    {
      MAT::Itskov* its = new MAT::Itskov();
      its->Unpack(data);
      return its;
    }
    case ParObject_MatList:
    {
      MAT::MatList* matlist = new MAT::MatList();
      matlist->Unpack(data);
      return matlist;
    }
    case ParObject_ElastHyper:
    {
      MAT::ElastHyper* elhy = new MAT::ElastHyper();
      elhy->Unpack(data);
      return elhy;
    }
    case ParObject_FourierIso:
    {
      MAT::FourierIso* fourieriso = new MAT::FourierIso();
      fourieriso->Unpack(data);
      return fourieriso;
    }
    case ParObject_HolzapfelCardio:
    {
      MAT::HolzapfelCardio* holzapfelcard = new MAT::HolzapfelCardio();
      holzapfelcard->Unpack(data);
      return holzapfelcard;
    }
    case ParObject_HumphreyCardio:
    {
      MAT::HumphreyCardio* humcard = new MAT::HumphreyCardio();
      humcard->Unpack(data);
      return humcard;
    }
    case ParObject_MortarNode:
    {
      double x[3];
      vector<int> dofs(0);
      MORTAR::MortarNode* node = new MORTAR::MortarNode(0,x,0,0,dofs,false);
      node->Unpack(data);
      return node;
    }
    case ParObject_MortarElement:
    {
      MORTAR::MortarElement* ele = new MORTAR::MortarElement(0,
                                                     DRT::Element::element_mortar,
                                                     0,DRT::Element::dis_none,
                                                     0,NULL,false);
      ele->Unpack(data);
      return ele;
    }
    case ParObject_CoNode:
    {
      double x[3];
      vector<int> dofs(0);
      CONTACT::CoNode* node = new CONTACT::CoNode(0,x,0,0,dofs,false,false);
      node->Unpack(data);
      return node;
    }
    case ParObject_FriNode:
    {
      double x[3];
      vector<int> dofs(0);
      CONTACT::FriNode* node = new CONTACT::FriNode(0,x,0,0,dofs,false,false);
      node->Unpack(data);
      return node;
    }
    case ParObject_CoElement:
    {
      CONTACT::CoElement* ele = new CONTACT::CoElement(0,
                                                       DRT::Element::element_contact,
                                                       0,DRT::Element::dis_none,
                                                       0,NULL,false);
      ele->Unpack(data);
      return ele;
    }
    case ParObject_ConstraintElement2:
    {
      DRT::ELEMENTS::ConstraintElement2* object = new DRT::ELEMENTS::ConstraintElement2(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_ConstraintElement2Register:
    {
      DRT::ELEMENTS::ConstraintElement2Register* object =
                      new DRT::ELEMENTS::ConstraintElement2Register(DRT::Element::element_constraintelement2);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_ConstraintElement3:
    {
      DRT::ELEMENTS::ConstraintElement3* object = new DRT::ELEMENTS::ConstraintElement3(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_ConstraintElement3Register:
    {
      DRT::ELEMENTS::ConstraintElement3Register* object =
                      new DRT::ELEMENTS::ConstraintElement3Register(DRT::Element::element_constraintelement3);
      object->Unpack(data);
      return object;
    }
    break;
#if defined(D_FLUID2) || defined(D_FLUID3)
    case ParObject_Transport:
    {
      DRT::ELEMENTS::Transport* object =
                      new DRT::ELEMENTS::Transport(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_THERMO
    case ParObject_Thermo:
    {
      DRT::ELEMENTS::Thermo* object =
                      new DRT::ELEMENTS::Thermo(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
#endif
    case ParObject_PlasticNeoHooke:
    {
      MAT::PlasticNeoHooke* plastic = new MAT::PlasticNeoHooke();
      plastic->Unpack(data);
      return plastic;
    }
#ifdef D_RED_AIRWAYS
    case ParObject_RedAirway:
    {
      DRT::ELEMENTS::RedAirway* object = new DRT::ELEMENTS::RedAirway(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_RedAirwayRegister:
    {
      DRT::ELEMENTS::RedAirwayRegister* object = new DRT::ELEMENTS::RedAirwayRegister(DRT::Element::element_red_airway);
      object->Unpack(data);
      return object;
    }
    break;
#endif //D_RED_AIRWAY
    default:
      dserror("Unknown type of ParObject instance: %d",type);
    break;
  }

  return NULL;
#endif
}

/*----------------------------------------------------------------------*
 |  allocate an element of a specific type (public)          mwgee 03|07|
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::UTILS::Factory(const string eletype,
                                              const string eledistype,
                                              const int id,
                                              const int owner)
{
#if 1
  ReferenceParObjectTypes();
  return ParObjectFactory::Instance().Create( eletype, eledistype, id, owner );
#else
  enum TypeofElement
  {
    none,
    shell8,
    wall1,
    fluid2,
    fluid3,
    xfluid3,
    xdiff3,
    combust3,
    ale2,
    ale3,
    bele2,
    bele3,
    vele3,
    so_hex8,
    so_sh8,
    so_sh8p8,
    so_nurbs27,
    so_hex27,
    so_hex20,
    so_tet4,
    ptet,
    nstet,
    so_tet10,
    so_weg6,
    so_shw6,
    sodisp,
    beam2,
    beam2r,
    beam3,
    beam3ii,
    truss3,
    truss2,
    torsion3,
    torsion2,
    constrele2,
    constrele3,
    transport,
    art_ele,
    so_hex8p1j1,
    thermo,
    red_airway_ele
  };

  TypeofElement type = none;
  if (eletype=="none"); // dummy
  else if (eletype=="SHELL8") type = shell8;
  else if (eletype=="WALL")  type = wall1;
  else if (eletype=="FLUID2") type = fluid2;
  else if (eletype=="CONDIF2") type = transport; //backward compatibility
  else if (eletype=="CONDIF3") type = transport; //backward compatibility
  else if (eletype=="FLUID3") type = fluid3;
  else if (eletype=="XFLUID3") type = xfluid3;
  else if (eletype=="XDIFF3") type = xdiff3;
  else if (eletype=="COMBUST3") type = combust3;
  else if (eletype=="ALE2") type = ale2;
  else if (eletype=="ALE3") type = ale3;
  else if (eletype=="BELE2") type = bele2;
  else if (eletype=="BELE3") type = bele3;
  else if (eletype=="VELE3") type = vele3;
  else if (eletype=="SOLIDH8") type = so_hex8;
  else if (eletype=="SOLIDSH8") type = so_sh8;
  else if (eletype=="SOLIDSH8P8") type = so_sh8p8;
  else if (eletype=="SONURBS27") type = so_nurbs27;
  else if (eletype=="SOLIDH27") type = so_hex27;
  else if (eletype=="SOLIDH20") type = so_hex20;
  else if (eletype=="SOLIDT4") type = so_tet4;
  else if (eletype=="PTET4") type = ptet;
  else if (eletype=="NSTET4") type = nstet;
  else if (eletype=="SOLIDT10") type = so_tet10;
  else if (eletype=="SOLIDW6") type = so_weg6;
  else if (eletype=="SOLIDSHW6") type = so_shw6;
  else if (eletype=="SOLID3") type = sodisp;
  else if (eletype=="BEAM2") type = beam2;
  else if (eletype=="BEAM2R") type = beam2r;
  else if (eletype=="BEAM3") type = beam3;
  else if (eletype=="BEAM3II") type = beam3ii;
  else if (eletype=="TRUSS3") type = truss3;
  else if (eletype=="TRUSS2") type = truss2;
  else if (eletype=="TORSION3") type = torsion3;
  else if (eletype=="TORSION2") type = torsion2;
  else if (eletype=="CONSTRELE2") type = constrele2;
  else if (eletype=="CONSTRELE3") type = constrele3;
  else if (eletype=="TRANSP") type = transport;
  else if (eletype=="ART")    type = art_ele;
  else if (eletype=="SOLIDH8P1J1") type = so_hex8p1j1;
  else if (eletype=="THERMO") type = thermo;
  else if (eletype=="RED_AIRWAY")    type = red_airway_ele;
  // continue to add elements here....
  else
  {
    cout << endl << eletype << endl;
    dserror("Unknown type of finite element");
  }


  switch (type)
  {
#ifdef D_BEAM2
    case beam2:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Beam2(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_BEAM2R
    case beam2r:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Beam2r(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_BEAM3
    case beam3:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Beam3(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_BEAM3II
    case beam3ii:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Beam3ii(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_TRUSS3
    case truss3:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Truss3(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_TRUSS2
    case truss2:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Truss2(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_TORSION3
    case torsion3:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Torsion3(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_TORSION2
    case torsion2:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Torsion2(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_SHELL8
    case shell8:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Shell8(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_FLUID2
    case fluid2:
    {
      Teuchos::RCP<DRT::Element> ele;

      if(eledistype=="NURBS4" || eledistype=="NURBS9")
      {
        ele = rcp(new DRT::ELEMENTS::NURBS::Fluid2Nurbs(id,owner));
      }
      else
      {
        ele = rcp(new DRT::ELEMENTS::Fluid2(id,owner));
      }
      return ele;
    }
    break;
#endif
#ifdef D_FLUID3
    case fluid3:
    {
      Teuchos::RCP<DRT::Element> ele;

      if(eledistype=="NURBS8" || eledistype=="NURBS27")
      {
        ele = rcp(new DRT::ELEMENTS::NURBS::Fluid3Nurbs(id,owner));
      }
      else
      {
        ele = rcp(new DRT::ELEMENTS::Fluid3(id,owner));
      }
      return ele;
    }
    break;
    case xfluid3:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::XFluid3(id,owner));
      return ele;
    }
    break;
    case xdiff3:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::XDiff3(id,owner));
      return ele;
    }
    break;
    case combust3:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Combust3(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_ALE
    case ale3:
    {
      Teuchos::RCP<DRT::Element> ele;

      if(eledistype=="NURBS27")
      {
        ele = rcp(new DRT::ELEMENTS::NURBS::Ale3Nurbs(id,owner));
      }
      else
      {
        ele = rcp(new DRT::ELEMENTS::Ale3(id,owner));
      }

      return ele;
    }
    break;
#endif
#ifdef D_ALE
    case ale2:
    {
      Teuchos::RCP<DRT::Element> ele;

      if(eledistype=="NURBS4" || eledistype=="NURBS9")
      {
        ele = rcp(new DRT::ELEMENTS::NURBS::Ale2Nurbs(id,owner));
      }
      else
      {
        ele = rcp(new DRT::ELEMENTS::Ale2(id,owner));
      }

      return ele;
    }
    break;
#endif
    case bele3:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Bele3(id,owner));
      return ele;
    }
    break;
    case vele3:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Vele3(id,owner));
      return ele;
    }
    break;
    case bele2:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Bele2(id,owner));
      return ele;
    }
    break;
    case constrele2:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::ConstraintElement2(id,owner));
      return ele;
    }
    break;
    case constrele3:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::ConstraintElement3(id,owner));
      return ele;
    }
    break;
#ifdef D_WALL1
    case wall1:
    {
      Teuchos::RCP<DRT::Element> ele;
      if(eledistype=="NURBS4" || eledistype=="NURBS9")
      {
        ele = rcp(new DRT::ELEMENTS::NURBS::Wall1Nurbs(id,owner));
      }
      else
      {
        ele = rcp(new DRT::ELEMENTS::Wall1(id,owner));
      }
      return ele;
    }
    break;
#endif
#ifdef D_SOLID3
    case so_hex8:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_hex8(id,owner));
      return ele;
    }
    break;
    case so_sh8:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_sh8(id,owner));
      return ele;
    }
    break;
    case so_sh8p8:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_sh8p8(id,owner));
      return ele;
    }
    break;
    case so_hex27:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_hex27(id,owner));
      return ele;
    }
    break;
    case so_nurbs27:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::NURBS::So_nurbs27(id,owner));
      return ele;
    }
    break;
    case so_hex20:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_hex20(id,owner));
      return ele;
    }
    break;
    case so_weg6:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_weg6(id,owner));
      return ele;
    }
    break;
    case so_shw6:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_shw6(id,owner));
      return ele;
    }
    break;
    case sodisp:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::SoDisp(id,owner));
      return ele;
    }
    break;
    case so_tet4:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_tet4(id,owner));
      return ele;
    }
    break;
    case ptet:
    {
      RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Ptet(id,owner));
      return ele;
    }
    break;
    case nstet:
    {
      RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::NStet(id,owner));
      return ele;
    }
    break;
    case so_tet10:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_tet10(id,owner));
      return ele;
    }
    break;
    case so_hex8p1j1:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_Hex8P1J1(id,owner));
      return ele;
    }
    break;
#endif
#if defined(D_FLUID2) || defined(D_FLUID3)
    case transport:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Transport(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_ARTNET
    case art_ele:
    {
       Teuchos::RCP<DRT::Element> ele =  rcp(new DRT::ELEMENTS::Artery(id,owner));
       return ele;
    }
    break;
#endif
#ifdef D_THERMO
    case thermo:
    {
      Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Thermo(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_RED_AIRWAYS
    case red_airway_ele:
    {
       Teuchos::RCP<DRT::Element> ele =  rcp(new DRT::ELEMENTS::RedAirway(id,owner));
       return ele;
    }
    break;
#endif
    // continue to add types of elements here....
    default:
      dserror("Unknown type '%s' of finite element", eletype.c_str());
    break;
  }

  return null;
#endif
}


/*----------------------------------------------------------------------*
 |  partition a graph using metis  (public)                  mwgee 11/06|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsGraph> DRT::UTILS::PartGraphUsingMetis(
                                             const Epetra_CrsGraph& graph,
                                             const Epetra_Vector& weights)
{
  const int myrank   = graph.Comm().MyPID();
  const int numproc  = graph.Comm().NumProc();

  if (numproc==1)
  {
    Teuchos::RCP<Epetra_CrsGraph> outgraph = rcp(new Epetra_CrsGraph(graph));
    return outgraph;
  }

  // proc that will do the serial partitioning
  // the graph is collapsed to this proc
  // Normally this would be proc 0 but 0 always has so much to do.... ;-)
  int workrank = 0;
  if (graph.Comm().NumProc()>1) workrank=1;

  // get rowmap of the graph
  const Epetra_BlockMap& tmp = graph.RowMap();
  Epetra_Map rowmap(tmp.NumGlobalElements(),tmp.NumMyElements(),
                    tmp.MyGlobalElements(),0,graph.Comm());

  // build a target map that stores everything on proc workrank
  // We have arbirtary gids here and we do not tell metis about
  // them. So we have to keep rowrecv until the redistributed map is
  // build.


  // rowrecv is a fully redundant vector (size of number of nodes)
  vector<int> rowrecv(rowmap.NumGlobalElements());

  // after Allreduce rowrecv contains
  //
  // *-+-+-    -+-+-*-+-+-    -+-+-*-           -*-+-+-    -+-+-*
  // * | | .... | | * | | .... | | * ..........  * | | .... | | *
  // *-+-+-    -+-+-*-+-+-    -+-+-*-           -*-+-+-    -+-+-*
  //   gids stored     gids stored                  gids stored
  //  on first proc  on second proc                 on last proc
  //
  // the ordering of the gids on the procs is arbitrary (as copied
  // from the map)
  LINALG::AllreduceEMap(rowrecv, rowmap);

  // construct an epetra map from the list of gids
  Epetra_Map tmap(rowmap.NumGlobalElements(),
                  // if ..........    then ............... else
                  (myrank == workrank) ? (int)rowrecv.size() : 0,
                  &rowrecv[0],
                  0,
                  rowmap.Comm());

  // export the graph to tmap
  Epetra_CrsGraph tgraph(Copy,tmap,108,false);
  Epetra_Export exporter(rowmap,tmap);
  int err = tgraph.Export(graph,exporter,Add);
  if (err<0) dserror("Graph export returned err=%d",err);
  tgraph.FillComplete();
  tgraph.OptimizeStorage();

  // export the weights to tmap
  Epetra_Vector tweights(tmap,false);
  err = tweights.Export(weights,exporter,Insert);
  if (err<0) dserror("Vector export returned err=%d",err);

  // do partitioning using metis on workrank
  vector<int> part(tmap.NumMyElements());
  if (myrank==workrank)
  {
    // metis requests indexes. So we need a reverse lookup from gids
    // to indexes.
    map<int,int> idxmap;
    for (unsigned i=0; i<rowrecv.size(); ++i)
    {
      idxmap[rowrecv[i]] = i;
    }

    // xadj points from index i to the index of the
    // first adjacent node
    vector<int> xadj(tmap.NumMyElements()+1);

    // a list of adjacent nodes, adressed using xadj
    vector<int> adjncy(tgraph.NumGlobalNonzeros()); // the size is an upper bound

    // rowrecv(i)       rowrecv(i+1)                      node gids
    //     ^                 ^
    //     |                 |
    //     | idxmap          | idxmap
    //     |                 |
    //     v                 v
    //     i                i+1                       equivalent indices
    //     |                 |
    //     | xadj            | xadj
    //     |                 |
    //     v                 v
    //    +-+-+-+-+-+-+-+-+-+-+                -+-+-+
    //    | | | | | | | | | | | ............... | | |      adjncy
    //    +-+-+-+-+-+-+-+-+-+-+                -+-+-+
    //
    //    |       i's       |    (i+1)'s
    //    |    neighbours   |   neighbours           (numbered by equivalent indices)
    //

    vector<int> vwgt(tweights.MyLength());
    for (int i=0; i<tweights.MyLength(); ++i) vwgt[i] = (int)tweights[i];

    int count=0;
    xadj[0] = 0;
    for (int row=0; row<tgraph.NumMyRows(); ++row)
    {
      //cout << "xadj[" << row << "] = " << xadj[row] << endl;
      int grid = tgraph.RowMap().GID(row);
      int numindices;
      int* lindices;
      int err = tgraph.ExtractMyRowView(row,numindices,lindices);
      if (err) dserror("Epetra_CrsGraph::ExtractMyRowView returned err=%d",err);
      //cout << "adjncy: ";
      for (int col=0; col<numindices; ++col)
      {
        int gcid = tgraph.ColMap().GID(lindices[col]);
        if (gcid==grid) continue;
        adjncy[count] = idxmap[gcid];
        //cout << adjncy[count] << " ";
        ++count;
      }
      //cout << endl;
      xadj[row+1] = count;
    }
    //cout << "xadj[" << xadj.size()-1 << "] = " << xadj[xadj.size()-1] << endl;
    //cout << "tgraph.NumGlobalNonzeros() " << tgraph.NumGlobalNonzeros() << endl
    //     << "tmap.NumMyElements()       " << tmap.NumMyElements() << endl
    //     << "count                      " << count << endl;

    idxmap.clear();

    if (numproc<8) // better for smaller no. of partitions
    {
#ifdef PARALLEL
      int wgtflag=2;
      int numflag=0;
      int npart=numproc;
      int options[5] = { 0,3,1,1,0 };
      int edgecut=0;
      int nummyele = tmap.NumMyElements();
      METIS_PartGraphRecursive(&nummyele,
                               &xadj[0],
                               &adjncy[0],
                               &vwgt[0],
                               NULL,
                               &wgtflag,
                               &numflag,
                               &npart,
                               options,
                               &edgecut,
                               &part[0]);
#endif
    }
    else
    {
#ifdef PARALLEL
      int wgtflag=2;
      int numflag=0;
      int npart=numproc;
      int options[5] = { 0,3,1,1,0 };
      int edgecut=0;
      int nummyele = tmap.NumMyElements();
      METIS_PartGraphKway(&nummyele,
                          &xadj[0],
                          &adjncy[0],
                          &vwgt[0],
                          NULL,
                          &wgtflag,
                          &numflag,
                          &npart,
                          options,
                          &edgecut,
                          &part[0]);
#endif
    }
  } // if (myrank==workrank)

  // broadcast partitioning result
  int size = tmap.NumMyElements();
  tmap.Comm().Broadcast(&size,1,workrank);
  part.resize(size);
  tmap.Comm().Broadcast(&part[0],size,workrank);

  // loop part and count no. of nodes belonging to me
  // (we reuse part to save on memory)
  int count=0;
  for (int i=0; i<size; ++i)
    if (part[i]==myrank)
    {
      part[count] = rowrecv[i];
      ++count;
    }

  rowrecv.clear();

  // create map with new layout
  Epetra_Map newmap(size,count,&part[0],0,graph.Comm());

  // create the output graph and export to it
  Teuchos::RCP<Epetra_CrsGraph> outgraph =
                           rcp(new Epetra_CrsGraph(Copy,newmap,108,false));
  Epetra_Export exporter2(graph.RowMap(),newmap);
  err = outgraph->Export(graph,exporter2,Add);
  if (err<0) dserror("Graph export returned err=%d",err);
  outgraph->FillComplete();
  outgraph->OptimizeStorage();

  return outgraph;
}



/*----------------------------------------------------------------------*
 |  locally extract a subset of values  (public)            mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyValues(const Epetra_Vector& global,
                                 vector<double>& local,
                                 const vector<int>& lm)
{
  const size_t ldim = lm.size();
  local.resize(ldim);
  for (size_t i=0; i<ldim; ++i)
  {
    const int lid = global.Map().LID(lm[i]);
    if (lid<0)
      dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",global.Comm().MyPID(),lm[i]);
    local[i] = global[lid];
  }
  return;
}


/*----------------------------------------------------------------------*
 |  locally extract a subset of values  (public)             henke 12/09|
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyValues(const Epetra_Vector&      global,
                                 Epetra_SerialDenseVector& local,
                                 const vector<int>&        lm)
{
  const size_t ldim = lm.size();
  local.Size(ldim);
  for (size_t i=0; i<ldim; ++i)
  {
    const int lid = global.Map().LID(lm[i]);
    if (lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",global.Comm().MyPID(),lm[i]);
    local[i] = global[lid];
  }
  return;
}


/*----------------------------------------------------------------------*
 | extract local values from global node-based (multi) vector           |
 |                                                          henke 06/09 |
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyNodeBasedValues(
    const DRT::Element* ele,
    std::vector<double>& local,
    const Epetra_MultiVector& global)
{
  const int numnode = ele->NumNode();
  const int numcol = global.NumVectors();
  local.resize(numnode*numcol);

  // loop over element nodes
  for (int i=0; i<numnode; ++i)
  {
    const int nodegid = (ele->Nodes()[i])->Id();
    const int lid = global.Map().LID(nodegid);
    if (lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",global.Comm().MyPID(),nodegid);

    // loop over multi vector columns (numcol=1 for Epetra_Vector)
    for (int col=0; col<numcol; col++)
    {
      double* globalcolumn = (global)[col];
      local[col+(numcol*i)] = globalcolumn[lid];
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | extract local values from global node-based multi vector   gjb 08/08 |
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyNodeBasedValues(
    const DRT::Element* ele,
    Epetra_SerialDenseVector& local,
    const RCP<Epetra_MultiVector>& global,
    const int nsd
    )
{
  if (global==null) dserror("received a TEUCHOS::null pointer");
  if (nsd > global->NumVectors())
    dserror("Requested %d of %d available columns", nsd,global->NumVectors());
  const int iel = ele->NumNode(); // number of nodes
  if (local.Length()!=(iel*nsd)) dserror("vector size mismatch.");

  for (int i=0; i<nsd; i++)
  {
    // access actual component column of multi-vector
    double* globalcolumn = (*global)[i];
    // loop over the element nodes
    for (int j=0;j<iel;j++)
    {
      const int nodegid = (ele->Nodes()[j])->Id();
      const int lid = global->Map().LID(nodegid);
      local(i+(nsd*j))=globalcolumn[lid];
    }
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::RedistributeWithNewNodalDistribution(
    DRT::Discretization&     dis,
    const Epetra_Map&        noderowmap,
    const Epetra_Map&        nodecolmap
    )
{
  // redistribute nodes to column (ghost) map
  dis.ExportColumnNodes(nodecolmap);

  Teuchos::RCP< Epetra_Map > elerowmap;
  Teuchos::RCP< Epetra_Map > elecolmap;

  // now we have all elements in a linear map roweles
  // build resonable maps for elements from the
  // already valid and final node maps
  dis.BuildElementRowColumn(noderowmap, nodecolmap, elerowmap, elecolmap);

  // we can now export elements to resonable row element distribution
  dis.ExportRowElements(*elerowmap);

  // export to the column map / create ghosting of elements
  dis.ExportColumnElements(*elecolmap);

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::PrintParallelDistribution(const DRT::Discretization& dis)
{
  const int numproc=dis.Comm().NumProc();

  if(numproc>1)
  {
    const int myrank=dis.Comm().MyPID();

    vector<int> my_n_nodes     (numproc,0);
    vector<int>    n_nodes     (numproc,0);
    vector<int> my_n_ghostnodes(numproc,0);
    vector<int>    n_ghostnodes(numproc,0);
    vector<int> my_n_elements  (numproc,0);
    vector<int>    n_elements  (numproc,0);
    vector<int> my_n_ghostele  (numproc,0);
    vector<int>    n_ghostele  (numproc,0);

    my_n_nodes     [myrank]=dis.NumMyRowNodes();
    my_n_ghostnodes[myrank]=dis.NumMyColNodes()-my_n_nodes[myrank];
    my_n_elements  [myrank]=dis.NumMyRowElements();
    my_n_ghostele  [myrank]=dis.NumMyColElements()-my_n_elements[myrank];

    dis.Comm().SumAll(&my_n_nodes     [0],&n_nodes     [0],numproc);
    dis.Comm().SumAll(&my_n_ghostnodes[0],&n_ghostnodes[0],numproc);
    dis.Comm().SumAll(&my_n_elements  [0],&n_elements  [0],numproc);
    dis.Comm().SumAll(&my_n_ghostele  [0],&n_ghostele  [0],numproc);

    if(myrank==0)
    {
      cout << endl;
      cout <<"   Discretization: " << dis.Name() << endl;
      printf("   +-----+---------------+--------------+-----------------+----------------+\n");
      printf("   | PID |  n_rownodes   | n_ghostnodes |  n_rowelements  |   n_ghostele   |\n");
      printf("   +-----+---------------+--------------+-----------------+----------------+\n");
      for(int npid=0;npid<numproc;++npid)
      {
        printf("   | %3d | %13d | %12d | %15d | %14d |\n",npid,n_nodes[npid],n_ghostnodes[npid],n_elements[npid],n_ghostele[npid]);
        printf("   +-----+---------------+--------------+-----------------+----------------+\n");
      }
      cout << endl;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> DRT::UTILS::GetColVersionOfRowVector(
    const Teuchos::RCP<const DRT::Discretization> dis,
    const Teuchos::RCP<const Epetra_Vector> state)
{
  // note that this routine has the same functionality as SetState,
  // although here we do not store the new vector anywhere
  // maybe this routine can be used in SetState or become a member function of the discretization class

  if (!dis->HaveDofs()) dserror("FillComplete() was not called");
  const Epetra_Map* colmap = dis->DofColMap();
  const Epetra_BlockMap& vecmap = state->Map();

  // if it's already in column map just set a reference
  // This is a rought test, but it might be ok at this place. It is an
  // error anyway to hand in a vector that is not related to our dof
  // maps.
  if (vecmap.PointSameAs(*colmap))
    return state;
  // if it's not in column map export and allocate
  else
  {
    Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*colmap,false);
    LINALG::Export(*state,*tmp);
    return tmp;
  }
}

/*----------------------------------------------------------------------*
 |(private)                                                   tk 06/10  |
 |recompute nodecolmap of standard discretization to include all        |
 |nodes as of subdicretization                                          |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> DRT::UTILS::ComputeNodeColMap(
         const RCP<DRT::Discretization> sourcedis,  ///< standard discretization we want to redistribute
         const RCP<DRT::Discretization> subdis ///< subdiscretization prescribing ghosting
         )
{
  const Epetra_Map* oldcolnodemap = sourcedis->NodeColMap();

  vector<int> mycolnodes(oldcolnodemap->NumMyElements());
  oldcolnodemap->MyGlobalElements (&mycolnodes[0]);
  for (int inode = 0; inode != subdis->NumMyColNodes(); ++inode)
  {
      const DRT::Node* newnode = subdis->lColNode(inode);
      const int gid = newnode->Id();
      if (!(sourcedis->HaveGlobalNode(gid)))
      {
          mycolnodes.push_back(gid);
      }
  }

  // now reconstruct the extended colmap
  RCP<Epetra_Map> newcolnodemap = rcp(new Epetra_Map(-1,
                                     mycolnodes.size(),
                                     &mycolnodes[0],
                                     0,
                                     sourcedis->Comm()));
  return newcolnodemap;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeStructure3DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(ns[i*lrows]);

#ifdef D_SHELL8
  // special for shell8
  Epetra_SerialDenseMatrix dir;
  DRT::Element* dwele = dis.lRowElement(0);
  if (dwele->ElementObjectType()==DRT::ELEMENTS::Shell8Type::Instance())
  {
    dir.Shape(dis.NumMyRowNodes(),3);
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      DRT::ELEMENTS::Shell8* s8 =
        dynamic_cast<DRT::ELEMENTS::Shell8*>(actnode->Elements()[0]);
      if (!s8) dserror("Cannot cast to Shell8");
      int j;
      for (j=0; j<s8->NumNode(); ++j)
        if (s8->Nodes()[j]->Id() == actnode->Id()) break;
      if (j==s8->NumNode()) dserror("Can't find matching node - weird!");
      double h2 = (*s8->GetThickness())[j]/2.0;
      // get director
      const Epetra_SerialDenseMatrix* a3ref = s8->GetDirectors();
      dir(i,0) = (*a3ref)(0,j)*h2;
      dir(i,1) = (*a3ref)(1,j)*h2;
      dir(i,2) = (*a3ref)(2,j)*h2;
    }
  }
#endif

  /* the rigid body modes for structures are:
        xtrans   ytrans  ztrans   xrot       yrot       zrot
        mode[0]  mode[1] mode[2]  mode[3]    mode[4]    mode[5]
  -----------------------------------------------------------
  x   |    1       0       0       0          z-z0      -y+y0
  y   |    0       1       0      -z+z0       0          x-x0
  z   |    0       0       1       y-y0      -x+x0       0
  dx  |    0       0       0       0          a3        -a2
  dy  |    0       0       0      -a3         0          a1
  dz  |    0       0       0       a2        -a1         0
  */

  // works straight for bricks as well
//   if (ele->Type() == DRT::Element::element_shell8 ||
//       ele->Type() == DRT::Element::element_ale3 ||
//       ele->Type() == DRT::Element::element_so_hex8 ||
//       ele->Type() == DRT::Element::element_so_hex20 ||
//       ele->Type() == DRT::Element::element_so_hex27 ||
//       ele->Type() == DRT::Element::element_sosh8 ||
//       ele->Type() == DRT::Element::element_so_tet4 ||
//       ele->Type() == DRT::Element::element_so_tet10 ||
//       ele->Type() == DRT::Element::element_so_weg6 ||
//       ele->Type() == DRT::Element::element_sodisp ||
//       ele->Type() == DRT::Element::element_so_shw6 ||
//       ele->Type() == DRT::Element::element_truss3 ||
//       ele->Type() == DRT::Element::element_torsion3)
  {
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      const double* x = actnode->X();
      vector<int> dofs = dis.Dof(actnode);
      for (unsigned j=0; j<dofs.size(); ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");
        switch (j) // j is degree of freedom
        {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] = x[2] - x0[2];
          mode[5][lid] = -x[1] + x0[1];
        break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = -x[2] + x0[2];
          mode[4][lid] = 0.0;
          mode[5][lid] = x[0] - x0[0];
        break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
          mode[3][lid] = x[1] - x0[1];
          mode[4][lid] = -x[0] + x0[0];
          mode[5][lid] = 0.0;
        break;
#ifdef D_SHELL8
        case 3:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] = dir(i,2);
          mode[5][lid] = -dir(i,1);
        break;
        case 4:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = -dir(i,2);
          mode[4][lid] = 0.0;
          mode[5][lid] = dir(i,0);
        break;
        case 5:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = dir(i,1);
          mode[4][lid] = -dir(i,0);
          mode[5][lid] = 0.0;
        break;
#endif
        default:
          dserror("Only dofs 0 - 5 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // if (ele->Type() == DRT::Element::element_shell8)
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeStructure2DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(ns[i*lrows]);

//   else if (ele->Type() == DRT::Element::element_wall1 ||
//            ele->Type() == DRT::Element::element_ale2 ||
//            ele->Type() == DRT::Element::element_torsion2)
  {
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      const double* x = actnode->X();
      vector<int> dofs = dis.Dof(actnode);
      for (unsigned j=0; j<dofs.size(); ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");
        switch (j) // j is degree of freedom
        {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = -x[1] + x0[1];
        break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = x[0] - x0[0];
        break;
        default:
          dserror("Only dofs 0 - 1 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // else if (ele->Type() == DRT::Element::element_wall1)
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeBeam2DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(ns[i*lrows]);

  /* the variable mode[i][j] describes how the position of a
   * node changes with respect to the j-th degree of freedom
   * in case that the i-th rigid body mode is applied to the
   * structure; the structure altogether always has 3 rigid body
   * modes in R^2 and 6 in R^3; these modes are translation in
   * each coordinate direction, respectively, and rotation
   * around each axis, respectively. This is summed up in the
   * following table where in the left column x,y,z denote
   * translations in x-, y- and z-direction of a node due to
   * the application of a rigid body mode, whereas dx,dy,dz
   * denote increments of the node's rotational degrees of
   * freedom, which relate to a rotation around the x-,y-
   * and z-axis.
   *
        xtrans   ytrans  ztrans   xrot       yrot       zrot
        mode[0]  mode[1] mode[2]  mode[3]    mode[4]    mode[5]
  -----------------------------------------------------------
  x   |    1       0       0       0          z-z0      -y+y0
  y   |    0       1       0      -z+z0       0          x-x0
  z   |    0       0       1       y-y0      -x+x0       0
  dx  |    0       0       0       1          0          0
  dy  |    0       0       0       0          1          0
  dz  |    0       0       0       0          0          1

  for example the first line means: a translation of a node in
  x-direction may be caused either by a x-translation of the whole
  structure (which is rigid body mode 0) or by a rotation either
  around the y-axis or the z-axis. In case of a rotation dtheta around the
  y-axis for example the resulting x-translation is dtheta times the
  lever arm z - z0. Here z0 represents the z-coordinate of the point
  around which the structure is rotated. This point may be chosen
  arbitrarily and by the algorithms underlying to this method it is
  chosen automatically according to some mathematical considerations.
  Note that this holds true for infinitesimal rotations dtehta, only, of
  course.
  On the other hand e.g. the fourth column means that a rigid body
  rotation dtheta around the x-axis entails translations (-z+z0)*dtheta
  in y-direction and (y-y0)*dtheta in z-direction and a rotation
  increment 1*dtheta of the rotational degree of freedom related to the
  x-axis.
  */

  /* for beam2 and beam2r elements the above table reduces to
   *
        xtrans   ytrans    zrot
        mode[0]  mode[1]   mode[2]
  -----------------------------------------------------------
  x   |    1       0       -y+y0
  y   |    0       1       x-x0
  dz  |    0       0       1
  note: for the here employed Timoshenko and Reissner beam elements a rigid
  body rotation entails also an increment of the rotation degree of freedom
  dz which makes the director of the beam move accordingly; only then
  a rotation does not lead to any shear stress and is truely a rigid
  body rotation
  */

  /*two dimensional beam beam elements, where each node has 2 translational
   * degrees of freedom and one rotational degree of freedom*/
//   else if (ele->Type() == DRT::Element::element_beam2 ||
//            ele->Type() == DRT::Element::element_beam2r)
  {
    //looping through all nodes
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      //getting pointer at current node
      DRT::Node* actnode = dis.lRowNode(i);

      //getting coordinates of current node
      const double* x = actnode->X();

      //getting number of degrees of freedom of current node
      vector<int> dofs = dis.Dof(actnode);

      //looping through all degrees of freedom of a node
      for (unsigned j=0; j<dofs.size(); ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");
        // j is degree of freedom; each case refers to one line in the above table
        switch (j)
        {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = -x[1] + x0[1];
        break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = x[0] - x0[0];
        break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
        break;
        default:
          dserror("Only dofs 0 - 2 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // else if (ele->Type() == DRT::Element::element_beam2 || ele->Type() == DRT::Element::element_beam2r)
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeBeam3DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(ns[i*lrows]);

  /* for beam3 elements the relation between rigid body modes and
   * increments on the degrees of freedom is non-trivial since
   * rotational increments in 3D are non-additive in general. In
   * general this relation may require calling all the elements.
   * However, in opposition to the SHELL8 element it is not
   * sufficient to just call a director saved in the element.
   * Rather to calculate proper increments for the rotational
   * degrees of freedom due to a rigid body rotation of the
   * complete structure, the triad at each node is required in
   * order to transform non-additive increments into additive ones.
   * However, the beam3 element currently does not save the nodal
   * triads as a class variable, but only the triads at each Gauss
   * point. In the following a wrong (!!!) dummy version is implemneted
   * but commented out. In this dummy version the rotational degrees of
   * freedom are treated identically to the additive translational
   * degrees of freedom. Activating and using this part of the code
   * quickly reveals the problems of such a naive implemnetation.
   * Usually the equation solver simply does not work with this
   * dummy code, i.e. the iterative solution process does not converge.
   * If Algebraic Multigrid methods should be really used for beam3
   * elements, one first has to develop efficient special methods for
   * these elements. Currently trying to use Algebraic multigrid methods
   * for beam3 elements just amounts to an error as no properly working
   * implementation has been available so far*/

//   else if (ele->Type() == DRT::Element::element_beam3 || ele->Type() == DRT::Element::element_beam3ii)
  {
    //looping through all nodes
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      //getting pointer at current node
      DRT::Node* actnode = dis.lRowNode(i);

      //getting coordinates of current node
      const double* x = actnode->X();

      //getting number of degrees of freedom of current node
      vector<int> dofs = dis.Dof(actnode);

      //looping through all degrees of freedom of a node
      for (unsigned j=0; j<dofs.size(); ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");
        // j is degree of freedom; each case refers to one line in the above table
        switch (j)
        {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] =  x[2] - x0[2];
          mode[5][lid] = -x[1] + x0[1];
        break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = -x[2] + x0[2];
          mode[4][lid] = 0.0;
          mode[5][lid] =  x[0] - x0[0];
        break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
          mode[3][lid] =  x[1] - x0[1];
          mode[4][lid] = -x[0] + x0[0];
          mode[5][lid] = 0.0;
        break;
        case 3:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 1.0;
          mode[4][lid] = 0.0;
          mode[5][lid] = 0.0;
        break;
        case 4:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] = 1.0;
          mode[5][lid] = 0.0;
        break;
        case 5:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] = 0.0;
          mode[5][lid] = 1.0;
        break;
        default:
          dserror("Only dofs 0 - 5 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)

  } // else if (ele->Type() == DRT::Element::element_beam3 || ele->Type() == DRT::Element::element_beam3ii)
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeXFluid3DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(ns[i*lrows]);

  /* the rigid body modes for fluids are:
        xtrans   ytrans  ztrans   pressure
        mode[0]  mode[1] mode[2]  mode[3]
  ----------------------------------------
  x   |    1       0       0       0
  y   |    0       1       0       0
  z   |    0       0       1       0
  p   |    0       0       0       1
  */

//   else if (ele->Type() == DRT::Element::element_xfluid3 ||
//            ele->Type() == DRT::Element::element_combust3 ||
//            ele->Type() == DRT::Element::element_smoothrod ||
//            ele->Type() == DRT::Element::element_sosh8p8)
  {
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      vector<int> dofs = dis.Dof(actnode);
      for (unsigned j=0; j<dofs.size(); ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");
        switch (j) // j is degree of freedom
        {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
          mode[3][lid] = 0.0;
        break;
        case 3:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 1.0;
        break;
        case 4:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        case 5:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        case 6:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        case 7:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        default:
          dserror("Only dofs 0 - 7 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // else if (ele->Type() == DRT::Element::element_fluid3 or ele->Type() == DRT::Element::element_xfluid3)
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeFluid2DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(ns[i*lrows]);

//   else if (ele->Type() == DRT::Element::element_fluid2)
  {
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      vector<int> dofs = dis.Dof(actnode);
      for (unsigned j=0; j<dofs.size(); ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");
        switch (j) // j is degree of freedom
        {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
        break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = 0.0;
        break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
        break;
        default:
          dserror("Only dofs 0 - 2 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // else if (ele->Type() == DRT::Element::element_fluid2)
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeFluid3DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(ns[i*lrows]);

//   else if (ele->Type() == DRT::Element::element_transport or
//       ele->Type() == DRT::Element::element_fluid3)
  {
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      vector<int> dofs = dis.Dof(actnode);
      const unsigned int ndof = dofs.size();
      for (unsigned j=0; j<ndof; ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");

        for (unsigned k=0; k<ndof; ++k)
        {
          if (k == j)
            mode[k][lid] = 1.0;
          else
            mode[k][lid] = 0.0;
        }

      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // else if (ele->Type() == DRT::Element::element_transport)
}


#endif  // #ifdef CCADISCRET
