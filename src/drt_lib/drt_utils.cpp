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

#include "linalg_utils.H"
#include "linalg_mapextractor.H"

#include "drt_utils.H"
#include "drt_node.H"
#include "../drt_nurbs_discret/drt_control_point.H"
#include "drt_dofset.H"
#include "drt_discret.H"
#include "../drt_beam2/beam2.H"
#include "../drt_beam2r/beam2r.H"
#include "../drt_beam3/beam3.H"
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
#include "../drt_bele3/bele3.H"
#include "../drt_bele3/vele3.H"
#include "../drt_constraint/constraint_element2.H"
#include "../drt_constraint/constraint_element3.H"
#include "../drt_w1/wall1.H"
#include "../drt_w1/wall1_nurbs.H"
#include "../drt_so3/so_hex8.H"
#include "../drt_so3/so_hex20.H"
#include "../drt_so3/so_hex27.H"
#include "../drt_so3/so_sh8.H"
#include "../drt_so3/so_sh8p8.H"
#include "../drt_so3/so_tet4.H"
#include "../drt_so3/so_ptet.H"
#include "../drt_so3/so_tet10.H"
#include "../drt_so3/so_weg6.H"
#include "../drt_so3/so_shw6.H"
#include "../drt_so3/so_disp.H"
#include "../drt_so3/so_hex8p1j1.H"
#include "../drt_thermo/thermo_element.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/micromaterial.H"
#include "../drt_mat/neohooke.H"
#include "../drt_mat/logneohooke.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/aaaraghavanvorp_damage.H"
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
#include "../drt_mat/mixfrac.H"
#include "../drt_mat/sutherland.H"
#include "../drt_mat/arrhenius_spec.H"
#include "../drt_mat/arrhenius_temp.H"
#include "../drt_mat/arrhenius_pv.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"
#include "../drt_mat/biocell.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/charmm.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_mat/cnst_1d_art.H"
#include "../drt_mat/fourieriso.H"
#include "../drt_mat/holzapfelcardiovascular.H"
#include "../drt_mat/itskov.H"
#include "../drt_mat/plasticneohooke.H"
#include "../drt_contact/drt_cnode.H"
#include "../drt_contact/drt_celement.H"
#include "../drt_art_net/artery.H"
#include "drt_dserror.H"
#include "standardtypes_cpp.H"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;



/*----------------------------------------------------------------------*
 |  allocate an instance of a specific impl. of ParObject (public) mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::UTILS::Factory(const vector<char>& data)
{
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
    case ParObject_ScatraMat:
    {
      MAT::ScatraMat* scatra_mat = new MAT::ScatraMat();
      scatra_mat->Unpack(data);
      return scatra_mat;
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
    case ParObject_CNode:
    {
      double x[3];
      vector<int> dofs(0);
      CONTACT::CNode* node = new CONTACT::CNode(0,x,0,0,dofs,false,false);
      node->Unpack(data);
      return node;
    }
    case ParObject_CElement:
    {
      CONTACT::CElement* ele = new CONTACT::CElement(0,
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
    default:
      dserror("Unknown type of ParObject instance: %d",type);
    break;
  }

  return NULL;
}

/*----------------------------------------------------------------------*
 |  allocate an element of a specific type (public)          mwgee 03|07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::Element> DRT::UTILS::Factory(const string eletype,
                                              const string eledistype,
                                              const int id,
                                              const int owner)
{
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
    bele3,
    vele3,
    so_hex8,
    so_sh8,
    so_sh8p8,
    so_hex27,
    so_hex20,
    so_tet4,
    ptet,
    so_tet10,
    so_weg6,
    so_shw6,
    sodisp,
    beam2,
    beam2r,
    beam3,
    truss3,
    truss2,
    torsion3,
    torsion2,
    constrele2,
    constrele3,
    transport,
    art_ele,
    so_hex8p1j1,
    thermo
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
  else if (eletype=="BELE3") type = bele3;
  else if (eletype=="VELE3") type = vele3;
  else if (eletype=="SOLIDH8") type = so_hex8;
  else if (eletype=="SOLIDSH8") type = so_sh8;
  else if (eletype=="SOLIDSH8P8") type = so_sh8p8;
  else if (eletype=="SOLIDH27") type = so_hex27;
  else if (eletype=="SOLIDH20") type = so_hex20;
  else if (eletype=="SOLIDT4") type = so_tet4;
  else if (eletype=="PTET4") type = ptet;
  else if (eletype=="SOLIDT10") type = so_tet10;
  else if (eletype=="SOLIDW6") type = so_weg6;
  else if (eletype=="SOLIDSHW6") type = so_shw6;
  else if (eletype=="SOLID3") type = sodisp;
  else if (eletype=="BEAM2") type = beam2;
  else if (eletype=="BEAM2R") type = beam2r;
  else if (eletype=="BEAM3") type = beam3;
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
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::Beam2(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_BEAM2R
    case beam2r:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::Beam2r(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_BEAM3
    case beam3:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::Beam3(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_TRUSS3
    case truss3:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::Truss3(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_TRUSS2
    case truss2:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::Truss2(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_TORSION3
    case torsion3:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::Torsion3(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_TORSION2
    case torsion2:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::Torsion2(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_SHELL8
    case shell8:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::Shell8(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_FLUID2
    case fluid2:
    {
      RefCountPtr<DRT::Element> ele;

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
      RefCountPtr<DRT::Element> ele;

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
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::XFluid3(id,owner));
      return ele;
    }
    break;
    case xdiff3:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::XDiff3(id,owner));
      return ele;
    }
    break;
    case combust3:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::Combust3(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_ALE
    case ale3:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::Ale3(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_ALE
    case ale2:
    {
      RefCountPtr<DRT::Element> ele;

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
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::Bele3(id,owner));
      return ele;
    }
    break;
    case vele3:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::Vele3(id,owner));
      return ele;
    }
    break;
    case constrele2:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::ConstraintElement2(id,owner));
      return ele;
    }
    break;
    case constrele3:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::ConstraintElement3(id,owner));
      return ele;
    }
    break;
#ifdef D_WALL1
    case wall1:
    {
      RefCountPtr<DRT::Element> ele;
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
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_hex8(id,owner));
      return ele;
    }
    break;
    case so_sh8:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_sh8(id,owner));
      return ele;
    }
    break;
    case so_sh8p8:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_sh8p8(id,owner));
      return ele;
    }
    break;
    case so_hex27:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_hex27(id,owner));
      return ele;
    }
    break;
    case so_hex20:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_hex20(id,owner));
      return ele;
    }
    break;
    case so_weg6:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_weg6(id,owner));
      return ele;
    }
    break;
    case so_shw6:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_shw6(id,owner));
      return ele;
    }
    break;
    case sodisp:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::SoDisp(id,owner));
      return ele;
    }
    break;
    case so_tet4:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_tet4(id,owner));
      return ele;
    }
    break;
    case ptet:
    {
      RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Ptet(id,owner));
      return ele;
    }
    break;
    case so_tet10:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_tet10(id,owner));
      return ele;
    }
    break;
    case so_hex8p1j1:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_Hex8P1J1(id,owner));
      return ele;
    }
    break;
#endif
#if defined(D_FLUID2) || defined(D_FLUID3)
    case transport:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::Transport(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_ARTNET
    case art_ele:
    {
       RefCountPtr<DRT::Element> ele =  rcp(new DRT::ELEMENTS::Artery(id,owner));
       return ele;
    }
    break;
#endif
#ifdef D_THERMO
    case thermo:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::ELEMENTS::Thermo(id,owner));
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
}


/*----------------------------------------------------------------------*
 |  partition a graph using metis  (public)                  mwgee 11/06|
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_CrsGraph> DRT::UTILS::PartGraphUsingMetis(
                                             const Epetra_CrsGraph& graph,
                                             const Epetra_Vector& weights)
{
  const int myrank   = graph.Comm().MyPID();
  const int numproc  = graph.Comm().NumProc();

  if (numproc==1)
  {
    RefCountPtr<Epetra_CrsGraph> outgraph = rcp(new Epetra_CrsGraph(graph));
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
  RefCountPtr<Epetra_CrsGraph> outgraph =
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

#endif  // #ifdef CCADISCRET
