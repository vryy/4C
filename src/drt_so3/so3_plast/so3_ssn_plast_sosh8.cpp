/*----------------------------------------------------------------------*/
/*!
\brief
This file contains everything additionally needed for plastic solid shell elements.
Basically, it is a copy of the elastic so_sh8 element. (Virtual inheritance in the
solid elements would help to reduce this code redundancy.)

The Solid-Shell element technology is based on the work of
(1) Vu-Quoc, Tan: "Optimal solid shells for non-linear analyses
                   of multilayer composites", CMAME 2003
(2) Klinkel, Gruttmann, Wagner: "A robust non-linear solid shell element
                                 based on a mixed variational fromulation"

Refer also to the Semesterarbeit of Alexander Popp, 2006

\maintainer Matthias Mayr
\level 2
*/

/*----------------------------------------------------------------------*/

#include "so3_ssn_plast_sosh8.H"
#include "so3_ssn_plast.H"
#include "../../drt_lib/drt_linedefinition.H"
#include "../../drt_mat/plasticelasthyper.H"
#include "Epetra_SerialDenseSolver.h"

#include "../../drt_structure_new/str_elements_paramsinterface.H"
#include "../../linalg/linalg_serialdensematrix.H"
#include "../../linalg/linalg_serialdensevector.H"
#include "../../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 | build an instance of plast type                         seitz 05/14 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh8PlastType DRT::ELEMENTS::So_sh8PlastType::instance_;
std::pair<bool, LINALG::Matrix<DRT::ELEMENTS::So3_Plast<DRT::Element::hex8>::nsd_,
                    DRT::ELEMENTS::So3_Plast<DRT::Element::hex8>::nsd_>>
    DRT::ELEMENTS::So_sh8Plast::jac_refe_;
std::pair<bool, LINALG::Matrix<DRT::ELEMENTS::So3_Plast<DRT::Element::hex8>::nsd_,
                    DRT::ELEMENTS::So3_Plast<DRT::Element::hex8>::nsd_>>
    DRT::ELEMENTS::So_sh8Plast::jac_curr_;
std::pair<bool,
    LINALG::Matrix<DRT::ELEMENTS::So_sh8Plast::num_ans * DRT::ELEMENTS::So_sh8Plast::num_sp,
        DRT::ELEMENTS::So_sh8Plast::numdofperelement_>>
    DRT::ELEMENTS::So_sh8Plast::B_ans_loc_;
std::pair<bool,
    LINALG::Matrix<DRT::ELEMENTS::So_sh8Plast::numstr_, DRT::ELEMENTS::So_sh8Plast::numstr_>>
    DRT::ELEMENTS::So_sh8Plast::TinvT_;



DRT::ELEMENTS::So_sh8PlastType& DRT::ELEMENTS::So_sh8PlastType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 05/14 |
| is called in ElementRegisterType                                     |
*----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::So_sh8PlastType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So_sh8Plast* object = new DRT::ELEMENTS::So_sh8Plast(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 05/14 |
| is called from ParObjectFactory                                      |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_sh8PlastType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDSH8PLAST")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_sh8Plast(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 05/14 |
| virtual method of ElementType                                        |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_sh8PlastType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_sh8Plast(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
| initialise the element (public)                          seitz 05/14 |
*----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_sh8PlastType::Initialize(DRT::Discretization& dis)
{
  // sosh8_gmshplotdis(dis);

  int num_morphed_so_hex8_easmild = 0;
  int num_morphed_so_hex8_easnone = 0;

  // Loop through all elements
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    // get the actual element
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So_sh8Plast* actele =
        dynamic_cast<DRT::ELEMENTS::So_sh8Plast*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_sh8* failed");

    if (!actele->nodes_rearranged_)
    {
      switch (actele->thickdir_)
      {
        // check for automatic definition of thickness direction
        case DRT::ELEMENTS::So_sh8Plast::autoj:
        {
          actele->thickdir_ = actele->findthickdir();
          break;
        }
        // check for enforced definition of thickness direction
        case DRT::ELEMENTS::So_sh8Plast::globx:
        {
          LINALG::Matrix<NUMDIM_SOH8, 1> thickdirglo(true);
          thickdirglo(0) = 1.0;
          actele->thickdir_ = actele->enfthickdir(thickdirglo);
          break;
        }
        case DRT::ELEMENTS::So_sh8Plast::globy:
        {
          LINALG::Matrix<NUMDIM_SOH8, 1> thickdirglo(true);
          thickdirglo(1) = 1.0;
          actele->thickdir_ = actele->enfthickdir(thickdirglo);
          break;
        }
        case DRT::ELEMENTS::So_sh8Plast::globz:
        {
          LINALG::Matrix<NUMDIM_SOH8, 1> thickdirglo(true);
          thickdirglo(2) = 1.0;
          actele->thickdir_ = actele->enfthickdir(thickdirglo);
          break;
        }
        default:
          break;
      }

      int new_nodeids[NUMNOD_SOH8];

      switch (actele->thickdir_)
      {
        case DRT::ELEMENTS::So_sh8Plast::globx:
        case DRT::ELEMENTS::So_sh8Plast::globy:
        case DRT::ELEMENTS::So_sh8Plast::globz:
        {
          dserror("This should have been replaced by auto(r|s|t)");
          break;
        }
        case DRT::ELEMENTS::So_sh8Plast::autor:
        case DRT::ELEMENTS::So_sh8Plast::enfor:
        {
          // resorting of nodes,
          // such that previous local r-dir is local t-dir afterwards
          new_nodeids[0] = actele->NodeIds()[7];
          new_nodeids[1] = actele->NodeIds()[4];
          new_nodeids[2] = actele->NodeIds()[0];
          new_nodeids[3] = actele->NodeIds()[3];
          new_nodeids[4] = actele->NodeIds()[6];
          new_nodeids[5] = actele->NodeIds()[5];
          new_nodeids[6] = actele->NodeIds()[1];
          new_nodeids[7] = actele->NodeIds()[2];
          actele->SetNodeIds(NUMNOD_SOH8, new_nodeids);
          actele->nodes_rearranged_ = true;
          break;
        }
        case DRT::ELEMENTS::So_sh8Plast::autos:
        case DRT::ELEMENTS::So_sh8Plast::enfos:
        {
          // resorting of nodes,
          // such that previous local s-dir is local t-dir afterwards
          new_nodeids[0] = actele->NodeIds()[4];
          new_nodeids[1] = actele->NodeIds()[5];
          new_nodeids[2] = actele->NodeIds()[1];
          new_nodeids[3] = actele->NodeIds()[0];
          new_nodeids[4] = actele->NodeIds()[7];
          new_nodeids[5] = actele->NodeIds()[6];
          new_nodeids[6] = actele->NodeIds()[2];
          new_nodeids[7] = actele->NodeIds()[3];
          actele->SetNodeIds(NUMNOD_SOH8, new_nodeids);
          actele->nodes_rearranged_ = true;
          break;
        }
        case DRT::ELEMENTS::So_sh8Plast::autot:
        case DRT::ELEMENTS::So_sh8Plast::enfot:
        {
          // no resorting necessary
          for (int node = 0; node < 8; ++node)
          {
            new_nodeids[node] = actele->NodeIds()[node];
          }
          actele->SetNodeIds(NUMNOD_SOH8, new_nodeids);
          actele->nodes_rearranged_ = true;
          break;
        }
        case DRT::ELEMENTS::So_sh8Plast::undefined:
        {
          if (actele->eastype_ == DRT::ELEMENTS::soh8p_eassosh8)
          {
            // here comes plan B: morph So_sh8 to So_hex8
            actele->ReInitEas(DRT::ELEMENTS::soh8p_easmild);
            actele->anstype_ = So_sh8Plast::ansnone_p;
            actele->InitJacobianMapping();
            num_morphed_so_hex8_easmild++;
          }
          else if (actele->eastype_ == DRT::ELEMENTS::soh8p_easnone)
          {
            // here comes plan B: morph So_sh8 to So_hex8
            actele->ReInitEas(DRT::ELEMENTS::soh8p_easnone);
            actele->anstype_ = So_sh8Plast::ansnone_p;
            actele->InitJacobianMapping();
            num_morphed_so_hex8_easnone++;
          }
          else if (actele->eastype_ == DRT::ELEMENTS::soh8p_easmild)
          {
            // this might happen in post filter (for morped sosh8->soh8)
            actele->ReInitEas(DRT::ELEMENTS::soh8p_easmild);
            actele->anstype_ = So_sh8Plast::ansnone_p;
            actele->InitJacobianMapping();
          }
          else if (actele->eastype_ == DRT::ELEMENTS::soh8p_easnone)
          {
            // this might happen in post filter (for morped sosh8->soh8)
            actele->anstype_ = So_sh8Plast::ansnone_p;
            actele->InitJacobianMapping();
          }
          else
            dserror("Undefined EAS type");
          break;
        }
        case DRT::ELEMENTS::So_sh8Plast::none:
          break;
        default:
          dserror("no thickness direction for So_sh8");
          break;
      }
    }
  }

  if (num_morphed_so_hex8_easmild > 0)
  {
    std::cout << std::endl
              << num_morphed_so_hex8_easmild
              << " Sosh8Plast-Elements have no clear 'thin' direction and have morphed to "
                 "So_hex8Plast with eas_mild"
              << std::endl;
  }
  if (num_morphed_so_hex8_easnone > 0)
  {
    std::cout << std::endl
              << num_morphed_so_hex8_easnone
              << " Sosh8Plast-Elements have no clear 'thin' direction and have morphed to "
                 "So_hex8Plast with eas_none"
              << std::endl;
  }

  // fill complete again to reconstruct element-node pointers,
  // but without element init, etc.
  dis.FillComplete(false, false, false);

  // loop again to init Jacobian for Sosh8's
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So_sh8Plast* actele =
        dynamic_cast<DRT::ELEMENTS::So_sh8Plast*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_sh8* failed");
    actele->InitJacobianMapping();
  }

  return 0;
}

/*---------------------------------------------------------------------*
|                                                          seitz 05/14 |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8PlastType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDSH8PLAST"];

  defs["HEX8"]
      .AddIntVector("HEX8", 8)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddNamedString("EAS")
      .AddNamedString("ANS")
      .AddNamedString("THICKDIR")
      .AddOptionalNamedDoubleVector("FIBER1", 3)
      .AddOptionalNamedDoubleVector("FIBER2", 3)
      .AddOptionalNamedDoubleVector("FIBER3", 3);

}  // SetupElementDefinition()

/*----------------------------------------------------------------------*
 | ctor (public)                                            seitz 05/14 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh8Plast::So_sh8Plast(int id, int owner)
    : So_base(id, owner), DRT::ELEMENTS::So3_Plast<DRT::Element::hex8>(id, owner)
{
  thickdir_ = globx;
  nodes_rearranged_ = false;
  thickvec_.resize(3, 0.);
  return;
}

/*----------------------------------------------------------------------*
 | copy-ctor (public)                                       seitz 05/14 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh8Plast::So_sh8Plast(const DRT::ELEMENTS::So_sh8Plast& old)
    : So_base(old), DRT::ELEMENTS::So3_Plast<DRT::Element::hex8>(old)
{
  return;
}

/*----------------------------------------------------------------------*
 | deep copy this instance of Solid3 and return pointer to              |
 | it (public)                                              seitz 05/14 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_sh8Plast::Clone() const
{
  DRT::ELEMENTS::So_sh8Plast* newelement = new DRT::ELEMENTS::So_sh8Plast(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 | pack data (public)                                       seitz 05/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8Plast::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class So3_Plast Element
  DRT::ELEMENTS::So3_Plast<DRT::Element::hex8>::Pack(data);
  // thickdir
  AddtoPack(data, thickdir_);
  AddtoPack(data, thickvec_);
  AddtoPack(data, anstype_);
  AddtoPack(data, nodes_rearranged_);

  return;
}

/*----------------------------------------------------------------------*
 | unpack data (public)                                     seitz 05/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8Plast::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class So_hex8 Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::ELEMENTS::So3_Plast<DRT::Element::hex8>::Unpack(basedata);
  // thickdir
  thickdir_ = static_cast<ThicknessDirection>(ExtractInt(position, data));
  ExtractfromPack(position, data, thickvec_);
  anstype_ = static_cast<ANSType>(ExtractInt(position, data));
  nodes_rearranged_ = ExtractInt(position, data);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

DRT::ELEMENTS::So_sh8Plast::~So_sh8Plast() { return; }

void DRT::ELEMENTS::So_sh8Plast::Print(std::ostream& os) const
{
  os << "So_sh8Plast ";
  Element::Print(os);
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 | read this element, get the material (public)             seitz 05/14 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_sh8Plast::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  std::string buffer;
  linedef->ExtractString("KINEM", buffer);

  // geometrically linear
  if (buffer == "linear")
  {
    dserror("no linear kinematics");
  }
  // geometrically non-linear with Total Lagrangean approach
  else if (buffer == "nonlinear")
  {
    kintype_ = INPAR::STR::kinem_nonlinearTotLag;
    // everything ok
  }
  else
    dserror("Reading of SO3_PLAST element failed! KINEM unknown");

  DRT::UTILS::GaussIntegration ip(DRT::Element::hex8, 3);
  numgpt_ = ip.NumPoints();
  xsi_.resize(numgpt_);
  wgt_.resize(numgpt_);
  for (int gp = 0; gp < numgpt_; ++gp)
  {
    wgt_[gp] = ip.Weight(gp);
    const double* gpcoord = ip.Point(gp);
    for (int idim = 0; idim < nsd_; idim++) xsi_[gp](idim) = gpcoord[idim];
  }

  // no fbar
  fbar_ = false;

  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);

  SetMaterial(material);

  Teuchos::RCP<MAT::So3Material> so3mat = SolidMaterial();
  so3mat->Setup(numgpt_, linedef);
  so3mat->ValidKinematics(INPAR::STR::kinem_nonlinearTotLag);
  if (HavePlasticSpin())
    plspintype_ = plspin;
  else
    plspintype_ = zerospin;

  // EAS
  linedef->ExtractString("EAS", buffer);

  if (buffer == "none")
    eastype_ = soh8p_easnone;
  else if (buffer == "sosh8")
    eastype_ = soh8p_eassosh8;
  else
    dserror("unknown EAS type for so3_plast");

  // initialize EAS data
  EasInit();

  // read ANS technology flag
  linedef->ExtractString("ANS", buffer);
  if (buffer == "sosh8")
  {
    anstype_ = anssosh8_p;
  }
  // no ANS technology
  else if (buffer == "none")
  {
    anstype_ = ansnone_p;
  }
  else
    dserror("Reading of SO_SH8 ANS technology failed");

  linedef->ExtractString("THICKDIR", buffer);
  nodes_rearranged_ = false;

  // global X
  if (buffer == "xdir") thickdir_ = globx;
  // global Y
  else if (buffer == "ydir")
    thickdir_ = globy;
  // global Z
  else if (buffer == "zdir")
    thickdir_ = globz;
  // find automatically through Jacobian of Xrefe
  else if (buffer == "auto")
    thickdir_ = autoj;
  // local r
  else if (buffer == "rdir")
    thickdir_ = enfor;
  // local s
  else if (buffer == "sdir")
    thickdir_ = enfos;
  // local t
  else if (buffer == "tdir")
    thickdir_ = enfot;
  // no noderearrangement
  else if (buffer == "none")
  {
    thickdir_ = none;
    nodes_rearranged_ = true;
  }
  else
    dserror("Reading of SO_SH8 thickness direction failed");

  // plasticity related stuff
  KbbInv_.resize(numgpt_, LINALG::SerialDenseMatrix(plspintype_, plspintype_, true));
  Kbd_.resize(numgpt_, LINALG::SerialDenseMatrix(plspintype_, numdofperelement_, true));
  fbeta_.resize(numgpt_, LINALG::SerialDenseVector(plspintype_, true));
  dDp_last_iter_.resize(numgpt_, LINALG::SerialDenseVector(plspintype_, true));
  dDp_inc_.resize(numgpt_, LINALG::SerialDenseVector(plspintype_, true));

  Teuchos::ParameterList plparams = DRT::Problem::Instance()->SemiSmoothPlastParams();
  plparams.set<PROBLEM_TYP>("PROBLEM_TYP", DRT::Problem::Instance()->ProblemType());
  ReadParameterList(Teuchos::rcpFromRef<Teuchos::ParameterList>(plparams));

  if (tsi_)
    dserror(
        "no tsi for solid shells.\n"
        "Go back to the revision introducing this error to find a version with TSI.\n"
        "Note: the strain rate used for the Gough Joule effect was not calculated\n"
        "consistently with the ANS modifications");

  return true;
}

/*----------------------------------------------------------------------*
 |                                                          seitz 05/14 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh8Plast::ThicknessDirection DRT::ELEMENTS::So_sh8Plast::findthickdir()
{
  // update element geometry
  LINALG::Matrix<nen_, nsd_> xrefe(false);  // material coord. of element
  for (int i = 0; i < nen_; ++i)
  {
    xrefe(i, 0) = this->Nodes()[i]->X()[0];
    xrefe(i, 1) = this->Nodes()[i]->X()[1];
    xrefe(i, 2) = this->Nodes()[i]->X()[2];
  }
  // vector of df(origin), ie parametric derivatives of shape functions
  // evaluated at the origin (r,s,t)=(0,0,0)
  const double df0_vector[] = {-0.125, -0.125, -0.125, +0.125, -0.125, -0.125, +0.125, +0.125,
      -0.125, -0.125, +0.125, -0.125, -0.125, -0.125, +0.125, +0.125, -0.125, +0.125, +0.125,
      +0.125, +0.125, -0.125, +0.125, +0.125};
  // shape function derivatives, evaluated at origin (r=s=t=0.0)
  LINALG::Matrix<nsd_, nen_> df0(df0_vector);

  // compute Jacobian, evaluated at element origin (r=s=t=0.0)
  // (J0_i^A) = (X^A_{,i})^T
  LINALG::Matrix<nsd_, nsd_> jac0;
  jac0.MultiplyNN(df0, xrefe);
  // compute inverse of Jacobian at element origin
  // (Jinv0_A^i) = (X^A_{,i})^{-T}
  LINALG::Matrix<nsd_, nsd_> iJ0(jac0);
  iJ0.Invert();

  // separate "stretch"-part of J-mapping between parameter and global space
  // (G0^ji) = (Jinv0^j_B) (krondelta^BA) (Jinv0_A^i)
  LINALG::Matrix<nsd_, nsd_> jac0stretch;
  jac0stretch.MultiplyTN(iJ0, iJ0);
  const double r_stretch = sqrt(jac0stretch(0, 0));
  const double s_stretch = sqrt(jac0stretch(1, 1));
  const double t_stretch = sqrt(jac0stretch(2, 2));

  // minimal stretch equivalents with "thinnest" direction
  // const double max_stretch = max(r_stretch, max(s_stretch, t_stretch));
  double max_stretch = -1.0;

  ThicknessDirection thickdir = none;  // of actual element
  int thick_index = -1;

  if (r_stretch >= s_stretch and r_stretch >= t_stretch)
  {
    max_stretch = r_stretch;
    if ((max_stretch / s_stretch <= 1.5) || (max_stretch / t_stretch <= 1.5))
    {
      // std::cout << "ID: " << this->Id() << ", has aspect ratio of: ";
      // std::cout << max_stretch / s_stretch << " , " << max_stretch / t_stretch << std::endl;
      // dserror("Solid-Shell element geometry has not a shell aspect ratio");
      return undefined;
    }
    thickdir = autor;
    thick_index = 0;
  }
  else if (s_stretch > r_stretch and s_stretch >= t_stretch)
  {
    max_stretch = s_stretch;
    if ((max_stretch / r_stretch <= 1.5) || (max_stretch / t_stretch <= 1.5))
    {
      // std::cout << "ID: " << this->Id() << ", has aspect ratio of: ";
      // std::cout << max_stretch / r_stretch << " , " << max_stretch / t_stretch << std::endl;
      // dserror("Solid-Shell element geometry has not a shell aspect ratio");
      return undefined;
    }
    thickdir = autos;
    thick_index = 1;
  }
  else if (t_stretch > r_stretch and t_stretch > s_stretch)
  {
    max_stretch = t_stretch;
    if ((max_stretch / r_stretch <= 1.5) || (max_stretch / s_stretch <= 1.5))
    {
      // std::cout << "ID: " << this->Id() << ", has aspect ratio of: ";
      // std::cout << max_stretch / r_stretch << " , " << max_stretch / s_stretch << std::endl;
      // dserror("Solid-Shell element geometry has not a shell aspect ratio");
      return undefined;
    }
    thickdir = autot;
    thick_index = 2;
  }

  if (thick_index == -1)
    dserror("Trouble with thick_index=%d %g,%g,%g,%g", thick_index, r_stretch, s_stretch, t_stretch,
        max_stretch);

  // thickness-vector in parameter-space, has 1.0 in thickness-coord
  LINALG::Matrix<nsd_, 1> loc_thickvec(true);
  loc_thickvec(thick_index) = 1.0;
  // thickness-vector in global coord is J times local thickness-vector
  // (X^A) = (J0_i^A)^T . (xi_i)
  LINALG::Matrix<nsd_, 1> glo_thickvec;
  glo_thickvec.MultiplyTN(jac0, loc_thickvec);
  // return doubles of thickness-vector
  thickvec_.resize(3);
  thickvec_[0] = glo_thickvec(0);
  thickvec_[1] = glo_thickvec(1);
  thickvec_[2] = glo_thickvec(2);

  return thickdir;
}

/*----------------------------------------------------------------------*
 |                                                          seitz 05/14 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh8Plast::ThicknessDirection DRT::ELEMENTS::So_sh8Plast::enfthickdir(
    LINALG::Matrix<nsd_, 1>& thickdirglo)
{
  // update element geometry
  LINALG::Matrix<nen_, nsd_> xrefe(false);  // material coord. of element
  for (int i = 0; i < nen_; ++i)
  {
    xrefe(i, 0) = this->Nodes()[i]->X()[0];
    xrefe(i, 1) = this->Nodes()[i]->X()[1];
    xrefe(i, 2) = this->Nodes()[i]->X()[2];
  }
  // vector of df(origin), ie parametric derivatives of shape functions
  // evaluated at the origin (r,s,t)=(0,0,0)
  const double df0_vector[numdofperelement_ * nen_] = {-0.125, -0.125, -0.125, +0.125, -0.125,
      -0.125, +0.125, +0.125, -0.125, -0.125, +0.125, -0.125, -0.125, -0.125, +0.125, +0.125,
      -0.125, +0.125, +0.125, +0.125, +0.125, -0.125, +0.125, +0.125};
  // shape function derivatives, evaluated at origin (r=s=t=0.0)
  LINALG::Matrix<nsd_, nen_> df0(df0_vector);

  // compute Jacobian, evaluated at element origin (r=s=t=0.0)
  // (J0_i^A) = (X^A_{,i})^T
  LINALG::Matrix<nsd_, nsd_> jac0(false);
  jac0.MultiplyNN(df0, xrefe);

  // compute inverse of Jacobian at element origin
  // (Jinv0_A^i) = (X^A_{,i})^{-T}
  LINALG::Matrix<nsd_, nsd_> iJ0(jac0);
  iJ0.Invert();

  // make enforced global thickness direction a unit vector
  const double thickdirglolength = thickdirglo.Norm2();
  thickdirglo.Scale(1.0 / thickdirglolength);

  // pull thickness direction from global to contra-variant local
  // (dxi^i) = (Jinv0_A^i)^T . (dX^A)
  LINALG::Matrix<nsd_, 1> thickdirlocsharp(false);
  thickdirlocsharp.MultiplyTN(iJ0, thickdirglo);

#if 0
  // metric tensor
  // (G0_ji) = (J0_j^B) (delta_BA) (J0^A_i)
  LINALG::Matrix<nsd_,nsd_> metrflat(false);
  metrflat.MultiplyNT(jac0,jac0);

  // co-variant local enforced thickness direction
  // (dxi_j) = (G0_ji) (dxi^i)
  LINALG::Matrix<nsd_,1> thickdirlocflat(false);
  thickdirlocflat.MultiplyNN(metrflat,thickdirlocsharp);

  // thickdirloclength
  const double thickdirloclength = thickdirlocsharp.Dot(thickdirlocflat);

  // check if transformation was successful
  if (fabs(thickdirglolength-thickdirloclength)>EPS6)
    dserror("Transformation erroneous: Vector length is not in-variant: %g!=%g",
            thickdirglolength, thickdirloclength);
#endif

  // identify parametric co-ordinate closest to enforced thickness direction
  int thick_index = -1;
  double thickdirlocmax = 0.0;
  for (int i = 0; i < nsd_; ++i)
  {
    if (fabs(thickdirlocsharp(i)) > thickdirlocmax)
    {
      thickdirlocmax = fabs(thickdirlocsharp(i));
      thick_index = i;
    }
  }
  const double tol = 0.9;  // should be larger than 1/sqrt(2)=0.707
  // check if parametric co-ordinate is clear
  if (thickdirlocmax < tol * thickdirlocsharp.Norm2())
    dserror(
        "could not clearly identify a parametric direction pointing along enforced thickness "
        "direction");

  ThicknessDirection thickdir = none;  // of actual element
  if (thick_index == 0)
  {
    thickdir = autor;
  }
  else if (thick_index == 1)
  {
    thickdir = autos;
  }
  else if (thick_index == 2)
  {
    thickdir = autot;
  }
  else
  {
    dserror("Trouble with thick_index=%g", thick_index);
  }

  // thickness-vector in parameter-space, has 1.0 in thickness-coord
  LINALG::Matrix<nsd_, 1> loc_thickvec(true);
  loc_thickvec(thick_index) = 1.0;
  // thickness-vector in global coord is J times local thickness-vector
  // (X^A) = (J0_i^A)^T . (xi_i)
  LINALG::Matrix<nsd_, 1> glo_thickvec;
  glo_thickvec.MultiplyTN(jac0, loc_thickvec);
  // return doubles of thickness-vector
  thickvec_.resize(3);
  thickvec_[0] = glo_thickvec(0);
  thickvec_[1] = glo_thickvec(1);
  thickvec_[2] = glo_thickvec(2);

  return thickdir;
}

/*----------------------------------------------------------------------*
 |  evaluate 'T'-transformation matrix )                       maf 05/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8Plast::EvaluateT(
    const LINALG::Matrix<nsd_, nsd_>& jac, LINALG::Matrix<numstr_, numstr_>& TinvT)
{
  // build T^T transformation matrix which maps
  // between global (r,s,t)-coordinates and local (x,y,z)-coords
  // later, invert the transposed to map from local to global
  // see literature for details (e.g. Andelfinger)
  // it is based on the voigt notation for strains: xx,yy,zz,xy,yz,xz
  TinvT(0, 0) = jac(0, 0) * jac(0, 0);
  TinvT(1, 0) = jac(1, 0) * jac(1, 0);
  TinvT(2, 0) = jac(2, 0) * jac(2, 0);
  TinvT(3, 0) = 2 * jac(0, 0) * jac(1, 0);
  TinvT(4, 0) = 2 * jac(1, 0) * jac(2, 0);
  TinvT(5, 0) = 2 * jac(0, 0) * jac(2, 0);

  TinvT(0, 1) = jac(0, 1) * jac(0, 1);
  TinvT(1, 1) = jac(1, 1) * jac(1, 1);
  TinvT(2, 1) = jac(2, 1) * jac(2, 1);
  TinvT(3, 1) = 2 * jac(0, 1) * jac(1, 1);
  TinvT(4, 1) = 2 * jac(1, 1) * jac(2, 1);
  TinvT(5, 1) = 2 * jac(0, 1) * jac(2, 1);

  TinvT(0, 2) = jac(0, 2) * jac(0, 2);
  TinvT(1, 2) = jac(1, 2) * jac(1, 2);
  TinvT(2, 2) = jac(2, 2) * jac(2, 2);
  TinvT(3, 2) = 2 * jac(0, 2) * jac(1, 2);
  TinvT(4, 2) = 2 * jac(1, 2) * jac(2, 2);
  TinvT(5, 2) = 2 * jac(0, 2) * jac(2, 2);

  TinvT(0, 3) = jac(0, 0) * jac(0, 1);
  TinvT(1, 3) = jac(1, 0) * jac(1, 1);
  TinvT(2, 3) = jac(2, 0) * jac(2, 1);
  TinvT(3, 3) = jac(0, 0) * jac(1, 1) + jac(1, 0) * jac(0, 1);
  TinvT(4, 3) = jac(1, 0) * jac(2, 1) + jac(2, 0) * jac(1, 1);
  TinvT(5, 3) = jac(0, 0) * jac(2, 1) + jac(2, 0) * jac(0, 1);


  TinvT(0, 4) = jac(0, 1) * jac(0, 2);
  TinvT(1, 4) = jac(1, 1) * jac(1, 2);
  TinvT(2, 4) = jac(2, 1) * jac(2, 2);
  TinvT(3, 4) = jac(0, 1) * jac(1, 2) + jac(1, 1) * jac(0, 2);
  TinvT(4, 4) = jac(1, 1) * jac(2, 2) + jac(2, 1) * jac(1, 2);
  TinvT(5, 4) = jac(0, 1) * jac(2, 2) + jac(2, 1) * jac(0, 2);

  TinvT(0, 5) = jac(0, 0) * jac(0, 2);
  TinvT(1, 5) = jac(1, 0) * jac(1, 2);
  TinvT(2, 5) = jac(2, 0) * jac(2, 2);
  TinvT(3, 5) = jac(0, 0) * jac(1, 2) + jac(1, 0) * jac(0, 2);
  TinvT(4, 5) = jac(1, 0) * jac(2, 2) + jac(2, 0) * jac(1, 2);
  TinvT(5, 5) = jac(0, 0) * jac(2, 2) + jac(2, 0) * jac(0, 2);

  // now evaluate T^{-T} with solver
  LINALG::FixedSizeSerialDenseSolver<numstr_, numstr_, 1> solve_for_inverseT;
  solve_for_inverseT.SetMatrix(TinvT);
  int err2 = solve_for_inverseT.Factor();
  int err = solve_for_inverseT.Invert();
  if ((err != 0) && (err2 != 0)) dserror("Inversion of Tinv (Jacobian) failed");
  return;
}


/*----------------------------------------------------------------------*
 |  setup of constant ANS data (private)                       maf 05/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8Plast::Anssetup(
    const LINALG::Matrix<nen_, nsd_>& xrefe,               // material element coords
    const LINALG::Matrix<nen_, nsd_>& xcurr,               // current element coords
    std::vector<LINALG::Matrix<nsd_, nen_>>** deriv_sp,    // derivs eval. at all sampling points
    std::vector<LINALG::Matrix<nsd_, nsd_>>& jac_sps,      // jac at all sampling points
    std::vector<LINALG::Matrix<nsd_, nsd_>>& jac_cur_sps,  // current jac at all sampling points
    LINALG::Matrix<num_ans * num_sp, numdofperelement_>& B_ans_loc)  // modified B
{
  // static matrix object of derivs at sampling points, kept in memory
  static std::vector<LINALG::Matrix<nsd_, nen_>> df_sp(num_sp);
  static bool dfsp_eval;  // flag for re-evaluate everything

  if (dfsp_eval != 0)
  {                      // if true f,df already evaluated
    *deriv_sp = &df_sp;  // return adress of static object to target of pointer
  }
  else
  {
    /*====================================================================*/
    /* 8-node hexhedra Solid-Shell node topology
     * and location of sampling points A to H                             */
    /*--------------------------------------------------------------------*/
    /*                      t
     *                      |
     *             4========|================7
     *          // |        |              //||
     *        //   |        |            //  ||
     *      //     |        |   D      //    ||
     *     5=======E=================6       H
     *    ||       |        |        ||      ||
     *    ||   A   |        o--------||-- C -------s
     *    ||       |       /         ||      ||
     *    F        0----- B ---------G ------3
     *    ||     //     /            ||    //
     *    ||   //     /              ||  //
     *    || //     r                ||//
     *     1=========================2
     *
     */
    /*====================================================================*/
    // (r,s,t) locations of sampling points A,B,C,D,E,F,G,H
    // numsp = 8 here set explicitly to allow direct initializing
    //                A,   B,   C,   D,   E,   F,   G,   H
    double r[8] = {0.0, 1.0, 0.0, -1.0, -1.0, 1.0, 1.0, -1.0};
    double s[8] = {-1.0, 0.0, 1.0, 0.0, -1.0, -1.0, 1.0, 1.0};
    double t[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // fill up df_sp w.r.t. rst directions (NUMDIM) at each sp
    for (int i = 0; i < num_sp; ++i)
    {
      // df wrt to r "+0" for each node(0..7) at each sp [i]
      df_sp[i](0, 0) = -(1.0 - s[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](0, 1) = (1.0 - s[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](0, 2) = (1.0 + s[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](0, 3) = -(1.0 + s[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](0, 4) = -(1.0 - s[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](0, 5) = (1.0 - s[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](0, 6) = (1.0 + s[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](0, 7) = -(1.0 + s[i]) * (1.0 + t[i]) * 0.125;

      // df wrt to s "+1" for each node(0..7) at each sp [i]
      df_sp[i](1, 0) = -(1.0 - r[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](1, 1) = -(1.0 + r[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](1, 2) = (1.0 + r[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](1, 3) = (1.0 - r[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](1, 4) = -(1.0 - r[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](1, 5) = -(1.0 + r[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](1, 6) = (1.0 + r[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](1, 7) = (1.0 - r[i]) * (1.0 + t[i]) * 0.125;

      // df wrt to t "+2" for each node(0..7) at each sp [i]
      df_sp[i](2, 0) = -(1.0 - r[i]) * (1.0 - s[i]) * 0.125;
      df_sp[i](2, 1) = -(1.0 + r[i]) * (1.0 - s[i]) * 0.125;
      df_sp[i](2, 2) = -(1.0 + r[i]) * (1.0 + s[i]) * 0.125;
      df_sp[i](2, 3) = -(1.0 - r[i]) * (1.0 + s[i]) * 0.125;
      df_sp[i](2, 4) = (1.0 - r[i]) * (1.0 - s[i]) * 0.125;
      df_sp[i](2, 5) = (1.0 + r[i]) * (1.0 - s[i]) * 0.125;
      df_sp[i](2, 6) = (1.0 + r[i]) * (1.0 + s[i]) * 0.125;
      df_sp[i](2, 7) = (1.0 - r[i]) * (1.0 + s[i]) * 0.125;
    }

    // return adresses of just evaluated matrices
    *deriv_sp = &df_sp;  // return adress of static object to target of pointer
    dfsp_eval = 1;       // now all arrays are filled statically
  }

  for (int sp = 0; sp < num_sp; ++sp)
  {
    // compute (REFERENCE) Jacobian matrix at all sampling points
    jac_sps[sp].Multiply(df_sp[sp], xrefe);
    // compute CURRENT Jacobian matrix at all sampling points
    jac_cur_sps[sp].Multiply(df_sp[sp], xcurr);
  }

  /*
  ** Compute modified B-operator in local(parametric) space,
  ** evaluated at all sampling points
  */
  // loop over each sampling point
  LINALG::Matrix<nsd_, nsd_> jac_cur;
  for (int sp = 0; sp < num_sp; ++sp)
  {
    /* compute the CURRENT Jacobian matrix at the sampling point:
    **         [ xcurr_,r  ycurr_,r  zcurr_,r ]
    **  Jcur = [ xcurr_,s  ycurr_,s  zcurr_,s ]
    **         [ xcurr_,t  ycurr_,t  zcurr_,t ]
    ** Used to transform the global displacements into parametric space
    */
    jac_cur.Multiply(df_sp[sp], xcurr);

    // fill up B-operator
    for (int inode = 0; inode < nen_; ++inode)
    {
      for (int dim = 0; dim < nsd_; ++dim)
      {
        // modify B_loc_tt = N_t.X_t
        B_ans_loc(sp * num_ans + 0, inode * 3 + dim) = df_sp[sp](2, inode) * jac_cur(2, dim);
        // modify B_loc_st = N_s.X_t + N_t.X_s
        B_ans_loc(sp * num_ans + 1, inode * 3 + dim) =
            df_sp[sp](1, inode) * jac_cur(2, dim) + df_sp[sp](2, inode) * jac_cur(1, dim);
        // modify B_loc_rt = N_r.X_t + N_t.X_r
        B_ans_loc(sp * num_ans + 2, inode * 3 + dim) =
            df_sp[sp](0, inode) * jac_cur(2, dim) + df_sp[sp](2, inode) * jac_cur(0, dim);
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------*
 |                                                          seitz 05/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8Plast::ReInitEas(const DRT::ELEMENTS::So3Plast_EASType EASType)
{
  neas_ = DRT::ELEMENTS::PlastEasTypeToNumEasV(EASType);
  eastype_ = EASType;

  if (eastype_ != soh8p_easnone)
  {
    KaaInv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(neas_, neas_, true));
    Kad_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(neas_, numdofperelement_, true));
    feas_ = Teuchos::rcp(new LINALG::SerialDenseVector(neas_, true));
    alpha_eas_ = Teuchos::rcp(new LINALG::SerialDenseVector(neas_, true));
    alpha_eas_last_timestep_ = Teuchos::rcp(new LINALG::SerialDenseVector(neas_, true));
    alpha_eas_delta_over_last_timestep_ = Teuchos::rcp(new LINALG::SerialDenseVector(neas_, true));
    alpha_eas_inc_ = Teuchos::rcp(new LINALG::SerialDenseVector(neas_, true));
    Kba_ = Teuchos::rcp(new std::vector<LINALG::SerialDenseMatrix>(
        numgpt_, LINALG::SerialDenseMatrix(5, neas_, true)));
  }

  return;
}

/*----------------------------------------------------------------------*
 |                                                          seitz 05/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8Plast::nln_stiffmass(std::vector<double>& disp,  // current displacements
    std::vector<double>& vel,                                              // current velocities
    std::vector<double>& temp,                                             // current temperatures
    LINALG::Matrix<numdofperelement_, numdofperelement_>* stiffmatrix,  // element stiffness matrix
    LINALG::Matrix<numdofperelement_, numdofperelement_>* massmatrix,   // element mass matrix
    LINALG::Matrix<numdofperelement_, 1>* force,      // element internal force vector
    LINALG::Matrix<numgpt_post, numstr_>* elestress,  // stresses at GP
    LINALG::Matrix<numgpt_post, numstr_>* elestrain,  // strains at GP
    Teuchos::ParameterList& params,                   // algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,            // stress output option
    const INPAR::STR::StrainType iostrain             // strain output option
)
{
  InvalidEleData();
  const bool is_tangDis = StrParamsInterface().GetPredictorType() == INPAR::STR::pred_tangdis;

  FillPositionArrays(disp, vel, temp);

  // get plastic hyperelastic material
  MAT::PlasticElastHyper* plmat = NULL;
  if (Material()->MaterialType() == INPAR::MAT::m_plelasthyper)
    plmat = static_cast<MAT::PlasticElastHyper*>(Material().get());
  else
    dserror("so3_ssn_plast elements only with PlasticElastHyper material");

  if (eastype_ != soh8p_easnone)
  {
    EvaluateCenter();
    EasSetup();
  }

  // EAS matrix block
  Epetra_SerialDenseMatrix Kda(numdofperelement_, neas_);

  // ANS modified rows of bop in local(parameter) coords
  LINALG::Matrix<num_ans * num_sp, numdofperelement_> B_ans_loc;
  // Jacobian evaluated at all ANS sampling points
  std::vector<LINALG::Matrix<nsd_, nsd_>> jac_sps(num_sp);
  // CURRENT Jacobian evaluated at all ANS sampling points
  std::vector<LINALG::Matrix<nsd_, nsd_>> jac_cur_sps(num_sp);
  // pointer to derivs evaluated at all sampling points
  std::vector<LINALG::Matrix<nsd_, nen_>>* deriv_sp = NULL;  // derivs eval. at all sampling points
  // evaluate all necessary variables for ANS
  Anssetup(Xrefe(), Xcurr(), &deriv_sp, jac_sps, jac_cur_sps, SetB_ans_loc());

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp = 0; gp < numgpt_; ++gp)
  {
    InvalidGpData();
    // shape functions (shapefunct) and their first derivatives (deriv)
    DRT::UTILS::shape_function<DRT::Element::hex8>(xsi_[gp], SetShapeFunction());
    DRT::UTILS::shape_function_deriv1<DRT::Element::hex8>(xsi_[gp], SetDerivShapeFunction());

    Kinematics(gp);

    /* get the inverse of the Jacobian matrix which looks like:
     **            [ x_,r  y_,r  z_,r ]^-1
     **     J^-1 = [ x_,s  y_,s  z_,s ]
     **            [ x_,t  y_,t  z_,t ]
     */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    SetDerivShapeFunctionXYZ().Multiply(InvJ(), DerivShapeFunction());  // (6.21)

    AnsStrains(gp, jac_sps, jac_cur_sps);
    if (eastype_ != soh8p_easnone)
    {
      EasShape(gp);
      EasEnhanceStrains();
    }

    // strain output *********************************
    OutputStrains(gp, iostrain, elestress);

    // material call *********************************************
    plmat->EvaluateElast(&DefgrdMod(), &DeltaLp(), &SetPK2(), &SetCmat(), gp, Id());
    // material call *********************************************

    // return gp stresses
    OutputStress(gp, iostress, elestress);

    // integrate usual internal force and stiffness matrix
    double detJ_w = DetJ() * wgt_[gp];
    // integrate elastic internal force vector **************************
    // update internal force vector
    if (force != NULL) IntegrateForce(gp, *force);

    // update stiffness matrix
    if (stiffmatrix != NULL)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<numstr_, numdofperelement_> cb;
      cb.Multiply(Cmat(), Bop());
      stiffmatrix->MultiplyTN(detJ_w, Bop(), cb, 1.0);

      // intergrate `geometric' stiffness matrix and add to keu *****************
      // here also the ANS interpolation comes into play
      double r = xsi_.at(gp)(0);
      double s = xsi_.at(gp)(1);
      for (int inod = 0; inod < nen_; ++inod)
      {
        for (int jnod = 0; jnod < nen_; ++jnod)
        {
          LINALG::Matrix<numstr_, 1> G_ij;
          G_ij(0) = DerivShapeFunction()(0, inod) * DerivShapeFunction()(0, jnod);  // rr-dir
          G_ij(1) = DerivShapeFunction()(1, inod) * DerivShapeFunction()(1, jnod);  // ss-dir
          G_ij(3) = DerivShapeFunction()(0, inod) * DerivShapeFunction()(1, jnod) +
                    DerivShapeFunction()(1, inod) * DerivShapeFunction()(0, jnod);  // rs-dir

          // do the ANS related stuff if wanted!
          if (anstype_ == anssosh8_p)
          {
            // ANS modification in tt-dir
            G_ij(2) = 0.25 * (1 - r) * (1 - s) * (*deriv_sp)[4](2, inod) * (*deriv_sp)[4](2, jnod) +
                      0.25 * (1 + r) * (1 - s) * (*deriv_sp)[5](2, inod) * (*deriv_sp)[5](2, jnod) +
                      0.25 * (1 + r) * (1 + s) * (*deriv_sp)[6](2, inod) * (*deriv_sp)[6](2, jnod) +
                      0.25 * (1 - r) * (1 + s) * (*deriv_sp)[7](2, inod) * (*deriv_sp)[7](2, jnod);
            // ANS modification in st-dir
            G_ij(4) = 0.5 * ((1 + r) * ((*deriv_sp)[1](1, inod) * (*deriv_sp)[1](2, jnod) +
                                           (*deriv_sp)[1](2, inod) * (*deriv_sp)[1](1, jnod)) +
                                (1 - r) * ((*deriv_sp)[3](1, inod) * (*deriv_sp)[3](2, jnod) +
                                              (*deriv_sp)[3](2, inod) * (*deriv_sp)[3](1, jnod)));
            // ANS modification in rt-dir
            G_ij(5) = 0.5 * ((1 - s) * ((*deriv_sp)[0](0, inod) * (*deriv_sp)[0](2, jnod) +
                                           (*deriv_sp)[0](2, inod) * (*deriv_sp)[0](0, jnod)) +
                                (1 + s) * ((*deriv_sp)[2](0, inod) * (*deriv_sp)[2](2, jnod) +
                                              (*deriv_sp)[2](2, inod) * (*deriv_sp)[2](0, jnod)));
          }
          else if (anstype_ == ansnone_p)
          {
            G_ij(2) = DerivShapeFunction()(2, inod) * DerivShapeFunction()(2, jnod);  // tt-dir
            G_ij(4) = DerivShapeFunction()(2, inod) * DerivShapeFunction()(1, jnod) +
                      DerivShapeFunction()(1, inod) * DerivShapeFunction()(2, jnod);  // st-dir
            G_ij(5) = DerivShapeFunction()(0, inod) * DerivShapeFunction()(2, jnod) +
                      DerivShapeFunction()(2, inod) * DerivShapeFunction()(0, jnod);  // rt-dir
          }
          else
            dserror("Cannot build geometric stiffness matrix on your ANS-choice!");

          // transformation of local(parameter) space 'back' to global(material) space
          LINALG::Matrix<MAT::NUM_STRESS_3D, 1> G_ij_glob;
          G_ij_glob.Multiply(TinvT(), G_ij);

          // Scalar Gij results from product of G_ij with stress, scaled with detJ*weights
          const double Gij = detJ_w * PK2().Dot(G_ij_glob);

          // add "geometric part" Gij times detJ*weights to stiffness matrix
          (*stiffmatrix)(nsd_ * inod + 0, nsd_ * jnod + 0) += Gij;
          (*stiffmatrix)(nsd_ * inod + 1, nsd_ * jnod + 1) += Gij;
          (*stiffmatrix)(nsd_ * inod + 2, nsd_ * jnod + 2) += Gij;
        }
      }  // end of intergrate `geometric' stiffness ******************************

      // EAS technology: integrate matrices --------------------------------- EAS
      if (eastype_ != soh8p_easnone)
      {
        // integrate Kaa: Kaa += (M^T . cmat . M) * detJ * w(gp)
        // integrate Kda: Kad += (M^T . cmat . B) * detJ * w(gp)
        // integrate feas: feas += (M^T . sigma) * detJ *wp(gp)
        LINALG::SerialDenseMatrix cM(numstr_, neas_, true);  // temporary c . M
        switch (eastype_)
        {
          case soh8p_eassosh8:
            LINALG::DENSEFUNCTIONS::multiply<double, numstr_, numstr_,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas>(
                cM.A(), Cmat().A(), M_eas().A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, numstr_,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas>(
                1.0, *KaaInv_, detJ_w, M_eas(), cM);
            LINALG::DENSEFUNCTIONS::multiplyTN<double,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, numstr_,
                numdofperelement_>(1.0, Kad_->A(), detJ_w, M_eas().A(), cb.A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double, numdofperelement_, numstr_,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas>(
                1.0, Kda.A(), detJ_w, cb.A(), M_eas().A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, numstr_, 1>(
                1.0, feas_->A(), detJ_w, M_eas().A(), PK2().A());
            break;
          case soh8p_easmild:
            LINALG::DENSEFUNCTIONS::multiply<double, numstr_, numstr_,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas>(
                cM.A(), Cmat().A(), M_eas().A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, numstr_,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas>(
                1.0, *KaaInv_, detJ_w, M_eas(), cM);
            LINALG::DENSEFUNCTIONS::multiplyTN<double,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, numstr_,
                numdofperelement_>(1.0, Kad_->A(), detJ_w, M_eas().A(), cb.A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double, numdofperelement_, numstr_,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas>(
                1.0, Kda.A(), detJ_w, cb.A(), M_eas().A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, numstr_, 1>(
                1.0, feas_->A(), detJ_w, M_eas().A(), PK2().A());
            break;
          case soh8p_easnone:
            break;
          default:
            dserror("Don't know what to do with EAS type %d", eastype_);
            break;
        }
      }  // ---------------------------------------------------------------- EAS
    }    // end of stiffness matrix

    if (massmatrix != NULL)  // evaluate mass matrix +++++++++++++++++++++++++
      IntegrateMassMatrix(gp, *massmatrix);

    // plastic modifications
    if ((stiffmatrix != NULL || force != NULL))
    {
      if (HavePlasticSpin())
      {
        if (eastype_ != soh8p_easnone)
          CondensePlasticity<plspin>(DefgrdMod(), DeltaLp(), Bop(), &DerivShapeFunctionXYZ(), NULL,
              detJ_w, gp, 0, params, force, stiffmatrix, &M_eas(), &Kda);
        else
          CondensePlasticity<plspin>(DefgrdMod(), DeltaLp(), Bop(), &DerivShapeFunctionXYZ(), NULL,
              detJ_w, gp, 0, params, force, stiffmatrix);
      }
      else
      {
        if (eastype_ != soh8p_easnone)
          CondensePlasticity<zerospin>(DefgrdMod(), DeltaLp(), Bop(), &DerivShapeFunctionXYZ(),
              NULL, detJ_w, gp, 0, params, force, stiffmatrix, &M_eas(), &Kda);
        else
          CondensePlasticity<zerospin>(DefgrdMod(), DeltaLp(), Bop(), &DerivShapeFunctionXYZ(),
              NULL, detJ_w, gp, 0, params, force, stiffmatrix);
      }
    }  // plastic modifications
  }    // gp loop

  // Static condensation EAS --> stiff ********************************
  if (stiffmatrix != NULL && !is_tangDis && eastype_ != soh8p_easnone)
  {
    Epetra_SerialDenseSolver solve_for_inverseKaa;
    solve_for_inverseKaa.SetMatrix(*KaaInv_);
    solve_for_inverseKaa.Invert();

    Epetra_SerialDenseMatrix kdakaai(numdofperelement_, neas_);
    switch (eastype_)
    {
      case soh8p_eassosh8:
        LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas>(
            0., kdakaai.A(), 1., Kda.A(), KaaInv_->A());
        if (stiffmatrix != NULL)
          LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, numdofperelement_>(
              1., stiffmatrix->A(), -1., kdakaai.A(), Kad_->A());
        if (force != NULL)
          LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, 1>(
              1., force->A(), -1., kdakaai.A(), feas_->A());
        break;
      case soh8p_easmild:
        LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas>(
            0., kdakaai.A(), 1., Kda.A(), KaaInv_->A());
        if (stiffmatrix != NULL)
          LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, numdofperelement_>(
              1., stiffmatrix->A(), -1., kdakaai.A(), Kad_->A());
        if (force != NULL)
          LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, 1>(
              1., force->A(), -1., kdakaai.A(), feas_->A());
        break;
      case soh8p_easnone:
        break;
      default:
        dserror("Don't know what to do with EAS type %d", eastype_);
        break;
    }
  }


  return;
}


/*----------------------------------------------------------------------*
 |                                                          seitz 05/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8Plast::CalculateBop(LINALG::Matrix<numstr_, numdofperelement_>* bop,
    const LINALG::Matrix<nsd_, nsd_>* defgrd, const LINALG::Matrix<nsd_, nen_>* N_XYZ, const int gp)
{
  SetJac_refe().Multiply(DerivShapeFunction(), Xrefe());
  SetJac_curr().Multiply(DerivShapeFunction(), Xcurr());

  if (gp < 0 || gp > 7) dserror("invalid gp number");

  // set up B-Operator in local(parameter) element space including ANS
  LINALG::Matrix<numstr_, numdofperelement_> bop_loc;
  for (int inode = 0; inode < NUMNOD_SOH8; ++inode)
  {
    for (int dim = 0; dim < NUMDIM_SOH8; ++dim)
    {
      // B_loc_rr = N_r.X_r
      bop_loc(0, inode * 3 + dim) = DerivShapeFunction()(0, inode) * Jac_curr()(0, dim);
      // B_loc_ss = N_s.X_s
      bop_loc(1, inode * 3 + dim) = DerivShapeFunction()(1, inode) * Jac_curr()(1, dim);
      // B_loc_rs = N_r.X_s + N_s.X_r
      bop_loc(3, inode * 3 + dim) = DerivShapeFunction()(0, inode) * Jac_curr()(1, dim) +
                                    DerivShapeFunction()(1, inode) * Jac_curr()(0, dim);

      // do the ANS related stuff
      if (anstype_ == anssosh8_p)
      {
        double r = xsi_.at(gp)(0);
        double s = xsi_.at(gp)(1);
        // B_loc_tt = interpolation along (r x s) of ANS B_loc_tt
        //          = (1-r)(1-s)/4 * B_ans(SP E) + (1+r)(1-s)/4 * B_ans(SP F)
        //           +(1+r)(1+s)/4 * B_ans(SP G) + (1-r)(1+s)/4 * B_ans(SP H)
        bop_loc(2, inode * 3 + dim) =
            0.25 * (1 - r) * (1 - s) * B_ans_loc()(0 + 4 * num_ans, inode * 3 + dim)     // E
            + 0.25 * (1 + r) * (1 - s) * B_ans_loc()(0 + 5 * num_ans, inode * 3 + dim)   // F
            + 0.25 * (1 + r) * (1 + s) * B_ans_loc()(0 + 6 * num_ans, inode * 3 + dim)   // G
            + 0.25 * (1 - r) * (1 + s) * B_ans_loc()(0 + 7 * num_ans, inode * 3 + dim);  // H
        // B_loc_st = interpolation along r of ANS B_loc_st
        //          = (1+r)/2 * B_ans(SP B) + (1-r)/2 * B_ans(SP D)
        bop_loc(4, inode * 3 + dim) =
            0.5 * (1.0 + r) * B_ans_loc()(1 + 1 * num_ans, inode * 3 + dim)     // B
            + 0.5 * (1.0 - r) * B_ans_loc()(1 + 3 * num_ans, inode * 3 + dim);  // D

        // B_loc_rt = interpolation along s of ANS B_loc_rt
        //          = (1-s)/2 * B_ans(SP A) + (1+s)/2 * B_ans(SP C)
        bop_loc(5, inode * 3 + dim) =
            0.5 * (1.0 - s) * B_ans_loc()(2 + 0 * num_ans, inode * 3 + dim)     // A
            + 0.5 * (1.0 + s) * B_ans_loc()(2 + 2 * num_ans, inode * 3 + dim);  // C
      }
      else if (anstype_ == ansnone_p)
      {
        // B_loc_tt = N_t.X_t
        bop_loc(2, inode * 3 + dim) = DerivShapeFunction()(2, inode) * Jac_curr()(2, dim);
        // B_loc_st = N_t.X_s + N_s.X_t
        bop_loc(4, inode * 3 + dim) = DerivShapeFunction()(2, inode) * Jac_curr()(1, dim) +
                                      DerivShapeFunction()(1, inode) * Jac_curr()(2, dim);

        // B_loc_rt = N_r.X_t + N_t.X_r
        bop_loc(5, inode * 3 + dim) = DerivShapeFunction()(0, inode) * Jac_curr()(2, dim) +
                                      DerivShapeFunction()(2, inode) * Jac_curr()(0, dim);
      }
      else
        dserror("Cannot build bop_loc based on your ANS-choice!");
    }
  }

  // transformation from local (parameter) element space to global(material) space
  // with famous 'T'-matrix already used for EAS but now evaluated at each gp
  EvaluateT(Jac_refe(), SetTinvT());
  SetBop().Multiply(TinvT(), bop_loc);
}


void DRT::ELEMENTS::So_sh8Plast::AnsStrains(
    const int gp, std::vector<LINALG::Matrix<nsd_, nsd_>>& jac_sps,  // jac at all sampling points
    std::vector<LINALG::Matrix<nsd_, nsd_>>& jac_cur_sps  // current jac at all sampling points
)
{
  // local GL strain vector lstrain={E11,E22,E33,2*E12,2*E23,2*E31}
  LINALG::Matrix<numstr_, 1> lstrain;
  // evaluate glstrains in local(parameter) coords
  // Err = 0.5 * (dx/dr * dx/dr^T - dX/dr * dX/dr^T)
  lstrain(0) =
      0.5 * (+(Jac_curr()(0, 0) * Jac_curr()(0, 0) + Jac_curr()(0, 1) * Jac_curr()(0, 1) +
                 Jac_curr()(0, 2) * Jac_curr()(0, 2)) -
                (Jac_refe()(0, 0) * Jac_refe()(0, 0) + Jac_refe()(0, 1) * Jac_refe()(0, 1) +
                    Jac_refe()(0, 2) * Jac_refe()(0, 2)));
  // Ess = 0.5 * (dy/ds * dy/ds^T - dY/ds * dY/ds^T)
  lstrain(1) =
      0.5 * (+(Jac_curr()(1, 0) * Jac_curr()(1, 0) + Jac_curr()(1, 1) * Jac_curr()(1, 1) +
                 Jac_curr()(1, 2) * Jac_curr()(1, 2)) -
                (Jac_refe()(1, 0) * Jac_refe()(1, 0) + Jac_refe()(1, 1) * Jac_refe()(1, 1) +
                    Jac_refe()(1, 2) * Jac_refe()(1, 2)));
  // Ers = (dx/ds * dy/dr^T - dX/ds * dY/dr^T)
  lstrain(3) = (+(Jac_curr()(0, 0) * Jac_curr()(1, 0) + Jac_curr()(0, 1) * Jac_curr()(1, 1) +
                    Jac_curr()(0, 2) * Jac_curr()(1, 2)) -
                (Jac_refe()(0, 0) * Jac_refe()(1, 0) + Jac_refe()(0, 1) * Jac_refe()(1, 1) +
                    Jac_refe()(0, 2) * Jac_refe()(1, 2)));


  // do the ANS related stuff if wanted!
  if (anstype_ == anssosh8_p)
  {
    double r = xsi_.at(gp)(0);
    double s = xsi_.at(gp)(1);
    // ANS modification of strains ************************************** ANS
    double dxdt_A = 0.0;
    double dXdt_A = 0.0;
    double dydt_B = 0.0;
    double dYdt_B = 0.0;
    double dxdt_C = 0.0;
    double dXdt_C = 0.0;
    double dydt_D = 0.0;
    double dYdt_D = 0.0;

    double dzdt_E = 0.0;
    double dZdt_E = 0.0;
    double dzdt_F = 0.0;
    double dZdt_F = 0.0;
    double dzdt_G = 0.0;
    double dZdt_G = 0.0;
    double dzdt_H = 0.0;
    double dZdt_H = 0.0;

    // vector product of rows of jacobians at corresponding sampling point    std::cout <<
    // jac_cur_sps;
    for (int dim = 0; dim < nsd_; ++dim)
    {
      dxdt_A += jac_cur_sps[0](0, dim) * jac_cur_sps[0](2, dim);  // g_13^A
      dXdt_A += jac_sps[0](0, dim) * jac_sps[0](2, dim);          // G_13^A
      dydt_B += jac_cur_sps[1](1, dim) * jac_cur_sps[1](2, dim);  // g_23^B
      dYdt_B += jac_sps[1](1, dim) * jac_sps[1](2, dim);          // G_23^B
      dxdt_C += jac_cur_sps[2](0, dim) * jac_cur_sps[2](2, dim);  // g_13^C
      dXdt_C += jac_sps[2](0, dim) * jac_sps[2](2, dim);          // G_13^C
      dydt_D += jac_cur_sps[3](1, dim) * jac_cur_sps[3](2, dim);  // g_23^D
      dYdt_D += jac_sps[3](1, dim) * jac_sps[3](2, dim);          // G_23^D

      dzdt_E += jac_cur_sps[4](2, dim) * jac_cur_sps[4](2, dim);
      dZdt_E += jac_sps[4](2, dim) * jac_sps[4](2, dim);
      dzdt_F += jac_cur_sps[5](2, dim) * jac_cur_sps[5](2, dim);
      dZdt_F += jac_sps[5](2, dim) * jac_sps[5](2, dim);
      dzdt_G += jac_cur_sps[6](2, dim) * jac_cur_sps[6](2, dim);
      dZdt_G += jac_sps[6](2, dim) * jac_sps[6](2, dim);
      dzdt_H += jac_cur_sps[7](2, dim) * jac_cur_sps[7](2, dim);
      dZdt_H += jac_sps[7](2, dim) * jac_sps[7](2, dim);
    }
    // E33: remedy of curvature thickness locking
    // Ett = 0.5* ( (1-r)(1-s)/4 * Ett(SP E) + ... + (1-r)(1+s)/4 * Ett(SP H) )
    lstrain(2) = 0.5 * (0.25 * (1 - r) * (1 - s) * (dzdt_E - dZdt_E) +
                           0.25 * (1 + r) * (1 - s) * (dzdt_F - dZdt_F) +
                           0.25 * (1 + r) * (1 + s) * (dzdt_G - dZdt_G) +
                           0.25 * (1 - r) * (1 + s) * (dzdt_H - dZdt_H));
    // E23: remedy of transverse shear locking
    // Est = (1+r)/2 * Est(SP B) + (1-r)/2 * Est(SP D)
    lstrain(4) = 0.5 * (1 + r) * (dydt_B - dYdt_B) + 0.5 * (1 - r) * (dydt_D - dYdt_D);
    // E13: remedy of transverse shear locking
    // Ert = (1-s)/2 * Ert(SP A) + (1+s)/2 * Ert(SP C)
    lstrain(5) = 0.5 * (1 - s) * (dxdt_A - dXdt_A) + 0.5 * (1 + s) * (dxdt_C - dXdt_C);
    // ANS modification of strains ************************************** ANS
  }
  else if (anstype_ == ansnone_p)
  {
    // No ANS!
    // Ett = 0.5 * (dz/dt * dz/dt^T - dZ/dt * dZ/dt^T)
    lstrain(2) =
        0.5 * (+(Jac_curr()(2, 0) * Jac_curr()(2, 0) + Jac_curr()(2, 1) * Jac_curr()(2, 1) +
                   Jac_curr()(2, 2) * Jac_curr()(2, 2)) -
                  (Jac_refe()(2, 0) * Jac_refe()(2, 0) + Jac_refe()(2, 1) * Jac_refe()(2, 1) +
                      Jac_refe()(2, 2) * Jac_refe()(2, 2)));
    // Est = (dz/ds * dy/dt^T - dZ/ds * dY/dt^T)
    lstrain(4) = (+(Jac_curr()(2, 0) * Jac_curr()(1, 0) + Jac_curr()(2, 1) * Jac_curr()(1, 1) +
                      Jac_curr()(2, 2) * Jac_curr()(1, 2)) -
                  (Jac_refe()(2, 0) * Jac_refe()(1, 0) + Jac_refe()(2, 1) * Jac_refe()(1, 1) +
                      Jac_refe()(2, 2) * Jac_refe()(1, 2)));
    // Est = (dz/dr * dx/dt^T - dZ/dr * dX/dt^T)
    lstrain(5) = (+(Jac_curr()(2, 0) * Jac_curr()(0, 0) + Jac_curr()(2, 1) * Jac_curr()(0, 1) +
                      Jac_curr()(2, 2) * Jac_curr()(0, 2)) -
                  (Jac_refe()(2, 0) * Jac_refe()(0, 0) + Jac_refe()(2, 1) * Jac_refe()(0, 1) +
                      Jac_refe()(2, 2) * Jac_refe()(0, 2)));
  }
  else
    dserror("Cannot build local strains based on your ANS-choice!");

  // transformation of local glstrains 'back' to global(material) space
  static LINALG::Matrix<numstr_, 1> glstrain(false);
  glstrain.Multiply(TinvT(), lstrain);

  for (int i = 0; i < nsd_; ++i) SetRCG()(i, i) = 2. * glstrain(i) + 1.;
  SetRCG()(0, 1) = SetRCG()(1, 0) = glstrain(3);
  SetRCG()(2, 1) = SetRCG()(1, 2) = glstrain(4);
  SetRCG()(0, 2) = SetRCG()(2, 0) = glstrain(5);

  // calculate deformation gradient consistent with modified GL strain tensor
  CalcConsistentDefgrd();
}
