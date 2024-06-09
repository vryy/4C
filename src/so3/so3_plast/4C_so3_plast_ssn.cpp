/*----------------------------------------------------------------------*/
/*! \file
\brief
\level 2

*/


/*----------------------------------------------------------------------*
 | headers                                                  seitz 07/13 |
 *----------------------------------------------------------------------*/
#include "4C_so3_plast_ssn.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_tsi.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_mat_plasticelasthyper.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_surface.hpp"
#include "4C_thermo_ele_impl_utils.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | ctor (public)                                            seitz 07/13 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::So3Plast<distype>::So3Plast(int id, int owner)
    : SoBase(id, owner),
      fbar_(false),
      KbbInv_(std::vector<Core::LinAlg::SerialDenseMatrix>(0)),
      Kbd_(std::vector<Core::LinAlg::SerialDenseMatrix>(0)),
      fbeta_(std::vector<Core::LinAlg::SerialDenseVector>(0)),
      dDp_last_iter_(std::vector<Core::LinAlg::SerialDenseVector>(0)),
      dDp_inc_(std::vector<Core::LinAlg::SerialDenseVector>(0)),
      plspintype_(plspin),
      KaaInv_(Teuchos::null),
      Kad_(Teuchos::null),
      KaT_(Teuchos::null),
      KdT_eas_(Teuchos::null),
      feas_(Teuchos::null),
      Kba_(Teuchos::null),
      alpha_eas_(Teuchos::null),
      alpha_eas_last_timestep_(Teuchos::null),
      alpha_eas_delta_over_last_timestep_(Teuchos::null),
      alpha_eas_inc_(Teuchos::null),
      eastype_(soh8p_easnone),
      neas_(0),
      tsi_(false),
      is_nitsche_contact_(false)
{
  if (distype == Core::FE::CellType::nurbs27)
    SetNurbsElement() = true;
  else
    SetNurbsElement() = false;
  return;
}


/*----------------------------------------------------------------------*
 | copy-ctor (public)                                       seitz 07/13 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::So3Plast<distype>::So3Plast(const Discret::ELEMENTS::So3Plast<distype>& old)
    : SoBase(old)
{
  if (distype == Core::FE::CellType::nurbs27)
    SetNurbsElement() = true;
  else
    SetNurbsElement() = false;
  return;
}


/*----------------------------------------------------------------------*
 | deep copy this instance of Solid3 and return pointer to  seitz 07/13 |
 | it (public)                                                          |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Core::Elements::Element* Discret::ELEMENTS::So3Plast<distype>::Clone() const
{
  auto* newelement = new Discret::ELEMENTS::So3Plast<distype>(*this);

  return newelement;
}


template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::nen_, 1>>
    Discret::ELEMENTS::So3Plast<distype>::shapefunct_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::nsd_,
                    Discret::ELEMENTS::So3Plast<distype>::nen_>>
    Discret::ELEMENTS::So3Plast<distype>::deriv_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::nsd_,
                    Discret::ELEMENTS::So3Plast<distype>::nsd_>>
    Discret::ELEMENTS::So3Plast<distype>::invJ_;
template <Core::FE::CellType distype>
std::pair<bool, double> Discret::ELEMENTS::So3Plast<distype>::detJ_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::nsd_,
                    Discret::ELEMENTS::So3Plast<distype>::nen_>>
    Discret::ELEMENTS::So3Plast<distype>::N_XYZ_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::nsd_,
                    Discret::ELEMENTS::So3Plast<distype>::nsd_>>
    Discret::ELEMENTS::So3Plast<distype>::defgrd_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::nsd_,
                    Discret::ELEMENTS::So3Plast<distype>::nsd_>>
    Discret::ELEMENTS::So3Plast<distype>::defgrd_mod_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::nsd_,
                    Discret::ELEMENTS::So3Plast<distype>::nsd_>>
    Discret::ELEMENTS::So3Plast<distype>::rcg_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::nsd_,
                    Discret::ELEMENTS::So3Plast<distype>::nsd_>>
    Discret::ELEMENTS::So3Plast<distype>::delta_Lp_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::numstr_,
                    Discret::ELEMENTS::So3Plast<distype>::numdofperelement_>>
    Discret::ELEMENTS::So3Plast<distype>::bop_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::numstr_, 1>>
    Discret::ELEMENTS::So3Plast<distype>::pk2_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::numstr_,
                    Discret::ELEMENTS::So3Plast<distype>::numstr_>>
    Discret::ELEMENTS::So3Plast<distype>::cmat_;

template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::nen_,
                    Discret::ELEMENTS::So3Plast<distype>::nsd_>>
    Discret::ELEMENTS::So3Plast<distype>::xrefe_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::nen_,
                    Discret::ELEMENTS::So3Plast<distype>::nsd_>>
    Discret::ELEMENTS::So3Plast<distype>::xcurr_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::nen_,
                    Discret::ELEMENTS::So3Plast<distype>::nsd_>>
    Discret::ELEMENTS::So3Plast<distype>::xcurr_rate_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::nen_, 1>>
    Discret::ELEMENTS::So3Plast<distype>::etemp_;

template <Core::FE::CellType distype>
std::pair<bool, double> Discret::ELEMENTS::So3Plast<distype>::detF_;
template <Core::FE::CellType distype>
std::pair<bool, double> Discret::ELEMENTS::So3Plast<distype>::detF_0_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::nsd_,
                    Discret::ELEMENTS::So3Plast<distype>::nsd_>>
    Discret::ELEMENTS::So3Plast<distype>::inv_defgrd_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::nsd_,
                    Discret::ELEMENTS::So3Plast<distype>::nsd_>>
    Discret::ELEMENTS::So3Plast<distype>::inv_defgrd_0_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::nsd_,
                    Discret::ELEMENTS::So3Plast<distype>::nen_>>
    Discret::ELEMENTS::So3Plast<distype>::N_XYZ_0_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::numstr_, 1>>
    Discret::ELEMENTS::So3Plast<distype>::rcg_vec_;
template <Core::FE::CellType distype>
std::pair<bool, double> Discret::ELEMENTS::So3Plast<distype>::f_bar_fac_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::numdofperelement_, 1>>
    Discret::ELEMENTS::So3Plast<distype>::htensor_;

template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::numstr_,
                    Discret::ELEMENTS::So3Plast<distype>::numstr_>>
    Discret::ELEMENTS::So3Plast<distype>::T0invT_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::nsd_,
                    Discret::ELEMENTS::So3Plast<distype>::nsd_>>
    Discret::ELEMENTS::So3Plast<distype>::jac_0_;
template <Core::FE::CellType distype>
std::pair<bool, double> Discret::ELEMENTS::So3Plast<distype>::det_jac_0_;
template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::SerialDenseMatrix> Discret::ELEMENTS::So3Plast<distype>::M_eas_;

template <Core::FE::CellType distype>
std::pair<bool, Core::LinAlg::Matrix<Discret::ELEMENTS::So3Plast<distype>::nen_, 1>>
    Discret::ELEMENTS::So3Plast<distype>::weights_;
template <Core::FE::CellType distype>
std::pair<bool, std::vector<Core::LinAlg::SerialDenseVector>>
    Discret::ELEMENTS::So3Plast<distype>::knots_;

/*----------------------------------------------------------------------*
 |                                                          seitz 05/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::So3Plast<distype>::NumVolume() const
{
  switch (distype)
  {
    case Core::FE::CellType::tet4:
    case Core::FE::CellType::hex8:
    case Core::FE::CellType::hex18:
    case Core::FE::CellType::hex27:
    case Core::FE::CellType::nurbs27:
      return 0;
      break;
    default:
      FOUR_C_THROW("unknown distpye for So3Plast");
      break;
      return 0;
  }
}

/*----------------------------------------------------------------------*
 |                                                          seitz 05/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::So3Plast<distype>::NumSurface() const
{
  switch (distype)
  {
    case Core::FE::CellType::hex8:
    case Core::FE::CellType::hex18:
    case Core::FE::CellType::hex27:
    case Core::FE::CellType::nurbs27:
      return 6;
      break;
    case Core::FE::CellType::tet4:
      return 4;
      break;
    default:
      FOUR_C_THROW("unknown distpye for So3Plast");
      break;
      return 0;
  }
}

/*----------------------------------------------------------------------*
 |                                                          seitz 05/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::So3Plast<distype>::NumLine() const
{
  switch (distype)
  {
    case Core::FE::CellType::hex8:
    case Core::FE::CellType::hex18:
    case Core::FE::CellType::hex27:
    case Core::FE::CellType::nurbs27:
      return 12;
      break;
    case Core::FE::CellType::tet4:
      return 6;
      break;
    default:
      FOUR_C_THROW("unknown distpye for So3Plast");
      break;
      return 0;
  }
}

/*----------------------------------------------------------------------*
 |                                                          seitz 05/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::So3Plast<distype>::Lines()
{
  return Core::Communication::ElementBoundaryFactory<StructuralLine, Core::Elements::Element>(
      Core::Communication::buildLines, *this);
}

/*----------------------------------------------------------------------*
 |                                                          seitz 05/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::So3Plast<distype>::Surfaces()
{
  return Core::Communication::ElementBoundaryFactory<StructuralSurface, Core::Elements::Element>(
      Core::Communication::buildSurfaces, *this);
}

/*----------------------------------------------------------------------*
 | pack data (public)                                       seitz 07/13 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // add base class Element
  SoBase::Pack(data);

  // Gauss points and weights
  const auto size2 = (int)xsi_.size();
  add_to_pack(data, size2);
  for (int i = 0; i < size2; ++i) add_to_pack(data, xsi_[i]);
  add_to_pack(data, wgt_);

  // parameters
  add_to_pack(data, (int)fbar_);

  // plastic spin type
  add_to_pack(data, (int)plspintype_);

  // tsi
  add_to_pack(data, (int)tsi_);
  if (tsi_)
  {
    add_to_pack(data, (int)KbT_->size());
    for (unsigned i = 0; i < KbT_->size(); i++)
    {
      add_to_pack(data, (*dFintdT_)[i]);
      add_to_pack(data, (*KbT_)[i]);
      add_to_pack(data, (*temp_last_)[i]);
    }
  }

  // EAS element technology
  add_to_pack(data, (int)eastype_);
  add_to_pack(data, neas_);
  if (eastype_ != soh8p_easnone)
  {
    add_to_pack(data, (*alpha_eas_));
    add_to_pack(data, (*alpha_eas_last_timestep_));
    add_to_pack(data, (*alpha_eas_delta_over_last_timestep_));
  }

  // history at each Gauss point
  int histsize = dDp_last_iter_.size();
  add_to_pack(data, histsize);
  if (histsize != 0)
    for (int i = 0; i < histsize; i++) add_to_pack(data, dDp_last_iter_[i]);

  // nitsche contact
  add_to_pack(data, (int)is_nitsche_contact_);
  if (is_nitsche_contact_)
  {
    add_to_pack(data, cauchy_);
    add_to_pack(data, cauchy_deriv_);
    if (tsi_) add_to_pack(data, cauchy_deriv_T_);
  }

  return;
}  // Pack()


/*----------------------------------------------------------------------*
 | unpack data (public)                                     seitz 07/13 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  SoBase::Unpack(basedata);

  // Gauss points and weights
  int size2 = ExtractInt(position, data);
  xsi_.resize(size2, Core::LinAlg::Matrix<nsd_, 1>(true));
  for (int i = 0; i < size2; ++i) extract_from_pack(position, data, xsi_[i]);
  extract_from_pack(position, data, wgt_);
  numgpt_ = wgt_.size();

  // paramters
  fbar_ = (bool)ExtractInt(position, data);

  // plastic spin type
  plspintype_ = static_cast<PlSpinType>(ExtractInt(position, data));

  // tsi
  tsi_ = (bool)ExtractInt(position, data);
  if (tsi_)
  {
    dFintdT_ = Teuchos::rcp(new std::vector<Core::LinAlg::Matrix<numdofperelement_, 1>>(numgpt_));
    KbT_ = Teuchos::rcp(new std::vector<Core::LinAlg::SerialDenseVector>(
        numgpt_, Core::LinAlg::SerialDenseVector(plspintype_, true)));
    temp_last_ = Teuchos::rcp(new std::vector<double>(numgpt_));
    int size = ExtractInt(position, data);
    for (int i = 0; i < size; i++)
    {
      extract_from_pack(position, data, (*dFintdT_)[i]);
      extract_from_pack(position, data, (*KbT_)[i]);
      extract_from_pack(position, data, (*temp_last_)[i]);
    }
  }

  // EAS element technology
  eastype_ = static_cast<Discret::ELEMENTS::So3PlastEasType>(ExtractInt(position, data));
  extract_from_pack(position, data, neas_);

  // no EAS
  if (eastype_ == soh8p_easnone)
  {
    KaaInv_ = Teuchos::null;
    Kad_ = Teuchos::null;
    KaT_ = Teuchos::null;
    KdT_eas_ = Teuchos::null;
    feas_ = Teuchos::null;
    Kba_ = Teuchos::null;
    alpha_eas_ = Teuchos::null;
    alpha_eas_last_timestep_ = Teuchos::null;
    alpha_eas_delta_over_last_timestep_ = Teuchos::null;
    alpha_eas_inc_ = Teuchos::null;
  }
  else
  {
    KaaInv_ = Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(neas_, neas_, true));
    Kad_ = Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(neas_, numdofperelement_, true));
    if (tsi_)
    {
      KaT_ = Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(neas_, nen_, true));
      KdT_eas_ = Teuchos::rcp(new Core::LinAlg::Matrix<numdofperelement_, nen_>);
    }
    feas_ = Teuchos::rcp(new Core::LinAlg::SerialDenseVector(neas_, true));
    Kba_ = Teuchos::rcp(new std::vector<Core::LinAlg::SerialDenseMatrix>(
        numgpt_, Core::LinAlg::SerialDenseMatrix(plspintype_, neas_, true)));
    alpha_eas_ = Teuchos::rcp(new Core::LinAlg::SerialDenseVector(neas_, true));
    alpha_eas_last_timestep_ = Teuchos::rcp(new Core::LinAlg::SerialDenseVector(neas_, true));
    alpha_eas_delta_over_last_timestep_ =
        Teuchos::rcp(new Core::LinAlg::SerialDenseVector(neas_, true));
    alpha_eas_inc_ = Teuchos::rcp(new Core::LinAlg::SerialDenseVector(neas_, true));
  }

  KbbInv_.resize(numgpt_, Core::LinAlg::SerialDenseMatrix(plspintype_, plspintype_, true));
  Kbd_.resize(numgpt_, Core::LinAlg::SerialDenseMatrix(plspintype_, numdofperelement_, true));
  fbeta_.resize(numgpt_, Core::LinAlg::SerialDenseVector(plspintype_, true));
  dDp_last_iter_.resize(numgpt_, Core::LinAlg::SerialDenseVector(plspintype_, true));
  dDp_inc_.resize(numgpt_, Core::LinAlg::SerialDenseVector(plspintype_, true));

  if (eastype_ != soh8p_easnone)
  {
    extract_from_pack(position, data, (*alpha_eas_));
    extract_from_pack(position, data, (*alpha_eas_last_timestep_));
    extract_from_pack(position, data, (*alpha_eas_delta_over_last_timestep_));
  }

  int size = ExtractInt(position, data);
  for (int i = 0; i < size; i++) extract_from_pack(position, data, dDp_last_iter_[i]);

  // Nitsche contact stuff
  is_nitsche_contact_ = (bool)ExtractInt(position, data);
  if (is_nitsche_contact_)
  {
    extract_from_pack(position, data, cauchy_);
    extract_from_pack(position, data, cauchy_deriv_);
    if (tsi_)
      extract_from_pack(position, data, cauchy_deriv_T_);
    else
      cauchy_deriv_T_.resize(0);
  }
  else
  {
    cauchy_.resize(0);
    cauchy_deriv_.resize(0);
    cauchy_deriv_T_.resize(0);
  }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;

}  // Unpack()


/*----------------------------------------------------------------------*
 | print this element (public)                              seitz 07/13 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::Print(std::ostream& os) const
{
  os << "So3Plast ";
  return;
}


/*----------------------------------------------------------------------*
 | read this element, get the material (public)             seitz 07/13 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
bool Discret::ELEMENTS::So3Plast<distype>::ReadElement(
    const std::string& eletype, const std::string& eledistype, Input::LineDefinition* linedef)
{
  std::string buffer;
  linedef->ExtractString("KINEM", buffer);

  // geometrically linear
  if (buffer == "linear")
  {
    FOUR_C_THROW("no linear kinematics");
  }
  // geometrically non-linear with Total Lagrangean approach
  else if (buffer == "nonlinear")
  {
    kintype_ = Inpar::STR::KinemType::nonlinearTotLag;
    // everything ok
  }
  else
    FOUR_C_THROW("Reading of SO3_PLAST element failed! KINEM unknown");

  // fbar
  if (linedef->HaveNamed("FBAR"))
  {
    std::string fb;
    linedef->ExtractString("FBAR", fb);
    if (fb == "yes")
      fbar_ = true;
    else if (fb == "no")
      fbar_ = false;
    else
      FOUR_C_THROW("unknown fbar option (valid: yes/no)");
  }

  // quadrature
  if (linedef->HaveNamed("NUMGP"))
  {
    if (distype != Core::FE::CellType::hex8)
      FOUR_C_THROW("You may only choose the Gauss point number for SOLIDH8PLAST");
    if (Global::Problem::Instance()->GetProblemType() == Core::ProblemType::tsi)
      FOUR_C_THROW("You may not choose the Gauss point number in TSI problems");

    int ngp = 0;
    linedef->ExtractInt("NUMGP", ngp);

    switch (ngp)
    {
      case 8:
      {
        Core::FE::IntPointsAndWeights<nsd_> intpoints(Core::FE::GaussRule3D::hex_8point);
        numgpt_ = intpoints.IP().nquad;
        xsi_.resize(numgpt_);
        wgt_.resize(numgpt_);
        for (int gp = 0; gp < numgpt_; ++gp)
        {
          wgt_[gp] = (intpoints.IP().qwgt)[gp];
          const double* gpcoord = (intpoints.IP().qxg)[gp];
          for (int idim = 0; idim < nsd_; idim++) xsi_[gp](idim) = gpcoord[idim];
        }
        break;
      }
      case 9:
      {
        Core::FE::GaussIntegration ip(distype, 3);
        numgpt_ = ip.NumPoints() + 1;
        xsi_.resize(numgpt_);
        wgt_.resize(numgpt_);
        for (int gp = 0; gp < numgpt_ - 1; ++gp)
        {
          wgt_[gp] = 5. / 9.;
          const double* gpcoord = ip.Point(gp);
          for (int idim = 0; idim < nsd_; idim++) xsi_[gp](idim) = gpcoord[idim];
        }
        // 9th quadrature point at element center
        xsi_[numgpt_ - 1](0) = 0.;
        xsi_[numgpt_ - 1](1) = 0.;
        xsi_[numgpt_ - 1](2) = 0.;
        wgt_[numgpt_ - 1] = 32. / 9.;
        break;
      }
      case 27:
      {
        Core::FE::IntPointsAndWeights<nsd_> intpoints(Core::FE::GaussRule3D::hex_27point);
        numgpt_ = intpoints.IP().nquad;
        xsi_.resize(numgpt_);
        wgt_.resize(numgpt_);
        for (int gp = 0; gp < numgpt_; ++gp)
        {
          wgt_[gp] = (intpoints.IP().qwgt)[gp];
          const double* gpcoord = (intpoints.IP().qxg)[gp];
          for (int idim = 0; idim < nsd_; idim++) xsi_[gp](idim) = gpcoord[idim];
        }
        break;
      }
      default:
        FOUR_C_THROW("so3_plast doesn't know what to do with %i Gauss points", ngp);
        break;
    }
  }
  else  // default integration
  {
    Core::FE::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);
    numgpt_ = intpoints.IP().nquad;
    xsi_.resize(numgpt_);
    wgt_.resize(numgpt_);
    for (int gp = 0; gp < numgpt_; ++gp)
    {
      wgt_[gp] = (intpoints.IP().qwgt)[gp];
      const double* gpcoord = (intpoints.IP().qxg)[gp];
      for (int idim = 0; idim < nsd_; idim++) xsi_[gp](idim) = gpcoord[idim];
    }
  }

  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);

  SetMaterial(0, Mat::Factory(material));

  Teuchos::RCP<Mat::So3Material> so3mat = SolidMaterial();
  so3mat->Setup(numgpt_, linedef);
  so3mat->ValidKinematics(Inpar::STR::KinemType::nonlinearTotLag);


  // Validate that materials doesn't use extended update call.
  if (SolidMaterial()->UsesExtendedUpdate())
    FOUR_C_THROW("This element currently does not support the extended update call.");

  if (so3mat->MaterialType() != Core::Materials::m_plelasthyper)
    std::cout << "*** warning *** so3plast used w/o PlasticElastHyper material. Better use "
                 "standard solid element!\n";
  if (have_plastic_spin())
    plspintype_ = plspin;
  else
    plspintype_ = zerospin;

  // EAS
  if (linedef->HaveNamed("EAS"))
  {
    if (distype != Core::FE::CellType::hex8)
      FOUR_C_THROW("EAS in so3 plast currently only for HEX8 elements");

    linedef->ExtractString("EAS", buffer);

    if (buffer == "none")
      eastype_ = soh8p_easnone;
    else if (buffer == "mild")
      eastype_ = soh8p_easmild;
    else if (buffer == "full")
      eastype_ = soh8p_easfull;
    else
      FOUR_C_THROW("unknown EAS type for so3_plast");

    if (fbar_ && eastype_ != soh8p_easnone) FOUR_C_THROW("no combination of Fbar and EAS");
  }
  else
    eastype_ = soh8p_easnone;

  // initialize EAS data
  eas_init();

  // plasticity related stuff
  KbbInv_.resize(numgpt_, Core::LinAlg::SerialDenseMatrix(plspintype_, plspintype_, true));
  Kbd_.resize(numgpt_, Core::LinAlg::SerialDenseMatrix(plspintype_, numdofperelement_, true));
  fbeta_.resize(numgpt_, Core::LinAlg::SerialDenseVector(plspintype_, true));
  dDp_last_iter_.resize(numgpt_, Core::LinAlg::SerialDenseVector(plspintype_, true));
  dDp_inc_.resize(numgpt_, Core::LinAlg::SerialDenseVector(plspintype_, true));

  Teuchos::ParameterList plparams = Global::Problem::Instance()->semi_smooth_plast_params();
  Core::UTILS::AddEnumClassToParameterList(
      "Core::ProblemType", Global::Problem::Instance()->GetProblemType(), plparams);
  ReadParameterList(Teuchos::rcpFromRef<Teuchos::ParameterList>(plparams));


  return true;

}  // ReadElement()

/*----------------------------------------------------------------------*
 | get the nodes from so3 (public)                          seitz 07/13 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::So3Plast<distype>::UniqueParObjectId() const
{
  switch (distype)
  {
    case Core::FE::CellType::hex8:
    {
      return SoHex8PlastType::Instance().UniqueParObjectId();
      break;
    }  // hex8
    case Core::FE::CellType::hex27:
      return SoHex27PlastType::Instance().UniqueParObjectId();
      break;
    case Core::FE::CellType::tet4:
      return SoTet4PlastType::Instance().UniqueParObjectId();
      break;
    case Core::FE::CellType::nurbs27:
      return SoNurbs27PlastType::Instance().UniqueParObjectId();
      break;
    default:
      FOUR_C_THROW("unknown element type!");
      break;
  }
  // Intel compiler needs a return
  return -1;

}  // UniqueParObjectId()


/*----------------------------------------------------------------------*
 | get the nodes from so3 (public)                          seitz 07/13 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Core::Elements::ElementType& Discret::ELEMENTS::So3Plast<distype>::ElementType() const
{
  switch (distype)
  {
    case Core::FE::CellType::hex8:
    {
      return SoHex8PlastType::Instance();
      break;
    }
    case Core::FE::CellType::hex27:
      return SoHex27PlastType::Instance();
      break;
    case Core::FE::CellType::tet4:
      return SoTet4PlastType::Instance();
      break;
    case Core::FE::CellType::nurbs27:
      return SoNurbs27PlastType::Instance();
      break;
    default:
      FOUR_C_THROW("unknown element type!");
      break;
  }
  // Intel compiler needs a return
  return SoHex8PlastType::Instance();

};  // ElementType()


/*----------------------------------------------------------------------*
 | return names of visualization data (public)              seitz 07/13 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::VisNames(std::map<std::string, int>& names)
{
  Core::Elements::Element::VisNames(names);
  SolidMaterial()->VisNames(names);

  return;
}  // VisNames()

/*----------------------------------------------------------------------*
 | return visualization data (public)                       seitz 07/13 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
bool Discret::ELEMENTS::So3Plast<distype>::VisData(
    const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (Core::Elements::Element::VisData(name, data)) return true;

  return SolidMaterial()->VisData(name, data, numgpt_, Id());

}  // VisData()

/*----------------------------------------------------------------------*
 | read relevant parameters from paramter list              seitz 01/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::ReadParameterList(
    Teuchos::RCP<Teuchos::ParameterList> plparams)
{
  double cpl = plparams->get<double>("SEMI_SMOOTH_CPL");
  double s = plparams->get<double>("STABILIZATION_S");
  if (Material()->MaterialType() == Core::Materials::m_plelasthyper)
    static_cast<Mat::PlasticElastHyper*>(Material().get())->GetParams(s, cpl);

  Core::ProblemType probtype =
      Teuchos::getIntegralValue<Core::ProblemType>(*plparams, "Core::ProblemType");
  if (probtype == Core::ProblemType::tsi)
    tsi_ = true;
  else
    tsi_ = false;
  if (tsi_)
  {
    // get plastic hyperelastic material
    Mat::PlasticElastHyper* plmat = nullptr;
    if (Material()->MaterialType() == Core::Materials::m_plelasthyper)
      plmat = static_cast<Mat::PlasticElastHyper*>(Material().get());
    else
      FOUR_C_THROW("so3_ssn_plast elements only with PlasticElastHyper material");

    // get dissipation mode
    auto mode =
        Core::UTILS::IntegralValue<Inpar::TSI::DissipationMode>(*plparams, "DISSIPATION_MODE");

    // prepare material for tsi
    plmat->SetupTSI(numgpt_, numdofperelement_, (eastype_ != soh8p_easnone), mode);

    // setup element data
    dFintdT_ = Teuchos::rcp(new std::vector<Core::LinAlg::Matrix<numdofperelement_, 1>>(numgpt_));
    temp_last_ = Teuchos::rcp(new std::vector<double>(numgpt_, plmat->InitTemp()));
    KbT_ = Teuchos::rcp(new std::vector<Core::LinAlg::SerialDenseVector>(
        numgpt_, Core::LinAlg::SerialDenseVector(plspintype_, true)));

    if (eastype_ != soh8p_easnone)
    {
      KaT_ = Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(neas_, nen_, true));
      KdT_eas_ = Teuchos::rcp(new Core::LinAlg::Matrix<numdofperelement_, nen_>);
    }
    else
    {
      KaT_ = Teuchos::null;
      KdT_eas_ = Teuchos::null;
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
template <unsigned int num_cols>
void Discret::ELEMENTS::So3Plast<distype>::soh8_expol(
    Core::LinAlg::Matrix<numgpt_post, num_cols>& data, Epetra_MultiVector& expolData)
{
  if (distype != Core::FE::CellType::hex8) FOUR_C_THROW("soh8_expol called from non-hex8 element");

  // static variables, that are the same for every element
  static Core::LinAlg::Matrix<nen_, numgpt_post> expolOperator;
  static bool isfilled;

  if (isfilled == false)
  {
    double sq3 = sqrt(3.0);

    expolOperator(0, 0) = 1.25 + 0.75 * sq3;
    expolOperator(0, 1) = -0.25 - 0.25 * sq3;
    expolOperator(0, 2) = -0.25 + 0.25 * sq3;
    expolOperator(0, 3) = -0.25 - 0.25 * sq3;
    expolOperator(0, 4) = -0.25 - 0.25 * sq3;
    expolOperator(0, 5) = -0.25 + 0.25 * sq3;
    expolOperator(0, 6) = 1.25 - 0.75 * sq3;
    expolOperator(0, 7) = -0.25 + 0.25 * sq3;
    expolOperator(1, 1) = 1.25 + 0.75 * sq3;
    expolOperator(1, 2) = -0.25 - 0.25 * sq3;
    expolOperator(1, 3) = -0.25 + 0.25 * sq3;
    expolOperator(1, 4) = -0.25 + 0.25 * sq3;
    expolOperator(1, 5) = -0.25 - 0.25 * sq3;
    expolOperator(1, 6) = -0.25 + 0.25 * sq3;
    expolOperator(1, 7) = 1.25 - 0.75 * sq3;
    expolOperator(2, 2) = 1.25 + 0.75 * sq3;
    expolOperator(2, 3) = -0.25 - 0.25 * sq3;
    expolOperator(2, 4) = 1.25 - 0.75 * sq3;
    expolOperator(2, 5) = -0.25 + 0.25 * sq3;
    expolOperator(2, 6) = -0.25 - 0.25 * sq3;
    expolOperator(2, 7) = -0.25 + 0.25 * sq3;
    expolOperator(3, 3) = 1.25 + 0.75 * sq3;
    expolOperator(3, 4) = -0.25 + 0.25 * sq3;
    expolOperator(3, 5) = 1.25 - 0.75 * sq3;
    expolOperator(3, 6) = -0.25 + 0.25 * sq3;
    expolOperator(3, 7) = -0.25 - 0.25 * sq3;
    expolOperator(4, 4) = 1.25 + 0.75 * sq3;
    expolOperator(4, 5) = -0.25 - 0.25 * sq3;
    expolOperator(4, 6) = -0.25 + 0.25 * sq3;
    expolOperator(4, 7) = -0.25 - 0.25 * sq3;
    expolOperator(5, 5) = 1.25 + 0.75 * sq3;
    expolOperator(5, 6) = -0.25 - 0.25 * sq3;
    expolOperator(5, 7) = -0.25 + 0.25 * sq3;
    expolOperator(6, 6) = 1.25 + 0.75 * sq3;
    expolOperator(6, 7) = -0.25 - 0.25 * sq3;
    expolOperator(7, 7) = 1.25 + 0.75 * sq3;

    for (int i = 0; i < NUMNOD_SOH8; ++i)
    {
      for (int j = 0; j < i; ++j)
      {
        expolOperator(i, j) = expolOperator(j, i);
      }
    }

    isfilled = true;
  }

  Core::LinAlg::Matrix<nen_, num_cols> nodalData;
  nodalData.Multiply(expolOperator, data);

  // "assembly" of extrapolated nodal data
  for (int i = 0; i < nen_; ++i)
  {
    const int lid = expolData.Map().LID(NodeIds()[i]);
    if (lid >= 0)  // rownode
    {
      const double invmyadjele = 1.0 / Nodes()[i]->NumElement();
      for (unsigned int j = 0; j < num_cols; ++j)
        (*(expolData(j)))[lid] += nodalData(i, j) * invmyadjele;
    }
  }
  return;
}

template void Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex8>::soh8_expol(
    Core::LinAlg::Matrix<numgpt_post, 1>&, Epetra_MultiVector&);
template void Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex8>::soh8_expol(
    Core::LinAlg::Matrix<numgpt_post, numstr_>&, Epetra_MultiVector&);
template void Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex8>::soh8_expol(
    Core::LinAlg::Matrix<numgpt_post, 9>&, Epetra_MultiVector&);

/*----------------------------------------------------------------------*
 | Have plastic spin                                        seitz 05/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
bool Discret::ELEMENTS::So3Plast<distype>::have_plastic_spin()
{
  // get plastic hyperelastic material
  Mat::PlasticElastHyper* plmat = nullptr;
  if (Material()->MaterialType() == Core::Materials::m_plelasthyper)
    plmat = static_cast<Mat::PlasticElastHyper*>(Material().get());

  if (plmat != nullptr) return plmat->have_plastic_spin();

  return false;
}

int Discret::ELEMENTS::PlastEasTypeToNumEasV(Discret::ELEMENTS::So3PlastEasType et)
{
  switch (et)
  {
    case soh8p_easnone:
      return PlastEasTypeToNumEas<soh8p_easnone>::neas;
      break;
    case soh8p_easmild:
      return PlastEasTypeToNumEas<soh8p_easmild>::neas;
      break;
    case soh8p_easfull:
      return PlastEasTypeToNumEas<soh8p_easfull>::neas;
      break;
    case soh8p_eassosh8:
      return PlastEasTypeToNumEas<soh8p_eassosh8>::neas;
      break;
    case soh18p_eassosh18:
      return PlastEasTypeToNumEas<soh18p_eassosh18>::neas;
      break;
    default:
      FOUR_C_THROW("EAS type not implemented");
  }
  return -1;
}

FOUR_C_NAMESPACE_CLOSE

#include "4C_so3_ssn_plast_fwd.hpp"
