/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for a truss element used as mechanical link
       between two beam elements

\level 3

*/
/*----------------------------------------------------------------------*/

#include "4C_beaminteraction_link_truss.hpp"

#include "4C_beaminteraction_link.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_truss3.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


BEAMINTERACTION::BeamLinkTrussType BEAMINTERACTION::BeamLinkTrussType::instance_;


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Communication::ParObject* BEAMINTERACTION::BeamLinkTrussType::create(
    const std::vector<char>& data)
{
  BEAMINTERACTION::BeamLinkTruss* my_truss_linker = new BEAMINTERACTION::BeamLinkTruss();
  my_truss_linker->unpack(data);
  return my_truss_linker;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamLinkTruss::BeamLinkTruss()
    : BeamLinkPinJointed(),
      linkele_(Teuchos::null),
      bspotforces_(2, Core::LinAlg::SerialDenseVector(true))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
BEAMINTERACTION::BeamLinkTruss::BeamLinkTruss(const BEAMINTERACTION::BeamLinkTruss& old)
    : BEAMINTERACTION::BeamLinkPinJointed(old),
      bspotforces_(2, Core::LinAlg::SerialDenseVector(true))
{
  if (linkele_ != Teuchos::null)
    linkele_ = Teuchos::rcp_dynamic_cast<Discret::ELEMENTS::Truss3>(
        Teuchos::rcp(old.linkele_->clone(), true));
  else
    linkele_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<BEAMINTERACTION::BeamLink> BEAMINTERACTION::BeamLinkTruss::clone() const
{
  Teuchos::RCP<BEAMINTERACTION::BeamLinkTruss> newlinker =
      Teuchos::rcp(new BEAMINTERACTION::BeamLinkTruss(*this));
  return newlinker;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::init(int id, const std::vector<std::pair<int, int>>& eleids,
    const std::vector<Core::LinAlg::Matrix<3, 1>>& initpos,
    const std::vector<Core::LinAlg::Matrix<3, 3>>& inittriad,
    Inpar::BEAMINTERACTION::CrosslinkerType linkertype, double timelinkwasset)
{
  issetup_ = false;

  BeamLinkPinJointed::init(id, eleids, initpos, inittriad, linkertype, timelinkwasset);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::setup(const int matnum)
{
  check_init();

  // call setup of base class first
  BeamLinkPinJointed::setup(matnum);

  /* the idea is to use a truss element as auxiliary object that provides us with a
   * response force depending on the position of the two material points on the
   *  parent elements (i.e. binding spots) it is connected to;
   *
   * note: the element instance created in this way can only be used in a limited way
   *       because it is not embedded in a discretization. For example,
   *       Nodes() and other methods are not functional because the
   *       pointers to nodes are not set. Same for reference position of nodes via X() ...
   *
   *       We really only use it as a calculation routine for a sophisticated
   *       (displacement-reaction force) relation here! */
  linkele_ = Teuchos::rcp(new Discret::ELEMENTS::Truss3(-1, 0));

  // set material
  linkele_->set_material(0, Mat::Factory(matnum));

  // set cross-section area Fixme hard-coded dummy value for now
  linkele_->set_cross_sec(1.0);

  // set dummy node Ids, in order to make NumNodes() method of element return the correct number of
  // nodes
  int nodeids[] = {-1, -1};
  linkele_->set_node_ids(2, nodeids);

  std::vector<double> refpos(6, 0.0);

  for (unsigned int i = 0; i < 3; ++i)
  {
    refpos[i] = get_bind_spot_pos1()(i);
    refpos[3 + i] = get_bind_spot_pos2()(i);
  }

  linkele_->set_up_reference_geometry(refpos);

  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::pack(Core::Communication::PackBuffer& data) const
{
  check_init_setup();

  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class
  BeamLinkPinJointed::pack(data);

  // pack linker element
  if (linkele_ != Teuchos::null) linkele_->pack(data);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // extract base class
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  BeamLinkPinJointed::unpack(basedata);

  // Unpack data of sub material (these lines are copied from element.cpp)
  std::vector<char> dataele;
  extract_from_pack(position, data, dataele);
  if (dataele.size() > 0)
  {
    Core::Communication::ParObject* object =
        Core::Communication::Factory(dataele);  // Unpack is done here
    Discret::ELEMENTS::Truss3* linkele = dynamic_cast<Discret::ELEMENTS::Truss3*>(object);
    if (linkele == nullptr) FOUR_C_THROW("failed to unpack Truss3 object within BeamLinkTruss");
    linkele_ = Teuchos::rcp(linkele);
  }
  else
    linkele_ = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::BeamLinkTruss::evaluate_force(
    Core::LinAlg::SerialDenseVector& forcevec1, Core::LinAlg::SerialDenseVector& forcevec2)
{
  check_init_setup();

  std::map<std::string, std::vector<double>> ele_state;
  get_disp_for_element_evaluation(ele_state);

  Core::LinAlg::SerialDenseVector force(6, true);
  Core::LinAlg::SerialDenseMatrix stiffmat(6, 6, true);

  linkele_->calc_internal_force_stiff_tot_lag(ele_state, force, stiffmat);

  std::copy(&force(0), &force(0) + 3, &forcevec1(0));
  std::copy(&force(0) + 3, &force(0) + 6, &forcevec2(0));

  bspotforces_[0] = forcevec1;
  bspotforces_[1] = forcevec2;

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::BeamLinkTruss::evaluate_stiff(Core::LinAlg::SerialDenseMatrix& stiffmat11,
    Core::LinAlg::SerialDenseMatrix& stiffmat12, Core::LinAlg::SerialDenseMatrix& stiffmat21,
    Core::LinAlg::SerialDenseMatrix& stiffmat22)
{
  check_init_setup();

  std::map<std::string, std::vector<double>> ele_state;
  get_disp_for_element_evaluation(ele_state);

  Core::LinAlg::SerialDenseVector force(6, true);
  Core::LinAlg::SerialDenseMatrix stiffmat(6, 6, true);

  linkele_->calc_internal_force_stiff_tot_lag(ele_state, force, stiffmat);

  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      stiffmat11(i, j) = stiffmat(i, j);
      stiffmat12(i, j) = stiffmat(i, 3 + j);
      stiffmat21(i, j) = stiffmat(3 + i, j);
      stiffmat22(i, j) = stiffmat(3 + i, 3 + j);
    }
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::BeamLinkTruss::evaluate_force_stiff(
    Core::LinAlg::SerialDenseVector& forcevec1, Core::LinAlg::SerialDenseVector& forcevec2,
    Core::LinAlg::SerialDenseMatrix& stiffmat11, Core::LinAlg::SerialDenseMatrix& stiffmat12,
    Core::LinAlg::SerialDenseMatrix& stiffmat21, Core::LinAlg::SerialDenseMatrix& stiffmat22)
{
  check_init_setup();

  std::map<std::string, std::vector<double>> ele_state;
  get_disp_for_element_evaluation(ele_state);

  Core::LinAlg::SerialDenseVector force(6, true);
  Core::LinAlg::SerialDenseMatrix stiffmat(6, 6, true);

  linkele_->calc_internal_force_stiff_tot_lag(ele_state, force, stiffmat);

  std::copy(&force(0), &force(0) + 3, &forcevec1(0));
  std::copy(&force(0) + 3, &force(0) + 6, &forcevec2(0));

  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      stiffmat11(i, j) = stiffmat(i, j);
      stiffmat12(i, j) = stiffmat(i, 3 + j);
      stiffmat21(i, j) = stiffmat(3 + i, j);
      stiffmat22(i, j) = stiffmat(3 + i, 3 + j);
    }
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::reset_state(std::vector<Core::LinAlg::Matrix<3, 1>>& bspotpos,
    std::vector<Core::LinAlg::Matrix<3, 3>>& bspottriad)
{
  check_init_setup();

  BeamLinkPinJointed::reset_state(bspotpos, bspottriad);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::fill_state_variables_for_element_evaluation(
    Core::LinAlg::Matrix<6, 1, double>& absolute_nodal_positions) const
{
  for (unsigned int i = 0; i < 3; ++i)
  {
    absolute_nodal_positions(i) = get_bind_spot_pos1()(i);
    absolute_nodal_positions(3 + i) = get_bind_spot_pos2()(i);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::get_disp_for_element_evaluation(
    std::map<std::string, std::vector<double>>& ele_state) const
{
  const auto ref_position = linkele_->x();
  std::vector<double> disp(6, 0);
  for (unsigned int i = 0; i < 3; ++i)
  {
    disp[i] = get_bind_spot_pos1()(i) - ref_position(i, 0);
    disp[3 + i] = get_bind_spot_pos2()(i) - ref_position(i + 3, 0);
  }

  ele_state.emplace(std::make_pair("disp", disp));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::scale_linker_reference_length(double scalefac)
{
  linkele_->scale_reference_length(scalefac);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkTruss::get_binding_spot_force(
    int bspotid, Core::LinAlg::SerialDenseVector& bspotforce) const
{
  bspotforce = bspotforces_[bspotid];
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double BEAMINTERACTION::BeamLinkTruss::get_current_linker_length() const
{
  Core::LinAlg::Matrix<6, 1> xcurr;

  fill_state_variables_for_element_evaluation(xcurr);

  Core::LinAlg::Matrix<6, 1> curr_nodal_coords;
  curr_nodal_coords(0) = xcurr(0) - xcurr(3);
  curr_nodal_coords(1) = xcurr(1) - xcurr(4);
  curr_nodal_coords(2) = xcurr(2) - xcurr(5);
  curr_nodal_coords(3) = curr_nodal_coords(0);
  curr_nodal_coords(4) = curr_nodal_coords(1);
  curr_nodal_coords(5) = curr_nodal_coords(2);

  return linkele_->curr_length(curr_nodal_coords);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double BEAMINTERACTION::BeamLinkTruss::get_internal_energy() const
{
  return linkele_->get_internal_energy();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double BEAMINTERACTION::BeamLinkTruss::get_kinetic_energy() const { return 0.0; }

FOUR_C_NAMESPACE_CLOSE
