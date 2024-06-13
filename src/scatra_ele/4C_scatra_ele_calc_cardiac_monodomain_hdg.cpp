/*----------------------------------------------------------------------*/
/*! \file

 \brief main file containing routines for calculation of HDG cardiac monodomain element

\level 3

 *----------------------------------------------------------------------*/

#include "4C_scatra_ele_calc_cardiac_monodomain_hdg.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_fiber_node.hpp"
#include "4C_fem_general_fiber_node_holder.hpp"
#include "4C_fem_general_fiber_node_utils.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_fem_general_utils_polynomial.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_densematrix_multiply.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_myocard.hpp"
#include "4C_scatra_ele_calc_hdg.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype,
    probdim>::ScaTraEleCalcHDGCardiacMonodomain(const int numdofpernode, const int numscal,
    const std::string& disname)
    : Discret::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::ScaTraEleCalcHDG(
          numdofpernode, numscal, disname),
      values_mat_gp_all_(0),
      gp_mat_alpha_(0)
{
}

/*----------------------------------------------------------------------*
 | singleton access method                               hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>*
Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname, bool create)
{
  static std::map<std::string, ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>*> instances;

  if (create)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] =
          new ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>(numdofpernode, numscal, disname);
  }

  else if (instances.find(disname) != instances.end())
  {
    for (typename std::map<std::string,
             ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>*>::iterator i = instances.begin();
         i != instances.end(); ++i)
    {
      delete i->second;
      i->second = nullptr;
    }

    instances.clear();
    return nullptr;
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 |  prepare material parameter                           hoermann 11/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::prepare_materials_all(
    Core::Elements::Element* ele,                            //!< the element we are dealing with
    const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    Teuchos::RCP<std::vector<Core::LinAlg::SerialDenseMatrix>> difftensor)
{
  const Teuchos::RCP<Mat::Myocard>& actmat =
      Teuchos::rcp_dynamic_cast<Mat::Myocard>(ele->Material());
  Discret::ELEMENTS::ScaTraHDG* hdgele =
      dynamic_cast<Discret::ELEMENTS::ScaTraHDG*>(const_cast<Core::Elements::Element*>(ele));

  if (actmat->diffusion_at_ele_center())
  {
    // get diffusivity at ele center
    Core::LinAlg::Matrix<probdim, probdim> diff(true);
    actmat->Diffusivity(diff, 0);
    Core::LinAlg::SerialDenseMatrix difftensortmp(this->nsd_, this->nsd_);
    for (unsigned int i = 0; i < this->nsd_; ++i)
      for (unsigned int j = 0; j < this->nsd_; ++j) difftensortmp(i, j) = diff(i, j);
    (*difftensor).push_back(difftensortmp);
    Core::Nodes::FiberNode* fnode = dynamic_cast<Core::Nodes::FiberNode*>(ele->Nodes()[0]);
    if (fnode) FOUR_C_THROW("Fiber direction defined twice (nodes and elements)");
  }
  else
  {
    actmat->reset_diffusion_tensor();

    Teuchos::RCP<Core::FE::ShapeValues<distype>> shapes =
        Teuchos::rcp(new Core::FE::ShapeValues<distype>(1, false, 2 * hdgele->Degree()));

    shapes->Evaluate(*ele);

    std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, 1>> shapefcns(shapes->nqpoints_);

    for (std::size_t q = 0; q < shapes->nqpoints_; ++q)
    {
      for (std::size_t i = 0; i < shapes->ndofs_; ++i)
      {
        shapefcns[q](i) = shapes->funct(i, q);
      }
    }

    Core::Nodes::NodalFiberHolder gpFiberHolder;
    Core::Nodes::ProjectFibersToGaussPoints<distype>(ele->Nodes(), shapefcns, gpFiberHolder);

    std::vector<Core::LinAlg::Matrix<probdim, 1>> fibergp(shapes->nqpoints_);
    setup_cardiac_fibers<probdim>(gpFiberHolder, fibergp);

    for (unsigned int q = 0; q < shapes->nqpoints_; ++q) actmat->setup_diffusion_tensor(fibergp[q]);

    for (unsigned int q = 0; q < shapes->nqpoints_; ++q)
    {
      Core::LinAlg::SerialDenseMatrix difftensortmp(this->nsd_, this->nsd_);
      Core::LinAlg::Matrix<probdim, probdim> diff(true);
      actmat->Diffusivity(diff, q);
      for (unsigned int i = 0; i < this->nsd_; ++i)
        for (unsigned int j = 0; j < this->nsd_; ++j) difftensortmp(i, j) = diff(i, j);
      (*difftensor).push_back(difftensortmp);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  prepare material parameter                           hoermann 01/11 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::prepare_materials(
    Core::Elements::Element* ele,                            //!< the element we are dealing with
    const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    Teuchos::RCP<std::vector<Core::LinAlg::SerialDenseMatrix>> difftensor)
{
  if (distype == Core::FE::CellType::tet4 or distype == Core::FE::CellType::tet10)
    prepare_materials_tet(ele, material, k, difftensor);
  else
    prepare_materials_all(ele, material, k, difftensor);

  return;
}

/*----------------------------------------------------------------------*
 |  prepare material parameter                           hoermann 01/11 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::prepare_materials_tet(
    Core::Elements::Element* ele,                            //!< the element we are dealing with
    const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    Teuchos::RCP<std::vector<Core::LinAlg::SerialDenseMatrix>> difftensor)
{
  const Teuchos::RCP<Mat::Myocard>& actmat =
      Teuchos::rcp_dynamic_cast<Mat::Myocard>(ele->Material());
  Discret::ELEMENTS::ScaTraHDG* hdgele =
      dynamic_cast<Discret::ELEMENTS::ScaTraHDG*>(const_cast<Core::Elements::Element*>(ele));

  if (actmat->diffusion_at_ele_center())
  {
    // get diffusivity at ele center
    Core::LinAlg::Matrix<probdim, probdim> diff(true);
    actmat->Diffusivity(diff, 0);
    Core::LinAlg::SerialDenseMatrix difftensortmp(this->nsd_, this->nsd_);
    for (unsigned int i = 0; i < this->nsd_; ++i)
      for (unsigned int j = 0; j < this->nsd_; ++j) difftensortmp(i, j) = diff(i, j);
    (*difftensor).push_back(difftensortmp);
    Core::Nodes::FiberNode* fnode = dynamic_cast<Core::Nodes::FiberNode*>(ele->Nodes()[0]);
    if (fnode) FOUR_C_THROW("Fiber direction defined twice (nodes and elements)");
  }
  else
  {
    actmat->reset_diffusion_tensor();

    const Core::FE::IntPointsAndWeights<Core::FE::dim<distype>> intpoints(
        ScaTra::DisTypeToMatGaussRule<distype>::get_gauss_rule(2 * hdgele->Degree()));
    const std::size_t numgp = intpoints.IP().nquad;

    std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, 1>> shapefcns(numgp);

    Core::LinAlg::Matrix<probdim, 1> gp_coord(true);
    for (std::size_t q = 0; q < numgp; ++q)
    {
      // gaussian points coordinates
      for (int idim = 0; idim < Core::FE::dim<distype>; ++idim)
        gp_coord(idim) = intpoints.IP().qxg[q][idim];

      Core::FE::shape_function<distype>(gp_coord, shapefcns[q]);
    }

    Core::Nodes::NodalFiberHolder gpFiberHolder;
    Core::Nodes::ProjectFibersToGaussPoints<distype>(ele->Nodes(), shapefcns, gpFiberHolder);

    std::vector<Core::LinAlg::Matrix<probdim, 1>> fibergp(numgp);
    setup_cardiac_fibers<probdim>(gpFiberHolder, fibergp);

    for (unsigned int q = 0; q < numgp; ++q) actmat->setup_diffusion_tensor(fibergp[q]);

    for (unsigned int q = 0; q < numgp; ++q)
    {
      Core::LinAlg::SerialDenseMatrix difftensortmp(this->nsd_, this->nsd_);
      Core::LinAlg::Matrix<probdim, probdim> diff(true);
      actmat->Diffusivity(diff, q);
      for (unsigned int i = 0; i < this->nsd_; ++i)
        for (unsigned int j = 0; j < this->nsd_; ++j) difftensortmp(i, j) = diff(i, j);
      (*difftensor).push_back(difftensortmp);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::materials(
    const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    Core::LinAlg::SerialDenseMatrix& difftensor, Core::LinAlg::SerialDenseVector& ivecn,
    Core::LinAlg::SerialDenseVector& ivecnp, Core::LinAlg::SerialDenseMatrix& ivecnpderiv)
{
  if (material->MaterialType() == Core::Materials::m_myocard)
    mat_myocard(material, k, difftensor, ivecn, ivecnp, ivecnpderiv);
  else
    FOUR_C_THROW("Material type is not supported");

  return;
}


/*----------------------------------------------------------------------*
 |  Material ScaTra                                      hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::mat_myocard(
    const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    Core::LinAlg::SerialDenseMatrix& difftensor, Core::LinAlg::SerialDenseVector& ivecn,
    Core::LinAlg::SerialDenseVector& ivecnp, Core::LinAlg::SerialDenseMatrix& ivecnpderiv)
{
  const Teuchos::RCP<const Mat::Myocard>& actmat =
      Teuchos::rcp_dynamic_cast<const Mat::Myocard>(material);

  // coordinate of material gauss points
  Core::LinAlg::Matrix<probdim, 1> mat_gp_coord(true);
  // values of shape function at material gauss points
  Core::LinAlg::SerialDenseVector values_mat_gp(this->shapes_->ndofs_);

  double imatgpnpderiv(0.);
  double imatgpnp(0.);
  double imatgpn(0.);

  ivecn.putScalar(0.0);
  ivecnp.putScalar(0.0);
  ivecnpderiv.putScalar(0.0);

  // polynomial space to get the value of the shape function at the material gauss points
  Core::FE::PolynomialSpaceParams params(distype, this->shapes_->degree_, this->usescompletepoly_);
  polySpace_ = Core::FE::PolynomialSpaceCache<probdim>::Instance().Create(params);

  int nqpoints;

  if (distype == Core::FE::CellType::tet4 or distype == Core::FE::CellType::tet10)
  {
    int deg = 0;
    if (this->shapes_->degree_ == 1)
      deg = 4 * this->shapes_->degree_;
    else
      deg = 3 * this->shapes_->degree_;
    const Core::FE::IntPointsAndWeights<Core::FE::dim<distype>> intpoints(
        ScaTra::DisTypeToMatGaussRule<distype>::get_gauss_rule(deg));
    nqpoints = intpoints.IP().nquad;

    if (nqpoints != actmat->GetNumberOfGP())
      FOUR_C_THROW(
          "Number of quadrature points (%d) does not match number of points in material (%d)!",
          nqpoints, actmat->GetNumberOfGP());

    if (values_mat_gp_all_.empty() or
        values_mat_gp_all_.size() != (unsigned)actmat->GetNumberOfGP())
    {
      values_mat_gp_all_.resize(actmat->GetNumberOfGP());
      gp_mat_alpha_.resize(actmat->GetNumberOfGP());
    }

    if (unsigned(values_mat_gp_all_[0].numRows()) != this->shapes_->ndofs_)
    {
      for (int q = 0; q < nqpoints; ++q)
      {
        values_mat_gp_all_[q].size(this->shapes_->ndofs_);

        gp_mat_alpha_[q] = intpoints.IP().qwgt[q];
        // gaussian points coordinates
        for (int idim = 0; idim < Core::FE::dim<distype>; ++idim)
          mat_gp_coord(idim) = intpoints.IP().qxg[q][idim];

        polySpace_->Evaluate(mat_gp_coord, values_mat_gp_all_[q]);
      }
    }
  }
  else
  {
    int deg = 0;
    if (this->shapes_->degree_ == 1)
      deg = 4 * this->shapes_->degree_;
    else
      deg = 3 * this->shapes_->degree_;

    Teuchos::RCP<Core::FE::GaussPoints> quadrature_(
        Core::FE::GaussPointCache::Instance().Create(distype, deg));
    nqpoints = quadrature_->NumPoints();

    if (nqpoints != actmat->GetNumberOfGP())
      FOUR_C_THROW(
          "Number of quadrature points (%d) does not match number of points in material (%d)!",
          nqpoints, actmat->GetNumberOfGP());

    if (values_mat_gp_all_.empty() or
        values_mat_gp_all_.size() != (unsigned)actmat->GetNumberOfGP())
    {
      values_mat_gp_all_.resize(actmat->GetNumberOfGP());
      gp_mat_alpha_.resize(actmat->GetNumberOfGP());
    }

    if (unsigned(values_mat_gp_all_[0].numRows()) != this->shapes_->ndofs_)
    {
      for (int q = 0; q < nqpoints; ++q)
      {
        values_mat_gp_all_[q].size(this->shapes_->ndofs_);

        gp_mat_alpha_[q] = quadrature_->Weight(q);
        // gaussian points coordinates
        for (int idim = 0; idim < Core::FE::dim<distype>; ++idim)
          mat_gp_coord(idim) = quadrature_->Point(q)[idim];

        polySpace_->Evaluate(mat_gp_coord, values_mat_gp_all_[q]);
      }
    }
  }

  // Jacobian determinant
  double jacdet = this->shapes_->xjm.Determinant();

  for (int q = 0; q < nqpoints; ++q)
  {
    double phinpgp = 0.0;
    double phingp = 0.0;

    // loop over shape functions
    for (unsigned int i = 0; i < this->shapes_->ndofs_; ++i)
    {
      phingp += values_mat_gp_all_[q](i) * this->interiorPhin_(i);
      phinpgp += values_mat_gp_all_[q](i) * this->interiorPhinp_(i);
    }

    // Reaction term at material gauss points

    if (!this->scatrapara_->SemiImplicit())
    {
      imatgpnpderiv = actmat->ReaCoeffDeriv(phinpgp, this->dt(), q);
    }

    imatgpn = actmat->ReaCoeffN(phingp, this->dt(), q);
    imatgpnp = actmat->ReaCoeff(phinpgp, this->dt(), q);

    // loop over shape functions
    for (unsigned int i = 0; i < this->shapes_->ndofs_; ++i)
    {
      ivecn(i) += imatgpn * values_mat_gp_all_[q](i) * jacdet * gp_mat_alpha_[q];
    }

    if (!this->scatrapara_->SemiImplicit())
      for (unsigned int i = 0; i < this->shapes_->ndofs_; ++i)
      {
        for (unsigned int j = 0; j < this->shapes_->ndofs_; ++j)
          ivecnpderiv(i, j) += imatgpnpderiv * values_mat_gp_all_[q](i) * values_mat_gp_all_[q](j) *
                               jacdet * gp_mat_alpha_[q];
        ivecnp(i) += imatgpnp * values_mat_gp_all_[q](i) * jacdet * gp_mat_alpha_[q];
      }
  }

  return;
}  // ScaTraEleCalcHDGCardiacMonodomain<distype>::MatMyocard


/*----------------------------------------------------------------------*
 |  Material Time Update                                 hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::time_update_material(
    const Core::Elements::Element* ele  //!< the element we are dealing with
)
{
  std::vector<Teuchos::RCP<Mat::Myocard>> updatemat;

  // access the general material
  Teuchos::RCP<Core::Mat::Material> material = ele->Material();

  // first, determine the materials which need a time update, i.e. myocard materials
  if (material->MaterialType() == Core::Materials::m_matlist)
  {
    const Teuchos::RCP<Mat::MatList> actmat = Teuchos::rcp_dynamic_cast<Mat::MatList>(material);
    if (actmat->NumMat() < this->numscal_) FOUR_C_THROW("Not enough materials in MatList.");

    for (int k = 0; k < this->numscal_; ++k)
    {
      const int matid = actmat->MatID(k);
      Teuchos::RCP<Core::Mat::Material> singlemat = actmat->MaterialById(matid);

      if (singlemat->MaterialType() == Core::Materials::m_myocard)
      {
        // reference to Teuchos::rcp not possible here, since the material
        // is required to be not const for this application
        updatemat.push_back(Teuchos::rcp_dynamic_cast<Mat::Myocard>(singlemat));
      }
    }
  }

  if (material->MaterialType() == Core::Materials::m_myocard)
  {
    // reference to Teuchos::rcp not possible here, since the material is required to be
    // not const for this application
    updatemat.push_back(Teuchos::rcp_dynamic_cast<Mat::Myocard>(material));
  }

  if (updatemat.size() > 0)  // found at least one material to be updated
  {
    for (unsigned i = 0; i < updatemat.size(); i++) updatemat[i]->Update(Teuchos::null, 0.0);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Get Material Internal State for Restart              hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype,
    probdim>::get_material_internal_state(const Core::Elements::Element*
                                              ele,  //!< the element we are dealing with
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization)
{
  // NOTE: add integral values only for elements which are NOT ghosted!
  if (ele->Owner() == discretization.Comm().MyPID())
  {
    // access the general material
    Teuchos::RCP<Core::Mat::Material> material = ele->Material();
    Teuchos::RCP<Epetra_MultiVector> material_internal_state =
        params.get<Teuchos::RCP<Epetra_MultiVector>>("material_internal_state");

    if (material->MaterialType() == Core::Materials::m_myocard)
    {
      Teuchos::RCP<Mat::Myocard> material =
          Teuchos::rcp_dynamic_cast<Mat::Myocard>(ele->Material());
      for (int k = 0; k < material->get_number_of_internal_state_variables(); ++k)
      {
        double material_state = 0;
        unsigned int nqpoints = material->GetNumberOfGP();
        for (unsigned int q = 0; q < nqpoints; ++q)
        {
          material_state += material->GetInternalState(k, q);
        }
        int err =
            material_internal_state->ReplaceGlobalValue(ele->Id(), k, material_state / nqpoints);
        if (err != 0) FOUR_C_THROW("%i", err);
      }
    }

    params.set<Teuchos::RCP<Epetra_MultiVector>>(
        "material_internal_state", material_internal_state);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Set Material Internal State after Restart            hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype,
    probdim>::set_material_internal_state(const Core::Elements::Element*
                                              ele,  //!< the element we are dealing with
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization)
{
  // NOTE: add integral values only for elements which are NOT ghosted!
  if (ele->Owner() == discretization.Comm().MyPID())
  {
    // access the general material
    Teuchos::RCP<Core::Mat::Material> material = ele->Material();
    Teuchos::RCP<Epetra_MultiVector> material_internal_state =
        params.get<Teuchos::RCP<Epetra_MultiVector>>("material_internal_state");

    if (material->MaterialType() == Core::Materials::m_myocard)
    {
      Teuchos::RCP<Mat::Myocard> material =
          Teuchos::rcp_dynamic_cast<Mat::Myocard>(ele->Material());
      for (int k = 0; k < material->get_number_of_internal_state_variables(); ++k)
      {
        int nqpoints = material->GetNumberOfGP();
        for (int q = 0; q < nqpoints; ++q)
        {
          Teuchos::RCP<Epetra_Vector> material_internal_state_component =
              Teuchos::rcp((*material_internal_state)(k * nqpoints + q), false);
          material->SetInternalState(k, (*material_internal_state_component)[ele->Id()], q);
        }
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------*
 |  Project Material Field                               hoermann 01/17 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
int Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::project_material_field(
    const Core::Elements::Element* ele  //!< the element we are dealing with
)
{
  if (distype == Core::FE::CellType::tet4 or distype == Core::FE::CellType::tet10)
    return project_material_field_tet(ele);
  else
    return project_material_field_all(ele);
}


/*----------------------------------------------------------------------*
 |  Project Material Field                               hoermann 12/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
int Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype,
    probdim>::project_material_field_all(const Core::Elements::Element*
        ele  //!< the element we are dealing with
)
{
  const Teuchos::RCP<Mat::Myocard>& actmat =
      Teuchos::rcp_dynamic_cast<Mat::Myocard>(ele->Material());

  Discret::ELEMENTS::ScaTraHDG* hdgele =
      dynamic_cast<Discret::ELEMENTS::ScaTraHDG*>(const_cast<Core::Elements::Element*>(ele));

  int deg = 0;
  if (hdgele->Degree() == 1)
    deg = 4 * hdgele->Degree();
  else
    deg = 3 * hdgele->Degree();
  int degold = 0;
  if (hdgele->DegreeOld() == 1)
    degold = 4 * hdgele->DegreeOld();
  else
    degold = 3 * hdgele->DegreeOld();

  Teuchos::RCP<Core::FE::ShapeValues<distype>> shapes = Teuchos::rcp(
      new Core::FE::ShapeValues<distype>(hdgele->DegreeOld(), this->usescompletepoly_, deg));

  Teuchos::RCP<Core::FE::ShapeValues<distype>> shapes_old = Teuchos::rcp(
      new Core::FE::ShapeValues<distype>(hdgele->DegreeOld(), this->usescompletepoly_, degold));

  shapes->Evaluate(*ele);
  shapes_old->Evaluate(*ele);

  Core::LinAlg::SerialDenseMatrix massPartOld(shapes->ndofs_, shapes_old->nqpoints_);
  Core::LinAlg::SerialDenseMatrix massPartOldW(shapes->ndofs_, shapes_old->nqpoints_);
  Core::LinAlg::SerialDenseMatrix massPart(shapes->ndofs_, shapes->nqpoints_);
  Core::LinAlg::SerialDenseMatrix massPartW(shapes->ndofs_, shapes->nqpoints_);
  Core::LinAlg::SerialDenseMatrix Mmat(shapes->ndofs_, shapes->ndofs_);

  Core::LinAlg::SerialDenseMatrix state_variables(
      shapes_old->nqpoints_, actmat->get_number_of_internal_state_variables());

  if (shapes->ndofs_ != shapes_old->ndofs_)
    FOUR_C_THROW("Number of shape functions not identical!");

  for (unsigned int i = 0; i < shapes->ndofs_; ++i)
  {
    for (unsigned int q = 0; q < shapes->nqpoints_; ++q)
    {
      massPart(i, q) = shapes->shfunct(i, q);
      massPartW(i, q) = shapes->shfunct(i, q) * shapes->jfac(q);
    }
    for (unsigned int q = 0; q < shapes_old->nqpoints_; ++q)
    {
      massPartOld(i, q) = shapes_old->shfunct(i, q);
      massPartOldW(i, q) = shapes_old->shfunct(i, q) * shapes_old->jfac(q);
    }
  }

  Core::LinAlg::multiplyNT(Mmat, massPartOld, massPartOldW);

  for (unsigned int q = 0; q < shapes_old->nqpoints_; ++q)
    for (int k = 0; k < actmat->get_number_of_internal_state_variables(); ++k)
      state_variables(q, k) = actmat->GetInternalState(k, q);

  Core::LinAlg::SerialDenseMatrix tempMat1(
      shapes->ndofs_, actmat->get_number_of_internal_state_variables());
  Core::LinAlg::multiply(tempMat1, massPartOldW, state_variables);

  using ordinalType = Core::LinAlg::SerialDenseMatrix::ordinalType;
  using scalarType = Core::LinAlg::SerialDenseMatrix::scalarType;
  Teuchos::SerialDenseSolver<ordinalType, scalarType> inverseMat;
  inverseMat.setMatrix(Teuchos::rcpFromRef(Mmat));
  inverseMat.setVectors(Teuchos::rcpFromRef(tempMat1), Teuchos::rcpFromRef(tempMat1));
  inverseMat.factorWithEquilibration(true);
  int err2 = inverseMat.factor();
  int err = inverseMat.solve();
  if (err != 0 || err2 != 0) FOUR_C_THROW("Inversion of matrix failed with errorcode %d", err);

  Core::LinAlg::SerialDenseMatrix tempMat2(
      shapes->nqpoints_, actmat->get_number_of_internal_state_variables());
  Core::LinAlg::multiplyTN(tempMat2, massPart, tempMat1);

  actmat->SetGP(shapes->nqpoints_);
  actmat->resize_internal_state_variables();


  for (unsigned int q = 0; q < shapes->nqpoints_; ++q)
    for (int k = 0; k < actmat->get_number_of_internal_state_variables(); ++k)
      actmat->SetInternalState(k, tempMat2(q, k), q);

  return 0;
}


/*----------------------------------------------------------------------*
 |  Project Material Field for Tet                       hoermann 01/17 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
int Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype,
    probdim>::project_material_field_tet(const Core::Elements::Element*
        ele  //!< the element we are dealing with
)
{
  const Teuchos::RCP<Mat::Myocard>& actmat =
      Teuchos::rcp_dynamic_cast<Mat::Myocard>(ele->Material());

  Discret::ELEMENTS::ScaTraHDG* hdgele =
      dynamic_cast<Discret::ELEMENTS::ScaTraHDG*>(const_cast<Core::Elements::Element*>(ele));

  // polynomial space to get the value of the shape function at the material gauss points
  Core::FE::PolynomialSpaceParams params(distype, hdgele->DegreeOld(), this->usescompletepoly_);
  Teuchos::RCP<Core::FE::PolynomialSpace<probdim>> polySpace =
      Core::FE::PolynomialSpaceCache<probdim>::Instance().Create(params);

  int deg = 0;
  int degold = 0;

  if (hdgele->Degree() == 1)
    deg = 4 * hdgele->Degree();
  else
    deg = 3 * hdgele->Degree();

  if (hdgele->DegreeOld() == 1)
    degold = 4 * hdgele->DegreeOld();
  else
    degold = 3 * hdgele->DegreeOld();

  const Core::FE::IntPointsAndWeights<Core::FE::dim<distype>> intpoints_old(
      ScaTra::DisTypeToMatGaussRule<distype>::get_gauss_rule(degold));
  const Core::FE::IntPointsAndWeights<Core::FE::dim<distype>> intpoints(
      ScaTra::DisTypeToMatGaussRule<distype>::get_gauss_rule(deg));


  std::vector<Core::LinAlg::SerialDenseVector> shape_gp_old(intpoints_old.IP().nquad);
  std::vector<Core::LinAlg::SerialDenseVector> shape_gp(intpoints.IP().nquad);

  // coordinate of material gauss points
  Core::LinAlg::Matrix<probdim, 1> mat_gp_coord(true);

  for (int q = 0; q < intpoints_old.IP().nquad; ++q)
  {
    shape_gp_old[q].size(polySpace->Size());

    // gaussian points coordinates
    for (int idim = 0; idim < Core::FE::dim<distype>; ++idim)
      mat_gp_coord(idim) = intpoints_old.IP().qxg[q][idim];
    polySpace->Evaluate(mat_gp_coord, shape_gp_old[q]);
  }

  for (int q = 0; q < intpoints.IP().nquad; ++q)
  {
    shape_gp[q].size(polySpace->Size());

    // gaussian points coordinates
    for (int idim = 0; idim < Core::FE::dim<distype>; ++idim)
      mat_gp_coord(idim) = intpoints.IP().qxg[q][idim];
    polySpace->Evaluate(mat_gp_coord, shape_gp[q]);
  }

  this->shapes_->Evaluate(*ele);
  // Jacobian determinant
  double jacdet = this->shapes_->xjm.Determinant();



  Core::LinAlg::SerialDenseMatrix massPartOld(polySpace->Size(), shape_gp_old.size());
  Core::LinAlg::SerialDenseMatrix massPartOldW(polySpace->Size(), shape_gp_old.size());
  Core::LinAlg::SerialDenseMatrix massPart(polySpace->Size(), shape_gp.size());
  Core::LinAlg::SerialDenseMatrix massPartW(polySpace->Size(), shape_gp.size());
  Core::LinAlg::SerialDenseMatrix Mmat(polySpace->Size(), polySpace->Size());

  Core::LinAlg::SerialDenseMatrix state_variables(
      shape_gp_old.size(), actmat->get_number_of_internal_state_variables());

  for (unsigned int i = 0; i < polySpace->Size(); ++i)
  {
    for (unsigned int q = 0; q < shape_gp.size(); ++q)
    {
      massPart(i, q) = shape_gp[q](i);
      massPartW(i, q) = shape_gp[q](i) * jacdet * intpoints.IP().qwgt[q];
    }
    for (unsigned int q = 0; q < shape_gp_old.size(); ++q)
    {
      massPartOld(i, q) = shape_gp_old[q](i);
      massPartOldW(i, q) = shape_gp_old[q](i) * jacdet * intpoints_old.IP().qwgt[q];
    }
  }

  Core::LinAlg::multiplyNT(Mmat, massPartOld, massPartOldW);

  for (unsigned int q = 0; q < shape_gp_old.size(); ++q)
    for (int k = 0; k < actmat->get_number_of_internal_state_variables(); ++k)
      state_variables(q, k) = actmat->GetInternalState(k, q);

  Core::LinAlg::SerialDenseMatrix tempMat1(
      polySpace->Size(), actmat->get_number_of_internal_state_variables());
  Core::LinAlg::multiply(tempMat1, massPartOldW, state_variables);

  using ordinalType = Core::LinAlg::SerialDenseMatrix::ordinalType;
  using scalarType = Core::LinAlg::SerialDenseMatrix::scalarType;
  Teuchos::SerialDenseSolver<ordinalType, scalarType> inverseMat;
  inverseMat.setMatrix(Teuchos::rcpFromRef(Mmat));
  inverseMat.setVectors(Teuchos::rcpFromRef(tempMat1), Teuchos::rcpFromRef(tempMat1));
  inverseMat.factorWithEquilibration(true);
  int err2 = inverseMat.factor();
  int err = inverseMat.solve();
  if (err != 0 || err2 != 0) FOUR_C_THROW("Inversion of matrix failed with errorcode %d", err);

  Core::LinAlg::SerialDenseMatrix tempMat2(
      shape_gp.size(), actmat->get_number_of_internal_state_variables());
  Core::LinAlg::multiplyTN(tempMat2, massPart, tempMat1);

  actmat->SetGP(shape_gp.size());
  actmat->resize_internal_state_variables();


  for (unsigned int q = 0; q < shape_gp.size(); ++q)
    for (int k = 0; k < actmat->get_number_of_internal_state_variables(); ++k)
      actmat->SetInternalState(k, tempMat2(q, k), q);

  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/

template <Core::FE::CellType distype, int probdim>
template <std::size_t dim>
void Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::setup_cardiac_fibers(
    const Core::Nodes::NodalFiberHolder& fibers, std::vector<Core::LinAlg::Matrix<dim, 1>>& f)
{
  if (fibers.fibers_size() > 0)
  {
    const std::vector<Core::LinAlg::Matrix<3, 1>>& fib = fibers.get_fiber(0);
    f.resize(fib.size());
    for (std::size_t gp = 0; gp < fib.size(); ++gp)
    {
      for (std::size_t i = 0; i < dim; ++i)
      {
        f[gp](i) = fib[gp](i);
      }
    }
  }
  else if (fibers.contains_coordinate_system_direction(
               Core::Nodes::CoordinateSystemDirection::Circular) &&
           fibers.contains_coordinate_system_direction(
               Core::Nodes::CoordinateSystemDirection::Tangential))
  {
    const std::vector<Core::LinAlg::Matrix<3, 1>>& cir =
        fibers.get_coordinate_system_direction(Core::Nodes::CoordinateSystemDirection::Circular);
    const std::vector<Core::LinAlg::Matrix<3, 1>>& tan =
        fibers.get_coordinate_system_direction(Core::Nodes::CoordinateSystemDirection::Tangential);
    const std::vector<double>& helix = fibers.get_angle(Core::Nodes::AngleType::Helix);
    const std::vector<double>& transverse = fibers.get_angle(Core::Nodes::AngleType::Transverse);
    f.resize(cir.size());

    double deg2rad = M_PI / 180.;
    for (unsigned int gp = 0; gp < cir.size(); ++gp)
    {
      Core::LinAlg::Matrix<3, 1> rad(false);
      rad.CrossProduct(cir[gp], tan[gp]);

      double tmp1 = cos(helix[gp] * deg2rad) * cos(transverse[gp] * deg2rad);
      double tmp2 = sin(helix[gp] * deg2rad) * cos(transverse[gp] * deg2rad);
      double tmp3 = sin(transverse[gp] * deg2rad);

      for (unsigned int i = 0; i < 3; ++i)
      {
        f[gp](i) = tmp1 * cir[gp](i, 0) + tmp2 * tan[gp](i, 0) + tmp3 * rad(i, 0);
      }
      f[gp].Scale(1.0 / f[gp].Norm2());
    }
  }
  else
  {
    FOUR_C_THROW("You have to specify either FIBER1 or CIR, TAN, HELIX and TRANS");
  }
}



// template classes
// 1D elements
// template class
// Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::line2,1>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::line2,2>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::line2,3>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::line3,1>;

// 2D elements
template class Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::tri3>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::quad4, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::quad4, 3>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::quad8>;
template class Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::quad9, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::nurbs9, 2>;

// 3D elements
template class Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::hex8, 3>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::hex27, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::tet4, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::tet10, 3>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::pyramid5,
    3>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
