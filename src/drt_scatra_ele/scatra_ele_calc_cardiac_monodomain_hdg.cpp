/*----------------------------------------------------------------------*/
/*! \file

 \brief main file containing routines for calculation of HDG cardiac monodomain element

\level 3

\maintainer Martin Kronbichler
 *----------------------------------------------------------------------*/

#include "scatra_ele_calc_cardiac_monodomain_hdg.H"

#include "scatra_ele_calc_hdg.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/standardtypes_cpp.H"  // for EPS13 and so on

#include "scatra_ele_parameter_timint.H"
#include "scatra_ele_parameter_std.H"

#include "../drt_fiber/drt_fiber_node.H"

#include "../drt_mat/myocard.H"
#include "../drt_mat/matlist.H"

#include "../drt_fem_general/drt_utils_polynomial.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype,
    probdim>::ScaTraEleCalcHDGCardiacMonodomain(const int numdofpernode, const int numscal,
    const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::ScaTraEleCalcHDG(
          numdofpernode, numscal, disname),
      values_mat_gp_all_(0),
      gp_mat_alpha_(0)
{
}

/*----------------------------------------------------------------------*
 | singleton access method                               hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>*
DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::Instance(
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
      i->second = NULL;
    }

    instances.clear();
    return NULL;
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                 hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(0, 0, "", false);
}

/*----------------------------------------------------------------------*
 |  prepare material parameter                           hoermann 11/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::PrepareMaterialsAll(
    DRT::Element* ele,                                 //!< the element we are dealing with
    const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
    const int k,                                       //!< id of current scalar
    Teuchos::RCP<std::vector<Epetra_SerialDenseMatrix>> difftensor)
{
  const Teuchos::RCP<MAT::Myocard>& actmat =
      Teuchos::rcp_dynamic_cast<MAT::Myocard>(ele->Material());
  DRT::ELEMENTS::ScaTraHDG* hdgele =
      dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<DRT::Element*>(ele));

  if (actmat->DiffusionAtEleCenter())
  {
    // get diffusivity at ele center
    LINALG::Matrix<probdim, probdim> diff(true);
    actmat->Diffusivity(diff, 0);
    Epetra_SerialDenseMatrix difftensortmp(this->nsd_, this->nsd_);
    for (unsigned int i = 0; i < this->nsd_; ++i)
      for (unsigned int j = 0; j < this->nsd_; ++j) difftensortmp(i, j) = diff(i, j);
    (*difftensor).push_back(difftensortmp);
    DRT::FIBER::FiberNode* fnode = dynamic_cast<DRT::FIBER::FiberNode*>(ele->Nodes()[0]);
    if (fnode) dserror("Fiber direction defined twice (nodes and elements)");
  }
  else
  {
    actmat->ResetDiffusionTensor();

    // get fiber information at corners of element
    DRT::Node** nodes = ele->Nodes();
    int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
    std::vector<Epetra_SerialDenseVector> fibernodes(numnodes, Epetra_SerialDenseVector(3));
    std::vector<Epetra_SerialDenseVector> cirnodes(numnodes, Epetra_SerialDenseVector(3));
    std::vector<Epetra_SerialDenseVector> tannodes(numnodes, Epetra_SerialDenseVector(3));
    std::vector<double> helixnodes(numnodes);
    std::vector<double> transversenodes(numnodes);

    for (int inode = 0; inode < numnodes; ++inode)
    {
      DRT::FIBER::FiberNode* fnode = dynamic_cast<DRT::FIBER::FiberNode*>(nodes[inode]);
      if (!fnode) dserror("No fiber direction defined on nodes or elements");
      for (unsigned int j = 0; j < this->nsd_; ++j)
      {
        fibernodes[inode](j) = fnode->Fiber1()[j];
        cirnodes[inode](j) = fnode->Cir()[j];
        tannodes[inode](j) = fnode->Tan()[j];
      }
      helixnodes[inode] = fnode->Helix();
      transversenodes[inode] = fnode->Transverse();
    }

    Teuchos::RCP<DRT::UTILS::ShapeValues<distype>> shapes =
        Teuchos::rcp(new DRT::UTILS::ShapeValues<distype>(1, false, 2 * hdgele->Degree()));

    shapes->Evaluate(*ele);

    std::vector<LINALG::Matrix<probdim, 1>> fibergp(shapes->nqpoints_);

    // case fibers are interpolated
    if ((fibernodes[0]).Norm2() > 1e-13)
    {
      // interpolate fibers to integration points
      for (unsigned int i = 0; i < this->nsd_; ++i)
        for (unsigned int q = 0; q < shapes->nqpoints_; ++q)
          for (unsigned int j = 0; j < shapes->ndofs_; ++j)
            fibergp[q](i, 0) += shapes->funct(j, q) * fibernodes[j](i);
    }
    else  // case coordinate system is interpolated
    {
      std::vector<LINALG::Matrix<probdim, 1>> cirgp(shapes->nqpoints_);
      std::vector<LINALG::Matrix<probdim, 1>> tangp(shapes->nqpoints_);
      std::vector<LINALG::Matrix<probdim, 1>> radgp(shapes->nqpoints_);
      std::vector<double> helixgp(shapes->nqpoints_);
      std::vector<double> transversegp(shapes->nqpoints_);
      Epetra_SerialDenseVector tmpvector(3);


      // interpolate circumferential and tangential directions to integration points
      for (unsigned int q = 0; q < shapes->nqpoints_; ++q)
      {
        LINALG::Matrix<probdim, 1> tmptan;
        for (unsigned int j = 0; j < shapes->ndofs_; ++j)
        {
          for (unsigned int i = 0; i < this->nsd_; ++i)
          {
            cirgp[q](i, 0) += shapes->funct(j, q) * cirnodes[j](i);
            tmptan(i, 0) += shapes->funct(j, q) * tannodes[j](i);
          }
          helixgp[q] += shapes->funct(j, q) * helixnodes[j];
          transversegp[q] += shapes->funct(j, q) * transversenodes[j];
        }
        cirgp[q].Scale(1 / cirgp[q].Norm2());
        tmptan.Scale(1 / tmptan.Norm2());
        // ensure that the circumferential direction is orthogonal to the tangential direction
        for (unsigned int i = 0; i < this->nsd_; ++i)
          tangp[q](i, 0) = tmptan(i, 0) - tmptan.Dot(cirgp[q]) * cirgp[q](i, 0);
        tangp[q].Scale(1 / tangp[q].Norm2());
      }
      for (unsigned int q = 0; q < shapes->nqpoints_; ++q)
      {
        // Cross product and normalize
        tmpvector(0) = cirgp[q](1, 0) * tangp[q](2, 0) - cirgp[q](2, 0) * tangp[q](1, 0);
        tmpvector(1) = cirgp[q](2, 0) * tangp[q](0, 0) - cirgp[q](0, 0) * tangp[q](2, 0);
        tmpvector(2) = cirgp[q](0, 0) * tangp[q](1, 0) - cirgp[q](1, 0) * tangp[q](0, 0);

        for (unsigned int i = 0; i < this->nsd_; ++i) radgp[q] = tmpvector(i) / tmpvector.Norm2();
      }

      double scale = PI / 180.;
      for (unsigned int q = 0; q < shapes->nqpoints_; ++q)
      {
        double tmp1 = cos(helixgp[q] * scale) * cos(transversegp[q] * scale);
        double tmp2 = sin(helixgp[q] * scale) * cos(transversegp[q] * scale);
        double tmp3 = sin(transversegp[q] * scale);

        for (unsigned int i = 0; i < 3; ++i)
        {
          tmpvector(i) = tmp1 * cirgp[q](i, 0) + tmp2 * tangp[q](i, 0) + tmp3 * radgp[q](i, 0);
          fibergp[q](i, 0) = tmpvector(i) / tmpvector.Norm2();
        }
      }
    }

    for (unsigned int q = 0; q < shapes->nqpoints_; ++q) actmat->SetupDiffusionTensor(fibergp[q]);

    for (unsigned int q = 0; q < shapes->nqpoints_; ++q)
    {
      Epetra_SerialDenseMatrix difftensortmp(this->nsd_, this->nsd_);
      LINALG::Matrix<probdim, probdim> diff(true);
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
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::PrepareMaterials(
    DRT::Element* ele,                                 //!< the element we are dealing with
    const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
    const int k,                                       //!< id of current scalar
    Teuchos::RCP<std::vector<Epetra_SerialDenseMatrix>> difftensor)
{
  if (distype == DRT::Element::tet4 or distype == DRT::Element::tet10)
    PrepareMaterialsTet(ele, material, k, difftensor);
  else
    PrepareMaterialsAll(ele, material, k, difftensor);

  return;
}

/*----------------------------------------------------------------------*
 |  prepare material parameter                           hoermann 01/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::PrepareMaterialsTet(
    DRT::Element* ele,                                 //!< the element we are dealing with
    const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
    const int k,                                       //!< id of current scalar
    Teuchos::RCP<std::vector<Epetra_SerialDenseMatrix>> difftensor)
{
  const Teuchos::RCP<MAT::Myocard>& actmat =
      Teuchos::rcp_dynamic_cast<MAT::Myocard>(ele->Material());
  DRT::ELEMENTS::ScaTraHDG* hdgele =
      dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<DRT::Element*>(ele));

  if (actmat->DiffusionAtEleCenter())
  {
    // get diffusivity at ele center
    LINALG::Matrix<probdim, probdim> diff(true);
    actmat->Diffusivity(diff, 0);
    Epetra_SerialDenseMatrix difftensortmp(this->nsd_, this->nsd_);
    for (unsigned int i = 0; i < this->nsd_; ++i)
      for (unsigned int j = 0; j < this->nsd_; ++j) difftensortmp(i, j) = diff(i, j);
    (*difftensor).push_back(difftensortmp);
    DRT::FIBER::FiberNode* fnode = dynamic_cast<DRT::FIBER::FiberNode*>(ele->Nodes()[0]);
    if (fnode) dserror("Fiber direction defined twice (nodes and elements)");
  }
  else
  {
    actmat->ResetDiffusionTensor();

    // get fiber information at corners of element
    DRT::Node** nodes = ele->Nodes();
    int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
    std::vector<Epetra_SerialDenseVector> fibernodes(numnodes, Epetra_SerialDenseVector(3));
    std::vector<Epetra_SerialDenseVector> cirnodes(numnodes, Epetra_SerialDenseVector(3));
    std::vector<Epetra_SerialDenseVector> tannodes(numnodes, Epetra_SerialDenseVector(3));
    std::vector<double> helixnodes(numnodes);
    std::vector<double> transversenodes(numnodes);

    for (int inode = 0; inode < numnodes; ++inode)
    {
      DRT::FIBER::FiberNode* fnode = dynamic_cast<DRT::FIBER::FiberNode*>(nodes[inode]);
      if (!fnode) dserror("No fiber direction defined on nodes or elements");
      for (unsigned int j = 0; j < this->nsd_; ++j)
      {
        fibernodes[inode](j) = fnode->Fiber1()[j];
        cirnodes[inode](j) = fnode->Cir()[j];
        tannodes[inode](j) = fnode->Tan()[j];
      }
      helixnodes[inode] = fnode->Helix();
      transversenodes[inode] = fnode->Transverse();
    }

    const DRT::UTILS::IntPointsAndWeights<DRT::UTILS::DisTypeToDim<distype>::dim> intpoints(
        SCATRA::DisTypeToMatGaussRule<distype>::GetGaussRule(2 * hdgele->Degree()));
    unsigned int nqpoints = intpoints.IP().nquad;

    // coordinate of gauss points
    LINALG::Matrix<probdim, 1> gp_coord(true);
    LINALG::SerialDenseMatrix funct(this->nen_, nqpoints);

    for (int q = 0; q < intpoints.IP().nquad; ++q)
    {
      // gaussian points coordinates
      for (int idim = 0; idim < DRT::UTILS::DisTypeToDim<distype>::dim; ++idim)
        gp_coord(idim) = intpoints.IP().qxg[q][idim];

      LINALG::Matrix<4, 1> myfunct(funct.A() + q * this->nen_, true);
      DRT::UTILS::shape_function<distype>(gp_coord, myfunct);
    }

    std::vector<LINALG::Matrix<probdim, 1>> fibergp(nqpoints);

    // case fibers are interpolated
    if ((fibernodes[0]).Norm2() > 1e-13)
    {
      // interpolate fibers to integration points
      for (unsigned int i = 0; i < this->nsd_; ++i)
        for (unsigned int q = 0; q < nqpoints; ++q)
          for (unsigned int j = 0; j < this->nen_; ++j)
            fibergp[q](i, 0) += funct(j, q) * fibernodes[j](i);
    }
    else  // case coordinate system is interpolated
    {
      std::vector<LINALG::Matrix<probdim, 1>> cirgp(nqpoints);
      std::vector<LINALG::Matrix<probdim, 1>> tangp(nqpoints);
      std::vector<LINALG::Matrix<probdim, 1>> radgp(nqpoints);
      std::vector<double> helixgp(nqpoints);
      std::vector<double> transversegp(nqpoints);
      Epetra_SerialDenseVector tmpvector(3);


      // interpolate circumferential and tangential directions to integration points
      for (unsigned int q = 0; q < nqpoints; ++q)
      {
        LINALG::Matrix<probdim, 1> tmptan;
        for (unsigned int j = 0; j < this->nen_; ++j)
        {
          for (unsigned int i = 0; i < this->nsd_; ++i)
          {
            cirgp[q](i, 0) += funct(j, q) * cirnodes[j](i);
            tmptan(i, 0) += funct(j, q) * tannodes[j](i);
          }
          helixgp[q] += funct(j, q) * helixnodes[j];
          transversegp[q] += funct(j, q) * transversenodes[j];
        }
        cirgp[q].Scale(1 / cirgp[q].Norm2());
        tmptan.Scale(1 / tmptan.Norm2());
        // ensure that the circumferential direction is orthogonal to the tangential direction
        for (unsigned int i = 0; i < this->nsd_; ++i)
          tangp[q](i, 0) = tmptan(i, 0) - tmptan.Dot(cirgp[q]) * cirgp[q](i, 0);
        tangp[q].Scale(1 / tangp[q].Norm2());
      }
      for (unsigned int q = 0; q < nqpoints; ++q)
      {
        // Cross product and normalize
        tmpvector(0) = cirgp[q](1, 0) * tangp[q](2, 0) - cirgp[q](2, 0) * tangp[q](1, 0);
        tmpvector(1) = cirgp[q](2, 0) * tangp[q](0, 0) - cirgp[q](0, 0) * tangp[q](2, 0);
        tmpvector(2) = cirgp[q](0, 0) * tangp[q](1, 0) - cirgp[q](1, 0) * tangp[q](0, 0);

        for (unsigned int i = 0; i < this->nsd_; ++i) radgp[q] = tmpvector(i) / tmpvector.Norm2();
      }

      double scale = PI / 180.;
      for (unsigned int q = 0; q < nqpoints; ++q)
      {
        double tmp1 = cos(helixgp[q] * scale) * cos(transversegp[q] * scale);
        double tmp2 = sin(helixgp[q] * scale) * cos(transversegp[q] * scale);
        double tmp3 = sin(transversegp[q] * scale);

        for (unsigned int i = 0; i < 3; ++i)
        {
          tmpvector(i) = tmp1 * cirgp[q](i, 0) + tmp2 * tangp[q](i, 0) + tmp3 * radgp[q](i, 0);
          fibergp[q](i, 0) = tmpvector(i) / tmpvector.Norm2();
        }
      }
    }


    for (unsigned int q = 0; q < nqpoints; ++q) actmat->SetupDiffusionTensor(fibergp[q]);

    for (unsigned int q = 0; q < nqpoints; ++q)
    {
      Epetra_SerialDenseMatrix difftensortmp(this->nsd_, this->nsd_);
      LINALG::Matrix<probdim, probdim> diff(true);
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
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::Materials(
    const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
    const int k,                                       //!< id of current scalar
    Epetra_SerialDenseMatrix& difftensor, Epetra_SerialDenseVector& ivecn,
    Epetra_SerialDenseVector& ivecnp, Epetra_SerialDenseMatrix& ivecnpderiv)
{
  if (material->MaterialType() == INPAR::MAT::m_myocard)
    MatMyocard(material, k, difftensor, ivecn, ivecnp, ivecnpderiv);
  else
    dserror("Material type is not supported");

  return;
}


/*----------------------------------------------------------------------*
 |  Material ScaTra                                      hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::MatMyocard(
    const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
    const int k,                                       //!< id of current scalar
    Epetra_SerialDenseMatrix& difftensor, Epetra_SerialDenseVector& ivecn,
    Epetra_SerialDenseVector& ivecnp, Epetra_SerialDenseMatrix& ivecnpderiv)
{
  const Teuchos::RCP<const MAT::Myocard>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::Myocard>(material);

  // coordinate of material gauss points
  LINALG::Matrix<probdim, 1> mat_gp_coord(true);
  // values of shape function at material gauss points
  Epetra_SerialDenseVector values_mat_gp(this->shapes_->ndofs_);

  double imatgpnpderiv(0.);
  double imatgpnp(0.);
  double imatgpn(0.);

  ivecn.Scale(0.);
  ivecnp.Scale(0.);
  ivecnpderiv.Scale(0.);

  // polynomial space to get the value of the shape function at the material gauss points
  DRT::UTILS::PolynomialSpaceParams params(
      distype, this->shapes_->degree_, this->usescompletepoly_);
  polySpace_ = DRT::UTILS::PolynomialSpaceCache<probdim>::Instance().Create(params);

  int nqpoints;

  if (distype == DRT::Element::tet4 or distype == DRT::Element::tet10)
  {
    int deg = 0;
    if (this->shapes_->degree_ == 1)
      deg = 4 * this->shapes_->degree_;
    else
      deg = 3 * this->shapes_->degree_;
    const DRT::UTILS::IntPointsAndWeights<DRT::UTILS::DisTypeToDim<distype>::dim> intpoints(
        SCATRA::DisTypeToMatGaussRule<distype>::GetGaussRule(deg));
    nqpoints = intpoints.IP().nquad;

    if (nqpoints != actmat->GetNumberOfGP())
      dserror("Number of quadrature points (%d) does not match number of points in material (%d)!",
          nqpoints, actmat->GetNumberOfGP());

    if (values_mat_gp_all_.empty() or
        values_mat_gp_all_.size() != (unsigned)actmat->GetNumberOfGP())
    {
      values_mat_gp_all_.resize(actmat->GetNumberOfGP());
      gp_mat_alpha_.resize(actmat->GetNumberOfGP());
    }

    if (unsigned(values_mat_gp_all_[0].M()) != this->shapes_->ndofs_)
    {
      for (int q = 0; q < nqpoints; ++q)
      {
        values_mat_gp_all_[q].Size(this->shapes_->ndofs_);

        gp_mat_alpha_[q] = intpoints.IP().qwgt[q];
        // gaussian points coordinates
        for (int idim = 0; idim < DRT::UTILS::DisTypeToDim<distype>::dim; ++idim)
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

    Teuchos::RCP<DRT::UTILS::GaussPoints> quadrature_(
        DRT::UTILS::GaussPointCache::Instance().Create(distype, deg));
    nqpoints = quadrature_->NumPoints();

    if (nqpoints != actmat->GetNumberOfGP())
      dserror("Number of quadrature points (%d) does not match number of points in material (%d)!",
          nqpoints, actmat->GetNumberOfGP());

    if (values_mat_gp_all_.empty() or
        values_mat_gp_all_.size() != (unsigned)actmat->GetNumberOfGP())
    {
      values_mat_gp_all_.resize(actmat->GetNumberOfGP());
      gp_mat_alpha_.resize(actmat->GetNumberOfGP());
    }

    if (unsigned(values_mat_gp_all_[0].M()) != this->shapes_->ndofs_)
    {
      for (int q = 0; q < nqpoints; ++q)
      {
        values_mat_gp_all_[q].Size(this->shapes_->ndofs_);

        gp_mat_alpha_[q] = quadrature_->Weight(q);
        // gaussian points coordinates
        for (int idim = 0; idim < DRT::UTILS::DisTypeToDim<distype>::dim; ++idim)
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
      imatgpnpderiv = actmat->ReaCoeffDeriv(phinpgp, this->Dt(), q);
    }

    imatgpn = actmat->ReaCoeffN(phingp, this->Dt(), q);
    imatgpnp = actmat->ReaCoeff(phinpgp, this->Dt(), q);

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
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::TimeUpdateMaterial(
    const DRT::Element* ele  //!< the element we are dealing with
)
{
  std::vector<Teuchos::RCP<MAT::Myocard>> updatemat;

  // access the general material
  Teuchos::RCP<MAT::Material> material = ele->Material();

  // first, determine the materials which need a time update, i.e. myocard materials
  if (material->MaterialType() == INPAR::MAT::m_matlist)
  {
    const Teuchos::RCP<MAT::MatList> actmat = Teuchos::rcp_dynamic_cast<MAT::MatList>(material);
    if (actmat->NumMat() < this->numscal_) dserror("Not enough materials in MatList.");

    for (int k = 0; k < this->numscal_; ++k)
    {
      const int matid = actmat->MatID(k);
      Teuchos::RCP<MAT::Material> singlemat = actmat->MaterialById(matid);

      if (singlemat->MaterialType() == INPAR::MAT::m_myocard)
      {
        // reference to Teuchos::rcp not possible here, since the material
        // is required to be not const for this application
        updatemat.push_back(Teuchos::rcp_dynamic_cast<MAT::Myocard>(singlemat));
      }
    }
  }

  if (material->MaterialType() == INPAR::MAT::m_myocard)
  {
    // reference to Teuchos::rcp not possible here, since the material is required to be
    // not const for this application
    updatemat.push_back(Teuchos::rcp_dynamic_cast<MAT::Myocard>(material));
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
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::GetMaterialInternalState(
    const DRT::Element* ele,  //!< the element we are dealing with
    Teuchos::ParameterList& params, DRT::Discretization& discretization)
{
  // NOTE: add integral values only for elements which are NOT ghosted!
  if (ele->Owner() == discretization.Comm().MyPID())
  {
    // access the general material
    Teuchos::RCP<MAT::Material> material = ele->Material();
    Teuchos::RCP<Epetra_MultiVector> material_internal_state =
        params.get<Teuchos::RCP<Epetra_MultiVector>>("material_internal_state");

    if (material->MaterialType() == INPAR::MAT::m_myocard)
    {
      Teuchos::RCP<MAT::Myocard> material =
          Teuchos::rcp_dynamic_cast<MAT::Myocard>(ele->Material());
      for (int k = 0; k < material->GetNumberOfInternalStateVariables(); ++k)
      {
        double material_state = 0;
        unsigned int nqpoints = material->GetNumberOfGP();
        for (unsigned int q = 0; q < nqpoints; ++q)
        {
          material_state += material->GetInternalState(k, q);
        }
        int err =
            material_internal_state->ReplaceGlobalValue(ele->Id(), k, material_state / nqpoints);
        if (err != 0) dserror("%i", err);
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
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::SetMaterialInternalState(
    const DRT::Element* ele,  //!< the element we are dealing with
    Teuchos::ParameterList& params, DRT::Discretization& discretization)
{
  // NOTE: add integral values only for elements which are NOT ghosted!
  if (ele->Owner() == discretization.Comm().MyPID())
  {
    // access the general material
    Teuchos::RCP<MAT::Material> material = ele->Material();
    Teuchos::RCP<Epetra_MultiVector> material_internal_state =
        params.get<Teuchos::RCP<Epetra_MultiVector>>("material_internal_state");

    if (material->MaterialType() == INPAR::MAT::m_myocard)
    {
      Teuchos::RCP<MAT::Myocard> material =
          Teuchos::rcp_dynamic_cast<MAT::Myocard>(ele->Material());
      for (int k = 0; k < material->GetNumberOfInternalStateVariables(); ++k)
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
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::ProjectMaterialField(
    const DRT::Element* ele  //!< the element we are dealing with
)
{
  if (distype == DRT::Element::tet4 or distype == DRT::Element::tet10)
    return ProjectMaterialFieldTet(ele);
  else
    return ProjectMaterialFieldAll(ele);
}


/*----------------------------------------------------------------------*
 |  Project Material Field                               hoermann 12/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::ProjectMaterialFieldAll(
    const DRT::Element* ele  //!< the element we are dealing with
)
{
  const Teuchos::RCP<MAT::Myocard>& actmat =
      Teuchos::rcp_dynamic_cast<MAT::Myocard>(ele->Material());

  DRT::ELEMENTS::ScaTraHDG* hdgele =
      dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<DRT::Element*>(ele));

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

  Teuchos::RCP<DRT::UTILS::ShapeValues<distype>> shapes = Teuchos::rcp(
      new DRT::UTILS::ShapeValues<distype>(hdgele->DegreeOld(), this->usescompletepoly_, deg));

  Teuchos::RCP<DRT::UTILS::ShapeValues<distype>> shapes_old = Teuchos::rcp(
      new DRT::UTILS::ShapeValues<distype>(hdgele->DegreeOld(), this->usescompletepoly_, degold));

  shapes->Evaluate(*ele);
  shapes_old->Evaluate(*ele);

  Epetra_SerialDenseMatrix massPartOld(shapes->ndofs_, shapes_old->nqpoints_);
  Epetra_SerialDenseMatrix massPartOldW(shapes->ndofs_, shapes_old->nqpoints_);
  Epetra_SerialDenseMatrix massPart(shapes->ndofs_, shapes->nqpoints_);
  Epetra_SerialDenseMatrix massPartW(shapes->ndofs_, shapes->nqpoints_);
  Epetra_SerialDenseMatrix Mmat(shapes->ndofs_, shapes->ndofs_);

  Epetra_SerialDenseMatrix state_variables(
      shapes_old->nqpoints_, actmat->GetNumberOfInternalStateVariables());

  if (shapes->ndofs_ != shapes_old->ndofs_) dserror("Number of shape functions not identical!");

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

  Mmat.Multiply('N', 'T', 1.0, massPartOld, massPartOldW, 0.0);

  for (unsigned int q = 0; q < shapes_old->nqpoints_; ++q)
    for (int k = 0; k < actmat->GetNumberOfInternalStateVariables(); ++k)
      state_variables(q, k) = actmat->GetInternalState(k, q);

  Epetra_SerialDenseMatrix tempMat1(shapes->ndofs_, actmat->GetNumberOfInternalStateVariables());
  tempMat1.Multiply('N', 'N', 1.0, massPartOldW, state_variables, 0.0);

  Epetra_SerialDenseSolver inverseMat;
  inverseMat.SetMatrix(Mmat);
  inverseMat.SetVectors(tempMat1, tempMat1);
  inverseMat.FactorWithEquilibration(true);
  int err2 = inverseMat.Factor();
  int err = inverseMat.Solve();
  if (err != 0 || err2 != 0) dserror("Inversion of matrix failed with errorcode %d", err);

  Epetra_SerialDenseMatrix tempMat2(shapes->nqpoints_, actmat->GetNumberOfInternalStateVariables());
  tempMat2.Multiply('T', 'N', 1.0, massPart, tempMat1, 0.0);

  actmat->SetGP(shapes->nqpoints_);
  actmat->ResizeInternalStateVariables();


  for (unsigned int q = 0; q < shapes->nqpoints_; ++q)
    for (int k = 0; k < actmat->GetNumberOfInternalStateVariables(); ++k)
      actmat->SetInternalState(k, tempMat2(q, k), q);

  return 0;
}


/*----------------------------------------------------------------------*
 |  Project Material Field for Tet                       hoermann 01/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::ProjectMaterialFieldTet(
    const DRT::Element* ele  //!< the element we are dealing with
)
{
  const Teuchos::RCP<MAT::Myocard>& actmat =
      Teuchos::rcp_dynamic_cast<MAT::Myocard>(ele->Material());

  DRT::ELEMENTS::ScaTraHDG* hdgele =
      dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<DRT::Element*>(ele));

  // polynomial space to get the value of the shape function at the material gauss points
  DRT::UTILS::PolynomialSpaceParams params(distype, hdgele->DegreeOld(), this->usescompletepoly_);
  Teuchos::RCP<DRT::UTILS::PolynomialSpace<probdim>> polySpace =
      DRT::UTILS::PolynomialSpaceCache<probdim>::Instance().Create(params);

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

  const DRT::UTILS::IntPointsAndWeights<DRT::UTILS::DisTypeToDim<distype>::dim> intpoints_old(
      SCATRA::DisTypeToMatGaussRule<distype>::GetGaussRule(degold));
  const DRT::UTILS::IntPointsAndWeights<DRT::UTILS::DisTypeToDim<distype>::dim> intpoints(
      SCATRA::DisTypeToMatGaussRule<distype>::GetGaussRule(deg));


  std::vector<Epetra_SerialDenseVector> shape_gp_old(intpoints_old.IP().nquad);
  std::vector<Epetra_SerialDenseVector> shape_gp(intpoints.IP().nquad);

  // coordinate of material gauss points
  LINALG::Matrix<probdim, 1> mat_gp_coord(true);

  for (int q = 0; q < intpoints_old.IP().nquad; ++q)
  {
    shape_gp_old[q].Size(polySpace->Size());

    // gaussian points coordinates
    for (int idim = 0; idim < DRT::UTILS::DisTypeToDim<distype>::dim; ++idim)
      mat_gp_coord(idim) = intpoints_old.IP().qxg[q][idim];
    polySpace->Evaluate(mat_gp_coord, shape_gp_old[q]);
  }

  for (int q = 0; q < intpoints.IP().nquad; ++q)
  {
    shape_gp[q].Size(polySpace->Size());

    // gaussian points coordinates
    for (int idim = 0; idim < DRT::UTILS::DisTypeToDim<distype>::dim; ++idim)
      mat_gp_coord(idim) = intpoints.IP().qxg[q][idim];
    polySpace->Evaluate(mat_gp_coord, shape_gp[q]);
  }

  this->shapes_->Evaluate(*ele);
  // Jacobian determinant
  double jacdet = this->shapes_->xjm.Determinant();



  Epetra_SerialDenseMatrix massPartOld(polySpace->Size(), shape_gp_old.size());
  Epetra_SerialDenseMatrix massPartOldW(polySpace->Size(), shape_gp_old.size());
  Epetra_SerialDenseMatrix massPart(polySpace->Size(), shape_gp.size());
  Epetra_SerialDenseMatrix massPartW(polySpace->Size(), shape_gp.size());
  Epetra_SerialDenseMatrix Mmat(polySpace->Size(), polySpace->Size());

  Epetra_SerialDenseMatrix state_variables(
      shape_gp_old.size(), actmat->GetNumberOfInternalStateVariables());

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

  Mmat.Multiply('N', 'T', 1.0, massPartOld, massPartOldW, 0.0);

  for (unsigned int q = 0; q < shape_gp_old.size(); ++q)
    for (int k = 0; k < actmat->GetNumberOfInternalStateVariables(); ++k)
      state_variables(q, k) = actmat->GetInternalState(k, q);

  Epetra_SerialDenseMatrix tempMat1(polySpace->Size(), actmat->GetNumberOfInternalStateVariables());
  tempMat1.Multiply('N', 'N', 1.0, massPartOldW, state_variables, 0.0);

  Epetra_SerialDenseSolver inverseMat;
  inverseMat.SetMatrix(Mmat);
  inverseMat.SetVectors(tempMat1, tempMat1);
  inverseMat.FactorWithEquilibration(true);
  int err2 = inverseMat.Factor();
  int err = inverseMat.Solve();
  if (err != 0 || err2 != 0) dserror("Inversion of matrix failed with errorcode %d", err);

  Epetra_SerialDenseMatrix tempMat2(shape_gp.size(), actmat->GetNumberOfInternalStateVariables());
  tempMat2.Multiply('T', 'N', 1.0, massPart, tempMat1, 0.0);

  actmat->SetGP(shape_gp.size());
  actmat->ResizeInternalStateVariables();


  for (unsigned int q = 0; q < shape_gp.size(); ++q)
    for (int k = 0; k < actmat->GetNumberOfInternalStateVariables(); ++k)
      actmat->SetInternalState(k, tempMat2(q, k), q);

  return 0;
}



// template classes
// 1D elements
// template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::line2,1>;
// template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::line2,2>;
// template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::line2,3>;
// template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::line3,1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::tri3>;
// template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::quad9, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::nurbs9, 2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::nurbs27>;
