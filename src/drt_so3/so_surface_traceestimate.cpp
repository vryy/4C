/*!----------------------------------------------------------------------
\file so_surface_traceestimate.cpp

\brief evaluation for nitsche trace inequality estimate

\maintainer Alexander Seitz

\level 3
*----------------------------------------------------------------------*/

#include "so_surface.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_element_integration_select.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/so3_material.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_mat/material_service.H"
#include "../drt_mat/fourieriso.H"


/*----------------------------------------------------------------------*
 |                                                           seitz 11/16|
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::StructuralSurface::EstimateNitscheTraceMaxEigenvalueCombined(
    std::vector<double>& parent_disp)
{
  switch(ParentElement()->Shape())
  {
  case DRT::Element::hex8:
    if (Shape()==DRT::Element::quad4)
      return EstimateNitscheTraceMaxEigenvalueCombined<DRT::Element::hex8,DRT::Element::quad4>(parent_disp);
    else
      dserror("how can an hex8 element have a surface that is not quad4 ???");
    break;
  case DRT::Element::hex27:
      return EstimateNitscheTraceMaxEigenvalueCombined<DRT::Element::hex27,DRT::Element::quad9>(parent_disp);
      break;
  case DRT::Element::tet4:
      return EstimateNitscheTraceMaxEigenvalueCombined<DRT::Element::tet4,DRT::Element::tri3>(parent_disp);
      break;
  default:
    dserror("parent shape not implemented");
  }

  return 0;
}


template<
DRT::Element::DiscretizationType dt_vol,
DRT::Element::DiscretizationType dt_surf>
double DRT::ELEMENTS::StructuralSurface::EstimateNitscheTraceMaxEigenvalueCombined(
    std::vector<double>& parent_disp)
{
  const int dim = DRT::UTILS::DisTypeToDim<dt_vol>::dim;
  const int num_dof =
      DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement*DRT::UTILS::DisTypeToDim<dt_vol>::dim;
  const int dim_image =
      DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement*DRT::UTILS::DisTypeToDim<dt_vol>::dim-
      DRT::UTILS::DisTypeToDim<dt_vol>::dim*(DRT::UTILS::DisTypeToDim<dt_vol>::dim+1)/2;

  LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement,3> xrefe;
  LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement,3> xcurr;

  for (int i=0;i<ParentElement()->NumNode();++i)
    for (int d=0;d<dim;++d)
    {
      xrefe(i,d) = ParentElement()->Nodes()[i]->X()[d];
      xcurr(i,d) = xrefe(i,d) + parent_disp[i*dim+d];
    }

  LINALG::Matrix<num_dof,num_dof> vol, surf;

  TraceEstimateVolMatrix<dt_vol>(xrefe,xcurr,vol);
  TraceEstimateSurfMatrix<dt_vol,dt_surf>(xrefe,xcurr,surf);

  LINALG::Matrix<num_dof,dim_image> proj, tmp;
  SubspaceProjector<dt_vol>(xcurr,proj);

  LINALG::Matrix<dim_image,dim_image> vol_red,surf_red;

  tmp.Multiply(vol,proj);
  vol_red.MultiplyTN(proj,tmp);
  tmp.Multiply(surf,proj);
  surf_red.MultiplyTN(proj,tmp);

  Epetra_SerialDenseMatrix vol_red_sd(::View,vol_red.A(),dim_image,dim_image,dim_image);
  Epetra_SerialDenseMatrix surf_red_sd(::View,surf_red.A(),dim_image,dim_image,dim_image);

  return LINALG::GeneralizedEigen(surf_red_sd,vol_red_sd);
}

template<DRT::Element::DiscretizationType dt_vol>
void DRT::ELEMENTS::StructuralSurface::TraceEstimateVolMatrix(
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement,3>& xrefe,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement,3>& xcurr,
    LINALG::Matrix<
    DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement*3,
    DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement*3>& vol)
{
  const int dim = DRT::UTILS::DisTypeToDim<dt_vol>::dim;

  double jac;
  LINALG::Matrix<3,3> defgrd;
  LINALG::Matrix<3,3> rcg;
  LINALG::Matrix<6,1> glstrain;
  LINALG::Matrix<6,DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement*3> bop;
  LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement*3,6> bc;
  LINALG::Matrix<dim,DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement> N_XYZ;

  DRT::UTILS::IntPointsAndWeights<dim> ip(DRT::ELEMENTS::DisTypeToOptGaussRule<dt_vol>::rule);

  for (int gp=0;gp<ip.IP().nquad;++gp)
  {
    const LINALG::Matrix<3,1> xi(ip.IP().qxg[gp],false);
    Strains<dt_vol>(xrefe,xcurr,xi,jac,defgrd,glstrain,rcg,bop,N_XYZ);

    LINALG::Matrix<6,6> cmat(true);
    LINALG::Matrix<6,1> stress(true);
    Teuchos::ParameterList params;
    Teuchos::rcp_dynamic_cast<MAT::So3Material>(ParentElement()->Material())->
        Evaluate(&defgrd,&glstrain,params,&stress,&cmat,ParentElement()->Id());
    bc.MultiplyTN(bop,cmat);
    vol.Multiply(ip.IP().qwgt[gp]*jac,bc,bop,1.);

  }

  return;
}


template<
DRT::Element::DiscretizationType dt_vol,
DRT::Element::DiscretizationType dt_surf>
void DRT::ELEMENTS::StructuralSurface::TraceEstimateSurfMatrix(
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement,3>& xrefe,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement,3>& xcurr,
    LINALG::Matrix<
    DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement*3,
    DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement*3>& surf)
{
  const int dim = DRT::UTILS::DisTypeToDim<dt_vol>::dim;

  LINALG::Matrix<6,6> id4;
  for (int i=0;i<3;++i)id4(i,i)=1.;
  for (int i=3;i<6;++i)id4(i,i)=2.;

  LINALG::SerialDenseMatrix xrefe_surf(DRT::UTILS::DisTypeToNumNodePerEle<dt_surf>::numNodePerElement,dim);
  MaterialConfiguration(xrefe_surf);

  std::vector<double> n(3);
  LINALG::Matrix<3,1> n_v(&n[0],true);
  LINALG::Matrix<3,3> nn;
  double detA;
  double jac;
  LINALG::Matrix<3,3> defgrd;
  LINALG::Matrix<3,3> rcg;
  LINALG::Matrix<6,1> glstrain;
  LINALG::Matrix<6,DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement*3> bop;
  LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement*3,6> bc;
  LINALG::Matrix<dim,DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement> N_XYZ;

  DRT::UTILS::IntPointsAndWeights<dim-1> ip(DRT::ELEMENTS::DisTypeToOptGaussRule<dt_surf>::rule);
  LINALG::SerialDenseMatrix  deriv_surf(2,DRT::UTILS::DisTypeToNumNodePerEle<dt_surf>::numNodePerElement);

  for (int gp=0;gp<ip.IP().nquad;++gp)
  {
    DRT::UTILS::CollectedGaussPoints intpoints = DRT::UTILS::CollectedGaussPoints(1); //reserve just for 1 entry ...
      intpoints.Append(ip.IP().qxg[gp][0], ip.IP().qxg[gp][1],0.0, ip.IP().qwgt[gp]);

    // get coordinates of gauss point w.r.t. local parent coordinate system
    LINALG::SerialDenseMatrix pqxg(1,3);
    LINALG::Matrix<3,3> derivtrafo;

    DRT::UTILS::BoundaryGPToParentGP<3>( pqxg,
        derivtrafo,
        intpoints ,
        ParentElement()->Shape(),
        Shape(),
        FaceParentNumber());

    LINALG::Matrix<3,1> xi; for (int i=0;i<3;++i) xi(i)=pqxg(0,i);
    Strains<dt_vol>(xrefe,xcurr,xi,jac,defgrd,glstrain,rcg,bop,N_XYZ);

    LINALG::Matrix<6,6> cmat(true);
    LINALG::Matrix<6,1> stress(true);
    Teuchos::ParameterList params;
    Teuchos::rcp_dynamic_cast<MAT::So3Material>(ParentElement()->Material())->
        Evaluate(&defgrd,&glstrain,params,&stress,&cmat,ParentElement()->Id());

    DRT::UTILS::shape_function_2D_deriv1(deriv_surf,ip.IP().qxg[gp][0], ip.IP().qxg[gp][1],Shape());
    SurfaceIntegration(detA,n,xrefe_surf,deriv_surf);
    n_v.Scale(1./n_v.Norm2());
    nn.MultiplyNT(n_v,n_v);

    LINALG::Matrix<6,6> cn;
    MAT::AddSymmetricHolzapfelProduct(cn,rcg,nn,.25);

    LINALG::Matrix<6,6>tmp1,tmp2;
    tmp1.Multiply(cmat,id4);
    tmp2.Multiply(tmp1,cn);
    tmp1.Multiply(tmp2,id4);
    tmp2.Multiply(tmp1,cmat);

    LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement*3,6> tmp3;
    tmp3.MultiplyTN(bop,tmp2);

    surf.Multiply(detA*ip.IP().qwgt[gp],tmp3,bop,1.);
  }

  return;
}

template<DRT::Element::DiscretizationType dt_vol>
void DRT::ELEMENTS::StructuralSurface::Strains(
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement,3>& xrefe,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement,3>& xcurr,
    const LINALG::Matrix<3,1>& xi,
    double& jac,
    LINALG::Matrix<3,3>& defgrd,
    LINALG::Matrix<6,1>& glstrain,
    LINALG::Matrix<3,3>& rcg,
    LINALG::Matrix<6,DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement*3>& bop,
    LINALG::Matrix<3,DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement>& N_XYZ)
{
  const int dim = DRT::UTILS::DisTypeToDim<dt_vol>::dim;
  const int num_node = DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement;
  LINALG::Matrix<dim,num_node> deriv;
  DRT::UTILS::shape_function_deriv1<dt_vol>(xi,deriv);

  LINALG::Matrix<dim,dim> invJ;
  invJ.Multiply(deriv,xrefe);
  jac=invJ.Invert();
  N_XYZ.Multiply(invJ,deriv);
  defgrd.MultiplyTT(xcurr,N_XYZ);

  rcg.MultiplyTN(defgrd,defgrd);
  glstrain(0) = 0.5 * (rcg(0,0) - 1.0);
  glstrain(1) = 0.5 * (rcg(1,1) - 1.0);
  glstrain(2) = 0.5 * (rcg(2,2) - 1.0);
  glstrain(3) = rcg(0,1);
  glstrain(4) = rcg(1,2);
  glstrain(5) = rcg(2,0);

  for (int i=0; i<num_node; ++i)
  {
    bop(0,dim*i+0) = defgrd(0,0)*N_XYZ(0,i);
    bop(0,dim*i+1) = defgrd(1,0)*N_XYZ(0,i);
    bop(0,dim*i+2) = defgrd(2,0)*N_XYZ(0,i);
    bop(1,dim*i+0) = defgrd(0,1)*N_XYZ(1,i);
    bop(1,dim*i+1) = defgrd(1,1)*N_XYZ(1,i);
    bop(1,dim*i+2) = defgrd(2,1)*N_XYZ(1,i);
    bop(2,dim*i+0) = defgrd(0,2)*N_XYZ(2,i);
    bop(2,dim*i+1) = defgrd(1,2)*N_XYZ(2,i);
    bop(2,dim*i+2) = defgrd(2,2)*N_XYZ(2,i);
    /* ~~~ */
    bop(3,dim*i+0) = defgrd(0,0)*N_XYZ(1,i) + defgrd(0,1)*N_XYZ(0,i);
    bop(3,dim*i+1) = defgrd(1,0)*N_XYZ(1,i) + defgrd(1,1)*N_XYZ(0,i);
    bop(3,dim*i+2) = defgrd(2,0)*N_XYZ(1,i) + defgrd(2,1)*N_XYZ(0,i);
    bop(4,dim*i+0) = defgrd(0,1)*N_XYZ(2,i) + defgrd(0,2)*N_XYZ(1,i);
    bop(4,dim*i+1) = defgrd(1,1)*N_XYZ(2,i) + defgrd(1,2)*N_XYZ(1,i);
    bop(4,dim*i+2) = defgrd(2,1)*N_XYZ(2,i) + defgrd(2,2)*N_XYZ(1,i);
    bop(5,dim*i+0) = defgrd(0,2)*N_XYZ(0,i) + defgrd(0,0)*N_XYZ(2,i);
    bop(5,dim*i+1) = defgrd(1,2)*N_XYZ(0,i) + defgrd(1,0)*N_XYZ(2,i);
    bop(5,dim*i+2) = defgrd(2,2)*N_XYZ(0,i) + defgrd(2,0)*N_XYZ(2,i);
  }

  return;
}


template<DRT::Element::DiscretizationType dt_vol>
void DRT::ELEMENTS::StructuralSurface::SubspaceProjector(
      const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement,3>& xcurr,
      LINALG::Matrix<
      DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement*DRT::UTILS::DisTypeToDim<dt_vol>::dim,
      DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement*DRT::UTILS::DisTypeToDim<dt_vol>::dim-
      DRT::UTILS::DisTypeToDim<dt_vol>::dim*(DRT::UTILS::DisTypeToDim<dt_vol>::dim+1)/2>& proj)
{
  const int dim = DRT::UTILS::DisTypeToDim<dt_vol>::dim;
  const int num_node = DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement;
  if (dim!=3)
    dserror("this should be 3D");

  LINALG::Matrix<3,1> c;
  for (int r=0;r<(int)xcurr.M();++r)
    for (int d=0;d<(int)xcurr.N();++d)
      c(d)+=xcurr(r,d);
  c.Scale(1./xcurr.M());

  LINALG::Matrix<dim,1> r[3];
  for (int i=0;i<3;++i)
    r[i](i)=1.;

  // basis, where the first six entries are the rigid body modes and the
  // remaining are constructed to be orthogonal to the rigid body modes
  LINALG::Matrix<dim*num_node,1> basis[dim*num_node];

  // rigid body translations
  for (int i=0;i<dim;++i)
    for (int j=0;j<num_node;++j)
      basis[i](j*dim+i)=1.;

  // rigid body rotations
  for (int i=0;i<dim;++i)
    for (int j=0;j<num_node;++j)
    {
      LINALG::Matrix<3,1> x;
      for (int d=0;d<3;++d) x(d)=xcurr(j,d);
      x.Update(-1.,c,1.);
      LINALG::Matrix<3,1> cross;
      cross.CrossProduct(r[i],x);
      for (int k=0;k<3;++k)
        basis[i+3](j*3+k)=cross(k);
    }
  for (int i=0;i<6;++i)
    basis[i].Scale(1./basis[i].Norm2());

  // build the remaining basis vectors by generalized cross products
  for (int i=6;i<dim*num_node;++i)
  {
    double sign=+1.;
    int off=0;
    bool new_basis_found=false;
    for (off=0;(off<dim*num_node-i) && !new_basis_found;++off)
    {
      for (int j=0;j<i+1;++j)
      {
        LINALG::SerialDenseMatrix det(i,i,true);
        for (int c=0;c<i;++c)
        {
          for (int k=0;k<j;++k)
            det(k,c)=basis[c](k+off);
          for (int k=j;k<i;++k)
            det(k,c)=basis[c](k+1+off);
        }
        basis[i](j+off)=LINALG::DeterminantLU(det)*sign;
        sign*=-1.;
      }
      if (basis[i].Norm2()>1.e-6)
      {
        basis[i].Scale(1./basis[i].Norm2());
        new_basis_found=true;
      }
    }
    if (!new_basis_found)
      dserror("no new basis vector found");
  }

  // at this point basis should already contain an ONB.
  // due to cut-off errors we do another sweep of Gram-Schmidt
  for (int i=0;i<dim*num_node;++i)
  {
    const LINALG::Matrix<dim*num_node,1> tmp(basis[i]);
    for (int j=0;j<i;++j)
      basis[i].Update(-tmp.Dot(basis[j]),basis[j],1.);

    basis[i].Scale(1./basis[i].Norm2());
  }

  // hand out the projection matrix, i.e. the ONB not containing rigid body modes
  for (int i=0;i<dim*num_node;++i)
    for (int j=6;j<dim*num_node;++j)
      proj(i,j-6)=basis[j](i);
}



/*----------------------------------------------------------------------*
 |                                                           seitz 11/16|
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::StructuralSurface::EstimateNitscheTraceMaxEigenvalueTSI(
    std::vector<double>& parent_disp)
{
  switch(ParentElement()->Shape())
  {
  case DRT::Element::hex8:
    if (Shape()==DRT::Element::quad4)
      return EstimateNitscheTraceMaxEigenvalueTSI<DRT::Element::hex8,DRT::Element::quad4>(parent_disp);
    else
      dserror("how can an hex8 element have a surface that is not quad4 ???");
    break;
  case DRT::Element::hex27:
    return EstimateNitscheTraceMaxEigenvalueTSI<DRT::Element::hex27,DRT::Element::quad9>(parent_disp);
  case DRT::Element::tet4:
    return EstimateNitscheTraceMaxEigenvalueTSI<DRT::Element::tet4,DRT::Element::tri3>(parent_disp);
  default:
    dserror("parent shape not implemented");
  }

  return 0;
}

template<
DRT::Element::DiscretizationType dt_vol,
DRT::Element::DiscretizationType dt_surf>
double DRT::ELEMENTS::StructuralSurface::EstimateNitscheTraceMaxEigenvalueTSI(
    std::vector<double>& parent_disp)
{
  const int dim = DRT::UTILS::DisTypeToDim<dt_vol>::dim;
  const int num_dof =
      DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement;
  const int dim_image =
      DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement-1;

  LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement,3> xrefe;
  LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement,3> xcurr;

  for (int i=0;i<ParentElement()->NumNode();++i)
    for (int d=0;d<dim;++d)
    {
      xrefe(i,d) = ParentElement()->Nodes()[i]->X()[d];
      xcurr(i,d) = xrefe(i,d) + parent_disp[i*dim+d];
    }

  LINALG::Matrix<num_dof,num_dof> vol, surf;

  TraceEstimateVolMatrixTSI<dt_vol>(xrefe,xcurr,vol);
  TraceEstimateSurfMatrixTSI<dt_vol,dt_surf>(xrefe,xcurr,surf);


  LINALG::Matrix<num_dof,dim_image> proj, tmp;
  SubspaceProjectorScalar<dt_vol>(proj);

  LINALG::Matrix<dim_image,dim_image> vol_red,surf_red;

  tmp.Multiply(vol,proj);
  vol_red.MultiplyTN(proj,tmp);
  tmp.Multiply(surf,proj);
  surf_red.MultiplyTN(proj,tmp);

  Epetra_SerialDenseMatrix vol_red_sd(::View,vol_red.A(),dim_image,dim_image,dim_image);
  Epetra_SerialDenseMatrix surf_red_sd(::View,surf_red.A(),dim_image,dim_image,dim_image);

  return LINALG::GeneralizedEigen(surf_red_sd,vol_red_sd);
}

template<DRT::Element::DiscretizationType dt_vol>
void DRT::ELEMENTS::StructuralSurface::TraceEstimateVolMatrixTSI(
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement,3>& xrefe,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement,3>& xcurr,
    LINALG::Matrix<
    DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement,
    DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement>& vol)
{
  const int dim = DRT::UTILS::DisTypeToDim<dt_vol>::dim;
  const int num_node = DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement;

  double jac;
  LINALG::Matrix<3,3> defgrd;
  LINALG::Matrix<3,3> rcg;
  LINALG::Matrix<6,1> glstrain;
  LINALG::Matrix<6,DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement*3> bop;
  LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement*3,6> bc;
  LINALG::Matrix<dim,num_node> N_XYZ, iC_N_XYZ;

  DRT::UTILS::IntPointsAndWeights<dim> ip(DRT::ELEMENTS::DisTypeToOptGaussRule<dt_vol>::rule);

  if (ParentElement()->NumMaterial()<2)
    dserror("where's my second material");
  Teuchos::RCP<MAT::FourierIso> mat_thr = Teuchos::rcp_dynamic_cast<MAT::FourierIso>(ParentElement()->Material(1),true);
  const double k0=mat_thr->Conductivity();

  for (int gp=0;gp<ip.IP().nquad;++gp)
  {
    const LINALG::Matrix<3,1> xi(ip.IP().qxg[gp],false);
    Strains<dt_vol>(xrefe,xcurr,xi,jac,defgrd,glstrain,rcg,bop,N_XYZ);

    LINALG::Matrix<3,3> iC;
    iC.MultiplyTN(defgrd,defgrd);
    iC.Invert();

    iC_N_XYZ.Multiply(iC,N_XYZ);
    iC_N_XYZ.Scale(k0);

    vol.MultiplyTN(ip.IP().qwgt[gp]*jac,N_XYZ,iC_N_XYZ,1.);
  }
}


template<
DRT::Element::DiscretizationType dt_vol,
DRT::Element::DiscretizationType dt_surf>
void DRT::ELEMENTS::StructuralSurface::TraceEstimateSurfMatrixTSI(
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement,3>& xrefe,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement,3>& xcurr,
    LINALG::Matrix<
    DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement,
    DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement>& surf)
{
  const int dim = DRT::UTILS::DisTypeToDim<dt_vol>::dim;
  const int num_node = DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement;

  double jac;
  LINALG::Matrix<3,3> defgrd;
  LINALG::Matrix<3,3> rcg;
  LINALG::Matrix<6,1> glstrain;
  LINALG::Matrix<6,DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement*3> bop;
  LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement*3,6> bc;
  LINALG::Matrix<dim,num_node> N_XYZ;
  LINALG::Matrix<1,num_node> iCn_N_XYZ;

  LINALG::SerialDenseMatrix xrefe_surf(DRT::UTILS::DisTypeToNumNodePerEle<dt_surf>::numNodePerElement,dim);
  MaterialConfiguration(xrefe_surf);

  std::vector<double> n(3);
  LINALG::Matrix<3,1> n_v(&n[0],true), iCn;
  double detA;

  DRT::UTILS::IntPointsAndWeights<dim-1> ip(DRT::ELEMENTS::DisTypeToOptGaussRule<dt_surf>::rule);
  LINALG::SerialDenseMatrix  deriv_surf(2,DRT::UTILS::DisTypeToNumNodePerEle<dt_surf>::numNodePerElement);

  if (ParentElement()->NumMaterial()<2)
    dserror("where's my second material");
  Teuchos::RCP<MAT::FourierIso> mat_thr = Teuchos::rcp_dynamic_cast<MAT::FourierIso>(ParentElement()->Material(1),true);
  const double k0=mat_thr->Conductivity();

  for (int gp=0;gp<ip.IP().nquad;++gp)
  {
    DRT::UTILS::shape_function_2D_deriv1(deriv_surf,ip.IP().qxg[gp][0], ip.IP().qxg[gp][1],Shape());
    SurfaceIntegration(detA,n,xrefe_surf,deriv_surf);
    n_v.Scale(1./n_v.Norm2());

    DRT::UTILS::CollectedGaussPoints intpoints = DRT::UTILS::CollectedGaussPoints(1); //reserve just for 1 entry ...
      intpoints.Append(ip.IP().qxg[gp][0], ip.IP().qxg[gp][1],0.0, ip.IP().qwgt[gp]);

    // get coordinates of gauss point w.r.t. local parent coordinate system
    LINALG::SerialDenseMatrix pqxg(1,3);
    LINALG::Matrix<3,3> derivtrafo;

    DRT::UTILS::BoundaryGPToParentGP<3>( pqxg,
        derivtrafo,
        intpoints ,
        ParentElement()->Shape(),
        Shape(),
        FaceParentNumber());

    LINALG::Matrix<3,1> xi; for (int i=0;i<3;++i) xi(i)=pqxg(0,i);

    Strains<dt_vol>(xrefe,xcurr,xi,jac,defgrd,glstrain,rcg,bop,N_XYZ);

    LINALG::Matrix<3,3> iC;
    iC.MultiplyTN(defgrd,defgrd);
    iC.Invert();
    iCn.Multiply(iC,n_v);

    iCn_N_XYZ.MultiplyTN(iCn,N_XYZ);
    iCn_N_XYZ.Scale(k0);

    surf.MultiplyTN(detA*ip.IP().qwgt[gp],iCn_N_XYZ,iCn_N_XYZ,1.);
  }
}



template<DRT::Element::DiscretizationType dt_vol>
void DRT::ELEMENTS::StructuralSurface::SubspaceProjectorScalar(
      LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement,
        DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement-1>& proj)
{
  const int num_node = DRT::UTILS::DisTypeToNumNodePerEle<dt_vol>::numNodePerElement;
  LINALG::Matrix<num_node,1> basis[num_node];

  for (int i=0;i<num_node;++i)
    basis[0](i)=1.;

  for (int i=1;i<num_node;++i)
  {
    double sign=+1.;
    for (int j=0;j<i+1;++j)
    {
      LINALG::SerialDenseMatrix det(i,i,true);
      for (int c=0;c<i;++c)
      {
        for (int k=0;k<j;++k)
          det(k,c)=basis[c](k);
        for (int k=j;k<i;++k)
          det(k,c)=basis[c](k+1);
      }
      basis[i](j)=LINALG::DeterminantLU(det)*sign;
      sign*=-1.;
    }
    basis[i].Scale(1./basis[i].Norm2());
  }

  // hand out the projection matrix, i.e. the ONB not containing rigid body modes
  for (int i=0;i<num_node;++i)
    for (int j=1;j<num_node;++j)
      proj(i,j-1)=basis[j](i);
}
