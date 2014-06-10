/*!----------------------------------------------------------------------
\file drt_utils_shapevalues_hdg.H
\brief Evaluation of scalar shape functions for HDG

<pre>
Maintainer: Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include "drt_utils_shapevalues_hdg.H"
#include "drt_utils_boundary_integration.H"
#include "../drt_geometry/position_array.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                              kronbichler 05/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::UTILS::ShapeValues<distype>::ShapeValues(const unsigned int degree,
                                              const bool         completepoly,
                                              const unsigned int quadratureDegree)
:
degree_ (degree),
polySpace_(distype, degree, completepoly),
polySpaceFace_(DRT::UTILS::DisTypeToFaceShapeType<distype>::shape, degree, completepoly),
ndofs_ (polySpace_.Size()),
nfdofs_ (polySpaceFace_.Size()),
quadrature_ (DRT::UTILS::GaussPointCache::Instance().Create(distype, quadratureDegree)),
fquadrature_ (DRT::UTILS::GaussPointCache::Instance().Create(DRT::UTILS::DisTypeToFaceShapeType<distype>::shape,
                                                             quadratureDegree)),
nqpoints_ (quadrature_->NumPoints()),
nfqpoints_ (fquadrature_->NumPoints())
{
  // TODO: allow polynomials separate for each face, use generic polynomials

  Epetra_SerialDenseVector values(ndofs_);
  Epetra_SerialDenseMatrix derivs(nsd_, ndofs_);
  Epetra_SerialDenseVector faceValues(nfdofs_);

  xyzreal.Shape(nsd_, nqpoints_);
  xyzFreal.Shape(nsd_, nfqpoints_);

  funct.Shape(nen_, nqpoints_);
  functF.Shape(nfn_, nfqpoints_);

  shfunct.Shape(ndofs_,nqpoints_);
  shfunctAvg.Resize(ndofs_);
  shderiv.Shape(ndofs_*nsd_,nqpoints_);
  shderxy.Shape(ndofs_*nsd_,nqpoints_);
  jfac.Resize(nqpoints_);

  for (unsigned int q=0; q<nqpoints_; ++q )
  {
    // gauss point in real coordinates
    const double* gpcoord = quadrature_->Point(q);
    for (unsigned int idim=0;idim<nsd_;idim++)
      xsi(idim) = gpcoord[idim];

    polySpace_.Evaluate(xsi,values);
    polySpace_.Evaluate_deriv1(xsi,derivs);

    for (unsigned int i=0; i<ndofs_; ++i)
    {
      shfunct(i,q) = values(i);
      for (unsigned int d=0; d<nsd_; ++d)
        shderiv(i*nsd_+d,q) = derivs(d,i);
    }

    LINALG::Matrix<nen_,1> myfunct(funct.A()+q*nen_,true);
    DRT::UTILS::shape_function<distype>(xsi,myfunct);
  }

  shfunctFNoPermute.Shape(nfdofs_, nfqpoints_);
  shfunctF.Shape(nfdofs_, nfqpoints_);
  shfunctI.resize(nfaces_);
  for (unsigned int f=0;f<nfaces_; ++f)
    shfunctI[f].Shape(ndofs_, nfqpoints_);
  normals.Shape(nsd_, nfqpoints_);
  jfacF.Resize(nfqpoints_);

  for (unsigned int q=0; q<nfqpoints_; ++q )
  {
    const double* gpcoord = fquadrature_->Point(q);

    const unsigned int codim = nsd_-1;
    for (unsigned int idim=0;idim<codim;idim++)
      xsiF(idim) = gpcoord[idim];

    polySpaceFace_.Evaluate(xsiF,faceValues);
    for (unsigned int i=0; i<nfdofs_; ++i)
      shfunctFNoPermute(i,q) = faceValues(i);

    LINALG::Matrix<nfn_,1> myfunct(functF.A()+q*nfn_,true);
    DRT::UTILS::shape_function<DRT::UTILS::DisTypeToFaceShapeType<distype>::shape>(xsiF,myfunct);
  }

  Epetra_SerialDenseMatrix quadrature(nfqpoints_,nsd_,false);
  Epetra_SerialDenseMatrix trafo(nsd_,nsd_,false);
  for (unsigned int f=0; f<nfaces_; ++f)
  {
    DRT::UTILS::BoundaryGPToParentGP<nsd_>(quadrature,trafo,*fquadrature_,distype,
                                           DRT::UTILS::DisTypeToFaceShapeType<distype>::shape, f);
    for (unsigned int q=0; q<nfqpoints_; ++q)
    {
      for (unsigned int d=0; d<nsd_; ++d)
        xsi(d) = quadrature(q,d);
      polySpace_.Evaluate(xsi,values);
      for (unsigned int i=0; i<ndofs_; ++i)
        shfunctI[f](i,q) = values(i);
    }
  }

  if (nsd_ == 2)
    faceNodeOrder = DRT::UTILS::getEleNodeNumberingLines(distype);
  else if (nsd_ == 3)
    faceNodeOrder = DRT::UTILS::getEleNodeNumberingSurfaces(distype);
  else
    dserror("Not implemented for dim != 2, 3");
}



/*----------------------------------------------------------------------*
 |  Evaluate element-dependent shape data (public)    kronbichler 05/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void
DRT::UTILS::ShapeValues<distype>::Evaluate (const DRT::Element &ele)
{
  dsassert(ele.Shape() == distype, "Internal error");
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(&ele,xyze);

  for (unsigned int i=0; i<ndofs_; ++i)
    shfunctAvg(i) = 0.;
  double faceVol = 0.;

  // evaluate geometry
  for (unsigned int q=0; q<nqpoints_; ++q)
  {
    const double* gpcoord = quadrature_->Point(q);
    for (unsigned int idim=0;idim<nsd_;idim++)
      xsi(idim) = gpcoord[idim];

    DRT::UTILS::shape_function_deriv1<distype>(xsi,deriv);
    xjm.MultiplyNT(deriv,xyze);
    jfac(q) = xji.Invert(xjm) * quadrature_->Weight(q);

    LINALG::Matrix<nen_,1> myfunct(funct.A()+q*nen_,true);
    LINALG::Matrix<nsd_,1> mypoint(xyzreal.A()+q*nsd_,true);
    mypoint.MultiplyNN(xyze,myfunct);

    // transform shape functions
    for (unsigned int i=0; i<ndofs_; ++i)
      for (unsigned int d=0; d<nsd_; ++d)
      {
        shderxy(i*nsd_+d,q) = xji(d,0) * shderiv(i*nsd_,q);
        for (unsigned int e=1; e<nsd_; ++e)
          shderxy(i*nsd_+d,q) += xji(d,e) * shderiv(i*nsd_+e,q);
      }

    for (unsigned int i=0; i<ndofs_; ++i)
      shfunctAvg(i) += shfunct(i,q) * jfac(q);
    faceVol += jfac(q);
  }
  faceVol = 1./faceVol;
  for (unsigned int i=0; i<ndofs_; ++i)
    shfunctAvg(i) *= faceVol;
}



/*----------------------------------------------------------------------*
 |  Evaluate face-dependent shape data (public)       kronbichler 05/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void
DRT::UTILS::ShapeValues<distype>::EvaluateFace (const DRT::Element &ele,
                                                const unsigned int  face)
{
  const DRT::Element::DiscretizationType facedis = DRT::UTILS::DisTypeToFaceShapeType<distype>::shape;

  // get face position array from element position array
  dsassert(faceNodeOrder[face].size() == nfn_,
           "Internal error");
  for (unsigned int i=0; i<nfn_; ++i)
    for (unsigned int d=0; d<nsd_; ++d)
      xyzeF(d,i) = xyze(d,faceNodeOrder[face][i]);

  // evaluate geometry
  for (unsigned int q=0; q<nfqpoints_; ++q) {
    const double* gpcoord = fquadrature_->Point(q);
    for (unsigned int idim=0;idim<nsd_-1;idim++)
      xsiF(idim) = gpcoord[idim];

    DRT::UTILS::shape_function_deriv1<facedis>(xsiF,derivF);
    double jacdet = 0;
    DRT::UTILS::ComputeMetricTensorForBoundaryEle<facedis>(xyzeF,derivF,metricTensor,jacdet,&normal);
    for (unsigned int d=0; d<nsd_; ++d)
      normals(d,q) = normal(d);
    jfacF(q) = jacdet * fquadrature_->Weight(q);

    LINALG::Matrix<nfn_,1> myfunct(functF.A()+q*nfn_,true);
    LINALG::Matrix<nsd_,1> mypoint(xyzFreal.A()+q*nsd_,true);
    mypoint.MultiplyNN(xyzeF,myfunct);
  }

  AdjustFaceOrientation(ele, face);
}



/*----------------------------------------------------------------------*
 |  Reorder evaluated face shape functions (private)  kronbichler 05/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void
DRT::UTILS::ShapeValues<distype>::AdjustFaceOrientation (const DRT::Element &ele,
                                                         const unsigned int  face)
{
  // figure out how to permute face indices by checking permutation of nodes.
  // In case there is some permutation, we need to change the order of quadrature
  // points in the shape functions of the trace in order to get the correct
  // contributions to the matrix.
  const int * nodeIds = ele.NodeIds();
  const int * fnodeIds = ele.Faces()[face]->NodeIds();
  const int nqpoints1d = std::pow(nfqpoints_+0.001,1./(nsd_-1));
  // easy case: standard orientation
  bool standard = true;
  for (unsigned int i=0; i<nfn_; ++i)
    if (nodeIds[faceNodeOrder[face][i]] != fnodeIds[i])
      standard = false;
  if (standard)
  {
    for (unsigned int q=0; q<nfqpoints_; ++q)
      for (unsigned int i=0; i<nfdofs_; ++i)
        shfunctF(i,q) = shfunctFNoPermute(i,q);
  }
  // OK, the orientation is different from what I expect. see if we can find it
  else switch (nsd_)
  {
  case 2:
    // face flipped is the only case
    {
      dsassert(nodeIds[faceNodeOrder[face][1]] == fnodeIds[0] &&
               nodeIds[faceNodeOrder[face][0]] == fnodeIds[1], "Unknown face orientation in 2D");
      for (unsigned int q=0; q<nfqpoints_; ++q)
      {
        for (unsigned int i=0; i<nfdofs_; ++i)
          shfunctF(i,q) = shfunctFNoPermute(i,nfqpoints_-1-q);
      }
    }
    break;
  case 3:
    if (distype == DRT::Element::hex8 ||
        distype == DRT::Element::hex20 ||
        distype == DRT::Element::hex27 ||
        distype == DRT::Element::nurbs8 ||
        distype == DRT::Element::nurbs27)
    {
      if (nodeIds[faceNodeOrder[face][0]] == fnodeIds[1] &&
          nodeIds[faceNodeOrder[face][1]] == fnodeIds[0] &&
          nodeIds[faceNodeOrder[face][2]] == fnodeIds[3] &&
          nodeIds[faceNodeOrder[face][3]] == fnodeIds[2])    // x-direction mirrored
      {
        for (unsigned int q=0; q<nfqpoints_; ++q)
        {
          const int ax = q%nqpoints1d;
          const int ay = q/nqpoints1d;
          int permute = nqpoints1d-1-ax + ay * nqpoints1d;
          for (unsigned int i=0; i<nfdofs_; ++i)
            shfunctF(i,q) = shfunctFNoPermute(i,permute);
        }
      }
      else if (nodeIds[faceNodeOrder[face][0]] == fnodeIds[0] &&
               nodeIds[faceNodeOrder[face][1]] == fnodeIds[3] &&
               nodeIds[faceNodeOrder[face][2]] == fnodeIds[2] &&
               nodeIds[faceNodeOrder[face][3]] == fnodeIds[1])    // permute x and y
      {
        for (unsigned int q=0; q<nfqpoints_; ++q)
        {
          const int ax = q%nqpoints1d;
          const int ay = q/nqpoints1d;
          int permute = ay + ax * nqpoints1d;
          for (unsigned int i=0; i<nfdofs_; ++i)
            shfunctF(i,q) = shfunctFNoPermute(i,permute);
        }
      }
      else if (nodeIds[faceNodeOrder[face][0]] == fnodeIds[3] &&
               nodeIds[faceNodeOrder[face][1]] == fnodeIds[2] &&
               nodeIds[faceNodeOrder[face][2]] == fnodeIds[1] &&
               nodeIds[faceNodeOrder[face][3]] == fnodeIds[0])    // y mirrored
      {
        for (unsigned int q=0; q<nfqpoints_; ++q)
        {
          const int ax = q%nqpoints1d;
          const int ay = q/nqpoints1d;
          int permute = ax + (nqpoints1d-1-ay) * nqpoints1d;
          for (unsigned int i=0; i<nfdofs_; ++i)
            shfunctF(i,q) = shfunctFNoPermute(i,permute);
        }
      }
      else if (nodeIds[faceNodeOrder[face][0]] == fnodeIds[2] &&
               nodeIds[faceNodeOrder[face][1]] == fnodeIds[3] &&
               nodeIds[faceNodeOrder[face][2]] == fnodeIds[0] &&
               nodeIds[faceNodeOrder[face][3]] == fnodeIds[1])    // x and y mirrored
      {
        for (unsigned int q=0; q<nfqpoints_; ++q)
        {
          const int ax = q%nqpoints1d;
          const int ay = q/nqpoints1d;
          int permute = (nqpoints1d-1-ax) + (nqpoints1d-1-ay) * nqpoints1d;
          for (unsigned int i=0; i<nfdofs_; ++i)
            shfunctF(i,q) = shfunctFNoPermute(i,permute);
        }
      }
      else if (nodeIds[faceNodeOrder[face][0]] == fnodeIds[2] &&
               nodeIds[faceNodeOrder[face][1]] == fnodeIds[1] &&
               nodeIds[faceNodeOrder[face][2]] == fnodeIds[0] &&
               nodeIds[faceNodeOrder[face][3]] == fnodeIds[3])    // x and y mirrored and permuted
      {
        for (unsigned int q=0; q<nfqpoints_; ++q)
        {
          const int ax = q%nqpoints1d;
          const int ay = q/nqpoints1d;
          int permute = (nqpoints1d-1-ay) + (nqpoints1d-1-ax) * nqpoints1d; // note that this is for lexicographic ordering
          for (unsigned int i=0; i<nfdofs_; ++i)
            shfunctF(i,q) = shfunctFNoPermute(i,permute);
        }
      }
      else
      {
        for (unsigned int i=0; i<4; ++i)
          std::cout << nodeIds[faceNodeOrder[face][i]] << " " << fnodeIds[i] << "   ";
        std::cout << std::endl << std::flush;
        dserror("Unknown HEX face orientation in 3D");
      }
    }
    else if (distype == DRT::Element::tet4 ||
             distype == DRT::Element::tet10)
    {
      if (nodeIds[faceNodeOrder[face][0]] == fnodeIds[1] &&
          nodeIds[faceNodeOrder[face][1]] == fnodeIds[0] &&
          nodeIds[faceNodeOrder[face][2]] == fnodeIds[2])
      {
        // zeroth and first node permuted
        for (unsigned int q=0; q<nfqpoints_; ++q)
        {
          double point[2] = {fquadrature_->Point(q)[0], fquadrature_->Point(q)[1]};
          point[0] = 1.-point[1]-point[0];
          unsigned int r=0;
          for (; r<nfqpoints_; ++r)
            if (std::abs(fquadrature_->Point(r)[0]-point[0]) < 1e-14 &&
                std::abs(fquadrature_->Point(r)[1]-point[1]) < 1e-14)
              break;
          dsassert(r<nfqpoints_, "Quadrature points seem to not be symmetric");
          for (unsigned int i=0; i<nfdofs_; ++i)
            shfunctF(i,q) = shfunctFNoPermute(i,r);
        }
      }
      else if (nodeIds[faceNodeOrder[face][0]] == fnodeIds[2] &&
               nodeIds[faceNodeOrder[face][1]] == fnodeIds[1] &&
               nodeIds[faceNodeOrder[face][2]] == fnodeIds[0])
      {
        // zeroth and second node permuted
        for (unsigned int q=0; q<nfqpoints_; ++q)
        {
          double point[2] = {fquadrature_->Point(q)[0], fquadrature_->Point(q)[1]};
          point[1] = 1.-point[1]-point[0];
          unsigned int r=0;
          for (; r<nfqpoints_; ++r)
            if (std::abs(fquadrature_->Point(r)[0]-point[0]) < 1e-14 &&
                std::abs(fquadrature_->Point(r)[1]-point[1]) < 1e-14)
              break;
          dsassert(r<nfqpoints_, "Quadrature points seem to not be symmetric");
          for (unsigned int i=0; i<nfdofs_; ++i)
            shfunctF(i,q) = shfunctFNoPermute(i,r);
        }
      }
      else if (nodeIds[faceNodeOrder[face][0]] == fnodeIds[0] &&
               nodeIds[faceNodeOrder[face][1]] == fnodeIds[2] &&
               nodeIds[faceNodeOrder[face][2]] == fnodeIds[1])
      {
        // first and second node permuted
        for (unsigned int q=0; q<nfqpoints_; ++q)
        {
          double point[2] = {fquadrature_->Point(q)[1], fquadrature_->Point(q)[0]};
          unsigned int r=0;
          for (; r<nfqpoints_; ++r)
            if (std::abs(fquadrature_->Point(r)[0]-point[0]) < 1e-14 &&
                std::abs(fquadrature_->Point(r)[1]-point[1]) < 1e-14)
              break;
          dsassert(r<nfqpoints_, "Quadrature points seem to not be symmetric");
          for (unsigned int i=0; i<nfdofs_; ++i)
            shfunctF(i,q) = shfunctFNoPermute(i,r);
        }
      }
      else
      {
        for (unsigned int i=0; i<3; ++i)
          std::cout << nodeIds[faceNodeOrder[face][i]] << " " << fnodeIds[i] << "   ";
        std::cout << std::endl << std::flush;
        // need to transform quadrature point coordinate 0 into 1-p[0]-p[1]
        dserror("Unknown TET face orientation in 3D");
      }
    }
    else
      dserror("Shape type %s not yet implemented", (DRT::DistypeToString(distype)).c_str());
    break;
  default:
    dserror("Only implemented in 2D and 3D");
    break;
  }
}


// explicit instantiation of template classes
template struct DRT::UTILS::ShapeValues<DRT::Element::hex8>;
template struct DRT::UTILS::ShapeValues<DRT::Element::hex20>;
template struct DRT::UTILS::ShapeValues<DRT::Element::hex27>;
template struct DRT::UTILS::ShapeValues<DRT::Element::tet4>;
template struct DRT::UTILS::ShapeValues<DRT::Element::tet10>;
template struct DRT::UTILS::ShapeValues<DRT::Element::wedge6>;
template struct DRT::UTILS::ShapeValues<DRT::Element::pyramid5>;
template struct DRT::UTILS::ShapeValues<DRT::Element::quad4>;
template struct DRT::UTILS::ShapeValues<DRT::Element::quad8>;
template struct DRT::UTILS::ShapeValues<DRT::Element::quad9>;
template struct DRT::UTILS::ShapeValues<DRT::Element::tri3>;
template struct DRT::UTILS::ShapeValues<DRT::Element::tri6>;
template struct DRT::UTILS::ShapeValues<DRT::Element::nurbs9>;
template struct DRT::UTILS::ShapeValues<DRT::Element::nurbs27>;
