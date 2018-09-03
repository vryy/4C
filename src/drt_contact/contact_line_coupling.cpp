/*---------------------------------------------------------------------*/
/*!
\file contact_line_coupling.cpp

\brief A class to perform line clipping for line to surface contact +
       call for numerical integration

\level 2

\maintainer Alexander Popp

*/
/*---------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  header                                                   farah 07/16|
 *----------------------------------------------------------------------*/
#include "contact_line_coupling.H"
#include "contact_node.H"
#include "contact_element.H"
#include "contact_defines.H"
#include "contact_integrator.H"
#include "contact_integrator_factory.H"
#include "friction_node.H"

#include "../drt_mortar/mortar_defines.H"
#include "../drt_mortar/mortar_projector.H"

#include "../drt_inpar/inpar_contact.H"

#include "../drt_lib/drt_discret.H"

#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 |  ctor for lts/stl (public)                                farah 07/16|
 *----------------------------------------------------------------------*/
CONTACT::LineToSurfaceCoupling3d::LineToSurfaceCoupling3d(DRT::Discretization& idiscret, int dim,
    Teuchos::ParameterList& params, CoElement& pEle, Teuchos::RCP<MORTAR::MortarElement>& lEle,
    std::vector<CoElement*> surfEles, LineToSurfaceCoupling3d::intType type)
    : idiscret_(idiscret),
      dim_(dim),
      pEle_(pEle),
      lEle_(lEle),
      surfEles_(surfEles),
      currEle_(-1),
      imortar_(params),
      intType_(type)
{
  // empty constructor

  return;
}

/*----------------------------------------------------------------------*
 |  eval (public)                                            farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineToSurfaceCoupling3d::EvaluateCoupling()
{
  // clear entries of master vertices
  DoneBefore().clear();

  // loop over all found master elements
  for (int nele = 0; nele < NumberSurfaceElements(); ++nele)
  {
    // set internal counter
    CurrEle() = nele;

    // 1. init internal data
    Initialize();

    // 2. create aux plane for master ele
    AuxiliaryPlane();  //--> build everything based on line element

    // 3. create aux line for slave ele
    AuxiliaryLine();

    // 4. check orientation
    if (!CheckOrientation()) return;

    // 5. project master nodes onto auxplane
    ProjectMaster();

    // 6. project slave line elements onto auxplane
    ProjectSlave();

    // 7. perform line clipping
    LineClipping();

    // 8. intersections found?
    if ((int)InterSections().size() == 0 or (int) InterSections().size() == 1) continue;

    // 9. check length of Integration Line
    bool check = CheckLength();
    if (check == false) continue;

    // create empty lin vector
    std::vector<std::vector<GEN::pairedvector<int, double>>> linvertex(
        2, std::vector<GEN::pairedvector<int, double>>(
               3, 3 * LineElement()->NumNode() + 3 * SurfaceElement().NumNode() + linsize_));

    // 10. linearize vertices
    LinearizeVertices(linvertex);

    // 11. create intlines
    CreateIntegrationLines(linvertex);

    // 12. consistent dual shape
    ConsistDualShape();

    // 13. integration
    IntegrateLine();
  }  // end loop

  return;
}

/*----------------------------------------------------------------------*
 |  init internal variables                                  farah 08/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineToSurfaceCoupling3d::Initialize()
{
  // reset auxplane normal, center and length
  Auxn()[0] = 0.0;
  Auxn()[1] = 0.0;
  Auxn()[2] = 0.0;

  Auxc()[0] = 0.0;
  Auxc()[1] = 0.0;
  Auxc()[2] = 0.0;

  Lauxn() = 0.0;
  GetDerivAuxn().clear();
  GetDerivAuxc().clear();

  // reset normal of aux line
  //  AuxnLine()[0] = 0.0;
  //  AuxnLine()[1] = 0.0;
  //  AuxnLine()[2] = 0.0;
  //  GetDerivAuxnLine().clear();

  // clear all slave and master vertices
  SlaveVertices().clear();
  MasterVertices().clear();

  // clear previously found intersections
  InterSections().clear();
  TempInterSections().clear();

  // clear integration line
  IntLine() = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------*
 |  check orientation of line and surface element            farah 07/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineToSurfaceCoupling3d::CheckOrientation()
{
  // check if surface normal and line ele are parallel!

  // tolerance for line clipping
  const double sminedge = ParentElement().MinEdgeSize();
  const double mminedge = SurfaceElement().MinEdgeSize();
  const double tol = 0.001 * std::min(sminedge, mminedge);

  // -------------------------------------------
  // CHECK LINE TO SURFACE ORIENTATION!
  // calculate line ele vector
  double lvec[3] = {0.0, 0.0, 0.0};
  CoNode* ns1 = dynamic_cast<CoNode*>(LineElement()()->Nodes()[0]);
  CoNode* ns2 = dynamic_cast<CoNode*>(LineElement()()->Nodes()[1]);
  lvec[0] = ns1->xspatial()[0] - ns2->xspatial()[0];
  lvec[1] = ns1->xspatial()[1] - ns2->xspatial()[1];
  lvec[2] = ns1->xspatial()[2] - ns2->xspatial()[2];

  // calculate lengths
  const double lengthS = sqrt(lvec[0] * lvec[0] + lvec[1] * lvec[1] + lvec[2] * lvec[2]);
  const double lengthA = sqrt(AuxnSurf()[0] * AuxnSurf()[0] + AuxnSurf()[1] * AuxnSurf()[1] +
                              AuxnSurf()[2] * AuxnSurf()[2]);
  const double prod = lengthS * lengthA;
  if (prod < 1e-12) return false;

  // calculate scalar product
  double scaprod = lvec[0] * AuxnSurf()[0] + lvec[1] * AuxnSurf()[1] + lvec[2] * AuxnSurf()[2];
  scaprod = scaprod / (prod);
  double diff = abs(scaprod) - 1.0;

  if (abs(diff) < tol) return false;

  return true;
}


/*----------------------------------------------------------------------*
 |  calculate dual shape functions                           farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineToSurfaceCoupling3d::ConsistDualShape()
{
  INPAR::MORTAR::ShapeFcn shapefcn =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(imortar_, "LM_SHAPEFCN");
  INPAR::MORTAR::ConsistentDualType consistent =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ConsistentDualType>(imortar_, "LM_DUAL_CONSISTENT");

  if (shapefcn != INPAR::MORTAR::shape_dual && shapefcn != INPAR::MORTAR::shape_petrovgalerkin)
    return;

  if (consistent == INPAR::MORTAR::consistent_none) return;

  if (IType() == LineToSurfaceCoupling3d::lts) return;
  dserror("ERROR: consistent dual shapes for stl is experimental!");

  // slave nodes and dofs
  const int max_nnodes = 9;
  const int nnodes = SurfaceElement().NumNode();
  if (nnodes > max_nnodes)
    dserror(
        "this function is not implemented to handle elements with that many nodes. Just adjust "
        "max_nnodes above");
  const int ndof = 3;

  // get number of master nodes
  int mnodes = LineElement()->NumNode();

  // Dual shape functions coefficient matrix and linearization
  LINALG::SerialDenseMatrix ae(nnodes, nnodes, true);
  SurfaceElement().MoData().DerivDualShape() =
      Teuchos::rcp(new GEN::pairedvector<int, Epetra_SerialDenseMatrix>(
          (nnodes + mnodes) * ndof, 0, Epetra_SerialDenseMatrix(nnodes, nnodes)));
  GEN::pairedvector<int, Epetra_SerialDenseMatrix>& derivae =
      *(SurfaceElement().MoData().DerivDualShape());

  // various variables
  double detg = 0.0;
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // initialize matrices de and me
  LINALG::SerialDenseMatrix me(nnodes, nnodes, true);
  LINALG::SerialDenseMatrix de(nnodes, nnodes, true);

  GEN::pairedvector<int, LINALG::Matrix<max_nnodes + 1, max_nnodes>> derivde_new(
      (nnodes + mnodes) * ndof);

  // two-dim arrays of maps for linearization of me/de
  std::vector<std::vector<GEN::pairedvector<int, double>>> derivme(
      nnodes, std::vector<GEN::pairedvector<int, double>>(nnodes, (nnodes + mnodes) * ndof));
  std::vector<std::vector<GEN::pairedvector<int, double>>> derivde(
      nnodes, std::vector<GEN::pairedvector<int, double>>(nnodes, (nnodes + mnodes) * ndof));

  double A_tot = 0.;

  // get number of master nodes
  const int ncol = LineElement()->NumNode();

  Teuchos::RCP<MORTAR::IntCell> currcell = IntLine();

  A_tot += currcell->Area();

  // create an integrator for this cell
  CONTACT::CoIntegrator integrator(imortar_, currcell->Shape(), Comm());

  // check if the cells are tri3
  // there's nothing wrong about other shapes, but as long as they are all
  // tri3 we can perform the jacobian calculation ( and its deriv) outside
  // the Gauss point loop
  if (currcell->Shape() != DRT::Element::line2)
    dserror("only line2 integration cells at the moment. See comment in the code");

  detg = currcell->Jacobian();
  // directional derivative of cell Jacobian
  GEN::pairedvector<int, double> derivjaccell((nnodes + ncol) * ndof);
  currcell->DerivJacobian(derivjaccell);

  for (int gp = 0; gp < integrator.nGP(); ++gp)
  {
    // coordinates and weight
    double eta[2] = {integrator.Coordinate(gp, 0), integrator.Coordinate(gp, 1)};
    const double wgt = integrator.Weight(gp);

    // get global Gauss point coordinates
    double globgp[3] = {0.0, 0.0, 0.0};
    currcell->LocalToGlobal(eta, globgp, 0);

    // project Gauss point onto slave integration element
    double sxi[2] = {0.0, 0.0};
    double sprojalpha = 0.0;
    MORTAR::MortarProjector::Impl(SurfaceElement())
        ->ProjectGaussPointAuxn3D(globgp, Auxn(), SurfaceElement(), sxi, sprojalpha);

    // project Gauss point onto slave (parent) element
    double psxi[2] = {0., 0.};
    for (int i = 0; i < 2; ++i) psxi[i] = sxi[i];

    // create vector for shape function evaluation
    LINALG::SerialDenseVector sval(nnodes);
    LINALG::SerialDenseMatrix sderiv(nnodes, 2, true);

    // evaluate trace space shape functions at Gauss point
    SurfaceElement().EvaluateShape(psxi, sval, sderiv, nnodes);

    // additional data for contact calculation (i.e. incl. derivative of dual shape functions
    // coefficient matrix) GP slave coordinate derivatives
    std::vector<GEN::pairedvector<int, double>> dsxigp(2, (nnodes + ncol) * ndof);
    // GP slave coordinate derivatives
    std::vector<GEN::pairedvector<int, double>> dpsxigp(2, (nnodes + ncol) * ndof);
    // global GP coordinate derivative on integration element
    GEN::pairedvector<int, LINALG::Matrix<3, 1>> lingp((nnodes + ncol) * ndof);

    // compute global GP coordinate derivative
    static LINALG::Matrix<3, 1> svalcell;
    static LINALG::Matrix<3, 2> sderivcell;
    currcell->EvaluateShape(eta, svalcell, sderivcell);

    for (int v = 0; v < 2; ++v)
      for (int d = 0; d < 3; ++d)
        for (_CI p = (currcell->GetDerivVertex(v))[d].begin();
             p != (currcell->GetDerivVertex(v))[d].end(); ++p)
          lingp[p->first](d) += svalcell(v) * (p->second);

    // compute GP slave coordinate derivatives
    integrator.DerivXiGP3DAuxPlane(SurfaceElement(), sxi, currcell->Auxn(), dsxigp, sprojalpha,
        currcell->GetDerivAuxn(), lingp);

    // compute GP slave coordinate derivatives (parent element)
    dpsxigp = dsxigp;

    double fac = 0.;
    for (_CI p = derivjaccell.begin(); p != derivjaccell.end(); ++p)
    {
      LINALG::Matrix<max_nnodes + 1, max_nnodes>& dtmp = derivde_new[p->first];
      const double& ps = p->second;
      for (int j = 0; j < nnodes; ++j)
      {
        fac = wgt * sval[j] * ps;
        dtmp(nnodes, j) += fac;
        for (int k = 0; k < nnodes; ++k) dtmp(k, j) += fac * sval[k];
      }
    }

    for (int i = 0; i < 2; ++i)
      for (_CI p = dpsxigp[i].begin(); p != dpsxigp[i].end(); ++p)
      {
        LINALG::Matrix<max_nnodes + 1, max_nnodes>& dtmp = derivde_new[p->first];
        const double& ps = p->second;
        for (int j = 0; j < nnodes; ++j)
        {
          fac = wgt * sderiv(j, i) * detg * ps;
          dtmp(nnodes, j) += fac;
          for (int k = 0; k < nnodes; ++k)
          {
            dtmp(k, j) += fac * sval[k];
            dtmp(j, k) += fac * sval[k];
          }
        }
      }

    // computing de, derivde and me, derivme and kappa, derivkappa
    for (int j = 0; j < nnodes; ++j)
    {
      double fac;
      fac = sval[j] * wgt;
      // computing de
      de(j, j) += fac * detg;

      for (int k = 0; k < nnodes; ++k)
      {
        // computing me
        fac = wgt * sval[j] * sval[k];
        me(j, k) += fac * detg;
      }
    }
  }

  // in case of no overlap just return, as there is no integration area
  // and therefore the consistent dual shape functions are not defined.
  // This doesn't matter, as there is no associated integration domain anyway
  if (A_tot < 1.e-12) return;

  // invert bi-ortho matrix me
  //  LINALG::SerialDenseMatrix meinv = LINALG::InvertAndMultiplyByCholesky(me,de,ae);

  LINALG::NonSymmetricInverse(me, 4);
  LINALG::SerialDenseMatrix meinv = me;

  // build linearization of ae and store in derivdual
  // (this is done according to a quite complex formula, which
  // we get from the linearization of the biorthogonality condition:
  // Lin (Me * Ae = De) -> Lin(Ae)=Lin(De)*Inv(Me)-Ae*Lin(Me)*Inv(Me) )
  typedef GEN::pairedvector<int, LINALG::Matrix<max_nnodes + 1, max_nnodes>>::const_iterator _CIM;
  for (_CIM p = derivde_new.begin(); p != derivde_new.end(); ++p)
  {
    LINALG::Matrix<max_nnodes + 1, max_nnodes>& dtmp = derivde_new[p->first];
    Epetra_SerialDenseMatrix& pt = derivae[p->first];
    for (int i = 0; i < nnodes; ++i)
      for (int j = 0; j < nnodes; ++j)
      {
        pt(i, j) += meinv(i, j) * dtmp(nnodes, i);

        for (int k = 0; k < nnodes; ++k)
          for (int l = 0; l < nnodes; ++l) pt(i, j) -= ae(i, k) * meinv(l, j) * dtmp(l, k);
      }
  }

  // store ae matrix in slave element data container
  SurfaceElement().MoData().DualShape() = Teuchos::rcp(new LINALG::SerialDenseMatrix(ae));

  return;
}

/*----------------------------------------------------------------------*
 |  integration for LTS (public)                             farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineToSurfaceCoupling3d::IntegrateLine()
{
  // get solution strategy
  INPAR::CONTACT::SolvingStrategy sol =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(imortar_, "STRATEGY");

  // create integrator object
  Teuchos::RCP<CONTACT::CoIntegrator> integrator =
      CONTACT::INTEGRATOR::BuildIntegrator(sol, imortar_, IntLine()->Shape(), Comm());

  // perform integration
  if (IType() == LineToSurfaceCoupling3d::lts)
  {
    integrator->IntegrateDerivCell3DAuxPlaneLTS(
        ParentElement(), *LineElement(), SurfaceElement(), IntLine(), Auxn(), Comm());
  }
  else if (IType() == LineToSurfaceCoupling3d::stl)
  {
    integrator->IntegrateDerivCell3DAuxPlaneSTL(
        ParentElement(), *LineElement(), SurfaceElement(), IntLine(), Auxn(), Comm());
  }
  else
    dserror("ERROR: wrong integration type for line coupling!");

  return;
}

/*----------------------------------------------------------------------*
 |  check if all vertices are along one line                 farah 09/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineToSurfaceCoupling3d::CheckLineOnLine(MORTAR::Vertex& edgeVertex1,
    MORTAR::Vertex& edgeVertex0, MORTAR::Vertex& lineVertex1, MORTAR::Vertex& lineVertex0)
{
  // tolerance for line clipping
  const double sminedge = ParentElement().MinEdgeSize();
  const double mminedge = SurfaceElement().MinEdgeSize();
  const double tol = MORTARCLIPTOL * std::min(sminedge, mminedge);

  // check if point of edge is on line
  bool lineOnLine = false;
  double line[3] = {0.0, 0.0, 0.0};
  double edgeLine[3] = {0.0, 0.0, 0.0};

  for (int k = 0; k < 3; ++k)
  {
    line[k] = lineVertex1.Coord()[k] - lineVertex0.Coord()[k];
    edgeLine[k] = edgeVertex1.Coord()[k] - lineVertex0.Coord()[k];
  }

  double lengthLine = sqrt(line[0] * line[0] + line[1] * line[1] + line[2] * line[2]);
  double lengthedge =
      sqrt(edgeLine[0] * edgeLine[0] + edgeLine[1] * edgeLine[1] + edgeLine[2] * edgeLine[2]);

  if (lengthLine < tol) dserror("ERROR: Line Element is of zero length!");

  if (lengthedge < tol)
  {
    lineOnLine = true;
  }
  else
  {
    // calc scalar product
    double scaprod = line[0] * edgeLine[0] + line[1] * edgeLine[1] + line[2] * edgeLine[2];
    scaprod /= (lengthLine * lengthedge);

    if ((abs(scaprod) - tol < 1.0) and (abs(scaprod) + tol > 1.0)) lineOnLine = true;
  }

  if (!lineOnLine) return false;

  return true;
}

/*----------------------------------------------------------------------*
 |  geometric stuff (private)                                farah 08/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineToSurfaceCoupling3d::LineToLineClipping(MORTAR::Vertex& edgeVertex1,
    MORTAR::Vertex& edgeVertex0, MORTAR::Vertex& lineVertex1, MORTAR::Vertex& lineVertex0)
{
  // output bool
  const bool out = false;

  // tolerance for line clipping
  const double sminedge = ParentElement().MinEdgeSize();
  const double mminedge = SurfaceElement().MinEdgeSize();
  const double tol = MORTARCLIPTOL * std::min(sminedge, mminedge);

  bool lineOnLine = CheckLineOnLine(edgeVertex1, edgeVertex0, lineVertex1, lineVertex0);

  if (!lineOnLine) dserror("ERROR: vertices not along a line, but already checked!");

  double line[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; ++k) line[k] = lineVertex1.Coord()[k] - lineVertex0.Coord()[k];

  // LINE ON LINE!!! go on with real line to line clipping
  bool e0v0 = false;
  bool e0v1 = false;
  bool e1v0 = false;
  bool e1v1 = false;
  double prod0 = 0.0;
  double prod1 = 0.0;
  double prod2 = 0.0;
  double prod3 = 0.0;
  // check of both master vertices are out of line in 0 direction
  double lineEdge0Vert0[3] = {0.0, 0.0, 0.0};
  double lineEdge1Vert0[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; ++k)
  {
    lineEdge0Vert0[k] = edgeVertex0.Coord()[k] - lineVertex0.Coord()[k];
    lineEdge1Vert0[k] = edgeVertex1.Coord()[k] - lineVertex0.Coord()[k];
  }

  for (int k = 0; k < 3; ++k)
  {
    prod0 += lineEdge0Vert0[k] * line[k];
    prod1 += lineEdge1Vert0[k] * line[k];
  }

  if (prod0 < 0.0)
    e0v0 = false;
  else
    e0v0 = true;
  if (prod1 < 0.0)
    e1v0 = false;
  else
    e1v0 = true;


  // check of both master vertices are out of line in 1 direction
  double lineEdge0Vert1[3] = {0.0, 0.0, 0.0};
  double lineEdge1Vert1[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; ++k)
  {
    lineEdge0Vert1[k] = edgeVertex0.Coord()[k] - lineVertex1.Coord()[k];
    lineEdge1Vert1[k] = edgeVertex1.Coord()[k] - lineVertex1.Coord()[k];
  }

  for (int k = 0; k < 3; ++k)
  {
    prod2 -= lineEdge0Vert1[k] * line[k];
    prod3 -= lineEdge1Vert1[k] * line[k];
  }

  if (prod2 < 0.0)
    e0v1 = false;
  else
    e0v1 = true;
  if (prod3 < 0.0)
    e1v1 = false;
  else
    e1v1 = true;

  // check if vertices are lying on each other
  bool e0isV0 = true;
  bool e0isV1 = true;
  bool e1isV0 = true;
  bool e1isV1 = true;

  double test0[3] = {0.0, 0.0, 0.0};
  double test1[3] = {0.0, 0.0, 0.0};
  double test2[3] = {0.0, 0.0, 0.0};
  double test3[3] = {0.0, 0.0, 0.0};

  for (int k = 0; k < 3; ++k)
  {
    test0[k] = edgeVertex0.Coord()[k] - lineVertex0.Coord()[k];
    test1[k] = edgeVertex0.Coord()[k] - lineVertex1.Coord()[k];
    test2[k] = edgeVertex1.Coord()[k] - lineVertex0.Coord()[k];
    test3[k] = edgeVertex1.Coord()[k] - lineVertex1.Coord()[k];
  }

  double l0 = sqrt(test0[0] * test0[0] + test0[1] * test0[1] + test0[2] * test0[2]);
  double l1 = sqrt(test1[0] * test1[0] + test1[1] * test1[1] + test1[2] * test1[2]);
  double l2 = sqrt(test2[0] * test2[0] + test2[1] * test2[1] + test2[2] * test2[2]);
  double l3 = sqrt(test3[0] * test3[0] + test3[1] * test3[1] + test3[2] * test3[2]);

  if (abs(l0) > tol) e0isV0 = false;
  if (abs(l1) > tol) e0isV1 = false;
  if (abs(l2) > tol) e1isV0 = false;
  if (abs(l3) > tol) e1isV1 = false;

  // ========================================================
  // 1.: nodes on each other
  if (e0isV0 and e1isV1)
  {
    if (out) std::cout << "CASE 1" << std::endl;
    InterSections().push_back(MORTAR::Vertex(edgeVertex0.Coord(), MORTAR::Vertex::master,
        edgeVertex0.Nodeids(), NULL, NULL, false, false, NULL, -1));

    InterSections().push_back(MORTAR::Vertex(edgeVertex1.Coord(), MORTAR::Vertex::master,
        edgeVertex1.Nodeids(), NULL, NULL, false, false, NULL, -1));
  }
  // ========================================================
  // 2.: nodes on each other
  else if (e0isV1 and e1isV0)
  {
    if (out) std::cout << "CASE 2" << std::endl;

    InterSections().push_back(MORTAR::Vertex(edgeVertex0.Coord(), MORTAR::Vertex::master,
        edgeVertex0.Nodeids(), NULL, NULL, false, false, NULL, -1));

    InterSections().push_back(MORTAR::Vertex(edgeVertex1.Coord(), MORTAR::Vertex::master,
        edgeVertex1.Nodeids(), NULL, NULL, false, false, NULL, -1));
  }
  // ========================================================
  // 3.: e0 on v0 and e1 valid
  else if (e0isV0 and e1v0 and e1v1)
  {
    if (out) std::cout << "CASE 3" << std::endl;

    InterSections().push_back(MORTAR::Vertex(edgeVertex0.Coord(), MORTAR::Vertex::master,
        edgeVertex0.Nodeids(), NULL, NULL, false, false, NULL, -1));

    InterSections().push_back(MORTAR::Vertex(edgeVertex1.Coord(), MORTAR::Vertex::master,
        edgeVertex1.Nodeids(), NULL, NULL, false, false, NULL, -1));
  }
  // ========================================================
  // 4.: e1 on v0 and e0 valid
  else if (e1isV0 and e0v0 and e0v1)
  {
    if (out) std::cout << "CASE 4" << std::endl;

    InterSections().push_back(MORTAR::Vertex(edgeVertex0.Coord(), MORTAR::Vertex::master,
        edgeVertex0.Nodeids(), NULL, NULL, false, false, NULL, -1));

    InterSections().push_back(MORTAR::Vertex(edgeVertex1.Coord(), MORTAR::Vertex::master,
        edgeVertex1.Nodeids(), NULL, NULL, false, false, NULL, -1));
  }
  // ========================================================
  // 5.: e0 on v1 and e1 valid
  else if (e0isV1 and e1v0 and e1v1)
  {
    if (out) std::cout << "CASE 5" << std::endl;

    InterSections().push_back(MORTAR::Vertex(edgeVertex0.Coord(), MORTAR::Vertex::master,
        edgeVertex0.Nodeids(), NULL, NULL, false, false, NULL, -1));

    InterSections().push_back(MORTAR::Vertex(edgeVertex1.Coord(), MORTAR::Vertex::master,
        edgeVertex1.Nodeids(), NULL, NULL, false, false, NULL, -1));
  }
  // ========================================================
  // 6.: e1 on v1 and e0 valid
  else if (e1isV1 and e0v0 and e0v1)
  {
    if (out) std::cout << "CASE 6" << std::endl;

    InterSections().push_back(MORTAR::Vertex(edgeVertex0.Coord(), MORTAR::Vertex::master,
        edgeVertex0.Nodeids(), NULL, NULL, false, false, NULL, -1));

    InterSections().push_back(MORTAR::Vertex(edgeVertex1.Coord(), MORTAR::Vertex::master,
        edgeVertex1.Nodeids(), NULL, NULL, false, false, NULL, -1));
  }
  // ========================================================
  // 7.: e0 on v0 and e1 out of v1 but in v0
  else if (e0isV0 and e1v0 and !e1v1)
  {
    if (out) std::cout << "CASE 7" << std::endl;

    InterSections().push_back(MORTAR::Vertex(edgeVertex0.Coord(), MORTAR::Vertex::master,
        edgeVertex0.Nodeids(), NULL, NULL, false, false, NULL, -1));

    InterSections().push_back(MORTAR::Vertex(lineVertex1.Coord(), MORTAR::Vertex::projslave,
        lineVertex1.Nodeids(), NULL, NULL, false, false, NULL, -1));
  }
  // ========================================================
  // 8.: e1 on v0 and e0 out of v1 but in v0
  else if (e1isV0 and e0v0 and !e0v1)
  {
    if (out) std::cout << "CASE 8" << std::endl;

    InterSections().push_back(MORTAR::Vertex(edgeVertex1.Coord(), MORTAR::Vertex::master,
        edgeVertex1.Nodeids(), NULL, NULL, false, false, NULL, -1));

    InterSections().push_back(MORTAR::Vertex(lineVertex1.Coord(), MORTAR::Vertex::projslave,
        lineVertex1.Nodeids(), NULL, NULL, false, false, NULL, -1));
  }
  // ========================================================
  // 9.: e1 on v1 and e0 out of v0 but in v1
  else if (e1isV1 and !e0v0 and e0v1)
  {
    if (out) std::cout << "CASE 9" << std::endl;

    InterSections().push_back(MORTAR::Vertex(edgeVertex1.Coord(), MORTAR::Vertex::master,
        edgeVertex1.Nodeids(), NULL, NULL, false, false, NULL, -1));

    InterSections().push_back(MORTAR::Vertex(lineVertex0.Coord(), MORTAR::Vertex::projslave,
        lineVertex0.Nodeids(), NULL, NULL, false, false, NULL, -1));
  }
  // ========================================================
  // 10.: e0 on v1 and e1 out of v0 but in v1
  else if (e0isV1 and !e1v0 and e1v1)
  {
    if (out) std::cout << "CASE 10" << std::endl;

    InterSections().push_back(MORTAR::Vertex(edgeVertex0.Coord(), MORTAR::Vertex::master,
        edgeVertex0.Nodeids(), NULL, NULL, false, false, NULL, -1));

    InterSections().push_back(MORTAR::Vertex(lineVertex0.Coord(), MORTAR::Vertex::projslave,
        lineVertex0.Nodeids(), NULL, NULL, false, false, NULL, -1));
  }
  // ========================================================
  // 11.: e0 on v0 and e1 out of v0 but in v1
  else if (e0isV0 and !e1v0 and e1v1)
  {
    if (out) std::cout << "CASE 11" << std::endl;

    // true because no more intersections expected
    return true;
  }
  // ========================================================
  // 12.: e1 on v0 and e0 out of v0 but in v1
  else if (e1isV0 and !e0v0 and e0v1)
  {
    if (out) std::cout << "CASE 12" << std::endl;

    // true because no more intersections expected
    return true;
  }
  // ========================================================
  // 13.: e0 on v1 and e1 out of v1 but in v0
  else if (e0isV1 and !e1v1 and e1v0)
  {
    if (out) std::cout << "CASE 13" << std::endl;

    // true because no more intersections expected
    return true;
  }
  // ========================================================
  // 14.: e1 on v1 and e0 out of v1 but in v0
  else if (e1isV1 and !e0v1 and e0v0)
  {
    if (out) std::cout << "CASE 14" << std::endl;

    // true because no more intersections expected
    return true;
  }
  // ========================================================
  // 15.: all true --> both intersections master nodes
  else if (e0v1 and e1v1 and e0v0 and e1v0)
  {
    if (out) std::cout << "CASE 15" << std::endl;

    InterSections().push_back(MORTAR::Vertex(edgeVertex0.Coord(), MORTAR::Vertex::master,
        edgeVertex0.Nodeids(), NULL, NULL, false, false, NULL, -1));

    InterSections().push_back(MORTAR::Vertex(edgeVertex1.Coord(), MORTAR::Vertex::master,
        edgeVertex1.Nodeids(), NULL, NULL, false, false, NULL, -1));
  }
  // ========================================================
  // 16.: all slave nodes are projected: E0 out of V0  and E1 out of V1
  else if (!e0v0 and e1v0 and e0v1 and !e1v1)
  {
    if (out) std::cout << "CASE 16" << std::endl;

    InterSections().push_back(MORTAR::Vertex(lineVertex0.Coord(), MORTAR::Vertex::projslave,
        lineVertex0.Nodeids(), NULL, NULL, false, false, NULL, -1));

    InterSections().push_back(MORTAR::Vertex(lineVertex1.Coord(), MORTAR::Vertex::projslave,
        lineVertex1.Nodeids(), NULL, NULL, false, false, NULL, -1));
  }
  // ========================================================
  // 17.: all slave nodes are projected: E1 out of V0  and E0 out of V1
  else if (e0v0 and !e1v0 and !e0v1 and e1v1)
  {
    if (out) std::cout << "CASE 17" << std::endl;

    InterSections().push_back(MORTAR::Vertex(lineVertex0.Coord(), MORTAR::Vertex::projslave,
        lineVertex0.Nodeids(), NULL, NULL, false, false, NULL, -1));

    InterSections().push_back(MORTAR::Vertex(lineVertex1.Coord(), MORTAR::Vertex::projslave,
        lineVertex1.Nodeids(), NULL, NULL, false, false, NULL, -1));
  }
  // ========================================================
  // 18.: mixed: E0 and E1 pos to V0 and E0 pos to V1
  else if (e0v0 and e1v0 and e0v1 and !e1v1)
  {
    if (out) std::cout << "CASE 18" << std::endl;

    InterSections().push_back(MORTAR::Vertex(edgeVertex0.Coord(), MORTAR::Vertex::master,
        edgeVertex0.Nodeids(), NULL, NULL, false, false, NULL, -1));

    InterSections().push_back(MORTAR::Vertex(lineVertex1.Coord(), MORTAR::Vertex::projslave,
        lineVertex1.Nodeids(), NULL, NULL, false, false, NULL, -1));
  }
  // ========================================================
  // 19.: mixed: E0 and E1 pos to V0 and E1 pos to V1
  else if (e0v0 and e1v0 and !e0v1 and e1v1)
  {
    if (out) std::cout << "CASE 19" << std::endl;

    InterSections().push_back(MORTAR::Vertex(edgeVertex1.Coord(), MORTAR::Vertex::master,
        edgeVertex1.Nodeids(), NULL, NULL, false, false, NULL, -1));

    InterSections().push_back(MORTAR::Vertex(lineVertex1.Coord(), MORTAR::Vertex::projslave,
        lineVertex1.Nodeids(), NULL, NULL, false, false, NULL, -1));
  }
  // ========================================================
  // 20.: mixed: E0 neg and E1 pos to V0 and E0 and E1 pos to V1
  else if (!e0v0 and e1v0 and e0v1 and e1v1)
  {
    if (out) std::cout << "CASE 20" << std::endl;

    InterSections().push_back(MORTAR::Vertex(edgeVertex1.Coord(), MORTAR::Vertex::master,
        edgeVertex1.Nodeids(), NULL, NULL, false, false, NULL, -1));

    InterSections().push_back(MORTAR::Vertex(lineVertex0.Coord(), MORTAR::Vertex::projslave,
        lineVertex0.Nodeids(), NULL, NULL, false, false, NULL, -1));
  }
  // ========================================================
  // 21.: mixed: E1 neg and E0 pos to V0 and E0 and E1 pos to V1
  else if (e0v0 and !e1v0 and e0v1 and e1v1)
  {
    if (out) std::cout << "CASE 21" << std::endl;

    InterSections().push_back(MORTAR::Vertex(edgeVertex0.Coord(), MORTAR::Vertex::master,
        edgeVertex0.Nodeids(), NULL, NULL, false, false, NULL, -1));

    InterSections().push_back(MORTAR::Vertex(lineVertex0.Coord(), MORTAR::Vertex::projslave,
        lineVertex0.Nodeids(), NULL, NULL, false, false, NULL, -1));
  }
  // ========================================================
  // 22.: out: E0 and E1 in V0 and out of V1
  else if (e0v0 and e1v0 and !e0v1 and !e1v1)
  {
    if (out) std::cout << "CASE 22" << std::endl;

    // true because no more intersections expected
    return true;
  }
  // ========================================================
  // 23.: out: E0 and E1 in V1 and out of V0
  else if (!e0v0 and !e1v0 and e0v1 and e1v1)
  {
    if (out) std::cout << "CASE 23" << std::endl;

    // true because no more intersections expected
    return true;
  }
  // ========================================================
  // no valid intersection
  else
  {
    std::cout << "e0isV0 = " << e0isV0 << std::endl;
    std::cout << "e0isV1 = " << e0isV1 << std::endl;
    std::cout << "e1isV0 = " << e1isV0 << std::endl;
    std::cout << "e1isV1 = " << e1isV1 << std::endl;

    std::cout << "e0v0 = " << e0v0 << std::endl;
    std::cout << "e1v0 = " << e1v0 << std::endl;
    std::cout << "e0v1 = " << e0v1 << std::endl;
    std::cout << "e1v1 = " << e1v1 << std::endl;

    dserror("ERROR: Something went terribly wrong!");
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  geometric stuff (private)                                farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineToSurfaceCoupling3d::LineClipping()
{
  // output variable
  const bool out = false;

  // tolerance for line clipping
  const double sminedge = ParentElement().MinEdgeSize();
  const double mminedge = SurfaceElement().MinEdgeSize();
  const double tol = MORTARCLIPTOL * std::min(sminedge, mminedge);

  // vector with vertices
  InterSections().clear();
  TempInterSections().clear();

  // safety
  if (MasterVertices().size() < 3) dserror("ERROR: Invalid number of Master Vertices!");
  if (SlaveVertices().size() != 2) dserror("ERROR: Invalid number of Slave Vertices!");

  // set previous and next Vertex pointer for all elements in lists
  for (int i = 0; i < (int)MasterVertices().size(); ++i)
  {
    // standard case
    if (i != 0 && i != (int)MasterVertices().size() - 1)
    {
      MasterVertices()[i].AssignNext(&MasterVertices()[i + 1]);
      MasterVertices()[i].AssignPrev(&MasterVertices()[i - 1]);
    }
    // first element in list
    else if (i == 0)
    {
      MasterVertices()[i].AssignNext(&MasterVertices()[i + 1]);
      MasterVertices()[i].AssignPrev(&MasterVertices()[(int)MasterVertices().size() - 1]);
    }
    // last element in list
    else
    {
      MasterVertices()[i].AssignNext(&MasterVertices()[0]);
      MasterVertices()[i].AssignPrev(&MasterVertices()[i - 1]);
    }
  }

  // flip ordering
  std::reverse(SlaveVertices().begin(), SlaveVertices().end());

  // create line from slave vertices
  double slaveLine[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; ++k)
    slaveLine[k] = SlaveVertices()[1].Coord()[k] - SlaveVertices()[0].Coord()[k];


  // check for parallelity of line and edges and perform line to line clipping
  bool foundValidParallelity = false;

  // loop over master vertices to create master polygon lines
  for (int j = 0; j < (int)MasterVertices().size(); ++j)
  {
    // we need one edge first
    double edge[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < 3; ++k)
      edge[k] = (MasterVertices()[j].Next())->Coord()[k] - MasterVertices()[j].Coord()[k];

    // outward edge normals of polygon and slave line
    double np[3] = {0.0, 0.0, 0.0};
    double nl[3] = {0.0, 0.0, 0.0};
    np[0] = edge[1] * AuxnSurf()[2] - edge[2] * AuxnSurf()[1];
    np[1] = edge[2] * AuxnSurf()[0] - edge[0] * AuxnSurf()[2];
    np[2] = edge[0] * AuxnSurf()[1] - edge[1] * AuxnSurf()[0];
    nl[0] = slaveLine[1] * AuxnSurf()[2] - slaveLine[2] * AuxnSurf()[1];
    nl[1] = slaveLine[2] * AuxnSurf()[0] - slaveLine[0] * AuxnSurf()[2];
    nl[2] = slaveLine[0] * AuxnSurf()[1] - slaveLine[1] * AuxnSurf()[0];

    if (out)
    {
      std::cout << "==============================================" << std::endl;
      std::cout << "SLine= " << slaveLine[0] << "  " << slaveLine[1] << "  " << slaveLine[2]
                << std::endl;
      std::cout << "Pos1= " << SlaveVertices()[0].Coord()[0] << "  "
                << SlaveVertices()[0].Coord()[1] << "  " << SlaveVertices()[0].Coord()[2]
                << std::endl;
      std::cout << "Pos2= " << SlaveVertices()[1].Coord()[0] << "  "
                << SlaveVertices()[1].Coord()[1] << "  " << SlaveVertices()[1].Coord()[2]
                << std::endl;
      std::cout << "N slave= " << nl[0] << "  " << nl[1] << "  " << nl[2] << std::endl;


      std::cout << "==============================================" << std::endl;
      std::cout << "MEdge= " << edge[0] << "  " << edge[1] << "  " << edge[2] << std::endl;
      std::cout << "Pos1= " << (MasterVertices()[j].Next())->Coord()[0] << "  "
                << (MasterVertices()[j].Next())->Coord()[1] << "  "
                << (MasterVertices()[j].Next())->Coord()[2] << std::endl;
      std::cout << "Pos2= " << MasterVertices()[j].Coord()[0] << "  "
                << MasterVertices()[j].Coord()[1] << "  " << MasterVertices()[j].Coord()[2]
                << std::endl;
      std::cout << "N master= " << np[0] << "  " << np[1] << "  " << np[2] << std::endl;
    }

    // check for parallelity of edges
    double parallel = edge[0] * nl[0] + edge[1] * nl[1] + edge[2] * nl[2];
    if (abs(parallel) < tol)
    {
      // safety checks
      if (MasterVertices()[j].Next()->Nodeids().size() > 1)
        dserror("ERROR: Only one node id per master vertex allowed!");
      if (MasterVertices()[j].Nodeids().size() > 1)
        dserror("ERROR: Only one node id per master vertex allowed!");

      // store master node ids in set to guarantee uniquness
      std::pair<int, int> actIDs = std::pair<int, int>(
          MasterVertices()[j].Next()->Nodeids()[0], MasterVertices()[j].Nodeids()[0]);
      std::pair<int, int> actIDsTw = std::pair<int, int>(
          MasterVertices()[j].Nodeids()[0], MasterVertices()[j].Next()->Nodeids()[0]);

      // check if edge on line element
      foundValidParallelity = CheckLineOnLine(
          *MasterVertices()[j].Next(), MasterVertices()[j], SlaveVertices()[1], SlaveVertices()[0]);

      // check if processed before
      std::set<std::pair<int, int>>::iterator iter = DoneBefore().find(actIDs);
      std::set<std::pair<int, int>>::iterator itertw = DoneBefore().find(actIDsTw);

      // if not perform clipping of lines
      if (iter == DoneBefore().end() and itertw == DoneBefore().end())
      {
        // add to set of processed nodes
        DoneBefore().insert(actIDs);
        DoneBefore().insert(actIDsTw);

        if (foundValidParallelity)
        {
          // perform line-line clipping
          LineToLineClipping(*MasterVertices()[j].Next(), MasterVertices()[j], SlaveVertices()[1],
              SlaveVertices()[0]);

          if (out)
            std::cout << "MASTER IDS = " << MasterVertices()[j].Next()->Nodeids()[0] << "  "
                      << MasterVertices()[j].Nodeids()[0] << std::endl;
          break;
        }
        else
          continue;
      }
    }
  }  // end master vertex loop

  // if there is a line to line setting --> jump to node check
  if (!foundValidParallelity)
  {
    // loop over master vertices to create master polygon lines
    for (int j = 0; j < (int)MasterVertices().size(); ++j)
    {
      // we need two edges first
      double edge[3] = {0.0, 0.0, 0.0};
      for (int k = 0; k < 3; ++k)
        edge[k] = (MasterVertices()[j].Next())->Coord()[k] - MasterVertices()[j].Coord()[k];

      // outward edge normals of polygon and slave line
      double np[3] = {0.0, 0.0, 0.0};
      double nl[3] = {0.0, 0.0, 0.0};
      np[0] = edge[1] * AuxnSurf()[2] - edge[2] * AuxnSurf()[1];
      np[1] = edge[2] * AuxnSurf()[0] - edge[0] * AuxnSurf()[2];
      np[2] = edge[0] * AuxnSurf()[1] - edge[1] * AuxnSurf()[0];
      nl[0] = slaveLine[1] * AuxnSurf()[2] - slaveLine[2] * AuxnSurf()[1];
      nl[1] = slaveLine[2] * AuxnSurf()[0] - slaveLine[0] * AuxnSurf()[2];
      nl[2] = slaveLine[0] * AuxnSurf()[1] - slaveLine[1] * AuxnSurf()[0];

      if (out)
      {
        std::cout << "==============================================" << std::endl;
        std::cout << "SLine= " << slaveLine[0] << "  " << slaveLine[1] << "  " << slaveLine[2]
                  << std::endl;
        std::cout << "Pos1= " << SlaveVertices()[0].Coord()[0] << "  "
                  << SlaveVertices()[0].Coord()[1] << "  " << SlaveVertices()[0].Coord()[2]
                  << std::endl;
        std::cout << "Pos2= " << SlaveVertices()[1].Coord()[0] << "  "
                  << SlaveVertices()[1].Coord()[1] << "  " << SlaveVertices()[1].Coord()[2]
                  << std::endl;
        std::cout << "N slave= " << nl[0] << "  " << nl[1] << "  " << nl[2] << std::endl;


        std::cout << "==============================================" << std::endl;
        std::cout << "MEdge= " << edge[0] << "  " << edge[1] << "  " << edge[2] << std::endl;
        std::cout << "Pos1= " << (MasterVertices()[j].Next())->Coord()[0] << "  "
                  << (MasterVertices()[j].Next())->Coord()[1] << "  "
                  << (MasterVertices()[j].Next())->Coord()[2] << std::endl;
        std::cout << "Pos2= " << MasterVertices()[j].Coord()[0] << "  "
                  << MasterVertices()[j].Coord()[1] << "  " << MasterVertices()[j].Coord()[2]
                  << std::endl;
        std::cout << "N master= " << np[0] << "  " << np[1] << "  " << np[2] << std::endl;
      }

      // check for parallelity of edges
      double parallel = edge[0] * nl[0] + edge[1] * nl[1] + edge[2] * nl[2];
      if (abs(parallel) < tol)
      {
        continue;
      }

      // check for intersection of non-parallel edges
      double wec_p1 = 0.0;
      double wec_p2 = 0.0;
      for (int k = 0; k < 3; ++k)
      {
        wec_p1 += (SlaveVertices()[0].Coord()[k] - MasterVertices()[j].Coord()[k]) * np[k];
        wec_p2 += (SlaveVertices()[1].Coord()[k] - MasterVertices()[j].Coord()[k]) * np[k];
      }

      if (out)
      {
        std::cout << "WecP1 = " << wec_p1 << std::endl;
        std::cout << "WecP2 = " << wec_p2 << std::endl;
      }


      // change of sign means we have an intersection!
      if (wec_p1 * wec_p2 <= 0.0)
      {
        double wec_q1 = 0.0;
        double wec_q2 = 0.0;
        for (int k = 0; k < 3; ++k)
        {
          wec_q1 += (MasterVertices()[j].Coord()[k] - SlaveVertices()[0].Coord()[k]) * nl[k];
          wec_q2 +=
              ((MasterVertices()[j].Next())->Coord()[k] - SlaveVertices()[0].Coord()[k]) * nl[k];
        }

        if (out)
        {
          std::cout << "WecQ1 = " << wec_q1 << std::endl;
          std::cout << "WecQ2 = " << wec_q2 << std::endl;
        }

        if (wec_q1 * wec_q2 <= 0.0)
        {
          double alpha = wec_p1 / (wec_p1 - wec_p2);
          double alphaq = wec_q1 / (wec_q1 - wec_q2);

          if (alpha < 0.0 or alpha > 1.0) continue;
          if (alphaq < 0.0 or alphaq > 1.0) continue;

          std::vector<double> coords(3);
          for (int k = 0; k < 3; ++k)
          {
            coords[k] =
                (1 - alpha) * SlaveVertices()[0].Coord()[k] + alpha * SlaveVertices()[1].Coord()[k];
            if (abs(coords[k]) < tol) coords[k] = 0.0;
          }

          if (out)
          {
            std::cout << "Found intersection! (" << j << ") " << alpha << std::endl;
            std::cout << "coords 1: " << coords[0] << " " << coords[1] << " " << coords[2]
                      << std::endl;
          }

          // generate vectors of underlying node ids for lineclip (2x slave, 2x master)
          std::vector<int> lcids(4);
          lcids[0] = (int)(SlaveVertices()[0].Nodeids()[0]);
          lcids[1] = (int)(SlaveVertices()[1].Nodeids()[0]);
          lcids[2] = (int)(MasterVertices()[j].Nodeids()[0]);
          lcids[3] = (int)((MasterVertices()[j].Next())->Nodeids()[0]);

          // store intersection points
          TempInterSections().push_back(MORTAR::Vertex(coords, MORTAR::Vertex::lineclip, lcids,
              &SlaveVertices()[1], &SlaveVertices()[0], true, false, NULL, alpha));
        }
      }
    }  // end vertex loop

    // ===================================================
    // find interior node intersections
    //    if((int)TempInterSections().size()!=2)
    {
      for (int i = 0; i < (int)SlaveVertices().size(); ++i)
      {
        // keep track of inside / outside status
        bool outside = false;

        // check against all poly1 (slave) edges
        for (int j = 0; j < (int)MasterVertices().size(); ++j)
        {
          // we need diff vector and edge2 first
          double diff[3] = {0.0, 0.0, 0.0};
          double edge[3] = {0.0, 0.0, 0.0};
          for (int k = 0; k < 3; ++k)
          {
            diff[k] = SlaveVertices()[i].Coord()[k] - MasterVertices()[j].Coord()[k];
            edge[k] = (MasterVertices()[j].Next())->Coord()[k] - MasterVertices()[j].Coord()[k];
          }

          // compute distance from point on poly1 to edge
          double n[3] = {0.0, 0.0, 0.0};
          n[0] = edge[1] * AuxnSurf()[2] - edge[2] * AuxnSurf()[1];
          n[1] = edge[2] * AuxnSurf()[0] - edge[0] * AuxnSurf()[2];
          n[2] = edge[0] * AuxnSurf()[1] - edge[1] * AuxnSurf()[0];
          double ln = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
          for (int k = 0; k < 3; ++k) n[k] /= ln;

          double dist = diff[0] * n[0] + diff[1] * n[1] + diff[2] * n[2];

          // only keep point if not in outside halfspace
          if (dist - tol > 0.0)  // tends to include nodes
          {
            outside = true;
            break;
          }
        }  // end master loop

        if (outside)
        {
          // next slave vertex
          continue;
        }
        else
        {
          TempInterSections().push_back(
              MORTAR::Vertex(SlaveVertices()[i].Coord(), MORTAR::Vertex::projslave,
                  SlaveVertices()[i].Nodeids(), NULL, NULL, false, false, NULL, -1));
        }
      }
    }  // if intersections != 2

    // check positions of all found intersections
    std::vector<int> redundantLocalIDs;
    for (int i = 0; i < (int)TempInterSections().size(); ++i)
    {
      for (int j = i; j < (int)TempInterSections().size(); ++j)
      {
        // do not check same intersections
        if (i == j) continue;

        // distance vector
        double diff[3] = {0.0, 0.0, 0.0};
        for (int k = 0; k < 3; ++k)
          diff[k] = TempInterSections()[i].Coord()[k] - TempInterSections()[j].Coord()[k];
        double dist = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
        if (dist < tol)
        {
          // store redundant id
          redundantLocalIDs.push_back(j);
        }
      }
    }

    std::vector<MORTAR::Vertex> aux;
    for (int i = 0; i < (int)TempInterSections().size(); ++i)
    {
      bool vanish = false;
      for (int j = 0; j < (int)redundantLocalIDs.size(); ++j)
      {
        if (i == redundantLocalIDs[j]) vanish = true;
      }

      if (!vanish) aux.push_back(TempInterSections()[i]);
    }

    // store right vector to TempIntersections
    TempInterSections().clear();
    for (int i = 0; i < (int)aux.size(); ++i) TempInterSections().push_back(aux[i]);

    // ===================================================
    // check if intersection is close to a node
    for (int i = 0; i < (int)TempInterSections().size(); ++i)
    {
      // keep track of comparisons
      bool close = false;

      // check against all poly1 (slave) points
      for (int j = 0; j < (int)SlaveVertices().size(); ++j)
      {
        // distance vector
        double diff[3] = {0.0, 0.0, 0.0};
        for (int k = 0; k < 3; ++k)
          diff[k] = TempInterSections()[i].Coord()[k] - SlaveVertices()[j].Coord()[k];
        double dist = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);

        // only keep intersection point if not close
        if (dist <= tol)
        {
          // intersection is close to slave vertex!
          close = true;

          // store slave vertex as intersection point
          InterSections().push_back(
              MORTAR::Vertex(SlaveVertices()[j].Coord(), MORTAR::Vertex::projslave,
                  SlaveVertices()[j].Nodeids(), NULL, NULL, false, false, NULL, -1));
          break;
        }
      }

      // do only if no close slave point found
      if (!close)
      {
        // check against all poly2 (master) points
        for (int j = 0; j < (int)MasterVertices().size(); ++j)
        {
          // distance vector
          double diff[3] = {0.0, 0.0, 0.0};
          for (int k = 0; k < 3; ++k)
            diff[k] = TempInterSections()[i].Coord()[k] - MasterVertices()[j].Coord()[k];
          double dist = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);

          // only keep intersection point if not close
          if (dist <= tol)
          {
            // intersection is close to master vertex!
            close = true;

            InterSections().push_back(
                MORTAR::Vertex(MasterVertices()[j].Coord(), MORTAR::Vertex::master,
                    MasterVertices()[j].Nodeids(), NULL, NULL, false, false, NULL, -1));
            break;
          }
        }
      }

      // keep intersection point only if not close to any Slave/Master point
      if (!close) InterSections().push_back(TempInterSections()[i]);
    }
  }  // end found valid parallity

  // 2. check plausibility
  if (InterSections().size() > 2)
  {
    std::cout << "Intersections= " << InterSections().size() << std::endl;
    dserror("ERROR: intersections not possible!!!");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  create integration lines                                 farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineToSurfaceCoupling3d::CreateIntegrationLines(
    std::vector<std::vector<GEN::pairedvector<int, double>>>& linvertex)
{
  // get coordinates
  LINALG::Matrix<3, 3> coords;

  for (int i = 0; i < 3; ++i)
  {
    coords(i, 0) = InterSections()[0].Coord()[i];
    coords(i, 1) = InterSections()[1].Coord()[i];
    coords(i, 2) = -1;  // dummy;
  }

  // create Integration Line
  IntLine() = Teuchos::rcp(new MORTAR::IntCell(ParentElement().Id(), 2, coords, Auxn(),
      DRT::Element::line2, linvertex[0], linvertex[1],
      linvertex[1],  // dummy
      GetDerivAuxn()));

  return;
}

/*----------------------------------------------------------------------*
 |  linearize vertices                                       farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineToSurfaceCoupling3d::LinearizeVertices(
    std::vector<std::vector<GEN::pairedvector<int, double>>>& linvertex)
{
  // linearize all aux.plane slave and master nodes only ONCE
  // and use these linearizations later during lineclip linearization
  // (this speeds up the vertex linearizations in most cases, as we
  // never linearize the SAME slave or master vertex more than once)

  // number of nodes
  const int nsrows = LineElement()->NumNode();
  const int nmrows = SurfaceElement().NumNode();

  // prepare storage for slave and master linearizations
  std::vector<std::vector<GEN::pairedvector<int, double>>> linsnodes(
      nsrows, std::vector<GEN::pairedvector<int, double>>(3,
                  100 + linsize_ + 3 * LineElement()->NumNode() + 3 * SurfaceElement().NumNode()));
  std::vector<std::vector<GEN::pairedvector<int, double>>> linmnodes(
      nmrows, std::vector<GEN::pairedvector<int, double>>(3,
                  100 + linsize_ + 3 * LineElement()->NumNode() + 3 * SurfaceElement().NumNode()));

  // compute slave linearizations (nsrows)
  SlaveVertexLinearization(linsnodes);

  // compute master linearizations (nmrows)
  MasterVertexLinearization(linmnodes);

  //**********************************************************************
  // Line vertex linearization
  //**********************************************************************
  // loop over all clip Intersections vertices
  for (int i = 0; i < (int)InterSections().size(); ++i)
  {
    // references to current vertex and its linearization
    MORTAR::Vertex& currv = InterSections()[i];
    std::vector<GEN::pairedvector<int, double>>& currlin = linvertex[i];

    // decision on vertex type (slave, projmaster, linclip)
    if (currv.VType() == MORTAR::Vertex::projslave)
    {
      // get corresponding slave id
      int sid = currv.Nodeids()[0];

      // find corresponding slave node linearization
      int k = 0;
      while (k < nsrows)
      {
        if (LineElement()->NodeIds()[k] == sid) break;
        ++k;
      }

      // dserror if not found
      if (k == nsrows) dserror("ERROR: Slave Id not found!");

      // get the correct slave node linearization
      currlin = linsnodes[k];
    }
    else if (currv.VType() == MORTAR::Vertex::master)
    {
      // get corresponding master id
      int mid = currv.Nodeids()[0];

      // find corresponding master node linearization
      int k = 0;
      while (k < nmrows)
      {
        if (SurfaceElement().NodeIds()[k] == mid) break;
        ++k;
      }

      // dserror if not found
      if (k == nmrows) dserror("ERROR: Master Id not found!");

      // get the correct master node linearization
      currlin = linmnodes[k];
    }
    else if (currv.VType() == MORTAR::Vertex::lineclip)
    {
      // get references to the two slave vertices
      int sindex1 = -1;
      int sindex2 = -1;
      for (int j = 0; j < (int)SlaveVertices().size(); ++j)
      {
        if (SlaveVertices()[j].Nodeids()[0] == currv.Nodeids()[0]) sindex1 = j;
        if (SlaveVertices()[j].Nodeids()[0] == currv.Nodeids()[1]) sindex2 = j;
      }
      if (sindex1 < 0 || sindex2 < 0 || sindex1 == sindex2)
        dserror("ERROR: Lineclip linearization: (S) Something went wrong!");

      MORTAR::Vertex* sv1 = &SlaveVertices()[sindex1];
      MORTAR::Vertex* sv2 = &SlaveVertices()[sindex2];

      // get references to the two master vertices
      int mindex1 = -1;
      int mindex2 = -1;
      for (int j = 0; j < (int)MasterVertices().size(); ++j)
      {
        if (MasterVertices()[j].Nodeids()[0] == currv.Nodeids()[2]) mindex1 = j;
        if (MasterVertices()[j].Nodeids()[0] == currv.Nodeids()[3]) mindex2 = j;
      }
      if (mindex1 < 0 || mindex2 < 0 || mindex1 == mindex2)
        dserror("ERROR: Lineclip linearization: (M) Something went wrong!");

      MORTAR::Vertex* mv1 = &MasterVertices()[mindex1];
      MORTAR::Vertex* mv2 = &MasterVertices()[mindex2];

      // do lineclip vertex linearization
      LineclipVertexLinearization(currv, currlin, sv1, sv2, mv1, mv2, linsnodes, linmnodes);
    }

    else
      dserror("ERROR: VertexLinearization: Invalid Vertex Type!");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Linearization of lineclip vertex (3D) AuxPlane            popp 03/09|
 *----------------------------------------------------------------------*/
void CONTACT::LineToSurfaceCoupling3d::LineclipVertexLinearization(MORTAR::Vertex& currv,
    std::vector<GEN::pairedvector<int, double>>& currlin, MORTAR::Vertex* sv1, MORTAR::Vertex* sv2,
    MORTAR::Vertex* mv1, MORTAR::Vertex* mv2,
    std::vector<std::vector<GEN::pairedvector<int, double>>>& linsnodes,
    std::vector<std::vector<GEN::pairedvector<int, double>>>& linmnodes)
{
  // number of nodes
  const int nsrows = LineElement()->NumNode();
  const int nmrows = SurfaceElement().NumNode();

  // iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // compute factor Z
  double crossZ[3] = {0.0, 0.0, 0.0};
  crossZ[0] = (sv1->Coord()[1] - mv1->Coord()[1]) * (mv2->Coord()[2] - mv1->Coord()[2]) -
              (sv1->Coord()[2] - mv1->Coord()[2]) * (mv2->Coord()[1] - mv1->Coord()[1]);
  crossZ[1] = (sv1->Coord()[2] - mv1->Coord()[2]) * (mv2->Coord()[0] - mv1->Coord()[0]) -
              (sv1->Coord()[0] - mv1->Coord()[0]) * (mv2->Coord()[2] - mv1->Coord()[2]);
  crossZ[2] = (sv1->Coord()[0] - mv1->Coord()[0]) * (mv2->Coord()[1] - mv1->Coord()[1]) -
              (sv1->Coord()[1] - mv1->Coord()[1]) * (mv2->Coord()[0] - mv1->Coord()[0]);
  double Zfac = crossZ[0] * Auxn()[0] + crossZ[1] * Auxn()[1] + crossZ[2] * Auxn()[2];

  // compute factor N
  double crossN[3] = {0.0, 0.0, 0.0};
  crossN[0] = (sv2->Coord()[1] - sv1->Coord()[1]) * (mv2->Coord()[2] - mv1->Coord()[2]) -
              (sv2->Coord()[2] - sv1->Coord()[2]) * (mv2->Coord()[1] - mv1->Coord()[1]);
  crossN[1] = (sv2->Coord()[2] - sv1->Coord()[2]) * (mv2->Coord()[0] - mv1->Coord()[0]) -
              (sv2->Coord()[0] - sv1->Coord()[0]) * (mv2->Coord()[2] - mv1->Coord()[2]);
  crossN[2] = (sv2->Coord()[0] - sv1->Coord()[0]) * (mv2->Coord()[1] - mv1->Coord()[1]) -
              (sv2->Coord()[1] - sv1->Coord()[1]) * (mv2->Coord()[0] - mv1->Coord()[0]);
  double Nfac = crossN[0] * Auxn()[0] + crossN[1] * Auxn()[1] + crossN[2] * Auxn()[2];

  // slave edge vector
  double sedge[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; ++k) sedge[k] = sv2->Coord()[k] - sv1->Coord()[k];

  // prepare linearization derivZ
  double crossdZ1[3] = {0.0, 0.0, 0.0};
  double crossdZ2[3] = {0.0, 0.0, 0.0};
  double crossdZ3[3] = {0.0, 0.0, 0.0};
  crossdZ1[0] = (mv2->Coord()[1] - mv1->Coord()[1]) * Auxn()[2] -
                (mv2->Coord()[2] - mv1->Coord()[2]) * Auxn()[1];
  crossdZ1[1] = (mv2->Coord()[2] - mv1->Coord()[2]) * Auxn()[0] -
                (mv2->Coord()[0] - mv1->Coord()[0]) * Auxn()[2];
  crossdZ1[2] = (mv2->Coord()[0] - mv1->Coord()[0]) * Auxn()[1] -
                (mv2->Coord()[1] - mv1->Coord()[1]) * Auxn()[0];
  crossdZ2[0] = Auxn()[1] * (sv1->Coord()[2] - mv1->Coord()[2]) -
                Auxn()[2] * (sv1->Coord()[1] - mv1->Coord()[1]);
  crossdZ2[1] = Auxn()[2] * (sv1->Coord()[0] - mv1->Coord()[0]) -
                Auxn()[0] * (sv1->Coord()[2] - mv1->Coord()[2]);
  crossdZ2[2] = Auxn()[0] * (sv1->Coord()[1] - mv1->Coord()[1]) -
                Auxn()[1] * (sv1->Coord()[0] - mv1->Coord()[0]);
  crossdZ3[0] = (sv1->Coord()[1] - mv1->Coord()[1]) * (mv2->Coord()[2] - mv1->Coord()[2]) -
                (sv1->Coord()[2] - mv1->Coord()[2]) * (mv2->Coord()[1] - mv1->Coord()[1]);
  crossdZ3[1] = (sv1->Coord()[2] - mv1->Coord()[2]) * (mv2->Coord()[0] - mv1->Coord()[0]) -
                (sv1->Coord()[0] - mv1->Coord()[0]) * (mv2->Coord()[2] - mv1->Coord()[2]);
  crossdZ3[2] = (sv1->Coord()[0] - mv1->Coord()[0]) * (mv2->Coord()[1] - mv1->Coord()[1]) -
                (sv1->Coord()[1] - mv1->Coord()[1]) * (mv2->Coord()[0] - mv1->Coord()[0]);

  // prepare linearization derivN
  double crossdN1[3] = {0.0, 0.0, 0.0};
  double crossdN2[3] = {0.0, 0.0, 0.0};
  double crossdN3[3] = {0.0, 0.0, 0.0};
  crossdN1[0] = (mv2->Coord()[1] - mv1->Coord()[1]) * Auxn()[2] -
                (mv2->Coord()[2] - mv1->Coord()[2]) * Auxn()[1];
  crossdN1[1] = (mv2->Coord()[2] - mv1->Coord()[2]) * Auxn()[0] -
                (mv2->Coord()[0] - mv1->Coord()[0]) * Auxn()[2];
  crossdN1[2] = (mv2->Coord()[0] - mv1->Coord()[0]) * Auxn()[1] -
                (mv2->Coord()[1] - mv1->Coord()[1]) * Auxn()[0];
  crossdN2[0] = Auxn()[1] * (sv2->Coord()[2] - sv1->Coord()[2]) -
                Auxn()[2] * (sv2->Coord()[1] - sv1->Coord()[1]);
  crossdN2[1] = Auxn()[2] * (sv2->Coord()[0] - sv1->Coord()[0]) -
                Auxn()[0] * (sv2->Coord()[2] - sv1->Coord()[2]);
  crossdN2[2] = Auxn()[0] * (sv2->Coord()[1] - sv1->Coord()[1]) -
                Auxn()[1] * (sv2->Coord()[0] - sv1->Coord()[0]);
  crossdN3[0] = (sv2->Coord()[1] - sv1->Coord()[1]) * (mv2->Coord()[2] - mv1->Coord()[2]) -
                (sv2->Coord()[2] - sv1->Coord()[2]) * (mv2->Coord()[1] - mv1->Coord()[1]);
  crossdN3[1] = (sv2->Coord()[2] - sv1->Coord()[2]) * (mv2->Coord()[0] - mv1->Coord()[0]) -
                (sv2->Coord()[0] - sv1->Coord()[0]) * (mv2->Coord()[2] - mv1->Coord()[2]);
  crossdN3[2] = (sv2->Coord()[0] - sv1->Coord()[0]) * (mv2->Coord()[1] - mv1->Coord()[1]) -
                (sv2->Coord()[1] - sv1->Coord()[1]) * (mv2->Coord()[0] - mv1->Coord()[0]);

  // slave vertex linearization (2x)
  int sid1 = currv.Nodeids()[0];
  int sid2 = currv.Nodeids()[1];

  // find corresponding slave node linearizations
  int k = 0;
  while (k < nsrows)
  {
    if (LineElement()->NodeIds()[k] == sid1) break;
    ++k;
  }

  // dserror if not found
  if (k == nsrows) dserror("ERROR: Slave Id1 not found!");

  // get the correct slave node linearization
  std::vector<GEN::pairedvector<int, double>>& slavelin0 = linsnodes[k];

  k = 0;
  while (k < nsrows)
  {
    if (LineElement()->NodeIds()[k] == sid2) break;
    ++k;
  }

  // dserror if not found
  if (k == nsrows) dserror("ERROR: Slave Id2 not found!");

  // get the correct slave node linearization
  std::vector<GEN::pairedvector<int, double>>& slavelin1 = linsnodes[k];

  // master vertex linearization (2x)
  int mid1 = currv.Nodeids()[2];
  int mid2 = currv.Nodeids()[3];

  // find corresponding master node linearizations
  k = 0;
  while (k < nmrows)
  {
    if (SurfaceElement().NodeIds()[k] == mid1) break;
    ++k;
  }

  // dserror if not found
  if (k == nmrows) dserror("ERROR: Master Id1 not found!");

  // get the correct master node linearization
  std::vector<GEN::pairedvector<int, double>>& masterlin0 = linmnodes[k];

  k = 0;
  while (k < nmrows)
  {
    if (SurfaceElement().NodeIds()[k] == mid2) break;
    ++k;
  }

  // dserror if not found
  if (k == nmrows) dserror("ERROR: Master Id2 not found!");

  // get the correct master node linearization
  std::vector<GEN::pairedvector<int, double>>& masterlin1 = linmnodes[k];

  // linearization of element normal Auxn()
  std::vector<GEN::pairedvector<int, double>>& linauxn = GetDerivAuxn();

  const double ZNfac = Zfac / Nfac;
  const double ZNNfac = Zfac / (Nfac * Nfac);
  const double Nfacinv = 1.0 / Nfac;

  // bring everything together -> lineclip vertex linearization
  for (int k = 0; k < 3; ++k)
  {
    for (_CI p = slavelin0[k].begin(); p != slavelin0[k].end(); ++p)
    {
      currlin[k][p->first] += (p->second);
      currlin[k][p->first] += ZNfac * (p->second);
      for (int dim = 0; dim < 3; ++dim)
      {
        currlin[dim][p->first] -= sedge[dim] * Nfacinv * crossdZ1[k] * (p->second);
        currlin[dim][p->first] -= sedge[dim] * ZNNfac * crossdN1[k] * (p->second);
      }
    }
    for (_CI p = slavelin1[k].begin(); p != slavelin1[k].end(); ++p)
    {
      currlin[k][p->first] -= ZNfac * (p->second);
      for (int dim = 0; dim < 3; ++dim)
      {
        currlin[dim][p->first] += sedge[dim] * ZNNfac * crossdN1[k] * (p->second);
      }
    }
    for (_CI p = masterlin0[k].begin(); p != masterlin0[k].end(); ++p)
    {
      for (int dim = 0; dim < 3; ++dim)
      {
        currlin[dim][p->first] += sedge[dim] * Nfacinv * crossdZ1[k] * (p->second);
        currlin[dim][p->first] += sedge[dim] * Nfacinv * crossdZ2[k] * (p->second);
        currlin[dim][p->first] -= sedge[dim] * ZNNfac * crossdN2[k] * (p->second);
      }
    }
    for (_CI p = masterlin1[k].begin(); p != masterlin1[k].end(); ++p)
    {
      for (int dim = 0; dim < 3; ++dim)
      {
        currlin[dim][p->first] -= sedge[dim] * Nfacinv * crossdZ2[k] * (p->second);
        currlin[dim][p->first] += sedge[dim] * ZNNfac * crossdN2[k] * (p->second);
      }
    }
    for (_CI p = linauxn[k].begin(); p != linauxn[k].end(); ++p)
    {
      for (int dim = 0; dim < 3; ++dim)
      {
        currlin[dim][p->first] -= sedge[dim] * Nfacinv * crossdZ3[k] * (p->second);
        currlin[dim][p->first] += sedge[dim] * ZNNfac * crossdN3[k] * (p->second);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  compute and check length of intLine                      farah 07/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineToSurfaceCoupling3d::CheckLength()
{
  // tolerance
  const double sminedge = ParentElement().MinEdgeSize();
  const double mminedge = SurfaceElement().MinEdgeSize();
  const double tol = MORTARCLIPTOL * std::min(sminedge, mminedge);

  // compute distance vector
  double v[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i) v[i] = InterSections()[0].Coord()[i] - InterSections()[1].Coord()[i];

  // compute length
  double length = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

  // check
  if (length < tol) return false;

  return true;
}

/*----------------------------------------------------------------------*
 |  eval (public)                                            farah 07/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineToSurfaceCoupling3d::AuxiliaryPlane()
{
  // we first need the element center:
  // for quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  double loccenter[2] = {0.0, 0.0};

  DRT::Element::DiscretizationType dt = SurfaceElement().Shape();
  if (dt == DRT::Element::tri3 || dt == DRT::Element::tri6)
  {
    loccenter[0] = 1.0 / 3.0;
    loccenter[1] = 1.0 / 3.0;
  }
  else if (dt == DRT::Element::quad4 || dt == DRT::Element::quad8 || dt == DRT::Element::quad9)
  {
    loccenter[0] = 0.0;
    loccenter[1] = 0.0;
  }
  else
    dserror("ERROR: AuxiliaryPlane called for unknown element type");

  // compute element center via shape fct. interpolation
  SurfaceElement().LocalToGlobal(loccenter, Auxc(), 0);

  // we then compute the unit normal vector at the element center
  Lauxn() = SurfaceElement().ComputeUnitNormalAtXi(loccenter, AuxnSurf());
  //
  //  // compute aux normal linearization
  //  SurfaceElement().DerivUnitNormalAtXi(loccenter, GetDerivAuxn());

  // bye
  return true;
}

/*----------------------------------------------------------------------*
 |  create auxiliary line + normal                           farah 08/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineToSurfaceCoupling3d::AuxiliaryLine()
{
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  int nnodes = LineElement()->NumNode();
  if (nnodes != 2) dserror("ERROR: Auxiliary line calculation only for line2 elements!");

  // average nodal normals of line element
  linsize_ = 0;
  for (int i = 0; i < nnodes; ++i)
  {
    DRT::Node* node = idiscret_.gNode(LineElement()->NodeIds()[i]);
    if (!node) dserror("ERROR: Cannot find slave element with gid %", LineElement()->NodeIds()[i]);
    CoNode* mycnode = dynamic_cast<CoNode*>(node);
    if (!mycnode) dserror("ERROR: ProjectSlave: Null pointer!");

    linsize_ += mycnode->GetLinsize();
  }

  // TODO: this is a safety scaling. Correct linsize
  //      should be predicted
  linsize_ *= 100;

  // auxiliary normal
  GetDerivAuxn().resize(3, linsize_ * 10);

  // auxiliary center
  GetDerivAuxc().resize(3, linsize_ * 10);

  Auxc()[0] = 0.0;
  Auxc()[1] = 0.0;
  Auxc()[2] = 0.0;

  std::vector<GEN::pairedvector<int, double>> dauxn(3, 100);

  // average nodal normals of line element
  for (int i = 0; i < nnodes; ++i)
  {
    DRT::Node* node = idiscret_.gNode(LineElement()->NodeIds()[i]);
    if (!node) dserror("ERROR: Cannot find slave element with gid %", LineElement()->NodeIds()[i]);
    CoNode* mycnode = dynamic_cast<CoNode*>(node);
    if (!mycnode) dserror("ERROR: ProjectSlave: Null pointer!");

    Auxn()[0] += 0.5 * mycnode->MoData().n()[0];
    Auxn()[1] += 0.5 * mycnode->MoData().n()[1];
    Auxn()[2] += 0.5 * mycnode->MoData().n()[2];

    for (_CI p = mycnode->CoData().GetDerivN()[0].begin();
         p != mycnode->CoData().GetDerivN()[0].end(); ++p)
      (dauxn[0])[p->first] += 0.5 * (p->second);
    for (_CI p = mycnode->CoData().GetDerivN()[1].begin();
         p != mycnode->CoData().GetDerivN()[1].end(); ++p)
      (dauxn[1])[p->first] += 0.5 * (p->second);
    for (_CI p = mycnode->CoData().GetDerivN()[2].begin();
         p != mycnode->CoData().GetDerivN()[2].end(); ++p)
      (dauxn[2])[p->first] += 0.5 * (p->second);

    // new aux center
    for (int dim = 0; dim < Dim(); ++dim) Auxc()[dim] += 0.5 * mycnode->xspatial()[dim];

    (GetDerivAuxc()[0])[mycnode->Dofs()[0]] += 0.5;
    (GetDerivAuxc()[1])[mycnode->Dofs()[1]] += 0.5;
    (GetDerivAuxc()[2])[mycnode->Dofs()[2]] += 0.5;
  }

  // create tangent of line element
  double tangent[3] = {0.0, 0.0, 0.0};
  DRT::Node* node = idiscret_.gNode(LineElement()->NodeIds()[0]);
  if (!node) dserror("ERROR: Cannot find slave element with gid %", LineElement()->NodeIds()[0]);
  CoNode* mycnode = dynamic_cast<CoNode*>(node);
  if (!mycnode) dserror("ERROR: ProjectSlave: Null pointer!");

  tangent[0] += mycnode->xspatial()[0];
  tangent[1] += mycnode->xspatial()[1];
  tangent[2] += mycnode->xspatial()[2];

  DRT::Node* node2 = idiscret_.gNode(LineElement()->NodeIds()[1]);
  if (!node2) dserror("ERROR: Cannot find slave element with gid %", LineElement()->NodeIds()[1]);
  CoNode* mycnode2 = dynamic_cast<CoNode*>(node2);
  if (!mycnode2) dserror("ERROR: ProjectSlave: Null pointer!");

  tangent[0] -= mycnode2->xspatial()[0];
  tangent[1] -= mycnode2->xspatial()[1];
  tangent[2] -= mycnode2->xspatial()[2];

  Epetra_SerialDenseMatrix tanplane(3, 3);
  tanplane(0, 0) = 1 - (tangent[0] * tangent[0]);
  tanplane(0, 1) = -(tangent[0] * tangent[1]);
  tanplane(0, 2) = -(tangent[0] * tangent[2]);
  tanplane(1, 0) = -(tangent[1] * tangent[0]);
  tanplane(1, 1) = 1 - (tangent[1] * tangent[1]);
  tanplane(1, 2) = -(tangent[1] * tangent[2]);
  tanplane(2, 0) = -(tangent[2] * tangent[0]);
  tanplane(2, 1) = -(tangent[2] * tangent[1]);
  tanplane(2, 2) = 1 - (tangent[2] * tangent[2]);

  double finalauxn[3] = {0.0, 0.0, 0.0};
  finalauxn[0] =
      tanplane(0, 0) * Auxn()[0] + tanplane(0, 1) * Auxn()[1] + tanplane(0, 2) * Auxn()[2];
  finalauxn[1] =
      tanplane(1, 0) * Auxn()[0] + tanplane(1, 1) * Auxn()[1] + tanplane(1, 2) * Auxn()[2];
  finalauxn[2] =
      tanplane(2, 0) * Auxn()[0] + tanplane(2, 1) * Auxn()[1] + tanplane(2, 2) * Auxn()[2];

  // lin tangent
  std::vector<GEN::pairedvector<int, double>> dnmap_unit(3, linsize_ * 10);
  for (int i = 0; i < Dim(); ++i)
  {
    dnmap_unit[i][mycnode->Dofs()[i]] += 1;
    dnmap_unit[i][mycnode2->Dofs()[i]] -= 1;
  }


  std::vector<GEN::pairedvector<int, double>> tplanex(3, linsize_ * 10);
  std::vector<GEN::pairedvector<int, double>> tplaney(3, linsize_ * 10);
  std::vector<GEN::pairedvector<int, double>> tplanez(3, linsize_ * 10);

  for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
    tplanex[0][p->first] -= tangent[0] * p->second;
  for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
    tplanex[1][p->first] -= tangent[1] * p->second;
  for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
    tplanex[2][p->first] -= tangent[2] * p->second;

  for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
    tplaney[0][p->first] -= tangent[0] * p->second;
  for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
    tplaney[1][p->first] -= tangent[1] * p->second;
  for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
    tplaney[2][p->first] -= tangent[2] * p->second;

  for (_CI p = dnmap_unit[2].begin(); p != dnmap_unit[2].end(); ++p)
    tplanez[0][p->first] -= tangent[0] * p->second;
  for (_CI p = dnmap_unit[2].begin(); p != dnmap_unit[2].end(); ++p)
    tplanez[1][p->first] -= tangent[1] * p->second;
  for (_CI p = dnmap_unit[2].begin(); p != dnmap_unit[2].end(); ++p)
    tplanez[2][p->first] -= tangent[2] * p->second;

  //------------
  for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
    tplanex[0][p->first] -= tangent[0] * p->second;
  for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
    tplanex[1][p->first] -= tangent[0] * p->second;
  for (_CI p = dnmap_unit[2].begin(); p != dnmap_unit[2].end(); ++p)
    tplanex[2][p->first] -= tangent[0] * p->second;

  for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
    tplaney[0][p->first] -= tangent[1] * p->second;
  for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
    tplaney[1][p->first] -= tangent[1] * p->second;
  for (_CI p = dnmap_unit[2].begin(); p != dnmap_unit[2].end(); ++p)
    tplaney[2][p->first] -= tangent[1] * p->second;

  for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
    tplanez[0][p->first] -= tangent[2] * p->second;
  for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
    tplanez[1][p->first] -= tangent[2] * p->second;
  for (_CI p = dnmap_unit[2].begin(); p != dnmap_unit[2].end(); ++p)
    tplanez[2][p->first] -= tangent[2] * p->second;



  for (_CI p = dauxn[0].begin(); p != dauxn[0].end(); ++p)
    GetDerivAuxn()[0][p->first] += tanplane(0, 0) * p->second;
  for (_CI p = dauxn[1].begin(); p != dauxn[1].end(); ++p)
    GetDerivAuxn()[0][p->first] += tanplane(0, 1) * p->second;
  for (_CI p = dauxn[2].begin(); p != dauxn[2].end(); ++p)
    GetDerivAuxn()[0][p->first] += tanplane(0, 2) * p->second;

  for (_CI p = dauxn[0].begin(); p != dauxn[0].end(); ++p)
    GetDerivAuxn()[1][p->first] += tanplane(1, 0) * p->second;
  for (_CI p = dauxn[1].begin(); p != dauxn[1].end(); ++p)
    GetDerivAuxn()[1][p->first] += tanplane(1, 1) * p->second;
  for (_CI p = dauxn[2].begin(); p != dauxn[2].end(); ++p)
    GetDerivAuxn()[1][p->first] += tanplane(1, 2) * p->second;

  for (_CI p = dauxn[0].begin(); p != dauxn[0].end(); ++p)
    GetDerivAuxn()[2][p->first] += tanplane(2, 0) * p->second;
  for (_CI p = dauxn[1].begin(); p != dauxn[1].end(); ++p)
    GetDerivAuxn()[2][p->first] += tanplane(2, 1) * p->second;
  for (_CI p = dauxn[2].begin(); p != dauxn[2].end(); ++p)
    GetDerivAuxn()[2][p->first] += tanplane(2, 2) * p->second;

  //-----------------------------
  for (_CI p = tplanex[0].begin(); p != tplanex[0].end(); ++p)
    GetDerivAuxn()[0][p->first] += Auxn()[0] * p->second;
  for (_CI p = tplanex[1].begin(); p != tplanex[1].end(); ++p)
    GetDerivAuxn()[0][p->first] += Auxn()[1] * p->second;
  for (_CI p = tplanex[2].begin(); p != tplanex[2].end(); ++p)
    GetDerivAuxn()[0][p->first] += Auxn()[2] * p->second;

  for (_CI p = tplaney[0].begin(); p != tplaney[0].end(); ++p)
    GetDerivAuxn()[1][p->first] += Auxn()[0] * p->second;
  for (_CI p = tplaney[1].begin(); p != tplaney[1].end(); ++p)
    GetDerivAuxn()[1][p->first] += Auxn()[1] * p->second;
  for (_CI p = tplaney[2].begin(); p != tplaney[2].end(); ++p)
    GetDerivAuxn()[1][p->first] += Auxn()[2] * p->second;

  for (_CI p = tplanez[0].begin(); p != tplanez[0].end(); ++p)
    GetDerivAuxn()[2][p->first] += Auxn()[0] * p->second;
  for (_CI p = tplanez[1].begin(); p != tplanez[1].end(); ++p)
    GetDerivAuxn()[2][p->first] += Auxn()[1] * p->second;
  for (_CI p = tplanez[2].begin(); p != tplanez[2].end(); ++p)
    GetDerivAuxn()[2][p->first] += Auxn()[2] * p->second;


  Auxn()[0] = finalauxn[0];
  Auxn()[1] = finalauxn[1];
  Auxn()[2] = finalauxn[2];

  AuxnSurf()[0] = -Auxn()[0];
  AuxnSurf()[1] = -Auxn()[1];
  AuxnSurf()[2] = -Auxn()[2];

  return true;
}

/*----------------------------------------------------------------------*
 |  eval (public)                                            farah 07/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineToSurfaceCoupling3d::HasProjStatus() { return true; }


/*----------------------------------------------------------------------*
 |  eval (public)                                            farah 07/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineToSurfaceCoupling3d::ProjectSlave()
{
  // project slave nodes onto auxiliary plane
  int nnodes = LineElement()->NumNode();

  // initialize storage for slave coords + their ids
  std::vector<double> vertices(3);
  std::vector<int> snodeids(1);

  for (int i = 0; i < nnodes; ++i)
  {
    DRT::Node* node = idiscret_.gNode(LineElement()->NodeIds()[i]);
    if (!node) dserror("ERROR: Cannot find slave element with gid %", LineElement()->NodeIds()[i]);
    CoNode* mycnode = dynamic_cast<CoNode*>(node);
    if (!mycnode) dserror("ERROR: ProjectSlave: Null pointer!");

    // first build difference of point and element center
    // and then dot product with unit normal at center
    const double dist = (mycnode->xspatial()[0] - Auxc()[0]) * Auxn()[0] +
                        (mycnode->xspatial()[1] - Auxc()[1]) * Auxn()[1] +
                        (mycnode->xspatial()[2] - Auxc()[2]) * Auxn()[2];

    // compute projection
    for (int k = 0; k < 3; ++k) vertices[k] = mycnode->xspatial()[k] - dist * Auxn()[k];

    // get node id, too
    snodeids[0] = mycnode->Id();

    // store into vertex data structure
    SlaveVertices().push_back(MORTAR::Vertex(
        vertices, MORTAR::Vertex::projslave, snodeids, NULL, NULL, false, false, NULL, -1.0));
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Linearization of slave vertex (3D) AuxPlane              farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineToSurfaceCoupling3d::SlaveVertexLinearization(
    std::vector<std::vector<GEN::pairedvector<int, double>>>& currlin)
{
  // we first need the slave element center:
  // for quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  double scxi[2];

  DRT::Element::DiscretizationType dt = SurfaceElement().Shape();
  if (dt == MORTAR::MortarElement::tri3 || dt == MORTAR::MortarElement::tri6)
  {
    scxi[0] = 1.0 / 3.0;
    scxi[1] = 1.0 / 3.0;
  }
  else if (dt == MORTAR::MortarElement::quad4 || dt == MORTAR::MortarElement::quad8 ||
           dt == MORTAR::MortarElement::quad9)
  {
    scxi[0] = 0.0;
    scxi[1] = 0.0;
  }
  else
    dserror("ERROR: MasterVertexLinearization called for unknown element type");

  // evlauate shape functions + derivatives at scxi
  int nrow = SurfaceElement().NumNode();
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow, 2, true);
  SurfaceElement().EvaluateShape(scxi, sval, sderiv, nrow);

  // we need all participating slave nodes
  DRT::Node** snodes = SurfaceElement().Nodes();
  std::vector<MORTAR::MortarNode*> smrtrnodes(nrow);

  for (int i = 0; i < nrow; ++i)
  {
    smrtrnodes[i] = dynamic_cast<MORTAR::MortarNode*>(snodes[i]);
    if (!smrtrnodes[i]) dserror("ERROR: MasterVertexLinearization: Null pointer!");
  }

  // linearization of the SlaveIntEle spatial coords
  std::vector<std::vector<GEN::pairedvector<int, double>>> snodelin(0);

  // resize the linearizations
  snodelin.resize(nrow, std::vector<GEN::pairedvector<int, double>>(3, 1));

  // loop over all intEle nodes
  for (int in = 0; in < nrow; ++in)
    for (int dim = 0; dim < 3; ++dim) snodelin[in][dim][smrtrnodes[in]->Dofs()[dim]] += 1.;

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator
      _CI;  // linearization of element center Auxc()
  //  std::vector<GEN::pairedvector<int  ,double> > linauxc(3,10*SurfaceElement().NumNode()); //
  //  assume 3 dofs per node

  //  for (int i = 0; i < nrow; ++i)
  //      for (int dim=0; dim<3; ++dim)
  //        for (_CI p=snodelin[i][dim].begin(); p!=snodelin[i][dim].end(); ++p)
  //          linauxc[dim][p->first] = sval[i]*p->second;

  // linearization of element normal Auxn()
  std::vector<GEN::pairedvector<int, double>>& linauxn = GetDerivAuxn();

  // linearization of the MasterIntEle spatial coords
  std::vector<std::vector<GEN::pairedvector<int, double>>> mnodelin(0);

  // resize the linearizations
  mnodelin.resize(LineElement()->NumNode(), std::vector<GEN::pairedvector<int, double>>(3, 1));

  // loop over all intEle nodes
  for (int in = 0; in < LineElement()->NumNode(); ++in)
  {
    MORTAR::MortarNode* mrtrmnode =
        dynamic_cast<MORTAR::MortarNode*>(idiscret_.gNode(LineElement()->NodeIds()[in]));
    if (mrtrmnode == NULL) dserror("dynamic cast to mortar node went wrong");

    for (int dim = 0; dim < 3; ++dim) mnodelin[in][dim][mrtrmnode->Dofs()[dim]] += 1.;
  }

  // put everything together for slave vertex linearization
  // loop over all vertices
  for (int i = 0; i < LineElement()->NumNode(); ++i)
  {
    MORTAR::MortarNode* mrtrmnode =
        dynamic_cast<MORTAR::MortarNode*>(idiscret_.gNode(LineElement()->NodeIds()[i]));
    if (!mrtrmnode) dserror("cast to mortar node failed");

    // (1) slave node coordinates part
    for (_CI p = mnodelin[i][0].begin(); p != mnodelin[i][0].end(); ++p)
    {
      currlin[i][0][p->first] += (1.0 - Auxn()[0] * Auxn()[0]) * p->second;
      currlin[i][1][p->first] -= (Auxn()[0] * Auxn()[1]) * p->second;
      currlin[i][2][p->first] -= (Auxn()[0] * Auxn()[2]) * p->second;
    }
    for (_CI p = mnodelin[i][1].begin(); p != mnodelin[i][1].end(); ++p)
    {
      currlin[i][0][p->first] -= (Auxn()[0] * Auxn()[1]) * p->second;
      currlin[i][1][p->first] += (1.0 - Auxn()[1] * Auxn()[1]) * p->second;
      currlin[i][2][p->first] -= (Auxn()[1] * Auxn()[2]) * p->second;
    }
    for (_CI p = mnodelin[i][2].begin(); p != mnodelin[i][2].end(); ++p)
    {
      currlin[i][0][p->first] -= (Auxn()[2] * Auxn()[0]) * p->second;
      currlin[i][1][p->first] -= (Auxn()[2] * Auxn()[1]) * p->second;
      currlin[i][2][p->first] += (1.0 - Auxn()[2] * Auxn()[2]) * p->second;
    }

    // (2) slave element center coordinates (Auxc()) part
    for (_CI p = GetDerivAuxc()[0].begin(); p != GetDerivAuxc()[0].end(); ++p)
      for (int k = 0; k < 3; ++k) currlin[i][k][p->first] += Auxn()[0] * Auxn()[k] * (p->second);

    for (_CI p = GetDerivAuxc()[1].begin(); p != GetDerivAuxc()[1].end(); ++p)
      for (int k = 0; k < 3; ++k) currlin[i][k][p->first] += Auxn()[1] * Auxn()[k] * (p->second);

    for (_CI p = GetDerivAuxc()[2].begin(); p != GetDerivAuxc()[2].end(); ++p)
      for (int k = 0; k < 3; ++k) currlin[i][k][p->first] += Auxn()[2] * Auxn()[k] * (p->second);

    // (3) slave element normal (Auxn()) part
    double xdotn = (mrtrmnode->xspatial()[0] - Auxc()[0]) * Auxn()[0] +
                   (mrtrmnode->xspatial()[1] - Auxc()[1]) * Auxn()[1] +
                   (mrtrmnode->xspatial()[2] - Auxc()[2]) * Auxn()[2];

    for (_CI p = linauxn[0].begin(); p != linauxn[0].end(); ++p)
    {
      currlin[i][0][p->first] -= xdotn * (p->second);
      for (int k = 0; k < 3; ++k)
        currlin[i][k][p->first] -= (mrtrmnode->xspatial()[0] - Auxc()[0]) * Auxn()[k] * (p->second);
    }

    for (_CI p = linauxn[1].begin(); p != linauxn[1].end(); ++p)
    {
      currlin[i][1][p->first] -= xdotn * (p->second);
      for (int k = 0; k < 3; ++k)
        currlin[i][k][p->first] -= (mrtrmnode->xspatial()[1] - Auxc()[1]) * Auxn()[k] * (p->second);
    }

    for (_CI p = linauxn[2].begin(); p != linauxn[2].end(); ++p)
    {
      currlin[i][2][p->first] -= xdotn * (p->second);
      for (int k = 0; k < 3; ++k)
        currlin[i][k][p->first] -= (mrtrmnode->xspatial()[2] - Auxc()[2]) * Auxn()[k] * (p->second);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  eval (public)                                            farah 07/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineToSurfaceCoupling3d::ProjectMaster()
{
  // project master nodes onto auxiliary plane
  int nnodes = SurfaceElement().NumNode();
  DRT::Node** mynodes = SurfaceElement().Nodes();
  if (!mynodes) dserror("ERROR: ProjectMaster: Null pointer!");

  // initialize storage for master coords + their ids
  std::vector<double> vertices(3);
  std::vector<int> mnodeids(1);

  for (int i = 0; i < nnodes; ++i)
  {
    CoNode* mycnode = dynamic_cast<CoNode*>(mynodes[i]);
    if (!mycnode) dserror("ERROR: ProjectMaster: Null pointer!");

    // first build difference of point and element center
    // and then dot product with unit normal at center
    const double dist = (mycnode->xspatial()[0] - Auxc()[0]) * Auxn()[0] +
                        (mycnode->xspatial()[1] - Auxc()[1]) * Auxn()[1] +
                        (mycnode->xspatial()[2] - Auxc()[2]) * Auxn()[2];

    // compute projection
    for (int k = 0; k < 3; ++k) vertices[k] = mycnode->xspatial()[k] - dist * Auxn()[k];

    // get node id, too
    mnodeids[0] = mycnode->Id();

    // store into vertex data structure
    MasterVertices().push_back(MORTAR::Vertex(
        vertices, MORTAR::Vertex::master, mnodeids, NULL, NULL, false, false, NULL, -1.0));

    // std::cout << "->RealNode(M) " << mycnode->Id() << ": " << mycnode->xspatial()[0] << " " <<
    // mycnode->xspatial()[1] << " " << mycnode->xspatial()[2] << std::endl; std::cout <<
    // "->ProjNode(M) " << mycnode->Id() << ": " << vertices[0] << " " << vertices[1] << " " <<
    // vertices[2] << std::endl;
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Linearization of slave vertex (3D) AuxPlane               farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineToSurfaceCoupling3d::MasterVertexLinearization(
    std::vector<std::vector<GEN::pairedvector<int, double>>>& currlin)
{
  // we first need the slave element center:
  // for quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  double scxi[2];

  DRT::Element::DiscretizationType dt = SurfaceElement().Shape();
  if (dt == MORTAR::MortarElement::tri3 || dt == MORTAR::MortarElement::tri6)
  {
    scxi[0] = 1.0 / 3.0;
    scxi[1] = 1.0 / 3.0;
  }
  else if (dt == MORTAR::MortarElement::quad4 || dt == MORTAR::MortarElement::quad8 ||
           dt == MORTAR::MortarElement::quad9)
  {
    scxi[0] = 0.0;
    scxi[1] = 0.0;
  }
  else
    dserror("ERROR: SlaveVertexLinearization called for unknown element type");

  // evlauate shape functions + derivatives at scxi
  const int nrow = SurfaceElement().NumNode();
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow, 2, true);
  SurfaceElement().EvaluateShape(scxi, sval, sderiv, nrow);

  // we need all participating slave nodes
  DRT::Node** snodes = SurfaceElement().Nodes();
  std::vector<MORTAR::MortarNode*> smrtrnodes(nrow);

  for (int i = 0; i < nrow; ++i)
  {
    smrtrnodes[i] = dynamic_cast<MORTAR::MortarNode*>(snodes[i]);
    if (!smrtrnodes[i]) dserror("ERROR: SlaveVertexLinearization: Null pointer!");
  }

  // linearization of the IntEle spatial coords
  std::vector<std::vector<GEN::pairedvector<int, double>>> nodelin(0);

  // resize the linearizations
  nodelin.resize(nrow, std::vector<GEN::pairedvector<int, double>>(3, 1));

  // loop over all intEle nodes
  for (int in = 0; in < nrow; ++in)
    for (int dim = 0; dim < 3; ++dim) nodelin[in][dim][smrtrnodes[in]->Dofs()[dim]] += 1.;

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator
      _CI;  // linearization of element center Auxc()
  //  std  ::vector<GEN::pairedvector<int  ,double> > linauxc(3,SurfaceElement().NumNode()); //
  //  assume 3 dofs per node
  //
  //  for (int i = 0; i < nrow; ++i)
  //      for (int dim=0; dim<3; ++dim)
  //        for (_CI p=nodelin[i][dim].begin(); p!=nodelin[i][dim].end(); ++p)
  //          linauxc[dim][p->first] = sval[i]*p->second;

  // linearization of element normal Auxn()
  std::vector<GEN::pairedvector<int, double>>& linauxn = GetDerivAuxn();

  // put everything together for slave vertex linearization
  // loop over all vertices
  for (int i = 0; i < SurfaceElement().NumNode(); ++i)
  {
    MORTAR::MortarNode* mrtrsnode = dynamic_cast<MORTAR::MortarNode*>(SurfaceElement().Nodes()[i]);
    if (!mrtrsnode) dserror("cast to mortar node failed");

    // (1) slave node coordinates part
    for (_CI p = nodelin[i][0].begin(); p != nodelin[i][0].end(); ++p)
    {
      currlin[i][0][p->first] += (1.0 - Auxn()[0] * Auxn()[0]) * p->second;
      currlin[i][1][p->first] -= (Auxn()[0] * Auxn()[1]) * p->second;
      currlin[i][2][p->first] -= (Auxn()[0] * Auxn()[2]) * p->second;
    }
    for (_CI p = nodelin[i][1].begin(); p != nodelin[i][1].end(); ++p)
    {
      currlin[i][0][p->first] -= (Auxn()[0] * Auxn()[1]) * p->second;
      currlin[i][1][p->first] += (1.0 - Auxn()[1] * Auxn()[1]) * p->second;
      currlin[i][2][p->first] -= (Auxn()[1] * Auxn()[2]) * p->second;
    }
    for (_CI p = nodelin[i][2].begin(); p != nodelin[i][2].end(); ++p)
    {
      currlin[i][0][p->first] -= (Auxn()[2] * Auxn()[0]) * p->second;
      currlin[i][1][p->first] -= (Auxn()[2] * Auxn()[1]) * p->second;
      currlin[i][2][p->first] += (1.0 - Auxn()[2] * Auxn()[2]) * p->second;
    }

    // (2) slave element center coordinates (Auxc()) part
    for (_CI p = GetDerivAuxc()[0].begin(); p != GetDerivAuxc()[0].end(); ++p)
      for (int k = 0; k < 3; ++k) currlin[i][k][p->first] += Auxn()[0] * Auxn()[k] * (p->second);

    for (_CI p = GetDerivAuxc()[1].begin(); p != GetDerivAuxc()[1].end(); ++p)
      for (int k = 0; k < 3; ++k) currlin[i][k][p->first] += Auxn()[1] * Auxn()[k] * (p->second);

    for (_CI p = GetDerivAuxc()[2].begin(); p != GetDerivAuxc()[2].end(); ++p)
      for (int k = 0; k < 3; ++k) currlin[i][k][p->first] += Auxn()[2] * Auxn()[k] * (p->second);

    // (3) slave element normal (Auxn()) part
    double xdotn = (mrtrsnode->xspatial()[0] - Auxc()[0]) * Auxn()[0] +
                   (mrtrsnode->xspatial()[1] - Auxc()[1]) * Auxn()[1] +
                   (mrtrsnode->xspatial()[2] - Auxc()[2]) * Auxn()[2];

    for (_CI p = linauxn[0].begin(); p != linauxn[0].end(); ++p)
    {
      currlin[i][0][p->first] -= xdotn * (p->second);
      for (int k = 0; k < 3; ++k)
        currlin[i][k][p->first] -= (mrtrsnode->xspatial()[0] - Auxc()[0]) * Auxn()[k] * (p->second);
    }

    for (_CI p = linauxn[1].begin(); p != linauxn[1].end(); ++p)
    {
      currlin[i][1][p->first] -= xdotn * (p->second);
      for (int k = 0; k < 3; ++k)
        currlin[i][k][p->first] -= (mrtrsnode->xspatial()[1] - Auxc()[1]) * Auxn()[k] * (p->second);
    }

    for (_CI p = linauxn[2].begin(); p != linauxn[2].end(); ++p)
    {
      currlin[i][2][p->first] -= xdotn * (p->second);
      for (int k = 0; k < 3; ++k)
        currlin[i][k][p->first] -= (mrtrsnode->xspatial()[2] - Auxc()[2]) * Auxn()[k] * (p->second);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  get communicator                                         farah 07/16|
 *----------------------------------------------------------------------*/
const Epetra_Comm& CONTACT::LineToSurfaceCoupling3d::Comm() const { return idiscret_.Comm(); }


/*----------------------------------------------------------------------*
 |  ctor for ltl (public)                                    farah 07/16|
 *----------------------------------------------------------------------*/
CONTACT::LineToLineCouplingPoint3d::LineToLineCouplingPoint3d(DRT::Discretization& idiscret,
    int dim, Teuchos::ParameterList& params, Teuchos::RCP<MORTAR::MortarElement>& lsele,
    Teuchos::RCP<MORTAR::MortarElement>& lmele)
    : idiscret_(idiscret), dim_(dim), imortar_(params), lSele_(lsele), lMele_(lmele)
{
  // empty constructor
}


/*----------------------------------------------------------------------*
 |  get communicator                                         farah 07/16|
 *----------------------------------------------------------------------*/
const Epetra_Comm& CONTACT::LineToLineCouplingPoint3d::Comm() const { return idiscret_.Comm(); }

/*----------------------------------------------------------------------*
 |  eval                                                     farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineToLineCouplingPoint3d::EvaluateCoupling()
{
  // 1. check parallelity
  bool parallel = CheckParallelity();
  if (parallel) return;

  // 2. calc intersection
  // create empty points
  double sxi[1] = {0.0};
  double mxi[1] = {0.0};

  // create empty lin vectors
  GEN::pairedvector<int, double> dsxi(
      100 + 3 * LineMasterElement()->NumNode() + 3 * LineSlaveElement()->NumNode());
  GEN::pairedvector<int, double> dmxi(
      100 + 3 * LineMasterElement()->NumNode() + 3 * LineSlaveElement()->NumNode());
  LineIntersection(sxi, mxi, dsxi, dmxi);

  // 3. check solution
  bool valid = CheckIntersection(sxi, mxi);
  if (!valid) return;

  // 4. check if intersection was already done!
  for (int i = 0; i < LineSlaveElement()->NumNode(); ++i)
  {
    if (dynamic_cast<CoNode*>(LineSlaveElement()->Nodes()[i])->MoData().GetDltl().size() > 0)
      return;
  }

  // 5. evaluate terms
  EvaluateTerms(sxi, mxi, dsxi, dmxi);

  return;
}


/*----------------------------------------------------------------------*
 |  line-line intersection                                   farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineToLineCouplingPoint3d::EvaluateTerms(double* sxi, double* mxi,
    GEN::pairedvector<int, double>& dsxi, GEN::pairedvector<int, double>& dmxi)
{
  bool friction = false;
  INPAR::CONTACT::FrictionType ftype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(IParams(), "FRICTION");
  if (ftype != INPAR::CONTACT::friction_none) friction = true;

  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = LineSlaveElement()->Nodes();
  if (!mynodes) dserror("ERROR: IntegrateDerivCell3DAuxPlaneLTS: Null pointer!");
  DRT::Node** mnodes = LineMasterElement()->Nodes();
  if (!mnodes) dserror("ERROR: IntegrateDerivCell3DAuxPlaneLTS: Null pointer!");

  int nnodes = 2;
  int ndof = 3;
  int nrow = LineSlaveElement()->NumNode();
  int ncol = LineMasterElement()->NumNode();

  // slave values
  LINALG::SerialDenseVector sval(nnodes);
  LINALG::SerialDenseMatrix sderiv(nnodes, 1);
  LineSlaveElement()->EvaluateShape(sxi, sval, sderiv, nnodes);

  // master values
  LINALG::SerialDenseVector mval(nnodes);
  LINALG::SerialDenseMatrix mderiv(nnodes, 1);
  LineMasterElement()->EvaluateShape(mxi, mval, mderiv, nnodes);

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;
  typedef std::map<int, double>::const_iterator CI;

  int linsize = 0;
  for (int i = 0; i < nrow; ++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*>(mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  // TODO: this is for safety. Correct linsize
  //       should be predicted
  linsize *= 100;

  //**********************************************************************
  // geometric quantities
  //**********************************************************************
  double gpn[3] = {0.0, 0.0, 0.0};
  GEN::pairedvector<int, double> dgapgp(
      (ncol * ndof) + 10 * linsize);  // gap lin. without lm and jac.
  double gap[1] = {0.0};
  std::vector<GEN::pairedvector<int, double>> dnmap_unit(
      3, 10 * linsize);  // deriv of x,y and z comp. of gpn (unit)

  //**********************************************************************
  // evaluate at GP and lin char. quantities
  //**********************************************************************
  double sgpx[3] = {0.0, 0.0, 0.0};
  double mgpx[3] = {0.0, 0.0, 0.0};

  for (int i = 0; i < nrow; ++i)
  {
    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[i]);
    gpn[0] += sval[i] * mymrtrnode->MoData().n()[0];
    gpn[1] += sval[i] * mymrtrnode->MoData().n()[1];
    gpn[2] += sval[i] * mymrtrnode->MoData().n()[2];

    sgpx[0] += sval[i] * LineSlaveElement()->GetNodalCoords(0, i);
    sgpx[1] += sval[i] * LineSlaveElement()->GetNodalCoords(1, i);
    sgpx[2] += sval[i] * LineSlaveElement()->GetNodalCoords(2, i);
  }

  // build interpolation of master GP coordinates
  for (int i = 0; i < ncol; ++i)
  {
    mgpx[0] += mval[i] * LineMasterElement()->GetNodalCoords(0, i);
    mgpx[1] += mval[i] * LineMasterElement()->GetNodalCoords(1, i);
    mgpx[2] += mval[i] * LineMasterElement()->GetNodalCoords(2, i);
  }

  // normalize interpolated GP normal back to length 1.0 !!!
  double lengthn = sqrt(gpn[0] * gpn[0] + gpn[1] * gpn[1] + gpn[2] * gpn[2]);
  if (lengthn < 1.0e-12) dserror("ERROR: IntegrateAndDerivSegment: Divide by zero!");

  for (int i = 0; i < 3; ++i) gpn[i] /= lengthn;

  // build gap function at current GP
  for (int i = 0; i < Dim(); ++i) gap[0] += (mgpx[i] - sgpx[i]) * gpn[i];

  // build directional derivative of slave GP normal (non-unit)
  GEN::pairedvector<int, double> dmap_nxsl_gp(linsize);
  GEN::pairedvector<int, double> dmap_nysl_gp(linsize);
  GEN::pairedvector<int, double> dmap_nzsl_gp(linsize);

  for (int i = 0; i < nrow; ++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*>(mynodes[i]);

    GEN::pairedvector<int, double>& dmap_nxsl_i = cnode->CoData().GetDerivN()[0];
    GEN::pairedvector<int, double>& dmap_nysl_i = cnode->CoData().GetDerivN()[1];
    GEN::pairedvector<int, double>& dmap_nzsl_i = cnode->CoData().GetDerivN()[2];

    for (_CI p = dmap_nxsl_i.begin(); p != dmap_nxsl_i.end(); ++p)
      dmap_nxsl_gp[p->first] += sval[i] * (p->second);
    for (_CI p = dmap_nysl_i.begin(); p != dmap_nysl_i.end(); ++p)
      dmap_nysl_gp[p->first] += sval[i] * (p->second);
    for (_CI p = dmap_nzsl_i.begin(); p != dmap_nzsl_i.end(); ++p)
      dmap_nzsl_gp[p->first] += sval[i] * (p->second);

    for (_CI p = dsxi.begin(); p != dsxi.end(); ++p)
    {
      double valx = sderiv(i, 0) * cnode->MoData().n()[0];
      dmap_nxsl_gp[p->first] += valx * (p->second);
      double valy = sderiv(i, 0) * cnode->MoData().n()[1];
      dmap_nysl_gp[p->first] += valy * (p->second);
      double valz = sderiv(i, 0) * cnode->MoData().n()[2];
      dmap_nzsl_gp[p->first] += valz * (p->second);
    }
  }

  const double ll = lengthn * lengthn;
  const double linv = 1.0 / (lengthn);
  const double lllinv = 1.0 / (lengthn * lengthn * lengthn);
  const double sxsx = gpn[0] * gpn[0] * ll;
  const double sxsy = gpn[0] * gpn[1] * ll;
  const double sxsz = gpn[0] * gpn[2] * ll;
  const double sysy = gpn[1] * gpn[1] * ll;
  const double sysz = gpn[1] * gpn[2] * ll;
  const double szsz = gpn[2] * gpn[2] * ll;

  for (_CI p = dmap_nxsl_gp.begin(); p != dmap_nxsl_gp.end(); ++p)
  {
    dnmap_unit[0][p->first] += linv * (p->second);
    dnmap_unit[0][p->first] -= lllinv * sxsx * (p->second);
    dnmap_unit[1][p->first] -= lllinv * sxsy * (p->second);
    dnmap_unit[2][p->first] -= lllinv * sxsz * (p->second);
  }

  for (_CI p = dmap_nysl_gp.begin(); p != dmap_nysl_gp.end(); ++p)
  {
    dnmap_unit[1][p->first] += linv * (p->second);
    dnmap_unit[1][p->first] -= lllinv * sysy * (p->second);
    dnmap_unit[0][p->first] -= lllinv * sxsy * (p->second);
    dnmap_unit[2][p->first] -= lllinv * sysz * (p->second);
  }

  for (_CI p = dmap_nzsl_gp.begin(); p != dmap_nzsl_gp.end(); ++p)
  {
    dnmap_unit[2][p->first] += linv * (p->second);
    dnmap_unit[2][p->first] -= lllinv * szsz * (p->second);
    dnmap_unit[0][p->first] -= lllinv * sxsz * (p->second);
    dnmap_unit[1][p->first] -= lllinv * sysz * (p->second);
  }

  // add everything to dgapgp
  for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
    dgapgp[p->first] += (mgpx[0] - sgpx[0]) * (p->second);

  for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
    dgapgp[p->first] += (mgpx[1] - sgpx[1]) * (p->second);

  for (_CI p = dnmap_unit[2].begin(); p != dnmap_unit[2].end(); ++p)
    dgapgp[p->first] += (mgpx[2] - sgpx[2]) * (p->second);

  // lin slave nodes
  for (int z = 0; z < nrow; ++z)
  {
    CoNode* cnode = dynamic_cast<CoNode*>(mynodes[z]);
    for (int k = 0; k < 3; ++k) dgapgp[cnode->Dofs()[k]] -= sval[z] * gpn[k];
  }

  for (_CI p = dsxi.begin(); p != dsxi.end(); ++p)
  {
    double& dg = dgapgp[p->first];
    const double& ps = p->second;
    for (int z = 0; z < nrow; ++z)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(mynodes[z]);
      for (int k = 0; k < 3; ++k) dg -= gpn[k] * sderiv(z, 0) * cnode->xspatial()[k] * ps;
    }
  }

  //        MASTER
  // lin master nodes
  for (int z = 0; z < ncol; ++z)
  {
    CoNode* cnode = dynamic_cast<CoNode*>(mnodes[z]);
    for (int k = 0; k < 3; ++k) dgapgp[cnode->Dofs()[k]] += mval[z] * gpn[k];
  }

  for (_CI p = dmxi.begin(); p != dmxi.end(); ++p)
  {
    double& dg = dgapgp[p->first];
    const double& ps = p->second;
    for (int z = 0; z < ncol; ++z)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(mnodes[z]);
      for (int k = 0; k < 3; ++k) dg += gpn[k] * mderiv(z, 0) * cnode->xspatial()[k] * ps;
    }
  }

  // gap
  CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(mynodes[0]);

  // do not process slave side boundary nodes
  // (their row entries would be zero anyway!)
  if (cnode->IsOnBound()) return;

  if (gap[0] >= 0.0) return;

  double value[3] = {0.0, 0.0, 0.0};
  value[0] = (mgpx[0] - sgpx[0]);  // gap[0]*gpn[0];
  value[1] = (mgpx[1] - sgpx[1]);  // gap[0]*gpn[1];
  value[2] = (mgpx[2] - sgpx[2]);  // gap[0]*gpn[2];

  // add current Gauss point's contribution to gseg
  cnode->AddltlGapValue(value);

  double lengthv = 0.0;
  lengthv = sqrt(value[0] * value[0] + value[1] * value[1] + value[2] * value[2]);
  if (lengthv < 1e-12) dserror("ERROR: zero length!");
  value[0] /= lengthv;
  value[1] /= lengthv;
  value[2] /= lengthv;

  std::vector<std::map<int, double>>& dgmap = cnode->CoData().GetDerivGltl();

  for (_CI p = dgapgp.begin(); p != dgapgp.end(); ++p)
  {
    dgmap[0][p->first] += gpn[0] * (p->second);
    dgmap[1][p->first] += gpn[1] * (p->second);
    dgmap[2][p->first] += gpn[2] * (p->second);
  }

  for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
    dgmap[0][p->first] += gap[0] * (p->second);
  for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
    dgmap[1][p->first] += gap[0] * (p->second);
  for (_CI p = dnmap_unit[2].begin(); p != dnmap_unit[2].end(); ++p)
    dgmap[2][p->first] += gap[0] * (p->second);

  //*****************************************
  // integrate D and M matrix
  // integrate dseg
  for (int k = 0; k < nrow; ++k)
  {
    CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mynodes[k]);

    // multiply the two shape functions
    double prod = sval[k];  // this reduces to sval[k]

    if (abs(prod) > MORTARINTTOL) cnode->AddDltlValue(mnode->Id(), prod);
    //    if(abs(prod)>MORTARINTTOL) cnode->AddSNode(mnode->Id()); // only for friction!
  }

  // integrate mseg
  for (int k = 0; k < ncol; ++k)
  {
    CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mnodes[k]);

    // multiply the two shape functions
    double prod = mval[k];  // this reduces to mval[k]

    if (abs(prod) > MORTARINTTOL) cnode->AddMltlValue(mnode->Id(), prod);
    //    if(abs(prod)>MORTARINTTOL) cnode->AddMNode(mnode->Id());  // only for friction!
  }

  // integrate LinD
  for (int k = 0; k < nrow; ++k)
  {
    // global master node ID
    int mgid = LineSlaveElement()->Nodes()[k]->Id();
    static double fac = 0.0;

    // get the correct map as a reference
    std::map<int, double>& ddmap_jk = cnode->CoData().GetDerivDltl()[mgid];

    // (3) Lin(NMaster) - master GP coordinates
    fac = sderiv(k, 0);
    for (_CI p = dsxi.begin(); p != dsxi.end(); ++p)
    {
      ddmap_jk[p->first] += fac * (p->second);
    }
  }  // loop over slave nodes

  // integrate LinM
  for (int k = 0; k < ncol; ++k)
  {
    // global master node ID
    int mgid = LineMasterElement()->Nodes()[k]->Id();
    static double fac = 0.0;

    // get the correct map as a reference
    std::map<int, double>& dmmap_jk = cnode->CoData().GetDerivMltl()[mgid];

    // (1) Lin(Phi) - dual shape functions
    // this vanishes here since there are no deformation-dependent dual functions

    // (3) Lin(NMaster) - master GP coordinates
    fac = mderiv(k, 0);
    for (_CI p = dmxi.begin(); p != dmxi.end(); ++p)
    {
      dmmap_jk[p->first] += fac * (p->second);
    }
  }  // loop over master nodes


  //  std::cout << "element is evaluated!" << std::endl;

  //***************************************************************************
  if (friction)
  {
    // tangent:
    // first jump:
    double jump[3] = {0.0, 0.0, 0.0};
    double sgpxold[3] = {0.0, 0.0, 0.0};
    double mgpxold[3] = {0.0, 0.0, 0.0};

    int oldID = -1;

    // loop over all slave nodes
    for (int i = 0; i < idiscret_.NodeColMap()->NumMyElements(); ++i)
    {
      int gid1 = idiscret_.NodeColMap()->GID(i);
      DRT::Node* node1 = idiscret_.gNode(gid1);
      if (!node1) dserror("ERROR: Cannot find node with gid %", gid1);
      CoNode* contactnode = dynamic_cast<CoNode*>(node1);

      // here only slave nodes
      if (!contactnode->IsSlave()) continue;

      // check if dold is present
      if (dynamic_cast<FriNode*>(contactnode)->FriData().GetDOldLTL().size() < 1) continue;

      // store id
      oldID = gid1;

      // if we are here, break
      break;
    }

    // linearizations
    GEN::pairedvector<int, double> sgpxoldlinx(linsize);
    GEN::pairedvector<int, double> sgpxoldliny(linsize);
    GEN::pairedvector<int, double> sgpxoldlinz(linsize);

    GEN::pairedvector<int, double> mgpxoldlinx(linsize);
    GEN::pairedvector<int, double> mgpxoldliny(linsize);
    GEN::pairedvector<int, double> mgpxoldlinz(linsize);

    if (oldID > -1)
    {
      DRT::Node* node1 = idiscret_.gNode(oldID);
      if (!node1) dserror("ERROR: Cannot find node with gid %", oldID);
      CoNode* contactnode = dynamic_cast<CoNode*>(node1);

      // check if we have dold
      if (dynamic_cast<FriNode*>(contactnode)->FriData().GetDOldLTL().size() > 0)
      {
        for (_CI p = dynamic_cast<FriNode*>(contactnode)->FriData().GetDOldLTL().begin();
             p != dynamic_cast<FriNode*>(contactnode)->FriData().GetDOldLTL().end(); ++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_.gNode(gid3);
          if (!snode) dserror("ERROR: Cannot find node with gid");
          CoNode* csnode = dynamic_cast<CoNode*>(snode);

          for (int dim = 0; dim < Dim(); ++dim)
          {
            sgpxold[dim] += p->second * csnode->xspatial()[dim];
          }
          sgpxoldlinx[csnode->Dofs()[0]] += p->second;
          sgpxoldliny[csnode->Dofs()[1]] += p->second;
          sgpxoldlinz[csnode->Dofs()[2]] += p->second;
        }

        // safety
        if (dynamic_cast<FriNode*>(contactnode)->FriData().GetMOldLTL().size() < 1)
          dserror("ERROR: something went wrong!");

        for (CI p = dynamic_cast<FriNode*>(contactnode)->FriData().GetMOldLTL().begin();
             p != dynamic_cast<FriNode*>(contactnode)->FriData().GetMOldLTL().end(); ++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* mnode = idiscret_.gNode(gid3);
          if (!mnode) dserror("ERROR: Cannot find node with gid");
          CoNode* cmnode = dynamic_cast<CoNode*>(mnode);

          for (int dim = 0; dim < Dim(); ++dim)
          {
            mgpxold[dim] += p->second * cmnode->xspatial()[dim];
          }
          mgpxoldlinx[cmnode->Dofs()[0]] += p->second;
          mgpxoldliny[cmnode->Dofs()[1]] += p->second;
          mgpxoldlinz[cmnode->Dofs()[2]] += p->second;
        }
      }
    }

    // create slip
    for (int dim = 0; dim < Dim(); ++dim)
      jump[dim] = mgpx[dim] - mgpxold[dim] - (sgpx[dim] - sgpxold[dim]);

    Epetra_SerialDenseMatrix tanplane(3, 3);
    tanplane(0, 0) = 1 - (value[0] * value[0]);
    tanplane(0, 1) = -(value[0] * value[1]);
    tanplane(0, 2) = -(value[0] * value[2]);
    tanplane(1, 0) = -(value[1] * value[0]);
    tanplane(1, 1) = 1 - (value[1] * value[1]);
    tanplane(1, 2) = -(value[1] * value[2]);
    tanplane(2, 0) = -(value[2] * value[0]);
    tanplane(2, 1) = -(value[2] * value[1]);
    tanplane(2, 2) = 1 - (value[2] * value[2]);

    double finaljump[3] = {0.0, 0.0, 0.0};
    finaljump[0] = tanplane(0, 0) * jump[0] + tanplane(0, 1) * jump[1] + tanplane(0, 2) * jump[2];
    finaljump[1] = tanplane(1, 0) * jump[0] + tanplane(1, 1) * jump[1] + tanplane(1, 2) * jump[2];
    finaljump[2] = tanplane(2, 0) * jump[0] + tanplane(2, 1) * jump[1] + tanplane(2, 2) * jump[2];

    cnode->AddltlJumpValue(finaljump);
    //    std::cout << "jump = " << sqrt(finaljump[0]*finaljump[0] + finaljump[1]*finaljump[1]
    //    +finaljump[2]*finaljump[2]) << std::endl;

    std::vector<std::map<int, double>>& djmapfinal = cnode->CoData().GetDerivJumpltl();

    std::vector<GEN::pairedvector<int, double>> djmap(3, 100);

    // lin slave nodes
    for (int z = 0; z < nrow; ++z)
    {
      CoNode* node = dynamic_cast<CoNode*>(mynodes[z]);
      for (int k = 0; k < 3; ++k) djmap[k][node->Dofs()[k]] -= sval[z];
    }


    for (int k = 0; k < nrow; ++k)
    {
      CoNode* node = dynamic_cast<CoNode*>(mynodes[k]);
      for (_CI p = dsxi.begin(); p != dsxi.end(); ++p)
      {
        for (int z = 0; z < 3; ++z)
          djmap[z][p->first] -= sderiv(k, 0) * (p->second) * node->xspatial()[z];
      }
    }  // loop over slave nodes


    // lin master nodes
    for (int z = 0; z < ncol; ++z)
    {
      CoNode* node = dynamic_cast<CoNode*>(mnodes[z]);
      for (int k = 0; k < 3; ++k) djmap[k][node->Dofs()[k]] += mval[z];
    }

    for (int k = 0; k < ncol; ++k)
    {
      CoNode* node = dynamic_cast<CoNode*>(mnodes[k]);
      for (_CI p = dmxi.begin(); p != dmxi.end(); ++p)
      {
        for (int z = 0; z < 3; ++z)
          djmap[z][p->first] += mderiv(k, 0) * (p->second) * node->xspatial()[z];
      }
    }  // loop over slave nodes

    // sgpx and mgpx old
    for (_CI p = mgpxoldlinx.begin(); p != mgpxoldlinx.end(); ++p) djmap[0][p->first] -= p->second;
    for (_CI p = mgpxoldliny.begin(); p != mgpxoldliny.end(); ++p) djmap[1][p->first] -= p->second;
    for (_CI p = mgpxoldlinz.begin(); p != mgpxoldlinz.end(); ++p) djmap[2][p->first] -= p->second;

    for (_CI p = sgpxoldlinx.begin(); p != sgpxoldlinx.end(); ++p) djmap[0][p->first] += p->second;
    for (_CI p = sgpxoldliny.begin(); p != sgpxoldliny.end(); ++p) djmap[1][p->first] += p->second;
    for (_CI p = sgpxoldlinz.begin(); p != sgpxoldlinz.end(); ++p) djmap[2][p->first] += p->second;



    std::vector<GEN::pairedvector<int, double>> tplanex(3, 100);
    std::vector<GEN::pairedvector<int, double>> tplaney(3, 100);
    std::vector<GEN::pairedvector<int, double>> tplanez(3, 100);

    for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
      tplanex[0][p->first] -= gpn[0] * p->second;
    for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
      tplanex[1][p->first] -= gpn[1] * p->second;
    for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
      tplanex[2][p->first] -= gpn[2] * p->second;

    for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
      tplaney[0][p->first] -= gpn[0] * p->second;
    for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
      tplaney[1][p->first] -= gpn[1] * p->second;
    for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
      tplaney[2][p->first] -= gpn[2] * p->second;

    for (_CI p = dnmap_unit[2].begin(); p != dnmap_unit[2].end(); ++p)
      tplanez[0][p->first] -= gpn[0] * p->second;
    for (_CI p = dnmap_unit[2].begin(); p != dnmap_unit[2].end(); ++p)
      tplanez[1][p->first] -= gpn[1] * p->second;
    for (_CI p = dnmap_unit[2].begin(); p != dnmap_unit[2].end(); ++p)
      tplanez[2][p->first] -= gpn[2] * p->second;

    //------------
    for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
      tplanex[0][p->first] -= gpn[0] * p->second;
    for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
      tplanex[1][p->first] -= gpn[0] * p->second;
    for (_CI p = dnmap_unit[2].begin(); p != dnmap_unit[2].end(); ++p)
      tplanex[2][p->first] -= gpn[0] * p->second;

    for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
      tplaney[0][p->first] -= gpn[1] * p->second;
    for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
      tplaney[1][p->first] -= gpn[1] * p->second;
    for (_CI p = dnmap_unit[2].begin(); p != dnmap_unit[2].end(); ++p)
      tplaney[2][p->first] -= gpn[1] * p->second;

    for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
      tplanez[0][p->first] -= gpn[2] * p->second;
    for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
      tplanez[1][p->first] -= gpn[2] * p->second;
    for (_CI p = dnmap_unit[2].begin(); p != dnmap_unit[2].end(); ++p)
      tplanez[2][p->first] -= gpn[2] * p->second;

    //-----------------------------

    for (_CI p = djmap[0].begin(); p != djmap[0].end(); ++p)
      djmapfinal[0][p->first] += tanplane(0, 0) * p->second;
    for (_CI p = djmap[1].begin(); p != djmap[1].end(); ++p)
      djmapfinal[0][p->first] += tanplane(0, 1) * p->second;
    for (_CI p = djmap[2].begin(); p != djmap[2].end(); ++p)
      djmapfinal[0][p->first] += tanplane(0, 2) * p->second;

    for (_CI p = djmap[0].begin(); p != djmap[0].end(); ++p)
      djmapfinal[1][p->first] += tanplane(1, 0) * p->second;
    for (_CI p = djmap[1].begin(); p != djmap[1].end(); ++p)
      djmapfinal[1][p->first] += tanplane(1, 1) * p->second;
    for (_CI p = djmap[2].begin(); p != djmap[2].end(); ++p)
      djmapfinal[1][p->first] += tanplane(1, 2) * p->second;

    for (_CI p = djmap[0].begin(); p != djmap[0].end(); ++p)
      djmapfinal[2][p->first] += tanplane(2, 0) * p->second;
    for (_CI p = djmap[1].begin(); p != djmap[1].end(); ++p)
      djmapfinal[2][p->first] += tanplane(2, 1) * p->second;
    for (_CI p = djmap[2].begin(); p != djmap[2].end(); ++p)
      djmapfinal[2][p->first] += tanplane(2, 2) * p->second;

    //-----------------------------
    for (_CI p = tplanex[0].begin(); p != tplanex[0].end(); ++p)
      djmapfinal[0][p->first] += jump[0] * p->second;
    for (_CI p = tplanex[1].begin(); p != tplanex[1].end(); ++p)
      djmapfinal[0][p->first] += jump[1] * p->second;
    for (_CI p = tplanex[2].begin(); p != tplanex[2].end(); ++p)
      djmapfinal[0][p->first] += jump[2] * p->second;

    for (_CI p = tplaney[0].begin(); p != tplaney[0].end(); ++p)
      djmapfinal[1][p->first] += jump[0] * p->second;
    for (_CI p = tplaney[1].begin(); p != tplaney[1].end(); ++p)
      djmapfinal[1][p->first] += jump[1] * p->second;
    for (_CI p = tplaney[2].begin(); p != tplaney[2].end(); ++p)
      djmapfinal[1][p->first] += jump[2] * p->second;

    for (_CI p = tplanez[0].begin(); p != tplanez[0].end(); ++p)
      djmapfinal[2][p->first] += jump[0] * p->second;
    for (_CI p = tplanez[1].begin(); p != tplanez[1].end(); ++p)
      djmapfinal[2][p->first] += jump[1] * p->second;
    for (_CI p = tplanez[2].begin(); p != tplanez[2].end(); ++p)
      djmapfinal[2][p->first] += jump[2] * p->second;

  }  // end friction

  return;
}

/*----------------------------------------------------------------------*
 |  line-line intersection                                   farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineToLineCouplingPoint3d::LineIntersection(double* sxi, double* mxi,
    GEN::pairedvector<int, double>& dsxi, GEN::pairedvector<int, double>& dmxi)
{
  // flag for debug output
  const bool out = false;

  // only for line 2
  const int nnodes = 2;

  // prepare linearizations
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // calculate slave vector
  CoNode* ns1 = dynamic_cast<CoNode*>(LineSlaveElement()->Nodes()[0]);
  CoNode* ns2 = dynamic_cast<CoNode*>(LineSlaveElement()->Nodes()[1]);
  ns1->BuildAveragedEdgeTangent();
  ns2->BuildAveragedEdgeTangent();

  // calculate slave vector
  CoNode* nm1 = dynamic_cast<CoNode*>(LineMasterElement()->Nodes()[0]);
  CoNode* nm2 = dynamic_cast<CoNode*>(LineMasterElement()->Nodes()[1]);
  nm1->BuildAveragedEdgeTangent();
  nm2->BuildAveragedEdgeTangent();

  double lengths1 = sqrt(ns1->MoData().EdgeTangent()[0] * ns1->MoData().EdgeTangent()[0] +
                         ns1->MoData().EdgeTangent()[1] * ns1->MoData().EdgeTangent()[1] +
                         ns1->MoData().EdgeTangent()[2] * ns1->MoData().EdgeTangent()[2]);
  double lengths2 = sqrt(ns2->MoData().EdgeTangent()[0] * ns2->MoData().EdgeTangent()[0] +
                         ns2->MoData().EdgeTangent()[1] * ns2->MoData().EdgeTangent()[1] +
                         ns2->MoData().EdgeTangent()[2] * ns2->MoData().EdgeTangent()[2]);
  double lengthm1 = sqrt(nm1->MoData().EdgeTangent()[0] * nm1->MoData().EdgeTangent()[0] +
                         nm1->MoData().EdgeTangent()[1] * nm1->MoData().EdgeTangent()[1] +
                         nm1->MoData().EdgeTangent()[2] * nm1->MoData().EdgeTangent()[2]);
  double lengthm2 = sqrt(nm2->MoData().EdgeTangent()[0] * nm2->MoData().EdgeTangent()[0] +
                         nm2->MoData().EdgeTangent()[1] * nm2->MoData().EdgeTangent()[1] +
                         nm2->MoData().EdgeTangent()[2] * nm2->MoData().EdgeTangent()[2]);
  if (lengths1 < 1e-12 or lengths2 < 1e-12 or lengthm1 < 1e-12 or lengthm2 < 1e-12)
    dserror("ERROR: tangents zero length");

  // calc angle between tangents
  double ts1[3] = {0.0, 0.0, 0.0};
  double ts2[3] = {0.0, 0.0, 0.0};
  double tm1[3] = {0.0, 0.0, 0.0};
  double tm2[3] = {0.0, 0.0, 0.0};

  ts1[0] = ns1->MoData().EdgeTangent()[0];
  ts1[1] = ns1->MoData().EdgeTangent()[1];
  ts1[2] = ns1->MoData().EdgeTangent()[2];

  ts2[0] = ns2->MoData().EdgeTangent()[0];
  ts2[1] = ns2->MoData().EdgeTangent()[1];
  ts2[2] = ns2->MoData().EdgeTangent()[2];

  tm1[0] = nm1->MoData().EdgeTangent()[0];
  tm1[1] = nm1->MoData().EdgeTangent()[1];
  tm1[2] = nm1->MoData().EdgeTangent()[2];

  tm2[0] = nm2->MoData().EdgeTangent()[0];
  tm2[1] = nm2->MoData().EdgeTangent()[1];
  tm2[2] = nm2->MoData().EdgeTangent()[2];

  if (out)
  {
    std::cout << "slave 1 = " << ts1[0] << "  " << ts1[1] << "  " << ts1[2] << std::endl;
    std::cout << "slave 2 = " << ts2[0] << "  " << ts2[1] << "  " << ts2[2] << std::endl;
  }

  double test = ts1[0] * ts2[0] + ts1[1] * ts2[1] + ts1[2] * ts2[2];
  if (test < 1e-8)
  {
    ns2->MoData().EdgeTangent()[0] *= -1.0;
    ns2->MoData().EdgeTangent()[1] *= -1.0;
    ns2->MoData().EdgeTangent()[2] *= -1.0;

    ts2[0] *= -1.0;
    ts2[1] *= -1.0;
    ts2[2] *= -1.0;

    for (_CI p = ns2->CoData().GetDerivTangent()[0].begin();
         p != ns2->CoData().GetDerivTangent()[0].end(); ++p)
      ns2->CoData().GetDerivTangent()[0][p->first] *= -1.0;
    for (_CI p = ns2->CoData().GetDerivTangent()[1].begin();
         p != ns2->CoData().GetDerivTangent()[1].end(); ++p)
      ns2->CoData().GetDerivTangent()[1][p->first] *= -1.0;
    for (_CI p = ns2->CoData().GetDerivTangent()[2].begin();
         p != ns2->CoData().GetDerivTangent()[2].end(); ++p)
      ns2->CoData().GetDerivTangent()[2][p->first] *= -1.0;
  }
  if (out)
  {
    std::cout << "----------------" << std::endl;
    std::cout << "slave 1 = " << ts1[0] << "  " << ts1[1] << "  " << ts1[2] << std::endl;
    std::cout << "slave 2 = " << ts2[0] << "  " << ts2[1] << "  " << ts2[2] << std::endl;

    std::cout << "master 1 = " << tm1[0] << "  " << tm1[1] << "  " << tm1[2] << std::endl;
    std::cout << "master 2 = " << tm2[0] << "  " << tm2[1] << "  " << tm2[2] << std::endl;
  }

  test = tm1[0] * tm2[0] + tm1[1] * tm2[1] + tm1[2] * tm2[2];
  if (test < 1e-8)
  {
    nm2->MoData().EdgeTangent()[0] *= -1.0;
    nm2->MoData().EdgeTangent()[1] *= -1.0;
    nm2->MoData().EdgeTangent()[2] *= -1.0;

    tm2[0] *= -1.0;
    tm2[1] *= -1.0;
    tm2[2] *= -1.0;

    for (_CI p = nm2->CoData().GetDerivTangent()[0].begin();
         p != nm2->CoData().GetDerivTangent()[0].end(); ++p)
      nm2->CoData().GetDerivTangent()[0][p->first] *= -1.0;
    for (_CI p = nm2->CoData().GetDerivTangent()[1].begin();
         p != nm2->CoData().GetDerivTangent()[1].end(); ++p)
      nm2->CoData().GetDerivTangent()[1][p->first] *= -1.0;
    for (_CI p = nm2->CoData().GetDerivTangent()[2].begin();
         p != nm2->CoData().GetDerivTangent()[2].end(); ++p)
      nm2->CoData().GetDerivTangent()[2][p->first] *= -1.0;
  }
  if (out)
  {
    std::cout << "----------------" << std::endl;
    std::cout << "master 1 = " << tm1[0] << "  " << tm1[1] << "  " << tm1[2] << std::endl;
    std::cout << "master 2 = " << tm2[0] << "  " << tm2[1] << "  " << tm2[2] << std::endl;
  }

  // res norm
  double conv = 0.0;

  // start in the element center
  double xiS[1] = {0.0};  // xi_slave
  double xiM[1] = {0.0};  // xi_master

  // function f (vector-valued)
  double f[2] = {0.0, 0.0};

  // gradient of f (df/deta[0], df/deta[1])
  LINALG::Matrix<2, 2> df;

  // Newton
  for (int k = 0; k < MORTARMAXITER; ++k)
  {
    //**********************************************
    //  F CALCULATION                             //
    //**********************************************
    // slave values
    LINALG::SerialDenseVector sval(nnodes);
    LINALG::SerialDenseMatrix sderiv(nnodes, 1);
    LineSlaveElement()->EvaluateShape(xiS, sval, sderiv, nnodes);

    // master values
    LINALG::SerialDenseVector mval(nnodes);
    LINALG::SerialDenseMatrix mderiv(nnodes, 1);
    LineMasterElement()->EvaluateShape(xiM, mval, mderiv, nnodes);

    double xs[3] = {0.0, 0.0, 0.0};
    double xm[3] = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; ++i)
      xs[i] += sval[0] * ns1->xspatial()[i] + sval[1] * ns2->xspatial()[i];
    for (int i = 0; i < 3; ++i)
      xm[i] += mval[0] * nm1->xspatial()[i] + mval[1] * nm2->xspatial()[i];

    double xdiff[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) xdiff[i] = xs[i] - xm[i];

    // calculate tangents:
    double vs[3] = {0.0, 0.0, 0.0};
    double vm[3] = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; ++i) vs[i] += sval[0] * ts1[i] + sval[1] * ts2[i];
    for (int i = 0; i < 3; ++i) vm[i] += mval[0] * tm1[i] + mval[1] * tm2[i];

    f[0] = xdiff[0] * vs[0] + xdiff[1] * vs[1] + xdiff[2] * vs[2];
    f[1] = xdiff[0] * vm[0] + xdiff[1] * vm[1] + xdiff[2] * vm[2];

    // check for convergence
    conv = sqrt(f[0] * f[0] + f[1] * f[1]);
    if (conv <= MORTARCONVTOL) break;

    //**********************************************
    //   F GRADIENT CALCULATION                   //
    //**********************************************

    double xsderiv[3] = {0.0, 0.0, 0.0};
    double xmderiv[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i)
      xsderiv[i] += sderiv(0, 0) * ns1->xspatial()[i] + sderiv(1, 0) * ns2->xspatial()[i];
    for (int i = 0; i < 3; ++i)
      xmderiv[i] += mderiv(0, 0) * nm1->xspatial()[i] + mderiv(1, 0) * nm2->xspatial()[i];

    double vsderiv[3] = {0.0, 0.0, 0.0};
    double vmderiv[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) vsderiv[i] += sderiv(0, 0) * ts1[i] + sderiv(1, 0) * ts2[i];
    for (int i = 0; i < 3; ++i) vmderiv[i] += mderiv(0, 0) * tm1[i] + mderiv(1, 0) * tm2[i];

    df(0, 0) = xsderiv[0] * vs[0] + xsderiv[1] * vs[1] + xsderiv[2] * vs[2] +
               vsderiv[0] * xdiff[0] + vsderiv[1] * xdiff[1] + vsderiv[2] * xdiff[2];

    df(0, 1) = -xmderiv[0] * vs[0] - xmderiv[1] * vs[1] - xmderiv[2] * vs[2];

    df(1, 0) = xsderiv[0] * vm[0] + xsderiv[1] * vm[1] + xsderiv[2] * vm[2];

    df(1, 1) = -xmderiv[0] * vm[0] - xmderiv[1] * vm[1] - xmderiv[2] * vm[2] +
               vmderiv[0] * xdiff[0] + vmderiv[1] * xdiff[1] + vmderiv[2] * xdiff[2];

    //**********************************************
    //   solve deta = - inv(dF) * F               //
    //**********************************************
    double jacdet = df.Invert();
    if (abs(jacdet) < 1.0e-12)
    {
      sxi[0] = 1e12;
      mxi[0] = 1e12;
      return;
      dserror("ERROR: Singular Jacobian for projection");
    }

    // update eta and alpha
    xiS[0] += -df(0, 0) * f[0] - df(0, 1) * f[1];
    xiM[0] += -df(1, 0) * f[0] - df(1, 1) * f[1];
  }

  // Newton iteration unconverged
  if (conv > MORTARCONVTOL) dserror("ERROR: LTL intersection not converged!");

  //**********************************************
  //  Linearization                             //
  //**********************************************
  // slave values
  LINALG::SerialDenseVector sval(nnodes);
  LINALG::SerialDenseMatrix sderiv(nnodes, 1);
  LineSlaveElement()->EvaluateShape(xiS, sval, sderiv, nnodes);

  // master values
  LINALG::SerialDenseVector mval(nnodes);
  LINALG::SerialDenseMatrix mderiv(nnodes, 1);
  LineMasterElement()->EvaluateShape(xiM, mval, mderiv, nnodes);

  double xs[3] = {0.0, 0.0, 0.0};
  double xm[3] = {0.0, 0.0, 0.0};

  for (int i = 0; i < 3; ++i) xs[i] += sval[0] * ns1->xspatial()[i] + sval[1] * ns2->xspatial()[i];
  for (int i = 0; i < 3; ++i) xm[i] += mval[0] * nm1->xspatial()[i] + mval[1] * nm2->xspatial()[i];

  double xdiff[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i) xdiff[i] = xs[i] - xm[i];

  // calculate tangents:
  double vs[3] = {0.0, 0.0, 0.0};
  double vm[3] = {0.0, 0.0, 0.0};

  for (int i = 0; i < 3; ++i) vs[i] += sval[0] * ts1[i] + sval[1] * ts2[i];
  for (int i = 0; i < 3; ++i) vm[i] += mval[0] * tm1[i] + mval[1] * tm2[i];

  std::vector<GEN::pairedvector<int, double>> xLin(3, 1000);
  std::vector<GEN::pairedvector<int, double>> vsLin(3, 1000);
  std::vector<GEN::pairedvector<int, double>> vmLin(3, 1000);

  // global position difference
  for (int i = 0; i < 3; ++i) (xLin[i])[ns1->Dofs()[i]] += sval(0);
  for (int i = 0; i < 3; ++i) (xLin[i])[ns2->Dofs()[i]] += sval(1);

  for (int i = 0; i < 3; ++i) (xLin[i])[nm1->Dofs()[i]] -= mval(0);
  for (int i = 0; i < 3; ++i) (xLin[i])[nm2->Dofs()[i]] -= mval(1);

  // TODO: this would be the correct linearization! however, the old one works better. no idea
  // why!?!?!? tangent vector slave
  for (int i = 0; i < 3; ++i)
  {
    for (_CI p = ns1->CoData().GetDerivTangent()[i].begin();
         p != ns1->CoData().GetDerivTangent()[i].end(); ++p)
      (vsLin[i])[p->first] += sval[0] * p->second;

    for (_CI p = ns2->CoData().GetDerivTangent()[i].begin();
         p != ns2->CoData().GetDerivTangent()[i].end(); ++p)
      (vsLin[i])[p->first] += sval[1] * p->second;
  }

  // tangent vector master
  for (int i = 0; i < 3; ++i)
  {
    for (_CI p = nm1->CoData().GetDerivTangent()[i].begin();
         p != nm1->CoData().GetDerivTangent()[i].end(); ++p)
      (vmLin[i])[p->first] += mval[0] * p->second;

    for (_CI p = nm2->CoData().GetDerivTangent()[i].begin();
         p != nm2->CoData().GetDerivTangent()[i].end(); ++p)
      (vmLin[i])[p->first] += mval[1] * p->second;
  }

  // TODO: this is the old linearization:
  // tangent vector slave
  //  for(int i=0; i<3;++i)
  //    (vsLin[i])[ns1->Dofs()[i]] += 1;
  //  for(int i=0; i<3;++i)
  //    (vsLin[i])[ns2->Dofs()[i]] -= 1;
  //
  //  // tangent vector master
  //  for(int i=0; i<3;++i)
  //    (vmLin[i])[nm1->Dofs()[i]] += 1;
  //  for(int i=0; i<3;++i)
  //    (vmLin[i])[nm2->Dofs()[i]] -= 1;

  GEN::pairedvector<int, double> f0(1000);
  GEN::pairedvector<int, double> f1(1000);

  // lin xdiff * tangent + xdiff * lin tangent
  for (_CI p = xLin[0].begin(); p != xLin[0].end(); ++p) f0[p->first] += (p->second) * vs[0];
  for (_CI p = xLin[1].begin(); p != xLin[1].end(); ++p) f0[p->first] += (p->second) * vs[1];
  for (_CI p = xLin[2].begin(); p != xLin[2].end(); ++p) f0[p->first] += (p->second) * vs[2];

  for (_CI p = vsLin[0].begin(); p != vsLin[0].end(); ++p) f0[p->first] += (p->second) * xdiff[0];
  for (_CI p = vsLin[1].begin(); p != vsLin[1].end(); ++p) f0[p->first] += (p->second) * xdiff[1];
  for (_CI p = vsLin[2].begin(); p != vsLin[2].end(); ++p) f0[p->first] += (p->second) * xdiff[2];

  // lin xdiff * tangent + xdiff * lin tangent
  for (_CI p = xLin[0].begin(); p != xLin[0].end(); ++p) f1[p->first] += (p->second) * vm[0];
  for (_CI p = xLin[1].begin(); p != xLin[1].end(); ++p) f1[p->first] += (p->second) * vm[1];
  for (_CI p = xLin[2].begin(); p != xLin[2].end(); ++p) f1[p->first] += (p->second) * vm[2];

  for (_CI p = vmLin[0].begin(); p != vmLin[0].end(); ++p) f1[p->first] += (p->second) * xdiff[0];
  for (_CI p = vmLin[1].begin(); p != vmLin[1].end(); ++p) f1[p->first] += (p->second) * xdiff[1];
  for (_CI p = vmLin[2].begin(); p != vmLin[2].end(); ++p) f1[p->first] += (p->second) * xdiff[2];

  // end
  for (_CI p = f0.begin(); p != f0.end(); ++p) dsxi[p->first] -= (p->second) * df(0, 0);
  for (_CI p = f1.begin(); p != f1.end(); ++p) dsxi[p->first] -= (p->second) * df(0, 1);

  for (_CI p = f0.begin(); p != f0.end(); ++p) dmxi[p->first] -= (p->second) * df(1, 0);
  for (_CI p = f1.begin(); p != f1.end(); ++p) dmxi[p->first] -= (p->second) * df(1, 1);

  sxi[0] = xiS[0];
  mxi[0] = xiM[0];

  return;
}

/*----------------------------------------------------------------------*
 |  check if intersection is in para space interval          farah 07/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineToLineCouplingPoint3d::CheckIntersection(double* sxi, double* mxi)
{
  if (sxi[0] >= -1.0 - 1e-12 and sxi[0] <= 1.0 + 1e-12 and mxi[0] >= -1.0 - 1e-12 and
      mxi[0] <= 1.0 + 1e-12)
    return true;

  return false;
}

/*----------------------------------------------------------------------*
 |  check if line eles parallel                              farah 07/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineToLineCouplingPoint3d::CheckParallelity()
{
  // tolerance for line clipping
  const double sminedge = LineSlaveElement()->MinEdgeSize();
  const double mminedge = LineMasterElement()->MinEdgeSize();
  const double tol = MORTARCLIPTOL * std::min(sminedge, mminedge);

  double vs[3] = {0.0, 0.0, 0.0};
  double vm[3] = {0.0, 0.0, 0.0};

  // calculate slave vector
  CoNode* ns1 = dynamic_cast<CoNode*>(LineSlaveElement()->Nodes()[0]);
  CoNode* ns2 = dynamic_cast<CoNode*>(LineSlaveElement()->Nodes()[1]);

  vs[0] = ns1->xspatial()[0] - ns2->xspatial()[0];
  vs[1] = ns1->xspatial()[1] - ns2->xspatial()[1];
  vs[2] = ns1->xspatial()[2] - ns2->xspatial()[2];

  // calculate slave vector
  CoNode* nm1 = dynamic_cast<CoNode*>(LineMasterElement()->Nodes()[0]);
  CoNode* nm2 = dynamic_cast<CoNode*>(LineMasterElement()->Nodes()[1]);

  vm[0] = nm1->xspatial()[0] - nm2->xspatial()[0];
  vm[1] = nm1->xspatial()[1] - nm2->xspatial()[1];
  vm[2] = nm1->xspatial()[2] - nm2->xspatial()[2];

  // calculate lengths
  const double lengthS = sqrt(vs[0] * vs[0] + vs[1] * vs[1] + vs[2] * vs[2]);
  const double lengthM = sqrt(vm[0] * vm[0] + vm[1] * vm[1] + vm[2] * vm[2]);

  // calculate scalar product
  const double scaprod = vs[0] * vm[0] + vs[1] * vm[1] + vs[2] * vm[2];

  // proof if scalar product equals length product --> parallelity
  const double diff = abs(scaprod) - (lengthS * lengthM);

  if (abs(diff) < tol) return true;

  return false;
}



/*----------------------------------------------------------------------*
 |  calc current angle (rad) between edges                   farah 09/16|
 *----------------------------------------------------------------------*/
double CONTACT::LineToLineCouplingPoint3d::CalcCurrentAngle(
    GEN::pairedvector<int, double>& lineAngle)
{
  // define iterator for linerization
  typedef GEN::pairedvector<int, double>::const_iterator CI;

  // slave edge vector and master vector edge
  double vs[3] = {0.0, 0.0, 0.0};
  double vm[3] = {0.0, 0.0, 0.0};

  // calculate slave vector
  CoNode* ns1 = dynamic_cast<CoNode*>(LineSlaveElement()->Nodes()[0]);
  CoNode* ns2 = dynamic_cast<CoNode*>(LineSlaveElement()->Nodes()[1]);

  vs[0] = ns1->xspatial()[0] - ns2->xspatial()[0];
  vs[1] = ns1->xspatial()[1] - ns2->xspatial()[1];
  vs[2] = ns1->xspatial()[2] - ns2->xspatial()[2];

  // calculate slave vector
  CoNode* nm1 = dynamic_cast<CoNode*>(LineMasterElement()->Nodes()[0]);
  CoNode* nm2 = dynamic_cast<CoNode*>(LineMasterElement()->Nodes()[1]);

  vm[0] = nm1->xspatial()[0] - nm2->xspatial()[0];
  vm[1] = nm1->xspatial()[1] - nm2->xspatial()[1];
  vm[2] = nm1->xspatial()[2] - nm2->xspatial()[2];

  // calculate lengths
  const double lengthS = sqrt(vs[0] * vs[0] + vs[1] * vs[1] + vs[2] * vs[2]);
  const double lengthM = sqrt(vm[0] * vm[0] + vm[1] * vm[1] + vm[2] * vm[2]);

  // safety
  if (lengthS < 1e-12 or lengthM < 1e-12) dserror("ERROR: line elements of zero length!");

  // calculate scalar product
  double scaprod = vs[0] * vm[0] + vs[1] * vm[1] + vs[2] * vm[2];
  double scaledScaprod = scaprod / (lengthS * lengthM);
  double angleRad = acos(scaledScaprod);

  // check if we used the right angle
  bool switchSign = false;
  if (angleRad > 0.5 * M_PI)  // if angle is > 90°
  {
    switchSign = true;

    // change sign of master vector
    vm[0] = -vm[0];
    vm[1] = -vm[1];
    vm[2] = -vm[2];

    scaprod = vs[0] * vm[0] + vs[1] * vm[1] + vs[2] * vm[2];
    scaledScaprod = scaprod / (lengthS * lengthM);
    angleRad = acos(scaledScaprod);
  }

  //===============================================================
  // linearization

  // delta lengthM
  std::vector<GEN::pairedvector<int, double>> DlM(3, 1000);
  GEN::pairedvector<int, double> DlengthM(1000);

  // change sign of master vectors linearization
  if (switchSign)
  {
    DlM[0][nm1->Dofs()[0]] -= 1;
    DlM[0][nm2->Dofs()[0]] += 1;
    DlM[1][nm1->Dofs()[1]] -= 1;
    DlM[1][nm2->Dofs()[1]] += 1;
    DlM[2][nm1->Dofs()[2]] -= 1;
    DlM[2][nm2->Dofs()[2]] += 1;
  }
  else
  {
    DlM[0][nm1->Dofs()[0]] += 1;
    DlM[0][nm2->Dofs()[0]] -= 1;
    DlM[1][nm1->Dofs()[1]] += 1;
    DlM[1][nm2->Dofs()[1]] -= 1;
    DlM[2][nm1->Dofs()[2]] += 1;
    DlM[2][nm2->Dofs()[2]] -= 1;
  }


  for (CI p = DlM[0].begin(); p != DlM[0].end(); ++p)
    (DlengthM)[p->first] += (p->second) * vm[0] * 1.0 / (lengthM);
  for (CI p = DlM[1].begin(); p != DlM[1].end(); ++p)
    (DlengthM)[p->first] += (p->second) * vm[1] * 1.0 / (lengthM);
  for (CI p = DlM[2].begin(); p != DlM[2].end(); ++p)
    (DlengthM)[p->first] += (p->second) * vm[2] * 1.0 / (lengthM);

  // delta lengthS
  std::vector<GEN::pairedvector<int, double>> DlS(3, 1000);
  GEN::pairedvector<int, double> DlengthS(1000);

  DlS[0][ns1->Dofs()[0]] += 1;
  DlS[0][ns2->Dofs()[0]] -= 1;
  DlS[1][ns1->Dofs()[1]] += 1;
  DlS[1][ns2->Dofs()[1]] -= 1;
  DlS[2][ns1->Dofs()[2]] += 1;
  DlS[2][ns2->Dofs()[2]] -= 1;

  for (CI p = DlS[0].begin(); p != DlS[0].end(); ++p)
    (DlengthS)[p->first] += (p->second) * vs[0] * 1.0 / (lengthS);
  for (CI p = DlS[1].begin(); p != DlS[1].end(); ++p)
    (DlengthS)[p->first] += (p->second) * vs[1] * 1.0 / (lengthS);
  for (CI p = DlS[2].begin(); p != DlS[2].end(); ++p)
    (DlengthS)[p->first] += (p->second) * vs[2] * 1.0 / (lengthS);

  // lin lengthS * lengthM
  GEN::pairedvector<int, double> prodLength(1000);

  for (CI p = DlengthS.begin(); p != DlengthS.end(); ++p)
    (prodLength)[p->first] += (p->second) * lengthM;
  for (CI p = DlengthM.begin(); p != DlengthM.end(); ++p)
    (prodLength)[p->first] += (p->second) * lengthS;

  // lin scaprod
  GEN::pairedvector<int, double> scaProdlin(1000);

  for (CI p = DlS[0].begin(); p != DlS[0].end(); ++p) (scaProdlin)[p->first] += (p->second) * vm[0];
  for (CI p = DlS[1].begin(); p != DlS[1].end(); ++p) (scaProdlin)[p->first] += (p->second) * vm[1];
  for (CI p = DlS[2].begin(); p != DlS[2].end(); ++p) (scaProdlin)[p->first] += (p->second) * vm[2];

  for (CI p = DlM[0].begin(); p != DlM[0].end(); ++p) (scaProdlin)[p->first] += (p->second) * vs[0];
  for (CI p = DlM[1].begin(); p != DlM[1].end(); ++p) (scaProdlin)[p->first] += (p->second) * vs[1];
  for (CI p = DlM[2].begin(); p != DlM[2].end(); ++p) (scaProdlin)[p->first] += (p->second) * vs[2];

  // lin scaprod/lengthprod
  GEN::pairedvector<int, double> scaProdnormalizedLin(1000);
  for (CI p = scaProdlin.begin(); p != scaProdlin.end(); ++p)
    (scaProdnormalizedLin)[p->first] += (p->second) * 1.0 / (lengthS * lengthM);
  for (CI p = prodLength.begin(); p != prodLength.end(); ++p)
    (scaProdnormalizedLin)[p->first] -=
        (p->second) * scaprod * 1.0 / (lengthS * lengthM * lengthS * lengthM);

  // lin acos(scaledscaprod)
  double fac = (-1.0 / (sqrt(1.0 - scaledScaprod * scaledScaprod)));
  //  if(switchSign)
  //    fac *= -1.0;

  for (CI p = scaProdnormalizedLin.begin(); p != scaProdnormalizedLin.end(); ++p)
    (lineAngle)[p->first] += (p->second) * fac;

  // bye bye
  return angleRad;
}
