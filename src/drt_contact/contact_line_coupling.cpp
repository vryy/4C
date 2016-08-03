/*---------------------------------------------------------------------*/
/*!
\file contact_line_coupling.cpp

\brief A class to perform line clipping for line to surface contact +
       call for numerical integration

\level 2

\maintainer Philipp Farah

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
CONTACT::LineCoupling3d::LineCoupling3d(
    DRT::Discretization& idiscret,
    int dim,
    Teuchos::ParameterList& params,
    CoElement& pEle,
    Teuchos::RCP<MORTAR::MortarElement>& lEle,
    CoElement& surfEle,
    LineCoupling3d::intType type) :
  idiscret_(idiscret),
  dim_(dim),
  pEle_(pEle),
  lEle_(lEle),
  surfEle_(surfEle),
  imortar_(params),
  intType_(type)
{
  // empty constructor

  return;
}

/*----------------------------------------------------------------------*
 |  eval (public)                                            farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineCoupling3d::EvaluateCoupling()
{
  // 1. create aux plane for master ele
  AuxiliaryPlane();

  // 2. check orientation
  if(!CheckOrientation())
    return;

  // 3. project master nodes onto auxplane
  ProjectMaster();

  // 4. project slave line elements onto auxplane
  ProjectSlave();

  // 5. perform line clipping
  LineClipping();

  // 6. intersections found?
  if((int)InterSections().size() == 0 or (int)InterSections().size() == 1 )
    return;

  // 7. check length of Integration Line
  bool check = CheckLength();
  if(check==false)
    return;

  // create empty lin vector
  std::vector<std::vector<GEN::pairedvector<int, double> > > linvertex(
      2,
      std::vector<GEN::pairedvector<int, double> >(
          3,
          3 * LineElement()->NumNode() + 3 * SurfaceElement().NumNode()));

  // 8. linearize vertices
  LinearizeVertices(linvertex);

  // 9. create intlines
  CreateIntegrationLines(linvertex);

  // 10. consistent dual shape
  ConsistDualShape();

  // 11. integration
  IntegrateLine();

  return;
}

/*----------------------------------------------------------------------*
 |  calculate dual shape functions                           farah 07/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineCoupling3d::CheckOrientation()
{
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
  const double lengthA = sqrt(Auxn()[0] * Auxn()[0] + Auxn()[1] * Auxn()[1] + Auxn()[2] * Auxn()[2]);
  const double prod = lengthS*lengthA;
  if(prod<1e-12)
    return false;

  // calculate scalar product
  double scaprod = lvec[0] * Auxn()[0] + lvec[1] * Auxn()[1] + lvec[2] * Auxn()[2];
  scaprod = scaprod/(prod);
  double diff = abs(scaprod) -1.0;

  if(abs(diff)<tol)
    return false;

  return true;
}


/*----------------------------------------------------------------------*
 |  calculate dual shape functions                           farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineCoupling3d::ConsistDualShape()
{
  INPAR::MORTAR::ShapeFcn shapefcn =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(imortar_,"LM_SHAPEFCN");
  INPAR::MORTAR::ConsistentDualType consistent =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ConsistentDualType>(imortar_,"LM_DUAL_CONSISTENT");

  if (shapefcn != INPAR::MORTAR::shape_dual &&
      shapefcn != INPAR::MORTAR::shape_petrovgalerkin)
    return;

  if (consistent==INPAR::MORTAR::consistent_none)
    return;

  if(IType() == LineCoupling3d::lts)
    return;
  dserror("ERROR: consistent dual shapes for stl is experimental!");

  // slave nodes and dofs
  const int max_nnodes=9;
  const int nnodes = SurfaceElement().NumNode();
  if (nnodes>max_nnodes)
    dserror("this function is not implemented to handle elements with that many nodes. Just adjust max_nnodes above");
  const int ndof = 3;

  // get number of master nodes
  int mnodes = LineElement()->NumNode();

  // Dual shape functions coefficient matrix and linearization
  LINALG::SerialDenseMatrix ae(nnodes,nnodes,true);
  SurfaceElement().MoData().DerivDualShape() =
      Teuchos::rcp(new GEN::pairedvector<int,Epetra_SerialDenseMatrix>((nnodes+mnodes)*ndof,0,Epetra_SerialDenseMatrix(nnodes,nnodes)));
  GEN::pairedvector<int,Epetra_SerialDenseMatrix>& derivae=*(SurfaceElement().MoData().DerivDualShape());

  // various variables
  double detg=0.0;
  typedef GEN::pairedvector<int,double>::const_iterator _CI;

  // initialize matrices de and me
  LINALG::SerialDenseMatrix me(nnodes,nnodes,true);
  LINALG::SerialDenseMatrix de(nnodes,nnodes,true);

  GEN::pairedvector<int,LINALG::Matrix<max_nnodes+1,max_nnodes> > derivde_new((nnodes+mnodes)*ndof);

  // two-dim arrays of maps for linearization of me/de
  std::vector<std::vector<GEN::pairedvector<int,double> > >
    derivme(nnodes,std::vector<GEN::pairedvector<int,double> >(nnodes,(nnodes+mnodes)*ndof));
  std::vector<std::vector<GEN::pairedvector<int,double> > >
    derivde(nnodes,std::vector<GEN::pairedvector<int,double> >(nnodes,(nnodes+mnodes)*ndof));

  double A_tot=0.;

  // get number of master nodes
  const int ncol = LineElement()->NumNode();

  Teuchos::RCP<MORTAR::IntCell> currcell = IntLine();

  A_tot+=currcell->Area();

  // create an integrator for this cell
  CONTACT::CoIntegrator integrator(imortar_,currcell->Shape(),Comm());

  // check if the cells are tri3
  // there's nothing wrong about other shapes, but as long as they are all
  // tri3 we can perform the jacobian calculation ( and its deriv) outside
  // the Gauss point loop
  if (currcell->Shape()!=DRT::Element::line2)
    dserror("only line2 integration cells at the moment. See comment in the code");
  double eta[2]={0.,0.};
  detg=currcell->Jacobian(eta);
  // directional derivative of cell Jacobian
  GEN::pairedvector<int,double> derivjaccell((nnodes+ncol)*ndof);
  currcell->DerivJacobian(eta, derivjaccell);

  for (int gp=0;gp<integrator.nGP(); ++gp)
  {
    // coordinates and weight
    double eta[2]    = { integrator.Coordinate(gp,0), integrator.Coordinate(gp,1)};
    const double wgt = integrator.Weight(gp);

    // get global Gauss point coordinates
    double globgp[3] = { 0.0, 0.0, 0.0};
    currcell->LocalToGlobal(eta, globgp,0);

    // project Gauss point onto slave integration element
    double sxi[2] = { 0.0, 0.0};
    double sprojalpha = 0.0;
    MORTAR::MortarProjector::Impl(SurfaceElement())->
        ProjectGaussPointAuxn3D(globgp, Auxn(), SurfaceElement(), sxi, sprojalpha);

    // project Gauss point onto slave (parent) element
    double psxi[2] = {0.,0.};
    for (int i=0; i<2; ++i)
      psxi[i]=sxi[i];

    // create vector for shape function evaluation
    LINALG::SerialDenseVector sval (nnodes);
    LINALG::SerialDenseMatrix sderiv(nnodes,2,true);

    // evaluate trace space shape functions at Gauss point
    SurfaceElement().EvaluateShape(psxi, sval, sderiv, nnodes);

    // additional data for contact calculation (i.e. incl. derivative of dual shape functions coefficient matrix)
    // GP slave coordinate derivatives
    std::vector<GEN::pairedvector<int,double> > dsxigp(2,(nnodes+ncol)*ndof);
    // GP slave coordinate derivatives
    std::vector<GEN::pairedvector<int,double> > dpsxigp(2,(nnodes+ncol)*ndof);
    // global GP coordinate derivative on integration element
    GEN::pairedvector<int,LINALG::Matrix<3,1> > lingp((nnodes+ncol)*ndof);

    // compute global GP coordinate derivative
    static LINALG::Matrix<3,1> svalcell;
    static LINALG::Matrix<3,2> sderivcell;
    currcell->EvaluateShape(eta,svalcell,sderivcell);

    for (int v=0;v<2;++v)
      for (int d=0; d<3; ++d)
        for (_CI p=(currcell->GetDerivVertex(v))[d].begin();p!=(currcell->GetDerivVertex(v))[d].end();++p)
          lingp[p->first](d) += svalcell(v) * (p->second);

    // compute GP slave coordinate derivatives
    integrator.DerivXiGP3DAuxPlane(SurfaceElement(),sxi,currcell->Auxn(),dsxigp,sprojalpha,currcell->GetDerivAuxn(),lingp);

    // compute GP slave coordinate derivatives (parent element)
    dpsxigp=dsxigp;

    double fac=0.;
    for (_CI p=derivjaccell.begin();p!=derivjaccell.end();++p)
    {
      LINALG::Matrix<max_nnodes+1,max_nnodes>& dtmp = derivde_new[p->first];
      const double& ps = p->second;
      for (int j=0;j<nnodes;++j)
      {
        fac=wgt*sval[j]*ps;
        dtmp(nnodes,j)+=fac;
        for (int k=0; k<nnodes; ++k)
          dtmp(k,j)+=fac*sval[k];
      }
    }

    for (int i=0; i<2; ++i)
      for (_CI p=dpsxigp[i].begin(); p!=dpsxigp[i].end(); ++p)
      {
        LINALG::Matrix<max_nnodes+1,max_nnodes>& dtmp = derivde_new[p->first];
        const double& ps = p->second;
        for (int j=0;j<nnodes;++j)
        {
          fac = wgt*sderiv(j,i)*detg*ps;
          dtmp(nnodes,j) += fac;
          for (int k=0; k<nnodes; ++k)
          {
            dtmp(k,j)+=fac*sval[k];
            dtmp(j,k)+=fac*sval[k];
          }
        }
      }

    // computing de, derivde and me, derivme and kappa, derivkappa
    for (int j=0; j<nnodes; ++j)
    {
      double fac;
      fac = sval[j]*wgt;
      // computing de
      de(j,j)+=fac*detg;

      for (int k=0; k<nnodes; ++k)
      {
        // computing me
        fac = wgt*sval[j]*sval[k];
        me(j,k)+=fac*detg;
      }
    }
  }

  // in case of no overlap just return, as there is no integration area
  // and therefore the consistent dual shape functions are not defined.
  // This doesn't matter, as there is no associated integration domain anyway
  if (A_tot<1.e-12)
    return;

  // invert bi-ortho matrix me
//  LINALG::SerialDenseMatrix meinv = LINALG::InvertAndMultiplyByCholesky(me,de,ae);

  LINALG::NonSymmetricInverse(me,4);
  LINALG::SerialDenseMatrix meinv = me;

  // build linearization of ae and store in derivdual
  // (this is done according to a quite complex formula, which
  // we get from the linearization of the biorthogonality condition:
  // Lin (Me * Ae = De) -> Lin(Ae)=Lin(De)*Inv(Me)-Ae*Lin(Me)*Inv(Me) )
  typedef GEN::pairedvector<int,LINALG::Matrix<max_nnodes+1,max_nnodes> >::const_iterator _CIM;
  for (_CIM p=derivde_new.begin();p!=derivde_new.end();++p)
  {
    LINALG::Matrix<max_nnodes+1,max_nnodes>& dtmp = derivde_new[p->first];
    Epetra_SerialDenseMatrix& pt = derivae[p->first];
    for (int i=0;i<nnodes;++i)
      for (int j=0;j<nnodes;++j)
      {
          pt(i,j) += meinv(i,j)*dtmp(nnodes,i);

        for (int k=0; k<nnodes; ++k)
          for (int l=0; l<nnodes; ++l)
            pt(i,j) -= ae(i,k)*meinv(l,j)*dtmp(l,k);
      }
  }

  // store ae matrix in slave element data container
  SurfaceElement().MoData().DualShape() = Teuchos::rcp(new LINALG::SerialDenseMatrix(ae));

  return;
}

/*----------------------------------------------------------------------*
 |  integration for LTS (public)                             farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineCoupling3d::IntegrateLine()
{
  // get solution strategy
  INPAR::CONTACT::SolvingStrategy sol =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(imortar_,"STRATEGY");

  // create integrator object
  Teuchos::RCP<CONTACT::CoIntegrator> integrator =
      CONTACT::INTEGRATOR::BuildIntegrator(
          sol,
          imortar_,
          IntLine()->Shape(),
          Comm());

  // perform integration
  if(IType() == LineCoupling3d::lts)
  {
    integrator->IntegrateDerivCell3DAuxPlaneLTS(
        ParentElement(),
        *LineElement(),
        SurfaceElement(),
        IntLine(),
        Auxn(),
        Comm());
  }
  else if(IType() == LineCoupling3d::stl)
  {
    integrator->IntegrateDerivCell3DAuxPlaneSTL(
        ParentElement(),
        *LineElement(),
        SurfaceElement(),
        IntLine(),
        Auxn(),
        Comm());
  }
  else
    dserror("ERROR: wrong integration type for line coupling!");

  return;
}

/*----------------------------------------------------------------------*
 |  geometric stuff (private)                                farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineCoupling3d::LineClipping()
{
  // output variable
  const bool out = false;

  // tolerance for line clipping
  const double sminedge = ParentElement().MinEdgeSize();
  const double mminedge = SurfaceElement().MinEdgeSize();
  const double tol = MORTARCLIPTOL * std::min(sminedge, mminedge);

  // vector with vertices
  InterSections().clear();

  // vector with temp intersections vertices
  std::vector<MORTAR::Vertex> tempintersec;

  // safety
  if(MasterVertices().size()<3)
    dserror("ERROR: Invalid number of Master Vertices!");
  if(SlaveVertices().size()!=2)
    dserror("ERROR: Invalid number of Slave Vertices!");

  // set previous and next Vertex pointer for all elements in lists
  for (int i = 0; i < (int) MasterVertices().size(); ++i)
  {
    // standard case
    if (i != 0 && i != (int) MasterVertices().size() - 1)
    {
      MasterVertices()[i].AssignNext(&MasterVertices()[i + 1]);
      MasterVertices()[i].AssignPrev(&MasterVertices()[i - 1]);
    }
    // first element in list
    else if (i == 0)
    {
      MasterVertices()[i].AssignNext(&MasterVertices()[i + 1]);
      MasterVertices()[i].AssignPrev(&MasterVertices()[(int) MasterVertices().size() - 1]);
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
  double slaveLine[3] = { 0.0, 0.0, 0.0 };
  for (int k = 0; k < 3; ++k)
    slaveLine[k] = SlaveVertices()[1].Coord()[k] - SlaveVertices()[0].Coord()[k];

  // loop over master vertices to create master polygon lines
  for (int j = 0; j < (int) MasterVertices().size(); ++j)
  {
    // we need two edges first
    double edge[3] = { 0.0, 0.0, 0.0 };
    for (int k = 0; k < 3; ++k)
      edge[k] = (MasterVertices()[j].Next())->Coord()[k] - MasterVertices()[j].Coord()[k];

    // outward edge normals of polygon and slave line
    double np[3] = { 0.0, 0.0, 0.0 };
    double nl[3] = { 0.0, 0.0, 0.0 };
    np[0] = edge[1] * Auxn()[2] - edge[2] * Auxn()[1];
    np[1] = edge[2] * Auxn()[0] - edge[0] * Auxn()[2];
    np[2] = edge[0] * Auxn()[1] - edge[1] * Auxn()[0];
    nl[0] = slaveLine[1] * Auxn()[2] - slaveLine[2] * Auxn()[1];
    nl[1] = slaveLine[2] * Auxn()[0] - slaveLine[0] * Auxn()[2];
    nl[2] = slaveLine[0] * Auxn()[1] - slaveLine[1] * Auxn()[0];

    if(out)
    {
      std::cout << "==============================================" << std::endl;
      std::cout << "SLine= " << slaveLine[0] <<"  " << slaveLine[1] << "  " << slaveLine[2] << std::endl;
      std::cout << "Pos1= " << SlaveVertices()[0].Coord()[0] <<"  " << SlaveVertices()[0].Coord()[1] << "  " << SlaveVertices()[0].Coord()[2] << std::endl;
      std::cout << "Pos2= " << SlaveVertices()[1].Coord()[0] <<"  " << SlaveVertices()[1].Coord()[1] << "  " << SlaveVertices()[1].Coord()[2] << std::endl;
      std::cout << "N slave= " << nl[0] <<"  " << nl[1] << "  " << nl[2] << std::endl;


      std::cout << "==============================================" << std::endl;
      std::cout << "MEdge= " << edge[0] <<"  " << edge[1] << "  " << edge[2] << std::endl;
      std::cout << "Pos1= " << (MasterVertices()[j].Next())->Coord()[0] <<"  " << (MasterVertices()[j].Next())->Coord()[1] << "  " << (MasterVertices()[j].Next())->Coord()[2] << std::endl;
      std::cout << "Pos2= " << MasterVertices()[j].Coord()[0] <<"  " << MasterVertices()[j].Coord()[1] << "  " << MasterVertices()[j].Coord()[2] << std::endl;
      std::cout << "N master= " << np[0] <<"  " << np[1] << "  " << np[2] << std::endl;
    }

    // check for parallelity of edges
    double parallel = edge[0] * nl[0] + edge[1] * nl[1] + edge[2] * nl[2];
    if (abs(parallel) < tol)
    {
      if (out)
        std::cout << "WARNING: Detected two parallel edges! (" << j << "," << j << ")" << std::endl;
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

    if(out)
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
        wec_q1 += (MasterVertices()[j].Coord()[k]           - SlaveVertices()[0].Coord()[k]) * nl[k];
        wec_q2 += ((MasterVertices()[j].Next())->Coord()[k] - SlaveVertices()[0].Coord()[k]) * nl[k];
      }

      if(out)
      {
        std::cout << "WecQ1 = " << wec_q1 << std::endl;
        std::cout << "WecQ2 = " << wec_q2 << std::endl;
      }

      if (wec_q1 * wec_q2 <= 0.0)
      {
        double alpha  = wec_p1 / (wec_p1 - wec_p2);
        double alphaq = wec_q1 / (wec_q1 - wec_q2);

        if(alpha<0.0 or alpha >1.0)
          continue;
        if(alphaq<0.0 or alphaq >1.0)
          continue;

        std::vector<double> coords(3);
        for (int k = 0; k < 3; ++k)
        {
          coords[k] = (1 - alpha) * SlaveVertices()[0].Coord()[k]  + alpha * SlaveVertices()[1].Coord()[k];
          if (abs(coords[k]) < tol)
            coords[k] = 0.0;
        }

        if (out)
        {
          std::cout << "Found intersection! (" << j << ") " << alpha << std::endl;
          std::cout << "coords 1: " << coords[0] << " " << coords[1] << " " << coords[2] << std::endl;
        }

        // generate vectors of underlying node ids for lineclip (2x slave, 2x master)
        std::vector<int> lcids(4);
        lcids[0] = (int) (SlaveVertices()[0].Nodeids()[0]);
        lcids[1] = (int) (SlaveVertices()[1].Nodeids()[0]);
        lcids[2] = (int) (MasterVertices()[j].Nodeids()[0]);
        lcids[3] = (int) ((MasterVertices()[j].Next())->Nodeids()[0]);

        // store intersection points
        tempintersec.push_back(
            MORTAR::Vertex(
                coords,
                MORTAR::Vertex::lineclip,
                lcids,
                &SlaveVertices()[1],
                &SlaveVertices()[0],
                true,
                false,
                NULL,
                alpha));
      }
    }
  }

  // check if intersection is close to a node
  for (int i = 0; i < (int) tempintersec.size(); ++i)
  {
    // keep track of comparisons
    bool close = false;

    // check against all poly1 (slave) points
    for (int j = 0; j < (int) SlaveVertices().size(); ++j)
    {
      // distance vector
      double diff[3] = { 0.0, 0.0, 0.0 };
      for (int k = 0; k < 3; ++k)
        diff[k] = tempintersec[i].Coord()[k] - SlaveVertices()[j].Coord()[k];
      double dist = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);

      // only keep intersection point if not close
      if (dist <= tol)
      {
        close = true;
        break;
      }
    }

    // do only if no close poly1 (slave) point found
    if (!close)
    {
      // check against all poly2 (master) points
      for (int j = 0; j < (int) MasterVertices().size(); ++j)
      {
        // distance vector
        double diff[3] = { 0.0, 0.0, 0.0 };
        for (int k = 0; k < 3; ++k)
          diff[k] = tempintersec[i].Coord()[k] - MasterVertices()[j].Coord()[k];
        double dist = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);

        // only keep intersection point if not close
        if (dist <= tol)
        {
          close = true;
          break;
        }
      }
    }

    // keep intersection point only if not close to any Slave/Master point
    if (!close)
      InterSections().push_back(tempintersec[i]);
  }

  // 1. intersection: find interior node
  if((int)InterSections().size()!=2)
  {
    for (int i = 0; i < (int) SlaveVertices().size(); ++i)
    {
      // keep track of inside / outside status
      bool outside = false;

      // check against all poly1 (slave) edges
      for (int j = 0; j < (int) MasterVertices().size(); ++j)
      {
        // we need diff vector and edge2 first
        double diff[3] = { 0.0, 0.0, 0.0 };
        double edge[3] = { 0.0, 0.0, 0.0 };
        for (int k = 0; k < 3; ++k)
        {
          diff[k] = SlaveVertices()[i].Coord()[k]            - MasterVertices()[j].Coord()[k];
          edge[k] = (MasterVertices()[j].Next())->Coord()[k] - MasterVertices()[j].Coord()[k];
        }

        // compute distance from point on poly1 to edge
        double n[3] = { 0.0, 0.0, 0.0 };
        n[0] = edge[1] * Auxn()[2] - edge[2] * Auxn()[1];
        n[1] = edge[2] * Auxn()[0] - edge[0] * Auxn()[2];
        n[2] = edge[0] * Auxn()[1] - edge[1] * Auxn()[0];
        double ln = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        for (int k = 0; k < 3; ++k)
          n[k] /= ln;

        double dist = diff[0] * n[0] + diff[1] * n[1] + diff[2] * n[2];

        // only keep point if not in outside halfspace
        if (dist > tol)
        {
          outside = true;
          break;
        }
      } // end master loop

      if(outside)
      {
        continue;
      }
      else
      {
        InterSections().push_back(
            MORTAR::Vertex(
                SlaveVertices()[i].Coord(),
                MORTAR::Vertex::projslave,
                SlaveVertices()[i].Nodeids(),
                NULL,
                NULL,
                false,
                false,
                NULL,
                -1));
      }
    }
  }

  // 2. check plausibility
  if(InterSections().size()>2)
  {
    std::cout << "Intersections= " << InterSections().size() << std::endl;
    dserror("ERROR: intersections not possible!!!");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  create integration lines                                 farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineCoupling3d::CreateIntegrationLines(std::vector<std::vector<GEN::pairedvector<int, double> > >& linvertex)
{
  // get coordinates
  LINALG::Matrix<3,3> coords;

  for(int i = 0; i< 3; ++i)
  {
    coords(i,0) = InterSections()[0].Coord()[i];
    coords(i,1) = InterSections()[1].Coord()[i];
    coords(i,2) = -1; // dummy;
  }

  // create Integration Line
  IntLine() = Teuchos::rcp(new MORTAR::IntCell(
      ParentElement().Id(),
      2,
      coords,
      Auxn(),
      DRT::Element::line2,
      linvertex[0],
      linvertex[1],
      linvertex[1], // dummy
      GetDerivAuxn()));

  return;
}

/*----------------------------------------------------------------------*
 |  linearize vertices                                       farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineCoupling3d::LinearizeVertices(std::vector<std::vector<GEN::pairedvector<int, double> > >& linvertex)
{
  // linearize all aux.plane slave and master nodes only ONCE
  // and use these linearizations later during lineclip linearization
  // (this speeds up the vertex linearizations in most cases, as we
  // never linearize the SAME slave or master vertex more than once)

  // number of nodes
  const int nsrows = LineElement()->NumNode();
  const int nmrows = SurfaceElement().NumNode();

  // prepare storage for slave and master linearizations
  std::vector<std::vector<GEN::pairedvector<int, double> > > linsnodes(nsrows,
      std::vector<GEN::pairedvector<int, double> >(3, 3 * LineElement()->NumNode() + 3 * SurfaceElement().NumNode()));
  std::vector<std::vector<GEN::pairedvector<int, double> > > linmnodes(nmrows,
      std::vector<GEN::pairedvector<int, double> >(3, 3 * LineElement()->NumNode() + 3 * SurfaceElement().NumNode()));

  // compute slave linearizations (nsrows)
  SlaveVertexLinearization(linsnodes);

  // compute master linearizations (nmrows)
  MasterVertexLinearization(linmnodes);

  //**********************************************************************
  // Line vertex linearization
  //**********************************************************************
  // loop over all clip Intersections vertices
  for (int i = 0; i < (int) InterSections().size(); ++i)
  {
    // references to current vertex and its linearization
    MORTAR::Vertex& currv = InterSections()[i];
    std::vector<GEN::pairedvector<int, double> >& currlin = linvertex[i];

    // decision on vertex type (slave, projmaster, linclip)
    if (currv.VType() == MORTAR::Vertex::projslave)
    {
      // get corresponding slave id
      int sid = currv.Nodeids()[0];

      // find corresponding slave node linearization
      int k = 0;
      while (k < nsrows)
      {
        if (LineElement()->NodeIds()[k] == sid)
          break;
        ++k;
      }

      // dserror if not found
      if (k == nsrows)
        dserror("ERROR: Slave Id not found!");

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
        if (SurfaceElement().NodeIds()[k] == mid)
          break;
        ++k;
      }

      // dserror if not found
      if (k == nmrows)
        dserror("ERROR: Master Id not found!");

      // get the correct master node linearization
      currlin = linmnodes[k];
    }
    else if (currv.VType() == MORTAR::Vertex::lineclip)
    {
      // get references to the two slave vertices
      int sindex1 = -1;
      int sindex2 = -1;
      for (int j = 0; j < (int) SlaveVertices().size(); ++j)
      {
        if (SlaveVertices()[j].Nodeids()[0] == currv.Nodeids()[0])
          sindex1 = j;
        if (SlaveVertices()[j].Nodeids()[0] == currv.Nodeids()[1])
          sindex2 = j;
      }
      if (sindex1 < 0 || sindex2 < 0 || sindex1 == sindex2)
        dserror("ERROR: Lineclip linearization: (S) Something went wrong!");

      MORTAR::Vertex* sv1 = &SlaveVertices()[sindex1];
      MORTAR::Vertex* sv2 = &SlaveVertices()[sindex2];

      // get references to the two master vertices
      int mindex1 = -1;
      int mindex2 = -1;
      for (int j = 0; j < (int) MasterVertices().size(); ++j)
      {
        if (MasterVertices()[j].Nodeids()[0] == currv.Nodeids()[2])
          mindex1 = j;
        if (MasterVertices()[j].Nodeids()[0] == currv.Nodeids()[3])
          mindex2 = j;
      }
      if (mindex1 < 0 || mindex2 < 0 || mindex1 == mindex2)
        dserror("ERROR: Lineclip linearization: (M) Something went wrong!");

      MORTAR::Vertex* mv1 = &MasterVertices()[mindex1];
      MORTAR::Vertex* mv2 = &MasterVertices()[mindex2];

      // do lineclip vertex linearization
      LineclipVertexLinearization(
          currv,
          currlin,
          sv1,
          sv2,
          mv1,
          mv2,
          linsnodes,
          linmnodes);
    }

    else
      dserror("ERROR: VertexLinearization: Invalid Vertex Type!");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Linearization of lineclip vertex (3D) AuxPlane            popp 03/09|
 *----------------------------------------------------------------------*/
void CONTACT::LineCoupling3d::LineclipVertexLinearization(MORTAR::Vertex& currv,
    std::vector<GEN::pairedvector<int, double> >& currlin,
    MORTAR::Vertex* sv1, MORTAR::Vertex* sv2,
    MORTAR::Vertex* mv1, MORTAR::Vertex* mv2,
    std::vector<std::vector<GEN::pairedvector<int, double> > >& linsnodes,
    std::vector<std::vector<GEN::pairedvector<int, double> > >& linmnodes)
{
  // number of nodes
  const int nsrows = LineElement()->NumNode();
  const int nmrows = SurfaceElement().NumNode();

  // iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // compute factor Z
  double crossZ[3] =
  { 0.0, 0.0, 0.0 };
  crossZ[0] =   (sv1->Coord()[1] - mv1->Coord()[1])
              * (mv2->Coord()[2] - mv1->Coord()[2])
              - (sv1->Coord()[2] - mv1->Coord()[2])
              * (mv2->Coord()[1] - mv1->Coord()[1]);
  crossZ[1] =   (sv1->Coord()[2] - mv1->Coord()[2])
              * (mv2->Coord()[0] - mv1->Coord()[0])
              - (sv1->Coord()[0] - mv1->Coord()[0])
              * (mv2->Coord()[2] - mv1->Coord()[2]);
  crossZ[2] =   (sv1->Coord()[0] - mv1->Coord()[0])
              * (mv2->Coord()[1] - mv1->Coord()[1])
              - (sv1->Coord()[1] - mv1->Coord()[1])
          * (mv2->Coord()[0] - mv1->Coord()[0]);
  double Zfac = crossZ[0] * Auxn()[0] + crossZ[1] * Auxn()[1]
      + crossZ[2] * Auxn()[2];

  // compute factor N
  double crossN[3] = { 0.0, 0.0, 0.0 };
  crossN[0] =   (sv2->Coord()[1] - sv1->Coord()[1])
              * (mv2->Coord()[2] - mv1->Coord()[2])
              - (sv2->Coord()[2] - sv1->Coord()[2])
              * (mv2->Coord()[1] - mv1->Coord()[1]);
  crossN[1] =   (sv2->Coord()[2] - sv1->Coord()[2])
              * (mv2->Coord()[0] - mv1->Coord()[0])
              - (sv2->Coord()[0] - sv1->Coord()[0])
              * (mv2->Coord()[2] - mv1->Coord()[2]);
  crossN[2] =   (sv2->Coord()[0] - sv1->Coord()[0])
              * (mv2->Coord()[1] - mv1->Coord()[1])
              - (sv2->Coord()[1] - sv1->Coord()[1])
              * (mv2->Coord()[0] - mv1->Coord()[0]);
  double Nfac =   crossN[0] * Auxn()[0] + crossN[1] * Auxn()[1]
                + crossN[2] * Auxn()[2];

  // slave edge vector
  double sedge[3] = { 0.0, 0.0, 0.0 };
  for (int k = 0; k < 3; ++k)
    sedge[k] = sv2->Coord()[k] - sv1->Coord()[k];

  // prepare linearization derivZ
  double crossdZ1[3] = { 0.0, 0.0, 0.0 };
  double crossdZ2[3] = { 0.0, 0.0, 0.0 };
  double crossdZ3[3] = { 0.0, 0.0, 0.0 };
  crossdZ1[0] = (mv2->Coord()[1] - mv1->Coord()[1]) * Auxn()[2]
      - (mv2->Coord()[2] - mv1->Coord()[2]) * Auxn()[1];
  crossdZ1[1] = (mv2->Coord()[2] - mv1->Coord()[2]) * Auxn()[0]
      - (mv2->Coord()[0] - mv1->Coord()[0]) * Auxn()[2];
  crossdZ1[2] = (mv2->Coord()[0] - mv1->Coord()[0]) * Auxn()[1]
      - (mv2->Coord()[1] - mv1->Coord()[1]) * Auxn()[0];
  crossdZ2[0] = Auxn()[1] * (sv1->Coord()[2] - mv1->Coord()[2])
      - Auxn()[2] * (sv1->Coord()[1] - mv1->Coord()[1]);
  crossdZ2[1] = Auxn()[2] * (sv1->Coord()[0] - mv1->Coord()[0])
      - Auxn()[0] * (sv1->Coord()[2] - mv1->Coord()[2]);
  crossdZ2[2] = Auxn()[0] * (sv1->Coord()[1] - mv1->Coord()[1])
      - Auxn()[1] * (sv1->Coord()[0] - mv1->Coord()[0]);
  crossdZ3[0] = (sv1->Coord()[1] - mv1->Coord()[1])
      * (mv2->Coord()[2] - mv1->Coord()[2])
      - (sv1->Coord()[2] - mv1->Coord()[2])
          * (mv2->Coord()[1] - mv1->Coord()[1]);
  crossdZ3[1] = (sv1->Coord()[2] - mv1->Coord()[2])
      * (mv2->Coord()[0] - mv1->Coord()[0])
      - (sv1->Coord()[0] - mv1->Coord()[0])
          * (mv2->Coord()[2] - mv1->Coord()[2]);
  crossdZ3[2] = (sv1->Coord()[0] - mv1->Coord()[0])
      * (mv2->Coord()[1] - mv1->Coord()[1])
      - (sv1->Coord()[1] - mv1->Coord()[1])
          * (mv2->Coord()[0] - mv1->Coord()[0]);

  // prepare linearization derivN
  double crossdN1[3] = { 0.0, 0.0, 0.0 };
  double crossdN2[3] = { 0.0, 0.0, 0.0 };
  double crossdN3[3] = { 0.0, 0.0, 0.0 };
  crossdN1[0] =   (mv2->Coord()[1] - mv1->Coord()[1]) * Auxn()[2]
                - (mv2->Coord()[2] - mv1->Coord()[2]) * Auxn()[1];
  crossdN1[1] =   (mv2->Coord()[2] - mv1->Coord()[2]) * Auxn()[0]
                - (mv2->Coord()[0] - mv1->Coord()[0]) * Auxn()[2];
  crossdN1[2] =   (mv2->Coord()[0] - mv1->Coord()[0]) * Auxn()[1]
                - (mv2->Coord()[1] - mv1->Coord()[1]) * Auxn()[0];
  crossdN2[0] =   Auxn()[1] * (sv2->Coord()[2] - sv1->Coord()[2])
                - Auxn()[2] * (sv2->Coord()[1] - sv1->Coord()[1]);
  crossdN2[1] =   Auxn()[2] * (sv2->Coord()[0] - sv1->Coord()[0])
                - Auxn()[0] * (sv2->Coord()[2] - sv1->Coord()[2]);
  crossdN2[2] =   Auxn()[0] * (sv2->Coord()[1] - sv1->Coord()[1])
                - Auxn()[1] * (sv2->Coord()[0] - sv1->Coord()[0]);
  crossdN3[0] =   (sv2->Coord()[1] - sv1->Coord()[1])
                * (mv2->Coord()[2] - mv1->Coord()[2])
                - (sv2->Coord()[2] - sv1->Coord()[2])
                * (mv2->Coord()[1] - mv1->Coord()[1]);
  crossdN3[1] =   (sv2->Coord()[2] - sv1->Coord()[2])
                * (mv2->Coord()[0] - mv1->Coord()[0])
                - (sv2->Coord()[0] - sv1->Coord()[0])
                * (mv2->Coord()[2] - mv1->Coord()[2]);
  crossdN3[2] =   (sv2->Coord()[0] - sv1->Coord()[0])
                * (mv2->Coord()[1] - mv1->Coord()[1])
                - (sv2->Coord()[1] - sv1->Coord()[1])
                * (mv2->Coord()[0] - mv1->Coord()[0]);

  // slave vertex linearization (2x)
  int sid1 = currv.Nodeids()[0];
  int sid2 = currv.Nodeids()[1];

  // find corresponding slave node linearizations
  int k = 0;
  while (k < nsrows)
  {
    if (LineElement()->NodeIds()[k] == sid1)
      break;
    ++k;
  }

  // dserror if not found
  if (k == nsrows)
    dserror("ERROR: Slave Id1 not found!");

  // get the correct slave node linearization
  std::vector<GEN::pairedvector<int, double> >& slavelin0 = linsnodes[k];

  k = 0;
  while (k < nsrows)
  {
    if (LineElement()->NodeIds()[k] == sid2)
      break;
    ++k;
  }

  // dserror if not found
  if (k == nsrows)
    dserror("ERROR: Slave Id2 not found!");

  // get the correct slave node linearization
  std::vector<GEN::pairedvector<int, double> >& slavelin1 = linsnodes[k];

  // master vertex linearization (2x)
  int mid1 = currv.Nodeids()[2];
  int mid2 = currv.Nodeids()[3];

  // find corresponding master node linearizations
  k = 0;
  while (k < nmrows)
  {
    if (SurfaceElement().NodeIds()[k] == mid1)
      break;
    ++k;
  }

  // dserror if not found
  if (k == nmrows)
    dserror("ERROR: Master Id1 not found!");

  // get the correct master node linearization
  std::vector<GEN::pairedvector<int,double> >& masterlin0 = linmnodes[k];

  k = 0;
  while (k < nmrows)
  {
    if (SurfaceElement().NodeIds()[k] == mid2)
      break;
    ++k;
  }

  // dserror if not found
  if (k == nmrows)
    dserror("ERROR: Master Id2 not found!");

  // get the correct master node linearization
  std::vector<GEN::pairedvector<int,double> >& masterlin1 = linmnodes[k];

  // linearization of element normal Auxn()
  std::vector<GEN::pairedvector<int,double> >& linauxn = GetDerivAuxn();

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
        currlin[dim][p->first] -= sedge[dim] * Nfacinv * crossdZ1[k]
            * (p->second);
        currlin[dim][p->first] -= sedge[dim] * ZNNfac * crossdN1[k]
            * (p->second);

      }
    }
    for (_CI p = slavelin1[k].begin(); p != slavelin1[k].end(); ++p)
    {
      currlin[k][p->first] -= ZNfac * (p->second);
      for (int dim = 0; dim < 3; ++dim)
      {
        currlin[dim][p->first] += sedge[dim] * ZNNfac * crossdN1[k]
            * (p->second);
      }
    }
    for (_CI p = masterlin0[k].begin(); p != masterlin0[k].end(); ++p)
    {
      for (int dim = 0; dim < 3; ++dim)
      {
        currlin[dim][p->first] += sedge[dim] * Nfacinv * crossdZ1[k]
            * (p->second);
        currlin[dim][p->first] += sedge[dim] * Nfacinv * crossdZ2[k]
            * (p->second);
        currlin[dim][p->first] -= sedge[dim] * ZNNfac * crossdN2[k]
            * (p->second);
      }
    }
    for (_CI p = masterlin1[k].begin(); p != masterlin1[k].end(); ++p)
    {
      for (int dim = 0; dim < 3; ++dim)
      {
        currlin[dim][p->first] -= sedge[dim] * Nfacinv * crossdZ2[k]
            * (p->second);
        currlin[dim][p->first] += sedge[dim] * ZNNfac * crossdN2[k]
            * (p->second);
      }
    }
    for (_CI p = linauxn[k].begin(); p != linauxn[k].end(); ++p)
    {
      for (int dim = 0; dim < 3; ++dim)
      {
        currlin[dim][p->first] -= sedge[dim] * Nfacinv * crossdZ3[k]
            * (p->second);
        currlin[dim][p->first] += sedge[dim] * ZNNfac * crossdN3[k]
            * (p->second);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  compute and check length of intLine                      farah 07/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineCoupling3d::CheckLength()
{
  // tolerance
  const double sminedge = ParentElement().MinEdgeSize();
  const double mminedge = SurfaceElement().MinEdgeSize();
  const double tol = MORTARCLIPTOL * std::min(sminedge, mminedge);

  // compute distance vector
  double v[3] = {0.0, 0.0, 0.0};
  for(int i = 0; i< 3; ++i)
    v[i]= InterSections()[0].Coord()[i] - InterSections()[1].Coord()[i];

  // compute length
  double length = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

  // check
  if(length<tol)
    return false;

  return true;
}

/*----------------------------------------------------------------------*
 |  eval (public)                                            farah 07/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineCoupling3d::AuxiliaryPlane()
{
  // we first need the element center:
  // for quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  double loccenter[2] = { 0.0, 0.0 };

  DRT::Element::DiscretizationType dt = SurfaceElement().Shape();
  if (dt == DRT::Element::tri3 || dt == DRT::Element::tri6)
  {
    loccenter[0] = 1.0 / 3.0;
    loccenter[1] = 1.0 / 3.0;
  }
  else if (dt == DRT::Element::quad4 || dt == DRT::Element::quad8
        || dt == DRT::Element::quad9)
  {
    loccenter[0] = 0.0;
    loccenter[1] = 0.0;
  }
  else
    dserror("ERROR: AuxiliaryPlane called for unknown element type");

  // compute element center via shape fct. interpolation
  SurfaceElement().LocalToGlobal(loccenter, Auxc(), 0);

  // we then compute the unit normal vector at the element center
  Lauxn() = SurfaceElement().ComputeUnitNormalAtXi(loccenter, Auxn());

  // compute aux normal linearization
  SurfaceElement().DerivUnitNormalAtXi(loccenter, GetDerivAuxn());

  // bye
  return true;
}


/*----------------------------------------------------------------------*
 |  eval (public)                                            farah 07/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineCoupling3d::HasProjStatus()
{

  return true;
}


/*----------------------------------------------------------------------*
 |  eval (public)                                            farah 07/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineCoupling3d::ProjectSlave()
{
  // project slave nodes onto auxiliary plane
  int nnodes = LineElement()->NumNode();

  // initialize storage for slave coords + their ids
  std::vector<double> vertices(3);
  std::vector<int> snodeids(1);

  for (int i = 0; i < nnodes; ++i)
  {
    DRT::Node* node = idiscret_.gNode(LineElement()->NodeIds()[i]);
    if (!node)
      dserror("ERROR: Cannot find slave element with gid %", LineElement()->NodeIds()[i]);
    CoNode* mycnode = dynamic_cast<CoNode*>(node);
    if (!mycnode)
      dserror("ERROR: ProjectSlave: Null pointer!");

    // first build difference of point and element center
    // and then dot product with unit normal at center
    const double dist = (mycnode->xspatial()[0] - Auxc()[0]) * Auxn()[0]
        + (mycnode->xspatial()[1] - Auxc()[1]) * Auxn()[1]
        + (mycnode->xspatial()[2] - Auxc()[2]) * Auxn()[2];

    // compute projection
    for (int k = 0; k < 3; ++k)
      vertices[k] = mycnode->xspatial()[k] - dist * Auxn()[k];

    // get node id, too
    snodeids[0] = mycnode->Id();

    // store into vertex data structure
    SlaveVertices().push_back(
        MORTAR::Vertex(
            vertices,
            MORTAR::Vertex::projslave,
            snodeids,
            NULL,
            NULL,
            false,
            false,
            NULL,
            -1.0));

    //std::cout << "->RealNode(S) " << mycnode->Id() << ": " << mycnode->xspatial()[0] << " " << mycnode->xspatial()[1] << " " << mycnode->xspatial()[2] << std::endl;
    //std::cout << "->ProjNode(S) " << mycnode->Id() << ": " << vertices[0] << " " << vertices[1] << " " << vertices[2] << std::endl;
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Linearization of slave vertex (3D) AuxPlane              farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineCoupling3d::SlaveVertexLinearization(
    std::vector<std::vector<GEN::pairedvector<int, double> > >& currlin)
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
  else if (dt == MORTAR::MortarElement::quad4
      || dt == MORTAR::MortarElement::quad8
      || dt == MORTAR::MortarElement::quad9)
  {
    scxi[0] = 0.0;
    scxi[1] = 0.0;
  }
  else
    dserror("ERROR: MasterVertexLinearization called for unknown element type");

  // evlauate shape functions + derivatives at scxi
  int nrow = SurfaceElement().NumNode();
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,2,true);
  SurfaceElement().EvaluateShape(scxi,sval,sderiv,nrow);

  // we need all participating slave nodes
  DRT::Node** snodes = SurfaceElement().Nodes();
  std::vector<MORTAR::MortarNode*> smrtrnodes(nrow);

  for (int i=0;i<nrow;++i)
  {
    smrtrnodes[i] = dynamic_cast<MORTAR::MortarNode*>(snodes[i]);
    if (!smrtrnodes[i]) dserror("ERROR: MasterVertexLinearization: Null pointer!");
  }

  // linearization of the SlaveIntEle spatial coords
  std::vector<std::vector<GEN::pairedvector<int, double> > > snodelin(0);

  // resize the linearizations
  snodelin.resize(nrow,std::vector<GEN::pairedvector<int,double> >(3,1));

  // loop over all intEle nodes
  for (int in=0; in<nrow; ++in)
    for (int dim=0; dim<3; ++dim)
      snodelin[in][dim][smrtrnodes[in]->Dofs()[dim]]+=1.;

  // map iterator
  typedef GEN::pairedvector<int, double>  :: const_iterator _CI;    // linearization of element center Auxc()
  std  ::vector<GEN::pairedvector<int  ,double> > linauxc(3,SurfaceElement().NumNode()); // assume 3 dofs per node

  for (int i = 0; i < nrow; ++i)
      for (int dim=0; dim<3; ++dim)
        for (_CI p=snodelin[i][dim].begin(); p!=snodelin[i][dim].end(); ++p)
          linauxc[dim][p->first] = sval[i]*p->second;

  // linearization of element normal Auxn()
  std::vector<GEN::pairedvector<int, double> >& linauxn = GetDerivAuxn();

  // linearization of the MasterIntEle spatial coords
  std::vector<std::vector<GEN::pairedvector<int, double> > > mnodelin(0);

  // resize the linearizations
  mnodelin.resize(LineElement()->NumNode(),std::vector<GEN::pairedvector<int,double> >(3,1));

  // loop over all intEle nodes
  for (int in=0; in<LineElement()->NumNode(); ++in)
  {
    MORTAR::MortarNode* mrtrmnode=dynamic_cast<MORTAR::MortarNode*>(idiscret_.gNode(LineElement()->NodeIds()[in]));
    if (mrtrmnode==NULL)
      dserror("dynamic cast to mortar node went wrong");

    for (int dim=0; dim<3; ++dim)
      mnodelin[in][dim][mrtrmnode->Dofs()[dim]]+=1.;
  }

  // put everything together for slave vertex linearization
  // loop over all vertices
  for (int i=0; i<LineElement()->NumNode(); ++i)
  {
    MORTAR::MortarNode* mrtrmnode=dynamic_cast<MORTAR::MortarNode*>(idiscret_.gNode(LineElement()->NodeIds()[i]));
    if (!mrtrmnode)
      dserror("cast to mortar node failed");

    // (1) slave node coordinates part
    for (_CI p=mnodelin[i][0].begin(); p!=mnodelin[i][0].end(); ++p)
    {
      currlin[i][0][p->first] += (1.0 - Auxn()[0] * Auxn()[0])*p->second;
      currlin[i][1][p->first] -= (Auxn()[0] * Auxn()[1])      *p->second;
      currlin[i][2][p->first] -= (Auxn()[0] * Auxn()[2])      *p->second;
    }
    for (_CI p=mnodelin[i][1].begin(); p!=mnodelin[i][1].end(); ++p)
    {
      currlin[i][0][p->first] -= (Auxn()[0] * Auxn()[1])      *p->second;
      currlin[i][1][p->first] += (1.0 - Auxn()[1] * Auxn()[1])*p->second;
      currlin[i][2][p->first] -= (Auxn()[1] * Auxn()[2])      *p->second;
    }
    for (_CI p=mnodelin[i][2].begin(); p!=mnodelin[i][2].end(); ++p)
    {
      currlin[i][0][p->first] -= (Auxn()[2] * Auxn()[0])      *p->second;
      currlin[i][1][p->first] -= (Auxn()[2] * Auxn()[1])      *p->second;
      currlin[i][2][p->first] += (1.0 - Auxn()[2] * Auxn()[2])*p->second;
    }

    // (2) slave element center coordinates (Auxc()) part
    for (_CI p = linauxc[0].begin(); p != linauxc[0].end(); ++p)
      for (int k = 0; k < 3; ++k)
        currlin[i][k][p->first] += Auxn()[0] * Auxn()[k] * (p->second);

    for (_CI p = linauxc[1].begin(); p != linauxc[1].end(); ++p)
      for (int k = 0; k < 3; ++k)
        currlin[i][k][p->first] += Auxn()[1] * Auxn()[k] * (p->second);

    for (_CI p = linauxc[2].begin(); p != linauxc[2].end(); ++p)
      for (int k = 0; k < 3; ++k)
        currlin[i][k][p->first] += Auxn()[2] * Auxn()[k] * (p->second);

    // (3) slave element normal (Auxn()) part
    double xdotn = (mrtrmnode->xspatial()[0] - Auxc()[0]) * Auxn()[0]
        + (mrtrmnode->xspatial()[1] - Auxc()[1]) * Auxn()[1]
        + (mrtrmnode->xspatial()[2] - Auxc()[2]) * Auxn()[2];

    for (_CI p = linauxn[0].begin(); p != linauxn[0].end(); ++p)
    {
      currlin[i][0][p->first] -= xdotn * (p->second);
      for (int k = 0; k < 3; ++k)
        currlin[i][k][p->first] -= (mrtrmnode->xspatial()[0] - Auxc()[0]) * Auxn()[k]
            * (p->second);
    }

    for (_CI p = linauxn[1].begin(); p != linauxn[1].end(); ++p)
    {
      currlin[i][1][p->first] -= xdotn * (p->second);
      for (int k = 0; k < 3; ++k)
        currlin[i][k][p->first] -= (mrtrmnode->xspatial()[1] - Auxc()[1]) * Auxn()[k]
            * (p->second);
    }

    for (_CI p = linauxn[2].begin(); p != linauxn[2].end(); ++p)
    {
      currlin[i][2][p->first] -= xdotn * (p->second);
      for (int k = 0; k < 3; ++k)
        currlin[i][k][p->first] -= (mrtrmnode->xspatial()[2] - Auxc()[2]) * Auxn()[k]
            * (p->second);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  eval (public)                                            farah 07/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineCoupling3d::ProjectMaster()
{
  // project master nodes onto auxiliary plane
  int nnodes = SurfaceElement().NumNode();
  DRT::Node** mynodes = SurfaceElement().Nodes();
  if (!mynodes)
    dserror("ERROR: ProjectMaster: Null pointer!");

  // initialize storage for master coords + their ids
  std::vector<double> vertices(3);
  std::vector<int> mnodeids(1);

  for (int i = 0; i < nnodes; ++i)
  {
    CoNode* mycnode = dynamic_cast<CoNode*>(mynodes[i]);
    if (!mycnode)
      dserror("ERROR: ProjectMaster: Null pointer!");

    // first build difference of point and element center
    // and then dot product with unit normal at center
    const double dist = (mycnode->xspatial()[0] - Auxc()[0]) * Auxn()[0]
        + (mycnode->xspatial()[1] - Auxc()[1]) * Auxn()[1]
        + (mycnode->xspatial()[2] - Auxc()[2]) * Auxn()[2];

    // compute projection
    for (int k = 0; k < 3; ++k)
      vertices[k] = mycnode->xspatial()[k] - dist * Auxn()[k];

    // get node id, too
    mnodeids[0] = mycnode->Id();

    // store into vertex data structure
    MasterVertices().push_back(
        MORTAR::Vertex(
            vertices,
            MORTAR::Vertex::master,
            mnodeids,
            NULL,
            NULL,
            false,
            false,
            NULL,
            -1.0));

    //std::cout << "->RealNode(M) " << mycnode->Id() << ": " << mycnode->xspatial()[0] << " " << mycnode->xspatial()[1] << " " << mycnode->xspatial()[2] << std::endl;
    //std::cout << "->ProjNode(M) " << mycnode->Id() << ": " << vertices[0] << " " << vertices[1] << " " << vertices[2] << std::endl;
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Linearization of slave vertex (3D) AuxPlane               farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineCoupling3d::MasterVertexLinearization(
    std::vector<std::vector<GEN::pairedvector<int, double> > >& currlin)
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
  else if (dt == MORTAR::MortarElement::quad4
      || dt == MORTAR::MortarElement::quad8
      || dt == MORTAR::MortarElement::quad9)
  {
    scxi[0] = 0.0;
    scxi[1] = 0.0;
  }
  else
    dserror("ERROR: SlaveVertexLinearization called for unknown element type");

  // evlauate shape functions + derivatives at scxi
  const int nrow = SurfaceElement().NumNode();
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,2,true);
  SurfaceElement().EvaluateShape(scxi,sval,sderiv,nrow);

  // we need all participating slave nodes
  DRT::Node** snodes = SurfaceElement().Nodes();
  std::vector<MORTAR::MortarNode*> smrtrnodes(nrow);

  for (int i=0;i<nrow;++i)
  {
    smrtrnodes[i] = dynamic_cast<MORTAR::MortarNode*>(snodes[i]);
    if (!smrtrnodes[i]) dserror("ERROR: SlaveVertexLinearization: Null pointer!");
  }

  // linearization of the IntEle spatial coords
  std::vector<std::vector<GEN::pairedvector<int, double> > > nodelin(0);

  // resize the linearizations
  nodelin.resize(nrow,std::vector<GEN::pairedvector<int,double> >(3,1));

  // loop over all intEle nodes
  for (int in=0; in<nrow; ++in)
    for (int dim=0; dim<3; ++dim)
      nodelin[in][dim][smrtrnodes[in]->Dofs()[dim]]+=1.;

  // map iterator
  typedef GEN::pairedvector<int, double>  :: const_iterator _CI;    // linearization of element center Auxc()
  std  ::vector<GEN::pairedvector<int  ,double> > linauxc(3,SurfaceElement().NumNode()); // assume 3 dofs per node

  for (int i = 0; i < nrow; ++i)
      for (int dim=0; dim<3; ++dim)
        for (_CI p=nodelin[i][dim].begin(); p!=nodelin[i][dim].end(); ++p)
          linauxc[dim][p->first] = sval[i]*p->second;

  // linearization of element normal Auxn()
  std::vector<GEN::pairedvector<int, double> >& linauxn = GetDerivAuxn();

  // put everything together for slave vertex linearization
  // loop over all vertices
  for (int i=0; i<SurfaceElement().NumNode(); ++i)
  {
    MORTAR::MortarNode* mrtrsnode=dynamic_cast<MORTAR::MortarNode*>(SurfaceElement().Nodes()[i]);
    if (!mrtrsnode)
      dserror("cast to mortar node failed");

    // (1) slave node coordinates part
    for (_CI p=nodelin[i][0].begin(); p!=nodelin[i][0].end(); ++p)
    {
      currlin[i][0][p->first] += (1.0 - Auxn()[0] * Auxn()[0])*p->second;
      currlin[i][1][p->first] -= (Auxn()[0] * Auxn()[1])      *p->second;
      currlin[i][2][p->first] -= (Auxn()[0] * Auxn()[2])      *p->second;
    }
    for (_CI p=nodelin[i][1].begin(); p!=nodelin[i][1].end(); ++p)
    {
      currlin[i][0][p->first] -= (Auxn()[0] * Auxn()[1])      *p->second;
      currlin[i][1][p->first] += (1.0 - Auxn()[1] * Auxn()[1])*p->second;
      currlin[i][2][p->first] -= (Auxn()[1] * Auxn()[2])      *p->second;
    }
    for (_CI p=nodelin[i][2].begin(); p!=nodelin[i][2].end(); ++p)
    {
      currlin[i][0][p->first] -= (Auxn()[2] * Auxn()[0])      *p->second;
      currlin[i][1][p->first] -= (Auxn()[2] * Auxn()[1])      *p->second;
      currlin[i][2][p->first] += (1.0 - Auxn()[2] * Auxn()[2])*p->second;
    }

    // (2) slave element center coordinates (Auxc()) part
    for (_CI p = linauxc[0].begin(); p != linauxc[0].end(); ++p)
      for (int k = 0; k < 3; ++k)
        currlin[i][k][p->first] += Auxn()[0] * Auxn()[k] * (p->second);

    for (_CI p = linauxc[1].begin(); p != linauxc[1].end(); ++p)
      for (int k = 0; k < 3; ++k)
        currlin[i][k][p->first] += Auxn()[1] * Auxn()[k] * (p->second);

    for (_CI p = linauxc[2].begin(); p != linauxc[2].end(); ++p)
      for (int k = 0; k < 3; ++k)
        currlin[i][k][p->first] += Auxn()[2] * Auxn()[k] * (p->second);

    // (3) slave element normal (Auxn()) part
    double xdotn = (mrtrsnode->xspatial()[0] - Auxc()[0]) * Auxn()[0]
                 + (mrtrsnode->xspatial()[1] - Auxc()[1]) * Auxn()[1]
                 + (mrtrsnode->xspatial()[2] - Auxc()[2]) * Auxn()[2];

    for (_CI p = linauxn[0].begin(); p != linauxn[0].end(); ++p)
    {
      currlin[i][0][p->first] -= xdotn * (p->second);
      for (int k = 0; k < 3; ++k)
        currlin[i][k][p->first] -= (mrtrsnode->xspatial()[0] - Auxc()[0]) * Auxn()[k]
                                  * (p->second);
    }

    for (_CI p = linauxn[1].begin(); p != linauxn[1].end(); ++p)
    {
      currlin[i][1][p->first] -= xdotn * (p->second);
      for (int k = 0; k < 3; ++k)
        currlin[i][k][p->first] -= (mrtrsnode->xspatial()[1] - Auxc()[1]) * Auxn()[k]
                                  * (p->second);
    }

    for (_CI p = linauxn[2].begin(); p != linauxn[2].end(); ++p)
    {
      currlin[i][2][p->first] -= xdotn * (p->second);
      for (int k = 0; k < 3; ++k)
        currlin[i][k][p->first] -= (mrtrsnode->xspatial()[2] - Auxc()[2]) * Auxn()[k]
                                  * (p->second);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  get communicator                                         farah 07/16|
 *----------------------------------------------------------------------*/
const Epetra_Comm& CONTACT::LineCoupling3d::Comm() const
{
  return idiscret_.Comm();
}


/*----------------------------------------------------------------------*
 |  ctor for ltl (public)                                    farah 07/16|
 *----------------------------------------------------------------------*/
CONTACT::LineToLineCoupling3d::LineToLineCoupling3d(
      DRT::Discretization& idiscret,
      int dim,
      Teuchos::ParameterList& params,
      Teuchos::RCP<MORTAR::MortarElement>& lsele,
      Teuchos::RCP<MORTAR::MortarElement>& lmele) :
  idiscret_(idiscret),
  dim_(dim),
  imortar_(params),
  lSele_(lsele),
  lMele_(lmele)
{
  // empty constructor
}


/*----------------------------------------------------------------------*
 |  get communicator                                         farah 07/16|
 *----------------------------------------------------------------------*/
const Epetra_Comm& CONTACT::LineToLineCoupling3d::Comm() const
{
  return idiscret_.Comm();
}

/*----------------------------------------------------------------------*
 |  eval                                                     farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineToLineCoupling3d::EvaluateCoupling()
{
  // 1. check parallelity
  bool parallel = CheckParallelity();
  if(parallel)
    return;

  // 2. calc intersection
  // create empty points
  double sxi[1] = {0.0};
  double mxi[1] = {0.0};

  // create empty lin vectors
  GEN::pairedvector<int,double> dsxi(3*LineMasterElement()->NumNode() + 3*LineSlaveElement()->NumNode());
  GEN::pairedvector<int,double> dmxi(3*LineMasterElement()->NumNode() + 3*LineSlaveElement()->NumNode());
  LineIntersection(sxi, mxi, dsxi, dmxi);

  // 3. check solution
  bool valid = CheckIntersection(sxi,mxi);
  if(!valid)
    return;

  // 4. evaluate terms
  EvaluateTerms(sxi, mxi, dsxi, dmxi);

  return;
}


/*----------------------------------------------------------------------*
 |  line-line intersection                                   farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineToLineCoupling3d::EvaluateTerms(
    double* sxi,
    double* mxi,
    GEN::pairedvector<int,double>& dsxi,
    GEN::pairedvector<int,double>& dmxi)
{
  INPAR::MORTAR::ShapeFcn shape =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(imortar_,"LM_SHAPEFCN");

  double jac = 1.0;
  double wgt = 1.0;

  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = LineSlaveElement()->Nodes();
  if(!mynodes)
    dserror("ERROR: IntegrateDerivCell3DAuxPlaneLTS: Null pointer!");
  DRT::Node** mnodes = LineMasterElement()->Nodes();
  if(!mnodes)
    dserror("ERROR: IntegrateDerivCell3DAuxPlaneLTS: Null pointer!");

  int nnodes = 2;
  int ndof = 3;
  int nrow = LineSlaveElement()->NumNode();
  int ncol = LineMasterElement()->NumNode();
  int nlm = 1;


  // slave values
  LINALG::SerialDenseVector sval(nnodes);
  LINALG::SerialDenseMatrix sderiv(nnodes,1);
  LINALG::SerialDenseVector lmval(nlm);
  LINALG::SerialDenseMatrix lmderiv(nlm,1);
  LineSlaveElement()->EvaluateShape(sxi, sval, sderiv, nnodes);
  LineSlaveElement()->EvaluateShape(sxi, lmval, lmderiv, nnodes);

  lmval(0)     = 1.0;
  lmderiv(0,0) = 0.0;

  // master values
  LINALG::SerialDenseVector mval(nnodes);
  LINALG::SerialDenseMatrix mderiv(nnodes,1);
  LineMasterElement()->EvaluateShape(mxi, mval, mderiv, nnodes);

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator _CI;

  int linsize = 0;
  for (int i=0;i<nrow;++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*> (mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  //**********************************************************************
  // geometric quantities
  //**********************************************************************
  double gpn[3]      = {0.0,0.0,0.0};
  GEN::pairedvector<int,double> dgapgp((ncol*ndof)+10*linsize);                   // gap lin. without lm and jac.
  double gap[1]      = {0.0};
  std::vector<GEN::pairedvector<int,double> > dnmap_unit(3,10*linsize);           // deriv of x,y and z comp. of gpn (unit)

  //**********************************************************************
  // evaluate at GP and lin char. quantities
 //**********************************************************************
  double sgpx[3] = {0.0, 0.0, 0.0};
  double mgpx[3] = {0.0, 0.0, 0.0};

  for (int i=0;i<nrow;++i)
  {
    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[i]);
    gpn[0]+=sval[i]*mymrtrnode->MoData().n()[0];
    gpn[1]+=sval[i]*mymrtrnode->MoData().n()[1];
    gpn[2]+=sval[i]*mymrtrnode->MoData().n()[2];

    sgpx[0]+=sval[i] * LineSlaveElement()->GetNodalCoords(0,i);
    sgpx[1]+=sval[i] * LineSlaveElement()->GetNodalCoords(1,i);
    sgpx[2]+=sval[i] * LineSlaveElement()->GetNodalCoords(2,i);
  }

  // build interpolation of master GP coordinates
  for (int i=0;i<ncol;++i)
  {
    mgpx[0]+=mval[i]*LineMasterElement()->GetNodalCoords(0,i);
    mgpx[1]+=mval[i]*LineMasterElement()->GetNodalCoords(1,i);
    mgpx[2]+=mval[i]*LineMasterElement()->GetNodalCoords(2,i);
  }

  // normalize interpolated GP normal back to length 1.0 !!!
    double lengthn = sqrt(gpn[0]*gpn[0]+gpn[1]*gpn[1]+gpn[2]*gpn[2]);
  if (lengthn<1.0e-12) dserror("ERROR: IntegrateAndDerivSegment: Divide by zero!");

  for (int i=0;i<3;++i)
    gpn[i]/=lengthn;

  // build gap function at current GP
  for (int i=0;i<Dim();++i)
    gap[0]+=(mgpx[i]-sgpx[i])*gpn[i];

  // build directional derivative of slave GP normal (non-unit)
  GEN::pairedvector<int,double> dmap_nxsl_gp(linsize);
  GEN::pairedvector<int,double> dmap_nysl_gp(linsize);
  GEN::pairedvector<int,double> dmap_nzsl_gp(linsize);

  for (int i=0;i<nrow;++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*> (mynodes[i]);

    GEN::pairedvector<int,double>& dmap_nxsl_i = cnode->CoData().GetDerivN()[0];
    GEN::pairedvector<int,double>& dmap_nysl_i = cnode->CoData().GetDerivN()[1];
    GEN::pairedvector<int,double>& dmap_nzsl_i = cnode->CoData().GetDerivN()[2];

    for (_CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
      dmap_nxsl_gp[p->first] += sval[i]*(p->second);
    for (_CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
      dmap_nysl_gp[p->first] += sval[i]*(p->second);
    for (_CI p=dmap_nzsl_i.begin();p!=dmap_nzsl_i.end();++p)
      dmap_nzsl_gp[p->first] += sval[i]*(p->second);

    for (_CI p=dsxi.begin();p!=dsxi.end();++p)
    {
      double valx =  sderiv(i,0)*cnode->MoData().n()[0];
      dmap_nxsl_gp[p->first] += valx*(p->second);
      double valy =  sderiv(i,0)*cnode->MoData().n()[1];
      dmap_nysl_gp[p->first] += valy*(p->second);
      double valz =  sderiv(i,0)*cnode->MoData().n()[2];
      dmap_nzsl_gp[p->first] += valz*(p->second);
    }
  }

  const double ll     = lengthn*lengthn;
  const double linv   = 1.0/(lengthn);
  const double lllinv = 1.0/(lengthn*lengthn*lengthn);
  const double sxsx   = gpn[0]*gpn[0]*ll;
  const double sxsy   = gpn[0]*gpn[1]*ll;
  const double sxsz   = gpn[0]*gpn[2]*ll;
  const double sysy   = gpn[1]*gpn[1]*ll;
  const double sysz   = gpn[1]*gpn[2]*ll;
  const double szsz   = gpn[2]*gpn[2]*ll;

  for (_CI p=dmap_nxsl_gp.begin();p!=dmap_nxsl_gp.end();++p)
  {
    dnmap_unit[0][p->first] += linv*(p->second);
    dnmap_unit[0][p->first] -= lllinv*sxsx*(p->second);
    dnmap_unit[1][p->first] -= lllinv*sxsy*(p->second);
    dnmap_unit[2][p->first] -= lllinv*sxsz*(p->second);
  }

  for (_CI p=dmap_nysl_gp.begin();p!=dmap_nysl_gp.end();++p)
  {
    dnmap_unit[1][p->first] += linv*(p->second);
    dnmap_unit[1][p->first] -= lllinv*sysy*(p->second);
    dnmap_unit[0][p->first] -= lllinv*sxsy*(p->second);
    dnmap_unit[2][p->first] -= lllinv*sysz*(p->second);
  }

  for (_CI p=dmap_nzsl_gp.begin();p!=dmap_nzsl_gp.end();++p)
  {
    dnmap_unit[2][p->first] += linv*(p->second);
    dnmap_unit[2][p->first] -= lllinv*szsz*(p->second);
    dnmap_unit[0][p->first] -= lllinv*sxsz*(p->second);
    dnmap_unit[1][p->first] -= lllinv*sysz*(p->second);
  }

  // add everything to dgapgp
  for (_CI p=dnmap_unit[0].begin();p!=dnmap_unit[0].end();++p)
    dgapgp[p->first] += (mgpx[0]-sgpx[0]) * (p->second);

  for (_CI p=dnmap_unit[1].begin();p!=dnmap_unit[1].end();++p)
    dgapgp[p->first] += (mgpx[1]-sgpx[1]) * (p->second);

  for (_CI p=dnmap_unit[2].begin();p!=dnmap_unit[2].end();++p)
    dgapgp[p->first] += (mgpx[2]-sgpx[2]) *(p->second);

  // lin slave nodes
  for (int z=0;z<nrow;++z)
  {
    CoNode* cnode = dynamic_cast<CoNode*> (mynodes[z]);
    for (int k=0;k<3;++k)
      dgapgp[cnode->Dofs()[k]] -= sval[z] * gpn[k];
  }

  for (_CI p=dsxi.begin();p!=dsxi.end();++p)
  {
    double& dg = dgapgp[p->first] ;
    const double& ps = p->second;
    for (int z=0;z<nrow;++z)
    {
      CoNode* cnode = dynamic_cast<CoNode*> (mynodes[z]);
      for (int k=0;k<3;++k)
        dg -= gpn[k] * sderiv(z,0) * cnode->xspatial()[k] * ps;
    }
  }

  //        MASTER
  // lin master nodes
  for (int z=0;z<ncol;++z)
  {
    CoNode* cnode = dynamic_cast<CoNode*> (mnodes[z]);
    for (int k=0;k<3;++k)
      dgapgp[cnode->Dofs()[k]] += mval[z] * gpn[k];
  }

  for (_CI p=dmxi.begin();p!=dmxi.end();++p)
  {
    double& dg = dgapgp[p->first] ;
    const double& ps = p->second;
    for (int z=0;z<ncol;++z)
    {
      CoNode* cnode = dynamic_cast<CoNode*> (mnodes[z]);
      for (int k=0;k<3;++k)
        dg+=gpn[k] * mderiv(z,0) * cnode->xspatial()[k] * ps;
    }
  }



  // weighted gap
  for (int j=0;j<nlm;++j)
  {
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(mynodes[j]);

    double prod = 0.0;
    // Petrov-Galerkin approach (dual LM for D/M but standard LM for gap)
    if (shape == INPAR::MORTAR::shape_petrovgalerkin)
      prod = sval[j]*gap[0]*jac*wgt;
    // usual standard or dual LM approach
    else
      prod = lmval[j]*gap[0]*jac*wgt;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (cnode->IsOnBound())
      continue;

    // add current Gauss point's contribution to gseg
    cnode->AddltlGapValue(prod);
  }

  for (int iter=0;iter<nlm;++iter)
  {
    // map iterator
    typedef GEN::pairedvector<int,double>::const_iterator _CI;

    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[iter]);
    if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

    static double fac = 0.0;

    // get the corresponding map as a reference
    std::map<int,double>& dgmap = dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivGltl();

    // switch if Petrov-Galerkin approach for LM is applied
    if (shape == INPAR::MORTAR::shape_petrovgalerkin)
    {
      // (1) Lin(Phi) - does not exist here for Petrov-Galerkin approach

      // (2) Lin(N) - slave GP coordinates
      fac = wgt*sderiv(iter,0)*gap[0]*jac;
      for (_CI p=dsxi.begin();p!=dsxi.end();++p)
        dgmap[p->first] += fac*(p->second);


      // (3) Lin(g) - gap function
      fac = wgt*sval[iter]*jac;
      for (_CI p=dgapgp.begin();p!=dgapgp.end();++p)
        dgmap[p->first] += fac*(p->second);

      // (4) Lin(dsxideta) - intcell GP Jacobian
//      fac = wgt*sval[iter]*gap[0];
//      for (_CI p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
//        dgmap[p->first] += fac*(p->second);
    }

    // the usual standard or dual LM approach
    else
    {
      // (1) Lin(Phi) - dual shape functions
      // only line2 --> not required!

      // (2) Lin(Phi) - slave GP coordinates
        fac = wgt*lmderiv(iter,0)*gap[0]*jac;
        for (_CI p=dsxi.begin();p!=dsxi.end();++p)
          dgmap[p->first] += fac*(p->second);

      // (3) Lin(g) - gap function
      fac = wgt*lmval[iter]*jac;
      for (_CI p=dgapgp.begin();p!=dgapgp.end();++p)
        dgmap[p->first] += fac*(p->second);

      // (4) Lin(dsxideta) - intcell GP Jacobian
//      fac = wgt*lmval[iter]*gap[0];
//      for (_CI p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
//        dgmap[p->first] += fac*(p->second);
    }
  }

  // integrate D and M matrix
  // dual shape functions
  if (shape == INPAR::MORTAR::shape_dual ||
      shape == INPAR::MORTAR::shape_petrovgalerkin)
  {
    for (int j=0;j<nlm;++j)
    {
      CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(mynodes[j]);

      // integrate mseg
      for (int k=0; k<ncol; ++k)
      {
        CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mnodes[k]);

        // multiply the two shape functions
        double prod = lmval[j]*mval[k]*jac*wgt;

        if(abs(prod)>MORTARINTTOL) cnode->AddDltlValue(cnode->Id(),prod);
        if(abs(prod)>MORTARINTTOL) cnode->AddSNode(cnode->Id()); // only for friction!

        if(abs(prod)>MORTARINTTOL) cnode->AddMltsValue(mnode->Id(),prod);
        if(abs(prod)>MORTARINTTOL) cnode->AddMNode(mnode->Id());  // only for friction!
      }
    }
  }
  else
  {
    for (int j=0;j<nlm;++j)
    {
      CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(mynodes[j]);

      // integrate dseg
      for (int k=0; k<nrow; ++k)
      {
        CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mynodes[k]);

        // multiply the two shape functions
        double prod = lmval[j]*sval[k]*jac*wgt;

        if(abs(prod)>MORTARINTTOL) cnode->AddDltlValue(mnode->Id(),prod);
        if(abs(prod)>MORTARINTTOL) cnode->AddSNode(mnode->Id()); // only for friction!
      }

      // integrate mseg
      for (int k=0; k<ncol; ++k)
      {
        CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mnodes[k]);

        // multiply the two shape functions
        double prod = lmval[j]*mval[k]*jac*wgt;

        if(abs(prod)>MORTARINTTOL) cnode->AddMltlValue(mnode->Id(),prod);
        if(abs(prod)>MORTARINTTOL) cnode->AddMNode(mnode->Id());  // only for friction!
      }
    }
  }

  if (shape == INPAR::MORTAR::shape_dual ||
      shape == INPAR::MORTAR::shape_petrovgalerkin)
  {
    // integrate LinD
    for (int j=0; j<nlm;++j)
    {
      MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[j]);
      if (!mymrtrnode)
        dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

      int sgid = mymrtrnode->Id();
      std::map<int,double>& ddmap_jk = dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivDltl()[sgid];


      // integrate LinM
      for (int k=0; k<ncol; ++k)
      {
        // global master node ID
        int mgid = LineMasterElement()->Nodes()[k]->Id();
        static double fac = 0.0;

        // get the correct map as a reference
        std::map<int,double>& dmmap_jk = dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivMltl()[mgid];

        // (1) Lin(Phi) - dual shape functions
        // this vanishes here since there are no deformation-dependent dual functions

        // (2) Lin(NSlave) - slave GP coordinates
        fac = wgt*lmderiv(j, 0)*mval[k]*jac;
        for (_CI p=dsxi.begin(); p!=dsxi.end(); ++p)
        {
          dmmap_jk[p->first] += fac*(p->second);
          ddmap_jk[p->first] += fac*(p->second);
        }

//        fac = wgt*lmderiv(j, 1)*mval[k]*jac;
//        for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
//        {
//          dmmap_jk[p->first] += fac*(p->second);
//          ddmap_jk[p->first] += fac*(p->second);
//        }

        // (3) Lin(NMaster) - master GP coordinates
        fac = wgt*lmval[j]*mderiv(k, 0)*jac;
        for (_CI p=dmxi.begin(); p!=dmxi.end(); ++p)
        {
          dmmap_jk[p->first] += fac*(p->second);
          ddmap_jk[p->first] += fac*(p->second);
        }

//        fac = wgt*lmval[j]*mderiv(k, 1)*jac;
//        for (_CI p=dmxigp[1].begin(); p!=dmxigp[1].end(); ++p)
//        {
//          dmmap_jk[p->first] += fac*(p->second);
//          ddmap_jk[p->first] += fac*(p->second);
//        }

        // (4) Lin(dsxideta) - intcell GP Jacobian
//        fac = wgt*lmval[j]*mval[k];
//        for (_CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
//        {
//          dmmap_jk[p->first] += fac*(p->second);
//          ddmap_jk[p->first] += fac*(p->second);
//        }
      } // loop over master nodes
    }
  }
  else // standard case
  {
    // integrate LinD
    for (int j=0; j<nlm;++j)
    {
      MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[j]);
      if (!mymrtrnode)
        dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");


      // integrate LinD
      for (int k=0; k<nrow; ++k)
      {
        // global master node ID
        int mgid = LineSlaveElement()->Nodes()[k]->Id();
        static double fac = 0.0;

        // get the correct map as a reference
        std::map<int,double>& ddmap_jk = dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivDltl()[mgid];

        // (1) Lin(Phi) - dual shape functions
        // this vanishes here since there are no deformation-dependent dual functions

        // (2) Lin(NSlave) - slave GP coordinates
        fac = wgt*lmderiv(j, 0)*sval[k]*jac;
        for (_CI p=dsxi.begin(); p!=dsxi.end(); ++p)
        {
          ddmap_jk[p->first] += fac*(p->second);
        }

//        fac = wgt*lmderiv(j, 1)*sval[k]*jac;
//        for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
//        {
//          ddmap_jk[p->first] += fac*(p->second);
//        }

        // (3) Lin(NMaster) - master GP coordinates
        fac = wgt*lmval[j]*sderiv(k, 0)*jac;
        for (_CI p=dsxi.begin(); p!=dsxi.end(); ++p)
        {
          ddmap_jk[p->first] += fac*(p->second);
        }

//        fac = wgt*lmval[j]*sderiv(k, 1)*jac;
//        for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
//        {
//          ddmap_jk[p->first] += fac*(p->second);
//        }

//        fac = wgt*lmval[j]*mderiv(k, 1)*jac;
//        for (_CI p=dmxigp[1].begin(); p!=dmxigp[1].end(); ++p)
//        {
//          dmmap_jk[p->first] += fac*(p->second);
//          ddmap_jk[p->first] += fac*(p->second);
//        }

        // (4) Lin(dsxideta) - intcell GP Jacobian
//        fac = wgt*lmval[j]*sval[k];
//        for (_CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
//        {
//          ddmap_jk[p->first] += fac*(p->second);
//        }
      } // loop over slave nodes

      // integrate LinM
      for (int k=0; k<ncol; ++k)
      {
        // global master node ID
        int mgid = LineMasterElement()->Nodes()[k]->Id();
        static double fac = 0.0;

        // get the correct map as a reference
        std::map<int,double>& dmmap_jk = dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivMltl()[mgid];

        // (1) Lin(Phi) - dual shape functions
        // this vanishes here since there are no deformation-dependent dual functions

        // (2) Lin(NSlave) - slave GP coordinates
        fac = wgt*lmderiv(j, 0)*mval[k]*jac;
        for (_CI p=dsxi.begin(); p!=dsxi.end(); ++p)
        {
          dmmap_jk[p->first] += fac*(p->second);
        }

//        fac = wgt*lmderiv(j, 1)*mval[k]*jac;
//        for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
//        {
//          dmmap_jk[p->first] += fac*(p->second);
//        }

        // (3) Lin(NMaster) - master GP coordinates
        fac = wgt*lmval[j]*mderiv(k, 0)*jac;
        for (_CI p=dmxi.begin(); p!=dmxi.end(); ++p)
        {
          dmmap_jk[p->first] += fac*(p->second);
        }

//        fac = wgt*lmval[j]*mderiv(k, 1)*jac;
//        for (_CI p=dmxigp[1].begin(); p!=dmxigp[1].end(); ++p)
//        {
//          dmmap_jk[p->first] += fac*(p->second);
//          ddmap_jk[p->first] += fac*(p->second);
//        }

        // (4) Lin(dsxideta) - intcell GP Jacobian
//        fac = wgt*lmval[j]*mval[k];
//        for (_CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
//        {
//          dmmap_jk[p->first] += fac*(p->second);
//        }
      } // loop over master nodes
    }
  }


  return;
}

/*----------------------------------------------------------------------*
 |  line-line intersection                                   farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::LineToLineCoupling3d::LineIntersection(
    double* sxi,
    double* mxi,
    GEN::pairedvector<int,double>& dsxi,
    GEN::pairedvector<int,double>& dmxi)
{
  // only for line 2
  const int nnodes = 2;

  // calculate tangents:
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

  // res norm
  double conv = 0.0;

  // start in the element center
  double xiS[1] =  {0.0}; // xi_slave
  double xiM[1] =  {0.0}; // xi_master

  // function f (vector-valued)
  double f[2] =  { 0.0, 0.0 };

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
    LINALG::SerialDenseMatrix sderiv(nnodes,1);
    LineSlaveElement()->EvaluateShape(xiS, sval, sderiv, nnodes);

    // master values
    LINALG::SerialDenseVector mval(nnodes);
    LINALG::SerialDenseMatrix mderiv(nnodes,1);
    LineMasterElement()->EvaluateShape(xiM, mval, mderiv, nnodes);

    double xs[3] = {0.0, 0.0, 0.0};
    double xm[3] = {0.0, 0.0, 0.0};

    for(int i = 0; i < 3; ++i)
      xs[i] += sval[0]*ns1->xspatial()[i] + sval[1]*ns2->xspatial()[i];
    for(int i = 0; i < 3; ++i)
      xm[i] += mval[0]*nm1->xspatial()[i] + mval[1]*nm2->xspatial()[i];

    double xdiff[3] = {0.0, 0.0, 0.0};
    for(int i = 0; i < 3; ++i)
      xdiff[i] = xs[i] - xm[i];

    f[0] = xdiff[0]*vs[0] + xdiff[1]*vs[1] + xdiff[2]*vs[2];
    f[1] = xdiff[0]*vm[0] + xdiff[1]*vm[1] + xdiff[2]*vm[2];

    // check for convergence
    conv = sqrt(f[0] * f[0] + f[1] * f[1]);
    if (conv <= MORTARCONVTOL)
      break;

    //**********************************************
    //   F GRADIENT CALCULATION                   //
    //**********************************************

    double xsderiv[3] = {0.0, 0.0, 0.0};
    double xmderiv[3] = {0.0, 0.0, 0.0};
    for(int i = 0; i < 3; ++i)
      xsderiv[i] += sderiv(0,0)*ns1->xspatial()[i] + sderiv(1,0)*ns2->xspatial()[i];
    for(int i = 0; i < 3; ++i)
      xmderiv[i] += mderiv(0,0)*nm1->xspatial()[i] + mderiv(1,0)*nm2->xspatial()[i];


    df(0,0) = xsderiv[0] * vs[0] + xsderiv[1] * vs[1] + xsderiv[2] * vs[2];
    df(1,0) = xsderiv[0] * vm[0] + xsderiv[1] * vm[1] + xsderiv[2] * vm[2];

    df(0,1) = -xmderiv[0] * vs[0] - xmderiv[1] * vs[1] - xmderiv[2] * vs[2];
    df(1,1) = -xmderiv[0] * vm[0] - xmderiv[1] * vm[1] - xmderiv[2] * vm[2];

    //**********************************************
    //   solve deta = - inv(dF) * F               //
    //**********************************************
    double jacdet = df.Invert();
    if (abs(jacdet) < 1.0e-12)
      dserror("ERROR: Singular Jacobian for projection");

    // update eta and alpha
    xiS[0] += -df(0, 0) * f[0] - df(0, 1) * f[1];
    xiM[0] += -df(1, 0) * f[0] - df(1, 1) * f[1];
  }

  // Newton iteration unconverged
  if (conv > MORTARCONVTOL)
    dserror("ERROR: LTL intersection not converged!");


  //**********************************************
  //  Linearization                             //
  //**********************************************
  // prepare linearizations
  typedef GEN::pairedvector<int,double>::const_iterator _CI;

  // slave values
  LINALG::SerialDenseVector sval(nnodes);
  LINALG::SerialDenseMatrix sderiv(nnodes,1);
  LineSlaveElement()->EvaluateShape(xiS, sval, sderiv, nnodes);

  // master values
  LINALG::SerialDenseVector mval(nnodes);
  LINALG::SerialDenseMatrix mderiv(nnodes,1);
  LineMasterElement()->EvaluateShape(xiM, mval, mderiv, nnodes);

  double xs[3] = {0.0, 0.0, 0.0};
  double xm[3] = {0.0, 0.0, 0.0};

  for(int i = 0; i < 3; ++i)
    xs[i] += sval[0]*ns1->xspatial()[i] + sval[1]*ns2->xspatial()[i];
  for(int i = 0; i < 3; ++i)
    xm[i] += mval[0]*nm1->xspatial()[i] + mval[1]*nm2->xspatial()[i];

  double xdiff[3] = {0.0, 0.0, 0.0};
  for(int i = 0; i < 3; ++i)
    xdiff[i] = xs[i] - xm[i];

  std::vector<GEN::pairedvector<int,double> > xLin(3,1000);
  std::vector<GEN::pairedvector<int,double> > vsLin(3,1000);
  std::vector<GEN::pairedvector<int,double> > vmLin(3,1000);

  // global position difference
  for(int i=0; i<3;++i)
    (xLin[i])[ns1->Dofs()[i]] += sval(0);
  for(int i=0; i<3;++i)
    (xLin[i])[ns2->Dofs()[i]] += sval(1);

  for(int i=0; i<3;++i)
    (xLin[i])[nm1->Dofs()[i]] -= mval(0);
  for(int i=0; i<3;++i)
    (xLin[i])[nm2->Dofs()[i]] -= mval(1);

  // tangent vector slave
  for(int i=0; i<3;++i)
    (vsLin[i])[ns1->Dofs()[i]] += 1;
  for(int i=0; i<3;++i)
    (vsLin[i])[ns2->Dofs()[i]] -= 1;

  // tangent vector master
  for(int i=0; i<3;++i)
    (vmLin[i])[nm1->Dofs()[i]] += 1;
  for(int i=0; i<3;++i)
    (vmLin[i])[nm2->Dofs()[i]] -= 1;


  GEN::pairedvector<int,double> f0(1000);
  GEN::pairedvector<int,double> f1(1000);

  // lin xdiff * tangent + xdiff * lin tangent
  for (_CI p=xLin[0].begin();p!=xLin[0].end();++p)
    f0[p->first] += (p->second) * vs[0];
  for (_CI p=xLin[1].begin();p!=xLin[1].end();++p)
    f0[p->first] += (p->second) * vs[1];
  for (_CI p=xLin[2].begin();p!=xLin[2].end();++p)
    f0[p->first] += (p->second) * vs[2];

  for (_CI p=vsLin[0].begin();p!=vsLin[0].end();++p)
    f0[p->first] += (p->second) * xdiff[0];
  for (_CI p=vsLin[1].begin();p!=vsLin[1].end();++p)
    f0[p->first] += (p->second) * xdiff[1];
  for (_CI p=vsLin[2].begin();p!=vsLin[2].end();++p)
    f0[p->first] += (p->second) * xdiff[2];

  // lin xdiff * tangent + xdiff * lin tangent
  for (_CI p=xLin[0].begin();p!=xLin[0].end();++p)
    f1[p->first] += (p->second) * vm[0];
  for (_CI p=xLin[1].begin();p!=xLin[1].end();++p)
    f1[p->first] += (p->second) * vm[1];
  for (_CI p=xLin[2].begin();p!=xLin[2].end();++p)
    f1[p->first] += (p->second) * vm[2];

  for (_CI p=vmLin[0].begin();p!=vmLin[0].end();++p)
    f1[p->first] += (p->second) * xdiff[0];
  for (_CI p=vmLin[1].begin();p!=vmLin[1].end();++p)
    f1[p->first] += (p->second) * xdiff[1];
  for (_CI p=vmLin[2].begin();p!=vmLin[2].end();++p)
    f1[p->first] += (p->second) * xdiff[2];

  // end
  for (_CI p=f0.begin();p!=f0.end();++p)
    dsxi[p->first] -= (p->second) * df(0,0);
  for (_CI p=f1.begin();p!=f1.end();++p)
    dsxi[p->first] -= (p->second) * df(0,1);

  for (_CI p=f0.begin();p!=f0.end();++p)
    dmxi[p->first] -= (p->second) * df(1,0);
  for (_CI p=f1.begin();p!=f1.end();++p)
    dmxi[p->first] -= (p->second) * df(1,1);

  sxi[0] = xiS[0];
  mxi[0] = xiM[0];

  return;
}

/*----------------------------------------------------------------------*
 |  check if intersection is in para space interval          farah 07/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineToLineCoupling3d::CheckIntersection(double* sxi,double* mxi)
{
  if(sxi[0]>=-1.0-1e-12 and sxi[0]<=1.0+1e-12 and
      mxi[0]>=-1.0-1e-12 and mxi[0]<=1.0+1e-12 )
    return true;

  return false;
}

/*----------------------------------------------------------------------*
 |  check if line eles parallel                              farah 07/16|
 *----------------------------------------------------------------------*/
bool CONTACT::LineToLineCoupling3d::CheckParallelity()
{
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
  const double diff = scaprod - (lengthS * lengthM);
  if(abs(diff)<1e-12)
    return true;

  return false;
}













