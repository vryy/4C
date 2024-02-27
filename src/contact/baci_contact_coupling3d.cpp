/*---------------------------------------------------------------------*/
/*! \file
\brief Classes for mortar contact coupling in 3D.

\level 2


*/
/*---------------------------------------------------------------------*/

#include "baci_contact_coupling3d.hpp"

#include "baci_contact_defines.hpp"
#include "baci_contact_element.hpp"
#include "baci_contact_integrator.hpp"
#include "baci_contact_integrator_factory.hpp"
#include "baci_contact_interpolator.hpp"
#include "baci_contact_node.hpp"
#include "baci_inpar_contact.hpp"
#include "baci_lib_discret.hpp"
#include "baci_linalg_serialdensematrix.hpp"
#include "baci_linalg_serialdensevector.hpp"
#include "baci_linalg_utils_densematrix_inverse.hpp"
#include "baci_linalg_utils_densematrix_multiply.hpp"
#include "baci_mortar_calc_utils.hpp"
#include "baci_mortar_coupling3d_classes.hpp"
#include "baci_mortar_defines.hpp"
#include "baci_mortar_projector.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::Coupling3d::Coupling3d(DRT::Discretization& idiscret, int dim, bool quad,
    Teuchos::ParameterList& params, MORTAR::Element& sele, MORTAR::Element& mele)
    : MORTAR::Coupling3d(idiscret, dim, quad, params, sele, mele),
      stype_(CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(params, "STRATEGY"))
{
  // empty constructor

  return;
}

/*----------------------------------------------------------------------*
 |  Build auxiliary plane from slave element (public)         popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::AuxiliaryPlane()
{
  // we first need the element center:
  // for quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  double loccenter[2];

  CORE::FE::CellType dt = SlaveIntElement().Shape();
  if (dt == CORE::FE::CellType::tri3 || dt == CORE::FE::CellType::tri6)
  {
    loccenter[0] = 1.0 / 3.0;
    loccenter[1] = 1.0 / 3.0;
  }
  else if (dt == CORE::FE::CellType::quad4 || dt == CORE::FE::CellType::quad8 ||
           dt == CORE::FE::CellType::quad9)
  {
    loccenter[0] = 0.0;
    loccenter[1] = 0.0;
  }
  else
    dserror("AuxiliaryPlane called for unknown element type");

  // compute element center via shape fct. interpolation
  SlaveIntElement().LocalToGlobal(loccenter, Auxc(), 0);

  // we then compute the unit normal vector at the element center
  Lauxn() = SlaveIntElement().ComputeUnitNormalAtXi(loccenter, Auxn());

  // THIS IS CONTACT-SPECIFIC!
  // also compute linearization of the unit normal vector
  SlaveIntElement().DerivUnitNormalAtXi(loccenter, GetDerivAuxn());

  // std::cout << "Slave Element: " << SlaveIntElement().Id() << std::endl;
  // std::cout << "->Center: " << Auxc()[0] << " " << Auxc()[1] << " " << Auxc()[2] << std::endl;
  // std::cout << "->Normal: " << Auxn()[0] << " " << Auxn()[1] << " " << Auxn()[2] << std::endl;

  return true;
}

/*----------------------------------------------------------------------*
 |  Integration of cells (3D)                                 popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::IntegrateCells(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  /**********************************************************************/
  /* INTEGRATION                                                        */
  /* Integrate the Mortar matrix M and the weighted gap function g~ on  */
  /* the current integration cell of the slave / master element pair    */
  /**********************************************************************/

  static const INPAR::MORTAR::AlgorithmType algo =
      CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(imortar_, "ALGORITHM");

  // do nothing if there are no cells
  if (Cells().size() == 0) return false;

  // create a CONTACT integrator instance with correct NumGP and Dim
  // it is sufficient to do this once as all IntCells are triangles
  Teuchos::RCP<CONTACT::Integrator> integrator =
      CONTACT::INTEGRATOR::BuildIntegrator(stype_, imortar_, Cells()[0]->Shape(), Comm());
  // loop over all integration cells
  for (int i = 0; i < (int)(Cells().size()); ++i)
  {
    // integrate cell only if it has a non-zero area
    if (Cells()[i]->Area() < MORTARINTLIM * SlaveElementArea()) continue;

      // debug output of integration cells in GMSH
#ifdef MORTARGMSHCELLS
    GmshOutputCells(i);
#endif  // #ifdef MORTARGMSHCELLS
    // set segmentation status of all slave nodes
    // (hassegment_ of a slave node is true if ANY segment/cell
    // is integrated that contributes to this slave node)
    int nnodes = SlaveIntElement().NumNode();
    DRT::Node** mynodes = SlaveIntElement().Nodes();
    if (!mynodes) dserror("Null pointer!");
    for (int k = 0; k < nnodes; ++k)
    {
      MORTAR::Node* mycnode = dynamic_cast<MORTAR::Node*>(mynodes[k]);
      if (!mycnode) dserror("Null pointer!");
      mycnode->HasSegment() = true;
    }

    // *******************************************************************
    // different options for mortar integration
    // *******************************************************************
    // (1) no quadratic element(s) involved -> linear LM interpolation
    // (2) quadratic element(s) involved -> quadratic LM interpolation
    // (3) quadratic element(s) involved -> linear LM interpolation
    // (4) quadratic element(s) involved -> piecew. linear LM interpolation
    // *******************************************************************
    INPAR::MORTAR::LagMultQuad lmtype = LagMultQuad();

    // *******************************************************************
    // case (1)
    // *******************************************************************
    if (!Quad())
    {
      integrator->IntegrateDerivCell3DAuxPlane(
          SlaveElement(), MasterElement(), Cells()[i], Auxn(), Comm(), mparams_ptr);
    }
    // *******************************************************************
    // cases (2) and (3)
    // *******************************************************************
    else if ((Quad() and
                 (lmtype == INPAR::MORTAR::lagmult_quad or lmtype == INPAR::MORTAR::lagmult_lin or
                     lmtype == INPAR::MORTAR::lagmult_const)) or
             algo == INPAR::MORTAR::algorithm_gpts)
    {
      // check for standard shape functions and quadratic LM interpolation
      if (ShapeFcn() == INPAR::MORTAR::shape_standard && lmtype == INPAR::MORTAR::lagmult_quad &&
          (SlaveElement().Shape() == CORE::FE::CellType::quad8 ||
              SlaveElement().Shape() == CORE::FE::CellType::tri6))
        dserror(
            "Quad. LM interpolation for STANDARD 3D quadratic contact only feasible for "
            "quad9");

      // dynamic_cast to make sure to pass in IntElement&
      MORTAR::IntElement& sintref = dynamic_cast<MORTAR::IntElement&>(SlaveIntElement());
      MORTAR::IntElement& mintref = dynamic_cast<MORTAR::IntElement&>(MasterIntElement());

      // call integrator
      integrator->IntegrateDerivCell3DAuxPlaneQuad(
          SlaveElement(), MasterElement(), sintref, mintref, Cells()[i], Auxn());
    }

    // *******************************************************************
    // case (4)
    // *******************************************************************
    else if (Quad() && lmtype == INPAR::MORTAR::lagmult_pwlin)
    {
      // check for dual shape functions
      if (ShapeFcn() == INPAR::MORTAR::shape_dual ||
          ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
        dserror(
            "Piecewise linear LM interpolation not yet implemented for DUAL 3D quadratic "
            "contact");

      // dynamic_cast to make sure to pass in IntElement&
      MORTAR::IntElement& sintref = dynamic_cast<MORTAR::IntElement&>(SlaveIntElement());
      MORTAR::IntElement& mintref = dynamic_cast<MORTAR::IntElement&>(MasterIntElement());

      // call integrator
      integrator->IntegrateDerivCell3DAuxPlaneQuad(
          SlaveElement(), MasterElement(), sintref, mintref, Cells()[i], Auxn());
    }

    // *******************************************************************
    // undefined case
    // *******************************************************************
    else if (Quad() && lmtype == INPAR::MORTAR::lagmult_undefined)
    {
      dserror(
          "Lagrange multiplier interpolation for quadratic elements undefined\n"
          "If you are using 2nd order mortar elements, you need to specify LM_QUAD in MORTAR "
          "COUPLING section");
    }

    // *******************************************************************
    // other cases
    // *******************************************************************
    else
      dserror("IntegrateCells: Invalid case for 3D mortar contact LM interpolation");
    // *******************************************************************
  }  // cell loop

  return true;
}

/*----------------------------------------------------------------------*
 |  Linearization of clip polygon vertices (3D)               popp 02/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::VertexLinearization(
    std::vector<std::vector<CORE::GEN::pairedvector<int, double>>>& linvertex,
    std::map<int, double>& projpar, bool printderiv)
{
  // linearize all aux.plane slave and master nodes only ONCE
  // and use these linearizations later during lineclip linearization
  // (this speeds up the vertex linearizations in most cases, as we
  // never linearize the SAME slave or master vertex more than once)

  // number of nodes
  const int nsrows = SlaveIntElement().NumNode();
  const int nmrows = MasterIntElement().NumNode();

  // prepare storage for slave and master linearizations
  std::vector<std::vector<CORE::GEN::pairedvector<int, double>>> linsnodes(
      nsrows, std::vector<CORE::GEN::pairedvector<int, double>>(3, 3 * SlaveElement().NumNode()));
  std::vector<std::vector<CORE::GEN::pairedvector<int, double>>> linmnodes(
      nmrows, std::vector<CORE::GEN::pairedvector<int, double>>(
                  3, 3 * SlaveElement().NumNode() + 3 * MasterElement().NumNode()));

  // compute slave linearizations (nsrows)
  SlaveVertexLinearization(linsnodes);

  // compute master linearizations (nmrows)
  MasterVertexLinearization(linmnodes);

  //**********************************************************************
  // Clip polygon vertex linearization
  //**********************************************************************
  // loop over all clip polygon vertices
  for (int i = 0; i < (int)Clip().size(); ++i)
  {
    // references to current vertex and its linearization
    MORTAR::Vertex& currv = Clip()[i];
    std::vector<CORE::GEN::pairedvector<int, double>>& currlin = linvertex[i];

    // decision on vertex type (slave, projmaster, linclip)
    if (currv.VType() == MORTAR::Vertex::slave)
    {
      // get corresponding slave id
      int sid = currv.Nodeids()[0];

      // find corresponding slave node linearization
      int k = 0;
      while (k < nsrows)
      {
        if (SlaveIntElement().NodeIds()[k] == sid) break;
        ++k;
      }

      // dserror if not found
      if (k == nsrows) dserror("Slave Id not found!");

      // get the correct slave node linearization
      currlin = linsnodes[k];
    }
    else if (currv.VType() == MORTAR::Vertex::projmaster)
    {
      // get corresponding master id
      int mid = currv.Nodeids()[0];

      // find corresponding master node linearization
      int k = 0;
      while (k < nmrows)
      {
        if (MasterIntElement().NodeIds()[k] == mid) break;
        ++k;
      }

      // dserror if not found
      if (k == nmrows) dserror("Master Id not found!");

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
        dserror("Lineclip linearization: (S) Something went wrong!");

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
        dserror("Lineclip linearization: (M) Something went wrong!");

      MORTAR::Vertex* mv1 = &MasterVertices()[mindex1];
      MORTAR::Vertex* mv2 = &MasterVertices()[mindex2];

      // do lineclip vertex linearization
      LineclipVertexLinearization(currv, currlin, sv1, sv2, mv1, mv2, linsnodes, linmnodes);
    }

    else
      dserror("VertexLinearization: Invalid Vertex Type!");
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Linearization of slave vertex (3D) AuxPlane               popp 03/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::SlaveVertexLinearization(
    std::vector<std::vector<CORE::GEN::pairedvector<int, double>>>& currlin)
{
  // we first need the slave element center:
  // for quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  double scxi[2];

  CORE::FE::CellType dt = SlaveIntElement().Shape();
  if (dt == CORE::FE::CellType::tri3 || dt == CORE::FE::CellType::tri6)
  {
    scxi[0] = 1.0 / 3.0;
    scxi[1] = 1.0 / 3.0;
  }
  else if (dt == CORE::FE::CellType::quad4 || dt == CORE::FE::CellType::quad8 ||
           dt == CORE::FE::CellType::quad9)
  {
    scxi[0] = 0.0;
    scxi[1] = 0.0;
  }
  else
    dserror("SlaveVertexLinearization called for unknown element type");

  // evlauate shape functions + derivatives at scxi
  const int nrow = SlaveIntElement().NumNode();
  CORE::LINALG::SerialDenseVector sval(nrow);
  CORE::LINALG::SerialDenseMatrix sderiv(nrow, 2, true);
  SlaveIntElement().EvaluateShape(scxi, sval, sderiv, nrow);

  // we need all participating slave nodes
  DRT::Node** snodes = SlaveIntElement().Nodes();
  std::vector<MORTAR::Node*> smrtrnodes(nrow);

  for (int i = 0; i < nrow; ++i)
  {
    smrtrnodes[i] = dynamic_cast<MORTAR::Node*>(snodes[i]);
    if (!smrtrnodes[i]) dserror("SlaveVertexLinearization: Null pointer!");
  }

  // linearization of the IntEle spatial coords
  std::vector<std::vector<CORE::GEN::pairedvector<int, double>>> nodelin(0);
  MORTAR::IntElement* sIntEle = dynamic_cast<MORTAR::IntElement*>(&SlaveIntElement());

  if (sIntEle == nullptr)
  {
    // resize the linearizations
    nodelin.resize(nrow, std::vector<CORE::GEN::pairedvector<int, double>>(3, 1));

    // loop over all intEle nodes
    for (int in = 0; in < nrow; ++in)
      for (int dim = 0; dim < 3; ++dim) nodelin[in][dim][smrtrnodes[in]->Dofs()[dim]] += 1.;
  }
  else
    sIntEle->NodeLinearization(nodelin);

  // map iterator
  typedef CORE::GEN::pairedvector<int, double>::const_iterator
      _CI;  // linearization of element center Auxc()
  std ::vector<CORE::GEN::pairedvector<int, double>> linauxc(
      3, SlaveElement().NumNode());  // assume 3 dofs per node

  for (int i = 0; i < nrow; ++i)
    for (int dim = 0; dim < 3; ++dim)
      for (_CI p = nodelin[i][dim].begin(); p != nodelin[i][dim].end(); ++p)
        linauxc[dim][p->first] = sval[i] * p->second;

  // linearization of element normal Auxn()
  std::vector<CORE::GEN::pairedvector<int, double>>& linauxn = GetDerivAuxn();

  // put everything together for slave vertex linearization
  // loop over all vertices
  for (int i = 0; i < SlaveIntElement().NumNode(); ++i)
  {
    MORTAR::Node* mrtrsnode = dynamic_cast<MORTAR::Node*>(SlaveIntElement().Nodes()[i]);
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
    for (_CI p = linauxc[0].begin(); p != linauxc[0].end(); ++p)
      for (int k = 0; k < 3; ++k) currlin[i][k][p->first] += Auxn()[0] * Auxn()[k] * (p->second);

    for (_CI p = linauxc[1].begin(); p != linauxc[1].end(); ++p)
      for (int k = 0; k < 3; ++k) currlin[i][k][p->first] += Auxn()[1] * Auxn()[k] * (p->second);

    for (_CI p = linauxc[2].begin(); p != linauxc[2].end(); ++p)
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

  return true;
}

/*----------------------------------------------------------------------*
 |  Linearization of slave vertex (3D) AuxPlane               popp 03/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::MasterVertexLinearization(
    std::vector<std::vector<CORE::GEN::pairedvector<int, double>>>& currlin)
{
  // we first need the slave element center:
  // for quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  double scxi[2];

  CORE::FE::CellType dt = SlaveIntElement().Shape();
  if (dt == CORE::FE::CellType::tri3 || dt == CORE::FE::CellType::tri6)
  {
    scxi[0] = 1.0 / 3.0;
    scxi[1] = 1.0 / 3.0;
  }
  else if (dt == CORE::FE::CellType::quad4 || dt == CORE::FE::CellType::quad8 ||
           dt == CORE::FE::CellType::quad9)
  {
    scxi[0] = 0.0;
    scxi[1] = 0.0;
  }
  else
    dserror("MasterVertexLinearization called for unknown element type");

  // evlauate shape functions + derivatives at scxi
  int nrow = SlaveIntElement().NumNode();
  CORE::LINALG::SerialDenseVector sval(nrow);
  CORE::LINALG::SerialDenseMatrix sderiv(nrow, 2, true);
  SlaveIntElement().EvaluateShape(scxi, sval, sderiv, nrow);

  // we need all participating slave nodes
  DRT::Node** snodes = SlaveIntElement().Nodes();
  std::vector<MORTAR::Node*> smrtrnodes(nrow);

  for (int i = 0; i < nrow; ++i)
  {
    smrtrnodes[i] = dynamic_cast<MORTAR::Node*>(snodes[i]);
    if (!smrtrnodes[i]) dserror("MasterVertexLinearization: Null pointer!");
  }

  // linearization of the SlaveIntEle spatial coords
  std::vector<std::vector<CORE::GEN::pairedvector<int, double>>> snodelin(0);
  MORTAR::IntElement* sIntEle = dynamic_cast<MORTAR::IntElement*>(&SlaveIntElement());

  if (sIntEle == nullptr)
  {
    // resize the linearizations
    snodelin.resize(nrow, std::vector<CORE::GEN::pairedvector<int, double>>(3, 1));

    // loop over all intEle nodes
    for (int in = 0; in < nrow; ++in)
      for (int dim = 0; dim < 3; ++dim) snodelin[in][dim][smrtrnodes[in]->Dofs()[dim]] += 1.;
  }
  else
    sIntEle->NodeLinearization(snodelin);

  // map iterator
  typedef CORE::GEN::pairedvector<int, double>::const_iterator
      _CI;  // linearization of element center Auxc()
  std ::vector<CORE::GEN::pairedvector<int, double>> linauxc(
      3, SlaveElement().NumNode());  // assume 3 dofs per node

  for (int i = 0; i < nrow; ++i)
    for (int dim = 0; dim < 3; ++dim)
      for (_CI p = snodelin[i][dim].begin(); p != snodelin[i][dim].end(); ++p)
        linauxc[dim][p->first] = sval[i] * p->second;

  // linearization of element normal Auxn()
  std::vector<CORE::GEN::pairedvector<int, double>>& linauxn = GetDerivAuxn();

  // linearization of the MasterIntEle spatial coords
  std::vector<std::vector<CORE::GEN::pairedvector<int, double>>> mnodelin(0);
  MORTAR::IntElement* mIntEle = dynamic_cast<MORTAR::IntElement*>(&MasterIntElement());

  if (mIntEle == nullptr)
  {
    // resize the linearizations
    mnodelin.resize(
        MasterIntElement().NumNode(), std::vector<CORE::GEN::pairedvector<int, double>>(3, 1));

    // loop over all intEle nodes
    for (int in = 0; in < MasterIntElement().NumNode(); ++in)
    {
      MORTAR::Node* mrtrmnode = dynamic_cast<MORTAR::Node*>(MasterIntElement().Nodes()[in]);
      if (mrtrmnode == nullptr) dserror("dynamic cast to mortar node went wrong");

      for (int dim = 0; dim < 3; ++dim) mnodelin[in][dim][mrtrmnode->Dofs()[dim]] += 1.;
    }
  }
  else
    mIntEle->NodeLinearization(mnodelin);

  // put everything together for slave vertex linearization
  // loop over all vertices
  for (int i = 0; i < MasterIntElement().NumNode(); ++i)
  {
    MORTAR::Node* mrtrmnode = dynamic_cast<MORTAR::Node*>(MasterIntElement().Nodes()[i]);
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
    for (_CI p = linauxc[0].begin(); p != linauxc[0].end(); ++p)
      for (int k = 0; k < 3; ++k) currlin[i][k][p->first] += Auxn()[0] * Auxn()[k] * (p->second);

    for (_CI p = linauxc[1].begin(); p != linauxc[1].end(); ++p)
      for (int k = 0; k < 3; ++k) currlin[i][k][p->first] += Auxn()[1] * Auxn()[k] * (p->second);

    for (_CI p = linauxc[2].begin(); p != linauxc[2].end(); ++p)
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

  return true;
}

/*----------------------------------------------------------------------*
 |  Linearization of lineclip vertex (3D) AuxPlane            popp 03/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::LineclipVertexLinearization(MORTAR::Vertex& currv,
    std::vector<CORE::GEN::pairedvector<int, double>>& currlin, MORTAR::Vertex* sv1,
    MORTAR::Vertex* sv2, MORTAR::Vertex* mv1, MORTAR::Vertex* mv2,
    std::vector<std::vector<CORE::GEN::pairedvector<int, double>>>& linsnodes,
    std::vector<std::vector<CORE::GEN::pairedvector<int, double>>>& linmnodes)
{
  // number of nodes
  const int nsrows = SlaveIntElement().NumNode();
  const int nmrows = MasterIntElement().NumNode();

  // iterator
  typedef CORE::GEN::pairedvector<int, double>::const_iterator _CI;

  // compute factor Z
  std::array<double, 3> crossZ = {0.0, 0.0, 0.0};
  crossZ[0] = (sv1->Coord()[1] - mv1->Coord()[1]) * (mv2->Coord()[2] - mv1->Coord()[2]) -
              (sv1->Coord()[2] - mv1->Coord()[2]) * (mv2->Coord()[1] - mv1->Coord()[1]);
  crossZ[1] = (sv1->Coord()[2] - mv1->Coord()[2]) * (mv2->Coord()[0] - mv1->Coord()[0]) -
              (sv1->Coord()[0] - mv1->Coord()[0]) * (mv2->Coord()[2] - mv1->Coord()[2]);
  crossZ[2] = (sv1->Coord()[0] - mv1->Coord()[0]) * (mv2->Coord()[1] - mv1->Coord()[1]) -
              (sv1->Coord()[1] - mv1->Coord()[1]) * (mv2->Coord()[0] - mv1->Coord()[0]);
  double Zfac = crossZ[0] * Auxn()[0] + crossZ[1] * Auxn()[1] + crossZ[2] * Auxn()[2];

  // compute factor N
  std::array<double, 3> crossN = {0.0, 0.0, 0.0};
  crossN[0] = (sv2->Coord()[1] - sv1->Coord()[1]) * (mv2->Coord()[2] - mv1->Coord()[2]) -
              (sv2->Coord()[2] - sv1->Coord()[2]) * (mv2->Coord()[1] - mv1->Coord()[1]);
  crossN[1] = (sv2->Coord()[2] - sv1->Coord()[2]) * (mv2->Coord()[0] - mv1->Coord()[0]) -
              (sv2->Coord()[0] - sv1->Coord()[0]) * (mv2->Coord()[2] - mv1->Coord()[2]);
  crossN[2] = (sv2->Coord()[0] - sv1->Coord()[0]) * (mv2->Coord()[1] - mv1->Coord()[1]) -
              (sv2->Coord()[1] - sv1->Coord()[1]) * (mv2->Coord()[0] - mv1->Coord()[0]);
  double Nfac = crossN[0] * Auxn()[0] + crossN[1] * Auxn()[1] + crossN[2] * Auxn()[2];

  // slave edge vector
  std::array<double, 3> sedge = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; ++k) sedge[k] = sv2->Coord()[k] - sv1->Coord()[k];

  // prepare linearization derivZ
  std::array<double, 3> crossdZ1 = {0.0, 0.0, 0.0};
  std::array<double, 3> crossdZ2 = {0.0, 0.0, 0.0};
  std::array<double, 3> crossdZ3 = {0.0, 0.0, 0.0};
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
  std::array<double, 3> crossdN1 = {0.0, 0.0, 0.0};
  std::array<double, 3> crossdN2 = {0.0, 0.0, 0.0};
  std::array<double, 3> crossdN3 = {0.0, 0.0, 0.0};
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
    if (SlaveIntElement().NodeIds()[k] == sid1) break;
    ++k;
  }

  // dserror if not found
  if (k == nsrows) dserror("Slave Id1 not found!");

  // get the correct slave node linearization
  std::vector<CORE::GEN::pairedvector<int, double>>& slavelin0 = linsnodes[k];

  k = 0;
  while (k < nsrows)
  {
    if (SlaveIntElement().NodeIds()[k] == sid2) break;
    ++k;
  }

  // dserror if not found
  if (k == nsrows) dserror("Slave Id2 not found!");

  // get the correct slave node linearization
  std::vector<CORE::GEN::pairedvector<int, double>>& slavelin1 = linsnodes[k];

  // master vertex linearization (2x)
  int mid1 = currv.Nodeids()[2];
  int mid2 = currv.Nodeids()[3];

  // find corresponding master node linearizations
  k = 0;
  while (k < nmrows)
  {
    if (MasterIntElement().NodeIds()[k] == mid1) break;
    ++k;
  }

  // dserror if not found
  if (k == nmrows) dserror("Master Id1 not found!");

  // get the correct master node linearization
  std::vector<CORE::GEN::pairedvector<int, double>>& masterlin0 = linmnodes[k];

  k = 0;
  while (k < nmrows)
  {
    if (MasterIntElement().NodeIds()[k] == mid2) break;
    ++k;
  }

  // dserror if not found
  if (k == nmrows) dserror("Master Id2 not found!");

  // get the correct master node linearization
  std::vector<CORE::GEN::pairedvector<int, double>>& masterlin1 = linmnodes[k];

  // linearization of element normal Auxn()
  std::vector<CORE::GEN::pairedvector<int, double>>& linauxn = GetDerivAuxn();

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

  return true;
}

/*----------------------------------------------------------------------*
 |  Linearization of clip polygon center (3D)                 popp 02/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::CenterLinearization(
    const std::vector<std::vector<CORE::GEN::pairedvector<int, double>>>& linvertex,
    std::vector<CORE::GEN::pairedvector<int, double>>& lincenter)
{
  // preparations
  int clipsize = (int)(Clip().size());
  typedef CORE::GEN::pairedvector<int, double>::const_iterator CI;

  // number of nodes
  const int nsrows = SlaveElement().NumNode();
  const int nmrows = MasterElement().NumNode();

  std::vector<double> clipcenter(3);
  for (int k = 0; k < 3; ++k) clipcenter[k] = 0.0;
  double fac = 0.0;

  // first we need node averaged center
  std::array<double, 3> nac = {0.0, 0.0, 0.0};
  for (int i = 0; i < clipsize; ++i)
    for (int k = 0; k < 3; ++k) nac[k] += (Clip()[i].Coord()[k] / clipsize);

  // loop over all triangles of polygon (1st round: preparations)
  for (int i = 0; i < clipsize; ++i)
  {
    std::array<double, 3> xi_i = {0.0, 0.0, 0.0};
    std::array<double, 3> xi_ip1 = {0.0, 0.0, 0.0};

    // standard case
    if (i < clipsize - 1)
    {
      for (int k = 0; k < 3; ++k) xi_i[k] = Clip()[i].Coord()[k];
      for (int k = 0; k < 3; ++k) xi_ip1[k] = Clip()[i + 1].Coord()[k];
    }
    // last vertex of clip polygon
    else
    {
      for (int k = 0; k < 3; ++k) xi_i[k] = Clip()[clipsize - 1].Coord()[k];
      for (int k = 0; k < 3; ++k) xi_ip1[k] = Clip()[0].Coord()[k];
    }

    // triangle area
    std::array<double, 3> diff1 = {0.0, 0.0, 0.0};
    std::array<double, 3> diff2 = {0.0, 0.0, 0.0};
    for (int k = 0; k < 3; ++k) diff1[k] = xi_ip1[k] - xi_i[k];
    for (int k = 0; k < 3; ++k) diff2[k] = xi_i[k] - nac[k];

    std::array<double, 3> cross = {0.0, 0.0, 0.0};
    cross[0] = diff1[1] * diff2[2] - diff1[2] * diff2[1];
    cross[1] = diff1[2] * diff2[0] - diff1[0] * diff2[2];
    cross[2] = diff1[0] * diff2[1] - diff1[1] * diff2[0];

    double Atri = 0.5 * sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);

    // add contributions to clipcenter and fac
    fac += Atri;
    for (int k = 0; k < 3; ++k) clipcenter[k] += 1.0 / 3.0 * (xi_i[k] + xi_ip1[k] + nac[k]) * Atri;
  }

  // build factors for linearization
  std::array<double, 3> z = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; ++k) z[k] = clipcenter[k];
  double n = fac;

  // first we need linearization of node averaged center
  std::vector<CORE::GEN::pairedvector<int, double>> linnac(3, 3 * (nsrows + nmrows));
  const double clipsizeinv = 1.0 / clipsize;

  for (int i = 0; i < clipsize; ++i)
    for (int k = 0; k < 3; ++k)
      for (CI p = linvertex[i][k].begin(); p != linvertex[i][k].end(); ++p)
        linnac[k][p->first] += clipsizeinv * (p->second);

  // loop over all triangles of polygon (2nd round: linearization)
  for (int i = 0; i < clipsize; ++i)
  {
    std::array<double, 3> xi_i = {0.0, 0.0, 0.0};
    std::array<double, 3> xi_ip1 = {0.0, 0.0, 0.0};
    int iplus1 = 0;

    // standard case
    if (i < clipsize - 1)
    {
      for (int k = 0; k < 3; ++k) xi_i[k] = Clip()[i].Coord()[k];
      for (int k = 0; k < 3; ++k) xi_ip1[k] = Clip()[i + 1].Coord()[k];
      iplus1 = i + 1;
    }
    // last vertex of clip polygon
    else
    {
      for (int k = 0; k < 3; ++k) xi_i[k] = Clip()[clipsize - 1].Coord()[k];
      for (int k = 0; k < 3; ++k) xi_ip1[k] = Clip()[0].Coord()[k];
      iplus1 = 0;
    }

    // triangle area
    std::array<double, 3> diff1 = {0.0, 0.0, 0.0};
    std::array<double, 3> diff2 = {0.0, 0.0, 0.0};
    for (int k = 0; k < 3; ++k) diff1[k] = xi_ip1[k] - xi_i[k];
    for (int k = 0; k < 3; ++k) diff2[k] = xi_i[k] - nac[k];

    std::array<double, 3> cross = {0.0, 0.0, 0.0};
    cross[0] = diff1[1] * diff2[2] - diff1[2] * diff2[1];
    cross[1] = diff1[2] * diff2[0] - diff1[0] * diff2[2];
    cross[2] = diff1[0] * diff2[1] - diff1[1] * diff2[0];

    double Atri = 0.5 * sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);

    // linearization of cross
    std::vector<CORE::GEN::pairedvector<int, double>> lincross(3, 3 * (nsrows + nmrows));

    for (CI p = linvertex[i][0].begin(); p != linvertex[i][0].end(); ++p)
    {
      lincross[1][p->first] += diff1[2] * (p->second);
      lincross[1][p->first] += diff2[2] * (p->second);
      lincross[2][p->first] -= diff1[1] * (p->second);
      lincross[2][p->first] -= diff2[1] * (p->second);
    }
    for (CI p = linvertex[i][1].begin(); p != linvertex[i][1].end(); ++p)
    {
      lincross[0][p->first] -= diff1[2] * (p->second);
      lincross[0][p->first] -= diff2[2] * (p->second);
      lincross[2][p->first] += diff1[0] * (p->second);
      lincross[2][p->first] += diff2[0] * (p->second);
    }
    for (CI p = linvertex[i][2].begin(); p != linvertex[i][2].end(); ++p)
    {
      lincross[0][p->first] += diff1[1] * (p->second);
      lincross[0][p->first] += diff2[1] * (p->second);
      lincross[1][p->first] -= diff1[0] * (p->second);
      lincross[1][p->first] -= diff2[0] * (p->second);
    }

    for (CI p = linvertex[iplus1][0].begin(); p != linvertex[iplus1][0].end(); ++p)
    {
      lincross[1][p->first] -= diff2[2] * (p->second);
      lincross[2][p->first] += diff2[1] * (p->second);
    }
    for (CI p = linvertex[iplus1][1].begin(); p != linvertex[iplus1][1].end(); ++p)
    {
      lincross[0][p->first] += diff2[2] * (p->second);
      lincross[2][p->first] -= diff2[0] * (p->second);
    }
    for (CI p = linvertex[iplus1][2].begin(); p != linvertex[iplus1][2].end(); ++p)
    {
      lincross[0][p->first] -= diff2[1] * (p->second);
      lincross[1][p->first] += diff2[0] * (p->second);
    }

    for (CI p = linnac[0].begin(); p != linnac[0].end(); ++p)
    {
      lincross[1][p->first] -= diff1[2] * (p->second);
      lincross[2][p->first] += diff1[1] * (p->second);
    }
    for (CI p = linnac[1].begin(); p != linnac[1].end(); ++p)
    {
      lincross[0][p->first] += diff1[2] * (p->second);
      lincross[2][p->first] -= diff1[0] * (p->second);
    }
    for (CI p = linnac[2].begin(); p != linnac[2].end(); ++p)
    {
      lincross[0][p->first] -= diff1[1] * (p->second);
      lincross[1][p->first] += diff1[0] * (p->second);
    }

    // linearization of triangle area
    CORE::GEN::pairedvector<int, double> linarea(3 * (nsrows + nmrows));
    for (int k = 0; k < 3; ++k)
      for (CI p = lincross[k].begin(); p != lincross[k].end(); ++p)
        linarea[p->first] += 0.25 / Atri * cross[k] * (p->second);

    const double fac1 = 1.0 / (3.0 * n);

    // put everything together
    for (int k = 0; k < 3; ++k)
    {
      for (CI p = linvertex[i][k].begin(); p != linvertex[i][k].end(); ++p)
        lincenter[k][p->first] += fac1 * Atri * (p->second);

      for (CI p = linvertex[iplus1][k].begin(); p != linvertex[iplus1][k].end(); ++p)
        lincenter[k][p->first] += fac1 * Atri * (p->second);

      for (CI p = linnac[k].begin(); p != linnac[k].end(); ++p)
        lincenter[k][p->first] += fac1 * Atri * (p->second);

      for (CI p = linarea.begin(); p != linarea.end(); ++p)
      {
        lincenter[k][p->first] += fac1 * (xi_i[k] + xi_ip1[k] + nac[k]) * (p->second);
        lincenter[k][p->first] -= z[k] / (n * n) * (p->second);
      }
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::Coupling3dQuad::Coupling3dQuad(DRT::Discretization& idiscret, int dim, bool quad,
    Teuchos::ParameterList& params, MORTAR::Element& sele, MORTAR::Element& mele,
    MORTAR::IntElement& sintele, MORTAR::IntElement& mintele)
    : CONTACT::Coupling3d(idiscret, dim, quad, params, sele, mele),
      sintele_(sintele),
      mintele_(mintele)
{
  //  3D quadratic coupling only for quadratic ansatz type
  if (!Quad()) dserror("Coupling3dQuad called for non-quadratic ansatz!");

  return;
}

/*----------------------------------------------------------------------*
 |  get communicator  (public)                               farah 01/13|
 *----------------------------------------------------------------------*/
const Epetra_Comm& CONTACT::Coupling3dManager::Comm() const { return idiscret_.Comm(); }

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::Coupling3dManager::Coupling3dManager(DRT::Discretization& idiscret, int dim, bool quad,
    Teuchos::ParameterList& params, MORTAR::Element* sele, std::vector<MORTAR::Element*> mele)
    : idiscret_(idiscret),
      dim_(dim),
      quad_(quad),
      imortar_(params),
      sele_(sele),
      mele_(mele),
      ncells_(0),
      stype_(CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(params, "STRATEGY"))
{
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 01/13|
 *----------------------------------------------------------------------*/
CONTACT::Coupling3dQuadManager::Coupling3dQuadManager(DRT::Discretization& idiscret, int dim,
    bool quad, Teuchos::ParameterList& params, MORTAR::Element* sele,
    std::vector<MORTAR::Element*> mele)
    : MORTAR::Coupling3dQuadManager(idiscret, dim, quad, params, sele, mele),
      CONTACT::Coupling3dManager(idiscret, dim, quad, params, sele, mele),
      smintpairs_(-1),
      intcells_(-1)
{
  return;
}


/*----------------------------------------------------------------------*
 |  Evaluate mortar coupling pairs                            popp 03/09|
 *----------------------------------------------------------------------*/
void CONTACT::Coupling3dManager::IntegrateCoupling(
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  // get algorithm
  INPAR::MORTAR::AlgorithmType algo =
      CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(imortar_, "ALGORITHM");

  // prepare linearizations
  if (algo == INPAR::MORTAR::algorithm_mortar)
    dynamic_cast<CONTACT::Element&>(SlaveElement()).PrepareDderiv(MasterElements());

  // decide which type of numerical integration scheme

  //**********************************************************************
  // STANDARD INTEGRATION (SEGMENTS)
  //**********************************************************************
  if (IntType() == INPAR::MORTAR::inttype_segments)
  {
    // loop over all master elements associated with this slave element
    for (int m = 0; m < (int)MasterElements().size(); ++m)
    {
      // create Coupling3d object and push back
      Coupling().push_back(Teuchos::rcp(
          new Coupling3d(idiscret_, dim_, false, imortar_, SlaveElement(), MasterElement(m))));

      // do coupling
      Coupling()[m]->EvaluateCoupling();

      // store number of intcells
      ncells_ += (int)(Coupling()[m]->Cells()).size();
    }

    // special treatment of boundary elements
    ConsistDualShape();

    // TODO modification of boundary shapes!

    // integrate cells
    for (int i = 0; i < (int)Coupling().size(); ++i)
    {
      // temporary m-matrix linearization of this slave/master pair
      if (algo == INPAR::MORTAR::algorithm_mortar)
        dynamic_cast<CONTACT::Element&>(SlaveElement()).PrepareMderiv(MasterElements(), i);

      // integrate cells
      Coupling()[i]->IntegrateCells(mparams_ptr);

      // assemble m-matrix for this slave/master pair
      if (algo == INPAR::MORTAR::algorithm_mortar)
        dynamic_cast<CONTACT::Element&>(SlaveElement())
            .AssembleMderivToNodes(Coupling()[i]->MasterElement());
    }
  }

  //**********************************************************************
  // ELEMENT-BASED INTEGRATION
  //**********************************************************************
  else if (IntType() == INPAR::MORTAR::inttype_elements ||
           IntType() == INPAR::MORTAR::inttype_elements_BS)
  {
    if ((int)MasterElements().size() == 0) return;

    if (!Quad())
    {
      bool boundary_ele = false;
      bool proj = false;

      /* find all feasible master elements (this check is inherent in the
       * segment based integration)                    hiermeier 04/16 */
      std::vector<MORTAR::Element*> feasible_ma_eles(MasterElements().size());
      FindFeasibleMasterElements(feasible_ma_eles);

      // create an integrator instance with correct NumGP and Dim
      Teuchos::RCP<CONTACT::Integrator> integrator =
          CONTACT::INTEGRATOR::BuildIntegrator(stype_, imortar_, SlaveElement().Shape(), Comm());

      // Perform integration and linearization
      integrator->IntegrateDerivEle3D(
          SlaveElement(), feasible_ma_eles, &boundary_ele, &proj, Comm(), mparams_ptr);


      if (IntType() == INPAR::MORTAR::inttype_elements_BS)
      {
        if (boundary_ele == true)
        {
          // loop over all master elements associated with this slave element
          for (int m = 0; m < (int)MasterElements().size(); ++m)
          {
            // create Coupling3d object and push back
            Coupling().push_back(Teuchos::rcp(new Coupling3d(
                idiscret_, dim_, false, imortar_, SlaveElement(), MasterElement(m))));

            // do coupling
            Coupling()[m]->EvaluateCoupling();

            // store number of intcells
            ncells_ += (int)(Coupling()[m]->Cells()).size();
          }

          // special treatment of boundary elements
          ConsistDualShape();

          // integrate cells
          for (int i = 0; i < (int)Coupling().size(); ++i)
          {
            // temporary m-matrix linearization of this slave/master pair
            if (algo == INPAR::MORTAR::algorithm_mortar)
              dynamic_cast<CONTACT::Element&>(SlaveElement()).PrepareMderiv(MasterElements(), i);

            // integrate cells
            Coupling()[i]->IntegrateCells(mparams_ptr);

            // assemble m-matrix for this slave/master pair
            if (algo == INPAR::MORTAR::algorithm_mortar)
              dynamic_cast<CONTACT::Element&>(SlaveElement())
                  .AssembleMderivToNodes(Coupling()[i]->MasterElement());
          }
        }
      }
    }
    else
    {
      dserror(
          "You should not be here! This coupling manager is not able to perform mortar coupling "
          "for high-order elements.");
    }
  }
  //**********************************************************************
  // INVALID TYPE OF NUMERICAL INTEGRATION
  //**********************************************************************
  else
  {
    dserror("Invalid type of numerical integration!");
  }

  // free memory of dual shape function coefficient matrix
  SlaveElement().MoData().ResetDualShape();
  SlaveElement().MoData().ResetDerivDualShape();

  // assemble element contribution to nodes
  if (algo == INPAR::MORTAR::algorithm_mortar)
  {
    bool dual = (ShapeFcn() == INPAR::MORTAR::shape_dual) ||
                (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin);
    dynamic_cast<CONTACT::Element&>(SlaveElement()).AssembleDderivToNodes(dual);
  }
  return;
}


/*----------------------------------------------------------------------*
 |  Evaluate coupling pairs                                 farah 09/14 |
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3dManager::EvaluateCoupling(
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  // decide which type of coupling should be evaluated
  INPAR::MORTAR::AlgorithmType algo =
      CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(imortar_, "ALGORITHM");

  //*********************************
  // Mortar Contact
  //*********************************
  if (algo == INPAR::MORTAR::algorithm_mortar || algo == INPAR::MORTAR::algorithm_gpts)
    IntegrateCoupling(mparams_ptr);

  //*********************************
  // Error
  //*********************************
  else
    dserror("chose contact algorithm not supported!");

  // interpolate temperatures in TSI case
  if (imortar_.get<int>("PROBTYPE") == INPAR::CONTACT::tsi)
    NTS::Interpolator(imortar_, dim_).InterpolateMasterTemp3D(SlaveElement(), MasterElements());

  return true;
}


/*----------------------------------------------------------------------*
 |  Evaluate mortar coupling pairs for Quad-coupling         farah 09/14|
 *----------------------------------------------------------------------*/
void CONTACT::Coupling3dQuadManager::IntegrateCoupling(
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  // get algorithm type
  INPAR::MORTAR::AlgorithmType algo = CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(
      MORTAR::Coupling3dQuadManager::imortar_, "ALGORITHM");

  // prepare linearizations
  if (algo == INPAR::MORTAR::algorithm_mortar)
    dynamic_cast<CONTACT::Element&>(SlaveElement()).PrepareDderiv(MasterElements());

  // decide which type of numerical integration scheme

  //**********************************************************************
  // STANDARD INTEGRATION (SEGMENTS)
  //**********************************************************************
  if (IntType() == INPAR::MORTAR::inttype_segments)
  {
    Coupling().resize(0);

    // build linear integration elements from quadratic MORTAR::Elements
    std::vector<Teuchos::RCP<MORTAR::IntElement>> sauxelements(0);
    std::vector<std::vector<Teuchos::RCP<MORTAR::IntElement>>> mauxelements(
        MasterElements().size());
    SplitIntElements(SlaveElement(), sauxelements);

    // loop over all master elements associated with this slave element
    for (int m = 0; m < (int)MasterElements().size(); ++m)
    {
      // build linear integration elements from quadratic MORTAR::Elements
      mauxelements[m].resize(0);
      SplitIntElements(*MasterElements()[m], mauxelements[m]);

      // loop over all IntElement pairs for coupling
      for (int i = 0; i < (int)sauxelements.size(); ++i)
      {
        for (int j = 0; j < (int)mauxelements[m].size(); ++j)
        {
          Coupling().push_back(Teuchos::rcp(new Coupling3dQuad(Discret(), Dim(), true, Params(),
              SlaveElement(), *MasterElements()[m], *sauxelements[i], *mauxelements[m][j])));

          Coupling()[Coupling().size() - 1]->EvaluateCoupling();

          // increase counter of slave/master integration pairs and intcells
          smintpairs_ += 1;
          intcells_ += (int)Coupling()[Coupling().size() - 1]->Cells().size();
        }  // for maux
      }    // for saux
    }      // for m

    ConsistDualShape();

    // integrate cells
    for (int i = 0; i < (int)Coupling().size(); ++i)
    {
      if (algo == INPAR::MORTAR::algorithm_mortar)
        dynamic_cast<CONTACT::Element&>(SlaveElement())
            .PrepareMderiv(MasterElements(), i % mauxelements.size());

      Coupling()[i]->IntegrateCells(mparams_ptr);

      if (algo == INPAR::MORTAR::algorithm_mortar)
        dynamic_cast<CONTACT::Element&>(SlaveElement())
            .AssembleMderivToNodes(Coupling()[i]->MasterElement());
    }
  }

  //**********************************************************************
  // FAST INTEGRATION (ELEMENTS)
  //**********************************************************************
  else if (IntType() == INPAR::MORTAR::inttype_elements ||
           IntType() == INPAR::MORTAR::inttype_elements_BS)
  {
    // check for standard shape functions and quadratic LM interpolation
    if (ShapeFcn() == INPAR::MORTAR::shape_standard &&
        LagMultQuad() == INPAR::MORTAR::lagmult_quad &&
        (SlaveElement().Shape() == CORE::FE::CellType::quad8 ||
            SlaveElement().Shape() == CORE::FE::CellType::tri6))
      dserror(
          "Quad. LM interpolation for STANDARD 3D quadratic contact only feasible for "
          "quad9");

    if ((int)MasterElements().size() == 0) return;

    // create an integrator instance with correct NumGP and Dim
    Teuchos::RCP<CONTACT::Integrator> integrator =
        CONTACT::INTEGRATOR::BuildIntegrator(stype_, Params(), SlaveElement().Shape(), Comm());

    bool boundary_ele = false;
    bool proj = false;

    // Perform integration and linearization
    integrator->IntegrateDerivEle3D(
        SlaveElement(), MasterElements(), &boundary_ele, &proj, Comm(), mparams_ptr);

    if (IntType() == INPAR::MORTAR::inttype_elements_BS)
    {
      if (boundary_ele == true)
      {
        Coupling().resize(0);

        // build linear integration elements from quadratic MORTAR::Elements
        std::vector<Teuchos::RCP<MORTAR::IntElement>> sauxelements(0);
        std::vector<std::vector<Teuchos::RCP<MORTAR::IntElement>>> mauxelements(
            MasterElements().size());
        SplitIntElements(SlaveElement(), sauxelements);

        // loop over all master elements associated with this slave element
        for (int m = 0; m < (int)MasterElements().size(); ++m)
        {
          // build linear integration elements from quadratic MORTAR::Elements
          mauxelements[m].resize(0);
          SplitIntElements(*MasterElements()[m], mauxelements[m]);

          // loop over all IntElement pairs for coupling
          for (int i = 0; i < (int)sauxelements.size(); ++i)
          {
            for (int j = 0; j < (int)mauxelements[m].size(); ++j)
            {
              Coupling().push_back(Teuchos::rcp(new Coupling3dQuad(Discret(), Dim(), true, Params(),
                  SlaveElement(), *MasterElements()[m], *sauxelements[i], *mauxelements[m][j])));

              Coupling()[Coupling().size() - 1]->EvaluateCoupling();

              // increase counter of slave/master integration pairs and intcells
              smintpairs_ += 1;
              intcells_ += (int)Coupling()[Coupling().size() - 1]->Cells().size();
            }  // for maux
          }    // for saux
        }      // for m

        ConsistDualShape();

        for (int i = 0; i < (int)Coupling().size(); ++i)
        {
          if (algo == INPAR::MORTAR::algorithm_mortar)
            dynamic_cast<CONTACT::Element&>(SlaveElement())
                .PrepareMderiv(MasterElements(), i % mauxelements.size());
          Coupling()[i]->IntegrateCells(mparams_ptr);
          if (algo == INPAR::MORTAR::algorithm_mortar)
            dynamic_cast<CONTACT::Element&>(SlaveElement())
                .AssembleMderivToNodes(Coupling()[i]->MasterElement());
        }
      }
    }
  }
  //**********************************************************************
  // INVALID
  //**********************************************************************
  else
  {
    dserror("Invalid type of numerical integration");
  }

  // free memory of consistent dual shape function coefficient matrix
  SlaveElement().MoData().ResetDualShape();
  SlaveElement().MoData().ResetDerivDualShape();

  if (algo == INPAR::MORTAR::algorithm_mortar)
    dynamic_cast<CONTACT::Element&>(SlaveElement())
        .AssembleDderivToNodes((ShapeFcn() == INPAR::MORTAR::shape_dual ||
                                ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin));

  return;
}


/*----------------------------------------------------------------------*
 |  Evaluate coupling pairs for Quad-coupling                farah 01/13|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3dQuadManager::EvaluateCoupling(
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  // decide which type of coupling should be evaluated
  INPAR::MORTAR::AlgorithmType algo =
      CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(Params(), "ALGORITHM");

  //*********************************
  // Mortar Contact
  //*********************************
  if (algo == INPAR::MORTAR::algorithm_mortar || algo == INPAR::MORTAR::algorithm_gpts)
    IntegrateCoupling(mparams_ptr);

  //*********************************
  // Error
  //*********************************
  else
    dserror("chosen contact algorithm not supported!");

  return true;
}


/*----------------------------------------------------------------------*
 |  Calculate dual shape functions                           seitz 07/13|
 *----------------------------------------------------------------------*/
void CONTACT::Coupling3dManager::ConsistDualShape()
{
  static const INPAR::MORTAR::AlgorithmType algo =
      CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(imortar_, "ALGORITHM");
  if (algo != INPAR::MORTAR::algorithm_mortar) return;

  // For standard shape functions no modification is necessary
  // A switch erlier in the process improves computational efficiency
  INPAR::MORTAR::ConsistentDualType consistent =
      CORE::UTILS::IntegralValue<INPAR::MORTAR::ConsistentDualType>(imortar_, "LM_DUAL_CONSISTENT");
  if (ShapeFcn() == INPAR::MORTAR::shape_standard || consistent == INPAR::MORTAR::consistent_none)
    return;

  // Consistent modification not yet checked for constant LM interpolation
  if (Quad() == true && LagMultQuad() == INPAR::MORTAR::lagmult_const &&
      consistent != INPAR::MORTAR::consistent_none)
    dserror("Consistent dual shape functions not yet checked for constant LM interpolation!");

  if (consistent == INPAR::MORTAR::consistent_all && IntType() != INPAR::MORTAR::inttype_segments)
    dserror(
        "Consistent dual shape functions on all elements only for segment-based "
        "integration");

  // do nothing if there are no coupling pairs
  if (Coupling().size() == 0) return;

  // check for boundary elements in segment-based integration
  // (fast integration already has this check, so that ConsistDualShape()
  // is only called for boundary elements)
  //
  // For NURBS elements, always compute consistent dual functions.
  // This improves robustness, since the duality is enforced at exactly
  // the same quadrature points, that the mortar integrals etc. are evaluated.
  // For Lagrange FE, the calculation of dual shape functions for fully
  // projecting elements is ok, since the integrands are polynomials (except
  // the jacobian)
  if (IntType() == INPAR::MORTAR::inttype_segments &&
      consistent == INPAR::MORTAR::consistent_boundary)
  {
    // check, if slave element is fully projecting
    // for convenience, we don't check each quadrature point
    // but only the element nodes. This usually does the job.
    bool boundary_ele = false;

    CORE::FE::CellType dt_s = SlaveElement().Shape();

    double sxi_test[2] = {0.0, 0.0};
    double alpha_test = 0.0;
    bool proj_test = false;

    DRT::Node** mynodes_test = SlaveElement().Nodes();
    if (!mynodes_test) dserror("HasProjStatus: Null pointer!");

    if (dt_s == CORE::FE::CellType::quad4 || dt_s == CORE::FE::CellType::quad8 ||
        dt_s == CORE::FE::CellType::nurbs9)
    {
      for (int s_test = 0; s_test < SlaveElement().NumNode(); ++s_test)
      {
        if (s_test == 0)
        {
          sxi_test[0] = -1.0;
          sxi_test[1] = -1.0;
        }
        else if (s_test == 1)
        {
          sxi_test[0] = -1.0;
          sxi_test[1] = 1.0;
        }
        else if (s_test == 2)
        {
          sxi_test[0] = 1.0;
          sxi_test[1] = -1.0;
        }
        else if (s_test == 3)
        {
          sxi_test[0] = 1.0;
          sxi_test[1] = 1.0;
        }
        else if (s_test == 4)
        {
          sxi_test[0] = 1.0;
          sxi_test[1] = 0.0;
        }
        else if (s_test == 5)
        {
          sxi_test[0] = 0.0;
          sxi_test[1] = 1.0;
        }
        else if (s_test == 6)
        {
          sxi_test[0] = -1.0;
          sxi_test[1] = 0.0;
        }
        else if (s_test == 7)
        {
          sxi_test[0] = 0.0;
          sxi_test[1] = -1.0;
        }

        proj_test = false;
        for (int bs_test = 0; bs_test < (int)Coupling().size(); ++bs_test)
        {
          double mxi_test[2] = {0.0, 0.0};
          MORTAR::Projector::Impl(SlaveElement(), Coupling()[bs_test]->MasterIntElement())
              ->ProjectGaussPoint3D(SlaveElement(), sxi_test,
                  Coupling()[bs_test]->MasterIntElement(), mxi_test, alpha_test);

          CORE::FE::CellType dt = Coupling()[bs_test]->MasterIntElement().Shape();
          if (dt == CORE::FE::CellType::quad4 || dt == CORE::FE::CellType::quad8 ||
              dt == CORE::FE::CellType::quad9)
          {
            if (mxi_test[0] >= -1.0 && mxi_test[1] >= -1.0 && mxi_test[0] <= 1.0 &&
                mxi_test[1] <= 1.0)
              proj_test = true;
          }
          else if (dt == CORE::FE::CellType::tri3 || dt == CORE::FE::CellType::tri6)
          {
            if (mxi_test[0] >= 0.0 && mxi_test[1] >= 0.0 && mxi_test[0] <= 1.0 &&
                mxi_test[1] <= 1.0 && mxi_test[0] + mxi_test[1] <= 1.0)
              proj_test = true;
          }
          else
          {
            dserror("Non valid element type for master discretization!");
          }
        }
        if (proj_test == false) boundary_ele = true;
      }
    }

    else if (dt_s == CORE::FE::CellType::tri3 || dt_s == CORE::FE::CellType::tri6)
    {
      for (int s_test = 0; s_test < SlaveElement().NumNode(); ++s_test)
      {
        if (s_test == 0)
        {
          sxi_test[0] = 0.0;
          sxi_test[1] = 0.0;
        }
        else if (s_test == 1)
        {
          sxi_test[0] = 1.0;
          sxi_test[1] = 0.0;
        }
        else if (s_test == 2)
        {
          sxi_test[0] = 0.0;
          sxi_test[1] = 1.0;
        }
        else if (s_test == 3)
        {
          sxi_test[0] = 0.5;
          sxi_test[1] = 0.0;
        }
        else if (s_test == 4)
        {
          sxi_test[0] = 0.5;
          sxi_test[1] = 0.5;
        }
        else if (s_test == 5)
        {
          sxi_test[0] = 0.0;
          sxi_test[1] = 0.5;
        }

        proj_test = false;
        for (int bs_test = 0; bs_test < (int)Coupling().size(); ++bs_test)
        {
          double mxi_test[2] = {0.0, 0.0};
          MORTAR::Projector::Impl(SlaveElement(), Coupling()[bs_test]->MasterElement())
              ->ProjectGaussPoint3D(SlaveElement(), sxi_test, Coupling()[bs_test]->MasterElement(),
                  mxi_test, alpha_test);

          CORE::FE::CellType dt = Coupling()[bs_test]->MasterElement().Shape();
          if (dt == CORE::FE::CellType::quad4 || dt == CORE::FE::CellType::quad8 ||
              dt == CORE::FE::CellType::quad9)
          {
            if (mxi_test[0] >= -1.0 && mxi_test[1] >= -1.0 && mxi_test[0] <= 1.0 &&
                mxi_test[1] <= 1.0)
              proj_test = true;
          }
          else if (dt == CORE::FE::CellType::tri3 || dt == CORE::FE::CellType::tri6)
          {
            if (mxi_test[0] >= 0.0 && mxi_test[1] >= 0.0 && mxi_test[0] <= 1.0 &&
                mxi_test[1] <= 1.0 && mxi_test[0] + mxi_test[1] <= 1.0)
              proj_test = true;
          }
          else
          {
            dserror("Non valid element type for master discretization!");
          }
        }
        if (proj_test == false) boundary_ele = true;
      }
    }

    else
      dserror(
          "Calculation of consistent dual shape functions called for non-valid slave element "
          "shape!");

    if (boundary_ele == false) return;
  }

  // slave nodes and dofs
  const int max_nnodes = 9;
  const int nnodes = SlaveElement().NumNode();
  if (nnodes > max_nnodes)
    dserror(
        "this function is not implemented to handle elements with that many nodes. Just adjust "
        "max_nnodes above");
  const int ndof = 3;
  const int msize = MasterElements().size();

  // get number of master nodes
  int mnodes = 0;
  for (int m = 0; m < msize; ++m) mnodes += MasterElements()[m]->NumNode();

  // Dual shape functions coefficient matrix and linearization
  SlaveElement().MoData().DerivDualShape() =
      Teuchos::rcp(new CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>(
          (nnodes + mnodes) * ndof, 0, CORE::LINALG::SerialDenseMatrix(nnodes, nnodes)));
  CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& derivae =
      *(SlaveElement().MoData().DerivDualShape());

  // various variables
  double detg = 0.0;
  typedef CORE::GEN::pairedvector<int, double>::const_iterator _CI;

  // initialize matrices de and me
  CORE::LINALG::SerialDenseMatrix me(nnodes, nnodes, true);
  CORE::LINALG::SerialDenseMatrix de(nnodes, nnodes, true);

  CORE::GEN::pairedvector<int, CORE::LINALG::Matrix<max_nnodes + 1, max_nnodes>> derivde_new(
      (nnodes + mnodes) * ndof);

  // two-dim arrays of maps for linearization of me/de
  std::vector<std::vector<CORE::GEN::pairedvector<int, double>>> derivme(
      nnodes, std::vector<CORE::GEN::pairedvector<int, double>>(nnodes, (nnodes + mnodes) * ndof));
  std::vector<std::vector<CORE::GEN::pairedvector<int, double>>> derivde(
      nnodes, std::vector<CORE::GEN::pairedvector<int, double>>(nnodes, (nnodes + mnodes) * ndof));

  double A_tot = 0.;
  // loop over all master elements associated with this slave element
  for (int m = 0; m < (int)Coupling().size(); ++m)
  {
    if (!Coupling()[m]->RoughCheckCenters()) continue;
    if (!Coupling()[m]->RoughCheckOrient()) continue;
    if (!Coupling()[m]->RoughCheckCenters()) continue;

    // get number of master nodes
    const int ncol = Coupling()[m]->MasterElement().NumNode();

    // loop over all integration cells
    for (int c = 0; c < (int)Coupling()[m]->Cells().size(); ++c)
    {
      Teuchos::RCP<MORTAR::IntCell> currcell = Coupling()[m]->Cells()[c];

      A_tot += currcell->Area();

      // create an integrator for this cell
      CONTACT::Integrator integrator(imortar_, currcell->Shape(), Comm());

      // check if the cells are tri3
      // there's nothing wrong about other shapes, but as long as they are all
      // tri3 we can perform the jacobian calculation ( and its deriv) outside
      // the Gauss point loop
      if (currcell->Shape() != CORE::FE::CellType::tri3)
        dserror("only tri3 integration cells at the moment. See comment in the code");

      detg = currcell->Jacobian();
      // directional derivative of cell Jacobian
      CORE::GEN::pairedvector<int, double> derivjaccell((nnodes + ncol) * ndof);
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
        MORTAR::Projector::Impl(Coupling()[m]->SlaveIntElement())
            ->ProjectGaussPointAuxn3D(
                globgp, Coupling()[m]->Auxn(), Coupling()[m]->SlaveIntElement(), sxi, sprojalpha);

        // project Gauss point onto slave (parent) element
        double psxi[2] = {0., 0.};
        double psprojalpha = 0.0;
        if (Quad())
        {
          MORTAR::IntElement* ie =
              dynamic_cast<MORTAR::IntElement*>(&(Coupling()[m]->SlaveIntElement()));
          if (ie == nullptr) dserror("nullptr pointer");
          MORTAR::Projector::Impl(SlaveElement())
              ->ProjectGaussPointAuxn3D(
                  globgp, Coupling()[m]->Auxn(), SlaveElement(), psxi, psprojalpha);
          // ie->MapToParent(sxi,psxi); // old way of doing it via affine map... wrong (popp
          // 05/2016)
        }
        else
          for (int i = 0; i < 2; ++i) psxi[i] = sxi[i];

        // create vector for shape function evaluation
        CORE::LINALG::SerialDenseVector sval(nnodes);
        CORE::LINALG::SerialDenseMatrix sderiv(nnodes, 2, true);

        // evaluate trace space shape functions at Gauss point
        if (LagMultQuad() == INPAR::MORTAR::lagmult_lin)
          SlaveElement().EvaluateShapeLagMultLin(
              INPAR::MORTAR::shape_standard, psxi, sval, sderiv, nnodes);
        else
          SlaveElement().EvaluateShape(psxi, sval, sderiv, nnodes);

        // additional data for contact calculation (i.e. incl. derivative of dual shape functions
        // coefficient matrix) GP slave coordinate derivatives
        std::vector<CORE::GEN::pairedvector<int, double>> dsxigp(2, (nnodes + ncol) * ndof);
        // GP slave coordinate derivatives
        std::vector<CORE::GEN::pairedvector<int, double>> dpsxigp(2, (nnodes + ncol) * ndof);
        // global GP coordinate derivative on integration element
        CORE::GEN::pairedvector<int, CORE::LINALG::Matrix<3, 1>> lingp((nnodes + ncol) * ndof);

        // compute global GP coordinate derivative
        static CORE::LINALG::Matrix<3, 1> svalcell;
        static CORE::LINALG::Matrix<3, 2> sderivcell;
        currcell->EvaluateShape(eta, svalcell, sderivcell);

        for (int v = 0; v < 3; ++v)
          for (int d = 0; d < 3; ++d)
            for (_CI p = (currcell->GetDerivVertex(v))[d].begin();
                 p != (currcell->GetDerivVertex(v))[d].end(); ++p)
              lingp[p->first](d) += svalcell(v) * (p->second);

        // compute GP slave coordinate derivatives
        integrator.DerivXiGP3DAuxPlane(Coupling()[m]->SlaveIntElement(), sxi, currcell->Auxn(),
            dsxigp, sprojalpha, currcell->GetDerivAuxn(), lingp);

        // compute GP slave coordinate derivatives (parent element)
        if (Quad())
        {
          MORTAR::IntElement* ie =
              dynamic_cast<MORTAR::IntElement*>(&(Coupling()[m]->SlaveIntElement()));
          if (ie == nullptr) dserror("wtf");
          integrator.DerivXiGP3DAuxPlane(SlaveElement(), psxi, currcell->Auxn(), dpsxigp,
              psprojalpha, currcell->GetDerivAuxn(), lingp);
          // ie->MapToParent(dsxigp,dpsxigp); // old way of doing it via affine map... wrong (popp
          // 05/2016)
        }
        else
          dpsxigp = dsxigp;

        double fac = 0.;
        for (_CI p = derivjaccell.begin(); p != derivjaccell.end(); ++p)
        {
          CORE::LINALG::Matrix<max_nnodes + 1, max_nnodes>& dtmp = derivde_new[p->first];
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
            CORE::LINALG::Matrix<max_nnodes + 1, max_nnodes>& dtmp = derivde_new[p->first];
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
    }  // cells
  }    // master elements

  // in case of no overlap just return, as there is no integration area
  // and therefore the consistent dual shape functions are not defined.
  // This doesn't matter, as there is no associated integration domain anyway
  if (A_tot < 1.e-12) return;

  // declare dual shape functions coefficient matrix and
  // inverse of matrix M_e
  CORE::LINALG::SerialDenseMatrix ae(nnodes, nnodes, true);
  CORE::LINALG::SerialDenseMatrix meinv(nnodes, nnodes, true);

  // compute matrix A_e and inverse of matrix M_e for
  // linear interpolation of quadratic element
  if (LagMultQuad() == INPAR::MORTAR::lagmult_lin)
  {
    // declare and initialize to zero inverse of Matrix M_e
    CORE::LINALG::SerialDenseMatrix meinv(nnodes, nnodes, true);

    if (SlaveElement().Shape() == CORE::FE::CellType::tri6)
    {
      // reduce me to non-zero nodes before inverting
      CORE::LINALG::Matrix<3, 3> melin;
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k) melin(j, k) = me(j, k);

      // invert bi-ortho matrix melin
      CORE::LINALG::Inverse(melin);

      // re-inflate inverse of melin to full size
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k) meinv(j, k) = melin(j, k);
    }
    else if (SlaveElement().Shape() == CORE::FE::CellType::quad8 ||
             SlaveElement().Shape() == CORE::FE::CellType::quad9)
    {
      // reduce me to non-zero nodes before inverting
      CORE::LINALG::Matrix<4, 4> melin;
      for (int j = 0; j < 4; ++j)
        for (int k = 0; k < 4; ++k) melin(j, k) = me(j, k);

      // invert bi-ortho matrix melin
      CORE::LINALG::Inverse(melin);

      // re-inflate inverse of melin to full size
      for (int j = 0; j < 4; ++j)
        for (int k = 0; k < 4; ++k) meinv(j, k) = melin(j, k);
    }
    else
      dserror("incorrect element shape for linear interpolation of quadratic element!");

    // get solution matrix with dual parameters
    CORE::LINALG::multiply(ae, de, meinv);
  }
  // compute matrix A_e and inverse of matrix M_e for all other cases
  else
    meinv = CORE::LINALG::InvertAndMultiplyByCholesky(me, de, ae);

  // build linearization of ae and store in derivdual
  // (this is done according to a quite complex formula, which
  // we get from the linearization of the biorthogonality condition:
  // Lin (Me * Ae = De) -> Lin(Ae)=Lin(De)*Inv(Me)-Ae*Lin(Me)*Inv(Me) )
  typedef CORE::GEN::pairedvector<int,
      CORE::LINALG::Matrix<max_nnodes + 1, max_nnodes>>::const_iterator _CIM;
  for (_CIM p = derivde_new.begin(); p != derivde_new.end(); ++p)
  {
    CORE::LINALG::Matrix<max_nnodes + 1, max_nnodes>& dtmp = derivde_new[p->first];
    CORE::LINALG::SerialDenseMatrix& pt = derivae[p->first];
    for (int i = 0; i < nnodes; ++i)
      for (int j = 0; j < nnodes; ++j)
      {
        pt(i, j) += meinv(i, j) * dtmp(nnodes, i);

        for (int k = 0; k < nnodes; ++k)
          for (int l = 0; l < nnodes; ++l) pt(i, j) -= ae(i, k) * meinv(l, j) * dtmp(l, k);
      }
  }

  // store ae matrix in slave element data container
  SlaveElement().MoData().DualShape() = Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix(ae));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Coupling3dManager::FindFeasibleMasterElements(
    std::vector<MORTAR::Element*>& feasible_ma_eles) const
{
  // feasibility counter
  std::size_t fcount = 0;
  for (std::size_t m = 0; m < MasterElements().size(); ++m)
  {
    // Build a instance of the MORTAR::Coupling3d object (no linearization needed).
    Teuchos::RCP<MORTAR::Coupling3d> coup = Teuchos::rcp(
        new MORTAR::Coupling3d(idiscret_, dim_, false, imortar_, SlaveElement(), MasterElement(m)));

    // Building the master element normals and check the angles.
    if (coup->RoughCheckOrient())
    {
      feasible_ma_eles[fcount] = &MasterElement(m);
      ++fcount;
    }
  }
  feasible_ma_eles.resize(fcount);

  return;
}

BACI_NAMESPACE_CLOSE
