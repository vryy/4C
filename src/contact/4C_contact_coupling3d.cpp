/*---------------------------------------------------------------------*/
/*! \file
\brief Classes for mortar contact coupling in 3D.

\level 2


*/
/*---------------------------------------------------------------------*/

#include "4C_contact_coupling3d.hpp"

#include "4C_contact_defines.hpp"
#include "4C_contact_element.hpp"
#include "4C_contact_integrator.hpp"
#include "4C_contact_integrator_factory.hpp"
#include "4C_contact_interpolator.hpp"
#include "4C_contact_node.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_densematrix_inverse.hpp"
#include "4C_linalg_utils_densematrix_multiply.hpp"
#include "4C_mortar_calc_utils.hpp"
#include "4C_mortar_coupling3d_classes.hpp"
#include "4C_mortar_defines.hpp"
#include "4C_mortar_projector.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::Coupling3d::Coupling3d(Discret::Discretization& idiscret, int dim, bool quad,
    Teuchos::ParameterList& params, Mortar::Element& sele, Mortar::Element& mele)
    : Mortar::Coupling3d(idiscret, dim, quad, params, sele, mele),
      stype_(Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(params, "STRATEGY"))
{
  // empty constructor

  return;
}

/*----------------------------------------------------------------------*
 |  Build auxiliary plane from slave element (public)         popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::auxiliary_plane()
{
  // we first need the element center:
  // for quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  double loccenter[2];

  Core::FE::CellType dt = SlaveIntElement().Shape();
  if (dt == Core::FE::CellType::tri3 || dt == Core::FE::CellType::tri6)
  {
    loccenter[0] = 1.0 / 3.0;
    loccenter[1] = 1.0 / 3.0;
  }
  else if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
           dt == Core::FE::CellType::quad9)
  {
    loccenter[0] = 0.0;
    loccenter[1] = 0.0;
  }
  else
    FOUR_C_THROW("auxiliary_plane called for unknown element type");

  // compute element center via shape fct. interpolation
  SlaveIntElement().LocalToGlobal(loccenter, Auxc(), 0);

  // we then compute the unit normal vector at the element center
  Lauxn() = SlaveIntElement().compute_unit_normal_at_xi(loccenter, Auxn());

  // THIS IS CONTACT-SPECIFIC!
  // also compute linearization of the unit normal vector
  SlaveIntElement().DerivUnitNormalAtXi(loccenter, get_deriv_auxn());

  // std::cout << "Slave Element: " << SlaveIntElement().Id() << std::endl;
  // std::cout << "->Center: " << Auxc()[0] << " " << Auxc()[1] << " " << Auxc()[2] << std::endl;
  // std::cout << "->Normal: " << Auxn()[0] << " " << Auxn()[1] << " " << Auxn()[2] << std::endl;

  return true;
}

/*----------------------------------------------------------------------*
 |  Integration of cells (3D)                                 popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::IntegrateCells(const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr)
{
  /**********************************************************************/
  /* INTEGRATION                                                        */
  /* Integrate the Mortar matrix M and the weighted gap function g~ on  */
  /* the current integration cell of the slave / master element pair    */
  /**********************************************************************/

  static const Inpar::Mortar::AlgorithmType algo =
      Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(imortar_, "ALGORITHM");

  // do nothing if there are no cells
  if (Cells().size() == 0) return false;

  // create a CONTACT integrator instance with correct num_gp and Dim
  // it is sufficient to do this once as all IntCells are triangles
  Teuchos::RCP<CONTACT::Integrator> integrator =
      CONTACT::INTEGRATOR::BuildIntegrator(stype_, imortar_, Cells()[0]->Shape(), Comm());
  // loop over all integration cells
  for (int i = 0; i < (int)(Cells().size()); ++i)
  {
    // integrate cell only if it has a non-zero area
    if (Cells()[i]->Area() < MORTARINTLIM * slave_element_area()) continue;

      // debug output of integration cells in GMSH
#ifdef MORTARGMSHCELLS
    GmshOutputCells(i);
#endif
    // set segmentation status of all slave nodes
    // (hassegment_ of a slave node is true if ANY segment/cell
    // is integrated that contributes to this slave node)
    int nnodes = SlaveIntElement().num_node();
    Core::Nodes::Node** mynodes = SlaveIntElement().Nodes();
    if (!mynodes) FOUR_C_THROW("Null pointer!");
    for (int k = 0; k < nnodes; ++k)
    {
      Mortar::Node* mycnode = dynamic_cast<Mortar::Node*>(mynodes[k]);
      if (!mycnode) FOUR_C_THROW("Null pointer!");
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
    Inpar::Mortar::LagMultQuad lmtype = LagMultQuad();

    // *******************************************************************
    // case (1)
    // *******************************************************************
    if (!Quad())
    {
      integrator->integrate_deriv_cell3_d_aux_plane(
          SlaveElement(), MasterElement(), Cells()[i], Auxn(), Comm(), mparams_ptr);
    }
    // *******************************************************************
    // cases (2) and (3)
    // *******************************************************************
    else if ((Quad() and
                 (lmtype == Inpar::Mortar::lagmult_quad or lmtype == Inpar::Mortar::lagmult_lin or
                     lmtype == Inpar::Mortar::lagmult_const)) or
             algo == Inpar::Mortar::algorithm_gpts)
    {
      // check for standard shape functions and quadratic LM interpolation
      if (ShapeFcn() == Inpar::Mortar::shape_standard && lmtype == Inpar::Mortar::lagmult_quad &&
          (SlaveElement().Shape() == Core::FE::CellType::quad8 ||
              SlaveElement().Shape() == Core::FE::CellType::tri6))
        FOUR_C_THROW(
            "Quad. LM interpolation for STANDARD 3D quadratic contact only feasible for "
            "quad9");

      // dynamic_cast to make sure to pass in IntElement&
      Mortar::IntElement& sintref = dynamic_cast<Mortar::IntElement&>(SlaveIntElement());
      Mortar::IntElement& mintref = dynamic_cast<Mortar::IntElement&>(MasterIntElement());

      // call integrator
      integrator->integrate_deriv_cell3_d_aux_plane_quad(
          SlaveElement(), MasterElement(), sintref, mintref, Cells()[i], Auxn());
    }

    // *******************************************************************
    // case (4)
    // *******************************************************************
    else if (Quad() && lmtype == Inpar::Mortar::lagmult_pwlin)
    {
      // check for dual shape functions
      if (ShapeFcn() == Inpar::Mortar::shape_dual ||
          ShapeFcn() == Inpar::Mortar::shape_petrovgalerkin)
        FOUR_C_THROW(
            "Piecewise linear LM interpolation not yet implemented for DUAL 3D quadratic "
            "contact");

      // dynamic_cast to make sure to pass in IntElement&
      Mortar::IntElement& sintref = dynamic_cast<Mortar::IntElement&>(SlaveIntElement());
      Mortar::IntElement& mintref = dynamic_cast<Mortar::IntElement&>(MasterIntElement());

      // call integrator
      integrator->integrate_deriv_cell3_d_aux_plane_quad(
          SlaveElement(), MasterElement(), sintref, mintref, Cells()[i], Auxn());
    }

    // *******************************************************************
    // undefined case
    // *******************************************************************
    else if (Quad() && lmtype == Inpar::Mortar::lagmult_undefined)
    {
      FOUR_C_THROW(
          "Lagrange multiplier interpolation for quadratic elements undefined\n"
          "If you are using 2nd order mortar elements, you need to specify LM_QUAD in MORTAR "
          "COUPLING section");
    }

    // *******************************************************************
    // other cases
    // *******************************************************************
    else
      FOUR_C_THROW("IntegrateCells: Invalid case for 3D mortar contact LM interpolation");
    // *******************************************************************
  }  // cell loop

  return true;
}

/*----------------------------------------------------------------------*
 |  Linearization of clip polygon vertices (3D)               popp 02/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::VertexLinearization(
    std::vector<std::vector<Core::Gen::Pairedvector<int, double>>>& linvertex,
    std::map<int, double>& projpar, bool printderiv)
{
  // linearize all aux.plane slave and master nodes only ONCE
  // and use these linearizations later during lineclip linearization
  // (this speeds up the vertex linearizations in most cases, as we
  // never linearize the SAME slave or master vertex more than once)

  // number of nodes
  const int nsrows = SlaveIntElement().num_node();
  const int nmrows = MasterIntElement().num_node();

  // prepare storage for slave and master linearizations
  std::vector<std::vector<Core::Gen::Pairedvector<int, double>>> linsnodes(
      nsrows, std::vector<Core::Gen::Pairedvector<int, double>>(3, 3 * SlaveElement().num_node()));
  std::vector<std::vector<Core::Gen::Pairedvector<int, double>>> linmnodes(
      nmrows, std::vector<Core::Gen::Pairedvector<int, double>>(
                  3, 3 * SlaveElement().num_node() + 3 * MasterElement().num_node()));

  // compute slave linearizations (nsrows)
  slave_vertex_linearization(linsnodes);

  // compute master linearizations (nmrows)
  master_vertex_linearization(linmnodes);

  //**********************************************************************
  // Clip polygon vertex linearization
  //**********************************************************************
  // loop over all clip polygon vertices
  for (int i = 0; i < (int)Clip().size(); ++i)
  {
    // references to current vertex and its linearization
    Mortar::Vertex& currv = Clip()[i];
    std::vector<Core::Gen::Pairedvector<int, double>>& currlin = linvertex[i];

    // decision on vertex type (slave, projmaster, linclip)
    if (currv.v_type() == Mortar::Vertex::slave)
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

      // FOUR_C_THROW if not found
      if (k == nsrows) FOUR_C_THROW("Slave Id not found!");

      // get the correct slave node linearization
      currlin = linsnodes[k];
    }
    else if (currv.v_type() == Mortar::Vertex::projmaster)
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

      // FOUR_C_THROW if not found
      if (k == nmrows) FOUR_C_THROW("Master Id not found!");

      // get the correct master node linearization
      currlin = linmnodes[k];
    }
    else if (currv.v_type() == Mortar::Vertex::lineclip)
    {
      // get references to the two slave vertices
      int sindex1 = -1;
      int sindex2 = -1;
      for (int j = 0; j < (int)slave_vertices().size(); ++j)
      {
        if (slave_vertices()[j].Nodeids()[0] == currv.Nodeids()[0]) sindex1 = j;
        if (slave_vertices()[j].Nodeids()[0] == currv.Nodeids()[1]) sindex2 = j;
      }
      if (sindex1 < 0 || sindex2 < 0 || sindex1 == sindex2)
        FOUR_C_THROW("Lineclip linearization: (S) Something went wrong!");

      Mortar::Vertex* sv1 = &slave_vertices()[sindex1];
      Mortar::Vertex* sv2 = &slave_vertices()[sindex2];

      // get references to the two master vertices
      int mindex1 = -1;
      int mindex2 = -1;
      for (int j = 0; j < (int)master_vertices().size(); ++j)
      {
        if (master_vertices()[j].Nodeids()[0] == currv.Nodeids()[2]) mindex1 = j;
        if (master_vertices()[j].Nodeids()[0] == currv.Nodeids()[3]) mindex2 = j;
      }
      if (mindex1 < 0 || mindex2 < 0 || mindex1 == mindex2)
        FOUR_C_THROW("Lineclip linearization: (M) Something went wrong!");

      Mortar::Vertex* mv1 = &master_vertices()[mindex1];
      Mortar::Vertex* mv2 = &master_vertices()[mindex2];

      // do lineclip vertex linearization
      lineclip_vertex_linearization(currv, currlin, sv1, sv2, mv1, mv2, linsnodes, linmnodes);
    }

    else
      FOUR_C_THROW("VertexLinearization: Invalid Vertex Type!");
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Linearization of slave vertex (3D) AuxPlane               popp 03/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::slave_vertex_linearization(
    std::vector<std::vector<Core::Gen::Pairedvector<int, double>>>& currlin)
{
  // we first need the slave element center:
  // for quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  double scxi[2];

  Core::FE::CellType dt = SlaveIntElement().Shape();
  if (dt == Core::FE::CellType::tri3 || dt == Core::FE::CellType::tri6)
  {
    scxi[0] = 1.0 / 3.0;
    scxi[1] = 1.0 / 3.0;
  }
  else if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
           dt == Core::FE::CellType::quad9)
  {
    scxi[0] = 0.0;
    scxi[1] = 0.0;
  }
  else
    FOUR_C_THROW("slave_vertex_linearization called for unknown element type");

  // evlauate shape functions + derivatives at scxi
  const int nrow = SlaveIntElement().num_node();
  Core::LinAlg::SerialDenseVector sval(nrow);
  Core::LinAlg::SerialDenseMatrix sderiv(nrow, 2, true);
  SlaveIntElement().evaluate_shape(scxi, sval, sderiv, nrow);

  // we need all participating slave nodes
  Core::Nodes::Node** snodes = SlaveIntElement().Nodes();
  std::vector<Mortar::Node*> smrtrnodes(nrow);

  for (int i = 0; i < nrow; ++i)
  {
    smrtrnodes[i] = dynamic_cast<Mortar::Node*>(snodes[i]);
    if (!smrtrnodes[i]) FOUR_C_THROW("slave_vertex_linearization: Null pointer!");
  }

  // linearization of the IntEle spatial coords
  std::vector<std::vector<Core::Gen::Pairedvector<int, double>>> nodelin(0);
  Mortar::IntElement* sIntEle = dynamic_cast<Mortar::IntElement*>(&SlaveIntElement());

  if (sIntEle == nullptr)
  {
    // resize the linearizations
    nodelin.resize(nrow, std::vector<Core::Gen::Pairedvector<int, double>>(3, 1));

    // loop over all intEle nodes
    for (int in = 0; in < nrow; ++in)
      for (int dim = 0; dim < 3; ++dim) nodelin[in][dim][smrtrnodes[in]->Dofs()[dim]] += 1.;
  }
  else
    sIntEle->NodeLinearization(nodelin);

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator
      _CI;  // linearization of element center Auxc()
  std ::vector<Core::Gen::Pairedvector<int, double>> linauxc(
      3, SlaveElement().num_node());  // assume 3 dofs per node

  for (int i = 0; i < nrow; ++i)
    for (int dim = 0; dim < 3; ++dim)
      for (_CI p = nodelin[i][dim].begin(); p != nodelin[i][dim].end(); ++p)
        linauxc[dim][p->first] = sval[i] * p->second;

  // linearization of element normal Auxn()
  std::vector<Core::Gen::Pairedvector<int, double>>& linauxn = get_deriv_auxn();

  // put everything together for slave vertex linearization
  // loop over all vertices
  for (int i = 0; i < SlaveIntElement().num_node(); ++i)
  {
    Mortar::Node* mrtrsnode = dynamic_cast<Mortar::Node*>(SlaveIntElement().Nodes()[i]);
    if (!mrtrsnode) FOUR_C_THROW("cast to mortar node failed");

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
bool CONTACT::Coupling3d::master_vertex_linearization(
    std::vector<std::vector<Core::Gen::Pairedvector<int, double>>>& currlin)
{
  // we first need the slave element center:
  // for quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  double scxi[2];

  Core::FE::CellType dt = SlaveIntElement().Shape();
  if (dt == Core::FE::CellType::tri3 || dt == Core::FE::CellType::tri6)
  {
    scxi[0] = 1.0 / 3.0;
    scxi[1] = 1.0 / 3.0;
  }
  else if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
           dt == Core::FE::CellType::quad9)
  {
    scxi[0] = 0.0;
    scxi[1] = 0.0;
  }
  else
    FOUR_C_THROW("master_vertex_linearization called for unknown element type");

  // evlauate shape functions + derivatives at scxi
  int nrow = SlaveIntElement().num_node();
  Core::LinAlg::SerialDenseVector sval(nrow);
  Core::LinAlg::SerialDenseMatrix sderiv(nrow, 2, true);
  SlaveIntElement().evaluate_shape(scxi, sval, sderiv, nrow);

  // we need all participating slave nodes
  Core::Nodes::Node** snodes = SlaveIntElement().Nodes();
  std::vector<Mortar::Node*> smrtrnodes(nrow);

  for (int i = 0; i < nrow; ++i)
  {
    smrtrnodes[i] = dynamic_cast<Mortar::Node*>(snodes[i]);
    if (!smrtrnodes[i]) FOUR_C_THROW("master_vertex_linearization: Null pointer!");
  }

  // linearization of the SlaveIntEle spatial coords
  std::vector<std::vector<Core::Gen::Pairedvector<int, double>>> snodelin(0);
  Mortar::IntElement* sIntEle = dynamic_cast<Mortar::IntElement*>(&SlaveIntElement());

  if (sIntEle == nullptr)
  {
    // resize the linearizations
    snodelin.resize(nrow, std::vector<Core::Gen::Pairedvector<int, double>>(3, 1));

    // loop over all intEle nodes
    for (int in = 0; in < nrow; ++in)
      for (int dim = 0; dim < 3; ++dim) snodelin[in][dim][smrtrnodes[in]->Dofs()[dim]] += 1.;
  }
  else
    sIntEle->NodeLinearization(snodelin);

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator
      _CI;  // linearization of element center Auxc()
  std ::vector<Core::Gen::Pairedvector<int, double>> linauxc(
      3, SlaveElement().num_node());  // assume 3 dofs per node

  for (int i = 0; i < nrow; ++i)
    for (int dim = 0; dim < 3; ++dim)
      for (_CI p = snodelin[i][dim].begin(); p != snodelin[i][dim].end(); ++p)
        linauxc[dim][p->first] = sval[i] * p->second;

  // linearization of element normal Auxn()
  std::vector<Core::Gen::Pairedvector<int, double>>& linauxn = get_deriv_auxn();

  // linearization of the MasterIntEle spatial coords
  std::vector<std::vector<Core::Gen::Pairedvector<int, double>>> mnodelin(0);
  Mortar::IntElement* mIntEle = dynamic_cast<Mortar::IntElement*>(&MasterIntElement());

  if (mIntEle == nullptr)
  {
    // resize the linearizations
    mnodelin.resize(
        MasterIntElement().num_node(), std::vector<Core::Gen::Pairedvector<int, double>>(3, 1));

    // loop over all intEle nodes
    for (int in = 0; in < MasterIntElement().num_node(); ++in)
    {
      Mortar::Node* mrtrmnode = dynamic_cast<Mortar::Node*>(MasterIntElement().Nodes()[in]);
      if (mrtrmnode == nullptr) FOUR_C_THROW("dynamic cast to mortar node went wrong");

      for (int dim = 0; dim < 3; ++dim) mnodelin[in][dim][mrtrmnode->Dofs()[dim]] += 1.;
    }
  }
  else
    mIntEle->NodeLinearization(mnodelin);

  // put everything together for slave vertex linearization
  // loop over all vertices
  for (int i = 0; i < MasterIntElement().num_node(); ++i)
  {
    Mortar::Node* mrtrmnode = dynamic_cast<Mortar::Node*>(MasterIntElement().Nodes()[i]);
    if (!mrtrmnode) FOUR_C_THROW("cast to mortar node failed");

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
bool CONTACT::Coupling3d::lineclip_vertex_linearization(Mortar::Vertex& currv,
    std::vector<Core::Gen::Pairedvector<int, double>>& currlin, Mortar::Vertex* sv1,
    Mortar::Vertex* sv2, Mortar::Vertex* mv1, Mortar::Vertex* mv2,
    std::vector<std::vector<Core::Gen::Pairedvector<int, double>>>& linsnodes,
    std::vector<std::vector<Core::Gen::Pairedvector<int, double>>>& linmnodes)
{
  // number of nodes
  const int nsrows = SlaveIntElement().num_node();
  const int nmrows = MasterIntElement().num_node();

  // iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

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

  // FOUR_C_THROW if not found
  if (k == nsrows) FOUR_C_THROW("Slave Id1 not found!");

  // get the correct slave node linearization
  std::vector<Core::Gen::Pairedvector<int, double>>& slavelin0 = linsnodes[k];

  k = 0;
  while (k < nsrows)
  {
    if (SlaveIntElement().NodeIds()[k] == sid2) break;
    ++k;
  }

  // FOUR_C_THROW if not found
  if (k == nsrows) FOUR_C_THROW("Slave Id2 not found!");

  // get the correct slave node linearization
  std::vector<Core::Gen::Pairedvector<int, double>>& slavelin1 = linsnodes[k];

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

  // FOUR_C_THROW if not found
  if (k == nmrows) FOUR_C_THROW("Master Id1 not found!");

  // get the correct master node linearization
  std::vector<Core::Gen::Pairedvector<int, double>>& masterlin0 = linmnodes[k];

  k = 0;
  while (k < nmrows)
  {
    if (MasterIntElement().NodeIds()[k] == mid2) break;
    ++k;
  }

  // FOUR_C_THROW if not found
  if (k == nmrows) FOUR_C_THROW("Master Id2 not found!");

  // get the correct master node linearization
  std::vector<Core::Gen::Pairedvector<int, double>>& masterlin1 = linmnodes[k];

  // linearization of element normal Auxn()
  std::vector<Core::Gen::Pairedvector<int, double>>& linauxn = get_deriv_auxn();

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
    const std::vector<std::vector<Core::Gen::Pairedvector<int, double>>>& linvertex,
    std::vector<Core::Gen::Pairedvector<int, double>>& lincenter)
{
  // preparations
  int clipsize = (int)(Clip().size());
  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;

  // number of nodes
  const int nsrows = SlaveElement().num_node();
  const int nmrows = MasterElement().num_node();

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
  std::vector<Core::Gen::Pairedvector<int, double>> linnac(3, 3 * (nsrows + nmrows));
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
    std::vector<Core::Gen::Pairedvector<int, double>> lincross(3, 3 * (nsrows + nmrows));

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
    Core::Gen::Pairedvector<int, double> linarea(3 * (nsrows + nmrows));
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
CONTACT::Coupling3dQuad::Coupling3dQuad(Discret::Discretization& idiscret, int dim, bool quad,
    Teuchos::ParameterList& params, Mortar::Element& sele, Mortar::Element& mele,
    Mortar::IntElement& sintele, Mortar::IntElement& mintele)
    : CONTACT::Coupling3d(idiscret, dim, quad, params, sele, mele),
      sintele_(sintele),
      mintele_(mintele)
{
  //  3D quadratic coupling only for quadratic ansatz type
  if (!Quad()) FOUR_C_THROW("Coupling3dQuad called for non-quadratic ansatz!");

  return;
}

/*----------------------------------------------------------------------*
 |  get communicator  (public)                               farah 01/13|
 *----------------------------------------------------------------------*/
const Epetra_Comm& CONTACT::Coupling3dManager::Comm() const { return idiscret_.Comm(); }

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::Coupling3dManager::Coupling3dManager(Discret::Discretization& idiscret, int dim, bool quad,
    Teuchos::ParameterList& params, Mortar::Element* sele, std::vector<Mortar::Element*> mele)
    : idiscret_(idiscret),
      dim_(dim),
      quad_(quad),
      imortar_(params),
      sele_(sele),
      mele_(mele),
      ncells_(0),
      stype_(Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(params, "STRATEGY"))
{
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 01/13|
 *----------------------------------------------------------------------*/
CONTACT::Coupling3dQuadManager::Coupling3dQuadManager(Discret::Discretization& idiscret, int dim,
    bool quad, Teuchos::ParameterList& params, Mortar::Element* sele,
    std::vector<Mortar::Element*> mele)
    : Mortar::Coupling3dQuadManager(idiscret, dim, quad, params, sele, mele),
      CONTACT::Coupling3dManager(idiscret, dim, quad, params, sele, mele),
      smintpairs_(-1),
      intcells_(-1)
{
  return;
}


/*----------------------------------------------------------------------*
 |  Evaluate mortar coupling pairs                            popp 03/09|
 *----------------------------------------------------------------------*/
void CONTACT::Coupling3dManager::integrate_coupling(
    const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr)
{
  // get algorithm
  Inpar::Mortar::AlgorithmType algo =
      Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(imortar_, "ALGORITHM");

  // prepare linearizations
  if (algo == Inpar::Mortar::algorithm_mortar)
    dynamic_cast<CONTACT::Element&>(SlaveElement()).PrepareDderiv(MasterElements());

  // decide which type of numerical integration scheme

  //**********************************************************************
  // STANDARD INTEGRATION (SEGMENTS)
  //**********************************************************************
  if (IntType() == Inpar::Mortar::inttype_segments)
  {
    // loop over all master elements associated with this slave element
    for (int m = 0; m < (int)MasterElements().size(); ++m)
    {
      // create Coupling3d object and push back
      Coupling().push_back(Teuchos::rcp(
          new Coupling3d(idiscret_, dim_, false, imortar_, SlaveElement(), MasterElement(m))));

      // do coupling
      Coupling()[m]->evaluate_coupling();

      // store number of intcells
      ncells_ += (int)(Coupling()[m]->Cells()).size();
    }

    // special treatment of boundary elements
    consist_dual_shape();

    // TODO modification of boundary shapes!

    // integrate cells
    for (int i = 0; i < (int)Coupling().size(); ++i)
    {
      // temporary m-matrix linearization of this slave/master pair
      if (algo == Inpar::Mortar::algorithm_mortar)
        dynamic_cast<CONTACT::Element&>(SlaveElement()).PrepareMderiv(MasterElements(), i);

      // integrate cells
      Coupling()[i]->IntegrateCells(mparams_ptr);

      // assemble m-matrix for this slave/master pair
      if (algo == Inpar::Mortar::algorithm_mortar)
        dynamic_cast<CONTACT::Element&>(SlaveElement())
            .assemble_mderiv_to_nodes(Coupling()[i]->MasterElement());
    }
  }

  //**********************************************************************
  // ELEMENT-BASED INTEGRATION
  //**********************************************************************
  else if (IntType() == Inpar::Mortar::inttype_elements ||
           IntType() == Inpar::Mortar::inttype_elements_BS)
  {
    if ((int)MasterElements().size() == 0) return;

    if (!Quad())
    {
      bool boundary_ele = false;
      bool proj = false;

      /* find all feasible master elements (this check is inherent in the
       * segment based integration)                    hiermeier 04/16 */
      std::vector<Mortar::Element*> feasible_ma_eles(MasterElements().size());
      find_feasible_master_elements(feasible_ma_eles);

      // create an integrator instance with correct num_gp and Dim
      Teuchos::RCP<CONTACT::Integrator> integrator =
          CONTACT::INTEGRATOR::BuildIntegrator(stype_, imortar_, SlaveElement().Shape(), Comm());

      // Perform integration and linearization
      integrator->IntegrateDerivEle3D(
          SlaveElement(), feasible_ma_eles, &boundary_ele, &proj, Comm(), mparams_ptr);


      if (IntType() == Inpar::Mortar::inttype_elements_BS)
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
            Coupling()[m]->evaluate_coupling();

            // store number of intcells
            ncells_ += (int)(Coupling()[m]->Cells()).size();
          }

          // special treatment of boundary elements
          consist_dual_shape();

          // integrate cells
          for (int i = 0; i < (int)Coupling().size(); ++i)
          {
            // temporary m-matrix linearization of this slave/master pair
            if (algo == Inpar::Mortar::algorithm_mortar)
              dynamic_cast<CONTACT::Element&>(SlaveElement()).PrepareMderiv(MasterElements(), i);

            // integrate cells
            Coupling()[i]->IntegrateCells(mparams_ptr);

            // assemble m-matrix for this slave/master pair
            if (algo == Inpar::Mortar::algorithm_mortar)
              dynamic_cast<CONTACT::Element&>(SlaveElement())
                  .assemble_mderiv_to_nodes(Coupling()[i]->MasterElement());
          }
        }
      }
    }
    else
    {
      FOUR_C_THROW(
          "You should not be here! This coupling manager is not able to perform mortar coupling "
          "for high-order elements.");
    }
  }
  //**********************************************************************
  // INVALID TYPE OF NUMERICAL INTEGRATION
  //**********************************************************************
  else
  {
    FOUR_C_THROW("Invalid type of numerical integration!");
  }

  // free memory of dual shape function coefficient matrix
  SlaveElement().MoData().ResetDualShape();
  SlaveElement().MoData().ResetDerivDualShape();

  // assemble element contribution to nodes
  if (algo == Inpar::Mortar::algorithm_mortar)
  {
    bool dual = (ShapeFcn() == Inpar::Mortar::shape_dual) ||
                (ShapeFcn() == Inpar::Mortar::shape_petrovgalerkin);
    dynamic_cast<CONTACT::Element&>(SlaveElement()).assemble_dderiv_to_nodes(dual);
  }
  return;
}


/*----------------------------------------------------------------------*
 |  Evaluate coupling pairs                                 farah 09/14 |
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3dManager::evaluate_coupling(
    const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr)
{
  // decide which type of coupling should be evaluated
  Inpar::Mortar::AlgorithmType algo =
      Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(imortar_, "ALGORITHM");

  //*********************************
  // Mortar Contact
  //*********************************
  if (algo == Inpar::Mortar::algorithm_mortar || algo == Inpar::Mortar::algorithm_gpts)
    integrate_coupling(mparams_ptr);

  //*********************************
  // Error
  //*********************************
  else
    FOUR_C_THROW("chose contact algorithm not supported!");

  // interpolate temperatures in TSI case
  if (imortar_.get<int>("PROBTYPE") == Inpar::CONTACT::tsi)
    NTS::Interpolator(imortar_, dim_).interpolate_master_temp3_d(SlaveElement(), MasterElements());

  return true;
}


/*----------------------------------------------------------------------*
 |  Evaluate mortar coupling pairs for Quad-coupling         farah 09/14|
 *----------------------------------------------------------------------*/
void CONTACT::Coupling3dQuadManager::integrate_coupling(
    const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr)
{
  // get algorithm type
  Inpar::Mortar::AlgorithmType algo = Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(
      Mortar::Coupling3dQuadManager::imortar_, "ALGORITHM");

  // prepare linearizations
  if (algo == Inpar::Mortar::algorithm_mortar)
    dynamic_cast<CONTACT::Element&>(SlaveElement()).PrepareDderiv(MasterElements());

  // decide which type of numerical integration scheme

  //**********************************************************************
  // STANDARD INTEGRATION (SEGMENTS)
  //**********************************************************************
  if (IntType() == Inpar::Mortar::inttype_segments)
  {
    Coupling().resize(0);

    // build linear integration elements from quadratic Mortar::Elements
    std::vector<Teuchos::RCP<Mortar::IntElement>> sauxelements(0);
    std::vector<std::vector<Teuchos::RCP<Mortar::IntElement>>> mauxelements(
        MasterElements().size());
    split_int_elements(SlaveElement(), sauxelements);

    // loop over all master elements associated with this slave element
    for (int m = 0; m < (int)MasterElements().size(); ++m)
    {
      // build linear integration elements from quadratic Mortar::Elements
      mauxelements[m].resize(0);
      split_int_elements(*MasterElements()[m], mauxelements[m]);

      // loop over all IntElement pairs for coupling
      for (int i = 0; i < (int)sauxelements.size(); ++i)
      {
        for (int j = 0; j < (int)mauxelements[m].size(); ++j)
        {
          Coupling().push_back(Teuchos::rcp(new Coupling3dQuad(Discret(), Dim(), true, Params(),
              SlaveElement(), *MasterElements()[m], *sauxelements[i], *mauxelements[m][j])));

          Coupling()[Coupling().size() - 1]->evaluate_coupling();

          // increase counter of slave/master integration pairs and intcells
          smintpairs_ += 1;
          intcells_ += (int)Coupling()[Coupling().size() - 1]->Cells().size();
        }  // for maux
      }    // for saux
    }      // for m

    consist_dual_shape();

    // integrate cells
    for (int i = 0; i < (int)Coupling().size(); ++i)
    {
      if (algo == Inpar::Mortar::algorithm_mortar)
        dynamic_cast<CONTACT::Element&>(SlaveElement())
            .PrepareMderiv(MasterElements(), i % mauxelements.size());

      Coupling()[i]->IntegrateCells(mparams_ptr);

      if (algo == Inpar::Mortar::algorithm_mortar)
        dynamic_cast<CONTACT::Element&>(SlaveElement())
            .assemble_mderiv_to_nodes(Coupling()[i]->MasterElement());
    }
  }

  //**********************************************************************
  // FAST INTEGRATION (ELEMENTS)
  //**********************************************************************
  else if (IntType() == Inpar::Mortar::inttype_elements ||
           IntType() == Inpar::Mortar::inttype_elements_BS)
  {
    // check for standard shape functions and quadratic LM interpolation
    if (ShapeFcn() == Inpar::Mortar::shape_standard &&
        LagMultQuad() == Inpar::Mortar::lagmult_quad &&
        (SlaveElement().Shape() == Core::FE::CellType::quad8 ||
            SlaveElement().Shape() == Core::FE::CellType::tri6))
      FOUR_C_THROW(
          "Quad. LM interpolation for STANDARD 3D quadratic contact only feasible for "
          "quad9");

    if ((int)MasterElements().size() == 0) return;

    // create an integrator instance with correct num_gp and Dim
    Teuchos::RCP<CONTACT::Integrator> integrator =
        CONTACT::INTEGRATOR::BuildIntegrator(stype_, Params(), SlaveElement().Shape(), Comm());

    bool boundary_ele = false;
    bool proj = false;

    // Perform integration and linearization
    integrator->IntegrateDerivEle3D(
        SlaveElement(), MasterElements(), &boundary_ele, &proj, Comm(), mparams_ptr);

    if (IntType() == Inpar::Mortar::inttype_elements_BS)
    {
      if (boundary_ele == true)
      {
        Coupling().resize(0);

        // build linear integration elements from quadratic Mortar::Elements
        std::vector<Teuchos::RCP<Mortar::IntElement>> sauxelements(0);
        std::vector<std::vector<Teuchos::RCP<Mortar::IntElement>>> mauxelements(
            MasterElements().size());
        split_int_elements(SlaveElement(), sauxelements);

        // loop over all master elements associated with this slave element
        for (int m = 0; m < (int)MasterElements().size(); ++m)
        {
          // build linear integration elements from quadratic Mortar::Elements
          mauxelements[m].resize(0);
          split_int_elements(*MasterElements()[m], mauxelements[m]);

          // loop over all IntElement pairs for coupling
          for (int i = 0; i < (int)sauxelements.size(); ++i)
          {
            for (int j = 0; j < (int)mauxelements[m].size(); ++j)
            {
              Coupling().push_back(Teuchos::rcp(new Coupling3dQuad(Discret(), Dim(), true, Params(),
                  SlaveElement(), *MasterElements()[m], *sauxelements[i], *mauxelements[m][j])));

              Coupling()[Coupling().size() - 1]->evaluate_coupling();

              // increase counter of slave/master integration pairs and intcells
              smintpairs_ += 1;
              intcells_ += (int)Coupling()[Coupling().size() - 1]->Cells().size();
            }  // for maux
          }    // for saux
        }      // for m

        consist_dual_shape();

        for (int i = 0; i < (int)Coupling().size(); ++i)
        {
          if (algo == Inpar::Mortar::algorithm_mortar)
            dynamic_cast<CONTACT::Element&>(SlaveElement())
                .PrepareMderiv(MasterElements(), i % mauxelements.size());
          Coupling()[i]->IntegrateCells(mparams_ptr);
          if (algo == Inpar::Mortar::algorithm_mortar)
            dynamic_cast<CONTACT::Element&>(SlaveElement())
                .assemble_mderiv_to_nodes(Coupling()[i]->MasterElement());
        }
      }
    }
  }
  //**********************************************************************
  // INVALID
  //**********************************************************************
  else
  {
    FOUR_C_THROW("Invalid type of numerical integration");
  }

  // free memory of consistent dual shape function coefficient matrix
  SlaveElement().MoData().ResetDualShape();
  SlaveElement().MoData().ResetDerivDualShape();

  if (algo == Inpar::Mortar::algorithm_mortar)
    dynamic_cast<CONTACT::Element&>(SlaveElement())
        .assemble_dderiv_to_nodes((ShapeFcn() == Inpar::Mortar::shape_dual ||
                                   ShapeFcn() == Inpar::Mortar::shape_petrovgalerkin));

  return;
}


/*----------------------------------------------------------------------*
 |  Evaluate coupling pairs for Quad-coupling                farah 01/13|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3dQuadManager::evaluate_coupling(
    const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr)
{
  // decide which type of coupling should be evaluated
  Inpar::Mortar::AlgorithmType algo =
      Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(Params(), "ALGORITHM");

  //*********************************
  // Mortar Contact
  //*********************************
  if (algo == Inpar::Mortar::algorithm_mortar || algo == Inpar::Mortar::algorithm_gpts)
    integrate_coupling(mparams_ptr);

  //*********************************
  // Error
  //*********************************
  else
    FOUR_C_THROW("chosen contact algorithm not supported!");

  return true;
}


/*----------------------------------------------------------------------*
 |  Calculate dual shape functions                           seitz 07/13|
 *----------------------------------------------------------------------*/
void CONTACT::Coupling3dManager::consist_dual_shape()
{
  static const Inpar::Mortar::AlgorithmType algo =
      Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(imortar_, "ALGORITHM");
  if (algo != Inpar::Mortar::algorithm_mortar) return;

  // For standard shape functions no modification is necessary
  // A switch erlier in the process improves computational efficiency
  Inpar::Mortar::ConsistentDualType consistent =
      Core::UTILS::IntegralValue<Inpar::Mortar::ConsistentDualType>(imortar_, "LM_DUAL_CONSISTENT");
  if (ShapeFcn() == Inpar::Mortar::shape_standard || consistent == Inpar::Mortar::consistent_none)
    return;

  // Consistent modification not yet checked for constant LM interpolation
  if (Quad() == true && LagMultQuad() == Inpar::Mortar::lagmult_const &&
      consistent != Inpar::Mortar::consistent_none)
    FOUR_C_THROW("Consistent dual shape functions not yet checked for constant LM interpolation!");

  if (consistent == Inpar::Mortar::consistent_all && IntType() != Inpar::Mortar::inttype_segments)
    FOUR_C_THROW(
        "Consistent dual shape functions on all elements only for segment-based "
        "integration");

  // do nothing if there are no coupling pairs
  if (Coupling().size() == 0) return;

  // check for boundary elements in segment-based integration
  // (fast integration already has this check, so that consist_dual_shape()
  // is only called for boundary elements)
  //
  // For NURBS elements, always compute consistent dual functions.
  // This improves robustness, since the duality is enforced at exactly
  // the same quadrature points, that the mortar integrals etc. are evaluated.
  // For Lagrange FE, the calculation of dual shape functions for fully
  // projecting elements is ok, since the integrands are polynomials (except
  // the jacobian)
  if (IntType() == Inpar::Mortar::inttype_segments &&
      consistent == Inpar::Mortar::consistent_boundary)
  {
    // check, if slave element is fully projecting
    // for convenience, we don't check each quadrature point
    // but only the element nodes. This usually does the job.
    bool boundary_ele = false;

    Core::FE::CellType dt_s = SlaveElement().Shape();

    double sxi_test[2] = {0.0, 0.0};
    double alpha_test = 0.0;
    bool proj_test = false;

    Core::Nodes::Node** mynodes_test = SlaveElement().Nodes();
    if (!mynodes_test) FOUR_C_THROW("has_proj_status: Null pointer!");

    if (dt_s == Core::FE::CellType::quad4 || dt_s == Core::FE::CellType::quad8 ||
        dt_s == Core::FE::CellType::nurbs9)
    {
      for (int s_test = 0; s_test < SlaveElement().num_node(); ++s_test)
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
          Mortar::Projector::Impl(SlaveElement(), Coupling()[bs_test]->MasterIntElement())
              ->ProjectGaussPoint3D(SlaveElement(), sxi_test,
                  Coupling()[bs_test]->MasterIntElement(), mxi_test, alpha_test);

          Core::FE::CellType dt = Coupling()[bs_test]->MasterIntElement().Shape();
          if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
              dt == Core::FE::CellType::quad9)
          {
            if (mxi_test[0] >= -1.0 && mxi_test[1] >= -1.0 && mxi_test[0] <= 1.0 &&
                mxi_test[1] <= 1.0)
              proj_test = true;
          }
          else if (dt == Core::FE::CellType::tri3 || dt == Core::FE::CellType::tri6)
          {
            if (mxi_test[0] >= 0.0 && mxi_test[1] >= 0.0 && mxi_test[0] <= 1.0 &&
                mxi_test[1] <= 1.0 && mxi_test[0] + mxi_test[1] <= 1.0)
              proj_test = true;
          }
          else
          {
            FOUR_C_THROW("Non valid element type for master discretization!");
          }
        }
        if (proj_test == false) boundary_ele = true;
      }
    }

    else if (dt_s == Core::FE::CellType::tri3 || dt_s == Core::FE::CellType::tri6)
    {
      for (int s_test = 0; s_test < SlaveElement().num_node(); ++s_test)
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
          Mortar::Projector::Impl(SlaveElement(), Coupling()[bs_test]->MasterElement())
              ->ProjectGaussPoint3D(SlaveElement(), sxi_test, Coupling()[bs_test]->MasterElement(),
                  mxi_test, alpha_test);

          Core::FE::CellType dt = Coupling()[bs_test]->MasterElement().Shape();
          if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
              dt == Core::FE::CellType::quad9)
          {
            if (mxi_test[0] >= -1.0 && mxi_test[1] >= -1.0 && mxi_test[0] <= 1.0 &&
                mxi_test[1] <= 1.0)
              proj_test = true;
          }
          else if (dt == Core::FE::CellType::tri3 || dt == Core::FE::CellType::tri6)
          {
            if (mxi_test[0] >= 0.0 && mxi_test[1] >= 0.0 && mxi_test[0] <= 1.0 &&
                mxi_test[1] <= 1.0 && mxi_test[0] + mxi_test[1] <= 1.0)
              proj_test = true;
          }
          else
          {
            FOUR_C_THROW("Non valid element type for master discretization!");
          }
        }
        if (proj_test == false) boundary_ele = true;
      }
    }

    else
      FOUR_C_THROW(
          "Calculation of consistent dual shape functions called for non-valid slave element "
          "shape!");

    if (boundary_ele == false) return;
  }

  // slave nodes and dofs
  const int max_nnodes = 9;
  const int nnodes = SlaveElement().num_node();
  if (nnodes > max_nnodes)
    FOUR_C_THROW(
        "this function is not implemented to handle elements with that many nodes. Just adjust "
        "max_nnodes above");
  const int ndof = 3;
  const int msize = MasterElements().size();

  // get number of master nodes
  int mnodes = 0;
  for (int m = 0; m < msize; ++m) mnodes += MasterElements()[m]->num_node();

  // Dual shape functions coefficient matrix and linearization
  SlaveElement().MoData().DerivDualShape() =
      Teuchos::rcp(new Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>(
          (nnodes + mnodes) * ndof, 0, Core::LinAlg::SerialDenseMatrix(nnodes, nnodes)));
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& derivae =
      *(SlaveElement().MoData().DerivDualShape());

  // various variables
  double detg = 0.0;
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // initialize matrices de and me
  Core::LinAlg::SerialDenseMatrix me(nnodes, nnodes, true);
  Core::LinAlg::SerialDenseMatrix de(nnodes, nnodes, true);

  Core::Gen::Pairedvector<int, Core::LinAlg::Matrix<max_nnodes + 1, max_nnodes>> derivde_new(
      (nnodes + mnodes) * ndof);

  // two-dim arrays of maps for linearization of me/de
  std::vector<std::vector<Core::Gen::Pairedvector<int, double>>> derivme(
      nnodes, std::vector<Core::Gen::Pairedvector<int, double>>(nnodes, (nnodes + mnodes) * ndof));
  std::vector<std::vector<Core::Gen::Pairedvector<int, double>>> derivde(
      nnodes, std::vector<Core::Gen::Pairedvector<int, double>>(nnodes, (nnodes + mnodes) * ndof));

  double A_tot = 0.;
  // loop over all master elements associated with this slave element
  for (int m = 0; m < (int)Coupling().size(); ++m)
  {
    if (!Coupling()[m]->RoughCheckCenters()) continue;
    if (!Coupling()[m]->rough_check_orient()) continue;
    if (!Coupling()[m]->RoughCheckCenters()) continue;

    // get number of master nodes
    const int ncol = Coupling()[m]->MasterElement().num_node();

    // loop over all integration cells
    for (int c = 0; c < (int)Coupling()[m]->Cells().size(); ++c)
    {
      Teuchos::RCP<Mortar::IntCell> currcell = Coupling()[m]->Cells()[c];

      A_tot += currcell->Area();

      // create an integrator for this cell
      CONTACT::Integrator integrator(imortar_, currcell->Shape(), Comm());

      // check if the cells are tri3
      // there's nothing wrong about other shapes, but as long as they are all
      // tri3 we can perform the jacobian calculation ( and its deriv) outside
      // the Gauss point loop
      if (currcell->Shape() != Core::FE::CellType::tri3)
        FOUR_C_THROW("only tri3 integration cells at the moment. See comment in the code");

      detg = currcell->Jacobian();
      // directional derivative of cell Jacobian
      Core::Gen::Pairedvector<int, double> derivjaccell((nnodes + ncol) * ndof);
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
        Mortar::Projector::Impl(Coupling()[m]->SlaveIntElement())
            ->project_gauss_point_auxn3_d(
                globgp, Coupling()[m]->Auxn(), Coupling()[m]->SlaveIntElement(), sxi, sprojalpha);

        // project Gauss point onto slave (parent) element
        double psxi[2] = {0., 0.};
        double psprojalpha = 0.0;
        if (Quad())
        {
          Mortar::IntElement* ie =
              dynamic_cast<Mortar::IntElement*>(&(Coupling()[m]->SlaveIntElement()));
          if (ie == nullptr) FOUR_C_THROW("nullptr pointer");
          Mortar::Projector::Impl(SlaveElement())
              ->project_gauss_point_auxn3_d(
                  globgp, Coupling()[m]->Auxn(), SlaveElement(), psxi, psprojalpha);
          // ie->MapToParent(sxi,psxi); // old way of doing it via affine map... wrong (popp
          // 05/2016)
        }
        else
          for (int i = 0; i < 2; ++i) psxi[i] = sxi[i];

        // create vector for shape function evaluation
        Core::LinAlg::SerialDenseVector sval(nnodes);
        Core::LinAlg::SerialDenseMatrix sderiv(nnodes, 2, true);

        // evaluate trace space shape functions at Gauss point
        if (LagMultQuad() == Inpar::Mortar::lagmult_lin)
          SlaveElement().evaluate_shape_lag_mult_lin(
              Inpar::Mortar::shape_standard, psxi, sval, sderiv, nnodes);
        else
          SlaveElement().evaluate_shape(psxi, sval, sderiv, nnodes);

        // additional data for contact calculation (i.e. incl. derivative of dual shape functions
        // coefficient matrix) GP slave coordinate derivatives
        std::vector<Core::Gen::Pairedvector<int, double>> dsxigp(2, (nnodes + ncol) * ndof);
        // GP slave coordinate derivatives
        std::vector<Core::Gen::Pairedvector<int, double>> dpsxigp(2, (nnodes + ncol) * ndof);
        // global GP coordinate derivative on integration element
        Core::Gen::Pairedvector<int, Core::LinAlg::Matrix<3, 1>> lingp((nnodes + ncol) * ndof);

        // compute global GP coordinate derivative
        static Core::LinAlg::Matrix<3, 1> svalcell;
        static Core::LinAlg::Matrix<3, 2> sderivcell;
        currcell->evaluate_shape(eta, svalcell, sderivcell);

        for (int v = 0; v < 3; ++v)
          for (int d = 0; d < 3; ++d)
            for (_CI p = (currcell->GetDerivVertex(v))[d].begin();
                 p != (currcell->GetDerivVertex(v))[d].end(); ++p)
              lingp[p->first](d) += svalcell(v) * (p->second);

        // compute GP slave coordinate derivatives
        integrator.DerivXiGP3DAuxPlane(Coupling()[m]->SlaveIntElement(), sxi, currcell->Auxn(),
            dsxigp, sprojalpha, currcell->get_deriv_auxn(), lingp);

        // compute GP slave coordinate derivatives (parent element)
        if (Quad())
        {
          Mortar::IntElement* ie =
              dynamic_cast<Mortar::IntElement*>(&(Coupling()[m]->SlaveIntElement()));
          if (ie == nullptr) FOUR_C_THROW("wtf");
          integrator.DerivXiGP3DAuxPlane(SlaveElement(), psxi, currcell->Auxn(), dpsxigp,
              psprojalpha, currcell->get_deriv_auxn(), lingp);
          // ie->MapToParent(dsxigp,dpsxigp); // old way of doing it via affine map... wrong (popp
          // 05/2016)
        }
        else
          dpsxigp = dsxigp;

        double fac = 0.;
        for (_CI p = derivjaccell.begin(); p != derivjaccell.end(); ++p)
        {
          Core::LinAlg::Matrix<max_nnodes + 1, max_nnodes>& dtmp = derivde_new[p->first];
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
            Core::LinAlg::Matrix<max_nnodes + 1, max_nnodes>& dtmp = derivde_new[p->first];
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
  Core::LinAlg::SerialDenseMatrix ae(nnodes, nnodes, true);
  Core::LinAlg::SerialDenseMatrix meinv(nnodes, nnodes, true);

  // compute matrix A_e and inverse of matrix M_e for
  // linear interpolation of quadratic element
  if (LagMultQuad() == Inpar::Mortar::lagmult_lin)
  {
    // declare and initialize to zero inverse of Matrix M_e
    Core::LinAlg::SerialDenseMatrix meinv(nnodes, nnodes, true);

    if (SlaveElement().Shape() == Core::FE::CellType::tri6)
    {
      // reduce me to non-zero nodes before inverting
      Core::LinAlg::Matrix<3, 3> melin;
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k) melin(j, k) = me(j, k);

      // invert bi-ortho matrix melin
      Core::LinAlg::Inverse(melin);

      // re-inflate inverse of melin to full size
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k) meinv(j, k) = melin(j, k);
    }
    else if (SlaveElement().Shape() == Core::FE::CellType::quad8 ||
             SlaveElement().Shape() == Core::FE::CellType::quad9)
    {
      // reduce me to non-zero nodes before inverting
      Core::LinAlg::Matrix<4, 4> melin;
      for (int j = 0; j < 4; ++j)
        for (int k = 0; k < 4; ++k) melin(j, k) = me(j, k);

      // invert bi-ortho matrix melin
      Core::LinAlg::Inverse(melin);

      // re-inflate inverse of melin to full size
      for (int j = 0; j < 4; ++j)
        for (int k = 0; k < 4; ++k) meinv(j, k) = melin(j, k);
    }
    else
      FOUR_C_THROW("incorrect element shape for linear interpolation of quadratic element!");

    // get solution matrix with dual parameters
    Core::LinAlg::multiply(ae, de, meinv);
  }
  // compute matrix A_e and inverse of matrix M_e for all other cases
  else
    meinv = Core::LinAlg::InvertAndMultiplyByCholesky(me, de, ae);

  // build linearization of ae and store in derivdual
  // (this is done according to a quite complex formula, which
  // we get from the linearization of the biorthogonality condition:
  // Lin (Me * Ae = De) -> Lin(Ae)=Lin(De)*Inv(Me)-Ae*Lin(Me)*Inv(Me) )
  typedef Core::Gen::Pairedvector<int,
      Core::LinAlg::Matrix<max_nnodes + 1, max_nnodes>>::const_iterator _CIM;
  for (_CIM p = derivde_new.begin(); p != derivde_new.end(); ++p)
  {
    Core::LinAlg::Matrix<max_nnodes + 1, max_nnodes>& dtmp = derivde_new[p->first];
    Core::LinAlg::SerialDenseMatrix& pt = derivae[p->first];
    for (int i = 0; i < nnodes; ++i)
      for (int j = 0; j < nnodes; ++j)
      {
        pt(i, j) += meinv(i, j) * dtmp(nnodes, i);

        for (int k = 0; k < nnodes; ++k)
          for (int l = 0; l < nnodes; ++l) pt(i, j) -= ae(i, k) * meinv(l, j) * dtmp(l, k);
      }
  }

  // store ae matrix in slave element data container
  SlaveElement().MoData().DualShape() = Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(ae));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Coupling3dManager::find_feasible_master_elements(
    std::vector<Mortar::Element*>& feasible_ma_eles) const
{
  // feasibility counter
  std::size_t fcount = 0;
  for (std::size_t m = 0; m < MasterElements().size(); ++m)
  {
    // Build a instance of the Mortar::Coupling3d object (no linearization needed).
    Teuchos::RCP<Mortar::Coupling3d> coup = Teuchos::rcp(
        new Mortar::Coupling3d(idiscret_, dim_, false, imortar_, SlaveElement(), MasterElement(m)));

    // Building the master element normals and check the angles.
    if (coup->rough_check_orient())
    {
      feasible_ma_eles[fcount] = &MasterElement(m);
      ++fcount;
    }
  }
  feasible_ma_eles.resize(fcount);

  return;
}

FOUR_C_NAMESPACE_CLOSE
