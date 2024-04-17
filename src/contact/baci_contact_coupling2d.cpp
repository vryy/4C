/*---------------------------------------------------------------------*/
/*! \file
\brief Classes for mortar contact coupling in 2D.

\level 2


*/
/*---------------------------------------------------------------------*/

#include "baci_contact_coupling2d.hpp"

#include "baci_contact_defines.hpp"
#include "baci_contact_element.hpp"
#include "baci_contact_integrator.hpp"
#include "baci_contact_integrator_factory.hpp"
#include "baci_contact_interpolator.hpp"
#include "baci_contact_node.hpp"
#include "baci_lib_discret.hpp"
#include "baci_linalg_serialdensevector.hpp"
#include "baci_linalg_utils_densematrix_inverse.hpp"
#include "baci_linalg_utils_densematrix_multiply.hpp"
#include "baci_mortar_defines.hpp"
#include "baci_mortar_element.hpp"
#include "baci_mortar_node.hpp"
#include "baci_mortar_projector.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 06/09|
 *----------------------------------------------------------------------*/
CONTACT::Coupling2d::Coupling2d(DRT::Discretization& idiscret, int dim, bool quad,
    Teuchos::ParameterList& params, MORTAR::Element& sele, MORTAR::Element& mele)
    : MORTAR::Coupling2d(idiscret, dim, quad, params, sele, mele),
      stype_(CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(params, "STRATEGY"))
{
  // empty constructor

  return;
}


/*----------------------------------------------------------------------*
 |  Integrate slave / master overlap (public)                 popp 04/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling2d::IntegrateOverlap(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  // explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("IntegrateOverlap called without specific shape function defined!");

  /**********************************************************************/
  /* INTEGRATION                                                        */
  /* Depending on overlap and the xiproj_ entries integrate the Mortar  */
  /* matrices D and M and the weighted gap function g~ on the overlap   */
  /* of the current sl / ma pair.                                       */
  /**********************************************************************/

  // no integration if no overlap
  if (!overlap_) return false;

  // set segmentation status of all slave nodes
  // (hassegment_ of a slave node is true if ANY segment/cell
  // is integrated that contributes to this slave node)
  int nnodes = SlaveElement().NumNode();
  DRT::Node** mynodes = SlaveElement().Nodes();
  if (!mynodes) dserror("Null pointer!");
  for (int k = 0; k < nnodes; ++k)
  {
    MORTAR::Node* mycnode = dynamic_cast<MORTAR::Node*>(mynodes[k]);
    if (!mycnode) dserror("Null pointer!");
    mycnode->HasSegment() = true;
  }

  // local working copies of input variables
  double sxia = xiproj_[0];
  double sxib = xiproj_[1];
  double mxia = xiproj_[2];
  double mxib = xiproj_[3];

  // create a CONTACT integrator instance with correct NumGP and Dim
  Teuchos::RCP<CONTACT::Integrator> integrator =
      CONTACT::INTEGRATOR::BuildIntegrator(stype_, imortar_, SlaveElement().Shape(), Comm());
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
  // cases (1), (2) and (3)
  // *******************************************************************
  if (!Quad() || (Quad() && lmtype == INPAR::MORTAR::lagmult_quad) ||
      (Quad() && lmtype == INPAR::MORTAR::lagmult_lin) ||
      (Quad() && lmtype == INPAR::MORTAR::lagmult_const))
  {
    // ***********************************************************
    //                   Integrate stuff !!!                    //
    // ***********************************************************
    integrator->IntegrateDerivSegment2D(
        SlaveElement(), sxia, sxib, MasterElement(), mxia, mxib, Comm(), mparams_ptr);
    // ***********************************************************
    //                   END INTEGRATION !!!                    //
    // ***********************************************************
  }

  // *******************************************************************
  // case (4)
  // *******************************************************************
  else if (Quad() && lmtype == INPAR::MORTAR::lagmult_pwlin)
  {
    dserror("Piecewise linear LM not (yet?) implemented in 2D");
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
  {
    dserror("IntegrateOverlap: Invalid case for 2D mortar coupling LM interpolation");
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 06/09|
 *----------------------------------------------------------------------*/
CONTACT::Coupling2dManager::Coupling2dManager(DRT::Discretization& idiscret, int dim, bool quad,
    Teuchos::ParameterList& params, MORTAR::Element* sele, std::vector<MORTAR::Element*> mele)
    : MORTAR::Coupling2dManager(idiscret, dim, quad, params, sele, mele),
      stype_(CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(params, "STRATEGY"))
{
  // empty constructor
  return;
}


/*----------------------------------------------------------------------*
 |  get communicator  (public)                               farah 01/13|
 *----------------------------------------------------------------------*/
const Epetra_Comm& CONTACT::Coupling2dManager::Comm() const { return idiscret_.Comm(); }

/*----------------------------------------------------------------------*
 |  Evaluate coupling pairs                                  farah 10/14|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling2dManager::EvaluateCoupling(
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  if (MasterElements().size() == 0) return false;

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

  return true;
}


/*----------------------------------------------------------------------*
 |  Evaluate mortar coupling pairs                           Popp 03/09 |
 *----------------------------------------------------------------------*/
void CONTACT::Coupling2dManager::IntegrateCoupling(
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  //**********************************************************************
  // STANDARD INTEGRATION (SEGMENTS)
  //**********************************************************************
  if (IntType() == INPAR::MORTAR::inttype_segments)
  {
    // loop over all master elements associated with this slave element
    for (int m = 0; m < (int)MasterElements().size(); ++m)
    {
      // create Coupling2d object and push back
      Coupling().push_back(Teuchos::rcp(
          new Coupling2d(idiscret_, dim_, quad_, imortar_, SlaveElement(), MasterElement(m))));

      // project the element pair
      Coupling()[m]->Project();

      // check for element overlap
      Coupling()[m]->DetectOverlap();
    }

    // calculate consistent dual shape functions for this element
    ConsistDualShape();

    // do mortar integration
    for (int m = 0; m < (int)MasterElements().size(); ++m)
      Coupling()[m]->IntegrateOverlap(mparams_ptr);

    // free memory of consistent dual shape function coefficient matrix
    SlaveElement().MoData().ResetDualShape();
    SlaveElement().MoData().ResetDerivDualShape();
  }
  //**********************************************************************
  // FAST INTEGRATION (ELEMENTS)
  //**********************************************************************
  else if (IntType() == INPAR::MORTAR::inttype_elements or
           IntType() == INPAR::MORTAR::inttype_elements_BS)
  {
    if ((int)MasterElements().size() == 0) return;

    // create an integrator instance with correct NumGP and Dim
    Teuchos::RCP<CONTACT::Integrator> integrator =
        CONTACT::INTEGRATOR::BuildIntegrator(stype_, imortar_, SlaveElement().Shape(), Comm());

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
    // cases (1), (2) and (3)
    // *******************************************************************
    if (!Quad() || (Quad() && lmtype == INPAR::MORTAR::lagmult_quad) ||
        (Quad() && lmtype == INPAR::MORTAR::lagmult_lin))
    {
      // Test whether projection from slave to master surface is feasible -->
      // important for dual LM Fnc.
      // Contact_interface.cpp --> AssembleG
      for (unsigned m = 0; m < MasterElements().size(); ++m)
      {
        // create Coupling2d object and push back
        Coupling().push_back(Teuchos::rcp(
            new Coupling2d(idiscret_, dim_, quad_, imortar_, SlaveElement(), MasterElement(m))));

        // project the element pair
        Coupling()[m]->Project();
      }

      // Bool for identification of boundary elements
      bool boundary_ele = false;

      // ***********************************************************
      //                  START INTEGRATION !!!                   //
      // ***********************************************************
      integrator->IntegrateDerivEle2D(SlaveElement(), MasterElements(), &boundary_ele, mparams_ptr);
      // ***********************************************************
      //                   END INTEGRATION !!!                    //
      // ***********************************************************

      if (IntType() == INPAR::MORTAR::inttype_elements_BS and boundary_ele == true)
      {
        // switch, if consistent boundary modification chosen
        if (CORE::UTILS::IntegralValue<int>(imortar_, "LM_DUAL_CONSISTENT") == true &&
            ShapeFcn() != INPAR::MORTAR::shape_standard  // so for petrov-Galerkin and dual
        )
        {
          // loop over all master elements associated with this slave element
          for (int m = 0; m < (int)MasterElements().size(); ++m)
          {
            // create Coupling2d object and push back
            Coupling().push_back(Teuchos::rcp(new Coupling2d(
                idiscret_, dim_, quad_, imortar_, SlaveElement(), MasterElement(m))));

            // project the element pair
            Coupling()[m]->Project();

            // check for element overlap
            Coupling()[m]->DetectOverlap();
          }

          // calculate consistent dual shape functions for this element
          ConsistDualShape();

          // do mortar integration
          for (int m = 0; m < (int)MasterElements().size(); ++m)
            Coupling()[m]->IntegrateOverlap(mparams_ptr);

          // free memory of consistent dual shape function coefficient matrix
          SlaveElement().MoData().ResetDualShape();
          SlaveElement().MoData().ResetDerivDualShape();
        }

        // segment-based integration for boundary elements
        else
        {
          // loop over all master elements associated with this slave element
          for (int m = 0; m < (int)MasterElements().size(); ++m)
          {
            // create Coupling2d object and push back
            Coupling().push_back(Teuchos::rcp(new Coupling2d(
                idiscret_, dim_, quad_, imortar_, SlaveElement(), MasterElement(m))));

            // project the element pair
            Coupling()[m]->Project();

            // check for element overlap
            Coupling()[m]->DetectOverlap();

            // integrate the element overlap
            Coupling()[m]->IntegrateOverlap(mparams_ptr);
          }
        }
      }
      else
      {
        // nothing to do...
      }
    }
    // *******************************************************************
    // case (4)
    // *******************************************************************
    else if (Quad() && lmtype == INPAR::MORTAR::lagmult_pwlin)
    {
      dserror("Piecewise linear LM not (yet?) implemented in 2D");
    }

    // *******************************************************************
    // undefined case
    // *******************************************************************
    else if (Quad() && lmtype == INPAR::MORTAR::lagmult_undefined)
    {
      dserror(
          "Lagrange multiplier interpolation for quadratic elements undefined\n"
          "If you are using 2nd order mortar elements, you need to specify LM_QUAD"
          " in MORTAR COUPLING section");
    }

    // *******************************************************************
    // other cases
    // *******************************************************************
    else
    {
      dserror("Integrate: Invalid case for 2D mortar coupling LM interpolation");
    }
  }
  //**********************************************************************
  // INVALID TYPE OF NUMERICAL INTEGRATION
  //**********************************************************************
  else
  {
    dserror("Invalid type of numerical integration");
  }
}


/*----------------------------------------------------------------------*
 |  Calculate dual shape functions                           seitz 07/13|
 *----------------------------------------------------------------------*/
void CONTACT::Coupling2dManager::ConsistDualShape()
{
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

  // do nothing if there are no coupling pairs
  if (Coupling().size() == 0) return;

  const int nnodes = SlaveElement().NumNode();
  const int ndof = 2;

  int linsize = 0;
  for (int i = 0; i < nnodes; ++i)
  {
    Node* cnode = dynamic_cast<Node*>(SlaveElement().Nodes()[i]);
    linsize += cnode->GetLinsize();
  }

  int mnodes = 0;
  for (int m = 0; m < (int)Coupling().size(); ++m) mnodes += MasterElements()[m]->NumNode();

  // detect entire overlap
  double ximin = 1.0;
  double ximax = -1.0;
  CORE::GEN::pairedvector<int, double> dximin(linsize + ndof * mnodes);
  CORE::GEN::pairedvector<int, double> dximax(linsize + ndof * mnodes);

  // loop over all master elements associated with this slave element
  for (int m = 0; m < (int)Coupling().size(); ++m)
  {
    double sxia = Coupling()[m]->XiProj()[0];
    double sxib = Coupling()[m]->XiProj()[1];
    double mxia = Coupling()[m]->XiProj()[2];
    double mxib = Coupling()[m]->XiProj()[3];

    // no overlap for this slave-master pair --> continue with next pair
    if (sxia == 0.0 && sxib == 0.0) continue;

    // for contact we need the derivatives as well
    bool startslave = false;
    bool endslave = false;

    if (sxia == -1.0)
      startslave = true;
    else
      startslave = false;
    if (sxib == 1.0)
      endslave = true;
    else
      endslave = false;

    // create an integrator for this segment
    CONTACT::Integrator integrator(imortar_, SlaveElement().Shape(), Comm());

    std::vector<CORE::GEN::pairedvector<int, double>> ximaps(4, linsize + ndof * mnodes);
    // get directional derivatives of sxia, sxib, mxia, mxib
    integrator.DerivXiAB2D(SlaveElement(), sxia, sxib, MasterElement(m), mxia, mxib, ximaps,
        startslave, endslave, linsize);

    // get element contact integration area
    // and for contact derivatives of beginning and end
    if ((sxia != 0.0 || sxib != 0.0) && (sxia >= -1.0 && sxia <= 1.0) &&
        (sxib >= -1.0 && sxib <= 1.0))
    {
      if (sxia < ximin && sxia >= -1. && sxia <= 1.)
      {
        ximin = sxia;
        dximin = ximaps[0];
      }
      if (sxib > ximax && sxib >= -1. && sxib <= 1.)
      {
        ximax = sxib;
        dximax = ximaps[1];
      }
    }
  }

  // map iterator
  typedef CORE::GEN::pairedvector<int, double>::const_iterator _CI;

  // no overlap: the applied dual shape functions don't matter, as the integration domain is void
  if ((ximax == -1.0 && ximin == 1.0) || (ximax - ximin < 4. * MORTARINTLIM)) return;

  // fully projecting element: no modification necessary
  if (ximin == -1.0 && ximax == 1.0) return;

  // calculate consistent dual schape functions (see e.g. Cichosz et.al.:
  // Consistent treatment of boundaries with mortar contact formulations, CMAME 2010

  // store derivae into element
  SlaveElement().MoData().DerivDualShape() =
      Teuchos::rcp(new CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>(
          linsize + 2 * ndof * mnodes, 0, CORE::LINALG::SerialDenseMatrix(nnodes, nnodes)));
  CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& derivae =
      *(SlaveElement().MoData().DerivDualShape());

  // compute entries to bi-ortho matrices me/de with Gauss quadrature
  MORTAR::ElementIntegrator integrator(SlaveElement().Shape());

  // prepare for calculation of dual shape functions
  CORE::LINALG::SerialDenseMatrix me(nnodes, nnodes, true);
  CORE::LINALG::SerialDenseMatrix de(nnodes, nnodes, true);
  // two-dim arrays of maps for linearization of me/de
  std::vector<std::vector<CORE::GEN::pairedvector<int, double>>> derivme(nnodes,
      std::vector<CORE::GEN::pairedvector<int, double>>(nnodes, linsize + 2 * ndof * mnodes));
  std::vector<std::vector<CORE::GEN::pairedvector<int, double>>> derivde(nnodes,
      std::vector<CORE::GEN::pairedvector<int, double>>(nnodes, linsize + 2 * ndof * mnodes));

  CORE::LINALG::SerialDenseVector sval(nnodes);
  CORE::LINALG::SerialDenseMatrix sderiv(nnodes, 1, true);
  CORE::LINALG::SerialDenseMatrix ssecderiv(nnodes, 1);

  for (int gp = 0; gp < integrator.nGP(); ++gp)
  {
    // coordinates and weight
    std::array<double, 2> eta = {integrator.Coordinate(gp, 0), 0.0};
    double wgt = integrator.Weight(gp);

    // coordinate transformation sxi->eta (slave MORTAR::Element->Overlap)
    double sxi[2] = {0.0, 0.0};
    sxi[0] = 0.5 * (1.0 - eta[0]) * ximin + 0.5 * (1.0 + eta[0]) * ximax;

    // evaluate trace space shape functions
    if (LagMultQuad() == INPAR::MORTAR::lagmult_lin)
      SlaveElement().EvaluateShapeLagMultLin(
          INPAR::MORTAR::shape_standard, sxi, sval, sderiv, nnodes);
    else
      SlaveElement().EvaluateShape(sxi, sval, sderiv, nnodes);
    SlaveElement().Evaluate2ndDerivShape(sxi, ssecderiv, nnodes);

    // evaluate the two slave side Jacobians
    double dxdsxi = SlaveElement().Jacobian(sxi);
    double dsxideta = -0.5 * ximin + 0.5 * ximax;

    // evaluate linearizations *******************************************
    // evaluate the derivative dxdsxidsxi = Jac,xi
    double djacdxi[2] = {0.0, 0.0};
    dynamic_cast<CONTACT::Element&>(SlaveElement()).DJacDXi(djacdxi, sxi, ssecderiv);
    double dxdsxidsxi = djacdxi[0];  // only 2D here

    // evalute the GP slave coordinate derivatives
    CORE::GEN::pairedvector<int, double> dsxigp(linsize + ndof * mnodes);
    for (_CI p = dximin.begin(); p != dximin.end(); ++p)
      dsxigp[p->first] += 0.5 * (1 - eta[0]) * (p->second);
    for (_CI p = dximax.begin(); p != dximax.end(); ++p)
      dsxigp[p->first] += 0.5 * (1 + eta[0]) * (p->second);

    // evaluate the Jacobian derivative
    CORE::GEN::pairedvector<int, double> derivjac(SlaveElement().NumNode() * Dim());
    SlaveElement().DerivJacobian(sxi, derivjac);

    // integrate dual shape matrices de, me and their linearizations
    for (int j = 0; j < nnodes; ++j)
    {
      double fac;
      // de and linearization
      de(j, j) += wgt * sval[j] * dxdsxi * dsxideta;

      // (1) linearization of slave gp coordinates in ansatz function j for derivative of de
      fac = wgt * sderiv(j, 0) * dxdsxi * dsxideta;
      for (_CI p = dsxigp.begin(); p != dsxigp.end(); ++p)
        derivde[j][j][p->first] += fac * (p->second);

      // (2) linearization dsxideta - segment end coordinates
      fac = 0.5 * wgt * sval[j] * dxdsxi;
      for (_CI p = dximin.begin(); p != dximin.end(); ++p)
        derivde[j][j][p->first] -= fac * (p->second);
      fac = 0.5 * wgt * sval[j] * dxdsxi;
      for (_CI p = dximax.begin(); p != dximax.end(); ++p)
        derivde[j][j][p->first] += fac * (p->second);

      // (3) linearization dxdsxi - slave GP jacobian
      fac = wgt * sval[j] * dsxideta;
      for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
        derivde[j][j][p->first] += fac * (p->second);

      // (4) linearization dxdsxi - slave GP coordinates
      fac = wgt * sval[j] * dsxideta * dxdsxidsxi;
      for (_CI p = dsxigp.begin(); p != dsxigp.end(); ++p)
        derivde[j][j][p->first] += fac * (p->second);

      // me and linearization
      for (int k = 0; k < nnodes; ++k)
      {
        me(j, k) += wgt * sval[j] * sval[k] * dxdsxi * dsxideta;

        // (1) linearization of slave gp coordinates in ansatz function for derivative of me
        fac = wgt * sval[k] * dxdsxi * dsxideta * sderiv(j, 0);
        for (_CI p = dsxigp.begin(); p != dsxigp.end(); ++p)
        {
          derivme[j][k][p->first] += fac * (p->second);
          derivme[k][j][p->first] += fac * (p->second);
        }

        // (2) linearization dsxideta - segment end coordinates
        fac = 0.5 * wgt * sval[j] * sval[k] * dxdsxi;
        for (_CI p = dximin.begin(); p != dximin.end(); ++p)
          derivme[j][k][p->first] -= fac * (p->second);
        fac = 0.5 * wgt * sval[j] * sval[k] * dxdsxi;
        for (_CI p = dximax.begin(); p != dximax.end(); ++p)
          derivme[j][k][p->first] += fac * (p->second);

        // (3) linearization dxdsxi - slave GP jacobian
        fac = wgt * sval[j] * sval[k] * dsxideta;
        for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
          derivme[j][k][p->first] += fac * (p->second);

        // (4) linearization dxdsxi - slave GP coordinates
        fac = wgt * sval[j] * sval[k] * dsxideta * dxdsxidsxi;
        for (_CI p = dsxigp.begin(); p != dsxigp.end(); ++p)
          derivme[j][k][p->first] += fac * (p->second);
      }
    }
  }

  // declare dual shape functions coefficient matrix and
  // inverse of matrix M_e
  CORE::LINALG::SerialDenseMatrix ae(nnodes, nnodes, true);
  CORE::LINALG::SerialDenseMatrix meinv(nnodes, nnodes, true);

  // compute matrix A_e and inverse of matrix M_e for
  // linear interpolation of quadratic element
  if (LagMultQuad() == INPAR::MORTAR::lagmult_lin)
  {
    // how many non-zero nodes
    const int nnodeslin = 2;

    // reduce me to non-zero nodes before inverting
    CORE::LINALG::Matrix<nnodeslin, nnodeslin> melin;
    for (int j = 0; j < nnodeslin; ++j)
      for (int k = 0; k < nnodeslin; ++k) melin(j, k) = me(j, k);

    // invert bi-ortho matrix melin
    CORE::LINALG::Inverse(melin);

    // re-inflate inverse of melin to full size
    for (int j = 0; j < nnodeslin; ++j)
      for (int k = 0; k < nnodeslin; ++k) meinv(j, k) = melin(j, k);

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

  // loop over all entries of ae (index i,j)
  for (int i = 0; i < nnodes; ++i)
  {
    for (int j = 0; j < nnodes; ++j)
    {
      // compute Lin(Ae) according to formula above
      for (int l = 0; l < nnodes; ++l)  // loop over sum l
      {
        // part1: Lin(De)*Inv(Me)
        for (_CI p = derivde[i][l].begin(); p != derivde[i][l].end(); ++p)
          derivae[p->first](i, j) += meinv(l, j) * (p->second);

        // part2: Ae*Lin(Me)*Inv(Me)
        for (int k = 0; k < nnodes; ++k)  // loop over sum k
          for (_CI p = derivme[k][l].begin(); p != derivme[k][l].end(); ++p)
            derivae[p->first](i, j) -= ae(i, k) * meinv(l, j) * (p->second);
      }
    }
  }

  // store ae matrix in slave element data container
  SlaveElement().MoData().DualShape() = Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix(ae));

  return;
}

FOUR_C_NAMESPACE_CLOSE
