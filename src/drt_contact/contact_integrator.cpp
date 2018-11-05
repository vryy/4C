/*---------------------------------------------------------------------*/
/*!
\file contact_integrator.cpp

\brief A class to perform integrations of Mortar matrices on the overlap
       of two MortarElements in 1D and 2D (derived version for contact)

\level 2

\maintainer Matthias Mayr

*/
/*---------------------------------------------------------------------*/
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "contact_integrator.H"
#include "contact_node.H"
#include "contact_element.H"
#include "contact_defines.H"
#include "contact_paramsinterface.H"
#include "friction_node.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_mortar/mortar_projector.H"
#include "../drt_mortar/mortar_coupling3d_classes.H"
#include "../drt_mortar/mortar_calc_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_wear.H"

// headers for poro contact integration
#include "../drt_mat/structporo.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 10/13|
 *----------------------------------------------------------------------*/
CONTACT::CoIntegrator::CoIntegrator(Teuchos::ParameterList& params,
    DRT::Element::DiscretizationType eletype, const Epetra_Comm& comm)
    : imortar_(params),
      Comm_(comm),
      dim_(imortar_.get<int>("DIMENSION")),
      shapefcn_(DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(imortar_, "LM_SHAPEFCN")),
      lagmultquad_(DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(imortar_, "LM_QUAD")),
      gpslip_(DRT::INPUT::IntegralValue<int>(imortar_, "GP_SLIP_INCR")),
      algo_(DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(imortar_, "ALGORITHM")),
      stype_(DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(imortar_, "STRATEGY")),
      cppnormal_(DRT::INPUT::IntegralValue<int>(imortar_, "CPP_NORMALS")),
      wearlaw_(DRT::INPUT::IntegralValue<INPAR::WEAR::WearLaw>(imortar_, "WEARLAW")),
      wearimpl_(false),
      wearside_(INPAR::WEAR::wear_slave),
      weartype_(INPAR::WEAR::wear_intstate),
      wearshapefcn_(INPAR::WEAR::wear_shape_standard),
      sswear_(DRT::INPUT::IntegralValue<int>(imortar_, "SSWEAR")),
      wearcoeff_(-1.0),
      wearcoeffm_(-1.0),
      ssslip_(imortar_.get<double>("SSSLIP")),
      nonsmooth_(false),
      integrationtype_(DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(imortar_, "INTTYPE"))
{
  // init gp
  InitializeGP(eletype);

  // wear specific
  if (wearlaw_ != INPAR::WEAR::wear_none)
  {
    // set wear contact status
    INPAR::WEAR::WearTimInt wtimint =
        DRT::INPUT::IntegralValue<INPAR::WEAR::WearTimInt>(params, "WEARTIMINT");
    if (wtimint == INPAR::WEAR::wear_impl) wearimpl_ = true;

    // wear surface
    wearside_ = DRT::INPUT::IntegralValue<INPAR::WEAR::WearSide>(imortar_, "WEAR_SIDE");

    // wear algorithm
    weartype_ = DRT::INPUT::IntegralValue<INPAR::WEAR::WearType>(imortar_, "WEARTYPE");

    // wear shape function
    wearshapefcn_ = DRT::INPUT::IntegralValue<INPAR::WEAR::WearShape>(imortar_, "WEAR_SHAPEFCN");

    // wear coefficient
    wearcoeff_ = imortar_.get<double>("WEARCOEFF");

    // wear coefficient
    wearcoeffm_ = imortar_.get<double>("WEARCOEFF_MASTER");
  }

  nonsmooth_ = DRT::INPUT::IntegralValue<int>(imortar_, "NONSMOOTH_GEOMETRIES");

  return;
}

/*----------------------------------------------------------------------*
 |  check for boundary elements                              farah 02/14|
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::BoundarySegmCheck2D(
    MORTAR::MortarElement& sele, std::vector<MORTAR::MortarElement*> meles)
{
  double sxi_test[2] = {0.0, 0.0};
  bool proj_test = false;
  bool boundary_ele = false;

  double glob_test[3] = {0.0, 0.0, 0.0};

  DRT::Node** mynodes_test = sele.Nodes();
  if (!mynodes_test) dserror("ERROR: HasProjStatus: Null pointer!");

  if (sele.Shape() == DRT::Element::line2 || sele.Shape() == DRT::Element::nurbs2)
  {
    for (int s_test = 0; s_test < 2; ++s_test)
    {
      if (s_test == 0)
        sxi_test[0] = -1.0;
      else
        sxi_test[0] = 1.0;

      proj_test = false;
      for (int bs_test = 0; bs_test < (int)meles.size(); ++bs_test)
      {
        double mxi_test[2] = {0.0, 0.0};
        MORTAR::MortarProjector::Impl(sele, *meles[bs_test])
            ->ProjectGaussPoint2D(sele, sxi_test, *meles[bs_test], mxi_test);

        if ((mxi_test[0] >= -1.0) && (mxi_test[0] <= 1.0))
        {
          // get hasproj
          sele.LocalToGlobal(sxi_test, glob_test, 0);
          for (int ii = 0; ii < sele.NumNode(); ++ii)
          {
            MORTAR::MortarNode* mycnode_test = dynamic_cast<MORTAR::MortarNode*>(mynodes_test[ii]);
            if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

            if (glob_test[0] == mycnode_test->xspatial()[0] &&
                glob_test[1] == mycnode_test->xspatial()[1] &&
                glob_test[2] == mycnode_test->xspatial()[2])
              mycnode_test->HasProj() = true;
          }

          glob_test[0] = 0.0;
          glob_test[1] = 0.0;
          glob_test[2] = 0.0;

          proj_test = true;
        }
      }
      if (proj_test == false) boundary_ele = true;
    }
  }
  else if (sele.Shape() == DRT::Element::line3 || sele.Shape() == DRT::Element::nurbs3)
  {
    for (int s_test = 0; s_test < 3; ++s_test)
    {
      if (s_test == 0)
        sxi_test[0] = -1.0;
      else if (s_test == 1)
        sxi_test[0] = 0.0;
      else if (s_test == 2)
        sxi_test[0] = 1.0;

      proj_test = false;
      for (int bs_test = 0; bs_test < (int)meles.size(); ++bs_test)
      {
        double mxi_test[2] = {0.0, 0.0};
        MORTAR::MortarProjector::Impl(sele, *meles[bs_test])
            ->ProjectGaussPoint2D(sele, sxi_test, *meles[bs_test], mxi_test);

        if ((mxi_test[0] >= -1.0) && (mxi_test[0] <= 1.0))
        {
          // get hasproj
          sele.LocalToGlobal(sxi_test, glob_test, 0);
          for (int ii = 0; ii < sele.NumNode(); ++ii)
          {
            MORTAR::MortarNode* mycnode_test = dynamic_cast<MORTAR::MortarNode*>(mynodes_test[ii]);
            if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

            if (glob_test[0] == mycnode_test->xspatial()[0] &&
                glob_test[1] == mycnode_test->xspatial()[1] &&
                glob_test[2] == mycnode_test->xspatial()[2])
              mycnode_test->HasProj() = true;
          }

          glob_test[0] = 0.0;
          glob_test[1] = 0.0;
          glob_test[2] = 0.0;

          proj_test = true;
        }
      }
      if (proj_test == false) boundary_ele = true;
    }
  }
  else
  {
    dserror("No valid element type for slave discretization!");
  }

  return boundary_ele;
}


/*----------------------------------------------------------------------*
 |  Initialize gauss points                                   popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::InitializeGP(DRT::Element::DiscretizationType eletype)
{
  //**********************************************************************
  // Create integration points according to eletype!
  //
  // For segment-based integration, we have pre-defined default
  // values for the Gauss rules according to the segment type.
  //
  // default for integrals on 1D lines:
  // --> 5 GP (degree of precision: 9)
  //
  // default for integrals on 2D triangles:
  // --> 7 GP (degree of precision: 5)
  //
  // default for integrals on 2D quadrilaterals:
  // --> 9 GP (degree of precision: 5)
  //
  // For element-based integration, we choose the Gauss rules according
  // to the user's wish (i.e. according to the parameter NUMGP_PER_DIM).
  //
  // possibilites for integrals on 1D lines:
  // --> 1,2,3,4,5,6,7,8,9,10,16,20,32 GPs
  //
  // possibilities for integrals on 2D triangles:
  // --> 1,3,6,7,12,37,64 GPs
  //
  // possibilities for integrals on 2D quadrilaterals
  // --> 1,4,9,16,25,36,49,64,81,100,256,400,1024 GPs
  //**********************************************************************

  // get numgp (for element-based integration)
  int numgp = imortar_.get<int>("NUMGP_PER_DIM");

  // get integration type
  INPAR::MORTAR::IntType integrationtype =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(imortar_, "INTTYPE");

  //**********************************************************************
  // choose Gauss rule according to (a) element type (b) input parameter
  //**********************************************************************
  switch (eletype)
  {
    case DRT::Element::line2:
    case DRT::Element::line3:
    case DRT::Element::nurbs2:
    case DRT::Element::nurbs3:
    {
      // set default value for segment-based version first
      DRT::UTILS::GaussRule1D mygaussrule = DRT::UTILS::intrule_line_5point;

      // GP switch if element-based version and non-zero value provided by user
      if (integrationtype == INPAR::MORTAR::inttype_elements ||
          integrationtype == INPAR::MORTAR::inttype_elements_BS)
      {
        if (numgp > 0)
        {
          switch (numgp)
          {
            case 1:
            {
              dserror("Our experience says that 1 GP per slave element is not enough.");
              break;
            }
            case 2:
            {
              mygaussrule = DRT::UTILS::intrule_line_2point;
              break;
            }
            case 3:
            {
              mygaussrule = DRT::UTILS::intrule_line_3point;
              break;
            }
            case 4:
            {
              mygaussrule = DRT::UTILS::intrule_line_4point;
              break;
            }
            case 5:
            {
              mygaussrule = DRT::UTILS::intrule_line_5point;
              break;
            }
            case 6:
            {
              mygaussrule = DRT::UTILS::intrule_line_6point;
              break;
            }
            case 7:
            {
              mygaussrule = DRT::UTILS::intrule_line_7point;
              break;
            }
            case 8:
            {
              mygaussrule = DRT::UTILS::intrule_line_8point;
              break;
            }
            case 9:
            {
              mygaussrule = DRT::UTILS::intrule_line_9point;
              break;
            }
            case 10:
            {
              mygaussrule = DRT::UTILS::intrule_line_10point;
              break;
            }
            case 16:
            {
              mygaussrule = DRT::UTILS::intrule_line_16point;
              break;
            }
            case 20:
            {
              mygaussrule = DRT::UTILS::intrule_line_20point;
              break;
            }
            case 32:
            {
              mygaussrule = DRT::UTILS::intrule_line_32point;
              break;
            }
            case 50:
            {
              mygaussrule = DRT::UTILS::intrule_line_50point;
              break;
            }
            default:
            {
              dserror("Requested GP-Number is not implemented!");
              break;
            }
          }
        }
      }

      const DRT::UTILS::IntegrationPoints1D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.Reshape(nGP(), 2);
      weights_.resize(nGP());
      for (int i = 0; i < nGP(); ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = 0.0;
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case DRT::Element::tri3:
    case DRT::Element::tri6:
    {
      // set default value for segment-based version first
      DRT::UTILS::GaussRule2D mygaussrule = DRT::UTILS::intrule_tri_7point;
      if (integrationtype == INPAR::MORTAR::inttype_segments)
      {
        if (numgp > 0) switch (numgp)
          {
            case 1:
              mygaussrule = DRT::UTILS::intrule_tri_1point;
              break;
            case 3:
              mygaussrule = DRT::UTILS::intrule_tri_3point;
              break;
            case 7:
              mygaussrule = DRT::UTILS::intrule_tri_7point;
              break;
            case 16:
              mygaussrule = DRT::UTILS::intrule_tri_16point;
              break;
            case 37:
              mygaussrule = DRT::UTILS::intrule_tri_37point;
              break;
            default:
              dserror("unknown tri gauss rule");
              break;
          }
      }

      // GP switch if element-based version and non-zero value provided by user
      if (integrationtype == INPAR::MORTAR::inttype_elements ||
          integrationtype == INPAR::MORTAR::inttype_elements_BS)
      {
        if (numgp > 0)
        {
          switch (numgp)
          {
            case 1:
            {
              mygaussrule = DRT::UTILS::intrule_tri_3point;
              break;
            }
            case 2:
            {
              mygaussrule = DRT::UTILS::intrule_tri_6point;
              break;
            }
            case 3:
            {
              mygaussrule = DRT::UTILS::intrule_tri_7point;
              break;
            }
            case 4:
            {
              mygaussrule = DRT::UTILS::intrule_tri_12point;
              break;
            }
            case 5:
            {
              mygaussrule = DRT::UTILS::intrule_tri_12point;
              break;
            }
            case 6:
            {
              mygaussrule = DRT::UTILS::intrule_tri_37point;
              break;
            }
            case 7:
            {
              mygaussrule = DRT::UTILS::intrule_tri_37point;
              break;
            }
            case 8:
            {
              mygaussrule = DRT::UTILS::intrule_tri_64point;
              break;
            }
            case 9:
            {
              mygaussrule = DRT::UTILS::intrule_tri_64point;
              break;
            }
            case 10:
            {
              mygaussrule = DRT::UTILS::intrule_tri_64point;
              break;
            }
            case 20:
            {
              mygaussrule = DRT::UTILS::intrule_tri_64point;
              break;
            }
            default:
            {
              dserror("Requested GP-Number is not implemented!");
              break;
            }
          }
        }
      }

      const DRT::UTILS::IntegrationPoints2D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.Reshape(nGP(), 2);
      weights_.resize(nGP());
      for (int i = 0; i < nGP(); ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case DRT::Element::quad4:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
    case DRT::Element::nurbs4:
    case DRT::Element::nurbs8:
    case DRT::Element::nurbs9:
    {
      // set default value for segment-based version first
      DRT::UTILS::GaussRule2D mygaussrule = DRT::UTILS::intrule_quad_9point;

      // GP switch if element-based version and non-zero value provided by user
      if (integrationtype == INPAR::MORTAR::inttype_elements ||
          integrationtype == INPAR::MORTAR::inttype_elements_BS)
      {
        if (numgp > 0)
        {
          switch (numgp)
          {
            case 1:
            {
              mygaussrule = DRT::UTILS::intrule_quad_1point;
              break;
            }
            case 2:
            {
              mygaussrule = DRT::UTILS::intrule_quad_4point;
              break;
            }
            case 3:
            {
              mygaussrule = DRT::UTILS::intrule_quad_9point;
              break;
            }
            case 4:
            {
              mygaussrule = DRT::UTILS::intrule_quad_16point;
              break;
            }
            case 5:
            {
              mygaussrule = DRT::UTILS::intrule_quad_25point;
              break;
            }
            case 6:
            {
              mygaussrule = DRT::UTILS::intrule_quad_36point;
              break;
            }
            case 7:
            {
              mygaussrule = DRT::UTILS::intrule_quad_49point;
              break;
            }
            case 8:
            {
              mygaussrule = DRT::UTILS::intrule_quad_64point;
              break;
            }
            case 9:
            {
              mygaussrule = DRT::UTILS::intrule_quad_81point;
              break;
            }
            case 10:
            {
              mygaussrule = DRT::UTILS::intrule_quad_100point;
              break;
            }
            case 16:
            {
              mygaussrule = DRT::UTILS::intrule_quad_256point;
              break;
            }
            case 20:
            {
              mygaussrule = DRT::UTILS::intrule_quad_400point;
              break;
            }
            case 32:
            {
              mygaussrule = DRT::UTILS::intrule_quad_1024point;
              break;
            }
            default:
            {
              dserror("Requested GP-Number is not implemented!");
              break;
            }
          }
        }
      }

      const DRT::UTILS::IntegrationPoints2D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.Reshape(nGP(), 2);
      weights_.resize(nGP());
      for (int i = 0; i < nGP(); ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    default:
    {
      dserror("ERROR: MortarIntegrator: This contact element type is not implemented!");
      break;
    }
  }  // switch(eletype)

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateDerivSegment2D(MORTAR::MortarElement& sele, double& sxia,
    double& sxib, MORTAR::MortarElement& mele, double& mxia, double& mxib, const Epetra_Comm& comm,
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr = Teuchos::null;
  if (not mparams_ptr.is_null())
  {
    cparams_ptr = Teuchos::rcp_dynamic_cast<CONTACT::ParamsInterface>(mparams_ptr, true);
  }
  IntegrateDerivSegment2D(sele, sxia, sxib, mele, mxia, mxib, comm, cparams_ptr);
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize a 1D slave / master overlap (2D)  popp 02/09|
 |  This method integrates the overlap M matrix and weighted gap g~     |
 |  and stores it in mseg and gseg respectively. Moreover, derivatives  |
 |  LinM and Ling are built and stored directly into the adjacent nodes.|
 |  (Thus this method combines EVERYTHING before done separately in     |
 |  IntegrateM, IntegrateG, DerivM and DerivG!)                         |
 |  Also wear is integrated.                                            |
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateDerivSegment2D(MORTAR::MortarElement& sele, double& sxia,
    double& sxib, MORTAR::MortarElement& mele, double& mxia, double& mxib, const Epetra_Comm& comm,
    const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  // skip this segment, if too small
  if (sxib - sxia < 4. * MORTARINTLIM) return;

  // *********************************************************************
  // Check integrator input for non-reasonable quantities
  // *********************************************************************

  // explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivSegment2D called without specific shape function defined!");

  // Petrov-Galerkin approach for LM not yet implemented for quadratic FE
  if (sele.Shape() == MORTAR::MortarElement::line3 &&
      ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
    dserror(
        "ERROR: Petrov-Galerkin approach not yet implemented for 2-D quadratic FE interpolation");

  // check for problem dimension
  if (Dim() != 2) dserror("ERROR: 2D integration method called for non-2D problem");

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateDerivSegment2D called on a wrong type of MortarElement pair!");
  if ((sxia < -1.0) || (sxib > 1.0))
    dserror("ERROR: IntegrateDerivSegment2D called with infeasible slave limits!");
  if ((mxia < -1.0) || (mxib > 1.0))
    dserror("ERROR: IntegrateDerivSegment2D called with infeasible master limits!");

  // *********************************************************************
  // Prepare integration
  // *********************************************************************
  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = Dim();

  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if (!mynodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k = 0; k < nrow; ++k)
  {
    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  // decide whether linear LM are used for quadratic FE here
  // and whether displacement shape fct. modification has to be considered or not
  // this is the case for dual linear Lagrange multipliers on line3 elements
  bool linlm = false;
  bool dualquad = false;
  if (LagMultQuad() == INPAR::MORTAR::lagmult_lin && sele.Shape() == DRT::Element::line3)
  {
    bound = false;  // crosspoints and linear LM NOT at the same time!!!!
    linlm = true;
    if (shapefcn_ == INPAR::MORTAR::shape_dual) dualquad = true;
  }

  // prepare directional derivative of dual shape functions
  GEN::pairedvector<int, Epetra_SerialDenseMatrix> dualmap(
      2 * nrow, 0, Epetra_SerialDenseMatrix(nrow, nrow));
  if ((ShapeFcn() == INPAR::MORTAR::shape_dual ||
          ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) &&
      (sele.Shape() == MORTAR::MortarElement::line3 ||
          sele.Shape() == MORTAR::MortarElement::nurbs3 ||
          sele.MoData().DerivDualShape() != Teuchos::null))
    sele.DerivShapeDual(dualmap);

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow, 1);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol, 1);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow, 1);

  LINALG::SerialDenseVector m2val(ncol);
  LINALG::SerialDenseMatrix m2deriv(ncol, 1);
  LINALG::SerialDenseVector lm2val(ncol);
  LINALG::SerialDenseMatrix lm2deriv(ncol, 1);

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseMatrix ssecderiv(nrow, 1);

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  int linsize = 0;
  for (int i = 0; i < nrow; ++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*>(mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  // safety
  linsize = linsize * 2;

  // *********************************************************************
  // Find out about whether start / end of overlap are slave or master!
  // CAUTION: be careful with positive rotation direction ("Umlaufsinn")
  // sxia -> belongs to sele.Nodes()[0]
  // sxib -> belongs to sele.Nodes()[1]
  // mxia -> belongs to mele.Nodes()[0]
  // mxib -> belongs to mele.Nodes()[1]
  // but slave and master have different positive rotation directions,
  // counter-clockwise for slave side, clockwise for master side!
  // this means that mxia belongs to sxib and vice versa!
  // *********************************************************************

  bool startslave = false;
  bool endslave = false;

  if (sele.NormalFac() * mele.NormalFac() > 0.)
  {
    if (sxia != -1.0 && mxib != 1.0)
      dserror("ERROR: First outer node is neither slave nor master node");
    if (sxib != 1.0 && mxia != -1.0)
      dserror("ERROR: Second outer node is neither slave nor master node");
  }
  else
  {
    if (sxia != -1. && mxia != -1.)
      dserror("ERROR: First outer node is neither slave nor master node");
    if (sxib != 1. && mxib != 1.)
      dserror("ERROR: Second outer node is neither slave nor master node");
  }
  if (sxia == -1.0)
    startslave = true;
  else
    startslave = false;
  if (sxib == 1.0)
    endslave = true;
  else
    endslave = false;

  // get directional derivatives of sxia, sxib, mxia, mxib
  std::vector<GEN::pairedvector<int, double>> ximaps(4, linsize + ndof * ncol);
  DerivXiAB2D(sele, sxia, sxib, mele, mxia, mxib, ximaps, startslave, endslave, linsize);

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < nGP(); ++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp, 0), 0.0};
    double wgt = Weight(gp);

    // coordinate transformation sxi->eta (slave MortarElement->Overlap)
    double sxi[2] = {0.0, 0.0};
    sxi[0] = 0.5 * (1.0 - eta[0]) * sxia + 0.5 * (1.0 + eta[0]) * sxib;

    // project Gauss point onto master element
    double mxi[2] = {0.0, 0.0};
    MORTAR::MortarProjector::Impl(sele, mele)->ProjectGaussPoint2D(sele, sxi, mele, mxi);

    // check GP projection
    if ((mxi[0] < mxia - 1e-4) || (mxi[0] > mxib + 1e-4))
    {
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Gauss point: " << sxi[0] << " " << sxi[1] << std::endl;
      std::cout << "Projection: " << mxi[0] << " " << mxi[1] << std::endl;
      std::cout << "Projection bounds: " << mxia << " " << mxib << std::endl;
      // dserror("ERROR: IntegrateAndDerivSegment: Gauss point projection failed! mxi=%d",mxi[0]);
    }

    // evaluate Lagrange multiplier shape functions (on slave element)
    if (LagMultQuad() == INPAR::MORTAR::lagmult_const)
      sele.EvaluateShapeLagMultConst(ShapeFcn(), sxi, lmval, lmderiv, nrow);
    else if (linlm)
      sele.EvaluateShapeLagMultLin(ShapeFcn(), sxi, lmval, lmderiv, nrow);
    else
      sele.EvaluateShapeLagMult(ShapeFcn(), sxi, lmval, lmderiv, nrow);

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi, sval, sderiv, nrow, dualquad);
    mele.EvaluateShape(mxi, mval, mderiv, ncol, false);

    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5 * sxia + 0.5 * sxib;

    // evaluate linearizations *******************************************
    // evaluate 2nd deriv of trace space shape functions (on slave element)
    sele.Evaluate2ndDerivShape(sxi, ssecderiv, nrow);

    // evaluate the derivative dxdsxidsxi = Jac,xi
    double djacdxi[2] = {0.0, 0.0};
    dynamic_cast<CONTACT::CoElement&>(sele).DJacDXi(djacdxi, sxi, ssecderiv);
    double dxdsxidsxi = djacdxi[0];  // only 2D here

    // evalute the GP slave coordinate derivatives
    std::vector<GEN::pairedvector<int, double>> dsxigp(
        1, GEN::pairedvector<int, double>(linsize + ndof * ncol));
    for (_CI p = ximaps[0].begin(); p != ximaps[0].end(); ++p)
      dsxigp[0][p->first] += 0.5 * (1.0 - eta[0]) * (p->second);
    for (_CI p = ximaps[1].begin(); p != ximaps[1].end(); ++p)
      dsxigp[0][p->first] += 0.5 * (1.0 + eta[0]) * (p->second);

    // evalute the GP master coordinate derivatives
    std::vector<GEN::pairedvector<int, double>> dmxigp(
        1, GEN::pairedvector<int, double>(linsize + ndof * ncol));
    DerivXiGP2D(sele, mele, sxi[0], mxi[0], dsxigp[0], dmxigp[0], linsize);

    // evaluate the Jacobian derivative
    GEN::pairedvector<int, double> derivjacSele(nrow * ndof);
    sele.DerivJacobian(sxi, derivjacSele);

    double jac = dsxideta * dxdsxi;
    GEN::pairedvector<int, double> derivjac(nrow * ndof + linsize);
    for (GEN::pairedvector<int, double>::const_iterator p = derivjacSele.begin();
         p != derivjacSele.end(); ++p)
      derivjac[p->first] += dsxideta * p->second;
    for (_CI p = ximaps[0].begin(); p != ximaps[0].end(); ++p)
      derivjac[p->first] -= .5 * dxdsxi * p->second;
    for (_CI p = ximaps[1].begin(); p != ximaps[1].end(); ++p)
      derivjac[p->first] += .5 * dxdsxi * p->second;
    for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
      derivjac[p->first] += dsxideta * dxdsxidsxi * p->second;

    //**********************************************************************
    // frequently reused quantities
    //**********************************************************************
    double gpn[2] = {0.0, 0.0};  // normalized normal at gp
    double gap[1] = {0.0};       // gap
    GEN::pairedvector<int, double> dgapgp(
        linsize + ndof * ncol);  // gap lin without weighting and jac
    std::vector<GEN::pairedvector<int, double>> dnmap_unit(
        2, (linsize + ndof * ncol));  // deriv of x and y comp. of gpn (unit)

    //**********************************************************************
    // evaluate at GP and lin char. quantities
    //**********************************************************************
    if (dualquad && linlm)
    {
      // declare and compute shape functions as well as derivatives for computation of gap
      // -> standard quadratic functions instead of only linear functions
      LINALG::SerialDenseVector svalgap(nrow);
      LINALG::SerialDenseMatrix sderivgap(nrow, 1);
      sele.EvaluateShape(sxi, svalgap, sderivgap, nrow);

      Gap_2D(sele, mele, svalgap, mval, sderivgap, mderiv, gap, gpn, dsxigp, dmxigp, dgapgp,
          dnmap_unit);
    }
    else
      Gap_2D(sele, mele, sval, mval, sderiv, mderiv, gap, gpn, dsxigp, dmxigp, dgapgp, dnmap_unit);

    IntegrateGP_2D(sele, mele, sval, lmval, mval, sderiv, mderiv, lmderiv, dualmap, wgt, jac,
        derivjac, gpn, dnmap_unit, *gap, dgapgp, sxi, mxi, dsxigp, dmxigp);

  }  // gp-loop

  return;
}

/*----------------------------------------------------------------------*
 |  check for boundary elements                              farah 07/14|
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::BoundarySegmCheck3D(
    MORTAR::MortarElement& sele, std::vector<MORTAR::MortarElement*> meles)
{
  double sxi_test[2] = {0.0, 0.0};
  bool proj_test = false;
  bool boundary_ele = false;
  double alpha_test = 0.0;
  double glob_test[3] = {0.0, 0.0, 0.0};
  const double tol = 1e-8;
  DRT::Node** mynodes_test = sele.Nodes();
  if (!mynodes_test) dserror("ERROR: HasProjStatus: Null pointer!");

  DRT::Element::DiscretizationType dt_s = sele.Shape();

  if (dt_s == DRT::Element::quad4)  //|| dt_s==DRT::Element::quad8 || dt_s==DRT::Element::quad9)
  {
    for (int s_test = 0; s_test < 4; ++s_test)
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

      proj_test = false;
      for (int bs_test = 0; bs_test < (int)meles.size(); ++bs_test)
      {
        double mxi_test[2] = {0.0, 0.0};
        MORTAR::MortarProjector::Impl(sele, *meles[bs_test])
            ->ProjectGaussPoint3D(sele, sxi_test, *meles[bs_test], mxi_test, alpha_test);
        DRT::Element::DiscretizationType dt = meles[bs_test]->Shape();

        if (dt == DRT::Element::quad4 || dt == DRT::Element::quad8 || dt == DRT::Element::quad9)
        {
          if (mxi_test[0] >= -1.0 && mxi_test[1] >= -1.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0)
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.NumNode(); ++ii)
            {
              MORTAR::MortarNode* mycnode_test =
                  dynamic_cast<MORTAR::MortarNode*>(mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0] > mycnode_test->xspatial()[0] - tol &&
                  glob_test[0] < mycnode_test->xspatial()[0] + tol &&
                  glob_test[1] > mycnode_test->xspatial()[1] - tol &&
                  glob_test[1] < mycnode_test->xspatial()[1] + tol &&
                  glob_test[2] > mycnode_test->xspatial()[2] - tol &&
                  glob_test[2] < mycnode_test->xspatial()[2] + tol)
                mycnode_test->HasProj() = true;
            }

            glob_test[0] = 0.0;
            glob_test[1] = 0.0;
            glob_test[2] = 0.0;

            proj_test = true;
          }
        }
        else if (dt == DRT::Element::tri3 || dt == DRT::Element::tri6)
        {
          if (mxi_test[0] >= 0.0 && mxi_test[1] >= 0.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0 && mxi_test[0] + mxi_test[1] <= 1.0)
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.NumNode(); ++ii)
            {
              MORTAR::MortarNode* mycnode_test =
                  dynamic_cast<MORTAR::MortarNode*>(mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0] > mycnode_test->xspatial()[0] - tol &&
                  glob_test[0] < mycnode_test->xspatial()[0] + tol &&
                  glob_test[1] > mycnode_test->xspatial()[1] - tol &&
                  glob_test[1] < mycnode_test->xspatial()[1] + tol &&
                  glob_test[2] > mycnode_test->xspatial()[2] - tol &&
                  glob_test[2] < mycnode_test->xspatial()[2] + tol)
                mycnode_test->HasProj() = true;
            }
            glob_test[0] = 0.0;
            glob_test[1] = 0.0;
            glob_test[2] = 0.0;

            proj_test = true;
          }
        }
        else
        {
          dserror("Non valid element type for master discretization!");
        }
      }
      if (proj_test == false) boundary_ele = true;
    }
  }
  else if (dt_s ==
           DRT::Element::quad9)  //|| dt_s==DRT::Element::quad8 || dt_s==DRT::Element::quad9)
  {
    for (int s_test = 0; s_test < 9; ++s_test)
    {
      if (s_test == 0)
      {
        sxi_test[0] = -1.0;
        sxi_test[1] = -1.0;
      }
      else if (s_test == 1)
      {
        sxi_test[0] = 0.0;
        sxi_test[1] = -1.0;
      }
      else if (s_test == 2)
      {
        sxi_test[0] = 1.0;
        sxi_test[1] = -1.0;
      }
      else if (s_test == 3)
      {
        sxi_test[0] = -1.0;
        sxi_test[1] = 0.0;
      }
      else if (s_test == 4)
      {
        sxi_test[0] = 0.0;
        sxi_test[1] = 0.0;
      }
      else if (s_test == 5)
      {
        sxi_test[0] = 1.0;
        sxi_test[1] = 0.0;
      }
      else if (s_test == 6)
      {
        sxi_test[0] = -1.0;
        sxi_test[1] = 1.0;
      }
      else if (s_test == 7)
      {
        sxi_test[0] = 0.0;
        sxi_test[1] = 1.0;
      }
      else if (s_test == 8)
      {
        sxi_test[0] = 1.0;
        sxi_test[1] = 1.0;
      }

      proj_test = false;
      for (int bs_test = 0; bs_test < (int)meles.size(); ++bs_test)
      {
        double mxi_test[2] = {0.0, 0.0};
        MORTAR::MortarProjector::Impl(sele, *meles[bs_test])
            ->ProjectGaussPoint3D(sele, sxi_test, *meles[bs_test], mxi_test, alpha_test);
        DRT::Element::DiscretizationType dt = meles[bs_test]->Shape();

        if (dt == DRT::Element::quad4 || dt == DRT::Element::quad8 || dt == DRT::Element::quad9)
        {
          if (mxi_test[0] >= -1.0 && mxi_test[1] >= -1.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0)
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.NumNode(); ++ii)
            {
              MORTAR::MortarNode* mycnode_test =
                  dynamic_cast<MORTAR::MortarNode*>(mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0] > mycnode_test->xspatial()[0] - tol &&
                  glob_test[0] < mycnode_test->xspatial()[0] + tol &&
                  glob_test[1] > mycnode_test->xspatial()[1] - tol &&
                  glob_test[1] < mycnode_test->xspatial()[1] + tol &&
                  glob_test[2] > mycnode_test->xspatial()[2] - tol &&
                  glob_test[2] < mycnode_test->xspatial()[2] + tol)
                mycnode_test->HasProj() = true;
            }

            glob_test[0] = 0.0;
            glob_test[1] = 0.0;
            glob_test[2] = 0.0;

            proj_test = true;
          }
        }
        else if (dt == DRT::Element::tri3 || dt == DRT::Element::tri6)
        {
          if (mxi_test[0] >= 0.0 && mxi_test[1] >= 0.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0 && mxi_test[0] + mxi_test[1] <= 1.0)
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.NumNode(); ++ii)
            {
              MORTAR::MortarNode* mycnode_test =
                  dynamic_cast<MORTAR::MortarNode*>(mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0] > mycnode_test->xspatial()[0] - tol &&
                  glob_test[0] < mycnode_test->xspatial()[0] + tol &&
                  glob_test[1] > mycnode_test->xspatial()[1] - tol &&
                  glob_test[1] < mycnode_test->xspatial()[1] + tol &&
                  glob_test[2] > mycnode_test->xspatial()[2] - tol &&
                  glob_test[2] < mycnode_test->xspatial()[2] + tol)
                mycnode_test->HasProj() = true;
            }
            glob_test[0] = 0.0;
            glob_test[1] = 0.0;
            glob_test[2] = 0.0;

            proj_test = true;
          }
        }
        else
        {
          dserror("Non valid element type for master discretization!");
        }
      }
      if (proj_test == false) boundary_ele = true;
    }
  }
  else if (dt_s ==
           DRT::Element::quad8)  //|| dt_s==DRT::Element::quad8 || dt_s==DRT::Element::quad9)
  {
    for (int s_test = 0; s_test < 8; ++s_test)
    {
      if (s_test == 0)
      {
        sxi_test[0] = -1.0;
        sxi_test[1] = -1.0;
      }
      else if (s_test == 1)
      {
        sxi_test[0] = 0.0;
        sxi_test[1] = -1.0;
      }
      else if (s_test == 2)
      {
        sxi_test[0] = 1.0;
        sxi_test[1] = -1.0;
      }
      else if (s_test == 3)
      {
        sxi_test[0] = -1.0;
        sxi_test[1] = 0.0;
      }
      else if (s_test == 4)
      {
        sxi_test[0] = 1.0;
        sxi_test[1] = 0.0;
      }
      else if (s_test == 5)
      {
        sxi_test[0] = -1.0;
        sxi_test[1] = 1.0;
      }
      else if (s_test == 6)
      {
        sxi_test[0] = 0.0;
        sxi_test[1] = 1.0;
      }
      else if (s_test == 7)
      {
        sxi_test[0] = 1.0;
        sxi_test[1] = 1.0;
      }

      proj_test = false;
      for (int bs_test = 0; bs_test < (int)meles.size(); ++bs_test)
      {
        double mxi_test[2] = {0.0, 0.0};
        MORTAR::MortarProjector::Impl(sele, *meles[bs_test])
            ->ProjectGaussPoint3D(sele, sxi_test, *meles[bs_test], mxi_test, alpha_test);
        DRT::Element::DiscretizationType dt = meles[bs_test]->Shape();

        if (dt == DRT::Element::quad4 || dt == DRT::Element::quad8 || dt == DRT::Element::quad9)
        {
          if (mxi_test[0] >= -1.0 && mxi_test[1] >= -1.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0)
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.NumNode(); ++ii)
            {
              MORTAR::MortarNode* mycnode_test =
                  dynamic_cast<MORTAR::MortarNode*>(mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0] > mycnode_test->xspatial()[0] - tol &&
                  glob_test[0] < mycnode_test->xspatial()[0] + tol &&
                  glob_test[1] > mycnode_test->xspatial()[1] - tol &&
                  glob_test[1] < mycnode_test->xspatial()[1] + tol &&
                  glob_test[2] > mycnode_test->xspatial()[2] - tol &&
                  glob_test[2] < mycnode_test->xspatial()[2] + tol)
                mycnode_test->HasProj() = true;
            }

            glob_test[0] = 0.0;
            glob_test[1] = 0.0;
            glob_test[2] = 0.0;

            proj_test = true;
          }
        }
        else if (dt == DRT::Element::tri3 || dt == DRT::Element::tri6)
        {
          if (mxi_test[0] >= 0.0 && mxi_test[1] >= 0.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0 && mxi_test[0] + mxi_test[1] <= 1.0)
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.NumNode(); ++ii)
            {
              MORTAR::MortarNode* mycnode_test =
                  dynamic_cast<MORTAR::MortarNode*>(mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0] > mycnode_test->xspatial()[0] - tol &&
                  glob_test[0] < mycnode_test->xspatial()[0] + tol &&
                  glob_test[1] > mycnode_test->xspatial()[1] - tol &&
                  glob_test[1] < mycnode_test->xspatial()[1] + tol &&
                  glob_test[2] > mycnode_test->xspatial()[2] - tol &&
                  glob_test[2] < mycnode_test->xspatial()[2] + tol)
                mycnode_test->HasProj() = true;
            }
            glob_test[0] = 0.0;
            glob_test[1] = 0.0;
            glob_test[2] = 0.0;

            proj_test = true;
          }
        }
        else
        {
          dserror("Non valid element type for master discretization!");
        }
      }
      if (proj_test == false) boundary_ele = true;
    }
  }
  // TRI-ELE
  else if (dt_s == DRT::Element::tri3)
  {
    for (int s_test = 0; s_test < 3; ++s_test)
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

      proj_test = false;
      for (int bs_test = 0; bs_test < (int)meles.size(); ++bs_test)
      {
        double mxi_test[2] = {0.0, 0.0};
        MORTAR::MortarProjector::Impl(sele, *meles[bs_test])
            ->ProjectGaussPoint3D(sele, sxi_test, *meles[bs_test], mxi_test, alpha_test);
        DRT::Element::DiscretizationType dt = meles[bs_test]->Shape();

        if (dt == DRT::Element::quad4 || dt == DRT::Element::quad8 || dt == DRT::Element::quad9)
        {
          if (mxi_test[0] >= -1.0 && mxi_test[1] >= -1.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0)  // Falls Position auf Element
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.NumNode(); ++ii)
            {
              MORTAR::MortarNode* mycnode_test =
                  dynamic_cast<MORTAR::MortarNode*>(mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0] > mycnode_test->xspatial()[0] - tol &&
                  glob_test[0] < mycnode_test->xspatial()[0] + tol &&
                  glob_test[1] > mycnode_test->xspatial()[1] - tol &&
                  glob_test[1] < mycnode_test->xspatial()[1] + tol &&
                  glob_test[2] > mycnode_test->xspatial()[2] - tol &&
                  glob_test[2] < mycnode_test->xspatial()[2] + tol)
                mycnode_test->HasProj() = true;
            }

            glob_test[0] = 0.0;
            glob_test[1] = 0.0;
            glob_test[2] = 0.0;

            proj_test = true;
          }
        }
        else if (dt == DRT::Element::tri3 || dt == DRT::Element::tri6)
        {
          if (mxi_test[0] >= 0.0 && mxi_test[1] >= 0.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0 && mxi_test[0] + mxi_test[1] <= 1.0)  // Falls Position auf Element
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.NumNode(); ++ii)
            {
              MORTAR::MortarNode* mycnode_test =
                  dynamic_cast<MORTAR::MortarNode*>(mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0] > mycnode_test->xspatial()[0] - tol &&
                  glob_test[0] < mycnode_test->xspatial()[0] + tol &&
                  glob_test[1] > mycnode_test->xspatial()[1] - tol &&
                  glob_test[1] < mycnode_test->xspatial()[1] + tol &&
                  glob_test[2] > mycnode_test->xspatial()[2] - tol &&
                  glob_test[2] < mycnode_test->xspatial()[2] + tol)
                mycnode_test->HasProj() = true;
            }

            glob_test[0] = 0.0;
            glob_test[1] = 0.0;
            glob_test[2] = 0.0;

            proj_test = true;
          }
        }
        else
        {
          dserror("Non valid element type for master discretization!");
        }
      }
      if (proj_test == false) boundary_ele = true;
    }
  }
  else if (dt_s == DRT::Element::tri6)
  {
    for (int s_test = 0; s_test < 6; ++s_test)
    {
      if (s_test == 0)
      {
        sxi_test[0] = 0.0;
        sxi_test[1] = 0.0;
      }
      else if (s_test == 1)
      {
        sxi_test[0] = 0.5;
        sxi_test[1] = 0.0;
      }
      else if (s_test == 2)
      {
        sxi_test[0] = 1.0;
        sxi_test[1] = 0.0;
      }
      else if (s_test == 3)
      {
        sxi_test[0] = 0.0;
        sxi_test[1] = 0.5;
      }
      else if (s_test == 4)
      {
        sxi_test[0] = 0.5;
        sxi_test[1] = 0.5;
      }
      else if (s_test == 5)
      {
        sxi_test[0] = 0.0;
        sxi_test[1] = 1.0;
      }

      proj_test = false;
      for (int bs_test = 0; bs_test < (int)meles.size(); ++bs_test)
      {
        double mxi_test[2] = {0.0, 0.0};
        MORTAR::MortarProjector::Impl(sele, *meles[bs_test])
            ->ProjectGaussPoint3D(sele, sxi_test, *meles[bs_test], mxi_test, alpha_test);
        DRT::Element::DiscretizationType dt = meles[bs_test]->Shape();

        if (dt == DRT::Element::quad4 || dt == DRT::Element::quad8 || dt == DRT::Element::quad9)
        {
          if (mxi_test[0] >= -1.0 && mxi_test[1] >= -1.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0)  // Falls Position auf Element
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.NumNode(); ++ii)
            {
              MORTAR::MortarNode* mycnode_test =
                  dynamic_cast<MORTAR::MortarNode*>(mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0] > mycnode_test->xspatial()[0] - tol &&
                  glob_test[0] < mycnode_test->xspatial()[0] + tol &&
                  glob_test[1] > mycnode_test->xspatial()[1] - tol &&
                  glob_test[1] < mycnode_test->xspatial()[1] + tol &&
                  glob_test[2] > mycnode_test->xspatial()[2] - tol &&
                  glob_test[2] < mycnode_test->xspatial()[2] + tol)
                mycnode_test->HasProj() = true;
            }

            glob_test[0] = 0.0;
            glob_test[1] = 0.0;
            glob_test[2] = 0.0;

            proj_test = true;
          }
        }
        else if (dt == DRT::Element::tri3 || dt == DRT::Element::tri6)
        {
          if (mxi_test[0] >= 0.0 && mxi_test[1] >= 0.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0 && mxi_test[0] + mxi_test[1] <= 1.0)  // Falls Position auf Element
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.NumNode(); ++ii)
            {
              MORTAR::MortarNode* mycnode_test =
                  dynamic_cast<MORTAR::MortarNode*>(mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0] > mycnode_test->xspatial()[0] - tol &&
                  glob_test[0] < mycnode_test->xspatial()[0] + tol &&
                  glob_test[1] > mycnode_test->xspatial()[1] - tol &&
                  glob_test[1] < mycnode_test->xspatial()[1] + tol &&
                  glob_test[2] > mycnode_test->xspatial()[2] - tol &&
                  glob_test[2] < mycnode_test->xspatial()[2] + tol)
                mycnode_test->HasProj() = true;
            }

            glob_test[0] = 0.0;
            glob_test[1] = 0.0;
            glob_test[2] = 0.0;

            proj_test = true;
          }
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
  {
    dserror("Non valid element type for slave discretization!");
  }

  return boundary_ele;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateDerivEle3D(MORTAR::MortarElement& sele,
    std::vector<MORTAR::MortarElement*> meles, bool* boundary_ele, bool* proj_,
    const Epetra_Comm& comm, const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr = Teuchos::null;
  if (not mparams_ptr.is_null())
  {
    cparams_ptr = Teuchos::rcp_dynamic_cast<CONTACT::ParamsInterface>(mparams_ptr, true);
  }
  IntegrateDerivEle3D(sele, meles, boundary_ele, proj_, comm, cparams_ptr);
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize for lin and quad elements        farah 01/13|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateDerivEle3D(MORTAR::MortarElement& sele,
    std::vector<MORTAR::MortarElement*> meles, bool* boundary_ele, bool* proj_,
    const Epetra_Comm& comm, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  // explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivEle3D called without specific shape function defined!");

  // Petrov-Galerkin approach for LM not yet implemented
  if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin && sele.Shape() != DRT::Element::nurbs9)
    dserror(
        "ERROR: Petrov-Galerkin approach not yet implemented for 3-D quadratic FE interpolation");

  // check for problem dimension
  if (Dim() != 3) dserror("ERROR: 3D integration method called for non-3D problem");

  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = sele.Nodes();
  if (!mynodes) dserror("ERROR: IntegrateDerivEle3D: Null pointer!");

  // check input data
  for (int test = 0; test < (int)meles.size(); ++test)
  {
    if ((!sele.IsSlave()) || (meles[test]->IsSlave()))
      dserror("ERROR: IntegrateDerivEle3D called on a wrong type of MortarElement pair!");
  }

  int msize = meles.size();
  int nrow = sele.NumNode();
  int ndof = dynamic_cast<MORTAR::MortarNode*>(sele.Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow, 2, true);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow, 2, true);

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  GEN::pairedvector<int, Epetra_SerialDenseMatrix> dualmap(
      nrow * ndof, 0, Epetra_SerialDenseMatrix(nrow, nrow));
  if ((ShapeFcn() == INPAR::MORTAR::shape_dual ||
          ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) &&
      (sele.Shape() != MORTAR::MortarElement::tri3 ||
          sele.MoData().DerivDualShape() != Teuchos::null) &&
      LagMultQuad() != INPAR::MORTAR::lagmult_const)
    sele.DerivShapeDual(dualmap);

  // check whether quadratic elements with linear interpolation are present
  bool linlm = false;
  if ((ShapeFcn() == INPAR::MORTAR::shape_dual ||
          ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) &&
      (LagMultQuad() == INPAR::MORTAR::lagmult_lin) &&
      (sele.Shape() == DRT::Element::quad8 || sele.Shape() == DRT::Element::tri6))
    linlm = true;

  //********************************************************************
  //  Boundary_segmentation test -- HasProj() check
  //  if a slave-node has no projection onto each master element
  //  --> Boundary_ele==true
  //********************************************************************
  INPAR::MORTAR::IntType integrationtype =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(imortar_, "INTTYPE");

  //************************************************************************
  // Boundary Segmentation check -- HasProj()-check
  //************************************************************************
  if (sele.Shape() != DRT::Element::nurbs4 and sele.Shape() != DRT::Element::nurbs9)
    *boundary_ele = BoundarySegmCheck3D(sele, meles);

  int linsize = 0;
  for (int i = 0; i < nrow; ++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*>(mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  // Start integration if fast integration should be used or if there is no boundary element
  // for the fast_BS integration
  if (*boundary_ele == false || integrationtype == INPAR::MORTAR::inttype_elements)
  {
    //**********************************************************************
    // loop over all Gauss points for integration
    //**********************************************************************
    for (int gp = 0; gp < nGP(); ++gp)
    {
      int iter_proj = 0;
      // coordinates and weight
      double eta[2] = {Coordinate(gp, 0), Coordinate(gp, 1)};
      double wgt = Weight(gp);

      // note that the third component of sxi is necessary!
      // (although it will always be 0.0 of course)
      double sxi[2] = {0.0, 0.0};
      double mxi[2] = {0.0, 0.0};
      double projalpha = 0.0;

      // get Gauss point in slave element coordinates
      sxi[0] = eta[0];
      sxi[1] = eta[1];

      bool is_on_mele = true;

      // evaluate Lagrange multiplier shape functions (on slave element)
      if (LagMultQuad() == INPAR::MORTAR::lagmult_const)
        sele.EvaluateShapeLagMultConst(ShapeFcn(), sxi, lmval, lmderiv, nrow);
      else if (LagMultQuad() == INPAR::MORTAR::lagmult_lin)
        sele.EvaluateShapeLagMultLin(ShapeFcn(), sxi, lmval, lmderiv, nrow);
      else
        sele.EvaluateShapeLagMult(ShapeFcn(), sxi, lmval, lmderiv, nrow);

      // evaluate trace space shape functions (on both elements)
      sele.EvaluateShape(sxi, sval, sderiv, nrow, linlm);

      // evaluate the two Jacobians (int. cell and slave element)
      double jacslave = sele.Jacobian(sxi);

      // evaluate linearizations *******************************************
      // evaluate the slave Jacobian derivative
      GEN::pairedvector<int, double> jacslavemap(nrow * ndof + linsize);
      sele.DerivJacobian(sxi, jacslavemap);

      //**********************************************************************
      // loop over all mele
      //**********************************************************************
      for (int nummaster = 0; nummaster < msize; ++nummaster)
      {
        DRT::Element::DiscretizationType dt = meles[nummaster]->Shape();

        int nmnode = meles[nummaster]->NumNode();
        LINALG::SerialDenseVector mval(nmnode);
        LINALG::SerialDenseMatrix mderiv(nmnode, 2, true);


        // project Gauss point onto master element
        MORTAR::MortarProjector::Impl(sele, *meles[nummaster])
            ->ProjectGaussPoint3D(sele, sxi, *meles[nummaster], mxi, projalpha);

        // check GP projection
        const double tol = 0.00;
        is_on_mele = true;
        if (dt == DRT::Element::quad4 || dt == DRT::Element::quad8 || dt == DRT::Element::quad9 ||
            dt == DRT::Element::nurbs9)
        {
          if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[0] > 1.0 + tol ||
              mxi[1] > 1.0 + tol)
          {
            is_on_mele = false;
          }
        }
        else
        {
          if (mxi[0] < -tol || mxi[1] < -tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol ||
              mxi[0] + mxi[1] > 1.0 + 2 * tol)
          {
            is_on_mele = false;
          }
        }

        // gp is valid
        if (is_on_mele == true)
        {
          *proj_ = true;
          iter_proj += 1;

          // get mval
          meles[nummaster]->EvaluateShape(mxi, mval, mderiv, nmnode);

          // evaluate the GP slave coordinate derivatives
          std::vector<GEN::pairedvector<int, double>> dsxigp(2, 0);
          std::vector<GEN::pairedvector<int, double>> dmxigp(2, 4 * linsize + nmnode * ndof);
          DerivXiGP3D(sele, *meles[nummaster], sxi, mxi, dsxigp, dmxigp, projalpha);

          //**********************************************************************
          // geometric quantities
          //**********************************************************************
          double gpn[3] = {0.0, 0.0, 0.0};
          GEN::pairedvector<int, double> dgapgp(
              (nmnode * ndof) + linsize);  // gap lin. without lm and jac.
          double gap[1] = {0.0};
          std::vector<GEN::pairedvector<int, double>> dnmap_unit(
              3, ((nmnode * ndof) + linsize));  // deriv of x,y and z comp. of gpn (unit)

          //**********************************************************************
          // evaluate at GP and lin char. quantities
          //**********************************************************************
          if (linlm)
          {
            // declare and compute shape functions as well as derivatives for computation of gap
            // -> standard quadratic functions instead of only linear functions
            LINALG::SerialDenseVector svalgap(nrow);
            LINALG::SerialDenseMatrix sderivgap(nrow, 2, true);
            sele.EvaluateShape(sxi, svalgap, sderivgap, nrow);

            Gap_3D(sele, *meles[nummaster], svalgap, mval, sderivgap, mderiv, gap, gpn, dsxigp,
                dmxigp, dgapgp, dnmap_unit);
          }
          else
            Gap_3D(sele, *meles[nummaster], sval, mval, sderiv, mderiv, gap, gpn, dsxigp, dmxigp,
                dgapgp, dnmap_unit);

          if (algo_ == INPAR::MORTAR::algorithm_mortar)
            dynamic_cast<CONTACT::CoElement&>(sele).PrepareMderiv(meles, nummaster);

          IntegrateGP_3D(sele, *meles[nummaster], sval, lmval, mval, sderiv, mderiv, lmderiv,
              dualmap, wgt, jacslave, jacslavemap, gpn, dnmap_unit, gap[0], dgapgp, sxi, mxi,
              dsxigp, dmxigp);

          // assemble m-matrix for this slave/master pair
          if (algo_ == INPAR::MORTAR::algorithm_mortar)
            dynamic_cast<CONTACT::CoElement&>(sele).AssembleMderivToNodes(*meles[nummaster]);

        }  // is_on_mele
        if (is_on_mele == true) break;
      }  // mele loop

      // warning, if an element which is declared not to be on the boundary by the above test
      // has non-projectable Gauss points
      if (is_on_mele == false && *boundary_ele == false)
      {
        std::cout << "*** warning *** Non-boundary element has non-projectable Gauss point \n";
        DRT::Node** snodes = sele.Nodes();
        for (int e = 0; e < sele.NumNode(); ++e)
        {
          CoNode* cnode = dynamic_cast<CoNode*>(snodes[e]);
          std::cout << "#eID " << sele.Id() << " | #nID " << cnode->Id() << ": " << cnode->X()[0]
                    << ", " << cnode->X()[1] << ", " << cnode->X()[2] << std::endl;
        }
      }

      // if one gp has counterparts on 2 elements --> non-uniqueness
      if (iter_proj > 1) dserror("Multiple feasible projections of one integration point!");
    }  // GP-loop
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateDerivCell3DAuxPlane(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, Teuchos::RCP<MORTAR::IntCell> cell, double* auxn,
    const Epetra_Comm& comm, const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr = Teuchos::null;
  if (not mparams_ptr.is_null())
  {
    cparams_ptr = Teuchos::rcp_dynamic_cast<CONTACT::ParamsInterface>(mparams_ptr, true);
  }
  IntegrateDerivCell3DAuxPlane(sele, mele, cell, auxn, comm, cparams_ptr);
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize a 2D slave / master cell (3D)     popp 03/09|
 |  This method integrates the cell M matrix and weighted gap g~        |
 |  and stores it in mseg and gseg respectively. Moreover, derivatives  |
 |  LinM and Ling are built and stored directly into the adjacent nodes.|
 |  (Thus this method combines EVERYTHING before done separately in     |
 |  IntegrateM3D, IntegrateG3D, DerivM3D and DerivG3D!)                 |
 |  This is the auxiliary plane coupling version!!!                     |
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateDerivCell3DAuxPlane(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, Teuchos::RCP<MORTAR::IntCell> cell, double* auxn,
    const Epetra_Comm& comm, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  // explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without specific shape function defined!");

  // check for problem dimension
  if (Dim() != 3) dserror("ERROR: 3D integration method called for non-3D problem");

  // discretization type of master element
  DRT::Element::DiscretizationType sdt = sele.Shape();
  DRT::Element::DiscretizationType mdt = mele.Shape();

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called on a wrong type of MortarElement pair!");
  if (cell == Teuchos::null)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without integration cell");

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = Dim();

  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = sele.Nodes();
  if (!mynodes) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow, 2, true);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol, 2, true);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow, 2, true);
  LINALG::SerialDenseVector lmvalScale(nrow);
  LINALG::SerialDenseMatrix lmderivScale(nrow, 2, true);



  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  int linsize = 0;
  if (cppnormal_)
  {
    for (int i = 0; i < ncol; ++i)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(mele.Nodes()[i]);
      linsize += cnode->GetLinsize();
    }
    linsize += 1 + mele.NumNode();
  }
  else
  {
    for (int i = 0; i < nrow; ++i)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(mynodes[i]);
      linsize += cnode->GetLinsize();
    }
  }

  // for cpp normal we have more lin entries
  linsize *= 10;

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  GEN::pairedvector<int, Epetra_SerialDenseMatrix> dualmap(
      (nrow + ncol) * ndof, 0, Epetra_SerialDenseMatrix(nrow, nrow));
  if ((ShapeFcn() == INPAR::MORTAR::shape_dual ||
          ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) &&
      (sele.Shape() != MORTAR::MortarElement::tri3 ||
          sele.MoData().DerivDualShape() != Teuchos::null))
    sele.DerivShapeDual(dualmap);

  // check if the cells are tri3
  // there's nothing wrong about other shapes, but as long as they are all
  // tri3 we can perform the jacobian calculation ( and its deriv) outside
  // the Gauss point loop
  if (cell->Shape() != DRT::Element::tri3)
    dserror("only tri3 integration cells at the moment. See comment in the code");

  double jac = cell->Jacobian();
  // directional derivative of cell Jacobian
  GEN::pairedvector<int, double> jacintcellmap((nrow + ncol) * ndof);
  cell->DerivJacobian(jacintcellmap);

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < nGP(); ++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp, 0), Coordinate(gp, 1)};
    double wgt = Weight(gp);

    // get global Gauss point coordinates
    double globgp[3] = {0.0, 0.0, 0.0};
    cell->LocalToGlobal(eta, globgp, 0);

    double sxi[2] = {0.0, 0.0};
    double mxi[2] = {0.0, 0.0};

    // project Gauss point onto slave element
    // project Gauss point onto master element
    double sprojalpha = 0.0;
    double mprojalpha = 0.0;
    MORTAR::MortarProjector::Impl(sele)->ProjectGaussPointAuxn3D(
        globgp, auxn, sele, sxi, sprojalpha);
    MORTAR::MortarProjector::Impl(mele)->ProjectGaussPointAuxn3D(
        globgp, auxn, mele, mxi, mprojalpha);

    // check GP projection (SLAVE)
    double tol = 0.01;
    if (sdt == DRT::Element::quad4 || sdt == DRT::Element::quad8 || sdt == DRT::Element::quad9)
    {
      if (sxi[0] < -1.0 - tol || sxi[1] < -1.0 - tol || sxi[0] > 1.0 + tol || sxi[1] > 1.0 + tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }
    else
    {
      if (sxi[0] < -tol || sxi[1] < -tol || sxi[0] > 1.0 + tol || sxi[1] > 1.0 + tol ||
          sxi[0] + sxi[1] > 1.0 + 2 * tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }

    // check GP projection (MASTER)
    if (mdt == DRT::Element::quad4 || mdt == DRT::Element::quad8 || mdt == DRT::Element::quad9)
    {
      if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }
    else
    {
      if (mxi[0] < -tol || mxi[1] < -tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol ||
          mxi[0] + mxi[1] > 1.0 + 2 * tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }

    // evaluate Lagrange multiplier shape functions (on slave element)
    sele.EvaluateShapeLagMult(ShapeFcn(), sxi, lmval, lmderiv, nrow);
    sele.EvaluateShapeLagMult(ShapeFcn(), sxi, lmvalScale, lmderivScale, nrow, false);

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi, sval, sderiv, nrow);
    mele.EvaluateShape(mxi, mval, mderiv, ncol);

    // evaluate linearizations *******************************************

    // evaluate global GP coordinate derivative
    static LINALG::Matrix<3, 1> svalcell;
    static LINALG::Matrix<3, 2> sderivcell;
    cell->EvaluateShape(eta, svalcell, sderivcell);

    GEN::pairedvector<int, LINALG::Matrix<3, 1>> lingp((nrow + ncol) * ndof);

    for (int v = 0; v < 3; ++v)
      for (int d = 0; d < 3; ++d)
        for (_CI p = (cell->GetDerivVertex(v))[d].begin(); p != (cell->GetDerivVertex(v))[d].end();
             ++p)
          lingp[p->first](d) += svalcell(v) * (p->second);

    // evalute the GP slave coordinate derivatives
    std::vector<GEN::pairedvector<int, double>> dsxigp(2, (nrow + ncol) * ndof);
    DerivXiGP3DAuxPlane(sele, sxi, cell->Auxn(), dsxigp, sprojalpha, cell->GetDerivAuxn(), lingp);

    // evalute the GP master coordinate derivatives
    std::vector<GEN::pairedvector<int, double>> dmxigp(2, (nrow + ncol) * ndof);
    DerivXiGP3DAuxPlane(mele, mxi, cell->Auxn(), dmxigp, mprojalpha, cell->GetDerivAuxn(), lingp);

    //**********************************************************************
    // geometric quantities
    //**********************************************************************
    double gpn[3] = {0.0, 0.0, 0.0};
    GEN::pairedvector<int, double> dgapgp((ncol * ndof) + linsize);  // gap lin. without lm and jac.
    double gap[1] = {0.0};
    std::vector<GEN::pairedvector<int, double>> dnmap_unit(
        3, linsize);  // deriv of x,y and z comp. of gpn (unit)
    Gap_3D(sele, mele, sval, mval, sderiv, mderiv, gap, gpn, dsxigp, dmxigp, dgapgp, dnmap_unit);

    // integrate for area scaling:
    for (int i = 0; i < sele.NumNode(); ++i)
    {
      dynamic_cast<CONTACT::CoNode*>(sele.Nodes()[i])->MoData().GetDscale() +=
          wgt * jac * sval[i] * lmvalScale[i];
    }

    //**********************************************************************
    // evaluate at GP and lin char. quantities
    //**********************************************************************
    IntegrateGP_3D(sele, mele, sval, lmval, mval, sderiv, mderiv, lmderiv, dualmap, wgt, jac,
        jacintcellmap, gpn, dnmap_unit, gap[0], dgapgp, sxi, mxi, dsxigp, dmxigp);

  }  // end gp loop
  //**********************************************************************

  return;
}



/*----------------------------------------------------------------------*
 |  Integrate and linearize a 1D slave / master cell (2D)    farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateDerivCell3DAuxPlaneSTL(MORTAR::MortarElement& mele,
    MORTAR::MortarElement& lele, MORTAR::MortarElement& sele, Teuchos::RCP<MORTAR::IntCell> cell,
    double* auxn, const Epetra_Comm& comm)
{
  // explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without specific shape function defined!");

  // check for problem dimension
  if (Dim() != 3) dserror("ERROR: 3D integration method called for non-3D problem");

  // discretization type of master element
  DRT::Element::DiscretizationType sdt = sele.Shape();
  DRT::Element::DiscretizationType mdt = mele.Shape();

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called on a wrong type of MortarElement pair!");
  if (cell == Teuchos::null)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without integration cell");

  // number of nodes (slave, master)
  //  int nrowL = lele.NumNode();
  int nrow = sele.NumNode();
  int ncolP = mele.NumNode();
  int ncolL = lele.NumNode();
  int ndof = Dim();

  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = sele.Nodes();
  if (!mynodes) dserror("ERROR: IntegrateDerivCell3DAuxPlaneLTS: Null pointer!");
  DRT::Node** mnodes = lele.Nodes();
  if (!mnodes) dserror("ERROR: IntegrateDerivCell3DAuxPlaneLTS: Null pointer!");

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow, 2, true);
  LINALG::SerialDenseVector mval(ncolL);
  LINALG::SerialDenseMatrix mderiv(ncolL, 1, true);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow, 2, true);

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  GEN::pairedvector<int, Epetra_SerialDenseMatrix> dualmap(
      (nrow + ncolL) * ndof, 0, Epetra_SerialDenseMatrix(nrow, nrow));
  if ((ShapeFcn() == INPAR::MORTAR::shape_dual ||
          ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) &&
      (lele.Shape() != MORTAR::MortarElement::line2 ||
          sele.MoData().DerivDualShape() != Teuchos::null))
    sele.DerivShapeDual(dualmap);

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  int linsize = 0;
  for (int i = 0; i < nrow; ++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*>(mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  // check if the cells are tri3
  // there's nothing wrong about other shapes, but as long as they are all
  // tri3 we can perform the jacobian calculation ( and its deriv) outside
  // the Gauss point loop
  if (cell->Shape() != DRT::Element::line2) dserror("only line2 integration cells for LTS");

  double jac = cell->Jacobian();

  // directional derivative of cell Jacobian
  GEN::pairedvector<int, double> jacintcellmap((nrow + ncolL) * ndof);
  cell->DerivJacobian(jacintcellmap);

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < nGP(); ++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp, 0), Coordinate(gp, 1)};
    double wgt = Weight(gp);

    // get global Gauss point coordinates
    double globgp[3] = {0.0, 0.0, 0.0};
    cell->LocalToGlobal(eta, globgp, 0);

    // para coordinate on surface slave ele
    double sxi[2] = {0.0, 0.0};

    // para coordinate on surface master ele
    double mxi[2] = {0.0, 0.0};

    // para coordinate on line (edge) slave ele
    double lxi[2] = {0.0, 0.0};

    // project Gauss point onto slave/master surface element
    double sprojalpha = 0.0;
    double mprojalpha = 0.0;
    MORTAR::MortarProjector::Impl(sele)->ProjectGaussPointAuxn3D(
        globgp, auxn, sele, sxi, sprojalpha);
    MORTAR::MortarProjector::Impl(mele)->ProjectGaussPointAuxn3D(
        globgp, auxn, mele, mxi, mprojalpha);

    // check GP projection (SLAVE)
    double tol = 0.01;
    if (sdt == DRT::Element::quad4 || sdt == DRT::Element::quad8 || sdt == DRT::Element::quad9)
    {
      if (sxi[0] < -1.0 - tol || sxi[1] < -1.0 - tol || sxi[0] > 1.0 + tol || sxi[1] > 1.0 + tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }
    else
    {
      if (sxi[0] < -tol || sxi[1] < -tol || sxi[0] > 1.0 + tol || sxi[1] > 1.0 + tol ||
          sxi[0] + sxi[1] > 1.0 + 2 * tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }

    // check GP projection (MASTER)
    if (mdt == DRT::Element::quad4 || mdt == DRT::Element::quad8 || mdt == DRT::Element::quad9)
    {
      if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }
    else
    {
      if (mxi[0] < -tol || mxi[1] < -tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol ||
          mxi[0] + mxi[1] > 1.0 + 2 * tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }

    // convert sxi to lsxi!!!
    if (lele.Id() == 0)
    {
      lxi[0] = mxi[0];
    }
    else if (lele.Id() == 1)
    {
      lxi[0] = mxi[1];
    }
    else if (lele.Id() == 2)
    {
      lxi[0] = -mxi[0];
    }
    else if (lele.Id() == 3)
    {
      lxi[0] = -mxi[1];
    }
    else
      dserror("ERROR: invalid local id of line ele");

    // evaluate Lagrange multiplier shape functions (on slave element)
    sele.EvaluateShapeLagMult(ShapeFcn(), sxi, lmval, lmderiv, nrow);

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi, sval, sderiv, nrow);
    lele.EvaluateShape(lxi, mval, mderiv, ncolL);

    // evaluate linearizations *******************************************

    // evaluate global GP coordinate derivative
    static LINALG::Matrix<3, 1> svalcell;
    static LINALG::Matrix<3, 2> sderivcell;
    cell->EvaluateShape(eta, svalcell, sderivcell);

    GEN::pairedvector<int, LINALG::Matrix<3, 1>> lingp(100 * (nrow + ncolP) * ndof);

    for (int v = 0; v < 2; ++v)
      for (int d = 0; d < 3; ++d)
        for (_CI p = (cell->GetDerivVertex(v))[d].begin(); p != (cell->GetDerivVertex(v))[d].end();
             ++p)
          lingp[p->first](d) += svalcell(v) * (p->second);

    // evalute the GP slave coordinate derivatives
    std::vector<GEN::pairedvector<int, double>> dsxigp(2, (nrow + ncolP) * ndof);
    DerivXiGP3DAuxPlane(sele, sxi, cell->Auxn(), dsxigp, sprojalpha, cell->GetDerivAuxn(), lingp);

    // evalute the GP master coordinate derivatives
    std::vector<GEN::pairedvector<int, double>> dmxigp(2, (nrow + ncolP) * 100 * ndof);
    DerivXiGP3DAuxPlane(mele, mxi, cell->Auxn(), dmxigp, mprojalpha, cell->GetDerivAuxn(), lingp);

    // convert deriv sxi
    std::vector<GEN::pairedvector<int, double>> dlxigp(1, 100 * (nrow + ncolL) * ndof);

    if (lele.Id() == 0)
    {
      for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p) dlxigp[0][p->first] += (p->second);
    }
    else if (lele.Id() == 1)
    {
      for (_CI p = dmxigp[1].begin(); p != dmxigp[1].end(); ++p) dlxigp[0][p->first] += (p->second);
    }
    else if (lele.Id() == 2)
    {
      for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p) dlxigp[0][p->first] -= (p->second);
    }
    else if (lele.Id() == 3)
    {
      for (_CI p = dmxigp[1].begin(); p != dmxigp[1].end(); ++p) dlxigp[0][p->first] -= (p->second);
    }
    else
      dserror("ERROR: invalid local id of line ele");

    //**********************************************************************
    // geometric quantities
    //**********************************************************************
    double gpn[3] = {0.0, 0.0, 0.0};
    GEN::pairedvector<int, double> dgapgp(
        (ncolL * ndof) + 10 * linsize);  // gap lin. without lm and jac.
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

      sgpx[0] += sval[i] * sele.GetNodalCoords(0, i);
      sgpx[1] += sval[i] * sele.GetNodalCoords(1, i);
      sgpx[2] += sval[i] * sele.GetNodalCoords(2, i);
    }

    // build interpolation of master GP coordinates
    for (int i = 0; i < ncolL; ++i)
    {
      mgpx[0] += mval[i] * lele.GetNodalCoords(0, i);
      mgpx[1] += mval[i] * lele.GetNodalCoords(1, i);
      mgpx[2] += mval[i] * lele.GetNodalCoords(2, i);
    }

    // normalize interpolated GP normal back to length 1.0 !!!
    double lengthn = sqrt(gpn[0] * gpn[0] + gpn[1] * gpn[1] + gpn[2] * gpn[2]);
    if (lengthn < 1.0e-12) dserror("ERROR: IntegrateAndDerivSegment: Divide by zero!");

    for (int i = 0; i < 3; ++i) gpn[i] /= lengthn;

    // build gap function at current GP
    for (int i = 0; i < Dim(); ++i) gap[0] += (mgpx[i] - sgpx[i]) * gpn[i];

    // **************************
    // Linearization
    // **************************
    int linsize = 0;
    for (int i = 0; i < nrow; ++i)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(mynodes[i]);
      linsize += 10 * cnode->GetLinsize();
    }

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

      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
      {
        double valx = sderiv(i, 0) * cnode->MoData().n()[0];
        dmap_nxsl_gp[p->first] += valx * (p->second);
        double valy = sderiv(i, 0) * cnode->MoData().n()[1];
        dmap_nysl_gp[p->first] += valy * (p->second);
        double valz = sderiv(i, 0) * cnode->MoData().n()[2];
        dmap_nzsl_gp[p->first] += valz * (p->second);
      }

      for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
      {
        double valx = sderiv(i, 1) * cnode->MoData().n()[0];
        dmap_nxsl_gp[p->first] += valx * (p->second);
        double valy = sderiv(i, 1) * cnode->MoData().n()[1];
        dmap_nysl_gp[p->first] += valy * (p->second);
        double valz = sderiv(i, 1) * cnode->MoData().n()[2];
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

    for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < nrow; ++z)
      {
        CoNode* cnode = dynamic_cast<CoNode*>(mynodes[z]);
        for (int k = 0; k < 3; ++k) dg -= gpn[k] * sderiv(z, 0) * cnode->xspatial()[k] * ps;
      }
    }
    for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < nrow; ++z)
      {
        CoNode* cnode = dynamic_cast<CoNode*>(mynodes[z]);
        for (int k = 0; k < 3; ++k) dg -= gpn[k] * sderiv(z, 1) * cnode->xspatial()[k] * ps;
      }
    }

    //        MASTER
    // lin master nodes
    for (int z = 0; z < ncolL; ++z)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(mnodes[z]);
      for (int k = 0; k < 3; ++k) dgapgp[cnode->Dofs()[k]] += mval[z] * gpn[k];
    }

    for (_CI p = dlxigp[0].begin(); p != dlxigp[0].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < ncolL; ++z)
      {
        CoNode* cnode = dynamic_cast<CoNode*>(mnodes[z]);
        for (int k = 0; k < 3; ++k) dg += gpn[k] * mderiv(z, 0) * cnode->xspatial()[k] * ps;
      }
    }



    // weighted gap
    for (int j = 0; j < nrow; ++j)
    {
      CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(mynodes[j]);

      double prod = 0.0;
      // Petrov-Galerkin approach (dual LM for D/M but standard LM for gap)
      if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) prod = sval[j] * gap[0] * jac * wgt;
      // usual standard or dual LM approach
      else
        prod = lmval[j] * gap[0] * jac * wgt;

      // do not process slave side boundary nodes
      // (their row entries would be zero anyway!)
      if (cnode->IsOnBound()) continue;

      // add current Gauss point's contribution to gseg
      cnode->AddltsGapValue(prod);
    }

    for (int iter = 0; iter < nrow; ++iter)
    {
      // map iterator
      typedef GEN::pairedvector<int, double>::const_iterator _CI;

      MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[iter]);
      if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

      static double fac = 0.0;

      // get the corresponding map as a reference
      std::map<int, double>& dgmap =
          dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivGlts();

      // switch if Petrov-Galerkin approach for LM is applied
      if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
      {
        // (1) Lin(Phi) - does not exist here for Petrov-Galerkin approach

        // (2) Lin(N) - slave GP coordinates
        for (int d = 0; d < Dim() - 1; ++d)
        {
          fac = wgt * sderiv(iter, d) * gap[0] * jac;
          for (_CI p = dsxigp[d].begin(); p != dsxigp[d].end(); ++p)
            dgmap[p->first] += fac * (p->second);
        }

        // (3) Lin(g) - gap function
        fac = wgt * sval[iter] * jac;
        for (_CI p = dgapgp.begin(); p != dgapgp.end(); ++p) dgmap[p->first] += fac * (p->second);

        // (4) Lin(dsxideta) - intcell GP Jacobian
        fac = wgt * sval[iter] * gap[0];
        for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
          dgmap[p->first] += fac * (p->second);
      }

      // the usual standard or dual LM approach
      else
      {
        // (1) Lin(Phi) - dual shape functions
        // only line2 --> not required!

        // (2) Lin(Phi) - slave GP coordinates
        for (int d = 0; d < Dim() - 1; ++d)
        {
          fac = wgt * lmderiv(iter, d) * gap[0] * jac;
          for (_CI p = dsxigp[d].begin(); p != dsxigp[d].end(); ++p)
            dgmap[p->first] += fac * (p->second);
        }

        // (3) Lin(g) - gap function
        fac = wgt * lmval[iter] * jac;
        for (_CI p = dgapgp.begin(); p != dgapgp.end(); ++p) dgmap[p->first] += fac * (p->second);

        // (4) Lin(dsxideta) - intcell GP Jacobian
        fac = wgt * lmval[iter] * gap[0];
        for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
          dgmap[p->first] += fac * (p->second);
      }
    }

    // integrate D and M matrix
    // dual shape functions
    if (ShapeFcn() == INPAR::MORTAR::shape_dual ||
        ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
    {
      for (int j = 0; j < nrow; ++j)
      {
        CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(mynodes[j]);

        // integrate mseg
        for (int k = 0; k < ncolL; ++k)
        {
          CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mnodes[k]);

          // multiply the two shape functions
          double prod = lmval[j] * mval[k] * jac * wgt;

          if (abs(prod) > MORTARINTTOL) cnode->AddDltsValue(cnode->Id(), prod);
          if (abs(prod) > MORTARINTTOL) cnode->AddSNode(cnode->Id());  // only for friction!

          if (abs(prod) > MORTARINTTOL) cnode->AddMltsValue(mnode->Id(), prod);
          if (abs(prod) > MORTARINTTOL) cnode->AddMNode(mnode->Id());  // only for friction!
        }
      }
    }
    else
    {
      for (int j = 0; j < nrow; ++j)
      {
        CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(mynodes[j]);

        // integrate dseg
        for (int k = 0; k < nrow; ++k)
        {
          CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mynodes[k]);

          // multiply the two shape functions
          double prod = lmval[j] * sval[k] * jac * wgt;

          if (abs(prod) > MORTARINTTOL) cnode->AddDltsValue(mnode->Id(), prod);
          if (abs(prod) > MORTARINTTOL) cnode->AddSNode(mnode->Id());  // only for friction!
        }

        // integrate mseg
        for (int k = 0; k < ncolL; ++k)
        {
          CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mnodes[k]);

          // multiply the two shape functions
          double prod = lmval[j] * mval[k] * jac * wgt;

          if (abs(prod) > MORTARINTTOL) cnode->AddMltsValue(mnode->Id(), prod);
          if (abs(prod) > MORTARINTTOL) cnode->AddMNode(mnode->Id());  // only for friction!
        }
      }
    }

    if (ShapeFcn() == INPAR::MORTAR::shape_dual ||
        ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
    {
      // integrate LinD
      for (int j = 0; j < nrow; ++j)
      {
        MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[j]);
        if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

        int sgid = mymrtrnode->Id();
        std::map<int, double>& ddmap_jk =
            dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivDlts()[sgid];


        // integrate LinM
        for (int k = 0; k < ncolL; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          static double fac = 0.0;

          // get the correct map as a reference
          std::map<int, double>& dmmap_jk =
              dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivMlts()[mgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmderiv(j, 0) * mval[k] * jac;
          for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac * (p->second);
            ddmap_jk[p->first] += fac * (p->second);
          }

          fac = wgt * lmderiv(j, 1) * mval[k] * jac;
          for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
          {
            dmmap_jk[p->first] += fac * (p->second);
            ddmap_jk[p->first] += fac * (p->second);
          }

          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt * lmval[j] * mderiv(k, 0) * jac;
          for (_CI p = dlxigp[0].begin(); p != dlxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac * (p->second);
            ddmap_jk[p->first] += fac * (p->second);
          }

          //        fac = wgt*lmval[j]*mderiv(k, 1)*jac;
          //        for (_CI p=dmxigp[1].begin(); p!=dmxigp[1].end(); ++p)
          //        {
          //          dmmap_jk[p->first] += fac*(p->second);
          //          ddmap_jk[p->first] += fac*(p->second);
          //        }

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt * lmval[j] * mval[k];
          for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
          {
            dmmap_jk[p->first] += fac * (p->second);
            ddmap_jk[p->first] += fac * (p->second);
          }
        }  // loop over master nodes
      }
    }
    else  // standard case
    {
      // integrate LinD
      for (int j = 0; j < nrow; ++j)
      {
        MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[j]);
        if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");


        // integrate LinD
        for (int k = 0; k < nrow; ++k)
        {
          // global master node ID
          int mgid = sele.Nodes()[k]->Id();
          static double fac = 0.0;

          // get the correct map as a reference
          std::map<int, double>& ddmap_jk =
              dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivDlts()[mgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmderiv(j, 0) * sval[k] * jac;
          for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
          {
            ddmap_jk[p->first] += fac * (p->second);
          }

          fac = wgt * lmderiv(j, 1) * sval[k] * jac;
          for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
          {
            ddmap_jk[p->first] += fac * (p->second);
          }

          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt * lmval[j] * sderiv(k, 0) * jac;
          for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
          {
            ddmap_jk[p->first] += fac * (p->second);
          }

          fac = wgt * lmval[j] * sderiv(k, 1) * jac;
          for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
          {
            ddmap_jk[p->first] += fac * (p->second);
          }

          //        fac = wgt*lmval[j]*mderiv(k, 1)*jac;
          //        for (_CI p=dmxigp[1].begin(); p!=dmxigp[1].end(); ++p)
          //        {
          //          dmmap_jk[p->first] += fac*(p->second);
          //          ddmap_jk[p->first] += fac*(p->second);
          //        }

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt * lmval[j] * sval[k];
          for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
          {
            ddmap_jk[p->first] += fac * (p->second);
          }
        }  // loop over slave nodes

        // integrate LinM
        for (int k = 0; k < ncolL; ++k)
        {
          // global master node ID
          int mgid = lele.Nodes()[k]->Id();
          static double fac = 0.0;

          // get the correct map as a reference
          std::map<int, double>& dmmap_jk =
              dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivMlts()[mgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmderiv(j, 0) * mval[k] * jac;
          for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac * (p->second);
          }

          fac = wgt * lmderiv(j, 1) * mval[k] * jac;
          for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
          {
            dmmap_jk[p->first] += fac * (p->second);
          }

          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt * lmval[j] * mderiv(k, 0) * jac;
          for (_CI p = dlxigp[0].begin(); p != dlxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac * (p->second);
          }

          //        fac = wgt*lmval[j]*mderiv(k, 1)*jac;
          //        for (_CI p=dmxigp[1].begin(); p!=dmxigp[1].end(); ++p)
          //        {
          //          dmmap_jk[p->first] += fac*(p->second);
          //          ddmap_jk[p->first] += fac*(p->second);
          //        }

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt * lmval[j] * mval[k];
          for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
          {
            dmmap_jk[p->first] += fac * (p->second);
          }
        }  // loop over master nodes
      }
    }


  }  // end gp loop
  //**********************************************************************

  return;
}


/*----------------------------------------------------------------------*
 |  Integrate and linearize a 1D slave / master cell (2D)    farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateDerivCell3DAuxPlaneLTS(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& lsele, MORTAR::MortarElement& mele, Teuchos::RCP<MORTAR::IntCell> cell,
    double* auxn, const Epetra_Comm& comm)
{
  // explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without specific shape function defined!");

  // check for problem dimension
  if (Dim() != 3) dserror("ERROR: 3D integration method called for non-3D problem");

  // discretization type of master element
  DRT::Element::DiscretizationType sdt = sele.Shape();
  DRT::Element::DiscretizationType mdt = mele.Shape();

  // check input data
  //  if ((!sele.IsSlave()) || (mele.IsSlave()))
  //    dserror("ERROR: IntegrateDerivCell3DAuxPlane called on a wrong type of MortarElement
  //    pair!");
  if (cell == Teuchos::null)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without integration cell");

  // number of nodes (slave, master)
  int nrowL = lsele.NumNode();
  int nrowS = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = Dim();

  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = lsele.Nodes();
  if (!mynodes)
    dserror("ERROR: IntegrateDerivCell3DAuxPlaneLTS: Cannot access slave ndoes. Null pointer!");
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes)
    dserror("ERROR: IntegrateDerivCell3DAuxPlaneLTS: Cannot access master nodes. Null pointer!");

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrowL);
  LINALG::SerialDenseMatrix sderiv(nrowL, 1, true);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol, 2, true);
  LINALG::SerialDenseVector lmval(nrowL);
  LINALG::SerialDenseMatrix lmderiv(nrowL, 1, true);

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  GEN::pairedvector<int, Epetra_SerialDenseMatrix> dualmap(
      (nrowL + ncol) * ndof, 0, Epetra_SerialDenseMatrix(nrowL, nrowL));
  if ((ShapeFcn() == INPAR::MORTAR::shape_dual ||
          ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) &&
      (lsele.Shape() != MORTAR::MortarElement::line2 ||
          sele.MoData().DerivDualShape() != Teuchos::null))
    sele.DerivShapeDual(dualmap);

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  int linsize = 0;
  for (int i = 0; i < nrowL; ++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*>(mynodes[i]);
    linsize += cnode->GetLinsize();
  }
  linsize *= 100;  // TODO: remove this safety scaling

  // check if the cells are tri3
  // there's nothing wrong about other shapes, but as long as they are all
  // tri3 we can perform the jacobian calculation ( and its deriv) outside
  // the Gauss point loop
  if (cell->Shape() != DRT::Element::line2) dserror("only line2 integration cells for LTS");

  const double jac = cell->Jacobian();

  // directional derivative of cell Jacobian
  GEN::pairedvector<int, double> jacintcellmap((nrowL + ncol) * ndof + linsize);
  cell->DerivJacobian(jacintcellmap);

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < nGP(); ++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp, 0), Coordinate(gp, 1)};
    double wgt = Weight(gp);

    // get global Gauss point coordinates
    double globgp[3] = {0.0, 0.0, 0.0};
    cell->LocalToGlobal(eta, globgp, 0);

    // coordinate in parameter space on surface slave ele
    double sxi[2] = {0.0, 0.0};

    // coordinate in parameter space on surface master ele
    double mxi[2] = {0.0, 0.0};

    // coordinate in parameter space on line (edge) slave ele
    double lxi[2] = {0.0, 0.0};

    // project Gauss point onto slave/master surface element
    double sprojalpha = 0.0;
    double mprojalpha = 0.0;
    bool success = MORTAR::MortarProjector::Impl(sele)->ProjectGaussPointAuxn3D(
        globgp, auxn, sele, sxi, sprojalpha);
    if (!success)
    {
      //      dserror("ERROR: projection not possible!");
      return;
    }

    success = MORTAR::MortarProjector::Impl(mele)->ProjectGaussPointAuxn3D(
        globgp, auxn, mele, mxi, mprojalpha);
    if (!success)
    {
      //      dserror("ERROR: projection not possible!");
      return;
    }

    // check
    double checkpos[3] = {0.0, 0.0, 0.0};
    mele.LocalToGlobal(mxi, checkpos, 0);

    // check GP projection (SLAVE)
    const double sminedge = sele.MinEdgeSize();
    const double mminedge = mele.MinEdgeSize();
    const double tol = MORTARCLIPTOL * std::min(sminedge, mminedge);
    if (sdt == DRT::Element::quad4 || sdt == DRT::Element::quad8 || sdt == DRT::Element::quad9)
    {
      if (sxi[0] < -1.0 - tol || sxi[1] < -1.0 - tol || sxi[0] > 1.0 + tol || sxi[1] > 1.0 + tol)
      {
        std::cout << "\n***Warning: LTS: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }
    else
    {
      if (sxi[0] < -tol || sxi[1] < -tol || sxi[0] > 1.0 + tol || sxi[1] > 1.0 + tol ||
          sxi[0] + sxi[1] > 1.0 + 2 * tol)
      {
        std::cout << "\n***Warning: LTS: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }

    // check GP projection (MASTER)
    if (mdt == DRT::Element::quad4 || mdt == DRT::Element::quad8 || mdt == DRT::Element::quad9)
    {
      if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol)
      {
        std::cout << "\n***Warning: LTS: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
        std::cout << "globgp=   " << globgp[0] << "  " << globgp[1] << "  " << globgp[2]
                  << std::endl;
        std::cout << "checkpos= " << checkpos[0] << "  " << checkpos[1] << "  " << checkpos[2]
                  << std::endl;
      }
    }
    else
    {
      if (mxi[0] < -tol || mxi[1] < -tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol ||
          mxi[0] + mxi[1] > 1.0 + 2 * tol)
      {
        std::cout << "\n***Warning: LTS: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
        std::cout << "globgp=   " << globgp[0] << "  " << globgp[1] << "  " << globgp[2]
                  << std::endl;
        std::cout << "checkpos= " << checkpos[0] << "  " << checkpos[1] << "  " << checkpos[2]
                  << std::endl;
      }
    }

    // convert sxi to lsxi!!!
    if (lsele.Id() == 0)
    {
      lxi[0] = sxi[0];
    }
    else if (lsele.Id() == 1)
    {
      lxi[0] = sxi[1];
    }
    else if (lsele.Id() == 2)
    {
      lxi[0] = -sxi[0];
    }
    else if (lsele.Id() == 3)
    {
      lxi[0] = -sxi[1];
    }
    else
      dserror("ERROR: invalid local id of line ele");

    // evaluate Lagrange multiplier shape functions (on slave element)
    lsele.EvaluateShapeLagMult(ShapeFcn(), lxi, lmval, lmderiv, nrowL);

    // evaluate trace space shape functions (on both elements)
    lsele.EvaluateShape(lxi, sval, sderiv, nrowL);
    mele.EvaluateShape(mxi, mval, mderiv, ncol);

    // evaluate linearizations *******************************************

    // evaluate global GP coordinate derivative
    static LINALG::Matrix<3, 1> svalcell;
    static LINALG::Matrix<3, 2> sderivcell;
    cell->EvaluateShape(eta, svalcell, sderivcell);

    GEN::pairedvector<int, LINALG::Matrix<3, 1>> lingp((nrowS + ncol) * ndof + linsize);

    for (int v = 0; v < 2; ++v)
      for (int d = 0; d < 3; ++d)
        for (_CI p = (cell->GetDerivVertex(v))[d].begin(); p != (cell->GetDerivVertex(v))[d].end();
             ++p)
          lingp[p->first](d) += svalcell(v) * (p->second);

    // evalute the GP slave coordinate derivatives
    std::vector<GEN::pairedvector<int, double>> dsxigp(2, (nrowS + ncol) * ndof + linsize);
    DerivXiGP3DAuxPlane(sele, sxi, cell->Auxn(), dsxigp, sprojalpha, cell->GetDerivAuxn(), lingp);

    // evalute the GP master coordinate derivatives
    std::vector<GEN::pairedvector<int, double>> dmxigp(2, (nrowS + ncol) * ndof + linsize);
    DerivXiGP3DAuxPlane(mele, mxi, cell->Auxn(), dmxigp, mprojalpha, cell->GetDerivAuxn(), lingp);

    // convert deriv sxi
    std::vector<GEN::pairedvector<int, double>> dlxigp(1, (nrowS + ncol) * ndof + linsize);

    if (lsele.Id() == 0)
    {
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p) dlxigp[0][p->first] += (p->second);
    }
    else if (lsele.Id() == 1)
    {
      for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p) dlxigp[0][p->first] += (p->second);
    }
    else if (lsele.Id() == 2)
    {
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p) dlxigp[0][p->first] -= (p->second);
    }
    else if (lsele.Id() == 3)
    {
      for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p) dlxigp[0][p->first] -= (p->second);
    }
    else
      dserror("ERROR: invalid local id of line ele");

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

    for (int i = 0; i < nrowL; ++i)
    {
      MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[i]);
      gpn[0] += sval[i] * mymrtrnode->MoData().n()[0];
      gpn[1] += sval[i] * mymrtrnode->MoData().n()[1];
      gpn[2] += sval[i] * mymrtrnode->MoData().n()[2];

      sgpx[0] += sval[i] * lsele.GetNodalCoords(0, i);
      sgpx[1] += sval[i] * lsele.GetNodalCoords(1, i);
      sgpx[2] += sval[i] * lsele.GetNodalCoords(2, i);
    }

    // build interpolation of master GP coordinates
    for (int i = 0; i < ncol; ++i)
    {
      mgpx[0] += mval[i] * mele.GetNodalCoords(0, i);
      mgpx[1] += mval[i] * mele.GetNodalCoords(1, i);
      mgpx[2] += mval[i] * mele.GetNodalCoords(2, i);
    }

    // normalize interpolated GP normal back to length 1.0 !!!
    double lengthn = sqrt(gpn[0] * gpn[0] + gpn[1] * gpn[1] + gpn[2] * gpn[2]);
    if (lengthn < 1.0e-12) dserror("ERROR: IntegrateAndDerivSegment: Divide by zero!");

    for (int i = 0; i < 3; ++i) gpn[i] /= lengthn;

    // build gap function at current GP
    for (int i = 0; i < Dim(); ++i) gap[0] += (mgpx[i] - sgpx[i]) * gpn[i];

    // **************************
    // Linearization
    // **************************
    int linsize = 0;
    for (int i = 0; i < nrowL; ++i)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(mynodes[i]);
      linsize += 10 * cnode->GetLinsize();
    }

    // build directional derivative of slave GP normal (non-unit)
    GEN::pairedvector<int, double> dmap_nxsl_gp(linsize);
    GEN::pairedvector<int, double> dmap_nysl_gp(linsize);
    GEN::pairedvector<int, double> dmap_nzsl_gp(linsize);

    for (int i = 0; i < nrowL; ++i)
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

      for (_CI p = dlxigp[0].begin(); p != dlxigp[0].end(); ++p)
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
    for (int z = 0; z < nrowL; ++z)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(mynodes[z]);
      for (int k = 0; k < 3; ++k) dgapgp[cnode->Dofs()[k]] -= sval[z] * gpn[k];
    }

    for (_CI p = dlxigp[0].begin(); p != dlxigp[0].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < nrowL; ++z)
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

    for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < ncol; ++z)
      {
        CoNode* cnode = dynamic_cast<CoNode*>(mnodes[z]);
        for (int k = 0; k < 3; ++k) dg += gpn[k] * mderiv(z, 0) * cnode->xspatial()[k] * ps;
      }
    }

    for (_CI p = dmxigp[1].begin(); p != dmxigp[1].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < ncol; ++z)
      {
        CoNode* cnode = dynamic_cast<CoNode*>(mnodes[z]);
        for (int k = 0; k < 3; ++k) dg += gpn[k] * mderiv(z, 1) * cnode->xspatial()[k] * ps;
      }
    }

    // weighted gap
    for (int j = 0; j < nrowL; ++j)
    {
      CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(mynodes[j]);

      double prod = 0.0;
      // Petrov-Galerkin approach (dual LM for D/M but standard LM for gap)
      if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) prod = sval[j] * gap[0] * jac * wgt;
      // usual standard or dual LM approach
      else
        prod = lmval[j] * gap[0] * jac * wgt;

      // do not process slave side boundary nodes
      // (their row entries would be zero anyway!)
      if (cnode->IsOnCorner()) continue;

      // add current Gauss point's contribution to gseg
      cnode->AddltsGapValue(prod);
    }

    for (int iter = 0; iter < nrowL; ++iter)
    {
      // map iterator
      typedef GEN::pairedvector<int, double>::const_iterator _CI;

      MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[iter]);
      if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

      if (mymrtrnode->IsOnCorner()) continue;

      static double fac = 0.0;

      // get the corresponding map as a reference
      std::map<int, double>& dgmap =
          dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivGlts();

      // switch if Petrov-Galerkin approach for LM is applied
      if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
      {
        // (1) Lin(Phi) - does not exist here for Petrov-Galerkin approach

        // (2) Lin(N) - slave GP coordinates
        for (int d = 0; d < Dim() - 2; ++d)
        {
          fac = wgt * sderiv(iter, d) * gap[0] * jac;
          for (_CI p = dlxigp[d].begin(); p != dlxigp[d].end(); ++p)
            dgmap[p->first] += fac * (p->second);
        }

        // (3) Lin(g) - gap function
        fac = wgt * sval[iter] * jac;
        for (_CI p = dgapgp.begin(); p != dgapgp.end(); ++p) dgmap[p->first] += fac * (p->second);

        // (4) Lin(dsxideta) - intcell GP Jacobian
        fac = wgt * sval[iter] * gap[0];
        for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
          dgmap[p->first] += fac * (p->second);
      }

      // the usual standard or dual LM approach
      else
      {
        // (1) Lin(Phi) - dual shape functions
        // only line2 --> not required!

        // (2) Lin(Phi) - slave GP coordinates
        for (int d = 0; d < Dim() - 2; ++d)
        {
          fac = wgt * lmderiv(iter, d) * gap[0] * jac;
          for (_CI p = dlxigp[d].begin(); p != dlxigp[d].end(); ++p)
            dgmap[p->first] += fac * (p->second);
        }

        // (3) Lin(g) - gap function
        fac = wgt * lmval[iter] * jac;
        for (_CI p = dgapgp.begin(); p != dgapgp.end(); ++p) dgmap[p->first] += fac * (p->second);

        // (4) Lin(dsxideta) - intcell GP Jacobian
        fac = wgt * lmval[iter] * gap[0];
        for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
          dgmap[p->first] += fac * (p->second);
      }
    }
    //--------------------------------------------------------
    // integrate D and M matrix
    // dual shape functions
    if (ShapeFcn() == INPAR::MORTAR::shape_dual ||
        ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
    {
      bool bound = false;
      for (int j = 0; j < nrowL; ++j)
      {
        CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(mynodes[j]);
        if (cnode->IsOnCorner()) bound = true;
      }

      for (int j = 0; j < nrowL; ++j)
      {
        CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(mynodes[j]);
        if (cnode->IsOnCorner()) continue;

        // integrate mseg
        for (int k = 0; k < ncol; ++k)
        {
          CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mnodes[k]);

          // multiply the two shape functions
          double prod = lmval[j] * mval[k] * jac * wgt;

          if (!bound and abs(prod) > MORTARINTTOL)
          {
            if (abs(prod) > MORTARINTTOL) cnode->AddDltsValue(cnode->Id(), prod);
            if (abs(prod) > MORTARINTTOL) cnode->AddSNode(cnode->Id());  // only for friction!
          }

          if (abs(prod) > MORTARINTTOL) cnode->AddMltsValue(mnode->Id(), prod);
          if (abs(prod) > MORTARINTTOL) cnode->AddMNode(mnode->Id());  // only for friction!
        }

        // integrate dseg (boundary modification)
        if (bound)
        {
          //        bool j_boundnode = cnode->IsOnBoundorCE();

          for (int k = 0; k < nrowL; ++k)
          {
            CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mynodes[k]);
            //          bool k_boundnode = mnode->IsOnBoundorCE();

            // do not assemble off-diagonal terms if j,k are both non-boundary nodes
            //          if (!j_boundnode && !k_boundnode && (j!=k))
            //            continue;

            // multiply the two shape functions
            double prod = lmval[j] * sval[k] * jac * wgt;

            // isolate the dseg entries to be filled
            // (both the main diagonal and every other secondary diagonal)
            // and add current Gauss point's contribution to dseg
            if (mnode->IsOnCorner())
            {
              if (abs(prod) > MORTARINTTOL) cnode->AddDltsValue(mnode->Id(), prod);
              if (abs(prod) > MORTARINTTOL) cnode->AddSNode(mnode->Id());  // only for friction!
            }
            else
            {
              if (abs(prod) > MORTARINTTOL) cnode->AddDltsValue(cnode->Id(), prod);
              if (abs(prod) > MORTARINTTOL) cnode->AddSNode(cnode->Id());  // only for friction!
            }
          }
        }
      }

      // OLD DM DUAL IMPLEMENTATION
      //      for (int j=0;j<nrowL;++j)
      //      {
      //        CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(mynodes[j]);
      //        if (cnode->IsOnCorner())
      //          continue;
      //
      //        // integrate mseg
      //        for (int k=0; k<ncol; ++k)
      //        {
      //          CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mnodes[k]);
      //
      //          // multiply the two shape functions
      //          double prod = lmval[j]*mval[k]*jac*wgt;
      //
      //          if(abs(prod)>MORTARINTTOL) cnode->AddDltsValue(cnode->Id(),prod);
      //          if(abs(prod)>MORTARINTTOL) cnode->AddSNode(cnode->Id()); // only for friction!
      //
      //          if(abs(prod)>MORTARINTTOL) cnode->AddMltsValue(mnode->Id(),prod);
      //          if(abs(prod)>MORTARINTTOL) cnode->AddMNode(mnode->Id());  // only for friction!
      //        }
      //      }
    }
    // STANDARD SHAPE FUNCTIONS
    else  // std
    {
      for (int j = 0; j < nrowL; ++j)
      {
        CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(mynodes[j]);

        if (cnode->IsOnCorner()) continue;

        // integrate mseg
        for (int k = 0; k < ncol; ++k)
        {
          CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mnodes[k]);

          // multiply the two shape functions
          double prod = lmval[j] * mval[k] * jac * wgt;

          if (abs(prod) > MORTARINTTOL) cnode->AddMltsValue(mnode->Id(), prod);
          if (abs(prod) > MORTARINTTOL) cnode->AddMNode(mnode->Id());  // only for friction!
        }

        // integrate dseg
        for (int k = 0; k < nrowL; ++k)
        {
          CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(mynodes[k]);

          if (snode->IsOnCorner())
          {
            // multiply the two shape functions
            double prod = lmval[j] * sval[k] * jac * wgt;

            if (abs(prod) > MORTARINTTOL) cnode->AddDltsValue(snode->Id(), prod);
            if (abs(prod) > MORTARINTTOL) cnode->AddSNode(snode->Id());  // only for friction!
          }
          else
          {
            // multiply the two shape functions
            double prod = lmval[j] * sval[k] * jac * wgt;

            if (abs(prod) > MORTARINTTOL) cnode->AddDltsValue(snode->Id(), prod);
            if (abs(prod) > MORTARINTTOL) cnode->AddSNode(snode->Id());  // only for friction!
          }
        }
      }
    }

    // dual shape functions
    if (ShapeFcn() == INPAR::MORTAR::shape_dual ||
        ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
    {
      bool bound = false;
      for (int j = 0; j < nrowL; ++j)
      {
        CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(mynodes[j]);
        if (cnode->IsOnCorner()) bound = true;
      }

      if (bound)
      {
        // integrate LinD
        for (int j = 0; j < nrowL; ++j)
        {
          MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[j]);
          if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

          // node j is boundary node
          if (mymrtrnode->IsOnCorner()) continue;

          // integrate LinM
          for (int k = 0; k < ncol; ++k)
          {
            // global master node ID
            int mgid = mele.Nodes()[k]->Id();
            static double fac = 0.0;

            // get the correct map as a reference
            std::map<int, double>& dmmap_jk =
                dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivMlts()[mgid];

            // (1) Lin(Phi) - dual shape functions
            //          if (dualmap.size()>0)
            //          {
            //            for (int m=0;m<nrowL;++m)
            //            {
            //              fac = wgt*sval[m]*mval[k]*jac;
            //              for (GEN::pairedvector<int,Epetra_SerialDenseMatrix>::const_iterator
            //              p=dualmap.begin();
            //                  p!=dualmap.end();++p)
            //                dmmap_jk[p->first] += fac*(p->second)(j,m);
            //            }
            //          }

            // (2) Lin(NSlave) - slave GP coordinates
            fac = wgt * lmderiv(j, 0) * mval[k] * jac;
            for (_CI p = dlxigp[0].begin(); p != dlxigp[0].end(); ++p)
              dmmap_jk[p->first] += fac * (p->second);

            //            fac = wgt*lmderiv(j, 1)*mval[k]*jac;
            //            for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
            //              dmmap_jk[p->first] += fac*(p->second);

            // (3) Lin(NMaster) - master GP coordinates
            fac = wgt * lmval[j] * mderiv(k, 0) * jac;
            for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
              dmmap_jk[p->first] += fac * (p->second);

            fac = wgt * lmval[j] * mderiv(k, 1) * jac;
            for (_CI p = dmxigp[1].begin(); p != dmxigp[1].end(); ++p)
              dmmap_jk[p->first] += fac * (p->second);

            // (4) Lin(dsxideta) - intcell GP Jacobian
            fac = wgt * lmval[j] * mval[k];
            for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
              dmmap_jk[p->first] += fac * (p->second);
          }  // loop over master nodes

          // loop over slave nodes
          for (int k = 0; k < nrowL; ++k)
          {
            MORTAR::MortarNode* mymrtrnode2 = dynamic_cast<MORTAR::MortarNode*>(mynodes[k]);
            if (!mymrtrnode2) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

            // global master node ID
            int sgid = mymrtrnode2->Id();
            static double fac = 0.0;

            // node k is boundary node
            if (mymrtrnode2->IsOnCorner())
            {
              // move entry to derivM (with minus sign)
              // get the correct map as a reference
              std::map<int, double>& dmmap_jk =
                  dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivDlts()[sgid];

              // (1) Lin(Phi) - dual shape functions
              //            if (dualmap.size()>0)
              //            {
              //              for (int m=0;m<nrow;++m)
              //              {
              //                fac = wgt*sval[m]*sval[k]*jac;
              //                for (GEN::pairedvector<int,Epetra_SerialDenseMatrix>::const_iterator
              //                p=dualmap.begin();
              //                    p!=dualmap.end();++p)
              //                  dmmap_jk[p->first] -= fac*(p->second)(j,m);
              //              }
              //            }

              // (2) Lin(NSlave) - slave GP coordinates
              fac = wgt * lmderiv(j, 0) * sval[k] * jac;
              for (_CI p = dlxigp[0].begin(); p != dlxigp[0].end(); ++p)
                dmmap_jk[p->first] += fac * (p->second);

              //              fac = wgt*lmderiv(j, 1)*sval[k]*jac;
              //              for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
              //                dmmap_jk[p->first] -= fac*(p->second);

              // (3) Lin(NSlave) - slave GP coordinates
              fac = wgt * lmval[j] * sderiv(k, 0) * jac;
              for (_CI p = dlxigp[0].begin(); p != dlxigp[0].end(); ++p)
                dmmap_jk[p->first] += fac * (p->second);

              //              fac = wgt*lmval[j]*sderiv(k, 1)*jac;
              //              for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
              //                dmmap_jk[p->first] -= fac*(p->second);

              // (4) Lin(dsxideta) - intcell GP Jacobian
              fac = wgt * lmval[j] * sval[k];
              for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
                dmmap_jk[p->first] += fac * (p->second);
            }

            // node k is NO boundary node
            else
            {
              // get the correct map as a reference
              std::map<int, double>& ddmap_jk = dynamic_cast<CONTACT::CoNode*>(mymrtrnode)
                                                    ->CoData()
                                                    .GetDerivDlts()[mymrtrnode->Id()];

              // (1) Lin(Phi) - dual shape functions
              //            if (dualmap.size()>0)
              //            {
              //              for (int m=0;m<nrow;++m)
              //              {
              //                fac = wgt*sval[m]*sval[k]*jac;
              //                for (GEN::pairedvector<int,Epetra_SerialDenseMatrix>::const_iterator
              //                p=dualmap.begin();
              //                    p!=dualmap.end();++p)
              //                  ddmap_jk[p->first] += fac*(p->second)(j,m);
              //              }
              //            }

              // (2) Lin(NSlave) - slave GP coordinates
              fac = wgt * lmderiv(j, 0) * sval[k] * jac;
              for (_CI p = dlxigp[0].begin(); p != dlxigp[0].end(); ++p)
                ddmap_jk[p->first] += fac * (p->second);

              //              fac = wgt*lmderiv(j, 1)*sval[k]*jac;
              //              for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
              //                ddmap_jk[p->first] += fac*(p->second);

              // (3) Lin(NSlave) - slave GP coordinates
              fac = wgt * lmval[j] * sderiv(k, 0) * jac;
              for (_CI p = dlxigp[0].begin(); p != dlxigp[0].end(); ++p)
                ddmap_jk[p->first] += fac * (p->second);

              //              fac = wgt*lmval[j]*sderiv(k, 1)*jac;
              //              for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
              //                ddmap_jk[p->first] += fac*(p->second);

              // (4) Lin(dsxideta) - intcell GP Jacobian
              fac = wgt * lmval[j] * sval[k];
              for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
                ddmap_jk[p->first] += fac * (p->second);
            }
          }  // loop over slave nodes
        }
      }
      else
      {
        // OLD DM DUAL IMPLEMENTATION
        // integrate LinD
        for (int j = 0; j < nrowL; ++j)
        {
          MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[j]);
          if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

          int sgid = mymrtrnode->Id();
          std::map<int, double>& ddmap_jk =
              dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivDlts()[sgid];

          // integrate LinM
          for (int k = 0; k < ncol; ++k)
          {
            // global master node ID
            int mgid = mele.Nodes()[k]->Id();
            static double fac = 0.0;

            // get the correct map as a reference
            std::map<int, double>& dmmap_jk =
                dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivMlts()[mgid];

            // (1) Lin(Phi) - dual shape functions
            // this vanishes here since there are no deformation-dependent dual functions

            // (2) Lin(NSlave) - slave GP coordinates
            fac = wgt * lmderiv(j, 0) * mval[k] * jac;
            for (_CI p = dlxigp[0].begin(); p != dlxigp[0].end(); ++p)
            {
              dmmap_jk[p->first] += fac * (p->second);
              ddmap_jk[p->first] += fac * (p->second);
            }

            // (3) Lin(NMaster) - master GP coordinates
            fac = wgt * lmval[j] * mderiv(k, 0) * jac;
            for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
            {
              dmmap_jk[p->first] += fac * (p->second);
              ddmap_jk[p->first] += fac * (p->second);
            }

            fac = wgt * lmval[j] * mderiv(k, 1) * jac;
            for (_CI p = dmxigp[1].begin(); p != dmxigp[1].end(); ++p)
            {
              dmmap_jk[p->first] += fac * (p->second);
              ddmap_jk[p->first] += fac * (p->second);
            }

            // (4) Lin(dsxideta) - intcell GP Jacobian
            fac = wgt * lmval[j] * mval[k];
            for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
            {
              dmmap_jk[p->first] += fac * (p->second);
              ddmap_jk[p->first] += fac * (p->second);
            }
          }  // loop over master nodes
        }
      }
    }
    // STANDARD SHAPE FUNCTIONS
    else  // std
    {
      // integrate LinD
      for (int j = 0; j < nrowL; ++j)
      {
        MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[j]);
        if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

        if (mymrtrnode->IsOnCorner()) continue;

        // integrate LinM
        for (int k = 0; k < ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          static double fac = 0.0;

          // get the correct map as a reference
          std::map<int, double>& dmmap_jk =
              dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivMlts()[mgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmderiv(j, 0) * mval[k] * jac;
          for (_CI p = dlxigp[0].begin(); p != dlxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac * (p->second);
          }

          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt * lmval[j] * mderiv(k, 0) * jac;
          for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac * (p->second);
          }

          fac = wgt * lmval[j] * mderiv(k, 1) * jac;
          for (_CI p = dmxigp[1].begin(); p != dmxigp[1].end(); ++p)
          {
            dmmap_jk[p->first] += fac * (p->second);
          }

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt * lmval[j] * mval[k];
          for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
          {
            dmmap_jk[p->first] += fac * (p->second);
          }
        }  // loop over master nodes

        // integrate LinD
        for (int k = 0; k < nrowL; ++k)
        {
          // global master node ID
          int mgid = mynodes[k]->Id();
          static double fac = 0.0;

          if (dynamic_cast<CoNode*>(mynodes[k])->IsOnCorner())
          {
            // get the correct map as a reference
            std::map<int, double>& dmmap_jk =
                dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivDlts()[mgid];

            // (1) Lin(Phi) - dual shape functions
            // this vanishes here since there are no deformation-dependent dual functions

            // (2) Lin(NSlave) - slave GP coordinates
            fac = wgt * lmderiv(j, 0) * sval[k] * jac;
            for (_CI p = dlxigp[0].begin(); p != dlxigp[0].end(); ++p)
            {
              dmmap_jk[p->first] += fac * (p->second);
            }

            // (3) Lin(NMaster) - master GP coordinates
            fac = wgt * lmval[j] * sderiv(k, 0) * jac;
            for (_CI p = dlxigp[0].begin(); p != dlxigp[0].end(); ++p)
            {
              dmmap_jk[p->first] += fac * (p->second);
            }

            // (4) Lin(dsxideta) - intcell GP Jacobian
            fac = wgt * lmval[j] * sval[k];
            for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
            {
              dmmap_jk[p->first] += fac * (p->second);
            }
          }
          else
          {
            // get the correct map as a reference
            std::map<int, double>& dmmap_jk =
                dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivDlts()[mgid];

            // (1) Lin(Phi) - dual shape functions
            // this vanishes here since there are no deformation-dependent dual functions

            // (2) Lin(NSlave) - slave GP coordinates
            fac = wgt * lmderiv(j, 0) * sval[k] * jac;
            for (_CI p = dlxigp[0].begin(); p != dlxigp[0].end(); ++p)
            {
              dmmap_jk[p->first] += fac * (p->second);
            }

            // (3) Lin(NMaster) - master GP coordinates
            fac = wgt * lmval[j] * sderiv(k, 0) * jac;
            for (_CI p = dlxigp[0].begin(); p != dlxigp[0].end(); ++p)
            {
              dmmap_jk[p->first] += fac * (p->second);
            }

            // (4) Lin(dsxideta) - intcell GP Jacobian
            fac = wgt * lmval[j] * sval[k];
            for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
            {
              dmmap_jk[p->first] += fac * (p->second);
            }
          }
        }  // loop over master nodes
      }
    }


  }  // end gp loop
  //**********************************************************************

  return;
}


/*----------------------------------------------------------------------*
 |  Integrate and linearize a 2D slave / master cell (3D)     popp 03/09|
 |  This method integrates the cell M matrix and weighted gap g~        |
 |  and stores it in mseg and gseg respectively. Moreover, derivatives  |
 |  LinM and Ling are built and stored directly into the adjacent nodes.|
 |  (Thus this method combines EVERYTHING before done separately in     |
 |  IntegrateM3D, IntegrateG3D, DerivM3D and DerivG3D!)                 |
 |  This is the QUADRATIC auxiliary plane coupling version!!!           |
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateDerivCell3DAuxPlaneQuad(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, MORTAR::IntElement& sintele, MORTAR::IntElement& mintele,
    Teuchos::RCP<MORTAR::IntCell> cell, double* auxn)
{
  // get LMtype
  INPAR::MORTAR::LagMultQuad lmtype = LagMultQuad();

  // explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror(
        "ERROR: IntegrateDerivCell3DAuxPlaneQuad called without specific shape function defined!");

  // Petrov-Galerkin approach for LM not yet implemented
  if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin && sele.Shape() != DRT::Element::nurbs9)
    dserror(
        "ERROR: Petrov-Galerkin approach not yet implemented for 3-D quadratic FE interpolation");

  // check for problem dimension
  if (Dim() != 3) dserror("ERROR: 3D integration method called for non-3D problem");

  // discretization type of slave and master IntElement
  DRT::Element::DiscretizationType sdt = sintele.Shape();
  DRT::Element::DiscretizationType mdt = mintele.Shape();

  // discretization type of slave and master Element
  DRT::Element::DiscretizationType psdt = sele.Shape();
  DRT::Element::DiscretizationType pmdt = mele.Shape();

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror(
        "ERROR: IntegrateDerivCell3DAuxPlaneQuad called on a wrong type of MortarElement pair!");
  if (cell == Teuchos::null)
    dserror("ERROR: IntegrateDerivCell3DAuxPlaneQuad called without integration cell");

  // contact with wear
  bool wear = false;
  if (wearlaw_ != INPAR::WEAR::wear_none) wear = true;

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = Dim();
  int nintrow = sintele.NumNode();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow, 2, true);
  LINALG::SerialDenseVector svalmod(nrow);
  LINALG::SerialDenseMatrix sderivmod(nrow, 2, true);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol, 2, true);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow, 2, true);
  LINALG::SerialDenseVector lmintval(nintrow);
  LINALG::SerialDenseMatrix lmintderiv(nintrow, 2, true);

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseMatrix ssecderiv(nrow, 3);

  // get slave and master nodal coords for Jacobian / GP evaluation
  LINALG::SerialDenseMatrix scoord(3, sele.NumNode());
  sele.GetNodalCoords(scoord);
  LINALG::SerialDenseMatrix mcoord(3, mele.NumNode());
  mele.GetNodalCoords(mcoord);

  // nodal coords from previous time step and lagrange mulitplier
  Teuchos::RCP<LINALG::SerialDenseMatrix> scoordold;
  Teuchos::RCP<LINALG::SerialDenseMatrix> mcoordold;
  Teuchos::RCP<LINALG::SerialDenseMatrix> lagmult;

  // get them in the case of wear
  if (wear or DRT::INPUT::IntegralValue<int>(imortar_, "GP_SLIP_INCR") == true)
  {
    scoordold = Teuchos::rcp(new LINALG::SerialDenseMatrix(3, sele.NumNode()));
    mcoordold = Teuchos::rcp(new LINALG::SerialDenseMatrix(3, mele.NumNode()));
    lagmult = Teuchos::rcp(new LINALG::SerialDenseMatrix(3, sele.NumNode()));
    sele.GetNodalCoordsOld(*scoordold);
    mele.GetNodalCoordsOld(*mcoordold);
    sele.GetNodalLagMult(*lagmult);
  }

  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = sele.Nodes();
  if (!mynodes) dserror("ERROR: IntegrateDerivCell3DAuxPlaneQuad: Null pointer!");
  DRT::Node** myintnodes = sintele.Nodes();
  if (!myintnodes) dserror("ERROR: IntegrateDerivCell3DAuxPlaneQuad: Null pointer!");

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  bool duallin = false;
  GEN::pairedvector<int, Epetra_SerialDenseMatrix> dualmap(
      (nrow + ncol) * ndof, 0, Epetra_SerialDenseMatrix(nrow, nrow));
  if ((ShapeFcn() == INPAR::MORTAR::shape_dual ||
          ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) &&
      (sele.Shape() != MORTAR::MortarElement::tri3 ||
          sele.MoData().DerivDualShape() != Teuchos::null) &&
      LagMultQuad() != INPAR::MORTAR::lagmult_const)
  {
    duallin = true;
    sele.DerivShapeDual(dualmap);
  }

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k = 0; k < nrow; ++k)
  {
    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlaneQuad: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  // decide whether displacement shape fct. modification has to be considered or not
  // this is the case for dual quadratic Lagrange multipliers on quad8 and tri6 elements
  // check whether quadratic element with linear interpolation is present
  bool dualquad3d = false;
  bool linlm = false;
  if ((ShapeFcn() == INPAR::MORTAR::shape_dual ||
          ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) &&
      (lmtype == INPAR::MORTAR::lagmult_quad || lmtype == INPAR::MORTAR::lagmult_lin) &&
      (sele.Shape() == DRT::Element::quad8 || sele.Shape() == DRT::Element::tri6))
  {
    if (lmtype == INPAR::MORTAR::lagmult_quad)
      dualquad3d = true;
    else
      linlm = true;
  }

  // get linearization size
  int linsize = 0;
  for (int i = 0; i < nrow; ++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*>(mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  // check if the cells are tri3
  // there's nothing wrong about other shapes, but as long as they are all
  // tri3 we can perform the jacobian calculation ( and its deriv) outside
  // the Gauss point loop
  if (cell->Shape() != DRT::Element::tri3)
    dserror("only tri3 integration cells at the moment. See comment in the code");

  double jac = cell->Jacobian();
  // directional derivative of cell Jacobian
  GEN::pairedvector<int, double> jacintcellmap((nrow + ncol) * ndof);
  cell->DerivJacobian(jacintcellmap);

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < nGP(); ++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp, 0), Coordinate(gp, 1)};
    double wgt = Weight(gp);

    // get global Gauss point coordinates
    double globgp[3] = {0.0, 0.0, 0.0};
    cell->LocalToGlobal(eta, globgp, 0);

    double sxi[2] = {0.0, 0.0};
    double mxi[2] = {0.0, 0.0};

    // project Gauss point onto slave integration element
    // project Gauss point onto master integration element
    double sprojalpha = 0.0;
    double mprojalpha = 0.0;
    MORTAR::MortarProjector::Impl(sintele)->ProjectGaussPointAuxn3D(
        globgp, auxn, sintele, sxi, sprojalpha);
    MORTAR::MortarProjector::Impl(mintele)->ProjectGaussPointAuxn3D(
        globgp, auxn, mintele, mxi, mprojalpha);

    // check GP projection (SLAVE)
    const double tol = 0.01;
    if (sdt == DRT::Element::quad4 || sdt == DRT::Element::quad8 || sdt == DRT::Element::quad9)
    {
      if (sxi[0] < -1.0 - tol || sxi[1] < -1.0 - tol || sxi[0] > 1.0 + tol || sxi[1] > 1.0 + tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlaneQuad: Slave Gauss point projection "
                     "outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP (IntElement) projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }
    else
    {
      if (sxi[0] < -tol || sxi[1] < -tol || sxi[0] > 1.0 + tol || sxi[1] > 1.0 + tol ||
          sxi[0] + sxi[1] > 1.0 + 2 * tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlaneQuad: Slave Gauss point projection "
                     "outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP (IntElement) projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }

    // check GP projection (MASTER)
    if (mdt == DRT::Element::quad4 || mdt == DRT::Element::quad8 || mdt == DRT::Element::quad9)
    {
      if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlaneQuad: Master Gauss point "
                     "projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP (IntElement) projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }
    else
    {
      if (mxi[0] < -tol || mxi[1] < -tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol ||
          mxi[0] + mxi[1] > 1.0 + 2 * tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlaneQuad: Master Gauss point "
                     "projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP (IntElement) projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }

    // project Gauss point back to slave (parent) element
    // project Gauss point back to master (parent) element
    double psxi[2] = {0.0, 0.0};
    double pmxi[2] = {0.0, 0.0};
    double psprojalpha = 0.0;
    double pmprojalpha = 0.0;
    MORTAR::MortarProjector::Impl(sele)->ProjectGaussPointAuxn3D(
        globgp, auxn, sele, psxi, psprojalpha);
    MORTAR::MortarProjector::Impl(mele)->ProjectGaussPointAuxn3D(
        globgp, auxn, mele, pmxi, pmprojalpha);
    // sintele.MapToParent(sxi, psxi); // old way of doing it via affine map... wrong (popp 05/2016)
    // mintele.MapToParent(mxi, pmxi); // old way of doing it via affine map... wrong (popp 05/2016)

    // check GP projection (SLAVE)
    if (psdt == DRT::Element::quad4 || psdt == DRT::Element::quad8 || psdt == DRT::Element::quad9 ||
        psdt == DRT::Element::nurbs9)
    {
      if (psxi[0] < -1.0 - tol || psxi[1] < -1.0 - tol || psxi[0] > 1.0 + tol ||
          psxi[1] > 1.0 + tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlaneQuad: Slave Gauss point projection "
                     "outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << psxi[0] << " " << psxi[1] << std::endl;
      }
    }
    else
    {
      if (psxi[0] < -tol || psxi[1] < -tol || psxi[0] > 1.0 + tol || psxi[1] > 1.0 + tol ||
          psxi[0] + psxi[1] > 1.0 + 2 * tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlaneQuad: Slave Gauss point projection "
                     "outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << psxi[0] << " " << psxi[1] << std::endl;
      }
    }

    // check GP projection (MASTER)
    if (pmdt == DRT::Element::quad4 || pmdt == DRT::Element::quad8 || pmdt == DRT::Element::quad9 ||
        pmdt == DRT::Element::nurbs9)
    {
      if (pmxi[0] < -1.0 - tol || pmxi[1] < -1.0 - tol || pmxi[0] > 1.0 + tol ||
          pmxi[1] > 1.0 + tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlaneQuad: Master Gauss point "
                     "projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << pmxi[0] << " " << pmxi[1] << std::endl;
      }
    }
    else
    {
      if (pmxi[0] < -tol || pmxi[1] < -tol || pmxi[0] > 1.0 + tol || pmxi[1] > 1.0 + tol ||
          pmxi[0] + pmxi[1] > 1.0 + 2 * tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlaneQuad: Master Gauss point "
                     "projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << pmxi[0] << " " << pmxi[1] << std::endl;
      }
    }

    // evaluate Lagrange multiplier shape functions (on slave element)
    if (LagMultQuad() == INPAR::MORTAR::lagmult_const)
      sele.EvaluateShapeLagMultConst(ShapeFcn(), sxi, lmval, lmderiv, nrow);
    else if (bound)
      sele.EvaluateShapeLagMult(ShapeFcn(), psxi, lmval, lmderiv, nrow);
    else
    {
      sele.EvaluateShapeLagMult(ShapeFcn(), psxi, lmval, lmderiv, nrow);
      sintele.EvaluateShapeLagMult(ShapeFcn(), sxi, lmintval, lmintderiv, nintrow, false);
    }

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(psxi, sval, sderiv, nrow, linlm);
    if (dualquad3d) sele.EvaluateShape(psxi, svalmod, sderivmod, nrow, true);
    mele.EvaluateShape(pmxi, mval, mderiv, ncol);

    // evaluate 2nd deriv of trace space shape functions (on slave element)
    sele.Evaluate2ndDerivShape(psxi, ssecderiv, nrow);

    // evaluate global GP coordinate derivative
    LINALG::Matrix<3, 1> svalcell;
    LINALG::Matrix<3, 2> sderivcell;
    cell->EvaluateShape(eta, svalcell, sderivcell);

    GEN::pairedvector<int, LINALG::Matrix<3, 1>> lingp((nrow + ncol) * ndof);

    for (int v = 0; v < 3; ++v)
      for (int d = 0; d < 3; ++d)
        for (_CI p = (cell->GetDerivVertex(v))[d].begin(); p != (cell->GetDerivVertex(v))[d].end();
             ++p)
          lingp[p->first](d) += svalcell(v) * (p->second);

    // evalute the GP slave coordinate derivatives
    std::vector<GEN::pairedvector<int, double>> dsxigp(2, (nrow + ncol) * ndof);
    DerivXiGP3DAuxPlane(
        sintele, sxi, cell->Auxn(), dsxigp, sprojalpha, cell->GetDerivAuxn(), lingp);

    // evalute the GP master coordinate derivatives
    std::vector<GEN::pairedvector<int, double>> dmxigp(2, (nrow + ncol) * ndof);
    DerivXiGP3DAuxPlane(
        mintele, mxi, cell->Auxn(), dmxigp, mprojalpha, cell->GetDerivAuxn(), lingp);

    // evaluate the GP slave coordinate derivatives (parent element)
    // evaluate the GP master coordinate derivatives (parent element)
    std::vector<GEN::pairedvector<int, double>> dpsxigp(2, (nrow + ncol) * ndof);
    std::vector<GEN::pairedvector<int, double>> dpmxigp(2, (nrow + ncol) * ndof);
    DerivXiGP3DAuxPlane(
        sele, psxi, cell->Auxn(), dpsxigp, psprojalpha, cell->GetDerivAuxn(), lingp);
    DerivXiGP3DAuxPlane(
        mele, pmxi, cell->Auxn(), dpmxigp, pmprojalpha, cell->GetDerivAuxn(), lingp);
    // sintele.MapToParent(dsxigp,dpsxigp); // old way of doing it via affine map... wrong (popp
    // 05/2016) mintele.MapToParent(dmxigp,dpmxigp); // old way of doing it via affine map... wrong
    // (popp 05/2016)

    //**********************************************************************
    // frequently reused quantities
    //**********************************************************************
    double gpn[3] = {0.0, 0.0, 0.0};
    double gap[1] = {0.0};
    double jumpval[2] = {0.0, 0.0};  // jump for wear
    // double jumpvalv[2] = {0.0,0.0};  // jump for slipincr --> equal to jumpval
    double wearval[1] = {0.0};  // wear value
    double lengthn[1] = {0.0};  // length of gp normal gpn
    GEN::pairedvector<int, double> dsliptmatrixgp(
        (ncol * ndof) + linsize);                                    // deriv. of slip for wear
    GEN::pairedvector<int, double> dgapgp((ncol * ndof) + linsize);  // gap lin. without lm and jac.
    GEN::pairedvector<int, double> dweargp(
        (ncol * ndof) + linsize);  // wear lin. without lm and jac.
    std::vector<GEN::pairedvector<int, double>> dslipgp(
        2, ((ncol * ndof) + linsize));  // deriv. of slip for slipincr (xi, eta)
    std::vector<GEN::pairedvector<int, double>> dnmap_unit(
        3, linsize);  // deriv of x,y and z comp. of gpn (unit)

    //**********************************************************************
    // evaluate at GP and lin char. quantities
    //**********************************************************************
    // these cases don't need any special treatment for quadratic elements
    if (algo_ == INPAR::MORTAR::algorithm_gpts ||
        (lmtype == INPAR::MORTAR::lagmult_const && algo_ == INPAR::MORTAR::algorithm_mortar))
    {
      Gap_3D(
          sele, mele, sval, mval, sderiv, mderiv, gap, gpn, dpsxigp, dpmxigp, dgapgp, dnmap_unit);
      IntegrateGP_3D(sele, mele, sval, lmval, mval, sderiv, mderiv, lmderiv, dualmap, wgt, jac,
          jacintcellmap, gpn, dnmap_unit, gap[0], dgapgp, psxi, pmxi, dpsxigp, dpmxigp);
    }
    // special treatment for quadratic elements
    else
    {
      // integrate and lin gp gap
      if (ShapeFcn() == INPAR::MORTAR::shape_standard && lmtype == INPAR::MORTAR::lagmult_pwlin)
      {
        GP_3D_G_Quad_pwlin(sele, sintele, mele, sval, mval, lmintval, scoord, mcoord, sderiv,
            mderiv, gap, gpn, lengthn, jac, wgt, dsxigp, dmxigp, dgapgp, dnmap_unit);
      }
      else
      {
        if (linlm)
        {
          // declare and compute shape functions as well as derivatives for computation of gap
          // -> standard quadratic functions instead of only linear functions
          LINALG::SerialDenseVector svalgap(nrow);
          LINALG::SerialDenseMatrix sderivgap(nrow, 2, true);
          sele.EvaluateShape(psxi, svalgap, sderivgap, nrow);

          Gap_3D(sele, mele, svalgap, mval, sderivgap, mderiv, gap, gpn, dpsxigp, dpmxigp, dgapgp,
              dnmap_unit);
        }
        else
          Gap_3D(sele, mele, sval, mval, sderiv, mderiv, gap, gpn, dpsxigp, dpmxigp, dgapgp,
              dnmap_unit);

        GP_3D_wGap(sele, sval, lmval, gap, jac, wgt, true);
      }

      // compute cell D/M matrix *******************************************
      GP_3D_DM_Quad(sele, mele, sintele, lmval, lmintval, sval, mval, jac, wgt, nrow, nintrow, ncol,
          ndof, bound);

      //*******************************
      // WEAR stuff
      //*******************************
      // std. wear for all wear-algorithm types
      if (wear)
        GP_3D_Wear(sele, mele, sval, sderiv, mval, mderiv, lmval, lmderiv, lagmult, gpn, jac, wgt,
            jumpval, wearval, dsliptmatrixgp, dweargp, dsxigp, dmxigp, dnmap_unit, dualmap);

      // integrate T and E matrix for discr. wear
      if (WearType() == INPAR::WEAR::wear_primvar) GP_TE(sele, lmval, sval, jac, wgt, jumpval);

      //********************************************************************
      // compute cell linearization
      //********************************************************************
      if (ShapeFcn() == INPAR::MORTAR::shape_standard && lmtype == INPAR::MORTAR::lagmult_pwlin)
      {
        for (int iter = 0; iter < nintrow; ++iter)
        {
          // Lin DM
          GP_3D_DM_Quad_pwlin_Lin(iter, sele, sintele, mele, sval, mval, lmintval, sderiv, mderiv,
              lmintderiv, wgt, jac, dsxigp, dpsxigp, dpmxigp, jacintcellmap);

          // Lin gap
          GP_3D_G_Quad_pwlin_Lin(iter, sintele, sval, lmintval, sderiv, lmintderiv, *gap, gpn, jac,
              wgt, dgapgp, jacintcellmap, dsxigp);
        }
      }
      else
      {
        // Lin DM
        GP_3D_DM_Quad_Lin(duallin, sele, mele, sval, svalmod, mval, lmval, sderiv, mderiv, lmderiv,
            wgt, jac, dpsxigp, dpmxigp, jacintcellmap, dualmap, dualquad3d);

        for (int iter = 0; iter < nrow; ++iter)
        {
          // Lin gap
          GP_3D_G_Quad_Lin(iter, sele, mele, sval, svalmod, lmval, sderiv, lmderiv, *gap, gpn, jac,
              wgt, duallin, dgapgp, jacintcellmap, dpsxigp, dualmap, dualquad3d);

          // Lin wear matrices T and E for discr. wear
          if (WearType() == INPAR::WEAR::wear_primvar)
            GP_3D_TE_Lin(iter, sele, sval, lmval, sderiv, lmderiv, jac, wgt, jumpval, dsxigp,
                jacintcellmap, dsliptmatrixgp, dualmap);
        }
      }
    }
  }  // gp loop

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateDerivEle2D(MORTAR::MortarElement& sele,
    std::vector<MORTAR::MortarElement*> meles, bool* boundary_ele,
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr = Teuchos::null;
  if (not mparams_ptr.is_null())
  {
    cparams_ptr = Teuchos::rcp_dynamic_cast<CONTACT::ParamsInterface>(mparams_ptr, true);
  }
  IntegrateDerivEle2D(sele, meles, boundary_ele, cparams_ptr);
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize mortar terms                     farah 01/13|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateDerivEle2D(MORTAR::MortarElement& sele,
    std::vector<MORTAR::MortarElement*> meles, bool* boundary_ele,
    const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  // ********************************************************************
  // Check integrator input for non-reasonable quantities
  // *********************************************************************

  // explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivEle2D called without specific shape function defined!");

  // Petrov-Galerkin approach for LM not yet implemented for quadratic FE
  if (sele.Shape() == MORTAR::MortarElement::line3 &&
      ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
    dserror(
        "ERROR: Petrov-Galerkin approach not yet implemented for 2-D quadratic FE interpolation");

  // check for problem dimension
  if (Dim() != 2) dserror("ERROR: 2D integration method called for non-2D problem");

  // check input data
  for (int i = 0; i < (int)meles.size(); ++i)
  {
    if ((!sele.IsSlave()) || (meles[i]->IsSlave()))
      dserror("ERROR: IntegrateDerivEle2D called on a wrong type of MortarElement pair!");
  }

  // *********************************************************************
  // Define slave quantities
  // *********************************************************************

  // consider entire slave element --> parameter space [-1,1]
  double sxia = -1.0;
  double sxib = 1.0;

  // number of nodes (slave, master)
  int ndof = Dim();
  int nrow = 0;

  DRT::Node** mynodes = NULL;
  nrow = sele.NumNode();
  mynodes = sele.Nodes();

  // get slave element nodes themselves
  if (!mynodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow, 1);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow, 1);

  // get slave nodal coords for Jacobian / GP evaluation
  LINALG::SerialDenseMatrix scoord(3, nrow);
  sele.GetNodalCoords(scoord);

  // nodal coords from previous time step and lagrange mulitplier
  Teuchos::RCP<LINALG::SerialDenseMatrix> scoordold;
  Teuchos::RCP<LINALG::SerialDenseMatrix> mcoordold;
  Teuchos::RCP<LINALG::SerialDenseMatrix> lagmult;

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // prepare directional derivative of dual shape functions
  // this is only necessary for quadratic dual shape functions in 2D
  GEN::pairedvector<int, Epetra_SerialDenseMatrix> dualmap(
      nrow * ndof, 0, Epetra_SerialDenseMatrix(nrow, nrow));
  if ((ShapeFcn() == INPAR::MORTAR::shape_dual ||
          ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) &&
      (sele.Shape() == MORTAR::MortarElement::line3 ||
          sele.MoData().DerivDualShape() != Teuchos::null ||
          sele.Shape() == MORTAR::MortarElement::nurbs3))
    sele.DerivShapeDual(dualmap);

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k = 0; k < nrow; ++k)
  {
    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  // decide whether linear LM are used for quadratic FE here
  // and whether displacement shape fct. modification has to be considered or not
  // this is the case for dual linear Lagrange multipliers on line3 elements
  bool linlm = false;
  bool dualquad = false;
  if (LagMultQuad() == INPAR::MORTAR::lagmult_lin && sele.Shape() == DRT::Element::line3)
  {
    bound = false;  // crosspoints and linear LM NOT at the same time!!!!
    linlm = true;
    if (ShapeFcn() == INPAR::MORTAR::shape_dual) dualquad = true;
  }

  // get numerical integration type
  INPAR::MORTAR::IntType inttype =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(imortar_, "INTTYPE");

  //************************************************************************
  // Boundary Segmentation check -- HasProj()-check
  //************************************************************************
  if (inttype == INPAR::MORTAR::inttype_elements_BS)
    *boundary_ele = BoundarySegmCheck2D(sele, meles);

  int linsize = 0;
  for (int i = 0; i < nrow; ++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*>(mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  if (*boundary_ele == false || inttype == INPAR::MORTAR::inttype_elements)
  {
    //*************************************************************************
    //                loop over all Gauss points for integration
    //*************************************************************************
    for (int gp = 0; gp < nGP(); ++gp)
    {
      bool kink_projection = false;

      // coordinates and weight
      double eta[2] = {Coordinate(gp, 0), 0.0};
      double wgt = Weight(gp);

      // coordinate transformation sxi->eta (slave MortarElement->Overlap)
      double sxi[2] = {0.0, 0.0};
      sxi[0] = eta[0];

      // evaluate the two slave side Jacobians
      double dxdsxi = sele.Jacobian(sxi);

      // evaluate Lagrange multiplier shape functions (on slave element)
      if (LagMultQuad() == INPAR::MORTAR::lagmult_const)
        sele.EvaluateShapeLagMultConst(ShapeFcn(), sxi, lmval, lmderiv, nrow);
      else if (linlm)
        sele.EvaluateShapeLagMultLin(ShapeFcn(), sxi, lmval, lmderiv, nrow);
      else
        sele.EvaluateShapeLagMult(ShapeFcn(), sxi, lmval, lmderiv, nrow);

      // evaluate trace space shape functions
      sele.EvaluateShape(sxi, sval, sderiv, nrow, dualquad);

      //****************************************************************************************************************
      //                loop over all Master Elements
      //****************************************************************************************************************
      for (int nummaster = 0; nummaster < (int)meles.size(); ++nummaster)
      {
        // project Gauss point onto master element
        double mxi[2] = {0.0, 0.0};
        MORTAR::MortarProjector::Impl(sele, *meles[nummaster])
            ->ProjectGaussPoint2D(sele, sxi, *meles[nummaster], mxi);

        // gp on mele?
        if ((mxi[0] >= -1.0) && (mxi[0] <= 1.0) && (kink_projection == false))
        {
          kink_projection = true;

          int ncol = 0;

          ncol = meles[nummaster]->NumNode();
          ;

          LINALG::SerialDenseVector mval(ncol);
          LINALG::SerialDenseMatrix mderiv(ncol, 1);

          // get master nodal coords for Jacobian / GP evaluation
          LINALG::SerialDenseMatrix mcoord(3, ncol);
          meles[nummaster]->GetNodalCoords(mcoord);

          // evaluate trace space shape functions
          meles[nummaster]->EvaluateShape(mxi, mval, mderiv, ncol, false);

          // get directional derivatives of sxia, sxib, mxia, mxib --> derivatives of mxia/mxib not
          // required
          std::vector<GEN::pairedvector<int, double>> ximaps(4, linsize + ndof * ncol);
          bool startslave = true;
          bool endslave = true;
          double mxia = -0.1;  //--> arbitrary value
          double mxib = 0.1;   //--> arbitrary value
          DerivXiAB2D(sele, sxia, sxib, *meles[nummaster], mxia, mxib, ximaps, startslave, endslave,
              linsize);

          // evalute the GP slave coordinate derivatives --> no entries
          std::vector<GEN::pairedvector<int, double>> dsxigp(
              1, GEN::pairedvector<int, double>(linsize + ndof * ncol));
          for (_CI p = ximaps[0].begin(); p != ximaps[0].end(); ++p) dsxigp[0][p->first] = 0.0;

          // evalute the GP master coordinate derivatives
          std::vector<GEN::pairedvector<int, double>> dmxigp(
              1, GEN::pairedvector<int, double>(linsize + ndof * ncol));
          DerivXiGP2D(sele, *meles[nummaster], sxi[0], mxi[0], dsxigp[0], dmxigp[0], linsize);

          // evaluate the Jacobian derivative
          GEN::pairedvector<int, double> derivjac(nrow * ndof);
          sele.DerivJacobian(sxi, derivjac);  // direct derivative if xi^1_g does not change

          //**********************************************************************
          // frequently reused quantities
          //**********************************************************************
          double gpn[2] = {0.0, 0.0};  // normalized normal at gp
          std::vector<GEN::pairedvector<int, double>> dnmap_unit(
              2, (linsize + ndof * ncol));  // deriv of x and y comp. of gpn (unit)
          double gap[1] = {0.0};            // gap
          GEN::pairedvector<int, double> dgapgp(
              linsize + ndof * ncol);  // gap lin without weighting and jac

          //**********************************************************************
          // evaluate at GP and lin char. quantities
          //**********************************************************************
          if (dualquad && linlm)
          {
            // declare and compute shape functions as well as derivatives for computation of gap
            // -> standard quadratic functions instead of only linear functions
            LINALG::SerialDenseVector svalgap(nrow);
            LINALG::SerialDenseMatrix sderivgap(nrow, 1);
            sele.EvaluateShape(sxi, svalgap, sderivgap, nrow);

            Gap_2D(sele, *meles[nummaster], svalgap, mval, sderivgap, mderiv, gap, gpn, dsxigp,
                dmxigp, dgapgp, dnmap_unit);
          }
          else
            Gap_2D(sele, *meles[nummaster], sval, mval, sderiv, mderiv, gap, gpn, dsxigp, dmxigp,
                dgapgp, dnmap_unit);

          // integrate and lin gp gap
          IntegrateGP_2D(sele, *meles[nummaster], sval, lmval, mval, sderiv, mderiv, lmderiv,
              dualmap, wgt, dxdsxi, derivjac, gpn, dnmap_unit, gap[0], dgapgp, sxi, mxi, dsxigp,
              dmxigp);
        }
      }  // End Loop over all Master Elements
    }    // End Loop over all GP
  }      // boundary_ele check

  return;
}

/*----------------------------------------------------------------------*
 |  Integrate D                                              farah 09/14|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateD(
    MORTAR::MortarElement& sele, const Epetra_Comm& comm, bool lin)
{
  // ********************************************************************
  // Check integrator input for non-reasonable quantities
  // *********************************************************************

  // explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateD called without specific shape function defined!");

  // *********************************************************************
  // Define slave quantities
  // *********************************************************************
  // number of nodes (slave, master)
  int ndof = Dim();
  int nrow = 0;

  DRT::Node** mynodes = NULL;
  nrow = sele.NumNode();
  mynodes = sele.Nodes();

  // get slave element nodes themselves
  if (!mynodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow, ndof - 1);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow, ndof - 1);

  int linsize = 0;
  for (int i = 0; i < nrow; ++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*>(mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k = 0; k < nrow; ++k)
  {
    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  bool duallin = false;
  GEN::pairedvector<int, Epetra_SerialDenseMatrix> dualmap(
      nrow * ndof, 0, Epetra_SerialDenseMatrix(nrow, nrow));
  if ((sele.Shape() != MORTAR::MortarElement::tri3 and
          sele.Shape() != MORTAR::MortarElement::line2) ||
      sele.MoData().DerivDualShape() != Teuchos::null)
  {
    duallin = true;
    sele.DerivShapeDual(dualmap);
  }

  // d-matrix derivative
  // local for this element to avoid direct assembly into the nodes for performance reasons
  GEN::pairedvector<int, Epetra_SerialDenseMatrix> dMatrixDeriv(
      sele.NumNode() * Dim(), 0, Epetra_SerialDenseMatrix(sele.NumNode(), sele.NumNode()));

  //*************************************************************************
  //                loop over all Gauss points for integration
  //*************************************************************************
  for (int gp = 0; gp < nGP(); ++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp, 0), 0.0};
    double wgt = Weight(gp);
    if (ndof == 3) eta[1] = Coordinate(gp, 1);

    // coordinate transformation sxi->eta (slave MortarElement->Overlap)
    double sxi[2] = {0.0, 0.0};
    sxi[0] = eta[0];
    sxi[1] = eta[1];

    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);

    // evaluate linearizations *******************************************
    // evaluate the slave Jacobian derivative
    GEN::pairedvector<int, double> jacslavemap(nrow * ndof + linsize);
    sele.DerivJacobian(sxi, jacslavemap);

    // evaluate Lagrange multiplier shape functions (on slave element)
    sele.EvaluateShapeLagMult(INPAR::MORTAR::shape_dual, sxi, lmval, lmderiv, nrow);

    // evaluate trace space shape functions
    sele.EvaluateShape(sxi, sval, sderiv, nrow, false);

    //**********************************************************************
    // evaluate at GP and lin char. quantities
    //**********************************************************************
    // integrate D and M
    for (int j = 0; j < nrow; ++j)
    {
      CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(mynodes[j]);

      // integrate dseg (boundary modification)
      if (bound)
      {
        bool j_boundnode = cnode->IsOnBound();

        for (int k = 0; k < nrow; ++k)
        {
          CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mynodes[k]);
          bool k_boundnode = mnode->IsOnBound();

          // do not assemble off-diagonal terms if j,k are both non-boundary nodes
          if (!j_boundnode && !k_boundnode && (j != k)) continue;

          // multiply the two shape functions
          double prod = lmval[j] * sval[k] * dxdsxi * wgt;

          // isolate the dseg entries to be filled
          // (both the main diagonal and every other secondary diagonal)
          // and add current Gauss point's contribution to dseg
          if (mnode->IsOnBound())
          {
            if (abs(prod) > MORTARINTTOL) cnode->AddMValue(mnode->Id(), -prod);
            if (abs(prod) > MORTARINTTOL) cnode->AddMNode(mnode->Id());  // only for friction!
          }
          else
          {
            if (abs(prod) > MORTARINTTOL) cnode->AddDValue(mnode->Id(), prod);
            if (abs(prod) > MORTARINTTOL) cnode->AddSNode(mnode->Id());  // only for friction!
          }
        }
      }
      else
      {
        // integrate dseg
        for (int k = 0; k < nrow; ++k)
        {
          CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(mynodes[k]);

          // multiply the two shape functions
          double prod = lmval[j] * sval[k] * dxdsxi * wgt;


          if (sele.IsSlave())
          {
            if (abs(prod) > MORTARINTTOL) cnode->AddDValue(snode->Id(), prod);
          }

          // loop over slave dofs
          for (int jdof = 0; jdof < ndof; ++jdof)
          {
            int col = snode->Dofs()[jdof];

            if (sele.IsSlave())
            {
            }
            else
            {
              if (sele.Owner() == comm.MyPID())
              {
                if (abs(prod) > MORTARINTTOL)
                  dynamic_cast<CONTACT::FriNode*>(cnode)->AddD2Value(jdof, col, prod);
              }
            }
          }
        }
      }
    }

    if (lin)
    {
      typedef GEN::pairedvector<int, double>::const_iterator _CI;

      // (1) Lin(Phi) - dual shape functions
      if (duallin)
      {
        for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
             p != dualmap.end(); ++p)
        {
          Epetra_SerialDenseMatrix& dderivtmp = dMatrixDeriv[p->first];
          for (int j = 0; j < nrow; ++j)
            for (int k = 0; k < nrow; ++k)
              for (int m = 0; m < nrow; ++m)
                dderivtmp(j, j) += wgt * sval[m] * sval[k] * dxdsxi * (p->second)(j, m);
        }
      }

      // (4) Lin(dsxideta) - intcell GP Jacobian
      for (_CI p = jacslavemap.begin(); p != jacslavemap.end(); ++p)
      {
        Epetra_SerialDenseMatrix& dderivtmp = dMatrixDeriv[p->first];
        for (int j = 0; j < nrow; ++j)
          for (int k = 0; k < nrow; ++k) dderivtmp(j, j) += wgt * lmval[j] * sval[k] * (p->second);
      }
    }
  }  // End Loop over all GP

  sele.MoData().ResetDualShape();
  sele.MoData().ResetDerivDualShape();

  return;
}


/*----------------------------------------------------------------------*
 |  Compute penalty scaling factor kappa                     farah 11/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateKappaPenaltyLTS(MORTAR::MortarElement& ele)
{
  // number of nodes (slave)
  int nrow = ele.NumNode();

  // create empty objects for shape fct. evaluation
  LINALG::SerialDenseVector val(nrow);
  LINALG::SerialDenseMatrix deriv(nrow, 2, true);

  // get slave element nodes themselves
  DRT::Node** mynodes = ele.Nodes();

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < nGP(); ++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp, 0), 0.0};
    if (Dim() == 3) eta[1] = Coordinate(gp, 1);
    double wgt = Weight(gp);

    // evaluate shape functions
    ele.EvaluateShapeLagMult(ShapeFcn(), eta, val, deriv, nrow);

    // evaluate the Jacobian det
    const double jac = ele.Jacobian(eta);

    // compute cell gap vector *******************************************
    // loop over all gseg vector entries
    // nrow represents the slave side dofs !!!  */
    for (int j = 0; j < nrow; ++j)
    {
      CoNode* mymrtrnode = dynamic_cast<CoNode*>(mynodes[j]);

      if (mymrtrnode->IsOnCorner()) continue;

      // add current Gauss point's contribution to kappa
      mymrtrnode->CoData().Kappa() += val[j] * jac * wgt;
    }
  }
  //**********************************************************************

  return;
}


/*----------------------------------------------------------------------*
 |  Compute penalty scaling factor kappa                      popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateKappaPenalty(MORTAR::MortarElement& sele, double* sxia,
    double* sxib, Teuchos::RCP<Epetra_SerialDenseVector> gseg)
{
  // explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateKappaPenalty called without specific shape function defined!");

  // check input data
  if (!sele.IsSlave()) dserror("ERROR: IntegrateKappaPenalty called on a non-slave MortarElement!");
  if ((sxia[0] < -1.0) || (sxia[1] < -1.0) || (sxib[0] > 1.0) || (sxib[1] > 1.0))
    dserror("ERROR: IntegrateKappaPenalty called with infeasible slave limits!");

  // number of nodes (slave)
  int nrow = sele.NumNode();

  // create empty objects for shape fct. evaluation
  LINALG::SerialDenseVector val(nrow);
  LINALG::SerialDenseMatrix deriv(nrow, 2, true);

  // map iterator
  // typedef std::map<int,double>::const_iterator CI;

  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if (!mynodes) dserror("ERROR: IntegrateKappaPenalty: Null pointer!");

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k = 0; k < nrow; ++k)
  {
    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode) dserror("ERROR: IntegrateKappaPenalty: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < nGP(); ++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp, 0), 0.0};
    if (Dim() == 3) eta[1] = Coordinate(gp, 1);
    double wgt = Weight(gp);

    // evaluate shape functions
    if (bound)
      sele.EvaluateShapeLagMultLin(ShapeFcn(), eta, val, deriv, nrow);
    else
      sele.EvaluateShapeLagMult(ShapeFcn(), eta, val, deriv, nrow);

    // evaluate the Jacobian det
    const double jac = sele.Jacobian(eta);

    // compute cell gap vector *******************************************
    // loop over all gseg vector entries
    // nrow represents the slave side dofs !!!  */
    for (int j = 0; j < nrow; ++j)
    {
      // add current Gauss point's contribution to gseg
      (*gseg)(j) += val[j] * jac * wgt;
    }
    // compute cell gap vector *******************************************
  }
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 |  Compute penalty scaling factor kappa (3D piecewise lin)   popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateKappaPenalty(MORTAR::MortarElement& sele,
    MORTAR::IntElement& sintele, double* sxia, double* sxib,
    Teuchos::RCP<Epetra_SerialDenseVector> gseg)
{
  // get LMtype
  INPAR::MORTAR::LagMultQuad lmtype = LagMultQuad();

  // explicitly defined shape function type needed
  if (ShapeFcn() != INPAR::MORTAR::shape_standard)
    dserror("ERROR: IntegrateKappaPenalty -> you should not be here!");

  // check input data
  if (!sele.IsSlave()) dserror("ERROR: IntegrateKappaPenalty called on a non-slave MortarElement!");
  if ((sxia[0] < -1.0) || (sxia[1] < -1.0) || (sxib[0] > 1.0) || (sxib[1] > 1.0))
    dserror("ERROR: IntegrateKappaPenalty called with infeasible slave limits!");

  // number of nodes (slave)
  int nrow = sele.NumNode();
  int nintrow = sintele.NumNode();

  // create empty objects for shape fct. evaluation
  LINALG::SerialDenseVector val(nrow);
  LINALG::SerialDenseMatrix deriv(nrow, 2, true);
  LINALG::SerialDenseVector intval(nintrow);
  LINALG::SerialDenseMatrix intderiv(nintrow, 2, true);

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < nGP(); ++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp, 0), 0.0};
    if (Dim() == 3) eta[1] = Coordinate(gp, 1);
    double wgt = Weight(gp);

    // get global Gauss point coordinates
    double globgp[3] = {0.0, 0.0, 0.0};
    sintele.LocalToGlobal(eta, globgp, 0);

    // get normal vector
    double auxn[3] = {0.0, 0.0, 0.0};
    sintele.ComputeUnitNormalAtXi(eta, auxn);

    // project Gauss point back to slave (parent) element
    double psxi[2] = {0.0, 0.0};
    double psprojalpha = 0.0;
    MORTAR::MortarProjector::Impl(sele)->ProjectGaussPointAuxn3D(
        globgp, auxn, sele, psxi, psprojalpha);
    // sintele.MapToParent(eta,psxi); //old way of doing it via affine map... wrong (popp 05/2016)

    // evaluate shape functions
    sele.EvaluateShape(psxi, val, deriv, nrow);
    sintele.EvaluateShape(eta, intval, intderiv, nintrow);

    // evaluate the Jacobian det
    const double jac = sintele.Jacobian(eta);

    // compute cell gap vector *******************************************
    if (lmtype == INPAR::MORTAR::lagmult_pwlin)
    {
      for (int j = 0; j < nintrow; ++j)
      {
        // add current Gauss point's contribution to gseg
        (*gseg)(j) += intval[j] * jac * wgt;
      }
    }
    else
    {
      dserror("ERROR: Invalid LM interpolation case!");
    }
    // compute cell gap vector *******************************************
  }
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 |  Compute directional derivative of XiAB (2D)               popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::DerivXiAB2D(MORTAR::MortarElement& sele, double& sxia, double& sxib,
    MORTAR::MortarElement& mele, double& mxia, double& mxib,
    std::vector<GEN::pairedvector<int, double>>& derivxi, bool& startslave, bool& endslave,
    int& linsize)
{
  // check for problem dimension
  if (Dim() != 2) dserror("ERROR: 2D integration method called for non-2D problem");

  // we need the participating slave and master nodes
  DRT::Node** snodes = NULL;
  DRT::Node** mnodes = NULL;
  int numsnode = sele.NumNode();
  int nummnode = mele.NumNode();

  int ndof = Dim();

  snodes = sele.Nodes();
  mnodes = mele.Nodes();

  std::vector<MORTAR::MortarNode*> smrtrnodes(numsnode);
  std::vector<MORTAR::MortarNode*> mmrtrnodes(nummnode);

  for (int i = 0; i < numsnode; ++i)
  {
    smrtrnodes[i] = dynamic_cast<MORTAR::MortarNode*>(snodes[i]);
    if (!smrtrnodes[i]) dserror("ERROR: DerivXiAB2D: Null pointer!");
  }

  for (int i = 0; i < nummnode; ++i)
  {
    mmrtrnodes[i] = dynamic_cast<MORTAR::MortarNode*>(mnodes[i]);
    if (!mmrtrnodes[i]) dserror("ERROR: DerivXiAB2D: Null pointer!");
  }

  // we also need shape function derivs in A and B
  double psxia[2] = {sxia, 0.0};
  double psxib[2] = {sxib, 0.0};
  double pmxia[2] = {mxia, 0.0};
  double pmxib[2] = {mxib, 0.0};
  LINALG::SerialDenseVector valsxia(numsnode);
  LINALG::SerialDenseVector valsxib(numsnode);
  LINALG::SerialDenseVector valmxia(nummnode);
  LINALG::SerialDenseVector valmxib(nummnode);
  LINALG::SerialDenseMatrix derivsxia(numsnode, 1);
  LINALG::SerialDenseMatrix derivsxib(numsnode, 1);
  LINALG::SerialDenseMatrix derivmxia(nummnode, 1);
  LINALG::SerialDenseMatrix derivmxib(nummnode, 1);

  sele.EvaluateShape(psxia, valsxia, derivsxia, numsnode, false);
  sele.EvaluateShape(psxib, valsxib, derivsxib, numsnode, false);
  mele.EvaluateShape(pmxia, valmxia, derivmxia, nummnode, false);
  mele.EvaluateShape(pmxib, valmxib, derivmxib, nummnode, false);

  // prepare linearizations
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // compute leading constant for DerivXiBMaster if start node = slave node
  if (startslave == true)
  {
    // compute factors and leading constants for master
    double cmxib = 0.0;
    double fac_dxm_b = 0.0;
    double fac_dym_b = 0.0;
    double fac_xmsl_b = 0.0;
    double fac_ymsl_b = 0.0;

    // there are only actual 2 components for the normal in 2D but ComputeUnitNormalAtXi expects 3
    // components (otherwise invalid memory access)
    double normal[3] = {0., 0., 0.};
    sele.ComputeUnitNormalAtXi(psxia, normal);
    std::vector<GEN::pairedvector<int, double>> derivN;
    dynamic_cast<CONTACT::CoElement*>(&sele)->DerivUnitNormalAtXi(psxia, derivN);

    LINALG::SerialDenseVector* mval = NULL;
    LINALG::SerialDenseMatrix* mderiv = NULL;
    if (sele.NormalFac() * mele.NormalFac() > 0.)
    {
      mval = &valmxib;
      mderiv = &derivmxib;
    }
    else
    {
      mval = &valmxia;
      mderiv = &derivmxia;
    }

    for (int i = 0; i < nummnode; ++i)
    {
      fac_dxm_b += (*mderiv)(i, 0) * (mmrtrnodes[i]->xspatial()[0]);
      fac_dym_b += (*mderiv)(i, 0) * (mmrtrnodes[i]->xspatial()[1]);
      fac_xmsl_b += (*mval)[i] * (mmrtrnodes[i]->xspatial()[0]);
      fac_ymsl_b += (*mval)[i] * (mmrtrnodes[i]->xspatial()[1]);
    }

    cmxib = -1 / (fac_dxm_b * (normal[1]) - fac_dym_b * (normal[0]));
    // std::cout << "cmxib: " << cmxib << std::endl;

    for (int i = 0; i < numsnode; ++i)
    {
      fac_xmsl_b -= valsxia[i] * (smrtrnodes[i]->xspatial()[0]);
      fac_ymsl_b -= valsxia[i] * (smrtrnodes[i]->xspatial()[1]);
    }

    GEN::pairedvector<int, double> dmap_mxib(nummnode * ndof + linsize);

    // add derivative of slave node coordinates
    for (int i = 0; i < numsnode; ++i)
    {
      dmap_mxib[smrtrnodes[i]->Dofs()[0]] -= valsxia[i] * normal[1];
      dmap_mxib[smrtrnodes[i]->Dofs()[1]] += valsxia[i] * normal[0];
    }
    // add derivatives of master node coordinates
    for (int i = 0; i < nummnode; ++i)
    {
      dmap_mxib[mmrtrnodes[i]->Dofs()[0]] += (*mval)[i] * (normal[1]);
      dmap_mxib[mmrtrnodes[i]->Dofs()[1]] -= (*mval)[i] * (normal[0]);
    }

    // add derivative of slave node normal
    for (_CI p = derivN[0].begin(); p != derivN[0].end(); ++p)
      dmap_mxib[p->first] -= fac_ymsl_b * (p->second);
    for (_CI p = derivN[1].begin(); p != derivN[1].end(); ++p)
      dmap_mxib[p->first] += fac_xmsl_b * (p->second);

    // multiply all entries with cmxib
    for (_CI p = dmap_mxib.begin(); p != dmap_mxib.end(); ++p)
      dmap_mxib[p->first] = cmxib * (p->second);

    // return map to DerivM() method
    if (sele.NormalFac() * mele.NormalFac() > 0.)
      derivxi[3] = dmap_mxib;
    else
      derivxi[2] = dmap_mxib;
  }

  // compute leading constant for DerivXiAMaster if end node = slave node
  if (endslave == true)
  {
    // compute factors and leading constants for master
    double cmxia = 0.0;
    double fac_dxm_a = 0.0;
    double fac_dym_a = 0.0;
    double fac_xmsl_a = 0.0;
    double fac_ymsl_a = 0.0;

    // there are only 2 actual components for the normal in 2D but ComputeUnitNormalAtXi expects 3
    // components (otherwise invalid memory access)
    double normal[3] = {0., 0., 0.};
    sele.ComputeUnitNormalAtXi(psxib, normal);
    std::vector<GEN::pairedvector<int, double>> derivN;
    dynamic_cast<CONTACT::CoElement*>(&sele)->DerivUnitNormalAtXi(psxib, derivN);

    LINALG::SerialDenseVector* mval = NULL;
    LINALG::SerialDenseMatrix* mderiv = NULL;
    if (sele.NormalFac() * mele.NormalFac() > 0.)
    {
      mval = &valmxia;
      mderiv = &derivmxia;
    }
    else
    {
      mval = &valmxib;
      mderiv = &derivmxib;
    }

    for (int i = 0; i < nummnode; ++i)
    {
      fac_dxm_a += (*mderiv)(i, 0) * (mmrtrnodes[i]->xspatial()[0]);
      fac_dym_a += (*mderiv)(i, 0) * (mmrtrnodes[i]->xspatial()[1]);
      fac_xmsl_a += (*mval)[i] * (mmrtrnodes[i]->xspatial()[0]);
      fac_ymsl_a += (*mval)[i] * (mmrtrnodes[i]->xspatial()[1]);
    }

    cmxia = -1 / (fac_dxm_a * (smrtrnodes[1]->MoData().n()[1]) -
                     fac_dym_a * (smrtrnodes[1]->MoData().n()[0]));
    // std::cout << "cmxia: " << cmxia << std::endl;

    for (int i = 0; i < numsnode; ++i)
    {
      fac_xmsl_a -= valsxib[i] * (smrtrnodes[i]->xspatial()[0]);
      fac_ymsl_a -= valsxib[i] * (smrtrnodes[i]->xspatial()[1]);
    }

    GEN::pairedvector<int, double> dmap_mxia(nummnode * ndof + linsize);

    // add derivative of slave node coordinates
    for (int i = 0; i < numsnode; ++i)
    {
      dmap_mxia[smrtrnodes[i]->Dofs()[0]] -= valsxib[i] * normal[1];
      dmap_mxia[smrtrnodes[i]->Dofs()[1]] += valsxib[i] * normal[0];
    }

    // add derivatives of master node coordinates
    for (int i = 0; i < nummnode; ++i)
    {
      dmap_mxia[mmrtrnodes[i]->Dofs()[0]] += (*mval)[i] * (normal[1]);
      dmap_mxia[mmrtrnodes[i]->Dofs()[1]] -= (*mval)[i] * (normal[0]);
    }

    // add derivative of slave node normal
    for (_CI p = derivN[0].begin(); p != derivN[0].end(); ++p)
      dmap_mxia[p->first] -= fac_ymsl_a * (p->second);
    for (_CI p = derivN[1].begin(); p != derivN[1].end(); ++p)
      dmap_mxia[p->first] += fac_xmsl_a * (p->second);

    // multiply all entries with cmxia
    for (_CI p = dmap_mxia.begin(); p != dmap_mxia.end(); ++p)
      dmap_mxia[p->first] = cmxia * (p->second);

    // return map to DerivM() method
    if (sele.NormalFac() * mele.NormalFac() > 0.)
      derivxi[2] = dmap_mxia;
    else
      derivxi[3] = dmap_mxia;
  }

  // compute leading constant for DerivXiASlave if start node = master node
  if (startslave == false)
  {
    // compute factors and leading constants for slave
    double csxia = 0.0;
    double fac_dxsl_a = 0.0;
    double fac_dysl_a = 0.0;
    double fac_xslm_a = 0.0;
    double fac_yslm_a = 0.0;
    double fac_dnx_a = 0.0;
    double fac_dny_a = 0.0;
    double fac_nx_a = 0.0;
    double fac_ny_a = 0.0;

    LINALG::SerialDenseVector* mval = NULL;
    if (sele.NormalFac() * mele.NormalFac() > 0.)
      mval = &valmxib;
    else
      mval = &valmxia;

    for (int i = 0; i < numsnode; ++i)
    {
      fac_dxsl_a += derivsxia(i, 0) * (smrtrnodes[i]->xspatial()[0]);
      fac_dysl_a += derivsxia(i, 0) * (smrtrnodes[i]->xspatial()[1]);
      fac_xslm_a += valsxia[i] * (smrtrnodes[i]->xspatial()[0]);
      fac_yslm_a += valsxia[i] * (smrtrnodes[i]->xspatial()[1]);
      fac_dnx_a += derivsxia(i, 0) * (smrtrnodes[i]->MoData().n()[0]);
      fac_dny_a += derivsxia(i, 0) * (smrtrnodes[i]->MoData().n()[1]);
      fac_nx_a += valsxia[i] * (smrtrnodes[i]->MoData().n()[0]);
      fac_ny_a += valsxia[i] * (smrtrnodes[i]->MoData().n()[1]);
    }

    for (int i = 0; i < nummnode; ++i)
    {
      fac_xslm_a -= (*mval)[i] * (mmrtrnodes[i]->xspatial()[0]);
      fac_yslm_a -= (*mval)[i] * (mmrtrnodes[i]->xspatial()[1]);
    }

    csxia = -1. / (fac_dxsl_a * fac_ny_a - fac_dysl_a * fac_nx_a + fac_xslm_a * fac_dny_a -
                      fac_yslm_a * fac_dnx_a);
    // std::cout << "csxia: " << csxia << std::endl;


    GEN::pairedvector<int, double> dmap_sxia(nummnode * ndof + linsize);

    // add derivative of master node coordinates
    for (int i = 0; i < nummnode; ++i)
    {
      dmap_sxia[mmrtrnodes[i]->Dofs()[0]] -= (*mval)[i] * fac_ny_a;
      dmap_sxia[mmrtrnodes[i]->Dofs()[1]] += (*mval)[i] * fac_nx_a;
    }

    // add derivatives of slave node coordinates
    for (int i = 0; i < numsnode; ++i)
    {
      dmap_sxia[smrtrnodes[i]->Dofs()[0]] += valsxia[i] * fac_ny_a;
      dmap_sxia[smrtrnodes[i]->Dofs()[1]] -= valsxia[i] * fac_nx_a;
    }

    // add derivatives of slave node normals
    for (int i = 0; i < numsnode; ++i)
    {
      GEN::pairedvector<int, double>& nxmap_curr =
          static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
      GEN::pairedvector<int, double>& nymap_curr =
          static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];

      for (_CI p = nxmap_curr.begin(); p != nxmap_curr.end(); ++p)
        dmap_sxia[p->first] -= valsxia[i] * fac_yslm_a * (p->second);
      for (_CI p = nymap_curr.begin(); p != nymap_curr.end(); ++p)
        dmap_sxia[p->first] += valsxia[i] * fac_xslm_a * (p->second);
    }

    // multiply all entries with csxia
    for (_CI p = dmap_sxia.begin(); p != dmap_sxia.end(); ++p)
      dmap_sxia[p->first] = csxia * (p->second);

    // return map to DerivM() method
    derivxi[0] = dmap_sxia;
  }

  // compute leading constant for DerivXiBSlave if end node = master node
  if (endslave == false)
  {
    // compute factors and leading constants for slave
    double csxib = 0.0;
    double fac_dxsl_b = 0.0;
    double fac_dysl_b = 0.0;
    double fac_xslm_b = 0.0;
    double fac_yslm_b = 0.0;
    double fac_dnx_b = 0.0;
    double fac_dny_b = 0.0;
    double fac_nx_b = 0.0;
    double fac_ny_b = 0.0;

    LINALG::SerialDenseVector* mval = NULL;
    if (sele.NormalFac() * mele.NormalFac() > 0.)
      mval = &valmxia;
    else
      mval = &valmxib;

    for (int i = 0; i < numsnode; ++i)
    {
      fac_dxsl_b += derivsxib(i, 0) * (smrtrnodes[i]->xspatial()[0]);
      fac_dysl_b += derivsxib(i, 0) * (smrtrnodes[i]->xspatial()[1]);
      fac_xslm_b += valsxib[i] * (smrtrnodes[i]->xspatial()[0]);
      fac_yslm_b += valsxib[i] * (smrtrnodes[i]->xspatial()[1]);
      fac_dnx_b += derivsxib(i, 0) * (smrtrnodes[i]->MoData().n()[0]);
      fac_dny_b += derivsxib(i, 0) * (smrtrnodes[i]->MoData().n()[1]);
      fac_nx_b += valsxib[i] * (smrtrnodes[i]->MoData().n()[0]);
      fac_ny_b += valsxib[i] * (smrtrnodes[i]->MoData().n()[1]);
    }

    for (int i = 0; i < nummnode; ++i)
    {
      fac_xslm_b -= (*mval)[i] * (mmrtrnodes[i]->xspatial()[0]);
      fac_yslm_b -= (*mval)[i] * (mmrtrnodes[i]->xspatial()[1]);
    }

    csxib = -1 / (fac_dxsl_b * fac_ny_b - fac_dysl_b * fac_nx_b + fac_xslm_b * fac_dny_b -
                     fac_yslm_b * fac_dnx_b);
    // std::cout << "csxib: " << csxib << std::endl;

    GEN::pairedvector<int, double> dmap_sxib(nummnode * ndof + linsize);

    // add derivative of master node coordinates
    for (int i = 0; i < nummnode; ++i)
    {
      dmap_sxib[mmrtrnodes[i]->Dofs()[0]] -= (*mval)[i] * fac_ny_b;
      dmap_sxib[mmrtrnodes[i]->Dofs()[1]] += (*mval)[i] * fac_nx_b;
    }

    // add derivatives of slave node coordinates
    for (int i = 0; i < numsnode; ++i)
    {
      dmap_sxib[smrtrnodes[i]->Dofs()[0]] += valsxib[i] * fac_ny_b;
      dmap_sxib[smrtrnodes[i]->Dofs()[1]] -= valsxib[i] * fac_nx_b;
    }

    // add derivatives of slave node normals
    for (int i = 0; i < numsnode; ++i)
    {
      GEN::pairedvector<int, double>& nxmap_curr =
          static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
      GEN::pairedvector<int, double>& nymap_curr =
          static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];

      for (_CI p = nxmap_curr.begin(); p != nxmap_curr.end(); ++p)
        dmap_sxib[p->first] -= valsxib[i] * fac_yslm_b * (p->second);
      for (_CI p = nymap_curr.begin(); p != nymap_curr.end(); ++p)
        dmap_sxib[p->first] += valsxib[i] * fac_xslm_b * (p->second);
    }

    // multiply all entries with csxib
    for (_CI p = dmap_sxib.begin(); p != dmap_sxib.end(); ++p)
      dmap_sxib[p->first] = csxib * (p->second);

    // return map to DerivM() method
    derivxi[1] = dmap_sxib;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Compute directional derivative of XiGP master (2D)        popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::DerivXiGP2D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    double sxigp, double mxigp, const GEN::pairedvector<int, double>& derivsxi,
    GEN::pairedvector<int, double>& derivmxi, int& linsize)
{
  // check for problem dimension
  if (Dim() != 2) dserror("ERROR: 2D integration method called for non-2D problem");

  // we need the participating slave and master nodes
  DRT::Node** snodes = NULL;
  DRT::Node** mnodes = NULL;
  int numsnode = sele.NumNode();
  int nummnode = mele.NumNode();

  int ndof = Dim();

  snodes = sele.Nodes();
  mnodes = mele.Nodes();

  std::vector<MORTAR::MortarNode*> smrtrnodes(numsnode);
  std::vector<MORTAR::MortarNode*> mmrtrnodes(nummnode);

  for (int i = 0; i < numsnode; ++i)
  {
    smrtrnodes[i] = dynamic_cast<MORTAR::MortarNode*>(snodes[i]);
    if (!smrtrnodes[i]) dserror("ERROR: DerivXiAB2D: Null pointer!");
  }

  for (int i = 0; i < nummnode; ++i)
  {
    mmrtrnodes[i] = dynamic_cast<MORTAR::MortarNode*>(mnodes[i]);
    if (!mmrtrnodes[i]) dserror("ERROR: DerivXiAB2D: Null pointer!");
  }

  // we also need shape function derivs in A and B
  double psxigp[2] = {sxigp, 0.0};
  double pmxigp[2] = {mxigp, 0.0};
  LINALG::SerialDenseVector valsxigp(numsnode);
  LINALG::SerialDenseVector valmxigp(nummnode);
  LINALG::SerialDenseMatrix derivsxigp(numsnode, 1);
  LINALG::SerialDenseMatrix derivmxigp(nummnode, 1);

  sele.EvaluateShape(psxigp, valsxigp, derivsxigp, numsnode, false);
  mele.EvaluateShape(pmxigp, valmxigp, derivmxigp, nummnode, false);

  // we also need the GP slave coordinates + normal
  double sgpn[3] = {0.0, 0.0, 0.0};
  double sgpx[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < numsnode; ++i)
  {
    sgpn[0] += valsxigp[i] * smrtrnodes[i]->MoData().n()[0];
    sgpn[1] += valsxigp[i] * smrtrnodes[i]->MoData().n()[1];
    sgpn[2] += valsxigp[i] * smrtrnodes[i]->MoData().n()[2];

    sgpx[0] += valsxigp[i] * smrtrnodes[i]->xspatial()[0];
    sgpx[1] += valsxigp[i] * smrtrnodes[i]->xspatial()[1];
    sgpx[2] += valsxigp[i] * smrtrnodes[i]->xspatial()[2];
  }

  // FIXME: This does not have to be the UNIT normal (see 3D)!
  // The reason for this is that we linearize the Gauss point
  // projection from slave to master side here and this condition
  // only includes the Gauss point normal in a cross product.
  // When looking at MortarProjector::ProjectGaussPoint, one can see
  // that we do NOT use a unit normal there, either. Thus, why here?
  // First results suggest that it really makes no difference!

  // normalize interpolated GP normal back to length 1.0 !!!
  const double length = sqrt(sgpn[0] * sgpn[0] + sgpn[1] * sgpn[1] + sgpn[2] * sgpn[2]);
  if (length < 1.0e-12) dserror("ERROR: DerivXiGP2D: Divide by zero!");
  for (int i = 0; i < 3; ++i) sgpn[i] /= length;

  // compute factors and leading constants for master
  double cmxigp = 0.0;
  double fac_dxm_gp = 0.0;
  double fac_dym_gp = 0.0;
  double fac_xmsl_gp = 0.0;
  double fac_ymsl_gp = 0.0;

  for (int i = 0; i < nummnode; ++i)
  {
    fac_dxm_gp += derivmxigp(i, 0) * (mmrtrnodes[i]->xspatial()[0]);
    fac_dym_gp += derivmxigp(i, 0) * (mmrtrnodes[i]->xspatial()[1]);

    fac_xmsl_gp += valmxigp[i] * (mmrtrnodes[i]->xspatial()[0]);
    fac_ymsl_gp += valmxigp[i] * (mmrtrnodes[i]->xspatial()[1]);
  }

  cmxigp = -1 / (fac_dxm_gp * sgpn[1] - fac_dym_gp * sgpn[0]);
  // std::cout << "cmxigp: " << cmxigp << std::endl;

  fac_xmsl_gp -= sgpx[0];
  fac_ymsl_gp -= sgpx[1];

  // prepare linearization
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // build directional derivative of slave GP coordinates
  GEN::pairedvector<int, double> dmap_xsl_gp(linsize + nummnode * ndof);
  GEN::pairedvector<int, double> dmap_ysl_gp(linsize + nummnode * ndof);

  for (int i = 0; i < numsnode; ++i)
  {
    dmap_xsl_gp[smrtrnodes[i]->Dofs()[0]] += valsxigp[i];
    dmap_ysl_gp[smrtrnodes[i]->Dofs()[1]] += valsxigp[i];

    for (_CI p = derivsxi.begin(); p != derivsxi.end(); ++p)
    {
      double facx = derivsxigp(i, 0) * (smrtrnodes[i]->xspatial()[0]);
      double facy = derivsxigp(i, 0) * (smrtrnodes[i]->xspatial()[1]);
      dmap_xsl_gp[p->first] += facx * (p->second);
      dmap_ysl_gp[p->first] += facy * (p->second);
    }
  }

  // build directional derivative of slave GP normal
  GEN::pairedvector<int, double> dmap_nxsl_gp(linsize + nummnode * ndof);
  GEN::pairedvector<int, double> dmap_nysl_gp(linsize + nummnode * ndof);

  double sgpnmod[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i) sgpnmod[i] = sgpn[i] * length;

  GEN::pairedvector<int, double> dmap_nxsl_gp_mod(linsize + nummnode * ndof);
  GEN::pairedvector<int, double> dmap_nysl_gp_mod(linsize + nummnode * ndof);

  for (int i = 0; i < numsnode; ++i)
  {
    GEN::pairedvector<int, double>& dmap_nxsl_i =
        static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
    GEN::pairedvector<int, double>& dmap_nysl_i =
        static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];

    for (_CI p = dmap_nxsl_i.begin(); p != dmap_nxsl_i.end(); ++p)
      dmap_nxsl_gp_mod[p->first] += valsxigp[i] * (p->second);
    for (_CI p = dmap_nysl_i.begin(); p != dmap_nysl_i.end(); ++p)
      dmap_nysl_gp_mod[p->first] += valsxigp[i] * (p->second);

    for (_CI p = derivsxi.begin(); p != derivsxi.end(); ++p)
    {
      double valx = derivsxigp(i, 0) * smrtrnodes[i]->MoData().n()[0];
      dmap_nxsl_gp_mod[p->first] += valx * (p->second);
      double valy = derivsxigp(i, 0) * smrtrnodes[i]->MoData().n()[1];
      dmap_nysl_gp_mod[p->first] += valy * (p->second);
    }
  }

  const double sxsx = sgpnmod[0] * sgpnmod[0];
  const double sxsy = sgpnmod[0] * sgpnmod[1];
  const double sysy = sgpnmod[1] * sgpnmod[1];
  const double linv = 1.0 / length;
  const double lllinv = 1.0 / (length * length * length);

  for (_CI p = dmap_nxsl_gp_mod.begin(); p != dmap_nxsl_gp_mod.end(); ++p)
  {
    dmap_nxsl_gp[p->first] += linv * (p->second);
    dmap_nxsl_gp[p->first] -= lllinv * sxsx * (p->second);
    dmap_nysl_gp[p->first] -= lllinv * sxsy * (p->second);
  }

  for (_CI p = dmap_nysl_gp_mod.begin(); p != dmap_nysl_gp_mod.end(); ++p)
  {
    dmap_nysl_gp[p->first] += linv * (p->second);
    dmap_nysl_gp[p->first] -= lllinv * sysy * (p->second);
    dmap_nxsl_gp[p->first] -= lllinv * sxsy * (p->second);
  }

  // *********************************************************************
  // finally compute Lin(XiGP_master)
  // *********************************************************************

  // add derivative of slave GP coordinates
  for (_CI p = dmap_xsl_gp.begin(); p != dmap_xsl_gp.end(); ++p)
    derivmxi[p->first] -= sgpn[1] * (p->second);
  for (_CI p = dmap_ysl_gp.begin(); p != dmap_ysl_gp.end(); ++p)
    derivmxi[p->first] += sgpn[0] * (p->second);

  // add derivatives of master node coordinates
  for (int i = 0; i < nummnode; ++i)
  {
    derivmxi[mmrtrnodes[i]->Dofs()[0]] += valmxigp[i] * sgpn[1];
    derivmxi[mmrtrnodes[i]->Dofs()[1]] -= valmxigp[i] * sgpn[0];
  }

  // add derivative of slave GP normal
  for (_CI p = dmap_nxsl_gp.begin(); p != dmap_nxsl_gp.end(); ++p)
    derivmxi[p->first] -= fac_ymsl_gp * (p->second);
  for (_CI p = dmap_nysl_gp.begin(); p != dmap_nysl_gp.end(); ++p)
    derivmxi[p->first] += fac_xmsl_gp * (p->second);

  // multiply all entries with cmxigp
  for (_CI p = derivmxi.begin(); p != derivmxi.end(); ++p)
    derivmxi[p->first] = cmxigp * (p->second);

  return;
}

/*----------------------------------------------------------------------*
 |  Compute directional derivative of XiGP master (3D)        popp 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::DerivXiGP3D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    const double* sxigp, const double* mxigp,
    const std::vector<GEN::pairedvector<int, double>>& derivsxi,
    std::vector<GEN::pairedvector<int, double>>& derivmxi, const double alpha)
{
  // check for problem dimension
  if (Dim() != 3) dserror("ERROR: 3D integration method called for non-3D problem");

  // we need the participating slave and master nodes
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();
  std::vector<MORTAR::MortarNode*> smrtrnodes(sele.NumNode());
  std::vector<MORTAR::MortarNode*> mmrtrnodes(mele.NumNode());
  const int numsnode = sele.NumNode();
  const int nummnode = mele.NumNode();

  for (int i = 0; i < numsnode; ++i)
  {
    smrtrnodes[i] = dynamic_cast<MORTAR::MortarNode*>(snodes[i]);
    if (!smrtrnodes[i]) dserror("ERROR: DerivXiGP3D: Null pointer!");
  }

  for (int i = 0; i < nummnode; ++i)
  {
    mmrtrnodes[i] = dynamic_cast<MORTAR::MortarNode*>(mnodes[i]);
    if (!mmrtrnodes[i]) dserror("ERROR: DerivXiGP3D: Null pointer!");
  }

  // we also need shape function derivs at the GP
  LINALG::SerialDenseVector valsxigp(numsnode);
  LINALG::SerialDenseVector valmxigp(nummnode);
  LINALG::SerialDenseMatrix derivsxigp(numsnode, 2, true);
  LINALG::SerialDenseMatrix derivmxigp(nummnode, 2, true);

  sele.EvaluateShape(sxigp, valsxigp, derivsxigp, numsnode);
  mele.EvaluateShape(mxigp, valmxigp, derivmxigp, nummnode);

  // we also need the GP slave coordinates + normal
  double sgpn[3] = {0.0, 0.0, 0.0};
  double sgpx[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < numsnode; ++i)
    for (int k = 0; k < 3; ++k)
    {
      sgpn[k] += valsxigp[i] * smrtrnodes[i]->MoData().n()[k];
      sgpx[k] += valsxigp[i] * smrtrnodes[i]->xspatial()[k];
    }

  // build 3x3 factor matrix L
  LINALG::Matrix<3, 3> lmatrix(true);
  for (int k = 0; k < 3; ++k) lmatrix(k, 2) = -sgpn[k];
  for (int z = 0; z < nummnode; ++z)
    for (int k = 0; k < 3; ++k)
    {
      lmatrix(k, 0) += derivmxigp(z, 0) * mmrtrnodes[z]->xspatial()[k];
      lmatrix(k, 1) += derivmxigp(z, 1) * mmrtrnodes[z]->xspatial()[k];
    }

  // get inverse of the 3x3 matrix L (in place)
  if (abs(lmatrix.Determinant()) < 1e-12) dserror("ERROR: Singular lmatrix for derivgp3d");

  lmatrix.Invert();

  // build directional derivative of slave GP normal
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  int linsize = 0;
  for (int i = 0; i < numsnode; ++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*>(snodes[i]);
    linsize += cnode->GetLinsize();
  }

  GEN::pairedvector<int, double> dmap_nxsl_gp(linsize);
  GEN::pairedvector<int, double> dmap_nysl_gp(linsize);
  GEN::pairedvector<int, double> dmap_nzsl_gp(linsize);

  for (int i = 0; i < numsnode; ++i)
  {
    GEN::pairedvector<int, double>& dmap_nxsl_i =
        static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
    GEN::pairedvector<int, double>& dmap_nysl_i =
        static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];
    GEN::pairedvector<int, double>& dmap_nzsl_i =
        static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[2];

    for (_CI p = dmap_nxsl_i.begin(); p != dmap_nxsl_i.end(); ++p)
      dmap_nxsl_gp[p->first] += valsxigp[i] * (p->second);
    for (_CI p = dmap_nysl_i.begin(); p != dmap_nysl_i.end(); ++p)
      dmap_nysl_gp[p->first] += valsxigp[i] * (p->second);
    for (_CI p = dmap_nzsl_i.begin(); p != dmap_nzsl_i.end(); ++p)
      dmap_nzsl_gp[p->first] += valsxigp[i] * (p->second);

    for (_CI p = derivsxi[0].begin(); p != derivsxi[0].end(); ++p)
    {
      double valx = derivsxigp(i, 0) * smrtrnodes[i]->MoData().n()[0];
      dmap_nxsl_gp[p->first] += valx * (p->second);
      double valy = derivsxigp(i, 0) * smrtrnodes[i]->MoData().n()[1];
      dmap_nysl_gp[p->first] += valy * (p->second);
      double valz = derivsxigp(i, 0) * smrtrnodes[i]->MoData().n()[2];
      dmap_nzsl_gp[p->first] += valz * (p->second);
    }

    for (_CI p = derivsxi[1].begin(); p != derivsxi[1].end(); ++p)
    {
      double valx = derivsxigp(i, 1) * smrtrnodes[i]->MoData().n()[0];
      dmap_nxsl_gp[p->first] += valx * (p->second);
      double valy = derivsxigp(i, 1) * smrtrnodes[i]->MoData().n()[1];
      dmap_nysl_gp[p->first] += valy * (p->second);
      double valz = derivsxigp(i, 1) * smrtrnodes[i]->MoData().n()[2];
      dmap_nzsl_gp[p->first] += valz * (p->second);
    }
  }

  // start to fill linearization maps for master GP
  // (1) all master nodes coordinates part
  for (int z = 0; z < nummnode; ++z)
  {
    for (int k = 0; k < 3; ++k)
    {
      derivmxi[0][mmrtrnodes[z]->Dofs()[k]] -= valmxigp[z] * lmatrix(0, k);
      derivmxi[1][mmrtrnodes[z]->Dofs()[k]] -= valmxigp[z] * lmatrix(1, k);
    }
  }

  // (2) slave Gauss point coordinates part
  for (int z = 0; z < numsnode; ++z)
  {
    for (int k = 0; k < 3; ++k)
    {
      derivmxi[0][smrtrnodes[z]->Dofs()[k]] += valsxigp[z] * lmatrix(0, k);
      derivmxi[1][smrtrnodes[z]->Dofs()[k]] += valsxigp[z] * lmatrix(1, k);

      for (_CI p = derivsxi[0].begin(); p != derivsxi[0].end(); ++p)
      {
        derivmxi[0][p->first] +=
            derivsxigp(z, 0) * smrtrnodes[z]->xspatial()[k] * lmatrix(0, k) * (p->second);
        derivmxi[1][p->first] +=
            derivsxigp(z, 0) * smrtrnodes[z]->xspatial()[k] * lmatrix(1, k) * (p->second);
      }

      for (_CI p = derivsxi[1].begin(); p != derivsxi[1].end(); ++p)
      {
        derivmxi[0][p->first] +=
            derivsxigp(z, 1) * smrtrnodes[z]->xspatial()[k] * lmatrix(0, k) * (p->second);
        derivmxi[1][p->first] +=
            derivsxigp(z, 1) * smrtrnodes[z]->xspatial()[k] * lmatrix(1, k) * (p->second);
      }
    }
  }

  // (3) slave Gauss point normal part
  for (_CI p = dmap_nxsl_gp.begin(); p != dmap_nxsl_gp.end(); ++p)
  {
    derivmxi[0][p->first] += alpha * lmatrix(0, 0) * (p->second);
    derivmxi[1][p->first] += alpha * lmatrix(1, 0) * (p->second);
  }
  for (_CI p = dmap_nysl_gp.begin(); p != dmap_nysl_gp.end(); ++p)
  {
    derivmxi[0][p->first] += alpha * lmatrix(0, 1) * (p->second);
    derivmxi[1][p->first] += alpha * lmatrix(1, 1) * (p->second);
  }
  for (_CI p = dmap_nzsl_gp.begin(); p != dmap_nzsl_gp.end(); ++p)
  {
    derivmxi[0][p->first] += alpha * lmatrix(0, 2) * (p->second);
    derivmxi[1][p->first] += alpha * lmatrix(1, 2) * (p->second);
  }

  /*
  // check linearization
  typedef std::map<int,double>::const_iterator CI;
  std::cout << "\nLinearization of current master GP:" << std::endl;
  std::cout << "-> Coordinate 1:" << std::endl;
  for (CI p=derivmxi[0].begin();p!=derivmxi[0].end();++p)
    std::cout << p->first << " " << p->second << std::endl;
  std::cout << "-> Coordinate 2:" << std::endl;
  for (CI p=derivmxi[1].begin();p!=derivmxi[1].end();++p)
      std::cout << p->first << " " << p->second << std::endl;
  */

  return;
}

/*----------------------------------------------------------------------*
 |  Compute deriv. of XiGP slave / master AuxPlane (3D)       popp 03/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::DerivXiGP3DAuxPlane(MORTAR::MortarElement& ele, double* xigp,
    double* auxn, std::vector<GEN::pairedvector<int, double>>& derivxi, double& alpha,
    std::vector<GEN::pairedvector<int, double>>& derivauxn,
    GEN::pairedvector<int, LINALG::Matrix<3, 1>>& derivgp)
{
  // check for problem dimension
  if (Dim() != 3) dserror("ERROR: 3D integration method called for non-3D problem");

  // we need the participating element nodes
  DRT::Node** nodes = ele.Nodes();
  std::vector<MORTAR::MortarNode*> mrtrnodes(ele.NumNode());
  const int numnode = ele.NumNode();

  for (int i = 0; i < numnode; ++i)
  {
    mrtrnodes[i] = dynamic_cast<MORTAR::MortarNode*>(nodes[i]);
    if (!mrtrnodes[i]) dserror("ERROR: DerivXiGP3DAuxPlane: Null pointer!");
  }

  // we also need shape function derivs at the GP
  LINALG::SerialDenseVector valxigp(numnode);
  LINALG::SerialDenseMatrix derivxigp(numnode, 2, true);
  ele.EvaluateShape(xigp, valxigp, derivxigp, numnode);

  // build 3x3 factor matrix L
  LINALG::Matrix<3, 3> lmatrix(true);
  for (int k = 0; k < 3; ++k) lmatrix(k, 2) = -auxn[k];
  for (int z = 0; z < numnode; ++z)
    for (int k = 0; k < 3; ++k)
    {
      lmatrix(k, 0) += derivxigp(z, 0) * mrtrnodes[z]->xspatial()[k];
      lmatrix(k, 1) += derivxigp(z, 1) * mrtrnodes[z]->xspatial()[k];
    }

  // get inverse of the 3x3 matrix L (in place)
  lmatrix.Invert();

  // start to fill linearization maps for element GP
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // see if this is an IntEle
  MORTAR::IntElement* ie = dynamic_cast<MORTAR::IntElement*>(&ele);

  // (1) all nodes coordinates part
  if (!ie)
    for (int z = 0; z < numnode; ++z)
      for (int k = 0; k < 3; ++k)
      {
        derivxi[0][mrtrnodes[z]->Dofs()[k]] -= valxigp[z] * lmatrix(0, k);
        derivxi[1][mrtrnodes[z]->Dofs()[k]] -= valxigp[z] * lmatrix(1, k);
      }
  else
  {
    std::vector<std::vector<GEN::pairedvector<int, double>>> nodelin(0);
    ie->NodeLinearization(nodelin);
    for (int z = 0; z < numnode; ++z)
      for (int k = 0; k < 3; ++k)
        for (_CI p = nodelin[z][k].begin(); p != nodelin[z][k].end(); ++p)
        {
          derivxi[0][p->first] -= valxigp[z] * lmatrix(0, k) * p->second;
          derivxi[1][p->first] -= valxigp[z] * lmatrix(1, k) * p->second;
        }
  }

  // (2) Gauss point coordinates part
  for (GEN::pairedvector<int, LINALG::Matrix<3, 1>>::const_iterator p = derivgp.begin();
       p != derivgp.end(); ++p)
  {
    const int pf = p->first;
    double& derivxi0 = derivxi[0][pf];
    double& derivxi1 = derivxi[1][pf];
    const LINALG::Matrix<3, 1>& tmp = p->second;
    for (int d = 0; d < 3; ++d)
    {
      derivxi0 += lmatrix(0, d) * tmp(d);
      derivxi1 += lmatrix(1, d) * tmp(d);
    }
  }

  // (3) AuxPlane normal part
  for (_CI p = derivauxn[0].begin(); p != derivauxn[0].end(); ++p)
  {
    derivxi[0][p->first] += alpha * lmatrix(0, 0) * (p->second);
    derivxi[1][p->first] += alpha * lmatrix(1, 0) * (p->second);
  }
  for (_CI p = derivauxn[1].begin(); p != derivauxn[1].end(); ++p)
  {
    derivxi[0][p->first] += alpha * lmatrix(0, 1) * (p->second);
    derivxi[1][p->first] += alpha * lmatrix(1, 1) * (p->second);
  }
  for (_CI p = derivauxn[2].begin(); p != derivauxn[2].end(); ++p)
  {
    derivxi[0][p->first] += alpha * lmatrix(0, 2) * (p->second);
    derivxi[1][p->first] += alpha * lmatrix(1, 2) * (p->second);
  }

  /*
  // check linearization
  typedef std::map<int,double>::const_iterator CI;
  std::cout << "\nLinearization of current slave / master GP:" << std::endl;
  std::cout << "-> Coordinate 1:" << std::endl;
  for (CI p=derivxi[0].begin();p!=derivxi[0].end();++p)
    std::cout << p->first << " " << p->second << std::endl;
  std::cout << "-> Coordinate 2:" << std::endl;
  for (CI p=derivxi[1].begin();p!=derivxi[1].end();++p)
    std::cout << p->first << " " << p->second << std::endl;
  std::cout << "-> Coordinate 3:" << std::endl;
  for (CI p=derivxi[2].begin();p!=derivxi[2].end();++p)
    std::cout << p->first << " " << p->second << std::endl;
  */

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateGP_3D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseVector& mval, LINALG::SerialDenseMatrix& sderiv,
    LINALG::SerialDenseMatrix& mderiv, LINALG::SerialDenseMatrix& lmderiv,
    GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap, double& wgt, double& jac,
    GEN::pairedvector<int, double>& derivjac, double* normal,
    std::vector<GEN::pairedvector<int, double>>& dnmap_unit, double& gap,
    GEN::pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
    std::vector<GEN::pairedvector<int, double>>& derivsxi,
    std::vector<GEN::pairedvector<int, double>>& derivmxi)
{
  // different algorithms integrate different sh*t
  switch (algo_)
  {
    case INPAR::MORTAR::algorithm_mortar:
    {
      // is quadratic case?
      bool quad = sele.IsQuad();

      // weighted gap
      GP_3D_wGap(sele, sval, lmval, &gap, jac, wgt, quad);
      for (int j = 0; j < sele.NumNode(); ++j)
        GP_G_Lin(j, sele, mele, sval, mval, lmval, sderiv, lmderiv, gap, normal, jac, wgt,
            deriv_gap, derivjac, derivsxi, dualmap);

      // check bound
      bool bound = false;
      for (int i = 0; i < sele.NumNode(); ++i)
        if (dynamic_cast<CONTACT::CoNode*>(sele.Nodes()[i])->IsOnBoundorCE()) bound = true;

      // integrate D and M matrix
      GP_DM(sele, mele, lmval, sval, mval, jac, wgt, bound);

      // compute segment D/M linearization
      if (bound)
      {
        GP_3D_DM_Lin_bound(sele, mele, sval, mval, lmval, sderiv, lmderiv, mderiv, jac, wgt,
            derivjac, derivsxi, derivmxi, dualmap);
      }
      else
      {
        GP_3D_DM_Lin(sele, mele, sval, mval, lmval, sderiv, mderiv, lmderiv, wgt, jac, derivsxi,
            derivmxi, derivjac, dualmap);
      }

      // Creating the WEIGHTED tangential relative slip increment (non-objective)
      if (gpslip_)
      {
        int linsize = 0;
        for (int i = 0; i < sele.NumNode(); ++i)
          linsize += dynamic_cast<CoNode*>(sele.Nodes()[i])->GetLinsize();

        double jumpvalv[2] = {0.0, 0.0};  // jump for slipincr --> equal to jumpval
        std::vector<GEN::pairedvector<int, double>> dslipgp(
            2, ((mele.NumNode() * Dim()) + linsize));  // deriv. of slip for slipincr (xi, eta)

        // Creating the WEIGHTED tangential relative slip increment (non-objective)
        GP_3D_SlipIncr(sele, mele, sval, mval, lmval, sderiv, mderiv, jac, wgt, jumpvalv, derivsxi,
            derivmxi, dslipgp);
        // Lin weighted slip
        for (int iter = 0; iter < sele.NumNode(); ++iter)
          GP_3D_SlipIncr_Lin(iter, sele, sval, lmval, sderiv, lmderiv, jac, wgt, jumpvalv, derivjac,
              dslipgp, derivsxi, dualmap);
      }

      //*******************************
      // WEAR stuff
      //*******************************
      // std. wear for all wear-algorithm types
      if (wearlaw_ != INPAR::WEAR::wear_none)
      {
        int linsize = 0;
        for (int i = 0; i < sele.NumNode(); ++i)
          linsize += dynamic_cast<CoNode*>(sele.Nodes()[i])->GetLinsize();

        double jumpval[2] = {0.0, 0.0};  // jump for wear
        double wearval[1] = {0.0};       // wear value
        GEN::pairedvector<int, double> dsliptmatrixgp(
            (mele.NumNode() * Dim()) + linsize);  // deriv. of slip for wear
        GEN::pairedvector<int, double> dweargp(
            (mele.NumNode() * Dim()) + linsize);  // wear lin. without lm and jac.
        Teuchos::RCP<LINALG::SerialDenseMatrix> lagmult =
            Teuchos::rcp(new LINALG::SerialDenseMatrix(3, sele.NumNode()));
        sele.GetNodalLagMult(*lagmult);

        GP_3D_Wear(sele, mele, sval, sderiv, mval, mderiv, lmval, lmderiv, lagmult, normal, jac,
            wgt, jumpval, wearval, dsliptmatrixgp, dweargp, derivsxi, derivmxi, dnmap_unit,
            dualmap);

        if (WearType() == INPAR::WEAR::wear_intstate and wearimpl_ == true)
          for (int j = 0; j < sele.NumNode(); ++j)
            GP_3D_Wear_Lin(j, sele, sval, lmval, sderiv, lmderiv, jac, normal, wgt, *wearval,
                jumpval, dweargp, derivjac, derivsxi, dualmap);

        // integrate T and E matrix for discr. wear
        if (WearType() == INPAR::WEAR::wear_primvar) GP_TE(sele, lmval, sval, jac, wgt, jumpval);

        // both-sided discr wear specific stuff
        if (WearSide() == INPAR::WEAR::wear_both and WearType() == INPAR::WEAR::wear_primvar)
        {
          LINALG::SerialDenseVector lm2val(mele.NumNode());
          LINALG::SerialDenseMatrix lm2deriv(mele.NumNode(), Dim() - 1, true);
          mele.EvaluateShapeLagMult(ShapeFcn(), mxi, lm2val, lm2deriv, mele.NumNode());

          GP_TE_Master(sele, mele, lmval, lm2val, mval, jac, wgt, jumpval, Comm_);
        }

        // both-sided wear specific stuff
        if (WearSide() == INPAR::WEAR::wear_both and WearType() == INPAR::WEAR::wear_intstate)
        {
          LINALG::SerialDenseVector lm2val(mele.NumNode());
          LINALG::SerialDenseMatrix lm2deriv(mele.NumNode(), 2, true);
          mele.EvaluateShapeLagMult(ShapeFcn(), mxi, lm2val, lm2deriv, mele.NumNode());

          GP_D2(sele, mele, lm2val, mval, jac, wgt, Comm_);
        }

        // Lin wear matrices T and E for discr. wear
        if (wearimpl_ == true and WearType() == INPAR::WEAR::wear_primvar)
          for (int j = 0; j < sele.NumNode(); ++j)
            GP_3D_TE_Lin(j, sele, sval, lmval, sderiv, lmderiv, jac, wgt, jumpval, derivsxi,
                derivjac, dsliptmatrixgp, dualmap);



        if (WearSide() == INPAR::WEAR::wear_both and WearType() == INPAR::WEAR::wear_primvar and
            wearimpl_ == true)
        {
          LINALG::SerialDenseVector lm2val(mele.NumNode());
          LINALG::SerialDenseMatrix lm2deriv(mele.NumNode(), Dim() - 1, true);
          GEN::pairedvector<int, Epetra_SerialDenseMatrix> dual2map(
              (sele.NumNode() + mele.NumNode()) * Dim(), 0,
              Epetra_SerialDenseMatrix(mele.NumNode(), mele.NumNode()));
          mele.EvaluateShapeLagMult(ShapeFcn(), mxi, lm2val, lm2deriv, mele.NumNode());
          mele.DerivShapeDual(dual2map);

          for (int iter = 0; iter < mele.NumNode(); ++iter)
            GP_3D_TE_Master_Lin(iter, sele, mele, sval, mval, lmval, lm2val, sderiv, mderiv,
                lmderiv, lm2deriv, jac, wgt, jumpval, derivsxi, derivmxi, derivjac, dsliptmatrixgp,
                dualmap, dual2map, Comm_);
        }
      }

      //*******************************
      // PORO stuff
      //*******************************
      bool poroprob = false;
      if (imortar_.get<int>("PROBTYPE") == INPAR::CONTACT::poro)
        if (imortar_.get<bool>("CONTACTNOPEN"))  // evaluate additional terms just in case of no
                                                 // penectration condition
          poroprob = true;
      if (poroprob)
      {
        double ncoup[1] = {0.0};
        std::map<int, double> dncoupgp;      // ncoup lin. without lm and jac.
        std::map<int, double> dvelncoupgp;   // velocity ncoup lin. without lm and jac.
        std::map<int, double> dpresncoupgp;  // pressure ncoup lin. without lm and jac.

        GP_NCOUP_DERIV(sele, mele, sval, mval, lmval, sderiv, mderiv, ncoup, normal, jac, wgt, sxi,
            derivsxi, derivmxi, dncoupgp, dvelncoupgp, dpresncoupgp, dnmap_unit, false);
        // Lin ncoup condition
        for (int j = 0; j < sele.NumNode(); ++j)
          GP_NCOUP_LIN(j, sele, mele, sval, mval, lmval, sderiv, lmderiv, *ncoup, normal, jac, wgt,
              dncoupgp, dvelncoupgp, dpresncoupgp, derivjac, derivsxi, derivmxi, dualmap);
      }
      break;
    }

    default:
      dserror("contact integrator doesn't know what to do with this algorithm");
  };

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateGP_2D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseVector& mval, LINALG::SerialDenseMatrix& sderiv,
    LINALG::SerialDenseMatrix& mderiv, LINALG::SerialDenseMatrix& lmderiv,
    GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap, double& wgt, double& jac,
    GEN::pairedvector<int, double>& derivjac, double* normal,
    std::vector<GEN::pairedvector<int, double>>& dnmap_unit, double& gap,
    GEN::pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
    std::vector<GEN::pairedvector<int, double>>& derivsxi,
    std::vector<GEN::pairedvector<int, double>>& derivmxi)
{
  // different algorithms integrate different sh*t
  switch (algo_)
  {
    case INPAR::MORTAR::algorithm_mortar:
    {
      // decide whether boundary modification has to be considered or not
      // this is element-specific (is there a boundary node in this element?)
      bool bound = false;
      for (int k = 0; k < sele.NumNode(); ++k)
      {
        MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(sele.Nodes()[k]);
        if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");

        if (mymrtrnode->IsOnBoundorCE())
        {
          bound = true;
          break;
        }
      }

      // decide whether linear LM are used for quadratic FE here
      bool linlm = false;
      if (LagMultQuad() == INPAR::MORTAR::lagmult_lin && sele.Shape() == DRT::Element::line3)
      {
        bound = false;  // crosspoints and linear LM NOT at the same time!!!!
        linlm = true;
      }

      // weighted gap
      GP_2D_wGap(sele, sval, lmval, &gap, jac, wgt);
      for (int j = 0; j < sval.M(); ++j)
        GP_2D_G_Lin(j, sele, mele, sval, mval, lmval, sderiv, lmderiv, gap, normal, jac, wgt,
            deriv_gap, derivjac, derivsxi, dualmap);

      // integrate D and M matrix
      GP_DM(sele, mele, lmval, sval, mval, jac, wgt, bound);

      // compute segment D/M linearization  -- bound
      if (bound)
      {
        GP_2D_DM_Lin_bound(sele, mele, sval, mval, lmval, sderiv, mderiv, lmderiv, jac, wgt,
            derivjac, derivsxi, derivmxi, dualmap);
      }
      else
      {
        // D/M linearization
        for (int j = 0; j < sval.M(); ++j)
        {
          GP_2D_DM_Lin(j, bound, linlm, sele, mele, sval, mval, lmval, sderiv, mderiv, lmderiv, jac,
              wgt, derivsxi, derivmxi, derivjac, dualmap);
        }
      }

      // Creating the tangential relative slip increment (non-objective)
      if (DRT::INPUT::IntegralValue<int>(imortar_, "GP_SLIP_INCR") == true)
      {
        int linsize = 0;
        for (int i = 0; i < sele.NumNode(); ++i)
          linsize += dynamic_cast<CoNode*>(sele.Nodes()[i])->GetLinsize();
        linsize = linsize * 2;

        double jumpvalv[1] = {0.0};  // jump for slipincr --> equal to jumpval
        GEN::pairedvector<int, double> dslipgp(
            linsize + Dim() * mele.NumNode());  // deriv. of slip for slipincr

        GP_2D_SlipIncr(sele, mele, sval, mval, lmval, sderiv, mderiv, jac, wgt, jumpvalv, derivsxi,
            derivmxi, dslipgp, linsize);

        for (int iter = 0; iter < sele.NumNode(); ++iter)
          GP_2D_SlipIncr_Lin(iter, sele, sval, lmval, sderiv, lmderiv, jac, wgt, jumpvalv, derivsxi,
              dslipgp, derivjac, dualmap);
      }

      // wear specific stuff
      if (wearlaw_ != INPAR::WEAR::wear_none)
      {
        // nodal Lagrange mulitplier
        Teuchos::RCP<LINALG::SerialDenseMatrix> lagmult =
            Teuchos::rcp(new LINALG::SerialDenseMatrix(3, sele.NumNode()));
        sele.GetNodalLagMult(*lagmult);

        int linsize = 0;
        for (int i = 0; i < sele.NumNode(); ++i)
          linsize += dynamic_cast<CoNode*>(sele.Nodes()[i])->GetLinsize();
        linsize = linsize * 2;

        double jumpval[1] = {0.0};  // jump for wear
        double wearval[1] = {0.0};  // wear value
        GEN::pairedvector<int, double> dsliptmatrixgp(
            linsize + Dim() * mele.NumNode());  // deriv. of slip for wear
        GEN::pairedvector<int, double> dweargp(
            linsize + Dim() * mele.NumNode());  // wear lin without weighting and jac

        // std. wear for all wear-algorithm types
        GP_2D_Wear(sele, mele, sval, sderiv, mval, mderiv, lmval, lmderiv, lagmult, normal, jac,
            wgt, jumpval, wearval, dsliptmatrixgp, dweargp, derivsxi, derivmxi, dnmap_unit,
            dualmap);

        // integrate T and E matrix for discr. wear
        if (WearType() == INPAR::WEAR::wear_primvar) GP_TE(sele, lmval, sval, jac, wgt, jumpval);

        // both-sided discr wear specific stuff
        if (WearSide() == INPAR::WEAR::wear_both and WearType() == INPAR::WEAR::wear_primvar)
        {
          LINALG::SerialDenseVector lm2val(mele.NumNode());
          LINALG::SerialDenseMatrix lm2deriv(mele.NumNode(), Dim() - 1, true);
          mele.EvaluateShapeLagMult(ShapeFcn(), mxi, lm2val, lm2deriv, mele.NumNode());

          GP_TE_Master(sele, mele, lmval, lm2val, mval, jac, wgt, jumpval, Comm_);
        }

        // both-sided map wear specific stuff
        if (WearSide() == INPAR::WEAR::wear_both and WearType() == INPAR::WEAR::wear_intstate)
        {
          LINALG::SerialDenseVector lm2val(mele.NumNode());
          LINALG::SerialDenseMatrix lm2deriv(mele.NumNode(), 2, true);
          mele.EvaluateShapeLagMult(ShapeFcn(), mxi, lm2val, lm2deriv, mele.NumNode());

          GP_D2(sele, mele, lm2val, mval, jac, wgt, Comm_);
        }

        // Lin wear for impl. alg.
        if (wearimpl_ == true and WearType() == INPAR::WEAR::wear_intstate)
          for (int j = 0; j < sele.NumNode(); ++j)
            GP_2D_Wear_Lin(j, sele, sval, lmval, sderiv, lmderiv, jac, normal, wgt, wearval[0],
                jumpval, dweargp, derivjac, derivsxi, dualmap);

        // Lin wear T and E matrix
        if (wearimpl_ == true and WearType() == INPAR::WEAR::wear_primvar)
          for (int j = 0; j < sele.NumNode(); ++j)
            GP_2D_TE_Lin(j, sele, sval, lmval, sderiv, lmderiv, jac, wgt, jumpval, derivsxi,
                derivjac, dsliptmatrixgp, dualmap);

        if (WearSide() == INPAR::WEAR::wear_both and WearType() == INPAR::WEAR::wear_primvar and
            wearimpl_ == true)
        {
          LINALG::SerialDenseVector lm2val(mele.NumNode());
          LINALG::SerialDenseMatrix lm2deriv(mele.NumNode(), Dim() - 1, true);
          GEN::pairedvector<int, Epetra_SerialDenseMatrix> dual2map(
              (sele.NumNode() + mele.NumNode()) * Dim(), 0,
              Epetra_SerialDenseMatrix(mele.NumNode(), mele.NumNode()));
          mele.EvaluateShapeLagMult(ShapeFcn(), mxi, lm2val, lm2deriv, mele.NumNode());
          mele.DerivShapeDual(dual2map);

          for (int iter = 0; iter < mele.NumNode(); ++iter)
            GP_3D_TE_Master_Lin(iter, sele, mele, sval, mval, lmval, lm2val, sderiv, mderiv,
                lmderiv, lm2deriv, jac, wgt, jumpval, derivsxi, derivmxi, derivjac, dsliptmatrixgp,
                dualmap, dual2map, Comm_);
        }
      }  // if wear

      //*******************************
      // PORO stuff
      //*******************************
      bool poroprob = false;
      if (imortar_.get<int>("PROBTYPE") == INPAR::CONTACT::poro)
        if (imortar_.get<bool>("CONTACTNOPEN"))  // evaluate additional terms just in case of no
                                                 // penectration condition
          poroprob = true;
      if (poroprob)
      {
        double ncoup[1] = {0.0};
        std::map<int, double> dncoupgp;      // ncoup lin. without lm and jac.
        std::map<int, double> dvelncoupgp;   // velocity ncoup lin. without lm and jac.
        std::map<int, double> dpresncoupgp;  // pressure ncoup lin. without lm and jac.

        GP_NCOUP_DERIV(sele, mele, sval, mval, lmval, sderiv, mderiv, ncoup, normal, jac, wgt, sxi,
            derivsxi, derivmxi, dncoupgp, dvelncoupgp, dpresncoupgp, dnmap_unit, false);
        // Lin ncoup condition
        for (int j = 0; j < sele.NumNode(); ++j)
          GP_NCOUP_LIN(j, sele, mele, sval, mval, lmval, sderiv, lmderiv, *ncoup, normal, jac, wgt,
              dncoupgp, dvelncoupgp, dpresncoupgp, derivjac, derivsxi, derivmxi, dualmap);
      }

      break;
    }
    default:
      dserror("contact integrator doesn't know what to do with this algorithm");
  };

  return;
}


/*----------------------------------------------------------------------*
 |  Compute entries for D and M matrix at GP                 farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::GP_DM(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    LINALG::SerialDenseVector& lmval, LINALG::SerialDenseVector& sval,
    LINALG::SerialDenseVector& mval, double& jac, double& wgt, bool& bound)
{
  int nrow;
  int ncol;
  // get slave element nodes themselves
  DRT::Node** snodes = NULL;
  DRT::Node** mnodes = NULL;
  nrow = sele.NumNode();
  snodes = sele.Nodes();
  ncol = mele.NumNode();
  mnodes = mele.Nodes();

  // BOUNDARY NODE MODIFICATION **********************************
  // We have modified their neighbors' dual shape functions, so we
  // now have a problem with off-diagonal entries occurring in D.
  // Of course we want to keep the diagonality property of the D
  // matrix, but still we may not modify the whole Mortar coupling
  // setting! We achieve both by applying a quite simple but very
  // effective trick: The boundary nodes have already been defined
  // as being master nodes, so all we have to do here, is to shift
  // the off-diagonal terms from D to the respective place in M,
  // which is not diagonal anyway! (Mind the MINUS sign!!!)
  // *************************************************************

  // compute segment D/M matrix ****************************************
  // dual shape functions without locally linear Lagrange multipliers
  if ((ShapeFcn() == INPAR::MORTAR::shape_dual ||
          ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) &&
      LagMultQuad() != INPAR::MORTAR::lagmult_lin)
  {
    for (int j = 0; j < nrow; ++j)
    {
      CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(snodes[j]);

      if (cnode->IsOnBoundorCE()) continue;

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mnodes[k]);

        // multiply the two shape functions
        double prod = lmval[j] * mval[k] * jac * wgt;
        if (abs(prod) > MORTARINTTOL)
        {
          cnode->AddMValue(mnode->Id(), prod);
          cnode->AddMNode(mnode->Id());  // only for friction!
          if (!bound)
          {
            cnode->AddDValue(cnode->Id(), prod);
            cnode->AddSNode(cnode->Id());  // only for friction!
          }
        }
      }

      // integrate dseg (boundary modification)
      if (bound)
      {
        //        bool j_boundnode = cnode->IsOnBoundorCE();

        for (int k = 0; k < nrow; ++k)
        {
          CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(snodes[k]);

          //          bool k_boundnode = mnode->IsOnBoundorCE();

          // do not assemble off-diagonal terms if j,k are both non-boundary nodes
          //          if (!j_boundnode && !k_boundnode && (j!=k))
          //            continue;

          // multiply the two shape functions
          double prod = lmval[j] * sval[k] * jac * wgt;
          if (abs(prod) > MORTARINTTOL)
          {
            // for the "bound" flag, the node on the slave side does not carry a LM value
            // therefore, we assemble the negative value to the M-matrix (which contains
            // all nodes, that do not have a LM).
            if (mnode->IsOnBound())
            {
              cnode->AddMValue(mnode->Id(), -prod);
              cnode->AddMNode(mnode->Id());  // only for friction!
            }
            // in the non-smooth case, the nodes at corners/edges will have a LM value,
            // however the corresponding integration is not performed here (here we only
            // do the surface integration). Since the nodes will have a LM value, the corresponding
            // integrals should be assembled in the D-matrix.
            else if (mnode->IsOnBoundorCE())
            {
              // isolate the dseg entries to be filled
              // (both the main diagonal and every other secondary diagonal)
              // and add current Gauss point's contribution to dseg
              // NONSMOOTH MODIFICATION
              cnode->AddDValue(mnode->Id(), prod);
              cnode->AddSNode(mnode->Id());  // only for friction!
            }
            else
            {
              cnode->AddDValue(cnode->Id(), prod);
              cnode->AddSNode(cnode->Id());  // only for friction!
            }
          }
        }
      }
    }
  }
  // standard shape functions or dual shape functions with locally linear Lagrange multipliers
  else
  {
    for (int j = 0; j < nrow; ++j)
    {
      CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(snodes[j]);

      if ((shapefcn_ == INPAR::MORTAR::shape_standard && cnode->IsOnBoundorCE()) ||
          ((shapefcn_ == INPAR::MORTAR::shape_dual ||
               shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin) &&
              cnode->IsOnBound()))
        continue;

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mnodes[k]);

        // multiply the two shape functions
        double prod = lmval[j] * mval[k] * jac * wgt;
        if (abs(prod) > MORTARINTTOL)
        {
          cnode->AddMValue(mnode->Id(), prod);
          cnode->AddMNode(mnode->Id());  // only for friction!
        }
      }

      // integrate dseg
      for (int k = 0; k < nrow; ++k)
      {
        CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(snodes[k]);

        // multiply the two shape functions
        double prod = lmval[j] * sval[k] * jac * wgt;
        if (abs(prod) > MORTARINTTOL)
        {
          if ((ShapeFcn() == INPAR::MORTAR::shape_standard && snode->IsOnBoundorCE()) ||
              ((ShapeFcn() == INPAR::MORTAR::shape_dual ||
                   ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) &&
                  snode->IsOnBound()))
          {
            if (ShapeFcn() == INPAR::MORTAR::shape_standard && nonsmooth_)
            {
              if (snode->IsOnBound() and !snode->IsOnCornerEdge())
              {
                cnode->AddMValue(snode->Id(), -prod);
                cnode->AddMNode(snode->Id());  // only for friction!
              }
              else
              {
                cnode->AddDValue(snode->Id(), prod);
                cnode->AddSNode(snode->Id());  // only for friction!
              }
            }
            // this still needs to be checked for dual shape functions with locally linear Lagrange
            // multipliers
            else if ((ShapeFcn() == INPAR::MORTAR::shape_dual ||
                         ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) &&
                     LagMultQuad() == INPAR::MORTAR::lagmult_lin && nonsmooth_)
              dserror(
                  "ERROR: dual shape functions with locally linear Lagrange multipliers in "
                  "combination non-smooth approach still needs to be checked!");
            else
            {
              cnode->AddMValue(snode->Id(), -prod);
              cnode->AddMNode(snode->Id());  // only for friction!
            }
          }
          else
          {
            cnode->AddDValue(snode->Id(), prod);
            cnode->AddSNode(snode->Id());  // only for friction!
          }
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Compute entries for D and M matrix at GP (3D Quad)       farah 12/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_3D_DM_Quad(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, MORTAR::IntElement& sintele, LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseVector& lmintval, LINALG::SerialDenseVector& sval,
    LINALG::SerialDenseVector& mval, const double& jac, double& wgt, const int& nrow,
    const int& nintrow, const int& ncol, const int& ndof, bool& bound)
{
  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: Null pointer!");
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes) dserror("ERROR: Null pointer!");
  DRT::Node** sintnodes = sintele.Nodes();
  if (!sintnodes) dserror("ERROR: Null pointer for sintnodes!");

  // CASE 1/2: Standard LM shape functions and quadratic, linear or constant interpolation
  // CASE 5: dual LM shape functions and linear interpolation
  if ((ShapeFcn() == INPAR::MORTAR::shape_standard &&
          (LagMultQuad() == INPAR::MORTAR::lagmult_quad ||
              LagMultQuad() == INPAR::MORTAR::lagmult_lin ||
              LagMultQuad() == INPAR::MORTAR::lagmult_const)) ||
      ((ShapeFcn() == INPAR::MORTAR::shape_dual ||
           ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) &&
          LagMultQuad() == INPAR::MORTAR::lagmult_lin))
  {
    // compute all mseg and dseg matrix entries
    // loop over Lagrange multiplier dofs j
    for (int j = 0; j < nrow; ++j)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(snodes[j]);

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        CoNode* mnode = dynamic_cast<CoNode*>(mnodes[k]);

        // multiply the two shape functions
        double prod = lmval[j] * mval[k] * jac * wgt;
        if (abs(prod) > MORTARINTTOL)
        {
          cnode->AddMValue(mnode->Id(), prod);
          cnode->AddMNode(mnode->Id());
        }
      }

      // integrate dseg
      for (int k = 0; k < nrow; ++k)
      {
        CoNode* snode = dynamic_cast<CoNode*>(snodes[k]);

        // multiply the two shape functions
        double prod = lmval[j] * sval[k] * jac * wgt;
        if (abs(prod) > MORTARINTTOL)
        {
          if (snode->IsOnBound())
          {
            cnode->AddMValue(snode->Id(), -prod);
            cnode->AddMNode(snode->Id());
          }
          else
          {
            cnode->AddDValue(snode->Id(), prod);
            cnode->AddSNode(snode->Id());
          }
        }
      }
    }
  }
  // CASE 3: Standard LM shape functions and piecewise linear interpolation
  else if (ShapeFcn() == INPAR::MORTAR::shape_standard &&
           LagMultQuad() == INPAR::MORTAR::lagmult_pwlin)
  {
    // compute all mseg and dseg matrix entries
    // loop over Lagrange multiplier dofs j
    for (int j = 0; j < nintrow; ++j)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(sintnodes[j]);

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        CoNode* mnode = dynamic_cast<CoNode*>(mnodes[k]);

        // multiply the two shape functions
        double prod = lmintval[j] * mval[k] * jac * wgt;
        if (abs(prod) > MORTARINTTOL)
        {
          cnode->AddMValue(mnode->Id(), prod);
          cnode->AddMNode(mnode->Id());
        }
      }

      // integrate dseg
      for (int k = 0; k < nrow; ++k)
      {
        CoNode* snode = dynamic_cast<CoNode*>(snodes[k]);

        // multiply the two shape functions
        double prod = lmintval[j] * sval[k] * jac * wgt;
        if (abs(prod) > MORTARINTTOL)
        {
          if (snode->IsOnBound())
          {
            cnode->AddMValue(snode->Id(), -prod);
            cnode->AddMNode(snode->Id());
          }
          else
          {
            cnode->AddDValue(snode->Id(), prod);
            cnode->AddSNode(snode->Id());
          }
        }
      }
    }
  }

  // CASE 4: Dual LM shape functions and quadratic or constant interpolation
  else if ((ShapeFcn() == INPAR::MORTAR::shape_dual ||
               ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) &&
           (LagMultQuad() == INPAR::MORTAR::lagmult_quad ||
               LagMultQuad() == INPAR::MORTAR::lagmult_const))
  {
    // compute all mseg and dseg matrix entries
    // loop over Lagrange multiplier dofs j
    for (int j = 0; j < nrow; ++j)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(snodes[j]);

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        CoNode* mnode = dynamic_cast<CoNode*>(mnodes[k]);

        // multiply the two shape functions
        double prod = lmval[j] * mval[k] * jac * wgt;
        if (abs(prod) > MORTARINTTOL)
        {
          cnode->AddDValue(cnode->Id(), prod);
          cnode->AddSNode(cnode->Id());
          cnode->AddMValue(mnode->Id(), prod);
          cnode->AddMNode(mnode->Id());
        }
      }
    }
  }
  // INVALID CASES
  else
  {
    dserror("ERROR: Invalid integration case for 3D quadratic mortar!");
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Compute entries for weighted Gap at GP                   farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_2D_wGap(MORTAR::MortarElement& sele,
    LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& lmval, double* gap, double& jac,
    double& wgt)
{
  DRT::Node** mynodes = NULL;
  int nrow = 0;
  nrow = sele.NumNode();
  mynodes = sele.Nodes();

  // **************************
  // add to node
  // **************************
  for (int j = 0; j < nrow; ++j)
  {
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(mynodes[j]);

    double prod = 0.0;
    // Petrov-Galerkin approach (dual LM for D/M but standard LM for gap)
    if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) prod = sval[j] * gap[0] * jac * wgt;
    // usual standard or dual LM approach
    else
      prod = lmval[j] * gap[0] * jac * wgt;


    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (cnode->IsOnBoundorCE()) continue;

    // add current Gauss point's contribution to gseg
    cnode->AddgValue(prod);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Compute entries for weighted Gap at GP                   farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::GP_3D_wGap(MORTAR::MortarElement& sele, LINALG::SerialDenseVector& sval,
    LINALG::SerialDenseVector& lmval, double* gap, double& jac, double& wgt, bool quadratic,
    int nintrow)
{
  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // number of nodes (slave, master)
  const int nrow = sele.NumNode();

  // **************************
  // add to node
  // **************************
  if (!quadratic)
  {
    for (int j = 0; j < nrow; ++j)
    {
      CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(snodes[j]);

      double prod = 0.0;
      // Petrov-Galerkin approach (dual LM for D/M but standard LM for gap)
      if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) prod = sval[j] * gap[0] * jac * wgt;
      // usual standard or dual LM approach
      else
        prod = lmval[j] * gap[0] * jac * wgt;

      // do not process slave side boundary nodes
      // (their row entries would be zero anyway!)
      if (cnode->IsOnBoundorCE()) continue;

      // add current Gauss point's contribution to gseg
      cnode->AddgValue(prod);
    }
  }
  else
  {
    // compute cell gap vector *******************************************
    // CASE 1/2: Standard LM shape functions and quadratic or linear interpolation
    if (ShapeFcn() == INPAR::MORTAR::shape_standard &&
        (LagMultQuad() == INPAR::MORTAR::lagmult_quad ||
            LagMultQuad() == INPAR::MORTAR::lagmult_lin ||
            LagMultQuad() == INPAR::MORTAR::lagmult_const))
    {
      for (int j = 0; j < nrow; ++j)
      {
        CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(snodes[j]);

        double prod = 0.0;
        prod = lmval[j] * gap[0] * jac * wgt;

        if (cnode->IsOnBoundorCE()) continue;

        // add current Gauss point's contribution to gseg
        cnode->AddgValue(prod);
      }
    }

    // CASE 3: Standard LM shape functions and piecewise linear interpolation
    // Attention:  for this case, lmval represents lmintval !!!
    else if (ShapeFcn() == INPAR::MORTAR::shape_standard &&
             LagMultQuad() == INPAR::MORTAR::lagmult_pwlin)
    {
      if (nintrow == 0) dserror("ERROR!");

      for (int j = 0; j < nintrow; ++j)
      {
        CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(snodes[j]);

        double prod = 0.0;
        prod = lmval[j] * gap[0] * jac * wgt;

        if (cnode->IsOnBound()) continue;

        // add current Gauss point's contribution to gseg
        cnode->AddgValue(prod);
      }
    }

    // CASE 4: Dual LM shape functions and quadratic, linear or constant interpolation
    else if ((ShapeFcn() == INPAR::MORTAR::shape_dual ||
                 ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) &&
             (LagMultQuad() == INPAR::MORTAR::lagmult_quad ||
                 LagMultQuad() == INPAR::MORTAR::lagmult_lin ||
                 LagMultQuad() == INPAR::MORTAR::lagmult_const))
    {
      for (int j = 0; j < nrow; ++j)
      {
        CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(snodes[j]);

        if (cnode->IsOnBound()) continue;

        double prod = 0.0;
        // Petrov-Galerkin approach (dual LM for D/M but standard LM for gap)
        if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) prod = sval[j] * gap[0] * jac * wgt;
        // usual standard or dual LM approach
        else
          prod = lmval[j] * gap[0] * jac * wgt;

        // add current Gauss point's contribution to gseg
        cnode->AddgValue(prod);
      }
    }

    // INVALID CASES
    else
    {
      dserror("ERROR: Invalid integration case for 3D quadratic contact!");
    }
  }

  return;
}

void CONTACT::CoIntegrator::Gap_3D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
    LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& mderiv, double* gap, double* gpn,
    std::vector<GEN::pairedvector<int, double>>& dsxigp,
    std::vector<GEN::pairedvector<int, double>>& dmxigp, GEN::pairedvector<int, double>& dgapgp,
    std::vector<GEN::pairedvector<int, double>>& dnmap_unit)
{
  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();
  if (!snodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");
  if (!mnodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // number of nodes (slave, master)
  const int nrow = sele.NumNode();
  const int ncol = mele.NumNode();

  double sgpx[3] = {0.0, 0.0, 0.0};
  double mgpx[3] = {0.0, 0.0, 0.0};

  for (int i = 0; i < nrow; ++i)
  {
    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(snodes[i]);
    gpn[0] += sval[i] * mymrtrnode->MoData().n()[0];
    gpn[1] += sval[i] * mymrtrnode->MoData().n()[1];
    gpn[2] += sval[i] * mymrtrnode->MoData().n()[2];

    if (WearType() == INPAR::WEAR::wear_primvar)
    {
      FriNode* myfricnode = dynamic_cast<FriNode*>(mymrtrnode);
      sgpx[0] += sval[i] * (sele.GetNodalCoords(0, i) -
                               (myfricnode->MoData().n()[0]) * myfricnode->WearData().wcurr()[0]);
      sgpx[1] += sval[i] * (sele.GetNodalCoords(1, i) -
                               (myfricnode->MoData().n()[1]) * myfricnode->WearData().wcurr()[0]);
      sgpx[2] += sval[i] * (sele.GetNodalCoords(2, i) -
                               (myfricnode->MoData().n()[2]) * myfricnode->WearData().wcurr()[0]);
    }
    else
    {
      sgpx[0] += sval[i] * sele.GetNodalCoords(0, i);
      sgpx[1] += sval[i] * sele.GetNodalCoords(1, i);
      sgpx[2] += sval[i] * sele.GetNodalCoords(2, i);
    }
  }

  // build interpolation of master GP coordinates
  for (int i = 0; i < ncol; ++i)
  {
    if (WearSide() == INPAR::WEAR::wear_both and WearType() == INPAR::WEAR::wear_primvar)
    {
      FriNode* masternode = dynamic_cast<FriNode*>(mnodes[i]);

      mgpx[0] += mval[i] * (mele.GetNodalCoords(0, i) -
                               (masternode->MoData().n()[0] * masternode->WearData().wcurr()[0]));
      mgpx[1] += mval[i] * (mele.GetNodalCoords(1, i) -
                               (masternode->MoData().n()[1] * masternode->WearData().wcurr()[0]));
      mgpx[2] += mval[i] * (mele.GetNodalCoords(2, i) -
                               (masternode->MoData().n()[2] * masternode->WearData().wcurr()[0]));
    }
    else
    {
      mgpx[0] += mval[i] * mele.GetNodalCoords(0, i);
      mgpx[1] += mval[i] * mele.GetNodalCoords(1, i);
      mgpx[2] += mval[i] * mele.GetNodalCoords(2, i);
    }
  }

  // normalize interpolated GP normal back to length 1.0 !!!
  double lengthn = sqrt(gpn[0] * gpn[0] + gpn[1] * gpn[1] + gpn[2] * gpn[2]);
  if (lengthn < 1.0e-12) dserror("ERROR: IntegrateAndDerivSegment: Divide by zero!");

  for (int i = 0; i < 3; ++i) gpn[i] /= lengthn;

  // build gap function at current GP
  for (int i = 0; i < Dim(); ++i) gap[0] += (mgpx[i] - sgpx[i]) * gpn[i];

  // **************************
  // Linearization
  // **************************
  int linsize = 0;
  if (cppnormal_)
  {
    for (int i = 0; i < ncol; ++i)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(mnodes[i]);
      linsize += cnode->GetLinsize();
    }
    linsize += 1 + mele.NumNode();
  }
  else
  {
    for (int i = 0; i < nrow; ++i)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(snodes[i]);
      linsize += cnode->GetLinsize();
    }
  }
  linsize *= 10;

  // build directional derivative of slave GP normal (non-unit)
  GEN::pairedvector<int, double> dmap_nxsl_gp(linsize);
  GEN::pairedvector<int, double> dmap_nysl_gp(linsize);
  GEN::pairedvector<int, double> dmap_nzsl_gp(linsize);

  for (int i = 0; i < nrow; ++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*>(snodes[i]);

    GEN::pairedvector<int, double>& dmap_nxsl_i = cnode->CoData().GetDerivN()[0];
    GEN::pairedvector<int, double>& dmap_nysl_i = cnode->CoData().GetDerivN()[1];
    GEN::pairedvector<int, double>& dmap_nzsl_i = cnode->CoData().GetDerivN()[2];

    for (_CI p = dmap_nxsl_i.begin(); p != dmap_nxsl_i.end(); ++p)
      dmap_nxsl_gp[p->first] += sval[i] * (p->second);
    for (_CI p = dmap_nysl_i.begin(); p != dmap_nysl_i.end(); ++p)
      dmap_nysl_gp[p->first] += sval[i] * (p->second);
    for (_CI p = dmap_nzsl_i.begin(); p != dmap_nzsl_i.end(); ++p)
      dmap_nzsl_gp[p->first] += sval[i] * (p->second);

    for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
    {
      double valx = sderiv(i, 0) * cnode->MoData().n()[0];
      dmap_nxsl_gp[p->first] += valx * (p->second);
      double valy = sderiv(i, 0) * cnode->MoData().n()[1];
      dmap_nysl_gp[p->first] += valy * (p->second);
      double valz = sderiv(i, 0) * cnode->MoData().n()[2];
      dmap_nzsl_gp[p->first] += valz * (p->second);
    }

    for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
    {
      double valx = sderiv(i, 1) * cnode->MoData().n()[0];
      dmap_nxsl_gp[p->first] += valx * (p->second);
      double valy = sderiv(i, 1) * cnode->MoData().n()[1];
      dmap_nysl_gp[p->first] += valy * (p->second);
      double valz = sderiv(i, 1) * cnode->MoData().n()[2];
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

  // for wear as own discretization
  // lin slave nodes
  if (WearType() == INPAR::WEAR::wear_primvar)
  {
    for (int z = 0; z < nrow; ++z)
    {
      for (int k = 0; k < 3; ++k)
      {
        FriNode* frinode = dynamic_cast<FriNode*>(snodes[z]);

        dgapgp[frinode->Dofs()[k]] -= sval[z] * gpn[k];

        for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
          dgapgp[p->first] -=
              gpn[k] * sderiv(z, 0) *
              (frinode->xspatial()[k] - frinode->MoData().n()[k] * frinode->WearData().wcurr()[0]) *
              (p->second);

        for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
          dgapgp[p->first] -=
              gpn[k] * sderiv(z, 1) *
              (frinode->xspatial()[k] - frinode->MoData().n()[k] * frinode->WearData().wcurr()[0]) *
              (p->second);

        for (_CI p = frinode->CoData().GetDerivN()[k].begin();
             p != frinode->CoData().GetDerivN()[k].end(); ++p)
          dgapgp[p->first] += gpn[k] * sval[z] * frinode->WearData().wcurr()[0] * (p->second);
      }
    }
  }
  else
  {
    for (int z = 0; z < nrow; ++z)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(snodes[z]);
      for (int k = 0; k < 3; ++k) dgapgp[cnode->Dofs()[k]] -= sval[z] * gpn[k];
    }

    for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < nrow; ++z)
      {
        CoNode* cnode = dynamic_cast<CoNode*>(snodes[z]);
        for (int k = 0; k < 3; ++k) dg -= gpn[k] * sderiv(z, 0) * cnode->xspatial()[k] * ps;
      }
    }

    for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < nrow; ++z)
      {
        CoNode* cnode = dynamic_cast<CoNode*>(snodes[z]);
        for (int k = 0; k < 3; ++k) dg -= gpn[k] * sderiv(z, 1) * cnode->xspatial()[k] * ps;
      }
    }
  }

  //        MASTER
  if (WearSide() == INPAR::WEAR::wear_both and WearType() == INPAR::WEAR::wear_primvar)
  {
    for (int z = 0; z < ncol; ++z)
    {
      FriNode* frinode = dynamic_cast<FriNode*>(mnodes[z]);

      for (int k = 0; k < 3; ++k)
      {
        dgapgp[frinode->Dofs()[k]] += mval[z] * gpn[k];

        for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
          dgapgp[p->first] +=
              gpn[k] * mderiv(z, 0) *
              (frinode->xspatial()[k] - frinode->MoData().n()[k] * frinode->WearData().wcurr()[0]) *
              (p->second);

        for (_CI p = dmxigp[1].begin(); p != dmxigp[1].end(); ++p)
          dgapgp[p->first] +=
              gpn[k] * mderiv(z, 1) *
              (frinode->xspatial()[k] - frinode->MoData().n()[k] * frinode->WearData().wcurr()[0]) *
              (p->second);

        for (_CI p = frinode->CoData().GetDerivN()[k].begin();
             p != frinode->CoData().GetDerivN()[k].end(); ++p)
          dgapgp[p->first] -= gpn[k] * mval[z] * frinode->WearData().wcurr()[0] * (p->second);
      }
    }
  }
  else
  {
    // lin master nodes
    for (int z = 0; z < ncol; ++z)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(mnodes[z]);
      for (int k = 0; k < 3; ++k) dgapgp[cnode->Dofs()[k]] += mval[z] * gpn[k];
    }

    for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < ncol; ++z)
      {
        CoNode* cnode = dynamic_cast<CoNode*>(mnodes[z]);
        for (int k = 0; k < 3; ++k) dg += gpn[k] * mderiv(z, 0) * cnode->xspatial()[k] * ps;
      }
    }

    for (_CI p = dmxigp[1].begin(); p != dmxigp[1].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < ncol; ++z)
      {
        CoNode* cnode = dynamic_cast<CoNode*>(mnodes[z]);
        for (int k = 0; k < 3; ++k) dg += gpn[k] * mderiv(z, 1) * cnode->xspatial()[k] * ps;
      }
    }
  }

  return;
}

void CONTACT::CoIntegrator::Gap_2D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
    LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& mderiv, double* gap, double* gpn,
    std::vector<GEN::pairedvector<int, double>>& dsxigp,
    std::vector<GEN::pairedvector<int, double>>& dmxigp, GEN::pairedvector<int, double>& dgapgp,
    std::vector<GEN::pairedvector<int, double>>& dnmap_unit)
{
  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // get slave element nodes themselves
  DRT::Node** snodes = NULL;
  DRT::Node** mnodes = NULL;
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();

  int linsize = 0;
  for (int i = 0; i < nrow; ++i) linsize += dynamic_cast<CoNode*>(sele.Nodes()[i])->GetLinsize();
  // safety
  linsize = linsize * 2;

  int ndof = Dim();

  snodes = sele.Nodes();
  mnodes = mele.Nodes();

  // get slave nodal coords for GP evaluation
  LINALG::SerialDenseMatrix scoord(3, nrow);
  sele.GetNodalCoords(scoord);

  // get master nodal coords for GP evaluation
  LINALG::SerialDenseMatrix mcoord(3, ncol);
  mele.GetNodalCoords(mcoord);

  double sgpx[2] = {0.0, 0.0};
  double mgpx[2] = {0.0, 0.0};

  for (int i = 0; i < nrow; ++i)
  {
    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(snodes[i]);
    gpn[0] += sval[i] * mymrtrnode->MoData().n()[0];
    gpn[1] += sval[i] * mymrtrnode->MoData().n()[1];

    if (WearType() == INPAR::WEAR::wear_primvar)
    {
      FriNode* myfricnode = dynamic_cast<FriNode*>(mymrtrnode);
      double w = myfricnode->WearData().wcurr()[0] + myfricnode->WearData().waccu()[0];
      sgpx[0] += sval[i] * (scoord(0, i) - (myfricnode->MoData().n()[0]) * w);
      sgpx[1] += sval[i] * (scoord(1, i) - (myfricnode->MoData().n()[1]) * w);
    }
    else
    {
      sgpx[0] += sval[i] * scoord(0, i);
      sgpx[1] += sval[i] * scoord(1, i);
    }
  }

  // build interpolation of master GP coordinates
  for (int i = 0; i < ncol; ++i)
  {
    if (WearSide() == INPAR::WEAR::wear_both and WearType() == INPAR::WEAR::wear_primvar)
    {
      FriNode* mymrtrnodeM = dynamic_cast<FriNode*>(mnodes[i]);
      double w = mymrtrnodeM->WearData().wcurr()[0] + mymrtrnodeM->WearData().waccu()[0];
      mgpx[0] += mval[i] * (mcoord(0, i) - (mymrtrnodeM->MoData().n()[0]) * w);
      mgpx[1] += mval[i] * (mcoord(1, i) - (mymrtrnodeM->MoData().n()[1]) * w);
    }
    else
    {
      mgpx[0] += mval[i] * mcoord(0, i);
      mgpx[1] += mval[i] * mcoord(1, i);
    }
  }

  // normalize interpolated GP normal back to length 1.0 !!!
  double lengthn = sqrt(gpn[0] * gpn[0] + gpn[1] * gpn[1]);
  if (lengthn < 1.0e-12) dserror("ERROR: IntegrateAndDerivSegment: Divide by zero!");

  for (int i = 0; i < Dim(); ++i) gpn[i] /= lengthn;

  // build gap function at current GP
  for (int i = 0; i < Dim(); ++i) gap[0] += (mgpx[i] - sgpx[i]) * gpn[i];

  // **************************
  // Linearization
  // **************************

  // build directional derivative of slave GP normal (non-unit)
  GEN::pairedvector<int, double> dmap_nxsl_gp(ncol * ndof + linsize);
  GEN::pairedvector<int, double> dmap_nysl_gp(ncol * ndof + linsize);

  for (int i = 0; i < nrow; ++i)
  {
    MORTAR::MortarNode* snode = dynamic_cast<MORTAR::MortarNode*>(snodes[i]);

    GEN::pairedvector<int, double>& dmap_nxsl_i =
        dynamic_cast<CONTACT::CoNode*>(snodes[i])->CoData().GetDerivN()[0];
    GEN::pairedvector<int, double>& dmap_nysl_i =
        dynamic_cast<CONTACT::CoNode*>(snodes[i])->CoData().GetDerivN()[1];

    for (_CI p = dmap_nxsl_i.begin(); p != dmap_nxsl_i.end(); ++p)
      dmap_nxsl_gp[p->first] += sval[i] * (p->second);
    for (_CI p = dmap_nysl_i.begin(); p != dmap_nysl_i.end(); ++p)
      dmap_nysl_gp[p->first] += sval[i] * (p->second);

    for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
    {
      double valx = sderiv(i, 0) * snode->MoData().n()[0];
      dmap_nxsl_gp[p->first] += valx * (p->second);
      double valy = sderiv(i, 0) * snode->MoData().n()[1];
      dmap_nysl_gp[p->first] += valy * (p->second);
    }
  }

  // build directional derivative of slave GP normal (unit)
  const double ll = lengthn * lengthn;
  const double linv = 1.0 / lengthn;
  const double lllinv = 1.0 / (lengthn * lengthn * lengthn);
  const double sxsx = gpn[0] * gpn[0] * ll;  // gpn is the unit normal --> multiplication with ll
  const double sxsy = gpn[0] * gpn[1] * ll;  // to get the non-unit normal
  const double sysy = gpn[1] * gpn[1] * ll;

  for (_CI p = dmap_nxsl_gp.begin(); p != dmap_nxsl_gp.end(); ++p)
  {
    dnmap_unit[0][p->first] += linv * (p->second);
    dnmap_unit[0][p->first] -= lllinv * sxsx * (p->second);
    dnmap_unit[1][p->first] -= lllinv * sxsy * (p->second);
  }

  for (_CI p = dmap_nysl_gp.begin(); p != dmap_nysl_gp.end(); ++p)
  {
    dnmap_unit[1][p->first] += linv * (p->second);
    dnmap_unit[1][p->first] -= lllinv * sysy * (p->second);
    dnmap_unit[0][p->first] -= lllinv * sxsy * (p->second);
  }

  // *****************************************************************************
  // add everything to dgapgp                                                    *
  // *****************************************************************************
  for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
    dgapgp[p->first] += (mgpx[0] - sgpx[0]) * (p->second);

  for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
    dgapgp[p->first] += (mgpx[1] - sgpx[1]) * (p->second);

  // for wear as own discretization
  // slave nodes
  if (WearType() == INPAR::WEAR::wear_primvar)
  {
    for (int z = 0; z < nrow; ++z)
    {
      for (int k = 0; k < 2; ++k)
      {
        FriNode* frinode = dynamic_cast<FriNode*>(snodes[z]);
        double w = frinode->WearData().wcurr()[0] + frinode->WearData().waccu()[0];

        dgapgp[frinode->Dofs()[k]] -= sval[z] * gpn[k];

        for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
          dgapgp[p->first] -= gpn[k] * sderiv(z, 0) *
                              (frinode->xspatial()[k] - frinode->MoData().n()[k] * w) * (p->second);

        for (_CI p = frinode->CoData().GetDerivN()[k].begin();
             p != frinode->CoData().GetDerivN()[k].end(); ++p)
          dgapgp[p->first] += gpn[k] * sval[z] * w * (p->second);
      }
    }
  }
  else
  {
    for (int z = 0; z < nrow; ++z)
    {
      MORTAR::MortarNode* snode = dynamic_cast<MORTAR::MortarNode*>(snodes[z]);

      for (int k = 0; k < Dim(); ++k)
      {
        dgapgp[snode->Dofs()[k]] -= sval[z] * (gpn[k]);

        for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
          dgapgp[p->first] -= gpn[k] * sderiv(z, 0) * snode->xspatial()[k] * (p->second);
      }
    }
  }

  // **************************************************
  // master nodes
  if (WearSide() == INPAR::WEAR::wear_both and WearType() == INPAR::WEAR::wear_primvar)
  {
    for (int z = 0; z < ncol; ++z)
    {
      for (int k = 0; k < 2; ++k)
      {
        FriNode* frinode = dynamic_cast<FriNode*>(mnodes[z]);
        const double w = frinode->WearData().wcurr()[0] + frinode->WearData().waccu()[0];

        dgapgp[frinode->Dofs()[k]] += mval[z] * gpn[k];

        for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
          dgapgp[p->first] += gpn[k] * mderiv(z, 0) *
                              (frinode->xspatial()[k] - frinode->MoData().n()[k] * w) * (p->second);

        for (_CI p = frinode->CoData().GetDerivN()[k].begin();
             p != frinode->CoData().GetDerivN()[k].end(); ++p)
          dgapgp[p->first] -= gpn[k] * mval[z] * w * (p->second);
      }
    }
  }
  else
  {
    for (int z = 0; z < ncol; ++z)
    {
      MORTAR::MortarNode* mnode = dynamic_cast<MORTAR::MortarNode*>(mnodes[z]);

      for (int k = 0; k < Dim(); ++k)
      {
        dgapgp[mnode->Dofs()[k]] += mval[z] * gpn[k];

        for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
          dgapgp[p->first] += gpn[k] * mderiv(z, 0) * mnode->xspatial()[k] * (p->second);
      }
    }
  }


  return;
}


/*----------------------------------------------------------------------*
 |  Compute entries for weighted Gap at GP                   farah 12/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_3D_G_Quad_pwlin(MORTAR::MortarElement& sele,
    MORTAR::IntElement& sintele, MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval,
    LINALG::SerialDenseVector& mval, LINALG::SerialDenseVector& lmintval,
    LINALG::SerialDenseMatrix& scoord, LINALG::SerialDenseMatrix& mcoord,
    LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& mderiv, double* gap, double* gpn,
    double* lengthn, double& jac, double& wgt,
    const std::vector<GEN::pairedvector<int, double>>& dsxigp,
    const std::vector<GEN::pairedvector<int, double>>& dmxigp,
    GEN::pairedvector<int, double>& dgapgp, std::vector<GEN::pairedvector<int, double>>& dnmap_unit)
{
  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** sintnodes = sintele.Nodes();
  DRT::Node** mnodes = mele.Nodes();
  if (!snodes) dserror("ERROR: Null pointer!");
  if (!sintnodes) dserror("ERROR: Null pointer!");
  if (!mnodes) dserror("ERROR: Null pointer!");

  // number of nodes (slave, master)
  const int nrow = sele.NumNode();
  const int nintrow = sintele.NumNode();
  const int ncol = mele.NumNode();

  double sgpx[3] = {0.0, 0.0, 0.0};
  double mgpx[3] = {0.0, 0.0, 0.0};

  for (int i = 0; i < nrow; ++i)
  {
    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(snodes[i]);
    gpn[0] += sval[i] * mymrtrnode->MoData().n()[0];
    gpn[1] += sval[i] * mymrtrnode->MoData().n()[1];
    gpn[2] += sval[i] * mymrtrnode->MoData().n()[2];

    if (WearType() == INPAR::WEAR::wear_primvar)
    {
      FriNode* myfrinode = dynamic_cast<FriNode*>(snodes[i]);
      sgpx[0] += sval[i] *
                 (scoord(0, i) - (myfrinode->MoData().n()[0]) * myfrinode->WearData().wcurr()[0]);
      sgpx[1] += sval[i] *
                 (scoord(1, i) - (myfrinode->MoData().n()[1]) * myfrinode->WearData().wcurr()[0]);
      sgpx[2] += sval[i] *
                 (scoord(2, i) - (myfrinode->MoData().n()[2]) * myfrinode->WearData().wcurr()[0]);
    }
    else
    {
      sgpx[0] += sval[i] * scoord(0, i);
      sgpx[1] += sval[i] * scoord(1, i);
      sgpx[2] += sval[i] * scoord(2, i);
    }
  }

  // build interpolation of master GP coordinates
  for (int i = 0; i < ncol; ++i)
  {
    if (WearSide() == INPAR::WEAR::wear_both and WearType() == INPAR::WEAR::wear_primvar)
    {
      FriNode* masternode = dynamic_cast<FriNode*>(mnodes[i]);

      mgpx[0] += mval[i] *
                 (mcoord(0, i) - (masternode->MoData().n()[0] * masternode->WearData().wcurr()[0]));
      mgpx[1] += mval[i] *
                 (mcoord(1, i) - (masternode->MoData().n()[1] * masternode->WearData().wcurr()[0]));
      mgpx[2] += mval[i] *
                 (mcoord(2, i) - (masternode->MoData().n()[2] * masternode->WearData().wcurr()[0]));
    }
    else
    {
      mgpx[0] += mval[i] * mcoord(0, i);
      mgpx[1] += mval[i] * mcoord(1, i);
      mgpx[2] += mval[i] * mcoord(2, i);
    }
  }

  // normalize interpolated GP normal back to length 1.0 !!!
  lengthn[0] = sqrt(gpn[0] * gpn[0] + gpn[1] * gpn[1] + gpn[2] * gpn[2]);
  if (lengthn[0] < 1.0e-12) dserror("ERROR: IntegrateAndDerivSegment: Divide by zero!");

  for (int i = 0; i < 3; ++i) gpn[i] /= lengthn[0];

  // build gap function at current GP
  for (int i = 0; i < Dim(); ++i) gap[0] += (mgpx[i] - sgpx[i]) * gpn[i];

  // **************************
  // add to node
  // **************************
  // CASE 3: Standard LM shape functions and piecewise linear interpolation
  // Attention:  for this case, lmval represents lmintval !!!
  if (ShapeFcn() == INPAR::MORTAR::shape_standard && LagMultQuad() == INPAR::MORTAR::lagmult_pwlin)
  {
    for (int j = 0; j < nintrow; ++j)
    {
      CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(sintnodes[j]);

      double prod = 0.0;
      prod = lmintval[j] * gap[0] * jac * wgt;

      if (cnode->IsOnBound()) continue;

      // add current Gauss point's contribution to gseg
      cnode->AddgValue(prod);
    }
  }
  // INVALID CASES
  else
  {
    dserror("ERROR: Invalid integration case for 3D quadratic contact!");
  }

  // **************************
  // Linearization
  // **************************
  int linsize = 0;
  for (int i = 0; i < nrow; ++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*>(snodes[i]);
    linsize += cnode->GetLinsize();
  }

  // build directional derivative of slave GP normal (non-unit)
  GEN::pairedvector<int, double> dmap_nxsl_gp(linsize);
  GEN::pairedvector<int, double> dmap_nysl_gp(linsize);
  GEN::pairedvector<int, double> dmap_nzsl_gp(linsize);

  for (int i = 0; i < nrow; ++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*>(snodes[i]);

    GEN::pairedvector<int, double>& dmap_nxsl_i = cnode->CoData().GetDerivN()[0];
    GEN::pairedvector<int, double>& dmap_nysl_i = cnode->CoData().GetDerivN()[1];
    GEN::pairedvector<int, double>& dmap_nzsl_i = cnode->CoData().GetDerivN()[2];

    for (_CI p = dmap_nxsl_i.begin(); p != dmap_nxsl_i.end(); ++p)
      dmap_nxsl_gp[p->first] += sval[i] * (p->second);
    for (_CI p = dmap_nysl_i.begin(); p != dmap_nysl_i.end(); ++p)
      dmap_nysl_gp[p->first] += sval[i] * (p->second);
    for (_CI p = dmap_nzsl_i.begin(); p != dmap_nzsl_i.end(); ++p)
      dmap_nzsl_gp[p->first] += sval[i] * (p->second);

    for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
    {
      double valx = sderiv(i, 0) * cnode->MoData().n()[0];
      dmap_nxsl_gp[p->first] += valx * (p->second);
      double valy = sderiv(i, 0) * cnode->MoData().n()[1];
      dmap_nysl_gp[p->first] += valy * (p->second);
      double valz = sderiv(i, 0) * cnode->MoData().n()[2];
      dmap_nzsl_gp[p->first] += valz * (p->second);
    }

    for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
    {
      double valx = sderiv(i, 1) * cnode->MoData().n()[0];
      dmap_nxsl_gp[p->first] += valx * (p->second);
      double valy = sderiv(i, 1) * cnode->MoData().n()[1];
      dmap_nysl_gp[p->first] += valy * (p->second);
      double valz = sderiv(i, 1) * cnode->MoData().n()[2];
      dmap_nzsl_gp[p->first] += valz * (p->second);
    }
  }

  const double ll = lengthn[0] * lengthn[0];
  const double linv = 1.0 / (lengthn[0]);
  const double lllinv = 1.0 / (lengthn[0] * lengthn[0] * lengthn[0]);
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

  // for wear as own discretization
  // lin slave nodes
  if (WearType() == INPAR::WEAR::wear_primvar)
  {
    for (int z = 0; z < nrow; ++z)
    {
      for (int k = 0; k < 3; ++k)
      {
        FriNode* frinode = dynamic_cast<FriNode*>(snodes[z]);

        dgapgp[frinode->Dofs()[k]] -= sval[z] * gpn[k];

        for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
          dgapgp[p->first] -=
              gpn[k] * sderiv(z, 0) *
              (frinode->xspatial()[k] - frinode->MoData().n()[k] * frinode->WearData().wcurr()[0]) *
              (p->second);

        for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
          dgapgp[p->first] -=
              gpn[k] * sderiv(z, 1) *
              (frinode->xspatial()[k] - frinode->MoData().n()[k] * frinode->WearData().wcurr()[0]) *
              (p->second);

        for (_CI p = frinode->CoData().GetDerivN()[k].begin();
             p != frinode->CoData().GetDerivN()[k].end(); ++p)
          dgapgp[p->first] += gpn[k] * sval[z] * frinode->WearData().wcurr()[0] * (p->second);
      }
    }
  }
  else
  {
    for (int z = 0; z < nrow; ++z)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(snodes[z]);

      for (int k = 0; k < 3; ++k)
      {
        dgapgp[cnode->Dofs()[k]] -= sval[z] * gpn[k];

        for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
          dgapgp[p->first] -= gpn[k] * sderiv(z, 0) * cnode->xspatial()[k] * (p->second);

        for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
          dgapgp[p->first] -= gpn[k] * sderiv(z, 1) * cnode->xspatial()[k] * (p->second);
      }
    }
  }

  //        MASTER
  if (WearSide() == INPAR::WEAR::wear_both and WearType() == INPAR::WEAR::wear_primvar)
  {
    for (int z = 0; z < ncol; ++z)
    {
      FriNode* frinode = dynamic_cast<FriNode*>(mnodes[z]);

      for (int k = 0; k < 3; ++k)
      {
        dgapgp[frinode->Dofs()[k]] += mval[z] * gpn[k];

        for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
          dgapgp[p->first] +=
              gpn[k] * mderiv(z, 0) *
              (frinode->xspatial()[k] - frinode->MoData().n()[k] * frinode->WearData().wcurr()[0]) *
              (p->second);

        for (_CI p = dmxigp[1].begin(); p != dmxigp[1].end(); ++p)
          dgapgp[p->first] +=
              gpn[k] * mderiv(z, 1) *
              (frinode->xspatial()[k] - frinode->MoData().n()[k] * frinode->WearData().wcurr()[0]) *
              (p->second);

        for (_CI p = frinode->CoData().GetDerivN()[k].begin();
             p != frinode->CoData().GetDerivN()[k].end(); ++p)
          dgapgp[p->first] -= gpn[k] * mval[z] * frinode->WearData().wcurr()[0] * (p->second);
      }
    }
  }
  else
  {
    // lin master nodes
    for (int z = 0; z < ncol; ++z)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(mnodes[z]);

      for (int k = 0; k < 3; ++k)
      {
        dgapgp[cnode->Dofs()[k]] += mval[z] * gpn[k];

        for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
          dgapgp[p->first] += gpn[k] * mderiv(z, 0) * cnode->xspatial()[k] * (p->second);

        for (_CI p = dmxigp[1].begin(); p != dmxigp[1].end(); ++p)
          dgapgp[p->first] += gpn[k] * mderiv(z, 1) * cnode->xspatial()[k] * (p->second);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Do lin. entries for weighted Gap at GP                   farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_2D_G_Lin(int& iter, MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
    LINALG::SerialDenseVector& lmval, LINALG::SerialDenseMatrix& sderiv,
    LINALG::SerialDenseMatrix& lmderiv, double& gap, double* gpn, double& jac, double& wgt,
    GEN::pairedvector<int, double>& dgapgp, GEN::pairedvector<int, double>& jacintcellmap,
    std::vector<GEN::pairedvector<int, double>>& dsxigp,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap)
{
  // get slave element nodes themselves
  DRT::Node** snodes = NULL;
  int nrow = sele.NumNode();
  snodes = sele.Nodes();

  MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(snodes[iter]);
  if (!mymrtrnode) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");
  if (mymrtrnode->IsOnBoundorCE()) return;

  double fac = 0.0;

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // get the corresponding map as a reference
  std::map<int, double>& dgmap = dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();

  // switch if Petrov-Galerkin approach for LM is applied
  if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
  {
    // (1) Lin(Phi) - does not exist in gap for Petrov-Galerkin interpolation
    // as std shape functions are used here

    // (2) Lin(Phi) - slave GP coordinates
    fac = wgt * sderiv(iter, 0) * gap * jac;
    for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p) dgmap[p->first] += fac * (p->second);

    // (3) Lin(g) - gap function
    fac = wgt * sval[iter] * jac;
    for (_CI p = dgapgp.begin(); p != dgapgp.end(); ++p) dgmap[p->first] += fac * (p->second);

    // (4) Lin(dxdsxi) - slave GP Jacobian
    fac = wgt * sval[iter] * gap;
    for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
      dgmap[p->first] += fac * (p->second);
  }

  // the usual standard or dual LM interpolation
  else
  {
    // (1) Lin(Phi) - dual shape functions
    if (ShapeFcn() == INPAR::MORTAR::shape_dual)
    {
      for (int m = 0; m < nrow; ++m)
      {
        fac = wgt * sval[m] * gap * jac;
        for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
             p != dualmap.end(); ++p)
          dgmap[p->first] += fac * (p->second)(iter, m);
      }
    }

    // (2) Lin(Phi) - slave GP coordinates
    fac = wgt * lmderiv(iter, 0) * gap * jac;
    for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p) dgmap[p->first] += fac * (p->second);

    // (3) Lin(g) - gap function
    fac = wgt * lmval[iter] * jac;
    for (_CI p = dgapgp.begin(); p != dgapgp.end(); ++p) dgmap[p->first] += fac * (p->second);

    // (4) Lin(dxdsxi) - slave GP Jacobian
    fac = wgt * lmval[iter] * gap;
    for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
      dgmap[p->first] += fac * (p->second);
  }

  //****************************************************************
  // LIN WEAR W.R.T. LM
  //****************************************************************
  if (WearType() == INPAR::WEAR::wear_primvar)
  {
    // get master element nodes themselves
    const int ncol = mele.NumNode();
    DRT::Node** mnodes = mele.Nodes();

    std::map<int, double>& dgwmmap =
        dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivGW();

    for (int bl = 0; bl < nrow; ++bl)
    {
      MORTAR::MortarNode* wearnode = dynamic_cast<MORTAR::MortarNode*>(snodes[bl]);
      for (int z = 0; z < Dim(); ++z)
        dgwmmap[wearnode->Dofs()[0]] +=
            jac * wgt * lmval[iter] * (gpn[z] * sval[bl] * wearnode->MoData().n()[z]);
    }

    if (WearSide() == INPAR::WEAR::wear_both)
    {
      for (int bl = 0; bl < ncol; ++bl)
      {
        MORTAR::MortarNode* wearnodeM = dynamic_cast<MORTAR::MortarNode*>(mnodes[bl]);
        for (int z = 0; z < Dim(); ++z)
          dgwmmap[wearnodeM->Dofs()[0]] -=
              jac * wgt * lmval[iter] * (gpn[z] * mval[bl] * wearnodeM->MoData().n()[z]);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Do lin. entries for weighted Gap at GP                   farah 12/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_3D_G_Quad_pwlin_Lin(int& iter, MORTAR::IntElement& sintele,
    LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& lmintval,
    LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& lmintderiv, double& gap,
    double* gpn, double& jac, double& wgt, const GEN::pairedvector<int, double>& dgapgp,
    const GEN::pairedvector<int, double>& jacintcellmap,
    const std::vector<GEN::pairedvector<int, double>>& dsxigp)
{
  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator CI;

  // get slave element nodes themselves
  DRT::Node** sintnodes = sintele.Nodes();

  MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(sintnodes[iter]);
  if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

  double fac = 0.0;

  // CASE 3: Standard LM shape functions and piecewise linear interpolation
  if (ShapeFcn() == INPAR::MORTAR::shape_standard && LagMultQuad() == INPAR::MORTAR::lagmult_pwlin)
  {
    // get the corresponding map as a reference
    std::map<int, double>& dgmap = dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();

    // (1) Lin(Phi) - dual shape functions
    // this vanishes here since there are no deformation-dependent dual functions

    // (2) Lin(Phi) - slave GP coordinates
    fac = wgt * lmintderiv(iter, 0) * gap * jac;
    for (CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p) dgmap[p->first] += fac * (p->second);

    fac = wgt * lmintderiv(iter, 1) * gap * jac;
    for (CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p) dgmap[p->first] += fac * (p->second);

    // (3) Lin(g) - gap function
    fac = wgt * lmintval[iter] * jac;
    for (CI p = dgapgp.begin(); p != dgapgp.end(); ++p) dgmap[p->first] += fac * (p->second);

    // (4) Lin(dsxideta) - intcell GP Jacobian
    fac = wgt * lmintval[iter] * gap;
    for (CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
      dgmap[p->first] += fac * (p->second);
  }
  else
    dserror("shapefcn-lagmult combination not supported!");


  return;
}

/*----------------------------------------------------------------------*
 |  Do lin. entries for weighted Gap at GP                   farah 12/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_3D_G_Quad_Lin(int& iter, MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval,
    LINALG::SerialDenseVector& svalmod, LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& lmderiv, double& gap, double* gpn,
    double& jac, double& wgt, bool& duallin, const GEN::pairedvector<int, double>& dgapgp,
    const GEN::pairedvector<int, double>& jacintcellmap,
    const std::vector<GEN::pairedvector<int, double>>& dpsxigp,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap, bool dualquad3d)
{
  const int nrow = sele.NumNode();

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();

  MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(snodes[iter]);
  if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

  double fac = 0.0;

  // compute cell gap linearization ************************************
  // CASE 1/2: Standard LM shape functions and quadratic or linear interpolation
  if (ShapeFcn() == INPAR::MORTAR::shape_standard &&
      (LagMultQuad() == INPAR::MORTAR::lagmult_quad || LagMultQuad() == INPAR::MORTAR::lagmult_lin))
  {
    // get the corresponding map as a reference
    std::map<int, double>& dgmap = dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();

    // (1) Lin(Phi) - dual shape functions
    // this vanishes here since there are no deformation-dependent dual functions

    // (2) Lin(Phi) - slave GP coordinates
    fac = wgt * lmderiv(iter, 0) * gap * jac;
    for (_CI p = dpsxigp[0].begin(); p != dpsxigp[0].end(); ++p)
      dgmap[p->first] += fac * (p->second);

    fac = wgt * lmderiv(iter, 1) * gap * jac;
    for (_CI p = dpsxigp[1].begin(); p != dpsxigp[1].end(); ++p)
      dgmap[p->first] += fac * (p->second);

    // (3) Lin(g) - gap function
    fac = wgt * lmval[iter] * jac;
    for (_CI p = dgapgp.begin(); p != dgapgp.end(); ++p) dgmap[p->first] += fac * (p->second);

    // (4) Lin(dsxideta) - intcell GP Jacobian
    fac = wgt * lmval[iter] * gap;
    for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
      dgmap[p->first] += fac * (p->second);
  }

  // CASE 4: Dual LM shape functions and quadratic or linear interpolation
  else if ((ShapeFcn() == INPAR::MORTAR::shape_dual ||
               ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) &&
           (LagMultQuad() == INPAR::MORTAR::lagmult_quad ||
               LagMultQuad() == INPAR::MORTAR::lagmult_lin))
  {
    // get the corresponding map as a reference
    std::map<int, double>& dgmap = dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();

    if (ShapeFcn() == INPAR::MORTAR::shape_dual)
    {
      // (1) Lin(Phi) - dual shape functions
      if (duallin)
      {
        for (int m = 0; m < nrow; ++m)
        {
          if (dualquad3d)
            fac = wgt * svalmod[m] * gap * jac;
          else
            fac = wgt * sval[m] * gap * jac;

          for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
               p != dualmap.end(); ++p)
            dgmap[p->first] += fac * (p->second)(iter, m);
        }
      }

      // (2) Lin(Phi) - slave GP coordinates
      fac = wgt * lmderiv(iter, 0) * gap * jac;
      for (_CI p = dpsxigp[0].begin(); p != dpsxigp[0].end(); ++p)
        dgmap[p->first] += fac * (p->second);

      fac = wgt * lmderiv(iter, 1) * gap * jac;
      for (_CI p = dpsxigp[1].begin(); p != dpsxigp[1].end(); ++p)
        dgmap[p->first] += fac * (p->second);

      // (3) Lin(g) - gap function
      fac = wgt * lmval[iter] * jac;
      for (_CI p = dgapgp.begin(); p != dgapgp.end(); ++p) dgmap[p->first] += fac * (p->second);

      // (4) Lin(dsxideta) - intcell GP Jacobian
      fac = wgt * lmval[iter] * gap;
      for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
        dgmap[p->first] += fac * (p->second);
    }
    else if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
    {
      // (1) Lin(Phi) - dual shape functions

      // (2) Lin(Phi) - slave GP coordinates
      fac = wgt * sderiv(iter, 0) * gap * jac;
      for (_CI p = dpsxigp[0].begin(); p != dpsxigp[0].end(); ++p)
        dgmap[p->first] += fac * (p->second);

      fac = wgt * sderiv(iter, 1) * gap * jac;
      for (_CI p = dpsxigp[1].begin(); p != dpsxigp[1].end(); ++p)
        dgmap[p->first] += fac * (p->second);

      // (3) Lin(g) - gap function
      fac = wgt * sval[iter] * jac;
      for (_CI p = dgapgp.begin(); p != dgapgp.end(); ++p) dgmap[p->first] += fac * (p->second);

      // (4) Lin(dsxideta) - intcell GP Jacobian
      fac = wgt * sval[iter] * gap;
      for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
        dgmap[p->first] += fac * (p->second);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Do lin. entries for weighted Gap at GP                   farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::GP_G_Lin(int& iter, MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
    LINALG::SerialDenseVector& lmval, LINALG::SerialDenseMatrix& sderiv,
    LINALG::SerialDenseMatrix& lmderiv, double& gap, double* gpn, double& jac, double& wgt,
    GEN::pairedvector<int, double>& dgapgp, GEN::pairedvector<int, double>& jacintcellmap,
    std::vector<GEN::pairedvector<int, double>>& dsxigp,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap)
{
  const int nrow = sele.NumNode();
  const int ncol = mele.NumNode();

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();

  MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(snodes[iter]);
  if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

  if (mymrtrnode->IsOnBoundorCE()) return;

  static double fac = 0.0;

  // get the corresponding map as a reference
  std::map<int, double>& dgmap = dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();

  // switch if Petrov-Galerkin approach for LM is applied
  if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
  {
    // (1) Lin(Phi) - does not exist here for Petrov-Galerkin approach

    // (2) Lin(N) - slave GP coordinates
    for (int d = 0; d < Dim() - 1; ++d)
    {
      fac = wgt * sderiv(iter, d) * gap * jac;
      for (_CI p = dsxigp[d].begin(); p != dsxigp[d].end(); ++p)
        dgmap[p->first] += fac * (p->second);
    }

    // (3) Lin(g) - gap function
    fac = wgt * sval[iter] * jac;
    for (_CI p = dgapgp.begin(); p != dgapgp.end(); ++p) dgmap[p->first] += fac * (p->second);

    // (4) Lin(dsxideta) - intcell GP Jacobian
    fac = wgt * sval[iter] * gap;
    for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
      dgmap[p->first] += fac * (p->second);
  }

  // the usual standard or dual LM approach
  else
  {
    // (1) Lin(Phi) - dual shape functions
    if (dualmap.size() > 0)
      for (int m = 0; m < nrow; ++m)
      {
        fac = wgt * sval[m] * gap * jac;
        for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
             p != dualmap.end(); ++p)
          dgmap[p->first] += fac * (p->second)(iter, m);
      }

    // (2) Lin(Phi) - slave GP coordinates
    for (int d = 0; d < Dim() - 1; ++d)
    {
      fac = wgt * lmderiv(iter, d) * gap * jac;
      for (_CI p = dsxigp[d].begin(); p != dsxigp[d].end(); ++p)
        dgmap[p->first] += fac * (p->second);
    }

    // (3) Lin(g) - gap function
    fac = wgt * lmval[iter] * jac;
    for (_CI p = dgapgp.begin(); p != dgapgp.end(); ++p) dgmap[p->first] += fac * (p->second);

    // (4) Lin(dsxideta) - intcell GP Jacobian
    fac = wgt * lmval[iter] * gap;
    for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
      dgmap[p->first] += fac * (p->second);
  }

  //****************************************************************
  // LIN WEAR W.R.T. W
  //****************************************************************
  if (WearType() == INPAR::WEAR::wear_primvar)
  {
    std::map<int, double>& dgwmmap =
        dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivGW();

    for (int bl = 0; bl < nrow; ++bl)
    {
      MORTAR::MortarNode* wearnode = dynamic_cast<MORTAR::MortarNode*>(snodes[bl]);
      for (int z = 0; z < Dim(); ++z)
        dgwmmap[wearnode->Dofs()[0]] +=
            jac * wgt * lmval[iter] * (gpn[z] * sval[bl] * wearnode->MoData().n()[z]);
    }

    if (WearSide() == INPAR::WEAR::wear_both)
    {
      for (int bl = 0; bl < ncol; ++bl)
      {
        MORTAR::MortarNode* wearnodeM = dynamic_cast<MORTAR::MortarNode*>(mnodes[bl]);
        for (int z = 0; z < Dim(); ++z)
          dgwmmap[wearnodeM->Dofs()[0]] -=
              jac * wgt * lmval[iter] * (gpn[z] * mval[bl] * wearnodeM->MoData().n()[z]);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Compute entries for bound DM                             farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::GP_3D_DM_Lin_bound(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
    LINALG::SerialDenseVector& lmval, LINALG::SerialDenseMatrix& sderiv,
    LINALG::SerialDenseMatrix& lmderiv, LINALG::SerialDenseMatrix& mderiv, double& jac, double& wgt,
    const GEN::pairedvector<int, double>& derivjac,
    std::vector<GEN::pairedvector<int, double>>& dsxigp,
    std::vector<GEN::pairedvector<int, double>>& dmxigp,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap)
{
  const int nrow = sele.NumNode();
  const int ncol = mele.NumNode();

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  if (ShapeFcn() == INPAR::MORTAR::shape_standard)
  {
    // integrate LinD
    for (int j = 0; j < nrow; ++j)
    {
      MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(snodes[j]);
      if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

      // node j is boundary node
      if (mymrtrnode->IsOnBoundorCE()) continue;

      // integrate LinM
      for (int k = 0; k < ncol; ++k)
      {
        // global master node ID
        int mgid = mele.Nodes()[k]->Id();
        static double fac = 0.0;

        // get the correct map as a reference
        std::map<int, double>& dmmap_jk =
            dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

        // (1) Lin(Phi) - dual shape functions
        // this vanishes here since there are no deformation-dependent dual functions

        // (2) Lin(NSlave) - slave GP coordinates
        fac = wgt * lmderiv(j, 0) * mval[k] * jac;
        for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);

        fac = wgt * lmderiv(j, 1) * mval[k] * jac;
        for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);

        // (3) Lin(NMaster) - master GP coordinates
        fac = wgt * lmval[j] * mderiv(k, 0) * jac;
        for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);

        fac = wgt * lmval[j] * mderiv(k, 1) * jac;
        for (_CI p = dmxigp[1].begin(); p != dmxigp[1].end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);

        // (4) Lin(dsxideta) - intcell GP Jacobian
        fac = wgt * lmval[j] * mval[k];
        for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);
      }  // loop over master nodes

      for (int k = 0; k < nrow; ++k)
      {
        MORTAR::MortarNode* mymrtrnode2 = dynamic_cast<MORTAR::MortarNode*>(snodes[k]);
        if (!mymrtrnode2) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

        // global master node ID
        int sgid = mymrtrnode2->Id();
        static double fac = 0.0;

        // node k is boundary node
        if (mymrtrnode2->IsOnBoundorCE())
        {
          // get the correct map pointer
          std::map<int, double>* dmmap_jk = NULL;
          double sign = 0.;

          // for boundary nodes, we assemble to M-matrix, since the corresponding nodes do not
          // carry a LM-dof
          if (mymrtrnode2->IsOnBound())
          {
            dmmap_jk = &(dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[sgid]);
            sign = -1.;
          }
          // for corners/edges, we assemble to D-matrix, since the corresponding nodes
          // carry a LM-dof, but are integrated some place else. Here we only do surface
          // integration.
          else if (mymrtrnode2->IsOnCornerEdge())
          {
            dmmap_jk = &(dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid]);
            sign = +1.;
          }
          else
            dserror("you should never be here");

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = sign * wgt * lmderiv(j, 0) * sval[k] * jac;
          for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
            (*dmmap_jk)[p->first] += fac * (p->second);

          fac = sign * wgt * lmderiv(j, 1) * sval[k] * jac;
          for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
            (*dmmap_jk)[p->first] += fac * (p->second);

          // (3) Lin(NSlave) - slave GP coordinates
          fac = sign * wgt * lmval[j] * sderiv(k, 0) * jac;
          for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
            (*dmmap_jk)[p->first] += fac * (p->second);

          fac = sign * wgt * lmval[j] * sderiv(k, 1) * jac;
          for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
            (*dmmap_jk)[p->first] += fac * (p->second);

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = sign * wgt * lmval[j] * sval[k];
          for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
            (*dmmap_jk)[p->first] += fac * (p->second);
        }

        // node k is NO boundary node
        else
        {
          // get the correct map as a reference
          std::map<int, double>& ddmap_jk =
              dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmderiv(j, 0) * sval[k] * jac;
          for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);

          fac = wgt * lmderiv(j, 1) * sval[k] * jac;
          for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);

          // (3) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmval[j] * sderiv(k, 0) * jac;
          for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);

          fac = wgt * lmval[j] * sderiv(k, 1) * jac;
          for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt * lmval[j] * sval[k];
          for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);
        }
      }  // loop over slave nodes
    }
  }
  //--------------------------------------------------------
  else if (ShapeFcn() == INPAR::MORTAR::shape_dual ||
           ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
  {
    // integrate LinD
    for (int j = 0; j < nrow; ++j)
    {
      MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(snodes[j]);
      if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

      // node j is boundary node
      if (mymrtrnode->IsOnBoundorCE()) continue;

      // integrate LinM
      for (int k = 0; k < ncol; ++k)
      {
        // global master node ID
        int mgid = mele.Nodes()[k]->Id();
        static double fac = 0.0;

        // get the correct map as a reference
        std::map<int, double>& dmmap_jk =
            dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

        // (1) Lin(Phi) - dual shape functions
        if (dualmap.size() > 0)
        {
          for (int m = 0; m < nrow; ++m)
          {
            fac = wgt * sval[m] * mval[k] * jac;
            for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p =
                     dualmap.begin();
                 p != dualmap.end(); ++p)
              dmmap_jk[p->first] += fac * (p->second)(j, m);
          }
        }

        // (2) Lin(NSlave) - slave GP coordinates
        fac = wgt * lmderiv(j, 0) * mval[k] * jac;
        for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);

        fac = wgt * lmderiv(j, 1) * mval[k] * jac;
        for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);

        // (3) Lin(NMaster) - master GP coordinates
        fac = wgt * lmval[j] * mderiv(k, 0) * jac;
        for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);

        fac = wgt * lmval[j] * mderiv(k, 1) * jac;
        for (_CI p = dmxigp[1].begin(); p != dmxigp[1].end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);

        // (4) Lin(dsxideta) - intcell GP Jacobian
        fac = wgt * lmval[j] * mval[k];
        for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);
      }  // loop over master nodes

      // loop over slave nodes
      for (int k = 0; k < nrow; ++k)
      {
        MORTAR::MortarNode* mymrtrnode2 = dynamic_cast<MORTAR::MortarNode*>(snodes[k]);
        if (!mymrtrnode2) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

        // global master node ID
        int sgid = mymrtrnode2->Id();
        static double fac = 0.0;

        // node k is boundary node
        if (mymrtrnode2->IsOnBoundorCE())
        {
          // get the correct map pointer
          std::map<int, double>* dmmap_jk = NULL;
          double sign = 0.;

          // for boundary nodes, we assemble to M-matrix, since the corresponding nodes do not
          // carry a LM-dof
          if (mymrtrnode2->IsOnBound())
          {
            dmmap_jk = &(dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[sgid]);
            sign = -1.;
          }
          // for corners/edges, we assemble to D-matrix, since the corresponding nodes
          // carry a LM-dof, but are integrated some place else. Here we only do surface
          // integration.
          else if (mymrtrnode2->IsOnCornerEdge())
          {
            dmmap_jk = &(dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid]);
            sign = +1.;
          }
          else
            dserror("you should never be here");

          // (1) Lin(Phi) - dual shape functions
          if (dualmap.size() > 0)
          {
            for (int m = 0; m < nrow; ++m)
            {
              fac = sign * wgt * sval[m] * sval[k] * jac;
              for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p =
                       dualmap.begin();
                   p != dualmap.end(); ++p)
                (*dmmap_jk)[p->first] += fac * (p->second)(j, m);
            }
          }

          // (2) Lin(NSlave) - slave GP coordinates
          fac = sign * wgt * lmderiv(j, 0) * sval[k] * jac;
          for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
            (*dmmap_jk)[p->first] += fac * (p->second);

          fac = sign * wgt * lmderiv(j, 1) * sval[k] * jac;
          for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
            (*dmmap_jk)[p->first] += fac * (p->second);

          // (3) Lin(NSlave) - slave GP coordinates
          fac = sign * wgt * lmval[j] * sderiv(k, 0) * jac;
          for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
            (*dmmap_jk)[p->first] += fac * (p->second);

          fac = sign * wgt * lmval[j] * sderiv(k, 1) * jac;
          for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
            (*dmmap_jk)[p->first] += fac * (p->second);

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = sign * wgt * lmval[j] * sval[k];
          for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
            (*dmmap_jk)[p->first] += fac * (p->second);
        }

        // node k is NO boundary node
        else
        {
          // get the correct map as a reference
          std::map<int, double>& ddmap_jk =
              dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[mymrtrnode->Id()];

          // (1) Lin(Phi) - dual shape functions
          if (dualmap.size() > 0)
          {
            for (int m = 0; m < nrow; ++m)
            {
              fac = wgt * sval[m] * sval[k] * jac;
              for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p =
                       dualmap.begin();
                   p != dualmap.end(); ++p)
                ddmap_jk[p->first] += fac * (p->second)(j, m);
            }
          }

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmderiv(j, 0) * sval[k] * jac;
          for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);

          fac = wgt * lmderiv(j, 1) * sval[k] * jac;
          for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);

          // (3) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmval[j] * sderiv(k, 0) * jac;
          for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);

          fac = wgt * lmval[j] * sderiv(k, 1) * jac;
          for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt * lmval[j] * sval[k];
          for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);
        }
      }  // loop over slave nodes
    }
  }
  else
    dserror("ERROR: unknwon shape function type!");

  return;
}

/*----------------------------------------------------------------------*
 |  Compute entries for weighted Gap at GP                   farah 07/16|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_2D_DM_Lin_bound(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
    LINALG::SerialDenseVector& lmval, LINALG::SerialDenseMatrix& sderiv,
    LINALG::SerialDenseMatrix& mderiv, LINALG::SerialDenseMatrix& lmderiv, double& jac, double& wgt,
    const GEN::pairedvector<int, double>& derivjac,
    std::vector<GEN::pairedvector<int, double>>& derivsxi,
    std::vector<GEN::pairedvector<int, double>>& derivmxi,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap)
{
  const int nrow = sele.NumNode();
  const int ncol = mele.NumNode();

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  if (ShapeFcn() == INPAR::MORTAR::shape_standard)
  {
    // integrate LinD
    for (int j = 0; j < nrow; ++j)
    {
      MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(snodes[j]);
      if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

      // node j is boundary node
      if (mymrtrnode->IsOnBoundorCE()) continue;

      // integrate LinM
      for (int k = 0; k < ncol; ++k)
      {
        // global master node ID
        int mgid = mele.Nodes()[k]->Id();
        static double fac = 0.0;

        // get the correct map as a reference
        std::map<int, double>& dmmap_jk =
            dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

        // (1) Lin(Phi) - dual shape functions
        // this vanishes here since there are no deformation-dependent dual functions

        // (2) Lin(NSlave) - slave GP coordinates
        fac = wgt * lmderiv(j, 0) * mval[k] * jac;
        for (_CI p = derivsxi[0].begin(); p != derivsxi[0].end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);

        //        fac = wgt*lmderiv(j, 1)*mval[k]*jac;
        //        for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
        //          dmmap_jk[p->first] += fac*(p->second);

        // (3) Lin(NMaster) - master GP coordinates
        fac = wgt * lmval[j] * mderiv(k, 0) * jac;
        for (_CI p = derivmxi[0].begin(); p != derivmxi[0].end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);

        //        fac = wgt*lmval[j]*mderiv(k, 1)*jac;
        //        for (_CI p=dmxigp[1].begin(); p!=dmxigp[1].end(); ++p)
        //          dmmap_jk[p->first] += fac*(p->second);

        // (4) Lin(dsxideta) - intcell GP Jacobian
        fac = wgt * lmval[j] * mval[k];
        for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);
      }  // loop over master nodes

      for (int k = 0; k < nrow; ++k)
      {
        MORTAR::MortarNode* mymrtrnode2 = dynamic_cast<MORTAR::MortarNode*>(snodes[k]);
        if (!mymrtrnode2) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

        // global master node ID
        int sgid = mymrtrnode2->Id();
        static double fac = 0.0;

        // node k is boundary node
        if (mymrtrnode2->IsOnBoundorCE())
        {
          // move entry to derivM (with minus sign)
          // get the correct map as a reference
          std::map<int, double>& dmmap_jk =
              dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[sgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmderiv(j, 0) * sval[k] * jac;
          for (_CI p = derivsxi[0].begin(); p != derivsxi[0].end(); ++p)
            dmmap_jk[p->first] -= fac * (p->second);

          //          fac = wgt*lmderiv(j, 1)*sval[k]*jac;
          //          for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          //            dmmap_jk[p->first] -= fac*(p->second);

          // (3) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmval[j] * sderiv(k, 0) * jac;
          for (_CI p = derivsxi[0].begin(); p != derivsxi[0].end(); ++p)
            dmmap_jk[p->first] -= fac * (p->second);

          //          fac = wgt*lmval[j]*sderiv(k, 1)*jac;
          //          for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          //            dmmap_jk[p->first] -= fac*(p->second);

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt * lmval[j] * sval[k];
          for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
            dmmap_jk[p->first] -= fac * (p->second);
        }

        // node k is NO boundary node
        else
        {
          // get the correct map as a reference
          std::map<int, double>& ddmap_jk =
              dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmderiv(j, 0) * sval[k] * jac;
          for (_CI p = derivsxi[0].begin(); p != derivsxi[0].end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);

          //          fac = wgt*lmderiv(j, 1)*sval[k]*jac;
          //          for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          //            ddmap_jk[p->first] += fac*(p->second);

          // (3) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmval[j] * sderiv(k, 0) * jac;
          for (_CI p = derivsxi[0].begin(); p != derivsxi[0].end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);

          //          fac = wgt*lmval[j]*sderiv(k, 1)*jac;
          //          for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          //            ddmap_jk[p->first] += fac*(p->second);

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt * lmval[j] * sval[k];
          for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);
        }
      }  // loop over slave nodes
    }
  }
  //--------------------------------------------------------
  else if (ShapeFcn() == INPAR::MORTAR::shape_dual ||
           ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
  {
    // integrate LinD
    for (int j = 0; j < nrow; ++j)
    {
      MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(snodes[j]);
      if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

      // node j is boundary node
      if (mymrtrnode->IsOnBoundorCE()) continue;

      // integrate LinM
      for (int k = 0; k < ncol; ++k)
      {
        // global master node ID
        int mgid = mele.Nodes()[k]->Id();
        static double fac = 0.0;

        // get the correct map as a reference
        std::map<int, double>& dmmap_jk =
            dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

        // (1) Lin(Phi) - dual shape functions
        if (dualmap.size() > 0)
        {
          for (int m = 0; m < nrow; ++m)
          {
            fac = wgt * sval[m] * mval[k] * jac;
            for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p =
                     dualmap.begin();
                 p != dualmap.end(); ++p)
              dmmap_jk[p->first] += fac * (p->second)(j, m);
          }
        }

        // (2) Lin(NSlave) - slave GP coordinates
        fac = wgt * lmderiv(j, 0) * mval[k] * jac;
        for (_CI p = derivsxi[0].begin(); p != derivsxi[0].end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);

        //        fac = wgt*lmderiv(j, 1)*mval[k]*jac;
        //        for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
        //          dmmap_jk[p->first] += fac*(p->second);

        // (3) Lin(NMaster) - master GP coordinates
        fac = wgt * lmval[j] * mderiv(k, 0) * jac;
        for (_CI p = derivmxi[0].begin(); p != derivmxi[0].end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);

        //        fac = wgt*lmval[j]*mderiv(k, 1)*jac;
        //        for (_CI p=dmxigp[1].begin(); p!=dmxigp[1].end(); ++p)
        //          dmmap_jk[p->first] += fac*(p->second);

        // (4) Lin(dsxideta) - intcell GP Jacobian
        fac = wgt * lmval[j] * mval[k];
        for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);
      }  // loop over master nodes

      // loop over slave nodes
      for (int k = 0; k < nrow; ++k)
      {
        MORTAR::MortarNode* mymrtrnode2 = dynamic_cast<MORTAR::MortarNode*>(snodes[k]);
        if (!mymrtrnode2) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

        // global master node ID
        int sgid = mymrtrnode2->Id();
        static double fac = 0.0;

        // node k is boundary node
        if (mymrtrnode2->IsOnBoundorCE())
        {
          // move entry to derivM (with minus sign)
          // get the correct map as a reference
          std::map<int, double>& dmmap_jk =
              dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[sgid];

          // (1) Lin(Phi) - dual shape functions
          if (dualmap.size() > 0)
          {
            for (int m = 0; m < nrow; ++m)
            {
              fac = wgt * sval[m] * sval[k] * jac;
              for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p =
                       dualmap.begin();
                   p != dualmap.end(); ++p)
                dmmap_jk[p->first] -= fac * (p->second)(j, m);
            }
          }

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmderiv(j, 0) * sval[k] * jac;
          for (_CI p = derivsxi[0].begin(); p != derivsxi[0].end(); ++p)
            dmmap_jk[p->first] -= fac * (p->second);

          //          fac = wgt*lmderiv(j, 1)*sval[k]*jac;
          //          for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          //            dmmap_jk[p->first] -= fac*(p->second);

          // (3) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmval[j] * sderiv(k, 0) * jac;
          for (_CI p = derivsxi[0].begin(); p != derivsxi[0].end(); ++p)
            dmmap_jk[p->first] -= fac * (p->second);

          //          fac = wgt*lmval[j]*sderiv(k, 1)*jac;
          //          for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          //            dmmap_jk[p->first] -= fac*(p->second);

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt * lmval[j] * sval[k];
          for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
            dmmap_jk[p->first] -= fac * (p->second);
        }

        // node k is NO boundary node
        else
        {
          // get the correct map as a reference
          std::map<int, double>& ddmap_jk =
              dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[mymrtrnode->Id()];

          // (1) Lin(Phi) - dual shape functions
          if (dualmap.size() > 0)
          {
            for (int m = 0; m < nrow; ++m)
            {
              fac = wgt * sval[m] * sval[k] * jac;
              for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p =
                       dualmap.begin();
                   p != dualmap.end(); ++p)
                ddmap_jk[p->first] += fac * (p->second)(j, m);
            }
          }

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmderiv(j, 0) * sval[k] * jac;
          for (_CI p = derivsxi[0].begin(); p != derivsxi[0].end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);

          //          fac = wgt*lmderiv(j, 1)*sval[k]*jac;
          //          for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          //            ddmap_jk[p->first] += fac*(p->second);

          // (3) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmval[j] * sderiv(k, 0) * jac;
          for (_CI p = derivsxi[0].begin(); p != derivsxi[0].end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);

          //          fac = wgt*lmval[j]*sderiv(k, 1)*jac;
          //          for (_CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          //            ddmap_jk[p->first] += fac*(p->second);

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt * lmval[j] * sval[k];
          for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);
        }
      }  // loop over slave nodes
    }
  }
  else
    dserror("ERROR: unknwon shape function type!");

  return;
}


/*----------------------------------------------------------------------*
 |  Lin D and M matrix entries at GP                         farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_2D_DM_Ele_Lin(int& iter, bool& bound,
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval,
    LINALG::SerialDenseVector& mval, LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseMatrix& mderiv, double& dxdsxi, double& wgt,
    const GEN::pairedvector<int, double>& dmxigp, const GEN::pairedvector<int, double>& derivjac,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap)
{
  DRT::Node** snodes = NULL;
  DRT::Node** mnodes = NULL;

  int nrow = sele.NumNode();
  int ncol = mele.NumNode();

  snodes = sele.Nodes();
  mnodes = mele.Nodes();

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(snodes[iter]);
  if (!mymrtrnode) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  const int sgid = mymrtrnode->Id();

  // standard shape functions
  if (ShapeFcn() == INPAR::MORTAR::shape_standard)
  {
    // integrate LinM
    for (int k = 0; k < ncol; ++k)
    {
      // global master node ID
      int mgid = mnodes[k]->Id();
      double fac = 0.0;


      // get the correct map as a reference
      std::map<int, double>& dmmap_jk =
          dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

      // (1) Lin(Phi) - dual shape functions    --> 0

      // (2) Lin(NSlave) - slave GP coordinates --> 0

      // (3) Lin(NMaster) - master GP coordinates
      fac = wgt * lmval[iter] * mderiv(k, 0) * dxdsxi;
      for (_CI p = dmxigp.begin(); p != dmxigp.end(); ++p) dmmap_jk[p->first] += fac * (p->second);

      // (4) Lin(dsxideta) - segment end coordinates --> 0

      // (5) Lin(dxdsxi) - slave GP Jacobian
      fac = wgt * lmval[iter] * mval[k];
      for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
        dmmap_jk[p->first] += fac * (p->second);

      // (6) Lin(dxdsxi) - slave GP coordinates --> 0
    }  // loop over master nodes

    // integrate LinD
    for (int k = 0; k < nrow; ++k)
    {
      // global slave node ID
      int sgid = snodes[k]->Id();
      double fac = 0.0;

      // get the correct map as a reference
      std::map<int, double>& ddmap_jk =
          dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

      // (1) Lin(Phi) - dual shape functions --> 0

      // (2) Lin(NSlave) - slave GP coordinates --> 0

      // (3) Lin(NSlave) - slave GP coordinates --> 0

      // (4) Lin(dsxideta) - segment end coordinates --> 0

      // (5) Lin(dxdsxi) - slave GP Jacobian
      fac = wgt * lmval[iter] * sval[k];
      for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
        ddmap_jk[p->first] += fac * (p->second);

      // (6) Lin(dxdsxi) - slave GP coordinates --> 0
    }  // loop over slave nodes
  }

  // dual shape functions
  else if (ShapeFcn() == INPAR::MORTAR::shape_dual ||
           ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
  {
    // get the D-map as a reference
    std::map<int, double>& ddmap_jk =
        dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

    // integrate LinM and LinD (NO boundary modification)
    for (int k = 0; k < ncol; ++k)
    {
      // global master node ID
      int mgid = mnodes[k]->Id();
      double fac = 0.0;

      // get the correct map as a reference
      std::map<int, double>& dmmap_jk =
          dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

      // (1) Lin(Phi) - dual shape functions
      for (int m = 0; m < nrow; ++m)
      {
        fac = wgt * sval[m] * mval[k] * dxdsxi;
        for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
             p != dualmap.end(); ++p)
        {
          dmmap_jk[p->first] += fac * (p->second)(iter, m);
          if (!bound) ddmap_jk[p->first] += fac * (p->second)(iter, m);
        }
      }

      // (2) Lin(Phi) - slave GP coordinates --> 0

      // (3) Lin(NMaster) - master GP coordinates
      fac = wgt * lmval(iter, 0) * mderiv(k, 0) * dxdsxi;
      for (_CI p = dmxigp.begin(); p != dmxigp.end(); ++p)
      {
        dmmap_jk[p->first] += fac * (p->second);
        if (!bound) ddmap_jk[p->first] += fac * (p->second);
      }

      // (4) Lin(dsxideta) - segment end coordinates --> 0

      // (5) Lin(dxdsxi) - slave GP Jacobian
      fac = wgt * lmval[iter] * mval[k];
      for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
      {
        dmmap_jk[p->first] += fac * (p->second);
        if (!bound) ddmap_jk[p->first] += fac * (p->second);
      }

      // (6) Lin(dxdsxi) - slave GP coordinates --> 0
    }  // loop over master nodes
  }    // ShapeFcn() switch

  return;
}

/*----------------------------------------------------------------------*
 |  Lin D and M matrix entries at GP                         farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_2D_DM_Lin(int& iter, bool& bound, bool& linlm,
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval,
    LINALG::SerialDenseVector& mval, LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& mderiv,
    LINALG::SerialDenseMatrix& lmderiv, double& jac, double& wgt,
    const std::vector<GEN::pairedvector<int, double>>& dsxigp,
    const std::vector<GEN::pairedvector<int, double>>& dmxigp,
    const GEN::pairedvector<int, double>& derivjac,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap)
{
  DRT::Node** snodes = NULL;
  DRT::Node** mnodes = NULL;
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  snodes = sele.Nodes();
  mnodes = mele.Nodes();

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // **************** no edge modification *****************************
  // (and LinM also for edge node modification case)

  // check for linear LM interpolation in quadratic FE
  if (linlm)
  {
    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(snodes[iter]);
    if (!mymrtrnode) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");
    bool jbound = mymrtrnode->IsOnBound();

    // node j is boundary node
    if (jbound)
    {
      // do nothing as respective D and M entries are zero anyway
    }
    // node j is NO boundary node
    else
    {
      // integrate LinM
      for (int k = 0; k < ncol; ++k)
      {
        // global master node ID
        int mgid = mnodes[k]->Id();
        double fac = 0.0;

        // get the correct map as a reference
        std::map<int, double>& dmmap_jk =
            dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

        // (1) Lin(Phi) - dual shape functions
        // this vanishes here since there are no deformation-dependent dual functions

        // (2) Lin(NSlave) - slave GP coordinates
        fac = wgt * lmderiv(iter, 0) * mval[k] * jac;
        for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);

        // (3) Lin(NMaster) - master GP coordinates
        fac = wgt * lmval[iter] * mderiv(k, 0) * jac;
        for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);

        // (4) Lin(dxdsxi) - slave GP Jacobian
        fac = wgt * lmval[iter] * mval[k];
        for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);
      }  // loop over master nodes

      // integrate LinD
      for (int k = 0; k < nrow; ++k)
      {
        MORTAR::MortarNode* mymrtrnode2 = dynamic_cast<MORTAR::MortarNode*>(snodes[k]);
        if (!mymrtrnode2) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");
        bool kbound = mymrtrnode2->IsOnBound();

        // global master node ID
        int sgid = mymrtrnode2->Id();
        double fac = 0.0;

        // node k is boundary node
        if (kbound)
        {
          // move entry to derivM (with minus sign)
          // get the correct map as a reference
          std::map<int, double>& dmmap_jk =
              dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[sgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmderiv(iter, 0) * sval[k] * jac;
          for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
            dmmap_jk[p->first] -= fac * (p->second);

          // (3) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmval[iter] * sderiv(k, 0) * jac;
          for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
            dmmap_jk[p->first] -= fac * (p->second);

          // (4) Lin(dxdsxi) - slave GP Jacobian
          fac = wgt * lmval[iter] * sval[k];
          for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
            dmmap_jk[p->first] -= fac * (p->second);
        }

        // node k is NO boundary node
        else
        {
          // get the correct map as a reference
          std::map<int, double>& ddmap_jk =
              dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmderiv(iter, 0) * sval[k] * jac;
          for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);

          // (3) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmval[iter] * sderiv(k, 0) * jac;
          for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);

          // (4) Lin(dxdsxi) - slave GP Jacobian
          fac = wgt * lmval[iter] * sval[k];
          for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
            ddmap_jk[p->first] += fac * (p->second);
        }
      }  // loop over slave nodes
    }
  }
  // no linear LM interpolation for quadratic FE
  else
  {
    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(snodes[iter]);
    if (!mymrtrnode) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

    int sgid = mymrtrnode->Id();

    // standard shape functions
    if (ShapeFcn() == INPAR::MORTAR::shape_standard)
    {
      // integrate LinM
      for (int k = 0; k < ncol; ++k)
      {
        // global master node ID
        int mgid = mnodes[k]->Id();
        double fac = 0.0;

        // get the correct map as a reference
        std::map<int, double>& dmmap_jk =
            dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

        // (1) Lin(Phi) - dual shape functions
        // this vanishes here since there are no deformation-dependent dual functions

        // (2) Lin(NSlave) - slave GP coordinates
        fac = wgt * lmderiv(iter, 0) * mval[k] * jac;
        for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);

        // (3) Lin(NMaster) - master GP coordinates
        fac = wgt * lmval[iter] * mderiv(k, 0) * jac;
        for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);

        // (4) Lin(dxdsxi) - slave GP Jacobian
        fac = wgt * lmval[iter] * mval[k];
        for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
          dmmap_jk[p->first] += fac * (p->second);
      }  // loop over master nodes

      // integrate LinD
      for (int k = 0; k < nrow; ++k)
      {
        // global slave node ID
        int sgid = snodes[k]->Id();
        double fac = 0.0;

        // get the correct map as a reference
        std::map<int, double>& ddmap_jk =
            dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

        // (1) Lin(Phi) - dual shape functions
        // this vanishes here since there are no deformation-dependent dual functions

        // (2) Lin(NSlave) - slave GP coordinates
        fac = wgt * lmderiv(iter, 0) * sval[k] * jac;
        for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
          ddmap_jk[p->first] += fac * (p->second);

        // (3) Lin(NSlave) - slave GP coordinates
        fac = wgt * lmval[iter] * sderiv(k, 0) * jac;
        for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
          ddmap_jk[p->first] += fac * (p->second);

        // (4) Lin(dxdsxi) - slave GP Jacobian
        fac = wgt * lmval[iter] * sval[k];
        for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
          ddmap_jk[p->first] += fac * (p->second);
      }  // loop over slave nodes
    }

    // dual shape functions
    else if (ShapeFcn() == INPAR::MORTAR::shape_dual ||
             ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
    {
      // get the D-map as a reference
      std::map<int, double>& ddmap_jk =
          dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

      // integrate LinM and LinD (NO boundary modification)
      for (int k = 0; k < ncol; ++k)
      {
        // global master node ID
        int mgid = mnodes[k]->Id();
        double fac = 0.0;

        // get the correct map as a reference
        std::map<int, double>& dmmap_jk =
            dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

        // (1) Lin(Phi) - dual shape functions
        for (int m = 0; m < nrow; ++m)
        {
          fac = wgt * sval[m] * mval[k] * jac;
          for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
               p != dualmap.end(); ++p)
          {
            dmmap_jk[p->first] += fac * (p->second)(iter, m);
            if (!bound) ddmap_jk[p->first] += fac * (p->second)(iter, m);
          }
        }

        // (2) Lin(Phi) - slave GP coordinates
        fac = wgt * lmderiv(iter, 0) * mval[k] * jac;

        for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
        {
          dmmap_jk[p->first] += fac * (p->second);
          if (!bound) ddmap_jk[p->first] += fac * (p->second);
        }

        // (3) Lin(NMaster) - master GP coordinates
        fac = wgt * lmval(iter, 0) * mderiv(k, 0) * jac;

        for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
        {
          dmmap_jk[p->first] += fac * (p->second);
          if (!bound) ddmap_jk[p->first] += fac * (p->second);
        }

        // (4) Lin(dxdsxi) - slave GP Jacobian
        fac = wgt * lmval[iter] * mval[k];

        for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
        {
          dmmap_jk[p->first] += fac * (p->second);
          if (!bound) ddmap_jk[p->first] += fac * (p->second);
        }
      }  // loop over master nodes
    }    // ShapeFcn() switch
  }
  // compute segment D/M linearization *********************************
  return;
}

/*----------------------------------------------------------------------*
 |  Lin D and M matrix entries at GP - pwlin                 farah 12/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_3D_DM_Quad_pwlin_Lin(int& iter, MORTAR::MortarElement& sele,
    MORTAR::MortarElement& sintele, MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval,
    LINALG::SerialDenseVector& mval, LINALG::SerialDenseVector& lmintval,
    LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& mderiv,
    LINALG::SerialDenseMatrix& lmintderiv, double& wgt, double& jac,
    const std::vector<GEN::pairedvector<int, double>>& dsxigp,
    const std::vector<GEN::pairedvector<int, double>>& dpsxigp,
    const std::vector<GEN::pairedvector<int, double>>& dpmxigp,
    const GEN::pairedvector<int, double>& jacintcellmap)
{
  const int nrow = sele.NumNode();
  const int ncol = mele.NumNode();

  // get slave element nodes themselves
  DRT::Node** sintnodes = sintele.Nodes();

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator CI;

  // **************** no edge modification *****************************
  // (and LinM also for edge node modification case)

  MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(sintnodes[iter]);
  if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

  // integrate LinM
  for (int k = 0; k < ncol; ++k)
  {
    // global master node ID
    int mgid = mele.Nodes()[k]->Id();
    double fac = 0.0;

    // get the correct map as a reference
    std::map<int, double>& dmmap_jk =
        dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

    // (1) Lin(Phi) - dual shape functions
    // this vanishes here since there are no deformation-dependent dual functions

    // (2) Lin(NSlave) - slave GP coordinates
    fac = wgt * lmintderiv(iter, 0) * mval[k] * jac;
    for (CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
      dmmap_jk[p->first] += fac * (p->second);

    fac = wgt * lmintderiv(iter, 1) * mval[k] * jac;
    for (CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
      dmmap_jk[p->first] += fac * (p->second);

    // (3) Lin(NMaster) - master GP coordinates
    fac = wgt * lmintval[iter] * mderiv(k, 0) * jac;
    for (CI p = dpmxigp[0].begin(); p != dpmxigp[0].end(); ++p)
      dmmap_jk[p->first] += fac * (p->second);

    fac = wgt * lmintval[iter] * mderiv(k, 1) * jac;
    for (CI p = dpmxigp[1].begin(); p != dpmxigp[1].end(); ++p)
      dmmap_jk[p->first] += fac * (p->second);

    // (4) Lin(dsxideta) - intcell GP Jacobian
    fac = wgt * lmintval[iter] * mval[k];
    for (CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
      dmmap_jk[p->first] += fac * (p->second);
  }  // loop over master nodes

  // integrate LinD
  for (int k = 0; k < nrow; ++k)
  {
    // global master node ID
    int sgid = sele.Nodes()[k]->Id();
    double fac = 0.0;

    // get the correct map as a reference
    std::map<int, double>& ddmap_jk =
        dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

    // (1) Lin(Phi) - dual shape functions
    // this vanishes here since there are no deformation-dependent dual functions

    // (2) Lin(NSlave) - slave GP coordinates
    fac = wgt * lmintderiv(iter, 0) * sval[k] * jac;
    for (CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
      ddmap_jk[p->first] += fac * (p->second);

    fac = wgt * lmintderiv(iter, 1) * sval[k] * jac;
    for (CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
      ddmap_jk[p->first] += fac * (p->second);

    // (3) Lin(NSlave) - slave GP coordinates
    fac = wgt * lmintval[iter] * sderiv(k, 0) * jac;
    for (CI p = dpsxigp[0].begin(); p != dpsxigp[0].end(); ++p)
      ddmap_jk[p->first] += fac * (p->second);

    fac = wgt * lmintval[iter] * sderiv(k, 1) * jac;
    for (CI p = dpsxigp[1].begin(); p != dpsxigp[1].end(); ++p)
      ddmap_jk[p->first] += fac * (p->second);

    // (4) Lin(dsxideta) - intcell GP Jacobian
    fac = wgt * lmintval[iter] * sval[k];
    for (CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
      ddmap_jk[p->first] += fac * (p->second);
  }  // loop over slave nodes

  return;
}

/*----------------------------------------------------------------------*
 |  Lin D and M matrix entries at GP                         farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_3D_DM_Quad_Lin(bool& duallin, MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval,
    LINALG::SerialDenseVector& svalmod, LINALG::SerialDenseVector& mval,
    LINALG::SerialDenseVector& lmval, LINALG::SerialDenseMatrix& sderiv,
    LINALG::SerialDenseMatrix& mderiv, LINALG::SerialDenseMatrix& lmderiv, double& wgt, double& jac,
    const std::vector<GEN::pairedvector<int, double>>& dpsxigp,
    const std::vector<GEN::pairedvector<int, double>>& dpmxigp,
    const GEN::pairedvector<int, double>& jacintcellmap,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap, bool dualquad3d)
{
  const int nrow = sele.NumNode();
  const int ncol = mele.NumNode();

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();

  std::vector<MORTAR::MortarNode*> smnodes(nrow);
  for (int i = 0; i < nrow; ++i) smnodes[i] = dynamic_cast<MORTAR::MortarNode*>(snodes[i]);

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // **************** no edge modification *****************************
  // (and LinM also for edge node modification case)

  // compute cell D/M linearization ************************************
  // CASE 1: Standard LM shape functions and quadratic interpolation
  if (ShapeFcn() == INPAR::MORTAR::shape_standard && LagMultQuad() == INPAR::MORTAR::lagmult_quad)
  {
    static double fac1 = 0.;
    static double fac2 = 0.;
    // (1) Lin(Phi) - dual shape functions
    // this vanishes here since there are no deformation-dependent dual functions

    // (2) Lin(NSlave) - slave GP coordinates
    // (3) Lin(NSlave) - slave GP coordinates
    for (int d = 0; d < 2; ++d)
      for (_CI p = dpsxigp[d].begin(); p != dpsxigp[d].end(); ++p)
      {
        Epetra_SerialDenseMatrix& dMderiv =
            dynamic_cast<CONTACT::CoElement&>(sele).GetDderiv()[p->first];
        for (int j = 0; j < nrow; ++j)
          for (int k = 0; k < nrow; ++k)
          {
            fac1 = wgt * lmderiv(j, d) * sval[k] * jac;
            fac2 = wgt * lmval[j] * sderiv(k, d) * jac;
            dMderiv(j, k) += (fac1 + fac2) * (p->second);
          }
      }
    // (4) Lin(dsxideta) - intcell GP Jacobian
    for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
    {
      Epetra_SerialDenseMatrix& dMderiv =
          dynamic_cast<CONTACT::CoElement&>(sele).GetDderiv()[p->first];
      for (int j = 0; j < nrow; ++j)
        for (int k = 0; k < nrow; ++k) dMderiv(j, k) += wgt * lmval[j] * sval[k] * (p->second);
    }

    // (1) Lin(Phi) - dual shape functions
    // this vanishes here since there are no deformation-dependent dual functions

    // (2) Lin(NSlave) - slave GP coordinates
    for (int d = 0; d < 2; ++d)
      for (_CI p = dpsxigp[d].begin(); p != dpsxigp[d].end(); ++p)
      {
        Epetra_SerialDenseMatrix& mMderiv =
            dynamic_cast<CONTACT::CoElement&>(sele).GetMderiv()[p->first];
        for (int j = 0; j < nrow; ++j)
          for (int k = 0; k < ncol; ++k)
            mMderiv(j, k) += wgt * lmderiv(j, d) * mval[k] * jac * (p->second);
      }

    // (3) Lin(NMaster) - master GP coordinates
    for (int d = 0; d < 2; ++d)
      for (_CI p = dpmxigp[d].begin(); p != dpmxigp[d].end(); ++p)
      {
        Epetra_SerialDenseMatrix& mMderiv =
            dynamic_cast<CONTACT::CoElement&>(sele).GetMderiv()[p->first];
        for (int j = 0; j < nrow; ++j)
          for (int k = 0; k < ncol; ++k)
            mMderiv(j, k) += wgt * lmval[j] * mderiv(k, d) * jac * (p->second);
      }

    // (4) Lin(dsxideta) - intcell GP Jacobian
    for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
    {
      Epetra_SerialDenseMatrix& mMderiv =
          dynamic_cast<CONTACT::CoElement&>(sele).GetMderiv()[p->first];
      for (int j = 0; j < nrow; ++j)
        for (int k = 0; k < ncol; ++k) mMderiv(j, k) += wgt * lmval[j] * mval[k] * (p->second);
    }
  }

  // CASE 2: Standard LM shape functions and linear interpolation
  // CASE 5: dual LM shape functions and linear interpolation
  // (this has to be treated seperately here for LinDM because of bound)
  else if ((ShapeFcn() == INPAR::MORTAR::shape_standard ||
               ShapeFcn() == INPAR::MORTAR::shape_dual ||
               ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) &&
           LagMultQuad() == INPAR::MORTAR::lagmult_lin)
  {
    // integrate LinD
    for (int j = 0; j < nrow; ++j)
    {
      MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(snodes[j]);
      if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

      // node j is boundary node
      if (mymrtrnode->SetBound())
      {
        // do nothing as respective D and M entries are zero anyway
      }

      // node j is NO boundary node
      else
      {
        // integrate LinM
        for (int k = 0; k < ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          static double fac = 0.0;

          // get the correct map as a reference
          std::map<int, double>& dmmap_jk =
              dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

          // (1) Lin(Phi) - dual shape functions
          if (duallin)
          {
            typedef GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator _CIM;
            for (_CIM p = dualmap.begin(); p != dualmap.end(); ++p)
            {
              double lmderiv_d = 0.0;
              for (int m = 0; m < nrow; ++m) lmderiv_d += (p->second)(j, m) * sval[m];
              fac = wgt * lmderiv_d * mval[k] * jac;
              dmmap_jk[p->first] += fac;
            }
          }

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt * lmderiv(j, 0) * mval[k] * jac;
          for (_CI p = dpsxigp[0].begin(); p != dpsxigp[0].end(); ++p)
            dmmap_jk[p->first] += fac * (p->second);

          fac = wgt * lmderiv(j, 1) * mval[k] * jac;
          for (_CI p = dpsxigp[1].begin(); p != dpsxigp[1].end(); ++p)
            dmmap_jk[p->first] += fac * (p->second);

          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt * lmval[j] * mderiv(k, 0) * jac;
          for (_CI p = dpmxigp[0].begin(); p != dpmxigp[0].end(); ++p)
            dmmap_jk[p->first] += fac * (p->second);

          fac = wgt * lmval[j] * mderiv(k, 1) * jac;
          for (_CI p = dpmxigp[1].begin(); p != dpmxigp[1].end(); ++p)
            dmmap_jk[p->first] += fac * (p->second);

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt * lmval[j] * mval[k];
          for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
            dmmap_jk[p->first] += fac * (p->second);
        }  // loop over master nodes

        for (int k = 0; k < nrow; ++k)
        {
          MORTAR::MortarNode* mymrtrnode2 = dynamic_cast<MORTAR::MortarNode*>(snodes[k]);
          if (!mymrtrnode2) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

          // global master node ID
          int sgid = mymrtrnode2->Id();
          static double fac = 0.0;

          // node k is boundary node
          if (mymrtrnode2->IsOnBound())
          {
            // move entry to derivM (with minus sign)
            // get the correct map as a reference
            std::map<int, double>& dmmap_jk =
                dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[sgid];

            // (1) Lin(Phi) - dual shape functions
            if (duallin)
            {
              typedef GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator _CIM;
              fac = wgt * sval[k] * jac;
              for (_CIM p = dualmap.begin(); p != dualmap.end(); ++p)
                for (int m = 0; m < nrow; ++m)
                  dmmap_jk[p->first] -= fac * (p->second)(j, m) * sval[m];

              for (_CIM p = dualmap.begin(); p != dualmap.end(); ++p)
              {
                double lmderiv_d = 0.0;
                for (int m = 0; m < nrow; ++m) lmderiv_d += (p->second)(j, m) * sval[m];
                fac = wgt * lmderiv_d * mval[k] * jac;
                dmmap_jk[p->first] -= fac;
              }
            }

            // (2) Lin(NSlave) - slave GP coordinates
            fac = wgt * lmderiv(j, 0) * sval[k] * jac;
            for (_CI p = dpsxigp[0].begin(); p != dpsxigp[0].end(); ++p)
              dmmap_jk[p->first] -= fac * (p->second);

            fac = wgt * lmderiv(j, 1) * sval[k] * jac;
            for (_CI p = dpsxigp[1].begin(); p != dpsxigp[1].end(); ++p)
              dmmap_jk[p->first] -= fac * (p->second);

            // (3) Lin(NSlave) - slave GP coordinates
            fac = wgt * lmval[j] * sderiv(k, 0) * jac;
            for (_CI p = dpsxigp[0].begin(); p != dpsxigp[0].end(); ++p)
              dmmap_jk[p->first] -= fac * (p->second);

            fac = wgt * lmval[j] * sderiv(k, 1) * jac;
            for (_CI p = dpsxigp[1].begin(); p != dpsxigp[1].end(); ++p)
              dmmap_jk[p->first] -= fac * (p->second);

            // (4) Lin(dsxideta) - intcell GP Jacobian
            fac = wgt * lmval[j] * sval[k];
            for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
              dmmap_jk[p->first] -= fac * (p->second);
          }

          // node k is NO boundary node
          else
          {
            // get the correct map as a reference
            std::map<int, double>& ddmap_jk =
                dynamic_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

            // (1) Lin(Phi) - dual shape functions
            if (duallin)
            {
              typedef GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator _CIM;
              fac = wgt * sval[k] * jac;
              for (_CIM p = dualmap.begin(); p != dualmap.end(); ++p)
                for (int m = 0; m < nrow; ++m)
                  ddmap_jk[p->first] += fac * (p->second)(j, m) * sval[m];

              for (_CIM p = dualmap.begin(); p != dualmap.end(); ++p)
              {
                double lmderiv_d = 0.0;
                for (int m = 0; m < nrow; ++m) lmderiv_d += (p->second)(j, m) * sval[m];
                fac = wgt * lmderiv_d * mval[k] * jac;
                ddmap_jk[p->first] += fac;
              }
            }

            // (2) Lin(NSlave) - slave GP coordinates
            fac = wgt * lmderiv(j, 0) * sval[k] * jac;
            for (_CI p = dpsxigp[0].begin(); p != dpsxigp[0].end(); ++p)
              ddmap_jk[p->first] += fac * (p->second);

            fac = wgt * lmderiv(j, 1) * sval[k] * jac;
            for (_CI p = dpsxigp[1].begin(); p != dpsxigp[1].end(); ++p)
              ddmap_jk[p->first] += fac * (p->second);

            // (3) Lin(NSlave) - slave GP coordinates
            fac = wgt * lmval[j] * sderiv(k, 0) * jac;
            for (_CI p = dpsxigp[0].begin(); p != dpsxigp[0].end(); ++p)
              ddmap_jk[p->first] += fac * (p->second);

            fac = wgt * lmval[j] * sderiv(k, 1) * jac;
            for (_CI p = dpsxigp[1].begin(); p != dpsxigp[1].end(); ++p)
              ddmap_jk[p->first] += fac * (p->second);

            // (4) Lin(dsxideta) - intcell GP Jacobian
            fac = wgt * lmval[j] * sval[k];
            for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
              ddmap_jk[p->first] += fac * (p->second);
          }
        }  // loop over slave nodes
      }
    }
  }
  // CASE 4: Dual LM shape functions and quadratic interpolation
  else if ((ShapeFcn() == INPAR::MORTAR::shape_dual ||
               ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) &&
           LagMultQuad() == INPAR::MORTAR::lagmult_quad)
  {
    double fac = 0.;
    // (1) Lin(Phi) - dual shape functions
    if (duallin)
    {
      for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
           p != dualmap.end(); ++p)
      {
        Epetra_SerialDenseMatrix& dDderiv =
            dynamic_cast<CONTACT::CoElement&>(sele).GetDderiv()[p->first];
        Epetra_SerialDenseMatrix& dMderiv =
            dynamic_cast<CONTACT::CoElement&>(sele).GetMderiv()[p->first];
        for (int m = 0; m < nrow; ++m)
          for (int k = 0; k < ncol; ++k)
          {
            if (dualquad3d)
              fac = wgt * svalmod[m] * mval[k] * jac;
            else
              fac = wgt * sval[m] * mval[k] * jac;
            for (int j = 0; j < nrow; ++j)
            {
              dMderiv(j, k) += fac * (p->second)(j, m);
              dDderiv(j, j) += fac * (p->second)(j, m);
            }
          }
      }
    }

    // (2) Lin(Phi) - slave GP coordinates
    for (int d = 0; d < 2; ++d)
      for (_CI p = dpsxigp[d].begin(); p != dpsxigp[d].end(); ++p)
      {
        Epetra_SerialDenseMatrix& dDderiv =
            dynamic_cast<CONTACT::CoElement&>(sele).GetDderiv()[p->first];
        Epetra_SerialDenseMatrix& dMderiv =
            dynamic_cast<CONTACT::CoElement&>(sele).GetMderiv()[p->first];
        for (int j = 0; j < nrow; ++j)
          for (int k = 0; k < ncol; ++k)
          {
            fac = wgt * lmderiv(j, d) * mval[k] * jac;
            dMderiv(j, k) += fac * (p->second);
            dDderiv(j, j) += fac * (p->second);
          }
      }

    // (3) Lin(NMaster) - master GP coordinates
    for (int d = 0; d < 2; ++d)
      for (_CI p = dpmxigp[d].begin(); p != dpmxigp[d].end(); ++p)
      {
        Epetra_SerialDenseMatrix& dDderiv =
            dynamic_cast<CONTACT::CoElement&>(sele).GetDderiv()[p->first];
        Epetra_SerialDenseMatrix& dMderiv =
            dynamic_cast<CONTACT::CoElement&>(sele).GetMderiv()[p->first];
        for (int j = 0; j < nrow; ++j)
          for (int k = 0; k < ncol; ++k)
          {
            fac = wgt * lmval[j] * mderiv(k, d) * jac;
            dMderiv(j, k) += fac * (p->second);
            dDderiv(j, j) += fac * (p->second);
          }
      }

    // (4) Lin(dsxideta) - intcell GP Jacobian
    for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
    {
      Epetra_SerialDenseMatrix& dDderiv =
          dynamic_cast<CONTACT::CoElement&>(sele).GetDderiv()[p->first];
      Epetra_SerialDenseMatrix& dMderiv =
          dynamic_cast<CONTACT::CoElement&>(sele).GetMderiv()[p->first];
      for (int j = 0; j < nrow; ++j)
        for (int k = 0; k < ncol; ++k)
        {
          fac = wgt * lmval[j] * mval[k];
          dMderiv(j, k) += fac * (p->second);
          dDderiv(j, j) += fac * (p->second);
        }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Lin D and M matrix entries at GP                         farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_3D_DM_Lin(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
    LINALG::SerialDenseVector& lmval, LINALG::SerialDenseMatrix& sderiv,
    LINALG::SerialDenseMatrix& mderiv, LINALG::SerialDenseMatrix& lmderiv, double& wgt, double& jac,
    std::vector<GEN::pairedvector<int, double>>& dsxigp,
    std::vector<GEN::pairedvector<int, double>>& dmxigp,
    GEN::pairedvector<int, double>& jacintcellmap,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap)
{
  const int nrow = sele.NumNode();
  const int ncol = mele.NumNode();

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // standard shape functions
  if (ShapeFcn() == INPAR::MORTAR::shape_standard)
  {
    // integrate LinM
    // (1) Lin(Phi) - dual shape functions
    // this vanishes here since there are no deformation-dependent dual functions

    // (2) Lin(NSlave) - slave GP coordinates
    for (int d = 0; d < 2; ++d)
      for (_CI p = dsxigp[d].begin(); p != dsxigp[d].end(); ++p)
      {
        Epetra_SerialDenseMatrix& dMderiv =
            dynamic_cast<CONTACT::CoElement&>(sele).GetMderiv()[p->first];
        for (int j = 0; j < nrow; ++j)
          for (int k = 0; k < ncol; ++k)
            dMderiv(j, k) += wgt * lmderiv(j, d) * mval[k] * jac * (p->second);
      }

    // (3) Lin(NMaster) - master GP coordinates
    for (int d = 0; d < 2; ++d)
      for (_CI p = dmxigp[d].begin(); p != dmxigp[d].end(); ++p)
      {
        Epetra_SerialDenseMatrix& dMderiv =
            dynamic_cast<CONTACT::CoElement&>(sele).GetMderiv()[p->first];
        for (int j = 0; j < nrow; ++j)
          for (int k = 0; k < ncol; ++k)
            dMderiv(j, k) += wgt * lmval[j] * mderiv(k, d) * jac * (p->second);
      }

    // (4) Lin(dsxideta) - intcell GP Jacobian
    for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
    {
      Epetra_SerialDenseMatrix& dMderiv =
          dynamic_cast<CONTACT::CoElement&>(sele).GetMderiv()[p->first];
      for (int j = 0; j < nrow; ++j)
        for (int k = 0; k < ncol; ++k) dMderiv(j, k) += wgt * lmval[j] * mval[k] * (p->second);
    }

    // integrate LinD
    // (1) Lin(Phi) - dual shape functions
    // this vanishes here since there are no deformation-dependent dual functions

    // (2) Lin(NSlave) - slave GP coordinates
    // (3) Lin(NSlave) - slave GP coordinates
    for (int d = 0; d < 2; ++d)
      for (_CI p = dsxigp[d].begin(); p != dsxigp[d].end(); ++p)
      {
        Epetra_SerialDenseMatrix& dDderiv =
            dynamic_cast<CONTACT::CoElement&>(sele).GetDderiv()[p->first];
        for (int j = 0; j < nrow; ++j)
          for (int k = 0; k < nrow; ++k)
            dDderiv(j, k) +=
                (wgt * lmderiv(j, d) * sval[k] * jac + wgt * lmval[j] * sderiv(k, d) * jac) *
                (p->second);
      }

    // (4) Lin(dsxideta) - intcell GP Jacobian
    for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
    {
      Epetra_SerialDenseMatrix& dDderiv =
          dynamic_cast<CONTACT::CoElement&>(sele).GetDderiv()[p->first];
      for (int j = 0; j < nrow; ++j)
        for (int k = 0; k < nrow; ++k) dDderiv(j, k) += wgt * lmval[j] * sval[k] * (p->second);
    }
  }

  // dual shape functions
  else if (ShapeFcn() == INPAR::MORTAR::shape_dual ||
           ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
  {
    double fac = 0.;

    // (1) Lin(Phi) - dual shape functions
    if (dualmap.size() > 0)
      for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
           p != dualmap.end(); ++p)
      {
        Epetra_SerialDenseMatrix& dDderiv =
            dynamic_cast<CONTACT::CoElement&>(sele).GetDderiv()[p->first];
        Epetra_SerialDenseMatrix& dMderiv =
            dynamic_cast<CONTACT::CoElement&>(sele).GetMderiv()[p->first];
        for (int j = 0; j < nrow; ++j)
          for (int k = 0; k < ncol; ++k)
            for (int m = 0; m < nrow; ++m)
            {
              fac = wgt * sval[m] * mval[k] * jac * (p->second)(j, m);
              dDderiv(j, j) += fac;
              dMderiv(j, k) += fac;
            }
      }

    // (2) Lin(Phi) - slave GP coordinates
    for (int d = 0; d < 2; ++d)
      for (_CI p = dsxigp[d].begin(); p != dsxigp[d].end(); ++p)
      {
        Epetra_SerialDenseMatrix& dDderiv =
            dynamic_cast<CONTACT::CoElement&>(sele).GetDderiv()[p->first];
        Epetra_SerialDenseMatrix& dMderiv =
            dynamic_cast<CONTACT::CoElement&>(sele).GetMderiv()[p->first];
        for (int j = 0; j < nrow; ++j)
          for (int k = 0; k < ncol; ++k)
          {
            fac = wgt * lmderiv(j, d) * mval[k] * jac * (p->second);
            dDderiv(j, j) += fac;
            dMderiv(j, k) += fac;
          }
      }

    // (3) Lin(NMaster) - master GP coordinates
    for (int d = 0; d < 2; ++d)
      for (_CI p = dmxigp[d].begin(); p != dmxigp[d].end(); ++p)
      {
        Epetra_SerialDenseMatrix& dDderiv =
            dynamic_cast<CONTACT::CoElement&>(sele).GetDderiv()[p->first];
        Epetra_SerialDenseMatrix& dMderiv =
            dynamic_cast<CONTACT::CoElement&>(sele).GetMderiv()[p->first];
        for (int j = 0; j < nrow; ++j)
          for (int k = 0; k < ncol; ++k)
          {
            fac = wgt * lmval[j] * mderiv(k, d) * jac * (p->second);
            dDderiv(j, j) += fac;
            dMderiv(j, k) += fac;
          }
      }

    // (4) Lin(dsxideta) - intcell GP Jacobian
    for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
    {
      Epetra_SerialDenseMatrix& dDderiv =
          dynamic_cast<CONTACT::CoElement&>(sele).GetDderiv()[p->first];
      Epetra_SerialDenseMatrix& dMderiv =
          dynamic_cast<CONTACT::CoElement&>(sele).GetMderiv()[p->first];
      for (int j = 0; j < nrow; ++j)
        for (int k = 0; k < ncol; ++k)
        {
          fac = wgt * lmval[j] * mval[k] * (p->second);
          dDderiv(j, j) += fac;
          dMderiv(j, k) += fac;
        }
    }
  }  // ShapeFcn() switch
  // compute segment D/M linearization *********************************
  return;
}

/*----------------------------------------------------------------------*
 |  Compute entries for D2 matrix at GP                      farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_D2(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    LINALG::SerialDenseVector& lm2val, LINALG::SerialDenseVector& m2val, double& jac, double& wgt,
    const Epetra_Comm& comm)
{
  int ncol = mele.NumNode();
  int ndof = Dim();

  // get slave element nodes themselves
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  if (ShapeFcn() == INPAR::MORTAR::shape_dual || ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
  {
    for (int j = 0; j < ncol; ++j)
    {
      CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(mnodes[j]);

      // IMPORTANT: assembling to node is only allowed for master nodes
      //            associated to owned slave elements. This results
      //            to an unique entry distribution!
      if (sele.Owner() == comm.MyPID())
      {
        for (int jdof = 0; jdof < ndof; ++jdof)
        {
          for (int k = 0; k < ncol; ++k)
          {
            CONTACT::FriNode* mnode = dynamic_cast<CONTACT::FriNode*>(mnodes[k]);

            for (int kdof = 0; kdof < ndof; ++kdof)
            {
              int col = mnode->Dofs()[kdof];

              // multiply the two shape functions
              double prod = lm2val[j] * m2val[k] * jac * wgt;

              if ((jdof == kdof) and (j == k))
              {
                if (abs(prod) > MORTARINTTOL) cnode->InvolvedM() = true;
                if (abs(prod) > MORTARINTTOL) cnode->AddD2Value(jdof, col, prod);
              }
            }
          }
        }
      }
    }
  }
  else
    dserror("Both-sided wear just for dual shape functions!");

  return;
}



/*----------------------------------------------------------------------*
 |  Compute wear at GP (for expl/impl algor.)                farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_2D_Wear(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval, LINALG::SerialDenseMatrix& sderiv,
    LINALG::SerialDenseVector& mval, LINALG::SerialDenseMatrix& mderiv,
    LINALG::SerialDenseVector& lmval, LINALG::SerialDenseMatrix& lmderiv,
    Teuchos::RCP<LINALG::SerialDenseMatrix> lagmult, double* gpn, double& jac, double& wgt,
    double* jumpval, double* wearval, GEN::pairedvector<int, double>& dsliptmatrixgp,
    GEN::pairedvector<int, double>& dweargp,
    const std::vector<GEN::pairedvector<int, double>>& dsxigp,
    const std::vector<GEN::pairedvector<int, double>>& dmxigp,
    const std::vector<GEN::pairedvector<int, double>>& dnmap_unit,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap)
{
  const int nrow = sele.NumNode();
  const int ncol = mele.NumNode();
  const int ndof = Dim();

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();

  int linsize = 0;
  for (int i = 0; i < sele.NumNode(); ++i)
    linsize += dynamic_cast<CoNode*>(sele.Nodes()[i])->GetLinsize();
  linsize = linsize * 2;

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  //***********************************************************************
  // Here, the tangential relative slip increment is used and NOT the
  // nodal weighted tangential relative slip increment !!!
  // The reason for that is the slip for wear calculation is written
  // within the integral --> no double weighting allowed !
  // The wearcoefficient is not included in this calculation
  //***********************************************************************
  // for wearval
  double gpt[2] = {0.0, 0.0};
  double gplm[2] = {0.0, 0.0};
  double sgpjump[2] = {0.0, 0.0};
  double mgpjump[2] = {0.0, 0.0};
  double jump[2] = {0.0, 0.0};

  // for linearization
  double lm_lin = 0.0;
  double lengtht = 0.0;

  for (int i = 0; i < nrow; ++i)
  {
    CONTACT::CoNode* myconode = dynamic_cast<CONTACT::CoNode*>(snodes[i]);

    // nodal tangent interpolation
    gpt[0] += sval[i] * myconode->CoData().txi()[0];
    gpt[1] += sval[i] * myconode->CoData().txi()[1];

    // delta D
    sgpjump[0] += sval[i] * (sele.GetNodalCoords(0, i) - (sele.GetNodalCoordsOld(0, i)));
    sgpjump[1] += sval[i] * (sele.GetNodalCoords(1, i) - (sele.GetNodalCoordsOld(1, i)));

    // LM interpolation
    gplm[0] += lmval[i] * ((*lagmult)(0, i));
    gplm[1] += lmval[i] * ((*lagmult)(1, i));
  }

  // normalize interpolated GP tangent back to length 1.0 !!!
  lengtht = sqrt(gpt[0] * gpt[0] + gpt[1] * gpt[1]);
  if (abs(lengtht) < 1.0e-12) dserror("ERROR: IntegrateAndDerivSegment: Divide by zero!");

  for (int i = 0; i < Dim(); i++) gpt[i] /= lengtht;

  // interpolation of master GP jumps (relative displacement increment)
  for (int i = 0; i < ncol; ++i)
  {
    mgpjump[0] += mval[i] * (mele.GetNodalCoords(0, i) - mele.GetNodalCoordsOld(0, i));
    mgpjump[1] += mval[i] * (mele.GetNodalCoords(1, i) - mele.GetNodalCoordsOld(1, i));
  }

  // jump
  jump[0] = sgpjump[0] - mgpjump[0];
  jump[1] = sgpjump[1] - mgpjump[1];

  // evaluate wear
  // normal contact stress -- normal LM value
  for (int i = 0; i < Dim(); ++i)
  {
    wearval[0] += gpn[i] * gplm[i];
    lm_lin += gpn[i] * gplm[i];  // required for linearization
  }

  // value of relative tangential jump
  for (int i = 0; i < Dim(); ++i) jumpval[0] += gpt[i] * jump[i];

  // steady state slip
  if (sswear_) jumpval[0] = ssslip_;

  // no jump --> no wear
  if (abs(jumpval[0]) < 1e-12) return;

  // product
  // use non-abs value for implicit-wear algorithm
  // just for simple linear. maybe we change this in future
  if (wearimpl_ and WearType() != INPAR::WEAR::wear_primvar)
    wearval[0] = (wearval[0]) * abs(jumpval[0]);
  else
    wearval[0] = abs(wearval[0]) * abs(jumpval[0]);

  // compute segment wear vector ***************************************
  // nrow represents the slave side dofs !!!
  for (int j = 0; j < nrow; ++j)
  {
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(snodes[j]);

    double prod = 0.0;
    if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
      prod = sval[j] * wearval[0] * jac * wgt;
    else
      prod = lmval[j] * wearval[0] * jac * wgt;

    // add current Gauss point's contribution to wseg
    cnode->AddDeltaWeightedWearValue(prod);
  }

  //****************************************************************
  //   linearization for implicit algorithms
  //****************************************************************
  if (wearimpl_ || WearType() == INPAR::WEAR::wear_primvar)
  {
    // evaluate the GP wear function derivatives
    GEN::pairedvector<int, double> ddualgp_x(ndof * ncol + linsize);
    GEN::pairedvector<int, double> ddualgp_y(ndof * ncol + linsize);

    GEN::pairedvector<int, double> ddualgp_x_sxi(ndof * ncol + linsize);
    GEN::pairedvector<int, double> ddualgp_y_sxi(ndof * ncol + linsize);

    // lin. abs(x) = x/abs(x) * lin x.
    double xabsx = (jumpval[0] / abs(jumpval[0])) * lm_lin;
    double xabsxT = (jumpval[0] / abs(jumpval[0]));

    // **********************************************************************
    // (1) Lin of normal for LM -- deriv normal maps from weighted gap lin.
    for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
      dweargp[p->first] += abs(jumpval[0]) * gplm[0] * (p->second);

    for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
      dweargp[p->first] += abs(jumpval[0]) * gplm[1] * (p->second);

    // **********************************************************************
    // (2) Lin. of dual shape function of LM mult.
    for (int i = 0; i < nrow; ++i)
    {
      for (int m = 0; m < nrow; ++m)
      {
        for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
             p != dualmap.end(); ++p)
        {
          ddualgp_x[p->first] += ((*lagmult)(0, m)) * sval[m] * (p->second)(i, m);
          ddualgp_y[p->first] += ((*lagmult)(1, m)) * sval[m] * (p->second)(i, m);
        }
      }
    }
    for (_CI p = ddualgp_x.begin(); p != ddualgp_x.end(); ++p)
      dweargp[p->first] += abs(jumpval[0]) * gpn[0] * (p->second);
    for (_CI p = ddualgp_y.begin(); p != ddualgp_y.end(); ++p)
      dweargp[p->first] += abs(jumpval[0]) * gpn[1] * (p->second);

    // LM deriv
    for (int i = 0; i < nrow; ++i)
    {
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
      {
        ddualgp_x_sxi[p->first] += ((*lagmult)(0, i)) * lmderiv(i, 0) * (p->second);
        ddualgp_y_sxi[p->first] += ((*lagmult)(1, i)) * lmderiv(i, 0) * (p->second);
      }
    }
    for (_CI p = ddualgp_x_sxi.begin(); p != ddualgp_x_sxi.end(); ++p)
      dweargp[p->first] += abs(jumpval[0]) * gpn[0] * (p->second);
    for (_CI p = ddualgp_y_sxi.begin(); p != ddualgp_y_sxi.end(); ++p)
      dweargp[p->first] += abs(jumpval[0]) * gpn[1] * (p->second);


    // **********************************************************************
    // (3) absolute incremental slip linearization:
    // (a) build directional derivative of slave GP tagent (non-unit)
    GEN::pairedvector<int, double> dmap_txsl_gp(ndof * ncol + linsize);
    GEN::pairedvector<int, double> dmap_tysl_gp(ndof * ncol + linsize);

    for (int i = 0; i < nrow; ++i)
    {
      CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(snodes[i]);

      GEN::pairedvector<int, double>& dmap_txsl_i =
          dynamic_cast<CONTACT::CoNode*>(cnode)->CoData().GetDerivTxi()[0];
      GEN::pairedvector<int, double>& dmap_tysl_i =
          dynamic_cast<CONTACT::CoNode*>(cnode)->CoData().GetDerivTxi()[1];

      for (_CI p = dmap_txsl_i.begin(); p != dmap_txsl_i.end(); ++p)
        dmap_txsl_gp[p->first] += sval[i] * (p->second);
      for (_CI p = dmap_tysl_i.begin(); p != dmap_tysl_i.end(); ++p)
        dmap_tysl_gp[p->first] += sval[i] * (p->second);

      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
      {
        double valx = sderiv(i, 0) * dynamic_cast<CONTACT::CoNode*>(cnode)->CoData().txi()[0];
        dmap_txsl_gp[p->first] += valx * (p->second);
        double valy = sderiv(i, 0) * dynamic_cast<CONTACT::CoNode*>(cnode)->CoData().txi()[1];
        dmap_tysl_gp[p->first] += valy * (p->second);
      }
    }

    // (b) build directional derivative of slave GP tagent (unit)
    GEN::pairedvector<int, double> dmap_txsl_gp_unit(ndof * ncol + linsize);
    GEN::pairedvector<int, double> dmap_tysl_gp_unit(ndof * ncol + linsize);

    const double ll = lengtht * lengtht;
    const double linv = 1.0 / lengtht;
    const double lllinv = 1.0 / (lengtht * lengtht * lengtht);
    const double sxsx = gpt[0] * gpt[0] * ll;
    const double sxsy = gpt[0] * gpt[1] * ll;
    const double sysy = gpt[1] * gpt[1] * ll;

    for (_CI p = dmap_txsl_gp.begin(); p != dmap_txsl_gp.end(); ++p)
    {
      dmap_txsl_gp_unit[p->first] += linv * (p->second);
      dmap_txsl_gp_unit[p->first] -= lllinv * sxsx * (p->second);
      dmap_tysl_gp_unit[p->first] -= lllinv * sxsy * (p->second);
    }

    for (_CI p = dmap_tysl_gp.begin(); p != dmap_tysl_gp.end(); ++p)
    {
      dmap_tysl_gp_unit[p->first] += linv * (p->second);
      dmap_tysl_gp_unit[p->first] -= lllinv * sysy * (p->second);
      dmap_txsl_gp_unit[p->first] -= lllinv * sxsy * (p->second);
    }

    // add tangent lin. to dweargp
    for (_CI p = dmap_txsl_gp_unit.begin(); p != dmap_txsl_gp_unit.end(); ++p)
      dweargp[p->first] += xabsx * jump[0] * (p->second);

    for (_CI p = dmap_tysl_gp_unit.begin(); p != dmap_tysl_gp_unit.end(); ++p)
      dweargp[p->first] += xabsx * jump[1] * (p->second);

    // add tangent lin. to slip linearization for wear Tmatrix
    for (_CI p = dmap_txsl_gp_unit.begin(); p != dmap_txsl_gp_unit.end(); ++p)
      dsliptmatrixgp[p->first] += xabsxT * jump[0] * (p->second);

    for (_CI p = dmap_tysl_gp_unit.begin(); p != dmap_tysl_gp_unit.end(); ++p)
      dsliptmatrixgp[p->first] += xabsxT * jump[1] * (p->second);

    // **********************************************************************
    // (c) build directional derivative of jump
    GEN::pairedvector<int, double> dmap_slcoord_gp_x(ndof * ncol + linsize);
    GEN::pairedvector<int, double> dmap_slcoord_gp_y(ndof * ncol + linsize);

    GEN::pairedvector<int, double> dmap_mcoord_gp_x(ndof * ncol + linsize);
    GEN::pairedvector<int, double> dmap_mcoord_gp_y(ndof * ncol + linsize);

    GEN::pairedvector<int, double> dmap_coord_x(ndof * ncol + linsize);
    GEN::pairedvector<int, double> dmap_coord_y(ndof * ncol + linsize);

    // lin slave part -- sxi
    for (int i = 0; i < nrow; ++i)
    {
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
      {
        double valx = sderiv(i, 0) * (sele.GetNodalCoords(0, i) - (sele.GetNodalCoordsOld(0, i)));
        dmap_slcoord_gp_x[p->first] += valx * (p->second);
        double valy = sderiv(i, 0) * (sele.GetNodalCoords(1, i) - (sele.GetNodalCoordsOld(1, i)));
        dmap_slcoord_gp_y[p->first] += valy * (p->second);
      }
    }

    // lin master part -- mxi
    for (int i = 0; i < ncol; ++i)
    {
      for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
      {
        double valx = mderiv(i, 0) * (mele.GetNodalCoords(0, i) - (mele.GetNodalCoordsOld(0, i)));
        dmap_mcoord_gp_x[p->first] += valx * (p->second);
        double valy = mderiv(i, 0) * (mele.GetNodalCoords(1, i) - (mele.GetNodalCoordsOld(1, i)));
        dmap_mcoord_gp_y[p->first] += valy * (p->second);
      }
    }

    // deriv slave x-coords
    for (int i = 0; i < nrow; ++i)
    {
      CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(snodes[i]);

      dmap_slcoord_gp_x[snode->Dofs()[0]] += sval[i];
      dmap_slcoord_gp_y[snode->Dofs()[1]] += sval[i];
    }
    // deriv master x-coords
    for (int i = 0; i < ncol; ++i)
    {
      CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mnodes[i]);

      dmap_mcoord_gp_x[mnode->Dofs()[0]] += mval[i];
      dmap_mcoord_gp_y[mnode->Dofs()[1]] += mval[i];
    }

    // slave: add to jumplin
    for (_CI p = dmap_slcoord_gp_x.begin(); p != dmap_slcoord_gp_x.end(); ++p)
      dmap_coord_x[p->first] += (p->second);
    for (_CI p = dmap_slcoord_gp_y.begin(); p != dmap_slcoord_gp_y.end(); ++p)
      dmap_coord_y[p->first] += (p->second);

    // master: add to jumplin
    for (_CI p = dmap_mcoord_gp_x.begin(); p != dmap_mcoord_gp_x.end(); ++p)
      dmap_coord_x[p->first] -= (p->second);
    for (_CI p = dmap_mcoord_gp_y.begin(); p != dmap_mcoord_gp_y.end(); ++p)
      dmap_coord_y[p->first] -= (p->second);

    // add to dweargp
    for (_CI p = dmap_coord_x.begin(); p != dmap_coord_x.end(); ++p)
      dweargp[p->first] += xabsx * gpt[0] * (p->second);

    for (_CI p = dmap_coord_y.begin(); p != dmap_coord_y.end(); ++p)
      dweargp[p->first] += xabsx * gpt[1] * (p->second);

    // add tangent lin. to slip linearization for wear Tmatrix
    for (_CI p = dmap_coord_x.begin(); p != dmap_coord_x.end(); ++p)
      dsliptmatrixgp[p->first] += xabsxT * gpt[0] * (p->second);

    for (_CI p = dmap_coord_y.begin(); p != dmap_coord_y.end(); ++p)
      dsliptmatrixgp[p->first] += xabsxT * gpt[1] * (p->second);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Compute wear at GP (for expl/impl algor.)                farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_3D_Wear(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval, LINALG::SerialDenseMatrix& sderiv,
    LINALG::SerialDenseVector& mval, LINALG::SerialDenseMatrix& mderiv,
    LINALG::SerialDenseVector& lmval, LINALG::SerialDenseMatrix& lmderiv,
    Teuchos::RCP<LINALG::SerialDenseMatrix> lagmult, double* gpn, double& jac, double& wgt,
    double* jumpval, double* wearval, GEN::pairedvector<int, double>& dsliptmatrixgp,
    GEN::pairedvector<int, double>& dweargp,
    const std::vector<GEN::pairedvector<int, double>>& dsxigp,
    const std::vector<GEN::pairedvector<int, double>>& dmxigp,
    const std::vector<GEN::pairedvector<int, double>>& dnmap_unit,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap)
{
  const int nrow = sele.NumNode();
  const int ncol = mele.NumNode();
  const int ndof = Dim();

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  //***********************************************************************
  // Here, the tangential relative slip increment is used and NOT the
  // nodal weighted tangential relative slip increment !!!
  // The reason for that is the slip for wear calculation is written
  // within the integral --> no double weighting allowed !
  // The wearcoefficient is not included in this calculation
  //***********************************************************************

  // wear-lin
  LINALG::Matrix<3, 1> jump;
  LINALG::Matrix<3, 1> jumptan;
  LINALG::Matrix<3, 1> lm;
  LINALG::Matrix<3, 1> lmtan;
  LINALG::Matrix<3, 3> tanplane;

  double gplm[3] = {0.0, 0.0, 0.0};
  double lm_lin = 0.0;

  // tangent plane
  tanplane(0, 0) = 1.0 - (gpn[0] * gpn[0]);
  tanplane(0, 1) = -(gpn[0] * gpn[1]);
  tanplane(0, 2) = -(gpn[0] * gpn[2]);
  tanplane(1, 0) = -(gpn[1] * gpn[0]);
  tanplane(1, 1) = 1.0 - (gpn[1] * gpn[1]);
  tanplane(1, 2) = -(gpn[1] * gpn[2]);
  tanplane(2, 0) = -(gpn[2] * gpn[0]);
  tanplane(2, 1) = -(gpn[2] * gpn[1]);
  tanplane(2, 2) = 1.0 - (gpn[2] * gpn[2]);

  // interpolation of slave GP jumps (relative displacement increment)
  double sgpjump[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < nrow; ++i)
  {
    sgpjump[0] += sval[i] * (sele.GetNodalCoords(0, i) - sele.GetNodalCoordsOld(0, i));
    sgpjump[1] += sval[i] * (sele.GetNodalCoords(1, i) - sele.GetNodalCoordsOld(1, i));
    sgpjump[2] += sval[i] * (sele.GetNodalCoords(2, i) - sele.GetNodalCoordsOld(2, i));
  }

  // interpolation of master GP jumps (relative displacement increment)
  double mgpjump[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < ncol; ++i)
  {
    mgpjump[0] += mval[i] * (mele.GetNodalCoords(0, i) - mele.GetNodalCoordsOld(0, i));
    mgpjump[1] += mval[i] * (mele.GetNodalCoords(1, i) - mele.GetNodalCoordsOld(1, i));
    mgpjump[2] += mval[i] * (mele.GetNodalCoords(2, i) - mele.GetNodalCoordsOld(2, i));
  }

  // build jump (relative displacement increment) at current GP
  jump(0, 0) = (sgpjump[0] - mgpjump[0]);
  jump(1, 0) = (sgpjump[1] - mgpjump[1]);
  jump(2, 0) = (sgpjump[2] - mgpjump[2]);

  // build lagrange multiplier at current GP
  for (int i = 0; i < nrow; ++i)
  {
    lm(0, 0) += lmval[i] * (*lagmult)(0, i);
    lm(1, 0) += lmval[i] * (*lagmult)(1, i);
    lm(2, 0) += lmval[i] * (*lagmult)(2, i);
  }

  // build tangential jump
  jumptan.Multiply(tanplane, jump);

  // build tangential lm
  lmtan.Multiply(tanplane, lm);

  // evaluate wear
  // not including wearcoefficient
  // normal contact stress
  for (int i = 0; i < 3; ++i) wearval[0] += lm(i, 0) * gpn[i];

  // absolute value of relative tangential jump
  jumpval[0] = sqrt(jumptan(0, 0) * jumptan(0, 0) + jumptan(1, 0) * jumptan(1, 0) +
                    jumptan(2, 0) * jumptan(2, 0));

  // steady state wear
  if (sswear_) jumpval[0] = ssslip_;

  // no jump --> no wear
  if (abs(jumpval[0]) < 1e-12) return;

  // product
  if (wearimpl_)
    wearval[0] = (wearval[0]) * jumpval[0];
  else
    wearval[0] = abs(wearval[0]) * jumpval[0];

  // normal contact stress
  for (int i = 0; i < 3; ++i) gplm[i] = lm(i, 0);

  for (int i = 0; i < 3; ++i) lm_lin += gpn[i] * gplm[i];

  // add to node
  for (int j = 0; j < nrow; ++j)
  {
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(snodes[j]);

    double prod = 0.0;
    if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
      prod = sval[j] * wearval[0] * jac * wgt;
    else
      prod = lmval[j] * wearval[0] * jac * wgt;

    // add current Gauss point's contribution to wseg
    cnode->AddDeltaWeightedWearValue(prod);
  }

  // linearization without lm weighting and jac.
  if (wearimpl_ or WearType() == INPAR::WEAR::wear_primvar)
  {
    int linsize = 0;
    for (int i = 0; i < nrow; ++i)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(snodes[i]);
      linsize += cnode->GetLinsize();
    }

    // evaluate the GP wear function derivatives
    GEN::pairedvector<int, double> ddualgp_x((nrow + ncol) * ndof + linsize);
    GEN::pairedvector<int, double> ddualgp_y((nrow + ncol) * ndof + linsize);
    GEN::pairedvector<int, double> ddualgp_z((nrow + ncol) * ndof + linsize);

    GEN::pairedvector<int, double> ddualgp_x_sxi((nrow + ncol) * ndof + linsize);
    GEN::pairedvector<int, double> ddualgp_y_sxi((nrow + ncol) * ndof + linsize);
    GEN::pairedvector<int, double> ddualgp_z_sxi((nrow + ncol) * ndof + linsize);

    std::vector<std::vector<GEN::pairedvector<int, double>>> tanggp(
        3, std::vector<GEN::pairedvector<int, double>>(3, (ncol * ndof) + linsize));

    // lin. abs(x) = x/abs(x) * lin x.
    double xabsx = (1.0 / abs(jumpval[0])) * lm_lin;
    double absx = (1.0 / abs(jumpval[0]));


    // **********************************************************************
    // (1) Lin of normal for LM -- deriv normal maps from weighted gap lin.
    for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
      dweargp[p->first] += abs(jumpval[0]) * gplm[0] * (p->second);

    for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
      dweargp[p->first] += abs(jumpval[0]) * gplm[1] * (p->second);

    for (_CI p = dnmap_unit[2].begin(); p != dnmap_unit[2].end(); ++p)
      dweargp[p->first] += abs(jumpval[0]) * gplm[2] * (p->second);

    // **********************************************************************
    // (2) Lin. of dual shape function of LM mult.
    for (int i = 0; i < nrow; ++i)
    {
      for (int m = 0; m < nrow; ++m)
      {
        for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
             p != dualmap.end(); ++p)
        {
          ddualgp_x[p->first] += (*lagmult)(0, m) * sval[m] * (p->second)(i, m);
          ddualgp_y[p->first] += (*lagmult)(1, m) * sval[m] * (p->second)(i, m);
          ddualgp_z[p->first] += (*lagmult)(2, m) * sval[m] * (p->second)(i, m);
        }
      }
    }
    for (_CI p = ddualgp_x.begin(); p != ddualgp_x.end(); ++p)
      dweargp[p->first] += abs(jumpval[0]) * gpn[0] * (p->second);
    for (_CI p = ddualgp_y.begin(); p != ddualgp_y.end(); ++p)
      dweargp[p->first] += abs(jumpval[0]) * gpn[1] * (p->second);
    for (_CI p = ddualgp_z.begin(); p != ddualgp_z.end(); ++p)
      dweargp[p->first] += abs(jumpval[0]) * gpn[2] * (p->second);

    // LM deriv
    for (int i = 0; i < nrow; ++i)
    {
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
      {
        ddualgp_x_sxi[p->first] += (*lagmult)(0, i) * lmderiv(i, 0) * (p->second);
        ddualgp_y_sxi[p->first] += (*lagmult)(1, i) * lmderiv(i, 0) * (p->second);
        ddualgp_z_sxi[p->first] += (*lagmult)(2, i) * lmderiv(i, 0) * (p->second);
      }
      for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
      {
        ddualgp_x_sxi[p->first] += (*lagmult)(0, i) * lmderiv(i, 1) * (p->second);
        ddualgp_y_sxi[p->first] += (*lagmult)(1, i) * lmderiv(i, 1) * (p->second);
        ddualgp_z_sxi[p->first] += (*lagmult)(2, i) * lmderiv(i, 1) * (p->second);
      }
    }
    for (_CI p = ddualgp_x_sxi.begin(); p != ddualgp_x_sxi.end(); ++p)
      dweargp[p->first] += abs(jumpval[0]) * gpn[0] * (p->second);
    for (_CI p = ddualgp_y_sxi.begin(); p != ddualgp_y_sxi.end(); ++p)
      dweargp[p->first] += abs(jumpval[0]) * gpn[1] * (p->second);
    for (_CI p = ddualgp_z_sxi.begin(); p != ddualgp_z_sxi.end(); ++p)
      dweargp[p->first] += abs(jumpval[0]) * gpn[2] * (p->second);

    // **********************************************************************
    // (3) absolute incremental slip linearization:
    // (a) build directional derivative of slave GP tagent

    // lin tangplane: 1-n x n --> - ( dn x n + n x dn )
    for (int i = 0; i < 3; ++i)
    {
      for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
        tanggp[0][i][p->first] -= gpn[i] * (p->second);

      for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
        tanggp[1][i][p->first] -= gpn[i] * (p->second);

      for (_CI p = dnmap_unit[2].begin(); p != dnmap_unit[2].end(); ++p)
        tanggp[2][i][p->first] -= gpn[i] * (p->second);
    }
    for (int i = 0; i < 3; ++i)
    {
      for (_CI p = dnmap_unit[0].begin(); p != dnmap_unit[0].end(); ++p)
        tanggp[i][0][p->first] -= gpn[i] * (p->second);

      for (_CI p = dnmap_unit[1].begin(); p != dnmap_unit[1].end(); ++p)
        tanggp[i][1][p->first] -= gpn[i] * (p->second);

      for (_CI p = dnmap_unit[2].begin(); p != dnmap_unit[2].end(); ++p)
        tanggp[i][2][p->first] -= gpn[i] * (p->second);
    }

    GEN::pairedvector<int, double> dt0((ncol * ndof) + linsize);
    GEN::pairedvector<int, double> dt1((ncol * ndof) + linsize);
    GEN::pairedvector<int, double> dt2((ncol * ndof) + linsize);

    // xccord from tang jump lin --> lin tangplane * jump
    for (_CI p = tanggp[0][0].begin(); p != tanggp[0][0].end(); ++p)
      dt0[p->first] += (p->second) * jump(0, 0);

    for (_CI p = tanggp[0][1].begin(); p != tanggp[0][1].end(); ++p)
      dt0[p->first] += (p->second) * jump(1, 0);

    for (_CI p = tanggp[0][2].begin(); p != tanggp[0][2].end(); ++p)
      dt0[p->first] += (p->second) * jump(2, 0);

    // yccord from tang jump lin
    for (_CI p = tanggp[1][0].begin(); p != tanggp[1][0].end(); ++p)
      dt1[p->first] += (p->second) * jump(0, 0);

    for (_CI p = tanggp[1][1].begin(); p != tanggp[1][1].end(); ++p)
      dt1[p->first] += (p->second) * jump(1, 0);

    for (_CI p = tanggp[1][2].begin(); p != tanggp[1][2].end(); ++p)
      dt1[p->first] += (p->second) * jump(2, 0);

    // zccord from tang jump lin
    for (_CI p = tanggp[2][0].begin(); p != tanggp[2][0].end(); ++p)
      dt2[p->first] += (p->second) * jump(0, 0);

    for (_CI p = tanggp[2][1].begin(); p != tanggp[2][1].end(); ++p)
      dt2[p->first] += (p->second) * jump(1, 0);

    for (_CI p = tanggp[2][2].begin(); p != tanggp[2][2].end(); ++p)
      dt2[p->first] += (p->second) * jump(2, 0);


    // add to weargp :  1/abs(tangplane*jump) * LM * tanjump^T * (lin. tangplane * jump)
    for (_CI p = dt0.begin(); p != dt0.end(); ++p)
      dweargp[p->first] += xabsx * (p->second) * jumptan(0, 0);
    for (_CI p = dt1.begin(); p != dt1.end(); ++p)
      dweargp[p->first] += xabsx * (p->second) * jumptan(1, 0);
    for (_CI p = dt2.begin(); p != dt2.end(); ++p)
      dweargp[p->first] += xabsx * (p->second) * jumptan(2, 0);

    // slip lin. for discrete wear
    // u/abs(u) * lin tang * jump
    if (WearType() == INPAR::WEAR::wear_primvar)
    {
      for (_CI p = dt0.begin(); p != dt0.end(); ++p)
        dsliptmatrixgp[p->first] += absx * (p->second) * jumptan(0, 0);

      for (_CI p = dt1.begin(); p != dt1.end(); ++p)
        dsliptmatrixgp[p->first] += absx * (p->second) * jumptan(1, 0);

      for (_CI p = dt2.begin(); p != dt2.end(); ++p)
        dsliptmatrixgp[p->first] += absx * (p->second) * jumptan(2, 0);
    }

    // **********************************************************************
    // (c) build directional derivative of jump
    GEN::pairedvector<int, double> dmap_slcoord_gp_x((nrow + ncol) * ndof + linsize);
    GEN::pairedvector<int, double> dmap_slcoord_gp_y((nrow + ncol) * ndof + linsize);
    GEN::pairedvector<int, double> dmap_slcoord_gp_z((nrow + ncol) * ndof + linsize);

    GEN::pairedvector<int, double> dmap_mcoord_gp_x((nrow + ncol) * ndof + linsize);
    GEN::pairedvector<int, double> dmap_mcoord_gp_y((nrow + ncol) * ndof + linsize);
    GEN::pairedvector<int, double> dmap_mcoord_gp_z((nrow + ncol) * ndof + linsize);

    GEN::pairedvector<int, double> dmap_coord_x((nrow + ncol) * ndof + linsize);
    GEN::pairedvector<int, double> dmap_coord_y((nrow + ncol) * ndof + linsize);
    GEN::pairedvector<int, double> dmap_coord_z((nrow + ncol) * ndof + linsize);

    // lin slave part -- sxi
    for (int i = 0; i < nrow; ++i)
    {
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
      {
        double valx = sderiv(i, 0) * (sele.GetNodalCoords(0, i) - sele.GetNodalCoordsOld(0, i));
        dmap_slcoord_gp_x[p->first] += valx * (p->second);
        double valy = sderiv(i, 0) * (sele.GetNodalCoords(1, i) - sele.GetNodalCoordsOld(1, i));
        dmap_slcoord_gp_y[p->first] += valy * (p->second);
        double valz = sderiv(i, 0) * (sele.GetNodalCoords(2, i) - sele.GetNodalCoordsOld(2, i));
        dmap_slcoord_gp_z[p->first] += valz * (p->second);
      }
      for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
      {
        double valx = sderiv(i, 1) * (sele.GetNodalCoords(0, i) - sele.GetNodalCoordsOld(0, i));
        dmap_slcoord_gp_x[p->first] += valx * (p->second);
        double valy = sderiv(i, 1) * (sele.GetNodalCoords(1, i) - sele.GetNodalCoordsOld(1, i));
        dmap_slcoord_gp_y[p->first] += valy * (p->second);
        double valz = sderiv(i, 1) * (sele.GetNodalCoords(2, i) - sele.GetNodalCoordsOld(2, i));
        dmap_slcoord_gp_z[p->first] += valz * (p->second);
      }
    }

    // lin master part -- mxi
    for (int i = 0; i < ncol; ++i)
    {
      for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
      {
        double valx = mderiv(i, 0) * (mele.GetNodalCoords(0, i) - mele.GetNodalCoordsOld(0, i));
        dmap_mcoord_gp_x[p->first] += valx * (p->second);
        double valy = mderiv(i, 0) * (mele.GetNodalCoords(1, i) - mele.GetNodalCoordsOld(1, i));
        dmap_mcoord_gp_y[p->first] += valy * (p->second);
        double valz = mderiv(i, 0) * (mele.GetNodalCoords(2, i) - mele.GetNodalCoordsOld(2, i));
        dmap_mcoord_gp_z[p->first] += valz * (p->second);
      }
      for (_CI p = dmxigp[1].begin(); p != dmxigp[1].end(); ++p)
      {
        double valx = mderiv(i, 1) * (mele.GetNodalCoords(0, i) - mele.GetNodalCoordsOld(0, i));
        dmap_mcoord_gp_x[p->first] += valx * (p->second);
        double valy = mderiv(i, 1) * (mele.GetNodalCoords(1, i) - mele.GetNodalCoordsOld(1, i));
        dmap_mcoord_gp_y[p->first] += valy * (p->second);
        double valz = mderiv(i, 1) * (mele.GetNodalCoords(2, i) - mele.GetNodalCoordsOld(2, i));
        dmap_mcoord_gp_z[p->first] += valz * (p->second);
      }
    }

    // deriv slave x-coords
    for (int i = 0; i < nrow; ++i)
    {
      CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(snodes[i]);

      dmap_slcoord_gp_x[snode->Dofs()[0]] += sval[i];
      dmap_slcoord_gp_y[snode->Dofs()[1]] += sval[i];
      dmap_slcoord_gp_z[snode->Dofs()[2]] += sval[i];
    }
    // deriv master x-coords
    for (int i = 0; i < ncol; ++i)
    {
      CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mnodes[i]);

      dmap_mcoord_gp_x[mnode->Dofs()[0]] += mval[i];
      dmap_mcoord_gp_y[mnode->Dofs()[1]] += mval[i];
      dmap_mcoord_gp_z[mnode->Dofs()[2]] += mval[i];
    }

    // slave: add to jumplin
    for (_CI p = dmap_slcoord_gp_x.begin(); p != dmap_slcoord_gp_x.end(); ++p)
      dmap_coord_x[p->first] += (p->second);
    for (_CI p = dmap_slcoord_gp_y.begin(); p != dmap_slcoord_gp_y.end(); ++p)
      dmap_coord_y[p->first] += (p->second);
    for (_CI p = dmap_slcoord_gp_z.begin(); p != dmap_slcoord_gp_z.end(); ++p)
      dmap_coord_z[p->first] += (p->second);

    // master: add to jumplin
    for (_CI p = dmap_mcoord_gp_x.begin(); p != dmap_mcoord_gp_x.end(); ++p)
      dmap_coord_x[p->first] -= (p->second);
    for (_CI p = dmap_mcoord_gp_y.begin(); p != dmap_mcoord_gp_y.end(); ++p)
      dmap_coord_y[p->first] -= (p->second);
    for (_CI p = dmap_mcoord_gp_z.begin(); p != dmap_mcoord_gp_z.end(); ++p)
      dmap_coord_z[p->first] -= (p->second);

    // matrix vector prod -- tan
    GEN::pairedvector<int, double> lintan0((nrow + ncol) * ndof + linsize);
    GEN::pairedvector<int, double> lintan1((nrow + ncol) * ndof + linsize);
    GEN::pairedvector<int, double> lintan2((nrow + ncol) * ndof + linsize);

    for (_CI p = dmap_coord_x.begin(); p != dmap_coord_x.end(); ++p)
      lintan0[p->first] += tanplane(0, 0) * (p->second);
    for (_CI p = dmap_coord_y.begin(); p != dmap_coord_y.end(); ++p)
      lintan0[p->first] += tanplane(0, 1) * (p->second);
    for (_CI p = dmap_coord_z.begin(); p != dmap_coord_z.end(); ++p)
      lintan0[p->first] += tanplane(0, 2) * (p->second);

    for (_CI p = dmap_coord_x.begin(); p != dmap_coord_x.end(); ++p)
      lintan1[p->first] += tanplane(1, 0) * (p->second);
    for (_CI p = dmap_coord_y.begin(); p != dmap_coord_y.end(); ++p)
      lintan1[p->first] += tanplane(1, 1) * (p->second);
    for (_CI p = dmap_coord_z.begin(); p != dmap_coord_z.end(); ++p)
      lintan1[p->first] += tanplane(1, 2) * (p->second);

    for (_CI p = dmap_coord_x.begin(); p != dmap_coord_x.end(); ++p)
      lintan2[p->first] += tanplane(2, 0) * (p->second);
    for (_CI p = dmap_coord_y.begin(); p != dmap_coord_y.end(); ++p)
      lintan2[p->first] += tanplane(2, 1) * (p->second);
    for (_CI p = dmap_coord_z.begin(); p != dmap_coord_z.end(); ++p)
      lintan2[p->first] += tanplane(2, 2) * (p->second);

    // add to dweargp
    for (_CI p = lintan0.begin(); p != lintan0.end(); ++p)
      dweargp[p->first] += xabsx * jumptan(0, 0) * (p->second);
    for (_CI p = lintan1.begin(); p != lintan1.end(); ++p)
      dweargp[p->first] += xabsx * jumptan(1, 0) * (p->second);
    for (_CI p = lintan2.begin(); p != lintan2.end(); ++p)
      dweargp[p->first] += xabsx * jumptan(2, 0) * (p->second);

    // slip lin. for discrete wear
    // u/abs(u) * tang * lin jump
    if (WearType() == INPAR::WEAR::wear_primvar)
    {
      for (_CI p = lintan0.begin(); p != lintan0.end(); ++p)
        dsliptmatrixgp[p->first] += absx * (p->second) * jumptan(0, 0);

      for (_CI p = lintan1.begin(); p != lintan1.end(); ++p)
        dsliptmatrixgp[p->first] += absx * (p->second) * jumptan(1, 0);

      for (_CI p = lintan2.begin(); p != lintan2.end(); ++p)
        dsliptmatrixgp[p->first] += absx * (p->second) * jumptan(2, 0);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Compute entries for E and T matrix at GP (Slave)         farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_TE(MORTAR::MortarElement& sele,
    LINALG::SerialDenseVector& lmval, LINALG::SerialDenseVector& sval, double& jac, double& wgt,
    double* jumpval)
{
  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: Null pointer!");

  int nrow = sele.NumNode();

  if (WearShapeFcn() == INPAR::WEAR::wear_shape_standard)
  {
    for (int k = 0; k < nrow; ++k)
    {
      CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(snodes[k]);

      for (int j = 0; j < nrow; ++j)
      {
        CONTACT::FriNode* snode = dynamic_cast<CONTACT::FriNode*>(snodes[j]);

        // multiply the two shape functions
        double prod1 = sval[k] * lmval[j] * abs(*jumpval) * jac * wgt;
        double prod2 = sval[k] * sval[j] * jac * wgt;

        int col = snode->Dofs()[0];
        int row = 0;

        if (abs(prod1) > MORTARINTTOL) cnode->AddTValue(row, col, prod1);
        if (abs(prod2) > MORTARINTTOL) cnode->AddEValue(row, col, prod2);
      }
    }
  }
  else if (WearShapeFcn() == INPAR::WEAR::wear_shape_dual)
  {
    for (int k = 0; k < nrow; ++k)
    {
      CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(snodes[k]);

      for (int j = 0; j < nrow; ++j)
      {
        CONTACT::FriNode* snode = dynamic_cast<CONTACT::FriNode*>(snodes[j]);

        // multiply the two shape functions
        double prod1 = lmval[k] * lmval[j] * abs(*jumpval) * jac * wgt;
        double prod2 = lmval[k] * sval[j] * jac * wgt;

        int col = snode->Dofs()[0];
        int row = 0;

        if (abs(prod1) > MORTARINTTOL) cnode->AddTValue(row, col, prod1);

        // diagonal E matrix
        if (j == k)
          if (abs(prod2) > MORTARINTTOL) cnode->AddEValue(row, col, prod2);
      }
    }
  }
  else
    dserror("Chosen wear shape function not supported!");


  return;
}

/*----------------------------------------------------------------------*
 |  Compute entries for E and T matrix at GP (Master)        farah 11/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_TE_Master(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseVector& lm2val, LINALG::SerialDenseVector& mval, double& jac, double& wgt,
    double* jumpval, const Epetra_Comm& comm)
{
  if (sele.Owner() != comm.MyPID()) return;

  // mele is involved for both-sided wear
  mele.SetAttached() = true;

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // get master element nodes themselves
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  int nrow = mele.NumNode();
  int ncol = sele.NumNode();

  if (WearShapeFcn() == INPAR::WEAR::wear_shape_standard)
  {
    for (int k = 0; k < nrow; ++k)
    {
      CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(mnodes[k]);
      int row = 0;

      for (int j = 0; j < nrow; ++j)
      {
        CONTACT::FriNode* snode = dynamic_cast<CONTACT::FriNode*>(mnodes[j]);

        // multiply the two shape functions
        double prod2 = mval[k] * mval[j] * jac * wgt;

        int col = snode->Dofs()[0];
        if (abs(prod2) > MORTARINTTOL) cnode->AddEValue(row, col, prod2);
      }
      for (int j = 0; j < ncol; ++j)
      {
        CONTACT::FriNode* snode = dynamic_cast<CONTACT::FriNode*>(snodes[j]);

        // multiply the two shape functions
        double prod1 = mval[k] * lmval[j] * abs(*jumpval) * jac * wgt;

        int col = snode->Dofs()[0];
        if (abs(prod1) > MORTARINTTOL) cnode->AddTValue(row, col, prod1);
      }
    }
  }
  else if (WearShapeFcn() == INPAR::WEAR::wear_shape_dual)
  {
    for (int k = 0; k < nrow; ++k)
    {
      CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(mnodes[k]);
      int row = 0;

      for (int j = 0; j < nrow; ++j)
      {
        CONTACT::FriNode* snode = dynamic_cast<CONTACT::FriNode*>(mnodes[j]);

        // multiply the two shape functions
        double prod2 = lm2val[k] * mval[j] * jac * wgt;

        int col = snode->Dofs()[0];
        if (abs(prod2) > MORTARINTTOL and j == k) cnode->AddEValue(row, col, prod2);
      }
      for (int j = 0; j < ncol; ++j)
      {
        CONTACT::FriNode* snode = dynamic_cast<CONTACT::FriNode*>(snodes[j]);

        // multiply the two shape functions
        double prod1 = lm2val[k] * lmval[j] * abs(*jumpval) * jac * wgt;

        int col = snode->Dofs()[0];
        if (abs(prod1) > MORTARINTTOL) cnode->AddTValue(row, col, prod1);
      }
    }
  }
  else
    dserror("Chosen wear shape function not supported!");

  return;
}

/*----------------------------------------------------------------------*
 |  Compute lin for T and E matrix -- discr. wear            farah 11/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_2D_TE_Master_Lin(int& iter,  // like k
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval,
    LINALG::SerialDenseVector& mval, LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseMatrix& mderiv, LINALG::SerialDenseMatrix& lmderiv, double& dsxideta,
    double& dxdsxi, double& dxdsxidsxi, double& wgt, double* jumpval,
    const GEN::pairedvector<int, double>& dsxigp, const GEN::pairedvector<int, double>& dmxigp,
    const GEN::pairedvector<int, double>& derivjac,
    const GEN::pairedvector<int, double>& dsliptmatrixgp,
    const std::vector<GEN::pairedvector<int, double>>& ximaps,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap, const Epetra_Comm& comm)
{
  if (sele.Owner() != comm.MyPID()) return;

  const int nrow = sele.NumNode();
  const int ncol = mele.NumNode();

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // get master element nodes themselves
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mnodes[iter]);
  if (!mymrtrnode) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  if (WearShapeFcn() == INPAR::WEAR::wear_shape_standard)
  {
    // integrate LinT
    for (int j = 0; j < nrow; ++j)
    {
      // global master node ID
      int mgid = sele.Nodes()[j]->Id();
      double fac = 0.0;

      // get the correct map as a reference
      std::map<int, double>& tmmap_jk =
          dynamic_cast<CONTACT::FriNode*>(mymrtrnode)->WearData().GetDerivTw()[mgid];

      // (1) Lin(Phi) - dual shape functions
      for (int m = 0; m < nrow; ++m)
      {
        fac = wgt * mval[m] * sval[j] * dsxideta * dxdsxi * abs(jumpval[0]);
        for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
             p != dualmap.end(); ++p)
          tmmap_jk[p->first] += fac * (p->second)(iter, m);
      }

      // (2) Lin(Phi) - slave GP coordinates
      fac = wgt * lmderiv(j, 0) * mval[iter] * dsxideta * dxdsxi * abs(jumpval[0]);
      for (_CI p = dsxigp.begin(); p != dsxigp.end(); ++p) tmmap_jk[p->first] += fac * (p->second);

      // (3) Lin(Phi) - slave GP coordinates
      fac = wgt * lmval[j] * mderiv(iter, 0) * dsxideta * dxdsxi * abs(jumpval[0]);
      for (_CI p = dmxigp.begin(); p != dmxigp.end(); ++p) tmmap_jk[p->first] += fac * (p->second);

      // (4) Lin(dsxideta) - segment end coordinates
      fac = wgt * lmval[j] * mval[iter] * dxdsxi * abs(jumpval[0]);
      for (_CI p = ximaps[0].begin(); p != ximaps[0].end(); ++p)
        tmmap_jk[p->first] -= 0.5 * fac * (p->second);
      for (_CI p = ximaps[1].begin(); p != ximaps[1].end(); ++p)
        tmmap_jk[p->first] += 0.5 * fac * (p->second);

      // (5) Lin(dxdsxi) - slave GP Jacobian
      fac = wgt * lmval[j] * mval[iter] * dsxideta * abs(jumpval[0]);
      for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
        tmmap_jk[p->first] += fac * (p->second);

      // (6) Lin(dxdsxi) - slave GP coordinates
      fac = wgt * lmval[j] * mval[iter] * dsxideta * dxdsxidsxi * abs(jumpval[0]);
      for (_CI p = dsxigp.begin(); p != dsxigp.end(); ++p) tmmap_jk[p->first] += fac * (p->second);

      // (7) Lin(wear)
      fac = wgt * lmval[j] * mval[iter] * dsxideta * dxdsxi;
      for (_CI p = dsliptmatrixgp.begin(); p != dsliptmatrixgp.end(); ++p)  // dsliptmatrixgp
      {
        tmmap_jk[p->first] += fac * (p->second);
      }
    }  // end integrate linT

    //********************************************************************************
    // integrate LinE
    for (int j = 0; j < ncol; ++j)
    {
      // global master node ID
      int mgid = mele.Nodes()[j]->Id();
      double fac = 0.0;

      // get the correct map as a reference
      std::map<int, double>& emmap_jk =
          dynamic_cast<CONTACT::FriNode*>(mymrtrnode)->WearData().GetDerivE()[mgid];

      // (1) Lin(Phi) - slave GP coordinates
      fac = wgt * mderiv(iter, 0) * mval[j] * dsxideta * dxdsxi;
      for (_CI p = dmxigp.begin(); p != dmxigp.end(); ++p) emmap_jk[p->first] += fac * (p->second);

      // (2) Lin(Phi) - slave GP coordinates
      fac = wgt * mval[iter] * mderiv(j, 0) * dsxideta * dxdsxi;
      for (_CI p = dmxigp.begin(); p != dmxigp.end(); ++p) emmap_jk[p->first] += fac * (p->second);

      // (3) Lin(dsxideta) - segment end coordinates
      fac = wgt * mval[iter] * mval[j] * dxdsxi;
      for (_CI p = ximaps[0].begin(); p != ximaps[0].end(); ++p)
        emmap_jk[p->first] -= 0.5 * fac * (p->second);
      for (_CI p = ximaps[1].begin(); p != ximaps[1].end(); ++p)
        emmap_jk[p->first] += 0.5 * fac * (p->second);

      // (4) Lin(dxdsxi) - slave GP Jacobian
      fac = wgt * mval[iter] * mval[j] * dsxideta;
      for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
        emmap_jk[p->first] += fac * (p->second);

      // (5) Lin(dxdsxi) - slave GP coordinates
      fac = wgt * mval[iter] * mval[j] * dsxideta * dxdsxidsxi;
      for (_CI p = dsxigp.begin(); p != dsxigp.end(); ++p) emmap_jk[p->first] += fac * (p->second);
    }  // end integrate linE
  }
  else if (WearShapeFcn() ==
           INPAR::WEAR::wear_shape_dual)  //******************************************
  {
    dserror("Chosen shapefunctions 'wear_shape_dual' not supported!");
  }
  else
    dserror("Chosen shapefunctions not supported!");

  return;
}

/*----------------------------------------------------------------------*
 |  Compute lin for T and E matrix -- discr. wear            farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_2D_TE_Lin(int& iter, MORTAR::MortarElement& sele,
    LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& lmderiv, double& jac, double& wgt,
    double* jumpval, const std::vector<GEN::pairedvector<int, double>>& dsxigp,
    const GEN::pairedvector<int, double>& derivjac,
    const GEN::pairedvector<int, double>& dsliptmatrixgp,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap)
{
  const int nrow = sele.NumNode();

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(snodes[iter]);
  if (!mymrtrnode) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  if (WearShapeFcn() == INPAR::WEAR::wear_shape_standard)
  {
    // integrate LinT
    for (int j = 0; j < nrow; ++j)
    {
      // global master node ID
      int mgid = sele.Nodes()[j]->Id();
      double fac = 0.0;

      // get the correct map as a reference
      std::map<int, double>& tmmap_jk =
          dynamic_cast<CONTACT::FriNode*>(mymrtrnode)->WearData().GetDerivTw()[mgid];

      // (1) Lin(Phi) - dual shape functions
      for (int m = 0; m < nrow; ++m)
      {
        fac = wgt * sval[m] * sval[j] * jac * abs(jumpval[0]);
        for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
             p != dualmap.end(); ++p)
          tmmap_jk[p->first] += fac * (p->second)(iter, m);
      }

      // (2) Lin(Phi) - slave GP coordinates
      fac = wgt * lmderiv(j, 0) * sval[iter] * jac * abs(jumpval[0]);
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
        tmmap_jk[p->first] += fac * (p->second);

      // (3) Lin(Phi) - slave GP coordinates
      fac = wgt * lmval[j] * sderiv(iter, 0) * jac * abs(jumpval[0]);
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
        tmmap_jk[p->first] += fac * (p->second);

      // (4) Lin(dxdsxi) - slave GP Jacobian
      fac = wgt * lmval[j] * sval[iter] * abs(jumpval[0]);
      for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
        tmmap_jk[p->first] += fac * (p->second);

      // (7) Lin(wear)
      if (!sswear_)
      {
        fac = wgt * lmval[j] * sval[iter] * jac;
        for (_CI p = dsliptmatrixgp.begin(); p != dsliptmatrixgp.end(); ++p)  // dsliptmatrixgp
        {
          tmmap_jk[p->first] += fac * (p->second);
        }
      }

    }  // end integrate linT

    //********************************************************************************
    // integrate LinE
    for (int j = 0; j < nrow; ++j)
    {
      // global master node ID
      int mgid = sele.Nodes()[j]->Id();
      double fac = 0.0;

      // get the correct map as a reference
      std::map<int, double>& emmap_jk =
          dynamic_cast<CONTACT::FriNode*>(mymrtrnode)->WearData().GetDerivE()[mgid];

      // (1) Lin(Phi) - slave GP coordinates
      fac = wgt * sderiv(iter, 0) * sval[j] * jac;
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
        emmap_jk[p->first] += fac * (p->second);

      // (2) Lin(Phi) - slave GP coordinates
      fac = wgt * sval[iter] * sderiv(j, 0) * jac;
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
        emmap_jk[p->first] += fac * (p->second);

      // (4) Lin(dxdsxi) - slave GP Jacobian
      fac = wgt * sval[iter] * sval[j];
      for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
        emmap_jk[p->first] += fac * (p->second);
    }  // end integrate linE
  }
  else if (WearShapeFcn() ==
           INPAR::WEAR::wear_shape_dual)  //******************************************
  {
    // integrate LinT
    for (int j = 0; j < nrow; ++j)
    {
      // global master node ID
      int mgid = sele.Nodes()[j]->Id();
      double fac = 0.0;

      // get the correct map as a reference
      std::map<int, double>& tmmap_jk =
          dynamic_cast<CONTACT::FriNode*>(mymrtrnode)->WearData().GetDerivTw()[mgid];

      // (1) Lin(Phi) - dual shape functions
      for (int m = 0; m < nrow; ++m)
      {
        fac = wgt * sval[m] * lmval[j] * jac * abs(jumpval[0]);
        for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
             p != dualmap.end(); ++p)
          tmmap_jk[p->first] += fac * (p->second)(iter, m);
      }

      // (2) Lin(Phi) - slave GP coordinates
      fac = wgt * lmderiv(j, 0) * lmval[iter] * jac * abs(jumpval[0]);
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
        tmmap_jk[p->first] += fac * (p->second);


      // (1) Lin(Phi) - dual shape functions
      for (int m = 0; m < nrow; ++m)
      {
        fac = wgt * sval[m] * sval[j] * jac * abs(jumpval[0]);
        for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
             p != dualmap.end(); ++p)
          tmmap_jk[p->first] += fac * (p->second)(iter, m);
      }

      // (2) Lin(Phi) - slave GP coordinates
      fac = wgt * lmderiv(iter, 0) * lmval[j] * jac * abs(jumpval[0]);
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
        tmmap_jk[p->first] += fac * (p->second);

      // (3) Lin(dxdsxi) - slave GP Jacobian
      fac = wgt * lmval[j] * lmval[iter] * abs(jumpval[0]);
      for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
        tmmap_jk[p->first] += fac * (p->second);

      // (4) Lin(wear)
      fac = wgt * lmval[j] * lmval[iter] * jac;
      for (_CI p = dsliptmatrixgp.begin(); p != dsliptmatrixgp.end(); ++p)  // dsliptmatrixgp
        tmmap_jk[p->first] += fac * (p->second);
    }  // end integrate linT

    //********************************************************************************
    // integrate LinE
    for (int j = 0; j < nrow; ++j)
    {
      // global master node ID
      int mgid = sele.Nodes()[j]->Id();
      double fac = 0.0;

      // get the correct map as a reference
      std::map<int, double>& emmap_jk =
          dynamic_cast<CONTACT::FriNode*>(mymrtrnode)->WearData().GetDerivE()[mgid];

      // (1) Lin(Phi) - dual shape functions
      for (int m = 0; m < nrow; ++m)
      {
        fac = wgt * sval[m] * sval[j] * jac;
        for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
             p != dualmap.end(); ++p)
          emmap_jk[p->first] += fac * (p->second)(iter, m);
      }

      // (1) Lin(Phi) - slave GP coordinates
      fac = wgt * lmderiv(iter, 0) * sval[j] * jac;
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
        emmap_jk[p->first] += fac * (p->second);

      // (2) Lin(Phi) - slave GP coordinates
      fac = wgt * lmval[iter] * sderiv(j, 0) * jac;
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
        emmap_jk[p->first] += fac * (p->second);

      // (3) Lin(dxdsxi) - slave GP Jacobian
      fac = wgt * lmval[iter] * sval[j];
      for (_CI p = derivjac.begin(); p != derivjac.end(); ++p)
        emmap_jk[p->first] += fac * (p->second);
    }  // end integrate linE
  }
  else
    dserror("Chosen shape functions not supported!");

  return;
}

/*----------------------------------------------------------------------*
 |  Compute lin for T and E matrix -- discr. wear            farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_3D_TE_Lin(int& iter, MORTAR::MortarElement& sele,
    LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& lmderiv, double& jac, double& wgt,
    double* jumpval, const std::vector<GEN::pairedvector<int, double>>& dsxigp,
    const GEN::pairedvector<int, double>& jacintcellmap,
    const GEN::pairedvector<int, double>& dsliptmatrixgp,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap)
{
  const int nrow = sele.NumNode();

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(snodes[iter]);
  if (!mymrtrnode) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  if (WearShapeFcn() == INPAR::WEAR::wear_shape_standard)
  {
    // integrate LinT
    for (int j = 0; j < nrow; ++j)
    {
      // global master node ID
      int mgid = sele.Nodes()[j]->Id();
      double fac = 0.0;

      // get the correct map as a reference
      std::map<int, double>& dtmap_jk =
          dynamic_cast<CONTACT::FriNode*>(mymrtrnode)->WearData().GetDerivTw()[mgid];

      // (1) Lin(Phi) - dual shape functions
      if (dualmap.size() > 0)
        for (int m = 0; m < nrow; ++m)
        {
          fac = wgt * sval[m] * sval[iter] * jac * abs(jumpval[0]);
          for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
               p != dualmap.end(); ++p)
            dtmap_jk[p->first] += fac * (p->second)(j, m);
        }

      // (2) Lin(Phi) - slave GP coordinates
      fac = wgt * lmderiv(j, 0) * sval[iter] * jac * abs(jumpval[0]);
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
        dtmap_jk[p->first] += fac * (p->second);

      fac = wgt * lmderiv(j, 1) * sval[iter] * jac * abs(jumpval[0]);
      for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
        dtmap_jk[p->first] += fac * (p->second);


      // (3) Lin(NMaster) - master GP coordinates
      fac = wgt * lmval[j] * sderiv(iter, 0) * jac * abs(jumpval[0]);
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
        dtmap_jk[p->first] += fac * (p->second);

      fac = wgt * lmval[j] * sderiv(iter, 1) * jac * abs(jumpval[0]);
      for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
        dtmap_jk[p->first] += fac * (p->second);


      // (4) Lin(dsxideta) - intcell GP Jacobian
      fac = wgt * lmval[j] * sval[iter] * abs(jumpval[0]);
      for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
        dtmap_jk[p->first] += fac * (p->second);

      if (!sswear_)
      {
        fac = wgt * lmval[j] * sval[iter] * jac;
        for (_CI p = dsliptmatrixgp.begin(); p != dsliptmatrixgp.end(); ++p)
          dtmap_jk[p->first] += fac * (p->second);
      }
    }

    //********************************************************************************
    // integrate LinE
    for (int j = 0; j < nrow; ++j)
    {
      // global master node ID
      int mgid = sele.Nodes()[j]->Id();
      double fac = 0.0;

      // get the correct map as a reference
      std::map<int, double>& emmap_jk =
          dynamic_cast<CONTACT::FriNode*>(mymrtrnode)->WearData().GetDerivE()[mgid];

      // (2) Lin(Phi) - slave GP coordinates
      fac = wgt * sderiv(j, 0) * sval[iter] * jac;
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
        emmap_jk[p->first] += fac * (p->second);

      fac = wgt * sderiv(j, 1) * sval[iter] * jac;
      for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
        emmap_jk[p->first] += fac * (p->second);

      // (3) Lin(NMaster) - master GP coordinates
      fac = wgt * sval[j] * sderiv(iter, 0) * jac;
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
        emmap_jk[p->first] += fac * (p->second);

      fac = wgt * sval[j] * sderiv(iter, 1) * jac;
      for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
        emmap_jk[p->first] += fac * (p->second);

      // (4) Lin(dsxideta) - intcell GP Jacobian
      fac = wgt * sval[j] * sval[iter];
      for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
        emmap_jk[p->first] += fac * (p->second);

    }  // end integrate linE
  }
  else if (WearShapeFcn() == INPAR::WEAR::wear_shape_dual)
  {
    // integrate LinT
    for (int j = 0; j < nrow; ++j)
    {
      // global master node ID
      int mgid = sele.Nodes()[j]->Id();
      double fac = 0.0;

      // get the correct map as a reference
      std::map<int, double>& dtmap_jk =
          dynamic_cast<CONTACT::FriNode*>(mymrtrnode)->WearData().GetDerivTw()[mgid];

      //**********************************************
      // LM-shape function lin...
      //**********************************************
      // (1) Lin(Phi) - dual shape functions
      if (dualmap.size() > 0)
        for (int m = 0; m < nrow; ++m)
        {
          fac = wgt * sval[m] * lmval[j] * jac * abs(jumpval[0]);
          for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
               p != dualmap.end(); ++p)
            dtmap_jk[p->first] += fac * (p->second)(iter, m);
        }

      // (2) Lin(Phi) - slave GP coordinates
      fac = wgt * lmderiv(j, 0) * lmval[iter] * jac * abs(jumpval[0]);
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
        dtmap_jk[p->first] += fac * (p->second);

      fac = wgt * lmderiv(j, 1) * lmval[iter] * jac * abs(jumpval[0]);
      for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
        dtmap_jk[p->first] += fac * (p->second);

      //**********************************************
      // wear weighting lin...
      //**********************************************
      // (1) Lin(Phi) - dual shape functions
      if (dualmap.size() > 0)
        for (int m = 0; m < nrow; ++m)
        {
          fac = wgt * lmval[iter] * sval[m] * jac * abs(jumpval[0]);
          for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
               p != dualmap.end(); ++p)
            dtmap_jk[p->first] += fac * (p->second)(j, m);
        }

      // (2) Lin(Phi) - slave GP coordinates
      fac = wgt * lmval[j] * lmderiv(iter, 0) * jac * abs(jumpval[0]);
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
        dtmap_jk[p->first] += fac * (p->second);

      fac = wgt * lmval[j] * lmderiv(iter, 1) * jac * abs(jumpval[0]);
      for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
        dtmap_jk[p->first] += fac * (p->second);

      //**********************************************
      // rest
      //**********************************************
      // (4) Lin(dsxideta) - intcell GP Jacobian
      fac = wgt * lmval[j] * lmval[iter] * abs(jumpval[0]);
      for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
        dtmap_jk[p->first] += fac * (p->second);

      fac = wgt * lmval[j] * lmval[iter] * jac;
      for (_CI p = dsliptmatrixgp.begin(); p != dsliptmatrixgp.end(); ++p)
        dtmap_jk[p->first] += fac * (p->second);
    }

    //********************************************************************************
    // integrate LinE
    for (int j = 0; j < nrow; ++j)
    {
      // global master node ID
      int mgid = sele.Nodes()[j]->Id();
      double fac = 0.0;

      // get the correct map as a reference
      std::map<int, double>& emmap_jk =
          dynamic_cast<CONTACT::FriNode*>(mymrtrnode)->WearData().GetDerivE()[mgid];

      //**********************************************
      // wear weighting lin...
      //**********************************************
      // (1) Lin(Phi) - dual shape functions
      if (dualmap.size() > 0)
        for (int m = 0; m < nrow; ++m)
        {
          fac = wgt * sval[m] * sval[iter] * jac;
          for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
               p != dualmap.end(); ++p)
            emmap_jk[p->first] += fac * (p->second)(j, m);
        }

      // (2) Lin(Phi) - slave GP coordinates
      fac = wgt * sderiv(j, 0) * lmval[iter] * jac;
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
        emmap_jk[p->first] += fac * (p->second);

      fac = wgt * sderiv(j, 1) * lmval[iter] * jac;
      for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
        emmap_jk[p->first] += fac * (p->second);

      // (3) Lin(NMaster) - master GP coordinates
      fac = wgt * sval[j] * lmderiv(iter, 0) * jac;
      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
        emmap_jk[p->first] += fac * (p->second);

      fac = wgt * sval[j] * lmderiv(iter, 1) * jac;
      for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
        emmap_jk[p->first] += fac * (p->second);

      // (4) Lin(dsxideta) - intcell GP Jacobian
      fac = wgt * sval[j] * lmval[iter];
      for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
        emmap_jk[p->first] += fac * (p->second);

    }  // end integrate linE
  }
  else
    dserror("Choosen shapefunctions not supported!");

  return;
}

/*----------------------------------------------------------------------*
 |  Compute lin for T and E matrix -- discr. wear (master)   farah 11/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_3D_TE_Master_Lin(int& iter, MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
    LINALG::SerialDenseVector& lmval, LINALG::SerialDenseVector& lm2val,
    LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& mderiv,
    LINALG::SerialDenseMatrix& lmderiv, LINALG::SerialDenseMatrix& lm2deriv, double& jac,
    double& wgt, double* jumpval, const std::vector<GEN::pairedvector<int, double>>& dsxigp,
    const std::vector<GEN::pairedvector<int, double>>& dmxigp,
    const GEN::pairedvector<int, double>& jacintcellmap,
    const GEN::pairedvector<int, double>& dsliptmatrixgp,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dual2map, const Epetra_Comm& comm)
{
  if (sele.Owner() != comm.MyPID()) return;

  const int ncol = mele.NumNode();
  const int nrow = sele.NumNode();

  // get slave element nodes themselves
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mnodes[iter]);
  if (!mymrtrnode) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  if (WearShapeFcn() == INPAR::WEAR::wear_shape_standard)
  {
    // integrate LinT
    for (int j = 0; j < nrow; ++j)
    {
      // global master node ID
      int mgid = sele.Nodes()[j]->Id();
      double fac = 0.0;

      // get the correct map as a reference
      std::map<int, double>& dtmap_jk =
          dynamic_cast<CONTACT::FriNode*>(mymrtrnode)->WearData().GetDerivTw()[mgid];

      // (1) Lin(Phi) - dual shape functions
      if (dualmap.size() > 0)
        for (int m = 0; m < nrow; ++m)
        {
          fac = wgt * sval[m] * mval[iter] * jac * abs(jumpval[0]);
          for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
               p != dualmap.end(); ++p)
            dtmap_jk[p->first] += fac * (p->second)(j, m);
        }

      // (2) Lin(Phi) - slave GP coordinates
      for (int d = 0; d < Dim() - 1; ++d)
      {
        fac = wgt * lmderiv(j, d) * mval[iter] * jac * abs(jumpval[0]);
        for (_CI p = dsxigp[d].begin(); p != dsxigp[d].end(); ++p)
          dtmap_jk[p->first] += fac * (p->second);
      }

      // (3) Lin(NMaster) - master GP coordinates
      for (int d = 0; d < Dim() - 1; ++d)
      {
        fac = wgt * lmval[j] * mderiv(iter, d) * jac * abs(jumpval[0]);
        for (_CI p = dmxigp[d].begin(); p != dmxigp[d].end(); ++p)
          dtmap_jk[p->first] += fac * (p->second);
      }

      // (4) Lin(dsxideta) - intcell GP Jacobian
      fac = wgt * lmval[j] * mval[iter] * abs(jumpval[0]);
      for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
        dtmap_jk[p->first] += fac * (p->second);

      fac = wgt * lmval[j] * mval[iter] * jac;
      for (_CI p = dsliptmatrixgp.begin(); p != dsliptmatrixgp.end(); ++p)
        dtmap_jk[p->first] += fac * (p->second);
    }

    //********************************************************************************
    // integrate LinE
    for (int j = 0; j < ncol; ++j)
    {
      // global master node ID
      int mgid = mele.Nodes()[j]->Id();
      double fac = 0.0;

      // get the correct map as a reference
      std::map<int, double>& emmap_jk =
          dynamic_cast<CONTACT::FriNode*>(mymrtrnode)->WearData().GetDerivE()[mgid];

      // (2) Lin(Phi) - slave GP coordinates
      for (int d = 0; d < Dim() - 1; ++d)
      {
        fac = wgt * mderiv(j, d) * mval[iter] * jac;
        for (_CI p = dmxigp[d].begin(); p != dmxigp[d].end(); ++p)
          emmap_jk[p->first] += fac * (p->second);
      }

      // (3) Lin(NMaster) - master GP coordinates
      for (int d = 0; d < Dim() - 1; ++d)
      {
        fac = wgt * mval[j] * mderiv(iter, d) * jac;
        for (_CI p = dmxigp[d].begin(); p != dmxigp[d].end(); ++p)
          emmap_jk[p->first] += fac * (p->second);
      }

      // (4) Lin(dsxideta) - intcell GP Jacobian
      fac = wgt * mval[j] * mval[iter];
      for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
        emmap_jk[p->first] += fac * (p->second);

    }  // end integrate linE
  }
  //------------------------------------------------------------
  //------------------------------------------------------------
  //------------------------------------------------------------
  else if (WearShapeFcn() == INPAR::WEAR::wear_shape_dual)
  {
    // integrate LinT
    for (int j = 0; j < nrow; ++j)
    {
      // global master node ID
      int mgid = sele.Nodes()[j]->Id();
      double fac = 0.0;

      // get the correct map as a reference
      std::map<int, double>& dtmap_jk =
          dynamic_cast<CONTACT::FriNode*>(mymrtrnode)->WearData().GetDerivTw()[mgid];

      // (1) Lin(Phi) - dual shape functions
      if (dualmap.size() > 0)
        for (int m = 0; m < nrow; ++m)
        {
          fac = wgt * sval[m] * lm2val[iter] * jac * abs(jumpval[0]);
          for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
               p != dualmap.end(); ++p)
            dtmap_jk[p->first] += fac * (p->second)(j, m);
        }

      // (2) Lin(Phi) - slave GP coordinates
      for (int d = 0; d < Dim() - 1; ++d)
      {
        fac = wgt * lmderiv(j, d) * lm2val[iter] * jac * abs(jumpval[0]);
        for (_CI p = dsxigp[d].begin(); p != dsxigp[d].end(); ++p)
          dtmap_jk[p->first] += fac * (p->second);
      }

      //----------------------

      // (1) Lin(Phi) - dual shape functions
      if (dual2map.size() > 0)
        for (int m = 0; m < ncol; ++m)
        {
          fac = wgt * lmval[m] * mval[iter] * jac * abs(jumpval[0]);
          for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p =
                   dual2map.begin();
               p != dual2map.end(); ++p)
            dtmap_jk[p->first] += fac * (p->second)(j, m);
        }

      // (3) Lin(NMaster) - master GP coordinates
      for (int d = 0; d < Dim() - 1; ++d)
      {
        fac = wgt * lmval[j] * lm2deriv(iter, d) * jac * abs(jumpval[0]);
        for (_CI p = dmxigp[d].begin(); p != dmxigp[d].end(); ++p)
          dtmap_jk[p->first] += fac * (p->second);
      }

      //----------------------

      // (4) Lin(dsxideta) - intcell GP Jacobian
      fac = wgt * lmval[j] * lm2val[iter] * abs(jumpval[0]);
      for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
        dtmap_jk[p->first] += fac * (p->second);

      fac = wgt * lmval[j] * lm2val[iter] * jac;
      for (_CI p = dsliptmatrixgp.begin(); p != dsliptmatrixgp.end(); ++p)
        dtmap_jk[p->first] += fac * (p->second);
    }

    //********************************************************************************
    // integrate LinE
    for (int j = 0; j < ncol; ++j)
    {
      // global master node ID
      int mgid = mele.Nodes()[j]->Id();
      double fac = 0.0;

      // get the correct map as a reference
      std::map<int, double>& emmap_jk =
          dynamic_cast<CONTACT::FriNode*>(mymrtrnode)->WearData().GetDerivE()[mgid];

      // (2) Lin(Phi) - slave GP coordinates
      for (int d = 0; d < Dim() - 1; ++d)
      {
        fac = wgt * mderiv(j, d) * lm2val[iter] * jac;
        for (_CI p = dmxigp[d].begin(); p != dmxigp[d].end(); ++p)
          emmap_jk[p->first] += fac * (p->second);
      }

      //----------------------

      // (1) Lin(Phi) - dual shape functions
      if (dual2map.size() > 0)
        for (int m = 0; m < ncol; ++m)
        {
          fac = wgt * mval[m] * mval[iter] * jac;
          for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p =
                   dual2map.begin();
               p != dual2map.end(); ++p)
            emmap_jk[p->first] += fac * (p->second)(j, m);
        }

      // (3) Lin(NMaster) - master GP coordinates
      for (int d = 0; d < Dim() - 1; ++d)
      {
        fac = wgt * mval[j] * lm2deriv(iter, d) * jac;
        for (_CI p = dmxigp[d].begin(); p != dmxigp[d].end(); ++p)
          emmap_jk[p->first] += fac * (p->second);
      }

      //----------------------

      // (4) Lin(dsxideta) - intcell GP Jacobian
      fac = wgt * mval[j] * mval[iter];
      for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
        emmap_jk[p->first] += fac * (p->second);

    }  // end integrate linE
  }
  else
    dserror("Chosen shape functions not supported!");

  return;
}

/*----------------------------------------------------------------------*
 |  Compute entries for weighted slip increment at GP        farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_2D_SlipIncr(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
    LINALG::SerialDenseVector& lmval, LINALG::SerialDenseMatrix& sderiv,
    LINALG::SerialDenseMatrix& mderiv, double& jac, double& wgt, double* jumpvalv,
    const std::vector<GEN::pairedvector<int, double>>& dsxigp,
    const std::vector<GEN::pairedvector<int, double>>& dmxigp,
    GEN::pairedvector<int, double>& dslipgp, int& linsize)
{
  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();
  if (!snodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");
  if (!mnodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // number of nodes (slave, master)
  const int nrow = sele.NumNode();
  const int ncol = mele.NumNode();
  const int ndof = Dim();

  // LIN OF TANGENT
  GEN::pairedvector<int, double> dmap_txsl_gp(ncol * ndof + linsize);
  GEN::pairedvector<int, double> dmap_tysl_gp(ncol * ndof + linsize);

  // build interpolation of slave GP normal and coordinates
  double sjumpv[2] = {0.0, 0.0};
  double mjumpv[2] = {0.0, 0.0};
  double jumpv[2] = {0.0, 0.0};
  double tanv[2] = {0.0, 0.0};

  double tanlength = 0.0;
  for (int i = 0; i < nrow; ++i)
  {
    CONTACT::CoNode* myconode = dynamic_cast<CONTACT::CoNode*>(snodes[i]);

    // nodal tangent interpolation
    tanv[0] += sval[i] * myconode->CoData().txi()[0];
    tanv[1] += sval[i] * myconode->CoData().txi()[1];

    // delta D
    sjumpv[0] += sval[i] * (sele.GetNodalCoords(0, i) - sele.GetNodalCoordsOld(0, i));
    sjumpv[1] += sval[i] * (sele.GetNodalCoords(1, i) - sele.GetNodalCoordsOld(1, i));
  }

  for (int i = 0; i < ncol; ++i)
  {
    mjumpv[0] += mval[i] * (mele.GetNodalCoords(0, i) - mele.GetNodalCoordsOld(0, i));
    mjumpv[1] += mval[i] * (mele.GetNodalCoords(1, i) - mele.GetNodalCoordsOld(1, i));
  }

  // normalize interpolated GP tangent back to length 1.0 !!!
  tanlength = sqrt(tanv[0] * tanv[0] + tanv[1] * tanv[1]);
  if (tanlength < 1.0e-12) dserror("ERROR: IntegrateAndDerivSegment: Divide by zero!");

  for (int i = 0; i < 2; i++) tanv[i] /= tanlength;

  // jump
  jumpv[0] = sjumpv[0] - mjumpv[0];
  jumpv[1] = sjumpv[1] - mjumpv[1];

  // multiply with tangent
  // value of relative tangential jump
  for (int i = 0; i < 2; ++i) jumpvalv[0] += tanv[i] * jumpv[i];


  // *****************************************************************************
  // add everything to dslipgp                                                   *
  // *****************************************************************************
  for (int i = 0; i < nrow; ++i)
  {
    GEN::pairedvector<int, double>& dmap_txsl_i =
        dynamic_cast<CONTACT::CoNode*>(snodes[i])->CoData().GetDerivTxi()[0];
    GEN::pairedvector<int, double>& dmap_tysl_i =
        dynamic_cast<CONTACT::CoNode*>(snodes[i])->CoData().GetDerivTxi()[1];

    for (_CI p = dmap_txsl_i.begin(); p != dmap_txsl_i.end(); ++p)
      dmap_txsl_gp[p->first] += sval[i] * (p->second);
    for (_CI p = dmap_tysl_i.begin(); p != dmap_tysl_i.end(); ++p)
      dmap_tysl_gp[p->first] += sval[i] * (p->second);

    for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
    {
      double valx = sderiv(i, 0) * dynamic_cast<CONTACT::CoNode*>(snodes[i])->CoData().txi()[0];
      dmap_txsl_gp[p->first] += valx * (p->second);
      double valy = sderiv(i, 0) * dynamic_cast<CONTACT::CoNode*>(snodes[i])->CoData().txi()[1];
      dmap_tysl_gp[p->first] += valy * (p->second);
    }
  }

  // build directional derivative of slave GP tagent (unit)
  GEN::pairedvector<int, double> dmap_txsl_gp_unit(ncol * ndof + linsize);
  GEN::pairedvector<int, double> dmap_tysl_gp_unit(ncol * ndof + linsize);

  const double llv = tanlength * tanlength;
  const double linv = 1.0 / tanlength;
  const double lllinv = 1.0 / (tanlength * tanlength * tanlength);
  const double sxsxv = tanv[0] * tanv[0] * llv;
  const double sxsyv = tanv[0] * tanv[1] * llv;
  const double sysyv = tanv[1] * tanv[1] * llv;

  for (_CI p = dmap_txsl_gp.begin(); p != dmap_txsl_gp.end(); ++p)
  {
    dmap_txsl_gp_unit[p->first] += linv * (p->second);
    dmap_txsl_gp_unit[p->first] -= lllinv * sxsxv * (p->second);
    dmap_tysl_gp_unit[p->first] -= lllinv * sxsyv * (p->second);
  }

  for (_CI p = dmap_tysl_gp.begin(); p != dmap_tysl_gp.end(); ++p)
  {
    dmap_tysl_gp_unit[p->first] += linv * (p->second);
    dmap_tysl_gp_unit[p->first] -= lllinv * sysyv * (p->second);
    dmap_txsl_gp_unit[p->first] -= lllinv * sxsyv * (p->second);
  }

  for (_CI p = dmap_txsl_gp_unit.begin(); p != dmap_txsl_gp_unit.end(); ++p)
    dslipgp[p->first] += jumpv[0] * (p->second);

  for (_CI p = dmap_tysl_gp_unit.begin(); p != dmap_tysl_gp_unit.end(); ++p)
    dslipgp[p->first] += jumpv[1] * (p->second);

  // coord lin
  for (int z = 0; z < nrow; ++z)
  {
    FriNode* snode = dynamic_cast<FriNode*>(snodes[z]);
    for (int k = 0; k < 2; ++k)
    {
      dslipgp[snode->Dofs()[k]] += sval[z] * tanv[k];

      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
        dslipgp[p->first] += tanv[k] * sderiv(z, 0) *
                             (sele.GetNodalCoords(k, z) - sele.GetNodalCoordsOld(k, z)) *
                             (p->second);
    }
  }

  for (int z = 0; z < ncol; ++z)
  {
    FriNode* mnode = dynamic_cast<FriNode*>(mnodes[z]);
    for (int k = 0; k < 2; ++k)
    {
      dslipgp[mnode->Dofs()[k]] -= mval[z] * tanv[k];

      for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
        dslipgp[p->first] -= tanv[k] * mderiv(z, 0) *
                             (mele.GetNodalCoords(k, z) - mele.GetNodalCoordsOld(k, z)) *
                             (p->second);
    }
  }

  // ***************************
  // Add to node!
  for (int j = 0; j < nrow; ++j)
  {
    FriNode* snode = dynamic_cast<FriNode*>(snodes[j]);
    if (snode->IsOnBoundorCE()) continue;

    double prod = lmval[j] * jumpvalv[0] * jac * wgt;

    // add current Gauss point's contribution to jump
    snode->AddJumpValue(prod, 0);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Compute entries for slip increment at GP                 farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_3D_SlipIncr(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
    LINALG::SerialDenseVector& lmval, LINALG::SerialDenseMatrix& sderiv,
    LINALG::SerialDenseMatrix& mderiv, double& jac, double& wgt, double* jumpvalv,
    const std::vector<GEN::pairedvector<int, double>>& dsxigp,
    const std::vector<GEN::pairedvector<int, double>>& dmxigp,
    std::vector<GEN::pairedvector<int, double>>& dslipgp)
{
  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();
  if (!snodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");
  if (!mnodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // number of nodes (slave, master)
  const int nrow = sele.NumNode();
  const int ncol = mele.NumNode();
  const int ndof = Dim();

  int linsize = 0;
  for (int i = 0; i < nrow; ++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*>(snodes[i]);
    linsize += cnode->GetLinsize();
  }

  // build interpolation of slave GP normal and coordinates
  double sjumpv[3] = {0.0, 0.0, 0.0};
  double mjumpv[3] = {0.0, 0.0, 0.0};
  double jumpv[3] = {0.0, 0.0, 0.0};
  double tanv1[3] = {0.0, 0.0, 0.0};
  double tanv2[3] = {0.0, 0.0, 0.0};

  double jumpvalv1 = 0.0;
  double jumpvalv2 = 0.0;
  double tanlength1 = 0.0;
  double tanlength2 = 0.0;

  // LIN OF TANGENT
  std::map<int, double> dmap_txsl_gp;
  std::map<int, double> dmap_tysl_gp;

  for (int i = 0; i < nrow; ++i)
  {
    CONTACT::CoNode* myconode = dynamic_cast<CONTACT::CoNode*>(snodes[i]);

    // nodal tangent interpolation
    tanv1[0] += sval[i] * myconode->CoData().txi()[0];
    tanv1[1] += sval[i] * myconode->CoData().txi()[1];
    tanv1[2] += sval[i] * myconode->CoData().txi()[2];

    tanv2[0] += sval[i] * myconode->CoData().teta()[0];
    tanv2[1] += sval[i] * myconode->CoData().teta()[1];
    tanv2[2] += sval[i] * myconode->CoData().teta()[2];
    // delta D
    sjumpv[0] += sval[i] * (sele.GetNodalCoords(0, i) - sele.GetNodalCoordsOld(0, i));
    sjumpv[1] += sval[i] * (sele.GetNodalCoords(1, i) - sele.GetNodalCoordsOld(1, i));
    sjumpv[2] += sval[i] * (sele.GetNodalCoords(2, i) - sele.GetNodalCoordsOld(2, i));
  }

  for (int i = 0; i < ncol; ++i)
  {
    mjumpv[0] += mval[i] * (mele.GetNodalCoords(0, i) - mele.GetNodalCoordsOld(0, i));
    mjumpv[1] += mval[i] * (mele.GetNodalCoords(1, i) - mele.GetNodalCoordsOld(1, i));
    mjumpv[2] += mval[i] * (mele.GetNodalCoords(2, i) - mele.GetNodalCoordsOld(2, i));
  }

  // normalize interpolated GP tangent back to length 1.0 !!!
  tanlength1 = sqrt(tanv1[0] * tanv1[0] + tanv1[1] * tanv1[1] + tanv1[2] * tanv1[2]);
  if (tanlength1 < 1.0e-12) dserror("ERROR: IntegrateAndDerivSegment: Divide by zero!");

  tanlength2 = sqrt(tanv2[0] * tanv2[0] + tanv2[1] * tanv2[1] + tanv2[2] * tanv2[2]);
  if (tanlength2 < 1.0e-12) dserror("ERROR: IntegrateAndDerivSegment: Divide by zero!");

  for (int i = 0; i < 3; i++)
  {
    tanv1[i] /= tanlength1;
    tanv2[i] /= tanlength2;
  }

  // jump
  jumpv[0] = sjumpv[0] - mjumpv[0];
  jumpv[1] = sjumpv[1] - mjumpv[1];
  jumpv[2] = sjumpv[2] - mjumpv[2];

  // multiply with tangent
  // value of relative tangential jump
  for (int i = 0; i < 3; ++i)
  {
    jumpvalv1 += tanv1[i] * jumpv[i];
    jumpvalv2 += tanv2[i] * jumpv[i];
  }

  jumpvalv[0] = jumpvalv1;
  jumpvalv[1] = jumpvalv2;

  // ***************************
  // Add to node!
  for (int j = 0; j < nrow; ++j)
  {
    FriNode* snode = dynamic_cast<FriNode*>(snodes[j]);

    if (snode->IsOnBoundorCE()) continue;

    double prod1 = lmval[j] * jumpvalv1 * jac * wgt;
    double prod2 = lmval[j] * jumpvalv2 * jac * wgt;

    // add current Gauss point's contribution to gseg
    snode->AddJumpValue(prod1, 0);
    snode->AddJumpValue(prod2, 1);
  }

  //************* LIN TANGENT TXI *********************
  // build directional derivative of slave GP txi (non-unit)
  GEN::pairedvector<int, double> dmap_txix_gp(linsize + ncol * ndof);
  GEN::pairedvector<int, double> dmap_txiy_gp(linsize + ncol * ndof);
  GEN::pairedvector<int, double> dmap_txiz_gp(linsize + ncol * ndof);

  // slave GP txi (non-unit)
  GEN::pairedvector<int, double> dmap_txix_gp_unit(linsize + ncol * ndof);
  GEN::pairedvector<int, double> dmap_txiy_gp_unit(linsize + ncol * ndof);
  GEN::pairedvector<int, double> dmap_txiz_gp_unit(linsize + ncol * ndof);

  for (int i = 0; i < nrow; ++i)
  {
    GEN::pairedvector<int, double>& dmap_txsl_i =
        dynamic_cast<CONTACT::CoNode*>(snodes[i])->CoData().GetDerivTxi()[0];
    GEN::pairedvector<int, double>& dmap_tysl_i =
        dynamic_cast<CONTACT::CoNode*>(snodes[i])->CoData().GetDerivTxi()[1];
    GEN::pairedvector<int, double>& dmap_tzsl_i =
        dynamic_cast<CONTACT::CoNode*>(snodes[i])->CoData().GetDerivTxi()[2];

    for (_CI p = dmap_txsl_i.begin(); p != dmap_txsl_i.end(); ++p)
      dmap_txix_gp[p->first] += sval[i] * (p->second);
    for (_CI p = dmap_tysl_i.begin(); p != dmap_tysl_i.end(); ++p)
      dmap_txiy_gp[p->first] += sval[i] * (p->second);
    for (_CI p = dmap_tzsl_i.begin(); p != dmap_tzsl_i.end(); ++p)
      dmap_txiz_gp[p->first] += sval[i] * (p->second);

    const double txi_x = dynamic_cast<CONTACT::CoNode*>(snodes[i])->CoData().txi()[0];
    const double txi_y = dynamic_cast<CONTACT::CoNode*>(snodes[i])->CoData().txi()[1];
    const double txi_z = dynamic_cast<CONTACT::CoNode*>(snodes[i])->CoData().txi()[2];

    for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
    {
      double valx = sderiv(i, 0) * txi_x;
      dmap_txix_gp[p->first] += valx * (p->second);
      double valy = sderiv(i, 0) * txi_y;
      dmap_txiy_gp[p->first] += valy * (p->second);
      double valz = sderiv(i, 0) * txi_z;
      dmap_txiz_gp[p->first] += valz * (p->second);
    }

    for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
    {
      double valx = sderiv(i, 1) * txi_x;
      dmap_txix_gp[p->first] += valx * (p->second);
      double valy = sderiv(i, 1) * txi_y;
      dmap_txiy_gp[p->first] += valy * (p->second);
      double valz = sderiv(i, 1) * txi_z;
      dmap_txiz_gp[p->first] += valz * (p->second);
    }
  }

  // build directional derivative of slave GP txi (unit)
  const double ll1 = tanlength1 * tanlength1;
  const double linv1 = 1.0 / tanlength1;
  const double lllinv1 = 1.0 / (tanlength1 * tanlength1 * tanlength1);
  const double sxsx1 = tanv1[0] * tanv1[0] * ll1;
  const double sxsy1 = tanv1[0] * tanv1[1] * ll1;
  const double sxsz1 = tanv1[0] * tanv1[2] * ll1;
  const double sysy1 = tanv1[1] * tanv1[1] * ll1;
  const double sysz1 = tanv1[1] * tanv1[2] * ll1;
  const double szsz1 = tanv1[2] * tanv1[2] * ll1;

  for (_CI p = dmap_txix_gp.begin(); p != dmap_txix_gp.end(); ++p)
  {
    dmap_txix_gp_unit[p->first] += linv1 * (p->second);
    dmap_txix_gp_unit[p->first] -= lllinv1 * sxsx1 * (p->second);
    dmap_txiy_gp_unit[p->first] -= lllinv1 * sxsy1 * (p->second);
    dmap_txiz_gp_unit[p->first] -= lllinv1 * sxsz1 * (p->second);
  }

  for (_CI p = dmap_txiy_gp.begin(); p != dmap_txiy_gp.end(); ++p)
  {
    dmap_txiy_gp_unit[p->first] += linv1 * (p->second);
    dmap_txiy_gp_unit[p->first] -= lllinv1 * sysy1 * (p->second);
    dmap_txix_gp_unit[p->first] -= lllinv1 * sxsy1 * (p->second);
    dmap_txiz_gp_unit[p->first] -= lllinv1 * sysz1 * (p->second);
  }

  for (_CI p = dmap_txiz_gp.begin(); p != dmap_txiz_gp.end(); ++p)
  {
    dmap_txiz_gp_unit[p->first] += linv1 * (p->second);
    dmap_txiz_gp_unit[p->first] -= lllinv1 * szsz1 * (p->second);
    dmap_txix_gp_unit[p->first] -= lllinv1 * sxsz1 * (p->second);
    dmap_txiy_gp_unit[p->first] -= lllinv1 * sysz1 * (p->second);
  }


  //************* LIN TANGENT TETA *********************
  // build directional derivative of slave GP teta (non-unit)
  GEN::pairedvector<int, double> dmap_tetax_gp(linsize + ncol * ndof);
  GEN::pairedvector<int, double> dmap_tetay_gp(linsize + ncol * ndof);
  GEN::pairedvector<int, double> dmap_tetaz_gp(linsize + ncol * ndof);

  // slave GP teta (unit)
  GEN::pairedvector<int, double> dmap_tetax_gp_unit(linsize + ncol * ndof);
  GEN::pairedvector<int, double> dmap_tetay_gp_unit(linsize + ncol * ndof);
  GEN::pairedvector<int, double> dmap_tetaz_gp_unit(linsize + ncol * ndof);

  for (int i = 0; i < nrow; ++i)
  {
    GEN::pairedvector<int, double>& dmap_txsl_i =
        dynamic_cast<CONTACT::CoNode*>(snodes[i])->CoData().GetDerivTeta()[0];
    GEN::pairedvector<int, double>& dmap_tysl_i =
        dynamic_cast<CONTACT::CoNode*>(snodes[i])->CoData().GetDerivTeta()[1];
    GEN::pairedvector<int, double>& dmap_tzsl_i =
        dynamic_cast<CONTACT::CoNode*>(snodes[i])->CoData().GetDerivTeta()[2];

    for (_CI p = dmap_txsl_i.begin(); p != dmap_txsl_i.end(); ++p)
      dmap_tetax_gp[p->first] += sval[i] * (p->second);
    for (_CI p = dmap_tysl_i.begin(); p != dmap_tysl_i.end(); ++p)
      dmap_tetay_gp[p->first] += sval[i] * (p->second);
    for (_CI p = dmap_tzsl_i.begin(); p != dmap_tzsl_i.end(); ++p)
      dmap_tetaz_gp[p->first] += sval[i] * (p->second);

    const double teta_x = dynamic_cast<CONTACT::CoNode*>(snodes[i])->CoData().teta()[0];
    const double teta_y = dynamic_cast<CONTACT::CoNode*>(snodes[i])->CoData().teta()[1];
    const double teta_z = dynamic_cast<CONTACT::CoNode*>(snodes[i])->CoData().teta()[2];

    for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
    {
      double valx = sderiv(i, 0) * teta_x;
      dmap_tetax_gp[p->first] += valx * (p->second);
      double valy = sderiv(i, 0) * teta_y;
      dmap_tetay_gp[p->first] += valy * (p->second);
      double valz = sderiv(i, 0) * teta_z;
      dmap_tetaz_gp[p->first] += valz * (p->second);
    }

    for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
    {
      double valx = sderiv(i, 1) * teta_x;
      dmap_tetax_gp[p->first] += valx * (p->second);
      double valy = sderiv(i, 1) * teta_y;
      dmap_tetay_gp[p->first] += valy * (p->second);
      double valz = sderiv(i, 1) * teta_z;
      dmap_tetaz_gp[p->first] += valz * (p->second);
    }
  }

  // build directional derivative of slave GP teta (unit)
  const double ll2 = tanlength2 * tanlength2;
  const double linv2 = 1.0 / tanlength2;
  const double lllinv2 = 1.0 / (tanlength2 * tanlength2 * tanlength2);
  const double sxsx2 = tanv2[0] * tanv2[0] * ll2;
  const double sxsy2 = tanv2[0] * tanv2[1] * ll2;
  const double sxsz2 = tanv2[0] * tanv2[2] * ll2;
  const double sysy2 = tanv2[1] * tanv2[1] * ll2;
  const double sysz2 = tanv2[1] * tanv2[2] * ll2;
  const double szsz2 = tanv2[2] * tanv2[2] * ll2;

  for (_CI p = dmap_tetax_gp.begin(); p != dmap_tetax_gp.end(); ++p)
  {
    dmap_tetax_gp_unit[p->first] += linv2 * (p->second);
    dmap_tetax_gp_unit[p->first] -= lllinv2 * sxsx2 * (p->second);
    dmap_tetay_gp_unit[p->first] -= lllinv2 * sxsy2 * (p->second);
    dmap_tetaz_gp_unit[p->first] -= lllinv2 * sxsz2 * (p->second);
  }

  for (_CI p = dmap_tetay_gp.begin(); p != dmap_tetay_gp.end(); ++p)
  {
    dmap_tetay_gp_unit[p->first] += linv2 * (p->second);
    dmap_tetay_gp_unit[p->first] -= lllinv2 * sysy2 * (p->second);
    dmap_tetax_gp_unit[p->first] -= lllinv2 * sxsy2 * (p->second);
    dmap_tetaz_gp_unit[p->first] -= lllinv2 * sysz2 * (p->second);
  }

  for (_CI p = dmap_tetaz_gp.begin(); p != dmap_tetaz_gp.end(); ++p)
  {
    dmap_tetaz_gp_unit[p->first] += linv2 * (p->second);
    dmap_tetaz_gp_unit[p->first] -= lllinv2 * szsz2 * (p->second);
    dmap_tetax_gp_unit[p->first] -= lllinv2 * sxsz2 * (p->second);
    dmap_tetay_gp_unit[p->first] -= lllinv2 * sysz2 * (p->second);
  }

  // TXI
  for (_CI p = dmap_txix_gp_unit.begin(); p != dmap_txix_gp_unit.end(); ++p)
    dslipgp[0][p->first] += jumpv[0] * (p->second);

  for (_CI p = dmap_txiy_gp_unit.begin(); p != dmap_txiy_gp_unit.end(); ++p)
    dslipgp[0][p->first] += jumpv[1] * (p->second);

  for (_CI p = dmap_txiz_gp_unit.begin(); p != dmap_txiz_gp_unit.end(); ++p)
    dslipgp[0][p->first] += jumpv[2] * (p->second);

  // TETA
  for (_CI p = dmap_tetax_gp_unit.begin(); p != dmap_tetax_gp_unit.end(); ++p)
    dslipgp[1][p->first] += jumpv[0] * (p->second);

  for (_CI p = dmap_tetay_gp_unit.begin(); p != dmap_tetay_gp_unit.end(); ++p)
    dslipgp[1][p->first] += jumpv[1] * (p->second);

  for (_CI p = dmap_tetaz_gp_unit.begin(); p != dmap_tetaz_gp_unit.end(); ++p)
    dslipgp[1][p->first] += jumpv[2] * (p->second);

  // coord lin
  for (int z = 0; z < nrow; ++z)
  {
    FriNode* snode = dynamic_cast<FriNode*>(snodes[z]);

    for (int k = 0; k < 3; ++k)
    {
      dslipgp[0][snode->Dofs()[k]] += sval[z] * tanv1[k];
      dslipgp[1][snode->Dofs()[k]] += sval[z] * tanv2[k];

      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
      {
        dslipgp[0][p->first] += tanv1[k] * sderiv(z, 0) *
                                (sele.GetNodalCoords(k, z) - sele.GetNodalCoordsOld(k, z)) *
                                (p->second);
        dslipgp[1][p->first] += tanv2[k] * sderiv(z, 0) *
                                (sele.GetNodalCoords(k, z) - sele.GetNodalCoordsOld(k, z)) *
                                (p->second);
      }

      for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
      {
        dslipgp[0][p->first] += tanv1[k] * sderiv(z, 1) *
                                (sele.GetNodalCoords(k, z) - sele.GetNodalCoordsOld(k, z)) *
                                (p->second);
        dslipgp[1][p->first] += tanv2[k] * sderiv(z, 1) *
                                (sele.GetNodalCoords(k, z) - sele.GetNodalCoordsOld(k, z)) *
                                (p->second);
      }
    }
  }

  for (int z = 0; z < ncol; ++z)
  {
    FriNode* mnode = dynamic_cast<FriNode*>(mnodes[z]);

    for (int k = 0; k < 3; ++k)
    {
      dslipgp[0][mnode->Dofs()[k]] -= mval[z] * tanv1[k];
      dslipgp[1][mnode->Dofs()[k]] -= mval[z] * tanv2[k];

      for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
      {
        dslipgp[0][p->first] -= tanv1[k] * mderiv(z, 0) *
                                (mele.GetNodalCoords(k, z) - mele.GetNodalCoordsOld(k, z)) *
                                (p->second);
        dslipgp[1][p->first] -= tanv2[k] * mderiv(z, 0) *
                                (mele.GetNodalCoords(k, z) - mele.GetNodalCoordsOld(k, z)) *
                                (p->second);
      }

      for (_CI p = dmxigp[1].begin(); p != dmxigp[1].end(); ++p)
      {
        dslipgp[0][p->first] -= tanv1[k] * mderiv(z, 1) *
                                (mele.GetNodalCoords(k, z) - mele.GetNodalCoordsOld(k, z)) *
                                (p->second);
        dslipgp[1][p->first] -= tanv2[k] * mderiv(z, 1) *
                                (mele.GetNodalCoords(k, z) - mele.GetNodalCoordsOld(k, z)) *
                                (p->second);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Compute slipincr lin at GP                               farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_2D_SlipIncr_Lin(int& iter, MORTAR::MortarElement& sele,
    LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& lmderiv, double& jac, double& wgt,
    double* jumpvalv, const std::vector<GEN::pairedvector<int, double>>& dsxigp,
    const GEN::pairedvector<int, double>& dslipgp, const GEN::pairedvector<int, double>& derivjac,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap)
{
  DRT::Node** snodes = sele.Nodes();

  const int nrow = sele.NumNode();
  double fac = 0.0;

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  FriNode* snode = dynamic_cast<FriNode*>(snodes[iter]);
  if (snode->IsOnBoundorCE()) return;
  ;

  // get the corresponding map as a reference
  std::map<int, double>& djumpmap = snode->FriData().GetDerivVarJump()[0];

  // (1) Lin(Phi) - dual shape functions
  if (ShapeFcn() == INPAR::MORTAR::shape_dual)
  {
    for (int m = 0; m < nrow; ++m)
    {
      fac = wgt * sval[m] * jumpvalv[0] * jac;
      for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
           p != dualmap.end(); ++p)
        djumpmap[p->first] += fac * (p->second)(iter, m);
    }
  }

  // (2) Lin(Phi) - slave GP coordinates
  fac = wgt * lmderiv(iter, 0) * jumpvalv[0] * jac;
  for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
    djumpmap[p->first] += fac * (p->second);

  // (3) Lin(g) - gap function
  fac = wgt * lmval[iter] * jac;
  for (_CI p = dslipgp.begin(); p != dslipgp.end(); ++p) djumpmap[p->first] += fac * (p->second);

  // (4) Lin(dxdsxi) - slave GP Jacobian
  fac = wgt * lmval[iter] * jumpvalv[0];
  for (_CI p = derivjac.begin(); p != derivjac.end(); ++p) djumpmap[p->first] += fac * (p->second);

  return;
}

/*----------------------------------------------------------------------*
 |  Compute slipincr lin at   GP                             farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_3D_SlipIncr_Lin(int& iter, MORTAR::MortarElement& sele,
    LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& lmderiv, double& jac, double& wgt,
    double* jumpvalv, const GEN::pairedvector<int, double>& jacintcellmap,
    const std::vector<GEN::pairedvector<int, double>>& dslipgp,
    const std::vector<GEN::pairedvector<int, double>>& dsxigp,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap)
{
  DRT::Node** snodes = sele.Nodes();

  double nrow = sele.NumNode();

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  FriNode* snode = dynamic_cast<FriNode*>(snodes[iter]);

  if (snode->IsOnBoundorCE()) return;

  // get the corresponding map as a reference
  std::map<int, double>& djumpmap1 = snode->FriData().GetDerivVarJump()[0];
  std::map<int, double>& djumpmap2 = snode->FriData().GetDerivVarJump()[1];

  double fac1 = 0.0;
  double fac2 = 0.0;

  // (1) Lin(Phi) - dual shape functions
  if (ShapeFcn() == INPAR::MORTAR::shape_dual)
  {
    for (int m = 0; m < nrow; ++m)
    {
      fac1 = wgt * sval[m] * jumpvalv[0] * jac;
      fac2 = wgt * sval[m] * jumpvalv[1] * jac;

      for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
           p != dualmap.end(); ++p)
      {
        djumpmap1[p->first] += fac1 * (p->second)(iter, m);
        djumpmap2[p->first] += fac2 * (p->second)(iter, m);
      }
    }
  }

  // (2) Lin(Phi) - slave GP coordinates --> because of duality
  fac1 = wgt * lmderiv(iter, 0) * jumpvalv[0] * jac;
  fac2 = wgt * lmderiv(iter, 0) * jumpvalv[1] * jac;
  for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
  {
    djumpmap1[p->first] += fac1 * (p->second);
    djumpmap2[p->first] += fac2 * (p->second);
  }

  fac1 = wgt * lmderiv(iter, 1) * jumpvalv[0] * jac;
  fac2 = wgt * lmderiv(iter, 1) * jumpvalv[1] * jac;
  for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
  {
    djumpmap1[p->first] += fac1 * (p->second);
    djumpmap2[p->first] += fac2 * (p->second);
  }

  // (3) Lin(w) - wear function
  fac1 = wgt * lmval[iter] * jac;
  for (_CI p = dslipgp[0].begin(); p != dslipgp[0].end(); ++p)
    djumpmap1[p->first] += fac1 * (p->second);
  for (_CI p = dslipgp[1].begin(); p != dslipgp[1].end(); ++p)
    djumpmap2[p->first] += fac1 * (p->second);

  // (5) Lin(dxdsxi) - slave GP Jacobian
  fac1 = wgt * lmval[iter] * jumpvalv[0];
  fac2 = wgt * lmval[iter] * jumpvalv[1];
  for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
  {
    djumpmap1[p->first] += fac1 * (p->second);
    djumpmap2[p->first] += fac2 * (p->second);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Lin wear for impl. algor.                                farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_2D_Wear_Lin(int& iter, MORTAR::MortarElement& sele,
    LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& lmderiv, double& jac, double* gpn,
    double& wgt, double& wearval, double* jumpval, const GEN::pairedvector<int, double>& dweargp,
    const GEN::pairedvector<int, double>& derivjac,
    const std::vector<GEN::pairedvector<int, double>>& dsxigp,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap)
{
  const double wcoeff = wearcoeff_ + wearcoeffm_;

  double facw = 0.0;
  const int nrow = sele.NumNode();

  DRT::Node** snodes = sele.Nodes();

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // get the corresponding map as a reference
  CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(snodes[iter]);

  std::map<int, double>& dwmap = cnode->CoData().GetDerivW();

  if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
  {
    // (1) Lin(Phi) - dual shape functions resulting from weighting the wear
    // we use std. shape functions for shape_petrovgalerkin

    // (2) Lin(Phi) - slave GP coordinates --> because of duality
    facw = wcoeff * wgt * sderiv(iter, 0) * wearval * jac;
    for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
      dwmap[p->first] += facw * (p->second);

    // (3) Lin(w) - wear function
    facw = wcoeff * wgt * sval[iter] * jac;
    for (_CI p = dweargp.begin(); p != dweargp.end(); ++p) dwmap[p->first] += facw * (p->second);

    // (4) Lin(dxdsxi) - slave GP Jacobian
    facw = wcoeff * wgt * sval[iter] * wearval;
    for (_CI p = derivjac.begin(); p != derivjac.end(); ++p) dwmap[p->first] += facw * (p->second);
  }
  else  // no petrov_galerkin
  {
    // (1) Lin(Phi) - dual shape functions resulting from weighting the wear
    if (ShapeFcn() == INPAR::MORTAR::shape_dual)
    {
      for (int m = 0; m < nrow; ++m)
      {
        facw = wcoeff * wgt * sval[m] * wearval * jac;
        for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
             p != dualmap.end(); ++p)
          dwmap[p->first] += facw * (p->second)(iter, m);
      }
    }

    // (2) Lin(Phi) - slave GP coordinates --> because of duality
    facw = wcoeff * wgt * lmderiv(iter, 0) * wearval * jac;
    for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
      dwmap[p->first] += facw * (p->second);

    // (3) Lin(w) - wear function
    facw = wcoeff * wgt * lmval[iter] * jac;
    for (_CI p = dweargp.begin(); p != dweargp.end(); ++p) dwmap[p->first] += facw * (p->second);

    // (4) Lin(dxdsxi) - slave GP Jacobian
    facw = wcoeff * wgt * lmval[iter] * wearval;
    for (_CI p = derivjac.begin(); p != derivjac.end(); ++p) dwmap[p->first] += facw * (p->second);
  }

  //****************************************************************
  // LIN WEAR W.R.T. LM
  //****************************************************************
  std::map<int, double>& dwlmmap = cnode->CoData().GetDerivWlm();

  for (int bl = 0; bl < nrow; ++bl)
  {
    MORTAR::MortarNode* wearnode = dynamic_cast<MORTAR::MortarNode*>(snodes[bl]);

    if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
    {
      dwlmmap[wearnode->Dofs()[0]] +=
          wcoeff * sval[iter] * jac * wgt * abs(jumpval[0]) * gpn[0] * lmval[bl];
      dwlmmap[wearnode->Dofs()[1]] +=
          wcoeff * sval[iter] * jac * wgt * abs(jumpval[0]) * gpn[1] * lmval[bl];
    }
    else
    {
      dwlmmap[wearnode->Dofs()[0]] +=
          wcoeff * lmval[iter] * jac * wgt * abs(jumpval[0]) * gpn[0] * lmval[bl];
      dwlmmap[wearnode->Dofs()[1]] +=
          wcoeff * lmval[iter] * jac * wgt * abs(jumpval[0]) * gpn[1] * lmval[bl];
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Lin wear for impl. algor.                                farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_3D_Wear_Lin(int& iter, MORTAR::MortarElement& sele,
    LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& lmderiv, double& jac, double* gpn,
    double& wgt, double& wearval, double* jumpval, const GEN::pairedvector<int, double>& dweargp,
    const GEN::pairedvector<int, double>& jacintcellmap,
    const std::vector<GEN::pairedvector<int, double>>& dsxigp,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap)
{
  double facw = 0.0;
  const int nrow = sele.NumNode();

  DRT::Node** snodes = sele.Nodes();

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator CI;

  // get the corresponding map as a reference
  CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(snodes[iter]);

  // get the corresponding map as a reference
  std::map<int, double>& dwmap = cnode->CoData().GetDerivW();

  if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
  {
    // (1) Lin(Phi) - dual shape functions resulting from weighting the wear
    // --

    // (2) Lin(Phi) - slave GP coordinates --> because of duality
    facw = wearcoeff_ * wgt * sderiv(iter, 0) * wearval * jac;
    for (CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p) dwmap[p->first] += facw * (p->second);

    facw = wearcoeff_ * wgt * sderiv(iter, 1) * wearval * jac;
    for (CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p) dwmap[p->first] += facw * (p->second);

    // (3) Lin(w) - wear function
    facw = wearcoeff_ * wgt * sval[iter] * jac;
    for (CI p = dweargp.begin(); p != dweargp.end(); ++p) dwmap[p->first] += facw * (p->second);

    // (5) Lin(dxdsxi) - slave GP Jacobian
    facw = wearcoeff_ * wgt * sval[iter] * wearval;
    for (CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
      dwmap[p->first] += facw * (p->second);
  }
  else  // no pg
  {
    // (1) Lin(Phi) - dual shape functions resulting from weighting the wear
    if (ShapeFcn() == INPAR::MORTAR::shape_dual)
    {
      for (int m = 0; m < nrow; ++m)
      {
        facw = wearcoeff_ * wgt * sval[m] * wearval * jac;
        for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
             p != dualmap.end(); ++p)
          dwmap[p->first] += facw * (p->second)(iter, m);
      }
    }

    // (2) Lin(Phi) - slave GP coordinates --> because of duality
    facw = wearcoeff_ * wgt * lmderiv(iter, 0) * wearval * jac;
    for (CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p) dwmap[p->first] += facw * (p->second);

    facw = wearcoeff_ * wgt * lmderiv(iter, 1) * wearval * jac;
    for (CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p) dwmap[p->first] += facw * (p->second);

    // (3) Lin(w) - wear function
    facw = wearcoeff_ * wgt * lmval[iter] * jac;
    for (CI p = dweargp.begin(); p != dweargp.end(); ++p) dwmap[p->first] += facw * (p->second);

    // (5) Lin(dxdsxi) - slave GP Jacobian
    facw = wearcoeff_ * wgt * lmval[iter] * wearval;
    for (CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
      dwmap[p->first] += facw * (p->second);
  }


  //****************************************************************
  // LIN WEAR W.R.T. LM
  //****************************************************************
  std::map<int, double>& dwlmmap = cnode->CoData().GetDerivWlm();

  for (int bl = 0; bl < nrow; ++bl)
  {
    MORTAR::MortarNode* wearnode = dynamic_cast<MORTAR::MortarNode*>(snodes[bl]);

    if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
    {
      dwlmmap[wearnode->Dofs()[0]] +=
          wearcoeff_ * sval[iter] * jac * wgt * abs(jumpval[0]) * gpn[0] * lmval[bl];
      dwlmmap[wearnode->Dofs()[1]] +=
          wearcoeff_ * sval[iter] * jac * wgt * abs(jumpval[0]) * gpn[1] * lmval[bl];
      dwlmmap[wearnode->Dofs()[2]] +=
          wearcoeff_ * sval[iter] * jac * wgt * abs(jumpval[0]) * gpn[2] * lmval[bl];
    }
    else
    {
      dwlmmap[wearnode->Dofs()[0]] +=
          wearcoeff_ * lmval[iter] * jac * wgt * abs(jumpval[0]) * gpn[0] * lmval[bl];
      dwlmmap[wearnode->Dofs()[1]] +=
          wearcoeff_ * lmval[iter] * jac * wgt * abs(jumpval[0]) * gpn[1] * lmval[bl];
      dwlmmap[wearnode->Dofs()[2]] +=
          wearcoeff_ * lmval[iter] * jac * wgt * abs(jumpval[0]) * gpn[2] * lmval[bl];
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Compute entries for poro normal coupling condition      07/14 ager  |
 *----------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_NCOUP_DERIV(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
    LINALG::SerialDenseVector& lmval, LINALG::SerialDenseMatrix& sderiv,
    LINALG::SerialDenseMatrix& mderiv, double* ncoup, double* gpn, double& jac, double& wgt,
    double* gpcoord, const std::vector<GEN::pairedvector<int, double>>& dsxigp,
    const std::vector<GEN::pairedvector<int, double>>& dmxigp, std::map<int, double>& dncoupgp,
    std::map<int, double>& dvelncoupgp, std::map<int, double>& dpresncoupgp,
    std::vector<GEN::pairedvector<int, double>>& dnmap_unit, bool quadratic, int nintrow)
{
  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();
  if (!snodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");
  if (!mnodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // bool to decide if fluid quantities and structural velocities
  // exist for slave side on CoNode level and contribute to porofluid meshtying
  bool slaveporo = (sele.PhysType() == MORTAR::MortarElement::poro);
  // bool to decide if fluid quantities and structural velocities
  // exist for master side on CoNode level and contribute to porofluid meshtying
  bool masterporo = (mele.PhysType() == MORTAR::MortarElement::poro);

  if (!slaveporo and masterporo)
    dserror(
        "poroelastic mesh tying method needs the slave side to be poroelastic (invert master/slave "
        "definition)");
  // see CONTACT::PoroLagrangeStrategy::PoroLagrangeStrategy for comment

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  // if(twosided) //needed only for master side
  int ncol = mele.NumNode();  // not used for onesided porocontact!!!

  // get fluid velocities in GP
  double sgpfvel[3] = {0.0, 0.0, 0.0};
  double sgpsvel[3] = {0.0, 0.0, 0.0};

  double mgpfvel[3] = {0.0, 0.0, 0.0};
  double mgpsvel[3] = {0.0, 0.0, 0.0};

  // fill local vectors for corresponding problem dimension
  for (int k = 0; k < Dim(); ++k)
  {
    if (slaveporo)
      for (int i = 0; i < nrow; ++i)
      {
        CoNode* mymrtrnode = dynamic_cast<CoNode*>(snodes[i]);
        sgpfvel[k] += sval[i] * mymrtrnode->CoPoroData().fvel()[k];
        sgpsvel[k] += sval[i] * mymrtrnode->CoPoroData().svel()[k];
      }

    if (masterporo)
    {
      for (int i = 0; i < ncol; ++i)
      {
        CoNode* mymrtrnode = dynamic_cast<CoNode*>(mnodes[i]);
        mgpfvel[k] += mval[i] * mymrtrnode->CoPoroData().fvel()[k];
        mgpsvel[k] += mval[i] * mymrtrnode->CoPoroData().svel()[k];
      }
    }
  }

  ////////////////////////////////!!!Calculate slave side Porosity!!!////////////////////////////
  // get J
  std::map<int, double>
      sJLin;  // slave map for linearizations of the Deformation Gradient Determinant
  // get fluid pressure in GP
  double sgpfpres = 0.0;
  // porosity to be written into
  double sporosity;
  double sdphi_dp;  // linearization terms to be written into
  double sdphi_dJ;  // linearization terms to be written into

  if (slaveporo)
  {
    const double sJ = DetDeformationGradient(sele, wgt, gpcoord, sJLin);

    for (int i = 0; i < nrow; ++i)
    {
      CoNode* mymrtrnode = dynamic_cast<CoNode*>(snodes[i]);
      sgpfpres += sval[i] * mymrtrnode->CoPoroData().fpres()[0];
    }
    Teuchos::ParameterList sparams;  // empty parameter list;

    Teuchos::RCP<MAT::StructPoro> sstructmat =
        Teuchos::rcp_dynamic_cast<MAT::StructPoro>(sele.ParentElement()->Material(0));
    if (sstructmat == Teuchos::null)
      sstructmat = Teuchos::rcp_dynamic_cast<MAT::StructPoro>(sele.ParentElement()->Material(1));
    if (sstructmat == Teuchos::null) dserror("Cast to StructPoro failed!");
    sstructmat->ComputeSurfPorosity(sparams, sgpfpres, sJ, sele.FaceParentNumber(),
        1,  // finally check what to do here Todo:
        sporosity,
        &sdphi_dp,  // linearization of phi w.r.t. pressure
        &sdphi_dJ,  // linearization of phi w.r.t. defo-grad determinant
        NULL,       // dphi_dJdp not needed
        NULL,       // dphi_dJJ not needed
        NULL,       // dphi_dpp not needed
        false);
  }

  ////////////////////////////////!!!Calculate master side Porosity!!!////////////////////////////
  // get J

  std::map<int, double>
      mJLin;  // master map for linearizations of the Deformation Gradient Determinant
  // get fluid pressure in GP
  double mgpfpres = 0.0;
  // porosity to be written into
  double mporosity;
  double mdphi_dp;  // linearization terms to be written into
  double mdphi_dJ;  // linearization terms to be written into

  if (masterporo)
  {
    const double mJ = DetDeformationGradient(mele, wgt, gpcoord, mJLin);

    for (int i = 0; i < ncol; ++i)
    {
      CoNode* mymrtrnode = dynamic_cast<CoNode*>(mnodes[i]);
      mgpfpres += mval[i] * mymrtrnode->CoPoroData().fpres()[0];
    }
    Teuchos::ParameterList mparams;  // empty parameter list;

    Teuchos::RCP<MAT::StructPoro> mstructmat =
        Teuchos::rcp_dynamic_cast<MAT::StructPoro>(mele.ParentElement()->Material(0));
    if (mstructmat == Teuchos::null)
      mstructmat = Teuchos::rcp_dynamic_cast<MAT::StructPoro>(mele.ParentElement()->Material(1));
    if (mstructmat == Teuchos::null) dserror("Cast to StructPoro failed!");

    mstructmat->ComputeSurfPorosity(mparams, mgpfpres, mJ,
        mele.FaceParentNumber(),  // may not work
        1,                        // finally check what to do here Todo:
        mporosity,
        &mdphi_dp,  // linearization of phi w.r.t. pressure
        &mdphi_dJ,  // linearization of phi w.r.t. defo-grad determinant
        NULL,       // dphi_dJdp not needed
        NULL,       // dphi_dJJ not needed
        NULL,       // dphi_dpp not needed
        false);
  }
  ////////////////////////////////!!!Calculate Porosity done!!!////////////////////////////
  // build normal coupling term at current GP
  for (int i = 0; i < Dim(); ++i)
  {
    if (slaveporo) ncoup[0] += (sgpsvel[i] - sgpfvel[i]) * gpn[i] * sporosity;
    if (masterporo) ncoup[0] -= (mgpsvel[i] - mgpfvel[i]) * gpn[i] * mporosity;
  }
  // **************************
  // add to node
  // **************************
  if (!quadratic)
  {
    for (int j = 0; j < nrow; ++j)
    {
      CONTACT::CoNode* mrtrnode = dynamic_cast<CONTACT::CoNode*>(snodes[j]);

      double prod = 0.0;
      // Petrov-Galerkin approach (dual LM for D/M but standard LM for normal coupling)
      if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin) prod = sval[j] * ncoup[0] * jac * wgt;
      // usual standard or dual LM approach
      else
        prod = lmval[j] * ncoup[0] * jac * wgt;

      // do not process slave side boundary nodes
      // (their row entries would be zero anyway!)
      if (mrtrnode->IsOnBound())
      {
        dserror("CHECKME: isonbound_ active; should not happen here");
        continue;
      }

      // if (cnode->Owner()!=Comm_.MyPID()) continue;

      // add current Gauss point's contribution to gseg
      mrtrnode->AddNcoupValue(prod);
    }
  }

  //    // CASE 4: Dual LM shape functions and quadratic interpolation
  //    else if ((ShapeFcn() == INPAR::MORTAR::shape_dual || ShapeFcn() ==
  //    INPAR::MORTAR::shape_petrovgalerkin) &&
  //        LagMultQuad() == INPAR::MORTAR::lagmult_quad)
  //    {
  //      for (int j=0;j<nrow;++j)
  //      {
  //        CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(snodes[j]);
  //
  //        double prod = 0.0;
  //        prod = lmval[j]*gap[0]*jac*wgt;
  //
  //        // add current Gauss point's contribution to gseg
  //        cnode->AddgValue(prod);
  //      }
  //    }

  // INVALID CASES
  else
  {
    dserror("ERROR: Invalid integration case for 3D quadratic normal coupling mortar!");
  }
  //}
  // **************************
  // Linearization w.r.t. displacements
  // **************************

  // add everything to dncoupgp
  // dncoupgp ... (v(struct)-v(fluid)) * delta n
  for (int k = 0; k < Dim(); ++k)
  {
    for (_CI p = dnmap_unit[k].begin(); p != dnmap_unit[k].end(); ++p)
    {
      double slaveside = 0.0;
      if (slaveporo)
        slaveside = sporosity * (sgpsvel[k] - sgpfvel[k]);  // h.Willmann slave terms preparation

      double masterside = 0.0;
      if (masterporo)
        masterside = mporosity * (mgpsvel[k] - mgpfvel[k]);  // h.Willmann master terms preparation

      dncoupgp[p->first] += (slaveside - masterside) * (p->second);  // linearization of n,
    }
  }

  double timefac =
      imortar_.get<double>("porotimefac");  // TODO: move in final version to other place ChrAg

  if (slaveporo)
  {
    for (int i = 0; i < nrow; ++i)
    {
      CoNode* mrtrnode = dynamic_cast<CoNode*>(snodes[i]);

      for (int k = 0; k < Dim(); ++k)
      {
        dncoupgp[mrtrnode->Dofs()[k]] +=
            sporosity * sval[i] * gpn[k] * timefac;  // linearisation of vs add for master side
        for (unsigned int d = 0; d < dsxigp.size(); ++d)
        {
          for (_CI p = dsxigp[d].begin(); p != dsxigp[d].end();
               ++p)  // linearization of shape function N
          {
            dncoupgp[p->first] -=
                sporosity * gpn[k] * sderiv(i, d) * mrtrnode->CoPoroData().fvel()[k] * (p->second);
            dncoupgp[p->first] +=
                sporosity * gpn[k] * sderiv(i, d) * mrtrnode->CoPoroData().svel()[k] * (p->second);
          }
        }
      }
    }
  }  // if(slaveporo)

  // h.Willmann master side is considered now
  // MASTER
  if (masterporo)
  {
    // lin master nodes
    for (int i = 0; i < ncol; ++i)
    {
      CoNode* mrtrnode = dynamic_cast<CoNode*>(mnodes[i]);

      for (int k = 0; k < Dim(); ++k)
      {
        dncoupgp[mrtrnode->Dofs()[k]] -= mporosity * mval[i] * gpn[k] * timefac;
        for (unsigned int d = 0; d < dmxigp.size(); ++d)
        {
          for (_CI p = dmxigp[d].begin(); p != dmxigp[d].end(); ++p)
          {
            dncoupgp[p->first] +=
                mporosity * gpn[k] * mderiv(i, d) * mrtrnode->CoPoroData().fvel()[k] * (p->second);
            dncoupgp[p->first] -=
                mporosity * gpn[k] * mderiv(i, d) * mrtrnode->CoPoroData().svel()[k] * (p->second);
          }
        }
      }
    }
  }

  // h.Willmann Write Deformation Gradient Determinant Linearization to map
  typedef std::map<int, double>::iterator I;

  if (slaveporo)
  {
    for (I p = sJLin.begin(); p != sJLin.end(); ++p)
    {
      for (int i = 0; i < Dim(); ++i)
      {
        dncoupgp[p->first] += (sgpsvel[i] - sgpfvel[i]) * sdphi_dJ * gpn[i] * p->second;
      }
    }
  }
  if (masterporo)
  {
    for (I p = mJLin.begin(); p != mJLin.end(); ++p)
    {
      for (int i = 0; i < Dim(); ++i)
      {
        dncoupgp[p->first] -= (mgpsvel[i] - mgpfvel[i]) * mdphi_dJ * gpn[i] * p->second;
      }
    }
  }

  // **************************
  // Linearization w.r.t. fluid velocities
  // **************************

  if (slaveporo)
  {
    for (int z = 0; z < nrow; ++z)
    {
      CoNode* mrtrnode = dynamic_cast<CoNode*>(snodes[z]);

      for (int k = 0; k < Dim(); ++k)
      {
        dvelncoupgp[mrtrnode->Dofs()[k]] -= sporosity * sval[z] * gpn[k];
        // Because Slave Discretisation should be the poro dis!!!  --- add master for two sided poro
        // contact!!!

        // h.Willmann linearization w.r.t. fluid pressure:
        dpresncoupgp[mrtrnode->Dofs()[0]] +=
            sdphi_dp * sval[z] * (sgpsvel[k] - sgpfvel[k]) * gpn[k];
        //      std::cout<<mrtrnode->Dofs()[0]<<"mrtrnode->Dofs()[0]"<<std::endl;
      }
    }
  }
  // h.Willmann terms added for master side
  if (masterporo)
  {
    for (int z = 0; z < ncol; ++z)
    {
      CoNode* mrtrnode = dynamic_cast<CoNode*>(mnodes[z]);

      for (int k = 0; k < Dim(); ++k)
      {
        dvelncoupgp[mrtrnode->Dofs()[k]] += mporosity * mval[z] * gpn[k];

        // h.Willmann linearization w.r.t. fluid pressure:
        dpresncoupgp[mrtrnode->Dofs()[0]] -=
            mdphi_dp * mval[z] * (mgpsvel[k] - mgpfvel[k]) * gpn[k];
        //       std::cout<<mrtrnode->Dofs()[0]<<"mrtrnode->Dofs()[0]"<<std::endl;
      }
    }
  }
  return;
}

/*-----------------------------------------------------------------------------*
 |  Do lin. entries for weighted normal coupling condition at GP     ager 06/14|
 |  modified by h.Willmann 2015                                                |
 *----------------------------------------------------------------------------*/
void inline CONTACT::CoIntegrator::GP_NCOUP_LIN(int& iter, MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
    LINALG::SerialDenseVector& lmval, LINALG::SerialDenseMatrix& sderiv,
    LINALG::SerialDenseMatrix& lmderiv, double& ncoup, double* gpn, double& jac, double& wgt,
    const std::map<int, double>& dncoupgp, const std::map<int, double>& dvelncoupgp,
    const std::map<int, double>& dpresncoupgp, const GEN::pairedvector<int, double>& jacintcellmap,
    const std::vector<GEN::pairedvector<int, double>>& dsxigp,
    const std::vector<GEN::pairedvector<int, double>>& dmxigp,
    const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap)
{
  int nrow = sele.NumNode();
  // only sele.NumNode() is needed because the integration of the interface constraint is
  // evaluated on the slave side interface surface

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator _CI;
  typedef std::map<int, double>::const_iterator CI;

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  // DRT::Node** mnodes = mele.Nodes();

  CONTACT::CoNode* mymrtrnode = dynamic_cast<CONTACT::CoNode*>(snodes[iter]);
  if (!mymrtrnode) dserror("ERROR: CONTACT::CoIntegrator::GP_NCOUP_LIN: mymrtnode: Null pointer!");

  double fac = 0.0;

  // get the corresponding map as a reference
  std::map<int, double>& dgmap = mymrtrnode->CoPoroData().GetDerivnCoup();
  // switch if Petrov-Galerkin approach for LM is applied
  // for meshtying shape_dual is the only allowed case as long as porofluid and skeleton
  // meshtying conditions share the same input lines
  //(for contact they still share the input lines but shape_petrovgalerkin is allowed for
  // skeleton/structural contact) shape_petrovgalerkin
  if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
  {
    // (1) Lin(Phi) - does not exist here for Petrov-Galerkin approach

    // (2) Lin(N) - slave GP coordinates
    for (unsigned int k = 0; k < dsxigp.size(); ++k)
    {
      fac = wgt * sderiv(iter, k) * ncoup * jac;
      for (_CI p = dsxigp[k].begin(); p != dsxigp[k].end(); ++p)
        dgmap[p->first] += fac * (p->second);
    }
    // (3) Lin(g) - normal coupling condition ncoup function
    fac = wgt * sval[iter] * jac;
    for (CI p = dncoupgp.begin(); p != dncoupgp.end(); ++p) dgmap[p->first] += fac * (p->second);

    // (4) Lin(dsxideta) - intcell GP Jacobian
    fac = wgt * sval[iter] * ncoup;
    for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
      dgmap[p->first] += fac * (p->second);
  }

  // the usual standard or dual LM approach
  else
  {
    // (1) Lin(Phi) - dual shape functions
    if (dualmap.size() > 0)
      for (int m = 0; m < nrow; ++m)
      {
        fac = wgt * sval[m] * ncoup * jac;
        for (GEN::pairedvector<int, Epetra_SerialDenseMatrix>::const_iterator p = dualmap.begin();
             p != dualmap.end(); ++p)
        {
          dgmap[p->first] += fac * (p->second)(iter, m);
        }
      }

    // (2) Lin(Phi) - slave GP coordinates
    for (unsigned int k = 0; k < dsxigp.size(); ++k)
    {
      fac = wgt * lmderiv(iter, k) * ncoup * jac;
      for (_CI p = dsxigp[k].begin(); p != dsxigp[k].end(); ++p)
      {
        dgmap[p->first] += fac * (p->second);
      }
    }

    // (3) Lin(g) - normal coupling condition ncoup function
    fac = wgt * lmval[iter] * jac;
    for (CI p = dncoupgp.begin(); p != dncoupgp.end(); ++p)
    {
      dgmap[p->first] += fac * (p->second);
    }

    // (4) Lin(dsxideta) - intcell GP Jacobian
    fac = wgt * lmval[iter] * ncoup;
    for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
    {
      dgmap[p->first] += fac * (p->second);
    }
  }

  // velocity linearisation of the ncoupling condition!
  // get the corresponding map as a reference
  std::map<int, double>& dvelncoupmap = mymrtrnode->CoPoroData().GetVelDerivnCoup();
  // switch if Petrov-Galerkin approach for LM is applied
  if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
  {
    // (3) Lin(g) - normal coupling condition function
    fac = wgt * sval[iter] * jac;
    for (CI p = dvelncoupgp.begin(); p != dvelncoupgp.end(); ++p)
      dvelncoupmap[p->first] += fac * (p->second);
  }
  // the usual standard or dual LM approach
  else
  {
    // (3) Lin(g) - normal coupling condition function
    fac = wgt * lmval[iter] * jac;
    for (CI p = dvelncoupgp.begin(); p != dvelncoupgp.end(); ++p)
    {
      dvelncoupmap[p->first] += fac * (p->second);
    }
  }

  // h.Willmann
  // pressure linearisation of the ncoupling condition!
  // get the corresponding map as a reference
  std::map<int, double>& dpresncoupmap = mymrtrnode->CoPoroData().GetPresDerivnCoup();
  // switch if Petrov-Galerkin approach for LM is applied
  if (ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
  {
    // (3) Lin(g) - normal coupling condition function
    fac = wgt * sval[iter] * jac;
    for (CI p = dpresncoupgp.begin(); p != dpresncoupgp.end(); ++p)
      dpresncoupmap[p->first] += fac * (p->second);
  }
  // the usual standard or dual LM approach
  else
  {
    // (3) Lin(g) - normal coupling condition function
    fac = wgt * lmval[iter] * jac;
    for (CI p = dpresncoupgp.begin(); p != dpresncoupgp.end(); ++p)
    {
      dpresncoupmap[p->first] += fac * (p->second);
    }
  }
  return;
}

/*-----------------------------------------------------------------------------*
 |  Calculate Determinate of the Deformation Gradient at GP          ager 10/14|
 *----------------------------------------------------------------------------*/
double CONTACT::CoIntegrator::DetDeformationGradient(
    MORTAR::MortarElement& sele, double& wgt, double* gpcoord, std::map<int, double>& JLin)
{
  double J;
  DRT::Element::DiscretizationType distype = sele.ParentElement()->Shape();
  switch (distype)
  {
    case DRT::Element::hex8:
      J = TDetDeformationGradient<DRT::Element::hex8, 3>(sele, wgt, gpcoord, JLin);
      break;
    case DRT::Element::quad4:
      J = TDetDeformationGradient<DRT::Element::quad4, 2>(sele, wgt, gpcoord, JLin);
      break;
    default:
      dserror("DetDeformationGradient3D: Parent Element Type not templated yet, just add it here!");
      J = 0.0;  // To avoid warning!
      break;
  }

  return J;
}
/*-------------------------------------------------------------------------------*
 |  Templated Calculate Determinant of the Deformation Gradient at GP  ager 10/14|
 |  modified by h.Willmann 2015                                                  |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType parentdistype, int dim>
double CONTACT::CoIntegrator::TDetDeformationGradient(
    MORTAR::MortarElement& sele, double& wgt, double* gpcoord, std::map<int, double>& JLin)
{
  //! nen_: number of element nodes (T. Hughes: The Finite Element Method)
  static const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<parentdistype>::numNodePerElement;

  DRT::UTILS::CollectedGaussPoints intpoints =
      DRT::UTILS::CollectedGaussPoints(1);  // reserve just for 1 entry ...
  intpoints.Append(gpcoord[0], gpcoord[1], 0.0, wgt);

  // get coordinates of gauss point w.r.t. local parent coordinate system
  LINALG::SerialDenseMatrix pqxg(1, dim);
  LINALG::Matrix<dim, dim> derivtrafo(true);

  DRT::UTILS::BoundaryGPToParentGP<dim>(pqxg, derivtrafo, intpoints, sele.ParentElement()->Shape(),
      sele.Shape(), sele.FaceParentNumber());

  LINALG::Matrix<dim, 1> pxsi(true);

  // coordinates of the current integration point in parent coordinate system
  for (int idim = 0; idim < dim; idim++)
  {
    pxsi(idim) = pqxg(0, idim);
  }

  LINALG::Matrix<dim, numnodes> pderiv_loc(
      true);  // derivatives of parent element shape functions in parent element coordinate system

  // evaluate derivatives of parent element shape functions at current integration point in parent
  // coordinate system
  DRT::UTILS::shape_function_deriv1<parentdistype>(pxsi, pderiv_loc);
  //
  // get Jacobian matrix and determinant w.r.t. spatial configuration
  //
  // |J| = det(xjm) * det(Jmat^-1) = det(xjm) * 1/det(Jmat)
  //
  //    _                     _
  //   |  x_1,1  x_2,1  x_3,1  |           d x_i
  //   |  x_1,2  x_2,2  x_3,2  | = xjm  = --------
  //   |_ x_1,3  x_2,3  x_3,3 _|           d s_j
  //    _                     _
  //   |  X_1,1  X_2,1  X_3,1  |           d X_i
  //   |  X_1,2  X_2,2  X_3,2  | = Jmat = --------
  //   |_ X_1,3  X_2,3  X_3,3 _|           d s_j
  //
  LINALG::Matrix<dim, dim> xjm;
  LINALG::Matrix<dim, dim> Jmat;

  LINALG::Matrix<dim, numnodes> xrefe(true);  // material coord. of parent element
  LINALG::Matrix<dim, numnodes> xcurr(true);  // current  coord. of parent element

  // update element geometry of parent element
  {
    DRT::Node** nodes = sele.ParentElement()->Nodes();
    for (int inode = 0; inode < numnodes; ++inode)
    {
      for (unsigned int idof = 0; idof < dim; ++idof)
      {
        const double* x = nodes[inode]->X();
        xrefe(idof, inode) = x[idof];
        xcurr(idof, inode) = xrefe(idof, inode) + sele.MoData().ParentDisp()[inode * dim + idof];
        //          std::cout<< sele.MoData().ParentDisp()[inode*dim+idof] <<",
        //          sele.MoData().ParentDisp()[inode*dim+idof]"<<std::endl;
      }
    }
  }
  // std::cout<<xcurr<<"xcurr, "<<xrefe<<"xrefe"<<std::endl;
  xjm.MultiplyNT(pderiv_loc, xcurr);
  Jmat.MultiplyNT(pderiv_loc, xrefe);
  double det = xjm.Determinant();
  double detJ = Jmat.Determinant();
  const double J = det / detJ;

  // Linearisation of Jacobian (atm missing!)
  // D J[d] = J div(d) = J div(N_i d_i) = J dN_i/dx_j d_ij = J dN_i/dzeta_k dzeta_k/dx_j d_ij

  //  LINALG::Matrix<dim,numnodes>  auxJLin (true);   // matrix to be filled with linearization
  //  scalars
  xjm.Invert();
  {
    //    DRT::Node** nodes = sele.ParentElement()->Nodes();
    for (int inode = 0; inode < numnodes; ++inode)
    {
      for (unsigned int i = 0; i < dim; ++i)
      {
        for (unsigned int j = 0; j < dim; ++j)
        {
          JLin[sele.MoData().ParentDof().at(inode * dim + i)] +=
              J * pderiv_loc(j, inode) * xjm(i, j);
          //          auxJLin(idof,inode)+=J*pderiv_loc(j,inode)*xjm(i,j); //matrix that could be
          //          filled with linearizations
        }
      }
    }
  }
  return J;
}
