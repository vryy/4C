/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of Mortar matrices on the overlap
       of two Mortar::Elements in 1D and 2D (derived version for contact)

\level 2


*/
/*---------------------------------------------------------------------*/
#include "4C_contact_integrator.hpp"

#include "4C_contact_defines.hpp"
#include "4C_contact_element.hpp"
#include "4C_contact_friction_node.hpp"
#include "4C_contact_node.hpp"
#include "4C_contact_paramsinterface.hpp"
#include "4C_discretization_fem_general_utils_boundary_integration.hpp"
#include "4C_discretization_fem_general_utils_integration.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_inpar_wear.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_mat_structporo.hpp"
#include "4C_mortar_calc_utils.hpp"
#include "4C_mortar_coupling3d_classes.hpp"
#include "4C_mortar_defines.hpp"
#include "4C_mortar_projector.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 10/13|
 *----------------------------------------------------------------------*/
CONTACT::Integrator::Integrator(
    Teuchos::ParameterList& params, Core::FE::CellType eletype, const Epetra_Comm& comm)
    : imortar_(params),
      Comm_(comm),
      dim_(imortar_.get<int>("DIMENSION")),
      shapefcn_(Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(imortar_, "LM_SHAPEFCN")),
      lagmultquad_(Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(imortar_, "LM_QUAD")),
      gpslip_(Core::UTILS::IntegralValue<int>(imortar_, "GP_SLIP_INCR")),
      algo_(Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(imortar_, "ALGORITHM")),
      stype_(Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(imortar_, "STRATEGY")),
      cppnormal_(Core::UTILS::IntegralValue<int>(imortar_, "CPP_NORMALS")),
      wearlaw_(Core::UTILS::IntegralValue<Inpar::Wear::WearLaw>(imortar_, "WEARLAW")),
      wearimpl_(false),
      wearside_(Inpar::Wear::wear_slave),
      weartype_(Inpar::Wear::wear_intstate),
      wearshapefcn_(Inpar::Wear::wear_shape_standard),
      sswear_(Core::UTILS::IntegralValue<int>(imortar_, "SSWEAR")),
      wearcoeff_(-1.0),
      wearcoeffm_(-1.0),
      ssslip_(imortar_.get<double>("SSSLIP")),
      nonsmooth_(false),
      nonsmoothselfcontactsurface_(
          Core::UTILS::IntegralValue<int>(imortar_, "NONSMOOTH_CONTACT_SURFACE")),
      integrationtype_(Core::UTILS::IntegralValue<Inpar::Mortar::IntType>(imortar_, "INTTYPE"))
{
  // init gp
  initialize_gp(eletype);

  // wear specific
  if (wearlaw_ != Inpar::Wear::wear_none)
  {
    // set wear contact status
    Inpar::Wear::WearTimInt wtimint =
        Core::UTILS::IntegralValue<Inpar::Wear::WearTimInt>(params, "WEARTIMINT");
    if (wtimint == Inpar::Wear::wear_impl) wearimpl_ = true;

    // wear surface
    wearside_ = Core::UTILS::IntegralValue<Inpar::Wear::WearSide>(imortar_, "WEAR_SIDE");

    // wear algorithm
    weartype_ = Core::UTILS::IntegralValue<Inpar::Wear::WearType>(imortar_, "WEARTYPE");

    // wear shape function
    wearshapefcn_ = Core::UTILS::IntegralValue<Inpar::Wear::WearShape>(imortar_, "WEAR_SHAPEFCN");

    // wear coefficient
    wearcoeff_ = imortar_.get<double>("WEARCOEFF");

    // wear coefficient
    wearcoeffm_ = imortar_.get<double>("WEARCOEFF_MASTER");
  }

  nonsmooth_ = Core::UTILS::IntegralValue<int>(imortar_, "NONSMOOTH_GEOMETRIES");

  return;
}

/*----------------------------------------------------------------------*
 |  check for boundary elements                              farah 02/14|
 *----------------------------------------------------------------------*/
bool CONTACT::Integrator::BoundarySegmCheck2D(
    Mortar::Element& sele, std::vector<Mortar::Element*> meles)
{
  double sxi_test[2] = {0.0, 0.0};
  bool proj_test = false;
  bool boundary_ele = false;

  double glob_test[3] = {0.0, 0.0, 0.0};

  Core::Nodes::Node** mynodes_test = sele.Nodes();
  if (!mynodes_test) FOUR_C_THROW("has_proj_status: Null pointer!");

  if (sele.Shape() == Core::FE::CellType::line2 || sele.Shape() == Core::FE::CellType::nurbs2)
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
        Mortar::Projector::Impl(sele, *meles[bs_test])
            ->ProjectGaussPoint2D(sele, sxi_test, *meles[bs_test], mxi_test);

        if ((mxi_test[0] >= -1.0) && (mxi_test[0] <= 1.0))
        {
          // get hasproj
          sele.LocalToGlobal(sxi_test, glob_test, 0);
          for (int ii = 0; ii < sele.num_node(); ++ii)
          {
            Mortar::Node* mycnode_test = dynamic_cast<Mortar::Node*>(mynodes_test[ii]);
            if (!mycnode_test) FOUR_C_THROW("has_proj_status: Null pointer!");

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
  else if (sele.Shape() == Core::FE::CellType::line3 || sele.Shape() == Core::FE::CellType::nurbs3)
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
        Mortar::Projector::Impl(sele, *meles[bs_test])
            ->ProjectGaussPoint2D(sele, sxi_test, *meles[bs_test], mxi_test);

        if ((mxi_test[0] >= -1.0) && (mxi_test[0] <= 1.0))
        {
          // get hasproj
          sele.LocalToGlobal(sxi_test, glob_test, 0);
          for (int ii = 0; ii < sele.num_node(); ++ii)
          {
            Mortar::Node* mycnode_test = dynamic_cast<Mortar::Node*>(mynodes_test[ii]);
            if (!mycnode_test) FOUR_C_THROW("has_proj_status: Null pointer!");

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
    FOUR_C_THROW("No valid element type for slave discretization!");
  }

  return boundary_ele;
}


/*----------------------------------------------------------------------*
 |  Initialize gauss points                                   popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::initialize_gp(Core::FE::CellType eletype)
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
  Inpar::Mortar::IntType integrationtype =
      Core::UTILS::IntegralValue<Inpar::Mortar::IntType>(imortar_, "INTTYPE");

  //**********************************************************************
  // choose Gauss rule according to (a) element type (b) input parameter
  //**********************************************************************
  switch (eletype)
  {
    case Core::FE::CellType::line2:
    case Core::FE::CellType::line3:
    case Core::FE::CellType::nurbs2:
    case Core::FE::CellType::nurbs3:
    {
      // set default value for segment-based version first
      Core::FE::GaussRule1D mygaussrule = Core::FE::GaussRule1D::line_5point;

      // GP switch if element-based version and non-zero value provided by user
      if (integrationtype == Inpar::Mortar::inttype_elements ||
          integrationtype == Inpar::Mortar::inttype_elements_BS)
      {
        if (numgp > 0)
        {
          switch (numgp)
          {
            case 1:
            {
              FOUR_C_THROW("Our experience says that 1 GP per slave element is not enough.");
              break;
            }
            case 2:
            {
              mygaussrule = Core::FE::GaussRule1D::line_2point;
              break;
            }
            case 3:
            {
              mygaussrule = Core::FE::GaussRule1D::line_3point;
              break;
            }
            case 4:
            {
              mygaussrule = Core::FE::GaussRule1D::line_4point;
              break;
            }
            case 5:
            {
              mygaussrule = Core::FE::GaussRule1D::line_5point;
              break;
            }
            case 6:
            {
              mygaussrule = Core::FE::GaussRule1D::line_6point;
              break;
            }
            case 7:
            {
              mygaussrule = Core::FE::GaussRule1D::line_7point;
              break;
            }
            case 8:
            {
              mygaussrule = Core::FE::GaussRule1D::line_8point;
              break;
            }
            case 9:
            {
              mygaussrule = Core::FE::GaussRule1D::line_9point;
              break;
            }
            case 10:
            {
              mygaussrule = Core::FE::GaussRule1D::line_10point;
              break;
            }
            case 16:
            {
              mygaussrule = Core::FE::GaussRule1D::line_16point;
              break;
            }
            case 20:
            {
              mygaussrule = Core::FE::GaussRule1D::line_20point;
              break;
            }
            case 32:
            {
              mygaussrule = Core::FE::GaussRule1D::line_32point;
              break;
            }
            case 50:
            {
              mygaussrule = Core::FE::GaussRule1D::line_50point;
              break;
            }
            default:
            {
              FOUR_C_THROW("Requested GP-Number is not implemented!");
              break;
            }
          }
        }
      }

      const Core::FE::IntegrationPoints1D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(nGP(), 2);
      weights_.resize(nGP());
      for (int i = 0; i < nGP(); ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = 0.0;
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case Core::FE::CellType::tri3:
    case Core::FE::CellType::tri6:
    {
      // set default value for segment-based version first
      Core::FE::GaussRule2D mygaussrule = Core::FE::GaussRule2D::tri_7point;
      if (integrationtype == Inpar::Mortar::inttype_segments)
      {
        if (numgp > 0) switch (numgp)
          {
            case 1:
              mygaussrule = Core::FE::GaussRule2D::tri_1point;
              break;
            case 3:
              mygaussrule = Core::FE::GaussRule2D::tri_3point;
              break;
            case 7:
              mygaussrule = Core::FE::GaussRule2D::tri_7point;
              break;
            case 16:
              mygaussrule = Core::FE::GaussRule2D::tri_16point;
              break;
            case 37:
              mygaussrule = Core::FE::GaussRule2D::tri_37point;
              break;
            default:
              FOUR_C_THROW("unknown tri gauss rule");
              break;
          }
      }

      // GP switch if element-based version and non-zero value provided by user
      if (integrationtype == Inpar::Mortar::inttype_elements ||
          integrationtype == Inpar::Mortar::inttype_elements_BS)
      {
        if (numgp > 0)
        {
          switch (numgp)
          {
            case 1:
            {
              mygaussrule = Core::FE::GaussRule2D::tri_3point;
              break;
            }
            case 2:
            {
              mygaussrule = Core::FE::GaussRule2D::tri_6point;
              break;
            }
            case 3:
            {
              mygaussrule = Core::FE::GaussRule2D::tri_7point;
              break;
            }
            case 4:
            {
              mygaussrule = Core::FE::GaussRule2D::tri_12point;
              break;
            }
            case 5:
            {
              mygaussrule = Core::FE::GaussRule2D::tri_12point;
              break;
            }
            case 6:
            {
              mygaussrule = Core::FE::GaussRule2D::tri_37point;
              break;
            }
            case 7:
            {
              mygaussrule = Core::FE::GaussRule2D::tri_37point;
              break;
            }
            case 8:
            {
              mygaussrule = Core::FE::GaussRule2D::tri_64point;
              break;
            }
            case 9:
            {
              mygaussrule = Core::FE::GaussRule2D::tri_64point;
              break;
            }
            case 10:
            {
              mygaussrule = Core::FE::GaussRule2D::tri_64point;
              break;
            }
            case 20:
            {
              mygaussrule = Core::FE::GaussRule2D::tri_64point;
              break;
            }
            default:
            {
              FOUR_C_THROW("Requested GP-Number is not implemented!");
              break;
            }
          }
        }
      }

      const Core::FE::IntegrationPoints2D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(nGP(), 2);
      weights_.resize(nGP());
      for (int i = 0; i < nGP(); ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case Core::FE::CellType::quad4:
    case Core::FE::CellType::quad8:
    case Core::FE::CellType::quad9:
    case Core::FE::CellType::nurbs4:
    case Core::FE::CellType::nurbs8:
    case Core::FE::CellType::nurbs9:
    {
      // set default value for segment-based version first
      Core::FE::GaussRule2D mygaussrule = Core::FE::GaussRule2D::quad_9point;

      // GP switch if element-based version and non-zero value provided by user
      if (integrationtype == Inpar::Mortar::inttype_elements ||
          integrationtype == Inpar::Mortar::inttype_elements_BS)
      {
        if (numgp > 0)
        {
          switch (numgp)
          {
            case 1:
            {
              mygaussrule = Core::FE::GaussRule2D::quad_1point;
              break;
            }
            case 2:
            {
              mygaussrule = Core::FE::GaussRule2D::quad_4point;
              break;
            }
            case 3:
            {
              mygaussrule = Core::FE::GaussRule2D::quad_9point;
              break;
            }
            case 4:
            {
              mygaussrule = Core::FE::GaussRule2D::quad_16point;
              break;
            }
            case 5:
            {
              mygaussrule = Core::FE::GaussRule2D::quad_25point;
              break;
            }
            case 6:
            {
              mygaussrule = Core::FE::GaussRule2D::quad_36point;
              break;
            }
            case 7:
            {
              mygaussrule = Core::FE::GaussRule2D::quad_49point;
              break;
            }
            case 8:
            {
              mygaussrule = Core::FE::GaussRule2D::quad_64point;
              break;
            }
            case 9:
            {
              mygaussrule = Core::FE::GaussRule2D::quad_81point;
              break;
            }
            case 10:
            {
              mygaussrule = Core::FE::GaussRule2D::quad_100point;
              break;
            }
            case 16:
            {
              mygaussrule = Core::FE::GaussRule2D::quad_256point;
              break;
            }
            case 20:
            {
              mygaussrule = Core::FE::GaussRule2D::quad_400point;
              break;
            }
            case 32:
            {
              mygaussrule = Core::FE::GaussRule2D::quad_1024point;
              break;
            }
            default:
            {
              FOUR_C_THROW("Requested GP-Number is not implemented!");
              break;
            }
          }
        }
      }

      const Core::FE::IntegrationPoints2D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(nGP(), 2);
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
      FOUR_C_THROW("Mortar::Integrator: This contact element type is not implemented!");
      break;
    }
  }  // switch(eletype)

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::integrate_deriv_segment2_d(Mortar::Element& sele, double& sxia,
    double& sxib, Mortar::Element& mele, double& mxia, double& mxib, const Epetra_Comm& comm,
    const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr)
{
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr = Teuchos::null;
  if (not mparams_ptr.is_null())
  {
    cparams_ptr = Teuchos::rcp_dynamic_cast<CONTACT::ParamsInterface>(mparams_ptr, true);
  }
  integrate_deriv_segment2_d(sele, sxia, sxib, mele, mxia, mxib, comm, cparams_ptr);
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
void CONTACT::Integrator::integrate_deriv_segment2_d(Mortar::Element& sele, double& sxia,
    double& sxib, Mortar::Element& mele, double& mxia, double& mxib, const Epetra_Comm& comm,
    const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  // skip this segment, if too small
  if (sxib - sxia < 4. * MORTARINTLIM) return;

  // *********************************************************************
  // Check integrator input for non-reasonable quantities
  // *********************************************************************

  // explicitly defined shape function type needed
  if (shape_fcn() == Inpar::Mortar::shape_undefined)
    FOUR_C_THROW("integrate_deriv_segment2_d called without specific shape function defined!");

  // Petrov-Galerkin approach for LM not yet implemented for quadratic FE
  if (sele.Shape() == Core::FE::CellType::line3 &&
      shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
    FOUR_C_THROW("Petrov-Galerkin approach not yet implemented for 2-D quadratic FE interpolation");

  // check for problem dimension
  if (Dim() != 2) FOUR_C_THROW("2D integration method called for non-2D problem");

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    FOUR_C_THROW("integrate_deriv_segment2_d called on a wrong type of Mortar::Element pair!");
  if ((sxia < -1.0) || (sxib > 1.0))
    FOUR_C_THROW("integrate_deriv_segment2_d called with infeasible slave limits!");
  if ((mxia < -1.0) || (mxib > 1.0))
    FOUR_C_THROW("integrate_deriv_segment2_d called with infeasible master limits!");

  // *********************************************************************
  // Prepare integration
  // *********************************************************************
  // number of nodes (slave, master)
  int nrow = sele.num_node();
  int ncol = mele.num_node();
  int ndof = Dim();

  // get slave element nodes themselves
  Core::Nodes::Node** mynodes = sele.Nodes();
  if (!mynodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  for (int k = 0; k < nrow; ++k)
  {
    Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(mynodes[k]);
    if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_segment2_d: Null pointer!");
  }

  // decide whether linear LM are used for quadratic FE here
  // and whether displacement shape fct. modification has to be considered or not
  // this is the case for dual linear Lagrange multipliers on line3 elements
  bool linlm = false;
  bool dualquad = false;
  if (lag_mult_quad() == Inpar::Mortar::lagmult_lin && sele.Shape() == Core::FE::CellType::line3)
  {
    linlm = true;
    if (shapefcn_ == Inpar::Mortar::shape_dual) dualquad = true;
  }

  // prepare directional derivative of dual shape functions
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix> dualmap(
      2 * nrow, 0, Core::LinAlg::SerialDenseMatrix(nrow, nrow));
  if ((shape_fcn() == Inpar::Mortar::shape_dual ||
          shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) &&
      (sele.Shape() == Core::FE::CellType::line3 || sele.Shape() == Core::FE::CellType::nurbs3 ||
          sele.MoData().DerivDualShape() != Teuchos::null))
    sele.DerivShapeDual(dualmap);

  // create empty vectors for shape fct. evaluation
  Core::LinAlg::SerialDenseVector sval(nrow);
  Core::LinAlg::SerialDenseMatrix sderiv(nrow, 1);
  Core::LinAlg::SerialDenseVector mval(ncol);
  Core::LinAlg::SerialDenseMatrix mderiv(ncol, 1);
  Core::LinAlg::SerialDenseVector lmval(nrow);
  Core::LinAlg::SerialDenseMatrix lmderiv(nrow, 1);

  Core::LinAlg::SerialDenseVector m2val(ncol);
  Core::LinAlg::SerialDenseMatrix m2deriv(ncol, 1);
  Core::LinAlg::SerialDenseVector lm2val(ncol);
  Core::LinAlg::SerialDenseMatrix lm2deriv(ncol, 1);

  // create empty vectors for shape fct. evaluation
  Core::LinAlg::SerialDenseMatrix ssecderiv(nrow, 1);

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  int linsize = 0;
  for (int i = 0; i < nrow; ++i)
  {
    Node* cnode = dynamic_cast<Node*>(mynodes[i]);
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
      FOUR_C_THROW("First outer node is neither slave nor master node");
    if (sxib != 1.0 && mxia != -1.0)
      FOUR_C_THROW("Second outer node is neither slave nor master node");
  }
  else
  {
    if (sxia != -1. && mxia != -1.)
      FOUR_C_THROW("First outer node is neither slave nor master node");
    if (sxib != 1. && mxib != 1.)
      FOUR_C_THROW("Second outer node is neither slave nor master node");
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
  std::vector<Core::Gen::Pairedvector<int, double>> ximaps(4, linsize + ndof * ncol);
  DerivXiAB2D(sele, sxia, sxib, mele, mxia, mxib, ximaps, startslave, endslave, linsize);

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < nGP(); ++gp)
  {
    // coordinates and weight
    std::array<double, 2> eta = {Coordinate(gp, 0), 0.0};
    double wgt = Weight(gp);

    // coordinate transformation sxi->eta (slave Mortar::Element->Overlap)
    double sxi[2] = {0.0, 0.0};
    sxi[0] = 0.5 * (1.0 - eta[0]) * sxia + 0.5 * (1.0 + eta[0]) * sxib;

    // project Gauss point onto master element
    double mxi[2] = {0.0, 0.0};
    Mortar::Projector::Impl(sele, mele)->ProjectGaussPoint2D(sele, sxi, mele, mxi);

    // check GP projection
    if ((mxi[0] < mxia - 1e-4) || (mxi[0] > mxib + 1e-4))
    {
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Gauss point: " << sxi[0] << " " << sxi[1] << std::endl;
      std::cout << "Projection: " << mxi[0] << " " << mxi[1] << std::endl;
      std::cout << "Projection bounds: " << mxia << " " << mxib << std::endl;
      // FOUR_C_THROW("IntegrateAndDerivSegment: Gauss point projection failed! mxi=%d",mxi[0]);
    }

    // evaluate Lagrange multiplier shape functions (on slave element)
    if (lag_mult_quad() == Inpar::Mortar::lagmult_const)
      sele.evaluate_shape_lag_mult_const(shape_fcn(), sxi, lmval, lmderiv, nrow);
    else if (linlm)
      sele.evaluate_shape_lag_mult_lin(shape_fcn(), sxi, lmval, lmderiv, nrow);
    else
      sele.evaluate_shape_lag_mult(shape_fcn(), sxi, lmval, lmderiv, nrow);

    // evaluate trace space shape functions (on both elements)
    sele.evaluate_shape(sxi, sval, sderiv, nrow, dualquad);
    mele.evaluate_shape(mxi, mval, mderiv, ncol, false);

    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5 * sxia + 0.5 * sxib;

    // evaluate linearizations *******************************************
    // evaluate 2nd deriv of trace space shape functions (on slave element)
    sele.evaluate2nd_deriv_shape(sxi, ssecderiv, nrow);

    // evaluate the derivative dxdsxidsxi = Jac,xi
    double djacdxi[2] = {0.0, 0.0};
    dynamic_cast<CONTACT::Element&>(sele).DJacDXi(djacdxi, sxi, ssecderiv);
    double dxdsxidsxi = djacdxi[0];  // only 2D here

    // evalute the GP slave coordinate derivatives
    std::vector<Core::Gen::Pairedvector<int, double>> dsxigp(
        1, Core::Gen::Pairedvector<int, double>(linsize + ndof * ncol));
    for (_CI p = ximaps[0].begin(); p != ximaps[0].end(); ++p)
      dsxigp[0][p->first] += 0.5 * (1.0 - eta[0]) * (p->second);
    for (_CI p = ximaps[1].begin(); p != ximaps[1].end(); ++p)
      dsxigp[0][p->first] += 0.5 * (1.0 + eta[0]) * (p->second);

    // evalute the GP master coordinate derivatives
    std::vector<Core::Gen::Pairedvector<int, double>> dmxigp(
        1, Core::Gen::Pairedvector<int, double>(linsize + ndof * ncol));
    deriv_xi_g_p2_d(sele, mele, sxi[0], mxi[0], dsxigp[0], dmxigp[0], linsize);

    // evaluate the Jacobian derivative
    Core::Gen::Pairedvector<int, double> derivjacSele(nrow * ndof);
    sele.DerivJacobian(sxi, derivjacSele);

    double jac = dsxideta * dxdsxi;
    Core::Gen::Pairedvector<int, double> derivjac(nrow * ndof + linsize);
    for (Core::Gen::Pairedvector<int, double>::const_iterator p = derivjacSele.begin();
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
    double gap = 0.0;            // gap
    Core::Gen::Pairedvector<int, double> dgapgp(
        linsize + ndof * ncol);  // gap lin without weighting and jac
    std::vector<Core::Gen::Pairedvector<int, double>> dnmap_unit(
        2, (linsize + ndof * ncol));  // deriv of x and y comp. of gpn (unit)

    //**********************************************************************
    // evaluate at GP and lin char. quantities
    //**********************************************************************
    if (dualquad && linlm)
    {
      // declare and compute shape functions as well as derivatives for computation of gap
      // -> standard quadratic functions instead of only linear functions
      Core::LinAlg::SerialDenseVector svalgap(nrow);
      Core::LinAlg::SerialDenseMatrix sderivgap(nrow, 1);
      sele.evaluate_shape(sxi, svalgap, sderivgap, nrow);

      gap_2_d(sele, mele, svalgap, mval, sderivgap, mderiv, &gap, gpn, dsxigp, dmxigp, dgapgp,
          dnmap_unit);
    }
    else
      gap_2_d(
          sele, mele, sval, mval, sderiv, mderiv, &gap, gpn, dsxigp, dmxigp, dgapgp, dnmap_unit);

    integrate_gp_2_d(sele, mele, sval, lmval, mval, sderiv, mderiv, lmderiv, dualmap, wgt, jac,
        derivjac, gpn, dnmap_unit, gap, dgapgp, sxi, mxi, dsxigp, dmxigp);

  }  // gp-loop

  return;
}

/*----------------------------------------------------------------------*
 |  check for boundary elements                              farah 07/14|
 *----------------------------------------------------------------------*/
bool CONTACT::Integrator::BoundarySegmCheck3D(
    Mortar::Element& sele, std::vector<Mortar::Element*> meles)
{
  double sxi_test[2] = {0.0, 0.0};
  bool proj_test = false;
  bool boundary_ele = false;
  double alpha_test = 0.0;
  double glob_test[3] = {0.0, 0.0, 0.0};
  const double tol = 1e-8;
  Core::Nodes::Node** mynodes_test = sele.Nodes();
  if (!mynodes_test) FOUR_C_THROW("has_proj_status: Null pointer!");

  Core::FE::CellType dt_s = sele.Shape();

  if (dt_s == Core::FE::CellType::quad4)  //|| dt_s==Core::FE::CellType::quad8
                                          //|| dt_s==Core::FE::CellType::quad9)
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
        Mortar::Projector::Impl(sele, *meles[bs_test])
            ->ProjectGaussPoint3D(sele, sxi_test, *meles[bs_test], mxi_test, alpha_test);
        Core::FE::CellType dt = meles[bs_test]->Shape();

        if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
            dt == Core::FE::CellType::quad9)
        {
          if (mxi_test[0] >= -1.0 && mxi_test[1] >= -1.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0)
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.num_node(); ++ii)
            {
              Mortar::Node* mycnode_test = dynamic_cast<Mortar::Node*>(mynodes_test[ii]);
              if (!mycnode_test) FOUR_C_THROW("has_proj_status: Null pointer!");

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
        else if (dt == Core::FE::CellType::tri3 || dt == Core::FE::CellType::tri6)
        {
          if (mxi_test[0] >= 0.0 && mxi_test[1] >= 0.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0 && mxi_test[0] + mxi_test[1] <= 1.0)
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.num_node(); ++ii)
            {
              Mortar::Node* mycnode_test = dynamic_cast<Mortar::Node*>(mynodes_test[ii]);
              if (!mycnode_test) FOUR_C_THROW("has_proj_status: Null pointer!");

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
          FOUR_C_THROW("Non valid element type for master discretization!");
        }
      }
      if (proj_test == false) boundary_ele = true;
    }
  }
  else if (dt_s == Core::FE::CellType::quad9)  //|| dt_s==Core::FE::CellType::quad8 ||
                                               // dt_s==Core::FE::CellType::quad9)
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
        Mortar::Projector::Impl(sele, *meles[bs_test])
            ->ProjectGaussPoint3D(sele, sxi_test, *meles[bs_test], mxi_test, alpha_test);
        Core::FE::CellType dt = meles[bs_test]->Shape();

        if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
            dt == Core::FE::CellType::quad9)
        {
          if (mxi_test[0] >= -1.0 && mxi_test[1] >= -1.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0)
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.num_node(); ++ii)
            {
              Mortar::Node* mycnode_test = dynamic_cast<Mortar::Node*>(mynodes_test[ii]);
              if (!mycnode_test) FOUR_C_THROW("has_proj_status: Null pointer!");

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
        else if (dt == Core::FE::CellType::tri3 || dt == Core::FE::CellType::tri6)
        {
          if (mxi_test[0] >= 0.0 && mxi_test[1] >= 0.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0 && mxi_test[0] + mxi_test[1] <= 1.0)
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.num_node(); ++ii)
            {
              Mortar::Node* mycnode_test = dynamic_cast<Mortar::Node*>(mynodes_test[ii]);
              if (!mycnode_test) FOUR_C_THROW("has_proj_status: Null pointer!");

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
          FOUR_C_THROW("Non valid element type for master discretization!");
        }
      }
      if (proj_test == false) boundary_ele = true;
    }
  }
  else if (dt_s == Core::FE::CellType::quad8)  //||
                                               // dt_s==Core::FE::CellType::quad8
                                               //|| dt_s==Core::FE::CellType::quad9)
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
        Mortar::Projector::Impl(sele, *meles[bs_test])
            ->ProjectGaussPoint3D(sele, sxi_test, *meles[bs_test], mxi_test, alpha_test);
        Core::FE::CellType dt = meles[bs_test]->Shape();

        if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
            dt == Core::FE::CellType::quad9)
        {
          if (mxi_test[0] >= -1.0 && mxi_test[1] >= -1.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0)
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.num_node(); ++ii)
            {
              Mortar::Node* mycnode_test = dynamic_cast<Mortar::Node*>(mynodes_test[ii]);
              if (!mycnode_test) FOUR_C_THROW("has_proj_status: Null pointer!");

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
        else if (dt == Core::FE::CellType::tri3 || dt == Core::FE::CellType::tri6)
        {
          if (mxi_test[0] >= 0.0 && mxi_test[1] >= 0.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0 && mxi_test[0] + mxi_test[1] <= 1.0)
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.num_node(); ++ii)
            {
              Mortar::Node* mycnode_test = dynamic_cast<Mortar::Node*>(mynodes_test[ii]);
              if (!mycnode_test) FOUR_C_THROW("has_proj_status: Null pointer!");

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
          FOUR_C_THROW("Non valid element type for master discretization!");
        }
      }
      if (proj_test == false) boundary_ele = true;
    }
  }
  // TRI-ELE
  else if (dt_s == Core::FE::CellType::tri3)
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
        Mortar::Projector::Impl(sele, *meles[bs_test])
            ->ProjectGaussPoint3D(sele, sxi_test, *meles[bs_test], mxi_test, alpha_test);
        Core::FE::CellType dt = meles[bs_test]->Shape();

        if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
            dt == Core::FE::CellType::quad9)
        {
          if (mxi_test[0] >= -1.0 && mxi_test[1] >= -1.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0)  // Falls Position auf Element
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.num_node(); ++ii)
            {
              Mortar::Node* mycnode_test = dynamic_cast<Mortar::Node*>(mynodes_test[ii]);
              if (!mycnode_test) FOUR_C_THROW("has_proj_status: Null pointer!");

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
        else if (dt == Core::FE::CellType::tri3 || dt == Core::FE::CellType::tri6)
        {
          if (mxi_test[0] >= 0.0 && mxi_test[1] >= 0.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0 && mxi_test[0] + mxi_test[1] <= 1.0)  // Falls Position auf Element
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.num_node(); ++ii)
            {
              Mortar::Node* mycnode_test = dynamic_cast<Mortar::Node*>(mynodes_test[ii]);
              if (!mycnode_test) FOUR_C_THROW("has_proj_status: Null pointer!");

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
          FOUR_C_THROW("Non valid element type for master discretization!");
        }
      }
      if (proj_test == false) boundary_ele = true;
    }
  }
  else if (dt_s == Core::FE::CellType::tri6)
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
        Mortar::Projector::Impl(sele, *meles[bs_test])
            ->ProjectGaussPoint3D(sele, sxi_test, *meles[bs_test], mxi_test, alpha_test);
        Core::FE::CellType dt = meles[bs_test]->Shape();

        if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
            dt == Core::FE::CellType::quad9)
        {
          if (mxi_test[0] >= -1.0 && mxi_test[1] >= -1.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0)  // Falls Position auf Element
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.num_node(); ++ii)
            {
              Mortar::Node* mycnode_test = dynamic_cast<Mortar::Node*>(mynodes_test[ii]);
              if (!mycnode_test) FOUR_C_THROW("has_proj_status: Null pointer!");

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
        else if (dt == Core::FE::CellType::tri3 || dt == Core::FE::CellType::tri6)
        {
          if (mxi_test[0] >= 0.0 && mxi_test[1] >= 0.0 && mxi_test[0] <= 1.0 &&
              mxi_test[1] <= 1.0 && mxi_test[0] + mxi_test[1] <= 1.0)  // Falls Position auf Element
          {
            // get hasproj
            sele.LocalToGlobal(sxi_test, glob_test, 0);
            for (int ii = 0; ii < sele.num_node(); ++ii)
            {
              Mortar::Node* mycnode_test = dynamic_cast<Mortar::Node*>(mynodes_test[ii]);
              if (!mycnode_test) FOUR_C_THROW("has_proj_status: Null pointer!");

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
          FOUR_C_THROW("Non valid element type for master discretization!");
        }
      }
      if (proj_test == false) boundary_ele = true;
    }
  }
  else
  {
    FOUR_C_THROW("Non valid element type for slave discretization!");
  }

  return boundary_ele;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::IntegrateDerivEle3D(Mortar::Element& sele,
    std::vector<Mortar::Element*> meles, bool* boundary_ele, bool* proj_, const Epetra_Comm& comm,
    const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr)
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
void CONTACT::Integrator::IntegrateDerivEle3D(Mortar::Element& sele,
    std::vector<Mortar::Element*> meles, bool* boundary_ele, bool* proj_, const Epetra_Comm& comm,
    const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  // explicitly defined shape function type needed
  if (shape_fcn() == Inpar::Mortar::shape_undefined)
    FOUR_C_THROW("IntegrateDerivEle3D called without specific shape function defined!");

  // Petrov-Galerkin approach for LM not yet implemented
  if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin &&
      sele.Shape() != Core::FE::CellType::nurbs9)
    FOUR_C_THROW("Petrov-Galerkin approach not yet implemented for 3-D quadratic FE interpolation");

  // check for problem dimension
  if (Dim() != 3) FOUR_C_THROW("3D integration method called for non-3D problem");

  // get slave element nodes themselves for normal evaluation
  Core::Nodes::Node** mynodes = sele.Nodes();
  if (!mynodes) FOUR_C_THROW("IntegrateDerivEle3D: Null pointer!");

  // check input data
  for (int test = 0; test < (int)meles.size(); ++test)
  {
    if (((!sele.IsSlave()) || (meles[test]->IsSlave())) and (!imortar_.get<bool>("Two_half_pass")))
      FOUR_C_THROW("IntegrateDerivEle3D called on a wrong type of Mortar::Element pair!");
  }

  int msize = meles.size();
  int nrow = sele.num_node();
  int ndof = dynamic_cast<Mortar::Node*>(sele.Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  Core::LinAlg::SerialDenseVector sval(nrow);
  Core::LinAlg::SerialDenseMatrix sderiv(nrow, 2, true);
  Core::LinAlg::SerialDenseVector lmval(nrow);
  Core::LinAlg::SerialDenseMatrix lmderiv(nrow, 2, true);

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix> dualmap(
      nrow * ndof, 0, Core::LinAlg::SerialDenseMatrix(nrow, nrow));
  if ((shape_fcn() == Inpar::Mortar::shape_dual ||
          shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) &&
      (sele.Shape() != Core::FE::CellType::tri3 ||
          sele.MoData().DerivDualShape() != Teuchos::null) &&
      lag_mult_quad() != Inpar::Mortar::lagmult_const)
    sele.DerivShapeDual(dualmap);

  // check whether quadratic elements with linear interpolation are present
  bool linlm = false;
  if ((shape_fcn() == Inpar::Mortar::shape_dual ||
          shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) &&
      (lag_mult_quad() == Inpar::Mortar::lagmult_lin) &&
      (sele.Shape() == Core::FE::CellType::quad8 || sele.Shape() == Core::FE::CellType::tri6))
    linlm = true;

  //********************************************************************
  //  Boundary_segmentation test -- HasProj() check
  //  if a slave-node has no projection onto each master element
  //  --> Boundary_ele==true
  //********************************************************************
  Inpar::Mortar::IntType integrationtype =
      Core::UTILS::IntegralValue<Inpar::Mortar::IntType>(imortar_, "INTTYPE");

  //************************************************************************
  // Boundary Segmentation check -- HasProj()-check
  //************************************************************************
  if (sele.Shape() != Core::FE::CellType::nurbs4 and sele.Shape() != Core::FE::CellType::nurbs9)
    *boundary_ele = BoundarySegmCheck3D(sele, meles);

  int linsize = 0;
  for (int i = 0; i < nrow; ++i)
  {
    Node* cnode = dynamic_cast<Node*>(mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  // Start integration if fast integration should be used or if there is no boundary element
  // for the fast_BS integration
  if (*boundary_ele == false || integrationtype == Inpar::Mortar::inttype_elements)
  {
    //**********************************************************************
    // loop over all Gauss points for integration
    //**********************************************************************
    for (int gp = 0; gp < nGP(); ++gp)
    {
      int iter_proj = 0;
      // coordinates and weight
      std::array<double, 2> eta = {Coordinate(gp, 0), Coordinate(gp, 1)};
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
      if (lag_mult_quad() == Inpar::Mortar::lagmult_const)
        sele.evaluate_shape_lag_mult_const(shape_fcn(), sxi, lmval, lmderiv, nrow);
      else if (lag_mult_quad() == Inpar::Mortar::lagmult_lin)
        sele.evaluate_shape_lag_mult_lin(shape_fcn(), sxi, lmval, lmderiv, nrow);
      else
        sele.evaluate_shape_lag_mult(shape_fcn(), sxi, lmval, lmderiv, nrow);

      // evaluate trace space shape functions (on both elements)
      sele.evaluate_shape(sxi, sval, sderiv, nrow, linlm);

      // evaluate the two Jacobians (int. cell and slave element)
      double jacslave = sele.Jacobian(sxi);

      // evaluate linearizations *******************************************
      // evaluate the slave Jacobian derivative
      Core::Gen::Pairedvector<int, double> jacslavemap(nrow * ndof + linsize);
      sele.DerivJacobian(sxi, jacslavemap);

      //**********************************************************************
      // loop over all mele
      //**********************************************************************
      for (int nummaster = 0; nummaster < msize; ++nummaster)
      {
        Core::FE::CellType dt = meles[nummaster]->Shape();

        int nmnode = meles[nummaster]->num_node();
        Core::LinAlg::SerialDenseVector mval(nmnode);
        Core::LinAlg::SerialDenseMatrix mderiv(nmnode, 2, true);


        // project Gauss point onto master element
        Mortar::Projector::Impl(sele, *meles[nummaster])
            ->ProjectGaussPoint3D(sele, sxi, *meles[nummaster], mxi, projalpha);

        // check GP projection
        const double tol = 0.00;
        is_on_mele = true;
        if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
            dt == Core::FE::CellType::quad9 || dt == Core::FE::CellType::nurbs9)
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
          meles[nummaster]->evaluate_shape(mxi, mval, mderiv, nmnode);

          // evaluate the GP slave coordinate derivatives
          std::vector<Core::Gen::Pairedvector<int, double>> dsxigp(2, 0);
          std::vector<Core::Gen::Pairedvector<int, double>> dmxigp(2, 4 * linsize + nmnode * ndof);
          DerivXiGP3D(sele, *meles[nummaster], sxi, mxi, dsxigp, dmxigp, projalpha);

          //**********************************************************************
          // geometric quantities
          //**********************************************************************
          double gpn[3] = {0.0, 0.0, 0.0};
          Core::Gen::Pairedvector<int, double> dgapgp(
              (nmnode * ndof) + linsize);  // gap lin. without lm and jac.
          double gap = 0.0;
          std::vector<Core::Gen::Pairedvector<int, double>> dnmap_unit(
              3, ((nmnode * ndof) + linsize));  // deriv of x,y and z comp. of gpn (unit)

          //**********************************************************************
          // evaluate at GP and lin char. quantities
          //**********************************************************************
          if (linlm)
          {
            // declare and compute shape functions as well as derivatives for computation of gap
            // -> standard quadratic functions instead of only linear functions
            Core::LinAlg::SerialDenseVector svalgap(nrow);
            Core::LinAlg::SerialDenseMatrix sderivgap(nrow, 2, true);
            sele.evaluate_shape(sxi, svalgap, sderivgap, nrow);

            gap_3_d(sele, *meles[nummaster], svalgap, mval, sderivgap, mderiv, &gap, gpn, dsxigp,
                dmxigp, dgapgp, dnmap_unit);
          }
          else
            gap_3_d(sele, *meles[nummaster], sval, mval, sderiv, mderiv, &gap, gpn, dsxigp, dmxigp,
                dgapgp, dnmap_unit);

          if (algo_ == Inpar::Mortar::algorithm_mortar)
            dynamic_cast<CONTACT::Element&>(sele).PrepareMderiv(meles, nummaster);

          integrate_gp_3_d(sele, *meles[nummaster], sval, lmval, mval, sderiv, mderiv, lmderiv,
              dualmap, wgt, jacslave, jacslavemap, gpn, dnmap_unit, gap, dgapgp, sxi, mxi, dsxigp,
              dmxigp);

          // assemble m-matrix for this slave/master pair
          if (algo_ == Inpar::Mortar::algorithm_mortar)
            dynamic_cast<CONTACT::Element&>(sele).assemble_mderiv_to_nodes(*meles[nummaster]);

        }  // is_on_mele
        if (is_on_mele == true) break;
      }  // mele loop

      // warning, if an element which is declared not to be on the boundary by the above test
      // has non-projectable Gauss points
      if (is_on_mele == false && *boundary_ele == false)
      {
        std::cout << "*** warning *** Non-boundary element has non-projectable Gauss point \n";
        Core::Nodes::Node** snodes = sele.Nodes();
        for (int e = 0; e < sele.num_node(); ++e)
        {
          Node* cnode = dynamic_cast<Node*>(snodes[e]);
          std::cout << "#eID " << sele.Id() << " | #nID " << cnode->Id() << ": " << cnode->X()[0]
                    << ", " << cnode->X()[1] << ", " << cnode->X()[2] << std::endl;
        }
      }

      // if one gp has counterparts on 2 elements --> non-uniqueness
      if (iter_proj > 1) FOUR_C_THROW("Multiple feasible projections of one integration point!");
    }  // GP-loop
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::integrate_deriv_cell3_d_aux_plane(Mortar::Element& sele,
    Mortar::Element& mele, Teuchos::RCP<Mortar::IntCell> cell, double* auxn,
    const Epetra_Comm& comm, const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr)
{
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr = Teuchos::null;
  if (not mparams_ptr.is_null())
  {
    cparams_ptr = Teuchos::rcp_dynamic_cast<CONTACT::ParamsInterface>(mparams_ptr, true);
  }
  integrate_deriv_cell3_d_aux_plane(sele, mele, cell, auxn, comm, cparams_ptr);
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
void CONTACT::Integrator::integrate_deriv_cell3_d_aux_plane(Mortar::Element& sele,
    Mortar::Element& mele, Teuchos::RCP<Mortar::IntCell> cell, double* auxn,
    const Epetra_Comm& comm, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  // explicitly defined shape function type needed
  if (shape_fcn() == Inpar::Mortar::shape_undefined)
    FOUR_C_THROW(
        "integrate_deriv_cell3_d_aux_plane called without specific shape function defined!");

  // check for problem dimension
  if (Dim() != 3) FOUR_C_THROW("3D integration method called for non-3D problem");

  // discretization type of master element
  Core::FE::CellType sdt = sele.Shape();
  Core::FE::CellType mdt = mele.Shape();

  // check input data
  if (((!sele.IsSlave()) || (mele.IsSlave())) and (!imortar_.get<bool>("Two_half_pass")))
    FOUR_C_THROW(
        "integrate_deriv_cell3_d_aux_plane called on a wrong type of Mortar::Element pair!");
  if (cell == Teuchos::null)
    FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane called without integration cell");

  // number of nodes (slave, master)
  int nrow = sele.num_node();
  int ncol = mele.num_node();
  int ndof = Dim();

  // get slave element nodes themselves for normal evaluation
  Core::Nodes::Node** mynodes = sele.Nodes();
  if (!mynodes) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

  // create empty vectors for shape fct. evaluation
  Core::LinAlg::SerialDenseVector sval(nrow);
  Core::LinAlg::SerialDenseMatrix sderiv(nrow, 2, true);
  Core::LinAlg::SerialDenseVector mval(ncol);
  Core::LinAlg::SerialDenseMatrix mderiv(ncol, 2, true);
  Core::LinAlg::SerialDenseVector lmval(nrow);
  Core::LinAlg::SerialDenseMatrix lmderiv(nrow, 2, true);
  Core::LinAlg::SerialDenseVector lmvalScale(nrow);
  Core::LinAlg::SerialDenseMatrix lmderivScale(nrow, 2, true);



  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  int linsize = 0;
  if (cppnormal_)
  {
    for (int i = 0; i < ncol; ++i)
    {
      Node* cnode = dynamic_cast<Node*>(mele.Nodes()[i]);
      linsize += cnode->GetLinsize();
    }
    linsize += 1 + mele.num_node();
  }
  else
  {
    for (int i = 0; i < nrow; ++i)
    {
      Node* cnode = dynamic_cast<Node*>(mynodes[i]);
      linsize += cnode->GetLinsize();
    }
  }

  // for cpp normal we have more lin entries
  linsize *= 10;

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix> dualmap(
      (nrow + ncol) * ndof, 0, Core::LinAlg::SerialDenseMatrix(nrow, nrow));
  if ((shape_fcn() == Inpar::Mortar::shape_dual ||
          shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) &&
      (sele.Shape() != Core::FE::CellType::tri3 || sele.MoData().DerivDualShape() != Teuchos::null))
    sele.DerivShapeDual(dualmap);

  // check if the cells are tri3
  // there's nothing wrong about other shapes, but as long as they are all
  // tri3 we can perform the jacobian calculation ( and its deriv) outside
  // the Gauss point loop
  if (cell->Shape() != Core::FE::CellType::tri3)
    FOUR_C_THROW("only tri3 integration cells at the moment. See comment in the code");

  double jac = cell->Jacobian();
  // directional derivative of cell Jacobian
  Core::Gen::Pairedvector<int, double> jacintcellmap((nrow + ncol) * ndof);
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
    Mortar::Projector::Impl(sele)->project_gauss_point_auxn3_d(globgp, auxn, sele, sxi, sprojalpha);
    Mortar::Projector::Impl(mele)->project_gauss_point_auxn3_d(globgp, auxn, mele, mxi, mprojalpha);

    // check GP projection (SLAVE)
    double tol = 0.01;
    if (sdt == Core::FE::CellType::quad4 || sdt == Core::FE::CellType::quad8 ||
        sdt == Core::FE::CellType::quad9)
    {
      if (sxi[0] < -1.0 - tol || sxi[1] < -1.0 - tol || sxi[0] > 1.0 + tol || sxi[1] > 1.0 + tol)
      {
        std::cout
            << "\n***Warning: integrate_deriv_cell3_d_aux_plane: Gauss point projection outside!";
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
        std::cout
            << "\n***Warning: integrate_deriv_cell3_d_aux_plane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }

    // check GP projection (MASTER)
    if (mdt == Core::FE::CellType::quad4 || mdt == Core::FE::CellType::quad8 ||
        mdt == Core::FE::CellType::quad9)
    {
      if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol)
      {
        std::cout
            << "\n***Warning: integrate_deriv_cell3_d_aux_plane: Gauss point projection outside!";
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
        std::cout
            << "\n***Warning: integrate_deriv_cell3_d_aux_plane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }

    // evaluate Lagrange multiplier shape functions (on slave element)
    sele.evaluate_shape_lag_mult(shape_fcn(), sxi, lmval, lmderiv, nrow);
    sele.evaluate_shape_lag_mult(shape_fcn(), sxi, lmvalScale, lmderivScale, nrow, false);

    // evaluate trace space shape functions (on both elements)
    sele.evaluate_shape(sxi, sval, sderiv, nrow);
    mele.evaluate_shape(mxi, mval, mderiv, ncol);

    // evaluate linearizations *******************************************

    // evaluate global GP coordinate derivative
    static Core::LinAlg::Matrix<3, 1> svalcell;
    static Core::LinAlg::Matrix<3, 2> sderivcell;
    cell->evaluate_shape(eta, svalcell, sderivcell);

    Core::Gen::Pairedvector<int, Core::LinAlg::Matrix<3, 1>> lingp((nrow + ncol) * ndof);

    for (int v = 0; v < 3; ++v)
      for (int d = 0; d < 3; ++d)
        for (_CI p = (cell->GetDerivVertex(v))[d].begin(); p != (cell->GetDerivVertex(v))[d].end();
             ++p)
          lingp[p->first](d) += svalcell(v) * (p->second);

    // evalute the GP slave coordinate derivatives
    std::vector<Core::Gen::Pairedvector<int, double>> dsxigp(2, (nrow + ncol) * ndof);
    DerivXiGP3DAuxPlane(sele, sxi, cell->Auxn(), dsxigp, sprojalpha, cell->get_deriv_auxn(), lingp);

    // evalute the GP master coordinate derivatives
    std::vector<Core::Gen::Pairedvector<int, double>> dmxigp(2, (nrow + ncol) * ndof);
    DerivXiGP3DAuxPlane(mele, mxi, cell->Auxn(), dmxigp, mprojalpha, cell->get_deriv_auxn(), lingp);

    //**********************************************************************
    // geometric quantities
    //**********************************************************************
    double gpn[3] = {0.0, 0.0, 0.0};
    Core::Gen::Pairedvector<int, double> dgapgp(
        (ncol * ndof) + linsize);  // gap lin. without lm and jac.
    double gap = 0.0;
    std::vector<Core::Gen::Pairedvector<int, double>> dnmap_unit(
        3, linsize);  // deriv of x,y and z comp. of gpn (unit)
    gap_3_d(sele, mele, sval, mval, sderiv, mderiv, &gap, gpn, dsxigp, dmxigp, dgapgp, dnmap_unit);

    // integrate for area scaling:
    for (int i = 0; i < sele.num_node(); ++i)
    {
      dynamic_cast<CONTACT::Node*>(sele.Nodes()[i])->MoData().GetDscale() +=
          wgt * jac * sval[i] * lmvalScale[i];
    }

    //**********************************************************************
    // evaluate at GP and lin char. quantities
    //**********************************************************************
    if (!nonsmoothselfcontactsurface_)
      integrate_gp_3_d(sele, mele, sval, lmval, mval, sderiv, mderiv, lmderiv, dualmap, wgt, jac,
          jacintcellmap, gpn, dnmap_unit, gap, dgapgp, sxi, mxi, dsxigp, dmxigp);
    // for non-smooth (self) contact surfaces we do not use the smoothed normal, but instead we use
    // the normal that is already used for the projection
    else
      integrate_gp_3_d(sele, mele, sval, lmval, mval, sderiv, mderiv, lmderiv, dualmap, wgt, jac,
          jacintcellmap, cell->Auxn(), cell->get_deriv_auxn(), gap, dgapgp, sxi, mxi, dsxigp,
          dmxigp);

  }  // end gp loop
  //**********************************************************************

  return;
}



/*----------------------------------------------------------------------*
 |  Integrate and linearize a 1D slave / master cell (2D)    farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::integrate_deriv_cell3_d_aux_plane_stl(Mortar::Element& mele,
    Mortar::Element& lele, Mortar::Element& sele, Teuchos::RCP<Mortar::IntCell> cell, double* auxn,
    const Epetra_Comm& comm)
{
  // explicitly defined shape function type needed
  if (shape_fcn() == Inpar::Mortar::shape_undefined)
    FOUR_C_THROW(
        "integrate_deriv_cell3_d_aux_plane called without specific shape function defined!");

  // check for problem dimension
  if (Dim() != 3) FOUR_C_THROW("3D integration method called for non-3D problem");

  // discretization type of master element
  Core::FE::CellType sdt = sele.Shape();
  Core::FE::CellType mdt = mele.Shape();

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    FOUR_C_THROW(
        "integrate_deriv_cell3_d_aux_plane called on a wrong type of Mortar::Element pair!");
  if (cell == Teuchos::null)
    FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane called without integration cell");

  // number of nodes (slave, master)
  //  int nrowL = lele.num_node();
  int nrow = sele.num_node();
  int ncolP = mele.num_node();
  int ncolL = lele.num_node();
  int ndof = Dim();

  // get slave element nodes themselves for normal evaluation
  Core::Nodes::Node** mynodes = sele.Nodes();
  if (!mynodes) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane_lts: Null pointer!");
  Core::Nodes::Node** mnodes = lele.Nodes();
  if (!mnodes) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane_lts: Null pointer!");

  // create empty vectors for shape fct. evaluation
  Core::LinAlg::SerialDenseVector sval(nrow);
  Core::LinAlg::SerialDenseMatrix sderiv(nrow, 2, true);
  Core::LinAlg::SerialDenseVector mval(ncolL);
  Core::LinAlg::SerialDenseMatrix mderiv(ncolL, 1, true);
  Core::LinAlg::SerialDenseVector lmval(nrow);
  Core::LinAlg::SerialDenseMatrix lmderiv(nrow, 2, true);

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix> dualmap(
      (nrow + ncolL) * ndof, 0, Core::LinAlg::SerialDenseMatrix(nrow, nrow));
  if ((shape_fcn() == Inpar::Mortar::shape_dual ||
          shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) &&
      (lele.Shape() != Core::FE::CellType::line2 ||
          sele.MoData().DerivDualShape() != Teuchos::null))
    sele.DerivShapeDual(dualmap);

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  int linsize = 0;
  for (int i = 0; i < nrow; ++i)
  {
    Node* cnode = dynamic_cast<Node*>(mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  // check if the cells are tri3
  // there's nothing wrong about other shapes, but as long as they are all
  // tri3 we can perform the jacobian calculation ( and its deriv) outside
  // the Gauss point loop
  if (cell->Shape() != Core::FE::CellType::line2)
    FOUR_C_THROW("only line2 integration cells for LTS");

  double jac = cell->Jacobian();

  // directional derivative of cell Jacobian
  Core::Gen::Pairedvector<int, double> jacintcellmap((nrow + ncolL) * ndof);
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
    Mortar::Projector::Impl(sele)->project_gauss_point_auxn3_d(globgp, auxn, sele, sxi, sprojalpha);
    Mortar::Projector::Impl(mele)->project_gauss_point_auxn3_d(globgp, auxn, mele, mxi, mprojalpha);

    // check GP projection (SLAVE)
    double tol = 0.01;
    if (sdt == Core::FE::CellType::quad4 || sdt == Core::FE::CellType::quad8 ||
        sdt == Core::FE::CellType::quad9)
    {
      if (sxi[0] < -1.0 - tol || sxi[1] < -1.0 - tol || sxi[0] > 1.0 + tol || sxi[1] > 1.0 + tol)
      {
        std::cout
            << "\n***Warning: integrate_deriv_cell3_d_aux_plane: Gauss point projection outside!";
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
        std::cout
            << "\n***Warning: integrate_deriv_cell3_d_aux_plane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }

    // check GP projection (MASTER)
    if (mdt == Core::FE::CellType::quad4 || mdt == Core::FE::CellType::quad8 ||
        mdt == Core::FE::CellType::quad9)
    {
      if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol)
      {
        std::cout
            << "\n***Warning: integrate_deriv_cell3_d_aux_plane: Gauss point projection outside!";
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
        std::cout
            << "\n***Warning: integrate_deriv_cell3_d_aux_plane: Gauss point projection outside!";
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
      FOUR_C_THROW("invalid local id of line ele");

    // evaluate Lagrange multiplier shape functions (on slave element)
    sele.evaluate_shape_lag_mult(shape_fcn(), sxi, lmval, lmderiv, nrow);

    // evaluate trace space shape functions (on both elements)
    sele.evaluate_shape(sxi, sval, sderiv, nrow);
    lele.evaluate_shape(lxi, mval, mderiv, ncolL);

    // evaluate linearizations *******************************************

    // evaluate global GP coordinate derivative
    static Core::LinAlg::Matrix<3, 1> svalcell;
    static Core::LinAlg::Matrix<3, 2> sderivcell;
    cell->evaluate_shape(eta, svalcell, sderivcell);

    Core::Gen::Pairedvector<int, Core::LinAlg::Matrix<3, 1>> lingp(100 * (nrow + ncolP) * ndof);

    for (int v = 0; v < 2; ++v)
      for (int d = 0; d < 3; ++d)
        for (_CI p = (cell->GetDerivVertex(v))[d].begin(); p != (cell->GetDerivVertex(v))[d].end();
             ++p)
          lingp[p->first](d) += svalcell(v) * (p->second);

    // evalute the GP slave coordinate derivatives
    std::vector<Core::Gen::Pairedvector<int, double>> dsxigp(2, (nrow + ncolP) * ndof);
    DerivXiGP3DAuxPlane(sele, sxi, cell->Auxn(), dsxigp, sprojalpha, cell->get_deriv_auxn(), lingp);

    // evalute the GP master coordinate derivatives
    std::vector<Core::Gen::Pairedvector<int, double>> dmxigp(2, (nrow + ncolP) * 100 * ndof);
    DerivXiGP3DAuxPlane(mele, mxi, cell->Auxn(), dmxigp, mprojalpha, cell->get_deriv_auxn(), lingp);

    // convert deriv sxi
    std::vector<Core::Gen::Pairedvector<int, double>> dlxigp(1, 100 * (nrow + ncolL) * ndof);

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
      FOUR_C_THROW("invalid local id of line ele");

    //**********************************************************************
    // geometric quantities
    //**********************************************************************
    std::array<double, 3> gpn = {0.0, 0.0, 0.0};
    Core::Gen::Pairedvector<int, double> dgapgp(
        (ncolL * ndof) + 10 * linsize);  // gap lin. without lm and jac.
    double gap = 0.0;
    std::vector<Core::Gen::Pairedvector<int, double>> dnmap_unit(
        3, 10 * linsize);  // deriv of x,y and z comp. of gpn (unit)

    //**********************************************************************
    // evaluate at GP and lin char. quantities
    //**********************************************************************
    std::array<double, 3> sgpx = {0.0, 0.0, 0.0};
    std::array<double, 3> mgpx = {0.0, 0.0, 0.0};

    for (int i = 0; i < nrow; ++i)
    {
      Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(mynodes[i]);
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
    if (lengthn < 1.0e-12) FOUR_C_THROW("IntegrateAndDerivSegment: Divide by zero!");

    for (int i = 0; i < 3; ++i) gpn[i] /= lengthn;

    // build gap function at current GP
    for (int i = 0; i < Dim(); ++i) gap += (mgpx[i] - sgpx[i]) * gpn[i];

    // **************************
    // Linearization
    // **************************
    int linsize = 0;
    for (int i = 0; i < nrow; ++i)
    {
      Node* cnode = dynamic_cast<Node*>(mynodes[i]);
      linsize += 10 * cnode->GetLinsize();
    }

    // build directional derivative of slave GP normal (non-unit)
    Core::Gen::Pairedvector<int, double> dmap_nxsl_gp(linsize);
    Core::Gen::Pairedvector<int, double> dmap_nysl_gp(linsize);
    Core::Gen::Pairedvector<int, double> dmap_nzsl_gp(linsize);

    for (int i = 0; i < nrow; ++i)
    {
      Node* cnode = dynamic_cast<Node*>(mynodes[i]);

      Core::Gen::Pairedvector<int, double>& dmap_nxsl_i = cnode->Data().GetDerivN()[0];
      Core::Gen::Pairedvector<int, double>& dmap_nysl_i = cnode->Data().GetDerivN()[1];
      Core::Gen::Pairedvector<int, double>& dmap_nzsl_i = cnode->Data().GetDerivN()[2];

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
      Node* cnode = dynamic_cast<Node*>(mynodes[z]);
      for (int k = 0; k < 3; ++k) dgapgp[cnode->Dofs()[k]] -= sval[z] * gpn[k];
    }

    for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < nrow; ++z)
      {
        Node* cnode = dynamic_cast<Node*>(mynodes[z]);
        for (int k = 0; k < 3; ++k) dg -= gpn[k] * sderiv(z, 0) * cnode->xspatial()[k] * ps;
      }
    }
    for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < nrow; ++z)
      {
        Node* cnode = dynamic_cast<Node*>(mynodes[z]);
        for (int k = 0; k < 3; ++k) dg -= gpn[k] * sderiv(z, 1) * cnode->xspatial()[k] * ps;
      }
    }

    //        MASTER
    // lin master nodes
    for (int z = 0; z < ncolL; ++z)
    {
      Node* cnode = dynamic_cast<Node*>(mnodes[z]);
      for (int k = 0; k < 3; ++k) dgapgp[cnode->Dofs()[k]] += mval[z] * gpn[k];
    }

    for (_CI p = dlxigp[0].begin(); p != dlxigp[0].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < ncolL; ++z)
      {
        Node* cnode = dynamic_cast<Node*>(mnodes[z]);
        for (int k = 0; k < 3; ++k) dg += gpn[k] * mderiv(z, 0) * cnode->xspatial()[k] * ps;
      }
    }



    // weighted gap
    for (int j = 0; j < nrow; ++j)
    {
      CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(mynodes[j]);

      double prod = 0.0;
      // Petrov-Galerkin approach (dual LM for D/M but standard LM for gap)
      if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) prod = sval[j] * gap * jac * wgt;
      // usual standard or dual LM approach
      else
        prod = lmval[j] * gap * jac * wgt;

      // do not process slave side boundary nodes
      // (their row entries would be zero anyway!)
      if (cnode->IsOnBound()) continue;

      // add current Gauss point's contribution to gseg
      cnode->AddltsGapValue(prod);
    }

    for (int iter = 0; iter < nrow; ++iter)
    {
      // map iterator
      typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

      Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(mynodes[iter]);
      if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

      static double fac = 0.0;

      // get the corresponding map as a reference
      std::map<int, double>& dgmap =
          dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivGlts();

      // switch if Petrov-Galerkin approach for LM is applied
      if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
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
        // only line2 --> not required!

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
    }

    // integrate D and M matrix
    // dual shape functions
    if (shape_fcn() == Inpar::Mortar::shape_dual ||
        shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
    {
      for (int j = 0; j < nrow; ++j)
      {
        CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(mynodes[j]);

        // integrate mseg
        for (int k = 0; k < ncolL; ++k)
        {
          CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(mnodes[k]);

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
        CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(mynodes[j]);

        // integrate dseg
        for (int k = 0; k < nrow; ++k)
        {
          CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(mynodes[k]);

          // multiply the two shape functions
          double prod = lmval[j] * sval[k] * jac * wgt;

          if (abs(prod) > MORTARINTTOL) cnode->AddDltsValue(mnode->Id(), prod);
          if (abs(prod) > MORTARINTTOL) cnode->AddSNode(mnode->Id());  // only for friction!
        }

        // integrate mseg
        for (int k = 0; k < ncolL; ++k)
        {
          CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(mnodes[k]);

          // multiply the two shape functions
          double prod = lmval[j] * mval[k] * jac * wgt;

          if (abs(prod) > MORTARINTTOL) cnode->AddMltsValue(mnode->Id(), prod);
          if (abs(prod) > MORTARINTTOL) cnode->AddMNode(mnode->Id());  // only for friction!
        }
      }
    }

    if (shape_fcn() == Inpar::Mortar::shape_dual ||
        shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
    {
      // integrate LinD
      for (int j = 0; j < nrow; ++j)
      {
        Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(mynodes[j]);
        if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

        int sgid = mymrtrnode->Id();
        std::map<int, double>& ddmap_jk =
            dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivDlts()[sgid];


        // integrate LinM
        for (int k = 0; k < ncolL; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          static double fac = 0.0;

          // get the correct map as a reference
          std::map<int, double>& dmmap_jk =
              dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivMlts()[mgid];

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
        Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(mynodes[j]);
        if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");


        // integrate LinD
        for (int k = 0; k < nrow; ++k)
        {
          // global master node ID
          int mgid = sele.Nodes()[k]->Id();
          static double fac = 0.0;

          // get the correct map as a reference
          std::map<int, double>& ddmap_jk =
              dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivDlts()[mgid];

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
              dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivMlts()[mgid];

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
void CONTACT::Integrator::integrate_deriv_cell3_d_aux_plane_lts(Mortar::Element& sele,
    Mortar::Element& lsele, Mortar::Element& mele, Teuchos::RCP<Mortar::IntCell> cell, double* auxn,
    const Epetra_Comm& comm)
{
  // explicitly defined shape function type needed
  if (shape_fcn() == Inpar::Mortar::shape_undefined)
    FOUR_C_THROW(
        "integrate_deriv_cell3_d_aux_plane called without specific shape function defined!");

  // check for problem dimension
  if (Dim() != 3) FOUR_C_THROW("3D integration method called for non-3D problem");

  // discretization type of master element
  Core::FE::CellType sdt = sele.Shape();
  Core::FE::CellType mdt = mele.Shape();

  // check input data
  //  if ((!sele.IsSlave()) || (mele.IsSlave()))
  //    FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane called on a wrong type of Mortar::Element
  //    pair!");
  if (cell == Teuchos::null)
    FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane called without integration cell");

  // number of nodes (slave, master)
  int nrowL = lsele.num_node();
  int nrowS = sele.num_node();
  int ncol = mele.num_node();
  int ndof = Dim();

  // get slave element nodes themselves for normal evaluation
  Core::Nodes::Node** mynodes = lsele.Nodes();
  if (!mynodes)
    FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane_lts: Cannot access slave ndoes. Null pointer!");
  Core::Nodes::Node** mnodes = mele.Nodes();
  if (!mnodes)
    FOUR_C_THROW(
        "integrate_deriv_cell3_d_aux_plane_lts: Cannot access master nodes. Null pointer!");

  // create empty vectors for shape fct. evaluation
  Core::LinAlg::SerialDenseVector sval(nrowL);
  Core::LinAlg::SerialDenseMatrix sderiv(nrowL, 1, true);
  Core::LinAlg::SerialDenseVector mval(ncol);
  Core::LinAlg::SerialDenseMatrix mderiv(ncol, 2, true);
  Core::LinAlg::SerialDenseVector lmval(nrowL);
  Core::LinAlg::SerialDenseMatrix lmderiv(nrowL, 1, true);

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix> dualmap(
      (nrowL + ncol) * ndof, 0, Core::LinAlg::SerialDenseMatrix(nrowL, nrowL));
  if ((shape_fcn() == Inpar::Mortar::shape_dual ||
          shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) &&
      (lsele.Shape() != Core::FE::CellType::line2 ||
          sele.MoData().DerivDualShape() != Teuchos::null))
    sele.DerivShapeDual(dualmap);

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  int linsize = 0;
  for (int i = 0; i < nrowL; ++i)
  {
    Node* cnode = dynamic_cast<Node*>(mynodes[i]);
    linsize += cnode->GetLinsize();
  }
  linsize *= 100;  // TODO: remove this safety scaling

  // check if the cells are tri3
  // there's nothing wrong about other shapes, but as long as they are all
  // tri3 we can perform the jacobian calculation ( and its deriv) outside
  // the Gauss point loop
  if (cell->Shape() != Core::FE::CellType::line2)
    FOUR_C_THROW("only line2 integration cells for LTS");

  const double jac = cell->Jacobian();

  // directional derivative of cell Jacobian
  Core::Gen::Pairedvector<int, double> jacintcellmap((nrowL + ncol) * ndof + linsize);
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
    bool success = Mortar::Projector::Impl(sele)->project_gauss_point_auxn3_d(
        globgp, auxn, sele, sxi, sprojalpha);
    if (!success)
    {
      //      FOUR_C_THROW("projection not possible!");
      return;
    }

    success = Mortar::Projector::Impl(mele)->project_gauss_point_auxn3_d(
        globgp, auxn, mele, mxi, mprojalpha);
    if (!success)
    {
      //      FOUR_C_THROW("projection not possible!");
      return;
    }

    // check
    double checkpos[3] = {0.0, 0.0, 0.0};
    mele.LocalToGlobal(mxi, checkpos, 0);

    // check GP projection (SLAVE)
    const double sminedge = sele.MinEdgeSize();
    const double mminedge = mele.MinEdgeSize();
    const double tol = MORTARCLIPTOL * std::min(sminedge, mminedge);
    if (sdt == Core::FE::CellType::quad4 || sdt == Core::FE::CellType::quad8 ||
        sdt == Core::FE::CellType::quad9)
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
    if (mdt == Core::FE::CellType::quad4 || mdt == Core::FE::CellType::quad8 ||
        mdt == Core::FE::CellType::quad9)
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
      FOUR_C_THROW("invalid local id of line ele");

    // evaluate Lagrange multiplier shape functions (on slave element)
    lsele.evaluate_shape_lag_mult(shape_fcn(), lxi, lmval, lmderiv, nrowL);

    // evaluate trace space shape functions (on both elements)
    lsele.evaluate_shape(lxi, sval, sderiv, nrowL);
    mele.evaluate_shape(mxi, mval, mderiv, ncol);

    // evaluate linearizations *******************************************

    // evaluate global GP coordinate derivative
    static Core::LinAlg::Matrix<3, 1> svalcell;
    static Core::LinAlg::Matrix<3, 2> sderivcell;
    cell->evaluate_shape(eta, svalcell, sderivcell);

    Core::Gen::Pairedvector<int, Core::LinAlg::Matrix<3, 1>> lingp((nrowS + ncol) * ndof + linsize);

    for (int v = 0; v < 2; ++v)
      for (int d = 0; d < 3; ++d)
        for (_CI p = (cell->GetDerivVertex(v))[d].begin(); p != (cell->GetDerivVertex(v))[d].end();
             ++p)
          lingp[p->first](d) += svalcell(v) * (p->second);

    // evalute the GP slave coordinate derivatives
    std::vector<Core::Gen::Pairedvector<int, double>> dsxigp(2, (nrowS + ncol) * ndof + linsize);
    DerivXiGP3DAuxPlane(sele, sxi, cell->Auxn(), dsxigp, sprojalpha, cell->get_deriv_auxn(), lingp);

    // evalute the GP master coordinate derivatives
    std::vector<Core::Gen::Pairedvector<int, double>> dmxigp(2, (nrowS + ncol) * ndof + linsize);
    DerivXiGP3DAuxPlane(mele, mxi, cell->Auxn(), dmxigp, mprojalpha, cell->get_deriv_auxn(), lingp);

    // convert deriv sxi
    std::vector<Core::Gen::Pairedvector<int, double>> dlxigp(1, (nrowS + ncol) * ndof + linsize);

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
      FOUR_C_THROW("invalid local id of line ele");

    //**********************************************************************
    // geometric quantities
    //**********************************************************************
    std::array<double, 3> gpn = {0.0, 0.0, 0.0};
    Core::Gen::Pairedvector<int, double> dgapgp(
        (ncol * ndof) + 10 * linsize);  // gap lin. without lm and jac.
    double gap = 0.0;
    std::vector<Core::Gen::Pairedvector<int, double>> dnmap_unit(
        3, 10 * linsize);  // deriv of x,y and z comp. of gpn (unit)

    //**********************************************************************
    // evaluate at GP and lin char. quantities
    //**********************************************************************
    std::array<double, 3> sgpx = {0.0, 0.0, 0.0};
    std::array<double, 3> mgpx = {0.0, 0.0, 0.0};

    for (int i = 0; i < nrowL; ++i)
    {
      Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(mynodes[i]);
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
    if (lengthn < 1.0e-12) FOUR_C_THROW("IntegrateAndDerivSegment: Divide by zero!");

    for (int i = 0; i < 3; ++i) gpn[i] /= lengthn;

    // build gap function at current GP
    for (int i = 0; i < Dim(); ++i) gap += (mgpx[i] - sgpx[i]) * gpn[i];

    // **************************
    // Linearization
    // **************************
    int linsize = 0;
    for (int i = 0; i < nrowL; ++i)
    {
      Node* cnode = dynamic_cast<Node*>(mynodes[i]);
      linsize += 10 * cnode->GetLinsize();
    }

    // build directional derivative of slave GP normal (non-unit)
    Core::Gen::Pairedvector<int, double> dmap_nxsl_gp(linsize);
    Core::Gen::Pairedvector<int, double> dmap_nysl_gp(linsize);
    Core::Gen::Pairedvector<int, double> dmap_nzsl_gp(linsize);

    for (int i = 0; i < nrowL; ++i)
    {
      Node* cnode = dynamic_cast<Node*>(mynodes[i]);

      Core::Gen::Pairedvector<int, double>& dmap_nxsl_i = cnode->Data().GetDerivN()[0];
      Core::Gen::Pairedvector<int, double>& dmap_nysl_i = cnode->Data().GetDerivN()[1];
      Core::Gen::Pairedvector<int, double>& dmap_nzsl_i = cnode->Data().GetDerivN()[2];

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
      Node* cnode = dynamic_cast<Node*>(mynodes[z]);
      for (int k = 0; k < 3; ++k) dgapgp[cnode->Dofs()[k]] -= sval[z] * gpn[k];
    }

    for (_CI p = dlxigp[0].begin(); p != dlxigp[0].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < nrowL; ++z)
      {
        Node* cnode = dynamic_cast<Node*>(mynodes[z]);
        for (int k = 0; k < 3; ++k) dg -= gpn[k] * sderiv(z, 0) * cnode->xspatial()[k] * ps;
      }
    }
    //        MASTER
    // lin master nodes
    for (int z = 0; z < ncol; ++z)
    {
      Node* cnode = dynamic_cast<Node*>(mnodes[z]);
      for (int k = 0; k < 3; ++k) dgapgp[cnode->Dofs()[k]] += mval[z] * gpn[k];
    }

    for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < ncol; ++z)
      {
        Node* cnode = dynamic_cast<Node*>(mnodes[z]);
        for (int k = 0; k < 3; ++k) dg += gpn[k] * mderiv(z, 0) * cnode->xspatial()[k] * ps;
      }
    }

    for (_CI p = dmxigp[1].begin(); p != dmxigp[1].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < ncol; ++z)
      {
        Node* cnode = dynamic_cast<Node*>(mnodes[z]);
        for (int k = 0; k < 3; ++k) dg += gpn[k] * mderiv(z, 1) * cnode->xspatial()[k] * ps;
      }
    }

    // weighted gap
    for (int j = 0; j < nrowL; ++j)
    {
      CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(mynodes[j]);

      double prod = 0.0;
      // Petrov-Galerkin approach (dual LM for D/M but standard LM for gap)
      if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) prod = sval[j] * gap * jac * wgt;
      // usual standard or dual LM approach
      else
        prod = lmval[j] * gap * jac * wgt;

      // do not process slave side boundary nodes
      // (their row entries would be zero anyway!)
      if (cnode->IsOnCorner()) continue;

      // add current Gauss point's contribution to gseg
      cnode->AddltsGapValue(prod);
    }

    for (int iter = 0; iter < nrowL; ++iter)
    {
      // map iterator
      typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

      Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(mynodes[iter]);
      if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

      if (mymrtrnode->IsOnCorner()) continue;

      static double fac = 0.0;

      // get the corresponding map as a reference
      std::map<int, double>& dgmap =
          dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivGlts();

      // switch if Petrov-Galerkin approach for LM is applied
      if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
      {
        // (1) Lin(Phi) - does not exist here for Petrov-Galerkin approach

        // (2) Lin(N) - slave GP coordinates
        for (int d = 0; d < Dim() - 2; ++d)
        {
          fac = wgt * sderiv(iter, d) * gap * jac;
          for (_CI p = dlxigp[d].begin(); p != dlxigp[d].end(); ++p)
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
        // only line2 --> not required!

        // (2) Lin(Phi) - slave GP coordinates
        for (int d = 0; d < Dim() - 2; ++d)
        {
          fac = wgt * lmderiv(iter, d) * gap * jac;
          for (_CI p = dlxigp[d].begin(); p != dlxigp[d].end(); ++p)
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
    }
    //--------------------------------------------------------
    // integrate D and M matrix
    // dual shape functions
    if (shape_fcn() == Inpar::Mortar::shape_dual ||
        shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
    {
      bool bound = false;
      for (int j = 0; j < nrowL; ++j)
      {
        CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(mynodes[j]);
        if (cnode->IsOnCorner()) bound = true;
      }

      for (int j = 0; j < nrowL; ++j)
      {
        CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(mynodes[j]);
        if (cnode->IsOnCorner()) continue;

        // integrate mseg
        for (int k = 0; k < ncol; ++k)
        {
          CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(mnodes[k]);

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
            CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(mynodes[k]);
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
      //        CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(mynodes[j]);
      //        if (cnode->IsOnCorner())
      //          continue;
      //
      //        // integrate mseg
      //        for (int k=0; k<ncol; ++k)
      //        {
      //          CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(mnodes[k]);
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
        CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(mynodes[j]);

        if (cnode->IsOnCorner()) continue;

        // integrate mseg
        for (int k = 0; k < ncol; ++k)
        {
          CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(mnodes[k]);

          // multiply the two shape functions
          double prod = lmval[j] * mval[k] * jac * wgt;

          if (abs(prod) > MORTARINTTOL) cnode->AddMltsValue(mnode->Id(), prod);
          if (abs(prod) > MORTARINTTOL) cnode->AddMNode(mnode->Id());  // only for friction!
        }

        // integrate dseg
        for (int k = 0; k < nrowL; ++k)
        {
          CONTACT::Node* snode = dynamic_cast<CONTACT::Node*>(mynodes[k]);

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
    if (shape_fcn() == Inpar::Mortar::shape_dual ||
        shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
    {
      bool bound = false;
      for (int j = 0; j < nrowL; ++j)
      {
        CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(mynodes[j]);
        if (cnode->IsOnCorner()) bound = true;
      }

      if (bound)
      {
        // integrate LinD
        for (int j = 0; j < nrowL; ++j)
        {
          Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(mynodes[j]);
          if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

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
                dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivMlts()[mgid];

            // (1) Lin(Phi) - dual shape functions
            //          if (dualmap.size()>0)
            //          {
            //            for (int m=0;m<nrowL;++m)
            //            {
            //              fac = wgt*sval[m]*mval[k]*jac;
            //              for
            //              (Core::Gen::Pairedvector<int,Core::LinAlg::SerialDenseMatrix>::const_iterator
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
            Mortar::Node* mymrtrnode2 = dynamic_cast<Mortar::Node*>(mynodes[k]);
            if (!mymrtrnode2) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

            // global master node ID
            int sgid = mymrtrnode2->Id();
            static double fac = 0.0;

            // node k is boundary node
            if (mymrtrnode2->IsOnCorner())
            {
              // move entry to derivM (with minus sign)
              // get the correct map as a reference
              std::map<int, double>& dmmap_jk =
                  dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivDlts()[sgid];

              // (1) Lin(Phi) - dual shape functions
              //            if (dualmap.size()>0)
              //            {
              //              for (int m=0;m<nrow;++m)
              //              {
              //                fac = wgt*sval[m]*sval[k]*jac;
              //                for
              //                (Core::Gen::Pairedvector<int,Core::LinAlg::SerialDenseMatrix>::const_iterator
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
              std::map<int, double>& ddmap_jk =
                  dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivDlts()[mymrtrnode->Id()];

              // (1) Lin(Phi) - dual shape functions
              //            if (dualmap.size()>0)
              //            {
              //              for (int m=0;m<nrow;++m)
              //              {
              //                fac = wgt*sval[m]*sval[k]*jac;
              //                for
              //                (Core::Gen::Pairedvector<int,Core::LinAlg::SerialDenseMatrix>::const_iterator
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
          Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(mynodes[j]);
          if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

          int sgid = mymrtrnode->Id();
          std::map<int, double>& ddmap_jk =
              dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivDlts()[sgid];

          // integrate LinM
          for (int k = 0; k < ncol; ++k)
          {
            // global master node ID
            int mgid = mele.Nodes()[k]->Id();
            static double fac = 0.0;

            // get the correct map as a reference
            std::map<int, double>& dmmap_jk =
                dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivMlts()[mgid];

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
        Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(mynodes[j]);
        if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

        if (mymrtrnode->IsOnCorner()) continue;

        // integrate LinM
        for (int k = 0; k < ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          static double fac = 0.0;

          // get the correct map as a reference
          std::map<int, double>& dmmap_jk =
              dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivMlts()[mgid];

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

          if (dynamic_cast<Node*>(mynodes[k])->IsOnCorner())
          {
            // get the correct map as a reference
            std::map<int, double>& dmmap_jk =
                dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivDlts()[mgid];

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
                dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivDlts()[mgid];

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
void CONTACT::Integrator::integrate_deriv_cell3_d_aux_plane_quad(Mortar::Element& sele,
    Mortar::Element& mele, Mortar::IntElement& sintele, Mortar::IntElement& mintele,
    Teuchos::RCP<Mortar::IntCell> cell, double* auxn)
{
  // get LMtype
  Inpar::Mortar::LagMultQuad lmtype = lag_mult_quad();

  // explicitly defined shape function type needed
  if (shape_fcn() == Inpar::Mortar::shape_undefined)
    FOUR_C_THROW(
        "integrate_deriv_cell3_d_aux_plane_quad called without specific shape function defined!");

  // Petrov-Galerkin approach for LM not yet implemented
  if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin &&
      sele.Shape() != Core::FE::CellType::nurbs9)
    FOUR_C_THROW("Petrov-Galerkin approach not yet implemented for 3-D quadratic FE interpolation");

  // check for problem dimension
  if (Dim() != 3) FOUR_C_THROW("3D integration method called for non-3D problem");

  // discretization type of slave and master IntElement
  Core::FE::CellType sdt = sintele.Shape();
  Core::FE::CellType mdt = mintele.Shape();

  // discretization type of slave and master Element
  Core::FE::CellType psdt = sele.Shape();
  Core::FE::CellType pmdt = mele.Shape();

  // check input data
  if (((!sele.IsSlave()) || (mele.IsSlave())) and (!imortar_.get<bool>("Two_half_pass")))
    FOUR_C_THROW(
        "integrate_deriv_cell3_d_aux_plane_quad called on a wrong type of Mortar::Element pair!");
  if (cell == Teuchos::null)
    FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane_quad called without integration cell");

  // contact with wear
  bool wear = false;
  if (wearlaw_ != Inpar::Wear::wear_none) wear = true;

  // number of nodes (slave, master)
  int nrow = sele.num_node();
  int ncol = mele.num_node();
  int ndof = Dim();
  int nintrow = sintele.num_node();

  // create empty vectors for shape fct. evaluation
  Core::LinAlg::SerialDenseVector sval(nrow);
  Core::LinAlg::SerialDenseMatrix sderiv(nrow, 2, true);
  Core::LinAlg::SerialDenseVector svalmod(nrow);
  Core::LinAlg::SerialDenseMatrix sderivmod(nrow, 2, true);
  Core::LinAlg::SerialDenseVector mval(ncol);
  Core::LinAlg::SerialDenseMatrix mderiv(ncol, 2, true);
  Core::LinAlg::SerialDenseVector lmval(nrow);
  Core::LinAlg::SerialDenseMatrix lmderiv(nrow, 2, true);
  Core::LinAlg::SerialDenseVector lmintval(nintrow);
  Core::LinAlg::SerialDenseMatrix lmintderiv(nintrow, 2, true);

  // create empty vectors for shape fct. evaluation
  Core::LinAlg::SerialDenseMatrix ssecderiv(nrow, 3);

  // get slave and master nodal coords for Jacobian / GP evaluation
  Core::LinAlg::SerialDenseMatrix scoord(3, sele.num_node());
  sele.GetNodalCoords(scoord);
  Core::LinAlg::SerialDenseMatrix mcoord(3, mele.num_node());
  mele.GetNodalCoords(mcoord);

  // nodal coords from previous time step and lagrange mulitplier
  Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> scoordold;
  Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> mcoordold;
  Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> lagmult;

  // get them in the case of wear
  if (wear or Core::UTILS::IntegralValue<int>(imortar_, "GP_SLIP_INCR") == true)
  {
    scoordold = Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(3, sele.num_node()));
    mcoordold = Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(3, mele.num_node()));
    lagmult = Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(3, sele.num_node()));
    sele.GetNodalCoordsOld(*scoordold);
    mele.GetNodalCoordsOld(*mcoordold);
    sele.GetNodalLagMult(*lagmult);
  }

  // get slave element nodes themselves for normal evaluation
  Core::Nodes::Node** mynodes = sele.Nodes();
  if (!mynodes) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane_quad: Null pointer!");
  Core::Nodes::Node** myintnodes = sintele.Nodes();
  if (!myintnodes) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane_quad: Null pointer!");

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  bool duallin = false;
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix> dualmap(
      (nrow + ncol) * ndof, 0, Core::LinAlg::SerialDenseMatrix(nrow, nrow));
  if ((shape_fcn() == Inpar::Mortar::shape_dual ||
          shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) &&
      (sele.Shape() != Core::FE::CellType::tri3 ||
          sele.MoData().DerivDualShape() != Teuchos::null) &&
      lag_mult_quad() != Inpar::Mortar::lagmult_const)
  {
    duallin = true;
    sele.DerivShapeDual(dualmap);
  }

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k = 0; k < nrow; ++k)
  {
    Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(mynodes[k]);
    if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane_quad: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  // decide whether displacement shape fct. modification has to be considered or not
  // this is the case for dual quadratic Lagrange multipliers on quad8 and tri6 elements
  // check whether quadratic element with linear interpolation is present
  bool dualquad3d = false;
  bool linlm = false;
  if ((shape_fcn() == Inpar::Mortar::shape_dual ||
          shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) &&
      (lmtype == Inpar::Mortar::lagmult_quad || lmtype == Inpar::Mortar::lagmult_lin) &&
      (sele.Shape() == Core::FE::CellType::quad8 || sele.Shape() == Core::FE::CellType::tri6))
  {
    if (lmtype == Inpar::Mortar::lagmult_quad)
      dualquad3d = true;
    else
      linlm = true;
  }

  // get linearization size
  int linsize = 0;
  for (int i = 0; i < nrow; ++i)
  {
    Node* cnode = dynamic_cast<Node*>(mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  // check if the cells are tri3
  // there's nothing wrong about other shapes, but as long as they are all
  // tri3 we can perform the jacobian calculation ( and its deriv) outside
  // the Gauss point loop
  if (cell->Shape() != Core::FE::CellType::tri3)
    FOUR_C_THROW("only tri3 integration cells at the moment. See comment in the code");

  double jac = cell->Jacobian();
  // directional derivative of cell Jacobian
  Core::Gen::Pairedvector<int, double> jacintcellmap((nrow + ncol) * ndof);
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
    Mortar::Projector::Impl(sintele)->project_gauss_point_auxn3_d(
        globgp, auxn, sintele, sxi, sprojalpha);
    Mortar::Projector::Impl(mintele)->project_gauss_point_auxn3_d(
        globgp, auxn, mintele, mxi, mprojalpha);

    // check GP projection (SLAVE)
    const double tol = 0.01;
    if (sdt == Core::FE::CellType::quad4 || sdt == Core::FE::CellType::quad8 ||
        sdt == Core::FE::CellType::quad9)
    {
      if (sxi[0] < -1.0 - tol || sxi[1] < -1.0 - tol || sxi[0] > 1.0 + tol || sxi[1] > 1.0 + tol)
      {
        std::cout
            << "\n***Warning: integrate_deriv_cell3_d_aux_plane_quad: Slave Gauss point projection "
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
        std::cout
            << "\n***Warning: integrate_deriv_cell3_d_aux_plane_quad: Slave Gauss point projection "
               "outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP (IntElement) projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }

    // check GP projection (MASTER)
    if (mdt == Core::FE::CellType::quad4 || mdt == Core::FE::CellType::quad8 ||
        mdt == Core::FE::CellType::quad9)
    {
      if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol)
      {
        std::cout << "\n***Warning: integrate_deriv_cell3_d_aux_plane_quad: Master Gauss point "
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
        std::cout << "\n***Warning: integrate_deriv_cell3_d_aux_plane_quad: Master Gauss point "
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
    Mortar::Projector::Impl(sele)->project_gauss_point_auxn3_d(
        globgp, auxn, sele, psxi, psprojalpha);
    Mortar::Projector::Impl(mele)->project_gauss_point_auxn3_d(
        globgp, auxn, mele, pmxi, pmprojalpha);
    // sintele.MapToParent(sxi, psxi); // old way of doing it via affine map... wrong (popp 05/2016)
    // mintele.MapToParent(mxi, pmxi); // old way of doing it via affine map... wrong (popp 05/2016)

    // check GP projection (SLAVE)
    if (psdt == Core::FE::CellType::quad4 || psdt == Core::FE::CellType::quad8 ||
        psdt == Core::FE::CellType::quad9 || psdt == Core::FE::CellType::nurbs9)
    {
      if (psxi[0] < -1.0 - tol || psxi[1] < -1.0 - tol || psxi[0] > 1.0 + tol ||
          psxi[1] > 1.0 + tol)
      {
        std::cout
            << "\n***Warning: integrate_deriv_cell3_d_aux_plane_quad: Slave Gauss point projection "
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
        std::cout
            << "\n***Warning: integrate_deriv_cell3_d_aux_plane_quad: Slave Gauss point projection "
               "outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << psxi[0] << " " << psxi[1] << std::endl;
      }
    }

    // check GP projection (MASTER)
    if (pmdt == Core::FE::CellType::quad4 || pmdt == Core::FE::CellType::quad8 ||
        pmdt == Core::FE::CellType::quad9 || pmdt == Core::FE::CellType::nurbs9)
    {
      if (pmxi[0] < -1.0 - tol || pmxi[1] < -1.0 - tol || pmxi[0] > 1.0 + tol ||
          pmxi[1] > 1.0 + tol)
      {
        std::cout << "\n***Warning: integrate_deriv_cell3_d_aux_plane_quad: Master Gauss point "
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
        std::cout << "\n***Warning: integrate_deriv_cell3_d_aux_plane_quad: Master Gauss point "
                     "projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << pmxi[0] << " " << pmxi[1] << std::endl;
      }
    }

    // evaluate Lagrange multiplier shape functions (on slave element)
    if (lag_mult_quad() == Inpar::Mortar::lagmult_const)
      sele.evaluate_shape_lag_mult_const(shape_fcn(), sxi, lmval, lmderiv, nrow);
    else if (bound)
      sele.evaluate_shape_lag_mult(shape_fcn(), psxi, lmval, lmderiv, nrow);
    else
    {
      sele.evaluate_shape_lag_mult(shape_fcn(), psxi, lmval, lmderiv, nrow);
      sintele.evaluate_shape_lag_mult(shape_fcn(), sxi, lmintval, lmintderiv, nintrow, false);
    }

    // evaluate trace space shape functions (on both elements)
    sele.evaluate_shape(psxi, sval, sderiv, nrow, linlm);
    if (dualquad3d) sele.evaluate_shape(psxi, svalmod, sderivmod, nrow, true);
    mele.evaluate_shape(pmxi, mval, mderiv, ncol);

    // evaluate 2nd deriv of trace space shape functions (on slave element)
    sele.evaluate2nd_deriv_shape(psxi, ssecderiv, nrow);

    // evaluate global GP coordinate derivative
    Core::LinAlg::Matrix<3, 1> svalcell;
    Core::LinAlg::Matrix<3, 2> sderivcell;
    cell->evaluate_shape(eta, svalcell, sderivcell);

    Core::Gen::Pairedvector<int, Core::LinAlg::Matrix<3, 1>> lingp((nrow + ncol) * ndof);

    for (int v = 0; v < 3; ++v)
      for (int d = 0; d < 3; ++d)
        for (_CI p = (cell->GetDerivVertex(v))[d].begin(); p != (cell->GetDerivVertex(v))[d].end();
             ++p)
          lingp[p->first](d) += svalcell(v) * (p->second);

    // evalute the GP slave coordinate derivatives
    std::vector<Core::Gen::Pairedvector<int, double>> dsxigp(2, (nrow + ncol) * ndof);
    DerivXiGP3DAuxPlane(
        sintele, sxi, cell->Auxn(), dsxigp, sprojalpha, cell->get_deriv_auxn(), lingp);

    // evalute the GP master coordinate derivatives
    std::vector<Core::Gen::Pairedvector<int, double>> dmxigp(2, (nrow + ncol) * ndof);
    DerivXiGP3DAuxPlane(
        mintele, mxi, cell->Auxn(), dmxigp, mprojalpha, cell->get_deriv_auxn(), lingp);

    // evaluate the GP slave coordinate derivatives (parent element)
    // evaluate the GP master coordinate derivatives (parent element)
    std::vector<Core::Gen::Pairedvector<int, double>> dpsxigp(2, (nrow + ncol) * ndof);
    std::vector<Core::Gen::Pairedvector<int, double>> dpmxigp(2, (nrow + ncol) * ndof);
    DerivXiGP3DAuxPlane(
        sele, psxi, cell->Auxn(), dpsxigp, psprojalpha, cell->get_deriv_auxn(), lingp);
    DerivXiGP3DAuxPlane(
        mele, pmxi, cell->Auxn(), dpmxigp, pmprojalpha, cell->get_deriv_auxn(), lingp);
    // sintele.MapToParent(dsxigp,dpsxigp); // old way of doing it via affine map... wrong (popp
    // 05/2016) mintele.MapToParent(dmxigp,dpmxigp); // old way of doing it via affine map... wrong
    // (popp 05/2016)

    //**********************************************************************
    // frequently reused quantities
    //**********************************************************************
    double gpn[3] = {0.0, 0.0, 0.0};
    double gap = 0.0;
    double jumpval[2] = {
        0.0, 0.0};  // jump for wear
                    // double jumpvalv[2] = {0.0,0.0};  // jump for slipincr --> equal to jumpval
    double wearval = 0.0;  // wear value
    double lengthn = 0.0;  // length of gp normal gpn
    Core::Gen::Pairedvector<int, double> dsliptmatrixgp(
        (ncol * ndof) + linsize);  // deriv. of slip for wear
    Core::Gen::Pairedvector<int, double> dgapgp(
        (ncol * ndof) + linsize);  // gap lin. without lm and jac.
    Core::Gen::Pairedvector<int, double> dweargp(
        (ncol * ndof) + linsize);  // wear lin. without lm and jac.
    std::vector<Core::Gen::Pairedvector<int, double>> dslipgp(
        2, ((ncol * ndof) + linsize));  // deriv. of slip for slipincr (xi, eta)
    std::vector<Core::Gen::Pairedvector<int, double>> dnmap_unit(
        3, linsize);  // deriv of x,y and z comp. of gpn (unit)

    //**********************************************************************
    // evaluate at GP and lin char. quantities
    //**********************************************************************
    // these cases don't need any special treatment for quadratic elements
    if (algo_ == Inpar::Mortar::algorithm_gpts ||
        (lmtype == Inpar::Mortar::lagmult_const && algo_ == Inpar::Mortar::algorithm_mortar))
    {
      gap_3_d(
          sele, mele, sval, mval, sderiv, mderiv, &gap, gpn, dpsxigp, dpmxigp, dgapgp, dnmap_unit);
      integrate_gp_3_d(sele, mele, sval, lmval, mval, sderiv, mderiv, lmderiv, dualmap, wgt, jac,
          jacintcellmap, gpn, dnmap_unit, gap, dgapgp, psxi, pmxi, dpsxigp, dpmxigp);
    }
    // special treatment for quadratic elements
    else
    {
      // integrate and lin gp gap
      if (shape_fcn() == Inpar::Mortar::shape_standard && lmtype == Inpar::Mortar::lagmult_pwlin)
      {
        gp_3_d_g_quad_pwlin(sele, sintele, mele, sval, mval, lmintval, scoord, mcoord, sderiv,
            mderiv, &gap, gpn, &lengthn, jac, wgt, dsxigp, dmxigp, dgapgp, dnmap_unit);
      }
      else
      {
        if (linlm)
        {
          // declare and compute shape functions as well as derivatives for computation of gap
          // -> standard quadratic functions instead of only linear functions
          Core::LinAlg::SerialDenseVector svalgap(nrow);
          Core::LinAlg::SerialDenseMatrix sderivgap(nrow, 2, true);
          sele.evaluate_shape(psxi, svalgap, sderivgap, nrow);

          gap_3_d(sele, mele, svalgap, mval, sderivgap, mderiv, &gap, gpn, dpsxigp, dpmxigp, dgapgp,
              dnmap_unit);
        }
        else
          gap_3_d(sele, mele, sval, mval, sderiv, mderiv, &gap, gpn, dpsxigp, dpmxigp, dgapgp,
              dnmap_unit);

        gp_3_d_w_gap(sele, sval, lmval, &gap, jac, wgt, true);
      }

      // compute cell D/M matrix *******************************************
      gp_3_d_dm_quad(sele, mele, sintele, lmval, lmintval, sval, mval, jac, wgt, nrow, nintrow,
          ncol, ndof, bound);

      //*******************************
      // WEAR stuff
      //*******************************
      // std. wear for all wear-algorithm types
      if (wear)
        gp_3_d_wear(sele, mele, sval, sderiv, mval, mderiv, lmval, lmderiv, lagmult, gpn, jac, wgt,
            jumpval, &wearval, dsliptmatrixgp, dweargp, dsxigp, dmxigp, dnmap_unit, dualmap);

      // integrate T and E matrix for discr. wear
      if (wear_type() == Inpar::Wear::wear_primvar) gp_te(sele, lmval, sval, jac, wgt, jumpval);

      //********************************************************************
      // compute cell linearization
      //********************************************************************
      if (shape_fcn() == Inpar::Mortar::shape_standard && lmtype == Inpar::Mortar::lagmult_pwlin)
      {
        for (int iter = 0; iter < nintrow; ++iter)
        {
          // Lin DM
          gp_3_d_dm_quad_pwlin_lin(iter, sele, sintele, mele, sval, mval, lmintval, sderiv, mderiv,
              lmintderiv, wgt, jac, dsxigp, dpsxigp, dpmxigp, jacintcellmap);

          // Lin gap
          gp_3_d_g_quad_pwlin_lin(iter, sintele, sval, lmintval, sderiv, lmintderiv, gap, gpn, jac,
              wgt, dgapgp, jacintcellmap, dsxigp);
        }
      }
      else
      {
        // Lin DM
        gp_3_d_dm_quad_lin(duallin, sele, mele, sval, svalmod, mval, lmval, sderiv, mderiv, lmderiv,
            wgt, jac, dpsxigp, dpmxigp, jacintcellmap, dualmap, dualquad3d);

        for (int iter = 0; iter < nrow; ++iter)
        {
          // Lin gap
          gp_3_d_g_quad_lin(iter, sele, mele, sval, svalmod, lmval, sderiv, lmderiv, gap, gpn, jac,
              wgt, duallin, dgapgp, jacintcellmap, dpsxigp, dualmap, dualquad3d);

          // Lin wear matrices T and E for discr. wear
          if (wear_type() == Inpar::Wear::wear_primvar)
            gp_3_d_te_lin(iter, sele, sval, lmval, sderiv, lmderiv, jac, wgt, jumpval, dsxigp,
                jacintcellmap, dsliptmatrixgp, dualmap);
        }
      }
    }
  }  // gp loop

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::IntegrateDerivEle2D(Mortar::Element& sele,
    std::vector<Mortar::Element*> meles, bool* boundary_ele,
    const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr)
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
void CONTACT::Integrator::IntegrateDerivEle2D(Mortar::Element& sele,
    std::vector<Mortar::Element*> meles, bool* boundary_ele,
    const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  // ********************************************************************
  // Check integrator input for non-reasonable quantities
  // *********************************************************************

  // explicitly defined shape function type needed
  if (shape_fcn() == Inpar::Mortar::shape_undefined)
    FOUR_C_THROW("IntegrateDerivEle2D called without specific shape function defined!");

  // Petrov-Galerkin approach for LM not yet implemented for quadratic FE
  if (sele.Shape() == Core::FE::CellType::line3 &&
      shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
    FOUR_C_THROW("Petrov-Galerkin approach not yet implemented for 2-D quadratic FE interpolation");

  // check for problem dimension
  if (Dim() != 2) FOUR_C_THROW("2D integration method called for non-2D problem");

  // check input data
  for (int i = 0; i < (int)meles.size(); ++i)
  {
    if ((!sele.IsSlave()) || (meles[i]->IsSlave()))
      FOUR_C_THROW("IntegrateDerivEle2D called on a wrong type of Mortar::Element pair!");
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

  Core::Nodes::Node** mynodes = nullptr;
  nrow = sele.num_node();
  mynodes = sele.Nodes();

  // get slave element nodes themselves
  if (!mynodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  // create empty vectors for shape fct. evaluation
  Core::LinAlg::SerialDenseVector sval(nrow);
  Core::LinAlg::SerialDenseMatrix sderiv(nrow, 1);
  Core::LinAlg::SerialDenseVector lmval(nrow);
  Core::LinAlg::SerialDenseMatrix lmderiv(nrow, 1);

  // get slave nodal coords for Jacobian / GP evaluation
  Core::LinAlg::SerialDenseMatrix scoord(3, nrow);
  sele.GetNodalCoords(scoord);

  // nodal coords from previous time step and lagrange mulitplier
  Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> scoordold;
  Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> mcoordold;
  Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> lagmult;

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // prepare directional derivative of dual shape functions
  // this is only necessary for quadratic dual shape functions in 2D
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix> dualmap(
      nrow * ndof, 0, Core::LinAlg::SerialDenseMatrix(nrow, nrow));
  if ((shape_fcn() == Inpar::Mortar::shape_dual ||
          shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) &&
      (sele.Shape() == Core::FE::CellType::line3 ||
          sele.MoData().DerivDualShape() != Teuchos::null ||
          sele.Shape() == Core::FE::CellType::nurbs3))
    sele.DerivShapeDual(dualmap);

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  for (int k = 0; k < nrow; ++k)
  {
    Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(mynodes[k]);
    if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_segment2_d: Null pointer!");
  }

  // decide whether linear LM are used for quadratic FE here
  // and whether displacement shape fct. modification has to be considered or not
  // this is the case for dual linear Lagrange multipliers on line3 elements
  bool linlm = false;
  bool dualquad = false;
  if (lag_mult_quad() == Inpar::Mortar::lagmult_lin && sele.Shape() == Core::FE::CellType::line3)
  {
    linlm = true;
    if (shape_fcn() == Inpar::Mortar::shape_dual) dualquad = true;
  }

  // get numerical integration type
  Inpar::Mortar::IntType inttype =
      Core::UTILS::IntegralValue<Inpar::Mortar::IntType>(imortar_, "INTTYPE");

  //************************************************************************
  // Boundary Segmentation check -- HasProj()-check
  //************************************************************************
  if (inttype == Inpar::Mortar::inttype_elements_BS)
    *boundary_ele = BoundarySegmCheck2D(sele, meles);

  int linsize = 0;
  for (int i = 0; i < nrow; ++i)
  {
    Node* cnode = dynamic_cast<Node*>(mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  if (*boundary_ele == false || inttype == Inpar::Mortar::inttype_elements)
  {
    //*************************************************************************
    //                loop over all Gauss points for integration
    //*************************************************************************
    for (int gp = 0; gp < nGP(); ++gp)
    {
      bool kink_projection = false;

      // coordinates and weight
      std::array<double, 2> eta = {Coordinate(gp, 0), 0.0};
      double wgt = Weight(gp);

      // coordinate transformation sxi->eta (slave Mortar::Element->Overlap)
      double sxi[2] = {0.0, 0.0};
      sxi[0] = eta[0];

      // evaluate the two slave side Jacobians
      double dxdsxi = sele.Jacobian(sxi);

      // evaluate Lagrange multiplier shape functions (on slave element)
      if (lag_mult_quad() == Inpar::Mortar::lagmult_const)
        sele.evaluate_shape_lag_mult_const(shape_fcn(), sxi, lmval, lmderiv, nrow);
      else if (linlm)
        sele.evaluate_shape_lag_mult_lin(shape_fcn(), sxi, lmval, lmderiv, nrow);
      else
        sele.evaluate_shape_lag_mult(shape_fcn(), sxi, lmval, lmderiv, nrow);

      // evaluate trace space shape functions
      sele.evaluate_shape(sxi, sval, sderiv, nrow, dualquad);

      //****************************************************************************************************************
      //                loop over all Master Elements
      //****************************************************************************************************************
      for (int nummaster = 0; nummaster < (int)meles.size(); ++nummaster)
      {
        // project Gauss point onto master element
        double mxi[2] = {0.0, 0.0};
        Mortar::Projector::Impl(sele, *meles[nummaster])
            ->ProjectGaussPoint2D(sele, sxi, *meles[nummaster], mxi);

        // gp on mele?
        if ((mxi[0] >= -1.0) && (mxi[0] <= 1.0) && (kink_projection == false))
        {
          kink_projection = true;

          int ncol = 0;

          ncol = meles[nummaster]->num_node();

          Core::LinAlg::SerialDenseVector mval(ncol);
          Core::LinAlg::SerialDenseMatrix mderiv(ncol, 1);

          // get master nodal coords for Jacobian / GP evaluation
          Core::LinAlg::SerialDenseMatrix mcoord(3, ncol);
          meles[nummaster]->GetNodalCoords(mcoord);

          // evaluate trace space shape functions
          meles[nummaster]->evaluate_shape(mxi, mval, mderiv, ncol, false);

          // get directional derivatives of sxia, sxib, mxia, mxib --> derivatives of mxia/mxib not
          // required
          std::vector<Core::Gen::Pairedvector<int, double>> ximaps(4, linsize + ndof * ncol);
          bool startslave = true;
          bool endslave = true;
          double mxia = -0.1;  //--> arbitrary value
          double mxib = 0.1;   //--> arbitrary value
          DerivXiAB2D(sele, sxia, sxib, *meles[nummaster], mxia, mxib, ximaps, startslave, endslave,
              linsize);

          // evalute the GP slave coordinate derivatives --> no entries
          std::vector<Core::Gen::Pairedvector<int, double>> dsxigp(
              1, Core::Gen::Pairedvector<int, double>(linsize + ndof * ncol));
          for (_CI p = ximaps[0].begin(); p != ximaps[0].end(); ++p) dsxigp[0][p->first] = 0.0;

          // evalute the GP master coordinate derivatives
          std::vector<Core::Gen::Pairedvector<int, double>> dmxigp(
              1, Core::Gen::Pairedvector<int, double>(linsize + ndof * ncol));
          deriv_xi_g_p2_d(sele, *meles[nummaster], sxi[0], mxi[0], dsxigp[0], dmxigp[0], linsize);

          // evaluate the Jacobian derivative
          Core::Gen::Pairedvector<int, double> derivjac(nrow * ndof);
          sele.DerivJacobian(sxi, derivjac);  // direct derivative if xi^1_g does not change

          //**********************************************************************
          // frequently reused quantities
          //**********************************************************************
          double gpn[2] = {0.0, 0.0};  // normalized normal at gp
          std::vector<Core::Gen::Pairedvector<int, double>> dnmap_unit(
              2, (linsize + ndof * ncol));  // deriv of x and y comp. of gpn (unit)
          double gap = 0.0;                 // gap
          Core::Gen::Pairedvector<int, double> dgapgp(
              linsize + ndof * ncol);  // gap lin without weighting and jac

          //**********************************************************************
          // evaluate at GP and lin char. quantities
          //**********************************************************************
          if (dualquad && linlm)
          {
            // declare and compute shape functions as well as derivatives for computation of gap
            // -> standard quadratic functions instead of only linear functions
            Core::LinAlg::SerialDenseVector svalgap(nrow);
            Core::LinAlg::SerialDenseMatrix sderivgap(nrow, 1);
            sele.evaluate_shape(sxi, svalgap, sderivgap, nrow);

            gap_2_d(sele, *meles[nummaster], svalgap, mval, sderivgap, mderiv, &gap, gpn, dsxigp,
                dmxigp, dgapgp, dnmap_unit);
          }
          else
            gap_2_d(sele, *meles[nummaster], sval, mval, sderiv, mderiv, &gap, gpn, dsxigp, dmxigp,
                dgapgp, dnmap_unit);

          // integrate and lin gp gap
          integrate_gp_2_d(sele, *meles[nummaster], sval, lmval, mval, sderiv, mderiv, lmderiv,
              dualmap, wgt, dxdsxi, derivjac, gpn, dnmap_unit, gap, dgapgp, sxi, mxi, dsxigp,
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
void CONTACT::Integrator::IntegrateD(Mortar::Element& sele, const Epetra_Comm& comm, bool lin)
{
  // ********************************************************************
  // Check integrator input for non-reasonable quantities
  // *********************************************************************

  // explicitly defined shape function type needed
  if (shape_fcn() == Inpar::Mortar::shape_undefined)
    FOUR_C_THROW("IntegrateD called without specific shape function defined!");

  // *********************************************************************
  // Define slave quantities
  // *********************************************************************
  // number of nodes (slave, master)
  int ndof = Dim();
  int nrow = 0;

  Core::Nodes::Node** mynodes = nullptr;
  nrow = sele.num_node();
  mynodes = sele.Nodes();

  // get slave element nodes themselves
  if (!mynodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  // create empty vectors for shape fct. evaluation
  Core::LinAlg::SerialDenseVector sval(nrow);
  Core::LinAlg::SerialDenseMatrix sderiv(nrow, ndof - 1);
  Core::LinAlg::SerialDenseVector lmval(nrow);
  Core::LinAlg::SerialDenseMatrix lmderiv(nrow, ndof - 1);

  int linsize = 0;
  for (int i = 0; i < nrow; ++i)
  {
    Node* cnode = dynamic_cast<Node*>(mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k = 0; k < nrow; ++k)
  {
    Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(mynodes[k]);
    if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_segment2_d: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  bool duallin = false;
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix> dualmap(
      nrow * ndof, 0, Core::LinAlg::SerialDenseMatrix(nrow, nrow));
  if ((sele.Shape() != Core::FE::CellType::tri3 and sele.Shape() != Core::FE::CellType::line2) ||
      sele.MoData().DerivDualShape() != Teuchos::null)
  {
    duallin = true;
    sele.DerivShapeDual(dualmap);
  }

  // d-matrix derivative
  // local for this element to avoid direct assembly into the nodes for performance reasons
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix> dMatrixDeriv(
      sele.num_node() * Dim(), 0,
      Core::LinAlg::SerialDenseMatrix(sele.num_node(), sele.num_node()));

  //*************************************************************************
  //                loop over all Gauss points for integration
  //*************************************************************************
  for (int gp = 0; gp < nGP(); ++gp)
  {
    // coordinates and weight
    std::array<double, 2> eta = {Coordinate(gp, 0), 0.0};
    double wgt = Weight(gp);
    if (ndof == 3) eta[1] = Coordinate(gp, 1);

    // coordinate transformation sxi->eta (slave Mortar::Element->Overlap)
    double sxi[2] = {0.0, 0.0};
    sxi[0] = eta[0];
    sxi[1] = eta[1];

    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);

    // evaluate linearizations *******************************************
    // evaluate the slave Jacobian derivative
    Core::Gen::Pairedvector<int, double> jacslavemap(nrow * ndof + linsize);
    sele.DerivJacobian(sxi, jacslavemap);

    // evaluate Lagrange multiplier shape functions (on slave element)
    sele.evaluate_shape_lag_mult(Inpar::Mortar::shape_dual, sxi, lmval, lmderiv, nrow);

    // evaluate trace space shape functions
    sele.evaluate_shape(sxi, sval, sderiv, nrow, false);

    //**********************************************************************
    // evaluate at GP and lin char. quantities
    //**********************************************************************
    // integrate D and M
    for (int j = 0; j < nrow; ++j)
    {
      CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(mynodes[j]);

      // integrate dseg (boundary modification)
      if (bound)
      {
        bool j_boundnode = cnode->IsOnBound();

        for (int k = 0; k < nrow; ++k)
        {
          CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(mynodes[k]);
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
          CONTACT::Node* snode = dynamic_cast<CONTACT::Node*>(mynodes[k]);

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
      typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

      // (1) Lin(Phi) - dual shape functions
      if (duallin)
      {
        for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                 dualmap.begin();
             p != dualmap.end(); ++p)
        {
          Core::LinAlg::SerialDenseMatrix& dderivtmp = dMatrixDeriv[p->first];
          for (int j = 0; j < nrow; ++j)
            for (int k = 0; k < nrow; ++k)
              for (int m = 0; m < nrow; ++m)
                dderivtmp(j, j) += wgt * sval[m] * sval[k] * dxdsxi * (p->second)(j, m);
        }
      }

      // (4) Lin(dsxideta) - intcell GP Jacobian
      for (_CI p = jacslavemap.begin(); p != jacslavemap.end(); ++p)
      {
        Core::LinAlg::SerialDenseMatrix& dderivtmp = dMatrixDeriv[p->first];
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
void CONTACT::Integrator::integrate_kappa_penalty_lts(Mortar::Element& ele)
{
  // number of nodes (slave)
  int nrow = ele.num_node();

  // create empty objects for shape fct. evaluation
  Core::LinAlg::SerialDenseVector val(nrow);
  Core::LinAlg::SerialDenseMatrix deriv(nrow, 2, true);

  // get slave element nodes themselves
  Core::Nodes::Node** mynodes = ele.Nodes();

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
    ele.evaluate_shape_lag_mult(shape_fcn(), eta, val, deriv, nrow);

    // evaluate the Jacobian det
    const double jac = ele.Jacobian(eta);

    // compute cell gap vector *******************************************
    // loop over all gseg vector entries
    // nrow represents the slave side dofs !!!  */
    for (int j = 0; j < nrow; ++j)
    {
      Node* mymrtrnode = dynamic_cast<Node*>(mynodes[j]);

      if (mymrtrnode->IsOnCorner()) continue;

      // add current Gauss point's contribution to kappa
      mymrtrnode->Data().Kappa() += val[j] * jac * wgt;
    }
  }
  //**********************************************************************

  return;
}


/*----------------------------------------------------------------------*
 |  Compute penalty scaling factor kappa                      popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::integrate_kappa_penalty(Mortar::Element& sele, double* sxia, double* sxib,
    Teuchos::RCP<Core::LinAlg::SerialDenseVector> gseg)
{
  // explicitly defined shape function type needed
  if (shape_fcn() == Inpar::Mortar::shape_undefined)
    FOUR_C_THROW("integrate_kappa_penalty called without specific shape function defined!");

  // check input data
  if (!sele.IsSlave())
    FOUR_C_THROW("integrate_kappa_penalty called on a non-slave Mortar::Element!");
  if ((sxia[0] < -1.0) || (sxia[1] < -1.0) || (sxib[0] > 1.0) || (sxib[1] > 1.0))
    FOUR_C_THROW("integrate_kappa_penalty called with infeasible slave limits!");

  // number of nodes (slave)
  int nrow = sele.num_node();

  // create empty objects for shape fct. evaluation
  Core::LinAlg::SerialDenseVector val(nrow);
  Core::LinAlg::SerialDenseMatrix deriv(nrow, 2, true);

  // map iterator
  // typedef std::map<int,double>::const_iterator CI;

  // get slave element nodes themselves
  Core::Nodes::Node** mynodes = sele.Nodes();
  if (!mynodes) FOUR_C_THROW("integrate_kappa_penalty: Null pointer!");

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k = 0; k < nrow; ++k)
  {
    Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(mynodes[k]);
    if (!mymrtrnode) FOUR_C_THROW("integrate_kappa_penalty: Null pointer!");
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
      sele.evaluate_shape_lag_mult_lin(shape_fcn(), eta, val, deriv, nrow);
    else
      sele.evaluate_shape_lag_mult(shape_fcn(), eta, val, deriv, nrow);

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
void CONTACT::Integrator::integrate_kappa_penalty(Mortar::Element& sele,
    Mortar::IntElement& sintele, double* sxia, double* sxib,
    Teuchos::RCP<Core::LinAlg::SerialDenseVector> gseg)
{
  // get LMtype
  Inpar::Mortar::LagMultQuad lmtype = lag_mult_quad();

  // explicitly defined shape function type needed
  if (shape_fcn() != Inpar::Mortar::shape_standard)
    FOUR_C_THROW("integrate_kappa_penalty -> you should not be here!");

  // check input data
  if (!sele.IsSlave())
    FOUR_C_THROW("integrate_kappa_penalty called on a non-slave Mortar::Element!");
  if ((sxia[0] < -1.0) || (sxia[1] < -1.0) || (sxib[0] > 1.0) || (sxib[1] > 1.0))
    FOUR_C_THROW("integrate_kappa_penalty called with infeasible slave limits!");

  // number of nodes (slave)
  int nrow = sele.num_node();
  int nintrow = sintele.num_node();

  // create empty objects for shape fct. evaluation
  Core::LinAlg::SerialDenseVector val(nrow);
  Core::LinAlg::SerialDenseMatrix deriv(nrow, 2, true);
  Core::LinAlg::SerialDenseVector intval(nintrow);
  Core::LinAlg::SerialDenseMatrix intderiv(nintrow, 2, true);

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
    sintele.compute_unit_normal_at_xi(eta, auxn);

    // project Gauss point back to slave (parent) element
    double psxi[2] = {0.0, 0.0};
    double psprojalpha = 0.0;
    Mortar::Projector::Impl(sele)->project_gauss_point_auxn3_d(
        globgp, auxn, sele, psxi, psprojalpha);
    // sintele.MapToParent(eta,psxi); //old way of doing it via affine map... wrong (popp 05/2016)

    // evaluate shape functions
    sele.evaluate_shape(psxi, val, deriv, nrow);
    sintele.evaluate_shape(eta, intval, intderiv, nintrow);

    // evaluate the Jacobian det
    const double jac = sintele.Jacobian(eta);

    // compute cell gap vector *******************************************
    if (lmtype == Inpar::Mortar::lagmult_pwlin)
    {
      for (int j = 0; j < nintrow; ++j)
      {
        // add current Gauss point's contribution to gseg
        (*gseg)(j) += intval[j] * jac * wgt;
      }
    }
    else
    {
      FOUR_C_THROW("Invalid LM interpolation case!");
    }
    // compute cell gap vector *******************************************
  }
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 |  Compute directional derivative of XiAB (2D)               popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::DerivXiAB2D(Mortar::Element& sele, double& sxia, double& sxib,
    Mortar::Element& mele, double& mxia, double& mxib,
    std::vector<Core::Gen::Pairedvector<int, double>>& derivxi, bool& startslave, bool& endslave,
    int& linsize)
{
  // check for problem dimension
  if (Dim() != 2) FOUR_C_THROW("2D integration method called for non-2D problem");

  // we need the participating slave and master nodes
  Core::Nodes::Node** snodes = nullptr;
  Core::Nodes::Node** mnodes = nullptr;
  int numsnode = sele.num_node();
  int nummnode = mele.num_node();

  int ndof = Dim();

  snodes = sele.Nodes();
  mnodes = mele.Nodes();

  std::vector<Mortar::Node*> smrtrnodes(numsnode);
  std::vector<Mortar::Node*> mmrtrnodes(nummnode);

  for (int i = 0; i < numsnode; ++i)
  {
    smrtrnodes[i] = dynamic_cast<Mortar::Node*>(snodes[i]);
    if (!smrtrnodes[i]) FOUR_C_THROW("DerivXiAB2D: Null pointer!");
  }

  for (int i = 0; i < nummnode; ++i)
  {
    mmrtrnodes[i] = dynamic_cast<Mortar::Node*>(mnodes[i]);
    if (!mmrtrnodes[i]) FOUR_C_THROW("DerivXiAB2D: Null pointer!");
  }

  // we also need shape function derivs in A and B
  double psxia[2] = {sxia, 0.0};
  double psxib[2] = {sxib, 0.0};
  double pmxia[2] = {mxia, 0.0};
  double pmxib[2] = {mxib, 0.0};
  Core::LinAlg::SerialDenseVector valsxia(numsnode);
  Core::LinAlg::SerialDenseVector valsxib(numsnode);
  Core::LinAlg::SerialDenseVector valmxia(nummnode);
  Core::LinAlg::SerialDenseVector valmxib(nummnode);
  Core::LinAlg::SerialDenseMatrix derivsxia(numsnode, 1);
  Core::LinAlg::SerialDenseMatrix derivsxib(numsnode, 1);
  Core::LinAlg::SerialDenseMatrix derivmxia(nummnode, 1);
  Core::LinAlg::SerialDenseMatrix derivmxib(nummnode, 1);

  sele.evaluate_shape(psxia, valsxia, derivsxia, numsnode, false);
  sele.evaluate_shape(psxib, valsxib, derivsxib, numsnode, false);
  mele.evaluate_shape(pmxia, valmxia, derivmxia, nummnode, false);
  mele.evaluate_shape(pmxib, valmxib, derivmxib, nummnode, false);

  // prepare linearizations
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // compute leading constant for DerivXiBMaster if start node = slave node
  if (startslave == true)
  {
    // compute factors and leading constants for master
    double cmxib = 0.0;
    double fac_dxm_b = 0.0;
    double fac_dym_b = 0.0;
    double fac_xmsl_b = 0.0;
    double fac_ymsl_b = 0.0;

    // there are only actual 2 components for the normal in 2D but compute_unit_normal_at_xi expects
    // 3 components (otherwise invalid memory access)
    double normal[3] = {0., 0., 0.};
    sele.compute_unit_normal_at_xi(psxia, normal);
    std::vector<Core::Gen::Pairedvector<int, double>> derivN;
    dynamic_cast<CONTACT::Element*>(&sele)->DerivUnitNormalAtXi(psxia, derivN);

    Core::LinAlg::SerialDenseVector* mval = nullptr;
    Core::LinAlg::SerialDenseMatrix* mderiv = nullptr;
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

    Core::Gen::Pairedvector<int, double> dmap_mxib(nummnode * ndof + linsize);

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

    // there are only 2 actual components for the normal in 2D but compute_unit_normal_at_xi expects
    // 3 components (otherwise invalid memory access)
    double normal[3] = {0., 0., 0.};
    sele.compute_unit_normal_at_xi(psxib, normal);
    std::vector<Core::Gen::Pairedvector<int, double>> derivN;
    dynamic_cast<CONTACT::Element*>(&sele)->DerivUnitNormalAtXi(psxib, derivN);

    Core::LinAlg::SerialDenseVector* mval = nullptr;
    Core::LinAlg::SerialDenseMatrix* mderiv = nullptr;
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

    Core::Gen::Pairedvector<int, double> dmap_mxia(nummnode * ndof + linsize);

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

    Core::LinAlg::SerialDenseVector* mval = nullptr;
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


    Core::Gen::Pairedvector<int, double> dmap_sxia(nummnode * ndof + linsize);

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
      Core::Gen::Pairedvector<int, double>& nxmap_curr =
          static_cast<CONTACT::Node*>(smrtrnodes[i])->Data().GetDerivN()[0];
      Core::Gen::Pairedvector<int, double>& nymap_curr =
          static_cast<CONTACT::Node*>(smrtrnodes[i])->Data().GetDerivN()[1];

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

    Core::LinAlg::SerialDenseVector* mval = nullptr;
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

    Core::Gen::Pairedvector<int, double> dmap_sxib(nummnode * ndof + linsize);

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
      Core::Gen::Pairedvector<int, double>& nxmap_curr =
          static_cast<CONTACT::Node*>(smrtrnodes[i])->Data().GetDerivN()[0];
      Core::Gen::Pairedvector<int, double>& nymap_curr =
          static_cast<CONTACT::Node*>(smrtrnodes[i])->Data().GetDerivN()[1];

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
void CONTACT::Integrator::deriv_xi_g_p2_d(Mortar::Element& sele, Mortar::Element& mele,
    double sxigp, double mxigp, const Core::Gen::Pairedvector<int, double>& derivsxi,
    Core::Gen::Pairedvector<int, double>& derivmxi, int& linsize)
{
  // check for problem dimension
  if (Dim() != 2) FOUR_C_THROW("2D integration method called for non-2D problem");

  // we need the participating slave and master nodes
  Core::Nodes::Node** snodes = nullptr;
  Core::Nodes::Node** mnodes = nullptr;
  int numsnode = sele.num_node();
  int nummnode = mele.num_node();

  int ndof = Dim();

  snodes = sele.Nodes();
  mnodes = mele.Nodes();

  std::vector<Mortar::Node*> smrtrnodes(numsnode);
  std::vector<Mortar::Node*> mmrtrnodes(nummnode);

  for (int i = 0; i < numsnode; ++i)
  {
    smrtrnodes[i] = dynamic_cast<Mortar::Node*>(snodes[i]);
    if (!smrtrnodes[i]) FOUR_C_THROW("DerivXiAB2D: Null pointer!");
  }

  for (int i = 0; i < nummnode; ++i)
  {
    mmrtrnodes[i] = dynamic_cast<Mortar::Node*>(mnodes[i]);
    if (!mmrtrnodes[i]) FOUR_C_THROW("DerivXiAB2D: Null pointer!");
  }

  // we also need shape function derivs in A and B
  double psxigp[2] = {sxigp, 0.0};
  double pmxigp[2] = {mxigp, 0.0};
  Core::LinAlg::SerialDenseVector valsxigp(numsnode);
  Core::LinAlg::SerialDenseVector valmxigp(nummnode);
  Core::LinAlg::SerialDenseMatrix derivsxigp(numsnode, 1);
  Core::LinAlg::SerialDenseMatrix derivmxigp(nummnode, 1);

  sele.evaluate_shape(psxigp, valsxigp, derivsxigp, numsnode, false);
  mele.evaluate_shape(pmxigp, valmxigp, derivmxigp, nummnode, false);

  // we also need the GP slave coordinates + normal
  std::array<double, 3> sgpn = {0.0, 0.0, 0.0};
  std::array<double, 3> sgpx = {0.0, 0.0, 0.0};
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
  // When looking at Mortar::Projector::ProjectGaussPoint, one can see
  // that we do NOT use a unit normal there, either. Thus, why here?
  // First results suggest that it really makes no difference!

  // normalize interpolated GP normal back to length 1.0 !!!
  const double length = sqrt(sgpn[0] * sgpn[0] + sgpn[1] * sgpn[1] + sgpn[2] * sgpn[2]);
  if (length < 1.0e-12) FOUR_C_THROW("deriv_xi_g_p2_d: Divide by zero!");
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
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // build directional derivative of slave GP coordinates
  Core::Gen::Pairedvector<int, double> dmap_xsl_gp(linsize + nummnode * ndof);
  Core::Gen::Pairedvector<int, double> dmap_ysl_gp(linsize + nummnode * ndof);

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
  Core::Gen::Pairedvector<int, double> dmap_nxsl_gp(linsize + nummnode * ndof);
  Core::Gen::Pairedvector<int, double> dmap_nysl_gp(linsize + nummnode * ndof);

  std::array<double, 3> sgpnmod = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i) sgpnmod[i] = sgpn[i] * length;

  Core::Gen::Pairedvector<int, double> dmap_nxsl_gp_mod(linsize + nummnode * ndof);
  Core::Gen::Pairedvector<int, double> dmap_nysl_gp_mod(linsize + nummnode * ndof);

  for (int i = 0; i < numsnode; ++i)
  {
    Core::Gen::Pairedvector<int, double>& dmap_nxsl_i =
        static_cast<CONTACT::Node*>(smrtrnodes[i])->Data().GetDerivN()[0];
    Core::Gen::Pairedvector<int, double>& dmap_nysl_i =
        static_cast<CONTACT::Node*>(smrtrnodes[i])->Data().GetDerivN()[1];

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
void CONTACT::Integrator::DerivXiGP3D(Mortar::Element& sele, Mortar::Element& mele,
    const double* sxigp, const double* mxigp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& derivsxi,
    std::vector<Core::Gen::Pairedvector<int, double>>& derivmxi, const double alpha)
{
  // check for problem dimension
  if (Dim() != 3) FOUR_C_THROW("3D integration method called for non-3D problem");

  // we need the participating slave and master nodes
  Core::Nodes::Node** snodes = sele.Nodes();
  Core::Nodes::Node** mnodes = mele.Nodes();
  std::vector<Mortar::Node*> smrtrnodes(sele.num_node());
  std::vector<Mortar::Node*> mmrtrnodes(mele.num_node());
  const int numsnode = sele.num_node();
  const int nummnode = mele.num_node();

  for (int i = 0; i < numsnode; ++i)
  {
    smrtrnodes[i] = dynamic_cast<Mortar::Node*>(snodes[i]);
    if (!smrtrnodes[i]) FOUR_C_THROW("DerivXiGP3D: Null pointer!");
  }

  for (int i = 0; i < nummnode; ++i)
  {
    mmrtrnodes[i] = dynamic_cast<Mortar::Node*>(mnodes[i]);
    if (!mmrtrnodes[i]) FOUR_C_THROW("DerivXiGP3D: Null pointer!");
  }

  // we also need shape function derivs at the GP
  Core::LinAlg::SerialDenseVector valsxigp(numsnode);
  Core::LinAlg::SerialDenseVector valmxigp(nummnode);
  Core::LinAlg::SerialDenseMatrix derivsxigp(numsnode, 2, true);
  Core::LinAlg::SerialDenseMatrix derivmxigp(nummnode, 2, true);

  sele.evaluate_shape(sxigp, valsxigp, derivsxigp, numsnode);
  mele.evaluate_shape(mxigp, valmxigp, derivmxigp, nummnode);

  // we also need the GP slave coordinates + normal
  std::array<double, 3> sgpn = {0.0, 0.0, 0.0};
  std::array<double, 3> sgpx = {0.0, 0.0, 0.0};
  for (int i = 0; i < numsnode; ++i)
    for (int k = 0; k < 3; ++k)
    {
      sgpn[k] += valsxigp[i] * smrtrnodes[i]->MoData().n()[k];
      sgpx[k] += valsxigp[i] * smrtrnodes[i]->xspatial()[k];
    }

  // build 3x3 factor matrix L
  Core::LinAlg::Matrix<3, 3> lmatrix(true);
  for (int k = 0; k < 3; ++k) lmatrix(k, 2) = -sgpn[k];
  for (int z = 0; z < nummnode; ++z)
    for (int k = 0; k < 3; ++k)
    {
      lmatrix(k, 0) += derivmxigp(z, 0) * mmrtrnodes[z]->xspatial()[k];
      lmatrix(k, 1) += derivmxigp(z, 1) * mmrtrnodes[z]->xspatial()[k];
    }

  // get inverse of the 3x3 matrix L (in place)
  if (abs(lmatrix.Determinant()) < 1e-12) FOUR_C_THROW("Singular lmatrix for derivgp3d");

  lmatrix.Invert();

  // build directional derivative of slave GP normal
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  int linsize = 0;
  for (int i = 0; i < numsnode; ++i)
  {
    Node* cnode = dynamic_cast<Node*>(snodes[i]);
    linsize += cnode->GetLinsize();
  }

  Core::Gen::Pairedvector<int, double> dmap_nxsl_gp(linsize);
  Core::Gen::Pairedvector<int, double> dmap_nysl_gp(linsize);
  Core::Gen::Pairedvector<int, double> dmap_nzsl_gp(linsize);

  for (int i = 0; i < numsnode; ++i)
  {
    Core::Gen::Pairedvector<int, double>& dmap_nxsl_i =
        static_cast<CONTACT::Node*>(smrtrnodes[i])->Data().GetDerivN()[0];
    Core::Gen::Pairedvector<int, double>& dmap_nysl_i =
        static_cast<CONTACT::Node*>(smrtrnodes[i])->Data().GetDerivN()[1];
    Core::Gen::Pairedvector<int, double>& dmap_nzsl_i =
        static_cast<CONTACT::Node*>(smrtrnodes[i])->Data().GetDerivN()[2];

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
void CONTACT::Integrator::DerivXiGP3DAuxPlane(Mortar::Element& ele, double* xigp, double* auxn,
    std::vector<Core::Gen::Pairedvector<int, double>>& derivxi, double& alpha,
    std::vector<Core::Gen::Pairedvector<int, double>>& derivauxn,
    Core::Gen::Pairedvector<int, Core::LinAlg::Matrix<3, 1>>& derivgp)
{
  // check for problem dimension
  if (Dim() != 3) FOUR_C_THROW("3D integration method called for non-3D problem");

  // we need the participating element nodes
  Core::Nodes::Node** nodes = ele.Nodes();
  std::vector<Mortar::Node*> mrtrnodes(ele.num_node());
  const int numnode = ele.num_node();

  for (int i = 0; i < numnode; ++i)
  {
    mrtrnodes[i] = dynamic_cast<Mortar::Node*>(nodes[i]);
    if (!mrtrnodes[i]) FOUR_C_THROW("DerivXiGP3DAuxPlane: Null pointer!");
  }

  // we also need shape function derivs at the GP
  Core::LinAlg::SerialDenseVector valxigp(numnode);
  Core::LinAlg::SerialDenseMatrix derivxigp(numnode, 2, true);
  ele.evaluate_shape(xigp, valxigp, derivxigp, numnode);

  // build 3x3 factor matrix L
  Core::LinAlg::Matrix<3, 3> lmatrix(true);
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
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // see if this is an IntEle
  Mortar::IntElement* ie = dynamic_cast<Mortar::IntElement*>(&ele);

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
    std::vector<std::vector<Core::Gen::Pairedvector<int, double>>> nodelin(0);
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
  for (Core::Gen::Pairedvector<int, Core::LinAlg::Matrix<3, 1>>::const_iterator p = derivgp.begin();
       p != derivgp.end(); ++p)
  {
    const int pf = p->first;
    double& derivxi0 = derivxi[0][pf];
    double& derivxi1 = derivxi[1][pf];
    const Core::LinAlg::Matrix<3, 1>& tmp = p->second;
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
void CONTACT::Integrator::integrate_gp_3_d(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseMatrix& sderiv,
    Core::LinAlg::SerialDenseMatrix& mderiv, Core::LinAlg::SerialDenseMatrix& lmderiv,
    Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap, double& wgt,
    double& jac, Core::Gen::Pairedvector<int, double>& derivjac, double* normal,
    std::vector<Core::Gen::Pairedvector<int, double>>& dnmap_unit, double& gap,
    Core::Gen::Pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
    std::vector<Core::Gen::Pairedvector<int, double>>& derivsxi,
    std::vector<Core::Gen::Pairedvector<int, double>>& derivmxi)
{
  // different algorithms integrate different sh*t
  switch (algo_)
  {
    case Inpar::Mortar::algorithm_mortar:
    {
      // is quadratic case?
      bool quad = sele.IsQuad();

      // weighted gap
      gp_3_d_w_gap(sele, sval, lmval, &gap, jac, wgt, quad);
      for (int j = 0; j < sele.num_node(); ++j)
        gp_g_lin(j, sele, mele, sval, mval, lmval, sderiv, lmderiv, gap, normal, jac, wgt,
            deriv_gap, derivjac, derivsxi, dualmap);

      // check bound
      bool bound = false;
      for (int i = 0; i < sele.num_node(); ++i)
        if (dynamic_cast<CONTACT::Node*>(sele.Nodes()[i])->IsOnBoundorCE()) bound = true;

      // integrate D and M matrix
      gp_dm(sele, mele, lmval, sval, mval, jac, wgt, bound);

      // compute segment D/M linearization
      if (bound)
      {
        gp_3_d_dm_lin_bound(sele, mele, sval, mval, lmval, sderiv, lmderiv, mderiv, jac, wgt,
            derivjac, derivsxi, derivmxi, dualmap);
      }
      else
      {
        gp_3_d_dm_lin(sele, mele, sval, mval, lmval, sderiv, mderiv, lmderiv, wgt, jac, derivsxi,
            derivmxi, derivjac, dualmap);
      }

      // Creating the WEIGHTED tangential relative slip increment (non-objective)
      if (gpslip_)
      {
        int linsize = 0;
        for (int i = 0; i < sele.num_node(); ++i)
          linsize += dynamic_cast<Node*>(sele.Nodes()[i])->GetLinsize();

        double jumpvalv[2] = {0.0, 0.0};  // jump for slipincr --> equal to jumpval
        std::vector<Core::Gen::Pairedvector<int, double>> dslipgp(
            2, ((mele.num_node() * Dim()) + linsize));  // deriv. of slip for slipincr (xi, eta)

        // Creating the WEIGHTED tangential relative slip increment (non-objective)
        gp_3_d_slip_incr(sele, mele, sval, mval, lmval, sderiv, mderiv, jac, wgt, jumpvalv,
            derivsxi, derivmxi, dslipgp);
        // Lin weighted slip
        for (int iter = 0; iter < sele.num_node(); ++iter)
          gp_3_d_slip_incr_lin(iter, sele, sval, lmval, sderiv, lmderiv, jac, wgt, jumpvalv,
              derivjac, dslipgp, derivsxi, dualmap);
      }

      //*******************************
      // WEAR stuff
      //*******************************
      // std. wear for all wear-algorithm types
      if (wearlaw_ != Inpar::Wear::wear_none)
      {
        int linsize = 0;
        for (int i = 0; i < sele.num_node(); ++i)
          linsize += dynamic_cast<Node*>(sele.Nodes()[i])->GetLinsize();

        double jumpval[2] = {0.0, 0.0};  // jump for wear
        double wearval = 0.0;            // wear value
        Core::Gen::Pairedvector<int, double> dsliptmatrixgp(
            (mele.num_node() * Dim()) + linsize);  // deriv. of slip for wear
        Core::Gen::Pairedvector<int, double> dweargp(
            (mele.num_node() * Dim()) + linsize);  // wear lin. without lm and jac.
        Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> lagmult =
            Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(3, sele.num_node()));
        sele.GetNodalLagMult(*lagmult);

        gp_3_d_wear(sele, mele, sval, sderiv, mval, mderiv, lmval, lmderiv, lagmult, normal, jac,
            wgt, jumpval, &wearval, dsliptmatrixgp, dweargp, derivsxi, derivmxi, dnmap_unit,
            dualmap);

        if (wear_type() == Inpar::Wear::wear_intstate and wearimpl_ == true)
          for (int j = 0; j < sele.num_node(); ++j)
            gp_3_d_wear_lin(j, sele, sval, lmval, sderiv, lmderiv, jac, normal, wgt, wearval,
                jumpval, dweargp, derivjac, derivsxi, dualmap);

        // integrate T and E matrix for discr. wear
        if (wear_type() == Inpar::Wear::wear_primvar) gp_te(sele, lmval, sval, jac, wgt, jumpval);

        // both-sided discr wear specific stuff
        if (wear_side() == Inpar::Wear::wear_both and wear_type() == Inpar::Wear::wear_primvar)
        {
          Core::LinAlg::SerialDenseVector lm2val(mele.num_node());
          Core::LinAlg::SerialDenseMatrix lm2deriv(mele.num_node(), Dim() - 1, true);
          mele.evaluate_shape_lag_mult(shape_fcn(), mxi, lm2val, lm2deriv, mele.num_node());

          gp_te_master(sele, mele, lmval, lm2val, mval, jac, wgt, jumpval, Comm_);
        }

        // both-sided wear specific stuff
        if (wear_side() == Inpar::Wear::wear_both and wear_type() == Inpar::Wear::wear_intstate)
        {
          Core::LinAlg::SerialDenseVector lm2val(mele.num_node());
          Core::LinAlg::SerialDenseMatrix lm2deriv(mele.num_node(), 2, true);
          mele.evaluate_shape_lag_mult(shape_fcn(), mxi, lm2val, lm2deriv, mele.num_node());

          gp_d2(sele, mele, lm2val, mval, jac, wgt, Comm_);
        }

        // Lin wear matrices T and E for discr. wear
        if (wearimpl_ == true and wear_type() == Inpar::Wear::wear_primvar)
          for (int j = 0; j < sele.num_node(); ++j)
            gp_3_d_te_lin(j, sele, sval, lmval, sderiv, lmderiv, jac, wgt, jumpval, derivsxi,
                derivjac, dsliptmatrixgp, dualmap);



        if (wear_side() == Inpar::Wear::wear_both and wear_type() == Inpar::Wear::wear_primvar and
            wearimpl_ == true)
        {
          Core::LinAlg::SerialDenseVector lm2val(mele.num_node());
          Core::LinAlg::SerialDenseMatrix lm2deriv(mele.num_node(), Dim() - 1, true);
          Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix> dual2map(
              (sele.num_node() + mele.num_node()) * Dim(), 0,
              Core::LinAlg::SerialDenseMatrix(mele.num_node(), mele.num_node()));
          mele.evaluate_shape_lag_mult(shape_fcn(), mxi, lm2val, lm2deriv, mele.num_node());
          mele.DerivShapeDual(dual2map);

          for (int iter = 0; iter < mele.num_node(); ++iter)
            gp_3_d_te_master_lin(iter, sele, mele, sval, mval, lmval, lm2val, sderiv, mderiv,
                lmderiv, lm2deriv, jac, wgt, jumpval, derivsxi, derivmxi, derivjac, dsliptmatrixgp,
                dualmap, dual2map, Comm_);
        }
      }

      //*******************************
      // PORO stuff
      //*******************************
      bool poroprob = false;
      if (imortar_.get<int>("PROBTYPE") == Inpar::CONTACT::poroelast ||
          imortar_.get<int>("PROBTYPE") == Inpar::CONTACT::poroscatra)
        if (imortar_.get<bool>("CONTACTNOPEN"))  // evaluate additional terms just in case of no
                                                 // penectration condition
          poroprob = true;
      if (poroprob)
      {
        double ncoup = 0.0;
        std::map<int, double> dncoupgp;      // ncoup lin. without lm and jac.
        std::map<int, double> dvelncoupgp;   // velocity ncoup lin. without lm and jac.
        std::map<int, double> dpresncoupgp;  // pressure ncoup lin. without lm and jac.

        gp_ncoup_deriv(sele, mele, sval, mval, lmval, sderiv, mderiv, &ncoup, normal, jac, wgt, sxi,
            derivsxi, derivmxi, dncoupgp, dvelncoupgp, dpresncoupgp, dnmap_unit, false);
        // Lin ncoup condition
        for (int j = 0; j < sele.num_node(); ++j)
          gp_ncoup_lin(j, sele, mele, sval, mval, lmval, sderiv, lmderiv, ncoup, normal, jac, wgt,
              dncoupgp, dvelncoupgp, dpresncoupgp, derivjac, derivsxi, derivmxi, dualmap);
      }
      break;
    }

    default:
      FOUR_C_THROW("contact integrator doesn't know what to do with this algorithm");
  };

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::integrate_gp_2_d(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseMatrix& sderiv,
    Core::LinAlg::SerialDenseMatrix& mderiv, Core::LinAlg::SerialDenseMatrix& lmderiv,
    Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap, double& wgt,
    double& jac, Core::Gen::Pairedvector<int, double>& derivjac, double* normal,
    std::vector<Core::Gen::Pairedvector<int, double>>& dnmap_unit, double& gap,
    Core::Gen::Pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
    std::vector<Core::Gen::Pairedvector<int, double>>& derivsxi,
    std::vector<Core::Gen::Pairedvector<int, double>>& derivmxi)
{
  // different algorithms integrate different sh*t
  switch (algo_)
  {
    case Inpar::Mortar::algorithm_mortar:
    {
      // decide whether boundary modification has to be considered or not
      // this is element-specific (is there a boundary node in this element?)
      bool bound = false;
      for (int k = 0; k < sele.num_node(); ++k)
      {
        Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(sele.Nodes()[k]);
        if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_segment2_d: Null pointer!");

        if (mymrtrnode->IsOnBoundorCE())
        {
          bound = true;
          break;
        }
      }

      // decide whether linear LM are used for quadratic FE here
      bool linlm = false;
      if (lag_mult_quad() == Inpar::Mortar::lagmult_lin &&
          sele.Shape() == Core::FE::CellType::line3)
      {
        bound = false;  // crosspoints and linear LM NOT at the same time!!!!
        linlm = true;
      }

      // weighted gap
      gp_2_d_w_gap(sele, sval, lmval, &gap, jac, wgt);
      for (int j = 0; j < sval.numRows(); ++j)
        gp_2_d_g_lin(j, sele, mele, sval, mval, lmval, sderiv, lmderiv, gap, normal, jac, wgt,
            deriv_gap, derivjac, derivsxi, dualmap);

      // integrate D and M matrix
      gp_dm(sele, mele, lmval, sval, mval, jac, wgt, bound);

      // compute segment D/M linearization  -- bound
      if (bound)
      {
        gp_2_d_dm_lin_bound(sele, mele, sval, mval, lmval, sderiv, mderiv, lmderiv, jac, wgt,
            derivjac, derivsxi, derivmxi, dualmap);
      }
      else
      {
        // D/M linearization
        for (int j = 0; j < sval.numRows(); ++j)
        {
          gp_2_d_dm_lin(j, bound, linlm, sele, mele, sval, mval, lmval, sderiv, mderiv, lmderiv,
              jac, wgt, derivsxi, derivmxi, derivjac, dualmap);
        }
      }

      // Creating the tangential relative slip increment (non-objective)
      if (Core::UTILS::IntegralValue<int>(imortar_, "GP_SLIP_INCR") == true)
      {
        int linsize = 0;
        for (int i = 0; i < sele.num_node(); ++i)
          linsize += dynamic_cast<Node*>(sele.Nodes()[i])->GetLinsize();
        linsize = linsize * 2;

        double jumpvalv = 0.0;  // jump for slipincr --> equal to jumpval
        Core::Gen::Pairedvector<int, double> dslipgp(
            linsize + Dim() * mele.num_node());  // deriv. of slip for slipincr

        gp_2_d_slip_incr(sele, mele, sval, mval, lmval, sderiv, mderiv, jac, wgt, &jumpvalv,
            derivsxi, derivmxi, dslipgp, linsize);

        for (int iter = 0; iter < sele.num_node(); ++iter)
          gp_2_d_slip_incr_lin(iter, sele, sval, lmval, sderiv, lmderiv, jac, wgt, &jumpvalv,
              derivsxi, dslipgp, derivjac, dualmap);
      }

      // wear specific stuff
      if (wearlaw_ != Inpar::Wear::wear_none)
      {
        // nodal Lagrange mulitplier
        Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> lagmult =
            Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(3, sele.num_node()));
        sele.GetNodalLagMult(*lagmult);

        int linsize = 0;
        for (int i = 0; i < sele.num_node(); ++i)
          linsize += dynamic_cast<Node*>(sele.Nodes()[i])->GetLinsize();
        linsize = linsize * 2;

        double jumpval = 0.0;  // jump for wear
        double wearval = 0.0;  // wear value
        Core::Gen::Pairedvector<int, double> dsliptmatrixgp(
            linsize + Dim() * mele.num_node());  // deriv. of slip for wear
        Core::Gen::Pairedvector<int, double> dweargp(
            linsize + Dim() * mele.num_node());  // wear lin without weighting and jac

        // std. wear for all wear-algorithm types
        gp_2_d_wear(sele, mele, sval, sderiv, mval, mderiv, lmval, lmderiv, lagmult, normal, jac,
            wgt, &jumpval, &wearval, dsliptmatrixgp, dweargp, derivsxi, derivmxi, dnmap_unit,
            dualmap);

        // integrate T and E matrix for discr. wear
        if (wear_type() == Inpar::Wear::wear_primvar) gp_te(sele, lmval, sval, jac, wgt, &jumpval);

        // both-sided discr wear specific stuff
        if (wear_side() == Inpar::Wear::wear_both and wear_type() == Inpar::Wear::wear_primvar)
        {
          Core::LinAlg::SerialDenseVector lm2val(mele.num_node());
          Core::LinAlg::SerialDenseMatrix lm2deriv(mele.num_node(), Dim() - 1, true);
          mele.evaluate_shape_lag_mult(shape_fcn(), mxi, lm2val, lm2deriv, mele.num_node());

          gp_te_master(sele, mele, lmval, lm2val, mval, jac, wgt, &jumpval, Comm_);
        }

        // both-sided map wear specific stuff
        if (wear_side() == Inpar::Wear::wear_both and wear_type() == Inpar::Wear::wear_intstate)
        {
          Core::LinAlg::SerialDenseVector lm2val(mele.num_node());
          Core::LinAlg::SerialDenseMatrix lm2deriv(mele.num_node(), 2, true);
          mele.evaluate_shape_lag_mult(shape_fcn(), mxi, lm2val, lm2deriv, mele.num_node());

          gp_d2(sele, mele, lm2val, mval, jac, wgt, Comm_);
        }

        // Lin wear for impl. alg.
        if (wearimpl_ == true and wear_type() == Inpar::Wear::wear_intstate)
          for (int j = 0; j < sele.num_node(); ++j)
            gp_2_d_wear_lin(j, sele, sval, lmval, sderiv, lmderiv, jac, normal, wgt, wearval,
                &jumpval, dweargp, derivjac, derivsxi, dualmap);

        // Lin wear T and E matrix
        if (wearimpl_ == true and wear_type() == Inpar::Wear::wear_primvar)
          for (int j = 0; j < sele.num_node(); ++j)
            gp_2_d_te_lin(j, sele, sval, lmval, sderiv, lmderiv, jac, wgt, &jumpval, derivsxi,
                derivjac, dsliptmatrixgp, dualmap);

        if (wear_side() == Inpar::Wear::wear_both and wear_type() == Inpar::Wear::wear_primvar and
            wearimpl_ == true)
        {
          Core::LinAlg::SerialDenseVector lm2val(mele.num_node());
          Core::LinAlg::SerialDenseMatrix lm2deriv(mele.num_node(), Dim() - 1, true);
          Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix> dual2map(
              (sele.num_node() + mele.num_node()) * Dim(), 0,
              Core::LinAlg::SerialDenseMatrix(mele.num_node(), mele.num_node()));
          mele.evaluate_shape_lag_mult(shape_fcn(), mxi, lm2val, lm2deriv, mele.num_node());
          mele.DerivShapeDual(dual2map);

          for (int iter = 0; iter < mele.num_node(); ++iter)
            gp_3_d_te_master_lin(iter, sele, mele, sval, mval, lmval, lm2val, sderiv, mderiv,
                lmderiv, lm2deriv, jac, wgt, &jumpval, derivsxi, derivmxi, derivjac, dsliptmatrixgp,
                dualmap, dual2map, Comm_);
        }
      }  // if wear

      //*******************************
      // PORO stuff
      //*******************************
      bool poroprob = false;
      if (imortar_.get<int>("PROBTYPE") == Inpar::CONTACT::poroelast ||
          imortar_.get<int>("PROBTYPE") == Inpar::CONTACT::poroscatra)
        if (imortar_.get<bool>("CONTACTNOPEN"))  // evaluate additional terms just in case of no
                                                 // penectration condition
          poroprob = true;
      if (poroprob)
      {
        double ncoup = 0.0;
        std::map<int, double> dncoupgp;      // ncoup lin. without lm and jac.
        std::map<int, double> dvelncoupgp;   // velocity ncoup lin. without lm and jac.
        std::map<int, double> dpresncoupgp;  // pressure ncoup lin. without lm and jac.

        gp_ncoup_deriv(sele, mele, sval, mval, lmval, sderiv, mderiv, &ncoup, normal, jac, wgt, sxi,
            derivsxi, derivmxi, dncoupgp, dvelncoupgp, dpresncoupgp, dnmap_unit, false);
        // Lin ncoup condition
        for (int j = 0; j < sele.num_node(); ++j)
          gp_ncoup_lin(j, sele, mele, sval, mval, lmval, sderiv, lmderiv, ncoup, normal, jac, wgt,
              dncoupgp, dvelncoupgp, dpresncoupgp, derivjac, derivsxi, derivmxi, dualmap);
      }

      break;
    }
    default:
      FOUR_C_THROW("contact integrator doesn't know what to do with this algorithm");
  };

  return;
}


/*----------------------------------------------------------------------*
 |  Compute entries for D and M matrix at GP                 farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::gp_dm(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& lmval, Core::LinAlg::SerialDenseVector& sval,
    Core::LinAlg::SerialDenseVector& mval, double& jac, double& wgt, bool& bound)
{
  int nrow;
  int ncol;
  // get slave element nodes themselves
  Core::Nodes::Node** snodes = nullptr;
  Core::Nodes::Node** mnodes = nullptr;
  nrow = sele.num_node();
  snodes = sele.Nodes();
  ncol = mele.num_node();
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
  if ((shape_fcn() == Inpar::Mortar::shape_dual ||
          shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) &&
      lag_mult_quad() != Inpar::Mortar::lagmult_lin)
  {
    for (int j = 0; j < nrow; ++j)
    {
      CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(snodes[j]);

      if (cnode->IsOnBoundorCE()) continue;

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(mnodes[k]);

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
          CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(snodes[k]);

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
      CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(snodes[j]);

      if ((shapefcn_ == Inpar::Mortar::shape_standard && cnode->IsOnBoundorCE()) ||
          ((shapefcn_ == Inpar::Mortar::shape_dual ||
               shapefcn_ == Inpar::Mortar::shape_petrovgalerkin) &&
              cnode->IsOnBound()))
        continue;

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(mnodes[k]);

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
        CONTACT::Node* snode = dynamic_cast<CONTACT::Node*>(snodes[k]);

        // multiply the two shape functions
        double prod = lmval[j] * sval[k] * jac * wgt;
        if (abs(prod) > MORTARINTTOL)
        {
          if ((shape_fcn() == Inpar::Mortar::shape_standard && snode->IsOnBoundorCE()) ||
              ((shape_fcn() == Inpar::Mortar::shape_dual ||
                   shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) &&
                  snode->IsOnBound()))
          {
            if (shape_fcn() == Inpar::Mortar::shape_standard && nonsmooth_)
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
            else if ((shape_fcn() == Inpar::Mortar::shape_dual ||
                         shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) &&
                     lag_mult_quad() == Inpar::Mortar::lagmult_lin && nonsmooth_)
              FOUR_C_THROW(
                  "dual shape functions with locally linear Lagrange multipliers in "
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
void inline CONTACT::Integrator::gp_3_d_dm_quad(Mortar::Element& sele, Mortar::Element& mele,
    Mortar::IntElement& sintele, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseVector& lmintval, Core::LinAlg::SerialDenseVector& sval,
    Core::LinAlg::SerialDenseVector& mval, const double& jac, double& wgt, const int& nrow,
    const int& nintrow, const int& ncol, const int& ndof, bool& bound)
{
  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();
  if (!snodes) FOUR_C_THROW("Null pointer!");
  Core::Nodes::Node** mnodes = mele.Nodes();
  if (!mnodes) FOUR_C_THROW("Null pointer!");
  Core::Nodes::Node** sintnodes = sintele.Nodes();
  if (!sintnodes) FOUR_C_THROW("Null pointer for sintnodes!");

  // CASE 1/2: Standard LM shape functions and quadratic, linear or constant interpolation
  // CASE 5: dual LM shape functions and linear interpolation
  if ((shape_fcn() == Inpar::Mortar::shape_standard &&
          (lag_mult_quad() == Inpar::Mortar::lagmult_quad ||
              lag_mult_quad() == Inpar::Mortar::lagmult_lin ||
              lag_mult_quad() == Inpar::Mortar::lagmult_const)) ||
      ((shape_fcn() == Inpar::Mortar::shape_dual ||
           shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) &&
          lag_mult_quad() == Inpar::Mortar::lagmult_lin))
  {
    // compute all mseg and dseg matrix entries
    // loop over Lagrange multiplier dofs j
    for (int j = 0; j < nrow; ++j)
    {
      Node* cnode = dynamic_cast<Node*>(snodes[j]);

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        Node* mnode = dynamic_cast<Node*>(mnodes[k]);

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
        Node* snode = dynamic_cast<Node*>(snodes[k]);

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
  else if (shape_fcn() == Inpar::Mortar::shape_standard &&
           lag_mult_quad() == Inpar::Mortar::lagmult_pwlin)
  {
    // compute all mseg and dseg matrix entries
    // loop over Lagrange multiplier dofs j
    for (int j = 0; j < nintrow; ++j)
    {
      Node* cnode = dynamic_cast<Node*>(sintnodes[j]);

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        Node* mnode = dynamic_cast<Node*>(mnodes[k]);

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
        Node* snode = dynamic_cast<Node*>(snodes[k]);

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
  else if ((shape_fcn() == Inpar::Mortar::shape_dual ||
               shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) &&
           (lag_mult_quad() == Inpar::Mortar::lagmult_quad ||
               lag_mult_quad() == Inpar::Mortar::lagmult_const))
  {
    // compute all mseg and dseg matrix entries
    // loop over Lagrange multiplier dofs j
    for (int j = 0; j < nrow; ++j)
    {
      Node* cnode = dynamic_cast<Node*>(snodes[j]);

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        Node* mnode = dynamic_cast<Node*>(mnodes[k]);

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
    FOUR_C_THROW("Invalid integration case for 3D quadratic mortar!");
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Compute entries for weighted Gap at GP                   farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::Integrator::gp_2_d_w_gap(Mortar::Element& sele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& lmval, double* gap,
    double& jac, double& wgt)
{
  Core::Nodes::Node** mynodes = nullptr;
  int nrow = 0;
  nrow = sele.num_node();
  mynodes = sele.Nodes();

  // **************************
  // add to node
  // **************************
  for (int j = 0; j < nrow; ++j)
  {
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(mynodes[j]);

    double prod = 0.0;
    // Petrov-Galerkin approach (dual LM for D/M but standard LM for gap)
    if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) prod = sval[j] * gap[0] * jac * wgt;
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
void CONTACT::Integrator::gp_3_d_w_gap(Mortar::Element& sele, Core::LinAlg::SerialDenseVector& sval,
    Core::LinAlg::SerialDenseVector& lmval, double* gap, double& jac, double& wgt, bool quadratic,
    int nintrow)
{
  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();
  if (!snodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  // number of nodes (slave, master)
  const int nrow = sele.num_node();

  // **************************
  // add to node
  // **************************
  if (!quadratic)
  {
    for (int j = 0; j < nrow; ++j)
    {
      CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(snodes[j]);

      double prod = 0.0;
      // Petrov-Galerkin approach (dual LM for D/M but standard LM for gap)
      if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) prod = sval[j] * gap[0] * jac * wgt;
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
    if (shape_fcn() == Inpar::Mortar::shape_standard &&
        (lag_mult_quad() == Inpar::Mortar::lagmult_quad ||
            lag_mult_quad() == Inpar::Mortar::lagmult_lin ||
            lag_mult_quad() == Inpar::Mortar::lagmult_const))
    {
      for (int j = 0; j < nrow; ++j)
      {
        CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(snodes[j]);

        double prod = 0.0;
        prod = lmval[j] * gap[0] * jac * wgt;

        if (cnode->IsOnBoundorCE()) continue;

        // add current Gauss point's contribution to gseg
        cnode->AddgValue(prod);
      }
    }

    // CASE 3: Standard LM shape functions and piecewise linear interpolation
    // Attention:  for this case, lmval represents lmintval !!!
    else if (shape_fcn() == Inpar::Mortar::shape_standard &&
             lag_mult_quad() == Inpar::Mortar::lagmult_pwlin)
    {
      if (nintrow == 0) FOUR_C_THROW("ERROR!");

      for (int j = 0; j < nintrow; ++j)
      {
        CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(snodes[j]);

        double prod = 0.0;
        prod = lmval[j] * gap[0] * jac * wgt;

        if (cnode->IsOnBound()) continue;

        // add current Gauss point's contribution to gseg
        cnode->AddgValue(prod);
      }
    }

    // CASE 4: Dual LM shape functions and quadratic, linear or constant interpolation
    else if ((shape_fcn() == Inpar::Mortar::shape_dual ||
                 shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) &&
             (lag_mult_quad() == Inpar::Mortar::lagmult_quad ||
                 lag_mult_quad() == Inpar::Mortar::lagmult_lin ||
                 lag_mult_quad() == Inpar::Mortar::lagmult_const))
    {
      for (int j = 0; j < nrow; ++j)
      {
        CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(snodes[j]);

        if (cnode->IsOnBound()) continue;

        double prod = 0.0;
        // Petrov-Galerkin approach (dual LM for D/M but standard LM for gap)
        if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) prod = sval[j] * gap[0] * jac * wgt;
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
      FOUR_C_THROW("Invalid integration case for 3D quadratic contact!");
    }
  }

  return;
}

void CONTACT::Integrator::gap_3_d(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& mval,
    Core::LinAlg::SerialDenseMatrix& sderiv, Core::LinAlg::SerialDenseMatrix& mderiv, double* gap,
    double* gpn, std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    std::vector<Core::Gen::Pairedvector<int, double>>& dmxigp,
    Core::Gen::Pairedvector<int, double>& dgapgp,
    std::vector<Core::Gen::Pairedvector<int, double>>& dnmap_unit)
{
  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();
  Core::Nodes::Node** mnodes = mele.Nodes();
  if (!snodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");
  if (!mnodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  // number of nodes (slave, master)
  const int nrow = sele.num_node();
  const int ncol = mele.num_node();

  std::array<double, 3> sgpx = {0.0, 0.0, 0.0};
  std::array<double, 3> mgpx = {0.0, 0.0, 0.0};

  for (int i = 0; i < nrow; ++i)
  {
    Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(snodes[i]);
    gpn[0] += sval[i] * mymrtrnode->MoData().n()[0];
    gpn[1] += sval[i] * mymrtrnode->MoData().n()[1];
    gpn[2] += sval[i] * mymrtrnode->MoData().n()[2];

    if (wear_type() == Inpar::Wear::wear_primvar)
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
    if (wear_side() == Inpar::Wear::wear_both and wear_type() == Inpar::Wear::wear_primvar)
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
  if (lengthn < 1.0e-12) FOUR_C_THROW("IntegrateAndDerivSegment: Divide by zero!");

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
      Node* cnode = dynamic_cast<Node*>(mnodes[i]);
      linsize += cnode->GetLinsize();
    }
    linsize += 1 + mele.num_node();
  }
  else
  {
    for (int i = 0; i < nrow; ++i)
    {
      Node* cnode = dynamic_cast<Node*>(snodes[i]);
      linsize += cnode->GetLinsize();
    }
  }
  linsize *= 10;

  // build directional derivative of slave GP normal (non-unit)
  Core::Gen::Pairedvector<int, double> dmap_nxsl_gp(linsize);
  Core::Gen::Pairedvector<int, double> dmap_nysl_gp(linsize);
  Core::Gen::Pairedvector<int, double> dmap_nzsl_gp(linsize);

  for (int i = 0; i < nrow; ++i)
  {
    Node* cnode = dynamic_cast<Node*>(snodes[i]);

    Core::Gen::Pairedvector<int, double>& dmap_nxsl_i = cnode->Data().GetDerivN()[0];
    Core::Gen::Pairedvector<int, double>& dmap_nysl_i = cnode->Data().GetDerivN()[1];
    Core::Gen::Pairedvector<int, double>& dmap_nzsl_i = cnode->Data().GetDerivN()[2];

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
  if (wear_type() == Inpar::Wear::wear_primvar)
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

        for (_CI p = frinode->Data().GetDerivN()[k].begin();
             p != frinode->Data().GetDerivN()[k].end(); ++p)
          dgapgp[p->first] += gpn[k] * sval[z] * frinode->WearData().wcurr()[0] * (p->second);
      }
    }
  }
  else
  {
    for (int z = 0; z < nrow; ++z)
    {
      Node* cnode = dynamic_cast<Node*>(snodes[z]);
      for (int k = 0; k < 3; ++k) dgapgp[cnode->Dofs()[k]] -= sval[z] * gpn[k];
    }

    for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < nrow; ++z)
      {
        Node* cnode = dynamic_cast<Node*>(snodes[z]);
        for (int k = 0; k < 3; ++k) dg -= gpn[k] * sderiv(z, 0) * cnode->xspatial()[k] * ps;
      }
    }

    for (_CI p = dsxigp[1].begin(); p != dsxigp[1].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < nrow; ++z)
      {
        Node* cnode = dynamic_cast<Node*>(snodes[z]);
        for (int k = 0; k < 3; ++k) dg -= gpn[k] * sderiv(z, 1) * cnode->xspatial()[k] * ps;
      }
    }
  }

  //        MASTER
  if (wear_side() == Inpar::Wear::wear_both and wear_type() == Inpar::Wear::wear_primvar)
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

        for (_CI p = frinode->Data().GetDerivN()[k].begin();
             p != frinode->Data().GetDerivN()[k].end(); ++p)
          dgapgp[p->first] -= gpn[k] * mval[z] * frinode->WearData().wcurr()[0] * (p->second);
      }
    }
  }
  else
  {
    // lin master nodes
    for (int z = 0; z < ncol; ++z)
    {
      Node* cnode = dynamic_cast<Node*>(mnodes[z]);
      for (int k = 0; k < 3; ++k) dgapgp[cnode->Dofs()[k]] += mval[z] * gpn[k];
    }

    for (_CI p = dmxigp[0].begin(); p != dmxigp[0].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < ncol; ++z)
      {
        Node* cnode = dynamic_cast<Node*>(mnodes[z]);
        for (int k = 0; k < 3; ++k) dg += gpn[k] * mderiv(z, 0) * cnode->xspatial()[k] * ps;
      }
    }

    for (_CI p = dmxigp[1].begin(); p != dmxigp[1].end(); ++p)
    {
      double& dg = dgapgp[p->first];
      const double& ps = p->second;
      for (int z = 0; z < ncol; ++z)
      {
        Node* cnode = dynamic_cast<Node*>(mnodes[z]);
        for (int k = 0; k < 3; ++k) dg += gpn[k] * mderiv(z, 1) * cnode->xspatial()[k] * ps;
      }
    }
  }

  return;
}

void CONTACT::Integrator::gap_2_d(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& mval,
    Core::LinAlg::SerialDenseMatrix& sderiv, Core::LinAlg::SerialDenseMatrix& mderiv, double* gap,
    double* gpn, std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    std::vector<Core::Gen::Pairedvector<int, double>>& dmxigp,
    Core::Gen::Pairedvector<int, double>& dgapgp,
    std::vector<Core::Gen::Pairedvector<int, double>>& dnmap_unit)
{
  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // get slave element nodes themselves
  Core::Nodes::Node** snodes = nullptr;
  Core::Nodes::Node** mnodes = nullptr;
  int nrow = sele.num_node();
  int ncol = mele.num_node();

  int linsize = 0;
  for (int i = 0; i < nrow; ++i) linsize += dynamic_cast<Node*>(sele.Nodes()[i])->GetLinsize();
  // safety
  linsize = linsize * 2;

  int ndof = Dim();

  snodes = sele.Nodes();
  mnodes = mele.Nodes();

  // get slave nodal coords for GP evaluation
  Core::LinAlg::SerialDenseMatrix scoord(3, nrow);
  sele.GetNodalCoords(scoord);

  // get master nodal coords for GP evaluation
  Core::LinAlg::SerialDenseMatrix mcoord(3, ncol);
  mele.GetNodalCoords(mcoord);

  std::array<double, 2> sgpx = {0.0, 0.0};
  std::array<double, 2> mgpx = {0.0, 0.0};

  for (int i = 0; i < nrow; ++i)
  {
    Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(snodes[i]);
    gpn[0] += sval[i] * mymrtrnode->MoData().n()[0];
    gpn[1] += sval[i] * mymrtrnode->MoData().n()[1];

    if (wear_type() == Inpar::Wear::wear_primvar)
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
    if (wear_side() == Inpar::Wear::wear_both and wear_type() == Inpar::Wear::wear_primvar)
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
  if (lengthn < 1.0e-12) FOUR_C_THROW("IntegrateAndDerivSegment: Divide by zero!");

  for (int i = 0; i < Dim(); ++i) gpn[i] /= lengthn;

  // build gap function at current GP
  for (int i = 0; i < Dim(); ++i) gap[0] += (mgpx[i] - sgpx[i]) * gpn[i];

  // **************************
  // Linearization
  // **************************

  // build directional derivative of slave GP normal (non-unit)
  Core::Gen::Pairedvector<int, double> dmap_nxsl_gp(ncol * ndof + linsize);
  Core::Gen::Pairedvector<int, double> dmap_nysl_gp(ncol * ndof + linsize);

  for (int i = 0; i < nrow; ++i)
  {
    Mortar::Node* snode = dynamic_cast<Mortar::Node*>(snodes[i]);

    Core::Gen::Pairedvector<int, double>& dmap_nxsl_i =
        dynamic_cast<CONTACT::Node*>(snodes[i])->Data().GetDerivN()[0];
    Core::Gen::Pairedvector<int, double>& dmap_nysl_i =
        dynamic_cast<CONTACT::Node*>(snodes[i])->Data().GetDerivN()[1];

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
  if (wear_type() == Inpar::Wear::wear_primvar)
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

        for (_CI p = frinode->Data().GetDerivN()[k].begin();
             p != frinode->Data().GetDerivN()[k].end(); ++p)
          dgapgp[p->first] += gpn[k] * sval[z] * w * (p->second);
      }
    }
  }
  else
  {
    for (int z = 0; z < nrow; ++z)
    {
      Mortar::Node* snode = dynamic_cast<Mortar::Node*>(snodes[z]);

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
  if (wear_side() == Inpar::Wear::wear_both and wear_type() == Inpar::Wear::wear_primvar)
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

        for (_CI p = frinode->Data().GetDerivN()[k].begin();
             p != frinode->Data().GetDerivN()[k].end(); ++p)
          dgapgp[p->first] -= gpn[k] * mval[z] * w * (p->second);
      }
    }
  }
  else
  {
    for (int z = 0; z < ncol; ++z)
    {
      Mortar::Node* mnode = dynamic_cast<Mortar::Node*>(mnodes[z]);

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
void inline CONTACT::Integrator::gp_3_d_g_quad_pwlin(Mortar::Element& sele,
    Mortar::IntElement& sintele, Mortar::Element& mele, Core::LinAlg::SerialDenseVector& sval,
    Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseVector& lmintval,
    Core::LinAlg::SerialDenseMatrix& scoord, Core::LinAlg::SerialDenseMatrix& mcoord,
    Core::LinAlg::SerialDenseMatrix& sderiv, Core::LinAlg::SerialDenseMatrix& mderiv, double* gap,
    double* gpn, double* lengthn, double& jac, double& wgt,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dmxigp,
    Core::Gen::Pairedvector<int, double>& dgapgp,
    std::vector<Core::Gen::Pairedvector<int, double>>& dnmap_unit)
{
  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();
  Core::Nodes::Node** sintnodes = sintele.Nodes();
  Core::Nodes::Node** mnodes = mele.Nodes();
  if (!snodes) FOUR_C_THROW("Null pointer!");
  if (!sintnodes) FOUR_C_THROW("Null pointer!");
  if (!mnodes) FOUR_C_THROW("Null pointer!");

  // number of nodes (slave, master)
  const int nrow = sele.num_node();
  const int nintrow = sintele.num_node();
  const int ncol = mele.num_node();

  std::array<double, 3> sgpx = {0.0, 0.0, 0.0};
  std::array<double, 3> mgpx = {0.0, 0.0, 0.0};

  for (int i = 0; i < nrow; ++i)
  {
    Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(snodes[i]);
    gpn[0] += sval[i] * mymrtrnode->MoData().n()[0];
    gpn[1] += sval[i] * mymrtrnode->MoData().n()[1];
    gpn[2] += sval[i] * mymrtrnode->MoData().n()[2];

    if (wear_type() == Inpar::Wear::wear_primvar)
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
    if (wear_side() == Inpar::Wear::wear_both and wear_type() == Inpar::Wear::wear_primvar)
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
  if (lengthn[0] < 1.0e-12) FOUR_C_THROW("IntegrateAndDerivSegment: Divide by zero!");

  for (int i = 0; i < 3; ++i) gpn[i] /= lengthn[0];

  // build gap function at current GP
  for (int i = 0; i < Dim(); ++i) gap[0] += (mgpx[i] - sgpx[i]) * gpn[i];

  // **************************
  // add to node
  // **************************
  // CASE 3: Standard LM shape functions and piecewise linear interpolation
  // Attention:  for this case, lmval represents lmintval !!!
  if (shape_fcn() == Inpar::Mortar::shape_standard &&
      lag_mult_quad() == Inpar::Mortar::lagmult_pwlin)
  {
    for (int j = 0; j < nintrow; ++j)
    {
      CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(sintnodes[j]);

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
    FOUR_C_THROW("Invalid integration case for 3D quadratic contact!");
  }

  // **************************
  // Linearization
  // **************************
  int linsize = 0;
  for (int i = 0; i < nrow; ++i)
  {
    Node* cnode = dynamic_cast<Node*>(snodes[i]);
    linsize += cnode->GetLinsize();
  }

  // build directional derivative of slave GP normal (non-unit)
  Core::Gen::Pairedvector<int, double> dmap_nxsl_gp(linsize);
  Core::Gen::Pairedvector<int, double> dmap_nysl_gp(linsize);
  Core::Gen::Pairedvector<int, double> dmap_nzsl_gp(linsize);

  for (int i = 0; i < nrow; ++i)
  {
    Node* cnode = dynamic_cast<Node*>(snodes[i]);

    Core::Gen::Pairedvector<int, double>& dmap_nxsl_i = cnode->Data().GetDerivN()[0];
    Core::Gen::Pairedvector<int, double>& dmap_nysl_i = cnode->Data().GetDerivN()[1];
    Core::Gen::Pairedvector<int, double>& dmap_nzsl_i = cnode->Data().GetDerivN()[2];

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
  if (wear_type() == Inpar::Wear::wear_primvar)
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

        for (_CI p = frinode->Data().GetDerivN()[k].begin();
             p != frinode->Data().GetDerivN()[k].end(); ++p)
          dgapgp[p->first] += gpn[k] * sval[z] * frinode->WearData().wcurr()[0] * (p->second);
      }
    }
  }
  else
  {
    for (int z = 0; z < nrow; ++z)
    {
      Node* cnode = dynamic_cast<Node*>(snodes[z]);

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
  if (wear_side() == Inpar::Wear::wear_both and wear_type() == Inpar::Wear::wear_primvar)
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

        for (_CI p = frinode->Data().GetDerivN()[k].begin();
             p != frinode->Data().GetDerivN()[k].end(); ++p)
          dgapgp[p->first] -= gpn[k] * mval[z] * frinode->WearData().wcurr()[0] * (p->second);
      }
    }
  }
  else
  {
    // lin master nodes
    for (int z = 0; z < ncol; ++z)
    {
      Node* cnode = dynamic_cast<Node*>(mnodes[z]);

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
void inline CONTACT::Integrator::gp_2_d_g_lin(int& iter, Mortar::Element& sele,
    Mortar::Element& mele, Core::LinAlg::SerialDenseVector& sval,
    Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseMatrix& sderiv, Core::LinAlg::SerialDenseMatrix& lmderiv, double& gap,
    double* gpn, double& jac, double& wgt, Core::Gen::Pairedvector<int, double>& dgapgp,
    Core::Gen::Pairedvector<int, double>& jacintcellmap,
    std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap)
{
  // get slave element nodes themselves
  Core::Nodes::Node** snodes = nullptr;
  int nrow = sele.num_node();
  snodes = sele.Nodes();

  Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(snodes[iter]);
  if (!mymrtrnode) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");
  if (mymrtrnode->IsOnBoundorCE()) return;

  double fac = 0.0;

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // get the corresponding map as a reference
  std::map<int, double>& dgmap = dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivG();

  // switch if Petrov-Galerkin approach for LM is applied
  if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
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
    if (shape_fcn() == Inpar::Mortar::shape_dual)
    {
      for (int m = 0; m < nrow; ++m)
      {
        fac = wgt * sval[m] * gap * jac;
        for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                 dualmap.begin();
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
  if (wear_type() == Inpar::Wear::wear_primvar)
  {
    // get master element nodes themselves
    const int ncol = mele.num_node();
    Core::Nodes::Node** mnodes = mele.Nodes();

    std::map<int, double>& dgwmmap = dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivGW();

    for (int bl = 0; bl < nrow; ++bl)
    {
      Mortar::Node* wearnode = dynamic_cast<Mortar::Node*>(snodes[bl]);
      for (int z = 0; z < Dim(); ++z)
        dgwmmap[wearnode->Dofs()[0]] +=
            jac * wgt * lmval[iter] * (gpn[z] * sval[bl] * wearnode->MoData().n()[z]);
    }

    if (wear_side() == Inpar::Wear::wear_both)
    {
      for (int bl = 0; bl < ncol; ++bl)
      {
        Mortar::Node* wearnodeM = dynamic_cast<Mortar::Node*>(mnodes[bl]);
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
void inline CONTACT::Integrator::gp_3_d_g_quad_pwlin_lin(int& iter, Mortar::IntElement& sintele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& lmintval,
    Core::LinAlg::SerialDenseMatrix& sderiv, Core::LinAlg::SerialDenseMatrix& lmintderiv,
    double& gap, double* gpn, double& jac, double& wgt,
    const Core::Gen::Pairedvector<int, double>& dgapgp,
    const Core::Gen::Pairedvector<int, double>& jacintcellmap,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp)
{
  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;

  // get slave element nodes themselves
  Core::Nodes::Node** sintnodes = sintele.Nodes();

  Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(sintnodes[iter]);
  if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

  double fac = 0.0;

  // CASE 3: Standard LM shape functions and piecewise linear interpolation
  if (shape_fcn() == Inpar::Mortar::shape_standard &&
      lag_mult_quad() == Inpar::Mortar::lagmult_pwlin)
  {
    // get the corresponding map as a reference
    std::map<int, double>& dgmap = dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivG();

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
    FOUR_C_THROW("shapefcn-lagmult combination not supported!");


  return;
}

/*----------------------------------------------------------------------*
 |  Do lin. entries for weighted Gap at GP                   farah 12/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::Integrator::gp_3_d_g_quad_lin(int& iter, Mortar::Element& sele,
    Mortar::Element& mele, Core::LinAlg::SerialDenseVector& sval,
    Core::LinAlg::SerialDenseVector& svalmod, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseMatrix& sderiv, Core::LinAlg::SerialDenseMatrix& lmderiv, double& gap,
    double* gpn, double& jac, double& wgt, bool& duallin,
    const Core::Gen::Pairedvector<int, double>& dgapgp,
    const Core::Gen::Pairedvector<int, double>& jacintcellmap,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dpsxigp,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap, bool dualquad3d)
{
  const int nrow = sele.num_node();

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();

  Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(snodes[iter]);
  if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

  double fac = 0.0;

  // compute cell gap linearization ************************************
  // CASE 1/2: Standard LM shape functions and quadratic or linear interpolation
  if (shape_fcn() == Inpar::Mortar::shape_standard &&
      (lag_mult_quad() == Inpar::Mortar::lagmult_quad ||
          lag_mult_quad() == Inpar::Mortar::lagmult_lin))
  {
    // get the corresponding map as a reference
    std::map<int, double>& dgmap = dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivG();

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
  else if ((shape_fcn() == Inpar::Mortar::shape_dual ||
               shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) &&
           (lag_mult_quad() == Inpar::Mortar::lagmult_quad ||
               lag_mult_quad() == Inpar::Mortar::lagmult_lin))
  {
    // get the corresponding map as a reference
    std::map<int, double>& dgmap = dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivG();

    if (shape_fcn() == Inpar::Mortar::shape_dual)
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

          for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                   dualmap.begin();
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
    else if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
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
void CONTACT::Integrator::gp_g_lin(int& iter, Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& mval,
    Core::LinAlg::SerialDenseVector& lmval, Core::LinAlg::SerialDenseMatrix& sderiv,
    Core::LinAlg::SerialDenseMatrix& lmderiv, double& gap, double* gpn, double& jac, double& wgt,
    Core::Gen::Pairedvector<int, double>& dgapgp,
    Core::Gen::Pairedvector<int, double>& jacintcellmap,
    std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap)
{
  const int nrow = sele.num_node();
  const int ncol = mele.num_node();

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();
  Core::Nodes::Node** mnodes = mele.Nodes();

  Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(snodes[iter]);
  if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

  if (mymrtrnode->IsOnBoundorCE()) return;

  static double fac = 0.0;

  // get the corresponding map as a reference
  std::map<int, double>& dgmap = dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivG();

  // switch if Petrov-Galerkin approach for LM is applied
  if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
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
        for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                 dualmap.begin();
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
  if (wear_type() == Inpar::Wear::wear_primvar)
  {
    std::map<int, double>& dgwmmap = dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivGW();

    for (int bl = 0; bl < nrow; ++bl)
    {
      Mortar::Node* wearnode = dynamic_cast<Mortar::Node*>(snodes[bl]);
      for (int z = 0; z < Dim(); ++z)
        dgwmmap[wearnode->Dofs()[0]] +=
            jac * wgt * lmval[iter] * (gpn[z] * sval[bl] * wearnode->MoData().n()[z]);
    }

    if (wear_side() == Inpar::Wear::wear_both)
    {
      for (int bl = 0; bl < ncol; ++bl)
      {
        Mortar::Node* wearnodeM = dynamic_cast<Mortar::Node*>(mnodes[bl]);
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
void CONTACT::Integrator::gp_3_d_dm_lin_bound(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& mval,
    Core::LinAlg::SerialDenseVector& lmval, Core::LinAlg::SerialDenseMatrix& sderiv,
    Core::LinAlg::SerialDenseMatrix& lmderiv, Core::LinAlg::SerialDenseMatrix& mderiv, double& jac,
    double& wgt, const Core::Gen::Pairedvector<int, double>& derivjac,
    std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    std::vector<Core::Gen::Pairedvector<int, double>>& dmxigp,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap)
{
  const int nrow = sele.num_node();
  const int ncol = mele.num_node();

  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  if (shape_fcn() == Inpar::Mortar::shape_standard)
  {
    // integrate LinD
    for (int j = 0; j < nrow; ++j)
    {
      Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(snodes[j]);
      if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

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
            dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivM()[mgid];

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
        Mortar::Node* mymrtrnode2 = dynamic_cast<Mortar::Node*>(snodes[k]);
        if (!mymrtrnode2) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

        // global master node ID
        int sgid = mymrtrnode2->Id();
        static double fac = 0.0;

        // node k is boundary node
        if (mymrtrnode2->IsOnBoundorCE())
        {
          // get the correct map pointer
          std::map<int, double>* dmmap_jk = nullptr;
          double sign = 0.;

          // for boundary nodes, we assemble to M-matrix, since the corresponding nodes do not
          // carry a LM-dof
          if (mymrtrnode2->IsOnBound())
          {
            dmmap_jk = &(dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivM()[sgid]);
            sign = -1.;
          }
          // for corners/edges, we assemble to D-matrix, since the corresponding nodes
          // carry a LM-dof, but are integrated some place else. Here we only do surface
          // integration.
          else if (mymrtrnode2->IsOnCornerEdge())
          {
            dmmap_jk = &(dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivD()[sgid]);
            sign = +1.;
          }
          else
            FOUR_C_THROW("you should never be here");

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
              dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivD()[sgid];

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
  else if (shape_fcn() == Inpar::Mortar::shape_dual ||
           shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
  {
    // integrate LinD
    for (int j = 0; j < nrow; ++j)
    {
      Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(snodes[j]);
      if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

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
            dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivM()[mgid];

        // (1) Lin(Phi) - dual shape functions
        if (dualmap.size() > 0)
        {
          for (int m = 0; m < nrow; ++m)
          {
            fac = wgt * sval[m] * mval[k] * jac;
            for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
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
        Mortar::Node* mymrtrnode2 = dynamic_cast<Mortar::Node*>(snodes[k]);
        if (!mymrtrnode2) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

        // global master node ID
        int sgid = mymrtrnode2->Id();
        static double fac = 0.0;

        // node k is boundary node
        if (mymrtrnode2->IsOnBoundorCE())
        {
          // get the correct map pointer
          std::map<int, double>* dmmap_jk = nullptr;
          double sign = 0.;

          // for boundary nodes, we assemble to M-matrix, since the corresponding nodes do not
          // carry a LM-dof
          if (mymrtrnode2->IsOnBound())
          {
            dmmap_jk = &(dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivM()[sgid]);
            sign = -1.;
          }
          // for corners/edges, we assemble to D-matrix, since the corresponding nodes
          // carry a LM-dof, but are integrated some place else. Here we only do surface
          // integration.
          else if (mymrtrnode2->IsOnCornerEdge())
          {
            dmmap_jk = &(dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivD()[sgid]);
            sign = +1.;
          }
          else
            FOUR_C_THROW("you should never be here");

          // (1) Lin(Phi) - dual shape functions
          if (dualmap.size() > 0)
          {
            for (int m = 0; m < nrow; ++m)
            {
              fac = sign * wgt * sval[m] * sval[k] * jac;
              for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
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
              dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivD()[mymrtrnode->Id()];

          // (1) Lin(Phi) - dual shape functions
          if (dualmap.size() > 0)
          {
            for (int m = 0; m < nrow; ++m)
            {
              fac = wgt * sval[m] * sval[k] * jac;
              for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
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
    FOUR_C_THROW("unknwon shape function type!");

  return;
}

/*----------------------------------------------------------------------*
 |  Compute entries for weighted Gap at GP                   farah 07/16|
 *----------------------------------------------------------------------*/
void inline CONTACT::Integrator::gp_2_d_dm_lin_bound(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& mval,
    Core::LinAlg::SerialDenseVector& lmval, Core::LinAlg::SerialDenseMatrix& sderiv,
    Core::LinAlg::SerialDenseMatrix& mderiv, Core::LinAlg::SerialDenseMatrix& lmderiv, double& jac,
    double& wgt, const Core::Gen::Pairedvector<int, double>& derivjac,
    std::vector<Core::Gen::Pairedvector<int, double>>& derivsxi,
    std::vector<Core::Gen::Pairedvector<int, double>>& derivmxi,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap)
{
  const int nrow = sele.num_node();
  const int ncol = mele.num_node();

  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  if (shape_fcn() == Inpar::Mortar::shape_standard)
  {
    // integrate LinD
    for (int j = 0; j < nrow; ++j)
    {
      Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(snodes[j]);
      if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

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
            dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivM()[mgid];

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
        Mortar::Node* mymrtrnode2 = dynamic_cast<Mortar::Node*>(snodes[k]);
        if (!mymrtrnode2) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

        // global master node ID
        int sgid = mymrtrnode2->Id();
        static double fac = 0.0;

        // node k is boundary node
        if (mymrtrnode2->IsOnBoundorCE())
        {
          // move entry to derivM (with minus sign)
          // get the correct map as a reference
          std::map<int, double>& dmmap_jk =
              dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivM()[sgid];

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
              dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivD()[sgid];

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
  else if (shape_fcn() == Inpar::Mortar::shape_dual ||
           shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
  {
    // integrate LinD
    for (int j = 0; j < nrow; ++j)
    {
      Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(snodes[j]);
      if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

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
            dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivM()[mgid];

        // (1) Lin(Phi) - dual shape functions
        if (dualmap.size() > 0)
        {
          for (int m = 0; m < nrow; ++m)
          {
            fac = wgt * sval[m] * mval[k] * jac;
            for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
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
        Mortar::Node* mymrtrnode2 = dynamic_cast<Mortar::Node*>(snodes[k]);
        if (!mymrtrnode2) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

        // global master node ID
        int sgid = mymrtrnode2->Id();
        static double fac = 0.0;

        // node k is boundary node
        if (mymrtrnode2->IsOnBoundorCE())
        {
          // move entry to derivM (with minus sign)
          // get the correct map as a reference
          std::map<int, double>& dmmap_jk =
              dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivM()[sgid];

          // (1) Lin(Phi) - dual shape functions
          if (dualmap.size() > 0)
          {
            for (int m = 0; m < nrow; ++m)
            {
              fac = wgt * sval[m] * sval[k] * jac;
              for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
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
              dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivD()[mymrtrnode->Id()];

          // (1) Lin(Phi) - dual shape functions
          if (dualmap.size() > 0)
          {
            for (int m = 0; m < nrow; ++m)
            {
              fac = wgt * sval[m] * sval[k] * jac;
              for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
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
    FOUR_C_THROW("unknwon shape function type!");

  return;
}


/*----------------------------------------------------------------------*
 |  Lin D and M matrix entries at GP                         farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::Integrator::gp_2_d_dm_ele_lin(int& iter, bool& bound, Mortar::Element& sele,
    Mortar::Element& mele, Core::LinAlg::SerialDenseVector& sval,
    Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseMatrix& mderiv, double& dxdsxi, double& wgt,
    const Core::Gen::Pairedvector<int, double>& dmxigp,
    const Core::Gen::Pairedvector<int, double>& derivjac,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap)
{
  Core::Nodes::Node** snodes = nullptr;
  Core::Nodes::Node** mnodes = nullptr;

  int nrow = sele.num_node();
  int ncol = mele.num_node();

  snodes = sele.Nodes();
  mnodes = mele.Nodes();

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(snodes[iter]);
  if (!mymrtrnode) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  const int sgid = mymrtrnode->Id();

  // standard shape functions
  if (shape_fcn() == Inpar::Mortar::shape_standard)
  {
    // integrate LinM
    for (int k = 0; k < ncol; ++k)
    {
      // global master node ID
      int mgid = mnodes[k]->Id();
      double fac = 0.0;


      // get the correct map as a reference
      std::map<int, double>& dmmap_jk =
          dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivM()[mgid];

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
          dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivD()[sgid];

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
  else if (shape_fcn() == Inpar::Mortar::shape_dual ||
           shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
  {
    // get the D-map as a reference
    std::map<int, double>& ddmap_jk =
        dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivD()[sgid];

    // integrate LinM and LinD (NO boundary modification)
    for (int k = 0; k < ncol; ++k)
    {
      // global master node ID
      int mgid = mnodes[k]->Id();
      double fac = 0.0;

      // get the correct map as a reference
      std::map<int, double>& dmmap_jk =
          dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivM()[mgid];

      // (1) Lin(Phi) - dual shape functions
      for (int m = 0; m < nrow; ++m)
      {
        fac = wgt * sval[m] * mval[k] * dxdsxi;
        for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                 dualmap.begin();
             p != dualmap.end(); ++p)
        {
          dmmap_jk[p->first] += fac * (p->second)(iter, m);
          if (!bound) ddmap_jk[p->first] += fac * (p->second)(iter, m);
        }
      }

      // (2) Lin(Phi) - slave GP coordinates --> 0

      // (3) Lin(NMaster) - master GP coordinates
      fac = wgt * lmval[iter] * mderiv(k, 0) * dxdsxi;
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
void inline CONTACT::Integrator::gp_2_d_dm_lin(int& iter, bool& bound, bool& linlm,
    Mortar::Element& sele, Mortar::Element& mele, Core::LinAlg::SerialDenseVector& sval,
    Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseMatrix& sderiv, Core::LinAlg::SerialDenseMatrix& mderiv,
    Core::LinAlg::SerialDenseMatrix& lmderiv, double& jac, double& wgt,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dmxigp,
    const Core::Gen::Pairedvector<int, double>& derivjac,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap)
{
  Core::Nodes::Node** snodes = nullptr;
  Core::Nodes::Node** mnodes = nullptr;
  int nrow = sele.num_node();
  int ncol = mele.num_node();
  snodes = sele.Nodes();
  mnodes = mele.Nodes();

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // **************** no edge modification *****************************
  // (and LinM also for edge node modification case)

  // check for linear LM interpolation in quadratic FE
  if (linlm)
  {
    Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(snodes[iter]);
    if (!mymrtrnode) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");
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
            dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivM()[mgid];

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
        Mortar::Node* mymrtrnode2 = dynamic_cast<Mortar::Node*>(snodes[k]);
        if (!mymrtrnode2) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");
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
              dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivM()[sgid];

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
              dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivD()[sgid];

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
    Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(snodes[iter]);
    if (!mymrtrnode) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

    int sgid = mymrtrnode->Id();

    // standard shape functions
    if (shape_fcn() == Inpar::Mortar::shape_standard)
    {
      // integrate LinM
      for (int k = 0; k < ncol; ++k)
      {
        // global master node ID
        int mgid = mnodes[k]->Id();
        double fac = 0.0;

        // get the correct map as a reference
        std::map<int, double>& dmmap_jk =
            dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivM()[mgid];

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
            dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivD()[sgid];

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
    else if (shape_fcn() == Inpar::Mortar::shape_dual ||
             shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
    {
      // get the D-map as a reference
      std::map<int, double>& ddmap_jk =
          dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivD()[sgid];

      // integrate LinM and LinD (NO boundary modification)
      for (int k = 0; k < ncol; ++k)
      {
        // global master node ID
        int mgid = mnodes[k]->Id();
        double fac = 0.0;

        // get the correct map as a reference
        std::map<int, double>& dmmap_jk =
            dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivM()[mgid];

        // (1) Lin(Phi) - dual shape functions
        for (int m = 0; m < nrow; ++m)
        {
          fac = wgt * sval[m] * mval[k] * jac;
          for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                   dualmap.begin();
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
        fac = wgt * lmval[iter] * mderiv(k, 0) * jac;

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
void inline CONTACT::Integrator::gp_3_d_dm_quad_pwlin_lin(int& iter, Mortar::Element& sele,
    Mortar::Element& sintele, Mortar::Element& mele, Core::LinAlg::SerialDenseVector& sval,
    Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseVector& lmintval,
    Core::LinAlg::SerialDenseMatrix& sderiv, Core::LinAlg::SerialDenseMatrix& mderiv,
    Core::LinAlg::SerialDenseMatrix& lmintderiv, double& wgt, double& jac,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dpsxigp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dpmxigp,
    const Core::Gen::Pairedvector<int, double>& jacintcellmap)
{
  const int nrow = sele.num_node();
  const int ncol = mele.num_node();

  // get slave element nodes themselves
  Core::Nodes::Node** sintnodes = sintele.Nodes();

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;

  // **************** no edge modification *****************************
  // (and LinM also for edge node modification case)

  Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(sintnodes[iter]);
  if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

  // integrate LinM
  for (int k = 0; k < ncol; ++k)
  {
    // global master node ID
    int mgid = mele.Nodes()[k]->Id();
    double fac = 0.0;

    // get the correct map as a reference
    std::map<int, double>& dmmap_jk =
        dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivM()[mgid];

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
        dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivD()[sgid];

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
void inline CONTACT::Integrator::gp_3_d_dm_quad_lin(bool& duallin, Mortar::Element& sele,
    Mortar::Element& mele, Core::LinAlg::SerialDenseVector& sval,
    Core::LinAlg::SerialDenseVector& svalmod, Core::LinAlg::SerialDenseVector& mval,
    Core::LinAlg::SerialDenseVector& lmval, Core::LinAlg::SerialDenseMatrix& sderiv,
    Core::LinAlg::SerialDenseMatrix& mderiv, Core::LinAlg::SerialDenseMatrix& lmderiv, double& wgt,
    double& jac, const std::vector<Core::Gen::Pairedvector<int, double>>& dpsxigp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dpmxigp,
    const Core::Gen::Pairedvector<int, double>& jacintcellmap,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap, bool dualquad3d)
{
  const int nrow = sele.num_node();
  const int ncol = mele.num_node();

  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();

  std::vector<Mortar::Node*> smnodes(nrow);
  for (int i = 0; i < nrow; ++i) smnodes[i] = dynamic_cast<Mortar::Node*>(snodes[i]);

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // **************** no edge modification *****************************
  // (and LinM also for edge node modification case)

  // compute cell D/M linearization ************************************
  // CASE 1: Standard LM shape functions and quadratic interpolation
  if (shape_fcn() == Inpar::Mortar::shape_standard &&
      lag_mult_quad() == Inpar::Mortar::lagmult_quad)
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
        Core::LinAlg::SerialDenseMatrix& dMderiv =
            dynamic_cast<CONTACT::Element&>(sele).GetDderiv()[p->first];
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
      Core::LinAlg::SerialDenseMatrix& dMderiv =
          dynamic_cast<CONTACT::Element&>(sele).GetDderiv()[p->first];
      for (int j = 0; j < nrow; ++j)
        for (int k = 0; k < nrow; ++k) dMderiv(j, k) += wgt * lmval[j] * sval[k] * (p->second);
    }

    // (1) Lin(Phi) - dual shape functions
    // this vanishes here since there are no deformation-dependent dual functions

    // (2) Lin(NSlave) - slave GP coordinates
    for (int d = 0; d < 2; ++d)
      for (_CI p = dpsxigp[d].begin(); p != dpsxigp[d].end(); ++p)
      {
        Core::LinAlg::SerialDenseMatrix& mMderiv =
            dynamic_cast<CONTACT::Element&>(sele).GetMderiv()[p->first];
        for (int j = 0; j < nrow; ++j)
          for (int k = 0; k < ncol; ++k)
            mMderiv(j, k) += wgt * lmderiv(j, d) * mval[k] * jac * (p->second);
      }

    // (3) Lin(NMaster) - master GP coordinates
    for (int d = 0; d < 2; ++d)
      for (_CI p = dpmxigp[d].begin(); p != dpmxigp[d].end(); ++p)
      {
        Core::LinAlg::SerialDenseMatrix& mMderiv =
            dynamic_cast<CONTACT::Element&>(sele).GetMderiv()[p->first];
        for (int j = 0; j < nrow; ++j)
          for (int k = 0; k < ncol; ++k)
            mMderiv(j, k) += wgt * lmval[j] * mderiv(k, d) * jac * (p->second);
      }

    // (4) Lin(dsxideta) - intcell GP Jacobian
    for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
    {
      Core::LinAlg::SerialDenseMatrix& mMderiv =
          dynamic_cast<CONTACT::Element&>(sele).GetMderiv()[p->first];
      for (int j = 0; j < nrow; ++j)
        for (int k = 0; k < ncol; ++k) mMderiv(j, k) += wgt * lmval[j] * mval[k] * (p->second);
    }
  }

  // CASE 2: Standard LM shape functions and linear interpolation
  // CASE 5: dual LM shape functions and linear interpolation
  // (this has to be treated seperately here for LinDM because of bound)
  else if ((shape_fcn() == Inpar::Mortar::shape_standard ||
               shape_fcn() == Inpar::Mortar::shape_dual ||
               shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) &&
           lag_mult_quad() == Inpar::Mortar::lagmult_lin)
  {
    // integrate LinD
    for (int j = 0; j < nrow; ++j)
    {
      Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(snodes[j]);
      if (!mymrtrnode) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

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
              dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivM()[mgid];

          // (1) Lin(Phi) - dual shape functions
          if (duallin)
          {
            typedef Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator
                _CIM;
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
          Mortar::Node* mymrtrnode2 = dynamic_cast<Mortar::Node*>(snodes[k]);
          if (!mymrtrnode2) FOUR_C_THROW("integrate_deriv_cell3_d_aux_plane: Null pointer!");

          // global master node ID
          int sgid = mymrtrnode2->Id();
          static double fac = 0.0;

          // node k is boundary node
          if (mymrtrnode2->IsOnBound())
          {
            // move entry to derivM (with minus sign)
            // get the correct map as a reference
            std::map<int, double>& dmmap_jk =
                dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivM()[sgid];

            // (1) Lin(Phi) - dual shape functions
            if (duallin)
            {
              typedef Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator
                  _CIM;
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
                dynamic_cast<CONTACT::Node*>(mymrtrnode)->Data().GetDerivD()[sgid];

            // (1) Lin(Phi) - dual shape functions
            if (duallin)
            {
              typedef Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator
                  _CIM;
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
  else if ((shape_fcn() == Inpar::Mortar::shape_dual ||
               shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) &&
           lag_mult_quad() == Inpar::Mortar::lagmult_quad)
  {
    double fac = 0.;
    // (1) Lin(Phi) - dual shape functions
    if (duallin)
    {
      for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
               dualmap.begin();
           p != dualmap.end(); ++p)
      {
        Core::LinAlg::SerialDenseMatrix& dDderiv =
            dynamic_cast<CONTACT::Element&>(sele).GetDderiv()[p->first];
        Core::LinAlg::SerialDenseMatrix& dMderiv =
            dynamic_cast<CONTACT::Element&>(sele).GetMderiv()[p->first];
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
        Core::LinAlg::SerialDenseMatrix& dDderiv =
            dynamic_cast<CONTACT::Element&>(sele).GetDderiv()[p->first];
        Core::LinAlg::SerialDenseMatrix& dMderiv =
            dynamic_cast<CONTACT::Element&>(sele).GetMderiv()[p->first];
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
        Core::LinAlg::SerialDenseMatrix& dDderiv =
            dynamic_cast<CONTACT::Element&>(sele).GetDderiv()[p->first];
        Core::LinAlg::SerialDenseMatrix& dMderiv =
            dynamic_cast<CONTACT::Element&>(sele).GetMderiv()[p->first];
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
      Core::LinAlg::SerialDenseMatrix& dDderiv =
          dynamic_cast<CONTACT::Element&>(sele).GetDderiv()[p->first];
      Core::LinAlg::SerialDenseMatrix& dMderiv =
          dynamic_cast<CONTACT::Element&>(sele).GetMderiv()[p->first];
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
void inline CONTACT::Integrator::gp_3_d_dm_lin(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& mval,
    Core::LinAlg::SerialDenseVector& lmval, Core::LinAlg::SerialDenseMatrix& sderiv,
    Core::LinAlg::SerialDenseMatrix& mderiv, Core::LinAlg::SerialDenseMatrix& lmderiv, double& wgt,
    double& jac, std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    std::vector<Core::Gen::Pairedvector<int, double>>& dmxigp,
    Core::Gen::Pairedvector<int, double>& jacintcellmap,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap)
{
  const int nrow = sele.num_node();
  const int ncol = mele.num_node();

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // standard shape functions
  if (shape_fcn() == Inpar::Mortar::shape_standard)
  {
    // integrate LinM
    // (1) Lin(Phi) - dual shape functions
    // this vanishes here since there are no deformation-dependent dual functions

    // (2) Lin(NSlave) - slave GP coordinates
    for (int d = 0; d < 2; ++d)
      for (_CI p = dsxigp[d].begin(); p != dsxigp[d].end(); ++p)
      {
        Core::LinAlg::SerialDenseMatrix& dMderiv =
            dynamic_cast<CONTACT::Element&>(sele).GetMderiv()[p->first];
        for (int j = 0; j < nrow; ++j)
          for (int k = 0; k < ncol; ++k)
            dMderiv(j, k) += wgt * lmderiv(j, d) * mval[k] * jac * (p->second);
      }

    // (3) Lin(NMaster) - master GP coordinates
    for (int d = 0; d < 2; ++d)
      for (_CI p = dmxigp[d].begin(); p != dmxigp[d].end(); ++p)
      {
        Core::LinAlg::SerialDenseMatrix& dMderiv =
            dynamic_cast<CONTACT::Element&>(sele).GetMderiv()[p->first];
        for (int j = 0; j < nrow; ++j)
          for (int k = 0; k < ncol; ++k)
            dMderiv(j, k) += wgt * lmval[j] * mderiv(k, d) * jac * (p->second);
      }

    // (4) Lin(dsxideta) - intcell GP Jacobian
    for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
    {
      Core::LinAlg::SerialDenseMatrix& dMderiv =
          dynamic_cast<CONTACT::Element&>(sele).GetMderiv()[p->first];
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
        Core::LinAlg::SerialDenseMatrix& dDderiv =
            dynamic_cast<CONTACT::Element&>(sele).GetDderiv()[p->first];
        for (int j = 0; j < nrow; ++j)
          for (int k = 0; k < nrow; ++k)
            dDderiv(j, k) +=
                (wgt * lmderiv(j, d) * sval[k] * jac + wgt * lmval[j] * sderiv(k, d) * jac) *
                (p->second);
      }

    // (4) Lin(dsxideta) - intcell GP Jacobian
    for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
    {
      Core::LinAlg::SerialDenseMatrix& dDderiv =
          dynamic_cast<CONTACT::Element&>(sele).GetDderiv()[p->first];
      for (int j = 0; j < nrow; ++j)
        for (int k = 0; k < nrow; ++k) dDderiv(j, k) += wgt * lmval[j] * sval[k] * (p->second);
    }
  }

  // dual shape functions
  else if (shape_fcn() == Inpar::Mortar::shape_dual ||
           shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
  {
    double fac = 0.;

    // (1) Lin(Phi) - dual shape functions
    if (dualmap.size() > 0)
      for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
               dualmap.begin();
           p != dualmap.end(); ++p)
      {
        Core::LinAlg::SerialDenseMatrix& dDderiv =
            dynamic_cast<CONTACT::Element&>(sele).GetDderiv()[p->first];
        Core::LinAlg::SerialDenseMatrix& dMderiv =
            dynamic_cast<CONTACT::Element&>(sele).GetMderiv()[p->first];
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
        Core::LinAlg::SerialDenseMatrix& dDderiv =
            dynamic_cast<CONTACT::Element&>(sele).GetDderiv()[p->first];
        Core::LinAlg::SerialDenseMatrix& dMderiv =
            dynamic_cast<CONTACT::Element&>(sele).GetMderiv()[p->first];
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
        Core::LinAlg::SerialDenseMatrix& dDderiv =
            dynamic_cast<CONTACT::Element&>(sele).GetDderiv()[p->first];
        Core::LinAlg::SerialDenseMatrix& dMderiv =
            dynamic_cast<CONTACT::Element&>(sele).GetMderiv()[p->first];
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
      Core::LinAlg::SerialDenseMatrix& dDderiv =
          dynamic_cast<CONTACT::Element&>(sele).GetDderiv()[p->first];
      Core::LinAlg::SerialDenseMatrix& dMderiv =
          dynamic_cast<CONTACT::Element&>(sele).GetMderiv()[p->first];
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
void inline CONTACT::Integrator::gp_d2(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& lm2val, Core::LinAlg::SerialDenseVector& m2val, double& jac,
    double& wgt, const Epetra_Comm& comm)
{
  int ncol = mele.num_node();
  int ndof = Dim();

  // get slave element nodes themselves
  Core::Nodes::Node** mnodes = mele.Nodes();
  if (!mnodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  if (shape_fcn() == Inpar::Mortar::shape_dual ||
      shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
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
    FOUR_C_THROW("Both-sided wear just for dual shape functions!");

  return;
}



/*----------------------------------------------------------------------*
 |  Compute wear at GP (for expl/impl algor.)                farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::Integrator::gp_2_d_wear(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseMatrix& sderiv,
    Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseMatrix& mderiv,
    Core::LinAlg::SerialDenseVector& lmval, Core::LinAlg::SerialDenseMatrix& lmderiv,
    Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> lagmult, double* gpn, double& jac, double& wgt,
    double* jumpval, double* wearval, Core::Gen::Pairedvector<int, double>& dsliptmatrixgp,
    Core::Gen::Pairedvector<int, double>& dweargp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dmxigp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dnmap_unit,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap)
{
  const int nrow = sele.num_node();
  const int ncol = mele.num_node();
  const int ndof = Dim();

  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();
  Core::Nodes::Node** mnodes = mele.Nodes();

  int linsize = 0;
  for (int i = 0; i < sele.num_node(); ++i)
    linsize += dynamic_cast<Node*>(sele.Nodes()[i])->GetLinsize();
  linsize = linsize * 2;

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  //***********************************************************************
  // Here, the tangential relative slip increment is used and NOT the
  // nodal weighted tangential relative slip increment !!!
  // The reason for that is the slip for wear calculation is written
  // within the integral --> no double weighting allowed !
  // The wearcoefficient is not included in this calculation
  //***********************************************************************
  // for wearval
  std::array<double, 2> gpt = {0.0, 0.0};
  std::array<double, 2> gplm = {0.0, 0.0};
  std::array<double, 2> sgpjump = {0.0, 0.0};
  std::array<double, 2> mgpjump = {0.0, 0.0};
  std::array<double, 2> jump = {0.0, 0.0};

  // for linearization
  double lm_lin = 0.0;
  double lengtht = 0.0;

  for (int i = 0; i < nrow; ++i)
  {
    CONTACT::Node* myconode = dynamic_cast<CONTACT::Node*>(snodes[i]);

    // nodal tangent interpolation
    gpt[0] += sval[i] * myconode->Data().txi()[0];
    gpt[1] += sval[i] * myconode->Data().txi()[1];

    // delta D
    sgpjump[0] += sval[i] * (sele.GetNodalCoords(0, i) - (sele.GetNodalCoordsOld(0, i)));
    sgpjump[1] += sval[i] * (sele.GetNodalCoords(1, i) - (sele.GetNodalCoordsOld(1, i)));

    // LM interpolation
    gplm[0] += lmval[i] * ((*lagmult)(0, i));
    gplm[1] += lmval[i] * ((*lagmult)(1, i));
  }

  // normalize interpolated GP tangent back to length 1.0 !!!
  lengtht = sqrt(gpt[0] * gpt[0] + gpt[1] * gpt[1]);
  if (abs(lengtht) < 1.0e-12) FOUR_C_THROW("IntegrateAndDerivSegment: Divide by zero!");

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
  if (wearimpl_ and wear_type() != Inpar::Wear::wear_primvar)
    wearval[0] = (wearval[0]) * abs(jumpval[0]);
  else
    wearval[0] = abs(wearval[0]) * abs(jumpval[0]);

  // compute segment wear vector ***************************************
  // nrow represents the slave side dofs !!!
  for (int j = 0; j < nrow; ++j)
  {
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(snodes[j]);

    double prod = 0.0;
    if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
      prod = sval[j] * wearval[0] * jac * wgt;
    else
      prod = lmval[j] * wearval[0] * jac * wgt;

    // add current Gauss point's contribution to wseg
    cnode->add_delta_weighted_wear_value(prod);
  }

  //****************************************************************
  //   linearization for implicit algorithms
  //****************************************************************
  if (wearimpl_ || wear_type() == Inpar::Wear::wear_primvar)
  {
    // evaluate the GP wear function derivatives
    Core::Gen::Pairedvector<int, double> ddualgp_x(ndof * ncol + linsize);
    Core::Gen::Pairedvector<int, double> ddualgp_y(ndof * ncol + linsize);

    Core::Gen::Pairedvector<int, double> ddualgp_x_sxi(ndof * ncol + linsize);
    Core::Gen::Pairedvector<int, double> ddualgp_y_sxi(ndof * ncol + linsize);

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
        for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                 dualmap.begin();
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
    Core::Gen::Pairedvector<int, double> dmap_txsl_gp(ndof * ncol + linsize);
    Core::Gen::Pairedvector<int, double> dmap_tysl_gp(ndof * ncol + linsize);

    for (int i = 0; i < nrow; ++i)
    {
      CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(snodes[i]);

      Core::Gen::Pairedvector<int, double>& dmap_txsl_i =
          dynamic_cast<CONTACT::Node*>(cnode)->Data().GetDerivTxi()[0];
      Core::Gen::Pairedvector<int, double>& dmap_tysl_i =
          dynamic_cast<CONTACT::Node*>(cnode)->Data().GetDerivTxi()[1];

      for (_CI p = dmap_txsl_i.begin(); p != dmap_txsl_i.end(); ++p)
        dmap_txsl_gp[p->first] += sval[i] * (p->second);
      for (_CI p = dmap_tysl_i.begin(); p != dmap_tysl_i.end(); ++p)
        dmap_tysl_gp[p->first] += sval[i] * (p->second);

      for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
      {
        double valx = sderiv(i, 0) * dynamic_cast<CONTACT::Node*>(cnode)->Data().txi()[0];
        dmap_txsl_gp[p->first] += valx * (p->second);
        double valy = sderiv(i, 0) * dynamic_cast<CONTACT::Node*>(cnode)->Data().txi()[1];
        dmap_tysl_gp[p->first] += valy * (p->second);
      }
    }

    // (b) build directional derivative of slave GP tagent (unit)
    Core::Gen::Pairedvector<int, double> dmap_txsl_gp_unit(ndof * ncol + linsize);
    Core::Gen::Pairedvector<int, double> dmap_tysl_gp_unit(ndof * ncol + linsize);

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
    Core::Gen::Pairedvector<int, double> dmap_slcoord_gp_x(ndof * ncol + linsize);
    Core::Gen::Pairedvector<int, double> dmap_slcoord_gp_y(ndof * ncol + linsize);

    Core::Gen::Pairedvector<int, double> dmap_mcoord_gp_x(ndof * ncol + linsize);
    Core::Gen::Pairedvector<int, double> dmap_mcoord_gp_y(ndof * ncol + linsize);

    Core::Gen::Pairedvector<int, double> dmap_coord_x(ndof * ncol + linsize);
    Core::Gen::Pairedvector<int, double> dmap_coord_y(ndof * ncol + linsize);

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
      CONTACT::Node* snode = dynamic_cast<CONTACT::Node*>(snodes[i]);

      dmap_slcoord_gp_x[snode->Dofs()[0]] += sval[i];
      dmap_slcoord_gp_y[snode->Dofs()[1]] += sval[i];
    }
    // deriv master x-coords
    for (int i = 0; i < ncol; ++i)
    {
      CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(mnodes[i]);

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
void inline CONTACT::Integrator::gp_3_d_wear(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseMatrix& sderiv,
    Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseMatrix& mderiv,
    Core::LinAlg::SerialDenseVector& lmval, Core::LinAlg::SerialDenseMatrix& lmderiv,
    Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> lagmult, double* gpn, double& jac, double& wgt,
    double* jumpval, double* wearval, Core::Gen::Pairedvector<int, double>& dsliptmatrixgp,
    Core::Gen::Pairedvector<int, double>& dweargp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dmxigp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dnmap_unit,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap)
{
  const int nrow = sele.num_node();
  const int ncol = mele.num_node();
  const int ndof = Dim();

  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();
  Core::Nodes::Node** mnodes = mele.Nodes();

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  //***********************************************************************
  // Here, the tangential relative slip increment is used and NOT the
  // nodal weighted tangential relative slip increment !!!
  // The reason for that is the slip for wear calculation is written
  // within the integral --> no double weighting allowed !
  // The wearcoefficient is not included in this calculation
  //***********************************************************************

  // wear-lin
  Core::LinAlg::Matrix<3, 1> jump;
  Core::LinAlg::Matrix<3, 1> jumptan;
  Core::LinAlg::Matrix<3, 1> lm;
  Core::LinAlg::Matrix<3, 1> lmtan;
  Core::LinAlg::Matrix<3, 3> tanplane;

  std::array<double, 3> gplm = {0.0, 0.0, 0.0};
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
  std::array<double, 3> sgpjump = {0.0, 0.0, 0.0};
  for (int i = 0; i < nrow; ++i)
  {
    sgpjump[0] += sval[i] * (sele.GetNodalCoords(0, i) - sele.GetNodalCoordsOld(0, i));
    sgpjump[1] += sval[i] * (sele.GetNodalCoords(1, i) - sele.GetNodalCoordsOld(1, i));
    sgpjump[2] += sval[i] * (sele.GetNodalCoords(2, i) - sele.GetNodalCoordsOld(2, i));
  }

  // interpolation of master GP jumps (relative displacement increment)
  std::array<double, 3> mgpjump = {0.0, 0.0, 0.0};
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
    if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
      prod = sval[j] * wearval[0] * jac * wgt;
    else
      prod = lmval[j] * wearval[0] * jac * wgt;

    // add current Gauss point's contribution to wseg
    cnode->add_delta_weighted_wear_value(prod);
  }

  // linearization without lm weighting and jac.
  if (wearimpl_ or wear_type() == Inpar::Wear::wear_primvar)
  {
    int linsize = 0;
    for (int i = 0; i < nrow; ++i)
    {
      Node* cnode = dynamic_cast<Node*>(snodes[i]);
      linsize += cnode->GetLinsize();
    }

    // evaluate the GP wear function derivatives
    Core::Gen::Pairedvector<int, double> ddualgp_x((nrow + ncol) * ndof + linsize);
    Core::Gen::Pairedvector<int, double> ddualgp_y((nrow + ncol) * ndof + linsize);
    Core::Gen::Pairedvector<int, double> ddualgp_z((nrow + ncol) * ndof + linsize);

    Core::Gen::Pairedvector<int, double> ddualgp_x_sxi((nrow + ncol) * ndof + linsize);
    Core::Gen::Pairedvector<int, double> ddualgp_y_sxi((nrow + ncol) * ndof + linsize);
    Core::Gen::Pairedvector<int, double> ddualgp_z_sxi((nrow + ncol) * ndof + linsize);

    std::vector<std::vector<Core::Gen::Pairedvector<int, double>>> tanggp(
        3, std::vector<Core::Gen::Pairedvector<int, double>>(3, (ncol * ndof) + linsize));

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
        for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                 dualmap.begin();
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

    Core::Gen::Pairedvector<int, double> dt0((ncol * ndof) + linsize);
    Core::Gen::Pairedvector<int, double> dt1((ncol * ndof) + linsize);
    Core::Gen::Pairedvector<int, double> dt2((ncol * ndof) + linsize);

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
    if (wear_type() == Inpar::Wear::wear_primvar)
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
    Core::Gen::Pairedvector<int, double> dmap_slcoord_gp_x((nrow + ncol) * ndof + linsize);
    Core::Gen::Pairedvector<int, double> dmap_slcoord_gp_y((nrow + ncol) * ndof + linsize);
    Core::Gen::Pairedvector<int, double> dmap_slcoord_gp_z((nrow + ncol) * ndof + linsize);

    Core::Gen::Pairedvector<int, double> dmap_mcoord_gp_x((nrow + ncol) * ndof + linsize);
    Core::Gen::Pairedvector<int, double> dmap_mcoord_gp_y((nrow + ncol) * ndof + linsize);
    Core::Gen::Pairedvector<int, double> dmap_mcoord_gp_z((nrow + ncol) * ndof + linsize);

    Core::Gen::Pairedvector<int, double> dmap_coord_x((nrow + ncol) * ndof + linsize);
    Core::Gen::Pairedvector<int, double> dmap_coord_y((nrow + ncol) * ndof + linsize);
    Core::Gen::Pairedvector<int, double> dmap_coord_z((nrow + ncol) * ndof + linsize);

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
      CONTACT::Node* snode = dynamic_cast<CONTACT::Node*>(snodes[i]);

      dmap_slcoord_gp_x[snode->Dofs()[0]] += sval[i];
      dmap_slcoord_gp_y[snode->Dofs()[1]] += sval[i];
      dmap_slcoord_gp_z[snode->Dofs()[2]] += sval[i];
    }
    // deriv master x-coords
    for (int i = 0; i < ncol; ++i)
    {
      CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(mnodes[i]);

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
    Core::Gen::Pairedvector<int, double> lintan0((nrow + ncol) * ndof + linsize);
    Core::Gen::Pairedvector<int, double> lintan1((nrow + ncol) * ndof + linsize);
    Core::Gen::Pairedvector<int, double> lintan2((nrow + ncol) * ndof + linsize);

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
    if (wear_type() == Inpar::Wear::wear_primvar)
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
void inline CONTACT::Integrator::gp_te(Mortar::Element& sele,
    Core::LinAlg::SerialDenseVector& lmval, Core::LinAlg::SerialDenseVector& sval, double& jac,
    double& wgt, double* jumpval)
{
  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();
  if (!snodes) FOUR_C_THROW("Null pointer!");

  int nrow = sele.num_node();

  if (wear_shape_fcn() == Inpar::Wear::wear_shape_standard)
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
  else if (wear_shape_fcn() == Inpar::Wear::wear_shape_dual)
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
    FOUR_C_THROW("Chosen wear shape function not supported!");


  return;
}

/*----------------------------------------------------------------------*
 |  Compute entries for E and T matrix at GP (Master)        farah 11/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::Integrator::gp_te_master(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& lmval, Core::LinAlg::SerialDenseVector& lm2val,
    Core::LinAlg::SerialDenseVector& mval, double& jac, double& wgt, double* jumpval,
    const Epetra_Comm& comm)
{
  if (sele.Owner() != comm.MyPID()) return;

  // mele is involved for both-sided wear
  mele.SetAttached() = true;

  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();
  if (!snodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  // get master element nodes themselves
  Core::Nodes::Node** mnodes = mele.Nodes();
  if (!mnodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  int nrow = mele.num_node();
  int ncol = sele.num_node();

  if (wear_shape_fcn() == Inpar::Wear::wear_shape_standard)
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
  else if (wear_shape_fcn() == Inpar::Wear::wear_shape_dual)
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
    FOUR_C_THROW("Chosen wear shape function not supported!");

  return;
}

/*----------------------------------------------------------------------*
 |  Compute lin for T and E matrix -- discr. wear            farah 11/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::Integrator::gp_2_d_te_master_lin(int& iter,  // like k
    Mortar::Element& sele, Mortar::Element& mele, Core::LinAlg::SerialDenseVector& sval,
    Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseMatrix& mderiv, Core::LinAlg::SerialDenseMatrix& lmderiv,
    double& dsxideta, double& dxdsxi, double& dxdsxidsxi, double& wgt, double* jumpval,
    const Core::Gen::Pairedvector<int, double>& dsxigp,
    const Core::Gen::Pairedvector<int, double>& dmxigp,
    const Core::Gen::Pairedvector<int, double>& derivjac,
    const Core::Gen::Pairedvector<int, double>& dsliptmatrixgp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& ximaps,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap,
    const Epetra_Comm& comm)
{
  if (sele.Owner() != comm.MyPID()) return;

  const int nrow = sele.num_node();
  const int ncol = mele.num_node();

  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();
  if (!snodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  // get master element nodes themselves
  Core::Nodes::Node** mnodes = mele.Nodes();
  if (!mnodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(mnodes[iter]);
  if (!mymrtrnode) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  if (wear_shape_fcn() == Inpar::Wear::wear_shape_standard)
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
        for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                 dualmap.begin();
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
  else if (wear_shape_fcn() ==
           Inpar::Wear::wear_shape_dual)  //******************************************
  {
    FOUR_C_THROW("Chosen shapefunctions 'wear_shape_dual' not supported!");
  }
  else
    FOUR_C_THROW("Chosen shapefunctions not supported!");

  return;
}

/*----------------------------------------------------------------------*
 |  Compute lin for T and E matrix -- discr. wear            farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::Integrator::gp_2_d_te_lin(int& iter, Mortar::Element& sele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseMatrix& sderiv, Core::LinAlg::SerialDenseMatrix& lmderiv, double& jac,
    double& wgt, double* jumpval, const std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const Core::Gen::Pairedvector<int, double>& derivjac,
    const Core::Gen::Pairedvector<int, double>& dsliptmatrixgp,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap)
{
  const int nrow = sele.num_node();

  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();
  if (!snodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(snodes[iter]);
  if (!mymrtrnode) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  if (wear_shape_fcn() == Inpar::Wear::wear_shape_standard)
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
        for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                 dualmap.begin();
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
  else if (wear_shape_fcn() ==
           Inpar::Wear::wear_shape_dual)  //******************************************
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
        for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                 dualmap.begin();
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
        for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                 dualmap.begin();
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
        for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                 dualmap.begin();
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
    FOUR_C_THROW("Chosen shape functions not supported!");

  return;
}

/*----------------------------------------------------------------------*
 |  Compute lin for T and E matrix -- discr. wear            farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::Integrator::gp_3_d_te_lin(int& iter, Mortar::Element& sele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseMatrix& sderiv, Core::LinAlg::SerialDenseMatrix& lmderiv, double& jac,
    double& wgt, double* jumpval, const std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const Core::Gen::Pairedvector<int, double>& jacintcellmap,
    const Core::Gen::Pairedvector<int, double>& dsliptmatrixgp,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap)
{
  const int nrow = sele.num_node();

  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();
  if (!snodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(snodes[iter]);
  if (!mymrtrnode) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  if (wear_shape_fcn() == Inpar::Wear::wear_shape_standard)
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
          for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                   dualmap.begin();
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
  else if (wear_shape_fcn() == Inpar::Wear::wear_shape_dual)
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
          for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                   dualmap.begin();
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
          for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                   dualmap.begin();
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
          for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                   dualmap.begin();
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
    FOUR_C_THROW("Choosen shapefunctions not supported!");

  return;
}

/*----------------------------------------------------------------------*
 |  Compute lin for T and E matrix -- discr. wear (master)   farah 11/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::Integrator::gp_3_d_te_master_lin(int& iter, Mortar::Element& sele,
    Mortar::Element& mele, Core::LinAlg::SerialDenseVector& sval,
    Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseVector& lm2val, Core::LinAlg::SerialDenseMatrix& sderiv,
    Core::LinAlg::SerialDenseMatrix& mderiv, Core::LinAlg::SerialDenseMatrix& lmderiv,
    Core::LinAlg::SerialDenseMatrix& lm2deriv, double& jac, double& wgt, double* jumpval,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dmxigp,
    const Core::Gen::Pairedvector<int, double>& jacintcellmap,
    const Core::Gen::Pairedvector<int, double>& dsliptmatrixgp,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dual2map,
    const Epetra_Comm& comm)
{
  if (sele.Owner() != comm.MyPID()) return;

  const int ncol = mele.num_node();
  const int nrow = sele.num_node();

  // get slave element nodes themselves
  Core::Nodes::Node** mnodes = mele.Nodes();
  if (!mnodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  Mortar::Node* mymrtrnode = dynamic_cast<Mortar::Node*>(mnodes[iter]);
  if (!mymrtrnode) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  if (wear_shape_fcn() == Inpar::Wear::wear_shape_standard)
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
          for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                   dualmap.begin();
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
  else if (wear_shape_fcn() == Inpar::Wear::wear_shape_dual)
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
          for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                   dualmap.begin();
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
          for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
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
          for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
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
    FOUR_C_THROW("Chosen shape functions not supported!");

  return;
}

/*----------------------------------------------------------------------*
 |  Compute entries for weighted slip increment at GP        farah 09/13|
 *----------------------------------------------------------------------*/
void inline CONTACT::Integrator::gp_2_d_slip_incr(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& mval,
    Core::LinAlg::SerialDenseVector& lmval, Core::LinAlg::SerialDenseMatrix& sderiv,
    Core::LinAlg::SerialDenseMatrix& mderiv, double& jac, double& wgt, double* jumpvalv,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dmxigp,
    Core::Gen::Pairedvector<int, double>& dslipgp, int& linsize)
{
  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();
  Core::Nodes::Node** mnodes = mele.Nodes();
  if (!snodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");
  if (!mnodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  // number of nodes (slave, master)
  const int nrow = sele.num_node();
  const int ncol = mele.num_node();
  const int ndof = Dim();

  // LIN OF TANGENT
  Core::Gen::Pairedvector<int, double> dmap_txsl_gp(ncol * ndof + linsize);
  Core::Gen::Pairedvector<int, double> dmap_tysl_gp(ncol * ndof + linsize);

  // build interpolation of slave GP normal and coordinates
  std::array<double, 2> sjumpv = {0.0, 0.0};
  std::array<double, 2> mjumpv = {0.0, 0.0};
  std::array<double, 2> jumpv = {0.0, 0.0};
  std::array<double, 2> tanv = {0.0, 0.0};

  double tanlength = 0.0;
  for (int i = 0; i < nrow; ++i)
  {
    CONTACT::Node* myconode = dynamic_cast<CONTACT::Node*>(snodes[i]);

    // nodal tangent interpolation
    tanv[0] += sval[i] * myconode->Data().txi()[0];
    tanv[1] += sval[i] * myconode->Data().txi()[1];

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
  if (tanlength < 1.0e-12) FOUR_C_THROW("IntegrateAndDerivSegment: Divide by zero!");

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
    Core::Gen::Pairedvector<int, double>& dmap_txsl_i =
        dynamic_cast<CONTACT::Node*>(snodes[i])->Data().GetDerivTxi()[0];
    Core::Gen::Pairedvector<int, double>& dmap_tysl_i =
        dynamic_cast<CONTACT::Node*>(snodes[i])->Data().GetDerivTxi()[1];

    for (_CI p = dmap_txsl_i.begin(); p != dmap_txsl_i.end(); ++p)
      dmap_txsl_gp[p->first] += sval[i] * (p->second);
    for (_CI p = dmap_tysl_i.begin(); p != dmap_tysl_i.end(); ++p)
      dmap_tysl_gp[p->first] += sval[i] * (p->second);

    for (_CI p = dsxigp[0].begin(); p != dsxigp[0].end(); ++p)
    {
      double valx = sderiv(i, 0) * dynamic_cast<CONTACT::Node*>(snodes[i])->Data().txi()[0];
      dmap_txsl_gp[p->first] += valx * (p->second);
      double valy = sderiv(i, 0) * dynamic_cast<CONTACT::Node*>(snodes[i])->Data().txi()[1];
      dmap_tysl_gp[p->first] += valy * (p->second);
    }
  }

  // build directional derivative of slave GP tagent (unit)
  Core::Gen::Pairedvector<int, double> dmap_txsl_gp_unit(ncol * ndof + linsize);
  Core::Gen::Pairedvector<int, double> dmap_tysl_gp_unit(ncol * ndof + linsize);

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
void inline CONTACT::Integrator::gp_3_d_slip_incr(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& mval,
    Core::LinAlg::SerialDenseVector& lmval, Core::LinAlg::SerialDenseMatrix& sderiv,
    Core::LinAlg::SerialDenseMatrix& mderiv, double& jac, double& wgt, double* jumpvalv,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dmxigp,
    std::vector<Core::Gen::Pairedvector<int, double>>& dslipgp)
{
  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();
  Core::Nodes::Node** mnodes = mele.Nodes();
  if (!snodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");
  if (!mnodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  // number of nodes (slave, master)
  const int nrow = sele.num_node();
  const int ncol = mele.num_node();
  const int ndof = Dim();

  int linsize = 0;
  for (int i = 0; i < nrow; ++i)
  {
    Node* cnode = dynamic_cast<Node*>(snodes[i]);
    linsize += cnode->GetLinsize();
  }

  // build interpolation of slave GP normal and coordinates
  std::array<double, 3> sjumpv = {0.0, 0.0, 0.0};
  std::array<double, 3> mjumpv = {0.0, 0.0, 0.0};
  std::array<double, 3> jumpv = {0.0, 0.0, 0.0};
  std::array<double, 3> tanv1 = {0.0, 0.0, 0.0};
  std::array<double, 3> tanv2 = {0.0, 0.0, 0.0};

  double jumpvalv1 = 0.0;
  double jumpvalv2 = 0.0;
  double tanlength1 = 0.0;
  double tanlength2 = 0.0;

  // LIN OF TANGENT
  std::map<int, double> dmap_txsl_gp;
  std::map<int, double> dmap_tysl_gp;

  for (int i = 0; i < nrow; ++i)
  {
    CONTACT::Node* myconode = dynamic_cast<CONTACT::Node*>(snodes[i]);

    // nodal tangent interpolation
    tanv1[0] += sval[i] * myconode->Data().txi()[0];
    tanv1[1] += sval[i] * myconode->Data().txi()[1];
    tanv1[2] += sval[i] * myconode->Data().txi()[2];

    tanv2[0] += sval[i] * myconode->Data().teta()[0];
    tanv2[1] += sval[i] * myconode->Data().teta()[1];
    tanv2[2] += sval[i] * myconode->Data().teta()[2];
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
  if (tanlength1 < 1.0e-12) FOUR_C_THROW("IntegrateAndDerivSegment: Divide by zero!");

  tanlength2 = sqrt(tanv2[0] * tanv2[0] + tanv2[1] * tanv2[1] + tanv2[2] * tanv2[2]);
  if (tanlength2 < 1.0e-12) FOUR_C_THROW("IntegrateAndDerivSegment: Divide by zero!");

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
  Core::Gen::Pairedvector<int, double> dmap_txix_gp(linsize + ncol * ndof);
  Core::Gen::Pairedvector<int, double> dmap_txiy_gp(linsize + ncol * ndof);
  Core::Gen::Pairedvector<int, double> dmap_txiz_gp(linsize + ncol * ndof);

  // slave GP txi (non-unit)
  Core::Gen::Pairedvector<int, double> dmap_txix_gp_unit(linsize + ncol * ndof);
  Core::Gen::Pairedvector<int, double> dmap_txiy_gp_unit(linsize + ncol * ndof);
  Core::Gen::Pairedvector<int, double> dmap_txiz_gp_unit(linsize + ncol * ndof);

  for (int i = 0; i < nrow; ++i)
  {
    Core::Gen::Pairedvector<int, double>& dmap_txsl_i =
        dynamic_cast<CONTACT::Node*>(snodes[i])->Data().GetDerivTxi()[0];
    Core::Gen::Pairedvector<int, double>& dmap_tysl_i =
        dynamic_cast<CONTACT::Node*>(snodes[i])->Data().GetDerivTxi()[1];
    Core::Gen::Pairedvector<int, double>& dmap_tzsl_i =
        dynamic_cast<CONTACT::Node*>(snodes[i])->Data().GetDerivTxi()[2];

    for (_CI p = dmap_txsl_i.begin(); p != dmap_txsl_i.end(); ++p)
      dmap_txix_gp[p->first] += sval[i] * (p->second);
    for (_CI p = dmap_tysl_i.begin(); p != dmap_tysl_i.end(); ++p)
      dmap_txiy_gp[p->first] += sval[i] * (p->second);
    for (_CI p = dmap_tzsl_i.begin(); p != dmap_tzsl_i.end(); ++p)
      dmap_txiz_gp[p->first] += sval[i] * (p->second);

    const double txi_x = dynamic_cast<CONTACT::Node*>(snodes[i])->Data().txi()[0];
    const double txi_y = dynamic_cast<CONTACT::Node*>(snodes[i])->Data().txi()[1];
    const double txi_z = dynamic_cast<CONTACT::Node*>(snodes[i])->Data().txi()[2];

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
  Core::Gen::Pairedvector<int, double> dmap_tetax_gp(linsize + ncol * ndof);
  Core::Gen::Pairedvector<int, double> dmap_tetay_gp(linsize + ncol * ndof);
  Core::Gen::Pairedvector<int, double> dmap_tetaz_gp(linsize + ncol * ndof);

  // slave GP teta (unit)
  Core::Gen::Pairedvector<int, double> dmap_tetax_gp_unit(linsize + ncol * ndof);
  Core::Gen::Pairedvector<int, double> dmap_tetay_gp_unit(linsize + ncol * ndof);
  Core::Gen::Pairedvector<int, double> dmap_tetaz_gp_unit(linsize + ncol * ndof);

  for (int i = 0; i < nrow; ++i)
  {
    Core::Gen::Pairedvector<int, double>& dmap_txsl_i =
        dynamic_cast<CONTACT::Node*>(snodes[i])->Data().GetDerivTeta()[0];
    Core::Gen::Pairedvector<int, double>& dmap_tysl_i =
        dynamic_cast<CONTACT::Node*>(snodes[i])->Data().GetDerivTeta()[1];
    Core::Gen::Pairedvector<int, double>& dmap_tzsl_i =
        dynamic_cast<CONTACT::Node*>(snodes[i])->Data().GetDerivTeta()[2];

    for (_CI p = dmap_txsl_i.begin(); p != dmap_txsl_i.end(); ++p)
      dmap_tetax_gp[p->first] += sval[i] * (p->second);
    for (_CI p = dmap_tysl_i.begin(); p != dmap_tysl_i.end(); ++p)
      dmap_tetay_gp[p->first] += sval[i] * (p->second);
    for (_CI p = dmap_tzsl_i.begin(); p != dmap_tzsl_i.end(); ++p)
      dmap_tetaz_gp[p->first] += sval[i] * (p->second);

    const double teta_x = dynamic_cast<CONTACT::Node*>(snodes[i])->Data().teta()[0];
    const double teta_y = dynamic_cast<CONTACT::Node*>(snodes[i])->Data().teta()[1];
    const double teta_z = dynamic_cast<CONTACT::Node*>(snodes[i])->Data().teta()[2];

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
void inline CONTACT::Integrator::gp_2_d_slip_incr_lin(int& iter, Mortar::Element& sele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseMatrix& sderiv, Core::LinAlg::SerialDenseMatrix& lmderiv, double& jac,
    double& wgt, double* jumpvalv, const std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const Core::Gen::Pairedvector<int, double>& dslipgp,
    const Core::Gen::Pairedvector<int, double>& derivjac,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap)
{
  Core::Nodes::Node** snodes = sele.Nodes();

  const int nrow = sele.num_node();
  double fac = 0.0;

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  FriNode* snode = dynamic_cast<FriNode*>(snodes[iter]);
  if (snode->IsOnBoundorCE()) return;

  // get the corresponding map as a reference
  std::map<int, double>& djumpmap = snode->FriData().GetDerivVarJump()[0];

  // (1) Lin(Phi) - dual shape functions
  if (shape_fcn() == Inpar::Mortar::shape_dual)
  {
    for (int m = 0; m < nrow; ++m)
    {
      fac = wgt * sval[m] * jumpvalv[0] * jac;
      for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
               dualmap.begin();
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
void inline CONTACT::Integrator::gp_3_d_slip_incr_lin(int& iter, Mortar::Element& sele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseMatrix& sderiv, Core::LinAlg::SerialDenseMatrix& lmderiv, double& jac,
    double& wgt, double* jumpvalv, const Core::Gen::Pairedvector<int, double>& jacintcellmap,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dslipgp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap)
{
  Core::Nodes::Node** snodes = sele.Nodes();

  double nrow = sele.num_node();

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  FriNode* snode = dynamic_cast<FriNode*>(snodes[iter]);

  if (snode->IsOnBoundorCE()) return;

  // get the corresponding map as a reference
  std::map<int, double>& djumpmap1 = snode->FriData().GetDerivVarJump()[0];
  std::map<int, double>& djumpmap2 = snode->FriData().GetDerivVarJump()[1];

  double fac1 = 0.0;
  double fac2 = 0.0;

  // (1) Lin(Phi) - dual shape functions
  if (shape_fcn() == Inpar::Mortar::shape_dual)
  {
    for (int m = 0; m < nrow; ++m)
    {
      fac1 = wgt * sval[m] * jumpvalv[0] * jac;
      fac2 = wgt * sval[m] * jumpvalv[1] * jac;

      for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
               dualmap.begin();
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
void inline CONTACT::Integrator::gp_2_d_wear_lin(int& iter, Mortar::Element& sele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseMatrix& sderiv, Core::LinAlg::SerialDenseMatrix& lmderiv, double& jac,
    double* gpn, double& wgt, double& wearval, double* jumpval,
    const Core::Gen::Pairedvector<int, double>& dweargp,
    const Core::Gen::Pairedvector<int, double>& derivjac,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap)
{
  const double wcoeff = wearcoeff_ + wearcoeffm_;

  double facw = 0.0;
  const int nrow = sele.num_node();

  Core::Nodes::Node** snodes = sele.Nodes();

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // get the corresponding map as a reference
  CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(snodes[iter]);

  std::map<int, double>& dwmap = cnode->Data().GetDerivW();

  if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
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
    if (shape_fcn() == Inpar::Mortar::shape_dual)
    {
      for (int m = 0; m < nrow; ++m)
      {
        facw = wcoeff * wgt * sval[m] * wearval * jac;
        for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                 dualmap.begin();
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
  std::map<int, double>& dwlmmap = cnode->Data().GetDerivWlm();

  for (int bl = 0; bl < nrow; ++bl)
  {
    Mortar::Node* wearnode = dynamic_cast<Mortar::Node*>(snodes[bl]);

    if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
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
void inline CONTACT::Integrator::gp_3_d_wear_lin(int& iter, Mortar::Element& sele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseMatrix& sderiv, Core::LinAlg::SerialDenseMatrix& lmderiv, double& jac,
    double* gpn, double& wgt, double& wearval, double* jumpval,
    const Core::Gen::Pairedvector<int, double>& dweargp,
    const Core::Gen::Pairedvector<int, double>& jacintcellmap,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap)
{
  double facw = 0.0;
  const int nrow = sele.num_node();

  Core::Nodes::Node** snodes = sele.Nodes();

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;

  // get the corresponding map as a reference
  CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(snodes[iter]);

  // get the corresponding map as a reference
  std::map<int, double>& dwmap = cnode->Data().GetDerivW();

  if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
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
    if (shape_fcn() == Inpar::Mortar::shape_dual)
    {
      for (int m = 0; m < nrow; ++m)
      {
        facw = wearcoeff_ * wgt * sval[m] * wearval * jac;
        for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                 dualmap.begin();
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
  std::map<int, double>& dwlmmap = cnode->Data().GetDerivWlm();

  for (int bl = 0; bl < nrow; ++bl)
  {
    Mortar::Node* wearnode = dynamic_cast<Mortar::Node*>(snodes[bl]);

    if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
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
void inline CONTACT::Integrator::gp_ncoup_deriv(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& mval,
    Core::LinAlg::SerialDenseVector& lmval, Core::LinAlg::SerialDenseMatrix& sderiv,
    Core::LinAlg::SerialDenseMatrix& mderiv, double* ncoup, double* gpn, double& jac, double& wgt,
    double* gpcoord, const std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dmxigp,
    std::map<int, double>& dncoupgp, std::map<int, double>& dvelncoupgp,
    std::map<int, double>& dpresncoupgp,
    std::vector<Core::Gen::Pairedvector<int, double>>& dnmap_unit, bool quadratic, int nintrow)
{
  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();
  Core::Nodes::Node** mnodes = mele.Nodes();
  if (!snodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");
  if (!mnodes) FOUR_C_THROW("IntegrateAndDerivSegment: Null pointer!");

  // bool to decide if fluid quantities and structural velocities
  // exist for slave side on Node level and contribute to porofluid meshtying
  bool slaveporo = (sele.PhysType() == Mortar::Element::poro);
  // bool to decide if fluid quantities and structural velocities
  // exist for master side on Node level and contribute to porofluid meshtying
  bool masterporo = (mele.PhysType() == Mortar::Element::poro);

  if (!slaveporo and masterporo)
    FOUR_C_THROW(
        "poroelastic mesh tying method needs the slave side to be poroelastic (invert master/slave "
        "definition)");
  // see CONTACT::LagrangeStrategyPoro::LagrangeStrategyPoro for comment

  // number of nodes (slave, master)
  int nrow = sele.num_node();
  // if(twosided) //needed only for master side
  int ncol = mele.num_node();  // not used for onesided porocontact!!!

  // get fluid velocities in GP
  std::array<double, 3> sgpfvel = {0.0, 0.0, 0.0};
  std::array<double, 3> sgpsvel = {0.0, 0.0, 0.0};

  std::array<double, 3> mgpfvel = {0.0, 0.0, 0.0};
  std::array<double, 3> mgpsvel = {0.0, 0.0, 0.0};

  // fill local vectors for corresponding problem dimension
  for (int k = 0; k < Dim(); ++k)
  {
    if (slaveporo)
      for (int i = 0; i < nrow; ++i)
      {
        Node* mymrtrnode = dynamic_cast<Node*>(snodes[i]);
        sgpfvel[k] += sval[i] * mymrtrnode->PoroData().fvel()[k];
        sgpsvel[k] += sval[i] * mymrtrnode->PoroData().svel()[k];
      }

    if (masterporo)
    {
      for (int i = 0; i < ncol; ++i)
      {
        Node* mymrtrnode = dynamic_cast<Node*>(mnodes[i]);
        mgpfvel[k] += mval[i] * mymrtrnode->PoroData().fvel()[k];
        mgpsvel[k] += mval[i] * mymrtrnode->PoroData().svel()[k];
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
    const double sJ = det_deformation_gradient(sele, wgt, gpcoord, sJLin);

    for (int i = 0; i < nrow; ++i)
    {
      Node* mymrtrnode = dynamic_cast<Node*>(snodes[i]);
      sgpfpres += sval[i] * mymrtrnode->PoroData().fpres()[0];
    }
    Teuchos::ParameterList sparams;  // empty parameter list;

    Teuchos::RCP<Mat::StructPoro> sstructmat =
        Teuchos::rcp_dynamic_cast<Mat::StructPoro>(sele.parent_element()->Material(0));
    if (sstructmat == Teuchos::null)
      sstructmat = Teuchos::rcp_dynamic_cast<Mat::StructPoro>(sele.parent_element()->Material(1));
    if (sstructmat == Teuchos::null) FOUR_C_THROW("Cast to StructPoro failed!");
    sstructmat->ComputeSurfPorosity(sparams, sgpfpres, sJ, sele.FaceParentNumber(),
        1,  // finally check what to do here Todo:
        sporosity,
        &sdphi_dp,  // linearization of phi w.r.t. pressure
        &sdphi_dJ,  // linearization of phi w.r.t. defo-grad determinant
        nullptr,    // dphi_dJdp not needed
        nullptr,    // dphi_dJJ not needed
        nullptr,    // dphi_dpp not needed
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
    const double mJ = det_deformation_gradient(mele, wgt, gpcoord, mJLin);

    for (int i = 0; i < ncol; ++i)
    {
      Node* mymrtrnode = dynamic_cast<Node*>(mnodes[i]);
      mgpfpres += mval[i] * mymrtrnode->PoroData().fpres()[0];
    }
    Teuchos::ParameterList mparams;  // empty parameter list;

    Teuchos::RCP<Mat::StructPoro> mstructmat =
        Teuchos::rcp_dynamic_cast<Mat::StructPoro>(mele.parent_element()->Material(0));
    if (mstructmat == Teuchos::null)
      mstructmat = Teuchos::rcp_dynamic_cast<Mat::StructPoro>(mele.parent_element()->Material(1));
    if (mstructmat == Teuchos::null) FOUR_C_THROW("Cast to StructPoro failed!");

    mstructmat->ComputeSurfPorosity(mparams, mgpfpres, mJ,
        mele.FaceParentNumber(),  // may not work
        1,                        // finally check what to do here Todo:
        mporosity,
        &mdphi_dp,  // linearization of phi w.r.t. pressure
        &mdphi_dJ,  // linearization of phi w.r.t. defo-grad determinant
        nullptr,    // dphi_dJdp not needed
        nullptr,    // dphi_dJJ not needed
        nullptr,    // dphi_dpp not needed
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
      CONTACT::Node* mrtrnode = dynamic_cast<CONTACT::Node*>(snodes[j]);

      double prod = 0.0;
      // Petrov-Galerkin approach (dual LM for D/M but standard LM for normal coupling)
      if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin) prod = sval[j] * ncoup[0] * jac * wgt;
      // usual standard or dual LM approach
      else
        prod = lmval[j] * ncoup[0] * jac * wgt;

      // do not process slave side boundary nodes
      // (their row entries would be zero anyway!)
      if (mrtrnode->IsOnBound())
      {
        FOUR_C_THROW("CHECKME: isonbound_ active; should not happen here");
        continue;
      }

      // add current Gauss point's contribution to gseg
      mrtrnode->AddNcoupValue(prod);
    }
  }

  //    // CASE 4: Dual LM shape functions and quadratic interpolation
  //    else if ((ShapeFcn() == Inpar::Mortar::shape_dual || ShapeFcn() ==
  //    Inpar::Mortar::shape_petrovgalerkin) &&
  //        LagMultQuad() == Inpar::Mortar::lagmult_quad)
  //    {
  //      for (int j=0;j<nrow;++j)
  //      {
  //        CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(snodes[j]);
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
    FOUR_C_THROW("Invalid integration case for 3D quadratic normal coupling mortar!");
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
      Node* mrtrnode = dynamic_cast<Node*>(snodes[i]);

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
                sporosity * gpn[k] * sderiv(i, d) * mrtrnode->PoroData().fvel()[k] * (p->second);
            dncoupgp[p->first] +=
                sporosity * gpn[k] * sderiv(i, d) * mrtrnode->PoroData().svel()[k] * (p->second);
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
      Node* mrtrnode = dynamic_cast<Node*>(mnodes[i]);

      for (int k = 0; k < Dim(); ++k)
      {
        dncoupgp[mrtrnode->Dofs()[k]] -= mporosity * mval[i] * gpn[k] * timefac;
        for (unsigned int d = 0; d < dmxigp.size(); ++d)
        {
          for (_CI p = dmxigp[d].begin(); p != dmxigp[d].end(); ++p)
          {
            dncoupgp[p->first] +=
                mporosity * gpn[k] * mderiv(i, d) * mrtrnode->PoroData().fvel()[k] * (p->second);
            dncoupgp[p->first] -=
                mporosity * gpn[k] * mderiv(i, d) * mrtrnode->PoroData().svel()[k] * (p->second);
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
      Node* mrtrnode = dynamic_cast<Node*>(snodes[z]);

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
      Node* mrtrnode = dynamic_cast<Node*>(mnodes[z]);

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
void inline CONTACT::Integrator::gp_ncoup_lin(int& iter, Mortar::Element& sele,
    Mortar::Element& mele, Core::LinAlg::SerialDenseVector& sval,
    Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseMatrix& sderiv, Core::LinAlg::SerialDenseMatrix& lmderiv,
    double& ncoup, double* gpn, double& jac, double& wgt, const std::map<int, double>& dncoupgp,
    const std::map<int, double>& dvelncoupgp, const std::map<int, double>& dpresncoupgp,
    const Core::Gen::Pairedvector<int, double>& jacintcellmap,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dmxigp,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap)
{
  int nrow = sele.num_node();
  // only sele.num_node() is needed because the integration of the interface constraint is
  // evaluated on the slave side interface surface

  // map iterator
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;
  typedef std::map<int, double>::const_iterator CI;

  // get slave element nodes themselves
  Core::Nodes::Node** snodes = sele.Nodes();
  // Core::Nodes::Node** mnodes = mele.Nodes();

  CONTACT::Node* mymrtrnode = dynamic_cast<CONTACT::Node*>(snodes[iter]);
  if (!mymrtrnode) FOUR_C_THROW("CONTACT::Integrator::gp_ncoup_lin: mymrtnode: Null pointer!");

  double fac = 0.0;

  // get the corresponding map as a reference
  std::map<int, double>& dgmap = mymrtrnode->PoroData().GetDerivnCoup();
  // switch if Petrov-Galerkin approach for LM is applied
  // for meshtying shape_dual is the only allowed case as long as porofluid and skeleton
  // meshtying conditions share the same input lines
  //(for contact they still share the input lines but shape_petrovgalerkin is allowed for
  // skeleton/structural contact) shape_petrovgalerkin
  if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
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
        for (Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>::const_iterator p =
                 dualmap.begin();
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
  std::map<int, double>& dvelncoupmap = mymrtrnode->PoroData().GetVelDerivnCoup();
  // switch if Petrov-Galerkin approach for LM is applied
  if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
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
  std::map<int, double>& dpresncoupmap = mymrtrnode->PoroData().GetPresDerivnCoup();
  // switch if Petrov-Galerkin approach for LM is applied
  if (shape_fcn() == Inpar::Mortar::shape_petrovgalerkin)
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
double CONTACT::Integrator::det_deformation_gradient(
    Mortar::Element& sele, double& wgt, double* gpcoord, std::map<int, double>& JLin)
{
  double J;
  Core::FE::CellType distype = sele.parent_element()->Shape();
  switch (distype)
  {
    case Core::FE::CellType::hex8:
      J = t_det_deformation_gradient<Core::FE::CellType::hex8, 3>(sele, wgt, gpcoord, JLin);
      break;
    case Core::FE::CellType::quad4:
      J = t_det_deformation_gradient<Core::FE::CellType::quad4, 2>(sele, wgt, gpcoord, JLin);
      break;
    default:
      FOUR_C_THROW(
          "DetDeformationGradient3D: Parent Element Type not templated yet, just add it here!");
      J = 0.0;  // To avoid warning!
      break;
  }

  return J;
}
/*-------------------------------------------------------------------------------*
 |  Templated Calculate Determinant of the Deformation Gradient at GP  ager 10/14|
 |  modified by h.Willmann 2015                                                  |
 *-------------------------------------------------------------------------------*/
template <Core::FE::CellType parentdistype, int dim>
double CONTACT::Integrator::t_det_deformation_gradient(
    Mortar::Element& sele, double& wgt, double* gpcoord, std::map<int, double>& JLin)
{
  //! nen_: number of element nodes (T. Hughes: The Finite Element Method)
  static const int numnodes = Core::FE::num_nodes<parentdistype>;

  Core::FE::CollectedGaussPoints intpoints =
      Core::FE::CollectedGaussPoints(1);  // reserve just for 1 entry ...
  intpoints.Append(gpcoord[0], gpcoord[1], 0.0, wgt);

  // get coordinates of gauss point w.r.t. local parent coordinate system
  Core::LinAlg::SerialDenseMatrix pqxg(1, dim);
  Core::LinAlg::Matrix<dim, dim> derivtrafo(true);

  Core::FE::BoundaryGPToParentGP<dim>(pqxg, derivtrafo, intpoints, sele.parent_element()->Shape(),
      sele.Shape(), sele.FaceParentNumber());

  Core::LinAlg::Matrix<dim, 1> pxsi(true);

  // coordinates of the current integration point in parent coordinate system
  for (int idim = 0; idim < dim; idim++)
  {
    pxsi(idim) = pqxg(0, idim);
  }

  Core::LinAlg::Matrix<dim, numnodes> pderiv_loc(
      true);  // derivatives of parent element shape functions in parent element coordinate system

  // evaluate derivatives of parent element shape functions at current integration point in parent
  // coordinate system
  Core::FE::shape_function_deriv1<parentdistype>(pxsi, pderiv_loc);
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
  Core::LinAlg::Matrix<dim, dim> xjm;
  Core::LinAlg::Matrix<dim, dim> Jmat;

  Core::LinAlg::Matrix<dim, numnodes> xrefe(true);  // material coord. of parent element
  Core::LinAlg::Matrix<dim, numnodes> xcurr(true);  // current  coord. of parent element

  // update element geometry of parent element
  {
    Core::Nodes::Node** nodes = sele.parent_element()->Nodes();
    for (int inode = 0; inode < numnodes; ++inode)
    {
      for (unsigned int idof = 0; idof < dim; ++idof)
      {
        const auto& x = nodes[inode]->X();
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

  //  Core::LinAlg::Matrix<dim,numnodes>  auxJLin (true);   // matrix to be filled with
  //  linearization scalars
  xjm.Invert();
  {
    //    Core::Nodes::Node** nodes = sele.parent_element()->Nodes();
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

FOUR_C_NAMESPACE_CLOSE
