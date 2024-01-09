/*-----------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of Mortar matrices on the overlap
of two MORTAR::Elements in 1D and 2D

\level 1

*/
/*----------------------------------------------------------------------*/

#include "baci_mortar_integrator.H"

#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_lib_element.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_linalg_serialdensevector.H"
#include "baci_linalg_sparsematrix.H"
#include "baci_mortar_calc_utils.H"
#include "baci_mortar_coupling3d_classes.H"
#include "baci_mortar_defines.H"
#include "baci_mortar_element.H"
#include "baci_mortar_node.H"
#include "baci_mortar_projector.H"
#include "baci_mortar_shape_utils.H"
#include "baci_utils_singleton_owner.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  impl...                                                  farah 01/14|
 *----------------------------------------------------------------------*/
MORTAR::Integrator* MORTAR::Integrator::Impl(
    MORTAR::Element& sele, MORTAR::Element& mele, Teuchos::ParameterList& params)
{
  switch (sele.Shape())
  {
    // 2D surface elements
    case CORE::FE::CellType::quad4:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::quad4:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::quad4,
              CORE::FE::CellType::quad4>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::quad8:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::quad4,
              CORE::FE::CellType::quad8>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::quad9:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::quad4,
              CORE::FE::CellType::quad9>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::tri3:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::quad4,
              CORE::FE::CellType::tri3>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::tri6:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::quad4,
              CORE::FE::CellType::tri6>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        default:
          dserror("Element combination not allowed!");
      }
      break;
    }
    case CORE::FE::CellType::quad8:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::quad4:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::quad8,
              CORE::FE::CellType::quad4>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::quad8:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::quad8,
              CORE::FE::CellType::quad8>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::quad9:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::quad8,
              CORE::FE::CellType::quad9>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::tri3:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::quad8,
              CORE::FE::CellType::tri3>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::tri6:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::quad8,
              CORE::FE::CellType::tri6>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        default:
          dserror("Element combination not allowed!");
      }
      break;
    }
    case CORE::FE::CellType::quad9:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::quad4:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::quad9,
              CORE::FE::CellType::quad4>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::quad8:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::quad9,
              CORE::FE::CellType::quad8>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::quad9:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::quad9,
              CORE::FE::CellType::quad9>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::tri3:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::quad9,
              CORE::FE::CellType::tri3>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::tri6:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::quad9,
              CORE::FE::CellType::tri6>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        default:
          dserror("Element combination not allowed!");
      }
      break;
    }
    case CORE::FE::CellType::tri3:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::quad4:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::tri3,
              CORE::FE::CellType::quad4>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::quad8:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::tri3,
              CORE::FE::CellType::quad8>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::quad9:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::tri3,
              CORE::FE::CellType::quad9>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::tri3:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::tri3,
              CORE::FE::CellType::tri3>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::tri6:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::tri3,
              CORE::FE::CellType::tri6>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        default:
          dserror("Element combination not allowed!");
      }
      break;
    }
    case CORE::FE::CellType::tri6:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::quad4:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::tri6,
              CORE::FE::CellType::quad4>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::quad8:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::tri6,
              CORE::FE::CellType::quad8>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::quad9:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::tri6,
              CORE::FE::CellType::quad9>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::tri3:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::tri6,
              CORE::FE::CellType::tri3>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::tri6:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::tri6,
              CORE::FE::CellType::tri6>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        default:
          dserror("Element combination not allowed!");
      }
      break;
    }
      // 1D surface elements
    case CORE::FE::CellType::line2:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::line2:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::line2,
              CORE::FE::CellType::line2>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::line3:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::line2,
              CORE::FE::CellType::line3>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        default:
          dserror("Element combination not allowed!");
      }
      break;
    }
    case CORE::FE::CellType::line3:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::line2:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::line3,
              CORE::FE::CellType::line2>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::line3:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::line3,
              CORE::FE::CellType::line3>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        default:
          dserror("Element combination not allowed!");
      }
      break;
    }

      //==================================================
      //                     NURBS
      //==================================================
      // 1D surface elements
    case CORE::FE::CellType::nurbs2:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::nurbs2:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::nurbs2,
              CORE::FE::CellType::nurbs2>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::nurbs3:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::nurbs2,
              CORE::FE::CellType::nurbs3>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        default:
          dserror("Element combination not allowed!");
      }
      break;
    }
    case CORE::FE::CellType::nurbs3:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::nurbs2:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::nurbs3,
              CORE::FE::CellType::nurbs2>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::nurbs3:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::nurbs3,
              CORE::FE::CellType::nurbs3>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        default:
          dserror("Element combination not allowed!");
      }
      break;
    }
    case CORE::FE::CellType::nurbs9:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::nurbs9:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::nurbs9,
              CORE::FE::CellType::nurbs9>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        case CORE::FE::CellType::nurbs4:
        {
          return MORTAR::IntegratorCalc<CORE::FE::CellType::nurbs9,
              CORE::FE::CellType::nurbs4>::Instance(CORE::UTILS::SingletonAction::create, params);
        }
        default:
          dserror("Element combination not allowed!");
      }
      break;
    }
    default:
      dserror("Error...");
      break;
  }
  return nullptr;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 01/14|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
MORTAR::IntegratorCalc<distypeS, distypeM>::IntegratorCalc(const Teuchos::ParameterList& params)
    : imortar_(params),
      shapefcn_(INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(params, "LM_SHAPEFCN")),
      lmquadtype_(INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(params, "LM_QUAD"))
{
  InitializeGP();
}

template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
MORTAR::IntegratorCalc<distypeS, distypeM>* MORTAR::IntegratorCalc<distypeS, distypeM>::Instance(
    CORE::UTILS::SingletonAction action, const Teuchos::ParameterList& params)
{
  static auto singleton_owner = CORE::UTILS::MakeSingletonOwner(
      [](const Teuchos::ParameterList& p)
      {
        return std::unique_ptr<MORTAR::IntegratorCalc<distypeS, distypeM>>(
            new MORTAR::IntegratorCalc<distypeS, distypeM>(p));
      });

  return singleton_owner.Instance(action, params);
}


/*----------------------------------------------------------------------*
 |  Initialize gauss points                                   popp 06/09|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
void MORTAR::IntegratorCalc<distypeS, distypeM>::InitializeGP()
{
  // get numgp (for element-based integration)
  int numgp = imortar_.get<int>("NUMGP_PER_DIM");

  // get integration type
  INPAR::MORTAR::IntType integrationtype =
      INPUT::IntegralValue<INPAR::MORTAR::IntType>(imortar_, "INTTYPE");

  // if we use segment-based integration, the shape of the cells has to be considered!
  CORE::FE::CellType intshape;
  if (integrationtype == INPAR::MORTAR::inttype_segments)
  {
    if (ndim_ == 2)
      intshape = CORE::FE::CellType::line2;
    else if (ndim_ == 3)
      intshape = CORE::FE::CellType::tri3;
    else
      dserror("wrong dimension!");
  }
  else
    intshape = distypeS;

  //**********************************************************************
  // choose Gauss rule according to (a) element type (b) input parameter
  //**********************************************************************
  switch (intshape)
  {
    case CORE::FE::CellType::line2:
    case CORE::FE::CellType::line3:
    case CORE::FE::CellType::nurbs2:
    case CORE::FE::CellType::nurbs3:
    {
      // set default value for segment-based version first
      CORE::FE::GaussRule1D mygaussrule = CORE::FE::GaussRule1D::line_5point;

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
              mygaussrule = CORE::FE::GaussRule1D::line_2point;
              break;
            }
            case 3:
            {
              mygaussrule = CORE::FE::GaussRule1D::line_3point;
              break;
            }
            case 4:
            {
              mygaussrule = CORE::FE::GaussRule1D::line_4point;
              break;
            }
            case 5:
            {
              mygaussrule = CORE::FE::GaussRule1D::line_5point;
              break;
            }
            case 6:
            {
              mygaussrule = CORE::FE::GaussRule1D::line_6point;
              break;
            }
            case 7:
            {
              mygaussrule = CORE::FE::GaussRule1D::line_7point;
              break;
            }
            case 8:
            {
              mygaussrule = CORE::FE::GaussRule1D::line_8point;
              break;
            }
            case 9:
            {
              mygaussrule = CORE::FE::GaussRule1D::line_9point;
              break;
            }
            case 10:
            {
              mygaussrule = CORE::FE::GaussRule1D::line_10point;
              break;
            }
            case 16:
            {
              mygaussrule = CORE::FE::GaussRule1D::line_16point;
              break;
            }
            case 20:
            {
              mygaussrule = CORE::FE::GaussRule1D::line_20point;
              break;
            }
            case 32:
            {
              mygaussrule = CORE::FE::GaussRule1D::line_32point;
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

      const CORE::FE::IntegrationPoints1D intpoints(mygaussrule);
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
    case CORE::FE::CellType::tri3:
    case CORE::FE::CellType::tri6:
    {
      // set default value for segment-based version first
      CORE::FE::GaussRule2D mygaussrule = CORE::FE::GaussRule2D::tri_7point;
      if (integrationtype == INPAR::MORTAR::inttype_segments)
      {
        if (numgp > 0) switch (numgp)
          {
            case 7:
              mygaussrule = CORE::FE::GaussRule2D::tri_7point;
              break;
            case 12:
              mygaussrule = CORE::FE::GaussRule2D::tri_12point;
              break;
            case 16:
              mygaussrule = CORE::FE::GaussRule2D::tri_16point;
              break;
            case 37:
              mygaussrule = CORE::FE::GaussRule2D::tri_37point;
              break;
            default:
              dserror("unknown tri gauss rule");
              break;
          }
      }
      // GP switch if element-based version and non-zero value provided by user
      else if (integrationtype == INPAR::MORTAR::inttype_elements ||
               integrationtype == INPAR::MORTAR::inttype_elements_BS)
      {
        if (numgp > 0)
        {
          switch (numgp)
          {
            case 1:
            {
              mygaussrule = CORE::FE::GaussRule2D::tri_3point;
              break;
            }
            case 2:
            {
              mygaussrule = CORE::FE::GaussRule2D::tri_6point;
              break;
            }
            case 3:
            {
              mygaussrule = CORE::FE::GaussRule2D::tri_7point;
              break;
            }
            case 4:
            {
              mygaussrule = CORE::FE::GaussRule2D::tri_12point;
              break;
            }
            case 5:
            {
              mygaussrule = CORE::FE::GaussRule2D::tri_12point;
              break;
            }
            case 6:
            {
              mygaussrule = CORE::FE::GaussRule2D::tri_37point;
              break;
            }
            case 7:
            {
              mygaussrule = CORE::FE::GaussRule2D::tri_37point;
              break;
            }
            case 8:
            {
              mygaussrule = CORE::FE::GaussRule2D::tri_64point;
              break;
            }
            case 9:
            {
              mygaussrule = CORE::FE::GaussRule2D::tri_64point;
              break;
            }
            case 10:
            {
              mygaussrule = CORE::FE::GaussRule2D::tri_64point;
              break;
            }
            case 20:
            {
              mygaussrule = CORE::FE::GaussRule2D::tri_64point;
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
      else
        dserror("unknown integration type!");

      const CORE::FE::IntegrationPoints2D intpoints(mygaussrule);
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
    case CORE::FE::CellType::quad4:
    case CORE::FE::CellType::quad8:
    case CORE::FE::CellType::quad9:
    case CORE::FE::CellType::nurbs4:
    case CORE::FE::CellType::nurbs8:
    case CORE::FE::CellType::nurbs9:
    {
      // set default value for segment-based version first
      CORE::FE::GaussRule2D mygaussrule = CORE::FE::GaussRule2D::quad_25point;

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
              mygaussrule = CORE::FE::GaussRule2D::quad_4point;
              break;
            }
            case 3:
            {
              mygaussrule = CORE::FE::GaussRule2D::quad_9point;
              break;
            }
            case 4:
            {
              mygaussrule = CORE::FE::GaussRule2D::quad_16point;
              break;
            }
            case 5:
            {
              mygaussrule = CORE::FE::GaussRule2D::quad_25point;
              break;
            }
            case 6:
            {
              mygaussrule = CORE::FE::GaussRule2D::quad_36point;
              break;
            }
            case 7:
            {
              mygaussrule = CORE::FE::GaussRule2D::quad_49point;
              break;
            }
            case 8:
            {
              mygaussrule = CORE::FE::GaussRule2D::quad_64point;
              break;
            }
            case 9:
            {
              mygaussrule = CORE::FE::GaussRule2D::quad_81point;
              break;
            }
            case 10:
            {
              mygaussrule = CORE::FE::GaussRule2D::quad_100point;
              break;
            }
            case 16:
            {
              mygaussrule = CORE::FE::GaussRule2D::quad_256point;
              break;
            }
            case 20:
            {
              mygaussrule = CORE::FE::GaussRule2D::quad_400point;
              break;
            }
            case 32:
            {
              mygaussrule = CORE::FE::GaussRule2D::quad_1024point;
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

      const CORE::FE::IntegrationPoints2D intpoints(mygaussrule);
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
      dserror("MORTAR::Integrator: This contact element type is not implemented!");
      break;
    }
  }  // switch(eletype)

  return;
}


/*--------------------------------------------------------------------------------------*
 | Integrate without segmentation --> more GP required                       farah 01/13|
 | Integration over the entire Slave-Element: no mapping sxi->eta                       |
 | required                                                                             |
 *--------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
void MORTAR::IntegratorCalc<distypeS, distypeM>::IntegrateEleBased2D(MORTAR::Element& sele,
    std::vector<MORTAR::Element*> meles, bool* boundary_ele, const Epetra_Comm& comm)
{
  // check for problem dimension
  if (ndim_ != 2) dserror("2D integration method called for non-2D problem");

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ndof = dynamic_cast<MORTAR::Node*>(sele.Nodes()[0])->NumDof();
  int nodemaster = meles[0]->NumNode();

  // create empty vectors for shape fct. evaluation
  static CORE::LINALG::Matrix<ns_, 1> sval;
  static CORE::LINALG::Matrix<nm_, 1> mval;
  static CORE::LINALG::Matrix<ns_, 1> lmval;

  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if (!mynodes) dserror("IntegrateAndDerivSegment: Null pointer!");

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k = 0; k < nrow; ++k)
  {
    MORTAR::Node* mymrtrnode = dynamic_cast<MORTAR::Node*>(mynodes[k]);
    if (!mymrtrnode) dserror("IntegrateDerivSegment2D: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  // decide whether linear LM are used for quadratic FE here
  bool linlm = false;
  bool dualquad = false;
  if (lmquadtype_ == INPAR::MORTAR::lagmult_lin && sele.Shape() == CORE::FE::CellType::line3)
  {
    bound = false;  // crosspoints and linear LM NOT at the same time!!!!
    linlm = true;
    if (shapefcn_ == INPAR::MORTAR::shape_dual) dualquad = true;
  }

  double sxia = -1;
  double sxib = 1;
  double sxi[2] = {0.0, 0.0};

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < nGP(); ++gp)  // loop to the end //Number of the GP
  {
    bool is_on_mele = false;
    bool kink_projection = false;
    // coordinates and weight of the GP
    double eta[2] = {Coordinate(gp, 0), 0.0};
    double wgt = Weight(gp);
    sxi[0] = 0.5 * (1 - eta[0]) * sxia + 0.5 * (1 + eta[0]) * sxib;

    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5 * sxia + 0.5 * sxib;

    // evaluate Lagrange multiplier shape functions (on slave element)
    if (lmquadtype_ == INPAR::MORTAR::lagmult_const)
      UTILS::EvaluateShape_LM_Const(shapefcn_, sxi, lmval, sele, nrow);
    else if (linlm)
      UTILS::EvaluateShape_LM_Lin(shapefcn_, sxi, lmval, sele, nrow);
    else
      UTILS::EvaluateShape_LM(shapefcn_, sxi, lmval, sele, nrow);

    // evaluate trace space shape functions (on both elements)
    UTILS::EvaluateShape_Displ(sxi, sval, sele, dualquad);

    //********************************************************************
    //  loop over all involved masterelements
    //********************************************************************
    for (int nummaster = 0; nummaster < (int)meles.size(); ++nummaster)
    {
      // get Master element nodes themselves
      DRT::Node** mnodes = meles[nummaster]->Nodes();
      if (!mnodes) dserror("EleBased_Integration: Null pointer!");

      // project Gauss point onto master element
      double mxi[2] = {0.0, 0.0};
      sxi[0] = eta[0];

      MORTAR::MortarProjector::Impl(sele, *meles[nummaster])
          ->ProjectGaussPoint2D(sele, eta, *meles[nummaster], mxi);

      // evaluate trace space shape functions (on both elements)
      UTILS::EvaluateShape_Displ(mxi, mval, *meles[nummaster], false);

      // check GP projection
      if ((mxi[0] >= -1) && (mxi[0] <= 1) && (kink_projection == false))
      {
        kink_projection = true;
        is_on_mele = true;

        // compute segment D/M matrix ****************************************
        double jac = dsxideta * dxdsxi;
        GP_DM(sele, *meles[nummaster], lmval, sval, mval, jac, wgt, nrow, nodemaster, ndof, bound,
            comm);

      }  // end - if Pojection on MasterElement
    }    // loop-end over all involved master elements

    // strong discontinuity --> Boundary Segmentation
    if (is_on_mele == false)
    {
      *boundary_ele = true;
      std::cout << "WARNING: strong disc. occurs !!!" << std::endl;
    }

  }  // loop-end over all gp

  return;
}


/*----------------------------------------------------------------------*
 |  Integrate and linearize a 1D slave / master overlap (2D)  popp 02/09|
 |  This method integrates the overlap D/M matrix and weighted gap g~   |
 |  and stores it in mseg and gseg respectively. Moreover, derivatives  |
 |  LinD/M and Ling are built and stored directly into adjacent nodes.  |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
void MORTAR::IntegratorCalc<distypeS, distypeM>::IntegrateSegment2D(MORTAR::Element& sele,
    double& sxia, double& sxib, MORTAR::Element& mele, double& mxia, double& mxib,
    const Epetra_Comm& comm)
{
  // get LMtype
  INPAR::MORTAR::LagMultQuad lmtype = lmquadtype_;

  // explicitly defined shape function type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror("IntegrateDerivSegment2D called without specific shape function defined!");

  // check for problem dimension
  if (ndim_ != 2) dserror("2D integration method called for non-2D problem");

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("IntegrateAndDerivSegment called on a wrong type of MORTAR::Element pair!");
  if ((sxia < -1.0) || (sxib > 1.0))
    dserror("IntegrateAndDerivSegment called with infeasible slave limits!");
  if ((mxia < -1.0) || (mxib > 1.0))
    dserror("IntegrateAndDerivSegment called with infeasible master limits!");

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = dynamic_cast<MORTAR::Node*>(sele.Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  static CORE::LINALG::Matrix<ns_, 1> sval;
  static CORE::LINALG::Matrix<nm_, 1> mval;
  static CORE::LINALG::Matrix<ns_, 1> lmval;

  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if (!mynodes) dserror("IntegrateAndDerivSegment: Null pointer!");

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  //---------------------------------
  // do trafo for bound elements
  CORE::LINALG::SerialDenseMatrix trafo(nrow, nrow, true);
  bool bound = false;

  // get number of bound nodes
  std::vector<int> ids;
  for (int i = 0; i < nrow; ++i)
  {
    Node* mymrtrnode = dynamic_cast<Node*>(sele.Nodes()[i]);
    if (mymrtrnode->IsOnBoundorCE())
    {
      // get local bound id
      ids.push_back(i);
      bound = true;
    }
  }

  int numbound = (int)ids.size();

  // if all bound: error
  if ((nrow - numbound) < 1e-12) return;

  const double factor = 1.0 / (nrow - numbound);
  // row loop
  for (int i = 0; i < nrow; ++i)
  {
    Node* mymrtrnode = dynamic_cast<Node*>(sele.Nodes()[i]);
    if (!mymrtrnode->IsOnBoundorCE())
    {
      trafo(i, i) = 1.0;
      for (int j = 0; j < (int)ids.size(); ++j) trafo(i, ids[j]) = factor;
    }
  }

  // decide whether linear LM are used for quadratic FE here
  bool linlm = false;
  bool dualquad = false;
  if (lmtype == INPAR::MORTAR::lagmult_lin && sele.Shape() == CORE::FE::CellType::line3)
  {
    bound = false;  // crosspoints and linear LM NOT at the same time!!!!
    linlm = true;
    if (shapefcn_ == INPAR::MORTAR::shape_dual) dualquad = true;
  }

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < nGP(); ++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp, 0), 0.0};
    double wgt = Weight(gp);

    // coordinate transformation sxi->eta (slave MORTAR::Element->Overlap)
    double sxi[2] = {0.0, 0.0};
    sxi[0] = 0.5 * (1 - eta[0]) * sxia + 0.5 * (1 + eta[0]) * sxib;

    // project Gauss point onto master element
    double mxi[2] = {0.0, 0.0};
    MORTAR::MortarProjector::Impl(sele, mele)->ProjectGaussPoint2D(sele, sxi, mele, mxi);

    // check GP projection
    if ((mxi[0] < mxia) || (mxi[0] > mxib))
    {
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Gauss point: " << sxi[0] << " " << sxi[1] << std::endl;
      std::cout << "Projection: " << mxi[0] << " " << mxi[1] << std::endl;
      dserror("IntegrateAndDerivSegment: Gauss point projection failed! mxi=%d", mxi[0]);
    }

    // evaluate Lagrange multiplier shape functions (on slave element)
    if (lmtype == INPAR::MORTAR::lagmult_const)
      UTILS::EvaluateShape_LM_Const(shapefcn_, sxi, lmval, sele, nrow);
    else if (linlm)
      UTILS::EvaluateShape_LM_Lin(shapefcn_, sxi, lmval, sele, nrow);
    else
      UTILS::EvaluateShape_LM(shapefcn_, sxi, lmval, sele, nrow);

    // transform shape functions for bound case
    if (bound)
    {
      CORE::LINALG::SerialDenseVector tempval(nrow, true);
      for (int i = 0; i < nrow; ++i)
        for (int j = 0; j < nrow; ++j) tempval(i) += trafo(i, j) * lmval(j);

      for (int i = 0; i < nrow; ++i) lmval(i) = tempval(i);
    }

    // evaluate trace space shape functions (on both elements)
    UTILS::EvaluateShape_Displ(sxi, sval, sele, dualquad);
    UTILS::EvaluateShape_Displ(mxi, mval, mele, false);

    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5 * sxia + 0.5 * sxib;

    // compute segment D/M matrix ****************************************
    double jac = dsxideta * dxdsxi;
    GP_DM(sele, mele, lmval, sval, mval, jac, wgt, nrow, ncol, ndof, bound, comm);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Compute entries for D and M matrix at GP                 farah 12/13|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
void inline MORTAR::IntegratorCalc<distypeS, distypeM>::GP_DM(MORTAR::Element& sele,
    MORTAR::Element& mele, CORE::LINALG::Matrix<ns_, 1>& lmval, CORE::LINALG::Matrix<ns_, 1>& sval,
    CORE::LINALG::Matrix<nm_, 1>& mval, double& jac, double& wgt, int& nrow, int& ncol, int& ndof,
    bool& bound, const Epetra_Comm& comm)
{
  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("IntegrateAndDerivSegment: Null pointer!");
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes) dserror("IntegrateAndDerivSegment: Null pointer!");

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
  if ((shapefcn_ == INPAR::MORTAR::shape_dual ||
          shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin) &&
      lmquadtype_ != INPAR::MORTAR::lagmult_lin)
  {
    for (int j = 0; j < nrow; ++j)
    {
      MORTAR::Node* cnode = dynamic_cast<MORTAR::Node*>(snodes[j]);

      if (cnode->Owner() != comm.MyPID()) continue;
      if (cnode->IsOnBoundorCE()) continue;

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        MORTAR::Node* mnode = dynamic_cast<MORTAR::Node*>(mnodes[k]);

        // multiply the two shape functions
        double prod = lmval(j) * mval(k) * jac * wgt;
        if (abs(prod) > MORTARINTTOL)
        {
          cnode->AddMValue(mnode->Id(), prod);
          if (!bound) cnode->AddDValue(cnode->Id(), prod);
        }
      }

      // integrate dseg (boundary modification)
      if (bound)
      {
        bool j_boundnode = cnode->IsOnBoundorCE();

        for (int k = 0; k < nrow; ++k)
        {
          MORTAR::Node* mnode = dynamic_cast<MORTAR::Node*>(snodes[k]);

          bool k_boundnode = mnode->IsOnBoundorCE();

          // do not assemble off-diagonal terms if j,k are both non-boundary nodes
          if (!j_boundnode && !k_boundnode && (j != k)) continue;

          // multiply the two shape functions
          double prod = lmval(j) * sval(k) * jac * wgt;
          if (abs(prod) > MORTARINTTOL)
          {
            // isolate the dseg entries to be filled
            // (both the main diagonal and every other secondary diagonal)
            // and add current Gauss point's contribution to dseg
            if (mnode->IsOnBoundorCE())
              cnode->AddMValue(mnode->Id(), -prod);
            else
              cnode->AddDValue(mnode->Id(), prod);
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
      MORTAR::Node* cnode = dynamic_cast<MORTAR::Node*>(snodes[j]);

      if (cnode->Owner() != comm.MyPID()) continue;
      if ((shapefcn_ == INPAR::MORTAR::shape_standard && cnode->IsOnBoundorCE()) ||
          ((shapefcn_ == INPAR::MORTAR::shape_dual ||
               shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin) &&
              cnode->IsOnBound()))
        continue;

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        MORTAR::Node* mnode = dynamic_cast<MORTAR::Node*>(mnodes[k]);

        // multiply the two shape functions
        double prod = lmval(j) * mval(k) * jac * wgt;
        if (abs(prod) > MORTARINTTOL) cnode->AddMValue(mnode->Id(), prod);
      }

      // integrate dseg
      for (int k = 0; k < nrow; ++k)
      {
        MORTAR::Node* snode = dynamic_cast<MORTAR::Node*>(snodes[k]);

        // multiply the two shape functions
        double prod = lmval(j) * sval(k) * jac * wgt;
        if (abs(prod) > MORTARINTTOL)
        {
          if ((shapefcn_ == INPAR::MORTAR::shape_standard && snode->IsOnBoundorCE()) ||
              ((shapefcn_ == INPAR::MORTAR::shape_dual ||
                   shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin) &&
                  snode->IsOnBound()))
            cnode->AddMValue(snode->Id(), -prod);
          else
            cnode->AddDValue(snode->Id(), prod);
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Compute entries for D and M matrix at GP (3D Quad)       farah 12/13|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
void inline MORTAR::IntegratorCalc<distypeS, distypeM>::GP_3D_DM_Quad(MORTAR::Element& sele,
    MORTAR::Element& mele, MORTAR::IntElement& sintele, CORE::LINALG::SerialDenseVector& lmval,
    CORE::LINALG::SerialDenseVector& lmintval, CORE::LINALG::Matrix<ns_, 1>& sval,
    CORE::LINALG::Matrix<nm_, 1>& mval, double& jac, double& wgt, int& nrow, int& nintrow,
    int& ncol, int& ndof, bool& bound)
{
  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("Null pointer!");
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes) dserror("Null pointer!");
  DRT::Node** sintnodes = sintele.Nodes();
  if (!sintnodes) dserror("Null pointer for sintnodes!");

  // CASES 1/2: standard LM shape functions and quadratic or linear interpolation
  // CASE 5: dual LM shape functions and linear interpolation
  // (here, we must NOT ignore the small off-diagonal terms for accurate convergence)
  if ((shapefcn_ == INPAR::MORTAR::shape_standard &&
          (lmquadtype_ == INPAR::MORTAR::lagmult_quad ||
              lmquadtype_ == INPAR::MORTAR::lagmult_lin)) ||
      ((shapefcn_ == INPAR::MORTAR::shape_dual ||
           shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin) &&
          (lmquadtype_ == INPAR::MORTAR::lagmult_lin ||
              lmquadtype_ == INPAR::MORTAR::lagmult_const)))
  {
    // compute all mseg and dseg matrix entries
    // loop over Lagrange multiplier dofs j
    for (int j = 0; j < nrow; ++j)
    {
      MORTAR::Node* cnode = dynamic_cast<MORTAR::Node*>(snodes[j]);

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        MORTAR::Node* mnode = dynamic_cast<MORTAR::Node*>(mnodes[k]);

        // multiply the two shape functions
        double prod = lmval[j] * mval(k) * jac * wgt;
        if (abs(prod) > MORTARINTTOL) cnode->AddMValue(mnode->Id(), prod);
      }

      // integrate dseg
      for (int k = 0; k < nrow; ++k)
      {
        MORTAR::Node* snode = dynamic_cast<MORTAR::Node*>(snodes[k]);

        // multiply the two shape functions
        double prod = lmval[j] * sval(k) * jac * wgt;
        if (abs(prod) > MORTARINTTOL)
        {
          if ((shapefcn_ == INPAR::MORTAR::shape_standard && snode->IsOnBoundorCE()) ||
              ((shapefcn_ == INPAR::MORTAR::shape_dual ||
                   shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin) &&
                  snode->IsOnBound()))
            cnode->AddMValue(snode->Id(), -prod);
          else
            cnode->AddDValue(snode->Id(), prod);
        }
      }
    }
  }
  // CASE 3: standard LM shape functions and piecewise linear interpolation
  else if (shapefcn_ == INPAR::MORTAR::shape_standard &&
           lmquadtype_ == INPAR::MORTAR::lagmult_pwlin)
  {
    // compute all mseg and dseg matrix entries
    // loop over Lagrange multiplier dofs j
    for (int j = 0; j < nintrow; ++j)
    {
      MORTAR::Node* cnode = dynamic_cast<MORTAR::Node*>(sintnodes[j]);

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        MORTAR::Node* mnode = dynamic_cast<MORTAR::Node*>(mnodes[k]);

        // multiply the two shape functions
        double prod = lmintval[j] * mval(k) * jac * wgt;
        if (abs(prod) > MORTARINTTOL) cnode->AddMValue(mnode->Id(), prod);
      }

      // integrate dseg
      for (int k = 0; k < nrow; ++k)
      {
        MORTAR::Node* snode = dynamic_cast<MORTAR::Node*>(snodes[k]);

        // multiply the two shape functions
        double prod = lmintval[j] * sval(k) * jac * wgt;
        if (abs(prod) > MORTARINTTOL)
        {
          if (snode->IsOnBound())
            cnode->AddMValue(snode->Id(), -prod);
          else
            cnode->AddDValue(snode->Id(), prod);
        }
      }
    }
  }
  // CASE 4: dual LM shape functions and quadratic interpolation
  else if ((shapefcn_ == INPAR::MORTAR::shape_dual ||
               shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin) &&
           lmquadtype_ == INPAR::MORTAR::lagmult_quad)
  {
    // compute all mseg and dseg matrix entries
    // loop over Lagrange multiplier dofs j
    for (int j = 0; j < nrow; ++j)
    {
      MORTAR::Node* cnode = dynamic_cast<MORTAR::Node*>(snodes[j]);

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        MORTAR::Node* mnode = dynamic_cast<MORTAR::Node*>(mnodes[k]);

        // multiply the two shape functions
        double prod = lmval[j] * mval(k) * jac * wgt;
        if (abs(prod) > MORTARINTTOL)
        {
          cnode->AddDValue(cnode->Id(), prod);
          cnode->AddMValue(mnode->Id(), prod);
        }
      }
    }
  }
  // INVALID CASES
  else
    dserror("Invalid integration case for 3D quadratic mortar!");

  return;
}


/*----------------------------------------------------------------------*
 |  Integrate Mmod on slave / master overlap (2D)             popp 01/08|
 |  This method integrates the modification to the Mortar matrix M      |
 |  for curved interface (Paper by Puso/Wohlmuth) from given local      |
 |  coordinates sxia to sxib. The corresponding master side local       |
 |  element coordinates given by mxia and mxib                          |
 |  Output is an LINALG::SerialDenseMatrix holding the int. values       |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>
MORTAR::IntegratorCalc<distypeS, distypeM>::IntegrateMmod2D(MORTAR::Element& sele, double& sxia,
    double& sxib, MORTAR::Element& mele, double& mxia, double& mxib)
{
  //**********************************************************************
  dserror("IntegrateMmod2D method is outdated!");
  //**********************************************************************

  // check for problem dimension
  if (ndim_ != 2) dserror("2D integration method called for non-2D problem");

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("IntegrateMmod2D called on a wrong type of MORTAR::Element pair!");
  if ((sxia < -1.0) || (sxib > 1.0))
    dserror("IntegrateMmod2D called with infeasible slave limits!");
  if ((mxia < -1.0) || (mxib > 1.0))
    dserror("IntegrateMmod2D called with infeasible master limits!");

  // create empty mmodseg object and wrap it with Teuchos::RCP
  int nrow = sele.NumNode();
  int nrowdof = ndim_;
  int ncol = mele.NumNode();
  int ncoldof = ndim_;
  Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> mmodseg =
      Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix(nrow * nrowdof, ncol * ncoldof));

  // create empty vectors for shape fct. evaluation
  CORE::LINALG::SerialDenseVector sval(nrow);
  CORE::LINALG::SerialDenseMatrix sderiv(nrow, 1);
  CORE::LINALG::SerialDenseVector mval(ncol);
  CORE::LINALG::SerialDenseMatrix mderiv(ncol, 1);
  CORE::LINALG::SerialDenseVector lmval(nrow);
  CORE::LINALG::SerialDenseMatrix lmderiv(nrow, 1);

  // loop over all Gauss points for integration
  for (int gp = 0; gp < nGP(); ++gp)
  {
    double eta[2] = {Coordinate(gp, 0), 0.0};
    double wgt = Weight(gp);

    double sxi[2] = {0.0, 0.0};
    double mxi[2] = {0.0, 0.0};

    // coordinate transformation sxi->eta (slave MORTAR::Element->Overlap)
    sxi[0] = 0.5 * (1 - eta[0]) * sxia + 0.5 * (1 + eta[0]) * sxib;

    // project Gauss point onto master element
    MORTAR::MortarProjector::Impl(sele, mele)->ProjectGaussPoint2D(sele, sxi, mele, mxi);

    // check GP projection
    if ((mxi[0] < mxia) || (mxi[0] > mxib))
    {
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Gauss point: " << sxi[0] << " " << sxi[1] << std::endl;
      std::cout << "Projection: " << mxi[0] << " " << mxi[1] << std::endl;
      dserror("IntegrateMmod2D: Gauss point projection failed!");
    }

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi, sval, sderiv, nrow);
    mele.EvaluateShape(mxi, mval, mderiv, ncol);

    // build the delta function of slave side shape functions
    double deltasval = sval[0] - sval[1];

    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5 * sxia + 0.5 * sxib;

    /* loop over all mmodseg matrix entries
     nrow represents the slave Lagrange multipliers !!!
     ncol represents the master dofs !!!
     (this DOES matter here for mmodseg, as it might
     sometimes be rectangular, not quadratic!)              */
    for (int j = 0; j < nrow * nrowdof; ++j)
    {
      for (int k = 0; k < ncol * ncoldof; ++k)
      {
        // multiply the two shape functions
        int mindex = (int)(k / ncoldof);
        double prod = 0.5 * deltasval * mval[mindex];
        // add current Gauss point's contribution to mmodseg
        (*mmodseg)(j, k) += prod * dxdsxi * dsxideta * wgt;
      }
    }
  }  // for (int gp=0;gp<nGP();++gp)

  // prepare computation of purely geometric part of Mmod entries
  Node* snode0 = dynamic_cast<Node*>(sele.Nodes()[0]);
  Node* snode1 = dynamic_cast<Node*>(sele.Nodes()[1]);

  // normals
  double n[2][2];
  n[0][0] = snode0->MoData().n()[0];
  n[0][1] = snode0->MoData().n()[1];
  n[1][0] = snode1->MoData().n()[0];
  n[1][1] = snode1->MoData().n()[1];

  // scalar product n1 * n2
  double n1n2 = 0.0;
  for (int i = 0; i < 2; ++i) n1n2 += n[0][i] * n[1][i];

  // vector product n1 x n2
  double n1xn2 = n[0][0] * n[1][1] - n[0][1] * n[1][0];

  // // multiply geometric part onto Mmod
  for (int i = 0; i < ncol; ++i)
  {
    (*mmodseg)(0, 0 + i * ncoldof) *= (1.0 - n1n2);
    (*mmodseg)(1, 0 + i * ncoldof) *= n1xn2;
    (*mmodseg)(0, 1 + i * ncoldof) *= -n1xn2;
    (*mmodseg)(1, 1 + i * ncoldof) *= (1.0 - n1n2);

    (*mmodseg)(2, 0 + i * ncoldof) *= (n1n2 - 1.0);
    (*mmodseg)(3, 0 + i * ncoldof) *= n1xn2;
    (*mmodseg)(2, 1 + i * ncoldof) *= -n1xn2;
    (*mmodseg)(3, 1 + i * ncoldof) *= (n1n2 - 1.0);
  }

  return mmodseg;
}


/*----------------------------------------------------------------------*
 |  Integrate and linearize without segmentation             farah 01/13|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
void MORTAR::IntegratorCalc<distypeS, distypeM>::IntegrateEleBased3D(MORTAR::Element& sele,
    std::vector<MORTAR::Element*> meles, bool* boundary_ele, const Epetra_Comm& comm)
{
  // explicitly defined shape function type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror("IntegrateDerivCell3DAuxPlane called without specific shape function defined!");

  // check for problem dimension
  if (ndim_ != 3) dserror("3D integration method called for non-3D problem");

  // discretization type of master element
  CORE::FE::CellType dt = meles[0]->Shape();

  // check input data
  for (int test = 0; test < (int)meles.size(); ++test)
  {
    if ((!sele.IsSlave()) || (meles[test]->IsSlave()))
      dserror("IntegrateDerivCell3D called on a wrong type of MORTAR::Element pair!");
  }

  int msize = meles.size();
  int nrow = sele.NumNode();
  int nmnode = meles[0]->NumNode();
  int ndof = dynamic_cast<MORTAR::Node*>(sele.Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  static CORE::LINALG::Matrix<ns_, 1> sval;
  static CORE::LINALG::Matrix<nm_, 1> mval;
  static CORE::LINALG::Matrix<ns_, 1> lmval;

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < nGP(); ++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp, 0), Coordinate(gp, 1)};
    double wgt = Weight(gp);

    // note that the third component of sxi is necessary!
    // (although it will always be 0.0 of course)
    // double tempsxi[3] = {0.0, 0.0, 0.0};
    double sxi[2] = {0.0, 0.0};
    double mxi[2] = {0.0, 0.0};
    double projalpha = 0.0;

    sxi[0] = eta[0];
    sxi[1] = eta[1];

    // evaluate the two Jacobians (int. cell and slave element)
    // double jaccell = cell->Jacobian(eta);
    double jacslave = sele.Jacobian(sxi);

    // evaluate Lagrange mutliplier shape functions (on slave element)
    if (lmquadtype_ == INPAR::MORTAR::lagmult_const)
      UTILS::EvaluateShape_LM_Const(shapefcn_, sxi, lmval, sele, nrow);
    else
      UTILS::EvaluateShape_LM(shapefcn_, sxi, lmval, sele, nrow);

    // evaluate trace space shape functions (on both elements)
    UTILS::EvaluateShape_Displ(sxi, sval, sele, false);

    // check for Boundary Segmentation
    bool projactable_gp = false;

    //*******************************************************************
    // loop over meles
    //*******************************************************************
    for (int nummaster = 0; nummaster < msize; ++nummaster)
    {
      // project Gauss point onto master element
      bool is_on_mele = MORTAR::MortarProjector::Impl(sele, *meles[nummaster])
                            ->ProjectGaussPoint3D(sele, sxi, *meles[nummaster], mxi, projalpha);
      if (not is_on_mele) continue;

      // check GP projection
      double tol = 0.00;
      if (dt == CORE::FE::CellType::quad4 || dt == CORE::FE::CellType::quad8 ||
          dt == CORE::FE::CellType::quad9 || dt == CORE::FE::CellType::nurbs8 ||
          dt == CORE::FE::CellType::nurbs9)
      {
        if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol)
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

      if (is_on_mele == true)
      {
        projactable_gp = true;

        // evaluate trace space shape functions (on both elements)
        UTILS::EvaluateShape_Displ(mxi, mval, *meles[nummaster], false);

        // compute cell D/M matrix *******************************************
        bool bound = false;
        GP_DM(sele, *meles[nummaster], lmval, sval, mval, jacslave, wgt, nrow, nmnode, ndof, bound,
            comm);
      }  // is_on_mele==true
    }    // loop over meles

    // strong discontinuity --> Boundary Segmentation
    if (projactable_gp == false) *boundary_ele = true;
  }

  return;
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
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
void MORTAR::IntegratorCalc<distypeS, distypeM>::IntegrateCell3DAuxPlane(MORTAR::Element& sele,
    MORTAR::Element& mele, Teuchos::RCP<MORTAR::IntCell> cell, double* auxn,
    const Epetra_Comm& comm)
{
  // explicitly defined shape function type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror("IntegrateDerivCell3DAuxPlane called without specific shape function defined!");

  // check for problem dimension
  if (ndim_ != 3) dserror("3D integration method called for non-3D problem");

  // discretization type of master element
  CORE::FE::CellType sdt = sele.Shape();
  CORE::FE::CellType mdt = mele.Shape();

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("IntegrateDerivCell3DAuxPlane called on a wrong type of MORTAR::Element pair!");
  if (cell == Teuchos::null)
    dserror("IntegrateDerivCell3DAuxPlane called without integration cell");

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = dynamic_cast<MORTAR::Node*>(sele.Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  static CORE::LINALG::Matrix<ns_, 1> sval;
  static CORE::LINALG::Matrix<nm_, 1> mval;
  static CORE::LINALG::Matrix<ns_, 1> lmval;

  //---------------------------------
  // do trafo for bound elements
  CORE::LINALG::SerialDenseMatrix trafo(nrow, nrow, true);
  bool bound = false;

  // get number of bound nodes
  std::vector<int> ids;
  for (int i = 0; i < nrow; ++i)
  {
    Node* mymrtrnode = dynamic_cast<Node*>(sele.Nodes()[i]);
    if (mymrtrnode->IsOnBoundorCE())
    {
      // get local bound id
      ids.push_back(i);
      bound = true;
    }
  }

  int numbound = (int)ids.size();

  // if all bound: error
  if ((nrow - numbound) < 1e-12) return;

  const double factor = 1.0 / (nrow - numbound);
  // row loop
  for (int i = 0; i < nrow; ++i)
  {
    Node* mymrtrnode = dynamic_cast<Node*>(sele.Nodes()[i]);
    if (!mymrtrnode->IsOnBoundorCE())
    {
      trafo(i, i) = 1.0;
      for (int j = 0; j < (int)ids.size(); ++j) trafo(i, ids[j]) = factor;
    }
  }

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
    if (sdt == CORE::FE::CellType::quad4 || sdt == CORE::FE::CellType::quad8 ||
        sdt == CORE::FE::CellType::quad9)
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
    if (mdt == CORE::FE::CellType::quad4 || mdt == CORE::FE::CellType::quad8 ||
        mdt == CORE::FE::CellType::quad9)
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

    // evaluate Lagrange mutliplier shape functions (on slave element)
    UTILS::EvaluateShape_LM(shapefcn_, sxi, lmval, sele, nrow);

    // transform shape functions for bound case
    if (bound)
    {
      CORE::LINALG::SerialDenseVector tempval(nrow, true);
      for (int i = 0; i < nrow; ++i)
        for (int j = 0; j < nrow; ++j) tempval(i) += trafo(i, j) * lmval(j);

      for (int i = 0; i < nrow; ++i) lmval(i) = tempval(i);
    }


    // evaluate trace space shape functions (on both elements)
    UTILS::EvaluateShape_Displ(sxi, sval, sele, false);
    UTILS::EvaluateShape_Displ(mxi, mval, mele, false);

    // evaluate the integration cell Jacobian
    double jac = cell->Jacobian();

    // compute cell D/M matrix *******************************************
    GP_DM(sele, mele, lmval, sval, mval, jac, wgt, nrow, ncol, ndof, bound, comm);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Integrate and linearize a 2D slave / master cell (3D)     popp 03/09|
 |  This method integrates the cell M matrix (and possibly D matrix)    |
 |  and stores it in mseg and dseg respectively.                        |
 |  This is the QUADRATIC auxiliary plane coupling version!!!           |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
void MORTAR::IntegratorCalc<distypeS, distypeM>::IntegrateCell3DAuxPlaneQuad(MORTAR::Element& sele,
    MORTAR::Element& mele, MORTAR::IntElement& sintele, MORTAR::IntElement& mintele,
    Teuchos::RCP<MORTAR::IntCell> cell, double* auxn)
{
  // get LMtype
  INPAR::MORTAR::LagMultQuad lmtype = lmquadtype_;

  // explicitly defined shape function type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror(
        "ERROR: IntegrateDerivCell3DAuxPlaneQuad called without specific shape function defined!");

  // check for problem dimension
  if (ndim_ != 3) dserror("3D integration method called for non-3D problem");

  if (cell->Shape() != CORE::FE::CellType::tri3) dserror("wrong cell shape!");

  // discretization type of slave and master IntElement
  CORE::FE::CellType sdt = sintele.Shape();
  CORE::FE::CellType mdt = mintele.Shape();

  // discretization type of slave and master Element
  CORE::FE::CellType psdt = sele.Shape();
  CORE::FE::CellType pmdt = mele.Shape();

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror(
        "ERROR: IntegrateDerivCell3DAuxPlaneQuad called on a wrong type of MORTAR::Element pair!");
  if (cell == Teuchos::null)
    dserror("IntegrateDerivCell3DAuxPlaneQuad called without integration cell");

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int nintrow = sintele.NumNode();
  int ndof = dynamic_cast<MORTAR::Node*>(sele.Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  CORE::LINALG::Matrix<ns_, 1> sval;
  CORE::LINALG::Matrix<nm_, 1> mval;
  CORE::LINALG::SerialDenseVector lmval(nrow);
  CORE::LINALG::SerialDenseMatrix lmderiv(nrow, 2, true);
  CORE::LINALG::SerialDenseVector lmintval(nintrow);
  CORE::LINALG::SerialDenseMatrix lmintderiv(nintrow, 2, true);

  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if (!mynodes) dserror("IntegrateDerivCell3DAuxPlaneQuad: Null pointer!");

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k = 0; k < nrow; ++k)
  {
    MORTAR::Node* mymrtrnode = dynamic_cast<MORTAR::Node*>(mynodes[k]);
    if (!mymrtrnode) dserror("IntegrateDerivSegment2D: Null pointer!");
    bound += mymrtrnode->IsOnBoundorCE();
  }

  // decide whether displacement shape fct. modification has to be considered or not
  // this is the case for dual quadratic and linear Lagrange multipliers on quad9/quad8/tri6
  // elements
  bool dualquad3d = false;
  if ((shapefcn_ == INPAR::MORTAR::shape_dual) &&
      (lmtype == INPAR::MORTAR::lagmult_quad || lmtype == INPAR::MORTAR::lagmult_lin) &&
      (sele.Shape() == CORE::FE::CellType::quad9 || sele.Shape() == CORE::FE::CellType::quad8 ||
          sele.Shape() == CORE::FE::CellType::tri6))
  {
    dualquad3d = true;
  }

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
    double tol = 0.01;
    if (sdt == CORE::FE::CellType::quad4 || sdt == CORE::FE::CellType::quad8 ||
        sdt == CORE::FE::CellType::quad9)
    {
      if (sxi[0] < -1.0 - tol || sxi[1] < -1.0 - tol || sxi[0] > 1.0 + tol || sxi[1] > 1.0 + tol)
      {
        std::cout
            << "\n***Warning: IntegrateDerivCell3DAuxPlane: Slave Gauss point projection outside!";
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
            << "\n***Warning: IntegrateDerivCell3DAuxPlane: Slave Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP (IntElement) projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }

    // check GP projection (MASTER)
    if (mdt == CORE::FE::CellType::quad4 || mdt == CORE::FE::CellType::quad8 ||
        mdt == CORE::FE::CellType::quad9)
    {
      if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol)
      {
        std::cout
            << "\n***Warning: IntegrateDerivCell3DAuxPlane: Master Gauss point projection outside!";
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
        std::cout
            << "\n***Warning: IntegrateDerivCell3DAuxPlane: Master Gauss point projection outside!";
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
    if (psdt == CORE::FE::CellType::quad4 || psdt == CORE::FE::CellType::quad8 ||
        psdt == CORE::FE::CellType::quad9 || psdt == CORE::FE::CellType::nurbs9)
    {
      if (psxi[0] < -1.0 - tol || psxi[1] < -1.0 - tol || psxi[0] > 1.0 + tol ||
          psxi[1] > 1.0 + tol)
      {
        std::cout
            << "\n***Warning: IntegrateDerivCell3DAuxPlane: Slave Gauss point projection outside!";
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
            << "\n***Warning: IntegrateDerivCell3DAuxPlane: Slave Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << psxi[0] << " " << psxi[1] << std::endl;
      }
    }

    // check GP projection (MASTER)
    if (pmdt == CORE::FE::CellType::quad4 || pmdt == CORE::FE::CellType::quad8 ||
        pmdt == CORE::FE::CellType::quad9 || pmdt == CORE::FE::CellType::nurbs9)
    {
      if (pmxi[0] < -1.0 - tol || pmxi[1] < -1.0 - tol || pmxi[0] > 1.0 + tol ||
          pmxi[1] > 1.0 + tol)
      {
        std::cout
            << "\n***Warning: IntegrateDerivCell3DAuxPlane: Master Gauss point projection outside!";
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
        std::cout
            << "\n***Warning: IntegrateDerivCell3DAuxPlane: Master Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << pmxi[0] << " " << pmxi[1] << std::endl;
      }
    }

    // evaluate Lagrange multiplier shape functions (on slave element)
    //    if (bound)
    //    {
    //      sele.EvaluateShapeLagMultLin(shapefcn_, psxi, lmval, lmderiv, nrow);
    //    }
    //    else
    //    {
    //      sele.EvaluateShapeLagMult(shapefcn_, psxi, lmval, lmderiv, nrow);
    //      sintele.EvaluateShapeLagMult(shapefcn_, sxi, lmintval, lmintderiv, nintrow);
    //    }

    if (lmtype == INPAR::MORTAR::lagmult_const)
      UTILS::EvaluateShape_LM_Const(shapefcn_, psxi, lmval, sele, nrow);
    else if (bound)
    {
      sele.EvaluateShapeLagMult(shapefcn_, psxi, lmval, lmderiv, nrow);
    }
    else
    {
      sele.EvaluateShapeLagMult(shapefcn_, psxi, lmval, lmderiv, nrow);
      sintele.EvaluateShapeLagMult(shapefcn_, sxi, lmintval, lmintderiv, nintrow, false);
    }

    // evaluate trace space shape functions (on both elements)
    UTILS::EvaluateShape_Displ(psxi, sval, sele, dualquad3d);
    UTILS::EvaluateShape_Displ(pmxi, mval, mele, false);

    // evaluate the integration cell Jacobian
    double jac = cell->Jacobian();

    // compute cell D/M matrix *******************************************
    GP_3D_DM_Quad(sele, mele, sintele, lmval, lmintval, sval, mval, jac, wgt, nrow, nintrow, ncol,
        ndof, bound);
  }
  //**********************************************************************

  return;
}


// line2 slave
template class MORTAR::IntegratorCalc<CORE::FE::CellType::line2, CORE::FE::CellType::line2>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::line2, CORE::FE::CellType::line3>;

// line3 slave
template class MORTAR::IntegratorCalc<CORE::FE::CellType::line3, CORE::FE::CellType::line2>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::line3, CORE::FE::CellType::line3>;

// quad4 slave
template class MORTAR::IntegratorCalc<CORE::FE::CellType::quad4, CORE::FE::CellType::quad4>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::quad4, CORE::FE::CellType::quad8>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::quad4, CORE::FE::CellType::quad9>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::quad4, CORE::FE::CellType::tri3>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::quad4, CORE::FE::CellType::tri6>;

// quad8 slave
template class MORTAR::IntegratorCalc<CORE::FE::CellType::quad8, CORE::FE::CellType::quad4>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::quad8, CORE::FE::CellType::quad8>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::quad8, CORE::FE::CellType::quad9>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::quad8, CORE::FE::CellType::tri3>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::quad8, CORE::FE::CellType::tri6>;

// quad9 slave
template class MORTAR::IntegratorCalc<CORE::FE::CellType::quad9, CORE::FE::CellType::quad4>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::quad9, CORE::FE::CellType::quad8>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::quad9, CORE::FE::CellType::quad9>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::quad9, CORE::FE::CellType::tri3>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::quad9, CORE::FE::CellType::tri6>;

// tri3 slave
template class MORTAR::IntegratorCalc<CORE::FE::CellType::tri3, CORE::FE::CellType::quad4>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::tri3, CORE::FE::CellType::quad8>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::tri3, CORE::FE::CellType::quad9>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::tri3, CORE::FE::CellType::tri3>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::tri3, CORE::FE::CellType::tri6>;

// tri6 slave
template class MORTAR::IntegratorCalc<CORE::FE::CellType::tri6, CORE::FE::CellType::quad4>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::tri6, CORE::FE::CellType::quad8>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::tri6, CORE::FE::CellType::quad9>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::tri6, CORE::FE::CellType::tri3>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::tri6, CORE::FE::CellType::tri6>;

//==================================================
//                     NURBS
//==================================================
// nurbs2 slave
template class MORTAR::IntegratorCalc<CORE::FE::CellType::nurbs2, CORE::FE::CellType::nurbs2>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::nurbs2, CORE::FE::CellType::nurbs3>;

// nurbs3 slave
template class MORTAR::IntegratorCalc<CORE::FE::CellType::nurbs3, CORE::FE::CellType::nurbs2>;
template class MORTAR::IntegratorCalc<CORE::FE::CellType::nurbs3, CORE::FE::CellType::nurbs3>;

BACI_NAMESPACE_CLOSE
