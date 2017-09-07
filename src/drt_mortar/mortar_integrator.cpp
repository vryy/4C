/*!----------------------------------------------------------------------
\file mortar_integrator.cpp

\brief A class to perform integrations of Mortar matrices on the overlap
of two MortarElements in 1D and 2D

\level 1

\maintainer Alexander Seitz

*-----------------------------------------------------------------------*/

#include "mortar_integrator.H"
#include "mortar_node.H"
#include "mortar_element.H"
#include "mortar_projector.H"
#include "mortar_coupling3d_classes.H"
#include "mortar_defines.H"
#include "mortar_calc_utils.H"
#include "mortar_shape_utils.H"

#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/drt_element.H"

#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_sparsematrix.H"

/*----------------------------------------------------------------------*
 |  impl...                                                  farah 01/14|
 *----------------------------------------------------------------------*/
MORTAR::MortarIntegrator* MORTAR::MortarIntegrator::Impl(MortarElement& sele,
    MortarElement& mele, Teuchos::ParameterList& params)
{
  switch (sele.Shape())
  {
  // 2D surface elements
  case DRT::Element::quad4:
  {
    switch (mele.Shape())
    {
    case DRT::Element::quad4:
    {
      return MortarIntegratorCalc<DRT::Element::quad4, DRT::Element::quad4>::Instance(
          true, params);
    }
    case DRT::Element::quad8:
    {
      return MortarIntegratorCalc<DRT::Element::quad4, DRT::Element::quad8>::Instance(
          true, params);
    }
    case DRT::Element::quad9:
    {
      return MortarIntegratorCalc<DRT::Element::quad4, DRT::Element::quad9>::Instance(
          true, params);
    }
    case DRT::Element::tri3:
    {
      return MortarIntegratorCalc<DRT::Element::quad4, DRT::Element::tri3>::Instance(
          true, params);
    }
    case DRT::Element::tri6:
    {
      return MortarIntegratorCalc<DRT::Element::quad4, DRT::Element::tri6>::Instance(
          true, params);
    }
    default:
      dserror("ERROR: Element combination not allowed!");
    }
    break;
  }
  case DRT::Element::quad8:
  {
    switch (mele.Shape())
    {
    case DRT::Element::quad4:
    {
      return MortarIntegratorCalc<DRT::Element::quad8, DRT::Element::quad4>::Instance(
          true, params);
    }
    case DRT::Element::quad8:
    {
      return MortarIntegratorCalc<DRT::Element::quad8, DRT::Element::quad8>::Instance(
          true, params);
    }
    case DRT::Element::quad9:
    {
      return MortarIntegratorCalc<DRT::Element::quad8, DRT::Element::quad9>::Instance(
          true, params);
    }
    case DRT::Element::tri3:
    {
      return MortarIntegratorCalc<DRT::Element::quad8, DRT::Element::tri3>::Instance(
          true, params);
    }
    case DRT::Element::tri6:
    {
      return MortarIntegratorCalc<DRT::Element::quad8, DRT::Element::tri6>::Instance(
          true, params);
    }
    default:
      dserror("ERROR: Element combination not allowed!");
    }
    break;
  }
  case DRT::Element::quad9:
  {
    switch (mele.Shape())
    {
    case DRT::Element::quad4:
    {
      return MortarIntegratorCalc<DRT::Element::quad9, DRT::Element::quad4>::Instance(
          true, params);
    }
    case DRT::Element::quad8:
    {
      return MortarIntegratorCalc<DRT::Element::quad9, DRT::Element::quad8>::Instance(
          true, params);
    }
    case DRT::Element::quad9:
    {
      return MortarIntegratorCalc<DRT::Element::quad9, DRT::Element::quad9>::Instance(
          true, params);
    }
    case DRT::Element::tri3:
    {
      return MortarIntegratorCalc<DRT::Element::quad9, DRT::Element::tri3>::Instance(
          true, params);
    }
    case DRT::Element::tri6:
    {
      return MortarIntegratorCalc<DRT::Element::quad9, DRT::Element::tri6>::Instance(
          true, params);
    }
    default:
      dserror("ERROR: Element combination not allowed!");
    }
    break;
  }
  case DRT::Element::tri3:
  {
    switch (mele.Shape())
    {
    case DRT::Element::quad4:
    {
      return MortarIntegratorCalc<DRT::Element::tri3, DRT::Element::quad4>::Instance(
          true, params);
    }
    case DRT::Element::quad8:
    {
      return MortarIntegratorCalc<DRT::Element::tri3, DRT::Element::quad8>::Instance(
          true, params);
    }
    case DRT::Element::quad9:
    {
      return MortarIntegratorCalc<DRT::Element::tri3, DRT::Element::quad9>::Instance(
          true, params);
    }
    case DRT::Element::tri3:
    {
      return MortarIntegratorCalc<DRT::Element::tri3, DRT::Element::tri3>::Instance(
          true, params);
    }
    case DRT::Element::tri6:
    {
      return MortarIntegratorCalc<DRT::Element::tri3, DRT::Element::tri6>::Instance(
          true, params);
    }
    default:
      dserror("ERROR: Element combination not allowed!");
    }
    break;
  }
  case DRT::Element::tri6:
  {
    switch (mele.Shape())
    {
    case DRT::Element::quad4:
    {
      return MortarIntegratorCalc<DRT::Element::tri6, DRT::Element::quad4>::Instance(
          true, params);
    }
    case DRT::Element::quad8:
    {
      return MortarIntegratorCalc<DRT::Element::tri6, DRT::Element::quad8>::Instance(
          true, params);
    }
    case DRT::Element::quad9:
    {
      return MortarIntegratorCalc<DRT::Element::tri6, DRT::Element::quad9>::Instance(
          true, params);
    }
    case DRT::Element::tri3:
    {
      return MortarIntegratorCalc<DRT::Element::tri6, DRT::Element::tri3>::Instance(
          true, params);
    }
    case DRT::Element::tri6:
    {
      return MortarIntegratorCalc<DRT::Element::tri6, DRT::Element::tri6>::Instance(
          true, params);
    }
    default:
      dserror("ERROR: Element combination not allowed!");
    }
    break;
  }
    //1D surface elements
  case DRT::Element::line2:
  {
    switch (mele.Shape())
    {
    case DRT::Element::line2:
    {
      return MortarIntegratorCalc<DRT::Element::line2, DRT::Element::line2>::Instance(
          true, params);
    }
    case DRT::Element::line3:
    {
      return MortarIntegratorCalc<DRT::Element::line2, DRT::Element::line3>::Instance(
          true, params);
    }
    default:
      dserror("ERROR: Element combination not allowed!");
    }
    break;
  }
  case DRT::Element::line3:
  {
    switch (mele.Shape())
    {
    case DRT::Element::line2:
    {
      return MortarIntegratorCalc<DRT::Element::line3, DRT::Element::line2>::Instance(
          true, params);
    }
    case DRT::Element::line3:
    {
      return MortarIntegratorCalc<DRT::Element::line3, DRT::Element::line3>::Instance(
          true, params);
    }
    default:
      dserror("ERROR: Element combination not allowed!");
    }
    break;
  }

    //==================================================
    //                     NURBS
    //==================================================
    //1D surface elements
  case DRT::Element::nurbs2:
  {
    switch (mele.Shape())
    {
    case DRT::Element::nurbs2:
    {
      return MortarIntegratorCalc<DRT::Element::nurbs2, DRT::Element::nurbs2>::Instance(
          true, params);
    }
    case DRT::Element::nurbs3:
    {
      return MortarIntegratorCalc<DRT::Element::nurbs2, DRT::Element::nurbs3>::Instance(
          true, params);
    }
    default:
      dserror("ERROR: Element combination not allowed!");
    }
    break;
  }
  case DRT::Element::nurbs3:
  {
    switch (mele.Shape())
    {
    case DRT::Element::nurbs2:
    {
      return MortarIntegratorCalc<DRT::Element::nurbs3, DRT::Element::nurbs2>::Instance(
          true, params);
    }
    case DRT::Element::nurbs3:
    {
      return MortarIntegratorCalc<DRT::Element::nurbs3, DRT::Element::nurbs3>::Instance(
          true, params);
    }
    default:
      dserror("ERROR: Element combination not allowed!");
    }
    break;
  }
  case DRT::Element::nurbs9:
  {
    switch (mele.Shape())
    {
    case DRT::Element::nurbs9:
    {
      return MortarIntegratorCalc<DRT::Element::nurbs9, DRT::Element::nurbs9>::Instance(
          true, params);
    }
    case DRT::Element::nurbs4:
    {
      return MortarIntegratorCalc<DRT::Element::nurbs9, DRT::Element::nurbs4>::Instance(
          true, params);
    }
    default:
      dserror("ERROR: Element combination not allowed!");
    }
    break;
  }
  default:
    dserror("Error...");
    break;
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
MORTAR::MortarIntegratorCalc<distypeS, distypeM>::MortarIntegratorCalc(
    Teuchos::ParameterList& params) :
    imortar_(params), shapefcn_(DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(params, "LM_SHAPEFCN")),
    lmquadtype_(DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(params,"LM_QUAD")),
    scale_(DRT::INPUT::IntegralValue<int>(imortar_, "LM_NODAL_SCALE"))
{
  InitializeGP();
}


/*----------------------------------------------------------------------*
 |  Instance (public)                                        farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
MORTAR::MortarIntegratorCalc<distypeS, distypeM> * MORTAR::MortarIntegratorCalc<distypeS, distypeM>::Instance(
    bool create, Teuchos::ParameterList& params)
{
  static MortarIntegratorCalc<distypeS, distypeM> * instance;
  if (create)
  {
    if (instance == NULL)
      instance = new MortarIntegratorCalc<distypeS, distypeM>(params);
  }
  else
  {
    if (instance != NULL)
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 |  Done (public)                                             farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void MORTAR::MortarIntegratorCalc<distypeS, distypeM>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false, imortar_);
}


/*----------------------------------------------------------------------*
 |  Initialize gauss points                                   popp 06/09|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
void MORTAR::MortarIntegratorCalc<distypeS, distypeM>::InitializeGP()
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
  INPAR::MORTAR::IntType integrationtype = DRT::INPUT::IntegralValue<
      INPAR::MORTAR::IntType>(imortar_, "INTTYPE");

  // if we use segment-based integration, the shape of the cells has to be considered!
  DRT::Element::DiscretizationType intshape;
  if (integrationtype == INPAR::MORTAR::inttype_segments)
  {
    if (ndim_ == 2)
      intshape = DRT::Element::line2;
    else if (ndim_ == 3)
      intshape = DRT::Element::tri3;
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
  case DRT::Element::line2:
  case DRT::Element::line3:
  case DRT::Element::nurbs2:
  case DRT::Element::nurbs3:
  {
    // set default value for segment-based version first
    DRT::UTILS::GaussRule1D mygaussrule = DRT::UTILS::intrule_line_5point;

    // GP switch if element-based version and non-zero value provided by user
    if (integrationtype == INPAR::MORTAR::inttype_elements
        || integrationtype == INPAR::MORTAR::inttype_elements_BS)
    {
      if (numgp > 0)
      {
        switch (numgp)
        {
        case 1:
        {
          dserror(
              "Our experience says that 1 GP per slave element is not enough.");
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
    if(integrationtype==INPAR::MORTAR::inttype_segments)
    {
      if (numgp>0)
      switch(numgp)
      {
      case 7 : mygaussrule=DRT::UTILS::intrule_tri_7point;  break;
      case 12: mygaussrule=DRT::UTILS::intrule_tri_12point; break;
      case 16: mygaussrule=DRT::UTILS::intrule_tri_16point; break;
      case 37: mygaussrule=DRT::UTILS::intrule_tri_37point; break;
      default: dserror("unknown tri gauss rule");           break;
      }
    }
    // GP switch if element-based version and non-zero value provided by user
    else if (integrationtype == INPAR::MORTAR::inttype_elements
        || integrationtype == INPAR::MORTAR::inttype_elements_BS)
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
    else
      dserror("ERROR: unknown integration type!");

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
    DRT::UTILS::GaussRule2D mygaussrule = DRT::UTILS::intrule_quad_25point;

    // GP switch if element-based version and non-zero value provided by user
    if (integrationtype == INPAR::MORTAR::inttype_elements
        || integrationtype == INPAR::MORTAR::inttype_elements_BS)
    {
      if (numgp > 0)
      {
        switch (numgp)
        {
        case 1:
        {
          dserror(
              "Our experience says that 1 GP per slave element is not enough.");
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
    dserror(
        "ERROR: MortarIntegrator: This contact element type is not implemented!");
    break;
  }
  } // switch(eletype)

  return;
}


/*--------------------------------------------------------------------------------------*
 | Integrate without segmentation --> more GP required                       farah 01/13|
 | Integration over the entire Slave-Element: no mapping sxi->eta                       |
 | required                                                                             |
 *--------------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void MORTAR::MortarIntegratorCalc<distypeS, distypeM>::IntegrateEleBased2D(
    MORTAR::MortarElement& sele, std::vector<MORTAR::MortarElement*> meles,
    bool *boundary_ele, const Epetra_Comm& comm)
{
  //check for problem dimension
  if (ndim_ != 2)
    dserror("ERROR: 2D integration method called for non-2D problem");

  bool scaling = false;
  if (scale_)
    scaling=true;

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ndof = dynamic_cast<MORTAR::MortarNode*>(sele.Nodes()[0])->NumDof();
  int nodemaster = meles[0]->NumNode();

  // create empty vectors for shape fct. evaluation
  static LINALG::Matrix<ns_, 1> sval;
  static LINALG::Matrix<nm_,1> mval;
  static LINALG::Matrix<ns_,1> lmval;

  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k=0;k<nrow;++k)
  {
    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  // decide whether linear LM are used for quadratic FE here
  bool linlm = false;
  if (lmquadtype_ == INPAR::MORTAR::lagmult_lin && sele.Shape() == DRT::Element::line3)
  {
    bound = false; // crosspoints and linear LM NOT at the same time!!!!
    linlm = true;
  }

  double sxia =-1;
  double sxib = 1;
  double sxi[2]= { 0.0 , 0.0};

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<nGP();++gp)//loop to the end //Number of the GP
  {
    bool is_on_mele      = false;
    bool kink_projection = false;
    // coordinates and weight of the GP
    double eta[2] = { Coordinate(gp,0), 0.0};
    double wgt = Weight(gp);
    sxi[0] = 0.5*(1-eta[0])*sxia + 0.5*(1+eta[0])*sxib;

    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5*sxia + 0.5*sxib;

    // evaluate Lagrange multiplier shape functions (on slave element)
    if (lmquadtype_ == INPAR::MORTAR::lagmult_const)
      UTILS::EvaluateShape_LM_Const(shapefcn_, sxi, lmval, sele, nrow);
    else if (linlm)
      UTILS::EvaluateShape_LM_Lin(shapefcn_,sxi,lmval,sele,nrow);
    else
      UTILS::EvaluateShape_LM(shapefcn_,sxi,lmval,sele,nrow);

    // evaluate trace space shape functions (on both elements)
    UTILS::EvaluateShape_Displ(sxi, sval,sele, false);

    //********************************************************************
    //  loop over all involved masterelements
    //********************************************************************
    for (int nummaster=0;nummaster<(int)meles.size();++nummaster)
    {
      // get Master element nodes themselves
      DRT::Node** mnodes = meles[nummaster]->Nodes();
      if(!mnodes) dserror("ERROR: EleBased_Integration: Null pointer!");

      // project Gauss point onto master element
      double mxi[2] = { 0.0, 0.0};
      sxi[0]=eta[0];

      MORTAR::MortarProjector::Impl(sele,*meles[nummaster])->ProjectGaussPoint2D(sele,eta,*meles[nummaster],mxi);

      // evaluate trace space shape functions (on both elements)
      UTILS::EvaluateShape_Displ(mxi, mval,*meles[nummaster], false);

      // check GP projection
      if ((mxi[0]>=-1) && (mxi[0]<=1) && (kink_projection==false))
      {
        kink_projection=true;
        is_on_mele=true;

        // compute segment D/M matrix ****************************************
        double jac = dsxideta*dxdsxi;
        GP_DM(sele,*meles[nummaster],lmval,sval,mval,jac, wgt,nrow,nodemaster,ndof,bound,comm);

        // compute nodal scaling factor **************************************
        if(scaling)
        GP_2D_Scaling(sele,sval,dsxideta,wgt);
      } // end - if Pojection on MasterElement
    } //loop-end over all involved master elements

    // strong discontinuity --> Boundary Segmentation
    if(is_on_mele==false)
    {
      *boundary_ele=true;
      std::cout << "WARNING: strong disc. occurs !!!"<< std::endl;
    }

  } //loop-end over all gp

  return;
}


/*----------------------------------------------------------------------*
 |  Integrate and linearize a 1D slave / master overlap (2D)  popp 02/09|
 |  This method integrates the overlap D/M matrix and weighted gap g~   |
 |  and stores it in mseg and gseg respectively. Moreover, derivatives  |
 |  LinD/M and Ling are built and stored directly into adjacent nodes.  |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void MORTAR::MortarIntegratorCalc<distypeS, distypeM>::IntegrateSegment2D(
    MORTAR::MortarElement& sele, double& sxia, double& sxib,
    MORTAR::MortarElement& mele, double& mxia, double& mxib,
    const Epetra_Comm& comm)
{
  // get LMtype
  INPAR::MORTAR::LagMultQuad lmtype = lmquadtype_;

  // explicitly defined shape function type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivSegment2D called without specific shape function defined!");

    //check for problem dimension
  if (ndim_!=2)
    dserror("ERROR: 2D integration method called for non-2D problem");

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateAndDerivSegment called on a wrong type of MortarElement pair!");
  if ((sxia<-1.0) || (sxib>1.0))
    dserror("ERROR: IntegrateAndDerivSegment called with infeasible slave limits!");
  if ((mxia<-1.0) || (mxib>1.0))
    dserror("ERROR: IntegrateAndDerivSegment called with infeasible master limits!");

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = dynamic_cast<MORTAR::MortarNode*>(sele.Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  static LINALG::Matrix<ns_, 1> sval;
  static LINALG::Matrix<nm_,1 > mval;
  static LINALG::Matrix<ns_, 1> lmval;

  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if (!mynodes)
    dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  //---------------------------------
  // do trafo for bound elements
  LINALG::SerialDenseMatrix trafo(nrow, nrow,true);
  bool bound = false;

  // get number of bound nodes
  std::vector<int> ids;
  for(int i = 0; i<nrow;++i)
  {
    MortarNode* mymrtrnode = dynamic_cast<MortarNode*>(sele.Nodes()[i]);
    if(mymrtrnode->IsOnBoundorCE())
    {
      // get local bound id
      ids.push_back(i);
      bound = true;
    }
  }

  int numbound = (int)ids.size();

  // if all bound: error
  if((nrow-numbound)<1e-12)
    return;

  const double factor = 1.0/(nrow-numbound);
  // row loop
  for(int i = 0; i<nrow;++i)
  {
    MortarNode* mymrtrnode = dynamic_cast<MortarNode*>(sele.Nodes()[i]);
    if (!mymrtrnode->IsOnBoundorCE())
    {
      trafo(i,i)=1.0;
      for(int j = 0; j<(int)ids.size();++j)
        trafo(i,ids[j]) = factor;
    }
  }

  // decide whether linear LM are used for quadratic FE here
  bool linlm = false;
  if (lmtype == INPAR::MORTAR::lagmult_lin
      && sele.Shape() == DRT::Element::line3)
  {
    bound = false; // crosspoints and linear LM NOT at the same time!!!!
    linlm = true;
  }

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < nGP(); ++gp)
  {
    // coordinates and weight
    double eta[2] = { Coordinate(gp, 0), 0.0 };
    double wgt = Weight(gp);

    // coordinate transformation sxi->eta (slave MortarElement->Overlap)
    double sxi[2] = { 0.0, 0.0 };
    sxi[0] = 0.5 * (1 - eta[0]) * sxia + 0.5 * (1 + eta[0]) * sxib;

    // project Gauss point onto master element
    double mxi[2] = { 0.0, 0.0 };
    MORTAR::MortarProjector::Impl(sele, mele)->ProjectGaussPoint2D(sele, sxi, mele, mxi);

    // check GP projection
    if ((mxi[0] < mxia) || (mxi[0] > mxib))
    {
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Gauss point: " << sxi[0] << " " << sxi[1] << std::endl;
      std::cout << "Projection: " << mxi[0] << " " << mxi[1] << std::endl;
      dserror("ERROR: IntegrateAndDerivSegment: Gauss point projection failed! mxi=%d",
          mxi[0]);
    }

    // evaluate Lagrange multiplier shape functions (on slave element)
    if (lmtype == INPAR::MORTAR::lagmult_const)
      UTILS::EvaluateShape_LM_Const(shapefcn_, sxi, lmval, sele, nrow);
    else if (linlm)
      UTILS::EvaluateShape_LM_Lin(shapefcn_, sxi, lmval, sele, nrow);
    else
      UTILS::EvaluateShape_LM(shapefcn_, sxi, lmval, sele, nrow);

    // transform shape functions for bound case
    if(bound)
    {
      LINALG::SerialDenseVector tempval(nrow,true);
      for(int i = 0; i < nrow; ++i)
        for(int j = 0; j < nrow; ++j)
          tempval(i) += trafo(i,j)*lmval(j);

      for(int i = 0; i < nrow; ++i)
        lmval(i) = tempval(i);
    }

    // evaluate trace space shape functions (on both elements)
    UTILS::EvaluateShape_Displ(sxi, sval, sele, false);
    UTILS::EvaluateShape_Displ(mxi, mval, mele, false);

    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5 * sxia + 0.5 * sxib;

    // compute segment D/M matrix ****************************************
    double jac = dsxideta * dxdsxi;
    GP_DM(sele, mele, lmval, sval, mval, jac, wgt, nrow, ncol, ndof, bound,
        comm);

    // compute nodal scaling factor **************************************
    if (scale_)
      GP_2D_Scaling(sele, sval, dsxideta, wgt);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Compute entries for D and M matrix at GP                 farah 12/13|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void inline MORTAR::MortarIntegratorCalc<distypeS, distypeM>::GP_DM(
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    LINALG::Matrix<ns_, 1>& lmval,
    LINALG::Matrix<ns_, 1>& sval,
    LINALG::Matrix<nm_, 1>& mval,
    double& jac,
    double& wgt, int& nrow, int& ncol,
    int& ndof, bool& bound,
    const Epetra_Comm& comm)
{
  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  if (!snodes)
    dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes)
    dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

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
  // standard shape functions
  if (shapefcn_ == INPAR::MORTAR::shape_standard)
  {
    for (int j = 0; j < nrow; ++j)
    {
      MORTAR::MortarNode* cnode = dynamic_cast<MORTAR::MortarNode*>(snodes[j]);

      if (cnode->Owner() != comm.MyPID())
        continue;

      if (cnode->IsOnBoundorCE())
        continue;

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        MORTAR::MortarNode* mnode =
            dynamic_cast<MORTAR::MortarNode*>(mnodes[k]);

        // multiply the two shape functions
        double prod = lmval(j) * mval(k) * jac * wgt;

        if (abs(prod) > MORTARINTTOL)
          cnode->AddMValue(mnode->Id(), prod);
      }

      // integrate dseg
      for (int k = 0; k < nrow; ++k)
      {
        MORTAR::MortarNode* snode =
            dynamic_cast<MORTAR::MortarNode*>(snodes[k]);

        // multiply the two shape functions
        double prod = lmval(j) * sval(k) * jac * wgt;

        if (snode->IsOnBoundorCE())
        {
          if (abs(prod) > MORTARINTTOL)
            cnode->AddMValue(snode->Id(), -prod);
        }
        else
        {
          if (abs(prod) > MORTARINTTOL)
            cnode->AddDValue(snode->Id(), prod);
        }
      }
    }
  }
  // dual shape functions
  else if (shapefcn_ == INPAR::MORTAR::shape_dual
        || shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
  {
    for (int j = 0; j < nrow; ++j)
    {
      MORTAR::MortarNode* cnode = dynamic_cast<MORTAR::MortarNode*>(snodes[j]);

      if (cnode->Owner() != comm.MyPID())
        continue;

      if (cnode->IsOnBoundorCE())
        continue;

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        MORTAR::MortarNode* mnode =
            dynamic_cast<MORTAR::MortarNode*>(mnodes[k]);

        // multiply the two shape functions
        double prod = lmval(j) * mval(k) * jac * wgt;

        if (!bound and abs(prod) > MORTARINTTOL)
          cnode->AddDValue(cnode->Id(), prod);

        if (abs(prod) > MORTARINTTOL)
          cnode->AddMValue(mnode->Id(), prod);
      }

      // integrate dseg (boundary modification)
      if (bound)
      {
        bool j_boundnode = cnode->IsOnBoundorCE();

        for (int k = 0; k < nrow; ++k)
        {
          MORTAR::MortarNode* mnode =
              dynamic_cast<MORTAR::MortarNode*>(snodes[k]);
          bool k_boundnode = mnode->IsOnBoundorCE();

          // do not assemble off-diagonal terms if j,k are both non-boundary nodes
          if (!j_boundnode && !k_boundnode && (j != k))
            continue;

          // multiply the two shape functions
          double prod = lmval(j) * sval(k) * jac * wgt;

          // isolate the dseg entries to be filled
          // (both the main diagonal and every other secondary diagonal)
          // and add current Gauss point's contribution to dseg
          if (mnode->IsOnBoundorCE())
          {
            if (abs(prod) > MORTARINTTOL)
              cnode->AddMValue(mnode->Id(), -prod);
          }
          else
          {
            if (abs(prod) > MORTARINTTOL)
              cnode->AddDValue(mnode->Id(), prod);
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
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void inline MORTAR::MortarIntegratorCalc<distypeS, distypeM>::GP_3D_DM_Quad(
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    MORTAR::IntElement& sintele, LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseVector& lmintval, LINALG::Matrix<ns_, 1>& sval,
    LINALG::Matrix<nm_, 1>& mval,
    double& jac,
    double& wgt, int& nrow, int& nintrow, int& ncol,
    int& ndof, bool& bound)
{
  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  if (!snodes)
    dserror("ERROR: Null pointer!");
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes)
    dserror("ERROR: Null pointer!");
  DRT::Node** sintnodes = sintele.Nodes();
  if (!sintnodes)
    dserror("ERROR: Null pointer for sintnodes!");

  // CASE 1/2: Standard LM shape functions and quadratic or linear interpolation
  if (shapefcn_ == INPAR::MORTAR::shape_standard
      && (lmquadtype_ == INPAR::MORTAR::lagmult_quad
          || lmquadtype_ == INPAR::MORTAR::lagmult_lin))
  {
    // compute all mseg and dseg matrix entries
    // loop over Lagrange multiplier dofs j
    for (int j = 0; j < nrow; ++j)
    {
      MORTAR::MortarNode* cnode = dynamic_cast<MORTAR::MortarNode*>(snodes[j]);
      if(cnode->IsOnBoundorCE())
        continue;

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        MORTAR::MortarNode* mnode =
            dynamic_cast<MORTAR::MortarNode*>(mnodes[k]);

        // multiply the two shape functions
        double prod = lmval[j] * mval(k) * jac * wgt;

        if (abs(prod) > MORTARINTTOL)
          cnode->AddMValue(mnode->Id(), prod);
      }

      // integrate dseg
      for (int k = 0; k < nrow; ++k)
      {
        MORTAR::MortarNode* snode =
            dynamic_cast<MORTAR::MortarNode*>(snodes[k]);

        // multiply the two shape functions
        double prod = lmval[j] * sval(k) * jac * wgt;

        if (snode->IsOnBoundorCE())
        {
          if (abs(prod) > MORTARINTTOL)
            cnode->AddMValue(snode->Id(), -prod);
        }
        else
        {
          if (abs(prod) > MORTARINTTOL)
            cnode->AddDValue(snode->Id(), prod);
        }
      }
    }
  }

  // CASE 3: Standard LM shape functions and piecewise linear interpolation
  else if (shapefcn_ == INPAR::MORTAR::shape_standard
      && lmquadtype_ == INPAR::MORTAR::lagmult_pwlin)
  {
    // compute all mseg and dseg matrix entries
    // loop over Lagrange multiplier dofs j
    for (int j = 0; j < nintrow; ++j)
    {
      MORTAR::MortarNode* cnode = dynamic_cast<MORTAR::MortarNode*>(sintnodes[j]);

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        MORTAR::MortarNode* mnode =
            dynamic_cast<MORTAR::MortarNode*>(mnodes[k]);

        // multiply the two shape functions
        double prod = lmintval[j] * mval(k) * jac * wgt;

        if (abs(prod) > MORTARINTTOL)
          cnode->AddMValue(mnode->Id(), prod);
      }

      // integrate dseg
      for (int k = 0; k < nrow; ++k)
      {
        MORTAR::MortarNode* snode =
            dynamic_cast<MORTAR::MortarNode*>(snodes[k]);

        // multiply the two shape functions
        double prod = lmintval[j] * sval(k) * jac * wgt;

        if (snode->IsOnBound())
        {
          if (abs(prod) > MORTARINTTOL)
            cnode->AddMValue(snode->Id(),-prod);
        }
        else
        {
          if (abs(prod) > MORTARINTTOL)
            cnode->AddDValue(snode->Id(), prod);
        }
      }
    }
  }

  // CASE 4: Dual LM shape functions and quadratic interpolation
  else if (shapefcn_ == INPAR::MORTAR::shape_dual
      && lmquadtype_ == INPAR::MORTAR::lagmult_quad)
  {
    // compute all mseg and dseg matrix entries
    // loop over Lagrange multiplier dofs j
    for (int j = 0; j < nrow; ++j)
    {
      MORTAR::MortarNode* cnode = dynamic_cast<MORTAR::MortarNode*>(snodes[j]);

      // integrate mseg
      for (int k = 0; k < ncol; ++k)
      {
        MORTAR::MortarNode* mnode =
            dynamic_cast<MORTAR::MortarNode*>(mnodes[k]);

        // multiply the two shape functions
        double prod = lmval[j] * mval(k) * jac * wgt;

        if (abs(prod) > MORTARINTTOL)
          cnode->AddDValue(cnode->Id(), prod);

        if (abs(prod) > MORTARINTTOL)
          cnode->AddMValue(mnode->Id(), prod);
      }
    }
  }

  // CASE 5: Dual LM shape functions and linear interpolation
  // (here, we must NOT ignore the small off-diagonal terms for accurate convergence)
  else if (shapefcn_ == INPAR::MORTAR::shape_dual
      && (   lmquadtype_ == INPAR::MORTAR::lagmult_lin
          || lmquadtype_ == INPAR::MORTAR::lagmult_const ) )
  {
    // compute all mseg and dseg matrix entries
    // loop over Lagrange multiplier dofs j
    for (int j = 0; j < nrow; ++j)
    {
      MORTAR::MortarNode* cnode = dynamic_cast<MORTAR::MortarNode*>(snodes[j]);

        // integrate mseg
        for (int k = 0; k < ncol; ++k)
        {
          MORTAR::MortarNode* mnode =
              dynamic_cast<MORTAR::MortarNode*>(mnodes[k]);

          // multiply the two shape functions
          double prod = lmval[j] * mval(k) * jac * wgt;

          if (abs(prod) > MORTARINTTOL)
            cnode->AddMValue(mnode->Id(), prod);
        }

      // integrate dseg
      for (int k = 0; k < nrow; ++k)
      {
        MORTAR::MortarNode* snode =
            dynamic_cast<MORTAR::MortarNode*>(snodes[k]);

        // multiply the two shape functions
        double prod = lmval[j] * sval(k) * jac * wgt;

        if (snode->IsOnBound())
        {
          if (abs(prod) > MORTARINTTOL)
            cnode->AddMValue(snode->Id(), -prod);
        }
        else
        {
          if (abs(prod) > MORTARINTTOL)
            cnode->AddDValue(snode->Id(), prod);
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
 |  Compute entries for scaling at GP                        farah 12/13|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void MORTAR::MortarIntegratorCalc<distypeS, distypeM>::GP_2D_Scaling(
    MORTAR::MortarElement& sele, LINALG::Matrix<ns_, 1>& sval,
    double& dsxideta, double& wgt)
{
  DRT::Node** snodes = sele.Nodes();

  for (int j = 0; j < ns_; ++j)
  {
    MORTAR::MortarNode* snode = dynamic_cast<MORTAR::MortarNode*>(snodes[j]);

    double prod = wgt * sval(j) * dsxideta / sele.Nodes()[j]->NumElement();
    snode->AddScValue(prod);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Compute entries for scaling at GP                        farah 12/13|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void MORTAR::MortarIntegratorCalc<distypeS, distypeM>::GP_3D_Scaling(
    MORTAR::MortarElement& sele, LINALG::Matrix<ns_, 1>& sval,
    double& jac, double& wgt,
    double* sxi)
{
  double jacsele = sele.Jacobian(sxi);

  DRT::Node** snodes = sele.Nodes();

  for (int j = 0; j < ns_; ++j)
  {
    MORTAR::MortarNode* snode = dynamic_cast<MORTAR::MortarNode*>(snodes[j]);

    double prod = (wgt * sval(j) * jac / jacsele)
        / (sele.Nodes()[j]->NumElement());
    prod *= 6.0;

    snode->AddScValue(prod);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Integrate Mmod on slave / master overlap (2D)             popp 01/08|
 |  This method integrates the modification to the Mortar matrix M      |
 |  for curved interface (Paper by Puso/Wohlmuth) from given local      |
 |  coordinates sxia to sxib. The corresponding master side local       |
 |  element coordinates given by mxia and mxib                          |
 |  Output is an Epetra_SerialDenseMatrix holding the int. values       |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
Teuchos::RCP<Epetra_SerialDenseMatrix> MORTAR::MortarIntegratorCalc<distypeS,distypeM>::IntegrateMmod2D(
    MORTAR::MortarElement& sele, double& sxia,
    double& sxib, MORTAR::MortarElement& mele, double& mxia, double& mxib)
{
  //**********************************************************************
  dserror("ERROR: IntegrateMmod2D method is outdated!");
  //**********************************************************************

  //check for problem dimension
  if (ndim_ != 2)
    dserror("ERROR: 2D integration method called for non-2D problem");

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror(
        "ERROR: IntegrateMmod2D called on a wrong type of MortarElement pair!");
  if ((sxia < -1.0) || (sxib > 1.0))
    dserror("ERROR: IntegrateMmod2D called with infeasible slave limits!");
  if ((mxia < -1.0) || (mxib > 1.0))
    dserror("ERROR: IntegrateMmod2D called with infeasible master limits!");

  // create empty mmodseg object and wrap it with Teuchos::RCP
  int nrow = sele.NumNode();
  int nrowdof = ndim_;
  int ncol = mele.NumNode();
  int ncoldof = ndim_;
  Teuchos::RCP<Epetra_SerialDenseMatrix> mmodseg = Teuchos::rcp(
      new Epetra_SerialDenseMatrix(nrow * nrowdof, ncol * ncoldof));

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow, 1);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol, 1);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow, 1);

  // loop over all Gauss points for integration
  for (int gp = 0; gp < nGP(); ++gp)
  {
    double eta[2] = { Coordinate(gp, 0), 0.0 };
    double wgt = Weight(gp);

    double sxi[2] = { 0.0, 0.0 };
    double mxi[2] = { 0.0, 0.0 };

    // coordinate transformation sxi->eta (slave MortarElement->Overlap)
    sxi[0] = 0.5 * (1 - eta[0]) * sxia + 0.5 * (1 + eta[0]) * sxib;

    // project Gauss point onto master element
    MORTAR::MortarProjector::Impl(sele, mele)->ProjectGaussPoint2D(sele, sxi,
        mele, mxi);

    // check GP projection
    if ((mxi[0] < mxia) || (mxi[0] > mxib))
    {
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Gauss point: " << sxi[0] << " " << sxi[1] << std::endl;
      std::cout << "Projection: " << mxi[0] << " " << mxi[1] << std::endl;
      dserror("ERROR: IntegrateMmod2D: Gauss point projection failed!");
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
        int mindex = (int) (k / ncoldof);
        double prod = 0.5 * deltasval * mval[mindex];
        // add current Gauss point's contribution to mmodseg
        (*mmodseg)(j, k) += prod * dxdsxi * dsxideta * wgt;
      }
    }
  } // for (int gp=0;gp<nGP();++gp)

  // prepare computation of purely geometric part of Mmod entries
  MortarNode* snode0 = dynamic_cast<MortarNode*>(sele.Nodes()[0]);
  MortarNode* snode1 = dynamic_cast<MortarNode*>(sele.Nodes()[1]);

  // normals
  double n[2][2];
  n[0][0] = snode0->MoData().n()[0];
  n[0][1] = snode0->MoData().n()[1];
  n[1][0] = snode1->MoData().n()[0];
  n[1][1] = snode1->MoData().n()[1];

  // scalar product n1 * n2
  double n1n2 = 0.0;
  for (int i = 0; i < 2; ++i)
    n1n2 += n[0][i] * n[1][i];

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
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void MORTAR::MortarIntegratorCalc<distypeS, distypeM>::IntegrateEleBased3D(
    MORTAR::MortarElement& sele, std::vector<MORTAR::MortarElement*> meles,
    bool *boundary_ele, const Epetra_Comm& comm)
{
  // explicitly defined shape function type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without specific shape function defined!");

  //check for problem dimension
  if (ndim_!=3)
    dserror("ERROR: 3D integration method called for non-3D problem");

  bool scaling = false;
  if (scale_)
    scaling=true;

  // discretization type of master element
  DRT::Element::DiscretizationType dt = meles[0]->Shape();

  // check input data
  for (int test=0; test<(int)meles.size();++test)
  {
    if ((!sele.IsSlave()) || (meles[test]->IsSlave()))
    dserror("ERROR: IntegrateDerivCell3D called on a wrong type of MortarElement pair!");
  }

  int msize = meles.size();
  int nrow = sele.NumNode();
  int nmnode = meles[0]->NumNode();
  int ndof = dynamic_cast<MORTAR::MortarNode*>(sele.Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  static LINALG::Matrix<ns_,1> sval;
  static LINALG::Matrix<nm_,1> mval;
  static LINALG::Matrix<ns_,1> lmval;

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<nGP();++gp)
  {
    // coordinates and weight
    double eta[2] = { Coordinate(gp,0), Coordinate(gp,1)};
    double wgt = Weight(gp);

    // note that the third component of sxi is necessary!
    // (although it will always be 0.0 of course)
    //double tempsxi[3] = {0.0, 0.0, 0.0};
    double sxi[2] = { 0.0, 0.0};
    double mxi[2] = { 0.0, 0.0};
    double projalpha = 0.0;

    sxi[0]=eta[0];
    sxi[1]=eta[1];

    // evaluate the two Jacobians (int. cell and slave element)
    //double jaccell = cell->Jacobian(eta);
    double jacslave = sele.Jacobian(sxi);

    // evaluate Lagrange mutliplier shape functions (on slave element)
    if (lmquadtype_ == INPAR::MORTAR::lagmult_const)
      UTILS::EvaluateShape_LM_Const(shapefcn_, sxi, lmval, sele, nrow);
    else
    UTILS::EvaluateShape_LM(shapefcn_,sxi,lmval,sele,nrow);

    // evaluate trace space shape functions (on both elements)
    UTILS::EvaluateShape_Displ(sxi, sval,sele, false);

    //check for Boundary Segmentation
    bool projactable_gp=false;

    //*******************************************************************
    // loop over meles
    //*******************************************************************
    for(int nummaster=0;nummaster<msize;++nummaster)
    {
      // project Gauss point onto master element
      bool is_on_mele = MORTAR::MortarProjector::Impl(sele,*meles[nummaster])->ProjectGaussPoint3D(sele,sxi,*meles[nummaster],mxi,projalpha);
      if(not is_on_mele)
      continue;

      // check GP projection
      double tol = 0.00;
      if (dt==DRT::Element::quad4 || dt==DRT::Element::quad8 || dt==DRT::Element::quad9 ||
          dt==DRT::Element::nurbs8 || dt==DRT::Element::nurbs9)
      {
        if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol)
        {
          is_on_mele=false;
        }
      }
      else
      {
        if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
        {
          is_on_mele=false;
        }
      }

      if (is_on_mele==true)
      {
        projactable_gp=true;

        // evaluate trace space shape functions (on both elements)
        UTILS::EvaluateShape_Displ(mxi, mval,*meles[nummaster], false);

        // compute cell D/M matrix *******************************************
        bool bound =false;
        GP_DM(sele,*meles[nummaster],lmval,sval,mval,jacslave,wgt,nrow,nmnode,ndof,bound,comm);

        // compute nodal scaling factor **************************************
        if (scaling)
          GP_3D_Scaling(sele,sval,jacslave,wgt,sxi);
        // compute nodal scaling factor **************************************

      } //is_on_mele==true
    } //loop over meles

    // strong discontinuity --> Boundary Segmentation
    if(projactable_gp==false)
    *boundary_ele=true;
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
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void MORTAR::MortarIntegratorCalc<distypeS, distypeM>::IntegrateCell3DAuxPlane(
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    Teuchos::RCP<MORTAR::IntCell> cell, double* auxn, const Epetra_Comm& comm)
{
  // explicitly defined shape function type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without specific shape function defined!");

  //check for problem dimension
  if (ndim_!=3)
    dserror("ERROR: 3D integration method called for non-3D problem");

  // discretization type of master element
  DRT::Element::DiscretizationType sdt = sele.Shape();
  DRT::Element::DiscretizationType mdt = mele.Shape();

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called on a wrong type of MortarElement pair!");
  if (cell==Teuchos::null)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without integration cell");

  bool scaling = false;
  if (scale_)
    scaling=true;

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = dynamic_cast<MORTAR::MortarNode*>(sele.Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  static LINALG::Matrix<ns_, 1> sval;
  static LINALG::Matrix<nm_,1> mval;
  static LINALG::Matrix<ns_,1> lmval;

  //---------------------------------
  // do trafo for bound elements
  LINALG::SerialDenseMatrix trafo(nrow, nrow,true);
  bool bound = false;

  // get number of bound nodes
  std::vector<int> ids;
  for(int i = 0; i<nrow;++i)
  {
    MortarNode* mymrtrnode = dynamic_cast<MortarNode*>(sele.Nodes()[i]);
    if(mymrtrnode->IsOnBoundorCE())
    {
      // get local bound id
      ids.push_back(i);
      bound = true;
    }
  }

  int numbound = (int)ids.size();

  // if all bound: error
  if((nrow-numbound)<1e-12)
    return;

  const double factor = 1.0/(nrow-numbound);
  // row loop
  for(int i = 0; i<nrow;++i)
  {
    MortarNode* mymrtrnode = dynamic_cast<MortarNode*>(sele.Nodes()[i]);
    if (!mymrtrnode->IsOnBoundorCE())
    {
      trafo(i,i)=1.0;
      for(int j = 0; j<(int)ids.size();++j)
        trafo(i,ids[j]) = factor;
    }
  }

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<nGP();++gp)
  {
    // coordinates and weight
    double eta[2] = { Coordinate(gp,0), Coordinate(gp,1)};
    double wgt = Weight(gp);

    // get global Gauss point coordinates
    double globgp[3] = { 0.0, 0.0, 0.0};
    cell->LocalToGlobal(eta,globgp,0);

    double sxi[2] = { 0.0, 0.0};
    double mxi[2] = { 0.0, 0.0};

    // project Gauss point onto slave element
    // project Gauss point onto master element
    double sprojalpha = 0.0;
    double mprojalpha = 0.0;
    MORTAR::MortarProjector::Impl(sele)->ProjectGaussPointAuxn3D(globgp,auxn,sele,sxi,sprojalpha);
    MORTAR::MortarProjector::Impl(mele)->ProjectGaussPointAuxn3D(globgp,auxn,mele,mxi,mprojalpha);

    // check GP projection (SLAVE)
    double tol = 0.01;
    if (sdt==DRT::Element::quad4 || sdt==DRT::Element::quad8 || sdt==DRT::Element::quad9)
    {
      if (sxi[0]<-1.0-tol || sxi[1]<-1.0-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }
    else
    {
      if (sxi[0]<-tol || sxi[1]<-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol || sxi[0]+sxi[1]>1.0+2*tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }

    // check GP projection (MASTER)
    if (mdt==DRT::Element::quad4 || mdt==DRT::Element::quad8 || mdt==DRT::Element::quad9)
    {
      if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }
    else
    {
      if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }

    // evaluate Lagrange mutliplier shape functions (on slave element)
    UTILS::EvaluateShape_LM(shapefcn_,sxi,lmval,sele,nrow);

    // transform shape functions for bound case
    if(bound)
    {
      LINALG::SerialDenseVector tempval(nrow,true);
      for(int i = 0; i < nrow; ++i)
        for(int j = 0; j < nrow; ++j)
          tempval(i) += trafo(i,j)*lmval(j);

      for(int i = 0; i < nrow; ++i)
        lmval(i) = tempval(i);
    }


    // evaluate trace space shape functions (on both elements)
    UTILS::EvaluateShape_Displ(sxi, sval,sele, false);
    UTILS::EvaluateShape_Displ(mxi, mval,mele, false);

    // evaluate the integration cell Jacobian
    double jac = cell->Jacobian();

    // compute cell D/M matrix *******************************************
    GP_DM(sele,mele,lmval,sval,mval,jac,wgt,nrow,ncol,ndof,bound,comm);

    // compute nodal scaling factor **************************************
    if (scaling)
      GP_3D_Scaling(sele,sval,jac,wgt,sxi);
    // compute nodal scaling factor **************************************
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Integrate and linearize a 2D slave / master cell (3D)     popp 03/09|
 |  This method integrates the cell M matrix (and possibly D matrix)    |
 |  and stores it in mseg and dseg respectively.                        |
 |  This is the QUADRATIC auxiliary plane coupling version!!!           |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void MORTAR::MortarIntegratorCalc<distypeS, distypeM>::IntegrateCell3DAuxPlaneQuad(
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    MORTAR::IntElement& sintele, MORTAR::IntElement& mintele,
    Teuchos::RCP<MORTAR::IntCell> cell, double* auxn)
{
  // get LMtype
  INPAR::MORTAR::LagMultQuad lmtype = lmquadtype_;

  // explicitly defined shape function type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivCell3DAuxPlaneQuad called without specific shape function defined!");

  //check for problem dimension
  if (ndim_!=3)
    dserror("ERROR: 3D integration method called for non-3D problem");

  if(cell->Shape()!=DRT::Element::tri3)
    dserror("ERROR: wrong cell shape!");

  // discretization type of slave and master IntElement
  DRT::Element::DiscretizationType sdt = sintele.Shape();
  DRT::Element::DiscretizationType mdt = mintele.Shape();

  // discretization type of slave and master Element
  DRT::Element::DiscretizationType psdt = sele.Shape();
  DRT::Element::DiscretizationType pmdt = mele.Shape();

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateDerivCell3DAuxPlaneQuad called on a wrong type of MortarElement pair!");
  if (cell==Teuchos::null)
    dserror("ERROR: IntegrateDerivCell3DAuxPlaneQuad called without integration cell");

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int nintrow = sintele.NumNode();
  int ndof = dynamic_cast<MORTAR::MortarNode*>(sele.Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  LINALG::Matrix<ns_, 1> sval;
  LINALG::Matrix<nm_, 1> mval;
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow, 2, true);
  LINALG::SerialDenseVector lmintval(nintrow);
  LINALG::SerialDenseMatrix lmintderiv(nintrow, 2, true);

  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if (!mynodes)
    dserror("ERROR: IntegrateDerivCell3DAuxPlaneQuad: Null pointer!");

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k = 0; k < nrow; ++k)
  {
    MORTAR::MortarNode* mymrtrnode =
        dynamic_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode)
      dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
    bound += mymrtrnode->IsOnBoundorCE();
  }

  // decide whether displacement shape fct. modification has to be considered or not
  // this is the case for dual quadratic and linear Lagrange multipliers on quad9/quad8/tri6 elements
  bool dualquad3d = false;
  if ((shapefcn_ == INPAR::MORTAR::shape_dual)
      && (lmtype == INPAR::MORTAR::lagmult_quad
          || lmtype == INPAR::MORTAR::lagmult_lin)
      && (sele.Shape() == DRT::Element::quad9
          || sele.Shape() == DRT::Element::quad8
          || sele.Shape() == DRT::Element::tri6))
  {
    dualquad3d = true;
  }

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < nGP(); ++gp)
  {
    // coordinates and weight
    double eta[2] = { Coordinate(gp, 0), Coordinate(gp, 1) };
    double wgt = Weight(gp);

    // get global Gauss point coordinates
    double globgp[3] = { 0.0, 0.0, 0.0 };
    cell->LocalToGlobal(eta, globgp, 0);

    double sxi[2] = { 0.0, 0.0 };
    double mxi[2] = { 0.0, 0.0 };

    // project Gauss point onto slave integration element
    // project Gauss point onto master integration element
    double sprojalpha = 0.0;
    double mprojalpha = 0.0;
    MORTAR::MortarProjector::Impl(sintele)->ProjectGaussPointAuxn3D(globgp,
        auxn, sintele, sxi, sprojalpha);
    MORTAR::MortarProjector::Impl(mintele)->ProjectGaussPointAuxn3D(globgp,
        auxn, mintele, mxi, mprojalpha);

    // check GP projection (SLAVE)
    double tol = 0.01;
    if (sdt == DRT::Element::quad4 || sdt == DRT::Element::quad8
        || sdt == DRT::Element::quad9)
    {
      if (sxi[0] < -1.0 - tol || sxi[1] < -1.0 - tol || sxi[0] > 1.0 + tol
          || sxi[1] > 1.0 + tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Slave Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1]                  << std::endl;
        std::cout << "Slave GP (IntElement) projection: " << sxi[0] << " " << sxi[1]       << std::endl;
      }
    }
    else
    {
      if (sxi[0] < -tol || sxi[1] < -tol || sxi[0] > 1.0 + tol
          || sxi[1] > 1.0 + tol || sxi[0] + sxi[1] > 1.0 + 2 * tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Slave Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1]                  << std::endl;
        std::cout << "Slave GP (IntElement) projection: " << sxi[0] << " " << sxi[1]       << std::endl;
      }
    }

    // check GP projection (MASTER)
    if (mdt == DRT::Element::quad4 || mdt == DRT::Element::quad8
        || mdt == DRT::Element::quad9)
    {
      if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[0] > 1.0 + tol
          || mxi[1] > 1.0 + tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Master Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1]                  << std::endl;
        std::cout << "Master GP (IntElement) projection: " << mxi[0] << " " << mxi[1]      << std::endl;
      }
    }
    else
    {
      if (mxi[0] < -tol || mxi[1] < -tol || mxi[0] > 1.0 + tol
          || mxi[1] > 1.0 + tol || mxi[0] + mxi[1] > 1.0 + 2 * tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Master Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1]                  << std::endl;
        std::cout << "Master GP (IntElement) projection: " << mxi[0] << " " << mxi[1]      << std::endl;
      }
    }

    // project Gauss point back to slave (parent) element
    // project Gauss point back to master (parent) element
    double psxi[2] = { 0.0, 0.0 };
    double pmxi[2] = { 0.0, 0.0 };
    double psprojalpha = 0.0;
    double pmprojalpha = 0.0;
    MORTAR::MortarProjector::Impl(sele)->ProjectGaussPointAuxn3D(globgp,
        auxn,sele,psxi,psprojalpha);
    MORTAR::MortarProjector::Impl(mele)->ProjectGaussPointAuxn3D(globgp,
        auxn,mele,pmxi,pmprojalpha);
    //sintele.MapToParent(sxi, psxi); // old way of doing it via affine map... wrong (popp 05/2016)
    //mintele.MapToParent(mxi, pmxi); // old way of doing it via affine map... wrong (popp 05/2016)

    // check GP projection (SLAVE)
    if (psdt == DRT::Element::quad4 || psdt == DRT::Element::quad8
        || psdt == DRT::Element::quad9 || psdt == DRT::Element::nurbs9)
    {
      if (psxi[0] < -1.0 - tol || psxi[1] < -1.0 - tol || psxi[0] > 1.0 + tol
          || psxi[1] > 1.0 + tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Slave Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1]                  << std::endl;
        std::cout << "Slave GP projection: " << psxi[0] << " " << psxi[1]       << std::endl;
      }
    }
    else
    {
      if (psxi[0] < -tol || psxi[1] < -tol || psxi[0] > 1.0 + tol
          || psxi[1] > 1.0 + tol || psxi[0] + psxi[1] > 1.0 + 2 * tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Slave Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1]                  << std::endl;
        std::cout << "Slave GP projection: " << psxi[0] << " " << psxi[1]       << std::endl;
      }
    }

    // check GP projection (MASTER)
    if (pmdt == DRT::Element::quad4 || pmdt == DRT::Element::quad8
        || pmdt == DRT::Element::quad9 || pmdt == DRT::Element::nurbs9)
    {
      if (pmxi[0] < -1.0 - tol || pmxi[1] < -1.0 - tol || pmxi[0] > 1.0 + tol
          || pmxi[1] > 1.0 + tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Master Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1]                  << std::endl;
        std::cout << "Master GP projection: " << pmxi[0] << " " << pmxi[1]      << std::endl;
      }
    }
    else
    {
      if (pmxi[0] < -tol || pmxi[1] < -tol || pmxi[0] > 1.0 + tol
          || pmxi[1] > 1.0 + tol || pmxi[0] + pmxi[1] > 1.0 + 2 * tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Master Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1]                  << std::endl;
        std::cout << "Master GP projection: " << pmxi[0] << " " << pmxi[1]      << std::endl;
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
    else if(bound)
    {
      sele.EvaluateShapeLagMult(shapefcn_, psxi, lmval, lmderiv, nrow);
    }
    else
    {
      sele.EvaluateShapeLagMult(shapefcn_, psxi, lmval, lmderiv, nrow);
      sintele.EvaluateShapeLagMult(shapefcn_, sxi, lmintval, lmintderiv, nintrow,false);
    }

    // evaluate trace space shape functions (on both elements)
    UTILS::EvaluateShape_Displ(psxi, sval, sele, dualquad3d);
    UTILS::EvaluateShape_Displ(pmxi, mval, mele, false);

    // evaluate the integration cell Jacobian
    double jac = cell->Jacobian();

    // compute cell D/M matrix *******************************************
    GP_3D_DM_Quad(sele, mele, sintele, lmval, lmintval, sval, mval, jac, wgt,
        nrow, nintrow, ncol, ndof, bound);
  }
  //**********************************************************************

  return;
}


//line2 slave
template class MORTAR::MortarIntegratorCalc<DRT::Element::line2, DRT::Element::line2>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::line2, DRT::Element::line3>;

//line3 slave
template class MORTAR::MortarIntegratorCalc<DRT::Element::line3, DRT::Element::line2>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::line3, DRT::Element::line3>;

//quad4 slave
template class MORTAR::MortarIntegratorCalc<DRT::Element::quad4, DRT::Element::quad4>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::quad4, DRT::Element::quad8>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::quad4, DRT::Element::quad9>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::quad4, DRT::Element::tri3>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::quad4, DRT::Element::tri6>;

//quad8 slave
template class MORTAR::MortarIntegratorCalc<DRT::Element::quad8, DRT::Element::quad4>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::quad8, DRT::Element::quad8>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::quad8, DRT::Element::quad9>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::quad8, DRT::Element::tri3>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::quad8, DRT::Element::tri6>;

//quad9 slave
template class MORTAR::MortarIntegratorCalc<DRT::Element::quad9, DRT::Element::quad4>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::quad9, DRT::Element::quad8>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::quad9, DRT::Element::quad9>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::quad9, DRT::Element::tri3>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::quad9, DRT::Element::tri6>;

//tri3 slave
template class MORTAR::MortarIntegratorCalc<DRT::Element::tri3, DRT::Element::quad4>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::tri3, DRT::Element::quad8>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::tri3, DRT::Element::quad9>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::tri3, DRT::Element::tri3>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::tri3, DRT::Element::tri6>;

//tri6 slave
template class MORTAR::MortarIntegratorCalc<DRT::Element::tri6, DRT::Element::quad4>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::tri6, DRT::Element::quad8>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::tri6, DRT::Element::quad9>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::tri6, DRT::Element::tri3>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::tri6, DRT::Element::tri6>;

//==================================================
//                     NURBS
//==================================================
//nurbs2 slave
template class MORTAR::MortarIntegratorCalc<DRT::Element::nurbs2, DRT::Element::nurbs2>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::nurbs2, DRT::Element::nurbs3>;

//nurbs3 slave
template class MORTAR::MortarIntegratorCalc<DRT::Element::nurbs3, DRT::Element::nurbs2>;
template class MORTAR::MortarIntegratorCalc<DRT::Element::nurbs3, DRT::Element::nurbs3>;
