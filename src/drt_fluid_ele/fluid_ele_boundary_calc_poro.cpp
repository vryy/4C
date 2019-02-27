/*!----------------------------------------------------------------------
\file fluid_ele_boundary_calc_poro.cpp

\brief evaluate boundary conditions not requiring parent-element evaluations

\level 2

\maintainer  Christoph Ager
*----------------------------------------------------------------------*/

#include "fluid_ele.H"
#include "../drt_lib/drt_element_integration_select.H"
#include "fluid_ele_action.H"
#include "fluid_ele_parameter_poro.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"

#include "../drt_geometry/position_array.H"
//
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_nurbs_discret/drt_nurbs_utils.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/fluidporo.H"
#include "../drt_mat/structporo.H"

#include "../drt_poroelast/poroelast_utils.H"

//#include "../drt_so3/so_poro_interface.H"

#include "../drt_lib/drt_discret.H"

// for dual shape functions
#include "../drt_volmortar/volmortar_shape.H"

#include "fluid_ele_boundary_calc_poro.H"

/*----------------------------------------------------------------------*
 * create/delete instance (public)                            vuong 06/11||
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>*
DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::Instance(bool create)
{
  static FluidEleBoundaryCalcPoro<distype>* instance;
  if (create)
  {
    if (instance == NULL)
    {
      instance = new FluidEleBoundaryCalcPoro<distype>();
    }
  }
  else
  {
    if (instance != NULL) delete instance;
    instance = NULL;
  }
  return instance;
}

/*----------------------------------------------------------------------*
 * destruction call (public)                                 vuong 06/11||
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since ele1 is a
  // cleanup call, we can do it this way.
  Instance(false);
}


/*----------------------------------------------------------------------*
 * constructor (protected)                                   vuong 06/11||
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::FluidEleBoundaryCalcPoro()
    : DRT::ELEMENTS::FluidBoundaryImpl<distype>::FluidBoundaryImpl()
{
  // pointer to class FluidImplParameterTimInt
  my::fldpara_ = DRT::ELEMENTS::FluidEleParameterPoro::Instance();
}

/*----------------------------------------------------------------------*
 * evaluate boundary action (public)                          vuong 06/11||
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::EvaluateAction(
    DRT::ELEMENTS::FluidBoundary* ele1, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::BoundaryAction act = DRT::INPUT::get<FLD::BoundaryAction>(params, "action");

  switch (act)
  {
    case FLD::no_penetration:
    {
      NoPenetration(ele1, params, discretization, lm, elemat1, elemat2, elevec1);
      break;
    }
    case FLD::no_penetrationIDs:
    {
      NoPenetrationIDs(ele1, params, discretization, elevec1, lm);
      break;
    }
    case FLD::poro_boundary:
    {
      PoroBoundary(ele1, params, discretization, lm, elemat1, elevec1);
      break;
    }
    case FLD::poro_prescoupl:
    {
      PressureCoupling(ele1, params, discretization, lm, elemat1, elevec1);
      break;
    }
    case FLD::fpsi_coupling:
    {
      // We skip all elements without any row nodes on this proc (will not contribute to the matrix
      // in the assembly of the matrix). Otherwise even fully ghosted Volume Elements would required
      // a ghosted Volume Element on the other side of the interface
      if (!ele1->HasOnlyGhostNodes(discretization.Comm().MyPID()))
        FPSICoupling(ele1, params, discretization, lm, elemat1, elevec1);
      break;
    }
    case FLD::calc_flowrate:
    {
      ComputeFlowRate(ele1, params, discretization, lm, elevec1);
      break;
    }
    case FLD::poro_splitnopenetration:
    {
      NoPenetrationMatAndRHS(ele1, params, discretization, lm, elemat1, elevec1);
      break;
    }
    case FLD::poro_splitnopenetration_OD:
    {
      NoPenetrationMatOD(ele1, params, discretization, lm, elemat1, elemat2);
      break;
    }
    case FLD::poro_splitnopenetration_ODpres:
    {
      NoPenetrationMatODPoroPres(ele1, params, discretization, lm, elemat1);
      break;
    }
    case FLD::poro_splitnopenetration_ODdisp:
    {
      NoPenetrationMatODPoroDisp(ele1, params, discretization, lm, elemat1);
      break;
    }
    default:
    {
      DRT::ELEMENTS::FluidBoundaryImpl<distype>::EvaluateAction(
          ele1, params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }
    break;
  }  // switch(act)
}

/*----------------------------------------------------------------------*
 * evaluate FPSI coupling condition (protected)               rauch 10/14|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::FPSICoupling(
    DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& plm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseVector& elevec1)
{
  switch (distype)
  {
    // 2D:
    case DRT::Element::line2:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::quad4)
      {
        this->FPSICoupling<DRT::Element::quad4>(ele, params, discretization, plm, elemat1, elevec1);
      }
      else
      {
        dserror(" expected combination line2/quad4 for surface/parent pair ");
      }
      break;
    }
    // 3D:
    case DRT::Element::quad4:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::hex8)
      {
        this->FPSICoupling<DRT::Element::hex8>(ele, params, discretization, plm, elemat1, elevec1);
      }
      else
      {
        dserror(" expected combination quad4/hex8 for surface/parent pair ");
      }
      break;
    }
    case DRT::Element::tri3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::tet4)
      {
        this->FPSICoupling<DRT::Element::tet4>(ele, params, discretization, plm, elemat1, elevec1);
      }
      else
      {
        dserror(" expected combination tri3/tet4 for surface/parent pair ");
      }
      break;
    }
    default:
    {
      dserror("surface/parent element pair not yet implemented. Just do it.\n");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 * evaluate FPSI coupling condition (protected)              rauch 10/14|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType pdistype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::FPSICoupling(
    DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& plm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseVector& elevec1)
{
  /*
   * Evaluate all Terms for the FPSI Boundary & Neumann Integration. The both conditions should be
   * splited into to methods later to avoid ifs in case of Neumann Integration!
   */

  /*
   rauch 01/2013

          /                  \
         |                    |
  (1)    |  (u - vs) o n , q  |             normal continuity of flux in porofluid equation
         |                    |
          \                  /  Gamma_Interface

          /                                                                \
         |                                                                  |
  (2)    |  J (tau - pf o I + gamma rho_f u dyadic u) o F^-T o N , delta d  |    equality of
  interface traction vector in structural equation | | \ /  Gamma_Interface

          /                                                          \
         |   1                                                        |
  (3)    | ------ n o (-pf o I - gamma rho_f u dyadic u) o n , w o n  |          equality of normal
  interface traction in fluid equation | rho_f | \ /  Gamma_Interface

          /                                                       \
         |  alphabj * mu_f                              I       I  |
  (4)    |  --------------- [u - (vs + phi(vf - vs))] o t , w o t  |             beavers-joseph
  condition in fluid equation |   rho_f sqrt(K)                                         | \ /
  Gamma_Interface


              nnod ->
             __ idof3 ->            __
     inod   |                         |
       idof2|                         |
        |   |                         |
      | V   |         elemat          |
      V     |                         |
            |                         |
            |                         |
            |__                     __|

   */


  // This function is only implemented for 3D
  if (my::bdrynsd_ != 2 and my::bdrynsd_ != 1)
  {
    dserror("Continuity boundary integral for FPSI coupling is only implemented for 3D and 2D!");
  }

  // number of parentnodes
  static const int nenparent = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

  // get the parent element
  DRT::ELEMENTS::Fluid* pele = ele->ParentElement();
  int currparenteleid = pele->Id();

  // get submatrix to fill
  const std::string block = params.get<std::string>("fillblock");

  // get map containing parent element facing current interface element
  const std::string tempstring("InterfaceFacingElementMap");
  Teuchos::RCP<std::map<int, int>> InterfaceFacingElementMap =
      params.get<Teuchos::RCP<std::map<int, int>>>(tempstring);
  std::map<int, int>::iterator it;

  // initialization of plenty of variables
  double fluiddynamicviscosity = 0.0;
  double permeability = 0.0;
  double beaversjosephcoefficient = 0.0;
  double normoftangential1 = 0.0;
  double normoftangential2 = 0.0;
  double normoftangential1_n = 0.0;
  double normoftangential2_n = 0.0;
  double scalarintegraltransformfac = 0.0;
  double tangentialfac = 0.0;

  LINALG::Matrix<my::nsd_, 1> neumannoverinflow(true);

  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;

  std::vector<double> my_displacements_np;
  std::vector<double> my_displacements_n;
  std::vector<double> my_parentdisp_np;
  std::vector<double> my_parentdisp_n;
  std::vector<double> porosity;

  LINALG::Matrix<my::nsd_, my::bdrynen_> evelnp(true);
  LINALG::Matrix<my::nsd_, my::bdrynen_> eveln(true);
  LINALG::Matrix<my::nsd_, nenparent> pevelnp(true);
  LINALG::Matrix<my::nsd_, nenparent> peveln(true);  // at previous time step n
  LINALG::Matrix<my::nsd_, my::bdrynen_> edispnp(true);
  LINALG::Matrix<my::nsd_, my::bdrynen_> egridvel(true);
  LINALG::Matrix<my::nsd_, my::bdrynen_> egridvel_n(true);
  LINALG::Matrix<1, my::bdrynen_> epressnp(true);
  LINALG::Matrix<1, my::bdrynen_> epressn(true);
  LINALG::Matrix<my::nsd_, 1> gridvelint(true);
  LINALG::Matrix<my::nsd_, 1> pxsi(true);
  LINALG::Matrix<1, 1> pressint(true);
  LINALG::Matrix<1, 1> pressint_n(true);  // at previous time step n
  LINALG::Matrix<my::nsd_, my::nsd_> dudxi(true);
  LINALG::Matrix<my::nsd_, my::nsd_> dudxi_n(true);  // at previous time step n
  LINALG::Matrix<my::nsd_, my::nsd_> dudxioJinv(true);
  LINALG::Matrix<my::nsd_, my::nsd_> dudxioJinv_n(true);  // at previous time step n
  LINALG::Matrix<1, 1> tangentialvelocity1(true);
  LINALG::Matrix<1, 1> tangentialvelocity2(true);
  LINALG::Matrix<1, 1> tangentialgridvelocity1(true);
  LINALG::Matrix<1, 1> tangentialgridvelocity2(true);
  LINALG::Matrix<1, 1> normalvelocity(true);

  LINALG::Matrix<my::nsd_, nenparent> xrefe;  // material coord. of parent element
  LINALG::Matrix<my::nsd_, nenparent> xcurr;  // current  coord. of parent element
  LINALG::Matrix<my::nsd_, nenparent>
      xcurr_n;  // current  coord. of parent element at previous time step n

  Teuchos::RCP<const Epetra_Vector> displacements_np = discretization.GetState("dispnp");
  Teuchos::RCP<const Epetra_Vector> displacements_n = discretization.GetState("dispn");
  Teuchos::RCP<const Epetra_Vector> fluidvelocity_np = discretization.GetState("velnp");
  Teuchos::RCP<const Epetra_Vector> fluidvelocity_n = discretization.GetState("veln");
  Teuchos::RCP<const Epetra_Vector> gridvelocity = discretization.GetState("gridv");

  if (fluidvelocity_np == Teuchos::null) dserror("Cannot get state vector 'fluidvelocity_np'");
  if (gridvelocity == Teuchos::null) dserror("Cannot get state vector 'gridvelocity'");
  if (displacements_np == Teuchos::null) dserror("Cannot get state vector 'displacements_np'");
  if (fluidvelocity_n == Teuchos::null) dserror("Cannot get state vector 'fluidvelocity_n'");
  if (displacements_n == Teuchos::null) dserror("Cannot get state vector 'displacements_n'");

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<my::bdrynsd_> intpoints(
      DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);
  //(const DRT::UTILS::IntPointsAndWeights<my::nsd_>
  // pintpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<pdistype>::rule);

  // get node coordinates
  // (we have a my::nsd_ dimensional domain, since my::nsd_ determines the dimension of
  // FluidBoundary element!)
  GEO::fillInitialPositionArray<distype, my::nsd_, LINALG::Matrix<my::nsd_, my::bdrynen_>>(
      ele, my::xyze_);
  GEO::fillInitialPositionArray<distype, my::nsd_, LINALG::Matrix<my::nsd_, my::bdrynen_>>(
      ele, my::xyze_n_);

  // get element location vector and ownerships
  ele->DRT::Element::LocationVector(discretization, lm, lmowner, lmstride);

  // get material parameters and constants needed to calculate matrix terms
  const Teuchos::ParameterList& fpsidynparams = DRT::Problem::Instance()->FPSIDynamicParams();

  Teuchos::RCP<MAT::Material> fluidmaterial;
  Teuchos::RCP<MAT::Material> generalmaterial;
  Teuchos::RCP<MAT::Material> currentmaterial;
  Teuchos::RCP<MAT::FluidPoro> porofluidmaterial;
  Teuchos::RCP<MAT::NewtonianFluid> newtonianfluidmaterial;

  currentmaterial = ele->ParentElement()->Material();

  if (discretization.Name() == "fluid")
  {
    if (block != "NeumannIntegration" && block != "NeumannIntegration_Ale")
    // InterfaceFacingElementMap in general does not have elements on NeumannIntegration
    //(just on the FPSI interface)
    {
      Teuchos::RCP<DRT::Discretization> porofluiddis =
          DRT::Problem::Instance()->GetDis("porofluid");
      it = InterfaceFacingElementMap->find(ele->Id());
      if (it == InterfaceFacingElementMap->end())
        dserror("Couldn't find ele %d in InterfaceFacingElementMap", ele->Id());

      DRT::Element* porofluidelement = porofluiddis->gElement(it->second);

      generalmaterial = porofluidelement->Material();
      porofluidmaterial = Teuchos::rcp_dynamic_cast<MAT::FluidPoro>(generalmaterial);
      permeability = porofluidmaterial->Permeability();
    }

    newtonianfluidmaterial = Teuchos::rcp_dynamic_cast<MAT::NewtonianFluid>(currentmaterial);

    fluiddynamicviscosity = newtonianfluidmaterial->Viscosity();
  }
  else if (discretization.Name() == "porofluid")
  {
    Teuchos::RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->GetDis("fluid");
    it = InterfaceFacingElementMap->find(ele->Id());
    if (it == InterfaceFacingElementMap->end())
      dserror("Couldn't find ele %d in InterfaceFacingElementMap", ele->Id());

    DRT::Element* fluidelement = fluiddis->gElement(it->second);

    fluidmaterial = fluidelement->Material();
    newtonianfluidmaterial = Teuchos::rcp_dynamic_cast<MAT::NewtonianFluid>(fluidmaterial);
    porofluidmaterial = Teuchos::rcp_dynamic_cast<MAT::FluidPoro>(currentmaterial);

    permeability = porofluidmaterial->Permeability();
    fluiddynamicviscosity = newtonianfluidmaterial->Viscosity();
  }

  if (block != "NeumannIntegration" && block != "NeumannIntegration_Ale")
  {
    // InterfaceFacingElementMap in general does not have elements on NeumannIntegration
    // calculate factor for the tangential interface condition on the free fluid field
    beaversjosephcoefficient = fpsidynparams.get<double>("ALPHABJ");
    tangentialfac = (beaversjosephcoefficient * fluiddynamicviscosity) / (sqrt(permeability));
  }

  const double timescale = params.get<double>("timescale", -1.0);
  if (timescale == -1.0) dserror("no timescale parameter in parameter list");

  if (displacements_np != Teuchos::null)
  {
    my_displacements_np.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*displacements_np, my_displacements_np, lm);
    DRT::UTILS::ExtractMyValues(*displacements_np, my_parentdisp_np, plm);
  }
  dsassert(my_displacements_np.size() != 0, "no displacement values for boundary element");
  dsassert(my_parentdisp_np.size() != 0, "no displacement values for parent element");

  if (displacements_n != Teuchos::null)
  {
    my_displacements_n.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*displacements_n, my_displacements_n, lm);
    DRT::UTILS::ExtractMyValues(*displacements_n, my_parentdisp_n, plm);
  }
  dsassert(
      my_displacements_n.size() != 0, "no displacement values for boundary element at time step n");
  dsassert(my_parentdisp_n.size() != 0, "no displacement values for parent element at time step n");

  // Add the deformation of the ALE mesh to the nodes coordinates
  for (int inode = 0; inode < my::bdrynen_; ++inode)
  {
    for (int idim = 0; idim < my::nsd_; ++idim)
    {
      my::xyze_(idim, inode) += my_displacements_np[my::numdofpernode_ * inode + idim];
      my::xyze_n_(idim, inode) += my_displacements_n[my::numdofpernode_ * inode + idim];
    }
  }

  // update element geometry of parent element
  {
    DRT::Node** nodes = pele->Nodes();
    for (int inode = 0; inode < nenparent; ++inode)
    {
      for (int idof = 0; idof < my::nsd_; ++idof)
      {
        const double* x = nodes[inode]->X();
        xrefe(idof, inode) = x[idof];
        xcurr(idof, inode) =
            xrefe(idof, inode) + my_parentdisp_np[inode * my::numdofpernode_ + idof];
        xcurr_n(idof, inode) =
            xrefe(idof, inode) + my_parentdisp_n[inode * my::numdofpernode_ + idof];
      }
    }
  }

  // extract local values from the global vectors
  std::vector<double> my_fluidvelocity_np(lm.size());
  DRT::UTILS::ExtractMyValues(*fluidvelocity_np, my_fluidvelocity_np, lm);
  std::vector<double> my_fluidvelocity_n(lm.size());  // at previous time step n
  DRT::UTILS::ExtractMyValues(*fluidvelocity_n, my_fluidvelocity_n, lm);
  std::vector<double> my_gridvelocity(lm.size());
  DRT::UTILS::ExtractMyValues(*gridvelocity, my_gridvelocity, lm);
  std::vector<double> my_parentfluidvelocity_np(plm.size());
  DRT::UTILS::ExtractMyValues(*fluidvelocity_np, my_parentfluidvelocity_np, plm);
  std::vector<double> my_parentfluidvelocity_n(plm.size());  // at previous time step n
  DRT::UTILS::ExtractMyValues(*fluidvelocity_n, my_parentfluidvelocity_n, plm);

  // split velocity and pressure, insert into element arrays
  for (int inode = 0; inode < my::bdrynen_; inode++)
  {
    for (int idim = 0; idim < my::nsd_; idim++)
    {
      evelnp(idim, inode) = my_fluidvelocity_np[idim + (inode * my::numdofpernode_)];
      eveln(idim, inode) = my_fluidvelocity_n[idim + (inode * my::numdofpernode_)];
      edispnp(idim, inode) = my_displacements_np[idim + (inode * my::numdofpernode_)];
      egridvel(idim, inode) = my_gridvelocity[idim + (inode * my::numdofpernode_)];
    }
    epressnp(inode) = my_fluidvelocity_np[my::nsd_ + (my::numdofpernode_ * inode)];
    epressn(inode) = my_fluidvelocity_n[my::nsd_ + (my::numdofpernode_ * inode)];
  }

  for (int inode = 0; inode < nenparent; inode++)
  {
    for (int idim = 0; idim < my::nsd_; idim++)
    {
      pevelnp(idim, inode) = my_parentfluidvelocity_np[idim + (inode * my::numdofpernode_)];
      peveln(idim, inode) = my_parentfluidvelocity_n[idim + (inode * my::numdofpernode_)];
    }
  }

  // get porosity values from parent element
  Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;

  // access structure discretization
  structdis = DRT::Problem::Instance()->GetDis("structure");
  // std::cout<<structdis<<endl;
  DRT::Element* structele = NULL;
  // get corresponding structure element (it has the same global ID as the porofluid element)
  if (discretization.Name() == "structure" or discretization.Name() == "porofluid")
  {
    structele = structdis->gElement(currparenteleid);
  }
  else if (discretization.Name() == "fluid" && block != "NeumannIntegration" &&
           block != "NeumannIntegration_Ale")
  {
    it = InterfaceFacingElementMap->find(ele->Id());
    structele = structdis->gElement(it->second);
  }

  if (structele == NULL && block != "NeumannIntegration" && block != "NeumannIntegration_Ale")
  {
    dserror("Structure element %i not on local processor", currparenteleid);
  }

  // get porous material
  Teuchos::RCP<MAT::StructPoro> structmat;
  if (block != "NeumannIntegration" && block != "NeumannIntegration_Ale")
  {
    structmat = Teuchos::rcp_dynamic_cast<MAT::StructPoro>(structele->Material());
    if (structmat->MaterialType() != INPAR::MAT::m_structporo)
    {
      dserror("invalid structure material for poroelasticity");
    }
  }

  // what's the current problem type?
  PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();
  double Lp = 0.0;
  if (probtype == prb_fps3i)
  {
    // get the conductivity of membrane at the interface
    Lp = params.get<double>("membrane conductivity");
  }

  // get coordinates of gauss points w.r.t. local parent coordinate system
  LINALG::SerialDenseMatrix pqxg(intpoints.IP().nquad, my::nsd_);
  LINALG::Matrix<my::nsd_, my::nsd_> derivtrafo(true);

  DRT::UTILS::BoundaryGPToParentGP<my::nsd_>(
      pqxg, derivtrafo, intpoints, pdistype, distype, ele->SurfaceNumber());

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////     Loop over Gauss-Points    /////////////////////
  ////////////////////////////////////////////////////////////////////////////
  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    // get shape functions and derivatives in the plane of the element
    LINALG::Matrix<nenparent, 1> pfunct(true);  // parent element shape function
    LINALG::Matrix<my::nsd_, nenparent> pderiv(
        true);  // derivatives of parent element shape functions in interface coordinate system
    LINALG::Matrix<my::nsd_, nenparent> pderiv_loc(
        true);  // derivatives of parent element shape functions in parent element coordinate system

    // coordinates of the current integration point in parent coordinate system
    for (int idim = 0; idim < my::nsd_; idim++)
    {
      pxsi(idim) = pqxg(gpid, idim);
    }

    // evalute parent element shape function at current integration point in parent coordinate
    // system
    DRT::UTILS::shape_function<pdistype>(pxsi, pfunct);
    // evaluate derivatives of parent element shape functions at current integration point in parent
    // coordinate system
    DRT::UTILS::shape_function_deriv1<pdistype>(pxsi, pderiv_loc);
    // transformation from parent element coordinate system to interface element coordinate system
    pderiv.MultiplyTN(derivtrafo, pderiv_loc);

    double dphi_dp = 0.0;
    double dphi_dJ = 0.0;
    double dphi_dJdp = 0.0;
    double dphi_dJJ = 0.0;
    double dphi_dpp = 0.0;
    double porosityint = 0.0;

    // get Jacobian matrix and determinant w.r.t. spatial configuration
    //
    // |J| = det(xjm) * det(Jmat^-1) = det(xjm) * 1/det(Jmat)
    //
    //    _                     _
    //   |  x_1,1  x_2,1  x_3,1  |           d x_i
    //   |  x_1,2  x_2,2  x_3,2  | = xjm  = --------
    //   |_ x_1,3  x_2,3  x_3,3 _|           d s_j
    //    _
    //   |  X_1,1  X_2,1  X_3,1  |           d X_i
    //   |  X_1,2  X_2,2  X_3,2  | = Jmat = --------
    //   |_ X_1,3  X_2,3  X_3,3 _|           d s_j
    //
    LINALG::Matrix<my::nsd_, my::nsd_> xjm;
    LINALG::Matrix<my::nsd_, my::nsd_> xjm_n;  // at previous time step n
    LINALG::Matrix<my::nsd_, my::nsd_> Jmat;
    xjm.MultiplyNT(pderiv_loc, xcurr);
    xjm_n.MultiplyNT(pderiv_loc, xcurr_n);
    Jmat.MultiplyNT(pderiv_loc, xrefe);
    double det = xjm.Determinant();
    double detJ = Jmat.Determinant();
    const double J = det / detJ;

    // inverse of transposed jacobian "ds/dx" (xjm)
    LINALG::Matrix<my::nsd_, my::nsd_> xji;
    LINALG::Matrix<my::nsd_, my::nsd_> xji_n;  // at previous time step n
    //    _                     _
    //   |  s_1,1  s_2,1  s_3,1  |           d s_i
    //   |  s_1,2  s_2,2  s_3,2  | = xji  = -------- ;  [xji] o [xjm] = I
    //   |_ s_1,3  s_2,3  s_3,3 _|           d x_j
    //    _
    xji.Invert(xjm);
    xji_n.Invert(xjm_n);

#ifdef DEBUG
    // check unitiy of  [xji] o [xjm]
    LINALG::Matrix<my::nsd_, my::nsd_> eye;
    eye.Multiply(xji, xjm);
    if (my::nsd_ == 3)
    {
      if (abs(eye(0, 0) - 1.0) > 1e-11 or abs(eye(1, 1) - 1.0) > 1e-11 or
          abs(eye(2, 2) - 1.0) > 1e-11)
      {
        std::cout << eye << std::endl;
        dserror("matrix times its inverse is not equal identity ... that sucks !!!");
      }
      if (abs(eye(0, 1)) > 1e-11 or abs(eye(0, 2)) > 1e-11 or abs(eye(1, 0)) > 1e-11 or
          abs(eye(1, 2)) > 1e-11 or abs(eye(2, 0)) > 1e-11 or abs(eye(2, 1)) > 1e-11)
      {
        std::cout << eye << std::endl;
        dserror("matrix times its inverse is not equal identity ... that sucks !!!");
      }
    }
    else if (my::nsd_ == 2)
    {
      if (abs(eye(0, 0) - 1.0) > 1e-11 or abs(eye(1, 1) - 1.0) > 1e-11)
      {
        std::cout << eye << std::endl;
        dserror("matrix times its inverse is not equal identity ... that sucks !!!");
      }
      if (abs(eye(0, 1)) > 1e-11 or abs(eye(1, 0)) > 1e-11)
      {
        std::cout << eye << std::endl;
        dserror("matrix times its inverse is not equal identity ... that sucks !!!");
      }
    }
#endif

    // evaluate my::unitnormal_ , my::deriv_, ...
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(my::funct_, my::deriv_, my::fac_,
        my::unitnormal_n_, my::drs_, my::xsi_, my::xyze_n_, intpoints, gpid, NULL, NULL,
        IsNurbs<distype>::isnurbs);

    // evaluate my::unitnormal_ , my::deriv_, ...
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(my::funct_, my::deriv_, my::fac_,
        my::unitnormal_, my::drs_, my::xsi_, my::xyze_, intpoints, gpid, NULL, NULL,
        IsNurbs<distype>::isnurbs);

    // my::fac_ = intpoints.IP().qwgt[gpid]*my::drs_ done in EvalShapeFuncAtBouIntPoint()
    const double timefac = my::fldparatimint_->TimeFac();
    const double timefacpre = my::fldparatimint_->TimeFacPre();
    const double timefacfacpre = my::fldparatimint_->TimeFacPre() * my::fac_;
    const double rhsfac = my::fldparatimint_->TimeFacRhs() * my::fac_;
    // const double rhsfacpre     = my::fldparatimint_->TimeFacRhsPre() * my::fac_;
    const double theta = my::fldparatimint_->Theta();

    // The integration factor is not multiplied with drs
    // since it is the same as the scaling factor for the unit normal derivatives
    // Therefore it cancels out!!
    const double fac = intpoints.IP().qwgt[gpid];

    // calculate variables at gausspoint
    my::velint_.Multiply(evelnp, my::funct_);
    gridvelint.Multiply(egridvel, my::funct_);
    pressint.Multiply(epressnp, my::funct_);
    pressint_n.Multiply(epressn, my::funct_);

    //                                         _              _
    //                                        | u1,1 u1,2 u1,3 |
    // dudxi = u_i,alhpa = N_A,alpha u^A_i =  | u2,1 u2,2 u2,3 |
    //                                        |_u3,1 u3,2 u3,3_|
    //
    dudxi.MultiplyNT(pevelnp, pderiv_loc);  // corrected: switched pevelnp and pderiv
    dudxi_n.MultiplyNT(peveln, pderiv_loc);

    //                                            l=_  1     2     3  _
    //         -1                               i=1| u1,x1 u1,x2 u1,x3 |
    // dudxi o J  = N_A,alpha u^A_i xi_alpha,l =  2| u2,x1 u2,x2 u2,x3 | = gradu
    //                                            3|_u3,x1 u3,x2 u3,x3_|
    //
    dudxioJinv.MultiplyNT(dudxi, xji);
    dudxioJinv_n.MultiplyNT(dudxi_n, xji_n);  // at previus time step n

    LINALG::Matrix<1, my::nsd_> graduon(true);
    LINALG::Matrix<1, my::nsd_> graduon_n(true);  // from previous time step
    //
    // l=  1     2     3
    // [  ...   ...   ...  ]
    //
    //
    for (int idof = 0; idof < my::nsd_; idof++)  // l Loop
    {
      for (int idof2 = 0; idof2 < my::nsd_; idof2++)
      {
        graduon(0, idof) += dudxioJinv(idof, idof2) * my::unitnormal_(idof2);
        graduon_n(0, idof) += dudxioJinv_n(idof, idof2) * my::unitnormal_n_(idof2);
      }
    }
    LINALG::Matrix<1, my::nsd_> graduTon(true);
    LINALG::Matrix<1, my::nsd_> graduTon_n(true);  // at previous time step n
    //
    // l=  1     2     3
    // [  ...   ...   ...  ]
    //
    //
    for (int idof = 0; idof < my::nsd_; idof++)  // l Loop
    {
      for (int idof2 = 0; idof2 < my::nsd_; idof2++)
      {
        graduTon(0, idof) += dudxioJinv(idof2, idof) * my::unitnormal_(idof2);
        graduTon_n(0, idof) += dudxioJinv_n(idof2, idof) * my::unitnormal_n_(idof2);
      }
    }

    if (discretization.Name() == "porofluid" or discretization.Name() == "structure")
      structmat->ComputeSurfPorosity(params, pressint(0, 0), J, ele->SurfaceNumber(), gpid,
          porosityint, &dphi_dp, &dphi_dJ, &dphi_dJdp, &dphi_dJJ, &dphi_dpp, false);
    else
      porosityint = 1.0;

#ifdef DEBUG
    if (porosityint < 0.00001)
    {
      std::cout << "Discretization: " << discretization.Name() << std::endl;
      std::cout << "SurfaceNumber:  " << ele->SurfaceNumber() << std::endl;
      std::cout << "Porosity:       " << porosityint << "  at gp: " << gpid << std::endl;
      std::cout << "Pressure at gp: " << pressint(0, 0) << std::endl;
      std::cout << "Jacobian:       " << J << std::endl;
      dserror("unreasonably low porosity for poro problem");
    }
#endif

    // dxyzdrs vector -> normal which is not normalized built from cross product of columns
    // of Jacobian matrix d(x,y,z)/d(r,s)
    LINALG::Matrix<my::bdrynsd_, my::nsd_> dxyzdrs(0.0);
    LINALG::Matrix<my::bdrynsd_, my::nsd_> dxyzdrs_n(0.0);
    dxyzdrs.MultiplyNT(my::deriv_, my::xyze_);
    dxyzdrs_n.MultiplyNT(my::deriv_, my::xyze_n_);

    // tangential surface vectors are columns of dxyzdrs
    LINALG::Matrix<my::nsd_, 1> tangential1(true);
    LINALG::Matrix<my::nsd_, 1> tangential2(true);
    LINALG::Matrix<my::nsd_, 1> tangential1_n(true);
    LINALG::Matrix<my::nsd_, 1> tangential2_n(true);

    for (int idof = 0; idof < my::nsd_; idof++)
    {
      tangential1(idof, 0) = dxyzdrs(0, idof);
      tangential1_n(idof, 0) = dxyzdrs_n(0, idof);
    }

    normoftangential1 = tangential1.Norm2();
    normoftangential1_n = tangential1_n.Norm2();

    // normalize tangential vectors
    tangential1.Scale(1 / normoftangential1);

    tangential1_n.Scale(1 / normoftangential1_n);

    if (my::nsd_ == 3)
    {
      for (int idof = 0; idof < my::nsd_; idof++)
      {
        tangential2(idof, 0) = dxyzdrs(1, idof);
        tangential2_n(idof, 0) = dxyzdrs_n(1, idof);
      }

      normoftangential2 = tangential2.Norm2();
      normoftangential2_n = tangential2_n.Norm2();

      // normalize tangential vectors
      tangential2.Scale(1 / normoftangential2);

      tangential2_n.Scale(1 / normoftangential2_n);
    }

    //                                                             I
    // calculate tangential structure velocity (gridvelocity) vs o t
    //
    // [my::nsd_ x 1] o [my::nsd_ x 1]
    //
    double tangentialvs1 = 0.0;
    double tangentialvs2 = 0.0;
    tangentialvs1 = gridvelint.Dot(tangential1);
    tangentialvs2 = gridvelint.Dot(tangential2);

    //                                          I
    // calculate tangential fluid velocity vf o t
    //
    // [my::nsd_ x 1] o [my::nsd_ x 1]
    //
    double tangentialvf1 = 0.0;
    double tangentialvf2 = 0.0;
    tangentialvf1 = my::velint_.Dot(tangential1);
    tangentialvf2 = my::velint_.Dot(tangential2);

    //  derivatives of surface tangentials with respect to mesh displacements
    //              I
    //            d t_i             I                               I   I
    //            -------- = 1/abs( t )* (N_L,(r,s) Kronecker^i_l - t_i t_l N_L,(r,s) )
    //            d d^L_l
    //
    //         _______________L=1_____________    ______________L=2_____________   ______ ...
    //     __ /l =  1         2         3     \  /l = 1          2        3     \ /       __
    //  i= |                                    |                                |          |
    //  t1 |  N_1,(r,s)-() -(...)      -(...)   |  N_2,(r,s)   ...       ...     |  ...     |
    //     |                                    |                                |          |
    //  t2 |  -(...)     N_1,(r,s)-()  -(...)   |    ...      N_2,(r,s)  ...     |  ...     |
    //     |                                    |                                |          |
    //  t3 |  -(...)     -(...)    N_1,(r,s)-() |    ...       ...     N_2,(r,s) |  ...     |
    //     |_                                                                              _|
    //
    LINALG::Matrix<my::nsd_, nenparent * my::nsd_> tangentialderiv1(true);
    LINALG::Matrix<my::nsd_, nenparent * my::nsd_> tangentialderiv2(true);

    for (int node = 0; node < nenparent; ++node)
    {
      // block diagonal entries
      for (int idof = 0; idof < my::nsd_; ++idof)
        tangentialderiv1(idof, (node * my::nsd_) + idof) = pderiv(0, node) / normoftangential1;

      // terms from linearization of norm
      for (int idof = 0; idof < my::nsd_; ++idof)
        for (int idof2 = 0; idof2 < my::nsd_; idof2++)
          tangentialderiv1(idof, (node * my::nsd_) + idof2) -=
              (tangential1(idof, 0) * tangential1(idof2, 0) * pderiv(0, node)) / normoftangential1;
    }
    if (my::nsd_ == 3)
    {
      for (int node = 0; node < nenparent; ++node)
      {
        // block diagonal entries
        for (int idof = 0; idof < my::nsd_; ++idof)
          tangentialderiv2(idof, (node * my::nsd_) + idof) = pderiv(1, node) / normoftangential2;

        // terms from linearization of norm
        for (int idof = 0; idof < my::nsd_; ++idof)
          for (int idof2 = 0; idof2 < my::nsd_; idof2++)
            tangentialderiv2(idof, (node * my::nsd_) + idof2) -=
                (tangential2(idof, 0) * tangential2(idof2, 0) * pderiv(1, node)) /
                normoftangential2;
      }
    }

    //          I        ___L=1___  __L=2___  ___ ...
    //        d t_j     /l=1 2 3  \/l=1 2 3 \/
    // vs_j --------- = [  x x x      x x x            ]
    //       d d^L_l
    //
    LINALG::Matrix<nenparent * my::nsd_, 1> vsotangentialderiv1(true);
    LINALG::Matrix<nenparent * my::nsd_, 1> vsotangentialderiv2(true);
    for (int inode = 0; inode < nenparent; inode++)
    {
      for (int idof = 0; idof < my::nsd_; idof++)
      {
        for (int idof2 = 0; idof2 < my::nsd_; idof2++)
        {
          vsotangentialderiv1((inode * my::nsd_) + idof, 0) +=
              gridvelint(idof2, 0) * tangentialderiv1(idof2, (inode * my::nsd_) + idof);
          vsotangentialderiv2((inode * my::nsd_) + idof, 0) +=
              gridvelint(idof2, 0) * tangentialderiv2(idof2, (inode * my::nsd_) + idof);
        }
      }
    }
    LINALG::Matrix<nenparent * my::nsd_, 1> vfotangentialderiv1(true);
    LINALG::Matrix<nenparent * my::nsd_, 1> vfotangentialderiv2(true);
    for (int inode = 0; inode < nenparent; inode++)
    {
      for (int idof = 0; idof < my::nsd_; idof++)
      {
        for (int idof2 = 0; idof2 < my::nsd_; idof2++)
        {
          vfotangentialderiv1((inode * my::nsd_) + idof, 0) +=
              my::velint_(idof2, 0) * tangentialderiv1(idof2, (inode * my::nsd_) + idof);
          vfotangentialderiv2((inode * my::nsd_) + idof, 0) +=
              my::velint_(idof2, 0) * tangentialderiv2(idof2, (inode * my::nsd_) + idof);
        }
      }
    }


    //  derivatives of surface normals with respect to mesh displacements:
    //                                 d n_i
    //                                --------
    //                                 d d^L_l
    //
    //  parent element shape functions are used because the matrix normalderiv
    //  must have the proper dimension to be compatible to the evaluation of
    //  the matrix terms. as built below the matrix normalderiv has more entries
    //  than needed to calculate the surface integrals since the derivatives of
    //  the parent element shape functions do not necessarily vanish at the boundary
    //  gauss points. later those additional entries are however multiplied by the
    //  weighting function in those gauss points which are only different from zero
    //  when they belong to an interface node. thus all terms not belonging to the
    //  interface and its corresponding basic functions become zero. this makes perfect
    //  sense for the normal and its linearization are well determined solely by the
    //  surface of the element.
    LINALG::Matrix<my::nsd_, nenparent * my::nsd_> normalderiv(true);

    if (my::nsd_ == 3)
      for (int node = 0; node < nenparent; ++node)
      {
        normalderiv(0, 3 * node) += 0.;
        normalderiv(0, 3 * node + 1) +=
            (pderiv(0, node) * dxyzdrs(1, 2) - pderiv(1, node) * dxyzdrs(0, 2));
        normalderiv(0, 3 * node + 2) +=
            (pderiv(1, node) * dxyzdrs(0, 1) - pderiv(0, node) * dxyzdrs(1, 1));

        normalderiv(1, 3 * node) +=
            (pderiv(1, node) * dxyzdrs(0, 2) - pderiv(0, node) * dxyzdrs(1, 2));
        normalderiv(1, 3 * node + 1) += 0.;
        normalderiv(1, 3 * node + 2) +=
            (pderiv(0, node) * dxyzdrs(1, 0) - pderiv(1, node) * dxyzdrs(0, 0));

        normalderiv(2, 3 * node) +=
            (pderiv(0, node) * dxyzdrs(1, 1) - pderiv(1, node) * dxyzdrs(0, 1));
        normalderiv(2, 3 * node + 1) +=
            (pderiv(1, node) * dxyzdrs(0, 0) - pderiv(0, node) * dxyzdrs(1, 0));
        normalderiv(2, 3 * node + 2) += 0.;
      }
    else
      for (int node = 0; node < nenparent; ++node)
      {
        normalderiv(0, my::nsd_ * node) += 0.;
        normalderiv(0, my::nsd_ * node + 1) += pderiv(0, node);

        normalderiv(1, my::nsd_ * node) += -pderiv(0, node);
        normalderiv(1, my::nsd_ * node + 1) += 0.;
      }


    // dxyzdrs(0,:) x dxyzdrs(1,:) non unit normal
    //           _     _       _     _
    //          |       |     |       |
    //          | x_1,r |     | x_1,s |
    //          |       |     |       |
    //          | x_2,r |  X  | x_2,s |
    //          |       |     |       |
    //          | x_3,r |     | x_3,s |
    //          |_     _|     |_     _|
    //
    LINALG::Matrix<my::nsd_, 1> normal(true);

    if (my::nsd_ == 3)
    {
      normal(0, 0) = dxyzdrs(0, 1) * dxyzdrs(1, 2) - dxyzdrs(0, 2) * dxyzdrs(1, 1);
      normal(1, 0) = dxyzdrs(0, 2) * dxyzdrs(1, 0) - dxyzdrs(0, 0) * dxyzdrs(1, 2);
      normal(2, 0) = dxyzdrs(0, 0) * dxyzdrs(1, 1) - dxyzdrs(0, 1) * dxyzdrs(1, 0);
    }
    else
    {
      normal(0, 0) = dxyzdrs(0, 1);
      normal(1, 0) = -dxyzdrs(0, 0);
    }
    // transformation factor for surface integrals without normal vector
    scalarintegraltransformfac = normal.Norm2();  // || x,r x x,s ||

    // linearization of || x,r x x,s || = ||n||
    //
    //                L=__                           1 2        ...     nenparent __
    //  d ||n||    l=  | |          |        |             |
    //  ------- :   1  |1/||n||*(n_2*(x_3,1 N_L,2 - x_3,2 N_L,1) + n_3*(x_2,2 N_L,1 - x_2,1 N_L,2))
    //  |          |        |             | d d^L_l     2  |1/||n||*(n_1*(x_3,2 N_L,1 - x_3,1 N_L,2)
    //  + n_3*(x_1,1 N_L,2 - x_1,2 N_L,1))    |          |        |             |
    //              3  |1/||n||*(n_1*(x_2,1 N_L,2 - x_2,2 N_L,1) + n_2*(x_1,2 N_L,1 - x_1,1 N_L,2))
    //              |          |        |             |
    //                 |_ |          |        |            _|
    //
    //
    LINALG::Matrix<my::nsd_, nenparent> linearizationofscalarintegraltransformfac(true);

    for (int node = 0; node < nenparent; ++node)
    {
      for (int ldof = 0; ldof < my::nsd_; ++ldof)
      {
        for (int idof = 0; idof < my::nsd_; ++idof)
        {
          linearizationofscalarintegraltransformfac(ldof, node) +=
              1 / scalarintegraltransformfac * normal(idof, 0) *
              normalderiv(idof, node * my::nsd_ + ldof);
        }
      }
    }


    //------------------------------------- d|J|/dd = d|J|/dF : dF/dd = |J| * F^-T . N_X = |J| * N_x
    //
    // linearization of jacobian determinant w.r.t. structural displacements
    LINALG::Matrix<1, my::nsd_ * nenparent> dJ_dds;
    // global derivatives of shape functions w.r.t x,y,z (material configuration)
    LINALG::Matrix<my::nsd_, nenparent> derxy;

    //                                        _                          _
    //            d  N_A      d xi_alpha     |  N1,1 N2,1 N3,1 N4,1 ...   |
    //  derxy  = ----------  ----------- =   |  N1,2 N2,2 N3,2 N4,2 ...   |
    //            d xi_alpha  d   x_j        |_ N1,3 N2,3 N3,3 N4,3 ...  _|
    //
    derxy.Multiply(xji, pderiv_loc);

    for (int i = 0; i < nenparent; i++)
      for (int j = 0; j < my::nsd_; j++) dJ_dds(j + i * my::nsd_) = J * derxy(j, i);

    //
    //
    //            d xi_beta
    //  N_L,beta  ---------- n^j = derxy o n
    //            d   x_j
    //
    LINALG::Matrix<1, nenparent> dNdxon(true);
    for (int inode = 0; inode < nenparent; inode++)
    {
      for (int idof = 0; idof < my::nsd_; idof++)
      {
        dNdxon(0, inode) += derxy(idof, inode) * my::unitnormal_(idof);
      }
    }

    // LINALG::Matrix<1,my::nsd_*nenparent> gradNon;
    LINALG::Matrix<1, nenparent> gradNon(true);
    LINALG::Matrix<1, my::nsd_ * nenparent> gradN(true);
    //              d xi_alpha
    //  N_L,alpha  ------------ [g_L x g_j]
    //              d  x_j
    //
    //      ___L=1___  __L=2___  ___ ...
    //     /j=1 2 3  \/j=1 2 3 \/
    //    [  x x x      x x x            ]
    //
    // gradN.MultiplyTT(pderiv,xji);
    for (int inode = 0; inode < nenparent; inode++)  // L     Loop
    {
      for (int idof = 0; idof < my::nsd_; idof++)  // j     Loop
      {
        for (int idof2 = 0; idof2 < my::nsd_; idof2++)  // alpha Loop
        {
          gradN(0, (inode * my::nsd_) + idof) += pderiv_loc(idof2, inode) * (xji(idof, idof2));
        }
        gradNon(0, inode) += gradN(0, inode * my::nsd_ + idof) * my::unitnormal_(idof);
      }
    }


    // gradient of u once contracted with linearization of normal
    //
    //                                L= 1 ... nenparent
    //                         i=   _ l= 1 ... my::nsd_        _
    //               d  n_j      1 |     ...                |
    //   N_A,j u^A_i -------- =  2 |     ...                |
    //               d d^L_l     3 |_    ...               _|
    //
    LINALG::Matrix<my::nsd_, my::nsd_ * nenparent> graduonormalderiv;
    graduonormalderiv.Multiply(dudxioJinv, normalderiv);

    // transposed gradient of u once contracted with linearization of normal
    //
    //                                L= 1 ... nenparent
    //                         i=   _ l= 1 ... my::nsd_        _
    //               d  n_j      1 |     ...                |
    //   N_A,i u^A_j -------- =  2 |     ...                |
    //               d d^L_l     3 |_    ...               _|
    //
    LINALG::Matrix<my::nsd_, my::nsd_ * nenparent> graduTonormalderiv;
    graduTonormalderiv.MultiplyTN(dudxioJinv, normalderiv);

    // Isn't that cool?
    LINALG::Matrix<1, nenparent> survivor;
    for (int inode = 0; inode < nenparent; inode++)
    {
      if (pfunct(inode) != 0)
      {
        survivor(0, inode) = 1.0;
      }
      else
      {
        survivor(0, inode) = 0.;
      }
    }

#ifdef DEBUG
    if (abs(scalarintegraltransformfac - my::drs_) > 1e-11)
    {
      std::cout << "my::drs_ = " << my::drs_ << std::endl;
      std::cout << "scalarintegraltransformfac = " << scalarintegraltransformfac << std::endl;
      dserror("scalarintegraltransformfac should be equal my::drs_ !");
    }
#endif

    normalvelocity.MultiplyTN(my::velint_, my::unitnormal_);

    ////////////////////////////////////////////////////////////////////////////
    //////////////////////////      Loop over Nodes       //////////////////////
    ////////////////////////////////////////////////////////////////////////////
    for (int inode = 0; inode < nenparent; inode++)
    {
      double normal_u_minus_vs = 0.0;
      LINALG::Matrix<1, my::nsd_> u_minus_vs(true);

      for (int idof = 0; idof < my::nsd_; idof++)
      {
        normal_u_minus_vs += my::unitnormal_(idof) * (my::velint_(idof) - gridvelint(idof));
        u_minus_vs(idof) = my::velint_(idof) - gridvelint(idof);
      }

      LINALG::Matrix<1, nenparent * my::nsd_> u_minus_vs_normalderiv(true);
      u_minus_vs_normalderiv.Multiply(u_minus_vs, normalderiv);


      ////////////////////////////////////////////////////////////////////////////
      ////////////////////////      Fill Element Matrix      /////////////////////
      ////////////////////////////////////////////////////////////////////////////
      for (int nnod = 0; nnod < nenparent; nnod++)
      {
        for (int idof2 = 0; idof2 < my::nsd_; idof2++)
        {
          if (block == "Porofluid_Freefluid")
          {
            /*
                    d(q,(u-vs) o n) / d(u)

                    evaluated on FluidField(): flip sign because my::unitnormal_ points in opposite
               direction
             */
            // double timefacfacpre =
            // DRT::ELEMENTS::FluidEleParameter::Instance(INPAR::FPSI::porofluid)->TimeFacPre();
            elemat1(inode * my::numdofpernode_ + my::nsd_, nnod * my::numdofpernode_ + idof2) -=
                ((timefacfacpre)*pfunct(inode) * my::unitnormal_(idof2) * pfunct(nnod));
          }
          else if (block == "Porofluid_Structure")
          {
            /*
                      d(q,(u-vs) o n) / d(ds)

                      evaluated on FluidField(): my::unitnormal_ points in wrong direction -> flip
               sign
             */
            // double timefacpre =
            // DRT::ELEMENTS::FluidEleParameter::Instance(INPAR::FPSI::porofluid)->TimeFacPre();
            elemat1(inode * my::numdofpernode_ + my::nsd_, nnod * my::numdofpernode_ + idof2) +=
                -u_minus_vs_normalderiv(0, nnod * my::nsd_ + idof2) * pfunct(inode) * timefacpre *
                    fac *
                    survivor(
                        nnod)  // no my::drs_ needed, since it is contained in the linearization
                               // w.r.t. nonunitnormal (normalderiv) -> timefacpre*fac instead of
                               // timefafacpre = timefacpre * my::fac_ (my::fac_ = fac*my::drs_)
                +
                pfunct(inode) * my::unitnormal_(idof2) * timescale * pfunct(nnod) * (timefacfacpre);
          }

          else if (block == "Fluid_Porofluid")
          {
            /*
                        d(w o n, pf_pm) / d(pf_pm) (3)

                        evaluated on PoroField(): flip sign because my::unitnormal_ points in
               opposite direction
             */
            elemat1((inode * my::numdofpernode_) + idof2, (nnod * my::numdofpernode_) + my::nsd_) -=
                (  // sign checked to be negative
                    pfunct(inode) * pfunct(nnod) * my::unitnormal_(idof2)

                        ) *
                my::fac_ * timefac;  // scalarintegraltransformfac;


            /*                              _                      _
                              I  alpha mu_f  |                        |   I  /
                        d(w o t,------------ | u - (vs + phi(vf -vs)) | o t / d(pfpm)
                                  rho_f K    |_           |          _|    /
                                 \_________/              V
                                tangentialfac         porosityint

                                evaluated on PoroField(): no sign flipping because there's no
               multiplication by my::unitnormal_

             */
            elemat1((inode * my::numdofpernode_) + idof2, (nnod * my::numdofpernode_) + my::nsd_) -=
                (  // sign checked to be negative
                    tangential1(idof2, 0) * (tangentialvf1 - tangentialvs1) +  // d phi / dpfpm
                    tangential2(idof2, 0) * (tangentialvf2 - tangentialvs2)

                        ) *
                pfunct(inode) * tangentialfac * dphi_dp * my::fac_ *
                timefac;  // scalarintegraltransformfac;

            for (int idof3 = 0; idof3 < my::nsd_; idof3++)
            {
              /*                              _                      _
                                I  alpha mu_f  |                        |   I  /
                          d(w o t,------------ | u - (vs + phi(vf -vs)) | o t / d(vf)
                                    rho_f K    |_           |          _|    /
                                   \_________/              V
                                  tangentialfac         porosityint

                                  evaluated on PoroField(): no sign flipping because there's no
                 multiplication by my::unitnormal_

               */
              elemat1((inode * my::numdofpernode_) + idof2, (nnod * my::numdofpernode_) + idof3) -=
                  (  // sign checked to be negative
                      tangential1(idof2, 0) * tangential1(idof3, 0) +
                      tangential2(idof2, 0) * tangential2(idof3, 0)

                          ) *
                  pfunct(inode) * pfunct(nnod) * porosityint * tangentialfac * my::fac_ * timefac;
            }
          }

          else if (block == "Fluid_Structure")
          {
            if (discretization.Name() == "porofluid")
            {
              /*
                        d(w o n, pf_pm * my::drs_) / d(ds)

                        evaluated on PoroField(): flip sign because my::unitnormal_ points in
                 opposite direction
               */
              for (int idof3 = 0; idof3 < my::nsd_; idof3++)
              {
                elemat1((inode * my::numdofpernode_) + idof2, (nnod * my::nsd_) + idof3) -=
                    (pfunct(inode) * normalderiv(idof2, (nnod * my::nsd_) + idof3)) *
                    pressint(0, 0) * fac * timefac *
                    survivor(
                        nnod);  // *my::fac_ since normalderiv is referring to the test function
              }                 // idof3

              /*                              _                      _
                              I  alpha mu_f  |                        |   I  /
                        d(w o t,------------ | u - (vs + phi(vf -vs)) | o t / d(ds)
                                  rho_f K    |_           |          _|    /
                                 \_________/              V
                                tangentialfac         porosityint

                                evaluated on PoroField():
               */
              for (int idof3 = 0; idof3 < my::nsd_; idof3++)
              {
                elemat1((inode * my::numdofpernode_) + idof2, (nnod * my::nsd_) + idof3) -=
                    ((tangential1(idof2, 0) *
                             (tangentialvs1 +
                                 porosityint *
                                     (tangentialvf1 - tangentialvs1)) +  // d ||n||/d d^L_l
                         tangential2(idof2, 0) *
                             (tangentialvs2 + porosityint * (tangentialvf2 - tangentialvs2))

                             ) *
                            (linearizationofscalarintegraltransformfac(idof3, nnod) / my::drs_) *
                            survivor(nnod)  // -> survivor(nnod) in order to filter the entries
                                            // which do not belong to the interface
                        +
                        (tangentialderiv1(idof2, (nnod * my::nsd_) + idof3) *
                                (porosityint * (tangentialvf1 - tangentialvs1)) +  // d t^i/d d^L_l
                            tangentialderiv2(idof2, (nnod * my::nsd_) + idof3) *
                                (porosityint * (tangentialvf2 - tangentialvs2))

                                ) *
                            survivor(nnod) +
                        (tangential1(idof2, 0) * (vfotangentialderiv1((nnod * my::nsd_) + idof3) -
                                                     vsotangentialderiv1((nnod * my::nsd_) +
                                                                         idof3)) +  // d t^j/d d^L_l
                            tangential2(idof2, 0) *
                                (vfotangentialderiv2((nnod * my::nsd_) + idof3) -
                                    vsotangentialderiv2((nnod * my::nsd_) + idof3))

                                ) *
                            porosityint * survivor(nnod) -
                        (tangential1(idof2, 0) *
                                tangential1(idof3, 0) +  // d vs / d d^L_l  (sign checked)
                            tangential2(idof2, 0) * tangential2(idof3, 0)

                                ) *
                            pfunct(nnod) * timescale * porosityint +
                        (tangential1(idof2, 0) *
                                (tangentialvf1 - tangentialvs1) +  // d phi / d d^L_l
                            tangential2(idof2, 0) * (tangentialvf2 - tangentialvs2)

                                ) *
                            dphi_dJ * dJ_dds((nnod * my::nsd_) + idof3) +
                        (tangential1(idof2, 0) *
                                tangential1(idof3,
                                    0) +  // d vs / d d^L_l (front term without phi) (sign checked)
                            tangential2(idof2, 0) * tangential2(idof3, 0)

                                ) *
                            pfunct(nnod) * timescale +
                        (tangentialderiv1(idof2, (nnod * my::nsd_) + idof3) *
                                tangentialvs1 +  // d t^i/d d^L_l (front term without phi)
                            tangentialderiv2(idof2, (nnod * my::nsd_) + idof3) * tangentialvs2

                            ) *
                            survivor(nnod) +
                        (tangential1(idof2, 0) *
                                vsotangentialderiv1(
                                    (nnod * my::nsd_) +
                                    idof3) +  // d t^j/d d^L_l (front term without phi)
                            tangential2(idof2, 0) * vsotangentialderiv2((nnod * my::nsd_) + idof3)

                                ) *
                            survivor(nnod)

                            ) *
                    pfunct(inode) * tangentialfac * my::fac_ * timefac;

                if (probtype == prb_fps3i)
                {
                  /*
                            d(w o n,(u-vs) o n) / d(ds)

                            evaluated on PoroField(): sign flip*/

                  elemat1((inode * my::numdofpernode_) + idof2, (nnod * my::nsd_) + idof3) +=
                      ((-u_minus_vs_normalderiv(0, nnod * my::nsd_ + idof2) * pfunct(inode) *
                               my::fac_ * timefac * survivor(nnod) +
                           pfunct(inode) * my::unitnormal_(idof2) * timescale * pfunct(nnod) *
                               my::fac_ * timefac) /
                          (Lp));
                }
              }  // idof3
            }
            else if (discretization.Name() == "fluid")
            {
              for (int idof3 = 0; idof3 < my::nsd_; idof3++)
              {
                elemat1(
                    (inode * my::numdofpernode_) + idof2, (nnod * my::numdofpernode_) + idof3) +=

                    ((tangential1(idof2, 0) * tangentialvf1 +  // d ||n||/d d^L_l
                         tangential2(idof2, 0) * tangentialvf2

                         ) * (linearizationofscalarintegraltransformfac(idof3, nnod) / my::drs_) *
                            survivor(nnod)  // -> survivor(nnod) in order to filter the entries
                                            // which do not belong to the interface
                        + (tangentialderiv1(idof2, (nnod * my::nsd_) + idof3) *
                                  tangentialvf1 +  // d t^i/d d^L_l
                              tangentialderiv2(idof2, (nnod * my::nsd_) + idof3) * tangentialvf2

                              ) *
                              survivor(nnod) +
                        (tangential1(idof2, 0) *
                                vfotangentialderiv1((nnod * my::nsd_) + idof3) +  // d t^j/d d^L_l
                            tangential2(idof2, 0) * vfotangentialderiv2((nnod * my::nsd_) + idof3)

                                ) *
                            survivor(nnod)) *
                    my::fac_ * timefac * pfunct(inode) * tangentialfac;
              }
            }
          }  // block Fluid_Structure

          else if (block == "Fluid_Fluid")
          {
            /*
                        d(w o t,tangentialfac * u o t) / d(du)
             */
            for (int idof3 = 0; idof3 < my::nsd_; idof3++)
            {
              elemat1((inode * my::numdofpernode_) + idof2, (nnod * my::numdofpernode_) + idof3) +=
                  (tangential1(idof2) * tangential1(idof3) +
                      tangential2(idof2) * tangential2(idof3)) *
                  pfunct(nnod) * pfunct(inode) * tangentialfac * my::fac_ * timefac;
              if (probtype == prb_fps3i)
              {
                /*
                     d(w o n,(u-vs) o n) / d(u)

                     evaluated on FluidField(): no sign flip */
                elemat1((inode * my::numdofpernode_) + idof2,
                    (nnod * my::numdofpernode_) + idof3) -= my::fac_ * timefac * pfunct(inode) *
                                                            my::unitnormal_(idof2) * pfunct(nnod) /
                                                            (Lp);
              }
            }
          }  // block Fluid_Fluid

          else if (block == "NeumannIntegration" and elemat1 != Teuchos::null)
          {
            if (discretization.Name() == "fluid")
            {
              /*
                        d (d,[tau - pf o I + gamma rho_f u dyadic u] o [x,1 x x,2]) / d(du)
                               |
                               V
                       2*mu*0.5*(u_i,j+u_j,i)

                       evaluated on FluidField()
               */
              elemat1((inode * my::numdofpernode_) + idof2, (nnod * my::numdofpernode_) + idof2) -=
                  (
                      // d (mu*(u_i,j+u_j,i)) / d u^L_l
                      pfunct(inode) * gradNon(0, nnod)  // d u_i,j / d u^L_l
                      ) *
                  fluiddynamicviscosity * my::fac_ * timefac;

              elemat1(
                  (inode * my::numdofpernode_) + idof2, (nnod * my::numdofpernode_) + my::nsd_) +=
                  (
                      // d (dd , pf o n) / d pf_B
                      // flip sign
                      pfunct(inode) * pfunct(nnod) * my::unitnormal_(idof2)) *
                  my::fac_ * timefac;

              for (int idof3 = 0; idof3 < my::nsd_; idof3++)
              {
                elemat1(
                    (inode * my::numdofpernode_) + idof2, (nnod * my::numdofpernode_) + idof3) -=
                    (
                        // d (2*mu*0.5*(u_i,j+u_j,i)) / d u^L_l
                        pfunct(inode) * gradN(0, (nnod * my::nsd_) + idof2) *
                        my::unitnormal_(idof3) * fluiddynamicviscosity  // d u_j,i / d u^L_l
                        ) *
                    my::fac_ * timefac;
              }
            }  // if dis=fluid
          }    // block NeumannIntegration

          else if (block == "NeumannIntegration_Ale")
          {
            for (int idof3 = 0; idof3 < my::nsd_; idof3++)
            {
              elemat1((inode * my::numdofpernode_) + idof2, (nnod * my::nsd_) + idof3) -=
                  (
                      // d (dd , - pf o n) / d d^L_l
                      -pfunct(inode) * pressint(0, 0) *
                          normalderiv(idof2, (nnod * my::nsd_) + idof3) * fac  // d n_j / d d^L_l

                      // d (dd, mu*u_i,j o n ) / d d^L_l
                      - fluiddynamicviscosity * pfunct(inode) * dudxioJinv(idof2, idof3) *
                            dNdxon(nnod) * my::fac_  // d ui,j / d d^L_l
                      + fluiddynamicviscosity * pfunct(inode) *
                            graduonormalderiv(idof2, (nnod * my::nsd_) + idof3) *
                            fac  // d n / d d^L_l

                      // d (dd, mu*u_j,i o n ) / d d^L_l
                      - fluiddynamicviscosity * pfunct(inode) * graduTon(0, idof3) *
                            derxy(idof2, nnod) * my::fac_  // d uj,i / d d^L,l
                      + fluiddynamicviscosity * pfunct(inode) *
                            graduTonormalderiv(idof2, (nnod * my::nsd_) + idof3) *
                            fac  // d n_j / d^L_l
                      ) *
                  timefac;  // split afterwards, as this is assembled into a blockmatrix
            }
          }  // block == "NeumannIntegration_Ale"
          else if (block == "Structure_Fluid")
          {
            /*
                                      d (d,[tau - pf o I + gamma rho_f u dyadic u] o [x,1 x x,2]) /
               d(du)
                                             |
                                             V
                                     2*mu*0.5*(u_i,j+u_j,i)

                                     evaluated on FluidField()
             */
            elemat1((inode * my::numdofpernode_) + idof2, (nnod * my::numdofpernode_) + idof2) +=
                ((
                     // d (mu*(u_i,j+u_j,i)) / d u^L_l

                     pfunct(inode) * gradNon(0, nnod)  // d u_i,j / d u^L_l

                     ) *
                    fluiddynamicviscosity * my::fac_ * theta);


            elemat1((inode * my::numdofpernode_) + idof2, (nnod * my::numdofpernode_) + my::nsd_) -=
                ((
                     // d (dd , pf o n) / d pf_B
                     // flip sign

                     pfunct(inode) * pfunct(nnod) * my::unitnormal_(idof2)

                         ) *
                    my::fac_ * theta);

            for (int idof3 = 0; idof3 < my::nsd_; idof3++)
            {
              elemat1((inode * my::numdofpernode_) + idof2, (nnod * my::numdofpernode_) + idof3) +=
                  (
                      // d (2*mu*0.5*(u_i,j+u_j,i)) / d u^L_l

                      pfunct(inode) * gradN(0, (nnod * my::nsd_) + idof2) *
                      my::unitnormal_(idof3)  // d u_j,i / d u^L_l
                      ) *
                  my::fac_ * theta * fluiddynamicviscosity;
            }
          }  // block structure_fluid

          else if (block == "Structure_Structure")
          {
            for (int idof3 = 0; idof3 < my::nsd_; idof3++)
            {
              elemat1((inode * my::numdofpernode_) + idof2, (nnod * my::numdofpernode_) + idof3) +=
                  (
                      // d (dd , - pf o n) / d d^L_l
                      -pfunct(inode) * pressint(0, 0) *
                          normalderiv(idof2, (nnod * my::nsd_) + idof3) * fac  // d n_j / d d^L_l

                      // d (dd, mu*u_i,j o n ) / d d^L_l
                      - fluiddynamicviscosity * pfunct(inode) * dudxioJinv(idof2, idof3) *
                            dNdxon(nnod) * my::fac_  // d ui,j / d d^L_l
                      + fluiddynamicviscosity * pfunct(inode) *
                            graduonormalderiv(idof2, (nnod * my::nsd_) + idof3) *
                            fac  // d n / d d^L_l

                      // d (dd, mu*u_j,i o n ) / d d^L_l
                      - fluiddynamicviscosity * pfunct(inode) * graduTon(0, idof3) *
                            derxy(idof2, nnod) * my::fac_  // d uj,i / d d^L,l
                      + fluiddynamicviscosity * pfunct(inode) *
                            graduTonormalderiv(idof2, (nnod * my::nsd_) + idof3) *
                            fac  // d n_j / d^L_l
                      ) *
                      survivor(nnod) * theta
                  // linearisation of the old timestep --> change of my::fac_
                  + (linearizationofscalarintegraltransformfac(idof3, nnod) * fac *
                        (pfunct(inode) *
                            (fluiddynamicviscosity * (graduon_n(idof2) + graduTon_n(idof2)) -
                                pressint_n(0, 0) * my::unitnormal_n_(idof2)))  // d (...)^n/d^L_l
                        ) *
                        survivor(nnod) * (1 - theta);  // <- only boundary dofs survive
            }
          }  // block structure_structure
          else if (block == "Structure_Ale")
          {
            for (int idof3 = 0; idof3 < my::nsd_; idof3++)
            {
              elemat1((inode * my::numdofpernode_) + idof2, (nnod * my::nsd_) + idof3) +=
                  (
                      // d (dd, mu*u_i,j o n ) / d d^L_l
                      -fluiddynamicviscosity * pfunct(inode) * dudxioJinv(idof2, idof3) *
                          dNdxon(nnod) * my::fac_  // d ui,j / d d^L_l

                      // d (dd, mu*u_j,i o n ) / d d^L_l
                      - fluiddynamicviscosity * pfunct(inode) * graduTon(0, idof3) *
                            derxy(idof2, nnod) * my::fac_  // d uj,i / d d^L,l

                      ) *
                  abs(survivor(0, nnod) - 1.0) * theta;  // <- only inner dofs survive
            }
          }  // block structure_ale

          else if (block == "defaultblock" && (block != "fluid" && block != "fluidfluid" &&
                                                  block != "structure" && block != "conti"))
          {
            dserror("no proper block specification available in parameterlist ...");
          }  // blocks
        }    // idof2
      }      // nnod
    }        // Loop over parent nodes (inode)

    tangentialvelocity1.MultiplyTN(my::velint_, tangential1);
    tangentialvelocity2.MultiplyTN(my::velint_, tangential2);
    tangentialgridvelocity1.MultiplyTN(gridvelint, tangential1);
    tangentialgridvelocity2.MultiplyTN(gridvelint, tangential2);


    ////////////////////////////////////////////////////////////////////////////
    //////////////////////////      Loop over Nodes       //////////////////////
    ////////////////////////////////////////////////////////////////////////////
    for (int inode = 0; inode < nenparent; inode++)
    {
      double normal_u_minus_vs = 0.0;
      LINALG::Matrix<1, my::nsd_> u_minus_vs(true);

      for (int idof = 0; idof < my::nsd_; idof++)
      {
        normal_u_minus_vs += my::unitnormal_(idof) * (my::velint_(idof) - gridvelint(idof));
        u_minus_vs(idof) = my::velint_(idof) - gridvelint(idof);
      }

      LINALG::Matrix<1, nenparent * my::nsd_> u_minus_vs_normalderiv(true);
      u_minus_vs_normalderiv.Multiply(u_minus_vs, normalderiv);

      ////////////////////////////////////////////////////////////////////////////
      ////////////////////////            Fill RHS           /////////////////////
      ////////////////////////////////////////////////////////////////////////////

      if (block == "conti")
      {
        /*
            Evaluated on FluidField() wears (+) in residual; multiplied by (-1) for RHS; switch sign
           because of opposite normal -> (+)
         */
        // double rhsfacpre =
        // DRT::ELEMENTS::FluidEleParameter::Instance(INPAR::FPSI::porofluid)->TimeFacRhsPre();
        elevec1(inode * my::numdofpernode_ + my::nsd_) +=
            rhsfac * pfunct(inode) * normal_u_minus_vs;

      }  // block conti

      else if (block == "structure")
      {
        /*
                    (2)  N * (tau - pf I) o n   << from last iteration at time n+1

                    evaluated on FluidField(); my::unitnormal_ opposite to strucutral unitnormal ->
           application of nanson's formula yields structural normal -> * (-1)
         */
        for (int idof2 = 0; idof2 < my::nsd_; idof2++)
        {
          elevec1(inode * my::numdofpernode_ + idof2) -=
              (theta * pfunct(inode) *
                      (fluiddynamicviscosity * (graduon(idof2) + graduTon(idof2)) -
                          pressint(0, 0) * my::unitnormal_(idof2)) +
                  (1.0 - theta) * pfunct(inode) *
                      (fluiddynamicviscosity * (graduon_n(idof2) + graduTon_n(idof2)) -
                          pressint_n(0, 0) * my::unitnormal_n_(idof2))) *
              survivor(inode) * my::fac_;
        }
      }  // block structure

      else if (block == "fluid")  // rhs of fluid evaluated on porofluidfield
      {
        /*
                  evaluated on PoroFluidField()

                  (3+4) - N*n * 1/rhof * (pf) + N*t*tangentialfac*[u- (vs + phi(vf-vs))]ot  << from
           last iteration at time n+1
         */
        for (int idof2 = 0; idof2 < my::nsd_; idof2++)
        {
          elevec1(inode * my::numdofpernode_ + idof2) +=
              (+(pfunct(inode) * my::unitnormal_(idof2) *
                   pressint(0, 0))  // pressure part // pressure part
                  + ((pfunct(inode) * tangential1(idof2) *
                         (tangentialgridvelocity1(0, 0) +
                             porosityint * (tangentialvelocity1(0, 0) -
                                               tangentialgridvelocity1(0, 0))))  // Beavers-Joseph
                        + (pfunct(inode) * tangential2(idof2) *
                              (tangentialgridvelocity2(0, 0) +
                                  porosityint * (tangentialvelocity2(0, 0) -
                                                    tangentialgridvelocity2(0, 0))))) *
                        tangentialfac) *
              rhsfac * survivor(inode);
        }
      }  // block fluid

      else if (block == "fluidfluid")  // rhs of fluid evaluated on fluidfield
      {
        /*
                    (4)  N*t*tangentialfac*[u]ot  << from last iteration at time n+1
         */
        for (int idof2 = 0; idof2 < my::nsd_; idof2++)
        {
          elevec1(inode * my::numdofpernode_ + idof2) -=
              (pfunct(inode) * tangential1(idof2) * tangentialvelocity1(0, 0) +
                  pfunct(inode) * tangential2(idof2) * tangentialvelocity2(0, 0)) *
              tangentialfac * rhsfac * survivor(inode);

          // in case of FPS3I we have to add the first Kedem-Katchalsky equation to prescribe the
          // volume flux see e.g. Kedem, O. T., and A. Katchalsky. "Thermodynamic analysis of the
          // permeability of biological membranes to non-electrolytes." Biochimica et biophysica
          // Acta 27 (1958): 229-246.
          // TODO: One could think of not using this equation, i.e. having L_p -> inf (Thon)
          if (probtype == prb_fps3i)
          {
            // evaluated on fluid field --> no sign flip
            elevec1(inode * my::numdofpernode_ + idof2) -=
                rhsfac * survivor(inode) * pfunct(inode) * normal_u_minus_vs / (Lp);
          }
        }
      }  // block fluidfluid

      else if (block == "NeumannIntegration")
      {
        if (discretization.Name() != "fluid")
        {
          dserror(
              "Tried to call NeumannIntegration on a discretization other than 'fluid'. \n"
              "You think that's funny, hu ?? Roundhouse-Kick !!!");
        }

        for (int idof2 = 0; idof2 < my::nsd_; idof2++)
        {
          elevec1(inode * my::numdofpernode_ + idof2) +=
              ((-pfunct(inode) * pressint(0, 0) * my::unitnormal_(idof2) * rhsfac) +
                  pfunct(inode) * fluiddynamicviscosity * (graduon(idof2) + graduTon(idof2)) *
                      rhsfac);
        }  // block NeumannIntegration

      }  // NeumannIntegration
    }    // Loop over interface nodes (inode)
  }      // Loop over integration points
  return;
}  // FPSI Coupling Terms

/*----------------------------------------------------------------------*
 * compute flow rate over boundary (protected)              vuong 06/11||
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::ComputeFlowRate(
    DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& plm, Epetra_SerialDenseVector& elevec1)
{
  switch (distype)
  {
    // 2D:
    case DRT::Element::line2:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::quad4)
      {
        this->ComputeFlowRate<DRT::Element::quad4>(ele, params, discretization, plm, elevec1);
      }
      else
      {
        dserror("expected combination line2/quad4 for line/parent pair");
      }
      break;
    }
    case DRT::Element::line3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::quad9)
      {
        this->ComputeFlowRate<DRT::Element::quad9>(ele, params, discretization, plm, elevec1);
      }
      else
      {
        dserror("expected combination line3/quad9 for line/parent pair");
      }
      break;
    }
    case DRT::Element::nurbs3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::nurbs9)
      {
        this->ComputeFlowRate<DRT::Element::nurbs9>(ele, params, discretization, plm, elevec1);
      }
      else
      {
        dserror("expected combination nurbs3/nurbs9 for line/parent pair");
      }
      break;
    }
    // 3D:
    case DRT::Element::quad4:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::hex8)
      {
        this->ComputeFlowRate<DRT::Element::hex8>(ele, params, discretization, plm, elevec1);
      }
      else
      {
        dserror("expected combination quad4/hex8 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::tri3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::tet4)
      {
        this->ComputeFlowRate<DRT::Element::tet4>(ele, params, discretization, plm, elevec1);
      }
      else
      {
        dserror("expected combination tri3/tet4 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::tri6:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::tet10)
      {
        this->ComputeFlowRate<DRT::Element::tet10>(ele, params, discretization, plm, elevec1);
      }
      else
      {
        dserror("expected combination tri6/tet10 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::quad9:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::hex27)
      {
        this->ComputeFlowRate<DRT::Element::hex27>(ele, params, discretization, plm, elevec1);
      }
      else
      {
        dserror("expected combination hex27/hex27 for surface/parent pair");
      }
      break;
    }
    default:
    {
      dserror("surface/parent element pair not yet implemented. Just do it.\n");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 * compute flow rate over boundary (protected)              vuong 06/11||
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType pdistype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::ComputeFlowRate(
    DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& plm, Epetra_SerialDenseVector& elevec1)
{
  // This function is only implemented for 3D and 2D
  if (my::bdrynsd_ != 2 and my::bdrynsd_ != 1)
    dserror("PoroBoundary is only implemented for 3D and 2D!");

  // get element location vector and ownerships
  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;
  ele->DRT::Element::LocationVector(discretization, lm, lmowner, lmstride);

  /// number of parentnodes
  static const int nenparent = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

  // get the parent element
  DRT::ELEMENTS::Fluid* pele = ele->ParentElement();

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<my::bdrynsd_> intpoints(
      DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a my::nsd_ dimensional domain, since my::nsd_ determines the dimension of
  // FluidBoundary element!)
  GEO::fillInitialPositionArray<distype, my::nsd_, LINALG::Matrix<my::nsd_, my::bdrynen_>>(
      ele, my::xyze_);

  // displacements
  Teuchos::RCP<const Epetra_Vector> dispnp;
  std::vector<double> mydispnp;
  std::vector<double> parentdispnp;

  dispnp = discretization.GetState("dispnp");
  if (dispnp != Teuchos::null)
  {
    mydispnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp, mydispnp, lm);
    DRT::UTILS::ExtractMyValues(*dispnp, parentdispnp, plm);
  }
  dsassert(mydispnp.size() != 0, "no displacement values for boundary element");
  dsassert(parentdispnp.size() != 0, "no displacement values for parent element");

  // Add the deformation of the ALE mesh to the nodes coordinates
  for (int inode = 0; inode < my::bdrynen_; ++inode)
    for (int idim = 0; idim < my::nsd_; ++idim)
      my::xyze_(idim, inode) += mydispnp[my::numdofpernode_ * inode + idim];

  // update element geometry of parent element
  LINALG::Matrix<my::nsd_, nenparent> xrefe;  // material coord. of parent element
  LINALG::Matrix<my::nsd_, nenparent> xcurr;  // current  coord. of parent element
  {
    DRT::Node** nodes = pele->Nodes();
    for (int i = 0; i < nenparent; ++i)
    {
      const double* x = nodes[i]->X();
      for (int j = 0; j < my::nsd_; ++j)
      {
        xrefe(j, i) = x[j];
        xcurr(j, i) = xrefe(j, i) + parentdispnp[i * my::numdofpernode_ + j];
      }
    }
  }

  // extract local values from the global vectors
  // renamed to "velaf" to be consistent in fluidimplicitintegration.cpp (krank 12/13)
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velaf");
  Teuchos::RCP<const Epetra_Vector> gridvel = discretization.GetState("gridv");

  if (velnp == Teuchos::null) dserror("Cannot get state vector 'velaf'");
  if (gridvel == Teuchos::null) dserror("Cannot get state vector 'gridv'");

  std::vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp, myvelnp, lm);
  std::vector<double> mygridvel(lm.size());
  DRT::UTILS::ExtractMyValues(*gridvel, mygridvel, lm);

  // allocate velocity vectors
  LINALG::Matrix<my::nsd_, my::bdrynen_> evelnp(true);
  LINALG::Matrix<my::bdrynen_, 1> epressnp(true);
  LINALG::Matrix<my::nsd_, my::bdrynen_> edispnp(true);
  LINALG::Matrix<my::nsd_, my::bdrynen_> egridvel(true);
  LINALG::Matrix<my::bdrynen_, 1> escaaf(true);
  LINALG::Matrix<my::bdrynen_, 1> eporosity(true);

  // split velocity and pressure, insert into element arrays
  for (int inode = 0; inode < my::bdrynen_; inode++)
  {
    for (int idim = 0; idim < my::nsd_; idim++)
    {
      evelnp(idim, inode) = myvelnp[idim + (inode * my::numdofpernode_)];
      edispnp(idim, inode) = mydispnp[idim + (inode * my::numdofpernode_)];
      egridvel(idim, inode) = mygridvel[idim + (inode * my::numdofpernode_)];
    }
    epressnp(inode) = myvelnp[my::nsd_ + (inode * my::numdofpernode_)];
  }

  ComputeNodalPorosity(ele, mydispnp, eporosity);

  // get coordinates of gauss points w.r.t. local parent coordinate system
  LINALG::SerialDenseMatrix pqxg(intpoints.IP().nquad, my::nsd_);
  LINALG::Matrix<my::nsd_, my::nsd_> derivtrafo(true);

  DRT::UTILS::BoundaryGPToParentGP<my::nsd_>(
      pqxg, derivtrafo, intpoints, pdistype, distype, ele->SurfaceNumber());


  // --------------------------------------------------
  // Now do the nurbs specific stuff
  // --------------------------------------------------

  // In the case of nurbs the normal vector is multiplied with normalfac
  double normalfac = 0.0;
  std::vector<Epetra_SerialDenseVector> mypknots(my::nsd_);
  std::vector<Epetra_SerialDenseVector> myknots(my::bdrynsd_);
  Epetra_SerialDenseVector weights(my::bdrynen_);
  Epetra_SerialDenseVector pweights(pele->NumNode());

  // for isogeometric elements --- get knotvectors for parent
  // element and surface element, get weights
  if (IsNurbs<distype>::isnurbs)
  {
    bool zero_size = DRT::NURBS::GetKnotVectorAndWeightsForNurbsBoundaryAndParent(pele, ele,
        ele->SurfaceNumber(), discretization, mypknots, myknots, pweights, weights, normalfac);

    if (zero_size)
    {
      return;
    }
  }
  // --------------------------------------------------

  // structure velocity at gausspoint
  LINALG::Matrix<my::nsd_, 1> gridvelint;

  // coordinates of gauss points of parent element
  LINALG::Matrix<my::nsd_, 1> pxsi(true);

  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    // get shape functions and derivatives in the plane of the element
    LINALG::Matrix<nenparent, 1> pfunct(true);
    LINALG::Matrix<my::nsd_, nenparent> pderiv_loc;

    // coordinates of the current integration point
    for (int idim = 0; idim < my::nsd_; idim++) pxsi(idim) = pqxg(gpid, idim);

    // get shape functions and derivatives of the parent element
    if (not IsNurbs<distype>::isnurbs)
    {
      // shape functions and their first derivatives of parent element
      DRT::UTILS::shape_function<pdistype>(pxsi, pfunct);
      DRT::UTILS::shape_function_deriv1<pdistype>(pxsi, pderiv_loc);
    }
    // only for NURBS!!!
    else
    {
      DRT::NURBS::UTILS::nurbs_get_funct_deriv(
          pfunct, pderiv_loc, pxsi, mypknots, pweights, pdistype);
    }

    // get Jacobian matrix and determinant w.r.t. spatial configuration
    // transposed jacobian "dx/ds"
    LINALG::Matrix<my::nsd_, my::nsd_> xjm;
    LINALG::Matrix<my::nsd_, my::nsd_> Jmat;
    xjm.MultiplyNT(pderiv_loc, xcurr);
    Jmat.MultiplyNT(pderiv_loc, xrefe);
    // jacobian determinant "det(dx/ds)"
    const double det = xjm.Determinant();
    // jacobian determinant "det(dX/ds)"
    const double detJ = Jmat.Determinant();
    // jacobian determinant "det(dx/dX) = det(dx/ds)/det(dX/ds)"
    const double J = det / detJ;

    // Computation of the integration factor & shape function at the Gauss point & derivative of the
    // shape function at the Gauss point Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(my::funct_, my::deriv_, my::fac_,
        my::unitnormal_, my::drs_, my::xsi_, my::xyze_, intpoints, gpid, &myknots, &weights,
        IsNurbs<distype>::isnurbs);

    // in the case of nurbs the normal vector must be scaled with a special factor
    if (IsNurbs<distype>::isnurbs) my::unitnormal_.Scale(normalfac);

    my::velint_.Multiply(evelnp, my::funct_);
    gridvelint.Multiply(egridvel, my::funct_);
    double press = epressnp.Dot(my::funct_);

    // double scalar = escaaf.Dot(my::funct_);

    double dphi_dp = 0.0;
    double dphi_dJ = 0.0;
    double porosity_gp = 0.0;

    // params.set<double>("scalar",scalar);

    ComputePorosityAtGP(
        params, ele, my::funct_, eporosity, press, J, gpid, porosity_gp, dphi_dp, dphi_dJ, false);

    // flowrate = uint o normal
    const double flowrate =
        (my::velint_.Dot(my::unitnormal_) - gridvelint.Dot(my::unitnormal_)) * porosity_gp;

    // store flowrate at first dof of each node
    // use negative value so that inflow is positiv
    for (int inode = 0; inode < my::bdrynen_; ++inode)
    {
      // see "A better consistency for low order stabilized finite element methods"
      // Jansen, Collis, Whiting, Shakib
      //
      // Here the principle is used to bring the flow rate to the outside world!!
      //
      // my::funct_ *  velint * n * fac
      //   |      |________________|
      //   |              |
      //   |         flow rate * fac  -> integral over Gamma
      //   |
      // flow rate is distributed to the single nodes of the element
      // = flow rate per node
      //
      // adding up all nodes (ghost elements are handled by the assembling strategy)
      // -> total flow rate at the desired boundary
      //
      // it can be interpreted as a rhs term
      //
      //  ( v , u o n)
      //               Gamma
      //
      elevec1[inode * my::numdofpernode_] += my::funct_(inode) * my::fac_ * flowrate;

      // alternative way:
      //
      //  velint * n * fac
      // |________________|
      //         |
      //    flow rate * fac  -> integral over Gamma
      //     = flow rate per element
      //
      //  adding up all elements (be aware of ghost elements!!)
      //  -> total flow rate at the desired boundary
      //     (is identical to the total flow rate computed above)
    }
  }
  return;
}  // DRT::ELEMENTS::FluidSurface::ComputeFlowRate

/*----------------------------------------------------------------------*
 * evaluate nodal no penetration condition (protected)       vuong 06/11||
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::NoPenetration(
    DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1)
{
  // This function is only implemented for 3D
  if (my::bdrynsd_ != 2 and my::bdrynsd_ != 1)
    dserror("NoPenetration is only implemented for 3D and 2D!");

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<my::bdrynsd_> intpoints(
      DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a my::nsd_ dimensional domain, since my::nsd_ determines the dimension of
  // FluidBoundary element!)
  GEO::fillInitialPositionArray<distype, my::nsd_, LINALG::Matrix<my::nsd_, my::bdrynen_>>(
      ele, my::xyze_);

  // displacements
  Teuchos::RCP<const Epetra_Vector> dispnp;
  std::vector<double> mydispnp;

  dispnp = discretization.GetState("dispnp");
  if (dispnp != Teuchos::null)
  {
    mydispnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp, mydispnp, lm);
  }
  dsassert(mydispnp.size() != 0, "no displacement values for boundary element");

  // Add the deformation of the ALE mesh to the nodes coordinates
  for (int inode = 0; inode < my::bdrynen_; ++inode)
    for (int idim = 0; idim < my::nsd_; ++idim)
      my::xyze_(idim, inode) += mydispnp[my::numdofpernode_ * inode + idim];

  Teuchos::RCP<const Epetra_Vector> condVector;
  std::vector<double> mycondVector;

  condVector = discretization.GetState("condVector");
  if (condVector == Teuchos::null)
    dserror("could not get state 'condVector'");
  else
  {
    mycondVector.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*condVector, mycondVector, lm);
  }
  dsassert(mycondVector.size() != 0, "no condition IDs values for boundary element");

  // calculate normal
  Epetra_SerialDenseVector normal;
  normal.Size(lm.size());

  // gauss point loop
  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the
    // shape function at the Gauss point Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(my::funct_, my::deriv_, my::fac_,
        my::unitnormal_, my::drs_, my::xsi_, my::xyze_, intpoints, gpid, NULL, NULL,
        IsNurbs<distype>::isnurbs);

    for (int inode = 0; inode < my::bdrynen_; ++inode)
    {
      for (int idim = 0; idim < my::nsd_; ++idim)
        normal(inode * my::numdofpernode_ + idim) +=
            my::unitnormal_(idim) * my::funct_(inode) * my::fac_;
      // pressure dof is set to zero
      normal(inode * my::numdofpernode_ + (my::nsd_)) = 0.0;
    }
  } /* end of loop over integration points gpid */

  LINALG::Matrix<my::numdofpernode_, 1> nodenormal(true);

  // check which matrix is to be filled
  POROELAST::coupltype coupling =
      params.get<POROELAST::coupltype>("coupling", POROELAST::undefined);

  if (coupling == POROELAST::fluidfluid)
  {
    // fill element matrix
    for (int inode = 0; inode < my::bdrynen_; inode++)
    {
      for (int i = 0; i < my::numdofpernode_; i++)
        nodenormal(i) = normal(inode * my::numdofpernode_ + i);
      double norm = nodenormal.Norm2();
      nodenormal.Scale(1 / norm);

      for (int idof = 0; idof < my::numdofpernode_; idof++)
      {
        if (mycondVector[inode * my::numdofpernode_ + idof] != 0.0)
        {
          for (int idof2 = 0; idof2 < my::numdofpernode_; idof2++)
            elemat1(inode * my::numdofpernode_ + idof, inode * my::numdofpernode_ + idof2) +=
                nodenormal(idof2);
        }
      }
    }
  }
  else if (coupling == POROELAST::fluidstructure)
  {
    // extract local values from the global vectors
    Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
    Teuchos::RCP<const Epetra_Vector> gridvel = discretization.GetState("gridv");

    if (velnp == Teuchos::null) dserror("Cannot get state vector 'velnp'");
    if (gridvel == Teuchos::null) dserror("Cannot get state vector 'gridv'");

    std::vector<double> myvelnp(lm.size());
    DRT::UTILS::ExtractMyValues(*velnp, myvelnp, lm);
    std::vector<double> mygridvel(lm.size());
    DRT::UTILS::ExtractMyValues(*gridvel, mygridvel, lm);

    // allocate velocity vectors
    LINALG::Matrix<my::nsd_, my::bdrynen_> evelnp(true);
    LINALG::Matrix<my::nsd_, my::bdrynen_> egridvel(true);

    // split velocity and pressure, insert into element arrays
    for (int inode = 0; inode < my::bdrynen_; inode++)
      for (int idim = 0; idim < my::nsd_; idim++)
      {
        evelnp(idim, inode) = myvelnp[idim + (inode * my::numdofpernode_)];
        egridvel(idim, inode) = mygridvel[idim + (inode * my::numdofpernode_)];
      }

    //  derivatives of surface normals wrt mesh displacements
    LINALG::Matrix<my::nsd_, my::bdrynen_ * my::nsd_> normalderiv(true);

    for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
    {
      // Computation of the integration factor & shape function at the Gauss point & derivative of
      // the shape function at the Gauss point Computation of the unit normal vector at the Gauss
      // points is not activated here Computation of nurb specific stuff is not activated here
      DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(my::funct_, my::deriv_, my::fac_,
          my::unitnormal_, my::drs_, my::xsi_, my::xyze_, intpoints, gpid, NULL, NULL,
          IsNurbs<distype>::isnurbs);

      // dxyzdrs vector -> normal which is not normalized
      LINALG::Matrix<my::bdrynsd_, my::nsd_> dxyzdrs(0.0);
      dxyzdrs.MultiplyNT(my::deriv_, my::xyze_);

      // The integration factor is not multiplied with drs
      // since it is the same as the scaling factor for the unit normal derivatives
      // Therefore it cancels out!!
      const double fac = intpoints.IP().qwgt[gpid];

      if (my::nsd_ == 3)
        for (int node = 0; node < my::bdrynen_; ++node)
        {
          normalderiv(0, my::nsd_ * node) += 0.;
          normalderiv(0, my::nsd_ * node + 1) +=
              (my::deriv_(0, node) * dxyzdrs(1, 2) - my::deriv_(1, node) * dxyzdrs(0, 2)) *
              my::funct_(node) * fac;
          normalderiv(0, my::nsd_ * node + 2) +=
              (my::deriv_(1, node) * dxyzdrs(0, 1) - my::deriv_(0, node) * dxyzdrs(1, 1)) *
              my::funct_(node) * fac;

          normalderiv(1, my::nsd_ * node) +=
              (my::deriv_(1, node) * dxyzdrs(0, 2) - my::deriv_(0, node) * dxyzdrs(1, 2)) *
              my::funct_(node) * fac;
          normalderiv(1, my::nsd_ * node + 1) += 0.;
          normalderiv(1, my::nsd_ * node + 2) +=
              (my::deriv_(0, node) * dxyzdrs(1, 0) - my::deriv_(1, node) * dxyzdrs(0, 0)) *
              my::funct_(node) * fac;

          normalderiv(2, my::nsd_ * node) +=
              (my::deriv_(0, node) * dxyzdrs(1, 1) - my::deriv_(1, node) * dxyzdrs(0, 1)) *
              my::funct_(node) * fac;
          normalderiv(2, my::nsd_ * node + 1) +=
              (my::deriv_(1, node) * dxyzdrs(0, 0) - my::deriv_(0, node) * dxyzdrs(1, 0)) *
              my::funct_(node) * fac;
          normalderiv(2, my::nsd_ * node + 2) += 0.;
        }
      else if (my::nsd_ == 2)
        for (int node = 0; node < my::bdrynen_; ++node)
        {
          normalderiv(0, my::nsd_ * node) += 0.;
          normalderiv(0, my::nsd_ * node + 1) += my::deriv_(0, node) * my::funct_(node) * fac;

          normalderiv(1, my::nsd_ * node) += -my::deriv_(0, node) * my::funct_(node) * fac;
          normalderiv(1, my::nsd_ * node + 1) += 0.;
        }
    }  // loop over gp

    // allocate auxiliary variable (= normalderiv^T * velocity)
    LINALG::Matrix<1, my::nsd_ * my::bdrynen_> temp(true);
    // allocate convective velocity at node
    LINALG::Matrix<1, my::nsd_> convvel(true);

    // elemat1.Shape(my::bdrynen_*my::numdofpernode_,my::bdrynen_*my::nsd_);
    // fill element matrix
    for (int inode = 0; inode < my::bdrynen_; inode++)
    {
      for (int i = 0; i < my::numdofpernode_; i++)
        nodenormal(i) = normal(inode * my::numdofpernode_ + i);

      double norm = nodenormal.Norm2();
      nodenormal.Scale(1 / norm);

      for (int idof = 0; idof < my::nsd_; idof++)
        convvel(idof) = evelnp(idof, inode) - egridvel(idof, inode);
      temp.Multiply(convvel, normalderiv);
      for (int idof = 0; idof < my::numdofpernode_; idof++)
      {
        // if(abs(nodenormal(idof)) > 0.5)
        if (mycondVector[inode * my::numdofpernode_ + idof] != 0.0)
        {
          for (int idof2 = 0; idof2 < my::nsd_; idof2++)
          {
            elemat1(inode * my::numdofpernode_ + idof, inode * my::nsd_ + idof2) +=
                temp(0, inode * my::nsd_ + idof2);
            elemat2(inode * my::numdofpernode_ + idof, inode * my::nsd_ + idof2) +=
                -nodenormal(idof2);
          }
          double normalconvvel = 0.0;
          for (int dim = 0; dim < my::nsd_; dim++) normalconvvel += convvel(dim) * nodenormal(dim);
          elevec1(inode * my::numdofpernode_ + idof) += -normalconvvel;
          break;
        }
      }
    }
  }  // coupling == "fluid structure"
  else
    dserror("unknown coupling type for no penetration boundary condition");

  return;
}

/*----------------------------------------------------------------------*
 * find dof IDs with no penetration condition  (protected)    vuong 06/11||
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::NoPenetrationIDs(
    DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, Epetra_SerialDenseVector& elevec1, std::vector<int>& lm)
{
  // This function is only implemented for 3D
  if (my::bdrynsd_ != 2 and my::bdrynsd_ != 1)
    dserror("NoPenetration is only implemented for 3D and 2D!");

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<my::bdrynsd_> intpoints(
      DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a my::nsd_ dimensional domain, since my::nsd_ determines the dimension of
  // FluidBoundary element!)
  // GEO::fillInitialPositionArray<distype,my::nsd_,Epetra_SerialDenseMatrix>(ele,my::xyze_);
  GEO::fillInitialPositionArray<distype, my::nsd_, LINALG::Matrix<my::nsd_, my::bdrynen_>>(
      ele, my::xyze_);

  // displacements
  Teuchos::RCP<const Epetra_Vector> dispnp;
  std::vector<double> mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp, mydispnp, lm);
    }
    dsassert(mydispnp.size() != 0, "no displacement values for boundary element");

    // Add the deformation of the ALE mesh to the nodes coordinates
    for (int inode = 0; inode < my::bdrynen_; ++inode)
      for (int idim = 0; idim < my::nsd_; ++idim)
        my::xyze_(idim, inode) += mydispnp[my::numdofpernode_ * inode + idim];
  }
  else
    dserror("fluid poro element not an ALE element!");

  // calculate normal
  Epetra_SerialDenseVector normal;
  normal.Size(lm.size());

  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the
    // shape function at the Gauss point Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(my::funct_, my::deriv_, my::fac_,
        my::unitnormal_, my::drs_, my::xsi_, my::xyze_, intpoints, gpid, NULL, NULL,
        IsNurbs<distype>::isnurbs);

    for (int inode = 0; inode < my::bdrynen_; ++inode)
    {
      for (int idim = 0; idim < my::nsd_; ++idim)
        normal(inode * my::numdofpernode_ + idim) +=
            my::unitnormal_(idim) * my::funct_(inode) * my::fac_;
      // pressure dof is set to zero
      normal(inode * my::numdofpernode_ + (my::nsd_)) = 0.0;
    }
  } /* end of loop over integration points gpid */

  LINALG::Matrix<my::numdofpernode_, 1> nodenormal(true);

  // fill element matrix
  for (int inode = 0; inode < my::bdrynen_; inode++)
  {
    for (int i = 0; i < my::numdofpernode_; i++)
      nodenormal(i) = normal(inode * my::numdofpernode_ + i);
    double norm = nodenormal.Norm2();
    nodenormal.Scale(1 / norm);

    bool isset = false;
    for (int idof = 0; idof < my::numdofpernode_; idof++)
    {
      if (isset == false and abs(nodenormal(idof)) > 0.5)
      {
        elevec1(inode * my::numdofpernode_ + idof) = 1.0;
        isset = true;
      }
      else  // no condition set on dof
        elevec1(inode * my::numdofpernode_ + idof) = 0.0;
    }
  }
  return;
}

/*------------------------------------------------------------------------------*
 * evaluate no penetration condition (weak substitution)  (protected) vuong 06/11||
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::PoroBoundary(
    DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& plm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseVector& elevec1)
{
  switch (distype)
  {
    // 2D:
    case DRT::Element::line2:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::quad4)
      {
        PoroBoundary<DRT::Element::quad4>(ele, params, discretization, plm, elemat1, elevec1);
      }
      else if (ele->ParentElement()->Shape() == DRT::Element::tri3)
      {
        PoroBoundary<DRT::Element::tri3>(ele, params, discretization, plm, elemat1, elevec1);
      }
      else
      {
        dserror("expected combination line2/quad4 or line2/tri3 for line/parent pair");
      }
      break;
    }
    case DRT::Element::line3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::quad9)
      {
        PoroBoundary<DRT::Element::quad9>(ele, params, discretization, plm, elemat1, elevec1);
      }
      else
      {
        dserror("expected combination line3/quad9 for line/parent pair");
      }
      break;
    }
    case DRT::Element::nurbs3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::nurbs9)
      {
        PoroBoundary<DRT::Element::nurbs9>(ele, params, discretization, plm, elemat1, elevec1);
      }
      else
      {
        dserror("expected combination nurbs3/nurbs9 for line/parent pair");
      }
      break;
    }
    // 3D:
    case DRT::Element::quad4:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::hex8)
      {
        PoroBoundary<DRT::Element::hex8>(ele, params, discretization, plm, elemat1, elevec1);
      }
      else
      {
        dserror("expected combination quad4/hex8 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::tri3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::tet4)
      {
        PoroBoundary<DRT::Element::tet4>(ele, params, discretization, plm, elemat1, elevec1);
      }
      else
      {
        dserror("expected combination tri3/tet4 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::tri6:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::tet10)
      {
        PoroBoundary<DRT::Element::tet10>(ele, params, discretization, plm, elemat1, elevec1);
      }
      else
      {
        dserror("expected combination tri6/tet10 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::quad9:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::hex27)
      {
        PoroBoundary<DRT::Element::hex27>(ele, params, discretization, plm, elemat1, elevec1);
      }
      else
      {
        dserror("expected combination hex27/hex27 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::nurbs9:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::nurbs27)
      {
        PoroBoundary<DRT::Element::nurbs27>(ele, params, discretization, plm, elemat1, elevec1);
      }
      else
      {
        dserror("expected combination nurbs9/nurbs27 for line/parent pair");
      }
      break;
    }
    default:
    {
      dserror("surface/parent element pair not yet implemented. Just do it.\n");
      break;
    }
  }
}

/*------------------------------------------------------------------------------*
 * evaluate no penetration condition (weak substitution)  (protected) vuong 06/11||
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType pdistype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::PoroBoundary(
    DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& plm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseVector& elevec1)
{
  // This function is only implemented for 3D and 2D
  if (my::bdrynsd_ != 2 and my::bdrynsd_ != 1)
    dserror("PoroBoundary is only implemented for 3D and 2D!");

  POROELAST::coupltype coupling =
      params.get<POROELAST::coupltype>("coupling", POROELAST::undefined);
  if (coupling == POROELAST::undefined) dserror("no coupling defined for poro-boundary condition");
  const bool offdiag(coupling == POROELAST::fluidstructure);

  // get timescale parameter from parameter list (depends on time integration scheme)
  double timescale = params.get<double>("timescale", -1.0);
  if (timescale == -1.0 and offdiag) dserror("no timescale parameter in parameter list");

  // reset timescale in stationary case
  if (my::fldparatimint_->IsStationary()) timescale = 0.0;

  // get element location vector and ownerships
  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;
  ele->DRT::Element::LocationVector(discretization, lm, lmowner, lmstride);

  /// number of parentnodes
  static const int nenparent = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

  // get the parent element
  DRT::ELEMENTS::Fluid* pele = ele->ParentElement();

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<my::bdrynsd_> intpoints(
      DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a my::nsd_ dimensional domain, since my::nsd_ determines the dimension of
  // FluidBoundary element!)
  GEO::fillInitialPositionArray<distype, my::nsd_, LINALG::Matrix<my::nsd_, my::bdrynen_>>(
      ele, my::xyze_);

  // displacements
  Teuchos::RCP<const Epetra_Vector> dispnp;
  std::vector<double> mydispnp;
  std::vector<double> parentdispnp;

  dispnp = discretization.GetState("dispnp");
  if (dispnp != Teuchos::null)
  {
    mydispnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp, mydispnp, lm);
    DRT::UTILS::ExtractMyValues(*dispnp, parentdispnp, plm);
  }
  dsassert(mydispnp.size() != 0, "no displacement values for boundary element");
  dsassert(parentdispnp.size() != 0, "no displacement values for parent element");

  // Add the deformation of the ALE mesh to the nodes coordinates
  for (int inode = 0; inode < my::bdrynen_; ++inode)
    for (int idim = 0; idim < my::nsd_; ++idim)
      my::xyze_(idim, inode) += mydispnp[my::numdofpernode_ * inode + idim];

  // update element geometry of parent element
  LINALG::Matrix<my::nsd_, nenparent> xrefe;  // material coord. of parent element
  LINALG::Matrix<my::nsd_, nenparent> xcurr;  // current  coord. of parent element
  {
    DRT::Node** nodes = pele->Nodes();
    for (int i = 0; i < nenparent; ++i)
    {
      for (int j = 0; j < my::nsd_; ++j)
      {
        const double* x = nodes[i]->X();
        xrefe(j, i) = x[j];
        xcurr(j, i) = xrefe(j, i) + parentdispnp[i * my::numdofpernode_ + j];
      }
    }
  }

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
  Teuchos::RCP<const Epetra_Vector> gridvel = discretization.GetState("gridv");
  Teuchos::RCP<const Epetra_Vector> scaaf = discretization.GetState("scaaf");

  if (velnp == Teuchos::null) dserror("Cannot get state vector 'velnp'");
  if (gridvel == Teuchos::null) dserror("Cannot get state vector 'gridv'");

  std::vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp, myvelnp, lm);
  std::vector<double> mygridvel(lm.size());
  DRT::UTILS::ExtractMyValues(*gridvel, mygridvel, lm);
  std::vector<double> myscaaf(lm.size());
  DRT::UTILS::ExtractMyValues(*scaaf, myscaaf, lm);

  // allocate velocity vectors
  LINALG::Matrix<my::nsd_, my::bdrynen_> evelnp(true);
  LINALG::Matrix<my::bdrynen_, 1> epressnp(true);
  LINALG::Matrix<my::nsd_, my::bdrynen_> edispnp(true);
  LINALG::Matrix<my::nsd_, my::bdrynen_> egridvel(true);
  LINALG::Matrix<my::bdrynen_, 1> escaaf(true);
  LINALG::Matrix<my::bdrynen_, 1> eporosity(true);

  // split velocity and pressure, insert into element arrays
  for (int inode = 0; inode < my::bdrynen_; inode++)
  {
    for (int idim = 0; idim < my::nsd_; idim++)
    {
      evelnp(idim, inode) = myvelnp[idim + (inode * my::numdofpernode_)];
      edispnp(idim, inode) = mydispnp[idim + (inode * my::numdofpernode_)];
      egridvel(idim, inode) = mygridvel[idim + (inode * my::numdofpernode_)];
    }
    epressnp(inode) = myvelnp[my::nsd_ + (inode * my::numdofpernode_)];
    escaaf(inode) = myscaaf[my::nsd_ + (inode * my::numdofpernode_)];
  }

  const bool porositydof = ComputeNodalPorosity(ele, mydispnp, eporosity);

  // get coordinates of gauss points w.r.t. local parent coordinate system
  LINALG::SerialDenseMatrix pqxg(intpoints.IP().nquad, my::nsd_);
  LINALG::Matrix<my::nsd_, my::nsd_> derivtrafo(true);

  DRT::UTILS::BoundaryGPToParentGP<my::nsd_>(
      pqxg, derivtrafo, intpoints, pdistype, distype, ele->SurfaceNumber());

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  // --------------------------------------------------

  // In the case of nurbs the normal vector is multiplied with normalfac
  double normalfac = 0.0;
  std::vector<Epetra_SerialDenseVector> mypknots(my::nsd_);
  std::vector<Epetra_SerialDenseVector> myknots(my::bdrynsd_);
  Epetra_SerialDenseVector weights(my::bdrynen_);
  Epetra_SerialDenseVector pweights(pele->NumNode());

  // for isogeometric elements --- get knotvectors for parent
  // element and surface element, get weights
  if (IsNurbs<distype>::isnurbs)
  {
    bool zero_size = DRT::NURBS::GetKnotVectorAndWeightsForNurbsBoundaryAndParent(pele, ele,
        ele->SurfaceNumber(), discretization, mypknots, myknots, pweights, weights, normalfac);

    if (zero_size)
    {
      return;
    }
  }
  // --------------------------------------------------
  // structure velocity at gausspoint
  LINALG::Matrix<my::nsd_, 1> gridvelint;

  // coordinates of gauss points of parent element
  LINALG::Matrix<my::nsd_, 1> pxsi(true);

  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    // get shape functions and derivatives in the plane of the element
    LINALG::Matrix<nenparent, 1> pfunct(true);
    LINALG::Matrix<my::nsd_, nenparent> pderiv;
    LINALG::Matrix<my::nsd_, nenparent> pderiv_loc;

    // coordinates of the current integration point
    for (int idim = 0; idim < my::nsd_; idim++) pxsi(idim) = pqxg(gpid, idim);

    // get shape functions and derivatives of the parent element
    if (not IsNurbs<distype>::isnurbs)
    {
      // shape functions and their first derivatives of parent element
      DRT::UTILS::shape_function<pdistype>(pxsi, pfunct);
      DRT::UTILS::shape_function_deriv1<pdistype>(pxsi, pderiv_loc);
    }
    // only for NURBS!!!
    else
    {
      DRT::NURBS::UTILS::nurbs_get_funct_deriv(
          pfunct, pderiv_loc, pxsi, mypknots, pweights, pdistype);
    }
    pderiv.MultiplyTN(derivtrafo, pderiv_loc);

    // get Jacobian matrix and determinant w.r.t. spatial configuration
    // transposed jacobian "dx/ds"
    LINALG::Matrix<my::nsd_, my::nsd_> xjm;
    LINALG::Matrix<my::nsd_, my::nsd_> Jmat;
    xjm.MultiplyNT(pderiv_loc, xcurr);
    Jmat.MultiplyNT(pderiv_loc, xrefe);
    // jacobian determinant "det(dx/ds)"
    const double det = xjm.Determinant();
    // jacobian determinant "det(dX/ds)"
    const double detJ = Jmat.Determinant();
    // jacobian determinant "det(dx/dX) = det(dx/ds)/det(dX/ds)"
    const double J = det / detJ;

    // Computation of the integration factor & shape function at the Gauss point & derivative of the
    // shape function at the Gauss point Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(my::funct_, my::deriv_, my::fac_,
        my::unitnormal_, my::drs_, my::xsi_, my::xyze_, intpoints, gpid, &myknots, &weights,
        IsNurbs<distype>::isnurbs);

    // in the case of nurbs the normal vector must be scaled with a special factor
    if (IsNurbs<distype>::isnurbs) my::unitnormal_.Scale(normalfac);

    const double timefacpre = my::fldparatimint_->TimeFacPre();
    const double timefacfacpre = my::fldparatimint_->TimeFacPre() * my::fac_;
    const double rhsfac = my::fldparatimint_->TimeFacRhs() * my::fac_;

    my::velint_.Multiply(evelnp, my::funct_);
    gridvelint.Multiply(egridvel, my::funct_);
    double press = epressnp.Dot(my::funct_);

    double scalar = escaaf.Dot(my::funct_);

    double dphi_dp = 0.0;
    double dphi_dJ = 0.0;
    double porosity_gp = 0.0;

    params.set<double>("scalar", scalar);

    ComputePorosityAtGP(
        params, ele, my::funct_, eporosity, press, J, gpid, porosity_gp, dphi_dp, dphi_dJ, false);

    // The integration factor is not multiplied with drs
    // since it is the same as the scaling factor for the unit normal derivatives
    // Therefore it cancels out!!
    const double fac = intpoints.IP().qwgt[gpid];

    //  derivatives of surface normals wrt mesh displacements
    LINALG::Matrix<my::nsd_, nenparent * my::nsd_> normalderiv(true);

    // dxyzdrs vector -> normal which is not normalized
    LINALG::Matrix<my::bdrynsd_, my::nsd_> dxyzdrs(0.0);
    dxyzdrs.MultiplyNT(my::deriv_, my::xyze_);

    if (my::nsd_ == 3)
      for (int node = 0; node < nenparent; ++node)
      {
        normalderiv(0, my::nsd_ * node) += 0.;
        normalderiv(0, my::nsd_ * node + 1) +=
            (pderiv(0, node) * dxyzdrs(1, 2) - pderiv(1, node) * dxyzdrs(0, 2));
        normalderiv(0, my::nsd_ * node + 2) +=
            (pderiv(1, node) * dxyzdrs(0, 1) - pderiv(0, node) * dxyzdrs(1, 1));

        normalderiv(1, my::nsd_ * node) +=
            (pderiv(1, node) * dxyzdrs(0, 2) - pderiv(0, node) * dxyzdrs(1, 2));
        normalderiv(1, my::nsd_ * node + 1) += 0.;
        normalderiv(1, my::nsd_ * node + 2) +=
            (pderiv(0, node) * dxyzdrs(1, 0) - pderiv(1, node) * dxyzdrs(0, 0));

        normalderiv(2, my::nsd_ * node) +=
            (pderiv(0, node) * dxyzdrs(1, 1) - pderiv(1, node) * dxyzdrs(0, 1));
        normalderiv(2, my::nsd_ * node + 1) +=
            (pderiv(1, node) * dxyzdrs(0, 0) - pderiv(0, node) * dxyzdrs(1, 0));
        normalderiv(2, my::nsd_ * node + 2) += 0.;
      }
    else  // if(my::nsd_==2)
      for (int node = 0; node < nenparent; ++node)
      {
        normalderiv(0, my::nsd_ * node) += 0.;
        normalderiv(0, my::nsd_ * node + 1) += pderiv(0, node);

        normalderiv(1, my::nsd_ * node) += -pderiv(0, node);
        normalderiv(1, my::nsd_ * node + 1) += 0.;
      }

    // in the case of nurbs the normal vector must be scaled with a special factor
    if (IsNurbs<distype>::isnurbs) normalderiv.Scale(normalfac);

    //------------------------------------------------dJ/dus = dJ/dF : dF/dus = J * F^-T . N_X = J *
    // N_x
    LINALG::Matrix<1, my::nsd_ * nenparent> dJ_dus;
    // global derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<my::nsd_, nenparent> derxy;
    // inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<my::nsd_, my::nsd_> xji;

    xji.Invert(xjm);
    derxy.Multiply(xji, pderiv_loc);

    for (int i = 0; i < nenparent; i++)
      for (int j = 0; j < my::nsd_; j++) dJ_dus(j + i * my::nsd_) = J * derxy(j, i);

    double normal_convel = 0.0;
    LINALG::Matrix<1, my::nsd_> convel;

    for (int idof = 0; idof < my::nsd_; idof++)
    {
      normal_convel += my::unitnormal_(idof) * my::velint_(idof);
      convel(idof) = my::velint_(idof);
    }

    if (not my::fldparatimint_->IsStationary())
      for (int idof = 0; idof < my::nsd_; idof++)
      {
        normal_convel += my::unitnormal_(idof) * (-gridvelint(idof));
        convel(idof) -= gridvelint(idof);
      }

    LINALG::Matrix<1, nenparent * my::nsd_> tmp;
    tmp.Multiply(convel, normalderiv);

    // fill element matrix
    {
      if (not offdiag)
      {
        for (int inode = 0; inode < nenparent; inode++)
          elevec1(inode * my::numdofpernode_ + my::nsd_) -=
              rhsfac * pfunct(inode) * porosity_gp * normal_convel;

        for (int inode = 0; inode < nenparent; inode++)
          for (int nnod = 0; nnod < nenparent; nnod++)
          {
            for (int idof2 = 0; idof2 < my::nsd_; idof2++)
              elemat1(inode * my::numdofpernode_ + my::nsd_, nnod * my::numdofpernode_ + idof2) +=
                  timefacfacpre * pfunct(inode) * porosity_gp * my::unitnormal_(idof2) *
                  pfunct(nnod);
            elemat1(inode * my::numdofpernode_ + my::nsd_, nnod * my::numdofpernode_ + my::nsd_) +=
                +timefacfacpre * pfunct(inode) * dphi_dp * normal_convel * pfunct(nnod);
          }
      }

      else if (not porositydof)
      {
        for (int inode = 0; inode < nenparent; inode++)
          for (int nnod = 0; nnod < nenparent; nnod++)
            for (int idof2 = 0; idof2 < my::nsd_; idof2++)
              elemat1(inode * my::numdofpernode_ + my::nsd_, nnod * my::nsd_ + idof2) +=
                  +tmp(0, nnod * my::nsd_ + idof2) * porosity_gp * pfunct(inode) * timefacpre *
                      fac -
                  pfunct(inode) * porosity_gp * my::unitnormal_(idof2) * timescale * pfunct(nnod) *
                      timefacfacpre +
                  pfunct(inode) * dphi_dJ * dJ_dus(nnod * my::nsd_ + idof2) * normal_convel *
                      timefacfacpre;
      }

      else  // offdiagonal and porositydof
        for (int inode = 0; inode < nenparent; inode++)
          for (int nnod = 0; nnod < nenparent; nnod++)
          {
            for (int idof2 = 0; idof2 < my::nsd_; idof2++)
              elemat1(inode * my::numdofpernode_ + my::nsd_, nnod * (my::nsd_ + 1) + idof2) +=
                  +tmp(0, nnod * my::nsd_ + idof2) * porosity_gp * pfunct(inode) * timefacpre *
                      fac -
                  pfunct(inode) * porosity_gp * my::unitnormal_(idof2) * timescale * pfunct(nnod) *
                      timefacfacpre +
                  pfunct(inode) * dphi_dJ * dJ_dus(nnod * my::nsd_ + idof2) * normal_convel *
                      timefacfacpre;
            elemat1(inode * my::numdofpernode_ + my::nsd_, nnod * (my::nsd_ + 1) + my::nsd_) +=
                pfunct(inode) * pfunct(nnod) * normal_convel * timefacfacpre;
          }
    }
  } /* end of loop over integration points gpid */
  return;
}

/*----------------------------------------------------------------------*
 * evaluate boundary pressure                     (protected) vuong 06/11||
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::PressureCoupling(
    DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseVector& elevec1)
{
  // This function is only implemented for 3D
  if (my::bdrynsd_ != 2 and my::bdrynsd_ != 1)
    dserror("PressureCoupling is only implemented for 2D and 3D!");

  POROELAST::coupltype coupling =
      params.get<POROELAST::coupltype>("coupling", POROELAST::undefined);
  if (coupling == POROELAST::undefined) dserror("no coupling defined for poro-boundary condition");
  const bool offdiag(coupling == POROELAST::fluidstructure);

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<my::bdrynsd_> intpoints(
      DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a my::nsd_ dimensional domain, since my::nsd_ determines the dimension of
  // FluidBoundary element!)
  GEO::fillInitialPositionArray<distype, my::nsd_, LINALG::Matrix<my::nsd_, my::bdrynen_>>(
      ele, my::xyze_);

  // displacements
  Teuchos::RCP<const Epetra_Vector> dispnp;
  std::vector<double> mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp, mydispnp, lm);
    }
    dsassert(mydispnp.size() != 0, "no displacement values for boundary element");

    // Add the deformation of the ALE mesh to the nodes coordinates
    for (int inode = 0; inode < my::bdrynen_; ++inode)
    {
      for (int idim = 0; idim < my::nsd_; ++idim)
      {
        my::xyze_(idim, inode) += mydispnp[my::numdofpernode_ * inode + idim];
      }
    }
  }

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");

  if (velnp == Teuchos::null) dserror("Cannot get state vector 'velnp'");

  std::vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp, myvelnp, lm);

  // allocate velocity vectors
  LINALG::Matrix<my::bdrynen_, 1> epressnp(true);

  // split velocity and pressure, insert into element arrays
  for (int inode = 0; inode < my::bdrynen_; inode++)
  {
    epressnp(inode) = myvelnp[my::nsd_ + (inode * my::numdofpernode_)];
  }

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  // --------------------------------------------------

  // In the case of nurbs the normal vector is multiplied with normalfac
  double normalfac = 0.0;
  std::vector<Epetra_SerialDenseVector> mypknots(my::nsd_);
  std::vector<Epetra_SerialDenseVector> myknots(my::bdrynsd_);
  Epetra_SerialDenseVector weights(my::bdrynen_);

  // for isogeometric elements --- get knotvectors for parent
  // element and surface element, get weights
  if (IsNurbs<distype>::isnurbs)
  {
    bool zero_size = DRT::NURBS::GetKnotVectorAndWeightsForNurbsBoundary(ele, ele->SurfaceNumber(),
        ele->ParentElement()->Id(), discretization, mypknots, myknots, weights, normalfac);
    if (zero_size)
    {
      return;
    }
  }
  // --------------------------------------------------

  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the
    // shape function at the Gauss point Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(my::funct_, my::deriv_, my::fac_,
        my::unitnormal_, my::drs_, my::xsi_, my::xyze_, intpoints, gpid, &myknots, &weights,
        IsNurbs<distype>::isnurbs);

    const double timefac = my::fldparatimint_->TimeFac();
    const double timefacfac = my::fldparatimint_->TimeFac() * my::fac_;
    const double rhsfac = my::fldparatimint_->TimeFacRhs() * my::fac_;

    // get pressure at integration point
    double press = my::funct_.Dot(epressnp);

    // dxyzdrs vector -> normal which is not normalized
    LINALG::Matrix<my::bdrynsd_, my::nsd_> dxyzdrs(0.0);
    dxyzdrs.MultiplyNT(my::deriv_, my::xyze_);

    // in the case of nurbs the normal vector must be scaled with a special factor
    if (IsNurbs<distype>::isnurbs) my::unitnormal_.Scale(normalfac);

    //  derivatives of surface normals wrt mesh displacements
    LINALG::Matrix<my::nsd_, my::bdrynen_ * my::nsd_> normalderiv(true);

    // The integration factor is not multiplied with drs
    // since it is the same as the scaling factor for the unit normal derivatives
    // Therefore it cancels out!!
    const double fac = intpoints.IP().qwgt[gpid];

    if (my::nsd_ == 3)
      for (int node = 0; node < my::bdrynen_; ++node)
      {
        normalderiv(0, 3 * node) += 0.;
        normalderiv(0, 3 * node + 1) +=
            (my::deriv_(0, node) * dxyzdrs(1, 2) - my::deriv_(1, node) * dxyzdrs(0, 2));
        normalderiv(0, 3 * node + 2) +=
            (my::deriv_(1, node) * dxyzdrs(0, 1) - my::deriv_(0, node) * dxyzdrs(1, 1));

        normalderiv(1, 3 * node) +=
            (my::deriv_(1, node) * dxyzdrs(0, 2) - my::deriv_(0, node) * dxyzdrs(1, 2));
        normalderiv(1, 3 * node + 1) += 0.;
        normalderiv(1, 3 * node + 2) +=
            (my::deriv_(0, node) * dxyzdrs(1, 0) - my::deriv_(1, node) * dxyzdrs(0, 0));

        normalderiv(2, 3 * node) +=
            (my::deriv_(0, node) * dxyzdrs(1, 1) - my::deriv_(1, node) * dxyzdrs(0, 1));
        normalderiv(2, 3 * node + 1) +=
            (my::deriv_(1, node) * dxyzdrs(0, 0) - my::deriv_(0, node) * dxyzdrs(1, 0));
        normalderiv(2, 3 * node + 2) += 0.;
      }
    else if (my::nsd_ == 2)
      for (int node = 0; node < my::bdrynen_; ++node)
      {
        normalderiv(0, my::nsd_ * node) += 0.;
        normalderiv(0, my::nsd_ * node + 1) += my::deriv_(0, node);

        normalderiv(1, my::nsd_ * node) += -my::deriv_(0, node);
        normalderiv(1, my::nsd_ * node + 1) += 0.;
      }

    // in the case of nurbs the normal vector must be scaled with a special factor
    if (IsNurbs<distype>::isnurbs) normalderiv.Scale(normalfac);

    // fill element matrix
    for (int inode = 0; inode < my::bdrynen_; inode++)
    {
      for (int idof = 0; idof < my::nsd_; idof++)
      {
        if (not offdiag)
          elevec1(inode * my::numdofpernode_ + idof) -=
              my::funct_(inode) * my::unitnormal_(idof) * press * rhsfac;
        for (int nnod = 0; nnod < my::bdrynen_; nnod++)
        {
          if (not offdiag)
            elemat1(inode * my::numdofpernode_ + idof, nnod * my::numdofpernode_ + my::nsd_) +=
                my::funct_(inode) * my::unitnormal_(idof) * my::funct_(nnod) * timefacfac;
          else
            for (int idof2 = 0; idof2 < my::nsd_; idof2++)
            {
              elemat1(inode * my::numdofpernode_ + idof, nnod * my::nsd_ + idof2) +=
                  normalderiv(idof, nnod * my::nsd_ + idof2) * press * my::funct_(inode) * timefac *
                  fac;
            }
        }
      }
    }
  } /* end of loop over integration points gpid */

  return;
}

/*----------------------------------------------------------------------*
 *  compute porosity at gauss point               (protected) vuong 06/11||
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::ComputePorosityAtGP(
    Teuchos::ParameterList& params, DRT::ELEMENTS::FluidBoundary* ele,
    const LINALG::Matrix<my::bdrynen_, 1>& funct, const LINALG::Matrix<my::bdrynen_, 1>& eporosity,
    double press, double J, int gp, double& porosity, double& dphi_dp, double& dphi_dJ, bool save)
{
  Teuchos::RCP<MAT::StructPoro> structmat =
      Teuchos::rcp_dynamic_cast<MAT::StructPoro>(ele->ParentElement()->Material(1));
  structmat->ComputeSurfPorosity(params, press, J, ele->SurfaceNumber(), gp, porosity, &dphi_dp,
      &dphi_dJ,
      NULL,  // dphi_dJdp not needed
      NULL,  // dphi_dJJ not needed
      NULL,  // dphi_dpp not needed
      save);
}

/*----------------------------------------------------------------------------------*
 * evaluate no penetration condition (lagrange multiplier)   (protected) vuong 11/14|
 *----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::NoPenetrationMatAndRHS(
    DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& k_fluid,
    Epetra_SerialDenseVector& rhs)
{
  switch (distype)
  {
    // 2D:
    case DRT::Element::line2:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::quad4)
      {
        NoPenetrationMatAndRHS<DRT::Element::quad4>(ele, params, discretization, lm, k_fluid, rhs);
      }
      else if (ele->ParentElement()->Shape() == DRT::Element::tri3)
      {
        NoPenetrationMatAndRHS<DRT::Element::tri3>(ele, params, discretization, lm, k_fluid, rhs);
      }
      else
      {
        dserror("expected combination line2/quad4 or line2/tri3 for line/parent pair");
      }
      break;
    }
    case DRT::Element::line3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::quad9)
      {
        NoPenetrationMatAndRHS<DRT::Element::quad9>(ele, params, discretization, lm, k_fluid, rhs);
      }
      else
      {
        dserror("expected combination line3/quad9 for line/parent pair");
      }
      break;
    }
    case DRT::Element::nurbs3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::nurbs9)
      {
        NoPenetrationMatAndRHS<DRT::Element::nurbs9>(ele, params, discretization, lm, k_fluid, rhs);
      }
      else
      {
        dserror("expected combination nurbs3/nurbs9 for line/parent pair");
      }
      break;
    }
    // 3D:
    case DRT::Element::quad4:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::hex8)
      {
        NoPenetrationMatAndRHS<DRT::Element::hex8>(ele, params, discretization, lm, k_fluid, rhs);
      }
      else
      {
        dserror("expected combination quad4/hex8 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::tri3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::tet4)
      {
        NoPenetrationMatAndRHS<DRT::Element::tet4>(ele, params, discretization, lm, k_fluid, rhs);
      }
      else
      {
        dserror("expected combination tri3/tet4 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::tri6:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::tet10)
      {
        NoPenetrationMatAndRHS<DRT::Element::tet10>(ele, params, discretization, lm, k_fluid, rhs);
      }
      else
      {
        dserror("expected combination tri6/tet10 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::quad9:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::hex27)
      {
        NoPenetrationMatAndRHS<DRT::Element::hex27>(ele, params, discretization, lm, k_fluid, rhs);
      }
      else
      {
        dserror("expected combination hex27/hex27 for surface/parent pair");
      }
      break;
    }
    default:
    {
      dserror("surface/parent element pair not yet implemented. Just do it.\n");
      break;
    }
  }
}

/*----------------------------------------------------------------------------------*
 * evaluate no penetration condition (lagrange multiplier)    (protected) vuong 06/11||
 *----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType pdistype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::NoPenetrationMatAndRHS(
    DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& k_fluid,
    Epetra_SerialDenseVector& rhs)
{
  // This function is only implemented for 3D
  if (my::bdrynsd_ != 2 and my::bdrynsd_ != 1)
    dserror("PressureCoupling is only implemented for 2D and 3D!");

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<my::bdrynsd_> intpoints(
      DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a my::nsd_ dimensional domain, since my::nsd_ determines the dimension of
  // FluidBoundary element!)
  GEO::fillInitialPositionArray<distype, my::nsd_, LINALG::Matrix<my::nsd_, my::bdrynen_>>(
      ele, my::xyze_);

  // displacements
  Teuchos::RCP<const Epetra_Vector> dispnp;
  std::vector<double> mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp, mydispnp, lm);
    }
    dsassert(mydispnp.size() != 0, "no displacement values for boundary element");

    // Add the deformation of the ALE mesh to the nodes coordinates
    for (int inode = 0; inode < my::bdrynen_; ++inode)
    {
      for (int idim = 0; idim < my::nsd_; ++idim)
      {
        my::xyze_(idim, inode) += mydispnp[my::nsd_ * inode + idim];
      }
    }
  }

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
  Teuchos::RCP<const Epetra_Vector> gridvel = discretization.GetState("gridv");

  if (velnp == Teuchos::null) dserror("Cannot get state vector 'velnp'");

  std::vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp, myvelnp, lm);

  std::vector<double> mygridvel(lm.size());
  DRT::UTILS::ExtractMyValues(*gridvel, mygridvel, lm);

  // allocate velocity vectors
  LINALG::Matrix<my::nsd_, my::bdrynen_> evelnp(true);
  LINALG::Matrix<my::nsd_, my::bdrynen_> egridvel(true);

  // split velocity and pressure, insert into element arrays
  for (int inode = 0; inode < my::bdrynen_; inode++)
    for (int idim = 0; idim < my::nsd_; idim++)
    {
      evelnp(idim, inode) = myvelnp[idim + (inode * my::nsd_)];
      egridvel(idim, inode) = mygridvel[idim + (inode * my::nsd_)];
    }

  // allocate convective velocity at node
  LINALG::Matrix<my::nsd_, my::bdrynen_> econvvel(true);
  econvvel += evelnp;
  if (not my::fldparatimint_->IsStationary()) econvvel -= egridvel;

  // --------------------------------------------------
  // parent element
  // --------------------------------------------------

  // get the parent element
  DRT::ELEMENTS::Fluid* pele = ele->ParentElement();

  /// number of parentnodes
  static const int nenparent = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

  // get element location vector and ownerships
  std::vector<int> plm;
  std::vector<int> plmowner;
  std::vector<int> plmstride;
  pele->DRT::Element::LocationVector(discretization, plm, plmowner, plmstride);

  std::vector<double> parentdispnp;
  DRT::UTILS::ExtractMyValues(*dispnp, parentdispnp, plm);

  // update element geometry of parent element
  LINALG::Matrix<my::nsd_, nenparent> xrefe;  // material coord. of parent element
  LINALG::Matrix<my::nsd_, nenparent> xcurr;  // current  coord. of parent element
  {
    DRT::Node** nodes = pele->Nodes();
    for (int i = 0; i < nenparent; ++i)
    {
      for (int j = 0; j < my::nsd_; ++j)
      {
        const double* x = nodes[i]->X();
        xrefe(j, i) = x[j];
        xcurr(j, i) = xrefe(j, i) + parentdispnp[i * my::numdofpernode_ + j];
      }
    }
  }


  std::vector<double> pvelnp(plm.size());
  DRT::UTILS::ExtractMyValues(*velnp, pvelnp, plm);

  // allocate vectors
  LINALG::Matrix<nenparent, 1> pepressnp(true);

  // split velocity and pressure, insert into element arrays
  for (int inode = 0; inode < nenparent; inode++)
  {
    pepressnp(inode) = pvelnp[my::nsd_ + (inode * my::numdofpernode_)];
  }

  // get coordinates of gauss points w.r.t. local parent coordinate system
  LINALG::SerialDenseMatrix pqxg(intpoints.IP().nquad, my::nsd_);
  LINALG::Matrix<my::nsd_, my::nsd_> derivtrafo(true);

  DRT::UTILS::BoundaryGPToParentGP<my::nsd_>(
      pqxg, derivtrafo, intpoints, pdistype, distype, ele->SurfaceNumber());

  // coordinates of gauss points of parent element
  LINALG::Matrix<my::nsd_, 1> pxsi(true);

  LINALG::Matrix<my::bdrynen_, 1> eporosity(true);

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  // --------------------------------------------------

  // In the case of nurbs the normal vector is multiplied with normalfac
  double normalfac = 0.0;
  std::vector<Epetra_SerialDenseVector> mypknots(my::nsd_);
  std::vector<Epetra_SerialDenseVector> myknots(my::bdrynsd_);
  Epetra_SerialDenseVector weights(my::bdrynen_);
  Epetra_SerialDenseVector pweights(pele->NumNode());

  // for isogeometric elements --- get knotvectors for parent
  // element and surface element, get weights
  if (IsNurbs<distype>::isnurbs)
  {
    bool zero_size = DRT::NURBS::GetKnotVectorAndWeightsForNurbsBoundaryAndParent(pele, ele,
        ele->SurfaceNumber(), discretization, mypknots, myknots, pweights, weights, normalfac);

    if (zero_size)
    {
      return;
    }
  }
  // --------------------------------------------------

  // --------------------------------------------------

  //! array for shape functions for boundary element
  // LINALG::Matrix<my::bdrynen_,1> dualfunct(true);

  // allocate convective velocity at gauss point
  LINALG::Matrix<my::nsd_, 1> convvel(true);

  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the
    // shape function at the Gauss point Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(my::funct_, my::deriv_, my::fac_,
        my::unitnormal_, my::drs_, my::xsi_, my::xyze_, intpoints, gpid, &myknots, &weights,
        IsNurbs<distype>::isnurbs);

    // --------------------------------------------------
    // parent element
    // --------------------------------------------------
    // get shape functions and derivatives in the plane of the element
    LINALG::Matrix<nenparent, 1> pfunct(true);
    LINALG::Matrix<my::nsd_, nenparent> pderiv_loc;

    // coordinates of the current integration point
    for (int idim = 0; idim < my::nsd_; idim++) pxsi(idim) = pqxg(gpid, idim);

    // get shape functions and derivatives of the parent element
    if (not IsNurbs<distype>::isnurbs)
    {
      // shape functions and their first derivatives of parent element
      DRT::UTILS::shape_function<pdistype>(pxsi, pfunct);
      DRT::UTILS::shape_function_deriv1<pdistype>(pxsi, pderiv_loc);
    }
    // only for NURBS!!!
    else
    {
      DRT::NURBS::UTILS::nurbs_get_funct_deriv(
          pfunct, pderiv_loc, pxsi, mypknots, pweights, pdistype);
    }

    // get Jacobian matrix and determinant w.r.t. spatial configuration
    // transposed jacobian "dx/ds"
    LINALG::Matrix<my::nsd_, my::nsd_> xjm;
    LINALG::Matrix<my::nsd_, my::nsd_> Jmat;
    xjm.MultiplyNT(pderiv_loc, xcurr);
    Jmat.MultiplyNT(pderiv_loc, xrefe);
    // jacobian determinant "det(dx/ds)"
    const double det = xjm.Determinant();
    // jacobian determinant "det(dX/ds)"
    const double detJ = Jmat.Determinant();
    // jacobian determinant "det(dx/dX) = det(dx/ds)/det(dX/ds)"
    const double J = det / detJ;

    double press = pepressnp.Dot(pfunct);

    double dphi_dp = 0.0;
    double dphi_dJ = 0.0;
    double porosity_gp = 0.0;

    ComputePorosityAtGP(
        params, ele, my::funct_, eporosity, press, J, gpid, porosity_gp, dphi_dp, dphi_dJ, false);

    // --------------------------------------------------
    // --------------------------------------------------

    // dxyzdrs vector -> normal which is not normalized
    LINALG::Matrix<my::bdrynsd_, my::nsd_> dxyzdrs(0.0);
    dxyzdrs.MultiplyNT(my::deriv_, my::xyze_);

    // in the case of nurbs the normal vector must be scaled with a special factor
    if (IsNurbs<distype>::isnurbs) my::unitnormal_.Scale(normalfac);

    convvel.Multiply(econvvel, my::funct_);

    // fill element matrix and rhs
    for (int inode = 0; inode < my::bdrynen_; inode++)
    {
      for (int idof = 0; idof < my::nsd_; idof++)
      {
        // residual for normal direction
        rhs(inode * my::nsd_) -=
            my::funct_(inode) * porosity_gp * my::unitnormal_(idof) * convvel(idof) * my::fac_;
      }

      for (int nnod = 0; nnod < my::bdrynen_; nnod++)
      {
        for (int idof2 = 0; idof2 < my::nsd_; idof2++)
        {
          k_fluid(inode * my::nsd_, nnod * my::nsd_ + idof2) += my::funct_(inode) * porosity_gp *
                                                                my::unitnormal_(idof2) *
                                                                my::funct_(nnod) * my::fac_;
        }
      }
    }

  } /* end of loop over integration points gpid */

  return;
}

/*----------------------------------------------------------------------------------*
 * evaluate no penetration condition (lagrange multiplier)   (protected) vuong 11/14|
 *----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::NoPenetrationMatOD(
    DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& k_struct,
    Epetra_SerialDenseMatrix& k_lambda)
{
  switch (distype)
  {
    // 2D:
    case DRT::Element::line2:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::quad4)
      {
        NoPenetrationMatOD<DRT::Element::quad4>(
            ele, params, discretization, lm, k_struct, k_lambda);
      }
      else if (ele->ParentElement()->Shape() == DRT::Element::tri3)
      {
        NoPenetrationMatOD<DRT::Element::tri3>(ele, params, discretization, lm, k_struct, k_lambda);
      }
      else
      {
        dserror("expected combination line2/quad4 or line2/tri3 for line/parent pair");
      }
      break;
    }
    case DRT::Element::line3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::quad9)
      {
        NoPenetrationMatOD<DRT::Element::quad9>(
            ele, params, discretization, lm, k_struct, k_lambda);
      }
      else
      {
        dserror("expected combination line3/quad9 for line/parent pair");
      }
      break;
    }
    case DRT::Element::nurbs3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::nurbs9)
      {
        NoPenetrationMatOD<DRT::Element::nurbs9>(
            ele, params, discretization, lm, k_struct, k_lambda);
      }
      else
      {
        dserror("expected combination nurbs3/nurbs9 for line/parent pair");
      }
      break;
    }
    // 3D:
    case DRT::Element::quad4:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::hex8)
      {
        NoPenetrationMatOD<DRT::Element::hex8>(ele, params, discretization, lm, k_struct, k_lambda);
      }
      else
      {
        dserror("expected combination quad4/hex8 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::tri3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::tet4)
      {
        NoPenetrationMatOD<DRT::Element::tet4>(ele, params, discretization, lm, k_struct, k_lambda);
      }
      else
      {
        dserror("expected combination tri3/tet4 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::tri6:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::tet10)
      {
        NoPenetrationMatOD<DRT::Element::tet10>(
            ele, params, discretization, lm, k_struct, k_lambda);
      }
      else
      {
        dserror("expected combination tri6/tet10 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::quad9:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::hex27)
      {
        NoPenetrationMatOD<DRT::Element::hex27>(
            ele, params, discretization, lm, k_struct, k_lambda);
      }
      else
      {
        dserror("expected combination hex27/hex27 for surface/parent pair");
      }
      break;
    }
    default:
    {
      dserror("surface/parent element pair not yet implemented. Just do it.\n");
      break;
    }
  }
}

/*-----------------------------------------------------------------------------------------------*
 * evaluate no penetration condition off diagonal (lagrange multiplier)    (protected) vuong 06/11|
 *------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType pdistype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::NoPenetrationMatOD(
    DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& k_struct,
    Epetra_SerialDenseMatrix& k_lambda)
{
  // This function is only implemented for 3D
  if (my::bdrynsd_ != 2 and my::bdrynsd_ != 1)
    dserror("PressureCoupling is only implemented for 2D and 3D!");

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<my::bdrynsd_> intpoints(
      DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a my::nsd_ dimensional domain, since my::nsd_ determines the dimension of
  // FluidBoundary element!)
  GEO::fillInitialPositionArray<distype, my::nsd_, LINALG::Matrix<my::nsd_, my::bdrynen_>>(
      ele, my::xyze_);

  // get timescale parameter from parameter list (depends on time integration scheme)
  double timescale = params.get<double>("timescale", -1.0);
  if (timescale == -1.0) dserror("no timescale parameter in parameter list");

  // reset timescale in stationary case
  if (my::fldparatimint_->IsStationary()) timescale = 0.0;

  // displacements
  Teuchos::RCP<const Epetra_Vector> dispnp;
  std::vector<double> mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp, mydispnp, lm);
    }
    dsassert(mydispnp.size() != 0, "no displacement values for boundary element");

    // Add the deformation of the ALE mesh to the nodes coordinates
    for (int inode = 0; inode < my::bdrynen_; ++inode)
    {
      for (int idim = 0; idim < my::nsd_; ++idim)
      {
        my::xyze_(idim, inode) += mydispnp[my::nsd_ * inode + idim];
      }
    }
  }

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
  Teuchos::RCP<const Epetra_Vector> gridvel = discretization.GetState("gridv");

  if (velnp == Teuchos::null) dserror("Cannot get state vector 'velnp'");

  std::vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp, myvelnp, lm);

  std::vector<double> mygridvel(lm.size());
  DRT::UTILS::ExtractMyValues(*gridvel, mygridvel, lm);

  // allocate velocity vectors
  LINALG::Matrix<my::nsd_, my::bdrynen_> evelnp(true);
  LINALG::Matrix<my::nsd_, my::bdrynen_> egridvel(true);

  // split velocity and pressure, insert into element arrays
  for (int inode = 0; inode < my::bdrynen_; inode++)
    for (int idim = 0; idim < my::nsd_; idim++)
    {
      evelnp(idim, inode) = myvelnp[idim + (inode * my::nsd_)];
      egridvel(idim, inode) = mygridvel[idim + (inode * my::nsd_)];
    }

  Teuchos::RCP<const Epetra_Vector> glambda = discretization.GetState("lambda");

  if (glambda == Teuchos::null) dserror("Cannot get state vector 'lambda'");

  std::vector<double> mylambda(lm.size());
  DRT::UTILS::ExtractMyValues(*glambda, mylambda, lm);

  LINALG::Matrix<my::nsd_, my::bdrynen_> elambda(true);

  // copy lagrange multiplier values into matrix
  for (int inode = 0; inode < my::bdrynen_; inode++)
    for (int idim = 0; idim < my::nsd_; idim++)
    {
      elambda(idim, inode) = mylambda[idim + (inode * my::nsd_)];
    }

  // allocate convective velocity at node
  LINALG::Matrix<my::nsd_, my::bdrynen_> econvvel(true);

  econvvel += evelnp;
  if (not my::fldparatimint_->IsStationary()) econvvel -= egridvel;

  // --------------------------------------------------
  // parent element
  // --------------------------------------------------

  // get the parent element
  DRT::ELEMENTS::Fluid* pele = ele->ParentElement();

  /// number of parentnodes
  static const int nenparent = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

  // get element location vector and ownerships
  std::vector<int> plm;
  std::vector<int> plmowner;
  std::vector<int> plmstride;
  pele->DRT::Element::LocationVector(discretization, plm, plmowner, plmstride);

  std::vector<double> parentdispnp;
  DRT::UTILS::ExtractMyValues(*dispnp, parentdispnp, plm);

  // update element geometry of parent element
  LINALG::Matrix<my::nsd_, nenparent> xrefe;  // material coord. of parent element
  LINALG::Matrix<my::nsd_, nenparent> xcurr;  // current  coord. of parent element
  {
    DRT::Node** nodes = pele->Nodes();
    for (int i = 0; i < nenparent; ++i)
    {
      for (int j = 0; j < my::nsd_; ++j)
      {
        const double* x = nodes[i]->X();
        xrefe(j, i) = x[j];
        xcurr(j, i) = xrefe(j, i) + parentdispnp[i * my::numdofpernode_ + j];
      }
    }
  }


  std::vector<double> pvelnp(plm.size());
  DRT::UTILS::ExtractMyValues(*velnp, pvelnp, plm);

  // allocate vectors
  LINALG::Matrix<nenparent, 1> pepressnp(true);

  // split velocity and pressure, insert into element arrays
  for (int inode = 0; inode < nenparent; inode++)
  {
    pepressnp(inode) = pvelnp[my::nsd_ + (inode * my::numdofpernode_)];
  }

  // get coordinates of gauss points w.r.t. local parent coordinate system
  LINALG::SerialDenseMatrix pqxg(intpoints.IP().nquad, my::nsd_);
  LINALG::Matrix<my::nsd_, my::nsd_> derivtrafo(true);

  DRT::UTILS::BoundaryGPToParentGP<my::nsd_>(
      pqxg, derivtrafo, intpoints, pdistype, distype, ele->SurfaceNumber());

  // coordinates of gauss points of parent element
  LINALG::Matrix<my::nsd_, 1> pxsi(true);

  LINALG::Matrix<my::bdrynen_, 1> eporosity(true);

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  // --------------------------------------------------

  // In the case of nurbs the normal vector is multiplied with normalfac
  double normalfac = 0.0;
  std::vector<Epetra_SerialDenseVector> mypknots(my::nsd_);
  std::vector<Epetra_SerialDenseVector> myknots(my::bdrynsd_);
  Epetra_SerialDenseVector weights(my::bdrynen_);
  Epetra_SerialDenseVector pweights(pele->NumNode());

  // for isogeometric elements --- get knotvectors for parent
  // element and surface element, get weights
  if (IsNurbs<distype>::isnurbs)
  {
    bool zero_size = DRT::NURBS::GetKnotVectorAndWeightsForNurbsBoundaryAndParent(pele, ele,
        ele->SurfaceNumber(), discretization, mypknots, myknots, pweights, weights, normalfac);

    if (zero_size)
    {
      return;
    }
  }
  // --------------------------------------------------
  // tangent vectors
  LINALG::Matrix<my::nsd_, 1> tangent1(true);
  LINALG::Matrix<my::nsd_, 1> tangent2(true);

  // allocate convective velocity at gauss point
  LINALG::Matrix<my::nsd_, 1> convvel(true);
  LINALG::Matrix<my::nsd_, 1> lambda(true);

  //! array for dual shape functions for boundary element
  LINALG::Matrix<my::bdrynen_, 1> dualfunct(true);

  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the
    // shape function at the Gauss point Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(my::funct_, my::deriv_, my::fac_,
        my::unitnormal_, my::drs_, my::xsi_, my::xyze_, intpoints, gpid, &myknots, &weights,
        IsNurbs<distype>::isnurbs);

    // --------------------------------------------------
    // parent element
    // --------------------------------------------------

    // get shape functions and derivatives in the plane of the element
    LINALG::Matrix<nenparent, 1> pfunct(true);
    LINALG::Matrix<my::nsd_, nenparent> pderiv_loc;

    // coordinates of the current integration point
    for (int idim = 0; idim < my::nsd_; idim++) pxsi(idim) = pqxg(gpid, idim);

    // get shape functions and derivatives of the parent element
    if (not IsNurbs<distype>::isnurbs)
    {
      // shape functions and their first derivatives of parent element
      DRT::UTILS::shape_function<pdistype>(pxsi, pfunct);
      DRT::UTILS::shape_function_deriv1<pdistype>(pxsi, pderiv_loc);
    }
    // only for NURBS!!!
    else
    {
      DRT::NURBS::UTILS::nurbs_get_funct_deriv(
          pfunct, pderiv_loc, pxsi, mypknots, pweights, pdistype);
    }

    // get Jacobian matrix and determinant w.r.t. spatial configuration
    // transposed jacobian "dx/ds"
    LINALG::Matrix<my::nsd_, my::nsd_> xjm;
    LINALG::Matrix<my::nsd_, my::nsd_> Jmat;
    xjm.MultiplyNT(pderiv_loc, xcurr);
    Jmat.MultiplyNT(pderiv_loc, xrefe);
    // jacobian determinant "det(dx/ds)"
    const double det = xjm.Determinant();
    // jacobian determinant "det(dX/ds)"
    const double detJ = Jmat.Determinant();
    // jacobian determinant "det(dx/dX) = det(dx/ds)/det(dX/ds)"
    const double J = det / detJ;

    double press = pepressnp.Dot(pfunct);

    double dphi_dp = 0.0;
    double dphi_dJ = 0.0;
    double porosity_gp = 0.0;

    ComputePorosityAtGP(
        params, ele, my::funct_, eporosity, press, J, gpid, porosity_gp, dphi_dp, dphi_dJ, false);

    // --------------------------------------------------
    // --------------------------------------------------

    double Axi[3];
    for (int i = 0; i < my::bdrynsd_; i++) Axi[i] = my::xsi_(i);
    for (int i = my::bdrynsd_; i < 3; i++) Axi[i] = 0.0;
    VOLMORTAR::UTILS::dual_shape_function<distype>(dualfunct, Axi, *ele);

    // dxyzdrs vector -> normal which is not normalized
    LINALG::Matrix<my::bdrynsd_, my::nsd_> dxyzdrs(0.0);
    dxyzdrs.MultiplyNT(my::deriv_, my::xyze_);

    // in the case of nurbs the normal vector must be scaled with a special factor
    if (IsNurbs<distype>::isnurbs) my::unitnormal_.Scale(normalfac);

    convvel.Multiply(econvvel, my::funct_);
    lambda.Multiply(elambda, dualfunct);

    //  derivatives of surface normals wrt mesh displacements
    LINALG::Matrix<my::nsd_, my::bdrynen_ * my::nsd_> normalderiv(true);
    LINALG::Matrix<my::nsd_, my::bdrynen_ * my::nsd_> tangent1deriv(true);
    LINALG::Matrix<my::nsd_, my::bdrynen_ * my::nsd_> tangent2deriv(true);

    // The integration factor is not multiplied with drs
    // since it is the same as the scaling factor for the unit normal derivatives
    // Therefore it cancels out!!
    const double fac = intpoints.IP().qwgt[gpid];

    if (my::nsd_ == 3)
    {
      for (int node = 0; node < my::bdrynen_; ++node)
      {
        normalderiv(0, 3 * node) += 0.;
        normalderiv(0, 3 * node + 1) +=
            (my::deriv_(0, node) * dxyzdrs(1, 2) - my::deriv_(1, node) * dxyzdrs(0, 2));
        normalderiv(0, 3 * node + 2) +=
            (my::deriv_(1, node) * dxyzdrs(0, 1) - my::deriv_(0, node) * dxyzdrs(1, 1));

        normalderiv(1, 3 * node) +=
            (my::deriv_(1, node) * dxyzdrs(0, 2) - my::deriv_(0, node) * dxyzdrs(1, 2));
        normalderiv(1, 3 * node + 1) += 0.;
        normalderiv(1, 3 * node + 2) +=
            (my::deriv_(0, node) * dxyzdrs(1, 0) - my::deriv_(1, node) * dxyzdrs(0, 0));

        normalderiv(2, 3 * node) +=
            (my::deriv_(0, node) * dxyzdrs(1, 1) - my::deriv_(1, node) * dxyzdrs(0, 1));
        normalderiv(2, 3 * node + 1) +=
            (my::deriv_(1, node) * dxyzdrs(0, 0) - my::deriv_(0, node) * dxyzdrs(1, 0));
        normalderiv(2, 3 * node + 2) += 0.;
      }

      // in the case of nurbs the normal vector must be scaled with a special factor
      if (IsNurbs<distype>::isnurbs) normalderiv.Scale(normalfac);

      if (abs(my::unitnormal_(0)) > 1.0e-6 || abs(my::unitnormal_(1)) > 1.0e-6)
      {
        tangent1(0) = -my::unitnormal_(1);
        tangent1(1) = my::unitnormal_(0);
        tangent1(2) = 0.0;

        for (int node = 0; node < my::bdrynen_; ++node)
        {
          tangent1deriv(0, 3 * node) = -normalderiv(1, 3 * node);
          tangent1deriv(0, 3 * node + 1) = -normalderiv(1, 3 * node + 1);
          tangent1deriv(0, 3 * node + 2) = -normalderiv(1, 3 * node + 2);

          tangent1deriv(1, 3 * node) = normalderiv(0, 3 * node);
          tangent1deriv(1, 3 * node + 1) = normalderiv(0, 3 * node + 1);
          tangent1deriv(1, 3 * node + 2) = normalderiv(0, 3 * node + 2);

          tangent1deriv(2, 3 * node) = 0.;
          tangent1deriv(2, 3 * node + 1) = 0.;
          tangent1deriv(2, 3 * node + 2) = 0.;
        }
      }
      else
      {
        tangent1(0) = 0.0;
        tangent1(1) = -my::unitnormal_(2);
        tangent1(2) = my::unitnormal_(1);

        for (int node = 0; node < my::bdrynen_; ++node)
        {
          tangent1deriv(0, 3 * node) = 0.0;
          tangent1deriv(0, 3 * node + 1) = 0.0;
          tangent1deriv(0, 3 * node + 2) = 0.0;

          tangent1deriv(1, 3 * node) = -normalderiv(2, 3 * node);
          tangent1deriv(1, 3 * node + 1) = -normalderiv(2, 3 * node + 1);
          tangent1deriv(1, 3 * node + 2) = -normalderiv(2, 3 * node + 2);

          tangent1deriv(2, 3 * node) = normalderiv(1, 3 * node);
          tangent1deriv(2, 3 * node + 1) = normalderiv(1, 3 * node + 1);
          tangent1deriv(2, 3 * node + 2) = normalderiv(1, 3 * node + 2);
        }
      }

      // teta follows from corkscrew rule (teta = n x txi)
      tangent2(0) = my::unitnormal_(1) * tangent1(2) - my::unitnormal_(2) * tangent1(1);
      tangent2(1) = my::unitnormal_(2) * tangent1(0) - my::unitnormal_(0) * tangent1(2);
      tangent2(2) = my::unitnormal_(0) * tangent1(1) - my::unitnormal_(1) * tangent1(0);

      for (int node = 0; node < my::bdrynen_; ++node)
      {
        for (int idim = 0; idim < 3; ++idim)
        {
          tangent2deriv(0, 3 * node + idim) =
              normalderiv(1, 3 * node + idim) * tangent1(2) +
              my::unitnormal_(1) * tangent1deriv(2, 3 * node + idim) -
              normalderiv(2, 3 * node + idim) * tangent1(1) -
              my::unitnormal_(2) * tangent1deriv(1, 3 * node + idim);

          tangent2deriv(1, 3 * node + idim) =
              normalderiv(2, 3 * node + idim) * tangent1(0) +
              my::unitnormal_(2) * tangent1deriv(0, 3 * node + idim) -
              normalderiv(0, 3 * node + idim) * tangent1(2) -
              my::unitnormal_(0) * tangent1deriv(2, 3 * node + idim);

          tangent2deriv(2, 3 * node + idim) =
              normalderiv(0, 3 * node + idim) * tangent1(1) +
              my::unitnormal_(0) * tangent1deriv(1, 3 * node + idim) -
              normalderiv(1, 3 * node + idim) * tangent1(0) -
              my::unitnormal_(1) * tangent1deriv(0, 3 * node + idim);
        }
      }
    }
    else if (my::nsd_ == 2)
    {
      for (int node = 0; node < my::bdrynen_; ++node)
      {
        normalderiv(0, my::nsd_ * node) += 0.;
        normalderiv(0, my::nsd_ * node + 1) += my::deriv_(0, node);

        normalderiv(1, my::nsd_ * node) += -my::deriv_(0, node);
        normalderiv(1, my::nsd_ * node + 1) += 0.;
      }

      // in the case of nurbs the normal vector must be scaled with a special factor
      if (IsNurbs<distype>::isnurbs) normalderiv.Scale(normalfac);

      // simple definition for txi
      tangent1(0) = -my::unitnormal_(1);
      tangent1(1) = my::unitnormal_(0);

      for (int node = 0; node < my::bdrynen_; ++node)
      {
        tangent1deriv(0, my::nsd_ * node) = -normalderiv(1, my::nsd_ * node);
        tangent1deriv(0, my::nsd_ * node + 1) = -normalderiv(1, my::nsd_ * node + 1);
        ;

        tangent1deriv(1, my::nsd_ * node) = normalderiv(0, my::nsd_ * node);
        tangent1deriv(1, my::nsd_ * node + 1) = normalderiv(0, my::nsd_ * node + 1);
      }
    }

    static LINALG::Matrix<1, my::bdrynen_ * my::nsd_> convvel_normalderiv(true);
    convvel_normalderiv.MultiplyTN(convvel, normalderiv);

    // fill element matrix
    for (int inode = 0; inode < my::bdrynen_; inode++)
    {
      const double funct_fac = my::funct_(inode) * porosity_gp * fac;
      for (int nnod = 0; nnod < my::bdrynen_; nnod++)
      {
        for (int idof = 0; idof < my::nsd_; idof++)
        {
          k_struct(inode * my::nsd_, nnod * my::nsd_ + idof) +=
              -my::unitnormal_(idof) * timescale * my::funct_(nnod) * my::funct_(inode) *
                  porosity_gp * my::fac_ +
              convvel_normalderiv(0, nnod * my::nsd_ + idof) * funct_fac;
        }
      }
    }

    if (my::nsd_ == 3)
    {
      static LINALG::Matrix<1, my::bdrynen_ * my::nsd_> lambda_tangent1deriv(true);
      lambda_tangent1deriv.MultiplyTN(lambda, tangent1deriv);
      static LINALG::Matrix<1, my::bdrynen_ * my::nsd_> lambda_tangent2deriv(true);
      lambda_tangent2deriv.MultiplyTN(lambda, tangent2deriv);

      for (int inode = 0; inode < my::bdrynen_; inode++)
      {
        const double funct_fac = my::funct_(inode) * fac;
        for (int nnod = 0; nnod < my::bdrynen_; nnod++)
        {
          for (int idof = 0; idof < my::nsd_; idof++)
          {
            k_struct(inode * my::nsd_ + 1, nnod * my::nsd_ + idof) +=
                lambda_tangent1deriv(0, nnod * my::nsd_ + idof) * funct_fac;
            k_struct(inode * my::nsd_ + 2, nnod * my::nsd_ + idof) +=
                lambda_tangent2deriv(0, nnod * my::nsd_ + idof) * funct_fac;
          }
        }
      }
    }
    else if (my::nsd_ == 2)
    {
      LINALG::Matrix<1, my::bdrynen_ * my::nsd_> lambda_tangent1deriv(true);
      lambda_tangent1deriv.MultiplyTN(lambda, tangent1deriv);

      for (int inode = 0; inode < my::bdrynen_; inode++)
      {
        const double funct_fac = my::funct_(inode) * fac;
        for (int nnod = 0; nnod < my::bdrynen_; nnod++)
        {
          for (int idof = 0; idof < my::nsd_; idof++)
          {
            k_struct(inode * my::nsd_ + 1, nnod * my::nsd_ + idof) +=
                lambda_tangent1deriv(0, nnod * my::nsd_ + idof) * funct_fac;
          }
        }
      }
    }

    if (my::nsd_ == 3)
    {
      for (int inode = 0; inode < my::bdrynen_; inode++)
      {
        const double funct_fac = my::funct_(inode) * my::fac_;
        for (int nnod = 0; nnod < my::bdrynen_; nnod++)
        {
          for (int idof = 0; idof < my::nsd_; idof++)
          {
            k_lambda(inode * my::nsd_ + 1, nnod * my::nsd_ + idof) +=
                tangent1(idof) * dualfunct(nnod) * funct_fac;
            k_lambda(inode * my::nsd_ + 2, nnod * my::nsd_ + idof) +=
                tangent2(idof) * dualfunct(nnod) * funct_fac;
          }
        }
      }
    }
    else if (my::nsd_ == 2)
    {
      for (int inode = 0; inode < my::bdrynen_; inode++)
      {
        const double funct_fac = my::funct_(inode) * my::fac_;
        for (int nnod = 0; nnod < my::bdrynen_; nnod++)
        {
          for (int idof = 0; idof < my::nsd_; idof++)
          {
            k_lambda(inode * my::nsd_ + 1, nnod * my::nsd_ + idof) +=
                tangent1(idof) * dualfunct(nnod) * funct_fac;
          }
        }
      }
    }
  } /* end of loop over integration points gpid */

  return;
}


/*----------------------------------------------------------------------------------*
 * evaluate no penetration condition (lagrange multiplier)   (protected) vuong 11/14|
 *----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::NoPenetrationMatODPoroPres(
    DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& k_pres)
{
  switch (distype)
  {
    // 2D:
    case DRT::Element::line2:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::quad4)
      {
        NoPenetrationMatODPoroPres<DRT::Element::quad4>(ele, params, discretization, lm, k_pres);
      }
      else if (ele->ParentElement()->Shape() == DRT::Element::tri3)
      {
        NoPenetrationMatODPoroPres<DRT::Element::tri3>(ele, params, discretization, lm, k_pres);
      }
      else
      {
        dserror("expected combination line2/quad4 or line2/tri3 for line/parent pair");
      }
      break;
    }
    case DRT::Element::line3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::quad9)
      {
        NoPenetrationMatODPoroPres<DRT::Element::quad9>(ele, params, discretization, lm, k_pres);
      }
      else
      {
        dserror("expected combination line3/quad9 for line/parent pair");
      }
      break;
    }
    case DRT::Element::nurbs3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::nurbs9)
      {
        NoPenetrationMatODPoroPres<DRT::Element::nurbs9>(ele, params, discretization, lm, k_pres);
      }
      else
      {
        dserror("expected combination nurbs3/nurbs9 for line/parent pair");
      }
      break;
    }
    // 3D:
    case DRT::Element::quad4:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::hex8)
      {
        NoPenetrationMatODPoroPres<DRT::Element::hex8>(ele, params, discretization, lm, k_pres);
      }
      else
      {
        dserror("expected combination quad4/hex8 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::tri3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::tet4)
      {
        NoPenetrationMatODPoroPres<DRT::Element::tet4>(ele, params, discretization, lm, k_pres);
      }
      else
      {
        dserror("expected combination tri3/tet4 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::tri6:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::tet10)
      {
        NoPenetrationMatODPoroPres<DRT::Element::tet10>(ele, params, discretization, lm, k_pres);
      }
      else
      {
        dserror("expected combination tri6/tet10 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::quad9:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::hex27)
      {
        NoPenetrationMatODPoroPres<DRT::Element::hex27>(ele, params, discretization, lm, k_pres);
      }
      else
      {
        dserror("expected combination hex27/hex27 for surface/parent pair");
      }
      break;
    }
    default:
    {
      dserror("surface/parent element pair not yet implemented. Just do it.\n");
      break;
    }
  }
}

/*-----------------------------------------------------------------------------------------------*
 * evaluate no penetration condition off diagonal (lagrange multiplier)    (protected) vuong 06/11|
 *------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType pdistype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::NoPenetrationMatODPoroPres(
    DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& k_pres)
{
  // This function is only implemented for 3D
  if (my::bdrynsd_ != 2 and my::bdrynsd_ != 1)
    dserror("PressureCoupling is only implemented for 2D and 3D!");

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<my::bdrynsd_> intpoints(
      DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a my::nsd_ dimensional domain, since my::nsd_ determines the dimension of
  // FluidBoundary element!)
  GEO::fillInitialPositionArray<distype, my::nsd_, LINALG::Matrix<my::nsd_, my::bdrynen_>>(
      ele, my::xyze_);

  // displacements
  Teuchos::RCP<const Epetra_Vector> dispnp;
  std::vector<double> mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp, mydispnp, lm);
    }
    dsassert(mydispnp.size() != 0, "no displacement values for boundary element");

    // Add the deformation of the ALE mesh to the nodes coordinates
    for (int inode = 0; inode < my::bdrynen_; ++inode)
    {
      for (int idim = 0; idim < my::nsd_; ++idim)
      {
        my::xyze_(idim, inode) += mydispnp[my::nsd_ * inode + idim];
      }
    }
  }

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
  Teuchos::RCP<const Epetra_Vector> gridvel = discretization.GetState("gridv");

  if (velnp == Teuchos::null) dserror("Cannot get state vector 'velnp'");

  std::vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp, myvelnp, lm);

  std::vector<double> mygridvel(lm.size());
  DRT::UTILS::ExtractMyValues(*gridvel, mygridvel, lm);

  // allocate velocity vectors
  LINALG::Matrix<my::nsd_, my::bdrynen_> evelnp(true);
  LINALG::Matrix<my::nsd_, my::bdrynen_> egridvel(true);

  // split velocity and pressure, insert into element arrays
  for (int inode = 0; inode < my::bdrynen_; inode++)
    for (int idim = 0; idim < my::nsd_; idim++)
    {
      evelnp(idim, inode) = myvelnp[idim + (inode * my::numdofpernode_)];
      egridvel(idim, inode) = mygridvel[idim + (inode * my::numdofpernode_)];
    }

  // allocate convective velocity at node
  LINALG::Matrix<my::nsd_, my::bdrynen_> econvvel(true);

  econvvel += evelnp;
  if (not my::fldparatimint_->IsStationary()) econvvel -= egridvel;

  // --------------------------------------------------
  // parent element
  // --------------------------------------------------

  // get the parent element
  DRT::ELEMENTS::Fluid* pele = ele->ParentElement();

  /// number of parentnodes
  static const int nenparent = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

  // get element location vector and ownerships
  std::vector<int> plm;
  std::vector<int> plmowner;
  std::vector<int> plmstride;
  pele->DRT::Element::LocationVector(discretization, plm, plmowner, plmstride);

  std::vector<double> parentdispnp;
  DRT::UTILS::ExtractMyValues(*dispnp, parentdispnp, plm);

  // update element geometry of parent element
  LINALG::Matrix<my::nsd_, nenparent> xrefe;  // material coord. of parent element
  LINALG::Matrix<my::nsd_, nenparent> xcurr;  // current  coord. of parent element
  {
    DRT::Node** nodes = pele->Nodes();
    for (int i = 0; i < nenparent; ++i)
    {
      for (int j = 0; j < my::nsd_; ++j)
      {
        const double* x = nodes[i]->X();
        xrefe(j, i) = x[j];
        xcurr(j, i) = xrefe(j, i) + parentdispnp[i * my::numdofpernode_ + j];
      }
    }
  }

  std::vector<double> pvelnp(plm.size());
  DRT::UTILS::ExtractMyValues(*velnp, pvelnp, plm);

  // allocate vectors
  LINALG::Matrix<nenparent, 1> pepressnp(true);

  // split velocity and pressure, insert into element arrays
  for (int inode = 0; inode < nenparent; inode++)
  {
    pepressnp(inode) = pvelnp[my::nsd_ + (inode * my::numdofpernode_)];
  }

  // get coordinates of gauss points w.r.t. local parent coordinate system
  LINALG::SerialDenseMatrix pqxg(intpoints.IP().nquad, my::nsd_);
  LINALG::Matrix<my::nsd_, my::nsd_> derivtrafo(true);

  DRT::UTILS::BoundaryGPToParentGP<my::nsd_>(
      pqxg, derivtrafo, intpoints, pdistype, distype, ele->SurfaceNumber());

  // coordinates of gauss points of parent element
  LINALG::Matrix<my::nsd_, 1> pxsi(true);

  LINALG::Matrix<my::bdrynen_, 1> eporosity(true);

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  // --------------------------------------------------

  // In the case of nurbs the normal vector is multiplied with normalfac
  double normalfac = 0.0;
  std::vector<Epetra_SerialDenseVector> mypknots(my::nsd_);
  std::vector<Epetra_SerialDenseVector> myknots(my::bdrynsd_);
  Epetra_SerialDenseVector weights(my::bdrynen_);
  Epetra_SerialDenseVector pweights(pele->NumNode());

  // for isogeometric elements --- get knotvectors for parent
  // element and surface element, get weights
  if (IsNurbs<distype>::isnurbs)
  {
    bool zero_size = DRT::NURBS::GetKnotVectorAndWeightsForNurbsBoundaryAndParent(pele, ele,
        ele->SurfaceNumber(), discretization, mypknots, myknots, pweights, weights, normalfac);

    if (zero_size)
    {
      return;
    }
  }

  // --------------------------------------------------
  // structure velocity at gausspoint
  // LINALG::Matrix<my::nsd_,1> gridvelint;
  LINALG::Matrix<my::nsd_, 1> convvel(true);

  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the
    // shape function at the Gauss point Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(my::funct_, my::deriv_, my::fac_,
        my::unitnormal_, my::drs_, my::xsi_, my::xyze_, intpoints, gpid, &myknots, &weights,
        IsNurbs<distype>::isnurbs);

    // --------------------------------------------------
    // parent element
    // --------------------------------------------------

    // get shape functions and derivatives in the plane of the element
    LINALG::Matrix<nenparent, 1> pfunct(true);
    LINALG::Matrix<my::nsd_, nenparent> pderiv_loc;

    // coordinates of the current integration point
    for (int idim = 0; idim < my::nsd_; idim++) pxsi(idim) = pqxg(gpid, idim);

    // get shape functions and derivatives of the parent element
    if (not IsNurbs<distype>::isnurbs)
    {
      // shape functions and their first derivatives of parent element
      DRT::UTILS::shape_function<pdistype>(pxsi, pfunct);
      DRT::UTILS::shape_function_deriv1<pdistype>(pxsi, pderiv_loc);
    }
    // only for NURBS!!!
    else
    {
      DRT::NURBS::UTILS::nurbs_get_funct_deriv(
          pfunct, pderiv_loc, pxsi, mypknots, pweights, pdistype);
    }

    // get Jacobian matrix and determinant w.r.t. spatial configuration
    // transposed jacobian "dx/ds"
    LINALG::Matrix<my::nsd_, my::nsd_> xjm;
    LINALG::Matrix<my::nsd_, my::nsd_> Jmat;
    xjm.MultiplyNT(pderiv_loc, xcurr);
    Jmat.MultiplyNT(pderiv_loc, xrefe);
    // jacobian determinant "det(dx/ds)"
    const double det = xjm.Determinant();
    // jacobian determinant "det(dX/ds)"
    const double detJ = Jmat.Determinant();
    // jacobian determinant "det(dx/dX) = det(dx/ds)/det(dX/ds)"
    const double J = det / detJ;

    double press = pepressnp.Dot(pfunct);

    double dphi_dp = 0.0;
    double dphi_dJ = 0.0;
    double porosity_gp = 0.0;

    // --------------------------------------------------
    // --------------------------------------------------

    ComputePorosityAtGP(
        params, ele, my::funct_, eporosity, press, J, gpid, porosity_gp, dphi_dp, dphi_dJ, false);

    // dxyzdrs vector -> normal which is not normalized
    LINALG::Matrix<my::bdrynsd_, my::nsd_> dxyzdrs(0.0);
    dxyzdrs.MultiplyNT(my::deriv_, my::xyze_);

    // in the case of nurbs the normal vector must be scaled with a special factor
    if (IsNurbs<distype>::isnurbs) my::unitnormal_.Scale(normalfac);

    convvel.Multiply(econvvel, my::funct_);
    double normal_convel = 0.0;

    for (int idof = 0; idof < my::nsd_; idof++)
      normal_convel += my::unitnormal_(idof) * convvel(idof);

    // fill element matrix
    for (int inode = 0; inode < my::bdrynen_; inode++)
    {
      const double funct_fac = my::funct_(inode) * my::fac_;
      for (int nnod = 0; nnod < my::bdrynen_; nnod++)
      {
        k_pres(inode * my::numdofpernode_, nnod * my::numdofpernode_ + my::nsd_) +=
            +normal_convel * dphi_dp * my::funct_(nnod) * funct_fac;
      }
    }

  }  // end of loop over integration points gpid

  return;
}

/*----------------------------------------------------------------------------------*
 * evaluate no penetration condition (lagrange multiplier)   (protected) vuong 11/14|
 *----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::NoPenetrationMatODPoroDisp(
    DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& plm, Epetra_SerialDenseMatrix& k_disp)
{
  switch (distype)
  {
    // 2D:
    case DRT::Element::line2:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::quad4)
      {
        NoPenetrationMatODPoroDisp<DRT::Element::quad4>(ele, params, discretization, plm, k_disp);
      }
      else if (ele->ParentElement()->Shape() == DRT::Element::tri3)
      {
        NoPenetrationMatODPoroDisp<DRT::Element::tri3>(ele, params, discretization, plm, k_disp);
      }
      else
      {
        dserror("expected combination line2/quad4 or line2/tri3 for line/parent pair");
      }
      break;
    }
    case DRT::Element::line3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::quad9)
      {
        NoPenetrationMatODPoroDisp<DRT::Element::quad9>(ele, params, discretization, plm, k_disp);
      }
      else
      {
        dserror("expected combination line3/quad9 for line/parent pair");
      }
      break;
    }
    case DRT::Element::nurbs3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::nurbs9)
      {
        NoPenetrationMatODPoroDisp<DRT::Element::nurbs9>(ele, params, discretization, plm, k_disp);
      }
      else
      {
        dserror("expected combination nurbs3/nurbs9 for line/parent pair");
      }
      break;
    }
    // 3D:
    case DRT::Element::quad4:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::hex8)
      {
        NoPenetrationMatODPoroDisp<DRT::Element::hex8>(ele, params, discretization, plm, k_disp);
      }
      else
      {
        dserror("expected combination quad4/hex8 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::tri3:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::tet4)
      {
        NoPenetrationMatODPoroDisp<DRT::Element::tet4>(ele, params, discretization, plm, k_disp);
      }
      else
      {
        dserror("expected combination tri3/tet4 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::tri6:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::tet10)
      {
        NoPenetrationMatODPoroDisp<DRT::Element::tet10>(ele, params, discretization, plm, k_disp);
      }
      else
      {
        dserror("expected combination tri6/tet10 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::quad9:
    {
      if (ele->ParentElement()->Shape() == DRT::Element::hex27)
      {
        NoPenetrationMatODPoroDisp<DRT::Element::hex27>(ele, params, discretization, plm, k_disp);
      }
      else
      {
        dserror("expected combination hex27/hex27 for surface/parent pair");
      }
      break;
    }
    default:
    {
      dserror("surface/parent element pair not yet implemented. Just do it.\n");
      break;
    }
  }
}

/*-----------------------------------------------------------------------------------------------*
 * evaluate no penetration condition off diagonal (lagrange multiplier)    (protected) vuong 06/11|
 *------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType pdistype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::NoPenetrationMatODPoroDisp(
    DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& plm, Epetra_SerialDenseMatrix& k_disp)
{
  // This function is only implemented for 3D
  if (my::bdrynsd_ != 2 and my::bdrynsd_ != 1)
    dserror("PressureCoupling is only implemented for 2D and 3D!");

  // get element location vector and ownerships
  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;
  ele->DRT::Element::LocationVector(discretization, lm, lmowner, lmstride);

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<my::bdrynsd_> intpoints(
      DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a my::nsd_ dimensional domain, since my::nsd_ determines the dimension of
  // FluidBoundary element!)
  GEO::fillInitialPositionArray<distype, my::nsd_, LINALG::Matrix<my::nsd_, my::bdrynen_>>(
      ele, my::xyze_);

  // displacements
  Teuchos::RCP<const Epetra_Vector> dispnp;
  std::vector<double> mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp, mydispnp, lm);
    }
    dsassert(mydispnp.size() != 0, "no displacement values for boundary element");

    // Add the deformation of the ALE mesh to the nodes coordinates
    for (int inode = 0; inode < my::bdrynen_; ++inode)
    {
      for (int idim = 0; idim < my::nsd_; ++idim)
      {
        my::xyze_(idim, inode) += mydispnp[my::numdofpernode_ * inode + idim];
      }
    }
  }

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
  Teuchos::RCP<const Epetra_Vector> gridvel = discretization.GetState("gridv");

  if (velnp == Teuchos::null) dserror("Cannot get state vector 'velnp'");

  std::vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp, myvelnp, lm);

  std::vector<double> mygridvel(lm.size());
  DRT::UTILS::ExtractMyValues(*gridvel, mygridvel, lm);

  // allocate velocity vectors
  LINALG::Matrix<my::nsd_, my::bdrynen_> evelnp(true);
  LINALG::Matrix<my::nsd_, my::bdrynen_> egridvel(true);

  // split velocity and pressure, insert into element arrays
  for (int inode = 0; inode < my::bdrynen_; inode++)
    for (int idim = 0; idim < my::nsd_; idim++)
    {
      evelnp(idim, inode) = myvelnp[idim + (inode * my::numdofpernode_)];
      egridvel(idim, inode) = mygridvel[idim + (inode * my::numdofpernode_)];
    }

  // allocate convective velocity at node
  LINALG::Matrix<my::nsd_, my::bdrynen_> econvvel(true);

  econvvel += evelnp;
  if (not my::fldparatimint_->IsStationary()) econvvel -= egridvel;

  // --------------------------------------------------
  // parent element
  // --------------------------------------------------

  // get the parent element
  DRT::ELEMENTS::Fluid* pele = ele->ParentElement();

  /// number of parentnodes
  static const int nenparent = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

  std::vector<double> parentdispnp;
  DRT::UTILS::ExtractMyValues(*dispnp, parentdispnp, plm);

  // update element geometry of parent element
  LINALG::Matrix<my::nsd_, nenparent> xrefe;  // material coord. of parent element
  LINALG::Matrix<my::nsd_, nenparent> xcurr;  // current  coord. of parent element
  {
    DRT::Node** nodes = pele->Nodes();
    for (int i = 0; i < nenparent; ++i)
    {
      for (int j = 0; j < my::nsd_; ++j)
      {
        const double* x = nodes[i]->X();
        xrefe(j, i) = x[j];
        xcurr(j, i) = xrefe(j, i) + parentdispnp[i * my::numdofpernode_ + j];
      }
    }
  }

  std::vector<double> pvelnp(plm.size());
  DRT::UTILS::ExtractMyValues(*velnp, pvelnp, plm);

  // allocate vectors
  LINALG::Matrix<nenparent, 1> pepressnp(true);

  // split velocity and pressure, insert into element arrays
  for (int inode = 0; inode < nenparent; inode++)
  {
    pepressnp(inode) = pvelnp[my::nsd_ + (inode * my::numdofpernode_)];
  }

  // get coordinates of gauss points w.r.t. local parent coordinate system
  LINALG::SerialDenseMatrix pqxg(intpoints.IP().nquad, my::nsd_);
  LINALG::Matrix<my::nsd_, my::nsd_> derivtrafo(true);

  DRT::UTILS::BoundaryGPToParentGP<my::nsd_>(
      pqxg, derivtrafo, intpoints, pdistype, distype, ele->SurfaceNumber());

  // coordinates of gauss points of parent element
  LINALG::Matrix<my::nsd_, 1> pxsi(true);

  LINALG::Matrix<my::bdrynen_, 1> eporosity(true);
  LINALG::Matrix<my::nsd_, 1> convvel(true);

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  // --------------------------------------------------

  // In the case of nurbs the normal vector is multiplied with normalfac
  double normalfac = 0.0;
  std::vector<Epetra_SerialDenseVector> mypknots(my::nsd_);
  std::vector<Epetra_SerialDenseVector> myknots(my::bdrynsd_);
  Epetra_SerialDenseVector weights(my::bdrynen_);
  Epetra_SerialDenseVector pweights(pele->NumNode());

  // for isogeometric elements --- get knotvectors for parent
  // element and surface element, get weights
  if (IsNurbs<distype>::isnurbs)
  {
    bool zero_size = DRT::NURBS::GetKnotVectorAndWeightsForNurbsBoundaryAndParent(pele, ele,
        ele->SurfaceNumber(), discretization, mypknots, myknots, pweights, weights, normalfac);

    if (zero_size)
    {
      return;
    }
  }

  // --------------------------------------------------
  // structure velocity at gausspoint
  // LINALG::Matrix<my::nsd_,1> gridvelint;

  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the
    // shape function at the Gauss point Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(my::funct_, my::deriv_, my::fac_,
        my::unitnormal_, my::drs_, my::xsi_, my::xyze_, intpoints, gpid, &myknots, &weights,
        IsNurbs<distype>::isnurbs);

    // --------------------------------------------------
    // parent element
    // --------------------------------------------------

    // get shape functions and derivatives in the plane of the element
    LINALG::Matrix<nenparent, 1> pfunct(true);
    LINALG::Matrix<my::nsd_, nenparent> pderiv_loc;

    // coordinates of the current integration point
    for (int idim = 0; idim < my::nsd_; idim++) pxsi(idim) = pqxg(gpid, idim);

    // get shape functions and derivatives of the parent element
    if (not IsNurbs<distype>::isnurbs)
    {
      // shape functions and their first derivatives of parent element
      DRT::UTILS::shape_function<pdistype>(pxsi, pfunct);
      DRT::UTILS::shape_function_deriv1<pdistype>(pxsi, pderiv_loc);
    }
    // only for NURBS!!!
    else
    {
      DRT::NURBS::UTILS::nurbs_get_funct_deriv(
          pfunct, pderiv_loc, pxsi, mypknots, pweights, pdistype);
    }

    // get Jacobian matrix and determinant w.r.t. spatial configuration
    // transposed jacobian "dx/ds"
    LINALG::Matrix<my::nsd_, my::nsd_> xjm;
    LINALG::Matrix<my::nsd_, my::nsd_> Jmat;
    xjm.MultiplyNT(pderiv_loc, xcurr);
    Jmat.MultiplyNT(pderiv_loc, xrefe);
    // jacobian determinant "det(dx/ds)"
    const double det = xjm.Determinant();
    // jacobian determinant "det(dX/ds)"
    const double detJ = Jmat.Determinant();
    // jacobian determinant "det(dx/dX) = det(dx/ds)/det(dX/ds)"
    const double J = det / detJ;

    double press = pepressnp.Dot(pfunct);

    //------------------------------------------------dJ/dus = dJ/dF : dF/dus = J * F^-T . N_X = J *
    // N_x
    LINALG::Matrix<1, my::nsd_ * nenparent> dJ_dus;
    // global derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<my::nsd_, nenparent> derxy;
    // inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<my::nsd_, my::nsd_> xji;

    xji.Invert(xjm);
    derxy.Multiply(xji, pderiv_loc);

    for (int i = 0; i < nenparent; i++)
      for (int j = 0; j < my::nsd_; j++) dJ_dus(j + i * my::nsd_) = J * derxy(j, i);

    // --------------------------------------------------
    // --------------------------------------------------

    double dphi_dp = 0.0;
    double dphi_dJ = 0.0;
    double porosity_gp = 0.0;


    ComputePorosityAtGP(
        params, ele, my::funct_, eporosity, press, J, gpid, porosity_gp, dphi_dp, dphi_dJ, false);

    // dxyzdrs vector -> normal which is not normalized
    LINALG::Matrix<my::bdrynsd_, my::nsd_> dxyzdrs(0.0);
    dxyzdrs.MultiplyNT(my::deriv_, my::xyze_);

    // in the case of nurbs the normal vector must be scaled with a special factor
    if (IsNurbs<distype>::isnurbs) my::unitnormal_.Scale(normalfac);

    convvel.Multiply(econvvel, my::funct_);
    double normal_convel = 0.0;

    for (int idof = 0; idof < my::nsd_; idof++)
      normal_convel += my::unitnormal_(idof) * convvel(idof);

    // fill element matrix
    for (int inode = 0; inode < nenparent; inode++)
    {
      const double funct_fac = pfunct(inode) * my::fac_;
      for (int nnod = 0; nnod < nenparent; nnod++)
      {
        for (int idof = 0; idof < my::nsd_; idof++)
        {
          k_disp(inode * my::numdofpernode_, nnod * my::nsd_ + idof) +=
              +normal_convel * dphi_dJ * dJ_dus(nnod * my::nsd_ + idof) * funct_fac;
        }
      }
    }

  }  // end of loop over integration points gpid

  return;
}


/*----------------------------------------------------------------------*
 * create/delete instance (public)                            vuong 06/11||
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleBoundaryCalcPoroP1<distype>*
DRT::ELEMENTS::FluidEleBoundaryCalcPoroP1<distype>::Instance(bool create)
{
  static FluidEleBoundaryCalcPoroP1<distype>* instance;
  if (create)
  {
    if (instance == NULL)
    {
      instance = new FluidEleBoundaryCalcPoroP1<distype>();
    }
  }
  else
  {
    if (instance != NULL) delete instance;
    instance = NULL;
  }
  return instance;
}

/*----------------------------------------------------------------------*
 * called upon destruction (public)                            vuong 06/11||
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoroP1<distype>::Done()
{
  // delete ele1 pointer! Afterwards we have to go! But since ele1 is a
  // cleanup call, we can do it ele1 way.
  Instance(false);
}

/*----------------------------------------------------------------------*
 * get porosity at nodes (protected)                        vuong 06/11||
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::FluidEleBoundaryCalcPoroP1<distype>::ComputeNodalPorosity(
    DRT::ELEMENTS::FluidBoundary* ele, const std::vector<double>& mydispnp,
    LINALG::Matrix<my::bdrynen_, 1>& eporosity)
{
  for (int inode = 0; inode < my::bdrynen_; inode++)
    eporosity(inode) = mydispnp[my::nsd_ + (inode * my::numdofpernode_)];

  return true;
}

/*----------------------------------------------------------------------*
 * compute porosity at gauss point (protected)                vuong 06/11||
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleBoundaryCalcPoroP1<distype>::ComputePorosityAtGP(
    Teuchos::ParameterList& params, DRT::ELEMENTS::FluidBoundary* ele,
    const LINALG::Matrix<my::bdrynen_, 1>& funct, const LINALG::Matrix<my::bdrynen_, 1>& eporosity,
    double press, double J, int gp, double& porosity, double& dphi_dp, double& dphi_dJ, bool save)
{
  porosity = eporosity.Dot(my::funct_);
  dphi_dp = 0.0;
  dphi_dJ = 0.0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoro<DRT::Element::quad4>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoro<DRT::Element::quad8>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoro<DRT::Element::quad9>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoro<DRT::Element::tri3>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoro<DRT::Element::tri6>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoro<DRT::Element::line2>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoro<DRT::Element::line3>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoro<DRT::Element::nurbs2>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoro<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoro<DRT::Element::nurbs4>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoro<DRT::Element::nurbs9>;

template class DRT::ELEMENTS::FluidEleBoundaryCalcPoroP1<DRT::Element::quad4>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoroP1<DRT::Element::quad8>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoroP1<DRT::Element::quad9>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoroP1<DRT::Element::tri3>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoroP1<DRT::Element::tri6>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoroP1<DRT::Element::line2>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoroP1<DRT::Element::line3>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoroP1<DRT::Element::nurbs2>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoroP1<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoroP1<DRT::Element::nurbs4>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcPoroP1<DRT::Element::nurbs9>;
