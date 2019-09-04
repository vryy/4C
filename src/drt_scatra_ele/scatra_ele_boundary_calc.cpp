/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra boundary terms at integration points

\level 1

\maintainer Anh-Tu Vuong
 */
/*----------------------------------------------------------------------*/
#include "scatra_ele_boundary_calc.H"
#include "scatra_ele.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"

#include <cstdlib>

#include "../drt_fem_general/drt_utils_boundary_integration.H"

#include "../drt_inpar/inpar_s2i.H"

#include "../drt_lib/drt_globalproblem.H"  // for curves and functions
#include "../drt_lib/standardtypes_cpp.H"  // for EPS12 and so on

#include "../drt_mat/fourieriso.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/scatra_mat.H"
#include "../drt_mat/thermostvenantkirchhoff.H"

#include "../drt_nurbs_discret/drt_nurbs_utils.H"

#include "../drt_fluid/fluid_rotsym_periodicbc.H"
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::ScaTraEleBoundaryCalc(
    const int numdofpernode, const int numscal, const std::string& disname)
    : scatraparamstimint_(DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance(
          disname)),  // params for time integration
      scatraparams_(DRT::ELEMENTS::ScaTraEleParameterStd::Instance(disname)),
      numdofpernode_(numdofpernode),
      numscal_(numscal),
      xyze_(true),  // initialize to zero
      weights_(true),
      myknots_(nsd_),
      mypknots_(nsd_ + 1),
      normalfac_(1.0),
      ephinp_(numdofpernode_, LINALG::Matrix<nen_, 1>(true)),
      edispnp_(true),
      diffus_(numscal_, 0),
      // valence_(numscal_,0),
      shcacp_(0.0),
      xsi_(true),
      funct_(true),
      deriv_(true),
      derxy_(true),
      normal_(true),
      velint_(true),
      metrictensor_(true),
      rotsymmpbc_(Teuchos::rcp(new FLD::RotationallySymmetricPeriodicBC<distype, nsd_ + 2,
          DRT::ELEMENTS::Fluid::none>()))
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::SetupCalc(
    DRT::FaceElement* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization)
{
  // get node coordinates (we have a nsd_+1 dimensional domain!)
  GEO::fillInitialPositionArray<distype, nsd_ + 1, LINALG::Matrix<nsd_ + 1, nen_>>(ele, xyze_);

  // Now do the nurbs specific stuff (for isogeometric elements)
  if (DRT::NURBS::IsNurbs(distype))
  {
    // for isogeometric elements --- get knotvectors for parent
    // element and boundary element, get weights
    bool zero_size =
        DRT::NURBS::GetKnotVectorAndWeightsForNurbsBoundary(ele, ele->FaceParentNumber(),
            ele->ParentElement()->Id(), discretization, mypknots_, myknots_, weights_, normalfac_);

    // if we have a zero sized element due to a interpolated point -> exit here
    if (zero_size) return -1;
  }  // Nurbs specific stuff

  // rotationally symmetric periodic bc's: do setup for current element
  rotsymmpbc_->Setup(ele);

  return 0;
}

/*----------------------------------------------------------------------*
 | evaluate element                                          fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::Evaluate(DRT::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra)
{
  //--------------------------------------------------------------------------------
  // preparations for element
  //--------------------------------------------------------------------------------
  if (SetupCalc(ele, params, discretization) == -1) return 0;

  //--------------------------------------------------------------------------------
  // extract element based or nodal values
  //--------------------------------------------------------------------------------
  ExtractElementAndNodeValues(ele, params, discretization, la);

  // check for the action parameter
  const SCATRA::BoundaryAction action = DRT::INPUT::get<SCATRA::BoundaryAction>(params, "action");
  // evaluate action
  EvaluateAction(ele, params, discretization, action, la, elemat1_epetra, elemat2_epetra,
      elevec1_epetra, elevec2_epetra, elevec3_epetra);

  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::ExtractElementAndNodeValues(
    DRT::FaceElement* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la)
{
  // get additional state vector for ALE case: grid displacement
  if (scatraparams_->IsAle())
  {
    // get number of dof-set associated with displacement related dofs
    const int ndsdisp = params.get<int>("ndsdisp");

    Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState(ndsdisp, "dispnp");
    if (dispnp == Teuchos::null) dserror("Cannot get state vector 'dispnp'");

    // determine number of displacement related dofs per node
    const int numdispdofpernode = la[ndsdisp].lm_.size() / nen_;

    // construct location vector for displacement related dofs
    std::vector<int> lmdisp((nsd_ + 1) * nen_, -1);
    for (int inode = 0; inode < nen_; ++inode)
      for (int idim = 0; idim < nsd_ + 1; ++idim)
        lmdisp[inode * (nsd_ + 1) + idim] = la[ndsdisp].lm_[inode * numdispdofpernode + idim];

    // extract local values of displacement field from global state vector
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd_ + 1, nen_>>(*dispnp, edispnp_, lmdisp);

    // add nodal displacements to point coordinates
    UpdateNodeCoordinates();
  }
  else
    edispnp_.Clear();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvaluateAction(DRT::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    SCATRA::BoundaryAction action, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  std::vector<int>& lm = la[0].lm_;

  switch (action)
  {
    case SCATRA::bd_calc_normal_vectors:
    {
      CalcNormalVectors(params, ele);
      break;
    }

    case SCATRA::bd_integrate_shape_functions:
    {
      // NOTE: add area value only for elements which are NOT ghosted!
      const bool addarea = (ele->Owner() == discretization.Comm().MyPID());
      IntegrateShapeFunctions(ele, params, elevec1_epetra, addarea);

      break;
    }

    case SCATRA::bd_calc_mass_matrix:
    {
      CalcMatMass(ele, elemat1_epetra);

      break;
    }

    case SCATRA::bd_calc_Neumann:
    {
      DRT::Condition* condition = params.get<DRT::Condition*>("condition");
      if (condition == NULL) dserror("Cannot access Neumann boundary condition!");

      EvaluateNeumann(ele, params, discretization, *condition, la, elevec1_epetra, 1.);

      break;
    }

    case SCATRA::bd_calc_Neumann_inflow:
    {
      NeumannInflow(ele, params, discretization, la, elemat1_epetra, elevec1_epetra);

      break;
    }

    case SCATRA::bd_calc_convective_heat_transfer:
    {
      // get the parent element including its material
      DRT::Element* parentele = ele->ParentElement();
      Teuchos::RCP<MAT::Material> mat = parentele->Material();

      // get values of scalar
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");

      // extract local values from global vector
      std::vector<LINALG::Matrix<nen_, 1>> ephinp(numdofpernode_, LINALG::Matrix<nen_, 1>(true));
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*phinp, ephinp, lm);

      // get condition
      Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition>>("condition");
      if (cond == Teuchos::null) dserror("Cannot access condition 'TransportThermoConvections'!");

      // get heat transfer coefficient and surrounding temperature
      const double heatranscoeff = cond->GetDouble("coeff");
      const double surtemp = cond->GetDouble("surtemp");

      ConvectiveHeatTransfer(
          ele, mat, ephinp, elemat1_epetra, elevec1_epetra, heatranscoeff, surtemp);

      break;
    }

    case SCATRA::bd_calc_weak_Dirichlet:
    {
      // get the parent element including its material
      DRT::Element* parentele = ele->ParentElement();
      Teuchos::RCP<MAT::Material> mat = parentele->Material();

      if (numscal_ > 1) dserror("not yet implemented for more than one scalar\n");

      switch (distype)
      {
        // 2D:
        case DRT::Element::line2:
        {
          if (ele->ParentElement()->Shape() == DRT::Element::quad4)
          {
            WeakDirichlet<DRT::Element::line2, DRT::Element::quad4>(
                ele, params, discretization, mat, elemat1_epetra, elevec1_epetra);
          }
          else
          {
            dserror("expected combination quad4/hex8 or line2/quad4 for surface/parent pair");
          }
          break;
        }

        // 3D:
        case DRT::Element::quad4:
        {
          if (ele->ParentElement()->Shape() == DRT::Element::hex8)
            WeakDirichlet<DRT::Element::quad4, DRT::Element::hex8>(
                ele, params, discretization, mat, elemat1_epetra, elevec1_epetra);

          else
            dserror("expected combination quad4/hex8 or line2/quad4 for surface/parent pair");

          break;
        }

        default:
        {
          dserror("not implemented yet\n");
          break;
        }
      }

      break;
    }

    case SCATRA::bd_calc_fs3i_surface_permeability:
    {
      EvaluateSurfacePermeability(ele, params, discretization, la, elemat1_epetra, elevec1_epetra);

      break;
    }

    case SCATRA::bd_calc_fps3i_surface_permeability:
    {
      EvaluateKedemKatchalsky(ele, params, discretization, la, elemat1_epetra, elevec1_epetra);

      break;
    }

    case SCATRA::bd_add_convective_mass_flux:
    {
      // calculate integral of convective mass/heat flux
      // NOTE: since results are added to a global vector via normal assembly
      //       it would be wrong to suppress results for a ghosted boundary!

      // get actual values of transported scalars
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");

      // extract local values from the global vector
      std::vector<LINALG::Matrix<nen_, 1>> ephinp(numdofpernode_, LINALG::Matrix<nen_, 1>(true));
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*phinp, ephinp, lm);

      // get number of dofset associated with velocity related dofs
      const int ndsvel = params.get<int>("ndsvel");

      // get convective (velocity - mesh displacement) velocity at nodes
      Teuchos::RCP<const Epetra_Vector> convel =
          discretization.GetState(ndsvel, "convective velocity field");
      if (convel == Teuchos::null) dserror("Cannot get state vector convective velocity");

      // determine number of velocity related dofs per node
      const int numveldofpernode = la[ndsvel].lm_.size() / nen_;

      // construct location vector for velocity related dofs
      std::vector<int> lmvel((nsd_ + 1) * nen_, -1);
      for (int inode = 0; inode < nen_; ++inode)
        for (int idim = 0; idim < nsd_ + 1; ++idim)
          lmvel[inode * (nsd_ + 1) + idim] = la[ndsvel].lm_[inode * numveldofpernode + idim];

      // we deal with a (nsd_+1)-dimensional flow field
      LINALG::Matrix<nsd_ + 1, nen_> econvel(true);

      // extract local values of convective velocity field from global state vector
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd_ + 1, nen_>>(*convel, econvel, lmvel);

      // rotate the vector field in the case of rotationally symmetric boundary conditions
      rotsymmpbc_->RotateMyValuesIfNecessary(econvel);

      // for the moment we ignore the return values of this method
      CalcConvectiveFlux(ele, ephinp, econvel, elevec1_epetra);
      // vector<double> locfluxintegral = CalcConvectiveFlux(ele,ephinp,evel,elevec1_epetra);
      // std::cout<<"locfluxintegral[0] = "<<locfluxintegral[0]<<std::endl;

      break;
    }

    case SCATRA::bd_calc_s2icoupling:
    {
      EvaluateS2ICoupling(
          ele, params, discretization, la, elemat1_epetra, elemat2_epetra, elevec1_epetra);

      break;
    }

    case SCATRA::bd_calc_s2icoupling_od:
    {
      EvaluateS2ICouplingOD(ele, params, discretization, la, elemat1_epetra);
      break;
    }

    case SCATRA::bd_calc_boundary_integral:
    {
      CalcBoundaryIntegral(ele, elevec1_epetra);
      break;
    }
    case SCATRA::bd_integrate_weighted_scalar:
    {
      IntegrateWeightedScalar(params, ele, elevec1_epetra);
      break;
    }
    case SCATRA::bd_calc_Robin:
    {
      CalcRobinBoundary(ele, params, discretization, la, elemat1_epetra, elevec1_epetra, 1.);
      break;
    }
    case SCATRA::bd_calc_mechanotransduction:
    {
      CalcMechanotransduction(ele, params, discretization, lm, elevec1_epetra);
      break;
    }


    default:
    {
      dserror("Not acting on this boundary action. Forgot implementation?");
      break;
    }
  }

  return 0;
}


/*----------------------------------------------------------------------*
 | evaluate Neumann boundary condition                        gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvaluateNeumann(DRT::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, DRT::Condition& condition,
    DRT::Element::LocationArray& la, Epetra_SerialDenseVector& elevec1, const double scalar)
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // find out whether we will use a time curve
  const double time = scatraparamstimint_->Time();

  // get values, switches and spatial functions from the condition
  // (assumed to be constant on element boundary)
  const int numdof = condition.GetInt("numdof");
  const std::vector<int>* onoff = condition.Get<std::vector<int>>("onoff");
  const std::vector<double>* val = condition.Get<std::vector<double>>("val");
  const std::vector<int>* func = condition.Get<std::vector<int>>("funct");

  if (numdofpernode_ != numdof)
    dserror(
        "The NUMDOF you have entered in your TRANSPORT NEUMANN CONDITION does not equal the number "
        "of scalars.");

  // integration loop
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    double fac = EvalShapeFuncAndIntFac(intpoints, iquad);

    // factor given by spatial function
    double functfac = 1.0;

    // determine global coordinates of current Gauss point
    double coordgp[3];  // we always need three coordinates for function evaluation!
    for (int i = 0; i < 3; ++i) coordgp[i] = 0.;
    for (int i = 0; i < nsd_; ++i)
    {
      coordgp[i] = 0.;
      for (int j = 0; j < nen_; ++j) coordgp[i] += xyze_(i, j) * funct_(j);
    }

    int functnum = -1;
    const double* coordgpref = &coordgp[0];  // needed for function evaluation

    for (int dof = 0; dof < numdofpernode_; ++dof)
    {
      if ((*onoff)[dof])  // is this dof activated?
      {
        // factor given by spatial function
        if (func) functnum = (*func)[dof];

        if (functnum > 0)
        {
          // evaluate function at current Gauss point (provide always 3D coordinates!)
          functfac = DRT::Problem::Instance()->Funct(functnum - 1).Evaluate(dof, coordgpref, time);
        }
        else
          functfac = 1.;

        const double val_fac_funct_fac = (*val)[dof] * fac * functfac;

        for (int node = 0; node < nen_; ++node)
          // TODO: with or without eps_
          elevec1[node * numdofpernode_ + dof] += scalar * funct_(node) * val_fac_funct_fac;
      }  // if ((*onoff)[dof])
    }    // loop over dofs
  }      // loop over integration points

  return 0;
}


/*----------------------------------------------------------------------*
 | calculate normals vectors                                   vg 03/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::CalcNormalVectors(
    Teuchos::ParameterList& params, DRT::FaceElement* ele)
{
  // access the global vector
  const Teuchos::RCP<Epetra_MultiVector> normals =
      params.get<Teuchos::RCP<Epetra_MultiVector>>("normal vectors", Teuchos::null);
  if (normals == Teuchos::null) dserror("Could not access vector 'normal vectors'");

  // determine constant outer normal to this element
  GetConstNormal(normal_, xyze_);

  // loop over the element nodes
  for (int j = 0; j < nen_; j++)
  {
    const int nodegid = (ele->Nodes()[j])->Id();
    if (normals->Map().MyGID(nodegid))
    {  // OK, the node belongs to this processor

      // scaling to a unit vector is performed on the global level after
      // assembly of nodal contributions since we have no reliable information
      // about the number of boundary elements adjacent to a node
      for (int dim = 0; dim < (nsd_ + 1); dim++)
      {
        normals->SumIntoGlobalValue(nodegid, dim, normal_(dim));
      }
    }
    // else: the node belongs to another processor; the ghosted
    //      element will contribute the right value on that proc
  }
}


/*----------------------------------------------------------------------*
 | calculate Neumann inflow boundary conditions                vg 03/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::NeumannInflow(const DRT::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& emat, Epetra_SerialDenseVector& erhs)
{
  // get location vector associated with primary dofset
  std::vector<int>& lm = la[0].lm_;

  // get parent element
  DRT::Element* parentele = ele->ParentElement();

  // get material of parent element
  Teuchos::RCP<MAT::Material> material = parentele->Material();

  // we don't know the parent element's lm vector; so we have to build it here
  const int nenparent = parentele->NumNode();
  std::vector<int> lmparent(nenparent);
  std::vector<int> lmparentowner;
  std::vector<int> lmparentstride;
  parentele->LocationVector(discretization, lmparent, lmparentowner, lmparentstride);

  // get values of scalar
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");

  // extract local values from global vector
  std::vector<LINALG::Matrix<nen_, 1>> ephinp(numdofpernode_, LINALG::Matrix<nen_, 1>(true));
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*phinp, ephinp, lm);

  // get number of dofset associated with velocity related dofs
  const int ndsvel = params.get<int>("ndsvel");

  // get convective (velocity - mesh displacement) velocity at nodes
  Teuchos::RCP<const Epetra_Vector> convel =
      discretization.GetState(ndsvel, "convective velocity field");
  if (convel == Teuchos::null) dserror("Cannot get state vector convective velocity");

  // determine number of velocity related dofs per node
  const int numveldofpernode = la[ndsvel].lm_.size() / nen_;

  // construct location vector for velocity related dofs
  std::vector<int> lmvel((nsd_ + 1) * nen_, -1);
  for (int inode = 0; inode < nen_; ++inode)
    for (int idim = 0; idim < nsd_ + 1; ++idim)
      lmvel[inode * (nsd_ + 1) + idim] = la[ndsvel].lm_[inode * numveldofpernode + idim];

  // we deal with a (nsd_+1)-dimensional flow field
  LINALG::Matrix<nsd_ + 1, nen_> econvel(true);

  // extract local values of convective velocity field from global state vector
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd_ + 1, nen_>>(*convel, econvel, lmvel);

  // rotate the vector field in the case of rotationally symmetric boundary conditions
  rotsymmpbc_->RotateMyValuesIfNecessary(econvel);

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over all scalars
  for (int k = 0; k < numdofpernode_; ++k)
  {
    // loop over all integration points
    for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
    {
      const double fac = EvalShapeFuncAndIntFac(intpoints, iquad, &normal_);

      // get velocity at integration point
      velint_.Multiply(econvel, funct_);

      // normal velocity
      const double normvel = velint_.Dot(normal_);

      if (normvel < -0.0001)
      {
        // set density to 1.0
        double dens = GetDensity(material, ephinp, k);

        // integration factor for left-hand side
        const double lhsfac = dens * normvel * scatraparamstimint_->TimeFac() * fac;

        // integration factor for right-hand side
        double rhsfac = 0.0;
        if (scatraparamstimint_->IsIncremental() and scatraparamstimint_->IsGenAlpha())
          rhsfac = lhsfac / scatraparamstimint_->AlphaF();
        else if (not scatraparamstimint_->IsIncremental() and scatraparamstimint_->IsGenAlpha())
          rhsfac = lhsfac * (1.0 - scatraparamstimint_->AlphaF()) / scatraparamstimint_->AlphaF();
        else if (scatraparamstimint_->IsIncremental() and not scatraparamstimint_->IsGenAlpha())
          rhsfac = lhsfac;

        // matrix
        for (int vi = 0; vi < nen_; ++vi)
        {
          const double vlhs = lhsfac * funct_(vi);

          const int fvi = vi * numdofpernode_ + k;

          for (int ui = 0; ui < nen_; ++ui)
          {
            const int fui = ui * numdofpernode_ + k;

            emat(fvi, fui) -= vlhs * funct_(ui);
          }
        }

        // scalar at integration point
        const double phi = funct_.Dot(ephinp[k]);

        // rhs
        const double vrhs = rhsfac * phi;
        for (int vi = 0; vi < nen_; ++vi)
        {
          const int fvi = vi * numdofpernode_ + k;

          erhs[fvi] += vrhs * funct_(vi);
        }
      }
    }
  }

  return;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::NeumannInflow


/*----------------------------------------------------------------------*
 | get density at integration point                          fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::GetDensity(
    Teuchos::RCP<const MAT::Material> material, const std::vector<LINALG::Matrix<nen_, 1>>& ephinp,
    const int k)
{
  // initialization
  double density(0.);

  // get density depending on material
  switch (material->MaterialType())
  {
    case INPAR::MAT::m_matlist:
    {
      const MAT::MatList* actmat = static_cast<const MAT::MatList*>(material.get());

      const int matid = actmat->MatID(0);

      if (actmat->MaterialById(matid)->MaterialType() == INPAR::MAT::m_scatra)
        // set density to unity
        density = 1.;
      else
        dserror("type of material found in material list is not supported");

      break;
    }

    case INPAR::MAT::m_matlist_reactions:
    case INPAR::MAT::m_scatra:
    {
      // set density to unity
      density = 1.;

      break;
    }

    default:
    {
      dserror("Invalid material type!");
      break;
    }
  }

  return density;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::GetDensity


/*----------------------------------------------------------------------*
 | calculate integral of convective flux across boundary      gjb 11/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
std::vector<double> DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::CalcConvectiveFlux(
    const DRT::FaceElement* ele, const std::vector<LINALG::Matrix<nen_, 1>>& ephinp,
    const LINALG::Matrix<nsd_ + 1, nen_>& evelnp, Epetra_SerialDenseVector& erhs)
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  std::vector<double> integralflux(numscal_);

  // loop over all scalars
  for (int k = 0; k < numscal_; ++k)
  {
    integralflux[k] = 0.0;

    // loop over all integration points
    for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
    {
      const double fac = EvalShapeFuncAndIntFac(intpoints, iquad, &normal_);

      // get velocity at integration point
      velint_.Multiply(evelnp, funct_);

      // normal velocity (note: normal_ is already a unit(!) normal)
      const double normvel = velint_.Dot(normal_);

      // scalar at integration point
      const double phi = funct_.Dot(ephinp[k]);

      const double val = phi * normvel * fac;
      integralflux[k] += val;
      // add contribution to provided vector (distribute over nodes using shape fct.)
      for (int vi = 0; vi < nen_; ++vi)
      {
        const int fvi = vi * numdofpernode_ + k;
        erhs[fvi] += val * funct_(vi);
      }
    }
  }

  return integralflux;

}  // ScaTraEleBoundaryCalc<distype>::ConvectiveFlux

/*----------------------------------------------------------------------*
 | calculate boundary cond. due to convective heat transfer    vg 10/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::ConvectiveHeatTransfer(
    const DRT::FaceElement* ele, Teuchos::RCP<const MAT::Material> material,
    const std::vector<LINALG::Matrix<nen_, 1>>& ephinp, Epetra_SerialDenseMatrix& emat,
    Epetra_SerialDenseVector& erhs, const double heatranscoeff, const double surtemp)
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over all scalars
  for (int k = 0; k < numdofpernode_; ++k)
  {
    // loop over all integration points
    for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
    {
      const double fac = EvalShapeFuncAndIntFac(intpoints, iquad, &normal_);

      // get specific heat capacity at constant volume
      double shc = 0.0;
      if (material->MaterialType() == INPAR::MAT::m_th_fourier_iso)
      {
        const MAT::FourierIso* actmat = static_cast<const MAT::FourierIso*>(material.get());

        shc = actmat->Capacity();
      }
      else if (material->MaterialType() == INPAR::MAT::m_thermostvenant)
      {
        const MAT::ThermoStVenantKirchhoff* actmat =
            static_cast<const MAT::ThermoStVenantKirchhoff*>(material.get());

        shc = actmat->Capacity();
      }
      else
        dserror("Material type is not supported for convective heat transfer!");

      // integration factor for left-hand side
      const double lhsfac = heatranscoeff * scatraparamstimint_->TimeFac() * fac / shc;

      // integration factor for right-hand side
      double rhsfac = 0.0;
      if (scatraparamstimint_->IsIncremental() and scatraparamstimint_->IsGenAlpha())
        rhsfac = lhsfac / scatraparamstimint_->AlphaF();
      else if (not scatraparamstimint_->IsIncremental() and scatraparamstimint_->IsGenAlpha())
        rhsfac = lhsfac * (1.0 - scatraparamstimint_->AlphaF()) / scatraparamstimint_->AlphaF();
      else if (scatraparamstimint_->IsIncremental() and not scatraparamstimint_->IsGenAlpha())
        rhsfac = lhsfac;

      // matrix
      for (int vi = 0; vi < nen_; ++vi)
      {
        const double vlhs = lhsfac * funct_(vi);

        const int fvi = vi * numdofpernode_ + k;

        for (int ui = 0; ui < nen_; ++ui)
        {
          const int fui = ui * numdofpernode_ + k;

          emat(fvi, fui) -= vlhs * funct_(ui);
        }
      }

      // scalar at integration point
      const double phi = funct_.Dot(ephinp[k]);

      // rhs
      const double vrhs = rhsfac * (phi - surtemp);
      for (int vi = 0; vi < nen_; ++vi)
      {
        const int fvi = vi * numdofpernode_ + k;

        erhs[fvi] += vrhs * funct_(vi);
      }
    }
  }

  return;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::ConvectiveHeatTransfer


/*-------------------------------------------------------------------------------------------------------------------------------------*
 | compute shape derivatives, i.e., derivatives of square root of determinant of metric tensor
 w.r.t. spatial coordinates   fang 11/17 |
 *-------------------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvalShapeDerivatives(
    LINALG::Matrix<nsd_ + 1, nen_>& shapederivatives  //!< shape derivatives to be computed
)
{
  // safety check
  if (nsd_ != 2) dserror("Computation of shape derivatives only implemented for 2D interfaces!");

  // compute derivatives of spatial coordinates w.r.t. reference coordinates
  static LINALG::Matrix<nsd_, nsd_ + 1> dxyz_drs;
  dxyz_drs.MultiplyNT(deriv_, xyze_);

  // compute basic components of shape derivatives
  const double xr(dxyz_drs(0, 0)), xs(dxyz_drs(1, 0)), yr(dxyz_drs(0, 1)), ys(dxyz_drs(1, 1)),
      zr(dxyz_drs(0, 2)), zs(dxyz_drs(1, 2));
  const double denominator_inv =
      1. / sqrt(xr * xr * ys * ys + xr * xr * zs * zs - 2 * xr * xs * yr * ys -
                2 * xr * xs * zr * zs + xs * xs * yr * yr + xs * xs * zr * zr + yr * yr * zs * zs -
                2 * yr * ys * zr * zs + ys * ys * zr * zr);
  const double numerator_xr = xr * ys * ys + xr * zs * zs - xs * yr * ys - xs * zr * zs;
  const double numerator_xs = -(xr * yr * ys + xr * zr * zs - xs * yr * yr - xs * zr * zr);
  const double numerator_yr = -(xr * xs * ys - xs * xs * yr - yr * zs * zs + ys * zr * zs);
  const double numerator_ys = xr * xr * ys - xr * xs * yr - yr * zr * zs + ys * zr * zr;
  const double numerator_zr = -(xr * xs * zs - xs * xs * zr + yr * ys * zs - ys * ys * zr);
  const double numerator_zs = xr * xr * zs - xr * xs * zr + yr * yr * zs - yr * ys * zr;

  // compute shape derivatives
  for (int ui = 0; ui < nen_; ++ui)
  {
    shapederivatives(0, ui) =
        denominator_inv * (numerator_xr * deriv_(0, ui) + numerator_xs * deriv_(1, ui));
    shapederivatives(1, ui) =
        denominator_inv * (numerator_yr * deriv_(0, ui) + numerator_ys * deriv_(1, ui));
    shapederivatives(2, ui) =
        denominator_inv * (numerator_zr * deriv_(0, ui) + numerator_zs * deriv_(1, ui));
  }

  return;
}


/*----------------------------------------------------------------------*
 | evaluate shape functions and int. factor at int. point     gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvalShapeFuncAndIntFac(
    const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
    const int iquad,                                         ///< id of current Gauss point
    LINALG::Matrix<1 + nsd_, 1>* normalvec  ///< normal vector at Gauss point(optional)
)
{
  // coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (int idim = 0; idim < nsd_; idim++)
  {
    xsi_(idim) = gpcoord[idim];
  }

  if (not DRT::NURBS::IsNurbs(distype))
  {
    // shape functions and their first derivatives
    DRT::UTILS::shape_function<distype>(xsi_, funct_);
    DRT::UTILS::shape_function_deriv1<distype>(xsi_, deriv_);
  }
  else  // nurbs elements are always somewhat special...
  {
    DRT::NURBS::UTILS::nurbs_get_funct_deriv(funct_, deriv_, xsi_, myknots_, weights_, distype);
  }

  // the metric tensor and the area of an infinitesimal surface/line element
  // optional: get normal at integration point as well
  double drs(0.0);
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>(
      xyze_, deriv_, metrictensor_, drs, normalvec);

  // for nurbs elements the normal vector must be scaled with a special orientation factor!!
  if (DRT::NURBS::IsNurbs(distype))
  {
    if (normalvec != NULL) normal_.Scale(normalfac_);
  }

  // return the integration factor
  return intpoints.IP().qwgt[iquad] * drs;
}


/*----------------------------------------------------------------------*
 | get constant normal                                        gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::GetConstNormal(
    LINALG::Matrix<nsd_ + 1, 1>& normal, const LINALG::Matrix<nsd_ + 1, nen_>& xyze)
{
  // determine normal to this element
  if (not DRT::NURBS::IsNurbs(distype))
  {
    switch (nsd_)
    {
      case 2:
      {
        LINALG::Matrix<3, 1> dist1(true), dist2(true);
        for (int i = 0; i < 3; i++)
        {
          dist1(i) = xyze(i, 1) - xyze(i, 0);
          dist2(i) = xyze(i, 2) - xyze(i, 0);
        }

        normal(0) = dist1(1) * dist2(2) - dist1(2) * dist2(1);
        normal(1) = dist1(2) * dist2(0) - dist1(0) * dist2(2);
        normal(2) = dist1(0) * dist2(1) - dist1(1) * dist2(0);
      }
      break;
      case 1:
      {
        normal(0) = xyze(1, 1) - xyze(1, 0);
        normal(1) = (-1.0) * (xyze(0, 1) - xyze(0, 0));
      }
      break;
      default:
        dserror("Illegal number of space dimensions: %d", nsd_);
        break;
    }  // switch(nsd)
  }
  else  // NURBS case
  {
    // ToDo: this is only a temporary solution in order to have something here.
    // Current handling of node-based normal vectors not applicable in NURBS case
#if 0
    // use one integration point at element center
    const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToStabGaussRule<distype>::rule);
    // hack: ele-id = -1
    // for nurbs elements the normal vector must be scaled with a special orientation factor!!
    // this is already part of this function call
    EvalShapeFuncAndIntFac(intpoints,0,&normal);
#endif
    normal(0) = 1.0;
  }

  // length of normal to this element
  const double length = normal.Norm2();
  // outward-pointing normal of length 1.0
  if (length > EPS10)
    normal.Scale(1 / length);
  else
    dserror("Zero length for element normal");

  return;
}  // ScaTraEleBoundaryCalc<distype>::GetConstNormal


/*----------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling condition       fang 10/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvaluateS2ICoupling(
    const DRT::FaceElement* ele,              ///< current boundary element
    Teuchos::ParameterList& params,           ///< parameter list
    DRT::Discretization& discretization,      ///< discretization
    DRT::Element::LocationArray& la,          ///< location array
    Epetra_SerialDenseMatrix& eslavematrix,   ///< element matrix for slave side
    Epetra_SerialDenseMatrix& emastermatrix,  ///< element matrix for master side
    Epetra_SerialDenseVector& eslaveresidual  ///< element residual for slave side
)
{
  // extract local nodal values on present and opposite sides of scatra-scatra interface
  ExtractNodeValues(discretization, la);
  std::vector<LINALG::Matrix<nen_, 1>> emasterphinp(numscal_);
  ExtractNodeValues(emasterphinp, discretization, la, "imasterphinp");

  // get current scatra-scatra interface coupling condition
  Teuchos::RCP<DRT::Condition> s2icondition = params.get<Teuchos::RCP<DRT::Condition>>("condition");
  if (s2icondition == Teuchos::null)
    dserror("Cannot access scatra-scatra interface coupling condition!");

  // dummy element matrix and vector
  Epetra_SerialDenseMatrix dummymatrix;
  Epetra_SerialDenseVector dummyvector;

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac =
        DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvalShapeFuncAndIntFac(intpoints, gpid);

    // evaluate overall integration factors
    const double timefacfac = scatraparamstimint_->TimeFac() * fac;
    const double timefacrhsfac = scatraparamstimint_->TimeFacRhs() * fac;
    if (timefacfac < 0. or timefacrhsfac < 0.) dserror("Integration factor is negative!");

    EvaluateS2ICouplingAtIntegrationPoint<distype>(*s2icondition, ephinp_, emasterphinp, funct_,
        funct_, funct_, funct_, numscal_, timefacfac, timefacrhsfac, eslavematrix, emastermatrix,
        dummymatrix, dummymatrix, eslaveresidual, dummyvector);
  }  // end of loop over integration points

  return;
}


/*---------------------------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling condition at integration point   fang 05/16 |
 *---------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType distype_master>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvaluateS2ICouplingAtIntegrationPoint(
    DRT::Condition& s2icondition,  //!< scatra-scatra interface coupling condition
    const std::vector<LINALG::Matrix<nen_, 1>>&
        eslavephinp,  //!< state variables at slave-side nodes
    const std::vector<
        LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>>&
        emasterphinp,                            //!< state variables at master-side nodes
    const LINALG::Matrix<nen_, 1>& funct_slave,  //!< slave-side shape function values
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        funct_master,                           //!< master-side shape function values
    const LINALG::Matrix<nen_, 1>& test_slave,  //!< slave-side test function values
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
        test_master,             //!< master-side test function values
    const int numscal,           //!< number of transported scalars
    const double timefacfac,     //!< time-integration factor times domain-integration factor
    const double timefacrhsfac,  //!< time-integration factor for right-hand side times
                                 //!< domain-integration factor
    Epetra_SerialDenseMatrix&
        k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseMatrix&
        k_sm,  //!< linearizations of slave-side residuals w.r.t. master-side dofs
    Epetra_SerialDenseMatrix&
        k_ms,  //!< linearizations of master-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseMatrix&
        k_mm,  //!< linearizations of master-side residuals w.r.t. master-side dofs
    Epetra_SerialDenseVector& r_s,  //!< slave-side residual vector
    Epetra_SerialDenseVector& r_m   //!< master-side residual vector
)
{
  // number of nodes of master-side mortar element
  const int nen_master = DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement;

  // loop over scalars
  for (int k = 0; k < numscal; ++k)
  {
    // evaluate dof values at current integration point on slave and master sides of scatra-scatra
    // interface
    const double slavephiint = funct_slave.Dot(eslavephinp[k]);
    const double masterphiint = funct_master.Dot(emasterphinp[k]);

    // compute matrix and vector contributions according to kinetic model for current scatra-scatra
    // interface coupling condition
    switch (s2icondition.GetInt("kinetic model"))
    {
      // constant permeability model
      case INPAR::S2I::kinetics_constperm:
      {
        // access real vector of constant permeabilities associated with current condition
        const std::vector<double>* permeabilities =
            s2icondition.GetMutable<std::vector<double>>("permeabilities");
        if (permeabilities == NULL)
          dserror("Cannot access vector of permeabilities for scatra-scatra interface coupling!");
        if (permeabilities->size() != (unsigned)numscal)
          dserror("Number of permeabilities does not match number of scalars!");

        // core residual
        const double N = timefacrhsfac * (*permeabilities)[k] * (slavephiint - masterphiint);

        // core linearizations
        const double dN_dc_slave = timefacfac * (*permeabilities)[k];
        const double dN_dc_master = -dN_dc_slave;

        if (k_ss.M() and k_sm.M() and r_s.Length())
        {
          for (int vi = 0; vi < nen_; ++vi)
          {
            const int fvi = vi * numscal + k;

            for (int ui = 0; ui < nen_; ++ui)
              k_ss(fvi, ui * numscal + k) += test_slave(vi) * dN_dc_slave * funct_slave(ui);

            for (int ui = 0; ui < nen_master; ++ui)
              k_sm(fvi, ui * numscal + k) += test_slave(vi) * dN_dc_master * funct_master(ui);

            r_s[fvi] -= test_slave(vi) * N;
          }
        }
        else if (k_ss.M() or k_sm.M() or r_s.Length())
          dserror("Must provide both slave-side matrices and slave-side vector or none of them!");

        if (k_ms.M() and k_mm.M() and r_m.Length())
        {
          for (int vi = 0; vi < nen_master; ++vi)
          {
            const int fvi = vi * numscal + k;

            for (int ui = 0; ui < nen_; ++ui)
              k_ms(fvi, ui * numscal + k) -= test_master(vi) * dN_dc_slave * funct_slave(ui);

            for (int ui = 0; ui < nen_master; ++ui)
              k_mm(fvi, ui * numscal + k) -= test_master(vi) * dN_dc_master * funct_master(ui);

            r_m[fvi] += test_master(vi) * N;
          }
        }
        else if (k_ms.M() or k_mm.M() or r_m.Length())
          dserror("Must provide both master-side matrices and master-side vector or none of them!");

        break;
      }

      default:
      {
        dserror("Kinetic model for scatra-scatra interface coupling not yet implemented!");
        break;
      }
    }
  }  // end of loop over scalars

  return;
}


/*---------------------------------------------------------------------------------------------------------------------------*
 | evaluate off-diagonal system matrix contributions associated with scatra-scatra interface
 coupling condition   fang 11/17 |
 *---------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvaluateS2ICouplingOD(
    const DRT::FaceElement* ele,            ///< current boundary element
    Teuchos::ParameterList& params,         ///< parameter list
    DRT::Discretization& discretization,    ///< discretization
    DRT::Element::LocationArray& la,        ///< location array
    Epetra_SerialDenseMatrix& eslavematrix  ///< element matrix for slave side
)
{
  // extract local nodal values on present and opposite side of scatra-scatra interface
  ExtractNodeValues(discretization, la);
  std::vector<LINALG::Matrix<nen_, 1>> emasterphinp(numscal_, LINALG::Matrix<nen_, 1>(true));
  ExtractNodeValues(emasterphinp, discretization, la, "imasterphinp");

  // get current scatra-scatra interface coupling condition
  Teuchos::RCP<DRT::Condition> s2icondition = params.get<Teuchos::RCP<DRT::Condition>>("condition");
  if (s2icondition == Teuchos::null)
    dserror("Cannot access scatra-scatra interface coupling condition!");

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions at current integration point
    EvalShapeFuncAndIntFac(intpoints, gpid);

    // evaluate shape derivatives
    static LINALG::Matrix<nsd_ + 1, nen_> shapederivatives;
    EvalShapeDerivatives(shapederivatives);

    // evaluate overall integration factor
    const double timefacwgt = scatraparamstimint_->TimeFac() * intpoints.IP().qwgt[gpid];
    if (timefacwgt < 0.) dserror("Integration factor is negative!");

    // loop over scalars
    for (int k = 0; k < numscal_; ++k)
    {
      // evaluate dof values at current integration point on slave and master sides of scatra-scatra
      // interface
      const double slavephiint = funct_.Dot(ephinp_[k]);
      const double masterphiint = funct_.Dot(emasterphinp[k]);

      // compute matrix contributions according to kinetic model for current scatra-scatra interface
      // coupling condition
      switch (s2icondition->GetInt("kinetic model"))
      {
        // constant permeability model
        case INPAR::S2I::kinetics_constperm:
        {
          // access real vector of constant permeabilities associated with current condition
          const std::vector<double>* permeabilities =
              s2icondition->GetMutable<std::vector<double>>("permeabilities");
          if (permeabilities == NULL)
            dserror("Cannot access vector of permeabilities for scatra-scatra interface coupling!");
          if (permeabilities->size() != (unsigned)numscal_)
            dserror("Number of permeabilities does not match number of scalars!");

          // core linearization
          const double dN_dd_slave =
              timefacwgt * (*permeabilities)[k] * (slavephiint - masterphiint);

          // loop over matrix columns
          for (int ui = 0; ui < nen_; ++ui)
          {
            const int fui = ui * 3;

            // loop over matrix rows
            for (int vi = 0; vi < nen_; ++vi)
            {
              const int fvi = vi * numscal_ + k;
              const double vi_dN_dd_slave = funct_(vi) * dN_dd_slave;

              // loop over spatial dimensions
              for (unsigned dim = 0; dim < 3; ++dim)
                // compute linearizations w.r.t. slave-side structural displacements
                eslavematrix(fvi, fui + dim) += vi_dN_dd_slave * shapederivatives(dim, ui);
            }
          }

          break;
        }

        default:
        {
          dserror("Kinetic model for scatra-scatra interface coupling not yet implemented!");
          break;
        }
      }  // selection of kinetic model
    }    // loop over scalars
  }      // loop over integration points

  return;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvaluateS2ICouplingOD


/*-----------------------------------------------------------------------------*
 | extract nodal state variables associated with boundary element   fang 01/17 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::ExtractNodeValues(
    const DRT::Discretization& discretization,  //!< discretization
    DRT::Element::LocationArray& la             //!< location array
)
{
  // extract nodal state variables associated with time t_{n+1} or t_{n+alpha_f}
  ExtractNodeValues(ephinp_, discretization, la);

  return;
}


/*-----------------------------------------------------------------------------*
 | extract nodal state variables associated with boundary element   fang 01/17 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::ExtractNodeValues(
    LINALG::Matrix<nen_, 1>& estate,            //!< nodal state variables
    const DRT::Discretization& discretization,  //!< discretization
    DRT::Element::LocationArray& la,            //!< location array
    const std::string& statename,               //!< name of relevant state
    const int& nds                              //!< number of relevant dofset
    ) const
{
  // initialize matrix vector
  std::vector<LINALG::Matrix<nen_, 1>> estate_temp(1, LINALG::Matrix<nen_, 1>(true));

  // call more general routine
  ExtractNodeValues(estate_temp, discretization, la, statename, nds);

  // copy extracted state variables
  estate = estate_temp[0];

  return;
}


/*-----------------------------------------------------------------------------*
 | extract nodal state variables associated with boundary element   fang 01/17 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::ExtractNodeValues(
    std::vector<LINALG::Matrix<nen_, 1>>& estate,  //!< nodal state variables
    const DRT::Discretization& discretization,     //!< discretization
    DRT::Element::LocationArray& la,               //!< location array
    const std::string& statename,                  //!< name of relevant state
    const int& nds                                 //!< number of relevant dofset
    ) const
{
  // extract global state vector from discretization
  const Teuchos::RCP<const Epetra_Vector> state = discretization.GetState(nds, statename);
  if (state == Teuchos::null)
    dserror("Cannot extract state vector \"" + statename + "\" from discretization!");

  // extract nodal state variables associated with boundary element
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*state, estate, la[nds].lm_);

  return;
}


/*----------------------------------------------------------------------------------*
 | calculate boundary integral, i.e., surface area of boundary element   fang 07/15 |
 *----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::CalcBoundaryIntegral(
    const DRT::FaceElement* ele,      //!< the element we are dealing with
    Epetra_SerialDenseVector& scalar  //!< result vector for scalar integral to be computed
)
{
  // initialize variable for boundary integral
  double boundaryintegral(0.);

  // get integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate values of shape functions and boundary integration factor at current integration
    // point
    const double fac =
        DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvalShapeFuncAndIntFac(intpoints, iquad);

    // add contribution from current integration point to boundary integral
    boundaryintegral += fac;
  }  // loop over integration points

  // write result into result vector
  scalar(0) = boundaryintegral;

  return;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::CalcBoundaryIntegral


/*----------------------------------------------------------------------------------*
 | calculate boundary mass matrix                                        fang 07/16 |
 *----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::CalcMatMass(
    const DRT::FaceElement* const element,  //!< boundary element
    Epetra_SerialDenseMatrix& massmatrix    //!< element mass matrix
)
{
  // get integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate values of shape functions and boundary integration factor at current integration
    // point
    const double fac =
        DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvalShapeFuncAndIntFac(intpoints, iquad);

    // add contribution from current integration point to element mass matrix
    for (int k = 0; k < numdofpernode_; ++k)
    {
      for (int vi = 0; vi < nen_; ++vi)
      {
        const int fvi = vi * numdofpernode_ + k;

        for (int ui = 0; ui < nen_; ++ui)
          massmatrix(fvi, ui * numdofpernode_ + k) += funct_(vi) * funct_(ui) * fac;
      }
    }
  }  // loop over integration points


  return;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::CalcBoundaryIntegral


/*----------------------------------------------------------------------------------*
 | calculate boundary integral of scalars                               rauch 08/16 |
 *----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::IntegrateWeightedScalar(
    Teuchos::ParameterList& params,   //!< parameter list
    const DRT::FaceElement* ele,      //!< the element we are dealing with
    Epetra_SerialDenseVector& result  //!< result vector for scalar integral to be computed
)
{
  // extract values from parameter list.
  // these values are set in scatra_timint_ost_endoexocytosis.
  const int scalarid =
      params.get<int>("ScalarID");  //!< dof id of integrated scalar (position in result vector)
  const double scalar = params.get<double>("scalar");  //!< scalar value to be integrated
  const double prefac =
      params.get<double>("user defined prefac");  //!< user defined factor to integral

  // get integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  //////////////////////////////////////
  // loop over integration points
  //////////////////////////////////////
  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac =
        DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvalShapeFuncAndIntFac(intpoints, gpid);

    // evaluate element right-hand side vector
    for (int vi = 0; vi < nen_; ++vi)
    {
      const int fvi = vi * numscal_ + scalarid;
      result[fvi] -= fac * prefac * funct_(vi) * scalar;
    }  // loop over nodes

  }  // loop over integration points

  return;

}  // DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::IntegrateWeightedScalar


/*----------------------------------------------------------------------*
 | evaluate Robin boundary condition                    schoeder 03/15  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::CalcRobinBoundary(DRT::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la,  ///< location array
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseVector& elevec1_epetra,
    const double scalar)
{
  //////////////////////////////////////////////////////////////////////
  //              get current condition and parameters
  //////////////////////////////////////////////////////////////////////

  // get current condition
  Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition>>("condition");
  if (cond == Teuchos::null) dserror("Cannot access condition 'TransportRobin'");

  // get on/off flags
  const std::vector<int>* onoff = cond->Get<std::vector<int>>("onoff");

  // safety check
  if ((int)(onoff->size()) != numscal_)
    dserror(
        "Mismatch in size for Robin boundary conditions, onoff has length %i, but you have %i "
        "scalars",
        onoff->size(), numscal_);

  // extract prefactor and reference value from condition
  const double prefac = cond->GetDouble("prefactor");
  const double refval = cond->GetDouble("refvalue");

  //////////////////////////////////////////////////////////////////////
  //                  read nodal values
  //////////////////////////////////////////////////////////////////////

  std::vector<int>& lm = la[0].lm_;

  // ------------get values of scalar transport------------------
  // extract global state vector from discretization
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (phinp == Teuchos::null) dserror("Cannot read state vector \"phinp\" from discretization!");

  // extract local nodal state variables from global state vector
  std::vector<LINALG::Matrix<nen_, 1>> ephinp(numdofpernode_, LINALG::Matrix<nen_, 1>(true));
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*phinp, ephinp, lm);

  //////////////////////////////////////////////////////////////////////
  //                  build RHS and StiffMat
  //////////////////////////////////////////////////////////////////////

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over all scalars
  for (int k = 0; k < numscal_; ++k)
  {
    // flag for dofs to be considered by robin conditions
    if ((*onoff)[k] == 1)
    {
      for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
      {
        // evaluate values of shape functions and domain integration factor at current integration
        // point
        const double intfac =
            DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvalShapeFuncAndIntFac(intpoints, gpid);

        // evaluate reference concentration factor
        const double refconcfac = FacForRefConc(gpid, ele, params, discretization);

        // evaluate overall integration factors
        const double fac_3 = prefac * intfac * refconcfac;

        // evaluate current scalar at current integration point
        const double phinp = funct_.Dot(ephinp[k]);

        // build RHS and matrix
        {
          //////////////////////////////////////////////////////////////////////
          //                  rhs
          //////////////////////////////////////////////////////////////////////
          const double vrhs = scatraparamstimint_->TimeFacRhs() * (phinp - refval) * fac_3;

          for (int vi = 0; vi < nen_; ++vi)
          {
            const int fvi = vi * numscal_ + k;

            elevec1_epetra[fvi] += vrhs * funct_(vi);
          }

          //////////////////////////////////////////////////////////////////////
          //                  matrix
          //////////////////////////////////////////////////////////////////////
          for (int vi = 0; vi < nen_; ++vi)
          {
            const double vlhs = scatraparamstimint_->TimeFac() * fac_3 * funct_(vi);
            const int fvi = vi * numscal_ + k;

            for (int ui = 0; ui < nen_; ++ui)
            {
              const int fui = ui * numdofpernode_ + k;

              elemat1_epetra(fvi, fui) -= vlhs * funct_(ui);
            }
          }
        }
      }  // loop over integration points
    }    // if((*onoff)[k]==1)
    // else //in the case of "OFF", a no flux condition is automatically applied

  }  // loop over scalars

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate surface/interface permeability                  Thon 11/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvaluateSurfacePermeability(
    const DRT::FaceElement* ele,          ///< current boundary element
    Teuchos::ParameterList& params,       ///< parameter list
    DRT::Discretization& discretization,  ///< discretization
    DRT::Element::LocationArray& la,      ///< location array
    Epetra_SerialDenseMatrix& elemat1,    ///< element matrix for slave side
    Epetra_SerialDenseVector& elevec1     ///< element residual for slave side
)
{
  //////////////////////////////////////////////////////////////////////
  //                  read nodal values
  //////////////////////////////////////////////////////////////////////

  if (scatraparamstimint_->IsGenAlpha() or not scatraparamstimint_->IsIncremental())
    dserror("Not a valid time integration scheme!");

  std::vector<int>& lm = la[0].lm_;

  // ------------get values of scalar transport------------------
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");
  // extract local values from global vector
  std::vector<LINALG::Matrix<nen_, 1>> ephinp(numdofpernode_, LINALG::Matrix<nen_, 1>(true));
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*phinp, ephinp, lm);

  //------------get membrane concentration at the interface (i.e. within the
  // membrane)------------------
  Teuchos::RCP<const Epetra_Vector> phibar = discretization.GetState("MembraneConcentration");
  if (phibar == Teuchos::null) dserror("Cannot get state vector 'MembraneConcentration'");
  // extract local values from global vector
  std::vector<LINALG::Matrix<nen_, 1>> ephibar(numdofpernode_, LINALG::Matrix<nen_, 1>(true));
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*phibar, ephibar, lm);

  // ------------get values of wall shear stress-----------------------
  // get number of dofset associated with pressure related dofs
  const int ndswss = params.get<int>("ndswss", -1);
  if (ndswss == -1) dserror("Cannot get number of dofset of wss vector");
  Teuchos::RCP<const Epetra_Vector> wss = discretization.GetState(ndswss, "WallShearStress");
  if (wss == Teuchos::null) dserror("Cannot get state vector 'WallShearStress'");

  // determine number of velocity (and pressure) related dofs per node
  const int numwssdofpernode = la[ndswss].lm_.size() / nen_;
  // construct location vector for wss related dofs
  std::vector<int> lmwss((nsd_ + 1) * nen_, -1);
  for (int inode = 0; inode < nen_; ++inode)
    for (int idim = 0; idim < nsd_ + 1; ++idim)
      lmwss[inode * (nsd_ + 1) + idim] = la[ndswss].lm_[inode * numwssdofpernode + idim];

  LINALG::Matrix<nsd_ + 1, nen_> ewss(true);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd_ + 1, nen_>>(*wss, ewss, lmwss);

  // rotate the vector field in the case of rotationally symmetric boundary conditions
  // rotsymmpbc_->RotateMyValuesIfNecessary(ewss);

  //////////////////////////////////////////////////////////////////////
  //                  get current condition
  //////////////////////////////////////////////////////////////////////

  Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition>>("condition");
  if (cond == Teuchos::null) dserror("Cannot access condition 'SurfacePermeability'");

  const std::vector<int>* onoff = cond->Get<std::vector<int>>("onoff");

  const double perm = cond->GetDouble("permeability coefficient");

  // get flag if concentration flux across membrane is affected by local wall shear stresses: 0->no
  // 1->yes
  const bool wss_onoff = (bool)cond->GetInt("wss onoff");

  const std::vector<double>* coeffs = cond->Get<std::vector<double>>("wss coeffs");

  //////////////////////////////////////////////////////////////////////
  //                  build RHS and StiffMat
  //////////////////////////////////////////////////////////////////////
  {
    // integration points and weights
    const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(
        SCATRA::DisTypeToOptGaussRule<distype>::rule);

    // define vector for wss concentration values at nodes
    //  LINALG::Matrix<nen_,1> fwssnod(true);

    // loop over all scalars
    for (int k = 0; k < numdofpernode_; ++k)
    {
      // flag for dofs to be considered by membrane equations of Kedem and Katchalsky
      if ((*onoff)[k] == 1)
      {
        // loop over all integration points
        for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
        {
          const double fac = EvalShapeFuncAndIntFac(intpoints, iquad, &normal_);
          const double refconcfac = FacForRefConc(iquad, ele, params, discretization);
          // integration factor for right-hand side
          double facfac = 0.0;
          if (scatraparamstimint_->IsIncremental() and not scatraparamstimint_->IsGenAlpha())
            facfac = scatraparamstimint_->TimeFac() * fac * refconcfac;
          else
            dserror("EvaluateSurfacePermeability: Requested scheme not yet implemented");

          // scalar at integration point
          const double phi = funct_.Dot(ephinp[k]);

          // permeabilty scaling factor (depending on the norm of the wss) at integration point
          const double facWSS = WSSinfluence(ewss, wss_onoff, coeffs);

          // matrix
          for (int vi = 0; vi < nen_; ++vi)
          {
            const double vlhs = facfac * facWSS * perm * funct_(vi);
            const int fvi = vi * numdofpernode_ + k;

            for (int ui = 0; ui < nen_; ++ui)
            {
              const int fui = ui * numdofpernode_ + k;

              elemat1(fvi, fui) += vlhs * funct_(ui);
            }
          }

          // rhs
          const double vrhs = facfac * facWSS * perm * phi;

          for (int vi = 0; vi < nen_; ++vi)
          {
            const int fvi = vi * numdofpernode_ + k;

            elevec1[fvi] -= vrhs * funct_(vi);
          }
        }
      }  // if((*onoff)[k]==1)
      // else //in the case of "OFF", a no flux condition is automatically applied
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate Kedem-Katchalsky interface                   Thon 11/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvaluateKedemKatchalsky(
    const DRT::FaceElement* ele,          ///< current boundary element
    Teuchos::ParameterList& params,       ///< parameter list
    DRT::Discretization& discretization,  ///< discretization
    DRT::Element::LocationArray& la,      ///< location array
    Epetra_SerialDenseMatrix& elemat1,    ///< element matrix for slave side
    Epetra_SerialDenseVector& elevec1     ///< element residual for slave side
)
{
  // safety checks
  if (scatraparamstimint_->IsGenAlpha() or not scatraparamstimint_->IsIncremental())
    dserror("Not a valid time integration scheme!");

  std::vector<int>& lm = la[0].lm_;

  // ------------get values of scalar transport------------------
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");
  // extract local values from global vector
  std::vector<LINALG::Matrix<nen_, 1>> ephinp(numdofpernode_, LINALG::Matrix<nen_, 1>(true));
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*phinp, ephinp, lm);


  //------------get membrane concentration at the interface (i.e. within the
  // membrane)------------------
  Teuchos::RCP<const Epetra_Vector> phibar = discretization.GetState("MembraneConcentration");
  if (phibar == Teuchos::null) dserror("Cannot get state vector 'MembraneConcentration'");
  // extract local values from global vector
  std::vector<LINALG::Matrix<nen_, 1>> ephibar(numdofpernode_, LINALG::Matrix<nen_, 1>(true));
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*phibar, ephibar, lm);


  //--------get values of pressure at the interface ----------------------
  // get number of dofset associated with pressure related dofs
  const int ndspres = params.get<int>("ndspres", -1);
  if (ndspres == -1) dserror("Cannot get number of dofset of pressure vector");
  Teuchos::RCP<const Epetra_Vector> pressure = discretization.GetState(ndspres, "Pressure");
  if (pressure == Teuchos::null) dserror("Cannot get state vector 'Pressure'");

  // determine number of velocity (and pressure) related dofs per node
  const int numveldofpernode = la[ndspres].lm_.size() / nen_;
  // construct location vector for pressure related dofs
  std::vector<int> lmpres(nen_, -1);
  for (int inode = 0; inode < nen_; ++inode)
    lmpres[inode] = la[ndspres].lm_[inode * numveldofpernode + nsd_ + 1];  // only pressure dofs

  LINALG::Matrix<nen_, 1> epressure(true);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*pressure, epressure, lmpres);

  // rotate the vector field in the case of rotationally symmetric boundary conditions
  // rotsymmpbc_->RotateMyValuesIfNecessary(epressure);


  // ------------get values of wall shear stress-----------------------
  // get number of dofset associated with pressure related dofs
  const int ndswss = params.get<int>("ndswss", -1);
  if (ndswss == -1) dserror("Cannot get number of dofset of wss vector");
  Teuchos::RCP<const Epetra_Vector> wss = discretization.GetState(ndswss, "WallShearStress");
  if (wss == Teuchos::null) dserror("Cannot get state vector 'WallShearStress'");

  // determine number of velocity (and pressure) related dofs per node
  const int numwssdofpernode = la[ndswss].lm_.size() / nen_;
  // construct location vector for wss related dofs
  std::vector<int> lmwss((nsd_ + 1) * nen_, -1);
  for (int inode = 0; inode < nen_; ++inode)
    for (int idim = 0; idim < nsd_ + 1; ++idim)
      lmwss[inode * (nsd_ + 1) + idim] = la[ndswss].lm_[inode * numwssdofpernode + idim];

  LINALG::Matrix<nsd_ + 1, nen_> ewss(true);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd_ + 1, nen_>>(*wss, ewss, lmwss);

  // rotate the vector field in the case of rotationally symmetric boundary conditions
  // rotsymmpbc_->RotateMyValuesIfNecessary(ewss);


  // ------------get current condition----------------------------------
  Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition>>("condition");
  if (cond == Teuchos::null)
    dserror("Cannot access condition 'DESIGN SCATRA COUPLING SURF CONDITIONS'");

  const std::vector<int>* onoff = cond->Get<std::vector<int>>("onoff");

  // get the standard permeability of the interface
  const double perm = cond->GetDouble("permeability coefficient");

  // get flag if concentration flux across membrane is affected by local wall shear stresses: 0->no
  // 1->yes
  const bool wss_onoff = (bool)cond->GetInt("wss onoff");
  const std::vector<double>* coeffs = cond->Get<std::vector<double>>("wss coeffs");

  // hydraulic conductivity at interface
  const double conductivity = cond->GetDouble("hydraulic conductivity");

  // Staverman filtration coefficient at interface
  const double sigma = cond->GetDouble("filtration coefficient");

  ///////////////////////////////////////////////////////////////////////////
  // ------------do the actual calculations----------------------------------
  ///////////////////////////////////////////////////////////////////////////

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over all scalars
  for (int k = 0; k < numdofpernode_; ++k)  // numdofpernode_//1
  {
    // flag for dofs to be considered by membrane equations of Kedem and Katchalsky
    if ((*onoff)[k] == 1)
    {
      // loop over all integration points
      for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
      {
        const double fac = EvalShapeFuncAndIntFac(intpoints, iquad, &normal_);

        // integration factor
        double facfac = 0.0;
        if (scatraparamstimint_->IsIncremental() and not scatraparamstimint_->IsGenAlpha())
          facfac = scatraparamstimint_->TimeFac() * fac;
        else
          dserror("Kedem-Katchalsky: Requested time integration scheme not yet implemented");

        // scalar at integration point
        const double phi = funct_.Dot(ephinp[k]);

        // pressure at integration point
        const double p = funct_.Dot(epressure);

        // mean concentration at integration point
        const double phibar = funct_.Dot(ephibar[k]);

        // mean concentration at integration point
        const double facWSS = WSSinfluence(ewss, wss_onoff, coeffs);


        // matrix
        for (int vi = 0; vi < nen_; ++vi)
        {
          const double vlhs = facfac * facWSS * perm * funct_(vi);

          const int fvi = vi * numdofpernode_ + k;

          for (int ui = 0; ui < nen_; ++ui)
          {
            const int fui = ui * numdofpernode_ + k;

            elemat1(fvi, fui) += vlhs * funct_(ui);
          }
        }

        // rhs
        const double vrhs =
            facfac * facWSS * (perm * phi + (1 - sigma) * phibar * conductivity * p);
        // J_s =f_WSS*[perm*(phi1-phi2)+(1-sigma)*phibar*conductivity*(p1-p2)]
        // J_s := solute flux through scalar scalar interface
        // perm:=membrane permeability
        // sigma:=Staverman filtration coefficient
        // phibar:= mean concentration within the membrane (for now: simply linear interpolated, but
        // other interpolations also possible) conductivity:=local hydraulic conductivity of
        // membrane

        for (int vi = 0; vi < nen_; ++vi)
        {
          const int fvi = vi * numdofpernode_ + k;

          elevec1[fvi] -= vrhs * funct_(vi);
        }
      }
    }  // if((*onoff)[k]==1)
    // else //in the case of "OFF", a no flux condition is automatically applied
  }
  return;
}


/*----------------------------------------------------------------------*
 |  Factor of WSS dependent interface transport           hemmler 07/14 |
 |  Calculated as in Calvez, V., "Mathematical and numerical modeling of early atherosclerotic
 lesions",ESAIM: Proceedings. Vol. 30. EDP Sciences, 2010
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::WSSinfluence(
    const LINALG::Matrix<nsd_ + 1, nen_>& ewss, const bool wss_onoff,
    const std::vector<double>* coeffs)
{
  // permeabilty scaling factor at integration point
  double facWSS = 1.0;

  if (wss_onoff)
  {
    LINALG::Matrix<nsd_ + 1, 1> wss(true);
    for (int ii = 0; ii < nsd_ + 1; ii++)
      for (int jj = 0; jj < nen_; jj++) wss(ii) += ewss(ii, jj) * funct_(jj);

    // euklidian norm of act node wss
    const double wss_norm = sqrt(wss(0) * wss(0) + wss(1) * wss(1) + wss(2) * wss(2));
    facWSS = log10(1 + coeffs->at(0) / (wss_norm + coeffs->at(1))) /
             log10(2);  // empirical function (log law) to account for influence of WSS;
  }
  // else //no WSS influence

  return facWSS;
}


/*----------------------------------------------------------------------*
 | computes mechanoresponsive scalar transport              rauch 01/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::CalcMechanotransduction(DRT::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1_epetra)
{
  const int RhoGEF_dof = DRT::Problem::Instance()
                             ->CellMigrationParams()
                             .sublist("SCALAR TRANSPORT DOF IDS")
                             .get<int>("RhoGEF");

  if (RhoGEF_dof > -1)
  {
    const double time = scatraparamstimint_->Time();
    params.set<double>("total time", time);
    params.set<int>("FromSactraBoundary", 1);
    const double dt = scatraparamstimint_->Dt();
    params.set<double>("delta time", dt);

    // cast FaceElement to TransportBoundary
    DRT::ELEMENTS::TransportBoundary* transp_brdyele =
        dynamic_cast<DRT::ELEMENTS::TransportBoundary*>(ele);
    if (transp_brdyele == NULL)
      dserror("cast from FaceElement to TransportBoundary element failed");

    // get number of dof-set associated with displacement related dofs
    const int ndsdisp = params.get<int>("ndsdisp");

    // get displacement state from scatra discretization
    Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState(ndsdisp, "dispnp");
    if (dispnp == Teuchos::null) dserror("Cannot get state vector 'dispnp'");

    // get parent element
    DRT::Element* parentele = ele->ParentElement();
    static const int nenparent = parentele->NumNode();
    if (nenparent != 8) dserror("only implemented for hex8 elements");

    DRT::Element::LocationArray parent_la(discretization.NumDofSets());
    parentele->LocationVector(discretization, parent_la, false);

    LINALG::Matrix<nsd_ + 1, 8> myeledispnp;
    // extract local values of displacement field from global state vector
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd_ + 1, 8>>(
        *dispnp, myeledispnp, parent_la[ndsdisp].lm_);

    // integration points and weights for boundary (!) gp
    const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(
        SCATRA::DisTypeToOptGaussRule<distype>::rule);

    // get coordinates of gauss points w.r.t. local parent coordinate system
    LINALG::SerialDenseMatrix pqxg(intpoints.IP().nquad, nsd_ + 1);
    LINALG::Matrix<nsd_ + 1, 1> pxsi(true);
    LINALG::Matrix<nsd_ + 1, nsd_ + 1> derivtrafo(true);
    DRT::UTILS::BoundaryGPToParentGP<nsd_ + 1>(
        pqxg, derivtrafo, intpoints, parentele->Shape(), distype, transp_brdyele->SurfaceNumber());

    // LINALG::Matrix<nsd_+1,nenparent>  xrefe; // material coord. of parent element
    LINALG::Matrix<3, 8> xrefe;
    LINALG::Matrix<3, 8> xcurr;  // current  coord. of parent element
    DRT::Node** nodes = parentele->Nodes();

    // update element geometry of parent element
    for (int i = 0; i < nenparent; ++i)
    {
      const double* x = nodes[i]->X();
      for (int j = 0; j < nsd_ + 1; ++j)
      {
        xrefe(j, i) = x[j];
        xcurr(j, i) = xrefe(j, i) + myeledispnp(j, i);
      }
    }
    Teuchos::RCP<MAT::So3Material> so3mat =
        Teuchos::rcp_dynamic_cast<MAT::So3Material>(parentele->Material(1));
    if (so3mat == Teuchos::null) dserror("cast to MAT::So3Material failed");

    params.set<int>("iostress", 0);  // needed for activefiber material; if output is requested only
                                     // active stresses are written

    double averagestress;
    averagestress = 0.0;

    // loop over all integration points
    for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
    {
      // get integration factor and normal at current boundary gp
      const double fac = EvalShapeFuncAndIntFac(intpoints, iquad, &normal_);

      // (Jo) coordinates of the current integration point in parent coordinate system
      for (int idim = 0; idim < nsd_ + 1; idim++)
      {
        pxsi(idim) = pqxg(iquad, idim);
      }
      // std::cout <<"pxsi: "<< pxsi<< '\n';
      LINALG::Matrix<3, 8> pderiv_loc(
          true);  // derivatives of parent element shape functions in parent element coordinate
                  // system //LINALG::Matrix<dim,numnodes> pderiv_loc(true);

      // evaluate derivatives of parent element shape functions at current integration point in
      // parent coordinate system
      DRT::UTILS::shape_function_deriv1<DRT::Element::hex8>(pxsi, pderiv_loc);

      const int NUMDIM_SOH8 = 3;
      const int NUMNOD_SOH8 = 8;

      /* get the inverse of the Jacobian matrix which looks like:
       **            [ x_,r  y_,r  z_,r ]^-1
       **     J^-1 = [ x_,s  y_,s  z_,s ]
       **            [ x_,t  y_,t  z_,t ]
       */
      //    std::cout <<"pderiv_loc"<< std::setprecision(3) << pderiv_loc << '\n';
      //    std::cout <<"xref: "<< std::setprecision(3) << xrefe << '\n';
      LINALG::Matrix<3, 3> invJ;
      // LINALG::Matrix<3,3>   Jmat;
      invJ.MultiplyNT(pderiv_loc, xrefe);
      invJ.Invert();
      // Jmat.MultiplyNT(pderiv_loc,xrefe);

      // compute derivatives N_XYZ at gp w.r.t. material coordinates
      // by N_XYZ = J^-1 * N_rst
      LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_XYZ;
      N_XYZ.Multiply(invJ, pderiv_loc);

      // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
      // build deformation gradient wrt to material configuration
      LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrd(false);
      defgrd.MultiplyNT(xcurr, N_XYZ);


      // Right Cauchy-Green tensor = F^T * F
      LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> cauchygreen;
      cauchygreen.MultiplyTN(defgrd, defgrd);

      // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
      Epetra_SerialDenseVector glstrain_epetra(MAT::NUM_STRESS_3D);
      LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(glstrain_epetra.A(), true);

      // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
      glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
      glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
      glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
      glstrain(3) = cauchygreen(0, 1);
      glstrain(4) = cauchygreen(1, 2);
      glstrain(5) = cauchygreen(2, 0);


      if (parentele->NumMaterial() < 2) dserror("only one material defined for scatra element");

      // evaluate the material to obtain stress
      LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat(true);
      LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress(true);
      params.set<LINALG::Matrix<nsd_ + 1, 1>>("xsi", pxsi);
      params.set<int>("gp", iquad);
      so3mat->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, parentele->Id());

      // transform PK2 to Cauchy
      LINALG::Matrix<3, 3> cauchystress(true);
      LINALG::Matrix<3, 3> PK2stress;

      // calculate the Jacobi-determinant
      double detF = defgrd.Determinant();
      if (abs(detF) < 1.0e-10)
      {
        std::cout << "Error: detF in PK2toCauchy =  \n" << detF << std::endl;
      }
      // Convert stress like 6x1-Voigt vector to 3x3 matrix
      PK2stress(0, 0) = stress(0);
      PK2stress(0, 1) = stress(3);
      PK2stress(0, 2) = stress(5);
      PK2stress(1, 0) = PK2stress(0, 1);
      PK2stress(1, 1) = stress(1);
      PK2stress(1, 2) = stress(4);
      PK2stress(2, 0) = PK2stress(0, 2);
      PK2stress(2, 1) = PK2stress(1, 2);
      PK2stress(2, 2) = stress(2);

      // sigma = 1/J * F * sigma * F^{T}
      LINALG::Matrix<3, 3> temp(true);
      temp.MultiplyNN(defgrd, PK2stress);
      cauchystress.MultiplyNT(temp, defgrd);
      cauchystress.Scale(1. / detF);

      LINALG::Matrix<3, 1> normalstress(true);
      LINALG::Matrix<3, 1> normal(true);
      normal(0) = normal_(0);
      normal(1) = normal_(1);
      normal(2) = normal_(2);
      normalstress.Multiply(cauchystress, normal);

      //    // DEBUG output
      //    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> stressmatrix;
      //    std::cout <<"stress: "<< std::setprecision(3) << stress << '\n';
      //    stressmatrix(0,0) = stress(0);
      //    stressmatrix(1,1) = stress(1);
      //    stressmatrix(2,2) = stress(2);
      //    stressmatrix(1,0) = stress(3);
      //    stressmatrix(0,1) = stress(3);
      //    stressmatrix(1,2) = stress(4);
      //    stressmatrix(2,1) = stress(4);
      //    stressmatrix(2,0) = stress(5);
      //    stressmatrix(0,2) = stress(5);
      //    std::cout <<"boundary stress: "<< std::setprecision(3) << normalstress << '\n';

      double forcemagnitude =
          1.0 * sqrt(normalstress(0) * normalstress(0) + normalstress(1) * normalstress(1) +
                     normalstress(2) * normalstress(2));

      const double timefac = scatraparamstimint_->TimeFacRhs();

      if (forcemagnitude > 1.0e-15)
      {
        for (int node = 0; node < nen_; ++node)
        {
          averagestress += fac * forcemagnitude * funct_(node) * timefac;
        }
      }

    }  // end gp loop

    double constant = params.get<double>("source const", -1);
    if (time > 0.0)
    {
      if (constant == -1) dserror("source const not available in scatra");
    }

    double source;  // = averagestress/area*constant;
    source = averagestress * constant;

    if (abs(source) > 1.0e-15)
    {
      for (int node = 0; node < nen_; ++node)
      {
        int dof = node * numdofpernode_ + RhoGEF_dof;
        elevec1_epetra[dof] += source;
      }
    }

  }  // if valid dof ids provided
  return;
}  // CalcMechanotransduction


/*----------------------------------------------------------------------*
 |  Integrate shapefunctions over surface (private)           gjb 02/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::IntegrateShapeFunctions(
    const DRT::FaceElement* ele, Teuchos::ParameterList& params, Epetra_SerialDenseVector& elevec1,
    const bool addarea)
{
  // access boundary area variable with its actual value
  double boundaryint = params.get<double>("area");

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    const double fac = EvalShapeFuncAndIntFac(intpoints, gpid);

    // compute integral of shape functions
    for (int node = 0; node < nen_; ++node)
      for (int k = 0; k < numdofpernode_; ++k)
        elevec1[node * numdofpernode_ + k] += funct_(node) * fac;

    // area calculation
    if (addarea) boundaryint += fac;
  }  // loop over integration points

  // add contribution to the global value
  params.set<double>("area", boundaryint);

  return;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::IntegrateShapeFunctions


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType bdistype, DRT::Element::DiscretizationType pdistype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::WeakDirichlet(DRT::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    Teuchos::RCP<const MAT::Material> material, Epetra_SerialDenseMatrix& elemat_epetra,
    Epetra_SerialDenseVector& elevec_epetra)
{
  //------------------------------------------------------------------------
  // Dirichlet boundary condition
  //------------------------------------------------------------------------
  Teuchos::RCP<DRT::Condition> dbc = params.get<Teuchos::RCP<DRT::Condition>>("condition");

  // check of total time
  const double time = scatraparamstimint_->Time();

  // get values and spatial functions from condition
  // (assumed to be constant on element boundary)
  const std::vector<double>* val = (*dbc).Get<std::vector<double>>("val");
  const std::vector<int>* func = (*dbc).Get<std::vector<int>>("funct");

  // assign boundary value multiplied by time-curve factor
  double dirichval = (*val)[0];

  // spatial function number
  const int funcnum = (*func)[0];

  //------------------------------------------------------------------------
  // preliminary definitions for (boundary) and parent element and
  // evaluation of nodal values of velocity and scalar based on parent
  // element nodes
  //------------------------------------------------------------------------
  // get the parent element
  DRT::Element* pele = ele->ParentElement();

  // number of spatial dimensions regarding (boundary) element
  static const int bnsd = DRT::UTILS::DisTypeToDim<bdistype>::dim;

  // number of spatial dimensions regarding parent element
  static const int pnsd = DRT::UTILS::DisTypeToDim<pdistype>::dim;

  // number of (boundary) element nodes
  static const int bnen = DRT::UTILS::DisTypeToNumNodePerEle<bdistype>::numNodePerElement;

  // number of parent element nodes
  static const int pnen = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

  // parent element location array
  DRT::Element::LocationArray pla(discretization.NumDofSets());
  pele->LocationVector(discretization, pla, false);

  // get number of dofset associated with velocity related dofs
  const int ndsvel = params.get<int>("ndsvel");

  // get convective (velocity - mesh displacement) velocity at nodes
  Teuchos::RCP<const Epetra_Vector> convel =
      discretization.GetState(ndsvel, "convective velocity field");
  if (convel == Teuchos::null) dserror("Cannot get state vector convective velocity");

  // determine number of velocity related dofs per node
  const int numveldofpernode = pla[ndsvel].lm_.size() / pnen;

  // construct location vector for velocity related dofs
  std::vector<int> plmvel(pnsd * pnen, -1);
  for (int inode = 0; inode < pnen; ++inode)
    for (int idim = 0; idim < pnsd; ++idim)
      plmvel[inode * pnsd + idim] = pla[ndsvel].lm_[inode * numveldofpernode + idim];

  // we deal with a (nsd_+1)-dimensional flow field
  LINALG::Matrix<pnsd, pnen> econvel(true);

  // extract local values of convective velocity field from global state vector
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<pnsd, pnen>>(*convel, econvel, plmvel);

  // rotate the vector field in the case of rotationally symmetric boundary conditions
  rotsymmpbc_->template RotateMyValuesIfNecessary<pnsd, pnen>(econvel);

  // get scalar values at parent element nodes
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");

  // extract local values from global vectors for parent element
  std::vector<LINALG::Matrix<pnen, 1>> ephinp(numscal_);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<pnen, 1>>(*phinp, ephinp, pla[0].lm_);

  //------------------------------------------------------------------------
  // preliminary definitions for integration loop
  //------------------------------------------------------------------------
  // reshape element matrices and vectors and init to zero, construct views
  elemat_epetra.Shape(pnen, pnen);
  elevec_epetra.Size(pnen);
  LINALG::Matrix<pnen, pnen> emat(elemat_epetra.A(), true);
  LINALG::Matrix<pnen, 1> erhs(elevec_epetra.A(), true);

  // (boundary) element local node coordinates
  LINALG::Matrix<pnsd, bnen> bxyze(true);
  GEO::fillInitialPositionArray<bdistype, pnsd, LINALG::Matrix<pnsd, bnen>>(ele, bxyze);

  // parent element local node coordinates
  LINALG::Matrix<pnsd, pnen> pxyze(true);
  GEO::fillInitialPositionArray<pdistype, pnsd, LINALG::Matrix<pnsd, pnen>>(pele, pxyze);

  // coordinates of integration points for (boundary) and parent element
  LINALG::Matrix<bnsd, 1> bxsi(true);
  LINALG::Matrix<pnsd, 1> pxsi(true);

  // transposed jacobian "dx/ds" and inverse of transposed jacobian "ds/dx"
  // for parent element
  LINALG::Matrix<pnsd, pnsd> pxjm(true);
  LINALG::Matrix<pnsd, pnsd> pxji(true);

  // metric tensor for (boundary) element
  LINALG::Matrix<bnsd, bnsd> bmetrictensor(true);

  // (outward-pointing) unit normal vector to (boundary) element
  LINALG::Matrix<pnsd, 1> bnormal(true);

  // velocity vector at integration point
  LINALG::Matrix<pnsd, 1> velint;

  // gradient of scalar value at integration point
  LINALG::Matrix<pnsd, 1> gradphi;

  // (boundary) element shape functions, local and global derivatives
  LINALG::Matrix<bnen, 1> bfunct(true);
  LINALG::Matrix<bnsd, bnen> bderiv(true);
  LINALG::Matrix<bnsd, bnen> bderxy(true);

  // parent element shape functions, local and global derivatives
  LINALG::Matrix<pnen, 1> pfunct(true);
  LINALG::Matrix<pnsd, pnen> pderiv(true);
  LINALG::Matrix<pnsd, pnen> pderxy(true);

  //------------------------------------------------------------------------
  // additional matrices and vectors for mixed-hybrid formulation
  //------------------------------------------------------------------------
  // for volume integrals
  LINALG::Matrix<pnsd * pnen, pnsd * pnen> mat_s_q(true);
  LINALG::Matrix<pnsd * pnen, pnen> mat_s_gradphi(true);

  LINALG::Matrix<pnsd * pnen, 1> vec_s_gradphi(true);

  // for boundary integrals
  LINALG::Matrix<pnen, pnsd * pnen> mat_w_q_o_n(true);
  LINALG::Matrix<pnsd * pnen, pnen> mat_s_o_n_phi(true);

  LINALG::Matrix<pnsd * pnen, 1> vec_s_o_n_phi_minus_g(true);

  // inverse matrix
  LINALG::Matrix<pnsd * pnen, pnsd * pnen> inv_s_q(true);

  //------------------------------------------------------------------------
  // check whether Nitsche (default) or mixed-hybrid formulation as well as
  // preliminary definitions and computations for Nitsche stabilization term
  //------------------------------------------------------------------------
  // default is Nitsche formulation
  bool mixhyb = false;

  // stabilization parameter for Nitsche term
  const double nitsche_stab_para = (*dbc).GetDouble("TauBscaling");

  // if stabilization parameter negative: mixed-hybrid formulation
  if (nitsche_stab_para < 0.0) mixhyb = true;

  // pre-factor for adjoint-consistency term:
  // either 1.0 (adjoint-consistent, default) or -1.0 (adjoint-inconsistent)
  double gamma = 1.0;
  const std::string* consistency = (*dbc).Get<std::string>("Choice of gamma parameter");
  if (*consistency == "adjoint-consistent")
    gamma = 1.0;
  else if (*consistency == "diffusive-optimal")
    gamma = -1.0;
  else
    dserror("unknown definition for gamma parameter: %s", (*consistency).c_str());

  // use one-point Gauss rule to do calculations at element center
  const DRT::UTILS::IntPointsAndWeights<bnsd> intpoints_tau(
      SCATRA::DisTypeToStabGaussRule<bdistype>::rule);

  // element surface area (1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double* gpcoord = (intpoints_tau.IP().qxg)[0];
  for (int idim = 0; idim < bnsd; idim++)
  {
    bxsi(idim) = gpcoord[idim];
  }
  DRT::UTILS::shape_function_deriv1<bdistype>(bxsi, bderiv);
  double drs = 0.0;
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<bdistype>(
      bxyze, bderiv, bmetrictensor, drs, &bnormal);
  const double area = intpoints_tau.IP().qwgt[0] * drs;

  // get number of dimensions for (boundary) element (convert from int to double)
  const double dim = (double)bnsd;

  // computation of characteristic length of (boundary) element
  // (2D: square root of element area, 1D: element length)
  const double h = std::pow(area, (1.0 / dim));

  //------------------------------------------------------------------------
  // preliminary computations for integration loop
  //------------------------------------------------------------------------
  // integration points and weights for (boundary) element and parent element
  const DRT::UTILS::IntPointsAndWeights<bnsd> bintpoints(
      SCATRA::DisTypeToOptGaussRule<bdistype>::rule);

  const DRT::UTILS::IntPointsAndWeights<pnsd> pintpoints(
      SCATRA::DisTypeToOptGaussRule<pdistype>::rule);

  // transfer integration-point coordinates of (boundary) element to parent element
  Epetra_SerialDenseMatrix pqxg(pintpoints.IP().nquad, pnsd);
  {
    Epetra_SerialDenseMatrix gps(bintpoints.IP().nquad, bnsd);

    for (int iquad = 0; iquad < bintpoints.IP().nquad; ++iquad)
    {
      const double* gpcoord = (bintpoints.IP().qxg)[iquad];

      for (int idim = 0; idim < bnsd; idim++)
      {
        gps(iquad, idim) = gpcoord[idim];
      }
    }
    if (pnsd == 2)
    {
      DRT::UTILS::BoundaryGPToParentGP2(pqxg, gps, pdistype, bdistype, ele->FaceParentNumber());
    }
    else if (pnsd == 3)
    {
      DRT::UTILS::BoundaryGPToParentGP3(pqxg, gps, pdistype, bdistype, ele->FaceParentNumber());
    }
  }

  //------------------------------------------------------------------------
  // integration loop 1: volume integrals (only for mixed-hybrid formulation)
  //------------------------------------------------------------------------
  if (mixhyb)
  {
    for (int iquad = 0; iquad < pintpoints.IP().nquad; ++iquad)
    {
      // reference coordinates of integration point from (boundary) element
      const double* gpcoord = (pintpoints.IP().qxg)[iquad];
      for (int idim = 0; idim < pnsd; idim++)
      {
        pxsi(idim) = gpcoord[idim];
      }

      // parent element shape functions and local derivatives
      DRT::UTILS::shape_function<pdistype>(pxsi, pfunct);
      DRT::UTILS::shape_function_deriv1<pdistype>(pxsi, pderiv);

      // Jacobian matrix and determinant of parent element (including check)
      pxjm.MultiplyNT(pderiv, pxyze);
      const double det = pxji.Invert(pxjm);
      if (det < 1E-16)
        dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", pele->Id(), det);

      // compute integration factor
      const double fac = pintpoints.IP().qwgt[iquad] * det;

      // compute global derivatives
      pderxy.Multiply(pxji, pderiv);

      //--------------------------------------------------------------------
      // loop over scalars (not yet implemented for more than one scalar)
      //--------------------------------------------------------------------
      // for(int k=0;k<numdofpernode_;++k)
      int k = 0;
      {
        // get viscosity
        if (material->MaterialType() == INPAR::MAT::m_scatra)
        {
          const MAT::ScatraMat* actmat = static_cast<const MAT::ScatraMat*>(material.get());

          dsassert(numdofpernode_ == 1, "more than 1 dof per node for SCATRA material");

          // get constant diffusivity
          diffus_[k] = actmat->Diffusivity();
        }
        else
          dserror("Material type is not supported");

        // gradient of current scalar value
        gradphi.Multiply(pderxy, ephinp[k]);

        // integration factor for left-hand side
        const double lhsfac = scatraparamstimint_->TimeFac() * fac;

        // integration factor for right-hand side
        double rhsfac = 0.0;
        if (scatraparamstimint_->IsIncremental() and scatraparamstimint_->IsGenAlpha())
          rhsfac = lhsfac / scatraparamstimint_->AlphaF();
        else if (not scatraparamstimint_->IsIncremental() and scatraparamstimint_->IsGenAlpha())
          rhsfac = lhsfac * (1.0 - scatraparamstimint_->AlphaF()) / scatraparamstimint_->AlphaF();
        else if (scatraparamstimint_->IsIncremental() and not scatraparamstimint_->IsGenAlpha())
          rhsfac = lhsfac;

        //--------------------------------------------------------------------
        //  matrix and vector additions due to mixed-hybrid formulation
        //--------------------------------------------------------------------
        /*
                       /         \
                  1   |   h   h  |
              - ----- |  s , q   |
                kappa |          |
                      \          / Omega
        */
        for (int vi = 0; vi < pnen; ++vi)
        {
          const int fvi = vi * numdofpernode_ + k;

          // const double vlhs = lhsfac*pfunct(vi);
          const double vlhs = lhsfac * (1.0 / diffus_[k]) * pfunct(vi);

          for (int ui = 0; ui < pnen; ++ui)
          {
            const int fui = ui * numdofpernode_ + k;

            for (int i = 0; i < pnsd; ++i)
            {
              mat_s_q(fvi * pnsd + i, fui * pnsd + i) -= vlhs * pfunct(ui);
            }
          }
        }

        /*
                       /                  \
                      |  h         /   h\  |
                    + | s  , grad | phi  | |
                      |            \    /  |
                       \                  / Omega
        */
        for (int vi = 0; vi < pnen; ++vi)
        {
          const int fvi = vi * numdofpernode_ + k;

          // const double vlhs = lhsfac*diffus_[k]*pfunct(vi);
          const double vlhs = lhsfac * pfunct(vi);

          for (int ui = 0; ui < pnen; ++ui)
          {
            const int fui = ui * numdofpernode_ + k;

            for (int i = 0; i < pnsd; ++i)
            {
              mat_s_gradphi(fvi * pnsd + i, fui) += vlhs * pderxy(i, ui);
            }
          }

          // const double vrhs = rhsfac*diffus_[k]*pfunct(vi);
          const double vrhs = rhsfac * pfunct(vi);

          for (int i = 0; i < pnsd; ++i)
          {
            vec_s_gradphi(fvi * pnsd + i) += vrhs * gradphi(i);
          }
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // integration loop 2: boundary integrals
  //------------------------------------------------------------------------
  for (int iquad = 0; iquad < bintpoints.IP().nquad; ++iquad)
  {
    // reference coordinates of integration point from (boundary) element
    const double* gpcoord = (bintpoints.IP().qxg)[iquad];
    for (int idim = 0; idim < bnsd; idim++)
    {
      bxsi(idim) = gpcoord[idim];
    }

    // (boundary) element shape functions
    DRT::UTILS::shape_function<bdistype>(bxsi, bfunct);
    DRT::UTILS::shape_function_deriv1<bdistype>(bxsi, bderiv);

    // global coordinates of current integration point from (boundary) element
    LINALG::Matrix<pnsd, 1> coordgp(true);
    for (int A = 0; A < bnen; ++A)
    {
      for (int j = 0; j < pnsd; ++j)
      {
        coordgp(j) += bxyze(j, A) * bfunct(A);
      }
    }

    // reference coordinates of integration point from parent element
    for (int idim = 0; idim < pnsd; idim++)
    {
      pxsi(idim) = pqxg(iquad, idim);
    }

    // parent element shape functions and local derivatives
    DRT::UTILS::shape_function<pdistype>(pxsi, pfunct);
    DRT::UTILS::shape_function_deriv1<pdistype>(pxsi, pderiv);

    // Jacobian matrix and determinant of parent element (including check)
    pxjm.MultiplyNT(pderiv, pxyze);
    const double det = pxji.Invert(pxjm);
    if (det < 1E-16)
      dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", pele->Id(), det);

    // compute measure tensor for surface element, infinitesimal area element drs
    // and (outward-pointing) unit normal vector
    DRT::UTILS::ComputeMetricTensorForBoundaryEle<bdistype>(
        bxyze, bderiv, bmetrictensor, drs, &bnormal);

    // for nurbs elements the normal vector must be scaled with a special orientation factor!!
    if (DRT::NURBS::IsNurbs(distype)) bnormal.Scale(normalfac_);

    // compute integration factor
    const double fac = bintpoints.IP().qwgt[iquad] * drs;

    // compute global derivatives
    pderxy.Multiply(pxji, pderiv);

#if 1
    //--------------------------------------------------------------------
    // check whether integration-point coordinates evaluated from
    // (boundary) and parent element match
    //--------------------------------------------------------------------
    LINALG::Matrix<pnsd, 1> check(true);
    LINALG::Matrix<pnsd, 1> diff(true);

    for (int A = 0; A < pnen; ++A)
    {
      for (int j = 0; j < pnsd; ++j)
      {
        check(j) += pxyze(j, A) * pfunct(A);
      }
    }

    diff = check;
    diff -= coordgp;

    const double norm = diff.Norm2();

    if (norm > 1e-9)
    {
      for (int j = 0; j < pnsd; ++j)
      {
        printf("%12.5e %12.5e\n", check(j), coordgp(j));
      }
      dserror("Gausspoint matching error %12.5e\n", norm);
    }
#endif

    //--------------------------------------------------------------------
    // factor for Dirichlet boundary condition given by spatial function
    //--------------------------------------------------------------------
    double functfac = 1.0;
    if (funcnum > 0)
    {
      // evaluate function at current integration point (important: a 3D position vector is
      // required)
      double coordgp3D[3];
      coordgp3D[0] = 0.0;
      coordgp3D[1] = 0.0;
      coordgp3D[2] = 0.0;
      for (int i = 0; i < pnsd; i++) coordgp3D[i] = coordgp(i);

      functfac = DRT::Problem::Instance()->Funct(funcnum - 1).Evaluate(0, &(coordgp3D[0]), time);
    }
    else
      functfac = 1.0;
    dirichval *= functfac;

    //--------------------------------------------------------------------
    // loop over scalars (not yet implemented for more than one scalar)
    //--------------------------------------------------------------------
    // for(int k=0;k<numdofpernode_;++k)
    int k = 0;
    {
      // get viscosity
      if (material->MaterialType() == INPAR::MAT::m_scatra)
      {
        const MAT::ScatraMat* actmat = static_cast<const MAT::ScatraMat*>(material.get());

        dsassert(numdofpernode_ == 1, "more than 1 dof per node for SCATRA material");

        // get constant diffusivity
        diffus_[k] = actmat->Diffusivity();
      }
      else
        dserror("Material type is not supported");

      // get scalar value at integration point
      const double phi = pfunct.Dot(ephinp[k]);

      // integration factor for left-hand side
      const double lhsfac = scatraparamstimint_->TimeFac() * fac;

      // integration factor for right-hand side
      double rhsfac = 0.0;
      if (scatraparamstimint_->IsIncremental() and scatraparamstimint_->IsGenAlpha())
        rhsfac = lhsfac / scatraparamstimint_->AlphaF();
      else if (not scatraparamstimint_->IsIncremental() and scatraparamstimint_->IsGenAlpha())
        rhsfac = lhsfac * (1.0 - scatraparamstimint_->AlphaF()) / scatraparamstimint_->AlphaF();
      else if (scatraparamstimint_->IsIncremental() and not scatraparamstimint_->IsGenAlpha())
        rhsfac = lhsfac;

      if (mixhyb)
      {
        //--------------------------------------------------------------------
        //  matrix and vector additions due to mixed-hybrid formulation
        //--------------------------------------------------------------------
        /*  consistency term
                    /           \
                   |  h   h     |
                 - | w , q  o n |
                   |            |
                   \            / Gamma
        */
        for (int vi = 0; vi < pnen; ++vi)
        {
          const int fvi = vi * numdofpernode_ + k;

          const double vlhs = lhsfac * pfunct(vi);

          for (int ui = 0; ui < pnen; ++ui)
          {
            const int fui = ui * numdofpernode_ + k;

            for (int i = 0; i < pnsd; ++i)
            {
              mat_w_q_o_n(fvi, fui * pnsd + i) -= vlhs * pfunct(ui) * bnormal(i);
            }
          }
        }

        /*  adjoint consistency term
                    /                 \
                   |  h          h    |
                 - | s  o n , phi - g |
                   |                  |
                   \                  / Gamma
        */
        for (int vi = 0; vi < pnen; ++vi)
        {
          const int fvi = vi * numdofpernode_ + k;

          const double vlhs = lhsfac * pfunct(vi);

          for (int ui = 0; ui < pnen; ++ui)
          {
            const int fui = ui * numdofpernode_ + k;

            for (int i = 0; i < pnsd; ++i)
            {
              mat_s_o_n_phi(fvi * pnsd + i, fui) -= vlhs * pfunct(ui) * bnormal(i);
            }
          }

          for (int i = 0; i < pnsd; ++i)
          {
            vec_s_o_n_phi_minus_g(fvi * pnsd + i) -=
                pfunct(vi) * bnormal(i) *
                (rhsfac * phi - scatraparamstimint_->TimeFac() * fac * dirichval);
          }
        }
      }
      else
      {
        // parameter alpha for Nitsche stabilization term
        const double alpha = nitsche_stab_para * diffus_[k] / h;

        // get velocity at integration point
        velint.Multiply(econvel, pfunct);

        // normal velocity
        const double normvel = velint.Dot(bnormal);

        // gradient of current scalar value
        gradphi.Multiply(pderxy, ephinp[k]);

        // gradient of current scalar value in normal direction
        const double gradphi_norm = bnormal.Dot(gradphi);

        //--------------------------------------------------------------------
        //  matrix and vector additions due to Nitsche formulation
        //--------------------------------------------------------------------
        /*  consistency term
                    /                           \
                   |  h                  h      |
                 - | w , kappa * grad(phi ) o n |
                   |                            |
                   \                            / Gamma
        */
        for (int vi = 0; vi < pnen; ++vi)
        {
          const int fvi = vi * numdofpernode_ + k;

          const double vlhs = lhsfac * pfunct(vi) * diffus_[k];

          for (int ui = 0; ui < pnen; ++ui)
          {
            const int fui = ui * numdofpernode_ + k;

            for (int i = 0; i < pnsd; ++i)
            {
              emat(fvi, fui) -= vlhs * pderxy(i, ui) * bnormal(i);
            }
          }

          const double vrhs = rhsfac * diffus_[k];

          erhs(fvi) += vrhs * pfunct(vi) * gradphi_norm;
        }

        /*  adjoint consistency term, inflow/outflow part
              / --          --                                        \
             |  |         h  |                      h           h     |
           - |  |(a o n) w  +| gamma * kappa *grad(w ) o n , phi - g  |
             |  |            |                                        |
             \  --          --                                        / Gamma_in/out
        */
        for (int vi = 0; vi < pnen; ++vi)
        {
          const int fvi = vi * numdofpernode_ + k;

          // compute diffusive part
          double prefac = 0.0;
          for (int i = 0; i < pnsd; ++i)
          {
            prefac += gamma * diffus_[k] * pderxy(i, vi) * bnormal(i);
          }

          // add convective part in case of inflow boundary
          if (normvel < -0.0001) prefac += normvel * pfunct(vi);

          const double vlhs = lhsfac * prefac;

          for (int ui = 0; ui < pnen; ++ui)
          {
            const int fui = ui * numdofpernode_ + k;

            emat(fvi, fui) -= vlhs * pfunct(ui);
          }

          erhs(fvi) += prefac * (rhsfac * phi - scatraparamstimint_->TimeFac() * fac * dirichval);
        }

        /*  stabilization term
                            /             \
                           |  h     h     |
                 + alpha * | w , phi - g  |
                           |              |
                           \              / Gamma
        */
        for (int vi = 0; vi < pnen; ++vi)
        {
          const int fvi = vi * numdofpernode_ + k;

          const double prefac = alpha * pfunct(vi);

          for (int ui = 0; ui < pnen; ++ui)
          {
            const int fui = ui * numdofpernode_ + k;

            emat(fvi, fui) += lhsfac * prefac * pfunct(ui);
          }

          erhs(fvi) -= prefac * (rhsfac * phi - scatraparamstimint_->TimeFac() * fac * dirichval);
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // local condensation (only for mixed-hybrid formulation)
  //------------------------------------------------------------------------
  if (mixhyb)
  {
    // matrix inversion of flux-flux block
    inv_s_q = mat_s_q;

    LINALG::FixedSizeSerialDenseSolver<pnsd * pnen, pnsd * pnen> solver;

    solver.SetMatrix(inv_s_q);
    solver.Invert();

    // computation of matrix-matrix and matrix vector products, local assembly
    for (int vi = 0; vi < pnen; ++vi)
    {
      for (int ui = 0; ui < pnen; ++ui)
      {
        for (int rr = 0; rr < pnsd * pnen; ++rr)
        {
          for (int mm = 0; mm < pnsd * pnen; ++mm)
          {
            emat(vi, ui) -= mat_w_q_o_n(vi, rr) * inv_s_q(rr, mm) *
                            (mat_s_gradphi(mm, ui) + mat_s_o_n_phi(mm, ui));
          }
        }
      }
    }

    for (int vi = 0; vi < pnen; ++vi)
    {
      for (int rr = 0; rr < pnsd * pnen; ++rr)
      {
        for (int mm = 0; mm < pnsd * pnen; ++mm)
        {
          erhs(vi) -= mat_w_q_o_n(vi, rr) * inv_s_q(rr, mm) *
                      (-vec_s_o_n_phi_minus_g(mm) - vec_s_gradphi(mm));
        }
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | calculate boundary conditions for                                    |
 | impl. Characteristic Galerkin (2nd order)                            |
 | time integration just for the reinitialization equation schott 04/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType bdistype,
    DRT::Element::DiscretizationType pdistype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::ReinitCharacteristicGalerkinBoundary(
    DRT::FaceElement* ele,                       //!< transport element
    Teuchos::ParameterList& params,              //!< parameter list
    DRT::Discretization& discretization,         //!< discretization
    Teuchos::RCP<const MAT::Material> material,  //!< material
    Epetra_SerialDenseMatrix& elemat_epetra,     //!< ele sysmat
    Epetra_SerialDenseVector& elevec_epetra      //!< ele rhs
)
{
  //------------------------------------------------------------------------
  // preliminary definitions for (boundary) and parent element and
  // evaluation of nodal values of velocity and scalar based on parent
  // element nodes
  //------------------------------------------------------------------------
  // get the parent element
  DRT::Element* pele = ele->ParentElement();

  // number of spatial dimensions regarding (boundary) element
  static const int bnsd = DRT::UTILS::DisTypeToDim<bdistype>::dim;

  // number of spatial dimensions regarding parent element
  static const int pnsd = DRT::UTILS::DisTypeToDim<pdistype>::dim;

  // number of (boundary) element nodes
  static const int bnen = DRT::UTILS::DisTypeToNumNodePerEle<bdistype>::numNodePerElement;

  // number of parent element nodes
  static const int pnen = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

  // parent element lm vector
  std::vector<int> plm;
  std::vector<int> plmowner;
  std::vector<int> plmstride;
  pele->LocationVector(discretization, plm, plmowner, plmstride);

  // get scalar values at parent element nodes
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");
  Teuchos::RCP<const Epetra_Vector> phin = discretization.GetState("phin");
  if (phinp == Teuchos::null) dserror("Cannot get state vector 'phin'");

  // extract local values from global vectors for parent element
  std::vector<LINALG::Matrix<pnen, 1>> ephinp(numscal_);
  std::vector<LINALG::Matrix<pnen, 1>> ephin(numscal_);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<pnen, 1>>(*phinp, ephinp, plm);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<pnen, 1>>(*phin, ephin, plm);

  //------------------------------------------------------------------------
  // preliminary definitions for integration loop
  //------------------------------------------------------------------------
  // reshape element matrices and vectors and init to zero, construct views
  elemat_epetra.Shape(pnen, pnen);
  elevec_epetra.Size(pnen);
  LINALG::Matrix<pnen, pnen> emat(elemat_epetra.A(), true);
  LINALG::Matrix<pnen, 1> erhs(elevec_epetra.A(), true);

  // (boundary) element local node coordinates
  LINALG::Matrix<pnsd, bnen> bxyze(true);
  GEO::fillInitialPositionArray<bdistype, pnsd, LINALG::Matrix<pnsd, bnen>>(ele, bxyze);

  // parent element local node coordinates
  LINALG::Matrix<pnsd, pnen> pxyze(true);
  GEO::fillInitialPositionArray<pdistype, pnsd, LINALG::Matrix<pnsd, pnen>>(pele, pxyze);

  // coordinates of integration points for (boundary) and parent element
  LINALG::Matrix<bnsd, 1> bxsi(true);
  LINALG::Matrix<pnsd, 1> pxsi(true);

  // transposed jacobian "dx/ds" and inverse of transposed jacobian "ds/dx"
  // for parent element
  LINALG::Matrix<pnsd, pnsd> pxjm(true);
  LINALG::Matrix<pnsd, pnsd> pxji(true);

  // metric tensor for (boundary) element
  LINALG::Matrix<bnsd, bnsd> bmetrictensor(true);

  // (outward-pointing) unit normal vector to (boundary) element
  LINALG::Matrix<pnsd, 1> bnormal(true);

  // velocity vector at integration point
  LINALG::Matrix<pnsd, 1> velint;

  // gradient of scalar value at integration point
  LINALG::Matrix<pnsd, 1> gradphi;

  // (boundary) element shape functions, local and global derivatives
  LINALG::Matrix<bnen, 1> bfunct(true);
  LINALG::Matrix<bnsd, bnen> bderiv(true);
  LINALG::Matrix<bnsd, bnen> bderxy(true);

  // parent element shape functions, local and global derivatives
  LINALG::Matrix<pnen, 1> pfunct(true);
  LINALG::Matrix<pnsd, pnen> pderiv(true);
  LINALG::Matrix<pnsd, pnen> pderxy(true);


  // use one-point Gauss rule to do calculations at element center
  const DRT::UTILS::IntPointsAndWeights<bnsd> intpoints_tau(
      SCATRA::DisTypeToStabGaussRule<bdistype>::rule);

  // element surface area (1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double* gpcoord = (intpoints_tau.IP().qxg)[0];
  for (int idim = 0; idim < bnsd; idim++)
  {
    bxsi(idim) = gpcoord[idim];
  }
  DRT::UTILS::shape_function_deriv1<bdistype>(bxsi, bderiv);
  double drs = 0.0;
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<bdistype>(
      bxyze, bderiv, bmetrictensor, drs, &bnormal);

  //------------------------------------------------------------------------
  // preliminary computations for integration loop
  //------------------------------------------------------------------------
  // integration points and weights for (boundary) element and parent element
  const DRT::UTILS::IntPointsAndWeights<bnsd> bintpoints(
      SCATRA::DisTypeToOptGaussRule<bdistype>::rule);

  const DRT::UTILS::IntPointsAndWeights<pnsd> pintpoints(
      SCATRA::DisTypeToOptGaussRule<pdistype>::rule);

  // transfer integration-point coordinates of (boundary) element to parent element
  Epetra_SerialDenseMatrix pqxg(pintpoints.IP().nquad, pnsd);
  {
    Epetra_SerialDenseMatrix gps(bintpoints.IP().nquad, bnsd);

    for (int iquad = 0; iquad < bintpoints.IP().nquad; ++iquad)
    {
      const double* gpcoord = (bintpoints.IP().qxg)[iquad];

      for (int idim = 0; idim < bnsd; idim++)
      {
        gps(iquad, idim) = gpcoord[idim];
      }
    }
    if (pnsd == 2)
    {
      DRT::UTILS::BoundaryGPToParentGP2(pqxg, gps, pdistype, bdistype, ele->FaceParentNumber());
    }
    else if (pnsd == 3)
    {
      DRT::UTILS::BoundaryGPToParentGP3(pqxg, gps, pdistype, bdistype, ele->FaceParentNumber());
    }
  }


  const double reinit_pseudo_timestepsize_factor = params.get<double>("pseudotimestepsize_factor");

  const double meshsize = getEleDiameter<pdistype>(pxyze);

  const double pseudo_timestep_size = meshsize * reinit_pseudo_timestepsize_factor;

  //------------------------------------------------------------------------
  // integration loop: boundary integrals
  //------------------------------------------------------------------------
  for (int iquad = 0; iquad < bintpoints.IP().nquad; ++iquad)
  {
    // reference coordinates of integration point from (boundary) element
    const double* gpcoord = (bintpoints.IP().qxg)[iquad];
    for (int idim = 0; idim < bnsd; idim++)
    {
      bxsi(idim) = gpcoord[idim];
    }

    // (boundary) element shape functions
    DRT::UTILS::shape_function<bdistype>(bxsi, bfunct);
    DRT::UTILS::shape_function_deriv1<bdistype>(bxsi, bderiv);

    // global coordinates of current integration point from (boundary) element
    LINALG::Matrix<pnsd, 1> coordgp(true);
    for (int A = 0; A < bnen; ++A)
    {
      for (int j = 0; j < pnsd; ++j)
      {
        coordgp(j) += bxyze(j, A) * bfunct(A);
      }
    }

    // reference coordinates of integration point from parent element
    for (int idim = 0; idim < pnsd; idim++)
    {
      pxsi(idim) = pqxg(iquad, idim);
    }

    // parent element shape functions and local derivatives
    DRT::UTILS::shape_function<pdistype>(pxsi, pfunct);
    DRT::UTILS::shape_function_deriv1<pdistype>(pxsi, pderiv);

    // Jacobian matrix and determinant of parent element (including check)
    pxjm.MultiplyNT(pderiv, pxyze);
    const double det = pxji.Invert(pxjm);
    if (det < 1E-16)
      dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", pele->Id(), det);

    // compute measure tensor for surface element, infinitesimal area element drs
    // and (outward-pointing) unit normal vector
    DRT::UTILS::ComputeMetricTensorForBoundaryEle<bdistype>(
        bxyze, bderiv, bmetrictensor, drs, &bnormal);

    // for nurbs elements the normal vector must be scaled with a special orientation factor!!
    if (DRT::NURBS::IsNurbs(distype)) bnormal.Scale(normalfac_);

    // compute integration factor
    const double fac_surface = bintpoints.IP().qwgt[iquad] * drs;

    // compute global derivatives
    pderxy.Multiply(pxji, pderiv);

    //--------------------------------------------------------------------
    // loop over scalars (not yet implemented for more than one scalar)
    //--------------------------------------------------------------------
    for (int dofindex = 0; dofindex < numdofpernode_; ++dofindex)
    {
      //----------  --------------      |                    |
      //  mat              -1/4* dtau^2 | w, n*grad(D(psi) ) |
      //--------------------------      |                    |

      LINALG::Matrix<1, pnen> derxy_normal;
      derxy_normal.Clear();
      derxy_normal.MultiplyTN(bnormal, pderxy);

      for (int vi = 0; vi < pnen; ++vi)
      {
        const int fvi = vi * numdofpernode_ + dofindex;

        for (int ui = 0; ui < pnen; ++ui)
        {
          const int fui = ui * numdofpernode_ + dofindex;

          emat(fvi, fui) -= pfunct(vi) *
                            (fac_surface * pseudo_timestep_size * pseudo_timestep_size / 4.0) *
                            derxy_normal(0, ui);
        }
      }

      //----------  --------------      |              m     |
      //  rhs               0.5* dtau^2 | w, n*grad(psi )    |
      //--------------------------      |                    |

      // update grad_dist_n
      LINALG::Matrix<pnsd, 1> grad_dist_n(true);
      grad_dist_n.Multiply(pderxy, ephin[dofindex]);

      LINALG::Matrix<1, 1> grad_dist_n_normal(true);
      grad_dist_n_normal.MultiplyTN(bnormal, grad_dist_n);

      for (int vi = 0; vi < pnen; ++vi)
      {
        const int fvi = vi * numdofpernode_ + dofindex;

        erhs(fvi) += pfunct(vi) * pseudo_timestep_size * pseudo_timestep_size * fac_surface / 2.0 *
                     grad_dist_n_normal(0, 0);
      }


      //                    |              m+1     m  |
      //    1/4*delta_tau^2 | w, n*grad(psi   - psi ) |
      //                    |              i          |
      // update grad_dist_n
      LINALG::Matrix<pnsd, 1> grad_dist_npi(true);
      grad_dist_npi.Multiply(pderxy, ephinp[dofindex]);

      LINALG::Matrix<1, 1> grad_dist_npi_normal;
      grad_dist_npi_normal.Clear();
      grad_dist_npi_normal.MultiplyTN(bnormal, grad_dist_npi);

      double Grad_Dpsi_normal = grad_dist_npi_normal(0, 0) - grad_dist_n_normal(0, 0);


      for (int vi = 0; vi < pnen; ++vi)
      {
        const int fvi = vi * numdofpernode_ + dofindex;

        erhs(fvi) += pfunct(vi) * Grad_Dpsi_normal * fac_surface * pseudo_timestep_size *
                     pseudo_timestep_size / 4.0;
      }

    }  // loop over scalars
  }    // loop over integration points

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::tri3>;
// explicit instantiation of template methods
template void DRT::ELEMENTS::ScaTraEleBoundaryCalc<
    DRT::Element::tri3>::EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::tri3>(DRT::Condition&,
    const std::vector<LINALG::Matrix<nen_, 1>>&,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
    const LINALG::Matrix<nen_, 1>&,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement,
        1>&,
    const LINALG::Matrix<nen_, 1>&,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement,
        1>&,
    const int, const double, const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&,
    Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&, Epetra_SerialDenseVector&,
    Epetra_SerialDenseVector&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalc<
    DRT::Element::tri3>::EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::quad4>(DRT::Condition&,
    const std::vector<LINALG::Matrix<nen_, 1>>&,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
    const LINALG::Matrix<nen_, 1>&,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement,
        1>&,
    const LINALG::Matrix<nen_, 1>&,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement,
        1>&,
    const int, const double, const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&,
    Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&, Epetra_SerialDenseVector&,
    Epetra_SerialDenseVector&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalc<
    DRT::Element::quad4>::EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::tri3>(DRT::Condition&,
    const std::vector<LINALG::Matrix<nen_, 1>>&,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
    const LINALG::Matrix<nen_, 1>&,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement,
        1>&,
    const LINALG::Matrix<nen_, 1>&,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement,
        1>&,
    const int, const double, const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&,
    Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&, Epetra_SerialDenseVector&,
    Epetra_SerialDenseVector&);
template void
DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::quad4>::EvaluateS2ICouplingAtIntegrationPoint<
    DRT::Element::quad4>(DRT::Condition&, const std::vector<LINALG::Matrix<nen_, 1>>&,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
    const LINALG::Matrix<nen_, 1>&,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement,
        1>&,
    const LINALG::Matrix<nen_, 1>&,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement,
        1>&,
    const int, const double, const double, Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&,
    Epetra_SerialDenseMatrix&, Epetra_SerialDenseMatrix&, Epetra_SerialDenseVector&,
    Epetra_SerialDenseVector&);
template class DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::nurbs9>;
