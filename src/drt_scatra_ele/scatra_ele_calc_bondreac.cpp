/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_bondreac.cpp

 \brief main file containing routines for calculation of scatra element with reactive scalars and bond dynamics.

 \level 2

 <pre>
   \maintainer Andreas Rauch
               rauch@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289 - 15240
 </pre>
 *----------------------------------------------------------------------*/


#include "scatra_ele_calc_bondreac.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_condition_utils.H"

#include "../drt_mat/matlist_bondreacs.H"
#include "../drt_mat/matlist_reactions.H"
#include "../drt_mat/structporo.H"
#include "../drt_mat/scatra_mat.H"
#include "../drt_mat/matlist.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"

#include "../drt_immersed_problem/immersed_field_exchange_manager.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcBondReac<distype,probdim>::ScaTraEleCalcBondReac(
    const int numdofpernode,
    const int numscal,
    const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::ScaTraEleCalc(
        numdofpernode,
        numscal,
        disname),
      DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype,probdim>::ScaTraEleCalcAdvReac(
        numdofpernode,
        numscal,
        disname)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype,int probdim>
DRT::ELEMENTS::ScaTraEleCalcBondReac<distype,probdim> *
DRT::ELEMENTS::ScaTraEleCalcBondReac<distype,probdim>::Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    const ScaTraEleCalcBondReac *delete_me )
{
  static std::map<std::pair<std::string,int>,ScaTraEleCalcBondReac<distype,probdim>* > instances;

  std::pair<std::string,int> key(disname,numdofpernode);

  if(delete_me == NULL)
  {
    if(instances.find(key) == instances.end())
      instances[key] = new ScaTraEleCalcBondReac<distype,probdim>(numdofpernode,numscal,disname);
  }

  else
  {
    // since we keep several instances around in the general case, we need to
    // find which of the instances to delete with this call. This is done by
    // letting the object to be deleted hand over the 'this' pointer, which is
    // located in the map and deleted
    for( typename std::map<std::pair<std::string,int>,ScaTraEleCalcBondReac<distype,probdim>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
      if ( i->second == delete_me )
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  return instances[key];
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcBondReac<distype,probdim>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, "", this );
}


/*----------------------------------------------------------------------*
 |  get the material constants  (private)                   rauch 12/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcBondReac<distype,probdim>::GetMaterialParams(
    const DRT::Element*  ele,       //!< the element we are dealing with
    std::vector<double>& densn,     //!< density at t_(n)
    std::vector<double>& densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
    std::vector<double>& densam,    //!< density at t_(n+alpha_M)
    double&              visc,      //!< fluid viscosity
    const int            iquad      //!< id of current gauss point
)
{
  // get the material
  Teuchos::RCP<MAT::Material> material = ele->Material();

  // We may have some reactive and some non-reactive elements in one discretization.
  // But since the calculation classes are singleton, we have to reset all reactive stuff in case
  // of non-reactive elements:
  advreac::ReaManager()->Clear(my::numscal_);

  if (material->MaterialType() == INPAR::MAT::m_matlist)
  {
    const Teuchos::RCP<const MAT::MatList> actmat = Teuchos::rcp_dynamic_cast<const MAT::MatList>(material);
    if (actmat->NumMat() != my::numscal_) dserror("Not enough materials in MatList.");

    for (int k = 0;k<my::numscal_;++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP< MAT::Material> singlemat = actmat->MaterialById(matid);

      advreac::Materials(singlemat,k,densn[k],densnp[k],densam[k],visc,iquad);
    }
  }

  else if (material->MaterialType() == INPAR::MAT::m_matlist_reactions)
  {
    const Teuchos::RCP<MAT::MatListReactions> actmat = Teuchos::rcp_dynamic_cast<MAT::MatListReactions>(material);
    if (actmat->NumMat() != my::numscal_) dserror("Not enough materials in MatList.");

    for (int k = 0;k<my::numscal_;++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP< MAT::Material> singlemat = actmat->MaterialById(matid);

      //Note: order is important here!!
      advreac::Materials(singlemat,k,densn[k],densnp[k],densam[k],visc,iquad);

      advreac::SetAdvancedReactionTerms(k,actmat,advreac::GetGpCoord()); //every reaction calculation stuff happens in here!!
    }
  }

  else if (material->MaterialType() == INPAR::MAT::m_matlist_bondreacs)
  {
    // get surface traction  and porosity at gauss point
    const double porosity  = GetPorosityFromBackgroundEle(ele,iquad);
    const double violation = GetAdhesionViolation(ele,iquad);

    const Teuchos::RCP<MAT::MatListBondReacs> actmat = Teuchos::rcp_dynamic_cast<MAT::MatListBondReacs>(material);
    if (actmat->NumMat() != my::numscal_) dserror("Not enough materials in MatList.");

    for (int k = 0;k<my::numscal_;++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP< MAT::Material> singlemat = actmat->MaterialById(matid);

      //Note: order is important here!!
      advreac::Materials(singlemat,k,densn[k],densnp[k],densam[k],visc,iquad);

      SetBondReactionTerms(k,actmat,violation,porosity,advreac::GetGpCoord()); //every reaction calculation stuff happens in here!!
    }
  }

  else
  {
    advreac::Materials(material,0,densn[0],densnp[0],densam[0],visc,iquad);
  }

  return;
} //ScaTraEleCalc::GetMaterialParams


/*-------------------------------------------------------------------------------*
 |  set reac. body force, reaction coefficient and derivatives       rauch 12/16 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcBondReac<distype,probdim>::SetBondReactionTerms(
    const int                                 k,            //!< index of current scalar
    const Teuchos::RCP<MAT::MatListBondReacs> matreaclist,  //!< index of current scalar
    const double                              violation,    //!< penalty violation at current gp
    const double                              porosity,     //!< average receptor-ligand distance
    const double* gpcoord
)
{
  const Teuchos::RCP<ScaTraEleReaManagerAdvReac> remanager = advreac::ReaManager();

  //! scalar values at t_(n+1) or t_(n+alpha_F)
  const std::vector<double>& phinp = my::scatravarmanager_->Phinp();

  //! scalar values at t_(n)
  const std::vector<double>& phin = my::scatravarmanager_->Phin();

  remanager->AddToReaBodyForce( matreaclist->CalcReaBodyForceTerm(k,phinp,phin,violation,porosity,gpcoord) ,k );

  matreaclist->CalcReaBodyForceDerivMatrix(k,remanager->GetReaBodyForceDerivVector(k),phinp,phin,violation,porosity,gpcoord);
}


/*-------------------------------------------------------------------------------*
 |  evaluate single bond traction at gauss point                     rauch 12/16 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
double DRT::ELEMENTS::ScaTraEleCalcBondReac<distype,probdim>::GetAdhesionViolation(
    const DRT::Element*  ele,       //!< the element we are dealing with
    const int            iquad      //!< id of current gauss point
) const
{
  // adhesion violation (return variable)
  double violation=-1234.0;

  // initialize result variable
  double traction=0.0;

  // get global problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // get RCP to auxdis which is evaluated in heterogeneous reaction strategy
  const Teuchos::RCP<const DRT::Discretization> auxdis =
      DRT::ImmersedFieldExchangeManager::Instance()->GetPointerToAuxDis();

  // get RCP to surface traction vector
  Teuchos::RCP<Epetra_Vector> surface_traction =
      DRT::ImmersedFieldExchangeManager::Instance()->GetPointerSurfaceTraction();

  // get discretization
  Teuchos::RCP<DRT::Discretization> dis = problem->GetDis("cellscatra");

  // get adhesion fixpoints
  Teuchos::RCP<const Epetra_Vector> cell_adhesion_fixpoints = dis->GetState(1,"AdhesionFixpoints");
  if(cell_adhesion_fixpoints == Teuchos::null)
    dserror("Cannot get AdhesionFixpoints vector from scatra discretization");

  // get displacement state from structure discretization
  Teuchos::RCP<const Epetra_Vector> dispnp = dis->GetState(1,"dispnp");
  if (dispnp == Teuchos::null)
    dserror("Cannot get displacement vector from scatra discretization");

  // is ele conditioned with CellFocalAdhesion condition?
  bool is_adhesion_element=true;

  // get element location vector
  DRT::Element::LocationArray la(auxdis->NumDofSets());
  ele->LocationVector(*auxdis,la,false);

  // get structure_lm from second dofset
  // the first dofset is the scatra surface and the second dofset the structure
  const std::vector<int>& cell_lm = la[1].lm_;

  // evaluate traction only for elements which wear a CellFocalAdhesion condition.
  // this condition is defined on the cell discretization.
  // the map of vector surface_traction_ contains only conditioned cell dofs.
  const size_t ldim = cell_lm.size();

  for (size_t i=0; i<ldim; ++i)
  {
    const int lid = (*surface_traction).Map().LID(cell_lm[i]);
    if (lid<0)
    {
      is_adhesion_element=false;
      break;
    }
  }

  // extract values if element is adhesion surface element
  if (is_adhesion_element)
  {
    // number of nodes and numdofs per node
    const int numNode = ele->NumNode();

    if(numNode!=4)
      dserror("only implemented for quad4 surface elements.");
    const int cell_numdofpernode = cell_lm.size()/numNode;

    // extract fixpoint values to helper variable
    std::vector<double> myadhesionfixpoints(ldim);
    DRT::UTILS::ExtractMyValues(*cell_adhesion_fixpoints,myadhesionfixpoints,cell_lm);

    // extract displacement values to helper variable
    std::vector<double> myvalues_displ(ldim);
    DRT::UTILS::ExtractMyValues(*dispnp,myvalues_displ,cell_lm);

    // calculate adhesion penalty violation
    std::vector<double> adhesion_violation_nd(numNode*cell_numdofpernode);

    for(int node=0;node<numNode;++node)
      for(int dof=0; dof<cell_numdofpernode;++dof)
        adhesion_violation_nd[node*cell_numdofpernode+dof] =
            myvalues_displ[node*cell_numdofpernode+dof] + (ele->Nodes())[node]->X()[dof] - myadhesionfixpoints[node*cell_numdofpernode+dof];

    // integration points and weights for boundary (!) gp --> quad4
    const DRT::UTILS::IntPointsAndWeights<2> intpoints (DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::quad4>::rule);

    // coordinates of current integration point in face element coordinate system --> QUAD4
    LINALG::Matrix<2,1> xsi(true);
    xsi(0) = intpoints.IP().qxg[iquad][0];
    xsi(1) = intpoints.IP().qxg[iquad][1];

    // shapefunct and derivates of face element in face element coordinate system
    LINALG::Matrix<4, 1> shapefunct;
    DRT::UTILS::shape_function<DRT::Element::quad4>(xsi,shapefunct);

    // get violation at gp from nodal violation
    std::vector<double> adhesion_violation_gp(cell_numdofpernode);

    for (int node=0; node<numNode; ++node)
      for (int dof=0; dof<cell_numdofpernode; ++dof)
        adhesion_violation_gp[dof] += shapefunct(node) * adhesion_violation_nd[node*cell_numdofpernode+dof];

    // extract values to helper variable mytraction
    std::vector<double> mytraction(cell_lm.size());
    DRT::UTILS::ExtractMyValues(*surface_traction,mytraction,cell_lm);

    // nodal traction vector
    LINALG::Matrix<4,1> traction_nd(true);

    // loop over all scatra element nodes to get surface traction at nodes (quad4)
    for (int node=0; node<numNode; node++)
      traction_nd(node,0) = mytraction[node*cell_numdofpernode];

    // interpolate traction at gp from nodal traction
    for (int node=0; node<numNode; node++)
      traction += shapefunct(node) * traction_nd(node,0);

    // write final traction value
    violation = sqrt( adhesion_violation_gp[0]*adhesion_violation_gp[0] +
                      adhesion_violation_gp[1]*adhesion_violation_gp[1] +
                      adhesion_violation_gp[2]*adhesion_violation_gp[2]    );
  }
  else
    violation=0.0;

  return violation;

}


/*-------------------------------------------------------------------------------*
 |  get background element porosity                                  rauch 12/16 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
double DRT::ELEMENTS::ScaTraEleCalcBondReac<distype,probdim>::GetPorosityFromBackgroundEle(
    const DRT::Element*  ele,       //!< the element we are dealing with
    const int            iquad      //!< id of current gauss point
) const
{
  // so far we assume spatial and temporal constant porosity
  Teuchos::RCP<DRT::Discretization> porostructdis = DRT::Problem::Instance()->GetDis("structure");
  const int elegid = porostructdis->ElementRowMap()->GID(0);
  Teuchos::RCP<MAT::Material> poromat = porostructdis->gElement(elegid)->Material();
  double porosity = Teuchos::rcp_dynamic_cast<MAT::StructPoro>(poromat)->Initporosity();

  return porosity;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcBondReac<DRT::Element::line2,1>;
template class DRT::ELEMENTS::ScaTraEleCalcBondReac<DRT::Element::line2,2>;
template class DRT::ELEMENTS::ScaTraEleCalcBondReac<DRT::Element::line2,3>;
template class DRT::ELEMENTS::ScaTraEleCalcBondReac<DRT::Element::line3,1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcBondReac<DRT::Element::tri3,2>;
template class DRT::ELEMENTS::ScaTraEleCalcBondReac<DRT::Element::tri3,3>;
template class DRT::ELEMENTS::ScaTraEleCalcBondReac<DRT::Element::tri6,2>;
template class DRT::ELEMENTS::ScaTraEleCalcBondReac<DRT::Element::quad4,2>;
template class DRT::ELEMENTS::ScaTraEleCalcBondReac<DRT::Element::quad4,3>;
template class DRT::ELEMENTS::ScaTraEleCalcBondReac<DRT::Element::quad9,2>;
template class DRT::ELEMENTS::ScaTraEleCalcBondReac<DRT::Element::nurbs9,2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcBondReac<DRT::Element::hex8,3>;
template class DRT::ELEMENTS::ScaTraEleCalcBondReac<DRT::Element::hex27,3>;
template class DRT::ELEMENTS::ScaTraEleCalcBondReac<DRT::Element::tet4,3>;
template class DRT::ELEMENTS::ScaTraEleCalcBondReac<DRT::Element::tet10,3>;
template class DRT::ELEMENTS::ScaTraEleCalcBondReac<DRT::Element::pyramid5,3>;
