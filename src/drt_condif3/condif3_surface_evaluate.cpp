/*!----------------------------------------------------------------------
\file condif3_surface_evaluate.cpp
\brief

Evaluate surface conditions

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>


 *----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "condif3.H"
#include "condif3_utils.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/convecdiffus.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/ion.H"

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             gjb 06/08 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Condif3Surface::Evaluate(ParameterList&            params,
                                            DRT::Discretization&      discretization,
                                            vector<int>&              lm,
                                            Epetra_SerialDenseMatrix& elemat1,
                                            Epetra_SerialDenseMatrix& elemat2,
                                            Epetra_SerialDenseVector& elevec1,
                                            Epetra_SerialDenseVector& elevec2,
                                            Epetra_SerialDenseVector& elevec3)
{
  // what kind of action do we have?
  DRT::ELEMENTS::Condif3Surface::ActionType act = Condif3Surface::none;
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action == "calc_condif_flux")
    act = Condif3Surface::calc_condif_flux;
  else if (action == "calc_elch_electrode_kinetics")
    act = Condif3Surface::calc_elch_electrode_kinetics;
  else dserror("Unknown type of action for Condif3_Surface");

  // get the material
  RefCountPtr<MAT::Material> mat = parent_->Material();

  MATERIAL* actmat = NULL;

  if(mat->MaterialType()== m_condif)
    actmat = static_cast<MAT::ConvecDiffus*>(mat.get())->MaterialData();
  else if (mat->MaterialType()== m_matlist)
    actmat = static_cast<MAT::MatList*>(mat.get())->MaterialData();
  else
    dserror("condif material expected but got type %d", mat->MaterialType());

  switch(act)
  {
  case Condif3Surface::none:
    dserror("action=none");
  case Condif3Surface::calc_condif_flux:
  {
    // get velocity values at the nodes (needed for total flux values)
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
    const int ielparent = parent_->NumNode();
    const int nsd=3;
    Epetra_SerialDenseVector evel(nsd*ielparent);
    DRT::UTILS::ExtractMyNodeBasedValues(parent_,evel,velocity);

    // get actual values of transported scalar
    RefCountPtr<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==null) dserror("Cannot get state vector 'phinp'");

    // we dont know the parent element's lm vector; so we have to build it here
    vector<int> lmparent(ielparent);
    vector<int> lmparentowner;
    parent_->LocationVector(discretization, lmparent, lmparentowner);

    // extract local values from the global vector for the parent(!) element 
    vector<double> myphinp(lmparent.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lmparent);

    // assure, that the values are in the same order as the parent element nodes
    for(int k=0;k<ielparent;++k)
    {
      Node* node = (parent_->Nodes())[k];
      vector<int> dof = discretization.Dof(node);
        // up to now, there's only one dof per node
        if (dof[0]!=lmparent[k])
        { cout<<"dof[0]= "<<dof[0]<<"  lmparent[j]="<<lmparent[k]<<endl;
          dserror("Dofs are not in the same order as the element nodes. Implement some resorting!");
        }
    }

    // access control parameter
    Condif3::FluxType fluxtype;
    string fluxtypestring = params.get<string>("fluxtype","noflux");
    if (fluxtypestring == "totalflux")
      fluxtype = Condif3::totalflux;
    else if (fluxtypestring == "diffusiveflux")
      fluxtype = Condif3::diffusiveflux;
    else
      fluxtype=Condif3::noflux;  //default value

    // do a loop for systems of transported scalars
    const int numdofpernode = parent_->numdofpernode_;
    for (int j = 0; j<numdofpernode; ++j)
    {
      // compute fluxes on each node of the parent element
      Epetra_SerialDenseMatrix eflux = parent_->CalculateFlux(myphinp,actmat,evel,fluxtype,j);

      // handle the result dofs in the right order (compare lm with lmparent)
      int dofcount = 0;
      for (int i=0; i<NumNode(); ++i)
      {
        for(int k = 0; k<ielparent;++k)
        {
          if (lm[i]==lmparent[k]) // dof ids match => assemble this value
          {
            dofcount++;
            // form arithmetic mean of assembled nodal flux vectors
            // => factor is the number of adjacent elements for each node
            double factor = (parent_->Nodes()[k])->NumElement();
            // here, we rely on the right order (see above!)
            //Node* node = (parent_->Nodes())[k];
            //vector<int> dof = discretization.Dof(node);
            //if (dof[0]!=lm[i]) dserror("mismatch");
            elevec1[i*numdofpernode+j]+=eflux(0,k)/factor;
            elevec2[i*numdofpernode+j]+=eflux(1,k)/factor;
            elevec3[i*numdofpernode+j]+=eflux(2,k)/factor;
          }
        }
      }
      if (dofcount != NumNode()) dserror("Expected dof for surface element is missing");

    } // loop over numdofpernode

  }
  break;
  case Condif3Surface::calc_elch_electrode_kinetics:
  {
    // get actual values of transported scalar
     RefCountPtr<const Epetra_Vector> phinp = discretization.GetState("phinp");
     if (phinp==null) dserror("Cannot get state vector 'phinp'");
     // extract local values from the global vector
     vector<double> ephinp(lm.size());
     DRT::UTILS::ExtractMyValues(*phinp,ephinp,lm);

    // get current condition
     Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition> >("condition");
     if (cond == Teuchos::null) dserror("Cannot access condition 'ElectrodeKinetics'");

    // access parameters of the condition
    double sign(1.0);
    const std::string* eltype = cond->Get<std::string>("electrode type");
    if ((*eltype)== "anode") sign = -1.0;
    const std::string* kinetics = cond->Get<std::string>("kinetic model");
    const int    reactantid = cond->Getint("reactant id");
    const double pot0 = cond->GetDouble("pot0");
    const double alphaa = cond->GetDouble("alpha_a");
    const double alphac = cond->GetDouble("alpha_c");
    const double i0 = cond->GetDouble("i0");
    const double frt = params.get<double>("frt"); // = F/RT

# if 0
    // print all parameters read from the current condition
    cout<<"electrode type = "<<*eltype<<endl;
    cout<<"sign           = "<<sign<<endl;
    cout<<"kinetic model  = "<<*kinetics<<endl;
    cout<<"reactant id    = "<<reactantid<<endl;
    cout<<"pot0           = "<<pot0<<endl;
    cout<<"alpha_a        = "<<alphaa<<endl;
    cout<<"alpha_c        = "<<alphac<<endl;
    cout<<"i0             = "<<i0<<endl<<endl;
    cout<<"F/RT           = "<<frt<<endl<<endl;
#endif

    EvaluateElectrodeKinetics(
        elemat1,
        elevec1,
        ephinp,
        actmat,
        sign,
        reactantid,
        kinetics,
        pot0,
        alphaa,
        alphac,
        i0,
        frt
        );
  }
  break;
  default:
    dserror("Unknown type of action for Condif3_Surface");
  } // end of switch(act)

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)   gjb 06/08 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Condif3Surface::EvaluateNeumann(
    ParameterList&            params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    vector<int>&              lm,
    Epetra_SerialDenseVector& elevec1)
{
  dserror("EvaluateNeumann not implemented.");
  return 0;
}


/*----------------------------------------------------------------------*
 | evaluate an electrode kinetics boundary condition (private) gjb 09/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif3Surface::EvaluateElectrodeKinetics(
    Epetra_SerialDenseMatrix& emat,
    Epetra_SerialDenseVector& erhs,
    const vector<double>&   ephinp,
    struct _MATERIAL*     material,
    const double&             sign,
    const int&               rctid,
    const std::string*    kinetics,
    const double&             pot0,
    const double&           alphaa,
    const double&           alphac,
    const double&               i0,
    const double&              frt
)
{
  if ((*kinetics) != "Butler-Volmer")
    dserror("Only Butler-Volmer model allowed. Got: %s",(*kinetics).c_str());

  // some parameters
  const int numdofpernode = parent_->numdofpernode_;
  const int numscal = numdofpernode-1;
  const DiscretizationType distype = this->Shape();
  const int iel   = this->NumNode();

  //pre-multiplication with 1/(F*z_1)
  double fz = 1.0/96485.3399;
  // get valence of the single(!) reactant
  if (material->mattyp == m_matlist)
  {
    if (material->m.matlist->matids[0] != rctid) 
      dserror("active species is not first scalar in material list!");
    // the active species is the FIRST material in the material list. ALWAYS!
    const _MATERIAL& singlemat =  DRT::Problem::Instance()->Material(rctid-1);
    if (singlemat.mattyp == m_ion)
      fz = fz/singlemat.m.ion->valence;
    else
      dserror("single material type is not 'ion'");
  }
  else
    dserror("material type is not a 'matlist' material");

  // Gaussian points
  GaussRule2D  gaussrule = intrule2D_undefined;
  switch(distype)
  {
  case quad4:
    gaussrule = intrule_quad_4point;
    break;
  case quad8: case quad9:
    gaussrule = intrule_quad_9point;
    break;
  case tri3 :
    gaussrule = intrule_tri_3point;
    break;
  case tri6:
    gaussrule = intrule_tri_6point;
    break;
  default:
    dserror("shape type unknown!\n");
  }
  const IntegrationPoints2D  intpoints = getIntegrationPoints2D(gaussrule);

  // allocate vector for shape functions and matrix for derivatives
  Epetra_SerialDenseVector  funct(iel);
  Epetra_SerialDenseMatrix  deriv(2,iel);

  // node coordinates
  Epetra_SerialDenseMatrix  xyze(3,iel);

  // the metric tensor and the area of an infinitesimal surface element
  Epetra_SerialDenseMatrix  metrictensor(2,2);
  double                    drs;

  // get node coordinates
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=this->Nodes()[i]->X()[0];
    xyze(1,i)=this->Nodes()[i]->X()[1];
    xyze(2,i)=this->Nodes()[i]->X()[2];
  }

  // concentration values of reactive species at element nodes
  Epetra_SerialDenseVector conreact(iel);

  // el. potential values at element nodes
  Epetra_SerialDenseVector pot(iel);
  for (int inode=0; inode< iel;++inode)
  {
    conreact[inode] += ephinp[inode*numdofpernode];
    pot[inode] += ephinp[inode*numdofpernode+numscal];
  }

  // concentration of active species at integration point
  static double conint;
  // el. potential at integration point
  static double potint;
  // surface overpotential eta at integration point
  static double eta;
  // a 'working variable'
  static double fac_fz_sign_i0_funct_vi;

  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
   *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.nquad; gpid++)
  {
    // get coordinates of integration point
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get shape functions and derivatives in the plane of the element
    shape_function_2D(funct, e0, e1, distype);
    shape_function_2D_deriv1(deriv, e0, e1, distype);

    // compute measure tensor for surface element and the infinitesimal
    // area element drs for the integration
    DRT::UTILS::ComputeMetricTensorForSurface(xyze,deriv,metrictensor,&drs);

    // values are multiplied by the product from inf. area element,
    // the gauss weight, the timecurve factor and the constant
    // belonging to the time integration algorithm (theta*dt for
    // one step theta, 2/3 for bdf with dt const.)
    const double fac = intpoints.qwgt[gpid] * drs;

    // elch-specific values at integration point:
    conint = 0.0;
    potint = 0.0;
    for (int node=0;node<iel;++node)
    {
      conint += funct[node]*conreact[node];
      potint += funct[node]*pot[node];
    }

    // anode:   eta= phi0 - phi
    // cathode: eta= phi - phi0
    eta = sign*(potint-pot0); 

    const double expterm = exp(alphaa*frt*eta)-exp((-alphac)*frt*eta);

    for (int vi=0; vi<iel; ++vi)
    {
      fac_fz_sign_i0_funct_vi = fac*fz*sign*i0*funct[vi];
      // ---------------------matrix
      for (int ui=0; ui<iel; ++ui)
      {
        emat(vi*numdofpernode,ui*numdofpernode) += fac_fz_sign_i0_funct_vi*funct[ui]*expterm; 
        emat(vi*numdofpernode,ui*numdofpernode+numscal) += fac_fz_sign_i0_funct_vi*conint*((alphaa*frt*sign*exp(alphaa*frt*eta))+(alphac*frt*sign*exp((-alphac)*frt*eta)))*funct[ui];
      }
      // ------------right-hand-side
      erhs[vi*numdofpernode] -= fac_fz_sign_i0_funct_vi*conint*expterm;
    }

  } // end of loop over integration points gpid

  return;

} // Condif3Surface::EvaluateElectrodeKinetics()


#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
