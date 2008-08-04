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
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_mat/convecdiffus.H"


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
  else dserror("Unknown type of action for Condif3_Surface");

  // get the material
  RefCountPtr<MAT::Material> mat = parent_->Material();
  if (mat->MaterialType() != m_condif)
    dserror("Material law is not a condif element");

  MATERIAL* actmat = NULL;
  if(mat->MaterialType()== m_condif)
    actmat = static_cast<MAT::ConvecDiffus*>(mat.get())->MaterialData();
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
    parent_->DRT::ELEMENTS::Condif3::ExtractMyNodeBasedValues(evel,velocity);

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


#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
