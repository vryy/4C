/*!----------------------------------------------------------------------
\file condif2_line_evaluate.cpp
\brief

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET

#include "condif2.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/convecdiffus.H"

using namespace DRT::UTILS;


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                               vg 08/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Condif2Line::Evaluate(ParameterList&            params,
                                         DRT::Discretization&      discretization,
                                         vector<int>&              lm,
                                         Epetra_SerialDenseMatrix& elemat1,
                                         Epetra_SerialDenseMatrix& elemat2,
                                         Epetra_SerialDenseVector& elevec1,
                                         Epetra_SerialDenseVector& elevec2,
                                         Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Condif2Line::ActionType act = Condif2Line::none;
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action == "calc_condif_flux")
    act = Condif2Line::calc_condif_flux;
  else dserror("Unknown type of action for Condif2Line");

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
    case Condif2Line::none:
      dserror("action=none");
    case Condif2Line::calc_condif_flux:
    {
      const int iel   = this->NumNode();

      // get velocity values at the nodes (needed for total flux values)
      const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
      const int ielparent = parent_->NumNode();
      const int nsd=2;
      Epetra_SerialDenseVector evel(nsd*ielparent);
      DRT::UTILS::ExtractMyNodeBasedValues2D(parent_,evel,velocity);

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
          {
            cout<<"dof[0]= "<<dof[0]<<"  lmparent[j]="<<lmparent[k]<<endl;
            dserror("Dofs are not in the same order as the element nodes. Implement some resorting!");
          }
      }

      // access control parameter
      Condif2::FluxType fluxtype;
      string fluxtypestring = params.get<string>("fluxtype","noflux");
      if (fluxtypestring == "totalflux")
        fluxtype = Condif2::totalflux;
      else if (fluxtypestring == "diffusiveflux")
        fluxtype = Condif2::diffusiveflux;
      else
        fluxtype=Condif2::noflux;  //default value

      // set flag for type of scalar
      string scaltypestr=params.get<string>("problem type");
      bool temperature = false;
      if (scaltypestr =="loma") temperature = true;

      // define vector for normal fluxes
      vector<double> mynormflux(lm.size());

      // node coordinates
      LINALG::SerialDenseMatrix xye(2,iel);
      for(int i=0;i<iel;i++)
      {
        xye(0,i)=this->Nodes()[i]->X()[0];
        xye(1,i)=this->Nodes()[i]->X()[1];
      }

      // determine normal to this element
      std::vector<double> normal(2);
      double length;

      normal[0] = xye(1,1) - xye(1,0);
      normal[1] = (-1.0)*(xye(0,1) - xye(0,0));

      // length of normal to this element
      length = sqrt(normal[0]*normal[0]+normal[1]*normal[1]);

      // outward-pointing normal of length 1.0
      for (int i=0; i<2; i++) normal[i] = normal[i] / length;

      // do a loop for systems of transported scalars
      const int numdofpernode = parent_->numdofpernode_;
      for (int j = 0; j<numdofpernode; ++j)
      {
        // compute fluxes on each node of the parent element
        Epetra_SerialDenseMatrix eflux = parent_->CalculateFlux(myphinp,actmat,temperature,evel,fluxtype,j);

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

              // calculate normal flux at present node
              mynormflux[i] = abs((eflux(0,k)*normal[0] + eflux(1,k)*normal[1])/factor);

              // calculate integral of normal flux
              // => only meaningful for one scalar, for the time being
              NormalFluxIntegral(params,discretization,lm,xye,mynormflux);

              // normal flux value is stored in elevec1, other elevecs set to 0.0
              elevec1[i*numdofpernode+j]+=mynormflux[i];
              elevec2[i*numdofpernode+j]+=0.0;
              elevec3[i*numdofpernode+j]+=0.0;
            }
          }
        }
        if (dofcount != NumNode()) dserror("Expected dof for surface element is missing");

      } // loop over numdofpernode

    }
    break;
    default:
      dserror("Unknown type of action for Condif2Line");
  } // end of switch(act)

  return 0;

} // DRT::ELEMENTS::Condif2Line::Evaluate


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif2Line::NormalFluxIntegral(
   ParameterList&                  params,
   DRT::Discretization&            discretization,
   const vector<int>&              lm,
   const Epetra_SerialDenseMatrix  xye,
   const vector<double>&           enormflux
)
{
  const int iel   = this->NumNode();
  const DiscretizationType distype = this->Shape();
  // allocate vector for shape functions and matrix for derivatives
  LINALG::SerialDenseVector funct       (iel);
  LINALG::SerialDenseVector deriv       (iel);

  double normfluxintegral = params.get<double>("normfluxintegral");

  const GaussRule1D gaussrule = getOptimalGaussrule(distype);

  const IntegrationPoints1D  intpoints = getIntegrationPoints1D(gaussrule);
  for (int iquad=0; iquad<intpoints.nquad; iquad++)
  {
    const double e1 = intpoints.qxg[iquad];

    // get shape functions and derivatives in the plane of the element
    DRT::UTILS::shape_function_1D(funct,e1,distype);
    // explicit determination of first derivatives, since respective function
    //DRT::UTILS::shape_function_1D_deriv1(deriv,e1,distype);
    // provides mismatch due to matrix/vector inconsistency
    switch (distype)
    {
      case DRT::Element::point1:
      {
        deriv(0) = 0.0;
        break;
      }
      case DRT::Element::line2:
      {
        deriv(0)= -0.5;
        deriv(1)= 0.5;
        break;
      }
      case DRT::Element::line3:
      {
        deriv(0)= e1 - 0.5;
        deriv(1)= e1 + 0.5;
        deriv(2)= -2.0*e1;
        break;
      }
      default:
         dserror("distype unknown\n");
    }

    // The Jacobian is computed using the formula
    //
    //            +-----
    //   dx_j(u)   \      dN_k(r)
    //   -------  = +     ------- * (x_j)_k
    //     du      /       du         |
    //            +-----    |         |
    //            node k    |         |
    //                  derivative    |
    //                   of shape     |
    //                   function     |
    //                           component of
    //                          node coordinate
    //
    LINALG::SerialDenseVector der_par(2);
    der_par[0]=0.0;
    der_par[1]=0.0;
    for (int node=0;node<iel;++node)
    {
      der_par[0] += deriv[node]*xye(0,node);
      der_par[1] += deriv[node]*xye(1,node);
    }

    // compute infinitesimal line element dr for integration along the line
    const double dr = sqrt(der_par[0]*der_par[0]+der_par[1]*der_par[1]);

    const double fac = intpoints.qwgt[iquad] * dr;

    // compute integral of normal flux
    for (int node=0;node<iel;++node)
    {
      for(int dim=0;dim<2;dim++)
      {
        normfluxintegral += funct[node] * enormflux[node] *fac;
      }
    }
  }

  params.set<double>("normfluxintegral",normfluxintegral);

}//DRT::ELEMENTS::Condif2Line::NormalFluxIntegral


/*----------------------------------------------------------------------*
 |  Integrate a Line Neumann boundary condition (public)        vg 09/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Condif2Line::EvaluateNeumann(
    ParameterList& params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    vector<int>&              lm,
    Epetra_SerialDenseVector& elevec1)
{
  const DiscretizationType distype = this->Shape();

  // set number of nodes
  const int iel   = this->NumNode();

  // Gaussian points
  const GaussRule1D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints1D  intpoints = getIntegrationPoints1D(gaussrule);

  // allocate vector for shape functions and matrix for derivatives
  Epetra_SerialDenseVector funct(iel);
  Epetra_SerialDenseVector deriv(iel);

  // node coordinates
  Epetra_SerialDenseMatrix xye(2,iel);
  for(int i=0;i<iel;i++)
  {
    xye(0,i)=this->Nodes()[i]->X()[0];
    xye(1,i)=this->Nodes()[i]->X()[1];
  }

  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve  = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

  // get values, switches and spatial functions from the condition
  // (assumed to be constant on element boundary)
  const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val"  );
  const vector<int>*    func  = condition.Get<vector<int> >   ("funct");

  // integration loop
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    const double e1 = intpoints.qxg[iquad];

    // get shape functions and derivatives for line element
    DRT::UTILS::shape_function_1D(funct,e1,distype);
    // explicit determination of first derivatives, since respective function
    //DRT::UTILS::shape_function_1D_deriv1(deriv,e1,distype);
    // provides mismatch due to matrix/vector inconsistency
    switch (distype)
    {
      case DRT::Element::point1:
      {
        deriv(0) = 0.0;
        break;
      }
      case DRT::Element::line2:
      {
        deriv(0)= -0.5;
        deriv(1)= 0.5;
        break;
      }
      case DRT::Element::line3:
      {
        deriv(0)= e1 - 0.5;
        deriv(1)= e1 + 0.5;
        deriv(2)= -2.0*e1;
        break;
      }
      default:
         dserror("distype unknown\n");
      }

    // The Jacobian is computed using the formula
    //
    //            +-----
    //   dx_j(u)   \      dN_k(r)
    //   -------  = +     ------- * (x_j)_k
    //     du      /       du         |
    //            +-----    |         |
    //            node k    |         |
    //                  derivative    |
    //                   of shape     |
    //                   function     |
    //                           component of
    //                          node coordinate
    //
    LINALG::SerialDenseVector der_par(2);
    der_par[0]=0.0;
    der_par[1]=0.0;
    for (int node=0;node<iel;++node)
    {
      der_par[0] += deriv[node]*xye(0,node);
      der_par[1] += deriv[node]*xye(1,node);
    }

    // compute infinitesimal line element dr for integration along the line
    const double dr = sqrt(der_par[0]*der_par[0]+der_par[1]*der_par[1]);

    // values are multiplied by the product from inf. area element,
    // the gauss weight and the timecurve factor
    const double fac = intpoints.qwgt[iquad] *dr* curvefac;

    // factor given by spatial function
    double functfac = 1.0;
    // determine coordinates of current Gauss point
    double coordgp[2];
    coordgp[0]=0.0;
    coordgp[1]=0.0;
    for (int i = 0; i< iel; i++)
    {
      coordgp[0]+=xye(0,i)*funct(i);
      coordgp[1]+=xye(1,i)*funct(i);
    }

    int functnum = -1;
    const double* coordgpref = &coordgp[0]; // needed for function evaluation

    // for the time being, only for single scalar transport equations
    numdofpernode_ = 1;

    for (int node=0;node<iel;++node)
    {
      for(int dof=0;dof<numdofpernode_;dof++)
      {
        // factor given by spatial function
        if (func) functnum = (*func)[dof];
        {
          if (functnum>0)
          {
            // evaluate function at current gauss point
            functfac = DRT::UTILS::FunctionManager::Instance().Funct(functnum-1).Evaluate(dof,coordgpref);
          }
          else
            functfac = 1.0;
        }

        elevec1[node]+= funct(node)*(*onoff)[dof]*(*val)[dof]*fac*functfac;
      }
    }
  } //end of loop over integration points

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
GaussRule1D DRT::ELEMENTS::Condif2Line::getOptimalGaussrule(const DiscretizationType& distype)
{
  GaussRule1D rule = intrule1D_undefined;
  switch (distype)
    {
        case line2:
          rule = intrule_line_2point;
          break;
        case line3:
          rule = intrule_line_3point;
          break;
    default:
    dserror("unknown number of nodes for gaussrule initialization");
    }
  return rule;
}


#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID2
