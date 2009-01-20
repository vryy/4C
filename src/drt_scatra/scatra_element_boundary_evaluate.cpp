/*!----------------------------------------------------------------------
\file scatra_element_boundary_evaluate.cpp
\brief

Evaluate boundary conditions for scalar transport problems

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 *----------------------------------------------------------------------*/
#if defined(D_FLUID2) || defined(D_FLUID3)
#ifdef CCADISCRET

#include "scatra_element.H"
#include "scatra_ele_impl.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/convecdiffus.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/ion.H"

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TransportBoundary::Evaluate(ParameterList&         params,
    DRT::Discretization&      discretization,
    vector<int>&              lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  // what kind of action do we have?
  DRT::ELEMENTS::TransportBoundary::ActionType act = TransportBoundary::none;
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action == "calc_condif_flux")
    act = TransportBoundary::calc_condif_flux;
  else if (action == "calc_elch_electrode_kinetics")
    act = TransportBoundary::calc_elch_electrode_kinetics;
  else if (action == "calc_therm_press")
    act = TransportBoundary::calc_therm_press;
  else 
    dserror("Unknown type of action for Transport boundary element");

  // get the material
  RefCountPtr<MAT::Material> mat = parent_->Material();

  _MATERIAL* actmat = NULL;

  if(mat->MaterialType()== m_condif)
    actmat = static_cast<MAT::ConvecDiffus*>(mat.get())->MaterialData();
    else if (mat->MaterialType()== m_matlist)
      actmat = static_cast<MAT::MatList*>(mat.get())->MaterialData();
      else
        dserror("condif or matlist material expected but got type %d", mat->MaterialType());

  switch(act)
  {
  case TransportBoundary::none:
    dserror("action==none");
  case TransportBoundary::calc_condif_flux:
  {
    const int iel = NumNode();

    // get velocity values at the nodes (needed for total flux values)
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
    const int ielparent = parent_->NumNode();
    // get spatial dimenion of the problem
    const int nsd = DRT::UTILS::getDimension(parent_->Shape());
    Epetra_SerialDenseVector evel(nsd*ielparent);
    DRT::UTILS::ExtractMyNodeBasedValues(parent_,evel,velocity,nsd);

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

    // set flag for type of scalar
    string scaltypestr=params.get<string>("problem type");
    const int numdofpernode = parent_->numdofpernode_;
    int numscal = parent_->numdofpernode_;
    bool temperature = false;
    if (scaltypestr =="loma") temperature = true;

    double frt(0.0);
    if (scaltypestr =="elch")
    {
      numscal -= 1; // ELCH case: last dof is for el. potential
      // get parameter F/RT
      frt = params.get<double>("frt");
    }

    // assure, that the values are in the same order as the parent element nodes
    for(int k=0;k<ielparent;++k)
    {
      Node* node = (parent_->Nodes())[k];
      vector<int> dof = discretization.Dof(node);
      for (unsigned int i=0 ; i< dof.size(); ++i)
      {
        if (dof[i]!=lmparent[k*numdofpernode + i])
        {
          cout<<"dof[i]= "<<dof[i]<<"  lmparent[k*numdofpernode + i]="<<lmparent[k*numdofpernode + i]<<endl;
          dserror("Dofs are not in the same order as the element nodes. Implement some resorting!");
        }
      }
    }

    // access control parameter
    string fluxtypestring = params.get<string>("fluxtype","noflux");

    // define vector for normal fluxes
    vector<double> mynormflux(lm.size());

    // node coordinates
    LINALG::SerialDenseMatrix xyze(nsd,iel);
    for(int i=0;i<nsd;i++)
    {
      for(int j=0;j<iel;j++)
      {
        xyze(i,j)=Nodes()[j]->X()[i];
      }
    }

    // determine normal to this element
    std::vector<double> normal(nsd);
    double length = 1.0;

    switch(nsd)
    {
    case 3:
    {
      std::vector<double> dist1(3), dist2(3);
      for (int i=0; i<3; i++)
      {
        dist1[i] = xyze(i,1)-xyze(i,0);
        dist2[i] = xyze(i,2)-xyze(i,0);
      }

      normal[0] = dist1[1]*dist2[2] - dist1[2]*dist2[1];
      normal[1] = dist1[2]*dist2[0] - dist1[0]*dist2[2];
      normal[2] = dist1[0]*dist2[1] - dist1[1]*dist2[0];

      // length of normal to this element
      length = sqrt( normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] );
    }
    break;
    case 2:
    {     normal[0] = xyze(1,1) - xyze(1,0);
    normal[1] = (-1.0)*(xyze(0,1) - xyze(0,0));

    // length of normal to this element
    length = sqrt(normal[0]*normal[0]+normal[1]*normal[1]);
    }
    break;
    default:
      dserror("Illegal number of space dimensions: %d",nsd);
    } // switch(nsd)

    // outward-pointing normal of length 1.0
    for (int i=0; i<nsd; i++)
    {normal[i] = normal[i] / length;}

    // do a loop for systems of transported scalars
    for (int j = 0; j<numscal; ++j)
    {
      // compute fluxes on each node of the parent element
      LINALG::SerialDenseMatrix eflux(3,ielparent);
      DRT::Element* parentele = (DRT::Element*) parent_;
      DRT::ELEMENTS::ScaTraImplInterface::Impl(parentele)->CalculateFluxSerialDense(
          eflux,
          parentele,
          myphinp,
          actmat,
          temperature,
          frt,
          evel,
          fluxtypestring,
          j);

      // handle the result dofs in the right order (compare lm with lmparent)
      int dofcount = 0;
      for (int i=0; i<NumNode(); ++i)
      {
        for(int k = 0; k<ielparent;++k)
        {
          if (lm[i*numdofpernode+j]==lmparent[k*numdofpernode+j]) // dof ids match => assemble this value
          {
            dofcount++;
            // form arithmetic mean of assembled nodal flux vectors
            // => factor is the number of adjacent elements for each node
            double factor = (parent_->Nodes()[k])->NumElement();

            // calculate normal flux at present node
            mynormflux[i] = 0.0;
            for (int l=0; l<nsd; l++)
            {
              mynormflux[i] += eflux(l,k)*normal[l];
            }

            // store normal flux vector for this node
            elevec1[i*numdofpernode+j]+=normal[0]*mynormflux[i]/factor;
            elevec2[i*numdofpernode+j]+=normal[1]*mynormflux[i]/factor;
            elevec3[i*numdofpernode+j]+=normal[2]*mynormflux[i]/factor;
          }
        }
      }
      if (dofcount != NumNode()) dserror("Expected dof for transport boundary element is missing");

      // calculate integral of normal flux
      // => only meaningful for one scalar, for the time being
      // NOTE: add integral value only for elements which are NOT ghosted!
      if(Owner() == discretization.Comm().MyPID())
      {
        NormalFluxIntegral(params,discretization,lm,xyze,mynormflux);
      }

    } // loop over numdofpernode

  }
  break;
  case TransportBoundary::calc_elch_electrode_kinetics:
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
    //const std::string* eltype = cond->Get<std::string>("electrode type");
    const std::string* kinetics = cond->Get<std::string>("kinetic model");
    const int    reactantid = cond->Getint("reactant id");
    double       pot0 = cond->GetDouble("pot0");
    const int    curvenum = cond->Getint("curve");
    const double alphaa = cond->GetDouble("alpha_a");
    const double alphac = cond->GetDouble("alpha_c");
    double       i0 = cond->GetDouble("i0");
    if (i0>0.0) dserror("i0 is positive, ergo not pointing INTO the domain: %f",i0);
    const double frt = params.get<double>("frt"); // = F/RT

    // get control parameter from parameter list
    const bool   iselch = params.get<bool>("iselch");
    const bool   is_stationary = params.get<bool>("using stationary formulation");
    const double time = params.get<double>("total time");
    double       timefac = 1.0;
    double       alphaF  = 1.0;
    if (not is_stationary)
    {
      // One-step-Theta:    timefac = theta*dt
      // BDF2:              timefac = 2/3 * dt
      // generalized-alpha: timefac = (gamma*alpha_F/alpha_M) * dt
      timefac = params.get<double>("time factor");
      alphaF = params.get<double>("alpha_F");
      timefac *= alphaF;
      if (timefac < 0.0) dserror("time factor is negative.");
      // we multiply i0 with timefac. So we do not have to give down this paramater
      // and perform autmatically the multiplication of matrix and rhs with timefac 
      i0 *= timefac;
    }
    // find out whether we shell use a time curve and get the factor
    // this feature can be also used for stationary "pseudo time loops"
    if (curvenum>=0)
    {
      const double curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
      // adjust potential at metal side accordingly
      pot0 *= curvefac;
    }

# if 0
    // print all parameters read from the current condition
    cout<<"electrode type = "<<*eltype<<endl;
    cout<<"sign           = "<<sign<<endl;
    cout<<"kinetic model  = "<<*kinetics<<endl;
    cout<<"reactant id    = "<<reactantid<<endl;
    cout<<"pot0(mod.)     = "<<pot0<<endl;
    cout<<"curvenum       = "<<curvenum<<endl;
    cout<<"alpha_a        = "<<alphaa<<endl;
    cout<<"alpha_c        = "<<alphac<<endl;
    cout<<"i0(mod.)       = "<<i0<<endl;
    cout<<"F/RT           = "<<frt<<endl<<endl;
#endif

    EvaluateElectrodeKinetics(
        elemat1,
        elevec1,
        ephinp,
        actmat,
        reactantid,
        kinetics,
        pot0,
        alphaa,
        alphac,
        i0,
        frt,
        iselch
    );
  }
  break;
  case TransportBoundary::calc_therm_press:
  {
    const int iel = NumNode();

    // get velocity values at the nodes (needed for total flux values)
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
    const int ielparent = parent_->NumNode();
    // get spatial dimenion of the problem
    const int nsd = DRT::UTILS::getDimension(parent_->Shape());
    Epetra_SerialDenseVector evel(nsd*ielparent);
    DRT::UTILS::ExtractMyNodeBasedValues(parent_,evel,velocity,nsd);

    // get actual values of transported scalar
    RefCountPtr<const Epetra_Vector> densnp = discretization.GetState("densnp");
    RefCountPtr<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (densnp==null || phinp==null)
      dserror("Cannot get state vector 'densnp' and/or 'phinp'");

    // we dont know the parent element's lm vector; so we have to build it here
    vector<int> lmparent(ielparent);
    vector<int> lmparentowner;
    parent_->LocationVector(discretization, lmparent, lmparentowner);

    // extract local values from the global vectors for the parent(!) element 
    vector<double> mydensnp(lmparent.size());
    vector<double> myphinp(lmparent.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lmparent);
    DRT::UTILS::ExtractMyValues(*densnp,mydensnp,lmparent);

    // this routine is only for low-mach-number flow
    bool temperature = true;

    // define vector for normal fluxes
    vector<double> mydiffflux(lm.size());
    vector<double> mydivu(lm.size());

    // node coordinates
    LINALG::SerialDenseMatrix xyze(nsd,iel);
    for(int i=0;i<nsd;i++)
    {
      for(int j=0;j<iel;j++)
      {
        xyze(i,j)=Nodes()[j]->X()[i];
      }
    }

    // determine normal to this element
    std::vector<double> normal(nsd);
    double length = 1.0;

    switch(nsd)
    {
    case 3:
    {
      std::vector<double> dist1(3), dist2(3);
      for (int i=0; i<3; i++)
      {
        dist1[i] = xyze(i,1)-xyze(i,0);
        dist2[i] = xyze(i,2)-xyze(i,0);
      }

      normal[0] = dist1[1]*dist2[2] - dist1[2]*dist2[1];
      normal[1] = dist1[2]*dist2[0] - dist1[0]*dist2[2];
      normal[2] = dist1[0]*dist2[1] - dist1[1]*dist2[0];

      // length of normal to this element
      length = sqrt( normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] );
    }
    break;
    case 2:
    {     normal[0] = xyze(1,1) - xyze(1,0);
    normal[1] = (-1.0)*(xyze(0,1) - xyze(0,0));

    // length of normal to this element
    length = sqrt(normal[0]*normal[0]+normal[1]*normal[1]);
    }
    break;
    default:
      dserror("Illegal number of space dimensions: %d",nsd);
    } // switch(nsd)

    // outward-pointing normal of length 1.0
    for (int i=0; i<nsd; i++)
    {normal[i] = normal[i] / length;}

    // compute fluxes on each node of the parent element
    LINALG::SerialDenseMatrix eflux(3,ielparent);
    DRT::Element* parentele = (DRT::Element*) parent_;
    // set some parameters
    double frt=0.0;
    string fluxtypestring("diffusiveflux");
    int j=0;
    DRT::ELEMENTS::ScaTraImplInterface::Impl(parentele)->CalculateFluxSerialDense(
          eflux,
          parentele,
          myphinp,
          actmat,
          temperature,
          frt,
          evel,
          fluxtypestring,
          j);

    for (int i=0; i<NumNode(); ++i)
    {
      for(int k = 0; k<ielparent;++k)
      {
        // calculate normal diffusive flux and velocity div. at present node
        mydiffflux[i] = 0.0;
        mydivu[i]     = 0.0;
        for (int l=0; l<nsd; l++)
        {
          mydiffflux[i] += eflux(l,k)*normal[l];
          mydivu[i]     += (evel[i*nsd+l]/mydensnp[i])*normal[l];
        }
      }
    }

    // calculate integral of normal diffusive flux and velocity divergence
    // NOTE: add integral value only for elements which are NOT ghosted!
    if(Owner() == discretization.Comm().MyPID())
    {
      DifffluxAndDivuIntegral(params,discretization,lm,xyze,mydiffflux,mydivu);
    }
  }
  break;
  default:
    dserror("Unknown type of action for Transport boundary element");
  } // end of switch(act)

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Surface/Line Neumann boundary condition       gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TransportBoundary::EvaluateNeumann(
    ParameterList&            params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    vector<int>&              lm,
    Epetra_SerialDenseVector& elevec1)
{
  const DiscretizationType distype = Shape();
  const int iel = NumNode();

  // switch for the number of space dimensions of the problem:
  const int nsd = DRT::UTILS::getDimension(parent_->Shape());
  if (nsd == 3)
  {
    DRT::UTILS::GaussRule2D  gaussrule = DRT::UTILS::intrule2D_undefined;
    switch(distype)
    {
    case DRT::Element::quad4:
      gaussrule = DRT::UTILS::intrule_quad_4point;
      break;
    case DRT::Element::quad8: case DRT::Element::quad9:
      gaussrule = DRT::UTILS::intrule_quad_9point;
      break;
    case DRT::Element::tri3 :
      gaussrule = DRT::UTILS::intrule_tri_3point;
      break;
    case DRT::Element::tri6:
      gaussrule = DRT::UTILS::intrule_tri_6point;
      break;
    default:
      dserror("shape type unknown!\n");
    }

    const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule);

    // allocate vector for shape functions and matrix for derivatives
    Epetra_SerialDenseVector  funct(iel);
    Epetra_SerialDenseMatrix  deriv(nsd,iel);

    // node coordinates
    Epetra_SerialDenseMatrix  xyze(nsd,iel);

    // the metric tensor and the area of an infinitesimal surface element
    Epetra_SerialDenseMatrix  metrictensor(2,2);
    double                    drs;

    // get node coordinates
    for(int i=0;i<nsd;i++)
    {
      for(int j=0;j<iel;j++)
      {
        xyze(i,j)=Nodes()[j]->X()[i];
      }
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
      // coordinates of the current integration point
      const double e1 = intpoints.qxg[iquad][0];
      const double e2 = intpoints.qxg[iquad][1];

      // shape functions and their first derivatives
      DRT::UTILS::shape_function_2D(funct,e1,e2,distype);
      DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);

      // compute measure tensor for surface element and the infinitesimal
      // area element drs for the integration
      DRT::UTILS::ComputeMetricTensorForSurface(xyze,deriv,metrictensor,&drs);

      // values are multiplied by the product from inf. area element,
      // the gauss weight and the timecurve factor
      const double fac = intpoints.qwgt[iquad] *drs* curvefac;

      // factor given by spatial function
      double functfac = 1.0;
      // determine global coordinates of current Gauss point
      double coordgp[3];
      coordgp[0]=0.0;
      coordgp[1]=0.0;
      coordgp[2]=0.0;
      for (int i = 0; i< nsd; i++)
      {
        for (int j = 0; j < iel; j++)
        {
          coordgp[i] += xyze(i,j) * funct[j];
        }
      }

      int functnum = -1;
      const double* coordgpref = &coordgp[0]; // needed for function evaluation

      numdofpernode_ = parent_->numdofpernode_;

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

  } //if (nsd == 3)
  else if (nsd==2)
  {
    // Gaussian points
    GaussRule1D intrule = intrule1D_undefined;
    switch (distype)
    {
    case line2:
      intrule = intrule_line_2point;
      break;
    case line3:
      intrule = intrule_line_3point;
      break;
    default:
      dserror("unknown number of nodes for gaussrule initialization");
    }

    const IntegrationPoints1D  intpoints(intrule);

    // allocate vector for shape functions and matrix for derivatives
    Epetra_SerialDenseVector funct(iel);
    Epetra_SerialDenseMatrix deriv(1,iel);

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
      const double e1 = intpoints.qxg[iquad][0];

      // get shape functions and derivatives for line element
      DRT::UTILS::shape_function_1D(funct,e1,distype);
      DRT::UTILS::shape_function_1D_deriv1(deriv,e1,distype);

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
        der_par[0] += deriv(0,node)*xye(0,node);
        der_par[1] += deriv(0,node)*xye(1,node);
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
  } //if (nsd == 2)
  else 
    dserror("Illegal number of space dimenions for problem: %d", nsd);

  return 0;
}


/*----------------------------------------------------------------------*
 | evaluate an electrode kinetics boundary condition (private) gjb 09/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TransportBoundary::EvaluateElectrodeKinetics(
    Epetra_SerialDenseMatrix& emat,
    Epetra_SerialDenseVector& erhs,
    const vector<double>&   ephinp,
    struct _MATERIAL*     material,
    const int&               rctid,
    const std::string*    kinetics,
    const double&             pot0,
    const double&           alphaa,
    const double&           alphac,
    const double&               i0,
    const double&              frt,
    const bool&             iselch
)
{
  if ((*kinetics) != "Butler-Volmer")
    dserror("Only Butler-Volmer model allowed. Got: %s",(*kinetics).c_str());

  // some parameters
  const int numdofpernode = parent_->numdofpernode_;
  const int numscal = numdofpernode-1;
  const DiscretizationType distype = Shape();
  const int iel   = NumNode();
  const int nsd = DRT::UTILS::getDimension(parent_->Shape());
  if (nsd != 3)
    dserror ("ElectrodeKinetics boundary condition only available for 3D. Have %dD",nsd);

  //pre-multiplication with 1/(F*z_1)
  double fz = 1.0/96485.3399;

  if (iselch)
  {
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
  }

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
  const IntegrationPoints2D  intpoints(gaussrule);

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
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
    xyze(2,i)=Nodes()[i]->X()[2];
  }

  // concentration values of reactive species at element nodes
  Epetra_SerialDenseVector conreact(iel);

  // el. potential values at element nodes
  Epetra_SerialDenseVector pot(iel);
  if(iselch)
  {
    for (int inode=0; inode< iel;++inode)
    {
      conreact[inode] += ephinp[inode*numdofpernode];
      pot[inode] += ephinp[inode*numdofpernode+numscal];
    }
  }
  else
  {
    for (int inode=0; inode< iel;++inode)
    {
      conreact[inode] = 1.0;
      pot[inode] += ephinp[inode*numdofpernode];
    }
  }

  // concentration of active species at integration point
  static double conint;
  // el. potential at integration point
  static double potint;
  // surface overpotential eta at integration point
  static double eta;
  // a 'working variable'
  static double fac_fz_i0_funct_vi;

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
    eta = (pot0 - potint);

    double gammak = 1.0;
    double pow_conint_gamma_k = pow(conint,gammak);

    if (iselch)
    {
      const double expterm = exp(alphaa*frt*eta)-exp((-alphac)*frt*eta);

      for (int vi=0; vi<iel; ++vi)
      {
        fac_fz_i0_funct_vi = fac*fz*i0*funct[vi];
        // ---------------------matrix
        for (int ui=0; ui<iel; ++ui)
        {
          emat(vi*numdofpernode,ui*numdofpernode) += fac_fz_i0_funct_vi*gammak*pow(conint,(gammak-1.0))*funct[ui]*expterm; 
          emat(vi*numdofpernode,ui*numdofpernode+numscal) += fac_fz_i0_funct_vi*pow_conint_gamma_k*(((-alphaa)*frt*exp(alphaa*frt*eta))+((-alphac)*frt*exp((-alphac)*frt*eta)))*funct[ui];
        }
        // ------------right-hand-side
        erhs[vi*numdofpernode] -= fac_fz_i0_funct_vi*pow_conint_gamma_k*expterm;
      }
    }
    else
    {
#if 1
      // Butler-Volmer kinetics
      const double expterm = exp(alphaa*eta)-exp((-alphac)*eta);
      const double exptermderiv = (((-alphaa)*exp(alphaa*eta))+((-alphac)*exp((-alphac)*eta)));

      for (int vi=0; vi<iel; ++vi)
      {
        const double fac_i0_funct_vi = fac*i0*funct[vi];
        // ---------------------matrix
        for (int ui=0; ui<iel; ++ui)
        {
          emat(vi*numdofpernode,ui*numdofpernode) += fac_i0_funct_vi*exptermderiv*funct[ui];
        }
        // ------------right-hand-side
        erhs[vi*numdofpernode] -= fac_i0_funct_vi*expterm;
      }
#else
      // Tafel kinetics
      const double expterm = -exp((-alphac)*eta);
      const double exptermderiv = alphac*expterm;

      for (int vi=0; vi<iel; ++vi)
      {
        const double fac_i0_funct_vi = fac*i0*funct[vi];
        // ---------------------matrix
        for (int ui=0; ui<iel; ++ui)
        {
          emat(vi*numdofpernode,ui*numdofpernode) += fac_i0_funct_vi*exptermderiv*funct[ui];
        }
        // ------------right-hand-side
        erhs[vi*numdofpernode] -= fac_i0_funct_vi*expterm;
      }
#endif
    } // if iselch

  } // end of loop over integration points gpid

  return;

} // TransportBoundary::EvaluateElectrodeKinetics()


/*----------------------------------------------------------------------*
 | calculate integral of normal flux (private)                  vg 09/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TransportBoundary::NormalFluxIntegral(
    ParameterList&                  params,
    DRT::Discretization&            discretization,
    const vector<int>&              lm,
    const Epetra_SerialDenseMatrix  xyze,
    const vector<double>&           enormflux
)
{
  const int iel   = NumNode();
  const DiscretizationType distype = Shape();
  // number space dimenions of the problem
  const int nsd = DRT::UTILS::getDimension(parent_->Shape());

  // allocate vector for shape functions and matrix for derivatives
  LINALG::SerialDenseVector funct       (iel);
  LINALG::SerialDenseMatrix deriv       (nsd,iel);

  // the metric tensor and the area of an infinitesimal surface element
  LINALG::SerialDenseMatrix metrictensor  (2,2);
  double                    drs;

  // get variable for integral of normal flux
  double normfluxintegral = params.get<double>("normfluxintegral");

  // integrations points and weights
  switch(nsd)
  {
  case 3:
  {
    DRT::UTILS::GaussRule2D  gaussrule = DRT::UTILS::intrule2D_undefined;
    switch(distype)
    {
    case DRT::Element::quad4:
      gaussrule = DRT::UTILS::intrule_quad_4point;
      break;
    case DRT::Element::quad8: case DRT::Element::quad9:
      gaussrule = DRT::UTILS::intrule_quad_9point;
      break;
    case DRT::Element::tri3 :
      gaussrule = DRT::UTILS::intrule_tri_3point;
      break;
    case DRT::Element::tri6:
      gaussrule = DRT::UTILS::intrule_tri_6point;
      break;
    default:
      dserror("shape type unknown!\n");
    }

    const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule);
    for (int gpid=0; gpid<intpoints.nquad; gpid++)
    {
      const double e0 = intpoints.qxg[gpid][0];
      const double e1 = intpoints.qxg[gpid][1];

      // get shape functions and derivatives in the plane of the element
      DRT::UTILS::shape_function_2D(funct, e0, e1, distype);
      DRT::UTILS::shape_function_2D_deriv1(deriv, e0, e1, distype);

      // compute measure tensor for surface element and the infinitesimal
      // area element drs for the integration
      DRT::UTILS::ComputeMetricTensorForSurface(xyze,deriv,metrictensor,&drs);

      const double fac = intpoints.qwgt[gpid] * drs;

      // compute integral of normal flux
      for (int node=0;node<iel;++node)
      {
        normfluxintegral += funct[node] * enormflux[node] * fac;
      }
    } // loop over integration points
  }
  break;
  case 2:
  {
    GaussRule1D intrule = intrule1D_undefined;
    switch (distype)
    {
    case line2:
      intrule = intrule_line_2point;
      break;
    case line3:
      intrule = intrule_line_3point;
      break;
    default:
      dserror("unknown number of nodes for gaussrule initialization");
    }

    const IntegrationPoints1D  intpoints(intrule);
    for (int iquad=0; iquad<intpoints.nquad; iquad++)
    {
      const double e1 = intpoints.qxg[iquad][0];

      // get shape functions and derivatives in the plane of the element
      DRT::UTILS::shape_function_1D(funct,e1,distype);
      DRT::UTILS::shape_function_1D_deriv1(deriv,e1,distype);

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
        der_par[0] += deriv(0,node)*xyze(0,node);
        der_par[1] += deriv(0,node)*xyze(1,node);
      }

      // compute infinitesimal line element dr for integration along the line
      const double dr = sqrt(der_par[0]*der_par[0]+der_par[1]*der_par[1]);

      const double fac = intpoints.qwgt[iquad] * dr;

      // compute integral of normal flux
      for (int node=0;node<iel;++node)
      {
          normfluxintegral += funct[node] * enormflux[node] *fac;
      }
    } // loop over integration points
  }
  break;
  default:
    dserror("Illegal number of space dimenions: %d",nsd);
  } // switch(nsd)

  // add contribution to the global value
  params.set<double>("normfluxintegral",normfluxintegral);

}//DRT::ELEMENTS::TransportBoundary::NormalFluxIntegral


/*----------------------------------------------------------------------*
 | calculate integral of normal diffusive flux + velocity div.  vg 09/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TransportBoundary::DifffluxAndDivuIntegral(
    ParameterList&                  params,
    DRT::Discretization&            discretization,
    const vector<int>&              lm,
    const Epetra_SerialDenseMatrix  xyze,
    const vector<double>&           ediffflux,
    const vector<double>&           edivu
)
{
  const int iel = NumNode();
  const DiscretizationType distype = Shape();
  // number space dimenions of the problem
  const int nsd = DRT::UTILS::getDimension(parent_->Shape());

  // allocate vector for shape functions and matrix for derivatives
  LINALG::SerialDenseVector funct(iel);
  LINALG::SerialDenseMatrix deriv(nsd,iel);

  // the metric tensor and the area of an infinitesimal surface element
  LINALG::SerialDenseMatrix metrictensor(2,2);
  double                    drs;

  // get variables for integrals of normal diffusive flux and velocity div.
  double difffluxintegral = params.get<double>("diffusive-flux integral");
  double divuintegral     = params.get<double>("velocity-divergence integral");

  // integrations points and weights
  switch(nsd)
  {
  case 3:
  {
    DRT::UTILS::GaussRule2D  gaussrule = DRT::UTILS::intrule2D_undefined;
    switch(distype)
    {
    case DRT::Element::quad4:
      gaussrule = DRT::UTILS::intrule_quad_4point;
      break;
    case DRT::Element::quad8: case DRT::Element::quad9:
      gaussrule = DRT::UTILS::intrule_quad_9point;
      break;
    case DRT::Element::tri3 :
      gaussrule = DRT::UTILS::intrule_tri_3point;
      break;
    case DRT::Element::tri6:
      gaussrule = DRT::UTILS::intrule_tri_6point;
      break;
    default:
      dserror("shape type unknown!\n");
    }

    const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule);
    for (int gpid=0; gpid<intpoints.nquad; gpid++)
    {
      const double e0 = intpoints.qxg[gpid][0];
      const double e1 = intpoints.qxg[gpid][1];

      // get shape functions and derivatives in the plane of the element
      DRT::UTILS::shape_function_2D(funct, e0, e1, distype);
      DRT::UTILS::shape_function_2D_deriv1(deriv, e0, e1, distype);

      // compute measure tensor for surface element and the infinitesimal
      // area element drs for the integration
      DRT::UTILS::ComputeMetricTensorForSurface(xyze,deriv,metrictensor,&drs);

      const double fac = intpoints.qwgt[gpid] * drs;

      // compute integral of normal flux
      for (int node=0;node<iel;++node)
      {
        difffluxintegral += funct[node] * ediffflux[node] * fac;
        divuintegral     += funct[node] * edivu[node] * fac;
      }
    } // loop over integration points
  }
  break;
  case 2:
  {
    GaussRule1D intrule = intrule1D_undefined;
    switch (distype)
    {
    case line2:
      intrule = intrule_line_2point;
      break;
    case line3:
      intrule = intrule_line_3point;
      break;
    default:
      dserror("unknown number of nodes for gaussrule initialization");
    }

    const IntegrationPoints1D  intpoints(intrule);
    for (int iquad=0; iquad<intpoints.nquad; iquad++)
    {
      const double e1 = intpoints.qxg[iquad][0];

      // get shape functions and derivatives in the plane of the element
      DRT::UTILS::shape_function_1D(funct,e1,distype);
      DRT::UTILS::shape_function_1D_deriv1(deriv,e1,distype);

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
        der_par[0] += deriv(0,node)*xyze(0,node);
        der_par[1] += deriv(0,node)*xyze(1,node);
      }

      // compute infinitesimal line element dr for integration along the line
      const double dr = sqrt(der_par[0]*der_par[0]+der_par[1]*der_par[1]);

      const double fac = intpoints.qwgt[iquad] * dr;

      // compute integral of normal flux
      for (int node=0;node<iel;++node)
      {
        difffluxintegral += funct[node] * ediffflux[node] * fac;
        divuintegral     += funct[node] * edivu[node] * fac;
      }
    } // loop over integration points
  }
  break;
  default:
    dserror("Illegal number of space dimenions: %d",nsd);
  } // switch(nsd)

  // add contributions to the global values
  params.set<double>("diffusive-flux integral",difffluxintegral);
  params.set<double>("velocity-divergence integral",divuintegral);

}//DRT::ELEMENTS::TransportBoundary::NormalFluxIntegral


#endif  // #ifdef CCADISCRET
#endif  // #if defined(D_FLUID2) || defined(D_FLUID3)
