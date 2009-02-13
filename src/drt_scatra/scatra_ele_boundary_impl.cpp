/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_impl.cpp

\brief Internal implementation of scalar transport boundary elements

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 */
/*----------------------------------------------------------------------*/

#if defined(D_FLUID2) || defined(D_FLUID3)
#ifdef CCADISCRET

#include "scatra_ele_boundary_impl.H"
#include "scatra_ele_impl.H"
#include "../drt_mat/convecdiffus.H"
#include "../drt_mat/sutherland_condif.H"
#include "../drt_mat/ion.H"
#include "../drt_mat/matlist.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_geometry/position_array.H"
#include "../drt_lib/linalg_serialdensematrix.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_lib/drt_function.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraBoundaryImplInterface* DRT::ELEMENTS::ScaTraBoundaryImplInterface::Impl(DRT::Element* ele)
{
  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the computations
  const int numdofpernode = ele->NumDofPerNode(*(ele->Nodes()[0]));
  int numscal = numdofpernode;
  if (DRT::Problem::Instance()->ProblemType() == "elch")
    numscal -= 1;

  switch (ele->Shape())
  {
  case DRT::Element::quad4:
  {
    static ScaTraBoundaryImpl<DRT::Element::quad4>* cp4;
    if (cp4==NULL)
      cp4 = new ScaTraBoundaryImpl<DRT::Element::quad4>(numdofpernode,numscal);
      return cp4;
  }
  /*  case DRT::Element::quad8:
  {
    static ScaTraBoundaryImpl<DRT::Element::quad8>* cp8;
    if (cp8==NULL)
      cp8 = new ScaTraImpl<DRT::Element::quad8>(numdofpernode,numscal);
    return cp8;
  }
  case DRT::Element::quad9:
  {
    static ScaTraBoundaryImpl<DRT::Element::quad9>* cp9;
    if (cp9==NULL)
      cp9 = new ScaTraImpl<DRT::Element::quad9>(numdofpernode,numscal);
    return cp9;
  }*/
  case DRT::Element::tri3:
  {
    static ScaTraBoundaryImpl<DRT::Element::tri3>* cp3;
    if (cp3==NULL)
      cp3 = new ScaTraBoundaryImpl<DRT::Element::tri3>(numdofpernode,numscal);
      return cp3;
  }
  /*  case DRT::Element::tri6:
  {
    static ScaTraBoundaryImpl<DRT::Element::tri6>* cp6;
    if (cp6==NULL)
      cp6 = new ScaTraBoundaryImpl<DRT::Element::tri6>(numdofpernode,numscal);
    return cp6;
  }*/
  case DRT::Element::line2:
  {
    static ScaTraBoundaryImpl<DRT::Element::line2>* cl2;
    if (cl2==NULL)
      cl2 = new ScaTraBoundaryImpl<DRT::Element::line2>(numdofpernode,numscal);
      return cl2;
  }/*
  case DRT::Element::line3:
  {
    static ScaTraBoundaryImpl<DRT::Element::line3>* cl3;
    if (cl3==NULL)
      cl3 = new ScaTraBoundaryImpl<DRT::Element::line3>(numdofpernode,numscal);
    return cl3;
  }*/
  default:
    dserror("Shape %d (%d nodes) not supported", ele->Shape(), ele->NumNode());
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::ScaTraBoundaryImpl
(int numdofpernode, 
    int numscal)
    : numdofpernode_(numdofpernode),
    numscal_(numscal),
    xyze_(true),  // initialize to zero
    bodyforce_(numdofpernode_,0),
    diffus_(numscal_,0),
    valence_(numscal_,0),
    shcacp_(0.0),
    xsi_(true),
    funct_(true),
    deriv_(true),
    derxy_(true),
    normal_(true),
    metrictensor_(true),
    fac_(0.0)
    {
        return;
    }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::Evaluate(
    DRT::ELEMENTS::TransportBoundary* ele,
    ParameterList&             params,
    DRT::Discretization&       discretization,
    vector<int>&               lm,
    Epetra_SerialDenseMatrix&  elemat1_epetra,
    Epetra_SerialDenseMatrix&  elemat2_epetra,
    Epetra_SerialDenseVector&  elevec1_epetra,
    Epetra_SerialDenseVector&  elevec2_epetra,
    Epetra_SerialDenseVector&  elevec3_epetra
)
{
  // First, do the things that are needed for all actions:
  // get the material (of the parent element)
  DRT::ELEMENTS::Transport* parentele = ele->ParentElement();
  RefCountPtr<MAT::Material> mat = parentele->Material();

  // Now, check for the action parameter
  const string action = params.get<string>("action","none");
  if (action=="calc_condif_flux")
  {
    // get actual values of transported scalars
    RefCountPtr<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==null) dserror("Cannot get state vector 'phinp'");

    // get velocity values at the nodes (needed for total flux values)
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
    const DRT::ELEMENTS::Transport* parentele = ele->ParentElement();
    const int ielparent = parentele->NumNode();
    // we deal with a (nsd_+1)-dimensional flow field 
    Epetra_SerialDenseVector evel((nsd_+1)*ielparent);
    DRT::UTILS::ExtractMyNodeBasedValues(parentele,evel,velocity,nsd_+1);

    // we dont know the parent element's lm vector; so we have to build it here
    vector<int> lmparent(ielparent);
    vector<int> lmparentowner;
    parentele->LocationVector(discretization, lmparent, lmparentowner);

    // extract local values from the global vector for the parent(!) element 
    vector<double> myphinp(lmparent.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lmparent);

    // get the averaged normal vectors at the nodes
    const RCP<Epetra_MultiVector> normals = params.get< RCP<Epetra_MultiVector> >("normal vectors",null);
    LINALG::SerialDenseVector enormals((nsd_+1)*iel);
    DRT::UTILS::ExtractMyNodeBasedValues(ele,enormals,normals,nsd_+1);

    // set flag for type of scalar
    string scaltypestr=params.get<string>("problem type");
    bool temperature = false;
    if (scaltypestr =="loma") temperature = true;

    double frt(0.0);
    if (scaltypestr =="elch")
    {
      // get parameter F/RT
      frt = params.get<double>("frt");
    }

    // assure, that the values are in the same order as the parent element nodes
    for(int k=0;k<ielparent;++k)
    {
      const Node* node = (parentele->Nodes())[k];
      vector<int> dof = discretization.Dof(node);
      for (unsigned int i=0 ; i< dof.size(); ++i)
      {
        if (dof[i]!=lmparent[k*numdofpernode_ + i])
        {
          cout<<"dof[i]= "<<dof[i]<<"  lmparent[k*numdofpernode + i]="<<lmparent[k*numdofpernode_ + i]<<endl;
          dserror("Dofs are not in the same order as the element nodes. Implement some resorting!");
        }
      }
    }

    // access control parameter
    string fluxtypestring = params.get<string>("fluxtype","noflux");

    // define vector for normal fluxes
    vector<double> mynormflux(lm.size());

    // get node coordinates (we have a nsd_+1 dimensional domain!)
    GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,iel> >(ele,xyze_);

    // do a loop for systems of transported scalars
    for (int j = 0; j<numscal_; ++j)
    {
      // compute fluxes on each node of the parent element
      LINALG::SerialDenseMatrix eflux(3,ielparent,true);
      DRT::Element* peleptr = (DRT::Element*) parentele;
      DRT::ELEMENTS::ScaTraImplInterface::Impl(peleptr)->CalculateFluxSerialDense(
          eflux,
          peleptr,
          myphinp,
          mat,
          temperature,
          frt,
          evel,
          fluxtypestring,
          j);

      // handle the result dofs in the right order (compare lm with lmparent)
      int dofcount = 0;
      for (int i=0; i<iel; ++i)
      {
        for(int k = 0; k<ielparent;++k)
        {
          if (lm[i*numdofpernode_+j]==lmparent[k*numdofpernode_+j]) // dof ids match => assemble this value
          {
            dofcount++;
            // form arithmetic mean of assembled nodal flux vectors
            // => factor is the number of adjacent elements for each node
            const double factor = (double) (ele->Nodes()[i])->NumElement();

            // calculate normal flux at present node
            mynormflux[i] = 0.0;
            for (int m=0; m<(nsd_+1); m++)
            {
              mynormflux[i] += eflux(m,k)*enormals(i*(nsd_+1)+m);
            }

            // store normal flux vector for this node
            elevec1_epetra[i*numdofpernode_+j]+=enormals(i*(nsd_+1))*mynormflux[i]/factor;
            elevec2_epetra[i*numdofpernode_+j]+=enormals(i*(nsd_+1)+1)*mynormflux[i]/factor;
            if (nsd_==2)
            elevec3_epetra[i*numdofpernode_+j]+=enormals(i*(nsd_+1)+2)*mynormflux[i]/factor;
          }
        }
      }
      if (dofcount != iel) dserror("Expected dof for transport boundary element is missing");

      // calculate integral of normal flux
      // => only meaningful for one scalar, for the time being (take first scalar!)
      // NOTE: add integral value only for elements which are NOT ghosted!
      if(j==0 && ele->Owner() == discretization.Comm().MyPID())
      {
        NormalFluxIntegral(ele,params,mynormflux);
      }

    } // loop over numscal
  }
  else if (action == "calc_normal_vectors")
  {
    // access the global vector
    const RCP<Epetra_MultiVector> normals = params.get< RCP<Epetra_MultiVector> >("normal vectors",null);
    if (normals == Teuchos::null) dserror("Could not access vector 'normal vectors'");

    // get node coordinates (we have a nsd_+1 dimensional domain!)
    GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,iel> >(ele,xyze_);

    // determine constant normal to this element
    GetConstNormal(normal_,xyze_);

    // loop over the element nodes
    for (int j=0;j<iel;j++)
    {
      const int nodegid = (ele->Nodes()[j])->Id();
      if (normals->Map().MyGID(nodegid) )
      { // OK, the node belongs to this processor

        // form arithmetic mean of normal vector
        // => numele is the number of adjacent elements for each node
        const double numele = (double) (ele->Nodes()[j])->NumElement();
        for (int dim=0; dim<(nsd_+1); dim++)
        {
          const double component = normal_(dim)/numele;
          normals->SumIntoGlobalValue(nodegid,dim,component);
        }
      }
      //else: the node belongs to another processor; the ghosted
      //      element will contribute the right value on that proc
    }
  }
  else if (action =="calc_elch_electrode_kinetics")
  {
    // get actual values of transported scalars
    RefCountPtr<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==null) dserror("Cannot get state vector 'phinp'");

    // extract local values from the global vector
    vector<double> ephinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,ephinp,lm);

    // get current condition
    Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition> >("condition");
    if (cond == Teuchos::null) dserror("Cannot access condition 'ElectrodeKinetics'");

    // access parameters of the condition
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
    // find out whether we shell use a time curve and get the factor
    // this feature can be also used for stationary "pseudo time loops"
    if (curvenum>=0)
    {
      const double curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
      // adjust potential at metal side accordingly
      pot0 *= curvefac;
    }

    const bool calc_status = params.get<bool>("calc_status",false);
    if (!calc_status)
    {
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


# if 0
      // print all parameters read from the current condition
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
          ele,
          elemat1_epetra,
          elevec1_epetra,
          ephinp,
          mat,
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
    else
    {
      // NOTE: add integral value only for elements which are NOT ghosted!
      if(ele->Owner() == discretization.Comm().MyPID())
      { ElectrodeStatus(
            ele,
            params,
            ephinp,
            kinetics,
            pot0,
            alphaa,
            alphac,
            i0,
            frt,
            iselch);} 
    }
  } 
  else if (action =="calc_therm_press")
  {
    // get actual values of transported scalars
    RefCountPtr<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==null) dserror("Cannot get state vector 'phinp'");

    // get velocity values at the nodes (needed for total flux values)
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
    DRT::ELEMENTS::Transport* parentele = ele->ParentElement();
    const int ielparent = parentele->NumNode();
    // we deal with a (nsd_+1)-dimensional flow field 
    Epetra_SerialDenseVector evel((nsd_+1)*ielparent);
    DRT::UTILS::ExtractMyNodeBasedValues(parentele,evel,velocity,nsd_+1);

    // get actual values of density
    RefCountPtr<const Epetra_Vector> densnp = discretization.GetState("densnp");
    if (densnp==null) dserror("Cannot get state vector 'densnp'");

    // we dont know the parent element's lm vector; so we have to build it here
    vector<int> lmparent(ielparent);
    vector<int> lmparentowner;
    parentele->LocationVector(discretization, lmparent, lmparentowner);

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

    // get node coordinates (we have a nsd_+1 dimensional domain!)
    GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,iel> >(ele,xyze_);

    // determine normal to this element
    GetConstNormal(normal_,xyze_);

    // compute fluxes on each node of the parent element
    LINALG::SerialDenseMatrix eflux(3,ielparent);
    DRT::Element* peleptr = (DRT::Element*) parentele;
    // set some parameters
    double frt=0.0;
    string fluxtypestring("diffusiveflux");
    int j=0;
    DRT::ELEMENTS::ScaTraImplInterface::Impl(parentele)->CalculateFluxSerialDense(
        eflux,
        peleptr,
        myphinp,
        mat,
        temperature,
        frt,
        evel,
        fluxtypestring,
        j);

    for (int i=0; i<iel; ++i)
    {
      for(int k = 0; k<ielparent;++k)
      {
        // calculate normal diffusive flux and velocity div. at present node
        mydiffflux[i] = 0.0;
        mydivu[i]     = 0.0;
        for (int l=0; l<nsd_+1; l++)
        {
          mydiffflux[i] += eflux(l,k)*normal_(l);
          mydivu[i]     += (evel[i*(nsd_+1)+l]/mydensnp[i])*normal_(l);
        }
      }
    }

    // calculate integral of normal diffusive flux and velocity divergence
    // NOTE: add integral value only for elements which are NOT ghosted!
    if(ele->Owner() == discretization.Comm().MyPID())
    {
      DifffluxAndDivuIntegral(ele,params,mydiffflux,mydivu);
    }
  }
  else
    dserror("Unknown type of action for Scatra Implementation: %s",action.c_str());

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Surface/Line Neumann boundary condition       gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::EvaluateNeumann(
    DRT::Element*             ele,
    ParameterList&            params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    vector<int>&              lm,
    Epetra_SerialDenseVector& elevec1)
    {
  // get node coordinates (we have a nsd_+1 dimensional domain!)
  GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,iel> >(ele,xyze_);

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

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
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    EvalShapeFuncAndIntFac(intpoints,iquad,ele->Id());

    // multiply integration factor with the timecurve factor
    fac_ *= curvefac;

    // factor given by spatial function
    double functfac = 1.0;
    // determine global coordinates of current Gauss point
    double coordgp[nsd_];
    for (int i = 0; i< nsd_; i++)
    {
      coordgp[i] = 0.0;
      for (int j = 0; j < iel; j++)
      {
        coordgp[i] += xyze_(i,j) * funct_(j);
      }
    }

    int functnum = -1;
    const double* coordgpref = &coordgp[0]; // needed for function evaluation


    for(int dof=0;dof<numdofpernode_;dof++)
    {
      if ((*onoff)[dof]) // is this dof activated?
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

        const double val_fac_functfac = (*val)[dof]*fac_*functfac;

        for (int node=0;node<iel;++node)
        {
          elevec1[node*numdofpernode_+dof] += funct_(node)*val_fac_functfac;
        }
      } // if ((*onoff)[dof])
    }
  } //end of loop over integration points

  return 0;
    }


/*----------------------------------------------------------------------*
 | evaluate shape functions and int. factor at int. point     gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::EvalShapeFuncAndIntFac(
    const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
    const int&                                   iquad,      ///< id of current Gauss point
    const int&                                   eleid       ///< the element id
)
{
  // coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (int idim=0;idim<nsd_;idim++)
  {xsi_(idim) = gpcoord[idim];}

  // shape functions and their first derivatives
  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);

  // the metric tensor and the area of an infinitesimal surface/line element
  double drs(0.0);
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>(xyze_,deriv_,metrictensor_,drs);

  // set the integration factor
  fac_ = intpoints.IP().qwgt[iquad] * drs;

  // say goodbye
  return;
}


/*----------------------------------------------------------------------*
 | evaluate an electrode kinetics boundary condition (private) gjb 01/09|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::EvaluateElectrodeKinetics(
    const DRT::Element*        ele,
    Epetra_SerialDenseMatrix& emat,
    Epetra_SerialDenseVector& erhs,
    const vector<double>&   ephinp,
    Teuchos::RCP<const MAT::Material> material,
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

  //pre-multiplication with 1/(F*z_1)
  double fz = 1.0/96485.3399;

  // get node coordinates (we have a nsd_+1 dimensional domain!)
  GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,iel> >(ele,xyze_);

  if (iselch)
  {
    // get valence of the single(!) reactant
    if (material->MaterialType() == INPAR::MAT::m_matlist)
    {
      const MAT::MatList* actmat = static_cast<const MAT::MatList*>(material.get());

      if (actmat->MatID(0) != rctid) 
        dserror("active species is not first scalar in material list!");
      // the active species is the FIRST material in the material list. ALWAYS!
      Teuchos::RCP<const MAT::Material> singlemat = actmat->MaterialById(rctid);
      if (singlemat->MaterialType() == INPAR::MAT::m_ion)
      {
        const MAT::Ion* actsinglemat = static_cast<const MAT::Ion*>(singlemat.get());
        fz = fz/actsinglemat->Valence();
      }
      else
        dserror("single material type is not 'ion'");
    }
    else
      dserror("material type is not a 'matlist' material");
  }

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // concentration values of reactive species at element nodes
  LINALG::Matrix<iel,1> conreact(true);

  // el. potential values at element nodes
  LINALG::Matrix<iel,1> pot(true);
  if(iselch)
  {
    for (int inode=0; inode< iel;++inode)
    {
      conreact(inode) = ephinp[inode*numdofpernode_];
      pot(inode) = ephinp[inode*numdofpernode_+numscal_];
    }
  }
  else
  {
    for (int inode=0; inode< iel;++inode)
    {
      conreact(inode) = 1.0;
      pot(inode) = ephinp[inode*numdofpernode_];
    }
  }

  // concentration of active species at integration point
  double conint;
  // el. potential at integration point
  double potint;
  // a 'working variable'
  double fac_fz_i0_funct_vi;

 /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    EvalShapeFuncAndIntFac(intpoints,gpid,ele->Id());

    // elch-specific values at integration point:
    conint = funct_.Dot(conreact);
    potint = funct_.Dot(pot);

    // surface overpotential eta at integration point
    const double eta = (pot0 - potint);

    double gammak = 1.0;
    double pow_conint_gamma_k = pow(conint,gammak);

    if (iselch)
    {
      const double expterm = exp(alphaa*frt*eta)-exp((-alphac)*frt*eta);

      for (int vi=0; vi<iel; ++vi)
      {
        fac_fz_i0_funct_vi = fac_*fz*i0*funct_(vi);
        // ---------------------matrix
        for (int ui=0; ui<iel; ++ui)
        {
          emat(vi*numdofpernode_,ui*numdofpernode_) += fac_fz_i0_funct_vi*gammak*pow(conint,(gammak-1.0))*funct_(ui)*expterm; 
          emat(vi*numdofpernode_,ui*numdofpernode_+numscal_) += fac_fz_i0_funct_vi*pow_conint_gamma_k*(((-alphaa)*frt*exp(alphaa*frt*eta))+((-alphac)*frt*exp((-alphac)*frt*eta)))*funct_(ui);
        }
        // ------------right-hand-side
        erhs[vi*numdofpernode_] -= fac_fz_i0_funct_vi*pow_conint_gamma_k*expterm;
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
        const double fac_i0_funct_vi = fac_*i0*funct_(vi);
        // ---------------------matrix
        for (int ui=0; ui<iel; ++ui)
        {
          emat(vi*numdofpernode_,ui*numdofpernode_) += fac_i0_funct_vi*exptermderiv*funct_(ui);
        }
        // ------------right-hand-side
        erhs[vi*numdofpernode_] -= fac_i0_funct_vi*expterm;
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

} // ScaTraBoundaryImpl<distype>::EvaluateElectrodeKinetics()


/*----------------------------------------------------------------------*
 | calculate electrode kinetics status information             gjb 01/09|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::ElectrodeStatus(
    const DRT::Element*        ele,
    ParameterList&          params,
    const vector<double>&   ephinp,
    const std::string*    kinetics,
    const double&             pot0,
    const double&           alphaa,
    const double&           alphac,
    const double&               i0,
    const double&              frt,
    const bool&             iselch
)
{
  // get node coordinates (we have a nsd_+1 dimensional domain!)
  GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,iel> >(ele,xyze_);

  // get variables with their current values
  double currentintegral  = params.get<double>("currentintegral");
  double boundaryint      = params.get<double>("boundaryintegral");
  double overpotentialint = params.get<double>("overpotentialintegral");
  double concentrationint = params.get<double>("concentrationintegral");

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // concentration values of reactive species at element nodes
  LINALG::Matrix<iel,1> conreact(true);

  // el. potential values at element nodes
  LINALG::Matrix<iel,1> pot(true);
  if(iselch)
  {
    for (int inode=0; inode< iel;++inode)
    {
      conreact(inode) = ephinp[inode*numdofpernode_];
      pot(inode) = ephinp[inode*numdofpernode_+numscal_];
    }
  }
  else
  {
    for (int inode=0; inode< iel;++inode)
    {
      conreact(inode) = 1.0;
      pot(inode) = ephinp[inode*numdofpernode_];
    }
  }

  // concentration of active species at integration point
  double conint;
  // el. potential at integration point
  double potint;

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    EvalShapeFuncAndIntFac(intpoints,gpid,ele->Id());

    // elch-specific values at integration point:
    conint = funct_.Dot(conreact);
    potint = funct_.Dot(pot);

    // surface overpotential eta at integration point
    const double eta = (pot0 - potint);
    //Butler-Volmer
    const double expterm = conint * (exp(alphaa*frt*eta)-exp((-alphac)*frt*eta));

    // compute integrals
    overpotentialint += eta * fac_;
    currentintegral += (-i0) * expterm * fac_; // the negative(!) normal flux density
    boundaryint += fac_;
    concentrationint += conint*fac_;

  } // loop over integration points

  // add contributions to the global values
  params.set<double>("currentintegral",currentintegral);
  params.set<double>("boundaryintegral",boundaryint);
  params.set<double>("overpotentialintegral",overpotentialint);
  params.set<double>("concentrationintegral",concentrationint);

} //ScaTraBoundaryImpl<distype>::ElectrodeStatus


/*----------------------------------------------------------------------*
 | get constant normal                                        gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::GetConstNormal(
    LINALG::Matrix<nsd_+1,1>&          normal,
    const LINALG::Matrix<nsd_+1,iel>&  xyze
    )
{
  // determine normal to this element
  switch(nsd_)
  {
  case 2:
  {
    LINALG::Matrix<3,1> dist1(true), dist2(true);
    for (int i=0; i<3; i++)
    {
      dist1(i) = xyze(i,1)-xyze(i,0);
      dist2(i) = xyze(i,2)-xyze(i,0);
    }

    normal(0) = dist1(1)*dist2(2) - dist1(2)*dist2(1);
    normal(1) = dist1(2)*dist2(0) - dist1(0)*dist2(2);
    normal(2) = dist1(0)*dist2(1) - dist1(1)*dist2(0);
  }
  break;
  case 1:
  {  
    normal(0) = xyze(1,1) - xyze(1,0);
    normal(1) = (-1.0)*(xyze(0,1) - xyze(0,0));
  }
  break;
  default:
    dserror("Illegal number of space dimensions: %d",nsd_);
  } // switch(nsd)

  // length of normal to this element
  const double length = normal.Norm2();
  // outward-pointing normal of length 1.0
  normal.Scale(1/length);

  return;
} // ScaTraBoundaryImpl<distype>::


/*----------------------------------------------------------------------*
 | calculate integral of normal flux (private)                  vg 09/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::NormalFluxIntegral(
    const DRT::Element*             ele,
    ParameterList&                  params,
    const vector<double>&           enormflux
)
{
  // get variables with their current values
  double normfluxintegral = params.get<double>("normfluxintegral");
  double boundaryint = params.get<double>("boundaryint");

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    EvalShapeFuncAndIntFac(intpoints,gpid,ele->Id());

    // compute integral of normal flux
    for (int node=0;node<iel;++node)
    {
      normfluxintegral += funct_(node) * enormflux[node] * fac_;
      boundaryint += funct_(node) * fac_;
    }
  } // loop over integration points

  // add contributions to the global values
  params.set<double>("normfluxintegral",normfluxintegral);
  params.set<double>("boundaryint",boundaryint);

} //ScaTraBoundaryImpl<distype>::NormalFluxIntegral


/*----------------------------------------------------------------------*
 | calculate integral of normal diffusive flux + velocity div.  vg 09/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::DifffluxAndDivuIntegral(
    const DRT::Element*             ele,
    ParameterList&                  params,
    const vector<double>&           ediffflux,
    const vector<double>&           edivu
)
{
  // get variables for integrals of normal diffusive flux and velocity div.
  double difffluxintegral = params.get<double>("diffusive-flux integral");
  double divuintegral     = params.get<double>("velocity-divergence integral");

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    EvalShapeFuncAndIntFac(intpoints,gpid,ele->Id());

    // compute integral of normal flux
    for (int node=0;node<iel;++node)
    {
      difffluxintegral += funct_(node) * ediffflux[node] * fac_;
      divuintegral     += funct_(node) * edivu[node] * fac_;
    }
  } // loop over integration points

  // add contributions to the global values
  params.set<double>("diffusive-flux integral",difffluxintegral);
  params.set<double>("velocity-divergence integral",divuintegral);

  return;

} //ScaTraBoundaryImpl<distype>::DifffluxAndDivuIntegral


#endif // CCADISCRET
#endif // D_FLUID3 or D_FLUID2
