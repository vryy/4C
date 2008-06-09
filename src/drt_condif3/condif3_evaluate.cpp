/*!
\file condif3_evaluate.cpp
\brief

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*/
#ifdef D_FLUID3
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

#include "condif3.H"
#include "../drt_mat/convecdiffus.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"


using namespace DRT::UTILS;

/*---------------------------------------------------------------------*
|  converts a string into an action for this element(private) gjb 06/08|
 *---------------------------------------------------------------------*/
DRT::ELEMENTS::Condif3::ActionType DRT::ELEMENTS::Condif3::convertStringToActionType(
  const string& action) const
{
  dsassert(action != "none", "No action supplied");

  DRT::ELEMENTS::Condif3::ActionType act = Condif3::none;
  if (action == "calc_condif_systemmat_and_residual")
    act = Condif3::calc_condif_systemmat_and_residual;
  else
    dserror("Unknown type of action for Condif3");
  return act;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              gjb 06/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Condif3::Evaluate(ParameterList& params,
    DRT::Discretization&      discretization,
    vector<int>&              lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  // get the action required
  const string action = params.get<string>("action","none");
  const DRT::ELEMENTS::Condif3::ActionType act = convertStringToActionType(action);

  // get the material
  RefCountPtr<MAT::Material> mat = Material();

  if (mat->MaterialType() != m_condif)
    dserror("Material law is not a condif element");

  MATERIAL* actmat = NULL;

  if(mat->MaterialType()== m_condif)
    actmat = static_cast<MAT::ConvecDiffus*>(mat.get())->MaterialData();
    else
      dserror("condif material expected but got type %d", mat->MaterialType());

  switch(act)
  {
  // the standard one-step-theta implementation
  case calc_condif_systemmat_and_residual:
  {
    // need current velocity and history vector
    RefCountPtr<const Epetra_Vector> hist = discretization.GetState("hist");
    if (hist==null) dserror("Cannot get state vector 'hist'");

    // extract local values from the global vector
    vector<double> myhist(lm.size());
    DRT::UTILS::ExtractMyValues(*hist,myhist,lm);

    // get control parameter
    const bool is_stationary = params.get<bool>("using stationary formulation",false);
    const double time = params.get<double>("total time",-1.0);

    // One-step-Theta: timefac = theta*dt
    // BDF2:           timefac = 2/3 * dt
    double timefac = 0.0;
    if (not is_stationary)
    {
      timefac = params.get<double>("thsl",-1.0);
      if (timefac < 0.0) dserror("No thsl supplied");
    }

    // get velocity values at the nodes
    // compare also with DRT::UTILS::ExtractMyValues()
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
    const int iel = NumNode();
    const int nsd=3;
    Epetra_SerialDenseVector evel(nsd*iel);
    ExtractMyNodeBasedValues(evel,velocity);

    // get flag for fine-scale subgrid diffusivity
    string fssgd = params.get<string>("fs subgrid diffusivity","No");
    if (fssgd != "No") dserror("fssgd not yet implemented!");

    // calculate element coefficient matrix and rhs
    Condif3SysMat(lm,myhist,&elemat1,&elemat2,&elevec1,elevec2,actmat,time,timefac,evel,fssgd,is_stationary);

  }
  break;
  default:
    dserror("Unknown type of action for Condif3");
  } // end of switch(act)

  return 0;
} //DRT::ELEMENTS::Condif3::Evaluate


/*----------------------------------------------------------------------*
 |  calculate system matrix and rhs for convec.-diff.(private) gjb 06/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif3::Condif3SysMat(
    vector<int>&              lm,
    vector<double>&           ehist,
    Epetra_SerialDenseMatrix* sys_mat,
    Epetra_SerialDenseMatrix* sys_mat_sd,
    Epetra_SerialDenseVector* residual,
    Epetra_SerialDenseVector& sugrvisc,
    struct _MATERIAL*         material,
    double                    time,
    double                    timefac,
    Epetra_SerialDenseVector& evel,
    string                    fssgd,
    bool                      is_stationary)
{

  /*------------------------------------------------- set element data */
  const int iel = NumNode();
  const DiscretizationType distype = this->Shape();
  const int nsd = 3;

  Epetra_SerialDenseMatrix xyze(nsd,iel);

  // get node coordinates
  for (int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
    xyze(2,i)=Nodes()[i]->X()[2];
  }

  // dead load in element nodes
  const Epetra_SerialDenseVector bodyforce = BodyForce(time);

  /*---------------------------------------------- get diffusivity ---*/
  if(material->mattyp != m_condif) dserror("Material law is not of type m_condif.");
  const double diffus = material->m.condif->diffusivity;

  /*----------------------------------------- declaration of variables ---*/
  Epetra_SerialDenseVector      funct(iel);
  Epetra_SerialDenseMatrix      deriv(nsd,iel);
  Epetra_SerialDenseMatrix      deriv2(6,iel);
  static Epetra_SerialDenseMatrix xjm(nsd,nsd);
  Epetra_SerialDenseMatrix      derxy(nsd,iel);
  Epetra_SerialDenseMatrix      derxy2(6,iel);
  static vector<double>         velint(nsd);
  double                        edeadng;
  double                        hist; /* history data at integration point      */
  double                        hk;  /* element length for calculation of tau  */
  double                        vel_norm, epe1, epe2, xi1, xi2;
  double                        mk=0.0;
  double                        tau; // stabilization parameter
  double                        kart; // artificial diffusivity

  /*----------------------------------------------------------------------*/
  // calculation of stabilization parameter
  /*----------------------------------------------------------------------*/

  //ToDo!
  
  /*----------------------------------------------------------------------*/
  // integration loop for one condif3 element
  /*----------------------------------------------------------------------*/

  // flag for higher order elements
  const bool higher_order_ele = isHigherOrderElement(distype);

  // gaussian points
  const GaussRule3D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints3D  intpoints = getIntegrationPoints3D(gaussrule);

  // integration loop
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    // coordiantes of the current integration point
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];
    const double e3 = intpoints.qxg[iquad][2];

    // shape functions and their derivatives
    DRT::UTILS::shape_function_3D(funct,e1,e2,e3,distype);
    DRT::UTILS::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

    if (higher_order_ele)
    {
       shape_function_3D_deriv2(deriv2,e1,e2,e3,distype);
    }

    /*----------------------------------------- compute Jacobian matrix */

    // get Jacobian matrix and determinant
    // actually compute its transpose....
    /*
      +-            -+ T      +-            -+
      | dx   dx   dx |        | dx   dy   dz |
      | --   --   -- |        | --   --   -- |
      | dr   ds   dt |        | dr   dr   dr |
      |              |        |              |
      | dy   dy   dy |        | dx   dy   dz |
      | --   --   -- |   =    | --   --   -- |
      | dr   ds   dt |        | ds   ds   ds |
      |              |        |              |
      | dz   dz   dz |        | dx   dy   dz |
      | --   --   -- |        | --   --   -- |
      | dr   ds   dt |        | dt   dt   dt |
      +-            -+        +-            -+
    */
    double dum;
     /*-------------------------------- determine jacobian at point r,s ---*/
     for (int i=0; i<nsd; i++)
     {
        for (int j=0; j<nsd; j++)
        {
           dum=0.0;
           for (int l=0; l<iel; l++)
           {
              dum += deriv(i,l)*xyze(j,l);
           }
           xjm(i,j)=dum;
        } /* end of loop j */
     } /* end of loop i */
     // ---------------------------------------- calculate determinant
     const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                        xjm(0,1)*xjm(1,2)*xjm(2,0)+
                        xjm(0,2)*xjm(1,0)*xjm(2,1)-
                        xjm(0,2)*xjm(1,1)*xjm(2,0)-
                        xjm(0,0)*xjm(1,2)*xjm(2,1)-
                        xjm(0,1)*xjm(1,0)*xjm(2,2);

     if (det < 0.0)
       dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", Id(), det);
     if (abs(det) < 1E-16)
       dserror("GLOBAL ELEMENT NO.%i\nZERO JACOBIAN DETERMINANT: %f", Id(), det);

     const double fac = intpoints.qwgt[iquad]*det; // Gauss weight * det(J)

    /*------------------------------------------------------------------*/
    /*                                         compute global derivates */
    /*------------------------------------------------------------------*/

     /*------------------------------------------------- initialization */
     for(int k=0;k<iel;k++)
     {
       derxy(0,k)=0.0;
       derxy(1,k)=0.0;
       derxy(2,k)=0.0;
     } /* end of loop over k */
     
     // ----------------------------------------------- inverse of jacobian
     static Epetra_SerialDenseMatrix       xij(nsd,nsd);
     double idet = 1./det;
     xij(0,0) = (  xjm(1,1)*xjm(2,2) - xjm(2,1)*xjm(1,2))*idet;
     xij(1,0) = (- xjm(1,0)*xjm(2,2) + xjm(2,0)*xjm(1,2))*idet;
     xij(2,0) = (  xjm(1,0)*xjm(2,1) - xjm(2,0)*xjm(1,1))*idet;
     xij(0,1) = (- xjm(0,1)*xjm(2,2) + xjm(2,1)*xjm(0,2))*idet;
     xij(1,1) = (  xjm(0,0)*xjm(2,2) - xjm(2,0)*xjm(0,2))*idet;
     xij(2,1) = (- xjm(0,0)*xjm(2,1) + xjm(2,0)*xjm(0,1))*idet;
     xij(0,2) = (  xjm(0,1)*xjm(1,2) - xjm(1,1)*xjm(0,2))*idet;
     xij(1,2) = (- xjm(0,0)*xjm(1,2) + xjm(1,0)*xjm(0,2))*idet;
     xij(2,2) = (  xjm(0,0)*xjm(1,1) - xjm(1,0)*xjm(0,1))*idet;
     /*---------------------------------------- calculate global derivatives */
     for (int k=0;k<iel;k++)
     {
       derxy(0,k) +=   xij(0,0) * deriv(0,k) + xij(0,1) * deriv(1,k) + xij(0,2) * deriv(2,k);
       derxy(1,k) +=   xij(1,0) * deriv(0,k) + xij(1,1) * deriv(1,k) + xij(1,2) * deriv(2,k);
       derxy(2,k) +=   xij(2,0) * deriv(0,k) + xij(2,1) * deriv(1,k) + xij(2,2) * deriv(2,k);
     } /* end of loop over k */


    /*--------------------------------- compute second global derivative */
    if (higher_order_ele) {cout<<"second derivative not yet supported"<<endl;}
    //condif2_gder2(xyze,xjm,derxy,derxy2,deriv2,iel);

    /*---------------------- get velocity at integration point */
    // use same shape functions for velocity as for unknown scalar field phi
    for (int i=0;i<nsd;i++)
    {
        velint[i]=0.0;
        for (int j=0;j<iel;j++)
        {
            velint[i] += funct[j]*evel[i+(nsd*j)];
        }
    } //end loop over i

    /*---------------- get history data (n,i) at integration point */
    hist=0.0;
    for (int j=0;j<iel;j++)
    {
      hist += funct[j]*ehist[j];
    }

    // get bodyforce in gausspoint
    edeadng = 0.0;
    for (int inode=0;inode<iel;inode++)
    {
      edeadng+= bodyforce[inode]*funct[inode];
    }

    /*-------------- perform integration for entire matrix and rhs ---*/
    if(is_stationary==false)
    {
      Condif3CalMat(*sys_mat,*sys_mat_sd,*residual,velint,hist,funct,derxy,derxy2,edeadng,tau,kart,fac,diffus,iel,fssgd,timefac);
      cout<<"calmat_instationary"<<endl;
    }
      else
      cout<<"calmat_stationary"<<endl;
      //condif2_calmat_stat(*sys_mat,*sys_mat_sd,*residual,velint,hist,funct,derxy,
      //                    derxy2,edeadng,tau,kart,fac,diffus,iel,fssgd);

#if 1
    // debug output
    cout<<"jacobi determinant = "<<det<<endl;
    cout<<"global derivatives: "<<derxy<<endl;
#endif

  } // integration loop
  

#if 1
  //debug output
  cout<<"diffusivity:"<<diffus<<endl;
  cout<<"node coordinates: "<<endl<<xyze<<endl;
  cout<<"sys_mat:"<<*sys_mat<<endl;
  cout<<"residual:"<<*residual<<endl;
  cout<<"velocity:"<<evel<<endl;
#endif

  return;
} // DRT::ELEMENTS::Condif3::Condif3SysMat


/*----------------------------------------------------------------------*
 |  do nothing (public)                                        gjb 06/08|
 |                                                                      |
 |  The function is just a dummy. For the condif elements, the          |
 |  integration of the volume neumann (body forces) loads takes place   |
 |  in the element. We need it there for the stabilisation terms!       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Condif3::EvaluateNeumann(ParameterList& params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    vector<int>&              lm,
    Epetra_SerialDenseVector& elevec1)
{
  return 0;
}


/*----------------------------------------------------------------------*
 |  get optimal gaussrule for discretization type  (private)   gjb 06/08|
 *----------------------------------------------------------------------*/
GaussRule3D DRT::ELEMENTS::Condif3::getOptimalGaussrule(const DiscretizationType& distype)
{
  GaussRule3D rule = intrule3D_undefined;
  switch (distype)
  {
  case hex8:
    rule = intrule_hex_8point;
    break;
  case hex20: case hex27:
    rule = intrule_hex_27point;
    break;
  case tet4:
    rule = intrule_tet_4point;
    break;
  case tet10:
    rule = intrule_tet_5point;
    break;
  default:
    dserror("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}


/*-------------------------------------------------------------------------*
 | check for higher order derivatives for shape functions(private)gjb 06/08|
 *-------------------------------------------------------------------------*/
bool DRT::ELEMENTS::Condif3::isHigherOrderElement(
  const DRT::Element::DiscretizationType& distype) const
{
  bool hoel = true;
  switch (distype)
  {
  case hex8: case hex20: case hex27: case tet10:
    hoel = true;
    break;
  case tet4: 
    hoel = false;
    break;
  case wedge6: case pyramid5: case wedge15: 
    //!!!TODO:  wedge und pyramid have 2nd derivatives!!!!!!!!!!!!!!!!!!!!!!!!
    dserror("wedge and pyramids have second derivatives!");
    break;
  default:
    dserror("distype unknown!");
  }
  return hoel;
}


/*----------------------------------------------------------------------*
 |  get the body force  (private)                              gjb 06/08|
 *----------------------------------------------------------------------*/
const Epetra_SerialDenseVector DRT::ELEMENTS::Condif3::BodyForce(const double time) const
{
   const int iel = NumNode();
   Epetra_SerialDenseVector bodyforce(iel); // initialization included
   vector<DRT::Condition*> myneumcond;

   // check whether all nodes have a unique VolumeNeumann condition
   DRT::UTILS::FindElementConditions(this, "VolumeNeumann", myneumcond);

   if (myneumcond.size()>1)
   {
     dserror("more than one VolumeNeumann cond on one node");
   }

   if (myneumcond.size()==1)
   {
    dserror("body force not implemented for condif3 element.");
    }

  return bodyforce;

} //Condif3::BodyForce


/*----------------------------------------------------------------------*
 |  extract my velocity dof from global vector (private)       gjb 06/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif3::ExtractMyNodeBasedValues(Epetra_SerialDenseVector& local,
     const RCP<Epetra_MultiVector>& global)
{
  const int nsd = global->NumVectors(); // get dimension
  const int iel = NumNode(); // number of nodes
  if (local.Length()!=(iel*nsd)) dserror("vector size mismatch.");

  for (int i=0; i<nsd; i++)
  {
    // get actual component column of velocity multi-vector
    double* globalcolumn = (*global)[i];
    // loop over the nodes
    for (int j=0;j<iel;j++)
    {
      const int nodegid = (Nodes()[j])->Id();
      const int lid = global->Map().LID(nodegid);
      local(i+(nsd*j))=globalcolumn[lid];
    }
  }
  return;
} //Condif3::ExtractMyNodeBasedValues


/*----------------------------------------------------------------------*
 |  evaluate instationary convection-diffusion matrix (private)gjb 06/08|
 *----------------------------------------------------------------------*/

/*
In this routine the Gauss point contributions to the elemental coefficient
matrix of a stabilized condif3 element are calculated for the instationary
case. The procedure is based on the Rothe method of first discretizing in
time. Hence the resulting terms include coefficients containing time
integration variables such as theta or delta t which are represented by
'timefac'.

The stabilization is based on the residuum:

R = phi + timefac * u * grad(phi) - timefac * diffus * laplace(phi) - rhsint

The corresponding weighting operators are
L = timefac * u * grad(w) +/- timefac * diffus * laplace(w)

'+': USFEM
'-': GLS


time-integration schemes:

one-step-theta:
rhsint = u_old + Theta dt f + (1-Theta) acc_old

BDF2:

generalised alpha:


The calculation proceeds as follows.
1) obtain single operators of R and L
2) build Galerkin terms from them
3) build stabilizing terms from them
4) build Galerkin and stabilizing terms of RHS

NOTE: Galerkin and stabilization matrices are calculated within one
      routine.


for further comments see comment lines within code.

</pre>
\param **estif      DOUBLE        (o)   ele stiffness matrix
\param  *eforce     DOUBLE        (o)   ele force vector
\param  *velint     DOUBLE        (i)   vel at INT point
\param  *hist       DOUBLE        (i)   rhs at INT point
\param  *funct      DOUBLE        (i)   nat. shape funcs
\param **derxy      DOUBLE        (i)   global coord. deriv
\param **derxy2     DOUBLE        (i)   2nd global coord. deriv.
\param   edeadng    DOUBLE        (i)   dead load
\param   tau        DOUBLE        (i)   stabilization parameter
\param   fac        DOUBLE        (i)   weighting factor
\param   diffus     DOUBLE        (i)   diffusivity
\param   iel        INT           (i)   number of nodes of act. ele
\return void
------------------------------------------------------------------------*/

void DRT::ELEMENTS::Condif3::Condif3CalMat(
    Epetra_SerialDenseMatrix& estif,
    Epetra_SerialDenseMatrix& esd,
    Epetra_SerialDenseVector& eforce,
    vector<double>&           velint,
    const double&             hist,
    Epetra_SerialDenseVector& funct,
    Epetra_SerialDenseMatrix& derxy,
    Epetra_SerialDenseMatrix& derxy2,
    const double&             edeadng,
    const double&             tau,
    const double&             kart,
    const double&             fac,
    const double&             diffus,
    const int&                iel,
    string                    fssgd,
    double                    timefac
    )
{
/*========================= further variables =========================*/

vector<double>            conv(iel);        /* convective part       */
vector<double>            diff(iel);        /* diffusive part        */
static double             rhsint;           /* rhs at int. point     */

// stabilization parameter
const double taufac = tau*fac;

// integration factors and coefficients of single terms
const double timefacfac  = timefac * fac;
const double timetaufac  = timefac * taufac;

/*-------------------------------- evaluate rhs at integration point ---*/
rhsint = hist + edeadng*timefac;

for (int i=0; i<iel; i++) /* loop over nodes of element */
{
   /* convective part */
   /* u_x * N,x  +  u_y * N,y          with  N .. form function matrix */
   conv[i] = velint[0] * derxy(0,i) + velint[1] * derxy(1,i) ;

   /* diffusive part */
   /* diffus * ( N,xx  +  N,yy  ) */
   diff[i] = diffus * ( derxy2(0,i) + derxy2(1,i) );

} // end of loop over nodes of element

/*--------------------------------- now build single stiffness terms ---*/

// -------------------------------------------System matrix
for (int vi=0; vi<iel; ++vi)
{
    for (int ui=0; ui<iel; ++ui)
    {
    /* Standard Galerkin terms: */
    /* transient term */
    //estif(vi, ui) += fac*funct_(vi)*funct_(ui) ;

    /* convective term */
    estif(vi, ui) += timefacfac*funct[vi]*conv[ui] ;

    /* diffusive term */
    estif(vi, ui) += timefacfac*diffus*(derxy(0, ui)*derxy(0, vi) + derxy(1, ui)*derxy(1, vi)) ;

    /* Stabilization terms: */
    /* 1) transient stabilization (USFEM assumed here, sign change necessary for GLS) */
    /* transient term */
    /*estif(vi, ui) += -taufac*funct_(vi)*funct_(ui) ;*/

    /* convective term */
    /*estif(vi, ui) += -timetaufac*funct_(vi)*conv_(ui) ;*/

    /* diffusive term */
    /*estif(vi, ui) += timetaufac*funct_(vi)*diff_(ui) ;*/

    /* 2) convective stabilization */
    /* transient term */
    estif(vi, ui) += taufac*conv[vi]*funct[ui];

    /* convective term */
    estif(vi, ui) += timetaufac*conv[vi]*conv[ui] ;

    /* diffusive term */
    estif(vi, ui) += -timetaufac*conv[vi]*diff[ui] ;

    /* 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */
    /* transient term */
    estif(vi, ui) += taufac*diff[vi]*funct[ui] ;

    /* convective term */
    estif(vi, ui) += timetaufac*diff[vi]*conv[ui] ; 

    /* diffusive term */
    estif(vi, ui) += -timetaufac*diff[vi]*diff[ui] ;
  }
}

// ----------------------------------------------RHS
for (int vi=0; vi<iel; ++vi)
{
  /* RHS source term */
  eforce[vi] += fac*funct[vi]*rhsint ;

  /* transient stabilization of RHS source term */
  /*eforce[vi] += -taufac*funct_(vi)*rhsint ;*/

  /* convective stabilization of RHS source term */
  eforce[vi] += taufac*conv[vi]*rhsint ;

  /* diffusive stabilization of RHS source term */
  eforce[vi] += taufac*diff[vi]*rhsint ;
}

return;
} //DRT:ELEMENTS:Condif3:Condif3CalMat


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*----------------------------------------------------------------------*
 |  init the element (public)                                mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Condif3Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
