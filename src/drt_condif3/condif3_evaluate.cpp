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
#include "condif3_utils.H"
#include "condif3_impl.H"
#include "../drt_mat/convecdiffus.H"
#include "../drt_mat/matlist.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H" // for CalculateFlux()
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include <Epetra_SerialDenseSolver.h>

#include "../drt_lib/drt_globalproblem.H"


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
  else if (action == "calc_initial_time_deriv")
    act = Condif3::calc_initial_time_deriv;
  else if (action == "calc_subgrid_diffusivity_matrix")
    act = Condif3::calc_subgrid_diffusivity_matrix;
  else if (action == "calc_condif_flux")
    act = Condif3::calc_condif_flux;
  else if (action == "calc_temp_and_dens")
    act = Condif3::calc_temp_and_dens;
  else if (action == "calc_elch_kwok_error")
    act = Condif3::calc_elch_kwok_error;
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

  MATERIAL* actmat = NULL;

  if(mat->MaterialType()== m_condif)
    actmat = static_cast<MAT::ConvecDiffus*>(mat.get())->MaterialData();
  else if (mat->MaterialType()== m_matlist)
    actmat = static_cast<MAT::MatList*>(mat.get())->MaterialData();
  else
    dserror("condif material expected but got type %d", mat->MaterialType());

  switch(act)
  {
  case DRT::ELEMENTS::Condif3::calc_condif_systemmat_and_residual:
  {
  return DRT::ELEMENTS::Condif3ImplInterface::Impl(this)->Evaluate(
           this,
           params,
           discretization,
           lm,
           elemat1,
           elemat2,
           elevec1,
           elevec2,
           elevec3,
           mat,
           actmat);
  }
  break;
  // calculate time derivative for time value t_0
  case DRT::ELEMENTS::Condif3::calc_initial_time_deriv:
  {
    return DRT::ELEMENTS::Condif3ImplInterface::Impl(this)->Evaluate(
             this,
             params,
             discretization,
             lm,
             elemat1,
             elemat2,
             elevec1,
             elevec2,
             elevec3,
             mat,
             actmat);
  }
  break;
    // calculate normalized subgrid-diffusivity matrix
  case DRT::ELEMENTS::Condif3::calc_subgrid_diffusivity_matrix:
  {
    return DRT::ELEMENTS::Condif3ImplInterface::Impl(this)->Evaluate(
             this,
             params,
             discretization,
             lm,
             elemat1,
             elemat2,
             elevec1,
             elevec2,
             elevec3,
             mat,
             actmat);
  }
  break;
  // calculate flux
  case calc_condif_flux:
  {
    // get velocity values at the nodes
    // compare also with DRT::UTILS::ExtractMyValues()
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
    const int iel = NumNode();
    const int nsd=3;
    Epetra_SerialDenseVector evel(nsd*iel);
    DRT::UTILS::ExtractMyNodeBasedValues(this,evel,velocity);

    // need current values of transported scalar
    RefCountPtr<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==null) dserror("Cannot get state vector 'phinp'");

    // extract local values from the global vectors
    vector<double> myphinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

    // assure, that the values are in the same order as the element nodes
    for(int k=0;k<iel;++k)
    {
      Node* node = (Nodes())[k];
      vector<int> dof = discretization.Dof(node);
      int numdof = dof.size();
        // up to now, there's only one dof per node
      for (int i=0;i<numdof;++i)
      {
        if (dof[i]!=lm[k*numdof+i])
        {
          cout<<"dof[i]= "<<dof[i]<<"  lm[k*numdof+i]="<<lm[k*numdof+i]<<endl;
          dserror("Dofs are not in the same order as the element nodes. Implement some resorting!");
        }
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

    // set flag for type of scalar
    string scaltypestr=params.get<string>("problem type");
    int numscal = numdofpernode_;
    bool temperature = false;
    if (scaltypestr =="loma") temperature = true;

    double frt(0.0);
    if (scaltypestr =="elch") 
    {
      numscal -= 1; // ELCH case: last dof is for el. potential
      // get parameter F/RT
      frt = params.get<double>("frt");
    }

    // do a loop for systems of transported scalars
    for (int i = 0; i<numscal; ++i)
    {
      Epetra_SerialDenseMatrix eflux = CalculateFlux(myphinp,actmat,temperature,frt,evel,fluxtype,i);

      for (int k=0;k<iel;k++)
      { // form arithmetic mean of assembled nodal flux vectors
        // => factor is the number of adjacent elements for each node
        double factor = (Nodes()[k])->NumElement();
        elevec1[k*numdofpernode_+i]+=eflux(0,k)/factor;
        elevec2[k*numdofpernode_+i]+=eflux(1,k)/factor;
        elevec3[k*numdofpernode_+i]+=eflux(2,k)/factor;
      }
    } // loop over numdofpernode
  }
  break;
  // calculate mean temperature and density
  case calc_temp_and_dens:
  {
    // need current scalar and density vector
    RefCountPtr<const Epetra_Vector> phinp = discretization.GetState("phinp");
    RefCountPtr<const Epetra_Vector> densnp = discretization.GetState("densnp");
    if (phinp==null || densnp==null)
      dserror("Cannot get state vector 'phinp' and/or 'densnp'");

    // extract local values from the global vectors
    vector<double> myphinp(lm.size());
    vector<double> mydensnp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);
    DRT::UTILS::ExtractMyValues(*densnp,mydensnp,lm);

    // calculate temperature, density and domain integral
    CalculateTempAndDens(params,myphinp,mydensnp);
  }
  break;
  case calc_elch_kwok_error:
  {
    // check if length suffices
    if (elevec1.Length() < 1) dserror("Result vector too short");
    // determine errors
    return DRT::ELEMENTS::Condif3ImplInterface::Impl(this)->Evaluate(
             this,
             params,
             discretization,
             lm,
             elemat1,
             elemat2,
             elevec1,
             elevec2,
             elevec3,
             mat,
             actmat);
  }
  break;
  default:
    dserror("Unknown type of action for Condif3");
  } // end of switch(act)

  return 0;
} //DRT::ELEMENTS::Condif3::Evaluate


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
 |  calculate mass flux                              (private) gjb 06/08|
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix DRT::ELEMENTS::Condif3::CalculateFlux(
    vector<double>&           ephinp,
    struct _MATERIAL*         material,
    const bool                temperature,
    const double              frt,
    Epetra_SerialDenseVector& evel,
    Condif3::FluxType         fluxtype,
    const int&                dofindex
)
{
  /*------------------------------------------------- set element data */
  const int iel = NumNode();
  const DiscretizationType distype = Shape();
  const int nsd = 3;

  Epetra_SerialDenseMatrix xyze(nsd,iel);
  Epetra_SerialDenseMatrix flux(nsd,iel);

  // get node coordinates
  for (int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
    xyze(2,i)=Nodes()[i]->X()[2];
  }

  // get diffusivity
  double diffus(0.0);
  double valence(0.0);
  double diffus_valence_frt(0.0);

  if (material->mattyp == m_matlist)
  {
    const int matid = material->m.matlist->matids[dofindex];
    const _MATERIAL& singlemat =  DRT::Problem::Instance()->Material(matid-1);

    if (singlemat.mattyp == m_condif)
      diffus = singlemat.m.condif->diffusivity;
    else if (singlemat.mattyp == m_ion)
    {
      diffus = singlemat.m.ion->diffusivity;
      valence = singlemat.m.ion->valence;
      diffus_valence_frt = diffus*valence*frt;
    }
    else
      dserror("type of material found in material list is not supported.");
  }
  else if (material->mattyp == m_condif)
  {
    dsassert(numdofpernode_==1,"more than 1 dof per node for condif material"); // paranoia?
    if (temperature)
      diffus = material->m.condif->diffusivity/material->m.condif->shc;
    else
      diffus = material->m.condif->diffusivity;
  }
  else
    dserror("Material type is not supported");

  /*----------------------------------------- declaration of variables ---*/
  Epetra_SerialDenseVector        funct(iel);
  Epetra_SerialDenseMatrix        deriv(nsd,iel);
  static Epetra_SerialDenseMatrix xjm(nsd,nsd);
  Epetra_SerialDenseMatrix        derxy(nsd,iel);

  vector< LINALG::Matrix<3,1> > nodecoords;
  nodecoords = DRT::UTILS::getEleNodeNumbering_nodes_reference(distype);

  if ((int) nodecoords.size() != iel) dserror("number of nodes does not match");

  // loop over all nodes
  for (int iquad=0; iquad<iel; ++iquad)
  {
    // coordinates of the current integration point
    const double e1 = nodecoords[iquad](0);
    const double e2 = nodecoords[iquad](1);
    const double e3 = nodecoords[iquad](2);

    // shape functions and their derivatives
    DRT::UTILS::shape_function_3D(funct,e1,e2,e3,distype);
    DRT::UTILS::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

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

    // ---------------------------------------inverse of transposed jacobian
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

    // gradient of electric potential
    Epetra_SerialDenseVector gradpot(3);
    if (frt > 0.0) // ELCH
    {
      for (int k=0;k<iel;k++)
      {
        gradpot[0] += derxy(0,k)*ephinp[k*numdofpernode_+(numdofpernode_-1)];
        gradpot[1] += derxy(1,k)*ephinp[k*numdofpernode_+(numdofpernode_-1)];
        gradpot[2] += derxy(2,k)*ephinp[k*numdofpernode_+(numdofpernode_-1)];
      } /* end of loop over k */
    }

    const double ephinpatnode = ephinp[iquad*numdofpernode_+dofindex];
    // add different flux contributions as specified by user input
    switch (fluxtype)
    {
      case Condif3::totalflux:
        if (frt > 0.0) // ELCH
        {
          //migration flux terms
          flux(0,iquad)-=diffus_valence_frt*gradpot[0]*ephinpatnode;
          flux(1,iquad)-=diffus_valence_frt*gradpot[1]*ephinpatnode;
          flux(2,iquad)-=diffus_valence_frt*gradpot[2]*ephinpatnode;
        }
        //convective flux terms
        flux(0,iquad)+=evel[iquad*nsd]*ephinpatnode;
        flux(1,iquad)+=evel[1+iquad*nsd]*ephinpatnode;
        flux(2,iquad)+=evel[2+iquad*nsd]*ephinpatnode;
        // no break statement here!
      case Condif3::diffusiveflux:
        //diffusive flux terms
        for (int k=0;k<iel;k++)
        {
          flux(0,iquad)+=-diffus*derxy(0,k)*ephinp[k*numdofpernode_+dofindex];
          flux(1,iquad)+=-diffus*derxy(1,k)*ephinp[k*numdofpernode_+dofindex];
          flux(2,iquad)+=-diffus*derxy(2,k)*ephinp[k*numdofpernode_+dofindex];
        }
        break;
      case Condif3::noflux:
        dserror("received noflux flag inside CONDIF3 flux evaluation");
    };

  } // loop over nodes

  return flux;
} // Condif3::CalculateFlux


/*----------------------------------------------------------------------*
 |  calculate temperature, density and domain integral          vg 09/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif3::CalculateTempAndDens(
    ParameterList&            params,
    vector<double>&           ephinp,
    vector<double>&           edensnp
)
{
  // get variables for integrals
  double tempint = params.get<double>("temperature integral");
  double densint = params.get<double>("density integral");
  double domint  = params.get<double>("domain integral");

  /*------------------------------------------------- set element data */
  const int iel = NumNode();
  const DiscretizationType distype = Shape();
  const int nsd = 3;

  // get node coordinates
  Epetra_SerialDenseMatrix xyze(nsd,iel);
  for (int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
    xyze(2,i)=Nodes()[i]->X()[2];
  }

  /*----------------------------------------- declaration of variables ---*/
  Epetra_SerialDenseVector        funct(iel);
  Epetra_SerialDenseMatrix        deriv(nsd,iel);
  static Epetra_SerialDenseMatrix xjm(nsd,nsd);

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(SCATRA::get3DOptimalGaussrule(distype));

  // integration loop
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    // coordinates of the current integration point
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];
    const double e3 = intpoints.qxg[iquad][2];

    // shape functions
    DRT::UTILS::shape_function_3D(funct,e1,e2,e3,distype);
    DRT::UTILS::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

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

    // calculate integrals of temperature, density and domain
    for (int i=0; i<iel; i++)
    {
      tempint += fac*funct[i]*ephinp[i];
      densint += fac*funct[i]*edensnp[i];
      domint  += fac*funct[i];
    }
  } // loop over nodes

  // return variables for integrals
  params.set<double>("temperature integral",tempint);
  params.set<double>("density integral",densint);
  params.set<double>("domain integral",domint);

  return;
} // Condif3::CalculateTempAndDens


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
