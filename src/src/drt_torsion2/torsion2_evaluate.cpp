/*!-----------------------------------------------------------------------------------------------------------
 \file torsion2_evaluate.cpp
 \brief three dimensional total Lagrange truss element (can be connected to beam3 elements and adapts assembly automatically according to the thereby changed number of nodal degrees of freedom)

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

 *-----------------------------------------------------------------------------------------------------------*/
#ifdef D_TORSION2
#ifdef CCADISCRET

#include "torsion2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/linalg_fixedsizematrix.H"

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                cyron 02/10|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Torsion2::Evaluate(ParameterList& params,
    DRT::Discretization& discretization,
    vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Torsion2::ActionType act = Torsion2::calc_none;
  // get the action required
  string action = params.get<string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff") act = Torsion2::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff") act = Torsion2::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Torsion2::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass") act = Torsion2::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass") act = Torsion2::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = Torsion2::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress") act = Torsion2::calc_struct_stress;
  else if (action=="calc_struct_eleload") act = Torsion2::calc_struct_eleload;
  else if (action=="calc_struct_fsiload") act = Torsion2::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep") act = Torsion2::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike") act = Torsion2::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep") act = Torsion2::calc_struct_reset_istep;
  else if (action=="postprocess_stress") act = Torsion2::postprocess_stress;
  else if (action=="calc_struct_ptcstiff") act = Torsion2::calc_struct_ptcstiff;
  else
    {
      cout<<action<<endl;
      dserror("Unknown type of action for Torsion2");
    }

  switch(act)
  {
    case Torsion2::calc_struct_ptcstiff:
    {
      EvaluatePTC(params, elemat1);
    }
    break;
    /*in case that only linear stiffness matrix is required b2_nlstiffmass is called with zero dispalcement and 
     residual values*/
    case Torsion2::calc_struct_linstiff:
    {
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");
    }
    break;

    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case Torsion2::calc_struct_nlnstiffmass:
    case Torsion2::calc_struct_nlnstifflmass:
    case Torsion2::calc_struct_nlnstiff:
    case Torsion2::calc_struct_internalforce:
    {
      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacemnet and residual values for each degree of freedom
      //
      // get element displacements
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      // get residual displacements
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (res==null) dserror("Cannot get state vectors 'residual displacement'");
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      
      // for engineering strains instead of total lagrange use t2_nlnstiffmass2
      if (act == Torsion2::calc_struct_nlnstiffmass)
        t2_nlnstiffmass(mydisp,&elemat1,&elemat2,&elevec1);
      else if (act == Torsion2::calc_struct_nlnstifflmass)
        t2_nlnstiffmass(mydisp,&elemat1,&elemat2,&elevec1);
      else if (act == Torsion2::calc_struct_nlnstiff)
        t2_nlnstiffmass(mydisp,&elemat1,NULL,&elevec1);
      else if (act == Torsion2::calc_struct_internalforce)
        t2_nlnstiffmass(mydisp,NULL,NULL,&elevec1);

      // for engineering strains instead of total lagrange use t2_nlnstiffmass2
      if (act == Torsion2::calc_struct_nlnstiffmass)
        t2_nlnstiffmass(mydisp,&elemat1,&elemat2,&elevec1);
      else if (act == Torsion2::calc_struct_nlnstifflmass)
        t2_nlnstiffmass(mydisp,&elemat1,&elemat2,&elevec1);
      else if (act == Torsion2::calc_struct_nlnstiff)
        t2_nlnstiffmass(mydisp,&elemat1,NULL,&elevec1);
      else if (act == Torsion2::calc_struct_internalforce)
        t2_nlnstiffmass(mydisp,NULL,NULL,&elevec1);
    
    }
    break;
    case calc_struct_update_istep:
    case calc_struct_update_imrlike:
    case calc_struct_reset_istep:
    case calc_struct_stress:
    break;
    case postprocess_stress:
    {
      //no stress calculation for postprocess. Does not really make sense!
      dserror("No stress output for Torsion2!");      
    }
    break;  
    default:
    dserror("Unknown type of action for Torsion2 %d", act);
  }
  return 0;

}
/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)                                       cyron 02/10|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Torsion2::EvaluateNeumann(ParameterList& params,
                                            DRT::Discretization&      discretization,
                                            DRT::Condition&           condition,
                                            vector<int>&              lm,
                                            Epetra_SerialDenseVector& elevec1,
                                            Epetra_SerialDenseMatrix* elemat1)
{
  /*torsion spring assumed to be infinitesimally small and thus no Neumann boundary conditions
   * can be assigned to this element*/
  return 0;
}


/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate PTC damping (public)                                                                  cyron 02/10|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Torsion2::EvaluatePTC(ParameterList& params,
                                      Epetra_SerialDenseMatrix& elemat1)
{
  std::cout<<"\nno PTC implemented for torsion2 element\n";

  return 0;
} //DRT::ELEMENTS::Torsion2::EvaluatePTC

/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 02/10|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion2::t2_nlnstiffmass( vector<double>& disp,
    Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector* force)
{
  //current node position (first entries 0 .. 1 for first node, 2 ..3 for second node , 4 .. 5 for third node)
  LINALG::Matrix<6,1> xcurr;
    
  //current nodal position
  for (int j=0; j<2; ++j)
  {
    xcurr(j  )   = Nodes()[0]->X()[j] + disp[    j];  //first node
    xcurr(j+2)   = Nodes()[1]->X()[j] + disp[2 + j];  //second node
    xcurr(j+4)   = Nodes()[2]->X()[j] + disp[4 + j];  //third node
  }
  
  //auxiliary vector for both internal force and stiffness matrix
  LINALG::Matrix<4,1> aux;
  for (int i=0; i<4; ++i) 
    aux(  i) = xcurr(2+i)-xcurr(  i);
  
  //current length of vectors 1-->2  und 2-->3
  LINALG::Matrix<2,1> lcurr;
  for (int j=0; j<2; ++j)
    lcurr(j) = sqrt( pow(aux(2*j),2) + pow(aux(2*j+1),2) );
  
  //unit vektor of the trusses
  LINALG::Matrix<4,1> unit;
  for (int i=0; i<2; ++i)
    for (int j=0; j<2; ++j)
      unit(j+i*2) = aux(j+i*2)/lcurr(i);
  
  //perpendicular vektor of the trusses
  LINALG::Matrix<4,1> perp;
  for (int i=0; i<2; ++i)
  {
    perp(  i*2)=-unit(1+i*2);
    perp(1+i*2)= unit(  i*2);
  }
     
  //computing the gradient of theta, which is equivalent to the variaton of theta (equation 2.10)
  LINALG::Matrix<6,1> grtheta;  
  for (int j=0; j<2; ++j)
  {
    grtheta(j)  =  perp(j)   / lcurr(0);    //virtual displacement of node 1
    grtheta(j+4)=  perp(j+2) / lcurr(1);    //virtual displacement of node 3
    grtheta(j+2)= -grtheta(j)-grtheta(j+4);     //virtual displacement of node 2
  }
  
  //current angle theta (equation 2.13)
  double deltatheta=0.0;
  double thetacurr=0.0;
  
  thetacurr=acos( ( aux(0)*aux(2) + aux(1)*aux(3) ) / lcurr(0) / lcurr(1) );
  if(( ( aux(0)*aux(2) + aux(1)*aux(3) )/ lcurr(0) / lcurr(1) ) > 1){
    if(( ( aux(0)*aux(2) + aux(1)*aux(3) )/ lcurr(0) / lcurr(1) -1)<10e-7){
      thetacurr=0;
    }
    else
      std::cout<<"\n calculation of the angle failed";
  }
         
  if (( aux(0)*aux(3) - aux(1)*aux(2) ) <0)
    thetacurr*=-1;
  
  deltatheta=thetacurr - theta_;
    
  //global internal forces (equation 2.12)
  if (force != NULL)
    for (int i=0; i<6; ++i)
      (*force)(i) =  springconstant_*deltatheta*grtheta(i);
  
  //stiffness matrix (equation 2.15)
  if (stiffmatrix != NULL){
    for (int i=0; i<6; ++i){
      for (int j=0; j<6; ++j)
        (*stiffmatrix)(i,j)=springconstant_*grtheta(i)*grtheta(j);    //first part of the stiffness matrix
    }
    
    LINALG::Matrix<2,2> tmp;
    tmp( 0, 0)= 0.0;
    tmp( 0, 1)=-1.0;
    tmp( 1, 0)= 1.0;
    tmp( 1, 1)= 0.0;
    
    LINALG::Matrix<2,2> A;    //equation 2.16
    LINALG::Matrix<2,2> B;    //equation 2.17
    for (int i=0;i<2;++i)
    {
      for (int j=0;j<2;++j)
      {
        A(i,j)=springconstant_*deltatheta/pow(lcurr(0),2)*( tmp(i,j)-2*perp(  i)*unit(  j) );
        B(i,j)=springconstant_*deltatheta/pow(lcurr(1),2)*( tmp(i,j)-2*perp(2+i)*unit(2+j) );
      }
    }
        
    //second part of the stiffness matrix (equation 2.18)
    for (int i=0; i<2; ++i)
    {
      for (int j=0; j<2; ++j)
      {
        (*stiffmatrix)(  i,  j)+=-A(i,j);
        (*stiffmatrix)(  i,2+j)+= A(i,j);
        (*stiffmatrix)(2+i,  j)+= A(i,j);
        (*stiffmatrix)(2+i,2+j)+=-A(i,j)+B(i,j);
        (*stiffmatrix)(2+i,4+j)+=-B(i,j);
        (*stiffmatrix)(4+i,2+j)+=-B(i,j);
        (*stiffmatrix)(4+i,4+j)+= B(i,j);
      }
    }
  }//stiffness matrix
   

  return;
} // DRT::ELEMENTS::Torsion2::t2_nlnstiffmass

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_TORSION2
