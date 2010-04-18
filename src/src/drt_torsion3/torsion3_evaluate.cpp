/*!-----------------------------------------------------------------------------------------------------------
 \file torsion3_evaluate.cpp
 \brief three dimensional total Lagrange truss element (can be connected to beam3 elements and adapts assembly automatically according to the thereby changed number of nodal degrees of freedom)

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

 *-----------------------------------------------------------------------------------------------------------*/
#ifdef D_TORSION3
#ifdef CCADISCRET

#include "torsion3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_fixedsizematrix.H"


/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 08/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Torsion3::Evaluate(ParameterList& params,
    DRT::Discretization& discretization,
    vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{ 
  DRT::ELEMENTS::Torsion3::ActionType act = Torsion3::calc_none;
  // get the action required
  string action = params.get<string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff") act = Torsion3::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff") act = Torsion3::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Torsion3::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass") act = Torsion3::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass") act = Torsion3::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = Torsion3::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress") act = Torsion3::calc_struct_stress;
  else if (action=="calc_struct_eleload") act = Torsion3::calc_struct_eleload;
  else if (action=="calc_struct_fsiload") act = Torsion3::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep") act = Torsion3::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike") act = Torsion3::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep") act = Torsion3::calc_struct_reset_istep;
  else if (action=="postprocess_stress") act = Torsion3::postprocess_stress;
  else if (action=="calc_struct_ptcstiff") act = Torsion3::calc_struct_ptcstiff;
  else 
  {
    cout<<action<<endl;
    dserror("Unknown type of action for Torsion3");
  }


  switch(act)
  {
    case Torsion3::calc_struct_ptcstiff:
    {
      //nothing to do
    }
    break;
    /*in case that only linear stiffness matrix is required b3_nlstiffmass is called with zero dispalcement and 
     residual values*/
    case Torsion3::calc_struct_linstiff:
    {
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");

    }
    break;

    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case Torsion3::calc_struct_nlnstiffmass:
    case Torsion3::calc_struct_nlnstifflmass:
    case Torsion3::calc_struct_nlnstiff:
    case Torsion3::calc_struct_internalforce:
    {
      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacemnet and residual values for each degree of freedom
      //
      // get element displcements
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      // get residual displacements
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (res==null) dserror("Cannot get state vectors 'residual displacement'");
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
           

      // for engineering strains instead of total lagrange use t3_nlnstiffmass2
      if (act == Torsion3::calc_struct_nlnstiffmass)
        t3_nlnstiffmass(params,mydisp,&elemat1,&elemat2,&elevec1);
      else if (act == Torsion3::calc_struct_nlnstifflmass)
        t3_nlnstiffmass(params,mydisp,&elemat1,&elemat2,&elevec1);
      else if (act == Torsion3::calc_struct_nlnstiff)
        t3_nlnstiffmass(params,mydisp,&elemat1,NULL,&elevec1);
      else if (act == Torsion3::calc_struct_internalforce)
        t3_nlnstiffmass(params,mydisp,NULL,NULL,&elevec1);   
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
      dserror("No stress output for Torsion3!");      
    }
    break;  
    default:
    dserror("Unknown type of action for Torsion3 %d", act);
  }
  return 0;

}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)                                       cyron 03/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Torsion3::EvaluateNeumann(ParameterList& params,
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

/*--------------------------------------------------------------------------------------*
 | evaluate nonlinear stiffness matrix and internal forces                    cyron 03/10|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3::t3_nlnstiffmass(ParameterList&            params,
                                              vector<double>&           disp,
                                              Epetra_SerialDenseMatrix* stiffmatrix,
                                              Epetra_SerialDenseMatrix* massmatrix,
                                              Epetra_SerialDenseVector* force)
{
  /*first displacement vector is modified for proper element evaluation in case of periodic boundary conditions; in case that
   *no periodic boundary conditions are to be applied the following code line may be ignored or deleted*/
  NodeShift<3,3>(params,disp);
  
  //current node position (first entries 0,1,2 for first node, 3,4,5 for second node , 6,7,8 for third node)
  LINALG::Matrix<9,1> xcurr;
    
  //current nodal position
  for (int j=0; j<3; ++j)
  {
    xcurr(j  )   = Nodes()[0]->X()[j] + disp[    j];  //first node
    xcurr(j+3)   = Nodes()[1]->X()[j] + disp[3 + j];  //second node
    xcurr(j+6)   = Nodes()[2]->X()[j] + disp[6 + j];  //third node
  }
    
  //auxiliary vector for both internal force and stiffness matrix
  LINALG::Matrix<6,1> aux;
  for (int j=0; j<6; ++j)
    aux(j) = xcurr(j+3)-xcurr(j);
  
  //current length of vectors 1-->2  and 2-->3
  LINALG::Matrix<2,1> lcurr;
  for (int j=0; j<2; ++j)
    lcurr(j) = sqrt( pow(aux(3*j),2) + pow(aux(3*j+1),2) + pow(aux(3*j+2),2) );
  
  //computing the change of angle theta (equation 3.3)
  double deltatheta=0.0;
  double dotprod=0.0;
  double s=0.0;
  
  for (int j=0; j<3; ++j)
    dotprod +=  aux(j) * aux(3+j);

  s = dotprod/lcurr(0)/lcurr(1);
  deltatheta=acos(s);
  

  //variation of theta (equation 3.4)
  LINALG::Matrix<9,1> grtheta;
  
  for (int j=0; j<3; ++j){    //virual displacement of node 1 and 3
    grtheta(  j)  = -aux(3+j)/lcurr(0)/lcurr(1)+dotprod*aux(  j)/pow(lcurr(0),3)/lcurr(1);    
    grtheta(6+j)  =  aux(  j)/lcurr(0)/lcurr(1)-dotprod*aux(3+j)/lcurr(0)/pow(lcurr(1),3);
  }
  
  for (int j=0; j<3; ++j)     //virtual displacement of node 2
    grtheta(3+j) = -grtheta(j)-grtheta(j+6);
  
  //auxiliary matrix for stiffness matrix (equation 3.9 and equation 3.10)
  LINALG::Matrix<6,3> A;
  LINALG::Matrix<6,3> B;
  
  for(int j=0; j<3;++j){
    for(int i=0; i<3; ++i){
      A(  j,i)= aux(3+j)*aux(  i)/pow(lcurr(0),3)/lcurr(1) + aux(3+i)*aux(  j)/pow(lcurr(0),3)/lcurr(1)-3*dotprod*aux(  j)*aux(  i)/pow(lcurr(0),5)/lcurr(1);
      A(3+j,i)= aux(  j)*aux(3+i)/lcurr(0)/pow(lcurr(1),3) + aux(  i)*aux(3+j)/lcurr(0)/pow(lcurr(1),3)-3*dotprod*aux(3+j)*aux(3+i)/lcurr(0)/pow(lcurr(1),5);
    }
    A(j,j)  += dotprod/pow(lcurr(0),3)/lcurr(1);
    A(3+j,j)+= dotprod/lcurr(0)/pow(lcurr(1),3);
  }
      
  for(int j=0; j<3; ++j)
  {
    for(int i=0; i<3;++i){
      B(  j,i)= -aux(  j)*aux(  i)/pow(lcurr(0),3)/lcurr(1) - aux(3+i)*aux(3+j)/lcurr(0)/pow(lcurr(1),3)+dotprod*aux(3+j)*aux(  i)/pow(lcurr(0),3)/pow(lcurr(1),3);
      B(3+j,i)= -aux(3+j)*aux(3+i)/lcurr(0)/pow(lcurr(1),3) - aux(  i)*aux(  j)/pow(lcurr(0),3)/lcurr(1)+dotprod*aux(  j)*aux(3+i)/pow(lcurr(0),3)/pow(lcurr(1),3);
    }
    B(j,j)  += 1/lcurr(0)/lcurr(1);
    B(3+j,j)+= 1/lcurr(0)/lcurr(1);
  }
  
  
  if((1-s)<=0.00000000001){
    //calculation of the linear stiffness matrix (if the angle between the trusses is almost zero)
    
    //internal forces (equation 3.17)
    if (force != NULL){
      for (int j=0; j<9; ++j)
        (*force)(j) = -springconstant_*grtheta(j);      
    }
    
    
    //stiffness matrix (equation 3.19)
    if (stiffmatrix != NULL) {
      for (int j=0; j<3; ++j) { 
        for (int i=0; i<3; ++i){
          (*stiffmatrix)(  j,  i) =-springconstant_*( -A(j  ,i) );
          (*stiffmatrix)(  j,3+i) =-springconstant_*(  A(j  ,i)+B(j+3,i) );
          (*stiffmatrix)(  j,6+i) =-springconstant_*( -B(j+3,i) );
          (*stiffmatrix)(3+j,  i) =-springconstant_*(  A(j  ,i)+B(j  ,i) );
          (*stiffmatrix)(3+j,3+i) =-springconstant_*( -A(j  ,i)-A(j+3,i)-B(j ,i)-B(j+3,i) );
          (*stiffmatrix)(3+j,6+i) =-springconstant_*(  A(j+3,i)+B(j+3,i) );
          (*stiffmatrix)(6+j,  i) =-springconstant_*( -B(j  ,i) );
          (*stiffmatrix)(6+j,3+i) =-springconstant_*(  A(j+3,i)+B(j  ,i) );
          (*stiffmatrix)(6+j,6+i) =-springconstant_*( -A(j+3,i) );
        }
      }
    } 
  } // if (linear stiffness matrix)
  
  else{
    //internal forces (equation 3.6)
    if (force != NULL){
      for (int j=0; j<9; ++j)
        (*force)(j) = -1/sqrt(1-s*s)*deltatheta*springconstant_*grtheta(j);     
    }
    
    
    //stiffness matrix (equation 3.8)
    if (stiffmatrix != NULL) {
      for (int i=0; i<9; ++i){
        for (int j=0; j<9; ++j)
          (*stiffmatrix)(i,j)=1/(1-s*s)*springconstant_*grtheta(i)*grtheta(j)
                              -s/pow(sqrt(1-s*s),3)*deltatheta*springconstant_*grtheta(i)*grtheta(j);
      }
    
      for (int j=0; j<3; ++j) { //equation 3.13
        for (int i=0; i<3; ++i){
          (*stiffmatrix)(  j,  i)+=springconstant_*deltatheta*(-1/sqrt(1-s*s))*( -A(j  ,i) );
          (*stiffmatrix)(  j,3+i)+=springconstant_*deltatheta*(-1/sqrt(1-s*s))*(  A(j  ,i)+B(j+3,i) );
          (*stiffmatrix)(  j,6+i)+=springconstant_*deltatheta*(-1/sqrt(1-s*s))*( -B(j+3,i) );
          (*stiffmatrix)(3+j,  i)+=springconstant_*deltatheta*(-1/sqrt(1-s*s))*(  A(j  ,i)+B(j  ,i) );
          (*stiffmatrix)(3+j,3+i)+=springconstant_*deltatheta*(-1/sqrt(1-s*s))*( -A(j  ,i)-A(j+3,i)-B(j ,i)-B(j+3,i) );
          (*stiffmatrix)(3+j,6+i)+=springconstant_*deltatheta*(-1/sqrt(1-s*s))*(  A(j+3,i)+B(j+3,i) );
          (*stiffmatrix)(6+j,  i)+=springconstant_*deltatheta*(-1/sqrt(1-s*s))*( -B(j  ,i) );
          (*stiffmatrix)(6+j,3+i)+=springconstant_*deltatheta*(-1/sqrt(1-s*s))*(  A(j+3,i)+B(j  ,i) );
          (*stiffmatrix)(6+j,6+i)+=springconstant_*deltatheta*(-1/sqrt(1-s*s))*( -A(j+3,i) );
        }
      }
    } //stiffness matrix
  
  } //else


}//stiffmatrix

/*-----------------------------------------------------------------------------------------------------------*
 | shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions; if two   |
 | nodes within one element are separated by a periodic boundary, one of them is shifted such that the final |
 | distance in R^3 is the same as the initial distance in the periodic space; the shift affects computation  |
 | on element level within that very iteration step, only (no change in global variables performed)          |                                 |
 |                                                                                       (public) cyron 10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim> //number of nodes, number of dimensions
inline void DRT::ELEMENTS::Torsion3::NodeShift(ParameterList& params,  //!<parameter list
                                            vector<double>& disp) //!<element disp vector
{
  /*get number of degrees of freedom per node; note: the following function assumes the same number of degrees
   *of freedom for each element node*/
  int numdof = NumDofPerNode(*(Nodes()[0]));



  /*only if periodic boundary conditions are in use, i.e. params.get<double>("PeriodLength",0.0) > 0.0, this
   * method has to change the displacement variables*/
  if(params.get<double>("PeriodLength",0.0) > 0.0)
    //loop through all nodes except for the first node which remains fixed as reference node
    for(int i=1;i<nnode;i++)
    {
      for(int dof= ndim - 1; dof > -1; dof--)
      {
        /*if the distance in some coordinate direction between some node and the first node becomes smaller by adding or subtracting
         * the period length, the respective node has obviously been shifted due to periodic boundary conditions and should be shifted
         * back for evaluation of element matrices and vectors; this way of detecting shifted nodes works as long as the element length
         * is smaller than half the periodic length*/
        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) + params.get<double>("PeriodLength",0.0) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
        {
          disp[numdof*i+dof] += params.get<double>("PeriodLength",0.0);

          /*the upper domain surface orthogonal to the z-direction is subject to shear Dirichlet boundary condition; the lower surface
           *is fixed by DBC. To avoid problmes when nodes exit the domain through the upper z-surface and reenter through the lower
           *z-surface, the shear has to be substracted from nodal coordinates in that case */
          if(dof == 2)
            disp[numdof*i+params.get<int>("OSCILLDIR",-1)] += params.get<double>("SHEARAMPLITUDE",0.0)*DRT::Problem::Instance()->Curve(params.get<int>("CURVENUMBER",-1)-1).f(params.get<double>("total time",0.0));
        }

        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - params.get<double>("PeriodLength",0.0) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
        {
          disp[numdof*i+dof] -= params.get<double>("PeriodLength",0.0);

          /*the upper domain surface orthogonal to the z-direction is subject to shear Dirichlet boundary condition; the lower surface
           *is fixed by DBC. To avoid problmes when nodes exit the domain through the lower z-surface and reenter through the upper
           *z-surface, the shear has to be added to nodal coordinates in that case */
          if(dof == 2)
            disp[numdof*i+params.get<int>("OSCILLDIR",-1)] -= params.get<double>("SHEARAMPLITUDE",0.0)*DRT::Problem::Instance()->Curve(params.get<int>("CURVENUMBER",-1)-1).f(params.get<double>("total time",0.0));
        }
      }
    }

return;

}//DRT::ELEMENTS::Torsion3::NodeShift

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_TORSION3
