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

#include "torsion2.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_mat/spring.H"

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                cyron 02/10|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Torsion2::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization,
    std::vector<int>& lm,
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
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      // get residual displacements
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (res==Teuchos::null) dserror("Cannot get state vectors 'residual displacement'");
      std::vector<double> myres(lm.size());
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
    
      
      /*
      //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
      //correctly or not by means of a numerically approximated stiffness matrix
      //The code block will work for all higher order elements.
      //if(Id() == 3) //limiting the following tests to certain element numbers
      {
        //assuming the same number of DOF for all nodes
        int numdof = NumDofPerNode(*(Nodes()[0]));
        int nnode  = NumNode();
        
        //variable to store numerically approximated stiffness matrix
        Epetra_SerialDenseMatrix stiff_approx;
        stiff_approx.Shape(numdof*nnode,numdof*nnode);


        //relative error of numerically approximated stiffness matrix
        Epetra_SerialDenseMatrix stiff_relerr;
        stiff_relerr.Shape(numdof*nnode,numdof*nnode);

        //characteristic length for numerical approximation of stiffness
        double h_rel = 1e-6;                                                       // evtl. Größenordnung anpassen

        //flag indicating whether approximation leads to significant relative error
        int outputflag = 0;

        //calculating strains in new configuration
        for(int i=0; i<numdof; i++) //for all dof
        {
          for(int k=0; k<nnode; k++)//for all nodes
          {

            Epetra_SerialDenseVector force_aux;
            force_aux.Size(numdof*nnode);

            //create new displacement and velocity vectors in order to store artificially modified displacements
            //vector<double> vel_aux(myvel);
            std::vector<double> disp_aux(mydisp);

            //modifying displacement artificially (for numerical derivative of internal forces):
            disp_aux[numdof*k + i] += h_rel;
            //vel_aux[numdof*k + i]  += h_rel / params.get<double>("delta time",0.01);

            t2_nlnstiffmass(disp_aux,NULL,NULL,&force_aux);  

            //computing derivative d(fint)/du numerically by finite difference
            for(int u = 0 ; u < numdof*nnode ; u++ )
              stiff_approx(u,k*numdof+i) = ( pow(force_aux[u],2) - pow(elevec1(u),2) )/ (h_rel * (force_aux[u] + elevec1(u) ) );

          } //for(int k=0; k<nnode; k++)//for all nodes
        } //for(int i=0; i<numdof; i++) //for all dof


        for(int line=0; line<numdof*nnode; line++)
        {
          for(int col=0; col<numdof*nnode; col++)
          {
            stiff_relerr(line,col)= fabs( ( pow(elemat1(line,col),2) - pow(stiff_approx(line,col),2) )/ ( (elemat1(line,col) + stiff_approx(line,col)) * elemat1(line,col) ));

            //suppressing small entries whose effect is only confusing and NaN entires (which arise due to zero entries)
            if ( fabs( stiff_relerr(line,col) ) < h_rel*500 || isnan( stiff_relerr(line,col)) || elemat1(line,col) == 0) //isnan = is not a number
              stiff_relerr(line,col) = 0;

            if ( stiff_relerr(line,col) > 0)
              outputflag = 1;
              
          } //for(int col=0; col<numdof*nnode; col++)
        } //for(int line=0; line<numdof*nnode; line++)

        if(outputflag ==1)
        {
          std::cout<<"\n\n acutally calculated stiffness matrix in Element "<<Id()<<": "<< elemat1;
          std::cout<<"\n\n approximated stiffness matrix in Element "<<Id()<<": "<< stiff_approx;
          std::cout<<"\n\n rel error stiffness matrix in Element "<<Id()<<": "<< stiff_relerr;
        }

      } //end of section in which numerical approximation for stiffness matrix is computed
      
      
      */      
      
      
    }
    break;
    case calc_struct_update_istep:
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
int DRT::ELEMENTS::Torsion2::EvaluateNeumann(Teuchos::ParameterList&           params,
                                            DRT::Discretization&      discretization,
                                            DRT::Condition&           condition,
                                            std::vector<int>&         lm,
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

int DRT::ELEMENTS::Torsion2::EvaluatePTC(Teuchos::ParameterList&            params,
                                         Epetra_SerialDenseMatrix& elemat1)
{
  std::cout<<"\nno PTC implemented for torsion2 element\n";

  return 0;
} //DRT::ELEMENTS::Torsion2::EvaluatePTC

/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 02/10|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion2::t2_nlnstiffmass(std::vector<double>&      disp,
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
  if(( ( aux(0)*aux(2) + aux(1)*aux(3) )/ lcurr(0) / lcurr(1) ) > 1)
  {
    if(( ( aux(0)*aux(2) + aux(1)*aux(3) )/ lcurr(0) / lcurr(1) -1)<10e-7)
    {
      thetacurr=0;
    }
    else
      std::cout<<"\n calculation of the angle failed";
  }
         
  if (( aux(0)*aux(3) - aux(1)*aux(2) ) <0)
    thetacurr*=-1;
  
  deltatheta=thetacurr - theta_;
  
  // spring constant from material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double spring = 0.0;

  // assignment of material parameters; only spring material is accepted for this element
  switch(currmat->MaterialType())
  {
    case INPAR::MAT::m_spring: // only elastic spring supported
    {
      const MAT::Spring* actmat = static_cast<const MAT::Spring*>(currmat.get());
      spring = actmat->Stiffness();
    }
    break;
    default:
    dserror("unknown or improper type of material law");
  }
  
  //__________________________________________________________________________________________
  // bending potential quadratic
  
  if (bendingpotential_==quadratic)
  {
  
	  //global internal forces (equation 2.12)
	  if (force != NULL)
	    for (int i=0; i<6; ++i)
	      (*force)(i) =  spring*deltatheta*grtheta(i);
	  
	  //stiffness matrix (equation 2.15)
	  if (stiffmatrix != NULL)
	  {
	    for (int i=0; i<6; ++i)
	    {
	      for (int j=0; j<6; ++j)
	        (*stiffmatrix)(i,j)=spring*grtheta(i)*grtheta(j);    //first part of the stiffness matrix
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
	        A(i,j)=spring*deltatheta/pow(lcurr(0),2)*( tmp(i,j)-2*perp(  i)*unit(  j) );
	        B(i,j)=spring*deltatheta/pow(lcurr(1),2)*( tmp(i,j)-2*perp(2+i)*unit(2+j) );
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
  
  } //bending potential quadratic
  
  //==================================================================================================================================================
  //____________________________________________________________________________________
  else
	//bending potetial cosine
	if (bendingpotential_==cosine) 
	{	    
	    double dotprod=0.0;
			  
			  for (int j=0; j<2; ++j)
			    dotprod +=  aux(j) * aux(2+j);
	  

			  //variation of theta (equation 3.4)
			  LINALG::Matrix<6,1> grtheta2;
			  
			  for (int j=0; j<2; ++j)
			  {    //virtual displacement of node 1 and 3, in vector b (3.4)
			    grtheta2(  j)  = -aux(2+j)/lcurr(0)/lcurr(1)+dotprod*aux(  j)/pow(lcurr(0),3)/lcurr(1);   
			    grtheta2(4+j)  =  aux(  j)/lcurr(0)/lcurr(1)-dotprod*aux(2+j)/lcurr(0)/pow(lcurr(1),3);
			  }
			  
			  for (int j=0; j<2; ++j)     //virtual displacement of node 2, in vector b (3.4)
			    grtheta2(2+j) = -grtheta2(j)-grtheta2(j+4);	
			  
			  //auxiliary matrix for stiffness matrix (equation 3.9 and equation 3.10)
			  LINALG::Matrix<4,2> A;
			  LINALG::Matrix<4,2> B;
			
			  for(int j=0; j<2;++j)
			  {
			    for(int i=0; i<2; ++i)  //auxiliary matrix A is A1 above A2
			    {
			      A(  j,i)= aux(2+j)*aux(  i)/pow(lcurr(0),3)/lcurr(1) + aux(2+i)*aux(  j)/pow(lcurr(0),3)/lcurr(1)-3*dotprod*aux(  j)*aux(  i)/pow(lcurr(0),5)/lcurr(1); //3->2? 3*dotprod
			      A(2+j,i)= aux(  j)*aux(2+i)/lcurr(0)/pow(lcurr(1),3) + aux(  i)*aux(2+j)/lcurr(0)/pow(lcurr(1),3)-3*dotprod*aux(2+j)*aux(2+i)/lcurr(0)/pow(lcurr(1),5);
			    }
			    A(j,j)  += dotprod/pow(lcurr(0),3)/lcurr(1);  //aus A1, Anteil auf Diagonale
			    A(2+j,j)+= dotprod/lcurr(0)/pow(lcurr(1),3);  //aus A2
			  }
			      
			  for(int j=0; j<2; ++j)  //auviliary matrix B is C2 above C1
			  {
			    for(int i=0; i<2;++i)
			    {
			      B(  j,i)= -aux(  j)*aux(  i)/pow(lcurr(0),3)/lcurr(1) - aux(2+i)*aux(2+j)/lcurr(0)/pow(lcurr(1),3)+dotprod*aux(2+j)*aux(  i)/pow(lcurr(0),3)/pow(lcurr(1),3);
			      B(2+j,i)= -aux(2+j)*aux(2+i)/lcurr(0)/pow(lcurr(1),3) - aux(  i)*aux(  j)/pow(lcurr(0),3)/lcurr(1)+dotprod*aux(  j)*aux(2+i)/pow(lcurr(0),3)/pow(lcurr(1),3);
			    }
			    B(j,j)  += 1/lcurr(0)/lcurr(1);
			    B(2+j,j)+= 1/lcurr(0)/lcurr(1);
			  }

        //internal forces (equation 3.16, same as for theta almost zero)
        if (force != NULL)
        {
          for (int j=0; j<6; ++j)
            (*force)(j) = -spring*grtheta2(j);      //-
        }


        //stiffness matrix (equation 3.17, same as for theta almost zero)
        if (stiffmatrix != NULL)
        {
          for (int j=0; j<2; ++j)
          {
            for (int i=0; i<2; ++i)
            {
              (*stiffmatrix)(  j,  i) =-spring*( -A(j  ,i) );
              (*stiffmatrix)(  j,2+i) =-spring*(  A(j  ,i)+B(j+2,i) );
              (*stiffmatrix)(  j,4+i) =-spring*( -B(j+2,i) );
              (*stiffmatrix)(2+j,  i) =-spring*(  A(j  ,i)+B(j  ,i) );
              (*stiffmatrix)(2+j,2+i) =-spring*( -A(j  ,i)-A(j+2,i)-B(j ,i)-B(j+2,i) );
              (*stiffmatrix)(2+j,4+i) =-spring*(  A(j+2,i)+B(j+2,i) );
              (*stiffmatrix)(4+j,  i) =-spring*( -B(j  ,i) );
              (*stiffmatrix)(4+j,2+i) =-spring*(  A(j+2,i)+B(j  ,i) );
              (*stiffmatrix)(4+j,4+i) =-spring*( -A(j+2,i) );
            }
          }
        }
		
		
	}//bending potetial cosine
  else 
	  std::cout<<"\n No such bending potential. Possible bending potentials: \n cosine \n quadratic"<<endl;
 
  return;
} // DRT::ELEMENTS::Torsion2::t2_nlnstiffmass

