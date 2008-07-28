/*!-----------------------------------------------------------------------------------------------------------
 \file truss3_evaluate.cpp
 \brief three dimensional total Lagrange truss element

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

 *-----------------------------------------------------------------------------------------------------------*/
#ifdef D_TRUSS3
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "truss3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

//externally defined structure for material data
extern struct _MATERIAL *mat;

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 08/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Truss3::Evaluate(ParameterList& params,
    DRT::Discretization& discretization,
    vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Truss3::ActionType act = Truss3::calc_none;
  // get the action required
  string action = params.get<string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff") act = Truss3::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff") act = Truss3::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Truss3::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass") act = Truss3::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass") act = Truss3::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = Truss3::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress") act = Truss3::calc_struct_stress;
  else if (action=="calc_struct_eleload") act = Truss3::calc_struct_eleload;
  else if (action=="calc_struct_fsiload") act = Truss3::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep") act = Truss3::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike") act = Truss3::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep") act = Truss3::calc_struct_reset_istep;
  else dserror("Unknown type of action for Truss3");

  switch(act)
  {
    /*in case that only linear stiffness matrix is required b3_nlstiffmass is called with zero dispalcement and 
     residual values*/
    case Truss3::calc_struct_linstiff:
    {
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");

    }
    break;

    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case Truss3::calc_struct_nlnstiffmass:
    case Truss3::calc_struct_nlnstifflmass:
    case Truss3::calc_struct_nlnstiff:
    case Truss3::calc_struct_internalforce:
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

      // get element velocities
      RefCountPtr<const Epetra_Vector> vel  = discretization.GetState("velocity");
      if (vel==null) dserror("Cannot get state vectors 'velocity'");
      vector<double> myvel(lm.size());
      DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

      if (act == Truss3::calc_struct_nlnstiffmass)
      t3_nlnstiffmass(mydisp,&elemat1,&elemat2,&elevec1);
      else if (act == Truss3::calc_struct_nlnstifflmass)
      {
        t3_nlnstiffmass(mydisp,&elemat1,&elemat2,&elevec1);
        // lump mass matrix (bborn 07/08)
        // the mass matrix is lumped anyway, cf #b3_nlnstiffmass
        //b3_lumpmass(&elemat2);
      }
      else if (act == Truss3::calc_struct_nlnstiff)
      t3_nlnstiffmass(mydisp,&elemat1,NULL,&elevec1);
      else if (act == Truss3::calc_struct_internalforce)
      t3_nlnstiffmass(mydisp,NULL,NULL,&elevec1);
    
    }
    break;
    case calc_struct_update_istep:
    case calc_struct_update_imrlike:
    {
      //nothing to do
    }
    break;
    case calc_struct_reset_istep:
    {
      //nothing to do
    }
    break;
    default:
    dserror("Unknown type of action for Truss3 %d", act);
  }
  return 0;

}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)                                       cyron 03/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Truss3::EvaluateNeumann(ParameterList& params,
    DRT::Discretization& discretization,
    DRT::Condition& condition,
    vector<int>& lm,
    Epetra_SerialDenseVector& elevec1)
{
  // get element displacements
  RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==null) dserror("Cannot get state vector 'displacement'");
  vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
  // get element velocities (UNCOMMENT IF NEEDED)
  /*
  RefCountPtr<const Epetra_Vector> vel  = discretization.GetState("velocity");
  if (vel==null) dserror("Cannot get state vectors 'velocity'");
  vector<double> myvel(lm.size());
  DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
  */

  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  // number of the load curve related with a specific line Neumann condition called
  if (curve) curvenum = (*curve)[0];
  // amplitude of load curve at current time called
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)//notation for this function similar to Crisfield, Volume 1;
  curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

  //jacobian determinant
  double det = lrefe_/2;

  // no. of nodes on this element; the following line is only valid for elements with constant number of 
  // degrees of freedom per node
  const int numdf = 6;
  const DiscretizationType distype = this->Shape();

  // gaussian points 
  const DRT::UTILS::IntegrationPoints1D intpoints = getIntegrationPoints1D(gaussrule_);

  //declaration of variable in order to store shape function
  Epetra_SerialDenseVector funct(NumNode());

  // get values and switches from the condition

  // onoff is related to the first 6 flags of a line Neumann condition in the input file;
  // value 1 for flag i says that condition is active for i-th degree of freedom
  const vector<int>* onoff = condition.Get<vector<int> >("onoff");
  // val is related to the 6 "val" fields after the onoff flags of the Neumann condition
  // in the input file; val gives the values of the force as a multiple of the prescribed load curve
  const vector<double>* val = condition.Get<vector<double> >("val");

  //integration loops	
  for (int ip=0; ip<intpoints.nquad; ++ip)
  {
    //integration points in parameter space and weights
    const double xi = intpoints.qxg[ip];
    const double wgt = intpoints.qwgt[ip];

    //evaluation of shape funcitons at Gauss points
    DRT::UTILS::shape_function_1D(funct,xi,distype);

    double fac=0;
    fac = wgt * det;

    // load vector ar
    double ar[numdf];
    // loop the dofs of a node

    for (int i=0; i<numdf; ++i)
    {
      ar[i] = fac * (*onoff)[i]*(*val)[i]*curvefac;
    }

    //sum up load components 
    for (int node=0; node<NumNode(); ++node)
    for (int dof=0; dof<numdf; ++dof)
    elevec1[node*numdf+dof] += funct[node] *ar[dof];

  } // for (int ip=0; ip<intpoints.nquad; ++ip)
  
  return 0;
}

/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 08/08|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_nlnstiffmass( vector<double>& disp,
    Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector* force)
{

  
  //current position of nodal degrees of freedom
  BlitzMat6x2 xcurr;
  
  //difference between coordinates of both nodes in current configuration, x21' Crisfield  Vol. 2 equ. (17.66a) and (17.72)
  BlitzVec3 x21;


  //nodal coordinates in current position
  for (int k=0; k<2; ++k) //looping over number of nodes
  {
    for (int j=0; j<3; ++j) 
    {
      xcurr(j,k)   = Nodes()[k]->X()[j] + disp[k*6+j]; //translational DOF
      xcurr(j+3,k) = disp[k*6+j+3]; //rotational DOF
    }
  }



  /* read material parameters using structure _MATERIAL which is defined by inclusion of      /
   / "../drt_lib/drt_timecurve.H"; note: material parameters have to be read in the evaluation /
   / function instead of e.g. Truss3_input.cpp or within the Truss3Register class since it is not/
   / sure that structure _MATERIAL is declared within those scopes properly whereas it is within/
   / the evaluation functions */

  // get the material law
  MATERIAL* currmat = &(mat[material_-1]);
  double ym;
  double sm;
  double density;

  //assignment of material parameters; only St.Venant material is accepted for this truss 
  switch(currmat->mattyp)
  {
    case m_stvenant:// only linear elastic material supported
    {
      ym = currmat->m.stvenant->youngs;
      sm = ym / (2*(1 + currmat->m.stvenant->possionratio));
      density = currmat->m.stvenant->density;
    }
    break;
    default:
    dserror("unknown or improper type of material law");
  }
  


  //computing global internal forces, Crisfield Vol. 2, equation (17.79)
  //note: X = [-I 0; -S -I; I 0; -S I] with -S = T^t; and S = S(x21)/2;
  if (force != NULL)
  {
    (*force).Size(6);
    for (int i=0; i<6; ++i)
    {
      (*force)(i)   += 0;
    }
  }

  //computing linear stiffness matrix
  if (stiffmatrix != NULL)
  {   
    //setting up basis of stiffness matrix according to Crisfield, Vol. 2, equation (17.81)   
    (*stiffmatrix).Shape(6,6);  


      //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
      //correctly or not: the function b3_nlnstiff_approx(mydisp) calculates the stiffness matrix approximated by
      //finite differences and finally the relative error is printed; in case that there is no significant error in any
      //element no printout is thrown
      //activating this part of code also the function b3_nlnstiff_approx(mydisp) has to be activated both in Truss3.H
      //and Truss3_evaluate.cpp
      /*
      if(Id() == 8) //limiting the following tests to certain element numbers
      {
       Epetra_SerialDenseMatrix stiff_approx;
       Epetra_SerialDenseMatrix stiff_relerr;
       stiff_approx.Shape(12,12);
       stiff_relerr.Shape(12,12);      
       double h_rel = 1e-7;
       int outputflag = 0;
       stiff_approx = b3_nlnstiff_approx(xcurr, h_rel, *force);
       
       for(int line=0; line<12; line++)
       {
         for(int col=0; col<12; col++)
         {
           if( fabs( (*stiffmatrix)(line,col) ) > h_rel)
             stiff_relerr(line,col)= abs( ((*stiffmatrix)(line,col) - stiff_approx(line,col))/(*stiffmatrix)(line,col) );
           else
           {
             if( fabs( stiff_approx(line,col) ) < h_rel*1000)
               stiff_relerr(line,col) = 0;
             else
               stiff_relerr(line,col)= abs( ((*stiffmatrix)(line,col) - stiff_approx(line,col))/(*stiffmatrix)(line,col) );
           }
           //suppressing small entries whose effect is only confusing
           if (stiff_relerr(line,col)<h_rel*100)
             stiff_relerr(line,col)=0;
           //there is no error if an entry is nan e.g. due to dirichlet boundary conditions
           if ( isnan( stiff_relerr(line,col) ) )
             stiff_relerr(line,col)=0;
           if (stiff_relerr(line,col)>0)
             outputflag = 1;  
          }
       
           if(outputflag ==1)
           {
             std::cout<<"\n\n acutally calculated stiffness matrix"<< *stiffmatrix;
             std::cout<<"\n\n approximated stiffness matrix"<< stiff_approx;    
             std::cout<<"\n\n rel error stiffness matrix"<< stiff_relerr;
           }    
         }
        } 
        */
   
  }
  
  /*calculating mass matrix; this truss3 element includes only a lumped mass matrix where for torsion and
   * bending the same moments of inertia are assumed; for slender trusss the influence of rotational moments
   * of inertia is dilute so that often rotational inertia is just set to zero; since within this code at one
   * point an LU-decomposition of the mass matrix is carried out this was avoided in the here implemented truss3
   * element and instead the above described simplification was assumed; Iyy_ = Izz_ was assumed for similar 
   * reasons */
  if (massmatrix != NULL)
  {
    (*massmatrix).Shape(6,6);
    for (int i=0; i<6; ++i)
    {
      (*massmatrix)(i,i) = 0.5*density*lrefe_*crosssec_;
    }
  }
  
  return;
} // DRT::ELEMENTS::Truss3::b3_nlnstiffmass


// lump mass matrix
void DRT::ELEMENTS::Truss3::t3_lumpmass(Epetra_SerialDenseMatrix* emass)
{
  // lump mass matrix
  if (emass != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (int c=0; c<(*emass).N(); ++c) // parse columns
    {
      double d = 0.0;
      for (int r=0; r<(*emass).M(); ++r) // parse rows
      {
        d += (*emass)(r,c); // accumulate row entries
        (*emass)(r,c) = 0.0;
      }
      (*emass)(c,c) = d; // apply sum of row entries on diagonal
    }
  }
}

//the following function can be activated in order to find bugs; it calculates a finite difference
//approximation of the nonlinear stiffness matrix; activate the follwing block for bug fixing only
/*
Epetra_SerialDenseMatrix DRT::ELEMENTS::Truss3::t3_nlnstiff_approx(BlitzMat6x2 xcurr, double h_rel, Epetra_SerialDenseVector force)
{
      //computing global internal forces, Crisfield Vol. 2, equation (17.79)
      //note: X = [-I 0; -S -I; I 0; -S I] with -S = T^t; and S = S(x21)/2;
      force_aux.Size(12);
      for (int u=0; u<3; ++u)
      {
        force_aux(u)   -= stressn_aux(u);
        force_aux(u+3) -= stressm_aux(u);
        force_aux(u+6) += stressn_aux(u);
        force_aux(u+9) += stressm_aux(u);
        
        for (int j=0; j<3; ++j)
        {      
          force_aux(u+3) -= stressn_aux(j)*spinx21_aux(u,j);      
          force_aux(u+9) -= stressn_aux(j)*spinx21_aux(u,j);       
        }
      }
      
      for(int u = 0;u<12;u++)
      {
        stiff_approx(u,i+k*6)= ( force_aux(u) - force(u) )/h_rel;
      }
   
    }
  }
   
  return stiff_approx;
} // DRT::ELEMENTS::Truss3::b3_nlnstiff_approx

*/

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_TRUSS3
