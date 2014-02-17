/*!-----------------------------------------------------------------------------------------------------------
 \file truss2_evaluate.cpp
 \brief two dimensional total Lagrange truss element 

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

 *-----------------------------------------------------------------------------------------------------------*/
#include "truss2.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 02/10|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Truss2::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization,
    std::vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{ 
  DRT::ELEMENTS::Truss2::ActionType act = Truss2::calc_none;
  // get the action required
  std::string action = params.get<std::string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff") act = Truss2::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff") act = Truss2::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Truss2::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass") act = Truss2::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass") act = Truss2::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = Truss2::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress") act = Truss2::calc_struct_stress;
  else if (action=="calc_struct_eleload") act = Truss2::calc_struct_eleload;
  else if (action=="calc_struct_fsiload") act = Truss2::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep") act = Truss2::calc_struct_update_istep;
  else if (action=="calc_struct_reset_istep") act = Truss2::calc_struct_reset_istep;
  else if (action=="postprocess_stress") act = Truss2::postprocess_stress;
  else if (action=="calc_struct_ptcstiff") act = Truss2::calc_struct_ptcstiff;
  else 
    {
      std::cout<<action<<std::endl;
      dserror("Unknown type of action for Truss2");
    }
  
  /*number of degrees of freedom actually assigned to the discretization by the first node
   *(allows connectin Truss2 and beam2 directly)*/
  int ActNumDof0 = discretization.NumDof(Nodes()[0]);

  switch(act)
  {
    case Truss2::calc_struct_ptcstiff:
    {
      EvaluatePTC(params, elemat1,ActNumDof0);
    }
    break;
    /*in case that only linear stiffness matrix is required b3_nlstiffmass is called with zero dispalcement and 
     residual values*/
    case Truss2::calc_struct_linstiff:
    {
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");

    }
    break;

    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case Truss2::calc_struct_nlnstiffmass:
    case Truss2::calc_struct_nlnstifflmass:
    case Truss2::calc_struct_nlnstiff:
    case Truss2::calc_struct_internalforce:
    {
      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacemnet and residual values for each degree of freedom
      //
      // get element displcements
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      // get residual displacements
      Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (res==Teuchos::null) dserror("Cannot get state vectors 'residual displacement'");
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
           
      // get element velocities (UNCOMMENT IF NEEDED)
      /*
      Teuchos::RCP<const Epetra_Vector> vel  = discretization.GetState("velocity");
      if (vel==Teuchos::null) dserror("Cannot get state vectors 'velocity'");
      std::vector<double> myvel(lm.size());
      DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
      */
      
      // for engineering strains instead of total lagrange use t2_nlnstiffmass2
      if (act == Truss2::calc_struct_nlnstiffmass)
      t2_nlnstiffmass(mydisp,&elemat1,&elemat2,&elevec1);
      else if (act == Truss2::calc_struct_nlnstifflmass)
      {
        t2_nlnstiffmass(mydisp,&elemat1,&elemat2,&elevec1);
        // lump mass matrix (bborn 07/08)
        // the mass matrix is lumped anyway, cf #b3_nlnstiffmass
        //b3_lumpmass(&elemat2);
      }
      else if (act == Truss2::calc_struct_nlnstiff)
      t2_nlnstiffmass(mydisp,&elemat1,NULL,&elevec1);
      else if (act == Truss2::calc_struct_internalforce)
      t2_nlnstiffmass(mydisp,NULL,NULL,&elevec1);
    
    }
    break;
    case calc_struct_update_istep:
    {
      //nothing to do
    }
    break;
    case calc_struct_reset_istep:
    {
      //nothing to do
    }
    break;
    case calc_struct_stress:
    {
      //no stress calculation implemented! Do not crash simulation and just keep quiet!
    }
    break;
    case postprocess_stress:
    {
      //no stress calculation for postprocess. Does not really make sense!
      dserror("No stress output for Truss2!");      
    }
    break;  
    default:
    dserror("Unknown type of action for Truss2 %d", act);
  }
  return 0;

}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)                                       cyron 02/10|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Truss2::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization,
    DRT::Condition& condition,
    std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
  // get element displacements
    Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
    if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement'");
    std::vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);


    // find out whether we will use a time curve
    bool usetime = true;
    const double time = params.get("total time",-1.0);
    if (time<0.0) usetime = false;

    // find out whether we will use a time curve and get the factor
    const std::vector<int>* curve = condition.Get<std::vector<int> >("curve");
    int curvenum = -1;
    // number of the load curve related with a specific line Neumann condition called
    if (curve) curvenum = (*curve)[0];
    // amplitude of load curve at current time called
    double curvefac = 1.0;
    if (curvenum>=0 && usetime)//notation for this function similar to Crisfield, Volume 1;
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

    //jacobian determinant
    double det = lrefe_/2;

    const DiscretizationType distype = this->Shape();

    // gaussian points
    const DRT::UTILS::IntegrationPoints1D intpoints(gaussrule_);

    //declaration of variable in order to store shape function
    Epetra_SerialDenseVector funct(NumNode());

    // get values and switches from the condition

    // onoff is related to the first 6 flags of a line Neumann condition in the input file;
    // value 1 for flag i says that condition is active for i-th degree of freedom
    const std::vector<int>* onoff = condition.Get<std::vector<int> >("onoff");
    // val is related to the 6 "val" fields after the onoff flags of the Neumann condition
    // in the input file; val gives the values of the force as a multiple of the prescribed load curve
    const std::vector<double>* val = condition.Get<std::vector<double> >("val");

    //integration loops
    for (int ip=0; ip<intpoints.nquad; ++ip)
    {
      //integration points in parameter space and weights
      const double xi = intpoints.qxg[ip][0];
      const double wgt = intpoints.qwgt[ip];

      //evaluation of shape funcitons at Gauss points
      DRT::UTILS::shape_function_1D(funct,xi,distype);

      double fac=0;
      fac = wgt * det;

      /*load vector ar; regardless of the actual number of degrees of freedom active with respect to this
       *element or certain nodes of it the vector val has always the lengths 6 and in order to deal with
       *possibly different numbers of acutally used DOF we always loop through all the 6*/
      double ar[6];
      // loop the dofs of a node

      for (int i = 0; i < 6; ++i)
        ar[i] = fac * (*onoff)[i]*(*val)[i]*curvefac;
    
      for (int dof=0; dof < 2; ++dof)
      {
        //computing entries for first node
        elevec1[dof] += funct[0] *ar[dof];
        //computing entries for first node
        elevec1[2 + dof] += funct[1] *ar[dof];
      }
        
    } // for (int ip=0; ip<intpoints.nquad; ++ip)

    return 0;
}



/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate PTC damping (public)                                                                  cyron 02/10|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Truss2::EvaluatePTC(Teuchos::ParameterList& params,
                                      Epetra_SerialDenseMatrix& elemat1,
                                      int& ActNumDof0)
{
  std::cout<<"\nno PTC implemented for truss2\n";

  return 0;
} //DRT::ELEMENTS::Truss2::EvaluatePTC

/*--------------------------------------------------------------------------------------*
 | switch between kintypes                                                      cyron 02/10|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss2::t2_nlnstiffmass( std::vector<double>& disp,
    Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector* force)
{
  switch(kintype_)
  {
  case tr2_totlag:
    t2_nlnstiffmass_totlag(disp,stiffmatrix,massmatrix,force);
  break;
  case tr2_engstrain:
    t2_nlnstiffmass_engstr(disp,stiffmatrix,massmatrix,force);
  break;
  default:
    dserror("Unknown type kintype_ for Truss2");
  }

  return;
}


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                  cyron 02/10|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss2::t2_nlnstiffmass_totlag( std::vector<double>& disp,
    Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector* force)
{     
  //current node position (first entries 0 ..1  for first node, 2 ..3 for second node)
  LINALG::Matrix<4,1> xcurr;
  
  /*current nodal displacement (first entries 0 .. 1 for first node, 2 ..3 for second node) compared
   * to reference configuration; note: in general this is not equal to the values in disp since the 
   * latter one referes to a nodal displacement compared to a reference configuration before the first
   * time step whereas the following variable referes to the displacement with respect to a reference
   * configuration which may have been set up at any point of time during the simulation (usually this
   * is only important if an element has been added to the discretization after the start of the simulation)*/
  LINALG::Matrix<4,1> ucurr; 
  
  //Green-Lagrange strain
  double epsilon;
  
  //auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  LINALG::Matrix<4,1> aux;

  //current nodal position
  for (int j=0; j<2; ++j) 
  {
    xcurr(j  )   = Nodes()[0]->X()[j] + disp[  j]; //first node
    xcurr(j+2)   = Nodes()[1]->X()[j] + disp[2 + j]; //second node
  }
  
  //current displacement = current position - reference position
  ucurr  = xcurr;
  ucurr -= X_;
  
   
  //computing auxiliary vector aux = N^T_{,xi} * N_{,xi} * xcurr
  aux(0) = 0.25 * (xcurr(0) - xcurr(2));
  aux(1) = 0.25 * (xcurr(1) - xcurr(3));
  aux(2) = 0.25 * (xcurr(2) - xcurr(0));
  aux(3) = 0.25 * (xcurr(3) - xcurr(1));
  
  
  //calculating strain epsilon from node position by scalar product:
  //epsilon = (xrefe + 0.5*ucurr)^T * N_{,s}^T * N_{,s} * d
  epsilon = 0;
  epsilon += (X_(0) + 0.5*ucurr(0)) * (ucurr(0) - ucurr(2));
  epsilon += (X_(1) + 0.5*ucurr(1)) * (ucurr(1) - ucurr(3));
  epsilon += (X_(2) + 0.5*ucurr(2)) * (ucurr(2) - ucurr(0));
  epsilon += (X_(3) + 0.5*ucurr(3)) * (ucurr(3) - ucurr(1));
  
  epsilon /= lrefe_*lrefe_;
  
 
  /* read material parameters using structure _MATERIAL which is defined by inclusion of      /
   / "../drt_lib/drt_timecurve.H"; note: material parameters have to be read in the evaluation /
   / function instead of e.g. Truss2_input.cpp or within the Truss2Register class since it is not/
   / sure that structure _MATERIAL is declared within those scopes properly whereas it is within/
   / the evaluation functions */

  // get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double ym = 0;
  double density = 0;

  //assignment of material parameters; only St.Venant material is accepted for this truss 
  switch(currmat->MaterialType())
  {
    case INPAR::MAT::m_stvenant:// only linear elastic material supported
    {
      const MAT::StVenantKirchhoff* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      ym = actmat->Youngs();
      density = actmat->Density();
    }
    break;
    default:
    dserror("unknown or improper type of material law");
  }
  

  //computing global internal forces
  if (force != NULL)
  {
    for (int i=0; i<2; ++i)
     (*force)(i) = (4*ym*crosssec_*epsilon/lrefe_) * aux(i);
    
    for (int i=0; i<2; ++i)
     (*force)(2 + i) = (4*ym*crosssec_*epsilon/lrefe_) * aux(i+2);    
  }
  

  //computing linear stiffness matrix
  if (stiffmatrix != NULL)
  {      
    for (int i=0; i<2; ++i)
    { 
        //stiffness entries for first node
        (*stiffmatrix)(i     ,i    )   =  (ym*crosssec_*epsilon/lrefe_);
        (*stiffmatrix)(i     ,2 + i)   = -(ym*crosssec_*epsilon/lrefe_);
        //stiffness entries for second node
        (*stiffmatrix)(i + 2 ,i + 2)   =  (ym*crosssec_*epsilon/lrefe_);
        (*stiffmatrix)(i + 2 ,i    )   = -(ym*crosssec_*epsilon/lrefe_);
    }
       
    for (int i=0; i<4; ++i)
      for (int j=0; j<4; ++j)
          (*stiffmatrix)(i,j) += (16*ym*crosssec_/pow(lrefe_,3))*aux(i)*aux(j);     
  }
  
  //calculating consistent mass matrix
  if (massmatrix != NULL)
  {
    for (int i=0; i<2; ++i)
    {
      (*massmatrix)(i    ,i    ) = density*lrefe_*crosssec_ / 3;
      (*massmatrix)(i + 2,i + 2) = density*lrefe_*crosssec_ / 3;
      (*massmatrix)(i    ,i + 2) = density*lrefe_*crosssec_ / 6;
      (*massmatrix)(i + 2,i    ) = density*lrefe_*crosssec_ / 6;
    }
  }
  
  return;
} // DRT::ELEMENTS::Truss2::t2_nlnstiffmass


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 02/10|
 | engineering strain measure, large displacements and rotations                                              |
  *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss2::t2_nlnstiffmass_engstr( std::vector<double>& disp,
    Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector* force
    )
{
  //current node position (first entries 0 .. 1 for first node, 2 ..3 for second node)
  LINALG::Matrix<4,1> xcurr;
  
  //Green-Lagrange strain
  double epsilon;
  
  //auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  LINALG::Matrix<4,1> aux;

  //current nodal position (first
  for (int j=0; j<2; ++j) 
  {
    xcurr(j  )   = Nodes()[0]->X()[j] + disp[  j]; //first node
    xcurr(j+2)   = Nodes()[1]->X()[j] + disp[2 + j]; //second node
  }
  
  //computing auxiliary vector aux = 4.0*N^T_{,xi} * N_{,xi} * xcurr
  aux(0) = (xcurr(0) - xcurr(2));
  aux(1) = (xcurr(1) - xcurr(3));
  aux(2) = (xcurr(2) - xcurr(0));
  aux(3) = (xcurr(3) - xcurr(1));
  
  
  double lcurr = sqrt(pow(aux(0),2)+pow(aux(1),2));

  //calculating strain epsilon from node position by scalar product:
  epsilon = (lcurr-lrefe_)/lrefe_;

  /* read material parameters using structure _MATERIAL which is defined by inclusion of      /
   / "../drt_lib/drt_timecurve.H"; note: material parameters have to be read in the evaluation /
   / function instead of e.g. Truss2_input.cpp or within the Truss2Register class since it is not/
   / sure that structure _MATERIAL is declared within those scopes properly whereas it is within/
   / the evaluation functions */

  // get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double ym = 0;
  double density = 0;

  //assignment of material parameters; only St.Venant material is accepted for this truss 
  switch(currmat->MaterialType())
  {
  case INPAR::MAT::m_stvenant:// only linear elastic material supported
    {
      const MAT::StVenantKirchhoff* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      ym = actmat->Youngs();
      density = actmat->Density();
    }
    break;
    default:
    dserror("unknown or improper type of material law");
  }
  
  // resulting force scaled by current length
  double forcescalar=(ym*crosssec_*epsilon)/lcurr;
  
  //computing global internal forces
  if (force != NULL)
  {  
    //node 1
    for (int i=0; i<2; ++i)
     (*force)(i) = forcescalar * aux(i);
    //node 2
    for (int i=0; i<2; ++i)
     (*force)(2 + i) = forcescalar * aux(i+2);
  }
  
  //computing linear stiffness matrix
  if (stiffmatrix != NULL)
  {        
    for (int i=0; i<2; ++i)
    { 
        //stiffness entries for first node
        (*stiffmatrix)(i    ,i    )   =  forcescalar;
        (*stiffmatrix)(i    ,2 + i)   = -forcescalar;
        //stiffness entries for second node
        (*stiffmatrix)(i + 2,i + 2)   =  forcescalar;
        (*stiffmatrix)(i + 2,i    )   = -forcescalar;
    }

    for (int i=0; i<4; ++i)
      for (int j=0; j<4; ++j)
        (*stiffmatrix)(i,j) += (ym*crosssec_/pow(lcurr,3))*aux(i)*aux(j);
 
  }
  
  //calculating consistent mass matrix
  if (massmatrix != NULL)
  {
    for (int i=0; i<2; ++i)
    {
      (*massmatrix)(i    ,i    ) = density*lrefe_*crosssec_ / 3;
      (*massmatrix)(i + 2,i + 2) = density*lrefe_*crosssec_ / 3;
      (*massmatrix)(i    ,i + 2) = density*lrefe_*crosssec_ / 6;
      (*massmatrix)(i + 2,i    ) = density*lrefe_*crosssec_ / 6;
    }
  }
  
  return;
} // DRT::ELEMENTS::Truss2::bt_nlnstiffmass2


// lump mass matrix
void DRT::ELEMENTS::Truss2::t2_lumpmass(Epetra_SerialDenseMatrix* emass)
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
        (*emass)(r,c) = 0;
      }
      (*emass)(c,c) = d; // apply sum of row entries on diagonal
    }
  }
}

