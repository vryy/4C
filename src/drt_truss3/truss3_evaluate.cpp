/*!-----------------------------------------------------------------------------------------------------------
 \file truss3_evaluate.cpp
 \brief three dimensional total Lagrange truss element (can be connected to beam3 elements and adapts assembly automatically according to the thereby changed number of nodal degrees of freedom)

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
  else if (action=="postprocess_stress") act = Truss3::postprocess_stress;
  else if (action=="calc_stat_force_damp") act = Truss3::calc_stat_force_damp;
  else if (action=="calc_struct_ptcstiff") act = Truss3::calc_struct_ptcstiff;
  else 
    {
      cout<<action<<endl;
      dserror("Unknown type of action for Truss3");
    }

  switch(act)
  {
    case Truss3::calc_stat_force_damp:
    case Truss3::calc_struct_ptcstiff:
    {   
      //Truss3 element does not provide any tool for calculating stochastical forces and related damping
      
      //Truss3 element does'nt need any special ptc tools to allow stable implicit dynamics with acceptable time step size
    }
    break;
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
      
      /*number of degrees of freedom actually assigned to the discretization by the first node
       *(allows flexible handling of the case that the truss3 element is connected directly to
       *a beam3 element)*/
      int ActNumDof0 = discretization.NumDof(Nodes()[0]);
      
      // get element velocities (UNCOMMENT IF NEEDED)
      /*
      RefCountPtr<const Epetra_Vector> vel  = discretization.GetState("velocity");
      if (vel==null) dserror("Cannot get state vectors 'velocity'");
      vector<double> myvel(lm.size());
      DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
      */
      
      // for engineering strains instead of total lagrange use t3_nlnstiffmass2
      if (act == Truss3::calc_struct_nlnstiffmass)
      t3_nlnstiffmass(mydisp,&elemat1,&elemat2,&elevec1,ActNumDof0);
      else if (act == Truss3::calc_struct_nlnstifflmass)
      {
        t3_nlnstiffmass(mydisp,&elemat1,&elemat2,&elevec1,ActNumDof0);
        // lump mass matrix (bborn 07/08)
        // the mass matrix is lumped anyway, cf #b3_nlnstiffmass
        //b3_lumpmass(&elemat2);
      }
      else if (act == Truss3::calc_struct_nlnstiff)
      t3_nlnstiffmass(mydisp,&elemat1,NULL,&elevec1,ActNumDof0);
      else if (act == Truss3::calc_struct_internalforce)
      t3_nlnstiffmass(mydisp,NULL,NULL,&elevec1,ActNumDof0);
    
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
    case calc_struct_stress:
    {
      //no stress calculation implemented! Do not crash simulation and just keep quiet!
    }
    break;
    case postprocess_stress:
    {
      //no stress calculation for postprocess. Does not really make sense!
      dserror("No stress output for Truss3!");      
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
  //first the actual number of DOF of each node is detected
  int ActNumDof0 = discretization.NumDof(Nodes()[0]);
  int ActNumDof1 = discretization.NumDof(Nodes()[1]);
  
  // get element displacements
  RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==null) dserror("Cannot get state vector 'displacement'");
  vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);


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

    /*load vector ar; regardless of the actual number of degrees of freedom active with respect to this
     *element or certain nodes of it the vector val has always the lengths 6 and in order to deal with
     *possibly different numbers of acutally used DOF we always loop through all the 6*/
    double ar[6];
    // loop the dofs of a node

    for (int i = 0; i < 6; ++i)
    {
      ar[i] = fac * (*onoff)[i]*(*val)[i]*curvefac;
    }
    
    //computing entries for first node
    for (int dof=0; dof < ActNumDof0; ++dof)
        elevec1[dof] += funct[0] *ar[dof];
    
    //computing entries for second node
    for (int dof=0; dof < ActNumDof1; ++dof)
      elevec1[ActNumDof0 + dof] += funct[1] *ar[dof];

  } // for (int ip=0; ip<intpoints.nquad; ++ip)
  
  return 0;
}

/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 08/08|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_nlnstiffmass( vector<double>& disp,
    Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector* force,
    int& ActNumDof0)
{     
  //current node position (first entries 0 .. 2 for first node, 3 ..5 for second node)
  BlitzVec6 xcurr;
  
  //Green-Lagrange strain
  double epsilon;
  
  //auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  BlitzVec6 aux;

  //current nodal position (first
  for (int j=0; j<3; ++j) 
  {
    xcurr(j  )   = Nodes()[0]->X()[j] + disp[  j]; //first node
    xcurr(j+3)   = Nodes()[1]->X()[j] + disp[ActNumDof0 + j]; //second node
  }
  
  //computing auxiliary vector aux = N^T_{,xi} * N_{,xi} * xcurr
  aux(0) = 0.25 * (xcurr[0] - xcurr[3]);
  aux(1) = 0.25 * (xcurr[1] - xcurr[4]);
  aux(2) = 0.25 * (xcurr[2] - xcurr[5]);
  aux(3) = 0.25 * (xcurr[3] - xcurr[0]);
  aux(4) = 0.25 * (xcurr[4] - xcurr[1]);
  aux(5) = 0.25 * (xcurr[5] - xcurr[2]);
  
  //calculating strain epsilon from node position by scalar product:
  //epsilon = (xrefe + 0.5*disp)^T * N_{,s}^T * N_{,s} * d
  epsilon = 0;
  epsilon += (Nodes()[0]->X()[0] + 0.5*disp[0]) * (disp[0] - disp[3]);
  epsilon += (Nodes()[0]->X()[1] + 0.5*disp[1]) * (disp[1] - disp[4]);
  epsilon += (Nodes()[0]->X()[2] + 0.5*disp[2]) * (disp[2] - disp[5]);
  epsilon += (Nodes()[1]->X()[0] + 0.5*disp[3]) * (disp[3] - disp[0]);
  epsilon += (Nodes()[1]->X()[1] + 0.5*disp[4]) * (disp[4] - disp[1]);
  epsilon += (Nodes()[1]->X()[2] + 0.5*disp[5]) * (disp[5] - disp[2]);
  epsilon /= lrefe_*lrefe_;


  /* read material parameters using structure _MATERIAL which is defined by inclusion of      /
   / "../drt_lib/drt_timecurve.H"; note: material parameters have to be read in the evaluation /
   / function instead of e.g. Truss3_input.cpp or within the Truss3Register class since it is not/
   / sure that structure _MATERIAL is declared within those scopes properly whereas it is within/
   / the evaluation functions */

  // get the material law
  MATERIAL* currmat = &(mat[material_-1]);
  double ym = 0;
  double sm = 0;
  double density = 0;

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

  //computing global internal forces
  if (force != NULL)
  {
    for (int i=0; i<3; ++i)
     (*force)(i) = (4*ym*crosssec_*epsilon/lrefe_) * aux(i);
    
    for (int i=0; i<3; ++i)
     (*force)(ActNumDof0 + i) = (4*ym*crosssec_*epsilon/lrefe_) * aux(i+3);
  }

  //computing linear stiffness matrix
  if (stiffmatrix != NULL)
  {      
    for (int i=0; i<3; ++i)
    { 
        //stiffness entries for first node
        (*stiffmatrix)(i              ,i             )   =  (ym*crosssec_*epsilon/lrefe_);
        (*stiffmatrix)(i              ,ActNumDof0 + i)   = -(ym*crosssec_*epsilon/lrefe_);
        //stiffness entries for second node
        (*stiffmatrix)(i + ActNumDof0 ,i + ActNumDof0)   =  (ym*crosssec_*epsilon/lrefe_);
        (*stiffmatrix)(i + ActNumDof0 ,i             )   = -(ym*crosssec_*epsilon/lrefe_);
    }

    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
      {
        //node 1
        (*stiffmatrix)(i              ,j            ) += (16*ym*crosssec_/pow(lrefe_,3))*aux(i)*aux(j); 
        //node 2
        (*stiffmatrix)(i + ActNumDof0 ,j +ActNumDof0) += (16*ym*crosssec_/pow(lrefe_,3))*aux(i+3)*aux(j+3);
      }
  }
  
  //calculating consistent mass matrix
  if (massmatrix != NULL)
  {
    for (int i=0; i<3; ++i)
    {
      (*massmatrix)(i             ,i             ) = density*lrefe_*crosssec_ / 3;
      (*massmatrix)(i + ActNumDof0,i + ActNumDof0) = density*lrefe_*crosssec_ / 3;
      (*massmatrix)(i             ,i + ActNumDof0) = density*lrefe_*crosssec_ / 6;
      (*massmatrix)(i + ActNumDof0,i             ) = density*lrefe_*crosssec_ / 6;
    }
  }
  
  return;
} // DRT::ELEMENTS::Truss3::t3_nlnstiffmass


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                      tk 10/08|
 | linear(!) strain measure, large displacements and rotations                                                |
  *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_nlnstiffmass2( vector<double>& disp,
    Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector* force,
    int& ActNumDof0)
{
  //current node position (first entries 0 .. 2 for first node, 3 ..5 for second node)
  BlitzVec6 xcurr;
  
  //Green-Lagrange strain
  double epsilon;
  
  //auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  BlitzVec6 aux;

  //current nodal position (first
  for (int j=0; j<3; ++j) 
  {
    xcurr(j  )   = Nodes()[0]->X()[j] + disp[  j]; //first node
    xcurr(j+3)   = Nodes()[1]->X()[j] + disp[ActNumDof0 + j]; //second node
  }
  
  //computing auxiliary vector aux = 4.0*N^T_{,xi} * N_{,xi} * xcurr
  aux(0) = (xcurr[0] - xcurr[3]);
  aux(1) = (xcurr[1] - xcurr[4]);
  aux(2) = (xcurr[2] - xcurr[5]);
  aux(3) = (xcurr[3] - xcurr[0]);
  aux(4) = (xcurr[4] - xcurr[1]);
  aux(5) = (xcurr[5] - xcurr[2]);
  
  double lcurr = sqrt(pow(aux(0),2)+pow(aux(1),2)+pow(aux(2),2));
  
  //calculating strain epsilon from node position by scalar product:
  epsilon = (lcurr-lrefe_)/lrefe_;

  /* read material parameters using structure _MATERIAL which is defined by inclusion of      /
   / "../drt_lib/drt_timecurve.H"; note: material parameters have to be read in the evaluation /
   / function instead of e.g. Truss3_input.cpp or within the Truss3Register class since it is not/
   / sure that structure _MATERIAL is declared within those scopes properly whereas it is within/
   / the evaluation functions */

  // get the material law
  MATERIAL* currmat = &(mat[material_-1]);
  double ym = 0;
  double density = 0;

  //assignment of material parameters; only St.Venant material is accepted for this truss 
  switch(currmat->mattyp)
  {
    case m_stvenant:// only linear elastic material supported
    {
      ym = currmat->m.stvenant->youngs;
      density = currmat->m.stvenant->density;
    }
    break;
    default:
    dserror("unknown or improper type of material law");
  }

  //computing global internal forces
  if (force != NULL)
  {  
    double forcescalar=(ym*crosssec_*epsilon)/lcurr;
    //node 1
    for (int i=0; i<3; ++i)
     (*force)(i) = forcescalar * aux(i);
    //node 2
    for (int i=0; i<3; ++i)
     (*force)(ActNumDof0 + i) = forcescalar * aux(i+3);
  }
  
  //computing linear stiffness matrix
  if (stiffmatrix != NULL)
  {      
    for (int i=0; i<3; ++i)
    { 
        //stiffness entries for first node
        (*stiffmatrix)(i              ,i             )   =  (ym*crosssec_*epsilon/lcurr);
        (*stiffmatrix)(i              ,i + ActNumDof0)   = -(ym*crosssec_*epsilon/lcurr);
        //stiffness entries for second node
        (*stiffmatrix)(i + ActNumDof0 ,i + ActNumDof0)   =  (ym*crosssec_*epsilon/lcurr);
        (*stiffmatrix)(i + ActNumDof0 ,i             )   = -(ym*crosssec_*epsilon/lcurr);
    }

    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
      {
        //node 1
        (*stiffmatrix)(i              ,j            ) += (ym*crosssec_/pow(lcurr,3))*aux(i)*aux(j); 
        //node 2
        (*stiffmatrix)(i + ActNumDof0 ,j +ActNumDof0) += (ym*crosssec_/pow(lcurr,3))*aux(i+3)*aux(j+3);
      }
  }
  
  //calculating consistent mass matrix
  if (massmatrix != NULL)
  {
    for (int i=0; i<3; ++i)
    {
      (*massmatrix)(i             ,i             ) = density*lrefe_*crosssec_ / 3;
      (*massmatrix)(i + ActNumDof0,i + ActNumDof0) = density*lrefe_*crosssec_ / 3;
      (*massmatrix)(i             ,i + ActNumDof0) = density*lrefe_*crosssec_ / 6;
      (*massmatrix)(i + ActNumDof0,i             ) = density*lrefe_*crosssec_ / 6;
    }
  }
  
  return;
} // DRT::ELEMENTS::Truss3::bt_nlnstiffmass2


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

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_TRUSS3
