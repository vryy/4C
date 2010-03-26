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

#include "truss3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/linalg_fixedsizematrix.H"
//including random number library of blitz for statistical forces
#include <random/normal.h>
#include "../drt_mat/stvenantkirchhoff.H"

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
  else if (action=="calc_brownian")        act = Truss3::calc_brownian;
  else if (action=="calc_struct_ptcstiff") act = Truss3::calc_struct_ptcstiff;
  else
    {
      cout<<action<<endl;
      dserror("Unknown type of action for Truss3");
    }

  switch(act)
  {
    case Truss3::calc_struct_ptcstiff:
    {
      EvaluatePTC(params, elemat1);
    }
    break;
    //action type for evaluating statistical forces
    case Truss3::calc_brownian:
    {
      /*in order to understand the way how in the following parallel computing is handled one has to
       * be aware of what usually happens: each processor evaluates forces on each of its elements no
       * matter whether it's a real element or only a ghost element; later on in the assembly of the
       * evaluation method each processor adds the thereby gained DOF forces only for the DOF of
       * which it is the row map owner; if one used this way for random forces the following problem
       * would arise: if two processors evaluate both the same element (one time as a real element, one time
       * as a ghost element) and assemble later only the random forces of their own DOF the assembly for
       * different DOF of the same element might take place by means of random forces evaluated during different
       * element calls (one time the real element was called, one time only the ghost element); as a consequence
       * prescribing any correlation function between random forces of different DOF is in general
       * impossible since the random forces of two different DOF may have been evaluated during different
       * element calls; the only solution to this problem is obviously to allow element evaluation
       * to one processor only (the row map owner processor) and to allow this processor to assemble the
       * thereby gained random forces later on also to DOF of which it is not the row map owner; so fist
       * we check in this method whether the current processor is the owner of the element (if not
       * evaluation is cancelled) and finally this processor assembles the evaluated forces for degrees of
       * freedom right no matter whether its the row map owner of these DOF (outside the element routine
       * in the evaluation method of the discretization this would be not possible by default*/

      /*
      //test whether current processor is row map owner of the element
      if(this->Owner() != discretization.Comm().MyPID()) return 0;

      // get element displacements (for use in shear flow fields)
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vector 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      //evaluation of statistical forces with the local statistical forces vector fstat
      Epetra_SerialDenseVector fstat(lm.size());
      EvaluateStatForceDamp(params,mydisp,fstat,elemat1);

      //all the above evaluated forces are used in assembly of the column map force vector no matter if the current processor if
      //the row map owner of these DOF
      //note carefully: a space between the two subsequal ">" signs is mandatory for the C++ parser in order to avoid confusion with ">>" for streams
      RCP<Epetra_Vector>    fstatcol = params.get<  RCP<Epetra_Vector> >("statistical force vector",Teuchos::null);

      for(unsigned int i = 0; i < lm.size(); i++)
      {
        //note: lm contains the global Ids of the degrees of freedom of this element
        //testing whether the fstatcol vector has really an element related with the i-th element of fstat by the i-the entry of lm
        if (!(fstatcol->Map()).MyGID(lm[i])) dserror("Sparse vector fstatcol does not have global row %d",lm[i]);

        //get local Id of the fstatcol vector related with a certain element of fstat
        int lid = (fstatcol->Map()).LID(lm[i]);

        //add to the related element of fstatcol the contribution of fstat
        (*fstatcol)[lid] += fstat[i];
      }
      */

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

      // get element velocities
      RefCountPtr<const Epetra_Vector> vel  = discretization.GetState("velocity");
      if (vel==null) dserror("Cannot get state vectors 'velocity'");
      vector<double> myvel(lm.size());
      DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

      // for engineering strains instead of total lagrange use t3_nlnstiffmass2
      if (act == Truss3::calc_struct_nlnstiffmass)
      t3_nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
      else if (act == Truss3::calc_struct_nlnstifflmass)
      {
        t3_nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
        // lump mass matrix (bborn 07/08)
        // the mass matrix is lumped anyway, cf #b3_nlnstiffmass
        //b3_lumpmass(&elemat2);
      }
      else if (act == Truss3::calc_struct_nlnstiff)
      t3_nlnstiffmass(params,myvel,mydisp,&elemat1,NULL,&elevec1);
      else if (act == Truss3::calc_struct_internalforce)
      t3_nlnstiffmass(params,myvel,mydisp,NULL,NULL,&elevec1);

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
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
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
  const vector<int>* onoff = condition.Get<vector<int> >("onoff");
  // val is related to the 6 "val" fields after the onoff flags of the Neumann condition
  // in the input file; val gives the values of the force as a multiple of the prescribed load curve
  const vector<double>* val = condition.Get<vector<double> >("val");

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
  
    for (int dof=0; dof < 3; ++dof)
    {
      //computing entries for first node
      elevec1[dof] += funct[0] *ar[dof];
      //computing entries for first node
      elevec1[3 + dof] += funct[1] *ar[dof];
    }
      
  } // for (int ip=0; ip<intpoints.nquad; ++ip)

  return 0;
}


/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate PTC damping (public)                                                                  cyron 01/09|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Truss3::EvaluatePTC(ParameterList& params,
                                      Epetra_SerialDenseMatrix& elemat1)
{
  double dti = params.get<double>("dti",0.0);

  for(int i= 0; i<6; i++)
    elemat1(i,i) += dti;

  return 0;
} //DRT::ELEMENTS::Truss3::EvaluatePTC

/*--------------------------------------------------------------------------------------*
 | switch between kintypes                                                      tk 11/08|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_nlnstiffmass(ParameterList& params,
    vector<double>&           vel,
    vector<double>&           disp,
    Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector* force)
{
  
  /*first displacement vector is modified for proper element evaluation in case of periodic boundary conditions; in case that
   *no periodic boundary conditions are to be applied the following code line may be ignored or deleted*/
  NodeShift<2,3>(params,disp);
  
  switch(kintype_)
  {
  case tr3_totlag:
    t3_nlnstiffmass_totlag(disp,stiffmatrix,massmatrix,force);
    return;
  case tr3_engstrain:
    t3_nlnstiffmass_engstr(disp,stiffmatrix,massmatrix,force);
    return;
  }
  
  /*the following function call applies statistical forces and damping matrix according to the fluctuation dissipation theorem;
   * it is dedicated to the application of truss3 elements in the frame of statistical mechanics problems; for these problems a
   * special vector has to be passed to the element packed in the params parameter list; in case that the control routine calling
   * the element does not attach this special vector to params the following method is just doing nothing, which means that for
   * any ordinary problem of structural mechanics it may be ignored*/
   CalcBrownian<2,3,3,3>(params,vel,disp,stiffmatrix,force);
  
}


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 08/08|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_nlnstiffmass_totlag( vector<double>& disp,
    Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector* force)
{
  //current node position (first entries 0 .. 2 for first node, 3 ..5 for second node)
  LINALG::Matrix<6,1> xcurr;

  /*current nodal displacement (first entries 0 .. 2 for first node, 3 ..5 for second node) compared
   * to reference configuration; note: in general this is not equal to the values in disp since the
   * latter one referes to a nodal displacement compared to a reference configuration before the first
   * time step whereas the following variable referes to the displacement with respect to a reference
   * configuration which may have been set up at any point of time during the simulation (usually this
   * is only important if an element has been added to the discretization after the start of the simulation)*/
  LINALG::Matrix<6,1> ucurr;

  //Green-Lagrange strain
  double epsilon;

  //auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  LINALG::Matrix<6,1> aux;

  //current nodal position
  for (int j=0; j<3; ++j)
  {
    xcurr(j  )   = Nodes()[0]->X()[j] + disp[  j]; //first node
    xcurr(j+3)   = Nodes()[1]->X()[j] + disp[3+j]; //second node
  }

  //current displacement = current position - reference position
  ucurr  = xcurr;
  ucurr -= X_;

  //computing auxiliary vector aux = N^T_{,xi} * N_{,xi} * xcurr
  aux(0) = 0.25 * (xcurr(0) - xcurr(3));
  aux(1) = 0.25 * (xcurr(1) - xcurr(4));
  aux(2) = 0.25 * (xcurr(2) - xcurr(5));
  aux(3) = 0.25 * (xcurr(3) - xcurr(0));
  aux(4) = 0.25 * (xcurr(4) - xcurr(1));
  aux(5) = 0.25 * (xcurr(5) - xcurr(2));

  //calculating strain epsilon from node position by scalar product:
  //epsilon = (xrefe + 0.5*ucurr)^T * N_{,s}^T * N_{,s} * d
  epsilon = 0;
  epsilon += (X_(0) + 0.5*ucurr(0)) * (ucurr(0) - ucurr(3));
  epsilon += (X_(1) + 0.5*ucurr(1)) * (ucurr(1) - ucurr(4));
  epsilon += (X_(2) + 0.5*ucurr(2)) * (ucurr(2) - ucurr(5));
  epsilon += (X_(3) + 0.5*ucurr(3)) * (ucurr(3) - ucurr(0));
  epsilon += (X_(4) + 0.5*ucurr(4)) * (ucurr(4) - ucurr(1));
  epsilon += (X_(5) + 0.5*ucurr(5)) * (ucurr(5) - ucurr(2));
  epsilon /= lrefe_*lrefe_;


  /* read material parameters using structure _MATERIAL which is defined by inclusion of      /
   / "../drt_lib/drt_timecurve.H"; note: material parameters have to be read in the evaluation /
   / function instead of e.g. Truss3_input.cpp or within the Truss3Register class since it is not/
   / sure that structure _MATERIAL is declared within those scopes properly whereas it is within/
   / the evaluation functions */

  // get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double ym = 0;
  double sm = 0;
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
    for (int i=0; i<6; ++i)
     (*force)(i) = (4*ym*crosssec_*epsilon/lrefe_) * aux(i);
  }


  //computing linear stiffness matrix
  if (stiffmatrix != NULL)
  {
    for (int i=0; i<3; ++i)
    {
        //stiffness entries for first node
        (*stiffmatrix)(i,i)     =  (ym*crosssec_*epsilon/lrefe_);
        (*stiffmatrix)(i,3+i)   = -(ym*crosssec_*epsilon/lrefe_);
        //stiffness entries for second node
        (*stiffmatrix)(i+3,i+3) =  (ym*crosssec_*epsilon/lrefe_);
        (*stiffmatrix)(i+3,i )  = -(ym*crosssec_*epsilon/lrefe_);
    }

    for (int i=0; i<6; ++i)
      for (int j=0; j<6; ++j)
        (*stiffmatrix)(i,j) += (16*ym*crosssec_/pow(lrefe_,3))*aux(i)*aux(j);
   }

  //calculating consistent mass matrix
  if (massmatrix != NULL)
  {
    for (int i=0; i<3; ++i)
    {
      (*massmatrix)(i,i) = density*lrefe_*crosssec_ / 3;
      (*massmatrix)(i+3,i+3) = density*lrefe_*crosssec_ / 3;
      (*massmatrix)(i,i+3) = density*lrefe_*crosssec_ / 6;
      (*massmatrix)(i+3,i) = density*lrefe_*crosssec_ / 6;
    }
  }

  return;
} // DRT::ELEMENTS::Truss3::t3_nlnstiffmass


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                      tk 10/08|
 | engineering strain measure, large displacements and rotations                                                |
  *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_nlnstiffmass_engstr( vector<double>& disp,
    Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector* force)
{
  //current node position (first entries 0 .. 2 for first node, 3 ..5 for second node)
  LINALG::Matrix<6,1> xcurr;

  //Green-Lagrange strain
  double epsilon;

  //auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  LINALG::Matrix<6,1> aux;

  //current nodal position (first
  for (int j=0; j<3; ++j)
  {
    xcurr(j  )   = Nodes()[0]->X()[j] + disp[  j]; //first node
    xcurr(j+3)   = Nodes()[1]->X()[j] + disp[3+j]; //second node
  }

  //computing auxiliary vector aux = 4.0*N^T_{,xi} * N_{,xi} * xcurr
  aux(0) = (xcurr(0) - xcurr(3));
  aux(1) = (xcurr(1) - xcurr(4));
  aux(2) = (xcurr(2) - xcurr(5));
  aux(3) = (xcurr(3) - xcurr(0));
  aux(4) = (xcurr(4) - xcurr(1));
  aux(5) = (xcurr(5) - xcurr(2));

  double lcurr = sqrt(pow(aux(0),2)+pow(aux(1),2)+pow(aux(2),2));

  //calculating strain epsilon from node position by scalar product:
  epsilon = (lcurr-lrefe_)/lrefe_;

  /* read material parameters using structure _MATERIAL which is defined by inclusion of      /
   / "../drt_lib/drt_timecurve.H"; note: material parameters have to be read in the evaluation /
   / function instead of e.g. Truss3_input.cpp or within the Truss3Register class since it is not/
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
    for (int i=0; i<6; ++i)
     (*force)(i) = forcescalar * aux(i);


  //computing linear stiffness matrix
  if (stiffmatrix != NULL)
  {
    for (int i=0; i<3; ++i)
    {
        //stiffness entries for first node
        (*stiffmatrix)(i,i)    =  forcescalar;
        (*stiffmatrix)(i,3+i)  = -forcescalar;
        //stiffness entries for second node
        (*stiffmatrix)(i+3,i+3)=  forcescalar;
        (*stiffmatrix)(i+3,i)  = -forcescalar;
    }

    for (int i=0; i<6; ++i)
      for (int j=0; j<6; ++j)
        (*stiffmatrix)(i,j) += (ym*crosssec_/pow(lrefe_,3))*aux(i)*aux(j);
  }

  //calculating consistent mass matrix
  if (massmatrix != NULL)
  {
    for (int i=0; i<3; ++i)
    {
      (*massmatrix)(i,i)     = density*lrefe_*crosssec_ / 3;
      (*massmatrix)(i+3,i+3) = density*lrefe_*crosssec_ / 3;
      (*massmatrix)(i,i+3)   = density*lrefe_*crosssec_ / 6;
      (*massmatrix)(i+3,i)   = density*lrefe_*crosssec_ / 6;
    }
  }
  
  return;
} // DRT::ELEMENTS::Truss3::bt_nlnstiffmass3


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

/*-----------------------------------------------------------------------------------------------------------*
 | computes damping coefficients per lengthand stores them in a matrix in the following order: damping of    |
 | translation parallel to filament axis, damping of translation orthogonal to filament axis, damping of     |
 | rotation around filament axis                                             (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Truss3::MyDampingConstants(ParameterList& params,LINALG::Matrix<3,1>& gamma, const INPAR::STATMECH::FrictionModel& frictionmodel)
{  
  //translational damping coefficients according to Howard, p. 107, table 6.2;
  gamma(0) = 2*PI*params.get<double>("ETA",0.0);
  gamma(1) = 4*PI*params.get<double>("ETA",0.0);
  //no rotational damping as no rotaional degrees of freedom
  gamma(2) = 0;
  

  //in case of an isotropic friction model the same damping coefficients are applied parallel to the polymer axis as perpendicular to it
  if(frictionmodel == INPAR::STATMECH::frictionmodel_isotropicconsistent || frictionmodel == INPAR::STATMECH::frictionmodel_isotropiclumped)
    gamma(0) = gamma(1);
 
}//DRT::ELEMENTS::Truss3::MyDampingConstants

/*-----------------------------------------------------------------------------------------------------------*
 |computes the number of different random numbers required in each time step for generation of stochastic    |
 |forces;                                                                    (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Truss3::HowManyRandomNumbersINeed()
{
  /*at each Gauss point one needs as many random numbers as randomly excited degrees of freedom, i.e. three
   *random numbers for the translational degrees of freedom*/
  return (3*2);

}//DRT::ELEMENTS::Beam3::HowManyRandomNumbersINeed

/*-----------------------------------------------------------------------------------------------------------*
 |computes velocity of background fluid and gradient of that velocity at a certain evaluation point in       |
 |the physical space                                                         (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int ndim> //number of dimensions of embedding space
void DRT::ELEMENTS::Truss3::MyBackgroundVelocity(ParameterList& params,  //!<parameter list
                                                const LINALG::Matrix<ndim,1>& evaluationpoint,  //!<point at which background velocity and its gradient has to be computed
                                                LINALG::Matrix<ndim,1>& velbackground,  //!< velocity of background fluid
                                                LINALG::Matrix<ndim,ndim>& velbackgroundgrad) //!<gradient of velocity of background fluid
{
  
  /*note: this function is not yet a general one, but always assumes a shear flow, where the velocity of the
   * background fluid is always directed in x-direction. In 3D the velocity increases linearly in z and equals zero for z = 0.
   * In 2D the velocity increases linearly in y and equals zero for y = 0. */
  
  velbackground.PutScalar(0);
  velbackground(0) = evaluationpoint(ndim-1) * params.get<double>("CURRENTSHEAR",0.0);
  
  velbackgroundgrad.PutScalar(0);
  velbackgroundgrad(0,ndim-1) = params.get<double>("CURRENTSHEAR",0.0);

}

/*-----------------------------------------------------------------------------------------------------------*
 | computes translational damping forces and stiffness (public)                                 cyron   03/10|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node
inline void DRT::ELEMENTS::Truss3::MyTranslationalDamping(ParameterList& params,  //!<parameter list
                                                  const vector<double>&     vel,  //!< element velocity vector
                                                  const vector<double>&     disp, //!<element disp vector
                                                  Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                                  Epetra_SerialDenseVector* force)//!< element internal force vector
{  
  //get time step size
  double dt = params.get<double>("delta time",0.0);
  
  //velocity and gradient of background velocity field
  LINALG::Matrix<ndim,1> velbackground;
  LINALG::Matrix<ndim,ndim> velbackgroundgrad;
  
  //evaluation point in physical space corresponding to a certain Gauss point in parameter space
  LINALG::Matrix<ndim,1> evaluationpoint;
  
  //get friction model according to which forces and damping are applied
  INPAR::STATMECH::FrictionModel frictionmodel = Teuchos::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");

  //damping coefficients for translational and rotatinal degrees of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel);
  
  //get vector jacobi with Jacobi determinants at each integration point (gets by default those values required for consistent damping matrix)
  vector<double> jacobi(jacobimass_);
  
  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  IntegrationType integrationtype = gaussexactintegration;
  if(frictionmodel == INPAR::STATMECH::frictionmodel_isotropiclumped)
  {
    integrationtype = lobattointegration;
    jacobi = jacobinode_;
  }
  
  //get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode,integrationtype));
  
  //matrix to store basis functions and their derivatives evaluated at a certain Gauss point
  LINALG::Matrix<1,nnode> funct;
  LINALG::Matrix<1,nnode> deriv;

  for(int gp=0; gp < gausspoints.nquad; gp++)
  {    
    //evaluate basis functions and their derivatives at current Gauss point
    DRT::UTILS::shape_function_1D(funct,gausspoints.qxg[gp][0],Shape());
    DRT::UTILS::shape_function_1D_deriv1(deriv,gausspoints.qxg[gp][0],Shape());
     
    //compute point in phyiscal space corresponding to Gauss point
    evaluationpoint.PutScalar(0);
    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop over all dimensions
      for(int j=0; j<ndim; j++)
        evaluationpoint(j) += funct(i)*(Nodes()[i]->X()[j]+disp[dof*i+j]);
    
    //compute velocity and gradient of background flow field at evaluationpoint
    MyBackgroundVelocity<ndim>(params,evaluationpoint,velbackground,velbackgroundgrad);

 
    //compute tangent vector t_{\par} at current Gauss point
    LINALG::Matrix<ndim,1> tpar(true);
    for(int i=0; i<nnode; i++)
      for(int k=0; k<ndim; k++)
        tpar(k) += deriv(i)*(Nodes()[i]->X()[k]+disp[dof*i+k]) / jacobi[gp];
    
    //compute velocity vector at this Gauss point
    LINALG::Matrix<ndim,1> velgp(true);
    for(int i=0; i<nnode; i++)
      for(int l=0; l<ndim; l++)
        velgp(l) += funct(i)*vel[dof*i+l]; 
    
    //compute matrix product (t_{\par} \otimes t_{\par}) \cdot velbackgroundgrad
    LINALG::Matrix<ndim,ndim> tpartparvelbackgroundgrad(true);
    for(int i=0; i<ndim; i++)
      for(int j=0; j<ndim; j++)
        for(int k=0; k<ndim; k++)
          tpartparvelbackgroundgrad(i,j) += tpar(i)*tpar(k)*velbackgroundgrad(k,j);
        
    //loop over all line nodes
    for(int i=0; i<nnode; i++)            
      //loop over lines of matrix t_{\par} \otimes t_{\par}
      for(int k=0; k<ndim; k++)
        //loop over columns of matrix t_{\par} \otimes t_{\par}
        for(int l=0; l<ndim; l++)           
        {               
          if(force != NULL)
            (*force)(i*dof+k)+= funct(i)*jacobi[gp]*gausspoints.qwgt[gp]*( (k==l)*gamma(1) + (gamma(0) - gamma(1))*tpar(k)*tpar(l) ) *(velgp(l)- velbackground(l));
          
          if(stiffmatrix != NULL)
            //loop over all column nodes
            for (int j=0; j<nnode; j++) 
            {
              (*stiffmatrix)(i*dof+k,j*dof+l) += gausspoints.qwgt[gp]*funct(i)*funct(j)*jacobi[gp]*(                 (k==l)*gamma(1) + (gamma(0) - gamma(1))*tpar(k)*tpar(l) ) / dt;
              (*stiffmatrix)(i*dof+k,j*dof+l) -= gausspoints.qwgt[gp]*funct(i)*funct(j)*jacobi[gp]*( velbackgroundgrad(k,l)*gamma(1) + (gamma(0) - gamma(1))*tpartparvelbackgroundgrad(k,l) ) ;             
              (*stiffmatrix)(i*dof+k,j*dof+k) += gausspoints.qwgt[gp]*funct(i)*deriv(j)*                                                   (gamma(0) - gamma(1))*tpar(l)*(velgp(l) - velbackground(l));
              (*stiffmatrix)(i*dof+k,j*dof+l) += gausspoints.qwgt[gp]*funct(i)*deriv(j)*                                                   (gamma(0) - gamma(1))*tpar(k)*(velgp(l) - velbackground(l));
            }    
        }   
  }
 
  return;
}//DRT::ELEMENTS::Truss3::MyTranslationalDamping(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes stochastic forces and resulting stiffness (public)                                  cyron   03/10|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::Truss3::MyStochasticForces(ParameterList& params,  //!<parameter list
                                              const vector<double>&     vel,  //!< element velocity vector
                                              const vector<double>&     disp, //!<element disp vector
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector* force)//!< element internal force vector
{
  //get friction model according to which forces and damping are applied
  INPAR::STATMECH::FrictionModel frictionmodel = Teuchos::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");
  
  //damping coefficients for three translational and one rotatinal degree of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel);
  

  //get vector jacobi with Jacobi determinants at each integration point (gets by default those values required for consistent damping matrix)
  vector<double> jacobi(jacobimass_);
  
  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  IntegrationType integrationtype = gaussexactintegration;
  if(frictionmodel == INPAR::STATMECH::frictionmodel_isotropiclumped)
  {
    integrationtype = lobattointegration;
    jacobi = jacobinode_;
  }
  
  //get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode,integrationtype));
  
  //matrix to store basis functions and their derivatives evaluated at a certain Gauss point
  LINALG::Matrix<1,nnode> funct;
  LINALG::Matrix<1,nnode> deriv;
  
  
  /*get pointer at Epetra multivector in parameter list linking to random numbers for stochastic forces with zero mean
   * and standard deviation (2*kT / dt)^0.5; note carefully: a space between the two subsequal ">" signs is mandatory
   * for the C++ parser in order to avoid confusion with ">>" for streams*/
   RCP<Epetra_MultiVector> randomnumbers = params.get<  RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null);
   


  for(int gp=0; gp < gausspoints.nquad; gp++)
  {
    //evaluate basis functions and their derivatives at current Gauss point
    DRT::UTILS::shape_function_1D(funct,gausspoints.qxg[gp][0],Shape());
    DRT::UTILS::shape_function_1D_deriv1(deriv,gausspoints.qxg[gp][0],Shape());
    
    //compute tangent vector t_{\par} at current Gauss point
    LINALG::Matrix<ndim,1> tpar(true);
    for(int i=0; i<nnode; i++)
      for(int k=0; k<ndim; k++)
        tpar(k) += deriv(i)*(Nodes()[i]->X()[k]+disp[dof*i+k]) / jacobi[gp];
     
    
    //loop over all line nodes
    for(int i=0; i<nnode; i++)             
      //loop dimensions with respect to lines
      for(int k=0; k<ndim; k++)
        //loop dimensions with respect to columns
        for(int l=0; l<ndim; l++)           
        {
          if(force != NULL)
            (*force)(i*dof+k) -= funct(i)*(sqrt(gamma(1))*(k==l) + (sqrt(gamma(0)) - sqrt(gamma(1)))*tpar(k)*tpar(l))*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(jacobi[gp]*gausspoints.qwgt[gp]);          

          if(stiffmatrix != NULL)
            //loop over all column nodes
            for (int j=0; j<nnode; j++) 
            {            
              (*stiffmatrix)(i*dof+k,j*dof+k) -= funct(i)*deriv(j)*tpar(l)*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(gausspoints.qwgt[gp]/ jacobi[gp])*(sqrt(gamma(0)) - sqrt(gamma(1)));   
              (*stiffmatrix)(i*dof+k,j*dof+l) -= funct(i)*deriv(j)*tpar(k)*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(gausspoints.qwgt[gp]/ jacobi[gp])*(sqrt(gamma(0)) - sqrt(gamma(1)));  
            }
        }  
  }

  return;
}//DRT::ELEMENTS::Truss3::MyStochasticForces(.)


/*-----------------------------------------------------------------------------------------------------------*
 | Assemble stochastic and viscous forces and respective stiffness according to fluctuation dissipation      |
 | theorem                                                                               (public) cyron 03/10|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::Truss3::CalcBrownian(ParameterList& params,
                                              const vector<double>&           vel,  //!< element velocity vector
                                              const vector<double>&           disp, //!< element displacement vector
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector* force) //!< element internal force vector
{   
  //if no random numbers for generation of stochastic forces are passed to the element no Brownian dynamics calculations are conducted
  if( params.get<  RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) == Teuchos::null)
    return;
  
  //add stiffness and forces due to translational damping effects
  MyTranslationalDamping<nnode,ndim,dof>(params,vel,disp,stiffmatrix,force); 

  //add stochastic forces and (if required) resulting stiffness
  MyStochasticForces<nnode,ndim,dof,randompergauss>(params,vel,disp,stiffmatrix,force);

return;

}//DRT::ELEMENTS::Truss3::CalcBrownian(.)

/*-----------------------------------------------------------------------------------------------------------*
 | shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions; if two   |
 | nodes within one element are separated by a periodic boundary, one of them is shifted such that the final |
 | distance in R^3 is the same as the initial distance in the periodic space; the shift affects computation  |
 | on element level within that very iteration step, only (no change in global variables performed)          |                                 |
 |                                                                                       (public) cyron 10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim> //number of nodes, number of dimensions
inline void DRT::ELEMENTS::Truss3::NodeShift(ParameterList& params,  //!<parameter list
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
      for(int dof=0; dof<ndim; dof++)
      {   
        /*if the distance in some coordinate direction between some node and the first node becomes smaller by adding or subtracting
         * the period length, the respective node has obviously been shifted due to periodic boundary conditions and should be shifted
         * back for evaluation of element matrices and vectors; this way of detecting shifted nodes works as long as the element length
         * is smaller than half the periodic length*/
        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) + params.get<double>("PeriodLength",0.0) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
          disp[numdof*i+dof] += params.get<double>("PeriodLength",0.0);
          
        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - params.get<double>("PeriodLength",0.0) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
          disp[numdof*i+dof] -= params.get<double>("PeriodLength",0.0);
      }
    }
return;

}//DRT::ELEMENTS::Truss3::NodeShift

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_TRUSS3
