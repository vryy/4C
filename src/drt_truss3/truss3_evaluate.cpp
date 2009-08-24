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
  else if (action=="calc_stat_force_damp") act = Truss3::calc_stat_force_damp;
  else if (action=="calc_struct_ptcstiff") act = Truss3::calc_struct_ptcstiff;
  else
    {
      cout<<action<<endl;
      dserror("Unknown type of action for Truss3");
    }

  /*number of degrees of freedom actually assigned to the discretization by the first node
   *(allows connectin truss3 and beam3 directly)*/
  int ActNumDof0 = discretization.NumDof(Nodes()[0]);

  switch(act)
  {
    case Truss3::calc_struct_ptcstiff:
    {
      EvaluatePTC(params, elemat1,ActNumDof0);
    }
    break;
    //action type for evaluating statistical forces
    case Truss3::calc_stat_force_damp:
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


      //test whether current processor is row map owner of the element
      if(this->Owner() != discretization.Comm().MyPID()) return 0;

      // get element displacements (for use in shear flow fields)
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vector 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      //evaluation of statistical forces with the local statistical forces vector fstat
      Epetra_SerialDenseVector fstat(lm.size());
      EvaluateStatForceDamp(params,mydisp,fstat,elemat1,ActNumDof0);

      /*all the above evaluated forces are used in assembly of the column map force vector no matter if the current processor if
       * the row map owner of these DOF*/
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
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
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

/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate Statistical forces and viscous damping (public)                                       cyron 01/09|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Truss3::EvaluateStatForceDamp(ParameterList& params,
                                                vector<double> mydisp,
                                                Epetra_SerialDenseVector& elevec1,
                                                Epetra_SerialDenseMatrix& elemat1,
                                                int& ActNumDof0)
{
  // frictional coefficient per unit length (approximated by the one of an infinitely long staff)
  double zeta = 4 * PI * lrefe_ * params.get<double>("ETA",0.0);

  /*the following funcitons are assuming linear interpolation of thermal and drag forces between the nodes
   * comparable to the case stoch_order == 1 for the beam3 element*/

  //computing damping matrix (by background fluid of thermal bath)

  //diagonal entries
  for(int i= 0; i<3; i++)
  {
    //first node
    elemat1(i,i) += zeta/3.0;
    //second node
    elemat1(ActNumDof0+i,ActNumDof0+i) += zeta/3.0;
  }

  //offdiagonal entries
  for(int i= 0; i<3; i++)
  {
    elemat1(i,ActNumDof0+i) += zeta/6.0;
    elemat1(ActNumDof0+i,i) += zeta/6.0;
  }

  //computing statistical forces due to fluctuation-dissipation theorem

  // thermal energy responsible for statistical forces
  double kT = params.get<double>("KT",0.0);

  //calculating standard deviation of statistical forces according to fluctuation dissipation theorem
  double stand_dev_trans = pow(2 * kT * (zeta/6) / params.get<double>("delta time",0.01),0.5);

  //creating a random generator object which creates random numbers with mean = 0 and standard deviation
  //stand_dev; using Blitz namespace "ranlib" for random number generation
  ranlib::Normal<double> normalGen(0,stand_dev_trans);

  //uncorrelated part of the statistical forces
  for(int i= 0; i<3; i++)
  {
    //first node
    elevec1[i] += normalGen.random();
    //second node
    elevec1[ActNumDof0+i] += normalGen.random();
  }

  //correlated part of the statistical forces
  for(int i= 0; i<3; i++)
  {
    double force = normalGen.random();

    //first node
    elevec1[i] += force;
    //second node
    elevec1[ActNumDof0+i] += force;
  }


  return 0;
} //DRT::ELEMENTS::Truss3::EvaluateStatisticalNeumann



/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate PTC damping (public)                                                                  cyron 01/09|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Truss3::EvaluatePTC(ParameterList& params,
                                      Epetra_SerialDenseMatrix& elemat1,
                                      int& ActNumDof0)
{
  double dti = params.get<double>("dti",0.0);

  //first node
  for(int i= 0; i<3; i++)
    elemat1(i,i) += dti;

  //second node
  for(int i= 0; i<3; i++)
    elemat1(i+ActNumDof0,i+ActNumDof0) += dti;

  return 0;
} //DRT::ELEMENTS::Truss3::EvaluatePTC

/*--------------------------------------------------------------------------------------*
 | switch between kintypes                                                      tk 11/08|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_nlnstiffmass( vector<double>& disp,
    Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector* force,
    int& ActNumDof0)
{
  switch(kintype_)
  {
  case tr3_totlag:
    t3_nlnstiffmass_totlag(disp,stiffmatrix,massmatrix,force,ActNumDof0);
    return;
  case tr3_engstrain:
    t3_nlnstiffmass_engstr(disp,stiffmatrix,massmatrix,force,ActNumDof0);
    return;
  }
}


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 08/08|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_nlnstiffmass_totlag( vector<double>& disp,
    Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector* force,
    int& ActNumDof0)
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
    xcurr(j+3)   = Nodes()[1]->X()[j] + disp[ActNumDof0 + j]; //second node
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
      sm = ym / (2*(1 + actmat->PoissonRatio()));
      density = actmat->Density();
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

    //auxiliary variables for handling indices:
    int id = 0;
    int jd = 0;

    for (int i=0; i<6; ++i)
    {
      for (int j=0; j<6; ++j)
      {
        if(i<3)
          id = i;
        else
          id = i + ActNumDof0 - 3;
        if(j<3)
          jd = j;
        else
          jd = j + ActNumDof0 - 3;

        (*stiffmatrix)(id,jd) += (16*ym*crosssec_/pow(lrefe_,3))*aux(i)*aux(j);
      }
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
 | engineering strain measure, large displacements and rotations                                                |
  *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_nlnstiffmass_engstr( vector<double>& disp,
    Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector* force,
    int& ActNumDof0)
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
    xcurr(j+3)   = Nodes()[1]->X()[j] + disp[ActNumDof0 + j]; //second node
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
  {
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
        (*stiffmatrix)(i              ,i             )   =  forcescalar;
        (*stiffmatrix)(i              ,ActNumDof0 + i)   = -forcescalar;
        //stiffness entries for second node
        (*stiffmatrix)(i + ActNumDof0 ,i + ActNumDof0)   =  forcescalar;
        (*stiffmatrix)(i + ActNumDof0 ,i             )   = -forcescalar;
    }

    //auxiliary variables for handling indices:
    int id = 0;
    int jd = 0;

    for (int i=0; i<6; ++i)
    {
      for (int j=0; j<6; ++j)
      {
        if(i<3)
          id = i;
        else
          id = i + ActNumDof0 - 3;
        if(j<3)
          jd = j;
        else
          jd = j + ActNumDof0 - 3;

        (*stiffmatrix)(id,jd) += (ym*crosssec_/pow(lrefe_,3))*aux(i)*aux(j);
      }
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
