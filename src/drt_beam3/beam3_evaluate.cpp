/*!-----------------------------------------------------------------------------------------------------------
 \file beam3_evaluate.cpp
 \brief

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

 *-----------------------------------------------------------------------------------------------------------*/
#ifdef D_BEAM3
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "beam3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

//externally defined structure for material data
extern struct _MATERIAL *mat;

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3::Evaluate(ParameterList& params,
    DRT::Discretization& discretization,
    vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Beam3::ActionType act = Beam3::calc_none;
  // get the action required
  string action = params.get<string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff") act = Beam3::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff") act = Beam3::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Beam3::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass") act = Beam3::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass") act = Beam3::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = Beam3::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress") act = Beam3::calc_struct_stress;
  else if (action=="calc_struct_eleload") act = Beam3::calc_struct_eleload;
  else if (action=="calc_struct_fsiload") act = Beam3::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep") act = Beam3::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike") act = Beam3::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep") act = Beam3::calc_struct_reset_istep;
  else dserror("Unknown type of action for Beam3");

  switch(act)
  {
    /*in case that only linear stiffness matrix is required b3_nlstiffmass is called with zero dispalcement and 
     residual values*/
    case Beam3::calc_struct_linstiff:
    {
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");

    }
    break;

    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case Beam3::calc_struct_nlnstiffmass:
    case Beam3::calc_struct_nlnstifflmass:
    case Beam3::calc_struct_nlnstiff:
    case Beam3::calc_struct_internalforce:
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

      if (act == Beam3::calc_struct_nlnstiffmass)
      b3_nlnstiffmass(mydisp,&elemat1,&elemat2,&elevec1);
      else if (act == Beam3::calc_struct_nlnstifflmass)
      {
        b3_nlnstiffmass(mydisp,&elemat1,&elemat2,&elevec1);
        // lump mass matrix (bborn 07/08)
        // the mass matrix is lumped anyway, cf #b3_nlnstiffmass
        //b3_lumpmass(&elemat2);
      }
      else if (act == Beam3::calc_struct_nlnstiff)
      b3_nlnstiffmass(mydisp,&elemat1,NULL,&elevec1);
      else if (act == Beam3::calc_struct_internalforce)
      b3_nlnstiffmass(mydisp,NULL,NULL,&elevec1);
    
    }
    break;
    case calc_struct_update_istep:
    case calc_struct_update_imrlike:
    {
      /*the action calc_struct_update_istep is called in the very end of a time step when the new dynamic
       * equilibrium has finally been found; this is the point where the variable representing the geomatric
       * status of the beam have to be updated; the geometric status is represented by means of the triad Tnew_,
       * the curvature curvnew_ and the angular values betaplusalphanew_ and betaminusalphanew_*/
      Told_ = Tnew_;
      curvold_ = curvnew_;
      betaplusalphaold_ = betaplusalphanew_;
      betaminusalphaold_ = betaminusalphanew_; 
    }
    break;
    case calc_struct_reset_istep:
    {
      /*the action calc_struct_reset_istep is called by the adaptive time step controller; carries out one test
       * step whose purpose is only figuring out a suitabel timestep; thus this step may be a very bad one in order
       * to iterated towards the new dynamic equilibrium and the thereby gained new geometric configuration should 
       * not be applied for any further iteration step; as a consequence the thereby generated change of the geometric
       * configuration should be canceled and the configuration should be reset to the value at the beginning of the 
       * time step*/
      Tnew_ = Told_;
      curvnew_ = curvold_;
      betaplusalphanew_ = betaplusalphaold_;
      betaminusalphanew_ = betaminusalphaold_; 
    }
    break;
    default:
    dserror("Unknown type of action for Beam3 %d", act);
  }
  return 0;

}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)                                       cyron 03/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Beam3::EvaluateNeumann(ParameterList& params,
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

  //stochastic forces due to fluctuation-dissipation theorem
  if (thermalenergy_ > 0)
  {
    extern struct _MATERIAL *mat;
    // get the material law and density
    MATERIAL* currmat = &(mat[material_-1]);
    double density;

    //assignment of material parameters; only St.Venant material is accepted for this beam 
    switch(currmat->mattyp)
    {
      case m_stvenant:// only linear elastic material supported
      {
        density = currmat->m.stvenant->density;
      }
      break;
      default:
      dserror("unknown or improper type of material law");
    }
        
    //calculating diagonal entry of damping matrix  
    double gammatrans = params.get<double>("damping factor M",0.0) * crosssec_ * density * lrefe_ / 2;
    double gammarot   = params.get<double>("damping factor M",0.0) * Iyy_      * density * lrefe_ / 2;
    
    
    //calculating standard deviation of statistical forces according to fluctuation dissipation theorem
    double standdevtrans = pow(2 * thermalenergy_ * gammatrans   / params.get<double>("delta time",0.01),0.5);
    double standdevrot   = pow(2 * thermalenergy_ * gammarot     / params.get<double>("delta time",0.01),0.5);
    
    //creating random generator objects which create random numbers with mean = 0 and standard deviation
    //standdevtrans and standdevrot; using Blitz namespace "ranlib" for random number generation
    ranlib::Normal<double> normalGenTrans(0,standdevtrans);
    ranlib::Normal<double> normalGenRot(0,standdevrot);
    
    //adding statistical forces 
    for (int i=0; i<3; ++i)
    {
      elevec1(i)   += normalGenTrans.random();
      //elevec1(i+3) += normalGenRot.random();
      elevec1(i+6) += normalGenTrans.random();
      //elevec1(i+9) += normalGenRot.random();
    }  
  
  }
  
  return 0;
}

/*-----------------------------------------------------------------------------------------------------------*
 | auxiliary functions for dealing with large rotations and nonlinear stiffness                    cyron 04/08|							     
 *----------------------------------------------------------------------------------------------------------*/
//computing basis of stiffness matrix of Crisfield, Vol. 2, equation (17.81)
inline void DRT::ELEMENTS::Beam3::computestiffbasis(const BlitzMat3x3& Tnew_, const BlitzVec3& Cm, const BlitzVec3& Cb, const BlitzMat3x3& spinx21, Epetra_SerialDenseMatrix& stiffmatrix)
{
  //calculating the first matrix of (17.81) directly involves multiplications of large matrices (e.g. with the 12x6-matrix X)
  //application of the definitions in (17.74) allows blockwise evaluation with multiplication and addition of 3x3-matrices only
  //in the follwoing all the blocks are treated separately; their name is related directly to their content:
  //e.g. TCmTt is the product of the 3 matrices T * C_m * T^t (with T and C_m according to (17.74) and (17.76)
  //for the blockwise calculation on which the following steps are based on the relation S^t = -S for spin matrices was applied
  
  BlitzMat3x3 TCmTt;
  BlitzMat3x3 TCbTt;
  BlitzMat3x3 STCmTt;
  BlitzMat3x3 STCmTtSt;
  
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      TCmTt(i,j) = 0.0;
      TCbTt(i,j) = 0.0;
      for (int k = 0; k < 3; ++k)
      {
        TCmTt(i,j) += Tnew_(i,k)*Cm(k)*Tnew_(j,k);
        TCbTt(i,j) += Tnew_(i,k)*Cb(k)*Tnew_(j,k);
      }
    }
  }
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      STCmTt(i,j) = 0.0;
      for (int k = 0; k < 3; ++k)
      {
        STCmTt(i,j) += spinx21(i,k)*TCmTt(k,j);
      }
    }
  }
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      STCmTtSt(i,j) = 0.0;
      for (int k = 0; k < 3; ++k)
      {
        STCmTtSt(i,j) += STCmTt(i,k)*spinx21(j,k);
      }
    }
  }
  //calculating basis of stiffness matrix by means of above blocks
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      stiffmatrix(i  ,j)   =  TCmTt(i,j);
      stiffmatrix(i  ,j+3) =  STCmTt(j,i);
      stiffmatrix(i  ,j+6) = -TCmTt(i,j);
      stiffmatrix(i  ,j+9) =  STCmTt(j,i);
      stiffmatrix(i+3,j)   =  STCmTt(i,j);
      stiffmatrix(i+3,j+3) =  STCmTtSt(i,j) + TCbTt(i,j);
      stiffmatrix(i+3,j+6) = -STCmTt(i,j);
      stiffmatrix(i+3,j+9) =  STCmTtSt(j,i) - TCbTt(i,j);
      stiffmatrix(i+6,j)   = -TCmTt(i,j);
      stiffmatrix(i+6,j+3) = -STCmTt(j,i);
      stiffmatrix(i+6,j+6) =  TCmTt(i,j);
      stiffmatrix(i+6,j+9) = -STCmTt(j,i);
      stiffmatrix(i+9,j)   =  STCmTt(i,j);
      stiffmatrix(i+9,j+3) =  STCmTtSt(i,j)-TCbTt(i,j);
      stiffmatrix(i+9,j+6) = -STCmTt(i,j);
      stiffmatrix(i+9,j+9) =  STCmTtSt(i,j)+TCbTt(i,j);
    }
  }
  
  //Checking variation of strains with respect to displacements
  
  return;
} // DRT::ELEMENTS::Beam3::computestiffbasis

//computing spin matrix out of a rotation vector
inline void DRT::ELEMENTS::Beam3::computespin(BlitzMat3x3& spin, BlitzVec3 rotationangle, const double& spinscale)
{
  XFEM::BLITZTINY::V_scale<3>(rotationangle,spinscale);
  spin(0,0) = 0;
  spin(0,1) = -rotationangle(2);
  spin(0,2) = rotationangle(1);
  spin(1,0) = rotationangle(2);
  spin(1,1) = 0;
  spin(1,2) = -rotationangle(0);
  spin(2,0) = -rotationangle(1);
  spin(2,1) = rotationangle(0);
  spin(2,2) = 0;
  return;
} // DRT::ELEMENTS::Beam3::computespin 

//computing rotation matrix out of a spin matrix
inline void DRT::ELEMENTS::Beam3::computerotation(BlitzMat3x3& rotationmatrix, const BlitzMat3x3& spin)
{
  rotationmatrix = spin;
  rotationmatrix(0,0) +=1;
  rotationmatrix(1,1) +=1;
  rotationmatrix(2,2) +=1;
  
  return;
} //DRT::ELEMENTS::Beam3::computerotation 

//computing stiffens matrix Ksigma1 according to Crisfield, Vol. 2, equation (17.83)
inline void DRT::ELEMENTS::Beam3::computeKsig1(Epetra_SerialDenseMatrix& Ksig1, const BlitzVec3& stressn, const BlitzVec3& stressm)
{
  Ksig1.Shape(12,12);
  BlitzMat3x3 Sn;
  BlitzMat3x3 Sm;
  computespin(Sn,stressn, 0.5);
  computespin(Sm,stressm, 0.5);

  for (int i=0; i<3; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      Ksig1(i ,j+3) = Sn(i,j);
      Ksig1(i ,j+9) = Sn(i,j);
      Ksig1(i+6,j+3) = -Sn(i,j);
      Ksig1(i+6,j+9) = -Sn(i,j);

      Ksig1(i+3,j+3) = Sm(i,j);
      Ksig1(i+3,j+9) = Sm(i,j);
      Ksig1(i+9,j+3) = -Sm(i,j);
      Ksig1(i+9,j+9) = -Sm(i,j);
    }
  }
  return;
} //DRT::ELEMENTS::Beam3::computeKsig1

//computing stiffens matrix Ksigma1 according to Crisfield, Vol. 2, equation (17.87) and (17.88)
inline void DRT::ELEMENTS::Beam3::computeKsig2(Epetra_SerialDenseMatrix& Ksig2, const BlitzVec3& stressn, const BlitzVec3& x21)
{
  Ksig2.Shape(12,12);
  BlitzMat3x3 Sn;
  BlitzMat3x3 Sx21;
  BlitzMat3x3 Y;
  
  computespin(Sn,stressn, 0.5); 
  computespin(Sx21,x21, 0.25);
  XFEM::BLITZTINY::MM_product<3,3,3>(Sx21,Sn,Y);

  for (int i=0; i<3; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      Ksig2(i+3,j ) = -Sn(i,j);
      Ksig2(i+3,j+6) = Sn(i,j);
      Ksig2(i+9,j ) = -Sn(i,j);
      Ksig2(i+9,j+6) = Sn(i,j);

      Ksig2(i+3,j+3) = Y(i,j);
      Ksig2(i+3,j+9) = Y(i,j);
      Ksig2(i+9,j+3) = Y(i,j);
      Ksig2(i+9,j+9) = Y(i,j);
    }
  }
  return;
} // DRT::ELEMENTS::Beam3::computeKsig2


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 01/08|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3::b3_nlnstiffmass( vector<double>& disp,
    Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector* force)
{
  //constitutive laws from Crisfield, Vol. 2, equation (17.76)
  BlitzVec3 Cm;
  BlitzVec3 Cb;
  
  //normal/shear strain and bending strain(curvature)
  BlitzVec3 epsilonn;
  BlitzVec3 epsilonm;
  
  //stress values n and m, Crisfield, Vol. 2, equation (17.78)
  BlitzVec3 stressn;
  BlitzVec3 stressm;
  
  //current position of nodal degrees of freedom
  BlitzMat6x2 xcurr;
  
  //difference between coordinates of both nodes in current configuration, x21' Crisfield  Vol. 2 equ. (17.66a) and (17.72)
  BlitzVec3 x21;
  
  //nonlinear parts of stiffness matrix, Crisfiel Vol. 2, equation (17.83) and (17.87)
  Epetra_SerialDenseMatrix Ksig1;
  Epetra_SerialDenseMatrix Ksig2;
  
  //auxiliary variables
  BlitzVec3 deltabetaplusalpha;
  BlitzVec3 deltabetaminusalpha;
  BlitzMat3x3 Saux;
  BlitzMat3x3 Raux;

  //nodal coordinates in current position
  for (int k=0; k<2; ++k) //looping over number of nodes
  {
    for (int j=0; j<3; ++j) 
    {
      xcurr(j,k)   = Nodes()[k]->X()[j] + disp[k*6+j]; //translational DOF
      xcurr(j+3,k) = disp[k*6+j+3]; //rotational DOF
    }
  }

  //first of all "new" variables have to be adopted to dispalcement passed in from BACI driver

  //difference between coordinates of both nodes, x21' Crisfield  Vol. 2 equ. (17.66a) and (17.72)
  for (int j=0; j<3; ++j)
  {
    x21(j) = xcurr(j,1) - xcurr(j,0);
    betaplusalphanew_(j) = xcurr(j+3,1) + xcurr(j+3,0);
    betaminusalphanew_(j) = xcurr(j+3,1) - xcurr(j+3,0);
  }

  deltabetaplusalpha  = betaplusalphanew_;
  deltabetaplusalpha -= betaplusalphaold_;
  
  deltabetaminusalpha  = betaminusalphanew_;
  deltabetaminusalpha -= betaminusalphaold_;

  //auxiliary matrix for update of Tnew_, Crisfield, Vol 2, equation (17.65)
  computespin(Saux,deltabetaplusalpha, 0.5);
  computerotation(Raux,Saux);
  XFEM::BLITZTINY::MM_product<3,3,3>(Raux,Told_,Tnew_);

  //computing triad for curvature update, Crisfield, Vol 2, equation (17.73)
  computespin(Saux,deltabetaplusalpha, 0.25);
  computerotation(Raux,Saux);
  XFEM::BLITZTINY::MM_product<3,3,3>(Raux,Told_,Tmid_);

  //updating curvature, Crisfield, Vol. 2, equation (17.72)
  XFEM::BLITZTINY::MtV_product<3,3>(Tmid_,deltabetaminusalpha,curvnew_);
  XFEM::BLITZTINY::V_scale<3>(curvnew_,1/lrefe_);
  curvnew_ += curvold_;

  //computing current axial and shear strain epsilon, Crisfield, Vol. 2, equation (17.67)
  XFEM::BLITZTINY::MtV_product<3,3>(Tnew_,x21,epsilonn);
  XFEM::BLITZTINY::V_scale<3>(epsilonn,1/lrefe_);
  epsilonn(0) -=  1;

  /* read material parameters using structure _MATERIAL which is defined by inclusion of      /
   / "../drt_lib/drt_timecurve.H"; note: material parameters have to be read in the evaluation /
   / function instead of e.g. Beam3_input.cpp or within the Beam3Register class since it is not/
   / sure that structure _MATERIAL is declared within those scopes properly whereas it is within/
   / the evaluation functions */

  // get the material law
  MATERIAL* currmat = &(mat[material_-1]);
  double ym;
  double sm;
  double density;

  //assignment of material parameters; only St.Venant material is accepted for this beam 
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
  
 
  //stress values n and m, Crisfield, Vol. 2, equation (17.76) and (17.78)
  epsilonn(0) *= ym*crosssec_;
  epsilonn(1) *= sm*crosssecshear_;
  epsilonn(2) *= sm*crosssecshear_;
  XFEM::BLITZTINY::MV_product<3,3>(Tnew_,epsilonn,stressn);  

  //turning bending strain epsilonm into bending stress stressm
  epsilonm = curvnew_;
  epsilonm(0) *= sm*Irr_;
  epsilonm(1) *= ym*Iyy_;
  epsilonm(2) *= ym*Izz_;
  XFEM::BLITZTINY::MV_product<3,3>(Tnew_,epsilonm,stressm); 
  
  //computing spin matrix S(x21)/2 according to Crisfield, Vol. 2, equation (17.74)
  BlitzMat3x3 spinx21;
  computespin(spinx21,x21,0.5);
  

  //computing global internal forces, Crisfield Vol. 2, equation (17.79)
  //note: X = [-I 0; -S -I; I 0; -S I] with -S = T^t; and S = S(x21)/2;
  if (force != NULL)
  {
    (*force).Size(12);
    for (int i=0; i<3; ++i)
    {
      (*force)(i)   -= stressn(i);
      (*force)(i+3) -= stressm(i);
      (*force)(i+6) += stressn(i);
      (*force)(i+9) += stressm(i);
      
      for (int j=0; j<3; ++j)
      {      
        (*force)(i+3) -= stressn(j)*spinx21(i,j);      
        (*force)(i+9) -= stressn(j)*spinx21(i,j);       
      }
    }
  }

  //computing linear stiffness matrix
  if (stiffmatrix != NULL)
  {
    //setting constitutive parameters , Crisfield, Vol. 2, equation (17.76)
    Cm(0) = ym*crosssec_/lrefe_;
    Cm(1) = sm*crosssecshear_/lrefe_;
    Cm(2) = sm*crosssecshear_/lrefe_;
    Cb(0) = sm*Irr_/lrefe_;
    Cb(1) = ym*Iyy_/lrefe_;
    Cb(2) = ym*Izz_/lrefe_;
    
    //setting up basis of stiffness matrix according to Crisfield, Vol. 2, equation (17.81)   
    (*stiffmatrix).Shape(12,12);  
    computestiffbasis(Tnew_,Cm,Cb,spinx21,(*stiffmatrix));
    
     
    //adding nonlinear (stress dependent) parts to tangent stiffness matrix, Crisfield, Vol. 2 equs. (17.83), (17.87), (17.89)
    computeKsig1(Ksig1,stressn,stressm);
    computeKsig2(Ksig2,stressn,x21);
    (*stiffmatrix) += Ksig1;
    (*stiffmatrix) += Ksig2;

      //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
      //correctly or not: the function b3_nlnstiff_approx(mydisp) calculates the stiffness matrix approximated by
      //finite differences and finally the relative error is printed; in case that there is no significant error in any
      //element no printout is thrown
      //activating this part of code also the function b3_nlnstiff_approx(mydisp) has to be activated both in Beam3.H
      //and Beam3_evaluate.cpp
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
  
  /*calculating mass matrix; this beam3 element includes only a lumped mass matrix where for torsion and
   * bending the same moments of inertia are assumed; for slender beams the influence of rotational moments
   * of inertia is dilute so that often rotational inertia is just set to zero; since within this code at one
   * point an LU-decomposition of the mass matrix is carried out this was avoided in the here implemented beam3
   * element and instead the above described simplification was assumed; Iyy_ = Izz_ was assumed for similar 
   * reasons */
  if (massmatrix != NULL)
  {
    (*massmatrix).Shape(12,12);
    for (int i=0; i<3; ++i)
    {
      (*massmatrix)(i  ,i  ) = 0.5*density*lrefe_*crosssec_;
      (*massmatrix)(i+3,i+3) = 0.5*density*lrefe_*Iyy_;
      (*massmatrix)(i+6,i+6) = 0.5*density*lrefe_*crosssec_;
      (*massmatrix)(i+9,i+9) = 0.5*density*lrefe_*Iyy_;
    }
  }
  
  return;
} // DRT::ELEMENTS::Beam3::b3_nlnstiffmass


// lump mass matrix
void DRT::ELEMENTS::Beam3::b3_lumpmass(Epetra_SerialDenseMatrix* emass)
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
Epetra_SerialDenseMatrix DRT::ELEMENTS::Beam3::b3_nlnstiff_approx(BlitzMat6x2 xcurr, double h_rel, Epetra_SerialDenseVector force)
{

  BlitzVec3 epsilonn_aux;
  BlitzVec3 epsilonm_aux;
  BlitzVec3 stressn_aux;
  BlitzVec3 stressm_aux;
  BlitzVec3 x21_aux;
  BlitzVec3 betaplusalphanew;
  BlitzVec3 betaminusalphanew;
  BlitzVec3 deltabetaplusalpha;
  BlitzVec3 deltabetaminusalpha;
  BlitzVec3 curvnew;
  BlitzMat3x3 Saux;
  BlitzMat3x3 Raux;
  BlitzMat3x3 Tnew;
  BlitzMat3x3 Tmid;
  BlitzMat6x2 xcurr_aux;
  Epetra_SerialDenseMatrix stiff_approx;
  Epetra_SerialDenseVector force_aux;
  stiff_approx.Shape(12,12);
    
  //calculating strains in new configuration
  for(int i=0; i<6; i++)
  {
    for(int k=0; k<2; k++)
    {
      xcurr_aux = xcurr;
      xcurr_aux(i,k) = xcurr_aux(i,k) + h_rel;
      
      
      //difference between coordinates of both nodes, x21_aux' Crisfield  Vol. 2 equ. (17.66a) and (17.72)
      for (int j=0; j<3; ++j)
      {
        x21_aux(j) = xcurr_aux(j,1) - xcurr_aux(j,0);
        betaplusalphanew(j) = xcurr_aux(j+3,1) + xcurr_aux(j+3,0);
        betaminusalphanew(j) = xcurr_aux(j+3,1) - xcurr_aux(j+3,0);
      }
      //the basis configuration is the new one of the actual calculation process; for numerical approximation
      //other configurations rotated out of this one by h_rel are considered
      deltabetaplusalpha  = betaplusalphanew;
      deltabetaplusalpha -= betaplusalphanew_;
    
      deltabetaminusalpha  = betaminusalphanew;
      deltabetaminusalpha -= betaminusalphanew_;

      //auxiliary matrix for update of Tnew_, Crisfield, Vol 2, equation (17.65)
      computespin(Saux,deltabetaplusalpha, 0.5);
      computerotation(Raux,Saux);
      XFEM::BLITZTINY::MM_product<3,3,3>(Raux,Tnew_,Tnew);

      //computing triad for curvature update, Crisfield, Vol 2, equation (17.73)
      computespin(Saux,deltabetaplusalpha, 0.25);
      computerotation(Raux,Saux);
      XFEM::BLITZTINY::MM_product<3,3,3>(Raux,Tnew_,Tmid);
      
      //updating curvature, Crisfield, Vol. 2, equation (17.72)
      XFEM::BLITZTINY::MtV_product<3,3>(Tmid,deltabetaminusalpha,curvnew);
      XFEM::BLITZTINY::V_scale<3>(curvnew,1/lrefe_);
    
      curvnew += curvnew_;
      epsilonm_aux = curvnew;

      //computing current axial and shear strain epsilon, Crisfield, Vol. 2, equation (17.67)
      XFEM::BLITZTINY::MtV_product<3,3>(Tnew,x21_aux,epsilonn_aux);
      XFEM::BLITZTINY::V_scale<3>(epsilonn_aux,1/lrefe_);
      epsilonn_aux(0) -=  1;

      // get the material law
      MATERIAL* currmat = &(mat[material_-1]);
      double ym;
      double sm;
      double density;

      //assignment of material parameters; only St.Venant material is accepted for this beam 
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

      //stress values n and m, Crisfield, Vol. 2, equation (17.76) and (17.78)
      epsilonn_aux(0) *= ym*crosssec_;
      epsilonn_aux(1) *= sm*crosssecshear_;
      epsilonn_aux(2) *= sm*crosssecshear_;
      XFEM::BLITZTINY::MV_product<3,3>(Tnew,epsilonn_aux,stressn_aux);  

      //turning bending strain epsilonm into bending stress stressm
      epsilonm_aux = curvnew;
      epsilonm_aux(0) *= sm*Irr_;
      epsilonm_aux(1) *= ym*Iyy_;
      epsilonm_aux(2) *= ym*Izz_;
      XFEM::BLITZTINY::MV_product<3,3>(Tnew,epsilonm_aux,stressm_aux); 
        
      //computing spin matrix S(x21_aux)/2 according to Crisfield, Vol. 2, equation (17.74)
      BlitzMat3x3 spinx21_aux;
      computespin(spinx21_aux,x21_aux,0.5);

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
} // DRT::ELEMENTS::Beam3::b3_nlnstiff_approx

*/

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM3
