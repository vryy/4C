/*!----------------------------------------------------------------------
\file so_surface_evaluate.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include <blitz/array.h>
#include "so_surface.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_surfstress/drt_surfstress_manager.H"
#include "../drt_surfstress/drt_potential_manager.H"
#include "../drt_statmech/bromotion_manager.H"

using UTILS::SurfStressManager;
using UTILS::PotentialManager;

/*----------------------------------------------------------------------*
 * Integrate a Surface Neumann boundary condition (public)     gee 04/08|
 * ---------------------------------------------------------------------*/
int DRT::ELEMENTS::StructuralSurface::EvaluateNeumann(ParameterList&           params,
                                                      DRT::Discretization&     discretization,
                                                      DRT::Condition&          condition,
                                                      vector<int>&             lm,
                                                      Epetra_SerialDenseVector& elevec1)
{
  // get type of condition
  enum LoadType
  {
    neum_none,
    neum_live,
    neum_orthopressure
  };
  // spatial or material configuration depends on the type of load
  enum Configuration
  {
    config_none,
    config_material,
    config_spatial,
    config_both
  };

  Configuration config = config_none;
  LoadType ltype       = neum_none;
  const string* type = condition.Get<string>("type");
  if      (*type == "neum_live")
  {
    ltype = neum_live;
    config = config_material;
  }
  else if (*type == "neum_orthopressure")
  {
    ltype = neum_orthopressure;
    config = config_spatial;
  }
  else
  {
    dserror("Unknown type of SurfaceNeumann condition");
  }

  // get values and switches from the condition
  const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val"  );

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve  = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

  // element geometry update
  const int numnode = NumNode();
  const int numdf=3;
  LINALG::SerialDenseMatrix x(numnode,numdf);
  LINALG::SerialDenseMatrix xc;
  switch (config)
  {
  case config_material:
    MaterialConfiguration(x);
  break;
  case config_spatial:
  {
    xc.LightShape(numnode,numdf);
#ifndef INVERSEDESIGNCREATE
    RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
    if (disp==null) dserror("Cannot get state vector 'displacement'");
    vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
    SpatialConfiguration(xc,mydisp);
#else
// in inverse design analysis, the current configuration is the reference
    MaterialConfiguration(xc);
#endif
  }
  break;
  case config_both:
  {
    xc.LightShape(numnode,3);
    RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
    if (disp==null) dserror("Cannot get state vector 'displacement'");
    vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
    MaterialConfiguration(x);
    SpatialConfiguration(xc,mydisp);
  }
  break;
  default: dserror("Unknown case of frame");
  break;
  }

  // allocate vector for shape functions and matrix for derivatives
  LINALG::SerialDenseVector  funct(numnode);
  LINALG::SerialDenseMatrix  deriv(2,numnode);

  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  const DRT::UTILS::IntegrationPoints2D  intpoints =
                                       getIntegrationPoints2D(gaussrule_);
  for (int gp=0; gp<intpoints.nquad; gp++)
  {
    const double e0 = intpoints.qxg[gp][0];
    const double e1 = intpoints.qxg[gp][1];

    // get shape functions and derivatives in the plane of the element
    DRT::UTILS::shape_function_2D(funct,e0,e1,Shape());
    DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());

    switch(ltype)
    {
    case neum_live:
    {
      LINALG::SerialDenseMatrix dxyzdrs(2,3);
      dxyzdrs.Multiply('N','N',1.0,deriv,x,0.0);
      LINALG::SerialDenseMatrix  metrictensor(2,2);
      metrictensor.Multiply('N','T',1.0,dxyzdrs,dxyzdrs,0.0);
      const double detA = sqrt( metrictensor(0,0)*metrictensor(1,1)
                                -metrictensor(0,1)*metrictensor(1,0));
      const double fac = intpoints.qwgt[gp] * detA * curvefac;
      for (int node=0; node < numnode; ++node)
        for(int dim=0 ; dim<3; ++dim)
          elevec1[node*numdf+dim]+= funct[node] * (*onoff)[dim] * (*val)[dim] * fac;
    }
    break;
    case neum_orthopressure:
    {
     if ((*onoff)[0] != 1) dserror("orthopressure on 1st dof only!");
      for (int checkdof = 1; checkdof < 3; ++checkdof)
        if ((*onoff)[checkdof] != 0) dserror("orthopressure on 1st dof only!");
      double ortho_value = (*val)[0];
      if (!ortho_value) dserror("no orthopressure value given!");
      vector<double> normal(3);
      SurfaceIntegration(normal,xc,deriv);
      const double fac = intpoints.qwgt[gp] * curvefac * ortho_value;
      for (int node=0; node < numnode; ++node)
        for(int dim=0 ; dim<3; dim++)
          elevec1[node*numdf+dim] += funct[node] * normal[dim] * fac;
    }
    break;
    default:
      dserror("Unknown type of SurfaceNeumann load");
    break;
    }

  } /* end of loop over integration points gp */

  return 0;
}

/*----------------------------------------------------------------------*
 * Evaluate normal at gp (private)                             gee 08/08|
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::SurfaceIntegration(vector<double>& normal,
                                                          const Epetra_SerialDenseMatrix& x,
                                                          const Epetra_SerialDenseMatrix& deriv)
{
  // note that the length of this normal is the area dA

  // compute dXYZ / drs
  LINALG::SerialDenseMatrix dxyzdrs(2,3);
  dxyzdrs.Multiply('N','N',1.0,deriv,x,0.0);

  normal[0] = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
  normal[1] = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
  normal[2] = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);

  return;
}

/*----------------------------------------------------------------------*
 * Evaluate sqrt of determinant of metric at gp (private)      gee 04/08|
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::SurfaceIntegration(double& detA,
                                                          vector<double>& normal,
                                                          const Epetra_SerialDenseMatrix& x,
                                                          const Epetra_SerialDenseMatrix& deriv)
{

  // compute dXYZ / drs
  LINALG::SerialDenseMatrix dxyzdrs(2,3);
  dxyzdrs.Multiply('N','N',1.0,deriv,x,0.0);

  /* compute covariant metric tensor G for surface element
  **                        | g11   g12 |
  **                    G = |           |
  **                        | g12   g22 |
  ** where (o denotes the inner product, xyz a vector)
  **
  **       dXYZ   dXYZ          dXYZ   dXYZ          dXYZ   dXYZ
  ** g11 = ---- o ----    g12 = ---- o ----    g22 = ---- o ----
  **        dr     dr            dr     ds            ds     ds
  */
  LINALG::SerialDenseMatrix metrictensor(2,2);
  metrictensor.Multiply('N','T',1.0,dxyzdrs,dxyzdrs,0.0);
  detA = sqrt( metrictensor(0,0)*metrictensor(1,1)-metrictensor(0,1)*metrictensor(1,0) );
  normal[0] = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
  normal[1] = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
  normal[2] = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);

  return;
}

/*----------------------------------------------------------------------*
 * Evaluate method for StructuralSurface-Elements               tk 10/07*
 * ---------------------------------------------------------------------*/
int DRT::ELEMENTS::StructuralSurface::Evaluate(ParameterList&            params,
                                               DRT::Discretization&      discretization,
                                               vector<int>&              lm,
                                               Epetra_SerialDenseMatrix& elematrix1,
                                               Epetra_SerialDenseMatrix& elematrix2,
                                               Epetra_SerialDenseVector& elevector1,
                                               Epetra_SerialDenseVector& elevector2,
                                               Epetra_SerialDenseVector& elevector3)
{
  const DiscretizationType distype = Shape();

  // start with "none"
  DRT::ELEMENTS::StructuralSurface::ActionType act = StructuralSurface::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_constrvol")        act = StructuralSurface::calc_struct_constrvol;
  else if (action=="calc_struct_volconstrstiff")   act = StructuralSurface::calc_struct_volconstrstiff;
  else if (action=="calc_struct_monitarea")        act = StructuralSurface::calc_struct_monitarea;
  else if (action=="calc_struct_constrarea")       act = StructuralSurface::calc_struct_constrarea;
  else if (action=="calc_struct_areaconstrstiff")  act = StructuralSurface::calc_struct_areaconstrstiff;
  else if (action=="calc_init_vol")                act = StructuralSurface::calc_init_vol;
  else if (action=="calc_surfstress_stiff")        act = StructuralSurface::calc_surfstress_stiff;
  else if (action=="calc_potential_stiff")         act = StructuralSurface::calc_potential_stiff;
  else if (action=="calc_brownian_motion")         act = StructuralSurface::calc_brownian_motion;
  else 
  {
    cout << action << endl;
    dserror("Unknown type of action for StructuralSurface");
  }

  //create communicator
  const Epetra_Comm& Comm = discretization.Comm();
  // what the element has to do
  switch(act)
  {
      //just compute the enclosed volume (e.g. for initialization)
      case calc_struct_constrvol:
      {
        if (distype!=quad4)
          dserror("Volume Constraint so far only works for quad4 surfaces!");
        //We are not interested in volume of ghosted elements
        if(Comm.MyPID()==Owner())
        {
          // element geometry update
          RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
          if (disp==null) dserror("Cannot get state vector 'displacement'");
          vector<double> mydisp(lm.size());
          DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
          const int numdim = 3;
          LINALG::SerialDenseMatrix xscurr(NumNode(),numdim);  // material coord. of element
          SpatialConfiguration(xscurr,mydisp);
          //call submethod for volume evaluation and store rseult in third systemvector
          double volumeele = ComputeConstrVols(xscurr);
          elevector3[0]=volumeele;
        }

      }
      break;
      case calc_struct_volconstrstiff:
      {
        if (distype!=quad4)
          dserror("Volume Constraint only works for quad4 surfaces!");
        // element geometry update
        RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp==null) dserror("Cannot get state vector 'displacement'");
        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
        const int numdim =3;
        LINALG::SerialDenseMatrix xscurr(NumNode(),numdim);  // material coord. of element
        SpatialConfiguration(xscurr,mydisp);
        //call submethods
        ComputeVolConstrStiff(xscurr,elematrix1);
        ComputeVolConstrDeriv(xscurr,elevector1);
        //update corresponding column in "constraint" matrix
        elevector2=elevector1;
        //call submethod for volume evaluation and store rseult in third systemvector
        double volumeele = ComputeConstrVols(xscurr);
        elevector3[0]=volumeele;
      }

      break;
      case calc_init_vol:
      {
        // the reference volume of the RVE (including inner
        // holes) is calculated by evaluating the following
        // surface integral:
        // V = 1/3*int(div(X))dV = 1/3*int(N*X)dA
        // with X being the reference coordinates and N the
        // normal vector of the surface element (exploiting the
        // fact that div(X)=1.0)

        // this is intended to be used in the serial case (microstructure)

        // NOTE: there must not be any holes penetrating the boundary!

        double V = params.get<double>("V0", 0.0);
        double dV = 0.0;
        const int numnode = NumNode();
        LINALG::SerialDenseMatrix x(numnode,3);
        MaterialConfiguration(x);

        // allocate vector for shape functions and matrix for derivatives
        LINALG::SerialDenseVector  funct(numnode);
        LINALG::SerialDenseMatrix  deriv(2,numnode);

        /*----------------------------------------------------------------------*
         |               start loop over integration points                     |
         *----------------------------------------------------------------------*/
        const DRT::UTILS::IntegrationPoints2D  intpoints =
          getIntegrationPoints2D(gaussrule_);

        for (int gp=0; gp<intpoints.nquad; gp++)
        {
          const double e0 = intpoints.qxg[gp][0];
          const double e1 = intpoints.qxg[gp][1];

          // get shape functions and derivatives in the plane of the element
          DRT::UTILS::shape_function_2D(funct,e0,e1,Shape());
          DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());

          vector<double> normal(3);
          double detA;
          SurfaceIntegration(detA,normal,x,deriv);
          const double fac = intpoints.qwgt[gp] * detA;

          double temp = 0.0;
          vector<double> X(3,0.);

          for (int i=0; i<numnode; i++)
          {
            X[0] += funct[i]*x(i,0);
            X[1] += funct[i]*x(i,1);
            X[2] += funct[i]*x(i,2);
          }

          for (int i=0;i<3;++i)
          {
            temp += normal[i]*normal[i];
          }

          if (temp<0.)
            dserror("calculation of initial volume failed in surface element");
          double absnorm = sqrt(temp);

          for (int i=0;i<3;++i)
          {
            normal[i] /= absnorm;
          }
          for (int i=0;i<3;++i)
          {
            dV += 1/3.0*fac*normal[i]*X[i];
          }
        }
        params.set("V0", V+dV);
      }
      break;

      case calc_surfstress_stiff:
      {
        RefCountPtr<SurfStressManager> surfstressman =
          params.get<RefCountPtr<SurfStressManager> >("surfstr_man", null);

        if (surfstressman==null)
          dserror("No SurfStressManager in Solid3 Surface available");

        RefCountPtr<DRT::Condition> cond = params.get<RefCountPtr<DRT::Condition> >("condition",null);
        if (cond==null)
          dserror("Condition not available in Solid3 Surface");

        if (distype!=quad4)
          cout << "Surface Stresses were only tested for quad4 surfaces! Use with caution!" << endl;

        double time = params.get<double>("total time",-1.0);
        double dt = params.get<double>("delta time",0.0);
        double alphaf = params.get<double>("alpha f",0.0);
        bool newstep = params.get<bool>("newstep", false);
        bool fintliketr = params.get<bool>("fintliketr", false);

        // element geometry update

        const int numnode = NumNode();
        LINALG::SerialDenseMatrix x(numnode,3);

        if (fintliketr)
        {
          RefCountPtr<const Epetra_Vector> disn = discretization.GetState("new displacement");
          if (disn==null) dserror("Cannot get state vector 'new displacement'");
          vector<double> mydisn(lm.size());
          DRT::UTILS::ExtractMyValues(*disn,mydisn,lm);
          SpatialConfiguration(x,mydisn);
        }
        else
        {
          RefCountPtr<const Epetra_Vector> dism = discretization.GetState("displacement");
          if (dism==null) dserror("Cannot get state vector 'displacement'");
          vector<double> mydism(lm.size());
          DRT::UTILS::ExtractMyValues(*dism,mydism,lm);
          SpatialConfiguration(x,mydism);
        }

        const DRT::UTILS::IntegrationPoints2D  intpoints =
          getIntegrationPoints2D(gaussrule_);

        // set up matrices and parameters needed for the evaluation of current
        // interfacial area and its derivatives w.r.t. the displacements

        int ndof = 3*numnode;                                     // overall number of surface dofs
        double A = 0.;                                            // interfacial area
        // we really want to zero out the following matrices -> no LINALG::SerialDenseMatrix
        // first partial derivatives
        RCP<Epetra_SerialDenseVector> Adiff = rcp(new Epetra_SerialDenseVector(ndof));
        // second partial derivatives
        RCP<Epetra_SerialDenseMatrix> Adiff2 = rcp(new Epetra_SerialDenseMatrix(ndof,ndof)); // second partial derivatives

        ComputeAreaDeriv(x, numnode, ndof, A, Adiff, Adiff2);

        if (cond->Type()==DRT::Condition::Surfactant)     // dynamic surfactant model
        {
          int curvenum = cond->Getint("curve");
          double k1xC = cond->GetDouble("k1xCbulk");
          double k2 = cond->GetDouble("k2");
          double m1 = cond->GetDouble("m1");
          double m2 = cond->GetDouble("m2");
          double gamma_0 = cond->GetDouble("gamma_0");
          double gamma_min = cond->GetDouble("gamma_min");
          double gamma_min_eq = cond->GetDouble("gamma_min_eq");
          double con_quot_max = (gamma_min_eq-gamma_min)/m2+1.;
          double con_quot_eq = (k1xC)/(k1xC+k2);

          if (fintliketr)
          {
            surfstressman->StiffnessAndInternalForces(curvenum, A, Adiff, Adiff2, A, Adiff, elevector1, elematrix1, this->Id(),
                                                      time, dt, 0, 0.0, k1xC, k2, m1, m2, gamma_0,
                                                      gamma_min, gamma_min_eq, con_quot_max,
                                                      con_quot_eq, alphaf, newstep, fintliketr);
          }
          else
          {
            // element geometry update (n+1)
            RefCountPtr<const Epetra_Vector> disn = discretization.GetState("new displacement");
            if (disn==null) dserror("Cannot get state vector 'new displacement'");
            vector<double> mydisn(lm.size());
            DRT::UTILS::ExtractMyValues(*disn,mydisn,lm);
            SpatialConfiguration(x,mydisn);

            // set up matrices and parameters needed for the evaluation of
            // interfacial area and its first derivative w.r.t. the displacements at (n+1)
            double Anew = 0.;                                            // interfacial area
            RCP<Epetra_SerialDenseVector> Adiffnew = rcp(new Epetra_SerialDenseVector(ndof));

            ComputeAreaDeriv(x, numnode, ndof, Anew, Adiffnew, null);


            surfstressman->StiffnessAndInternalForces(curvenum, A, Adiff, Adiff2, Anew, Adiffnew, elevector1, elematrix1, this->Id(),
                                                      time, dt, 0, 0.0, k1xC, k2, m1, m2, gamma_0,
                                                      gamma_min, gamma_min_eq, con_quot_max,
                                                      con_quot_eq, alphaf, newstep, fintliketr);
          }
        }
        else if (cond->Type()==DRT::Condition::SurfaceTension) // ideal liquid
        {
          int curvenum = cond->Getint("curve");
          double const_gamma = cond->GetDouble("gamma");
          RCP<Epetra_SerialDenseVector> Adiffnew = rcp(new Epetra_SerialDenseVector(ndof));
          surfstressman->StiffnessAndInternalForces(curvenum, A, Adiff, Adiff2, 0., Adiffnew, elevector1, elematrix1, this->Id(),
                                                    time, dt, 1, const_gamma, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0, 0.0, alphaf, newstep);
        }
        else
          dserror("Unknown condition type %d",cond->Type());
      }
      break;

      // compute additional stresses due to intermolecular potential forces
      case calc_potential_stiff:
      {
        RefCountPtr<PotentialManager> potentialmanager =
          params.get<RefCountPtr<PotentialManager> >("pot_man", null);
        if (potentialmanager==null)
          dserror("No PotentialManager in Solid3 Surface available");

        RefCountPtr<DRT::Condition> cond = params.get<RefCountPtr<DRT::Condition> >("condition",null);
        if (cond==null)
          dserror("Condition not available in Solid3 Surface");
     
        if (cond->Type()==DRT::Condition::LJ_Potential) // Lennard-Jones potential
        { 
          potentialmanager->StiffnessAndInternalForcesLJPotential(this, gaussrule_, params,lm, elematrix1, elevector1);
        }
        else
          dserror("Unknown condition type %d",cond->Type());
      }
      break;
      
      // compute stochastical forces due to Brownian Motion
      case calc_brownian_motion:
      {
        RefCountPtr<BroMotion_Manager> bromotion_manager =
          params.get<RefCountPtr<BroMotion_Manager> >("bromo_man", null);
        if (bromotion_manager==null)
          dserror("No Brownian Manager in Solid3 Surface available");
         
        RefCountPtr<DRT::Condition> cond = params.get<RefCountPtr<DRT::Condition> >("condition",null);
         if (cond==null)
           dserror("Condition not available in Solid3 Surface");
      
        if (cond->Type()==DRT::Condition::Brownian_Motion) // Brownian Motion
        { 
          bromotion_manager->StochasticalForces(this, gaussrule_, params,lm, elematrix1, elevector1);
        }
        else
          dserror("Unknown condition type %d",cond->Type());
      }
      break;

      //compute the area (e.g. for initialization)
      case calc_struct_monitarea:
      {
        if (distype!=quad4)
          dserror("Area monitor only works for quad4 surfaces!");
        //We are not interested in volume of ghosted elements
        if(Comm.MyPID()==Owner())
        {
          // element geometry update
          RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
          if (disp==null) dserror("Cannot get state vector 'displacement'");
          vector<double> mydisp(lm.size());
          DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
          const int numdim = 3;
          LINALG::SerialDenseMatrix xscurr(NumNode(),numdim);  // material coord. of element
          SpatialConfiguration(xscurr,mydisp);


          //get required projection method
          enum ProjType
          {
            none,
            xy,
            yz,
            xz
          };
          ProjType protype;

          RCP<DRT::Condition> condition = params.get<RefCountPtr<DRT::Condition> >("condition");
          const string* type = condition->Get<string>("projection");

          protype = none;
          if (!type);
          else if (*type == "xy") protype = xy;
          else if (*type == "yz") protype = yz;
          else if (*type == "xz") protype = xz;

          // To compute monitored area consider required projection method
          // and set according coordinates to zero
          if (protype == yz)
          {
            xscurr(0,0)=0;
            xscurr(1,0)=0;
            xscurr(2,0)=0;
            xscurr(3,0)=0;
          }
          else if (protype == xz)
          {
            xscurr(0,1)=0;
            xscurr(1,1)=0;
            xscurr(2,1)=0;
            xscurr(3,1)=0;
          }
          else if (protype == xy)
          {
            xscurr(0,2)=0;
            xscurr(1,2)=0;
            xscurr(2,2)=0;
            xscurr(3,2)=0;
          }

          double areaele=0.0;
          const DRT::UTILS::IntegrationPoints2D  intpoints =
                                               getIntegrationPoints2D(gaussrule_);
          // allocate matrix for derivatives of shape functions
          LINALG::SerialDenseMatrix  deriv(2,NumNode());

          //Compute area
          for (int gp=0; gp<intpoints.nquad; gp++)
          {
            const double e0 = intpoints.qxg[gp][0];
            const double e1 = intpoints.qxg[gp][1];

            // get shape functions and derivatives in the plane of the element
            DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());

            vector<double> normal(3);
            double detA;
            SurfaceIntegration(detA,normal,xscurr,deriv);
            const double fac = intpoints.qwgt[gp] * detA;
            areaele += fac;

          }

          //store result in third systemvector
          elevector3[0]=areaele;
        }

      }
      break;
      case calc_struct_constrarea:
      {
        dserror("Element routines for area constraint in 3D not implemented yet!");
      }
      break;
      case calc_struct_areaconstrstiff:
      {
        dserror("Element routines for area constraint in 3D not implemented yet!");
      }
      break;
      default:
        dserror("Unimplemented type of action for StructuralSurface");

  }
  return 0;
  }

/*----------------------------------------------------------------------*
 * Compute Volume between surface and xy-plane.                 tk 10/07*
 * Yields to the enclosed volume when summed up over all elements       *
 * ---------------------------------------------------------------------*/
double DRT::ELEMENTS::StructuralSurface::ComputeConstrVols(const LINALG::SerialDenseMatrix& xc)
{
  double volume =0;
  //Formula for volume computation based on calculation of Ulrich done
  //within the old code
  volume =-(1.0/12.0)*(2*xc(0,0)*xc(1,1)*xc(0,2) -
      2*xc(1,0)*xc(0,1)*xc(0,2) +
      2*xc(0,0)*xc(1,1)*xc(1,2) -
      2*xc(1,0)*xc(0,1)*xc(1,2) +
      xc(0,0)*xc(1,1)*xc(2,2) -
      2*xc(0,0)*xc(0,2)*xc(3,1) -
      xc(0,0)*xc(2,1)*xc(1,2) -
      xc(1,0)*xc(0,1)*xc(2,2) +
      xc(1,0)*xc(0,2)*xc(2,1) +
      xc(0,1)*xc(2,0)*xc(1,2) +
      2*xc(0,1)*xc(0,2)*xc(3,0) -
      xc(2,0)*xc(1,1)*xc(0,2) +
      xc(0,0)*xc(1,1)*xc(3,2) -
      xc(0,0)*xc(1,2)*xc(3,1) -
      xc(1,0)*xc(0,1)*xc(3,2) +
      xc(1,0)*xc(0,2)*xc(3,1) +
      2*xc(1,0)*xc(2,1)*xc(1,2) +
      xc(0,1)*xc(3,0)*xc(1,2) -
      2*xc(2,0)*xc(1,1)*xc(1,2) -
      xc(1,1)*xc(0,2)*xc(3,0) +
      xc(0,0)*xc(2,1)*xc(3,2) -
      xc(0,0)*xc(3,1)*xc(2,2) +
      2*xc(1,0)*xc(2,1)*xc(2,2) -
      xc(0,1)*xc(2,0)*xc(3,2) +
      xc(0,1)*xc(3,0)*xc(2,2) -
      2*xc(2,0)*xc(1,1)*xc(2,2) +
      xc(2,0)*xc(0,2)*xc(3,1) -
      xc(0,2)*xc(3,0)*xc(2,1) -
      2*xc(0,0)*xc(3,1)*xc(3,2) +
      xc(1,0)*xc(2,1)*xc(3,2) -
      xc(1,0)*xc(3,1)*xc(2,2) +
      2*xc(0,1)*xc(3,0)*xc(3,2) -
      xc(2,0)*xc(1,1)*xc(3,2) +
      xc(2,0)*xc(1,2)*xc(3,1) +
      xc(1,1)*xc(3,0)*xc(2,2) -
      xc(3,0)*xc(2,1)*xc(1,2) +
      2*xc(2,0)*xc(3,1)*xc(2,2) -
      2*xc(3,0)*xc(2,1)*xc(2,2) +
      2*xc(2,0)*xc(3,1)*xc(3,2) -
      2*xc(3,0)*xc(2,1)*xc(3,2));
  return volume;
}

/*----------------------------------------------------------------------*
 * Compute influence of volume constraint on stiffness matrix.  tk 10/07*
 * Second derivatives of volume with respect to the displacements       *
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::ComputeVolConstrStiff(const LINALG::SerialDenseMatrix& xc,
    Epetra_SerialDenseMatrix& elematrix)
{
  //Second derivatives of volume with respect to the displacements.
  //Only suitable for quad4 surfaces, since based on a symbolic calculation with mupad
  elematrix(0,0) = 0.0 ;
  elematrix(0,1) = 0.0 ;
  elematrix(0,2) = xc(1,1)*(1.0/6.0) - xc(3,1)*(1.0/6.0) ;
  elematrix(0,3) = 0.0 ;
  elematrix(0,4) = xc(0,2)*(1.0/6.0) + xc(1,2)*(1.0/6.0) + xc(2,2)*(1.0/12.0) + xc(3,2)*(1.0/12.0) ;
  elematrix(0,5) = xc(1,1)*(1.0/6.0) - xc(2,1)*(1.0/12.0) - xc(3,1)*(1.0/12.0) ;
  elematrix(0,6) = 0.0 ;
  elematrix(0,7) = xc(1,2)*(-1.0/12.0) + xc(3,2)*(1.0/12.0) ;
  elematrix(0,8) = xc(1,1)*(1.0/12.0) - xc(3,1)*(1.0/12.0) ;
  elematrix(0,9) = 0.0 ;
  elematrix(0,10) = xc(0,2)*(-1.0/6.0) - xc(1,2)*(1.0/12.0) - xc(2,2)*(1.0/12.0) - xc(3,2)*(1.0/6.0) ;
  elematrix(0,11) = xc(1,1)*(1.0/12.0) + xc(2,1)*(1.0/12.0) - xc(3,1)*(1.0/6.0) ;

  elematrix(1,0) = 0.0 ;
  elematrix(1,1) = 0.0 ;
  elematrix(1,2) = xc(1,0)*(-1.0/6.0) + xc(3,0)*(1.0/6.0) ;
  elematrix(1,3) = xc(0,2)*(-1.0/6.0) - xc(1,2)*(1.0/6.0) - xc(2,2)*(1.0/12.0) - xc(3,2)*(1.0/12.0) ;
  elematrix(1,4) = 0.0 ;
  elematrix(1,5) = xc(1,0)*(-1.0/6.0) + xc(2,0)*(1.0/12.0) + xc(3,0)*(1.0/12.0) ;
  elematrix(1,6) = xc(1,2)*(1.0/12.0) - xc(3,2)*(1.0/12.0) ;
  elematrix(1,7) = 0.0 ;
  elematrix(1,8) = xc(1,0)*(-1.0/12.0) + xc(3,0)*(1.0/12.0) ;
  elematrix(1,9) = xc(0,2)*(1.0/6.0) + xc(1,2)*(1.0/12.0) + xc(2,2)*(1.0/12.0) + xc(3,2)*(1.0/6.0) ;
  elematrix(1,10) = 0.0 ;
  elematrix(1,11) = xc(1,0)*(-1.0/12.0) - xc(2,0)*(1.0/12.0) + xc(3,0)*(1.0/6.0) ;

  elematrix(2,0) = xc(1,1)*(1.0/6.0) - xc(3,1)*(1.0/6.0) ;
  elematrix(2,1) = xc(1,0)*(-1.0/6.0) + xc(3,0)*(1.0/6.0) ;
  elematrix(2,2) = 0.0 ;
  elematrix(2,3) = xc(0,1)*(-1.0/6.0) + xc(2,1)*(1.0/12.0) + xc(3,1)*(1.0/12.0) ;
  elematrix(2,4) = xc(0,0)*(1.0/6.0) - xc(2,0)*(1.0/12.0) - xc(3,0)*(1.0/12.0) ;
  elematrix(2,5) = 0.0 ;
  elematrix(2,6) = xc(1,1)*(-1.0/12.0) + xc(3,1)*(1.0/12.0) ;
  elematrix(2,7) = xc(1,0)*(1.0/12.0) - xc(3,0)*(1.0/12.0) ;
  elematrix(2,8) = 0.0 ;
  elematrix(2,9) = xc(0,1)*(1.0/6.0) - xc(1,1)*(1.0/12.0) - xc(2,1)*(1.0/12.0) ;
  elematrix(2,10) = xc(0,0)*(-1.0/6.0) + xc(1,0)*(1.0/12.0) + xc(2,0)*(1.0/12.0) ;
  elematrix(2,11) = 0.0 ;

  elematrix(3,0) = 0.0 ;
  elematrix(3,1) = xc(0,2)*(-1.0/6.0) - xc(1,2)*(1.0/6.0) - xc(2,2)*(1.0/12.0) - xc(3,2)*(1.0/12.0) ;
  elematrix(3,2) = xc(0,1)*(-1.0/6.0) + xc(2,1)*(1.0/12.0) + xc(3,1)*(1.0/12.0) ;
  elematrix(3,3) = 0.0 ;
  elematrix(3,4) = 0.0 ;
  elematrix(3,5) = xc(0,1)*(-1.0/6.0) + xc(2,1)*(1.0/6.0) ;
  elematrix(3,6) = 0.0 ;
  elematrix(3,7) = xc(0,2)*(1.0/12.0) + xc(1,2)*(1.0/6.0) + xc(2,2)*(1.0/6.0) + xc(3,2)*(1.0/12.0) ;
  elematrix(3,8) = xc(0,1)*(-1.0/12.0) + xc(2,1)*(1.0/6.0) - xc(3,1)*(1.0/12.0) ;
  elematrix(3,9) = 0.0 ;
  elematrix(3,10) = xc(0,2)*(1.0/12.0) - xc(2,2)*(1.0/12.0) ;
  elematrix(3,11) = xc(0,1)*(-1.0/12.0) + xc(2,1)*(1.0/12.0) ;

  elematrix(4,0) = xc(0,2)*(1.0/6.0) + xc(1,2)*(1.0/6.0) + xc(2,2)*(1.0/12.0) + xc(3,2)*(1.0/12.0) ;
  elematrix(4,1) = 0.0 ;
  elematrix(4,2) = xc(0,0)*(1.0/6.0) - xc(2,0)*(1.0/12.0) - xc(3,0)*(1.0/12.0) ;
  elematrix(4,3) = 0.0 ;
  elematrix(4,4) = 0.0 ;
  elematrix(4,5) = xc(0,0)*(1.0/6.0) - xc(2,0)*(1.0/6.0) ;
  elematrix(4,6) = xc(0,2)*(-1.0/12.0) - xc(1,2)*(1.0/6.0) - xc(2,2)*(1.0/6.0) - xc(3,2)*(1.0/12.0) ;
  elematrix(4,7) = 0.0 ;
  elematrix(4,8) = xc(0,0)*(1.0/12.0) - xc(2,0)*(1.0/6.0) + xc(3,0)*(1.0/12.0) ;
  elematrix(4,9) = xc(0,2)*(-1.0/12.0) + xc(2,2)*(1.0/12.0) ;
  elematrix(4,10) = 0.0 ;
  elematrix(4,11) = xc(0,0)*(1.0/12.0) - xc(2,0)*(1.0/12.0) ;

  elematrix(5,0) = xc(1,1)*(1.0/6.0) - xc(2,1)*(1.0/12.0) - xc(3,1)*(1.0/12.0) ;
  elematrix(5,1) = xc(1,0)*(-1.0/6.0) + xc(2,0)*(1.0/12.0) + xc(3,0)*(1.0/12.0) ;
  elematrix(5,2) = 0.0 ;
  elematrix(5,3) = xc(0,1)*(-1.0/6.0) + xc(2,1)*(1.0/6.0) ;
  elematrix(5,4) = xc(0,0)*(1.0/6.0) - xc(2,0)*(1.0/6.0) ;
  elematrix(5,5) = 0.0 ;
  elematrix(5,6) = xc(0,1)*(1.0/12.0) - xc(1,1)*(1.0/6.0) + xc(3,1)*(1.0/12.0) ;
  elematrix(5,7) = xc(0,0)*(-1.0/12.0) + xc(1,0)*(1.0/6.0) - xc(3,0)*(1.0/12.0) ;
  elematrix(5,8) = 0.0 ;
  elematrix(5,9) = xc(0,1)*(1.0/12.0) - xc(2,1)*(1.0/12.0) ;
  elematrix(5,10) = xc(0,0)*(-1.0/12.0) + xc(2,0)*(1.0/12.0) ;
  elematrix(5,11) = 0.0 ;

  elematrix(6,0) = 0.0 ;
  elematrix(6,1) = xc(1,2)*(1.0/12.0) - xc(3,2)*(1.0/12.0) ;
  elematrix(6,2) = xc(1,1)*(-1.0/12.0) + xc(3,1)*(1.0/12.0) ;
  elematrix(6,3) = 0.0 ;
  elematrix(6,4) = xc(0,2)*(-1.0/12.0) - xc(1,2)*(1.0/6.0) - xc(2,2)*(1.0/6.0) - xc(3,2)*(1.0/12.0) ;
  elematrix(6,5) = xc(0,1)*(1.0/12.0) - xc(1,1)*(1.0/6.0) + xc(3,1)*(1.0/12.0) ;
  elematrix(6,6) = 0.0 ;
  elematrix(6,7) = 0.0 ;
  elematrix(6,8) = xc(1,1)*(-1.0/6.0) + xc(3,1)*(1.0/6.0) ;
  elematrix(6,9) = 0.0 ;
  elematrix(6,10) = xc(0,2)*(1.0/12.0) + xc(1,2)*(1.0/12.0) + xc(2,2)*(1.0/6.0) + xc(3,2)*(1.0/6.0) ;
  elematrix(6,11) = xc(0,1)*(-1.0/12.0) - xc(1,1)*(1.0/12.0) + xc(3,1)*(1.0/6.0) ;

  elematrix(7,0) = xc(1,2)*(-1.0/12.0) + xc(3,2)*(1.0/12.0) ;
  elematrix(7,1) = 0.0 ;
  elematrix(7,2) = xc(1,0)*(1.0/12.0) - xc(3,0)*(1.0/12.0) ;
  elematrix(7,3) = xc(0,2)*(1.0/12.0) + xc(1,2)*(1.0/6.0) + xc(2,2)*(1.0/6.0) + xc(3,2)*(1.0/12.0) ;
  elematrix(7,4) = 0.0 ;
  elematrix(7,5) = xc(0,0)*(-1.0/12.0) + xc(1,0)*(1.0/6.0) - xc(3,0)*(1.0/12.0) ;
  elematrix(7,6) = 0.0 ;
  elematrix(7,7) = 0.0 ;
  elematrix(7,8) = xc(1,0)*(1.0/6.0) - xc(3,0)*(1.0/6.0) ;
  elematrix(7,9) = xc(0,2)*(-1.0/12.0) - xc(1,2)*(1.0/12.0) - xc(2,2)*(1.0/6.0) - xc(3,2)*(1.0/6.0) ;
  elematrix(7,10) = 0.0 ;
  elematrix(7,11) = xc(0,0)*(1.0/12.0) + xc(1,0)*(1.0/12.0) - xc(3,0)*(1.0/6.0) ;

  elematrix(8,0) = xc(1,1)*(1.0/12.0) - xc(3,1)*(1.0/12.0) ;
  elematrix(8,1) = xc(1,0)*(-1.0/12.0) + xc(3,0)*(1.0/12.0) ;
  elematrix(8,2) = 0.0 ;
  elematrix(8,3) = xc(0,1)*(-1.0/12.0) + xc(2,1)*(1.0/6.0) - xc(3,1)*(1.0/12.0) ;
  elematrix(8,4) = xc(0,0)*(1.0/12.0) - xc(2,0)*(1.0/6.0) + xc(3,0)*(1.0/12.0) ;
  elematrix(8,5) = 0.0 ;
  elematrix(8,6) = xc(1,1)*(-1.0/6.0) + xc(3,1)*(1.0/6.0) ;
  elematrix(8,7) = xc(1,0)*(1.0/6.0) - xc(3,0)*(1.0/6.0) ;
  elematrix(8,8) = 0.0 ;
  elematrix(8,9) = xc(0,1)*(1.0/12.0) + xc(1,1)*(1.0/12.0) - xc(2,1)*(1.0/6.0) ;
  elematrix(8,10) = xc(0,0)*(-1.0/12.0) - xc(1,0)*(1.0/12.0) + xc(2,0)*(1.0/6.0) ;
  elematrix(8,11) = 0.0 ;

  elematrix(9,0) = 0.0 ;
  elematrix(9,1) = xc(0,2)*(1.0/6.0) + xc(1,2)*(1.0/12.0) + xc(2,2)*(1.0/12.0) + xc(3,2)*(1.0/6.0) ;
  elematrix(9,2) = xc(0,1)*(1.0/6.0) - xc(1,1)*(1.0/12.0) - xc(2,1)*(1.0/12.0) ;
  elematrix(9,3) = 0.0 ;
  elematrix(9,4) = xc(0,2)*(-1.0/12.0) + xc(2,2)*(1.0/12.0) ;
  elematrix(9,5) = xc(0,1)*(1.0/12.0) - xc(2,1)*(1.0/12.0) ;
  elematrix(9,6) = 0.0 ;
  elematrix(9,7) = xc(0,2)*(-1.0/12.0) - xc(1,2)*(1.0/12.0) - xc(2,2)*(1.0/6.0) - xc(3,2)*(1.0/6.0) ;
  elematrix(9,8) = xc(0,1)*(1.0/12.0) + xc(1,1)*(1.0/12.0) - xc(2,1)*(1.0/6.0) ;
  elematrix(9,9) = 0.0 ;
  elematrix(9,10) = 0.0 ;
  elematrix(9,11) = xc(0,1)*(1.0/6.0) - xc(2,1)*(1.0/6.0) ;

  elematrix(10,0) = xc(0,2)*(-1.0/6.0) - xc(1,2)*(1.0/12.0) - xc(2,2)*(1.0/12.0) - xc(3,2)*(1.0/6.0) ;
  elematrix(10,1) = 0.0 ;
  elematrix(10,2) = xc(0,0)*(-1.0/6.0) + xc(1,0)*(1.0/12.0) + xc(2,0)*(1.0/12.0) ;
  elematrix(10,3) = xc(0,2)*(1.0/12.0) - xc(2,2)*(1.0/12.0) ;
  elematrix(10,4) = 0.0 ;
  elematrix(10,5) = xc(0,0)*(-1.0/12.0) + xc(2,0)*(1.0/12.0) ;
  elematrix(10,6) = xc(0,2)*(1.0/12.0) + xc(1,2)*(1.0/12.0) + xc(2,2)*(1.0/6.0) + xc(3,2)*(1.0/6.0) ;
  elematrix(10,7) = 0.0 ;
  elematrix(10,8) = xc(0,0)*(-1.0/12.0) - xc(1,0)*(1.0/12.0) + xc(2,0)*(1.0/6.0) ;
  elematrix(10,9) = 0.0 ;
  elematrix(10,10) = 0.0 ;
  elematrix(10,11) = xc(0,0)*(-1.0/6.0) + xc(2,0)*(1.0/6.0) ;

  elematrix(11,0) = xc(1,1)*(1.0/12.0) + xc(2,1)*(1.0/12.0) - xc(3,1)*(1.0/6.0) ;
  elematrix(11,1) = xc(1,0)*(-1.0/12.0) - xc(2,0)*(1.0/12.0) + xc(3,0)*(1.0/6.0) ;
  elematrix(11,2) = 0.0 ;
  elematrix(11,3) = xc(0,1)*(-1.0/12.0) + xc(2,1)*(1.0/12.0) ;
  elematrix(11,4) = xc(0,0)*(1.0/12.0) - xc(2,0)*(1.0/12.0) ;
  elematrix(11,5) = 0.0 ;
  elematrix(11,6) = xc(0,1)*(-1.0/12.0) - xc(1,1)*(1.0/12.0) + xc(3,1)*(1.0/6.0) ;
  elematrix(11,7) = xc(0,0)*(1.0/12.0) + xc(1,0)*(1.0/12.0) - xc(3,0)*(1.0/6.0) ;
  elematrix(11,8) = 0.0 ;
  elematrix(11,9) = xc(0,1)*(1.0/6.0) - xc(2,1)*(1.0/6.0) ;
  elematrix(11,10) = xc(0,0)*(-1.0/6.0) + xc(2,0)*(1.0/6.0) ;
  elematrix(11,11) = 0.0 ;
  return;
}

/*----------------------------------------------------------------------*
 * Compute first derivatives of volume                          tk 10/07*
 * with respect to the displacements                                    *
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::ComputeVolConstrDeriv(const LINALG::SerialDenseMatrix& xc,
    Epetra_SerialDenseVector& elevector)
{

  //implementation based on symbolic calculation with mupad
  elevector[0] = (1.0/6.0)*xc(1,1)*xc(0,2) + (1.0/6.0)*xc(1,1)*xc(1,2) + (1.0/12.0)*xc(1,1)*xc(2,2) -
    (1.0/6.0)*xc(0,2)*xc(3,1) - (1.0/12.0)*xc(2,1)*xc(1,2) + (1.0/12.0)*xc(1,1)*xc(3,2) -
    (1.0/12.0)*xc(1,2)*xc(3,1) + (1.0/12.0)*xc(2,1)*xc(3,2) - (1.0/12.0)*xc(3,1)*xc(2,2) -
    (1.0/6.0)*xc(3,1)*xc(3,2) ;
  elevector[1] = (-1.0/6.0)*xc(1,0)*xc(0,2) - (1.0/6.0)*xc(1,0)*xc(1,2) - (1.0/12.0)*xc(1,0)*xc(2,2) +
    (1.0/12.0)*xc(2,0)*xc(1,2) + (1.0/6.0)*xc(0,2)*xc(3,0) - (1.0/12.0)*xc(1,0)*xc(3,2) +
    (1.0/12.0)*xc(3,0)*xc(1,2) - (1.0/12.0)*xc(2,0)*xc(3,2) + (1.0/12.0)*xc(3,0)*xc(2,2) +
    (1.0/6.0)*xc(3,0)*xc(3,2) ;
  elevector[2] = (1.0/6.0)*xc(0,0)*xc(1,1) - (1.0/6.0)*xc(1,0)*xc(0,1) - (1.0/6.0)*xc(0,0)*xc(3,1) +
    (1.0/12.0)*xc(1,0)*xc(2,1) + (1.0/6.0)*xc(0,1)*xc(3,0) - (1.0/12.0)*xc(2,0)*xc(1,1) +
    (1.0/12.0)*xc(1,0)*xc(3,1) - (1.0/12.0)*xc(1,1)*xc(3,0) + (1.0/12.0)*xc(2,0)*xc(3,1) -
    (1.0/12.0)*xc(3,0)*xc(2,1) ;
  elevector[3] = (-1.0/6.0)*xc(0,1)*xc(0,2) - (1.0/6.0)*xc(0,1)*xc(1,2) - (1.0/12.0)*xc(0,1)*xc(2,2) +
    (1.0/12.0)*xc(0,2)*xc(2,1) - (1.0/12.0)*xc(0,1)*xc(3,2) + (1.0/12.0)*xc(0,2)*xc(3,1) +
    (1.0/6.0)*xc(2,1)*xc(1,2) + (1.0/6.0)*xc(2,1)*xc(2,2) + (1.0/12.0)*xc(2,1)*xc(3,2) -
    (1.0/12.0)*xc(3,1)*xc(2,2) ;
  elevector[4] = (1.0/6.0)*xc(0,0)*xc(0,2) + (1.0/6.0)*xc(0,0)*xc(1,2) + (1.0/12.0)*xc(0,0)*xc(2,2) -
    (1.0/12.0)*xc(2,0)*xc(0,2) + (1.0/12.0)*xc(0,0)*xc(3,2) - (1.0/6.0)*xc(2,0)*xc(1,2) -
    (1.0/12.0)*xc(0,2)*xc(3,0) - (1.0/6.0)*xc(2,0)*xc(2,2) - (1.0/12.0)*xc(2,0)*xc(3,2) +
    (1.0/12.0)*xc(3,0)*xc(2,2) ;
  elevector[5] = (1.0/6.0)*xc(0,0)*xc(1,1) - (1.0/6.0)*xc(1,0)*xc(0,1) - (1.0/12.0)*xc(0,0)*xc(2,1) +
    (1.0/12.0)*xc(0,1)*xc(2,0) - (1.0/12.0)*xc(0,0)*xc(3,1) + (1.0/6.0)*xc(1,0)*xc(2,1) +
    (1.0/12.0)*xc(0,1)*xc(3,0) - (1.0/6.0)*xc(2,0)*xc(1,1) + (1.0/12.0)*xc(2,0)*xc(3,1) -
    (1.0/12.0)*xc(3,0)*xc(2,1) ;
  elevector[6] = (1.0/12.0)*xc(0,1)*xc(1,2) - (1.0/12.0)*xc(1,1)*xc(0,2) - (1.0/6.0)*xc(1,1)*xc(1,2) -
    (1.0/12.0)*xc(0,1)*xc(3,2) - (1.0/6.0)*xc(1,1)*xc(2,2) + (1.0/12.0)*xc(0,2)*xc(3,1) -
    (1.0/12.0)*xc(1,1)*xc(3,2) + (1.0/12.0)*xc(1,2)*xc(3,1) + (1.0/6.0)*xc(3,1)*xc(2,2) +
    (1.0/6.0)*xc(3,1)*xc(3,2) ;
  elevector[7] = (-1.0/12.0)*xc(0,0)*xc(1,2) + (1.0/12.0)*xc(1,0)*xc(0,2) + (1.0/6.0)*xc(1,0)*xc(1,2) +
    (1.0/12.0)*xc(0,0)*xc(3,2) + (1.0/6.0)*xc(1,0)*xc(2,2) - (1.0/12.0)*xc(0,2)*xc(3,0) +
    (1.0/12.0)*xc(1,0)*xc(3,2) - (1.0/12.0)*xc(3,0)*xc(1,2) - (1.0/6.0)*xc(3,0)*xc(2,2) -
    (1.0/6.0)*xc(3,0)*xc(3,2) ;
  elevector[8] = (1.0/12.0)*xc(0,0)*xc(1,1) - (1.0/12.0)*xc(1,0)*xc(0,1) - (1.0/12.0)*xc(0,0)*xc(3,1) +
    (1.0/6.0)*xc(1,0)*xc(2,1) + (1.0/12.0)*xc(0,1)*xc(3,0) - (1.0/6.0)*xc(2,0)*xc(1,1) -
    (1.0/12.0)*xc(1,0)*xc(3,1) + (1.0/12.0)*xc(1,1)*xc(3,0) + (1.0/6.0)*xc(2,0)*xc(3,1) -
    (1.0/6.0)*xc(3,0)*xc(2,1) ;
  elevector[9] = (1.0/6.0)*xc(0,1)*xc(0,2) + (1.0/12.0)*xc(0,1)*xc(1,2) - (1.0/12.0)*xc(1,1)*xc(0,2) +
    (1.0/12.0)*xc(0,1)*xc(2,2) - (1.0/12.0)*xc(0,2)*xc(2,1) + (1.0/6.0)*xc(0,1)*xc(3,2) +
    (1.0/12.0)*xc(1,1)*xc(2,2) - (1.0/12.0)*xc(2,1)*xc(1,2) - (1.0/6.0)*xc(2,1)*xc(2,2) -
    (1.0/6.0)*xc(2,1)*xc(3,2) ;
  elevector[10] = (-1.0/6.0)*xc(0,0)*xc(0,2) - (1.0/12.0)*xc(0,0)*xc(1,2) + (1.0/12.0)*xc(1,0)*xc(0,2) -
    (1.0/12.0)*xc(0,0)*xc(2,2) + (1.0/12.0)*xc(2,0)*xc(0,2) - (1.0/6.0)*xc(0,0)*xc(3,2) -
    (1.0/12.0)*xc(1,0)*xc(2,2) + (1.0/12.0)*xc(2,0)*xc(1,2) + (1.0/6.0)*xc(2,0)*xc(2,2) +
    (1.0/6.0)*xc(2,0)*xc(3,2) ;
  elevector[11] = (1.0/12.0)*xc(0,0)*xc(1,1) - (1.0/12.0)*xc(1,0)*xc(0,1) + (1.0/12.0)*xc(0,0)*xc(2,1) -
    (1.0/12.0)*xc(0,1)*xc(2,0) - (1.0/6.0)*xc(0,0)*xc(3,1) + (1.0/12.0)*xc(1,0)*xc(2,1) +
    (1.0/6.0)*xc(0,1)*xc(3,0) - (1.0/12.0)*xc(2,0)*xc(1,1) + (1.0/6.0)*xc(2,0)*xc(3,1) -
    (1.0/6.0)*xc(3,0)*xc(2,1) ;
  return;
}

/*----------------------------------------------------------------------*
 * Compute surface area and its first and second derivatives    lw 05/08*
 * with respect to the displacements                                    *
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::ComputeAreaDeriv(const LINALG::SerialDenseMatrix& x,
                                                        const int numnode,
                                                        const int ndof,
                                                        double& A,
                                                        RCP<Epetra_SerialDenseVector> Adiff,
                                                        RCP<Epetra_SerialDenseMatrix> Adiff2)
{
  const DRT::UTILS::IntegrationPoints2D  intpoints =
    getIntegrationPoints2D(gaussrule_);

  int ngp = intpoints.nquad;

  // allocate vector for shape functions and matrix for derivatives
  LINALG::SerialDenseMatrix  deriv(2,numnode);
  LINALG::SerialDenseMatrix  dxyzdrs(2,3);

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/

  for (int gpid = 0; gpid < ngp; ++gpid)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get derivatives of shape functions in the plane of the element
    DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());

    vector<double> normal(3);
    double detA;
    double Jac;
    SurfaceIntegration(detA,normal,x,deriv);
    Jac = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
    A += Jac*intpoints.qwgt[gpid];

    blitz::Array<double,2> ddet(3,ndof);
    blitz::Array<double,3> ddet2(3,ndof,ndof);
    ddet2 = 0.;
    blitz::Array<double,1> jacobi_deriv(ndof);

    dxyzdrs.Multiply('N','N',1.0,deriv,x,0.0);

    /*--------------- derivation of minor determiants of the Jacobian
     *----------------------------- with respect to the displacements */
    for (int i=0;i<numnode;++i)
    {
      ddet(0,3*i)   = 0.;
      ddet(0,3*i+1) = deriv(0,i)*dxyzdrs(1,2)-deriv(1,i)*dxyzdrs(0,2);
      ddet(0,3*i+2) = deriv(1,i)*dxyzdrs(0,1)-deriv(0,i)*dxyzdrs(1,1);

      ddet(1,3*i)   = deriv(1,i)*dxyzdrs(0,2)-deriv(0,i)*dxyzdrs(1,2);
      ddet(1,3*i+1) = 0.;
      ddet(1,3*i+2) = deriv(0,i)*dxyzdrs(1,0)-deriv(1,i)*dxyzdrs(0,0);

      ddet(2,3*i)   = deriv(0,i)*dxyzdrs(1,1)-deriv(1,i)*dxyzdrs(0,1);
      ddet(2,3*i+1) = deriv(1,i)*dxyzdrs(0,0)-deriv(0,i)*dxyzdrs(1,0);
      ddet(2,3*i+2) = 0.;

      jacobi_deriv(i*3)   = 1/Jac*(normal[2]*ddet(2,3*i  )+normal[1]*ddet(1,3*i  ));
      jacobi_deriv(i*3+1) = 1/Jac*(normal[2]*ddet(2,3*i+1)+normal[0]*ddet(0,3*i+1));
      jacobi_deriv(i*3+2) = 1/Jac*(normal[0]*ddet(0,3*i+2)+normal[1]*ddet(1,3*i+2));
    }

    /*--- calculation of first derivatives of current interfacial area
     *----------------------------- with respect to the displacements */
    for (int i=0;i<ndof;++i)
    {
      (*Adiff)[i] += jacobi_deriv(i)*intpoints.qwgt[gpid];
    }

    if (Adiff2!=null)
    {
      /*--------- second derivates of minor determiants of the Jacobian
       *----------------------------- with respect to the displacements */
      for (int n=0;n<numnode;++n)
      {
        for (int o=0;o<numnode;++o)
        {
          ddet2(0,n*3+1,o*3+2) = deriv(0,n)*deriv(1,o)-deriv(1,n)*deriv(0,o);
          ddet2(0,n*3+2,o*3+1) = - ddet2(0,n*3+1,o*3+2);

          ddet2(1,n*3  ,o*3+2) = deriv(1,n)*deriv(0,o)-deriv(0,n)*deriv(1,o);
          ddet2(1,n*3+2,o*3  ) = - ddet2(1,n*3,o*3+2);

          ddet2(2,n*3  ,o*3+1) = ddet2(0,n*3+1,o*3+2);
          ddet2(2,n*3+1,o*3  ) = - ddet2(2,n*3,o*3+1);
        }
      }

      /*- calculation of second derivatives of current interfacial areas
       *----------------------------- with respect to the displacements */
      for (int i=0;i<ndof;++i)
      {
        int var1, var2;

        if (i%3==0)           // displacement in x-direction
        {
          var1 = 1;
          var2 = 2;
        }
        else if ((i-1)%3==0)  // displacement in y-direction
        {
          var1 = 0;
          var2 = 2;
        }
        else if ((i-2)%3==0)  // displacement in z-direction
        {
          var1 = 0;
          var2 = 1;
        }
        else
        {
          dserror("calculation of second derivatives of interfacial area failed");
          exit(1);
        }

        for (int j=0;j<ndof;++j)
        {
          (*Adiff2)(i,j) += (-1/Jac*jacobi_deriv(j)*jacobi_deriv(i)+1/Jac*
                             (ddet(var1,i)*ddet(var1,j)+normal[var1]*ddet2(var1,i,j)+
                              ddet(var2,i)*ddet(var2,j)+normal[var2]*ddet2(var2,i,j)))*intpoints.qwgt[gpid];
        }
      }
    }
  }

  return;
}

#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOLID3
