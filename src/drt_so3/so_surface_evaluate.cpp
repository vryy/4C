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
//  cout << "SoSurface LM=" << lm.size() << " : ";
//  for (int i=0; i<(int)lm.size(); ++i) cout << lm[i] << " ";
//  cout << endl;

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
  const int numdf = 3;
  const int actnumdf = ActualNumDofPerNode(lm.size());
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
  const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
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
          elevec1[node*actnumdf+dim]+= funct[node] * (*onoff)[dim] * (*val)[dim] * fac;
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
          elevec1[node*actnumdf+dim] += funct[node] * normal[dim] * fac;
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
          double volumeele = ComputeConstrVols(xscurr,NumNode());
          elevector3[0]= volumeele;
        }
      }
      break;
      case calc_struct_volconstrstiff:
      {
        // element geometry update
        RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp==null) dserror("Cannot get state vector 'displacement'");
        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
        const int numdim =3;
        LINALG::SerialDenseMatrix xscurr(NumNode(),numdim);  // material coord. of element
        SpatialConfiguration(xscurr,mydisp);
        double volumeele;
        // first partial derivatives
        RCP<Epetra_SerialDenseVector> Vdiff1 = rcp(new Epetra_SerialDenseVector);
        // second partial derivatives
        RCP<Epetra_SerialDenseMatrix> Vdiff2 = rcp(new Epetra_SerialDenseMatrix);

        //call submethod to compute volume and its derivatives w.r.t. to current displ.
        ComputeVolDeriv(xscurr, NumNode(),numdim*NumNode(), volumeele, Vdiff1, Vdiff2);
        //update rhs vector and corresponding column in "constraint" matrix
        elevector1 = *Vdiff1;
        elevector2 = *Vdiff1;
        elematrix1 = *Vdiff2;
        //call submethod for volume evaluation and store result in third systemvector
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
        const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);

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
        bool newstep = params.get<bool>("newstep", false);

        // element geometry update

        const int numnode = NumNode();
        LINALG::SerialDenseMatrix x(numnode,3);

        RefCountPtr<const Epetra_Vector> dism = discretization.GetState("displacement");
        if (dism==null) dserror("Cannot get state vector 'displacement'");
        vector<double> mydism(lm.size());
        DRT::UTILS::ExtractMyValues(*dism,mydism,lm);
        SpatialConfiguration(x,mydism);

        const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);

        // set up matrices and parameters needed for the evaluation of current
        // interfacial area and its derivatives w.r.t. the displacements

        int ndof = 3*numnode;                                     // overall number of surface dofs
        double A;                                                 // interfacial area
        // first partial derivatives
        RCP<Epetra_SerialDenseVector> Adiff = rcp(new Epetra_SerialDenseVector);
        // second partial derivatives
        RCP<Epetra_SerialDenseMatrix> Adiff2 = rcp(new Epetra_SerialDenseMatrix);

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

          // element geometry update (n+1)
          RefCountPtr<const Epetra_Vector> disn = discretization.GetState("new displacement");
          if (disn==null) dserror("Cannot get state vector 'new displacement'");
          vector<double> mydisn(lm.size());
          DRT::UTILS::ExtractMyValues(*disn,mydisn,lm);
          SpatialConfiguration(x,mydisn);

          // set up matrices and parameters needed for the evaluation of
          // interfacial area and its first derivative w.r.t. the displacements at (n+1)
          double Anew = 0.;                                            // interfacial area
          RCP<Epetra_SerialDenseVector> Adiffnew = rcp(new Epetra_SerialDenseVector);

          ComputeAreaDeriv(x, numnode, ndof, Anew, Adiffnew, null);

          surfstressman->StiffnessAndInternalForces(curvenum, A, Adiff, Adiff2, Anew, Adiffnew, elevector1, elematrix1, this->Id(),
                                                    time, dt, 0, 0.0, k1xC, k2, m1, m2, gamma_0,
                                                    gamma_min, gamma_min_eq, con_quot_max,
                                                    con_quot_eq, newstep);
        }
        else if (cond->Type()==DRT::Condition::SurfaceTension) // ideal liquid
        {
          int curvenum = cond->Getint("curve");
          double const_gamma = cond->GetDouble("gamma");
          RCP<Epetra_SerialDenseVector> Adiffnew = rcp(new Epetra_SerialDenseVector);
          surfstressman->StiffnessAndInternalForces(curvenum, A, Adiff, Adiff2, 0., Adiffnew, elevector1, elematrix1, this->Id(),
                                                    time, dt, 1, const_gamma, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0, 0.0, newstep);
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
          const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
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
        // element geometry update
        RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp==null) dserror("Cannot get state vector 'displacement'");
        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
        const int numdim =3;
        LINALG::SerialDenseMatrix xscurr(NumNode(),numdim);  // material coord. of element
        SpatialConfiguration(xscurr,mydisp);
        // initialize variables
        double elearea;
        // first partial derivatives
        RCP<Epetra_SerialDenseVector> Adiff = rcp(new Epetra_SerialDenseVector);
        // second partial derivatives
        RCP<Epetra_SerialDenseMatrix> Adiff2 = rcp(new Epetra_SerialDenseMatrix);

        //call submethod
        ComputeAreaDeriv(xscurr, NumNode(),numdim*NumNode(), elearea, Adiff, Adiff2);
        // store result
        elevector3[0] = elearea;

      }
      break;
      case calc_struct_areaconstrstiff:
      {
        // element geometry update
        RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp==null) dserror("Cannot get state vector 'displacement'");
        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
        const int numdim =3;
        LINALG::SerialDenseMatrix xscurr(NumNode(),numdim);  // material coord. of element
        SpatialConfiguration(xscurr,mydisp);
        // initialize variables
        double elearea;
        // first partial derivatives
        RCP<Epetra_SerialDenseVector> Adiff = rcp(new Epetra_SerialDenseVector);
        // second partial derivatives
        RCP<Epetra_SerialDenseMatrix> Adiff2 = rcp(new Epetra_SerialDenseMatrix);

        //call submethod
        ComputeAreaDeriv(xscurr, NumNode(),numdim*NumNode(), elearea, Adiff, Adiff2);
        //update elematrices and elevectors
        elevector1 = *Adiff;
        elevector1.Scale(-1.0);
        elevector2 = elevector1;
        elematrix1 = *Adiff2;
        elematrix1.Scale(-1.0);
        elevector3[0] = elearea;
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
double DRT::ELEMENTS::StructuralSurface::ComputeConstrVols
(
    const LINALG::SerialDenseMatrix& xc,
    const int numnode
)
{
  double volume = 0.0;
  
  //Volume is calculated by evaluating the integral of z(r,s) over dA*,
  //where dA* is the projection of dA on the xy-plane.
  //Therefore separate current configuration between xy and z
  LINALG::SerialDenseMatrix xy= xc;
  LINALG::SerialDenseVector z (numnode);
  for (int i = 0; i < numnode; i++)
  {
    xy(i,2) = 0.0; // project by z_i = 0.0
    z(i) = xc(i,2); // extract z coordinate
  }
  
  // get gaussrule
  const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
  int ngp = intpoints.nquad;
  
  // allocate vector for shape functions and matrix for derivatives
  LINALG::SerialDenseVector  funct(numnode);
  LINALG::SerialDenseMatrix  deriv(2,numnode);

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/
  for (int gpid = 0; gpid < ngp; ++gpid)
  { 
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get shape functions and derivatives of shape functions in the plane of the element
    DRT::UTILS::shape_function_2D(funct,e0,e1,Shape());
    DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());
    
    double detA;
    // compute "metric tensor" deriv*xy, which is a 2x3 matrix with zero 3rd column 
    LINALG::SerialDenseMatrix metrictensor(2,3);
    metrictensor.Multiply('N','N',1.0,deriv,xy,0.0);
    //LINALG::SerialDenseMatrix metrictensor(2,2);
    //metrictensor.Multiply('N','T',1.0,dxyzdrs,dxyzdrs,0.0);
    detA =  metrictensor(0,0)*metrictensor(1,1)-metrictensor(0,1)*metrictensor(1,0);
    // add weighted volume at gausspoint
    volume -= funct.Dot(z)*detA*intpoints.qwgt[gpid];
    
  }  
  return volume;
}

/*----------------------------------------------------------------------*
 * Compute volume and its first and second derivatives          tk 02/09*
 * with respect to the displacements                                    *
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::ComputeVolDeriv
(
  const LINALG::SerialDenseMatrix& xc,
  const int numnode,
  const int ndof,
  double& V,
  RCP<Epetra_SerialDenseVector> Vdiff1,
  RCP<Epetra_SerialDenseMatrix> Vdiff2
)
{
  // necessary constants
  const int numdim = 3;
  
  // initialize
  V = 0.0;
  Vdiff1->Size(ndof);
  if (Vdiff2!=null) Vdiff2->Shape(ndof, ndof);
  
  //Volume is calculated by evaluating the integral of z(r,s) over dA*,
  //where dA* is the projection of dA on the xy-plane.
  //Therefore separate current configuration between xy and z
  LINALG::SerialDenseMatrix xy= xc;
  LINALG::SerialDenseVector z (numnode);
  for (int i = 0; i < numnode; i++)
  {
    xy(i,2) = 0.0; // project by z_i = 0.0
    z(i) = xc(i,2); // extract z coordinate
  }
  
  // get gaussrule
  const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
  int ngp = intpoints.nquad;
  
  // allocate vector for shape functions and matrix for derivatives
  LINALG::SerialDenseVector  funct(numnode);
  LINALG::SerialDenseMatrix  deriv(2,numnode);

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/
  for (int gpid = 0; gpid < ngp; ++gpid)
  { 
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get shape functions and derivatives of shape functions in the plane of the element
    DRT::UTILS::shape_function_2D(funct,e0,e1,Shape());
    DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());
    
    // evaluate Jacobi determinant, for projected dA*
    vector<double> normal(numdim);
    double detA;
    // compute "metric tensor" deriv*xy, which is a 2x3 matrix with zero 3rd column 
    LINALG::SerialDenseMatrix metrictensor(2,numdim);
    metrictensor.Multiply('N','N',1.0,deriv,xy,0.0);
    //metrictensor.Multiply('N','T',1.0,dxyzdrs,dxyzdrs,0.0);
    detA =  metrictensor(0,0)*metrictensor(1,1)-metrictensor(0,1)*metrictensor(1,0);
    const double dotprodz = funct.Dot(z);
    // add weighted volume at gausspoint
    V -= dotprodz*detA*intpoints.qwgt[gpid];
    
    //-------- compute first derivative
    for (int i = 0; i < numnode ; i++)
    {
      (*Vdiff1)[3*i]   += intpoints.qwgt[gpid]*dotprodz*(deriv(0,i)*metrictensor(1,1)-metrictensor(0,1)*deriv(1,i));
      (*Vdiff1)[3*i+1] += intpoints.qwgt[gpid]*dotprodz*(deriv(1,i)*metrictensor(0,0)-metrictensor(1,0)*deriv(0,i));
      (*Vdiff1)[3*i+2] += intpoints.qwgt[gpid]*funct[i]*detA;
    }
    
    //-------- compute second derivative
    if (Vdiff2!=null)
    {
      for (int i = 0; i < numnode ; i++)
      {
        for (int j = 0; j < numnode ; j++)
        {
          //"diagonal" (dV)^2/(dx_i dx_j) = 0, therefore only six entries have to be specified
          (*Vdiff2)(3*i,3*j+1) += intpoints.qwgt[gpid]*dotprodz*(deriv(0,i)*deriv(1,j)-deriv(1,i)*deriv(0,j));
          (*Vdiff2)(3*i+1,3*j) += intpoints.qwgt[gpid]*dotprodz*(deriv(0,j)*deriv(1,i)-deriv(1,j)*deriv(0,i));
          (*Vdiff2)(3*i,3*j+2) += intpoints.qwgt[gpid]*funct[j]*(deriv(0,i)*metrictensor(1,1)-metrictensor(0,1)*deriv(1,i));
          (*Vdiff2)(3*i+2,3*j) += intpoints.qwgt[gpid]*funct[i]*(deriv(0,j)*metrictensor(1,1)-metrictensor(0,1)*deriv(1,j));
          (*Vdiff2)(3*i+1,3*j+2) += intpoints.qwgt[gpid]*funct[j]*(deriv(1,i)*metrictensor(0,0)-metrictensor(1,0)*deriv(0,i));
          (*Vdiff2)(3*i+2,3*j+1) += intpoints.qwgt[gpid]*funct[i]*(deriv(1,j)*metrictensor(0,0)-metrictensor(1,0)*deriv(0,j));
        }
      }
    }
  
  }   
  
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
  // initialization
  A = 0.;
  Adiff->Size(ndof);

  if (Adiff2!=null) Adiff2->Shape(ndof, ndof);

  const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);

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
    SurfaceIntegration(detA,normal,x,deriv);
    A += detA*intpoints.qwgt[gpid];

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

      jacobi_deriv(i*3)   = 1/detA*(normal[2]*ddet(2,3*i  )+normal[1]*ddet(1,3*i  ));
      jacobi_deriv(i*3+1) = 1/detA*(normal[2]*ddet(2,3*i+1)+normal[0]*ddet(0,3*i+1));
      jacobi_deriv(i*3+2) = 1/detA*(normal[0]*ddet(0,3*i+2)+normal[1]*ddet(1,3*i+2));
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
          (*Adiff2)(i,j) += (-1/detA*jacobi_deriv(j)*jacobi_deriv(i)+1/detA*
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
