/*!----------------------------------------------------------------------
\file so_hex8_surface.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_hex8.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_surfstress/drt_surfstress_manager.H"
#include "../drt_lib/drt_container.H"

/*----------------------------------------------------------------------*
 * Integrate a Surface Neumann boundary condition (public)     maf 04/07*
 * ---------------------------------------------------------------------*/
int DRT::ELEMENTS::Soh8Surface::EvaluateNeumann(ParameterList&           params,
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
    neum_orthopressure,
    neum_consthydro_z,
    neum_increhydro_z,
    neum_live_FSI,
    neum_opres_FSI
  };
  LoadType ltype;
  const string* type = condition.Get<string>("type");
  if      (*type == "neum_live")          ltype = neum_live;
  //else if (*type == "neum_live_FSI")      ltype = neum_live_FSI;
  else if (*type == "neum_orthopressure") ltype = neum_orthopressure;
  else dserror("Unknown type of SurfaceNeumann condition");

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
  // **

  // element geometry update
  RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==null) dserror("Cannot get state vector 'displacement'");
  vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
  const int numnod = 4;
  Epetra_SerialDenseMatrix xsrefe(numnod,NUMDIM_SOH8);  // material coord. of element
  Epetra_SerialDenseMatrix xscurr(numnod,NUMDIM_SOH8);  // material coord. of element
  for (int i=0; i<numnod; ++i){
    xsrefe(i,0) = Nodes()[i]->X()[0];
    xsrefe(i,1) = Nodes()[i]->X()[1];
    xsrefe(i,2) = Nodes()[i]->X()[2];

    xscurr(i,0) = xsrefe(i,0) + mydisp[i*NODDOF_SOH8+0];
    xscurr(i,1) = xsrefe(i,1) + mydisp[i*NODDOF_SOH8+1];
    xscurr(i,2) = xsrefe(i,2) + mydisp[i*NODDOF_SOH8+2];
  }

  /*
  ** Here, we integrate a 4-node surface with 2xs(1,0) Gauss Points
  */
  const int ngp = 4;

  // gauss parameters
  const double gpweight = 1.0;
  const double gploc    = 1.0/sqrt(3.0);
  Epetra_SerialDenseMatrix gpcoord (ngp,2);
  gpcoord(0,0) = - gploc;
  gpcoord(0,1) = - gploc;
  gpcoord(1,0) =   gploc;
  gpcoord(1,1) = - gploc;
  gpcoord(2,0) = - gploc;
  gpcoord(2,1) =   gploc;
  gpcoord(3,0) =   gploc;
  gpcoord(3,1) =   gploc;

  for (int gpid = 0; gpid < 4; ++gpid) {    // loop over intergration points
    // get shape functions and derivatives of element surface
    vector<double> funct(ngp);                // 4 shape function values
    double drs;                               // surface area factor
    switch(ltype){
      case neum_live:{            // uniform load on reference configuration
        soh8_surface_integ(&funct,&drs,NULL,&xsrefe,gpcoord(gpid,0),gpcoord(gpid,1));
        double fac = gpweight * drs * curvefac;   // integration factor
        // distribute over element load vector
        for (int nodid=0; nodid < 4; ++nodid) {
          for(int dim=0; dim < NUMDIM_SOH8; ++dim) {
            elevec1[nodid*NUMDIM_SOH8 + dim] += funct[nodid] * (*onoff)[dim] * (*val)[dim] * fac;
          }
        }
      }
      break;
      case neum_orthopressure:{   // orthogonal pressure on deformed config.
        if ((*onoff)[0] != 1) dserror("orthopressure on 1st dof only!");
        for (int checkdof = 1; checkdof < 3; ++checkdof) {
          if ((*onoff)[checkdof] != 0) dserror("orthopressure on 1st dof only!");
        }
        double ortho_value = (*val)[0];
        if (!ortho_value) dserror("no orthopressure value given!");
        vector<double> unrm(NUMDIM_SOH8);
        soh8_surface_integ(&funct,&drs,&unrm,&xscurr,gpcoord(gpid,0),gpcoord(gpid,1));
        double fac = (-1.0) * gpweight * drs * curvefac;   // integration factor
        // distribute over element load vector
        for (int nodid=0; nodid < 4; ++nodid) {
          for(int dim=0; dim < NUMDIM_SOH8; ++dim) {
            elevec1[nodid*NUMDIM_SOH8 + dim] +=
                                funct[nodid] * ortho_value * unrm[dim] * fac;
          }
        }
      }
      break;
    default:
      dserror("Unknown type of SurfaceNeumann load");
    break;
    }

  }

//    cout << elevec1 << endl;
  return 0;
}

/*----------------------------------------------------------------------*
 * Evaluate sqrt of determinant of metric at gp (private)      maf 05/07*
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::Soh8Surface::soh8_surface_integ(
      vector<double>* funct,                 // (o) shape functions
      double* sqrtdetg,                      // (o) pointer to sqrt of det(g)
      vector<double>* unrm,                  // (o) unit normal
      const Epetra_SerialDenseMatrix* xs,    // (i) element coords
      const double r,                        // (i) coord in r-direction
      const double s)                        // (i) coord in s-direction
{
  // shape functions for 4 nodes
  (*funct)[0] = 0.25 * (1.0-r) * (1.0-s);
  (*funct)[1] = 0.25 * (1.0+r) * (1.0-s);
  (*funct)[2] = 0.25 * (1.0+r) * (1.0+s);
  (*funct)[3] = 0.25 * (1.0-r) * (1.0+s);
  // derivatives of 4 shape functions wrt 2 directions
  Epetra_SerialDenseMatrix deriv(4,2);
  deriv(0,0) = -0.25 * (1.0-s);
  deriv(0,1) = -0.25 * (1.0-r);
  deriv(1,0) =  0.25 * (1.0-s);
  deriv(1,1) = -0.25 * (1.0+r);
  deriv(2,0) =  0.25 * (1.0+s);
  deriv(2,1) =  0.25 * (1.0+r);
  deriv(3,0) = -0.25 * (1.0+s);
  deriv(3,1) =  0.25 * (1.0-r);

  // compute dXYZ / drs
  Epetra_SerialDenseMatrix dxyzdrs(2,3);
  dxyzdrs.Multiply('T','N',1.0,deriv,(*xs),1.0);

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
  Epetra_SerialDenseMatrix metrictensor(2,2);
  metrictensor.Multiply('N','T',1.0,dxyzdrs,dxyzdrs,1.0);
  (*sqrtdetg) = sqrt( metrictensor(0,0)*metrictensor(1,1)
                     -metrictensor(0,1)*metrictensor(1,0));
  if (unrm != NULL){
    (*unrm)[0] = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
    (*unrm)[1] = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
    (*unrm)[2] = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);
  }

  return;
}

/*----------------------------------------------------------------------*
 * Evaluate method for Soh8Surface-Elements                     tk 10/07*
 * ---------------------------------------------------------------------*/
int DRT::ELEMENTS::Soh8Surface::Evaluate(ParameterList& params,
        DRT::Discretization&      discretization,
        vector<int>&              lm,
        Epetra_SerialDenseMatrix& elematrix1,
        Epetra_SerialDenseMatrix& elematrix2,
        Epetra_SerialDenseVector& elevector1,
        Epetra_SerialDenseVector& elevector2,
        Epetra_SerialDenseVector& elevector3)
{
  const DiscretizationType distype = this->Shape();

  // start with "none"
  DRT::ELEMENTS::Soh8Surface::ActionType act = Soh8Surface::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_constrvol")       act = Soh8Surface::calc_struct_constrvol;
  else if (action=="calc_struct_volconstrstiff")  act= Soh8Surface::calc_struct_volconstrstiff;
  else if (action=="calc_struct_monitarea")      act=Soh8Surface::calc_struct_monitarea;
  else if (action=="calc_struct_areaconstrstiff") act=Soh8Surface::calc_struct_areaconstrstiff;
  else if (action=="calc_init_vol")               act= Soh8Surface::calc_init_vol;
  else if (action=="calc_surfstress_stiff")       act= Soh8Surface::calc_surfstress_stiff;
  else dserror("Unknown type of action for Soh8Surface");
  //create communicator
  const Epetra_Comm& Comm = discretization.Comm();
  // what the element has to do
  switch(act)
  {
      //just compute the enclosed volume (e.g. for initialization)
      case calc_struct_constrvol:
      {
        if (distype!=quad4)
          dserror("Volume Constraint online works for quad4 surfaces!");
        //We are not interested in volume of ghosted elements
        if(Comm.MyPID()==this->Owner())
        {
          // element geometry update
          RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
          if (disp==null) dserror("Cannot get state vector 'displacement'");
          vector<double> mydisp(lm.size());
          DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
          const int numnod = 4;
          Epetra_SerialDenseMatrix xsrefe(numnod,NUMDIM_SOH8);  // material coord. of element
          Epetra_SerialDenseMatrix xscurr(numnod,NUMDIM_SOH8);  // material coord. of element
          for (int i=0; i<numnod; ++i)
          {
            xsrefe(i,0) = Nodes()[i]->X()[0];
            xsrefe(i,1) = Nodes()[i]->X()[1];
            xsrefe(i,2) = Nodes()[i]->X()[2];

            xscurr(i,0) = xsrefe(i,0) + mydisp[i*NODDOF_SOH8+0];
            xscurr(i,1) = xsrefe(i,1) + mydisp[i*NODDOF_SOH8+1];
            xscurr(i,2) = xsrefe(i,2) + mydisp[i*NODDOF_SOH8+2];
          }
          //call submethod
          double volumeele =  ComputeConstrVols(xscurr);

          // get RIGHT volume out of parameterlist and maximum ConditionID
          char volname[30];
          const int ID =params.get("ConditionID",-1);
          const int maxID=params.get("MaxID",0);
          const int minID=params.get("MinID",1000000);
          if (ID<0)
          {
            dserror("Condition ID for volume constraint missing!");
          }
          if (maxID<ID)
          {
            params.set("MaxID",ID);
          }
          if (minID>ID)
          {
            params.set("MinID",ID);
          }
          sprintf(volname,"computed volume %d",ID);
          double volumecond = params.get(volname,0.0);
          //update volume in parameter list
          params.set(volname, volumecond+volumeele);
        }

      }
      break;
      case calc_struct_volconstrstiff:
      {
        if (distype!=quad4)
          dserror("Volume Constraint online works for quad4 surfaces!");
        // element geometry update
        RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp==null) dserror("Cannot get state vector 'displacement'");
        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
        const int numnod = 4;
        Epetra_SerialDenseMatrix xsrefe(numnod,NUMDIM_SOH8);  // material coord. of element
        Epetra_SerialDenseMatrix xscurr(numnod,NUMDIM_SOH8);  // material coord. of element
        for (int i=0; i<numnod; ++i)
        {
          xsrefe(i,0) = Nodes()[i]->X()[0];
          xsrefe(i,1) = Nodes()[i]->X()[1];
          xsrefe(i,2) = Nodes()[i]->X()[2];

          xscurr(i,0) = xsrefe(i,0) + mydisp[i*NODDOF_SOH8+0];
          xscurr(i,1) = xsrefe(i,1) + mydisp[i*NODDOF_SOH8+1];
          xscurr(i,2) = xsrefe(i,2) + mydisp[i*NODDOF_SOH8+2];
        }
        //call submethods
        ComputeVolconstrStiff(xscurr,elematrix1);
        ComputeVolconstrVolDeriv(xscurr,elevector1);
        //apply the right lagrange multiplier and right signs to matrix and vectors
        const int ID =params.get("ConditionID",-1);
        int numID=params.get("NumberofID",0);
        RCP<Epetra_Map> reducedmap = rcp(new Epetra_Map(*(params.get<RCP<Epetra_Map> >("ReducedMap",reducedmap))));
        RCP<Epetra_Vector> lambdav=rcp(new Epetra_Vector(*reducedmap));
        lambdav= params.get("LagrMultVector",lambdav);
        if (ID<0)
        {
          dserror("Condition ID for volume constraint missing!");
        }
        const int minID =params.get("MinID",0);
        //update corresponding column in "constraint" matrix
        elevector2=elevector1;
        elevector1.Scale(1*(*lambdav)[ID-minID]);
        elematrix1.Scale(1*(*lambdav)[ID-minID]);
        //call submethod for volume evaluation
        if(Comm.MyPID()==this->Owner())
        {
          double volumeele = ComputeConstrVols(xscurr);
          // get RIGHT volume out of parameterlist and maximum ConditionID
          char volname[30];
          sprintf(volname,"computed volume %d",ID);
          double volumecond = params.get(volname,0.0);
          //update volume in parameter list
          params.set(volname, volumecond+volumeele);
        }
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

        double V = params.get<double>("V0", 0.0);
        double dV = 0.0;
        const int numnod = 4;
        Epetra_SerialDenseMatrix xsrefe(numnod,NUMDIM_SOH8);  // material coord. of element
        Epetra_SerialDenseMatrix xscurr(numnod,NUMDIM_SOH8);  // material coord. of element
        for (int i=0; i<numnod; ++i){
          xsrefe(i,0) = Nodes()[i]->X()[0];
          xsrefe(i,1) = Nodes()[i]->X()[1];
          xsrefe(i,2) = Nodes()[i]->X()[2];
        }
        const int ngp = 4;

        // gauss parameters
        const double gpweight = 1.0;
        const double gploc    = 1.0/sqrt(3.0);
        Epetra_SerialDenseMatrix gpcoord (ngp,2);
        gpcoord(0,0) = - gploc;
        gpcoord(0,1) = - gploc;
        gpcoord(1,0) =   gploc;
        gpcoord(1,1) = - gploc;
        gpcoord(2,0) = - gploc;
        gpcoord(2,1) =   gploc;
        gpcoord(3,0) =   gploc;
        gpcoord(3,1) =   gploc;

        for (int gpid = 0; gpid < 4; ++gpid) {    // loop over integration points
          vector<double> funct(ngp);              // 4 shape function values
          double drs;                             // surface area factor
          vector<double> normal(NUMDIM_SOH8);
          double temp = 0.0;
          vector<double> X(NUMDIM_SOH8);

          soh8_surface_integ(&funct,&drs,&normal,&xsrefe,gpcoord(gpid,0),gpcoord(gpid,1));

          X[0] = funct[0]*xsrefe(0,0) + funct[1]*xsrefe(1,0) + funct[2]*xsrefe(2,0) + funct[3]*xsrefe(3,0);
          X[1] = funct[0]*xsrefe(0,1) + funct[1]*xsrefe(1,1) + funct[2]*xsrefe(2,1) + funct[3]*xsrefe(3,1);
          X[2] = funct[0]*xsrefe(0,2) + funct[1]*xsrefe(1,2) + funct[2]*xsrefe(2,2) + funct[3]*xsrefe(3,2);

          for (int i=0;i<3;++i){
            temp += normal[i]*normal[i];
          }
          double absnorm = sqrt(temp);
          for (int i=0;i<3;++i){
            normal[i] /= absnorm;
          }
          double fac = gpweight * drs;
          for (int i=0;i<3;++i){
            dV += 1/3.0*fac*normal[i]*X[i];
          }
        }
        params.set("V0", V+dV);
      }
      break;

      case calc_surfstress_stiff:
      {
        // element geometry update

        RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp==null) dserror("Cannot get state vector 'displacement'");
        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
        const int numnod = 4;
        Epetra_SerialDenseMatrix xsrefe(numnod,NUMDIM_SOH8);  // material coord. of element
        Epetra_SerialDenseMatrix xscurr(numnod,NUMDIM_SOH8);  // material coord. of element

        for (int i=0; i<numnod; ++i){
          xsrefe(i,0) = Nodes()[i]->X()[0];
          xsrefe(i,1) = Nodes()[i]->X()[1];
          xsrefe(i,2) = Nodes()[i]->X()[2];

          xscurr(i,0) = xsrefe(i,0) + mydisp[i*NODDOF_SOH8+0];
          xscurr(i,1) = xsrefe(i,1) + mydisp[i*NODDOF_SOH8+1];
          xscurr(i,2) = xsrefe(i,2) + mydisp[i*NODDOF_SOH8+2];
        }

        RefCountPtr<SurfStressManager> surfstressman =
          params.get<RefCountPtr<SurfStressManager> >("surfstr_man", null);

        if (surfstressman==null)
          dserror("No SurfStressManager in SoHex8Surface available");

        double time = params.get<double>("total time",-1.0);
        double dt = params.get<double>("delta time",0.0);
        RefCountPtr<DRT::Condition> cond = params.get<RefCountPtr<DRT::Condition> >("condition",null);
        if (cond==null)
          dserror("Condition not available in SoHex8Surface");

        // set up matrices and parameters needed for the evaluation of current
        // interfacial area and its derivatives w.r.t. the displacements

        int ngp=4;                                // number of Gauss points
        int nnode=4;                              // number of surface nodes
        double A;                                 // interfacial area
        Epetra_SerialDenseVector Adiff(12);       // first partial derivatives
        Epetra_SerialDenseMatrix Adiff2(12,12);   // second partial derivatives
        Epetra_SerialDenseVector gpweight(ngp);   // Gauss point weights
        for (int i=0;i<ngp;++i)
          gpweight[i]=1.;
        Epetra_SerialDenseMatrix gpcoord(ngp,2);  // Gauss point coordinates
        const double gploc    = 1.0/sqrt(3.0);
        gpcoord(0,0) = - gploc;
        gpcoord(0,1) = - gploc;
        gpcoord(1,0) =   gploc;
        gpcoord(1,1) = - gploc;
        gpcoord(2,0) = - gploc;
        gpcoord(2,1) =   gploc;
        gpcoord(3,0) =   gploc;
        gpcoord(3,1) =   gploc;
        vector<Epetra_SerialDenseMatrix> deriv(ngp);      // derivatives of shape functions
        vector<Epetra_SerialDenseMatrix> dxyzdrs(ngp);    // dXYZ / drs
        for (int gpid = 0; gpid < 4; ++gpid)
        {
          double r = gpcoord(gpid,0);
          double s = gpcoord(gpid,1);
          Epetra_SerialDenseMatrix der(4,2);
          der(0,0) = -0.25 * (1.0-s);
          der(0,1) = -0.25 * (1.0-r);
          der(1,0) =  0.25 * (1.0-s);
          der(1,1) = -0.25 * (1.0+r);
          der(2,0) =  0.25 * (1.0+s);
          der(2,1) =  0.25 * (1.0+r);
          der(3,0) = -0.25 * (1.0+s);
          der(3,1) =  0.25 * (1.0-r);
          Epetra_SerialDenseMatrix dx(2,3);
          dx.Multiply('T','N',1.0,der,xscurr,0.0);
          deriv[gpid]=der;
          dxyzdrs[gpid]=dx;
        }

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

          surfstressman->SurfaceCalc(dxyzdrs, deriv, gpcoord, gpweight, ngp, nnode, A, Adiff, Adiff2);
          surfstressman->StiffnessAndInternalForces(curvenum, A, Adiff, Adiff2, elevector1, elematrix1, this->Id(),
                                                    time, dt, 0, 0.0, k1xC, k2, m1, m2, gamma_0,
                                                    gamma_min, gamma_min_eq, con_quot_max,
                                                    con_quot_eq);
        }
        else if (cond->Type()==DRT::Condition::SurfaceTension) // ideal liquid
        {
          int curvenum = cond->Getint("curve");
          double const_gamma = cond->GetDouble("gamma");

          surfstressman->SurfaceCalc(dxyzdrs, deriv, gpcoord, gpweight, ngp, nnode, A, Adiff, Adiff2);
          surfstressman->StiffnessAndInternalForces(curvenum, A, Adiff, Adiff2, elevector1, elematrix1, this->Id(),
                                                    time, dt, 1, const_gamma, 0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
        else
          dserror("Unknown condition type %d",cond->Type());
      }
      break;

      //compute the area (e.g. for initialization)
      case calc_struct_monitarea:
        {
          if (distype!=quad4)
            dserror("Area Constraint online works for quad4 surfaces!");
          //We are not interested in volume of ghosted elements
          if(Comm.MyPID()==this->Owner())
          {
            // element geometry update
            RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
            if (disp==null) dserror("Cannot get state vector 'displacement'");
            vector<double> mydisp(lm.size());
            DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
            const int numnod = 4;
            Epetra_SerialDenseMatrix xsrefe(numnod,NUMDIM_SOH8);  // material coord. of element
            Epetra_SerialDenseMatrix xscurr(numnod,NUMDIM_SOH8);  // material coord. of element

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

            if (!type) protype =none;
            else if (*type == "xy") protype = xy;
            else if (*type == "yz") protype = yz;
            else if (*type == "xz") protype = xz;
            else protype = none;

            for (int i=0; i<numnod; ++i)
            {
              xsrefe(i,0) = Nodes()[i]->X()[0];
              xsrefe(i,1) = Nodes()[i]->X()[1];
              xsrefe(i,2) = Nodes()[i]->X()[2];

              // To compute monitored area consider required projection method
              // and set according coordinates to zero
              if (protype == yz)
              {
                xscurr(i,0)=0;
              }
              else
              {
                xscurr(i,0) = xsrefe(i,0) + mydisp[i*NODDOF_SOH8+0];
              }
              if (protype == xz)
              {
                xscurr(i,1)=0;
              }
              else
              {
                xscurr(i,1) = xsrefe(i,1) + mydisp[i*NODDOF_SOH8+1];
              }
              if (protype == xy)
              {
                xscurr(i,2)=0;
              }
              else
              {
                xscurr(i,2) = xsrefe(i,2) + mydisp[i*NODDOF_SOH8+2];
              }
            }

            double areaele=0.0;
            const int ngp = 4;

            // gauss parameters
            const double gpweight = 1.0;
            const double gploc    = 1.0/sqrt(3.0);
            Epetra_SerialDenseMatrix gpcoord (ngp,2);
            gpcoord(0,0) = - gploc;
            gpcoord(0,1) = - gploc;
            gpcoord(1,0) =   gploc;
            gpcoord(1,1) = - gploc;
            gpcoord(2,0) = - gploc;
            gpcoord(2,1) =   gploc;
            gpcoord(3,0) =   gploc;
            gpcoord(3,1) =   gploc;

            for (int gpid = 0; gpid < 4; ++gpid) {    // loop over intergration points
              vector<double> funct(ngp);              // 4 shape function values
              double drs;                             // surface area factor
              vector<double> normal(NUMDIM_SOH8);

              soh8_surface_integ(&funct,&drs,&normal,&xscurr,gpcoord(gpid,0),gpcoord(gpid,1));

              double fac = gpweight * drs;
              areaele += fac;

            }

            // get RIGHT volume out of parameterlist and maximum ConditionID
            char areaname[30];
            const int ID =params.get("ConditionID",-1);
            const int maxID=params.get("MaxID",0);
            const int minID=params.get("MinID",1000000);
            if (ID<0)
            {
              dserror("Condition ID for area constraint missing!");
            }
            if (maxID<ID)
            {
              params.set("MaxID",ID);
            }
            if (minID>ID)
            {
              params.set("MinID",ID);
            }
            sprintf(areaname,"computed area %d",ID);
            double areacond = params.get(areaname,0.0);
            //update area in parameter list
            params.set(areaname, areacond+areaele);
          }

        }
        break;
        case calc_struct_areaconstrstiff:
        {
          dserror("Element routines for area constraint not implemented yet!");
        }
        default:
          dserror("Unimplemented type of action for Soh8Surface");

  }
  return 0;
  }

/*----------------------------------------------------------------------*
 * Compute Volume between surface and xy-plane.                 tk 10/07*
 * Yields to the enclosed volume when summed up over all elements       *
 * ---------------------------------------------------------------------*/
double DRT::ELEMENTS::Soh8Surface::ComputeConstrVols(Epetra_SerialDenseMatrix xc)
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
void DRT::ELEMENTS::Soh8Surface::ComputeVolconstrStiff(Epetra_SerialDenseMatrix xc,
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
 * Compute derivative of volume                                 tk 10/07*
 * with respect to the displacements                                    *
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::Soh8Surface::ComputeVolconstrVolDeriv(Epetra_SerialDenseMatrix xc,
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

#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOLID3
