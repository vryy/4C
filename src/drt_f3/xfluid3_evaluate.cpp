/*!
\file xfluid3_evaluate.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef D_FLUID3
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

#include "xfluid3.H"
#include "xfluid3_impl.H"
#include "xfluid3_stationary.H"
#include "xfluid3_interpolation.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_xfem/dof_management.H"

#include <blitz/array.h>
#include <Epetra_SerialDenseSolver.h>

using namespace std;
using namespace DRT::UTILS;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


/*
  Depending on the type of the algorithm (the implementation) and the
  element type (tet, hex etc.), the elements allocate common static
  arrays.

  That means that for example all hex8 fluid elements of the stationary
  implementation have a pointer f8 to the same 'implementation class'
  containing all the element arrays for eight noded elements, and all
  wedge15 fluid elements of the same problem have a pointer f15 to
  the 'implementation class' containing all the element arrays for the
  15 noded element.
  
  */

DRT::ELEMENTS::XFluid3Impl* DRT::ELEMENTS::XFluid3::Impl(
        const int numnode,
        const int numparamvelx,
        const int numparamvely,
        const int numparamvelz,
        const int numparampres)
{
  switch (numnode)
  {
  case 8:
  {
    static XFluid3Impl* f8;
    if (f8==NULL)
      f8 = new XFluid3Impl(numnode,numnode,numnode,numnode,numnode);
    return f8;
  }
  case 20:
  {
    static XFluid3Impl* f20;
    if (f20==NULL)
      f20 = new XFluid3Impl(numnode,numnode,numnode,numnode,numnode);
    return f20;
  }
  case 27:
  {
    static XFluid3Impl* f27;
    if (f27==NULL)
      f27 = new XFluid3Impl(numnode,numnode,numnode,numnode,numnode);
    return f27;
  }
  case 4:
  {
    static XFluid3Impl* f4;
    if (f4==NULL)
      f4 = new XFluid3Impl(numnode,numnode,numnode,numnode,numnode);
    return f4;
  }
  case 10:
  {
    static XFluid3Impl* f10;
    if (f10==NULL)
      f10 = new XFluid3Impl(numnode,numnode,numnode,numnode,numnode);
    return f10;
  }
  case 6:
  {
    static XFluid3Impl* f6;
    if (f6==NULL)
      f6 = new XFluid3Impl(numnode,numnode,numnode,numnode,numnode);
    return f6;
  }
  case 15:
  {
    static XFluid3Impl* f15;
    if (f15==NULL)
      f15 = new XFluid3Impl(numnode,numnode,numnode,numnode,numnode);
    return f15;
  }
  case 5:
  {
    static XFluid3Impl* f5;
    if (f5==NULL)
      f5 = new XFluid3Impl(numnode,numnode,numnode,numnode,numnode);
    return f5;
  }
  default:
  {
    static XFluid3Impl* fx;
    if (fx!=NULL)
        delete fx;
    fx = new XFluid3Impl(numnode,numparamvelx,numparamvely,numparamvelz,numparampres);
    return fx;
  }
  }
  return NULL;
}


// converts a string into an Action for this element
DRT::ELEMENTS::XFluid3::ActionType DRT::ELEMENTS::XFluid3::convertStringToActionType(
              const string& action) const
{
  dsassert(action != "none", "No action supplied");

  DRT::ELEMENTS::XFluid3::ActionType act = XFluid3::none;
  if (action == "calc_fluid_systemmat_and_residual")
    act = XFluid3::calc_fluid_systemmat_and_residual;
  else if (action == "calc_fluid_genalpha_sysmat_and_residual")
    act = XFluid3::calc_fluid_genalpha_sysmat_and_residual;
  else if (action == "time update for subscales")
    act = XFluid3::calc_fluid_genalpha_update_for_subscales;
  else if (action == "time average for subscales and residual")
    act = XFluid3::calc_fluid_genalpha_average_for_subscales_and_residual;
  else if (action == "calc_fluid_stationary_systemmat_and_residual")
    act = XFluid3::calc_fluid_stationary_systemmat_and_residual;  
  else if (action == "calc_fluid_beltrami_error")
    act = XFluid3::calc_fluid_beltrami_error;
  else if (action == "calc_turbulence_statistics")
    act = XFluid3::calc_turbulence_statistics;
  else if (action == "calc_fluid_box_filter")
    act = XFluid3::calc_fluid_box_filter;
  else if (action == "calc_smagorinsky_const")
    act = XFluid3::calc_smagorinsky_const;
  else if (action == "store_xfem_info")
    act = XFluid3::store_xfem_info;
  else
    dserror("Unknown type of action for XFluid3");
  return act;
}

/*----------------------------------------------------------------------*
 // converts a string into an stabilisation action for this element
 //                                                          gammi 08/07
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XFluid3::StabilisationAction DRT::ELEMENTS::XFluid3::ConvertStringToStabAction(
  const string& action) const
{
  DRT::ELEMENTS::XFluid3::StabilisationAction act = stabaction_unspecified;

  map<string,StabilisationAction>::const_iterator iter=stabstrtoact_.find(action);

  if (iter != stabstrtoact_.end())
  {
    act = (*iter).second;
  }
  else
  {
    char errorout[200];
    sprintf(errorout,"looking for stab action (%s) not contained in map",action.c_str());
    dserror(errorout);
  }
  return act;
}


 /*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::XFluid3::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  // get the action required
  const string action = params.get<string>("action","none");
  const DRT::ELEMENTS::XFluid3::ActionType act = convertStringToActionType(action);

  // get the material
  RCP<MAT::Material> mat = Material();
  if (mat->MaterialType()!=m_fluid)
    dserror("newtonian fluid material expected but got type %d", mat->MaterialType());

  MATERIAL* actmat = static_cast<MAT::NewtonianFluid*>(mat.get())->MaterialData();

  switch(act)
  {
      case calc_fluid_systemmat_and_residual:
      {
        // first task is to compare that the dofs fit to the current state of the dofmanager
        // get access to global dofman
        const RCP<XFEM::DofManager> globaldofman = params.get< RCP< XFEM::DofManager > >("dofmanager",null);
        if (globaldofman == null) dserror("hey, you did not give me a globaldofmanager (says xfluid3 evaluate)!!!");
        
        // check for outdated dof information
        globaldofman->checkForConsistency((*this), eleDofManager_);
          
        // need current velocity and history vector
        RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
        RCP<const Epetra_Vector> hist  = discretization.GetState("hist");
        if (velnp==null || hist==null)
          dserror("Cannot get state vectors 'velnp' and/or 'hist'");

        // extract local values from the global vectors
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);
        vector<double> myhist(lm.size());
        DRT::UTILS::ExtractMyValues(*hist,myhist,lm);

        RCP<const Epetra_Vector> dispnp;
        vector<double> mydispnp;
        RCP<const Epetra_Vector> gridv;
        vector<double> mygridv;

        if (is_ale_)
        {
          dispnp = discretization.GetState("dispnp");
          if (dispnp==null) dserror("Cannot get state vectors 'dispnp'");
          mydispnp.resize(lm.size());
          DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);

          gridv = discretization.GetState("gridv");
          if (gridv==null) dserror("Cannot get state vectors 'gridv'");
          mygridv.resize(lm.size());
          DRT::UTILS::ExtractMyValues(*gridv,mygridv,lm);
        }

        const int numparampres = eleDofManager_.NumDofPerField(XFEM::PHYSICS::Pres);
        const int numparamvelx = eleDofManager_.NumDofPerField(XFEM::PHYSICS::Velx);
        const int numparamvely = eleDofManager_.NumDofPerField(XFEM::PHYSICS::Vely);
        const int numparamvelz = eleDofManager_.NumDofPerField(XFEM::PHYSICS::Velz);
        dsassert((numparamvelx == numparamvely and numparamvelx == numparamvelz and numparamvelx == numparampres),
                "for now, we enrich velocity and pressure together");
        
        // create blitz objects for element arrays
        blitz::Array<double, 1> eprenp(numparampres);
        blitz::Array<double, 2> evelnp(3,numparamvelx,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> evhist(3,numparamvelx,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> edispnp(3,numparamvelx,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> egridv(3,numparamvelx,blitz::ColumnMajorArray<2>());

        const vector<int> velxdof = eleDofManager_.LocalDofPosPerField(XFEM::PHYSICS::Velx);
        const vector<int> velydof = eleDofManager_.LocalDofPosPerField(XFEM::PHYSICS::Vely);
        const vector<int> velzdof = eleDofManager_.LocalDofPosPerField(XFEM::PHYSICS::Velz);
        const vector<int> presdof = eleDofManager_.LocalDofPosPerField(XFEM::PHYSICS::Pres);
        
        // split velocity and pressure, insert into element arrays
        for (int iparam=0; iparam<numparamvelx; ++iparam)   evelnp(0,iparam) = myvelnp[velxdof[iparam]];
        for (int iparam=0; iparam<numparamvely; ++iparam)   evelnp(1,iparam) = myvelnp[velydof[iparam]];
        for (int iparam=0; iparam<numparamvelz; ++iparam)   evelnp(2,iparam) = myvelnp[velzdof[iparam]];
        for (int iparam=0; iparam<numparampres; ++iparam)   eprenp(iparam) = myvelnp[presdof[iparam]];
        
        // the history vector contains the information of time step t_n (mass rhs!)
        for (int iparam=0; iparam<numparamvelx; ++iparam)   evhist(0,iparam) = myhist[velxdof[iparam]];
        for (int iparam=0; iparam<numparamvely; ++iparam)   evhist(1,iparam) = myhist[velydof[iparam]];
        for (int iparam=0; iparam<numparamvelz; ++iparam)   evhist(2,iparam) = myhist[velzdof[iparam]];

        if (is_ale_)
        {
          // assign grid velocity and grid displacement to element arrays
          for (int iparam=0; iparam<numparamvelx; ++iparam)   edispnp(0,iparam) = mydispnp[velxdof[iparam]];
          for (int iparam=0; iparam<numparamvely; ++iparam)   edispnp(1,iparam) = mydispnp[velydof[iparam]];
          for (int iparam=0; iparam<numparamvelz; ++iparam)   edispnp(2,iparam) = mydispnp[velzdof[iparam]];
          
          for (int iparam=0; iparam<numparamvelx; ++iparam)   egridv(0,iparam) = mygridv[velxdof[iparam]];
          for (int iparam=0; iparam<numparamvely; ++iparam)   egridv(1,iparam) = mygridv[velydof[iparam]];
          for (int iparam=0; iparam<numparamvelz; ++iparam)   egridv(2,iparam) = mygridv[velzdof[iparam]];
        }

        const RCP<XFEM::InterfaceHandle> ih = params.get< RCP< XFEM::InterfaceHandle > >("interfacehandle",null);
        dsassert(ih!=null, "you did not give the InterfaceHandle");
        
        //! information about domain integration cells
        const XFEM::DomainIntCells   domainIntCells   = ih->GetDomainIntCells(this->Id(),this->Shape());
        //! information about boundary integration cells
        const XFEM::BoundaryIntCells boundaryIntCells = ih->GetBoundaryIntCells(this->Id());
        
        // get control parameter
        const double time = params.get<double>("total time",-1.0);

        const bool newton = params.get<bool>("include reactive terms for linearisation",false);

        // the stabilisation scheme is hardcoded up to now --- maybe it's worth taking
        // this into the input or to choose a standard implementation and drop all ifs
        // on the element level -- if so, I would recommend to drop vstab...
        
        const bool pstab  = true;
        const bool supg   = true;
        const bool vstab  = true;
        const bool cstab  = true;

        // One-step-Theta: timefac = theta*dt
        // BDF2:           timefac = 2/3 * dt
        const double timefac = params.get<double>("thsl",-1.0);
        if (timefac < 0.0) dserror("No thsl supplied");

        // get flag for (fine-scale) subgrid viscosity (1=yes, 0=no)
        const int fssgv = params.get<int>("fs subgrid viscosity",0);

        // wrap epetra serial dense objects in blitz objects
        blitz::Array<double, 2> estif(elemat1.A(),
                                      blitz::shape(elemat1.M(),elemat1.N()),
                                      blitz::neverDeleteData,
                                      blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> esv(elemat2.A(),
                                    blitz::shape(elemat2.M(),elemat2.N()),
                                    blitz::neverDeleteData,
                                    blitz::ColumnMajorArray<2>());
        blitz::Array<double, 1> eforce(elevec1.Values(),
                                       blitz::shape(elevec1.Length()),
                                       blitz::neverDeleteData);

        // calculate element coefficient matrix and rhs     
        Impl(
                NumNode(),
                numparamvelx, 
                numparamvely, 
                numparamvelz, 
                numparampres)->Sysmat(this,
                       domainIntCells,
                       boundaryIntCells,
                       evelnp,
                       eprenp,
                       evhist,
                       edispnp,
                       egridv,
                       estif,
                       esv,
                       eforce,
                       actmat,
                       time,
                       timefac,
                       newton ,
                       fssgv  ,
                       pstab  ,
                       supg   ,
                       vstab  ,
                       cstab  );

        // This is a very poor way to transport the density to the
        // outside world. Is there a better one?
        params.set("density", actmat->m.fluid->density);

      }
      break;
      case calc_fluid_beltrami_error:
      {
        // add error only for elements which are not ghosted
        if(this->Owner() == discretization.Comm().MyPID())
        {

          // need current velocity and history vector
          RCP<const Epetra_Vector> vel_pre_np = discretization.GetState("u and p at time n+1 (converged)");
          if (vel_pre_np==null) dserror("Cannot get state vectors 'velnp'");

          // extract local values from the global vectors
          vector<double> my_vel_pre_np(lm.size());
          DRT::UTILS::ExtractMyValues(*vel_pre_np,my_vel_pre_np,lm);

          // split "my_vel_pre_np" into velocity part "myvelnp" and pressure part "myprenp"
          const int numnode = NumNode();
          vector<double> myprenp(numnode);
          vector<double> myvelnp(3*numnode);

          for (int i=0;i<numnode;++i)
          {
            myvelnp[0+(i*3)]=my_vel_pre_np[0+(i*4)];
            myvelnp[1+(i*3)]=my_vel_pre_np[1+(i*4)];
            myvelnp[2+(i*3)]=my_vel_pre_np[2+(i*4)];

            myprenp[i]=my_vel_pre_np[3+(i*4)];
          }

          // integrate beltrami error
          f3_int_beltrami_err(myvelnp,myprenp,actmat,params);
        }
      }
      break;
      case calc_turbulence_statistics:
      {
        // do nothing if you do not own this element
        if(this->Owner() == discretization.Comm().MyPID())
        {
          // --------------------------------------------------
          // extract velocities and pressure from the global distributed vectors

          // velocity and pressure values (n+1)
          RCP<const Epetra_Vector> velnp
            = discretization.GetState("u and p (n+1,converged)");


          if (velnp==null)
          {
            dserror("Cannot get state vectors 'velnp'");
          }

          // extract local values from the global vectors
          vector<double> mysol  (lm.size());
          DRT::UTILS::ExtractMyValues(*velnp,mysol,lm);

          // integrate mean values
          f3_calc_means(mysol,params);

        }
      }
      break;
      case calc_fluid_genalpha_update_for_subscales:
      {
        // most recent subscale pressure becomes the old subscale pressure
        // for the next timestep
        //
        //  ~n   ~n+1
        //  p <- p
        //
        sub_pre_old_ = sub_pre_;

        // the old subscale acceleration for the next timestep is calculated
        // on the fly, not stored on the element
        /*
                     ~n+1   ~n
             ~ n     u    - u     ~ n   / 1.0-gamma \
            acc  <-  --------- - acc * |  ---------  |
                     gamma*dt           \   gamma   /
        */
        
        const double dt     = params.get<double>("dt");
        const double gamma  = params.get<double>("gamma");

        sub_acc_old_ = (sub_vel_-sub_vel_old_)/(gamma*dt)
                       -
                       sub_acc_old_*(1.0-gamma)/gamma;

        // most recent subscale velocity becomes the old subscale velocity
        // for the next timestep
        //
        //  ~n   ~n+1
        //  u <- u
        //
        sub_vel_old_=sub_vel_;
      }
      break;
      case calc_fluid_stationary_systemmat_and_residual:
      {
          //cout << endl << "Processing element with GiD id = " << (this->Id()+1) <<":" << endl;
          
          // first task is to compare that the dofs fit to the current state of the dofmanager
          // get access to global dofman
          const RCP<XFEM::DofManager> globaldofman = params.get< RCP< XFEM::DofManager > >("dofmanager",null);
          if (globaldofman == null) dserror("hey, you did not give me a globaldofmanager (says xfluid3 evaluate)!!!");
          
          // check for outdated dof information
          globaldofman->checkForConsistency((*this), eleDofManager_);
          
          // get access to interface information
          const RCP<XFEM::InterfaceHandle> ih = params.get< RCP< XFEM::InterfaceHandle > >("interfacehandle",null);
          dsassert(ih!=null, "hey, you did not give the InterfaceHandle");
          
          // need current velocity/pressure 
          RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
          if (velnp==null)
            dserror("Cannot get state vector 'velnp'");

          // extract local values from the global vector
          vector<double> locval(lm.size());
          DRT::UTILS::ExtractMyValues(*velnp,locval,lm);
          //cout << "number of unknowns (node + element): " << lm.size() << endl;
          //cout << "number of unknowns (node):  " << (eleDofManager_.NumDofPerField(XFEM::PHYSICS::Velx)*4) << endl;
          //cout << "number of unknowns (tauxx): " << (eleDofManager_.NumDofPerField(XFEM::PHYSICS::Tauxx)) << endl;
          //cout << "NumDofPerNode(0) " << this->NumDofPerNode(*(this->Nodes()[0])) << endl;
          //cout << "NumDofPerNode(1) " << this->NumDofPerNode(*(this->Nodes()[1])) << endl;
          //cout << "NumDofPerNode(2) " << this->NumDofPerNode(*(this->Nodes()[2])) << endl;
          //cout << "NumDofPerElement " << this->NumDofPerElement() << endl;
          
          // do no calculation, if not needed
          if (lm.size() == 0)
              break;

          if (is_ale_)
          {
        	  dserror("No ALE support within stationary fluid solver.");
          }
          
          // get control parameter
          const double pseudotime = params.get<double>("total time",-1.0);
          if (pseudotime < 0.0)
        	  dserror("no value for total (pseudo-)time in the parameter list");

          const bool newton = params.get<bool>("include reactive terms for linearisation",false);
          const bool pstab  = true;
          const bool supg   = true;
          const bool vstab  = false;  // viscous stabilisation part switched off !!
          const bool cstab  = true;        

          // wrap epetra serial dense objects in blitz objects
          blitz::Array<double, 2> estif(elemat1.A(),
                                        blitz::shape(elemat1.M(),elemat1.N()),
                                        blitz::neverDeleteData,
                                        blitz::ColumnMajorArray<2>());
          blitz::Array<double, 1> eforce(elevec1.Values(),
                                         blitz::shape(elevec1.Length()),
                                         blitz::neverDeleteData);

          const XFEM::AssemblyType assembly_type = CheckForStandardEnrichmentsOnly(
                  eleDofManager_, NumNode(), NodeIds());
          
          // calculate element coefficient matrix and rhs
          XFLUID::callSysmat(assembly_type,
                  this, ih, eleDofManager_, locval, estif, eforce,
                  actmat, pseudotime, newton, pstab, supg, vstab, cstab);

          // This is a very poor way to transport the density to the
          // outside world. Is there a better one?
          params.set("density", actmat->m.fluid->density);	  
          break;
      }
      case store_xfem_info:
      {
    	// get access to global dofman
    	const RCP<XFEM::DofManager> globaldofman = params.get< RCP< XFEM::DofManager > >("dofmanager",null);
    	
    	const DRT::Element::DiscretizationType stressdistype = XFLUID::getStressInterpolationType3D(this->Shape());
    	const int numvirtualnodes = DRT::UTILS::getNumberOfElementNodes(stressdistype);
    	
    	// create local copy of information about dofs
    	eleDofManager_ = globaldofman->constructElementDofManager((*this), numvirtualnodes);
    	break;
      }
      default:
        dserror("Unknown type of action for XFluid3");
  } // end of switch(act)

  return 0;
} // end of DRT::ELEMENTS::Fluid3::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                      gammi 04/07|
 |                                                                      |
 |  The function is just a dummy. For the fluid elements, the           |
 |  integration of the volume neumann loads takes place in the element. |
 |  We need it there for the stabilisation terms!                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::XFluid3::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  return 0;
}

// get optimal gaussrule for discretization type
GaussRule3D DRT::ELEMENTS::XFluid3::getOptimalGaussrule(const DiscretizationType& distype)
{
    GaussRule3D rule = intrule3D_undefined;
    switch (distype)
    {
    case hex8:
        rule = intrule_hex_8point;
        break;
    case hex20: case hex27:
        rule = intrule_hex_27point;
        break;
    case tet4:
        rule = intrule_tet_4point;
        break;
    case tet10:
        rule = intrule_tet_5point;
        break;
    default:
        dserror("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}


/*---------------------------------------------------------------------*
 |  calculate error for beltrami test problem               gammi 04/07|
 *---------------------------------------------------------------------*/
void DRT::ELEMENTS::XFluid3::f3_int_beltrami_err(
  vector<double>&           evelnp,
  vector<double>&           eprenp,
  struct _MATERIAL*         material,
  ParameterList& 	    params
  )
{
  const int NSD = 3;

  // add element error to "integrated" error
  double velerr = params.get<double>("L2 integrated velocity error");
  double preerr = params.get<double>("L2 integrated pressure error");

  // set element data
  const int iel = NumNode();
  const DiscretizationType distype = this->Shape();

  Epetra_SerialDenseVector  funct(iel);
  Epetra_SerialDenseMatrix  xjm(3,3);
  Epetra_SerialDenseMatrix  deriv(3,iel);

  // get node coordinates of element
  Epetra_SerialDenseMatrix xyze(3,iel);
  for(int inode=0;inode<iel;inode++)
  {
    xyze(0,inode)=Nodes()[inode]->X()[0];
    xyze(1,inode)=Nodes()[inode]->X()[1];
    xyze(2,inode)=Nodes()[inode]->X()[2];
  }

  // set constants for analytical solution
  const double t = params.get("total time",-1.0);
  dsassert (t >= 0.0, "beltrami: no total time for error calculation");

  const double a      = PI/4.0;
  const double d      = PI/2.0;

  // get viscosity
  const double  visc = material->m.fluid->viscosity;

  double         preint;
  vector<double> velint  (3);
  vector<double> xint    (3);

  vector<double> u       (3);

  double         deltap;
  vector<double> deltavel(3);

  // gaussian points
  const GaussRule3D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints3D  intpoints = getIntegrationPoints3D(gaussrule);

  // start loop over integration points
  for (int iquad=0;iquad<intpoints.nquad;iquad++)
  {
    // declaration of gauss point variables
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];
    const double e3 = intpoints.qxg[iquad][2];
    shape_function_3D(funct,e1,e2,e3,distype);
    shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

    /*----------------------------------------------------------------------*
      | calculate Jacobian matrix and it's determinant (private) gammi  07/07|
      | Well, I think we actually compute its transpose....
      |
      |     +-            -+ T      +-            -+
      |     | dx   dx   dx |        | dx   dy   dz |
      |     | --   --   -- |        | --   --   -- |
      |     | dr   ds   dt |        | dr   dr   dr |
      |     |              |        |              |
      |     | dy   dy   dy |        | dx   dy   dz |
      |     | --   --   -- |   =    | --   --   -- |
      |     | dr   ds   dt |        | ds   ds   ds |
      |     |              |        |              |
      |     | dz   dz   dz |        | dx   dy   dz |
      |     | --   --   -- |        | --   --   -- |
      |     | dr   ds   dt |        | dt   dt   dt |
      |     +-            -+        +-            -+
      |
      *----------------------------------------------------------------------*/
    Epetra_SerialDenseMatrix    xjm(NSD,NSD);

    for (int isd=0; isd<NSD; isd++)
    {
      for (int jsd=0; jsd<NSD; jsd++)
      {
        double dum = 0.0;
        for (int inode=0; inode<iel; inode++)
        {
          dum += deriv(isd,inode)*xyze(jsd,inode);
        }
        xjm(isd,jsd) = dum;
      }
    }

    // determinant of jacobian matrix
    const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                       xjm(0,1)*xjm(1,2)*xjm(2,0)+
                       xjm(0,2)*xjm(1,0)*xjm(2,1)-
                       xjm(0,2)*xjm(1,1)*xjm(2,0)-
                       xjm(0,0)*xjm(1,2)*xjm(2,1)-
                       xjm(0,1)*xjm(1,0)*xjm(2,2);

    if(det < 0.0)
    {
        printf("\n");
        printf("GLOBAL ELEMENT NO.%i\n",Id());
        printf("NEGATIVE JACOBIAN DETERMINANT: %lf\n", det);
        dserror("Stopped not regulary!\n");
    }

    const double fac = intpoints.qwgt[iquad]*det;

    // get velocity sol at integration point
    for (int i=0;i<3;i++)
    {
      velint[i]=0.0;
      for (int j=0;j<iel;j++)
      {
        velint[i] += funct[j]*evelnp[i+(3*j)];
      }
    }

    // get pressure sol at integration point
    preint = 0;
    for (int inode=0;inode<iel;inode++)
    {
      preint += funct[inode]*eprenp[inode];
    }

    // get velocity sol at integration point
    for (int isd=0;isd<3;isd++)
    {
      xint[isd]=0.0;
      for (int inode=0;inode<iel;inode++)
      {
        xint[isd] += funct[inode]*xyze(isd,inode);
      }
    }

    // compute analytical pressure
    const double p = -a*a/2.0 *
        ( exp(2.0*a*xint[0])
        + exp(2.0*a*xint[1])
        + exp(2.0*a*xint[2])
        + 2.0 * sin(a*xint[0] + d*xint[1]) * cos(a*xint[2] + d*xint[0]) * exp(a*(xint[1]+xint[2]))
        + 2.0 * sin(a*xint[1] + d*xint[2]) * cos(a*xint[0] + d*xint[1]) * exp(a*(xint[2]+xint[0]))
        + 2.0 * sin(a*xint[2] + d*xint[0]) * cos(a*xint[1] + d*xint[2]) * exp(a*(xint[0]+xint[1]))
        )* exp(-2.0*visc*d*d*t);

    // compute analytical velocities
    u[0] = -a * ( exp(a*xint[0]) * sin(a*xint[1] + d*xint[2]) +
                  exp(a*xint[2]) * cos(a*xint[0] + d*xint[1]) ) * exp(-visc*d*d*t);
    u[1] = -a * ( exp(a*xint[1]) * sin(a*xint[2] + d*xint[0]) +
                  exp(a*xint[0]) * cos(a*xint[1] + d*xint[2]) ) * exp(-visc*d*d*t);
    u[2] = -a * ( exp(a*xint[2]) * sin(a*xint[0] + d*xint[1]) +
                  exp(a*xint[1]) * cos(a*xint[2] + d*xint[0]) ) * exp(-visc*d*d*t);

    // compute difference between analytical solution and numerical solution
    deltap = preint - p;

    for (int isd=0;isd<NSD;isd++)
    {
      deltavel[isd] = velint[isd]-u[isd];
    }

    // add square to L2 error
    for (int isd=0;isd<NSD;isd++)
    {
      velerr += deltavel[isd]*deltavel[isd]*fac;
    }
    preerr += deltap*deltap*fac;

  } // end of loop over integration points


  // we use the parameterlist as a container to transport the calculated
  // errors from the elements to the dynamic routine

  params.set<double>("L2 integrated velocity error",velerr);
  params.set<double>("L2 integrated pressure error",preerr);

  return;
}

/*---------------------------------------------------------------------*
 | Calculate spatial mean values for channel flow (cartesian mesh)
 |                                                           gammi 07/07
 |
 | The necessary element integration is performed in here. The element
 | is cut into at least two (HEX8) or three (quadratic elements) planes,
 | the spatial functions (velocity, pressure etc.) are integrated over
 | this plane and this element contribution is added to a processor local
 | vector (see formulas below for a exact description of the output).
 | The method assumes, that all elements are of the same rectangular
 | shape in the "inplanedirection". In addition, it is assumed that
 | the sampling planes are distributed equidistant in the element.
 |
 |
 |                      ^ normdirect       integration plane
 |                      |                /
 |                      |               /
 |                      |
 |                +-----|-------------+
 |               /|     |            /|
 |              / |     |           / |
 |             /  |     |          /  |
 |            /   |     |         /   |
 |           /    +-----|--------/----+ ---- additional integration
 |          /    /|     |       /    /|      plane (for quadratic elements)
 |         /    / |     |      /    / |
 |        +-------------------+    /  |
 |        |   /   |     *-----|---+------------>
 |        |  /    +----/------|--/----+         inplanedirect[1]
 |        | /    /    /       | /    /
 |        |/    /    /        |/    /   \
 |        +---------+---------+    /     \
 |        |   /    /          |   /       integration plane
 |        |  /    /           |  /
 |        | /    /            | /
 |        |/    /             |/
 |        +----/--------------+
 |            /
 |           /   inplanedirect[0]
 |
 |
 |  Example for a mean value evaluation:
 |
 |         1.0       /
 |  _               |
 |  u = -------- *  | u(x,y,z) dx dy dz =
 |      +---        |
 |       \         / A
 |       / area
 |      +---
 |
 |
 |        1.0      /
 |                |            area
 |  =  -------- * | u(r,s,t) * ---- dr ds dt
 |     +---       |              4
 |      \        /  [-1:1]^2
 |      / area
 |     +---
 |
 |
 |
 |         1.0      /
 |                 |            1
 |  =   -------- * | u(r,s,t) * - dr ds dt
 |                 |            4
 |       numele   /  [-1:1]^2
 |
 |                |                        |
 |                +------------------------+
 |             this is the integral we compute!
 |
 | The factor 1/4 is necessary since we use a reference element of
 | size 2x2
 |
 | The method computes:
 |                      _             _             _             _
 |             numele * u  , numele * v  , numele * w  , numele * p
 |                      ___           ___           ___           ___
 |                       ^2            ^2            ^2            ^2
 | and         numele * u  , numele * v  , numele * w  , numele * p
 |                      _ _           _ _           _ _  
 | as well as  numele * u*v, numele * u*w, numele * v*w
 | 
 | as well as numele.
 | All results are communicated via the parameter list!
 |
 *---------------------------------------------------------------------*/
void DRT::ELEMENTS::XFluid3::f3_calc_means(
  vector<double>&           sol  ,
  ParameterList& 	    params
  )
{

  // set element data
  const int iel = NumNode();
  const DiscretizationType distype = this->Shape();


  // the plane normal tells you in which plane the integration takes place
  const int normdirect = params.get<int>("normal direction to homogeneous plane");


  // the vector planes contains the coordinates of the homogeneous planes (in
  // wall normal direction)
  RCP<vector<double> > planes = params.get<RCP<vector<double> > >("coordinate vector for hom. planes");

  // get the pointers to the solution vectors
  RCP<vector<double> > sumu   = params.get<RCP<vector<double> > >("mean velocity u");
  RCP<vector<double> > sumv   = params.get<RCP<vector<double> > >("mean velocity v");
  RCP<vector<double> > sumw   = params.get<RCP<vector<double> > >("mean velocity w");
  RCP<vector<double> > sump   = params.get<RCP<vector<double> > >("mean pressure p");

  RCP<vector<double> > sumsqu = params.get<RCP<vector<double> > >("mean value u^2");
  RCP<vector<double> > sumsqv = params.get<RCP<vector<double> > >("mean value v^2");
  RCP<vector<double> > sumsqw = params.get<RCP<vector<double> > >("mean value w^2");
  RCP<vector<double> > sumuv  = params.get<RCP<vector<double> > >("mean value uv");
  RCP<vector<double> > sumuw  = params.get<RCP<vector<double> > >("mean value uw");
  RCP<vector<double> > sumvw  = params.get<RCP<vector<double> > >("mean value vw");
  RCP<vector<double> > sumsqp = params.get<RCP<vector<double> > >("mean value p^2");


  // get node coordinates of element
  Epetra_SerialDenseMatrix xyze(3,iel);
  for(int inode=0;inode<iel;inode++)
  {
    xyze(0,inode)=Nodes()[inode]->X()[0];
    xyze(1,inode)=Nodes()[inode]->X()[1];
    xyze(2,inode)=Nodes()[inode]->X()[2];
  }

  double min = xyze(normdirect,0);
  double max = xyze(normdirect,0);

  // set maximum and minimum value in wall normal direction
  for(int inode=0;inode<iel;inode++)
  {
    if(min > xyze(normdirect,inode))
    {
      min=xyze(normdirect,inode);
    }
    if(max < xyze(normdirect,inode))
    {
      max=xyze(normdirect,inode);
    }
  }

  // determine the ids of the homogeneous planes intersecting this element
  set<int> planesinele;
  for(unsigned nplane=0;nplane<planes->size();++nplane)
  {
    // get all available wall normal coordinates
    for(int nn=0;nn<iel;++nn)
    {
      if (min-2e-9 < (*planes)[nplane] && max+2e-9 > (*planes)[nplane])
      {
        planesinele.insert(nplane);
      }
    }
  }

  // remove lowest layer from planesinele to avoid double calculations. This is not done
  // for the first level (index 0) --- if deleted, shift the first integration point in
  // wall normal direction
  // the shift depends on the number of sampling planes in the element
  double shift=0;

  // set the number of planes which cut the element
  const int numplanesinele = planesinele.size();

  if(*planesinele.begin() != 0)
  {
    // this is not an element of the lowest element layer
    planesinele.erase(planesinele.begin());

    shift=2.0/((double) numplanesinele - 1.0);
  }
  else
  {
    // this is an element of the lowest element layer. Increase the counter
    // in order to compute the total number of elements in one layer
    int* count = params.get<int*>("count processed elements");

    (*count)++;
  }

  // determine the orientation of the rst system compared to the xyz system
  int elenormdirect=-1;
  bool upsidedown =false;
  // the only thing of interest is how normdirect is oriented in the
  // element coordinate system
  if(xyze(normdirect,4)-xyze(normdirect,0)>2e-9)
  {
    // t aligned
    elenormdirect =2;
    cout << "upsidedown false" <<&endl;
  }
  else if (xyze(normdirect,3)-xyze(normdirect,0)>2e-9)
  {
    // s aligned
    elenormdirect =1;
  }
  else if (xyze(normdirect,1)-xyze(normdirect,0)>2e-9)
  {
    // r aligned
    elenormdirect =0;
  }
  else if(xyze(normdirect,4)-xyze(normdirect,0)<-2e-9)
  {
    cout << xyze(normdirect,4)-xyze(normdirect,0) << &endl;
    // -t aligned
    elenormdirect =2;
    upsidedown =true;
    cout << "upsidedown true" <<&endl;
  }
  else if (xyze(normdirect,3)-xyze(normdirect,0)<-2e-9)
  {
    // -s aligned
    elenormdirect =1;
    upsidedown =true;
  }
  else if (xyze(normdirect,1)-xyze(normdirect,0)<-2e-9)
  {
    // -r aligned
    elenormdirect =0;
    upsidedown =true;
  }
  else
  {
    dserror("cannot determine orientation of plane normal in local coordinate system of element");
  }
  vector<int> inplanedirect;
  {
    set <int> inplanedirectset;
    for(int i=0;i<3;++i)
    {
      inplanedirectset.insert(i);
    }
    inplanedirectset.erase(elenormdirect);

    for(set<int>::iterator id = inplanedirectset.begin();id!=inplanedirectset.end() ;++id)
    {
      inplanedirect.push_back(*id);
    }
  }

  // allocate vector for shapefunctions
  Epetra_SerialDenseVector  funct(iel);

  // get the quad9 gaussrule for the in plane integration
  const IntegrationPoints2D  intpoints = getIntegrationPoints2D(intrule_quad_9point);

  // a hex8 element has two levels, the hex20 and hex27 element have three layers to sample
  // (now we allow even more)
  double layershift=0;

  // loop all levels in element
  for(set<int>::const_iterator id = planesinele.begin();id!=planesinele.end() ;++id)
  {
    // reset temporary values
    double ubar=0;
    double vbar=0;
    double wbar=0;
    double pbar=0;

    double usqbar=0;
    double vsqbar=0;
    double wsqbar=0;
    double uvbar =0;
    double uwbar =0;
    double vwbar =0;
    double psqbar=0;

    // get the intgration point in wall normal direction
    double e[3];

    e[elenormdirect]=-1.0+shift+layershift;
    if(upsidedown)
    {
      e[elenormdirect]*=-1;
    }

    // start loop over integration points in layer
    for (int iquad=0;iquad<intpoints.nquad;iquad++)
    {
      // get the other gauss point coordinates
      for(int i=0;i<2;++i)
      {
        e[inplanedirect[i]]=intpoints.qxg[iquad][i];
      }

      // compute the shape function values
      shape_function_3D(funct,e[0],e[1],e[2],distype);

      // check whether this gausspoint is really inside the desired plane
      {
        double x[3];
        x[0]=0;
        x[1]=0;
        x[2]=0;
        for(int inode=0;inode<iel;inode++)
        {
          x[0]+=funct[inode]*xyze(0,inode);
          x[1]+=funct[inode]*xyze(1,inode);
          x[2]+=funct[inode]*xyze(2,inode);
        }

        if(abs(x[normdirect]-(*planes)[*id])>2e-9)
        {
          dserror("Mixing up element cut planes during integration");
        }
      }

      //interpolated values at gausspoints
      double ugp=0;
      double vgp=0;
      double wgp=0;
      double pgp=0;

      // we assume that every 2d element we are integrating here is of
      // rectangular shape and every element is of the same size.
      // 1/4 is necessary since we use a reference element of size 2x2
      // the factor fac is omitting the element area up to now
      double fac=0.25*intpoints.qwgt[iquad];

      for(int inode=0;inode<iel;inode++)
      {
        ugp += funct[inode]*sol[inode*4  ];
        vgp += funct[inode]*sol[inode*4+1];
        wgp += funct[inode]*sol[inode*4+2];
        pgp += funct[inode]*sol[inode*4+3];
      }

      // add contribution to integral
      ubar   += ugp*fac;
      vbar   += vgp*fac;
      wbar   += wgp*fac;
      pbar   += pgp*fac;

      usqbar += ugp*ugp*fac;
      vsqbar += vgp*vgp*fac;
      wsqbar += wgp*wgp*fac;
      uvbar  += ugp*vgp*fac;
      uwbar  += ugp*wgp*fac;
      vwbar  += vgp*wgp*fac;
      psqbar += pgp*pgp*fac;
    } // end loop integration points

    // add increments from this layer to processor local vectors
    (*sumu  )[*id] += ubar;
    (*sumv  )[*id] += vbar;
    (*sumw  )[*id] += wbar;
    (*sump  )[*id] += pbar;

    (*sumsqu)[*id] += usqbar;
    (*sumsqv)[*id] += vsqbar;
    (*sumsqw)[*id] += wsqbar;
    (*sumuv) [*id] += uvbar;
    (*sumuw) [*id] += uwbar;
    (*sumvw) [*id] += vwbar;
    (*sumsqp)[*id] += psqbar;

    // jump to the next layer in the element.
    // in case of an hex8 element, the two coordinates are -1 and 1(+2)
    // for quadratic elements with three sample planes, we have -1,0(+1),1(+2)

    layershift+=2.0/((double) numplanesinele - 1.0);
  }


  return;
} // DRT::ELEMENTS::Fluid3::f3_calc_means

//----------------------------------------------------------------------
//
//----------------------------------------------------------------------
void DRT::ELEMENTS::XFluid3::f3_apply_box_filter(
    blitz::Array<double, 2>&  evelaf,
    blitz::Array<double, 1>&  vel_hat,
    blitz::Array<double, 2>&  reystr_hat,
    blitz::Array<double, 2>&  modeled_stress_grid_scale_hat,
    double&                   volume
    )
{

  //------------------------------------------------------------------
  //                     BLITZ CONFIGURATION
  //------------------------------------------------------------------
  //
  // We define the variables i,j,k to be indices to blitz arrays.
  // These are used for array expressions, that is matrix-vector
  // products in the following.

  blitz::firstIndex  i;   // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex  k;   // Placeholder for the third index

  // set element data
  const int iel = NumNode();
  const DiscretizationType distype = this->Shape();

  // allocate arrays for shapefunctions, derivatives and the transposed jacobian
  blitz::Array<double,1>  funct(iel);
  blitz::Array<double,2>  xjm  (3,3);
  blitz::Array<double,2>  deriv(3,iel,blitz::ColumnMajorArray<2>());


  // get node coordinates of element
  blitz::Array<double,2>  xyze(3,iel);
  for(int inode=0;inode<iel;inode++)
  {
    xyze(0,inode)=Nodes()[inode]->X()[0];
    xyze(1,inode)=Nodes()[inode]->X()[1];
    xyze(2,inode)=Nodes()[inode]->X()[2];
  }

  // use one point gauss rule to calculate tau at element center
  DRT::UTILS::GaussRule3D integrationrule_filter=DRT::UTILS::intrule_hex_1point;
  switch (distype)
  {
      case DRT::Element::hex8:
        integrationrule_filter = DRT::UTILS::intrule_hex_1point;
        break;
      case DRT::Element::tet4:
        integrationrule_filter = DRT::UTILS::intrule_tet_1point;
        break;
      case DRT::Element::tet10:
      case DRT::Element::hex20:
      case DRT::Element::hex27:
        dserror("the box filtering operation is only permitted for linear elements\n");
        break;
      default:
        dserror("invalid discretization type for fluid3");
  }

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints_onepoint(integrationrule_filter);
  
  // shape functions and derivs at element center
  const double e1    = intpoints_onepoint.qxg[0][0];
  const double e2    = intpoints_onepoint.qxg[0][1];
  const double e3    = intpoints_onepoint.qxg[0][2];
  const double wquad = intpoints_onepoint.qwgt[0];
  
  DRT::UTILS::shape_function_3D       (funct,e1,e2,e3,distype);
  DRT::UTILS::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

  // get Jacobian matrix and determinant
  xjm = blitz::sum(deriv(i,k)*xyze(j,k),k);
  const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                     xjm(0,1)*xjm(1,2)*xjm(2,0)+
                     xjm(0,2)*xjm(1,0)*xjm(2,1)-
                     xjm(0,2)*xjm(1,1)*xjm(2,0)-
                     xjm(0,0)*xjm(1,2)*xjm(2,1)-
                     xjm(0,1)*xjm(1,0)*xjm(2,2);

  //
  //             compute global first derivates
  //
  blitz::Array<double,2>  derxy(3,iel,blitz::ColumnMajorArray<2>());
  /*
    Use the Jacobian and the known derivatives in element coordinate
    directions on the right hand side to compute the derivatives in
    global coordinate directions

          +-                 -+     +-    -+      +-    -+
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dr    dr    dr   |     |  dx  |      |  dr  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |  *  | ---- |   =  | ---- | for all k
          |  ds    ds    ds   |     |  dy  |      |  ds  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dt    dt    dt   |     |  dz  |      |  dt  |
          +-                 -+     +-    -+      +-    -+

          Do one LU factorisation, everything else is backward substitution!

  */
  
  {
    // LAPACK solver
    Epetra_LAPACK          solver;
    
    // this copy of xjm will be used to calculate a in place factorisation
    blitz::Array<double,2> factorU(3,3,blitz::ColumnMajorArray<2>());
    factorU=xjm.copy();
    
    // a vector specifying the pivots (reordering)
    int pivot[3];

    // error code
    int ierr = 0;

    // Perform LU factorisation
    solver.GETRF(3,3,factorU.data(),3,&(pivot[0]),&ierr);
    
    if (ierr!=0)
    {
      dserror("Unable to perform LU factorisation during computation of derxy");
    }
    
    // backward substitution. The copy is required since GETRS replaces
    // the input with the result
    derxy =deriv.copy();
    solver.GETRS('N',3,iel,factorU.data(),3,&(pivot[0]),derxy.data(),3,&ierr);
    
    if (ierr!=0)
    {
      dserror("Unable to perform backward substitution after factorisation of jacobian");
    }
  }
  
  // get velocities (n+alpha_F,i) at integration point
  //
  //                 +-----
  //       n+af       \                  n+af
  //    vel    (x) =   +      N (x) * vel
  //                  /        j         j
  //                 +-----
  //                 node j
  //
  blitz::Array<double,1> velintaf (3);
  velintaf = blitz::sum(funct(j)*evelaf(i,j),j);


  // get velocity (n+alpha_F,i) derivatives at integration point
  //
  //       n+af      +-----  dN (x)
  //   dvel    (x)    \        k         n+af
  //   ----------- =   +     ------ * vel
  //       dx         /        dx        k
  //         j       +-----      j
  //                 node k
  //
  // j : direction of derivative x/y/z
  //
  blitz::Array<double,2>  vderxyaf(3,3);
  vderxyaf = blitz::sum(derxy(j,k)*evelaf(i,k),k);

  /*                            
                            +-     n+af          n+af    -+
          / h \       1.0   |  dvel_i  (x)   dvel_j  (x)  |
     eps | u   |    = --- * |  ----------- + -----------  |
          \   / ij    2.0   |      dx            dx       |
                            +-       j             i     -+
  */
  blitz::Array<double,2> epsilon(3,3,blitz::ColumnMajorArray<2>());
  epsilon = 0.5 * ( vderxyaf(i,j) + vderxyaf(j,i) );
  
  //
  // modeled part of subgrid scale stresses
  // 
  /*    +-                                 -+ 1                 
        |          / h \           / h \    | -         / h \   
        | 2 * eps | u   |   * eps | u   |   | 2  * eps | u   |  
        |          \   / kl        \   / kl |           \   / ij
        +-                                 -+                   
      
        |                                   |
        +-----------------------------------+
             'resolved' rate of strain
  */
  
  double rateofstrain = 0;
  
  for(int rr=0;rr<3;rr++)
  {
    for(int mm=0;mm<3;mm++)
    {
      rateofstrain += epsilon(rr,mm)*epsilon(rr,mm);
    }
  }
  rateofstrain *= 2.0;
  rateofstrain = sqrt(rateofstrain);

  
  //--------------------------------------------------
  // one point integrations

  // determine contribution to patch volume
  volume = wquad*det;

  // add contribution to integral over velocities
  vel_hat += velintaf*wquad*det;
  
  // add contribution to integral over reynolds stresses
  reystr_hat += velintaf(i)*velintaf(j)*wquad*det;
  
  // add contribution to integral over the modeled part of subgrid
  // scale stresses 
  modeled_stress_grid_scale_hat += rateofstrain * epsilon * wquad*det;

  return;
} // DRT::ELEMENTS::Fluid3::f3_apply_box_filter


void DRT::ELEMENTS::XFluid3::f3_calc_smag_const_LijMij_and_MijMij(
  blitz::Array<double, 2> evelaf_hat,
  blitz::Array<double, 3> ereynoldsstress_hat,
  blitz::Array<double, 3> efiltered_modeled_subgrid_stress_hat,
  double&                 LijMij,
  double&                 MijMij,
  double&                 center)
{

  //------------------------------------------------------------------
  //                     BLITZ CONFIGURATION
  //------------------------------------------------------------------
  //
  // We define the variables i,j,k to be indices to blitz arrays.
  // These are used for array expressions, that is matrix-vector
  // products in the following.

  blitz::firstIndex  i;   // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex  k;   // Placeholder for the third index

  // set element data
  const int iel = NumNode();
  const DiscretizationType distype = this->Shape();

  // allocate arrays for shapefunctions, derivatives and the transposed jacobian
  blitz::Array<double,1>  funct(iel);
  blitz::Array<double,2>  xjm  (3,3);
  blitz::Array<double,2>  deriv(3,iel,blitz::ColumnMajorArray<2>());


  //this will be the y-coordinate of a point in the element interior
  center = 0;
  
  // get node coordinates of element
  blitz::Array<double,2>  xyze(3,iel);
  for(int inode=0;inode<iel;inode++)
  {
    xyze(0,inode)=Nodes()[inode]->X()[0];
    xyze(1,inode)=Nodes()[inode]->X()[1];
    xyze(2,inode)=Nodes()[inode]->X()[2];

    center+=xyze(1,inode);
  }
  center/=iel;
  

  // use one point gauss rule to calculate tau at element center
  DRT::UTILS::GaussRule3D integrationrule_filter=DRT::UTILS::intrule_hex_1point;
  switch (distype)
  {
      case DRT::Element::hex8:
      case DRT::Element::hex20:
      case DRT::Element::hex27:
        integrationrule_filter = DRT::UTILS::intrule_hex_1point;
        break;
      case DRT::Element::tet4:
      case DRT::Element::tet10:
        integrationrule_filter = DRT::UTILS::intrule_tet_1point;
        break;
      default:
        dserror("invalid discretization type for fluid3");
  }

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints_onepoint(integrationrule_filter);
  const double e1    = intpoints_onepoint.qxg[0][0];
  const double e2    = intpoints_onepoint.qxg[0][1];
  const double e3    = intpoints_onepoint.qxg[0][2];
  
  // shape functions and derivs at element center
  DRT::UTILS::shape_function_3D       (funct,e1,e2,e3,distype);
  DRT::UTILS::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);
  
  // get element type constant for tau
  double mk=0.0;
  switch (distype)
  {
      case DRT::Element::tet4:
      case DRT::Element::hex8:
        mk = 0.333333333333333333333;
        break;
      case DRT::Element::hex20:
      case DRT::Element::hex27:
      case DRT::Element::tet10:
        mk = 0.083333333333333333333;
        break;
      default:
        dserror("type unknown!\n");
  }

  // get Jacobian matrix 
  xjm = blitz::sum(deriv(i,k)*xyze(j,k),k);

  //
  //             compute global first derivates
  //
  blitz::Array<double,2>  derxy(3,iel,blitz::ColumnMajorArray<2>());
  /*
    Use the Jacobian and the known derivatives in element coordinate
    directions on the right hand side to compute the derivatives in
    global coordinate directions

          +-                 -+     +-    -+      +-    -+
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dr    dr    dr   |     |  dx  |      |  dr  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |  *  | ---- |   =  | ---- | for all k
          |  ds    ds    ds   |     |  dy  |      |  ds  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dt    dt    dt   |     |  dz  |      |  dt  |
          +-                 -+     +-    -+      +-    -+

          Do one LU factorisation, everything else is backward substitution!

  */
  
  {
    // LAPACK solver
    Epetra_LAPACK          solver;
    
    // this copy of xjm will be used to calculate a in place factorisation
    blitz::Array<double,2> factorU(3,3,blitz::ColumnMajorArray<2>());
    factorU=xjm.copy();
    
    // a vector specifying the pivots (reordering)
    int pivot[3];

    // error code
    int ierr = 0;

    // Perform LU factorisation
    solver.GETRF(3,3,factorU.data(),3,&(pivot[0]),&ierr);
    
    if (ierr!=0)
    {
      dserror("Unable to perform LU factorisation during computation of derxy");
    }
    
    // backward substitution. The copy is required since GETRS replaces
    // the input with the result
    derxy =deriv.copy();
    solver.GETRS('N',3,iel,factorU.data(),3,&(pivot[0]),derxy.data(),3,&ierr);
    
    if (ierr!=0)
    {
      dserror("Unable to perform backward substitution after factorisation of jacobian");
    }
  }
  
  // get velocities (n+alpha_F,i) at integration point
  //
  //                 +-----
  //     ^ n+af       \                ^ n+af
  //    vel    (x) =   +      N (x) * vel
  //                  /        j         j
  //                 +-----
  //                 node j
  //
  blitz::Array<double,1> velintaf_hat (3);
  velintaf_hat = blitz::sum(funct(j)*evelaf_hat(i,j),j);


  // get velocity (n+alpha_F,i) derivatives at integration point
  //
  //     ^ n+af      +-----  dN (x)
  //   dvel    (x)    \        k       ^ n+af
  //   ----------- =   +     ------ * vel
  //       dx         /        dx        k
  //         j       +-----      j
  //                 node k
  //
  // j : direction of derivative x/y/z
  //
  blitz::Array<double,2>  vderxyaf_hat(3,3);
  vderxyaf_hat = blitz::sum(derxy(j,k)*evelaf_hat(i,k),k);

  // get filtered reynolds stress (n+alpha_F,i) at integration point
  //
  //                      +-----
  //        ^   n+af       \                   ^   n+af
  //    restress    (x) =   +      N (x) * restress
  //            ij         /        k              k, ij
  //                      +-----
  //                      node k
  //
  blitz::Array<double,2>  restress_hat(3,3);
  restress_hat = blitz::sum(funct(k)*ereynoldsstress_hat(i,j,k),k);

  // get filtered modeled subgrid stress (n+alpha_F,i) at integration point
  //
  //                                                 +-----
  //                   ^                   n+af       \                              ^                   n+af
  //    filtered_modeled_subgrid_stress_hat    (x) =   +      N (x) * filtered_modeled_subgrid_stress_hat
  //                                       ij         /        k                                         k, ij
  //                                                 +-----
  //                                                 node k
  //
  blitz::Array<double,2>  filtered_modeled_subgrid_stress_hat(3,3);
  filtered_modeled_subgrid_stress_hat = blitz::sum(funct(k)*efiltered_modeled_subgrid_stress_hat(i,j,k),k);

  
  /*                            
                            +-   ^ n+af        ^ n+af    -+
      ^   / h \       1.0   |  dvel_i  (x)   dvel_j  (x)  |
     eps | u   |    = --- * |  ----------- + -----------  |
          \   / ij    2.0   |      dx            dx       |
                            +-       j             i     -+
  */
  blitz::Array<double,2> epsilon_hat(3,3,blitz::ColumnMajorArray<2>());
  epsilon_hat = 0.5 * ( vderxyaf_hat(i,j) + vderxyaf_hat(j,i) );
  
  //
  // modeled part of subtestfilter scale stresses
  // 
  /*    +-                                 -+ 1                 
        |      ^   / h \       ^   / h \    | -     ^   / h \   
        | 2 * eps | u   |   * eps | u   |   | 2  * eps | u   |  
        |          \   / kl        \   / kl |           \   / ij
        +-                                 -+                   
      
        |                                   |
        +-----------------------------------+
             'resolved' rate of strain
  */
  
  double rateofstrain_hat = 0;
  
  for(int rr=0;rr<3;rr++)
  {
    for(int mm=0;mm<3;mm++)
    {
      rateofstrain_hat += epsilon_hat(rr,mm)*epsilon_hat(rr,mm);
    }
  }
  rateofstrain_hat *= 2.0;
  rateofstrain_hat = sqrt(rateofstrain_hat);
  
  blitz::Array<double, 2> L_ij (3,3,blitz::ColumnMajorArray<2>());
  blitz::Array<double, 2> M_ij (3,3,blitz::ColumnMajorArray<2>());

  L_ij = restress_hat(i,j) - velintaf_hat(i)*velintaf_hat(j);

  // this is sqrt(3)
  double filterwidthratio = 1.73;

  M_ij = filtered_modeled_subgrid_stress_hat(i,j)
         -
         filterwidthratio*filterwidthratio*rateofstrain_hat*epsilon_hat(i,j);

  LijMij =0;
  MijMij =0;
  for(int rr=0;rr<3;rr++)
  {
    for(int mm=0;mm<3;mm++)
    {
      LijMij += L_ij(rr,mm)*M_ij(rr,mm);
      MijMij += M_ij(rr,mm)*M_ij(rr,mm);
    }
  }
  
  return;
} // DRT::ELEMENTS::Fluid3::f3_calc_smag_const_LijMij_and_MijMij

//
// check for higher order derivatives for shape functions
//
bool DRT::ELEMENTS::XFluid3::isHigherOrderElement(
  const DRT::Element::DiscretizationType  distype) const
{
  bool hoel = true;
  switch (distype)
  {
  case hex8: case hex20: case hex27: case tet10: case wedge15:
    hoel = true;
    break;
  case tet4: case wedge6: case pyramid5:
    hoel = false;
    break;
  default:
    dserror("distype unknown!");
  }
  return hoel;
}

//
// check for element rewinding based on Jacobian determinant
//
bool DRT::ELEMENTS::XFluid3::checkRewinding() const
{
  const DiscretizationType distype = this->Shape();
  const int iel = NumNode();
  // use one point gauss rule to calculate tau at element center
  GaussRule3D integrationrule_1point = intrule3D_undefined;
  switch(distype)
  {
  case hex8: case hex20: case hex27:
      integrationrule_1point = intrule_hex_1point;
      break;
  case tet4: case tet10:
      integrationrule_1point = intrule_tet_1point;
      break;
  case wedge6: case wedge15:
      integrationrule_1point = intrule_wedge_1point;
      break;
  case pyramid5:
      integrationrule_1point = intrule_pyramid_1point;
      break;
  default:
      dserror("invalid discretization type for fluid3");
  }
  const IntegrationPoints3D  intpoints = getIntegrationPoints3D(integrationrule_1point);

  // shape functions derivatives
  const int NSD = 3;
  Epetra_SerialDenseMatrix    deriv(NSD, iel);
  Epetra_SerialDenseMatrix    xyze(NSD,iel);
  DRT::UTILS::shape_function_3D_deriv1(deriv,intpoints.qxg[0][0],intpoints.qxg[0][1],intpoints.qxg[0][2],distype);
  // get node coordinates
  const DRT::Node** nodes = this->Nodes();
  for (int inode=0; inode<iel; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze(0,inode) = x[0];
    xyze(1,inode) = x[1];
    xyze(2,inode) = x[2];
  }

  // get Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
  */
  Epetra_SerialDenseMatrix xjm(NSD,NSD);

  xjm.Multiply('N','T',1.0,deriv,xyze,0.0);


  
  const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                     xjm(0,1)*xjm(1,2)*xjm(2,0)+
                     xjm(0,2)*xjm(1,0)*xjm(2,1)-
                     xjm(0,2)*xjm(1,1)*xjm(2,0)-
                     xjm(0,0)*xjm(1,2)*xjm(2,1)-
                     xjm(0,1)*xjm(1,0)*xjm(2,2);
  if (det < 0.0) return true;

  return false;
}

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*----------------------------------------------------------------------*
 |  init the element (public)                                mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::XFluid3Register::Initialize(DRT::Discretization& dis)
{
  bool dofillcompleteagain = false;
  //-------------------- loop all my column elements and check rewinding
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    // get the actual element
    if (dis.lColElement(i)->Type() != DRT::Element::element_fluid3) continue;
    DRT::ELEMENTS::XFluid3* actele = dynamic_cast<DRT::ELEMENTS::XFluid3*>(dis.lColElement(i));
    if (!actele) dserror("cast to XFluid3* failed");
    
    const DRT::Element::DiscretizationType distype = actele->Shape();
    bool possiblytorewind = false;
    switch(distype)
    {
    case DRT::Element::hex8: case DRT::Element::hex20: case DRT::Element::hex27:
        possiblytorewind = true;
        break;
    case DRT::Element::tet4: case DRT::Element::tet10:
        possiblytorewind = true;
        break;
    case DRT::Element::wedge6: case DRT::Element::wedge15:
        possiblytorewind = true;
        break;
    case DRT::Element::pyramid5:
        possiblytorewind = true;
        break;
    default:
        dserror("invalid discretization type for fluid3");
    }
    
    if ( (possiblytorewind) && (!actele->donerewinding_) ) {
      actele->rewind_ = actele->checkRewinding();

      if (actele->rewind_) {
        if (distype==DRT::Element::tet4){
          int iel = actele->NumNode();
          int new_nodeids[iel];
          const int* old_nodeids;
          old_nodeids = actele->NodeIds();
          // rewinding of nodes to arrive at mathematically positive element
          new_nodeids[0] = old_nodeids[0];
          new_nodeids[1] = old_nodeids[2];
          new_nodeids[2] = old_nodeids[1];
          new_nodeids[3] = old_nodeids[3];
          actele->SetNodeIds(iel, new_nodeids);
        }
        else if (distype==DRT::Element::hex8){
          int iel = actele->NumNode();
          int new_nodeids[iel];
          const int* old_nodeids;
          old_nodeids = actele->NodeIds();
          // rewinding of nodes to arrive at mathematically positive element
          new_nodeids[0] = old_nodeids[4];
          new_nodeids[1] = old_nodeids[5];
          new_nodeids[2] = old_nodeids[6];
          new_nodeids[3] = old_nodeids[7];
          new_nodeids[4] = old_nodeids[0];
          new_nodeids[5] = old_nodeids[1];
          new_nodeids[6] = old_nodeids[2];
          new_nodeids[7] = old_nodeids[3];
          actele->SetNodeIds(iel, new_nodeids);
        }
        else if (distype==DRT::Element::wedge6){
          int iel = actele->NumNode();
          int new_nodeids[iel];
          const int* old_nodeids;
          old_nodeids = actele->NodeIds();
          // rewinding of nodes to arrive at mathematically positive element
          new_nodeids[0] = old_nodeids[3];
          new_nodeids[1] = old_nodeids[4];
          new_nodeids[2] = old_nodeids[5];
          new_nodeids[3] = old_nodeids[0];
          new_nodeids[4] = old_nodeids[1];
          new_nodeids[5] = old_nodeids[2];
          actele->SetNodeIds(iel, new_nodeids);
        }
        else if (distype == DRT::Element::pyramid5){
          int iel = actele->NumNode();
          int new_nodeids[iel];
          const int* old_nodeids;
          old_nodeids = actele->NodeIds();
          // rewinding of nodes to arrive at mathematically positive element
          new_nodeids[1] = old_nodeids[3];
          new_nodeids[3] = old_nodeids[1];
          // the other nodes can stay the same
          new_nodeids[0] = old_nodeids[0];
          new_nodeids[2] = old_nodeids[2];
          new_nodeids[4] = old_nodeids[4];
          actele->SetNodeIds(iel, new_nodeids);
        }
        else dserror("no rewinding scheme for this type of fluid3");
      }
      // process of rewinding done
      actele->donerewinding_ = true;
      dofillcompleteagain = true;
    }
  }
  // fill complete again to reconstruct element-node pointers,
  // but without element init, etc.
  if(dofillcompleteagain) dis.FillComplete(false,false,false);
  
  return 0;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
