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
#include "xfluid3_sysmat3.H"
#include "xfluid3_sysmat4.H"
#include "xfluid3_interpolation.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_xfem/dof_management.H"

#include <blitz/array.h>


// converts a string into an Action for this element
DRT::ELEMENTS::XFluid3::ActionType DRT::ELEMENTS::XFluid3::convertStringToActionType(
              const string& action) const
{
  dsassert(action != "none", "No action supplied");

  DRT::ELEMENTS::XFluid3::ActionType act = XFluid3::none;
  if (action == "calc_fluid_systemmat_and_residual")
    act = XFluid3::calc_fluid_systemmat_and_residual;
  else if (action == "calc_linear_fluid")
    act = XFluid3::calc_linear_fluid;
  else if (action == "calc_fluid_stationary_systemmat_and_residual")
    act = XFluid3::calc_fluid_stationary_systemmat_and_residual;  
  else if (action == "calc_fluid_beltrami_error")
    act = XFluid3::calc_fluid_beltrami_error;
  else if (action == "store_xfem_info")
    act = XFluid3::store_xfem_info;
  else if (action == "get_density")
    act = XFluid3::get_density;
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
    dserror("looking for stab action (%s) not contained in map",action.c_str());
  }
  return act;
}


 /*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::XFluid3::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    std::vector<int>&         lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix&,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector&,
                                    Epetra_SerialDenseVector&)
{
  // get the action required
  const string action = params.get<string>("action","none");
  const DRT::ELEMENTS::XFluid3::ActionType act = convertStringToActionType(action);

  // get the material
  const Teuchos::RCP<MAT::Material> mat = Material();
  if (mat->MaterialType()!=m_fluid)
    dserror("newtonian fluid material expected but got type %d", mat->MaterialType());

  const MATERIAL* actmat = static_cast<MAT::NewtonianFluid*>(mat.get())->MaterialData();

  switch(act)
  {
      case get_density:
      {
        // This is a very poor way to transport the density to the
        // outside world. Is there a better one?
        params.set("density", actmat->m.fluid->density);
        break;
      }
  
    //--------------------------------------------------
    //--------------------------------------------------
    // the standard one-step-theta implementation
    //--------------------------------------------------
    //--------------------------------------------------
      case calc_fluid_systemmat_and_residual:
      {
        // do no calculation, if not needed
        if (lm.empty())
            break;
        
        // need current velocity/pressure and history vector
        Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
        if (velnp==null)
            dserror("Cannot get state vector 'velnp'");
        Teuchos::RCP<const Epetra_Vector> hist  = discretization.GetState("hist");
        if (hist==null)
            dserror("Cannot get state vectors 'hist'");

        // extract local values from the global vectors
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);
        vector<double> myhist(lm.size());
        DRT::UTILS::ExtractMyValues(*hist,myhist,lm);

        if (is_ale_)
        {
            dserror("No ALE support within instationary fluid solver.");
        }

        // get control parameter
        const double time = params.get<double>("total time",-1.0);

        const bool newton = params.get<bool>("include reactive terms for linearisation",false);

        const bool pstab  = true;
        const bool supg   = true;
        const bool cstab  = true;

        // One-step-Theta: timefac = theta*dt
        // BDF2:           timefac = 2/3 * dt
        const double timefac = params.get<double>("thsl",-1.0);
        if (timefac < 0.0)
            dserror("No thsl supplied");
        
        const Teuchos::RCP<Epetra_Vector> ivelcol = params.get<Teuchos::RCP<Epetra_Vector> >("interface velocity",null);
        if (ivelcol==null)
            dserror("Cannot get interface velocity from parameters");
        
        const Teuchos::RCP<Epetra_Vector> iforcecol = params.get<Teuchos::RCP<Epetra_Vector> >("interface force",null);
        if (iforcecol==null)
            dserror("Cannot get interface force from parameters");

        const XFEM::AssemblyType assembly_type = CheckForStandardEnrichmentsOnly(
                eleDofManager_, NumNode(), NodeIds());
        
        //--------------------------------------------------
        // calculate element coefficient matrix and rhs
        //--------------------------------------------------
        XFLUID::callSysmat4(assembly_type,
                this, ih_, eleDofManager_, myvelnp, myhist, ivelcol, iforcecol, elemat1, elevec1,
                actmat, time, timefac, newton, pstab, supg, cstab, true);

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
          RefCountPtr<const Epetra_Vector> vel_pre_np = discretization.GetState("u and p at time n+1 (converged)");
          if (vel_pre_np==null)
              dserror("Cannot get state vectors 'velnp'");

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
      case calc_fluid_stationary_systemmat_and_residual:
      {
          // do no calculation, if not needed
          if (lm.empty())
              break;
          
          // need current velocity/pressure 
          Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
          if (velnp==null)
              dserror("Cannot get state vector 'velnp'");

          // extract local values from the global vector
          vector<double> locval(lm.size());
          DRT::UTILS::ExtractMyValues(*velnp,locval,lm);
          vector<double> locval_hist(lm.size(),0.0); // zero history vector
          
          const Teuchos::RCP<Epetra_Vector> ivelcol = params.get<Teuchos::RCP<Epetra_Vector> >("interface velocity",null);
          if (ivelcol==null)
              dserror("Cannot get interface velocity from parameters");
          
          const Teuchos::RCP<Epetra_Vector> iforcecol = params.get<Teuchos::RCP<Epetra_Vector> >("interface force",null);
          if (iforcecol==null)
              dserror("Cannot get interface force from parameters");
          
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
          const bool cstab  = true;

          const XFEM::AssemblyType assembly_type = CheckForStandardEnrichmentsOnly(
                  eleDofManager_, NumNode(), NodeIds());
          
#if 0
          const XFEM::BoundaryIntCells&  boundaryIntCells(ih_->GetBoundaryIntCells(this->Id()));
          if ((assembly_type == XFEM::xfem_assembly) and (boundaryIntCells.size() > 0))
          {
              const int entry = 4; // line in stiffness matrix to compare
              const double disturbance = 1.0e-4;

              // initialize locval
              for (unsigned i = 0;i < locval.size(); ++i)
              {
                  locval[i] = 0.0;
                  locval_hist[i] = 0.0;
              }
              // R_0
              // calculate element coefficient matrix and rhs
              XFLUID::callSysmat4(assembly_type,
                      this, ih_, eleDofManager_, locval, locval_hist, ivelcol, iforcecol, estif, eforce,
                      actmat, pseudotime, 1.0, newton, pstab, supg, cstab, false);

              blitz::Array<double, 1> eforce_0(locval.size());
              for (unsigned i = 0;i < locval.size(); ++i)
              {
                  eforce_0(i) = eforce(i);
              }
              
              // create disturbed vector
              vector<double> locval_disturbed(locval.size());
              for (unsigned i = 0;i < locval.size(); ++i)
              {
                  if (i == entry)
                  {
                      locval_disturbed[i] = locval[i] + disturbance;
                  }
                  else
                  {
                      locval_disturbed[i] = locval[i];
                  }
                  std::cout << locval[i] <<  " " << locval_disturbed[i] << endl;
              }
              

              // R_0+dx
              // calculate element coefficient matrix and rhs
              XFLUID::callSysmat4(assembly_type,
                      this, ih_, eleDofManager_, locval_disturbed, locval_hist, ivelcol, iforcecol, estif, eforce,
                      actmat, pseudotime, 1.0, newton, pstab, supg, cstab, false);

              
              
              // compare
              std::cout << "sekante" << endl;
              for (int i = 0;i < locval.size(); ++i)
              {
                  //cout << i << endl;
                  const double matrixentry = (eforce_0(i) - eforce(i))/disturbance;
                  printf("should be %+12.8E, is %+12.8E, factor = %5.2f, is %+12.8E, factor = %5.2f\n", matrixentry, estif(i, entry), estif(i, entry)/matrixentry, estif(entry,i), estif(entry,i)/matrixentry);
                  //cout << "should be: " << std::scientific << matrixentry << ", is: " << estif(entry, i) << " " << estif(i, entry) << endl;                
              }
              
              exit(0);
          }
          else
#endif
          {
          // calculate element coefficient matrix and rhs
          XFLUID::callSysmat4(assembly_type,
                  this, ih_, eleDofManager_, locval, locval_hist, ivelcol, iforcecol, elemat1, elevec1,
                  actmat, pseudotime, 1.0, newton, pstab, supg, cstab, false);
          }
          break;
      }
      case store_xfem_info:
      {
          // get access to global dofman
          const Teuchos::RCP<XFEM::DofManager> globaldofman = params.get< Teuchos::RCP< XFEM::DofManager > >("dofmanager",null);
          if (globaldofman == null)
            dserror("nope, I need a DofManager!");
          
          // create local copy of information about dofs
          const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz(XFLUID::getElementAnsatz(this->Shape()));
          
          eleDofManager_ = globaldofman->constructElementDofManager(*this, element_ansatz);
          
          // store pointer to interface handle
          ih_ = params.get< Teuchos::RCP< XFEM::InterfaceHandle > >("interfacehandle",null);
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
                                           std::vector<int>&         lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  return 0;
}

// get optimal gaussrule for discretization type
DRT::UTILS::GaussRule3D DRT::ELEMENTS::XFluid3::getOptimalGaussrule(const DiscretizationType& distype)
{
  DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule3D_undefined;
    switch (distype)
    {
    case hex8:
        rule = DRT::UTILS::intrule_hex_8point;
        break;
    case hex20: case hex27:
        rule = DRT::UTILS::intrule_hex_27point;
        break;
    case tet4:
        rule = DRT::UTILS::intrule_tet_4point;
        break;
    case tet10:
        rule = DRT::UTILS::intrule_tet_5point;
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
    std::vector<double>&      evelnp,
    std::vector<double>&      eprenp,
    const struct _MATERIAL*   material,
    ParameterList&            params
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
  const DRT::UTILS::GaussRule3D gaussrule = getOptimalGaussrule(distype);
  const DRT::UTILS::IntegrationPoints3D  intpoints = getIntegrationPoints3D(gaussrule);

  // start loop over integration points
  for (int iquad=0;iquad<intpoints.nquad;iquad++)
  {
    // declaration of gauss point variables
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];
    const double e3 = intpoints.qxg[iquad][2];
    DRT::UTILS::shape_function_3D(funct,e1,e2,e3,distype);
    DRT::UTILS::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

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
        printf("NEGATIVE JACOBIAN DETERMINANT: %f\n", det);
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


/*----------------------------------------------------------------------*
 |  init the element (public)                                mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::XFluid3Register::Initialize(DRT::Discretization&)
{
  return 0;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
