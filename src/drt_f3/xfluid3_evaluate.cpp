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

#include <Epetra_SerialDenseSolver.h>

#include "xfluid3.H"
#include "xfluid3_sysmat.H"
#include "xfluid3_interpolation.H"

#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/xdofmapcreation.H"
#include "../drt_xfem/enrichment_utils.H"


// converts a string into an Action for this element
DRT::ELEMENTS::XFluid3::ActionType DRT::ELEMENTS::XFluid3::convertStringToActionType(
              const string& action) const
{
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
  const string action(params.get<string>("action","none"));
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
      case calc_fluid_systemmat_and_residual:
      {
        // do no calculation, if not needed
        if (lm.empty())
            break;
        
        // extract local values from the global vectors
        DRT::ELEMENTS::XFluid3::MyState mystate;
        mystate.instationary = true;
        DRT::UTILS::ExtractMyValues(*discretization.GetState("velnp"),mystate.velnp,lm);
        DRT::UTILS::ExtractMyValues(*discretization.GetState("veln") ,mystate.veln ,lm);
        DRT::UTILS::ExtractMyValues(*discretization.GetState("velnm"),mystate.velnm,lm);
        DRT::UTILS::ExtractMyValues(*discretization.GetState("accn") ,mystate.accn ,lm);

        if (is_ale_)
        {
            dserror("No ALE support within instationary fluid solver.");
        }

        const bool newton = params.get<bool>("include reactive terms for linearisation",false);

        const bool pstab  = true;
        const bool supg   = true;
        const bool cstab  = true;

        // time integration factors
        const FLUID_TIMEINTTYPE timealgo = params.get<FLUID_TIMEINTTYPE>("timealgo");
        const double            dt       = params.get<double>("dt");
        const double            theta    = params.get<double>("theta");
        
        const Teuchos::RCP<const Epetra_Vector> ivelcol = params.get<Teuchos::RCP<const Epetra_Vector> >("interface velocity");
        const Teuchos::RCP<Epetra_Vector> iforcecol = params.get<Teuchos::RCP<Epetra_Vector> >("interface force");

        const bool ifaceForceContribution = discretization.ElementRowMap()->MyGID(this->Id());
        
        if (not params.get<bool>("DLM_condensation") or not ih_->ElementIntersected(Id())) // integrate and assemble all unknowns
        {
          const XFEM::AssemblyType assembly_type = CheckForStandardEnrichmentsOnly(
                  *eleDofManager_, NumNode(), NodeIds());
          
          // calculate element coefficient matrix and rhs
          XFLUID::callSysmat4(assembly_type,
                  this, ih_, *eleDofManager_, mystate, ivelcol, iforcecol, elemat1, elevec1,
                  actmat, timealgo, dt, theta, newton, pstab, supg, cstab, true, ifaceForceContribution);
        }
        else // create bigger element matrix and vector, assemble, condense and copy to small matrix provided by discretization
        {
          // sanity checks
          if (eleDofManager_->NumNodeDof() != eleDofManager_uncondensed_->NumNodeDof())
            dserror("NumNodeDof mismatch");
          if (eleDofManager_->NumElemDof() != 0)
            dserror("NumElemDof not 0");
          if (eleDofManager_uncondensed_->NumElemDof() == 0)
            dserror("uncondensed NumElemDof == 0");
          
          // stress update
          UpdateOldDLMAndDLMRHS(discretization, lm, mystate);
          
          // create uncondensed element matrix and vector
          const int numdof_uncond = eleDofManager_uncondensed_->NumDofElemAndNode();
          Epetra_SerialDenseMatrix elemat1_uncond(numdof_uncond,numdof_uncond);
          Epetra_SerialDenseVector elevec1_uncond(numdof_uncond);
          
          const XFEM::AssemblyType assembly_type = CheckForStandardEnrichmentsOnly(
                  *eleDofManager_uncondensed_, NumNode(), NodeIds());
          
          // calculate element coefficient matrix and rhs
          XFLUID::callSysmat4(assembly_type,
                  this, ih_, *eleDofManager_uncondensed_, mystate, ivelcol, iforcecol, elemat1_uncond, elevec1_uncond,
                  actmat, timealgo, dt, theta, newton, pstab, supg, cstab, true, ifaceForceContribution);
          
          // condensation
          CondenseDLMAndStoreOldIterationStep(elemat1_uncond, elevec1_uncond, elemat1, elevec1);
        }
        
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

          // extract local values from the global vector
          DRT::ELEMENTS::XFluid3::MyState mystate;
          mystate.instationary = false;
          DRT::UTILS::ExtractMyValues(*discretization.GetState("velnp"),mystate.velnp,lm);
          
          const Teuchos::RCP<const Epetra_Vector> ivelcol = params.get<Teuchos::RCP<const Epetra_Vector> >("interface velocity");
          const Teuchos::RCP<Epetra_Vector> iforcecol = params.get<Teuchos::RCP<Epetra_Vector> >("interface force");
          
          if (is_ale_)
          {
              dserror("No ALE support within stationary fluid solver.");
          }
          
          // time integration factors
          const FLUID_TIMEINTTYPE timealgo = params.get<FLUID_TIMEINTTYPE>("timealgo");
          const double            dt       = 1.0;
          const double            theta    = 1.0;

          const bool newton = params.get<bool>("include reactive terms for linearisation",false);
          const bool pstab  = true;
          const bool supg   = true;
          const bool cstab  = true;

          const bool ifaceForceContribution = discretization.ElementRowMap()->MyGID(this->Id());
          
          if (not params.get<bool>("DLM_condensation") or not ih_->ElementIntersected(Id())) // integrate and assemble all unknowns
          {
            const XFEM::AssemblyType assembly_type = CheckForStandardEnrichmentsOnly(
                    *eleDofManager_, NumNode(), NodeIds());
            
            // calculate element coefficient matrix and rhs
            XFLUID::callSysmat4(assembly_type,
                    this, ih_, *eleDofManager_, mystate, ivelcol, iforcecol, elemat1, elevec1,
                    actmat, timealgo, dt, theta, newton, pstab, supg, cstab, mystate.instationary, ifaceForceContribution);
            
          }
          else // create bigger element matrix and vector, assemble, condense and copy to small matrix provided by discretization
          {
            // sanity checks
            if (eleDofManager_->NumNodeDof() != eleDofManager_uncondensed_->NumNodeDof())
              dserror("NumNodeDof mismatch");
            if (eleDofManager_->NumElemDof() != 0)
              dserror("NumElemDof not 0");
//            if (eleDofManager_uncondensed_->NumElemDof() == 0)
//            {
//              const double boundarysize = XFEM::BoundaryCoverageRatio(*this,*ih_);
//              cout << "boundarysize = " << boundarysize << endl;
////              dserror("NumElemDof uncondensed == 0");
//            }
            
            // stress update
            UpdateOldDLMAndDLMRHS(discretization, lm, mystate);
            
            // create uncondensed element matrix and vector
            const int numdof_uncond = eleDofManager_uncondensed_->NumDofElemAndNode();
            Epetra_SerialDenseMatrix elemat1_uncond(numdof_uncond,numdof_uncond);
            Epetra_SerialDenseVector elevec1_uncond(numdof_uncond);
            
            const XFEM::AssemblyType assembly_type = CheckForStandardEnrichmentsOnly(
                    *eleDofManager_uncondensed_, NumNode(), NodeIds());
            
            // calculate element coefficient matrix and rhs
            XFLUID::callSysmat4(assembly_type,
                    this, ih_, *eleDofManager_uncondensed_, mystate, ivelcol, iforcecol, elemat1_uncond, elevec1_uncond,
                    actmat, timealgo, dt, theta, newton, pstab, supg, cstab, mystate.instationary, ifaceForceContribution);
            
            // condensation
            CondenseDLMAndStoreOldIterationStep(elemat1_uncond, elevec1_uncond, elemat1, elevec1);
          }
          
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
          break;
      }
      case store_xfem_info:
      {
          // store pointer to interface handle
          ih_ = params.get< Teuchos::RCP< XFEM::InterfaceHandleXFSI > >("interfacehandle",null);
        
          // get access to global dofman
          const Teuchos::RCP<XFEM::DofManager> globaldofman = params.get< Teuchos::RCP< XFEM::DofManager > >("dofmanager");
          
          const bool DLM_condensation = params.get<bool>("DLM_condensation");
          {
            // create local copy of information about dofs
            map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz;
            if (not DLM_condensation)
              element_ansatz = XFLUID::getElementAnsatz(this->Shape());
          
            eleDofManager_ = rcp(new XFEM::ElementDofManager(*this, element_ansatz, *globaldofman));
          }
          
          if (ih_->ElementIntersected(Id()))
          {
            std::set<XFEM::FieldEnr> enrfieldset;
            
            const std::set<int> xlabelset(eleDofManager_->getUniqueEnrichmentLabels());
            // loop condition labels
            for(std::set<int>::const_iterator labeliter = xlabelset.begin(); labeliter!=xlabelset.end(); ++labeliter)
            {
              const int label = *labeliter;

              // for surface with label, loop my col elements and add void enrichments to each elements member nodes
              if (ih_->ElementHasLabel(this->Id(), label))
              {
                XFEM::ApplyElementEnrichments(this, *ih_, label, false, enrfieldset);
              }
            };
            
            // create local copy of information about dofs
            const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz(XFLUID::getElementAnsatz(this->Shape()));
            
            // nodal dofs for ele
            eleDofManager_uncondensed_ = 
              rcp(new XFEM::ElementDofManager(*this, *eleDofManager_->getNodalDofSet(), enrfieldset, element_ansatz));
            
            const int nd = eleDofManager_uncondensed_->NumNodeDof();
            const int na = eleDofManager_uncondensed_->NumElemDof();
            
            DLM_info_ = Teuchos::rcp(new DLMInfo(nd,na));
          }
          else
          {
            eleDofManager_uncondensed_ = Teuchos::null;
            DLM_info_ = Teuchos::null;
          }
          
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


/*---------------------------------------------------------------------*
 *---------------------------------------------------------------------*/
void DRT::ELEMENTS::XFluid3::UpdateOldDLMAndDLMRHS(
    const DRT::Discretization&      discretization,
    const std::vector<int>&         lm,
    MyState&                        mystate
    ) const
{
  const int nd = eleDofManager_uncondensed_->NumNodeDof();
  const int na = eleDofManager_uncondensed_->NumElemDof();
  const int numdof_uncond = eleDofManager_uncondensed_->NumDofElemAndNode();
  
  // add Kda . res_d to feas
  // new alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas
  
  vector<double> res_d(lm.size());
  DRT::UTILS::ExtractMyValues(*discretization.GetState("nodal residual"),res_d,lm);
  
  //update old iteration residuum of the stresses
  for (int i=0;i<na;i++)
    for (int j=0;j<nd;j++)
      DLM_info_->oldfa_(i) += DLM_info_->oldKad_(i,j)*res_d[j];
  
  // compute element stresses
  for (int i=0;i<na;i++)
    for (int j=0;j<na;j++)
      DLM_info_->stressdofs_(i) -= DLM_info_->oldKaainv_(i,j)*DLM_info_->oldfa_(j);
  
  // increase size of element vector (old values stay and zeros are added)
  mystate.velnp.resize(numdof_uncond,0.0);
  mystate.veln .resize(numdof_uncond,0.0);
  mystate.velnm.resize(numdof_uncond,0.0);
  mystate.accn .resize(numdof_uncond,0.0);
  for (int i=0;i<na;i++)
    mystate.velnp[i+nd] = DLM_info_->stressdofs_(i);
}

/*---------------------------------------------------------------------*
 *---------------------------------------------------------------------*/
void DRT::ELEMENTS::XFluid3::CondenseDLMAndStoreOldIterationStep(
    const Epetra_SerialDenseMatrix& elemat1_uncond,
    const Epetra_SerialDenseVector& elevec1_uncond,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseVector& elevec1    
) const
{

  const int nd = eleDofManager_uncondensed_->NumNodeDof();
  const int na = eleDofManager_uncondensed_->NumElemDof();
  
  // copy nodal dof entries
  for (int i = 0; i < nd; ++i)
  {
    elevec1(i) = elevec1_uncond(i);
    for (int j = 0; j < nd; ++j)
    {
      elemat1(i,j) = elemat1_uncond(i,j);
    }
  }
  
  if (na > 0)
  {
    // note: the full (u,p,sigma) matrix is unsymmetric, hence we need both, rectangular matrizes Kda and Kad
    Epetra_SerialDenseMatrix Kda(nd,na);
    Epetra_SerialDenseMatrix Kaa(na,na);
    Epetra_SerialDenseMatrix Kad(na,nd);
    Epetra_SerialDenseVector fa(na);

    for (int i=0;i<nd;i++)
      for (int j=0;j<na;j++)
        Kda(i,j) = elemat1_uncond(i   ,j+nd);
    
    for (int i=0;i<na;i++)
      for (int j=0;j<na;j++)
        Kaa(i,j) = elemat1_uncond(i+nd,j+nd);

    for (int i=0;i<na;i++)
      for (int j=0;j<nd;j++)
        Kad(i,j) = elemat1_uncond(i+nd,j   );
    
    for (int i=0;i<na;i++)
      fa(i) = elevec1_uncond(i+nd);
    
    
    // DLM-stiffness matrix is: Kdd - Kda . Kaa^-1 . Kad
    // DLM-internal force is: fint - Kda . Kaa^-1 . feas
    
    // we need the inverse of Kaa
    Epetra_SerialDenseSolver solve_for_inverseKaa;
    solve_for_inverseKaa.SetMatrix(Kaa);
    solve_for_inverseKaa.Invert();

    {
      LINALG::SerialDenseMatrix KdaKaainv(nd,na); // temporary Kda.Kaa^{-1}
      KdaKaainv.Zero();
      
      //LINALG::DENSEFUNCTIONS::multiply<NUMDOF_SOH8,soh8_easfull,soh8_easfull>(KdaKaa, Kda, Kaa);
      for (int i=0;i<nd;i++)
        for (int j=0;j<na;j++)
          for (int k=0;k<na;k++)
            KdaKaainv(i,j) += Kda(i,k)*Kaa(k,j);
//                cout << endl << "KdaKaainv" << endl;
//                cout << KdaKaainv << endl;
      
      //LINALG::DENSEFUNCTIONS::multiply<NUMDOF_SOH8,soh8_easfull,NUMDOF_SOH8>(1.0, stiffmatrix->A(), -1.0, KdaKaa.A(), Kda.A());
      for (int i=0;i<nd;i++)
        for (int j=0;j<nd;j++)
          for (int k=0;k<na;k++)
            elemat1(i,j) -= KdaKaainv(i,k)*Kad(k,j);
      
      //LINALG::DENSEFUNCTIONS::multiply<NUMDOF_SOH8,soh8_easfull,1>(1.0, force->A(), -1.0, KdaKaa.A(), feas.A());
      for (int i=0;i<nd;i++)
        for (int j=0;j<na;j++)
          elevec1(i) -= KdaKaainv(i,j)*fa(j);
    }
   
    {
      // store current DLM data in iteration history
      for (int i=0; i<na; ++i)
        for (int j=0; j<na; ++j)
          DLM_info_->oldKaainv_(i,j) = Kaa(i,j);
      
      for (int i=0; i<na; ++i)
        for (int j=0; j<nd; ++j)
          DLM_info_->oldKad_(i,j) = Kad(i,j);

      for (int i=0; i<na; ++i)
        DLM_info_->oldfa_(i) = fa(i);
    }
  }
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
