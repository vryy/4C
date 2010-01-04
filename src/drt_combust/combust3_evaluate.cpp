/*!
\file combust3_evaluate.cpp
\brief

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include <Teuchos_TimeMonitor.hpp>

#include "combust3.H"
#include "combust3_sysmat.H"
#include "combust3_interpolation.H"

#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_mat/matlist.H"
#include "../drt_f3/fluid3_stabilization.H"
#include "../drt_inpar/inpar_fluid.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>


// converts a string into an Action for this element
DRT::ELEMENTS::Combust3::ActionType DRT::ELEMENTS::Combust3::convertStringToActionType(
              const string& action) const
{
  DRT::ELEMENTS::Combust3::ActionType act = Combust3::none;
  if (action == "calc_fluid_systemmat_and_residual")
    act = Combust3::calc_fluid_systemmat_and_residual;
  else if (action == "calc_linear_fluid")
    act = Combust3::calc_linear_fluid;
  else if (action == "calc_fluid_stationary_systemmat_and_residual")
    act = Combust3::calc_fluid_stationary_systemmat_and_residual;
  else if (action == "calc_fluid_beltrami_error")
    act = Combust3::calc_fluid_beltrami_error;
  else if (action == "calc_turbulence_statistics")
    act = Combust3::calc_turbulence_statistics;
  else if (action == "calc_fluid_box_filter")
    act = Combust3::calc_fluid_box_filter;
  else if (action == "calc_smagorinsky_const")
    act = Combust3::calc_smagorinsky_const;
  else if (action == "store_xfem_info")
    act = Combust3::store_xfem_info;
  else if (action == "get_density")
    act = Combust3::get_density;
  else if (action == "reset")
    act = Combust3::reset;
  else if (action == "set_output_mode")
    act = Combust3::set_output_mode;
  else
    dserror("Unknown type of action for Combust3");
  return act;
}

/*----------------------------------------------------------------------*
 // converts a string into an stabilisation action for this element
 //                                                          gammi 08/07
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Combust3::StabilisationAction DRT::ELEMENTS::Combust3::ConvertStringToStabAction(
  const string& action) const
{
  DRT::ELEMENTS::Combust3::StabilisationAction act = stabaction_unspecified;

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
 |  evaluate the element (public)                           g.bau 03/07 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Combust3::Evaluate(ParameterList& params,
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
  const DRT::ELEMENTS::Combust3::ActionType act = convertStringToActionType(action);

  // get the list of materials
  const Teuchos::RCP<MAT::Material> material = Material();

  switch(act)
  {
    case get_density:
    {
      std::cout << "Warning! The density is set to 1.0!" << std::endl;
      params.set("density", 1.0);
      break;
    }
    case reset:
    {
      // reset all information and make element unusable (e.g. it can't answer the numdof question anymore)
      // this way, one can see, if all information are generated correctly or whether something is left
      // from the last nonlinear iteration
      eleDofManager_ = Teuchos::null;
      ih_ = Teuchos::null;
      break;
    }
    case set_output_mode:
    {
      output_mode_ = true;
      // reset dof managers if present
      eleDofManager_ = Teuchos::null;
      ih_ = Teuchos::null;
      break;
    }
    case calc_fluid_systemmat_and_residual:
    {
      TEUCHOS_FUNC_TIME_MONITOR("COMBUST3 - evaluate - calc_fluid_systemmat_and_residual");

      // do no calculation, if not needed
      if (lm.empty())
        break;

      DRT::ELEMENTS::Combust3::MyState mystate;
      mystate.instationary = true;
      DRT::UTILS::ExtractMyValues(*discretization.GetState("velnp"),mystate.velnp,lm);
      DRT::UTILS::ExtractMyValues(*discretization.GetState("veln") ,mystate.veln ,lm);
      DRT::UTILS::ExtractMyValues(*discretization.GetState("velnm"),mystate.velnm,lm);
      DRT::UTILS::ExtractMyValues(*discretization.GetState("accn") ,mystate.accn ,lm);

      // get pointer to vector holding G-function values at the fluid nodes
      const Teuchos::RCP<Epetra_Vector> phinp = ih_->FlameFront()->Phinp();
#ifdef DEBUG
      // check if this element is the first element on this processor
      // remark:
      // The SameAs-operation requires MPI communication between processors. Therefore it can only
      // be performed once (at the beginning) on each processor. Otherwise some processors would
      // wait to receive MPI information, but would never get it, because some processores are
      // already done with their element loop. This will cause a mean parallel bug!   henke 11.08.09
      if(this->Id() == discretization.lRowElement(0)->Id())
      {
        // get map of this vector
        const Epetra_BlockMap& phimap = phinp->Map();
        // check, whether this map is still identical with the current node map in the discretization
        if (not phimap.SameAs(*discretization.NodeColMap())) dserror("node column map has changed!");
      }
#endif

      // extract G-function values to element level
      DRT::UTILS::ExtractMyNodeBasedValues(this, mystate.phinp, *phinp);

      const bool newton = params.get<bool>("include reactive terms for linearisation",false);

      const INPAR::COMBUST::CombustionType combusttype = params.get<INPAR::COMBUST::CombustionType>("combusttype");
      const double flamespeed = params.get<double>("flamespeed");
      const double nitschevel = params.get<double>("nitschevel");
      const double nitschepres = params.get<double>("nitschepres");

      // stabilization terms
      const bool pstab = true;
      const bool supg  = true;
      const bool cstab = true;
      // stabilization parameters
      const INPAR::FLUID::TauType tautype = Teuchos::getIntegralValue<INPAR::FLUID::TauType>(params.sublist("STABILIZATION"),"TAUTYPE");
      // check if stabilization parameter definition can be handled by combust3 element
      if (!(tautype == INPAR::FLUID::tautype_franca_barrenechea_valentin_wall or
          tautype == INPAR::FLUID::tautype_bazilevs))
        dserror("unknown type of stabilization parameter definition");

      // time integration parameters
      const FLUID_TIMEINTTYPE timealgo = params.get<FLUID_TIMEINTTYPE>("timealgo");
      const double            dt       = params.get<double>("dt");
      const double            theta    = params.get<double>("theta");

      const XFEM::AssemblyType assembly_type = XFEM::ComputeAssemblyType(
              *eleDofManager_, NumNode(), NodeIds());

      // calculate element coefficient matrix and rhs
      COMBUST::callSysmat(assembly_type,
        this, ih_, *eleDofManager_, mystate, elemat1, elevec1,
        material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, mystate.instationary,
        combusttype, flamespeed, nitschevel, nitschepres);
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
        std::vector<double> my_vel_pre_np(lm.size());
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
        f3_int_beltrami_err(myvelnp,myprenp,material,params);
      }
    }
    break;
    case calc_fluid_stationary_systemmat_and_residual:
    {
      TEUCHOS_FUNC_TIME_MONITOR("COMBUST3 - evaluate - calc_fluid_stationary_systemmat_and_residual");
      // do no calculation, if not needed
      if (lm.empty())
        break;

      // extract local values from the global vector
      DRT::ELEMENTS::Combust3::MyState mystate;
      mystate.instationary = false;
      DRT::UTILS::ExtractMyValues(*discretization.GetState("velnp"),mystate.velnp,lm);

      // get pointer to vector holding G-function values at the fluid nodes
      const Teuchos::RCP<Epetra_Vector> phinp = ih_->FlameFront()->Phinp();
#ifdef DEBUG
      // check if this element is the first element on this processor
      // remark:
      // The SameAs-operation requires MPI communication between processors. Therefore it can only
      // be performed once (at the beginning) on each processor. Otherwise some processors would
      // wait to receive MPI information, but would never get it, because some processores are
      // already done with their element loop. This will cause a mean parallel bug!   henke 11.08.09
      if(this->Id() == discretization.lRowElement(0)->Id())
      {
        // get map of this vector
        const Epetra_BlockMap& phimap = phinp->Map();
        // check, whether this map is still identical with the current node map in the discretization
        if (not phimap.SameAs(*discretization.NodeColMap())) dserror("node column map has changed!");
      }
#endif

      // extract G-function values to element level (used kink enrichment)
      DRT::UTILS::ExtractMyNodeBasedValues(this, mystate.phinp, *phinp);

      const INPAR::COMBUST::CombustionType combusttype = params.get<INPAR::COMBUST::CombustionType>("combusttype");
      const double flamespeed = params.get<double>("flamespeed");
      const double nitschevel = params.get<double>("nitschevel");
      const double nitschepres = params.get<double>("nitschepres");

      // time integration factors
      const FLUID_TIMEINTTYPE timealgo = params.get<FLUID_TIMEINTTYPE>("timealgo");
      const double            dt       = 1.0;
      const double            theta    = 1.0;

      const bool newton = params.get<bool>("include reactive terms for linearisation",false);
      // stabilization parameters
      const INPAR::FLUID::TauType tautype = Teuchos::getIntegralValue<INPAR::FLUID::TauType>(params.sublist("STABILIZATION"),"TAUTYPE");
      // check if stabilization parameter definition can be handled by combust3 element
      if (!(tautype == INPAR::FLUID::tautype_franca_barrenechea_valentin_wall or
          tautype == INPAR::FLUID::tautype_bazilevs))
        dserror("unknown type of stabilization parameter definition");
      const bool pstab  = true;
      const bool supg   = true;
      const bool cstab  = true;

      const XFEM::AssemblyType assembly_type = XFEM::ComputeAssemblyType(
              *eleDofManager_, NumNode(), NodeIds());

      {
        // calculate element coefficient matrix and rhs
        COMBUST::callSysmat(assembly_type,
                 this, ih_, *eleDofManager_, mystate, elemat1, elevec1,
                 material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, mystate.instationary,
                 combusttype, flamespeed, nitschevel, nitschepres);
      }
#if 0
          const XFEM::BoundaryIntCells&  boundaryIntCells(ih_->GetBoundaryIntCells(this->Id()));
          if ((assembly_type == XFEM::xfem_assembly) and (not boundaryIntCells.empty()))
          {
              const int entry = 4; // line in stiffness matrix to compare
              const double disturbance = 1.0e-4;

              // initialize locval
              for (std::size_t i = 0;i < locval.size(); ++i)
              {
                  locval[i] = 0.0;
                  locval_hist[i] = 0.0;
              }
              // R_0
              // calculate element coefficient matrix and rhs
              XFLUID::callSysmat4(assembly_type,
                      this, ih_, eleDofManager_, locval, locval_hist, ivelcol, iforcecol, estif, eforce,
                      mat, pseudotime, 1.0, newton, pstab, supg, cstab, false);

              LINALG::SerialDensevector eforce_0(locval.size());
              for (std::size_t i = 0;i < locval.size(); ++i)
              {
                  eforce_0(i) = eforce(i);
              }

              // create disturbed vector
              vector<double> locval_disturbed(locval.size());
              for (std::size_t i = 0;i < locval.size(); ++i)
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
                      mat, pseudotime, 1.0, newton, pstab, supg, cstab, false);



              // compare
              std::cout << "sekante" << endl;
              for (std::size_t i = 0;i < locval.size(); ++i)
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
    }
    break;
    case store_xfem_info:
    {
      TEUCHOS_FUNC_TIME_MONITOR("COMBUST3 - evaluate - store_xfem_info");

      // now the element can answer how many (XFEM) dofs it has
      output_mode_ = false;

      // store pointer to interface handle
      ih_ = params.get< Teuchos::RCP< COMBUST::InterfaceHandleCombust > >("interfacehandle",Teuchos::null);

      //---------------------------------------------------
      // find out whether an element is intersected or not
      //---------------------------------------------------
      // initialization call of fluid time integration scheme will end up here:
      // The initial flame front has not been incorporated into the fluid field -> no XFEM dofs, yet!
      if (ih_->FlameFront() == Teuchos::null)
      {
        this->intersected_ = false;
      }
      else // regular call
      {
        // get vector of integration cells for this element
        std::size_t numcells= ih_->GetNumDomainIntCells(this);
        if (numcells > 1)
        {
          // more than one integration cell -> element intersected
          this->intersected_ = true;
          // std::cout << "Mehrere DomainIntCells für Element "<<  this->Id() << " gefunden " << std::endl;
        }
        else if (numcells == 1)
        {
          // only one integration cell -> element not intersected
          this->intersected_ = false;
          // std::cout << "Eine DomainIntCell für Element "<<  this->Id() << " gefunden " << std::endl;
        }
        else // numcells = 0 or negative number
        {
          // no integration cell -> impossible, something went wrong!
          dserror ("There are no DomainIntCells for element %d ", this->Id());
        }
      }

      // get access to global dofman
      const Teuchos::RCP<XFEM::DofManager> globaldofman = params.get< Teuchos::RCP< XFEM::DofManager > >("dofmanager");

      // create local copy of information about dofs
      const COMBUST::CombustElementAnsatz elementAnsatz;
      eleDofManager_ = rcp(new XFEM::ElementDofManager(*this, elementAnsatz.getElementAnsatz(Shape()), *globaldofman));
    }
    break;
    default:
      dserror("Unknown type of action for Combust3");
  } // end of switch(act)

    return 0;
} // end of DRT::ELEMENTS::Fluid3::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                      gammi 04/07|
 |                                                                      |
 |  The function is just a dummy. For the fluid elements, the           |
 |  integration of the volume neumann (body forces) loads takes place   |
 |  in the element. We need it there for the stabilisation terms!       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Combust3::EvaluateNeumann(ParameterList& params,
                                            DRT::Discretization&      discretization,
                                            DRT::Condition&           condition,
                                            std::vector<int>&         lm,
                                            Epetra_SerialDenseVector& elevec1,
                                            Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}

// get optimal gaussrule for discretization type
DRT::UTILS::GaussRule3D DRT::ELEMENTS::Combust3::getOptimalGaussrule(const DiscretizationType& distype)
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
void DRT::ELEMENTS::Combust3::f3_int_beltrami_err(
    std::vector<double>&      evelnp,
    std::vector<double>&      eprenp,
    Teuchos::RCP<const MAT::Material> material,
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

  const double a      = M_PI/4.0;
  const double d      = M_PI/2.0;

  // get viscosity
  double  visc = 0.0;
  // just to be sure - actually this has already been checked in Combust3::Evaluate, before
  dsassert(material->MaterialType() == INPAR::MAT::m_matlist, "material is not of type m_matlist");
  const MAT::MatList* matlist = static_cast<const MAT::MatList*>(material.get());
  // use material MAT 3 (first in the material list) for Beltrami flow
  Teuchos::RCP<const MAT::Material> matptr = matlist->MaterialById(3);
  dsassert(matptr->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
  std::cout << "material MAT 3 is used for Beltrami flow" << std::endl;
  const MAT::NewtonianFluid* mat = static_cast<const MAT::NewtonianFluid*>(matptr.get());
  visc = mat->Viscosity();

  double         preint;
  vector<double> velint  (3);
  vector<double> xint    (3);

  vector<double> u       (3);

  double         deltap;
  vector<double> deltavel(3);

  // gaussian points
  const DRT::UTILS::GaussRule3D gaussrule = getOptimalGaussrule(distype);
  const DRT::UTILS::IntegrationPoints3D  intpoints(gaussrule);

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
    LINALG::Matrix<NSD,NSD>    xjm;

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
    const double det = xjm.Determinant();

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
int DRT::ELEMENTS::Combust3Register::Initialize(DRT::Discretization&)
{
  return 0;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
