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
#include <Teuchos_TimeMonitor.hpp>

#include "xfluid3.H"
#include "xfluid3_sysmat.H"
#include "xfluid3_interpolation.H"

#include "xfluid3_local_assembler.H"
#include "xfluid3_interpolation.H"
#include "xfluid3_utils.H"
#include "../drt_geometry/integrationcell_coordtrafo.H"
#include "../drt_xfem/physics.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/xdofmapcreation.H"
#include "../drt_xfem/enrichment_utils.H"


/*---------------------------------------------------------------------*
|  converts a string into an Action for this element                   |
*----------------------------------------------------------------------*/
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
  else if (action == "calc_fluid_projection_systemmat_and_residual")
    act = XFluid3::calc_fluid_projection_systemmat_and_residual;
  else if (action == "calc_fluid_beltrami_error")
    act = XFluid3::calc_fluid_beltrami_error;
  else if (action == "store_xfem_info")
    act = XFluid3::store_xfem_info;
  else if (action == "get_density")
    act = XFluid3::get_density;
  else if (action == "reset")
    act = XFluid3::reset;
  else if (action == "set_output_mode")
    act = XFluid3::set_output_mode;
  else if (action == "integrate_shape")
    act = XFluid3::integrate_shape;
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

static void SanityChecks(
    Teuchos::RCP<XFEM::ElementDofManager> eleDofManager,
    Teuchos::RCP<XFEM::ElementDofManager> eleDofManager_uncondensed
    )
{
  // sanity checks
  if (eleDofManager->NumNodeDof() != eleDofManager_uncondensed->NumNodeDof())
    dserror("NumNodeDof mismatch");
  if (eleDofManager->NumElemDof() != 0)
    dserror("NumElemDof not 0");
}


 /*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::XFluid3::Evaluate(ParameterList& params,
                                     DRT::Discretization&      discretization,
                                     std::vector<int>&         lm,
                                     Epetra_SerialDenseMatrix& elemat1,
                                     Epetra_SerialDenseMatrix& elemat2,
                                     Epetra_SerialDenseVector& elevec1,
                                     Epetra_SerialDenseVector& elevec2,
                                     Epetra_SerialDenseVector&)
{
  // get the action required
  const DRT::ELEMENTS::XFluid3::ActionType act =
      convertStringToActionType(params.get<std::string>("action"));

  // get the material
  const Teuchos::RCP<MAT::Material> mat = Material();
  if (mat->MaterialType()!=INPAR::MAT::m_fluid)
    dserror("newtonian fluid material expected but got type %d", mat->MaterialType());

  const MAT::NewtonianFluid* actmat = dynamic_cast<const MAT::NewtonianFluid*>(mat.get());

  switch(act)
  {
    case reset:
    {
      // reset all information and make element unusable (e.g. it can't answer the numdof question anymore)
      // this way, one can see, if all information are generated correctly or whether something is left
      // from the last nonlinear iteration
      eleDofManager_ = Teuchos::null;
      eleDofManager_uncondensed_ = Teuchos::null;
      ih_ = Teuchos::null;
      DLM_info_ = Teuchos::null;
      break;
    }
    case set_output_mode:
    {
      output_mode_ = params.get<bool>("output_mode");
      // reset dof managers if present
      eleDofManager_ = Teuchos::null;
      eleDofManager_uncondensed_ = Teuchos::null;
      ih_ = Teuchos::null;
      DLM_info_ = Teuchos::null;
      break;
    }
    case store_xfem_info:
    {
      // after this part the element can answer, how many DOFs it has
      output_mode_ = false;

      // store pointer to interface handle
      ih_ = params.get< Teuchos::RCP< XFEM::InterfaceHandleXFSI > >("interfacehandle");

      // get access to global dofman
      const Teuchos::RCP<XFEM::DofManager> globaldofman = params.get< Teuchos::RCP< XFEM::DofManager > >("dofmanager");

      const XFLUID::FluidElementAnsatz elementAnsatz;
      const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz_empty;
      const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz_filled(elementAnsatz.getElementAnsatz(this->Shape()));

      // always build the eledofman that fits to the global dofs
      // problem: tight connectivity to xdofmapcreation
      if (params.get<bool>("DLM_condensation"))
      {
        // assume no stress unknowns for the element
        eleDofManager_ = rcp(new XFEM::ElementDofManager(*this, element_ansatz_empty, *globaldofman));
      }
      else
      {
        // assume stress unknowns for the element
        eleDofManager_ = rcp(new XFEM::ElementDofManager(*this, element_ansatz_filled, *globaldofman));
      }

      // create an eledofman that has stress unknowns only for intersected elements
      // Note: condensation for unintersected elements is not handled, but also not needed
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
            const bool anothervoidenrichment_in_set = XFEM::EnrichmentInDofSet(XFEM::Enrichment::typeVoid, enrfieldset);
            if (not anothervoidenrichment_in_set)
            {
              XFEM::ApplyElementEnrichments(this, element_ansatz_filled, *ih_, label, XFEM::Enrichment::typeVoid, params.get<double>("boundaryRatioLimit"), enrfieldset);
            }
          }
        };

        int nd_old = -1;
        int na_old = -1;
        if (eleDofManager_uncondensed_ != Teuchos::null)
        {
          nd_old = eleDofManager_uncondensed_->NumNodeDof();
          na_old = eleDofManager_uncondensed_->NumElemDof();
        }

        // nodal dofs for ele
        eleDofManager_uncondensed_ =
          rcp(new XFEM::ElementDofManager(*this, eleDofManager_->getNodalDofSet(), enrfieldset, element_ansatz_filled));

        const int nd = eleDofManager_uncondensed_->NumNodeDof();
        const int na = eleDofManager_uncondensed_->NumElemDof();

        if (nd != nd_old or na != na_old)
        {

          DLM_info_ = Teuchos::rcp(new DLMInfo(nd,na));
        }
        else
        {
          cout << "saved steps" << endl;
        }
      }
      else
      {
        eleDofManager_uncondensed_ = Teuchos::null;
        DLM_info_ = Teuchos::null;
      }
      break;
    }
    case calc_fluid_systemmat_and_residual:
    {
      // do no calculation, if not needed
      if (lm.empty())
        break;

      // extract local values from the global vectors
      DRT::ELEMENTS::XFluid3::MyState mystate(discretization,lm,true);

      const Teuchos::RCP<Epetra_Vector> iforcecol = params.get<Teuchos::RCP<Epetra_Vector> >("interface force");

      double L2 = params.get<double>("L2");

      // time integration factors
      const FLUID_TIMEINTTYPE timealgo = params.get<FLUID_TIMEINTTYPE>("timealgo");
      const double            dt       = params.get<double>("dt");
      const double            theta    = params.get<double>("theta");

      const bool newton = params.get<bool>("include reactive terms for linearisation");
      const bool pstab  = true;
      const bool supg   = true;
      const bool cstab  = true;

      const bool ifaceForceContribution = discretization.ElementRowMap()->MyGID(this->Id());
      const bool monolithic_FSI = params.get<bool>("monolithic_FSI");

      RCP<Epetra_SerialDenseMatrix> Cuu;
      RCP<Epetra_SerialDenseMatrix> Mud;
      RCP<Epetra_SerialDenseMatrix> Mdu;
      RCP<Epetra_SerialDenseMatrix> Cdd;
      RCP<Epetra_SerialDenseVector> rhsd;

      RCP<Epetra_SerialDenseMatrix> Gds_uncond;
      RCP<Epetra_SerialDenseVector> rhsd_uncond;

      if (not params.get<bool>("DLM_condensation") or not ih_->ElementIntersected(Id())) // integrate and assemble all unknowns
      {
        const XFEM::AssemblyType assembly_type = CheckForStandardEnrichmentsOnly(
                *eleDofManager_, NumNode(), NodeIds());

        // calculate element coefficient matrix and rhs
        XFLUID::callSysmat4(assembly_type,
                this, ih_, *eleDofManager_, mystate, iforcecol, elemat1, elevec1, *Gds_uncond, *rhsd_uncond,
                mat, timealgo, dt, theta, newton, pstab, supg, cstab, mystate.instationary, ifaceForceContribution, monolithic_FSI, L2);
      }
      else // create bigger element matrix and vector, assemble, condense and copy to small matrix provided by discretization
      {
        // sanity checks
        SanityChecks(eleDofManager_, eleDofManager_uncondensed_);

        // stress update
        UpdateOldDLMAndDLMRHS(discretization, lm, mystate);

        // create uncondensed element matrix and vector
        const int numdof_uncond = eleDofManager_uncondensed_->NumDofElemAndNode();
        Epetra_SerialDenseMatrix elemat1_uncond(numdof_uncond,numdof_uncond);
        Epetra_SerialDenseVector elevec1_uncond(numdof_uncond);

        const XFEM::AssemblyType assembly_type = CheckForStandardEnrichmentsOnly(
                *eleDofManager_uncondensed_, NumNode(), NodeIds());

        if (ih_->ElementIntersected(Id()) and monolithic_FSI)
        {
          Cuu  = params.get<RCP<Epetra_SerialDenseMatrix> >("Cuu");
          Mud  = params.get<RCP<Epetra_SerialDenseMatrix> >("Mud");
          Mdu  = params.get<RCP<Epetra_SerialDenseMatrix> >("Mdu");
          Cdd  = params.get<RCP<Epetra_SerialDenseMatrix> >("Cdd");
          rhsd = params.get<RCP<Epetra_SerialDenseVector> >("rhsd");

          const std::set<int> begids = ih_->GetIntersectingBoundaryElementsGID(this->Id());
          const int numnode_b = 4;
          const size_t numpatchdof = 3*numnode_b*begids.size();
          Gds_uncond  = rcp(new Epetra_SerialDenseMatrix(numpatchdof, eleDofManager_uncondensed_->NumDofElemAndNode()));
          rhsd_uncond = rcp(new Epetra_SerialDenseVector(numpatchdof));
        }

        // calculate element coefficient matrix and rhs
        XFLUID::callSysmat4(assembly_type,
                this, ih_, *eleDofManager_uncondensed_, mystate, iforcecol, elemat1_uncond, elevec1_uncond, *Gds_uncond, *rhsd_uncond,
                mat, timealgo, dt, theta, newton, pstab, supg, cstab, mystate.instationary, ifaceForceContribution, monolithic_FSI, L2);

        // condensation
        CondenseDLMAndStoreOldIterationStep(
            elemat1_uncond, elevec1_uncond,
            *Gds_uncond, *rhsd_uncond,
            elemat1, elevec1,
            *Cuu, *Mud, *Mdu, *Cdd, *rhsd,
            monolithic_FSI
            );
      }
      params.set<double>("L2",L2);
      break;
    }
    case calc_fluid_stationary_systemmat_and_residual:
    {
      // do no calculation, if not needed
      if (lm.empty())
        break;

      // extract local values from the global vector
      DRT::ELEMENTS::XFluid3::MyState mystate(discretization,lm,false);

      const Teuchos::RCP<Epetra_Vector> iforcecol = params.get<Teuchos::RCP<Epetra_Vector> >("interface force");

      double L2 = params.get<double>("L2");

      // time integration factors
      const FLUID_TIMEINTTYPE timealgo = params.get<FLUID_TIMEINTTYPE>("timealgo");
      const double            dt       = 1.0;
      const double            theta    = 1.0;

      const bool newton = params.get<bool>("include reactive terms for linearisation");
      const bool pstab  = true;
      const bool supg   = true;
      const bool cstab  = true;

      const bool ifaceForceContribution = discretization.ElementRowMap()->MyGID(this->Id());
      const bool monolithic_FSI = params.get<bool>("monolithic_FSI");

      RCP<Epetra_SerialDenseMatrix> Cuu;
      RCP<Epetra_SerialDenseMatrix> Mud;
      RCP<Epetra_SerialDenseMatrix> Mdu;
      RCP<Epetra_SerialDenseMatrix> Cdd;
      RCP<Epetra_SerialDenseVector> rhsd;

      RCP<Epetra_SerialDenseMatrix> Gds_uncond;
      RCP<Epetra_SerialDenseVector> rhsd_uncond;

      if (not params.get<bool>("DLM_condensation") or not ih_->ElementIntersected(Id())) // integrate and assemble all unknowns
      {
        const XFEM::AssemblyType assembly_type = CheckForStandardEnrichmentsOnly(
                *eleDofManager_, NumNode(), NodeIds());

        // calculate element coefficient matrix and rhs
        XFLUID::callSysmat4(assembly_type,
                this, ih_, *eleDofManager_, mystate, iforcecol, elemat1, elevec1, *Gds_uncond, *rhsd_uncond,
                mat, timealgo, dt, theta, newton, pstab, supg, cstab, mystate.instationary, ifaceForceContribution, monolithic_FSI, L2);
      }
      else // create bigger element matrix and vector, assemble, condense and copy to small matrix provided by discretization
      {
        // sanity checks
        SanityChecks(eleDofManager_, eleDofManager_uncondensed_);

        // stress update
        UpdateOldDLMAndDLMRHS(discretization, lm, mystate);

        // create uncondensed element matrix and vector
        const int numdof_uncond = eleDofManager_uncondensed_->NumDofElemAndNode();
        Epetra_SerialDenseMatrix elemat1_uncond(numdof_uncond,numdof_uncond);
        Epetra_SerialDenseVector elevec1_uncond(numdof_uncond);

        const XFEM::AssemblyType assembly_type = CheckForStandardEnrichmentsOnly(
                *eleDofManager_uncondensed_, NumNode(), NodeIds());

        if (ih_->ElementIntersected(Id()) and monolithic_FSI)
        {
          Cuu  = params.get<RCP<Epetra_SerialDenseMatrix> >("Cuu");
          Mud  = params.get<RCP<Epetra_SerialDenseMatrix> >("Mud");
          Mdu  = params.get<RCP<Epetra_SerialDenseMatrix> >("Mdu");
          Cdd  = params.get<RCP<Epetra_SerialDenseMatrix> >("Cdd");
          rhsd = params.get<RCP<Epetra_SerialDenseVector> >("rhsd");

          const std::set<int> begids = ih_->GetIntersectingBoundaryElementsGID(this->Id());
          const int numnode_b = 4;
          const size_t numpatchdof = 3*numnode_b*begids.size();
          Gds_uncond  = rcp(new Epetra_SerialDenseMatrix(numpatchdof, eleDofManager_uncondensed_->NumDofElemAndNode()));
          rhsd_uncond = rcp(new Epetra_SerialDenseVector(numpatchdof));
        }

        // calculate element coefficient matrix and rhs
        XFLUID::callSysmat4(assembly_type,
                this, ih_, *eleDofManager_uncondensed_, mystate, iforcecol, elemat1_uncond, elevec1_uncond, *Gds_uncond, *rhsd_uncond,
                mat, timealgo, dt, theta, newton, pstab, supg, cstab, mystate.instationary, ifaceForceContribution, monolithic_FSI, L2);

        // condensation
        CondenseDLMAndStoreOldIterationStep(
            elemat1_uncond, elevec1_uncond,
            *Gds_uncond, *rhsd_uncond,
            elemat1, elevec1,
            *Cuu, *Mud, *Mdu, *Cdd, *rhsd,
            monolithic_FSI
            );
      }

      params.set<double>("L2",L2);

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
      break;
    }
    case calc_fluid_projection_systemmat_and_residual:
    {
      // do no calculation, if not needed
      if (lm.empty())
        break;

      // extract local values from the global vector
      DRT::ELEMENTS::XFluid3::MyState mystate(discretization,lm,true);

      const bool pstab  = true;


      const bool ifaceForceContribution = discretization.ElementRowMap()->MyGID(this->Id());

      const XFEM::AssemblyType assembly_type = CheckForStandardEnrichmentsOnly(
              *eleDofManager_, NumNode(), NodeIds());

      // calculate element coefficient matrix and rhs
      XFLUID::callSysmatProjection(assembly_type,
              this, ih_, *eleDofManager_, mystate, elemat1, elemat2, elevec1, elevec2,
              pstab, ifaceForceContribution);


      break;
    }
    case get_density:
    {
      // This is a very poor way to transport the density to the
      // outside world. Is there a better one?
      params.set("density", actmat->Density());
      break;
    }
    case integrate_shape:
    {
      // do no calculation, if not needed
      if (lm.empty())
        break;

      const XFEM::AssemblyType assembly_type = CheckForStandardEnrichmentsOnly(
          *eleDofManager_, NumNode(), NodeIds());

      // calculate element coefficient matrix and rhs
      integrateShapefunction(assembly_type, this, ih_, *eleDofManager_, elemat1, elevec1);

      break;
    }
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
        f3_int_beltrami_err(myvelnp,myprenp,mat,params);
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
 |  integration of the volume neumann (body forces) loads takes place   |
 |  in the element. We need it there for the stabilisation terms!       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::XFluid3::EvaluateNeumann(ParameterList& params,
                                            DRT::Discretization&      discretization,
                                            DRT::Condition&           condition,
                                            std::vector<int>&         lm,
                                            Epetra_SerialDenseVector& elevec1,
                                            Epetra_SerialDenseMatrix* elemat1)
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

  const double a      = PI/4.0;
  const double d      = PI/2.0;

  // get viscosity
  double  visc = 0.0;
  if(material->MaterialType() == INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = dynamic_cast<const MAT::NewtonianFluid*>(material.get());
    visc = actmat->Viscosity();
  }
  else
    dserror("Cannot handle material of type %d", material->MaterialType());

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

/*---------------------------------------------------------------------*
 *---------------------------------------------------------------------*/
void DRT::ELEMENTS::XFluid3::UpdateOldDLMAndDLMRHS(
    const DRT::Discretization&      discretization,
    const std::vector<int>&         lm,
    MyState&                        mystate
    ) const
{
  const int nu = eleDofManager_uncondensed_->NumNodeDof();
  const int ns = eleDofManager_uncondensed_->NumElemDof();

  if (ns > 0)
  {
    // add Kda . inc_velnp to feas
    // new alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas

    vector<double> inc_velnp(lm.size());
    DRT::UTILS::ExtractMyValues(*discretization.GetState("nodal increment"),inc_velnp,lm);

    static const Epetra_BLAS blas;

    // update old iteration residual of the stresses
    // DLM_info_->oldfa_(i) += DLM_info_->oldKad_(i,j)*inc_velnp[j];
    blas.GEMV('N', ns, nu,-1.0, DLM_info_->oldKad_.A(), DLM_info_->oldKad_.LDA(), &inc_velnp[0], 1.0, DLM_info_->oldfa_.A());

    // compute element stresses
    // DLM_info_->stressdofs_(i) -= DLM_info_->oldKaainv_(i,j)*DLM_info_->oldfa_(j);
    blas.GEMV('N', ns, ns,1.0, DLM_info_->oldKaainv_.A(), DLM_info_->oldKaainv_.LDA(), DLM_info_->oldfa_.A(), 1.0, DLM_info_->stressdofs_.A());

    // increase size of element vector (old values stay and zeros are added)
    const int numdof_uncond = eleDofManager_uncondensed_->NumDofElemAndNode();
    mystate.velnp.resize(numdof_uncond,0.0);
    mystate.veln .resize(numdof_uncond,0.0);
    mystate.velnm.resize(numdof_uncond,0.0);
    mystate.accn .resize(numdof_uncond,0.0);
    for (int i=0;i<ns;i++)
    {
      mystate.velnp[nu+i] = DLM_info_->stressdofs_(i);
    }
  }
}

/*---------------------------------------------------------------------*
 *---------------------------------------------------------------------*/
void DRT::ELEMENTS::XFluid3::CondenseDLMAndStoreOldIterationStep(
    const Epetra_SerialDenseMatrix& elemat1_uncond,
    const Epetra_SerialDenseVector& elevec1_uncond,
    const Epetra_SerialDenseMatrix& Gds_uncond,
    const Epetra_SerialDenseVector& rhsd_uncond,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix& Cuu,
    Epetra_SerialDenseMatrix& Mud,
    Epetra_SerialDenseMatrix& Mdu,
    Epetra_SerialDenseMatrix& Cdd,
    Epetra_SerialDenseVector& rhsd,
    const bool monolithic_FSI
) const
{

  const size_t nu = eleDofManager_uncondensed_->NumNodeDof();
  const size_t ns = eleDofManager_uncondensed_->NumElemDof();

  // copy nodal dof entries
  for (size_t i = 0; i < nu; ++i)
  {
    elevec1(i) = elevec1_uncond(i);
    for (size_t j = 0; j < nu; ++j)
    {
      elemat1(i,j) = elemat1_uncond(i,j);
    }
  }

  if (ns > 0)
  {
    // note: the full (u,p,sigma) matrix is asymmetric,
    // hence we need both rectangular matrices Kda and Kad
    LINALG::SerialDenseMatrix Gus(nu,ns);
    LINALG::SerialDenseMatrix Kssinv(ns,ns);
    LINALG::SerialDenseMatrix KGsu(ns,nu);
    LINALG::SerialDenseVector fs(ns);

//    cout << elemat1_uncond << endl;

    // copy data of uncondensed matrix into submatrices
    for (size_t i=0;i<nu;i++)
      for (size_t j=0;j<ns;j++)
        Gus(i,j) = elemat1_uncond(   i,nu+j);

    for (size_t i=0;i<ns;i++)
      for (size_t j=0;j<ns;j++)
        Kssinv(i,j) = elemat1_uncond(nu+i,nu+j);

    for (size_t i=0;i<ns;i++)
      for (size_t j=0;j<nu;j++)
        KGsu(i,j) = elemat1_uncond(nu+i,   j);

    for (size_t i=0;i<ns;i++)
      fs(i) = elevec1_uncond(nu+i);

    // DLM-stiffness matrix is: Kdd - Kda . Kaa^-1 . Kad
    // DLM-internal force is: fint - Kda . Kaa^-1 . feas

    // we need the inverse of Kaa
    Epetra_SerialDenseSolver solve_for_inverseKaa;
    solve_for_inverseKaa.SetMatrix(Kssinv);
    solve_for_inverseKaa.Invert();

    static const Epetra_BLAS blas;
    {
      LINALG::SerialDenseMatrix GusKssinv(nu,ns); // temporary Gus.Kss^{-1}

      // KusKssinv(i,j) = Kus(i,k)*Kssinv(k,j);
      blas.GEMM('N','N',nu,ns,ns,1.0,Gus.A(),Gus.LDA(),Kssinv.A(),Kssinv.LDA(),0.0,GusKssinv.A(),GusKssinv.LDA());

      // elemat1(i,j) += - KusKssinv(i,k)*Ksu(k,j);
      blas.GEMM('N','N',nu,nu,ns,-1.0,GusKssinv.A(),GusKssinv.LDA(),KGsu.A(),KGsu.LDA(),1.0,elemat1.A(),elemat1.LDA());

      // elevec1(i) += - KusKssinv(i,j)*fs(j);
      blas.GEMV('N', nu, ns,-1.0, GusKssinv.A(), GusKssinv.LDA(), fs.A(), 1.0, elevec1.A());

      if (monolithic_FSI)
      {
        const std::set<int> begids = ih_->GetIntersectingBoundaryElementsGID(this->Id());
        const int numnode_b = 4;
        const size_t nd = 3*numnode_b*begids.size();

        Epetra_SerialDenseMatrix Gds(nd, ns);
        Epetra_SerialDenseMatrix Gsd(ns, nd);

        for (size_t i=0;i<nd;i++)
        {
          rhsd(i) = rhsd_uncond(i);
          for (size_t j=0;j<ns;j++)
          {
            Gds(i,j) = Gds_uncond(i, nu+j);
            Gsd(j,i) = Gds_uncond(i, nu+j); //only for transient
          }
        }

        LINALG::SerialDenseMatrix GdsKssinv(nd,ns); // temporary Kds.Kss^{-1}

        // KdsKssinv(i,j) = Kds(i,k)*Kssinv(k,j);
        blas.GEMM('N','N',nd,ns,ns,1.0,Gds.A(),Gds.LDA(),Kssinv.A(),Kssinv.LDA(),0.0,GdsKssinv.A(),GdsKssinv.LDA());

//        for (size_t i=0;i<nu;i++)
//          for (size_t j=0;j<nu;j++)
//            for (size_t k=0;k<ns;k++)
//              Cuu(i,j) += KusKssinv(i,k)*Ksu(k,j);
        for (size_t i=0;i<nu;i++)
          for (size_t j=0;j<nd;j++)
            for (size_t k=0;k<ns;k++)
              Mud(i,j) += -GusKssinv(i,k)*Gsd(k,j);
        for (size_t i=0;i<nd;i++)
          for (size_t j=0;j<nu;j++)
            for (size_t k=0;k<ns;k++)
              Mdu(i,j) += -GdsKssinv(i,k)*KGsu(k,j);
        for (size_t i=0;i<nd;i++)
          for (size_t j=0;j<nd;j++)
            for (size_t k=0;k<ns;k++)
              Cdd(i,j) += -GdsKssinv(i,k)*Gsd(k,j);
        for (size_t i=0;i<nd;i++)
          for (size_t j=0;j<ns;j++)
            rhsd(i) += - GdsKssinv(i,j)*fs(j);
      }


    }

    // store current DLM data in iteration history
    //DLM_info_->oldKaainv_.Update(1.0,Kaa,0.0);
    blas.COPY(DLM_info_->oldKaainv_.M()*DLM_info_->oldKaainv_.N(), Kssinv.A(), DLM_info_->oldKaainv_.A());
    //DLM_info_->oldKad_.Update(1.0,Kad,0.0);
    blas.COPY(DLM_info_->oldKad_.M()*DLM_info_->oldKad_.N(), KGsu.A(), DLM_info_->oldKad_.A());
    //DLM_info_->oldfa_.Update(1.0,fa,0.0);
    blas.COPY(DLM_info_->oldfa_.M()*DLM_info_->oldfa_.N(), fs.A(), DLM_info_->oldfa_.A());
  }
}


//! size factor to allow fixed size arrays
///
/// to allow fixed size arrays for a unknown number of unknowns, we make them bigger than necessary
/// this factor is multiplied times numnode(distype) to get the size of many arrays
template<XFEM::AssemblyType ASSTYPE>
struct SizeFac {};
/// specialization of SizeFac for XFEM::standard_assembly
template<> struct SizeFac<XFEM::standard_assembly> {static const std::size_t fac = 1;};
/// specialization of SizeFac for XFEM::xfem_assembly
template<> struct SizeFac<XFEM::xfem_assembly>     {static const std::size_t fac = 3;};

/*!
  Calculate matrix and rhs for stationary problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE>
void integrateShapefunctionT(
        const DRT::Element*                             ele,     ///< the element those matrix is calculated
        const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,      ///< connection to the interface handler
        const XFEM::ElementDofManager&                  dofman,  ///< dofmanager of the current element
        Epetra_SerialDenseMatrix&                       estif,   ///< element matrix to calculate
        Epetra_SerialDenseVector&                       eforce   ///< element rhs to calculate
        )
{
  // initialize arrays
  estif.Scale(0.0);
  eforce.Scale(0.0);

  const int NUMDOF = 4;

  LocalAssembler<DISTYPE, ASSTYPE, NUMDOF> assembler(dofman, estif, eforce);

  // number of nodes for element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // space dimension for 3d fluid element
  const size_t nsd = 3;

  // get node coordinates of the current element
  static LINALG::Matrix<nsd,numnode> xyze;
  GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

  // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);

  // information about domain integration cells
  const GEO::DomainIntCells&  domainIntCells(ih->GetDomainIntCells(ele));

  // loop over integration cells
  for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
  {
    const LINALG::Matrix<nsd,1> cellcenter_xyz(cell->GetPhysicalCenterPosition());

    int labelnp = 0;

    if (ASSTYPE == XFEM::xfem_assembly)
    {
      // integrate only in fluid integration cells (works well only with void enrichments!!!)
      labelnp = ih->PositionWithinConditionNP(cellcenter_xyz);
      const std::set<int> xlabelset(dofman.getUniqueEnrichmentLabels());
      bool compute = false;
      if (labelnp == 0) // fluid
      {
        compute = true;
      }
      if (not compute)
      {
        continue; // next integration cell
      }
    }

    const XFEM::ElementEnrichmentValues enrvals(
        *ele,
        ih,
        dofman,
        cellcenter_xyz,
        XFEM::Enrichment::approachUnknown);

    const DRT::UTILS::GaussRule3D gaussrule = XFLUID::getXFEMGaussrule<DISTYPE>(ele, xyze, ih->ElementIntersected(ele->Id()),cell->Shape());

    // gaussian points
    const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

    // integration loop
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
      // coordinates of the current integration point in cell coordinates \eta
      LINALG::Matrix<nsd,1> pos_eta_domain;
      pos_eta_domain(0) = intpoints.qxg[iquad][0];
      pos_eta_domain(1) = intpoints.qxg[iquad][1];
      pos_eta_domain(2) = intpoints.qxg[iquad][2];

      // coordinates of the current integration point in element coordinates \xi
      LINALG::Matrix<nsd,1> posXiDomain;
      GEO::mapEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain, posXiDomain);

      const double detcell = GEO::detEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain);

      // shape functions and their first derivatives
      static LINALG::Matrix<numnode,1> funct;
      static LINALG::Matrix<nsd,numnode> deriv;
      DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

      // get transposed of the jacobian matrix d x / d \xi
      // xjm(i,j) = deriv(i,k)*xyze(j,k)
      static LINALG::Matrix<nsd,nsd> xjm;
      xjm.MultiplyNT(deriv,xyze);

      const double det = xjm.Determinant();
      const double fac = intpoints.qwgt[iquad]*det*detcell;

      if (det < 0.0)
      {
        dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
      }

      // inverse of jacobian
      static LINALG::Matrix<nsd,nsd> xji;
      xji.Invert(xjm);

      // compute global derivates: derxy(i,j) = xji(i,k) * deriv(k,j)
      static LINALG::Matrix<3,numnode> derxy;
      derxy.Multiply(xji,deriv);

      const size_t shpVecSize       = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

      static XFLUID::ApproxFunc<shpVecSize> shp;

      if (ASSTYPE == XFEM::xfem_assembly)
      {
        // temporary arrays
        static LINALG::Matrix<shpVecSize,1> enr_funct;
        static LINALG::Matrix<3,shpVecSize> enr_derxy;

        // shape function for nodal dofs
        enrvals.ComputeEnrichedNodalShapefunction(
            XFEM::PHYSICS::Velx,
            funct,
            derxy,
            enr_funct,
            enr_derxy);

        for (size_t iparam = 0; iparam != numparampres; ++iparam)
        {
          shp.d0(iparam) = enr_funct(iparam);
          shp.dx(iparam) = enr_derxy(0,iparam);
          shp.dy(iparam) = enr_derxy(1,iparam);
          shp.dz(iparam) = enr_derxy(2,iparam);
        }
      }
      else // standard assembly
      {
        // -> numparamvelx == numnode
        for (size_t iparam = 0; iparam < numnode; ++iparam)
        {
          shp.d0(iparam) = funct(iparam);
          shp.dx(iparam) = derxy(0,iparam);
          shp.dy(iparam) = derxy(1,iparam);
          shp.dz(iparam) = derxy(2,iparam);
        }
      }

      //////////////////////////////////////
      // now build single stiffness terms //
      //////////////////////////////////////

      assembler.template Vector<XFEM::PHYSICS::Pres>(shp.d0, fac);
      assembler.template Vector<XFEM::PHYSICS::Pres>(shp.d0, fac);
      assembler.template Vector<XFEM::PHYSICS::Pres>(shp.d0, fac);

    } // end loop over gauss points
  } // end loop over integration cells

  return;
}


void DRT::ELEMENTS::XFluid3::integrateShapefunction(
        const XFEM::AssemblyType&                       assembly_type,
        const DRT::Element*                             ele,     ///< the element those matrix is calculated
        const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,      ///< connection to the interface handler
        const XFEM::ElementDofManager&                  dofman,  ///< dofmanager of the current element
        Epetra_SerialDenseMatrix&                       estif,   ///< element matrix to calculate
        Epetra_SerialDenseVector&                       eforce   ///< element rhs to calculate
        )
{
// avoid repeating the same thing again and again
#define M_UNROLL_INTEGRATESHAPEFUNCTION_CASE(distype,asstype) {case distype : integrateShapefunctionT<distype, asstype>(ele, ih, dofman, estif, eforce); break;}

  if (assembly_type == XFEM::standard_assembly)
  {
    switch (ele->Shape())
    {
      M_UNROLL_INTEGRATESHAPEFUNCTION_CASE(DRT::Element::hex8 ,XFEM::standard_assembly)
      M_UNROLL_INTEGRATESHAPEFUNCTION_CASE(DRT::Element::hex20,XFEM::standard_assembly)
      M_UNROLL_INTEGRATESHAPEFUNCTION_CASE(DRT::Element::hex27,XFEM::standard_assembly)
      M_UNROLL_INTEGRATESHAPEFUNCTION_CASE(DRT::Element::tet4 ,XFEM::standard_assembly)
      M_UNROLL_INTEGRATESHAPEFUNCTION_CASE(DRT::Element::tet10,XFEM::standard_assembly)
      default:
        dserror("standard_assembly integrateShapefunctionT not templated yet");
    };
  }
  else
  {
    switch (ele->Shape())
    {
      M_UNROLL_INTEGRATESHAPEFUNCTION_CASE(DRT::Element::hex8 ,XFEM::xfem_assembly)
      M_UNROLL_INTEGRATESHAPEFUNCTION_CASE(DRT::Element::hex20,XFEM::xfem_assembly)
      M_UNROLL_INTEGRATESHAPEFUNCTION_CASE(DRT::Element::hex27,XFEM::xfem_assembly)
      M_UNROLL_INTEGRATESHAPEFUNCTION_CASE(DRT::Element::tet4 ,XFEM::xfem_assembly)
      M_UNROLL_INTEGRATESHAPEFUNCTION_CASE(DRT::Element::tet10,XFEM::xfem_assembly)
      default:
        dserror("xfem_assembly integrateShapefunctionT not templated yet");
    };
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
