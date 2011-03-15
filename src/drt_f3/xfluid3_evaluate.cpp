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
#include "../drt_xfem/xfem_element_utils.H"
#include "../drt_geometry/integrationcell_coordtrafo.H"
#include "../drt_xfem/physics.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/xdofmapcreation_fsi.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_xfem/interfacexfsi.H"
#include "../drt_fem_general/debug_nan.H"


/*---------------------------------------------------------------------*
|  converts a string into an Action for this element                   |
*----------------------------------------------------------------------*/
DRT::ELEMENTS::XFluid3::ActionType DRT::ELEMENTS::XFluid3::convertStringToActionType(
              const string& action) const
{
  DRT::ELEMENTS::XFluid3::ActionType act = XFluid3::none;
  if (action == "calc_fluid_systemmat_and_residual")
    act = XFluid3::calc_fluid_systemmat_and_residual;
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
  else if (action == "fluidfluidCoupling")
    act = XFluid3::fluidfluidCoupling;
  else if (action == "fluidxfluidCoupling")
      act = XFluid3::fluidxfluidCoupling;
  else
  {
    cout << "Unknown action: " << action << endl;
    dserror("Unknown type of action for XFluid3");
  }
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
      ih_ = NULL;
      DLM_info_ = Teuchos::null;
      break;
    }
    case set_output_mode:
    {
      output_mode_ = params.get<bool>("output_mode");
      // reset dof managers if present
      eleDofManager_ = Teuchos::null;
      eleDofManager_uncondensed_ = Teuchos::null;
      ih_ = NULL;
      DLM_info_ = Teuchos::null;
      break;
    }
    case store_xfem_info:
    {
      // after this part the element can answer, how many DOFs it has
      output_mode_ = false;

      // store pointer to interface handle
      ih_ = &*params.get< Teuchos::RCP< XFEM::InterfaceHandleXFSI > >("interfacehandle");

      // get access to global dofman
      const Teuchos::RCP<const XFEM::DofManager> globaldofman = params.get< Teuchos::RCP< XFEM::DofManager > >("dofmanager");

      Teuchos::RCP<XFEM::ElementAnsatz> elementAnsatz;
      switch (DRT::INPUT::get<INPAR::XFEM::BoundaryIntegralType>(params, "EMBEDDED_BOUNDARY"))
      {
      case INPAR::XFEM::BoundaryTypeSigma:
        elementAnsatz = rcp<XFLUID::FluidElementAnsatz>(new XFLUID::FluidElementAnsatz());
        break;
      case INPAR::XFEM::BoundaryTypeTauPressure:
        elementAnsatz = rcp<XFLUID::FluidElementAnsatzWithExtraElementPressure>(new XFLUID::FluidElementAnsatzWithExtraElementPressure());
        break;
      default:
        dserror("unknown boundary type");
      }
      const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz_empty;
      const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz_filled(elementAnsatz->getElementAnsatz(this->Shape()));

      // always build the eledofman that fits to the global dofs
      // problem: tight connectivity to xdofmapcreation
      if (params.get<bool>("DLM_condensation"))
      {
//        globaldofman->  dserror check
        // assume no stress unknowns for the element
        eleDofManager_ = rcp(new XFEM::ElementDofManager(*this, element_ansatz_empty, *globaldofman));
      }
      else
      {
        // assume stress unknowns for the element
        eleDofManager_ = rcp(new XFEM::ElementDofManager(*this, element_ansatz_filled, *globaldofman));
      }

      // create an eledofman that has stress unknowns only for intersected elements
      // Note: condensation for uncut elements is not possible/needed
      if (ih_->ElementIntersected(Id()))
      {
        std::set<XFEM::FieldEnr> enrfieldset;

        bool skipped_elem_enr;
        XFEM::processVoidEnrichmentForElement(
            this, element_ansatz_filled, *ih_, eleDofManager_->getUniqueEnrichmentLabels(),
            params.get<double>("boundaryRatioLimit"),
            params.get<double>("volumeRatioLimit"),
            enrfieldset,
            skipped_elem_enr);

        // nodal dofs for ele
        eleDofManager_uncondensed_ =
          rcp(new XFEM::ElementDofManager(
              *this, eleDofManager_->getNodalDofSet(), enrfieldset, element_ansatz_filled
              ));

        const RCP<vector<int> > ifacepatchlm = rcp(new vector<int>());
        const RCP<vector<int> > ifacepatchlmowner = rcp(new vector<int>());
        ih_->GetInterfacepatchLocationVectors(*this, ifacepatchlm, ifacepatchlmowner);

        if (DLM_info_ == Teuchos::null)
        {
          DLM_info_ = Teuchos::rcp(
              new DLMInfo(
                  eleDofManager_uncondensed_->NumNodeDof(),
                  eleDofManager_uncondensed_->NumElemDof(),
                  ifacepatchlm->size()
                  )
              );
        }
        else if (DLM_info_->oldGsui_.N() != (int)ifacepatchlm->size()) // different number of intersecting elements -> coupling matrices change size!
        {
          // todo: rescue stress instead of deleting
          cout << "stress rescue" << endl;
//          cout << "DLM_info_->oldGsui_.N() != ifacepatchlm->size() -> reset" << DLM_info_->oldGsui_.N() << "  " << ifacepatchlm->size() << endl;
          Teuchos::RCP<DLMInfo> DLM_info = Teuchos::rcp(
              new DLMInfo(
                  eleDofManager_uncondensed_->NumNodeDof(),
                  eleDofManager_uncondensed_->NumElemDof(),
                  ifacepatchlm->size()
                  ,*DLM_info_
                  )
              );
          DLM_info_ = DLM_info;
        }
        else
        {
//          cout << DLM_info_->oldGsui_.N() << "  " << DLM_info_->oldGsui_.M() << "  " << ifacepatchlm->size() << endl;
        }
      }
      else
      {
//        eleDofManager_uncondensed_ = Teuchos::null;
//        DLM_info_ = Teuchos::null;
      }
      break;
    }
    case stress_update:
    {
      // do no calculation, if not needed
      if (lm.empty())
        break;

      const Teuchos::RCP<Epetra_Vector> iforcecol = params.get<Teuchos::RCP<Epetra_Vector> >("interface force");

      // time integration factors
      const INPAR::FLUID::TimeIntegrationScheme timealgo = DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(params, "timealgo");

      // extract local values from the global vectors
      const bool instationary = (timealgo != INPAR::FLUID::timeint_stationary);

      DRT::ELEMENTS::XFluid3::MyState mystate(discretization,lm,instationary);

      const bool monolithic_FSI = params.get<bool>("monolithic_FSI");

      const RCP<const vector<int> > ifacepatchlm = params.get<RCP<vector<int> > >("ifacepatchlm");

      if (not params.get<bool>("DLM_condensation") or not ih_->ElementIntersected(Id())) // integrate and assemble all unknowns
      {

      }
      else // create bigger element matrix and vector, assemble, condense and copy to small matrix provided by discretization
      {
        // sanity checks
        SanityChecks(eleDofManager_, eleDofManager_uncondensed_);

        const RCP<const Epetra_Vector>  iterincxdomain = discretization.GetState("velpres nodal iterinc");
        const RCP<const Epetra_Vector>  iterinciface   = ih_->cutterdis()->GetState("veliface nodal iterinc");

        // stress update
        UpdateOldDLMAndDLMRHS(iterincxdomain, iterinciface, lm, *ifacepatchlm, mystate, monolithic_FSI);

      }
      break;
    }
    case calc_fluid_systemmat_and_residual:
    {
      // do no calculation, if not needed
      if (lm.empty())
        break;

      const Teuchos::RCP<Epetra_Vector> iforcecol = params.get<Teuchos::RCP<Epetra_Vector> >("interface force");

      double L2 = params.get<double>("L2");

      // time integration factors
      const INPAR::FLUID::TimeIntegrationScheme timealgo = DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(params, "timealgo");
      const double            dt       = params.get<double>("dt");
      const double            theta    = params.get<double>("theta");

      // extract local values from the global vectors
      const bool instationary = (timealgo != INPAR::FLUID::timeint_stationary);

      DRT::ELEMENTS::XFluid3::MyState mystate(discretization,lm,instationary);
      DRT::DEBUGGING::NaNChecker(mystate.velnp);

      const bool newton = params.get<bool>("include reactive terms for linearisation");
      const bool pstab  = true;
      const bool supg   = true;
      const bool cstab  = true;

      const bool ifaceForceContribution = discretization.ElementRowMap()->MyGID(this->Id());
      const bool monolithic_FSI = params.get<bool>("monolithic_FSI");

      const RCP<const vector<int> > ifacepatchlm = params.get<RCP<vector<int> > >("ifacepatchlm");

      if (not params.get<bool>("DLM_condensation") or not ih_->ElementIntersected(Id())) // integrate and assemble all unknowns
      {
        if (ih_->ElementIntersected(Id()))
        {
          const size_t nus = eleDofManager_uncondensed_->NumDofElemAndNode();
          const size_t nui = ifacepatchlm->size();
          fluidfluidmatrices_.Guis_uncond  = rcp(new Epetra_SerialDenseMatrix(nui, nus));
          fluidfluidmatrices_.Gsui_uncond  = rcp(new Epetra_SerialDenseMatrix(nus, nui));
          fluidfluidmatrices_.rhsui_uncond = rcp(new Epetra_SerialDenseVector(nui));
        }

        const XFEM::AssemblyType assembly_type = XFEM::ComputeAssemblyType(
                *eleDofManager_, NumNode(), NodeIds());

        // calculate element coefficient matrix and rhs
        XFLUID::callSysmat(DRT::INPUT::get<INPAR::XFEM::BoundaryIntegralType>(params, "EMBEDDED_BOUNDARY"), params, assembly_type,
                this, ih_, *eleDofManager_, mystate, iforcecol, elemat1, elevec1,
                mat, timealgo, dt, theta, newton, pstab, supg, cstab, ifaceForceContribution, monolithic_FSI, L2,fluidfluidmatrices_);

        if (ih_->ElementIntersected(Id()))
        {
          const size_t nui = ifacepatchlm->size();
          RCP<Epetra_SerialDenseMatrix> Cdd = rcp(new Epetra_SerialDenseMatrix(nui, nui));
          params.set("Cdu",fluidfluidmatrices_.Guis_uncond);
          params.set("Cud",fluidfluidmatrices_.Gsui_uncond);
          params.set("Cdd",Cdd);
          params.set("rhsd",fluidfluidmatrices_.rhsui_uncond);
        }
      }
      else // create bigger element matrix and vector, assemble, condense and copy to small matrix provided by discretization
      {
        // sanity checks
        SanityChecks(eleDofManager_, eleDofManager_uncondensed_);

        const RCP<const Epetra_Vector>  iterincxdomain = discretization.GetState("velpres nodal iterinc");
        const RCP<const Epetra_Vector>  iterinciface   = ih_->cutterdis()->GetState("veliface nodal iterinc");

        // stress update
        UpdateOldDLMAndDLMRHS(iterincxdomain, iterinciface, lm, *ifacepatchlm, mystate, monolithic_FSI);

        // create uncondensed element matrix and vector
        const int numdof_uncond = eleDofManager_uncondensed_->NumDofElemAndNode();
        Epetra_SerialDenseMatrix elemat1_uncond(numdof_uncond,numdof_uncond);
        Epetra_SerialDenseVector elevec1_uncond(numdof_uncond);

        RCP<Epetra_SerialDenseMatrix> Cfi;
        RCP<Epetra_SerialDenseMatrix> Cif;
        RCP<Epetra_SerialDenseMatrix> Cii;
        RCP<Epetra_SerialDenseVector> rhsi;

        if (ih_->ElementIntersected(Id()))
        {
          const size_t ndof_i = ifacepatchlm->size();
          const size_t ndof_ups = eleDofManager_uncondensed_->NumDofElemAndNode();

          fluidfluidmatrices_.Gsui_uncond  = rcp(new Epetra_SerialDenseMatrix(ndof_ups, ndof_i));
          fluidfluidmatrices_.Guis_uncond  = rcp(new Epetra_SerialDenseMatrix(ndof_i, ndof_ups));
          fluidfluidmatrices_.rhsui_uncond = rcp(new Epetra_SerialDenseVector(ndof_i));
          if (monolithic_FSI)
          {
            fluidfluidmatrices_.GNudi_uncond  = rcp(new Epetra_SerialDenseMatrix(ndof_ups, ndof_i));
            fluidfluidmatrices_.GNsdi_uncond  = rcp(new Epetra_SerialDenseMatrix(ndof_ups, ndof_i));
            fluidfluidmatrices_.GNdidi_uncond = rcp(new Epetra_SerialDenseMatrix(ndof_i, ndof_i));
          }

          const size_t ndof_up = lm.size();
          if (ndof_up != eleDofManager_uncondensed_->NumNodeDof() or
              ndof_up != eleDofManager_->NumNodeDof())
            dserror("this is a inconsistency and considered a BUG!");

          Cfi  = rcp(new Epetra_SerialDenseMatrix(ndof_up, ndof_i));
          Cif  = rcp(new Epetra_SerialDenseMatrix(ndof_i, ndof_up));
          Cii  = rcp(new Epetra_SerialDenseMatrix(ndof_i, ndof_i));
          rhsi = rcp(new Epetra_SerialDenseVector(ndof_i));
        }

        const XFEM::AssemblyType assembly_type = XFEM::ComputeAssemblyType(
                *eleDofManager_uncondensed_, NumNode(), NodeIds());

        // calculate element coefficient matrix and rhs
        XFLUID::callSysmat(DRT::INPUT::get<INPAR::XFEM::BoundaryIntegralType>(params, "EMBEDDED_BOUNDARY"), params, assembly_type,
                this, ih_, *eleDofManager_uncondensed_, mystate, iforcecol, elemat1_uncond, elevec1_uncond,
                mat, timealgo, dt, theta, newton, pstab, supg, cstab, ifaceForceContribution, monolithic_FSI, L2, fluidfluidmatrices_);

        const bool stationary_monolithic_FSI = (monolithic_FSI and (timealgo == INPAR::FLUID::timeint_stationary));
        const bool instationary_monolithic_FSI = (monolithic_FSI and (timealgo != INPAR::FLUID::timeint_stationary));

        if (stationary_monolithic_FSI)
        {
          fluidfluidmatrices_.Gsui_uncond->Scale(0.0);
        }
        else if (instationary_monolithic_FSI)
        {
          double theta_iface;
          if (params.get<bool>("interface second order"))
          { theta_iface = 0.5; }
          else
          { theta_iface = 1.0; }
          fluidfluidmatrices_.Gsui_uncond->Scale(1.0/(dt*theta_iface));
        }
        else
        {
          // fluid fluid coupling: everything is ok.
        }

        // condensation
        CondenseElementStressAndStoreOldIterationStep(
            elemat1_uncond, elevec1_uncond,
            *fluidfluidmatrices_.Gsui_uncond,
            *fluidfluidmatrices_.Guis_uncond,
            *fluidfluidmatrices_.rhsui_uncond,
            fluidfluidmatrices_.GNudi_uncond,
            fluidfluidmatrices_.GNsdi_uncond,
            fluidfluidmatrices_.GNdidi_uncond,
            elemat1, elevec1,
            *Cfi, *Cif, *Cii, *rhsi,
            *ifacepatchlm,
            monolithic_FSI,
            monolithic_FSI
            );

        if (ih_->ElementIntersected(Id()))
        {
          DRT::DEBUGGING::NaNChecker(*fluidfluidmatrices_.Guis_uncond);
          DRT::DEBUGGING::NaNChecker(*fluidfluidmatrices_.Gsui_uncond);
          DRT::DEBUGGING::NaNChecker(*fluidfluidmatrices_.rhsui_uncond);
          if (monolithic_FSI)
          {
            DRT::DEBUGGING::NaNChecker(*fluidfluidmatrices_.GNudi_uncond);
            DRT::DEBUGGING::NaNChecker(*fluidfluidmatrices_.GNsdi_uncond);
            DRT::DEBUGGING::NaNChecker(*fluidfluidmatrices_.GNdidi_uncond);
          }
          params.set("Cud",Cfi);
          params.set("Cdu",Cif);
          params.set("Cdd",Cii);
          params.set("rhsd",rhsi);
        }
      }

      DRT::DEBUGGING::NaNChecker(elevec1);
      DRT::DEBUGGING::NaNChecker(elemat1);

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

              std::exit(0);
          }
          else
#endif
	    break;
	  }
    case fluidfluidCoupling:
    {
      // do no calculation, if not needed
      if (lm.empty())
        break;

      const Teuchos::RCP<Epetra_Vector> iforcecol = params.get<Teuchos::RCP<Epetra_Vector> >("interface force");

      double L2 = params.get<double>("L2");

      // time integration factors
      const INPAR::FLUID::TimeIntegrationScheme timealgo = DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(params, "timealgo");
      const double            dt       = params.get<double>("dt");
      const double            theta    = params.get<double>("theta");

      // extract local values from the global vectors
      const bool instationary = (timealgo != INPAR::FLUID::timeint_stationary);

      DRT::ELEMENTS::XFluid3::MyState mystate(discretization,lm,instationary);

      const bool newton = params.get<bool>("include reactive terms for linearisation");
      const bool pstab  = true;
      const bool supg   = true;
      const bool cstab  = true;

      const bool monolithic_FSI = params.get<bool>("monolithic_FSI");

      const RCP<const vector<int> > ifacepatchlm = params.get<RCP<vector<int> > >("ifacepatchlm");

      if (not params.get<bool>("DLM_condensation") or not ih_->ElementIntersected(Id())) // integrate and assemble all unknowns
      {
        if (ih_->ElementIntersected(Id()))
        {
          const size_t nui = ifacepatchlm->size();
          fluidfluidmatrices_.Guis_uncond   = rcp(new Epetra_SerialDenseMatrix(nui, eleDofManager_uncondensed_->NumDofElemAndNode()));
          fluidfluidmatrices_.Gsui_uncond   = rcp(new Epetra_SerialDenseMatrix(eleDofManager_uncondensed_->NumDofElemAndNode(), nui));
          fluidfluidmatrices_.rhsui_uncond = rcp(new Epetra_SerialDenseVector(nui));
        }

        const XFEM::AssemblyType assembly_type = XFEM::ComputeAssemblyType(
            *eleDofManager_, NumNode(), NodeIds());

        // calculate element coefficient matrix and rhs
        XFLUID::callSysmat(DRT::INPUT::get<INPAR::XFEM::BoundaryIntegralType>(params, "EMBEDDED_BOUNDARY"), params, assembly_type,
            this, ih_, *eleDofManager_, mystate, iforcecol, elemat1, elevec1,
            mat, timealgo, dt, theta, newton, pstab, supg, cstab, false, monolithic_FSI, L2, fluidfluidmatrices_);

        if (ih_->ElementIntersected(Id()))
        {
          const size_t nui = ifacepatchlm->size();
          RCP<Epetra_SerialDenseMatrix> Cdd = rcp(new Epetra_SerialDenseMatrix(nui, nui));
          params.set("Cdu",fluidfluidmatrices_.Guis_uncond);
          params.set("Cud",fluidfluidmatrices_.Gsui_uncond);
          params.set("Cdd",Cdd);
          params.set("rhsd",fluidfluidmatrices_.rhsui_uncond);
        }
      }

      else // create bigger element matrix and vector, assemble, condense and copy to small matrix provided by discretization
      {
        // sanity checks
        SanityChecks(eleDofManager_, eleDofManager_uncondensed_);

        const RCP<const Epetra_Vector>  iterincxdomain = discretization.GetState("velpres nodal iterinc");
        const RCP<const Epetra_Vector>  iterinciface   = discretization.GetState("interface nodal iterinc");

        // stress update
        UpdateOldDLMAndDLMRHS(iterincxdomain, iterinciface, lm, *ifacepatchlm, mystate, true);

        // create uncondensed element matrix and vector
        const int numdof_uncond = eleDofManager_uncondensed_->NumDofElemAndNode();
        Epetra_SerialDenseMatrix elemat1_uncond(numdof_uncond,numdof_uncond);
        Epetra_SerialDenseVector elevec1_uncond(numdof_uncond);

        RCP<Epetra_SerialDenseMatrix> Cud;
        RCP<Epetra_SerialDenseMatrix> Cdu;
        RCP<Epetra_SerialDenseMatrix> Cdd;
        RCP<Epetra_SerialDenseVector> rhsd;

        if (ih_->ElementIntersected(Id()))
        {
          const size_t nui = ifacepatchlm->size();
          const size_t nus = eleDofManager_uncondensed_->NumDofElemAndNode();

          fluidfluidmatrices_.Guis_uncond  = rcp(new Epetra_SerialDenseMatrix(nui, nus));
          fluidfluidmatrices_.Gsui_uncond  = rcp(new Epetra_SerialDenseMatrix(nus, nui));
          fluidfluidmatrices_.rhsui_uncond = rcp(new Epetra_SerialDenseVector(nui));
          if (monolithic_FSI)
          {
            fluidfluidmatrices_.GNudi_uncond  = rcp(new Epetra_SerialDenseMatrix(nus, nui));
            fluidfluidmatrices_.GNsdi_uncond  = rcp(new Epetra_SerialDenseMatrix(nus, nui));
            fluidfluidmatrices_.GNdidi_uncond = rcp(new Epetra_SerialDenseMatrix(nui, nui));
          }
          Cud = rcp(new Epetra_SerialDenseMatrix(lm.size(), nui));
          Cdu = rcp(new Epetra_SerialDenseMatrix(nui, lm.size()));
          Cdd = rcp(new Epetra_SerialDenseMatrix(nui, nui));
          rhsd = rcp(new Epetra_SerialDenseVector(nui));
        }

        const XFEM::AssemblyType assembly_type = XFEM::ComputeAssemblyType(
                           *eleDofManager_uncondensed_, NumNode(), NodeIds());

        // calculate element coefficient matrix and rhs
        XFLUID::callSysmat(DRT::INPUT::get<INPAR::XFEM::BoundaryIntegralType>(params, "EMBEDDED_BOUNDARY"), params, assembly_type,
            this, ih_, *eleDofManager_uncondensed_, mystate, iforcecol, elemat1_uncond, elevec1_uncond,
            mat, timealgo, dt, theta, newton, pstab, supg, cstab, false, monolithic_FSI, L2, fluidfluidmatrices_);

        // condensation
        CondenseElementStressAndStoreOldIterationStep(
            elemat1_uncond, elevec1_uncond,
            *fluidfluidmatrices_.Gsui_uncond,
            *fluidfluidmatrices_.Guis_uncond,
            *fluidfluidmatrices_.rhsui_uncond,
            fluidfluidmatrices_.GNudi_uncond,
            fluidfluidmatrices_.GNsdi_uncond,
            fluidfluidmatrices_.GNdidi_uncond,
            elemat1, elevec1,
            *Cud, *Cdu, *Cdd, *rhsd,
            *ifacepatchlm,
            true,
            monolithic_FSI
        );

        if (ih_->ElementIntersected(Id()))
        {
          params.set("Cdu",Cdu);
          params.set("Cud",Cud);
          params.set("Cdd",Cdd);
          params.set("rhsd",rhsd);
        }
       }

      params.set<double>("L2",L2);
      break;
    }
    case fluidxfluidCoupling:
        {
          // do no calculation, if not needed
          if (lm.empty())
            break;

          const Teuchos::RCP<Epetra_Vector> iforcecol = params.get<Teuchos::RCP<Epetra_Vector> >("interface force");

          double L2 = params.get<double>("L2");

          // time integration factors
          const INPAR::FLUID::TimeIntegrationScheme timealgo = DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(params, "timealgo");
          const double            dt       = params.get<double>("dt");
          const double            theta    = params.get<double>("theta");

          // extract local values from the global vectors
          const bool instationary = (timealgo != INPAR::FLUID::timeint_stationary);

          DRT::ELEMENTS::XFluid3::MyState mystate(discretization,lm,instationary);

          const bool newton = params.get<bool>("include reactive terms for linearisation");
          const bool pstab  = true;
          const bool supg   = true;
          const bool cstab  = true;

          const bool monolithic_FSI = params.get<bool>("monolithic_FSI");

          const RCP<const vector<int> > ifacepatchlm = params.get<RCP<vector<int> > >("ifacepatchlm");

          if (not params.get<bool>("DLM_condensation") or not ih_->ElementIntersected(Id())) // integrate and assemble all unknowns
          {
            if (ih_->ElementIntersected(Id()))
            {
              const size_t nui = ifacepatchlm->size();
              fluidfluidmatrices_.Guis_uncond   = rcp(new Epetra_SerialDenseMatrix(nui, eleDofManager_uncondensed_->NumDofElemAndNode()));
              fluidfluidmatrices_.Gsui_uncond   = rcp(new Epetra_SerialDenseMatrix(eleDofManager_uncondensed_->NumDofElemAndNode(), nui));
              fluidfluidmatrices_.rhsui_uncond = rcp(new Epetra_SerialDenseVector(nui));
            }

            const XFEM::AssemblyType assembly_type = XFEM::ComputeAssemblyType(
                *eleDofManager_, NumNode(), NodeIds());

            // calculate element coefficient matrix and rhs
            XFLUID::callSysmat(DRT::INPUT::get<INPAR::XFEM::BoundaryIntegralType>(params, "EMBEDDED_BOUNDARY"), params, assembly_type,
                this, ih_, *eleDofManager_, mystate, iforcecol, elemat1, elevec1,
                mat, timealgo, dt, theta, newton, pstab, supg, cstab, false, monolithic_FSI, L2, fluidfluidmatrices_);

            if (ih_->ElementIntersected(Id()))
            {
              const size_t nui = ifacepatchlm->size();
              RCP<Epetra_SerialDenseMatrix> Cdd = rcp(new Epetra_SerialDenseMatrix(nui, nui));
              params.set("Cdu",fluidfluidmatrices_.Guis_uncond);
              params.set("Cud",fluidfluidmatrices_.Gsui_uncond);
              params.set("Cdd",Cdd);
              params.set("rhsd",fluidfluidmatrices_.rhsui_uncond);
            }
          }

          else // create bigger element matrix and vector, assemble, condense and copy to small matrix provided by discretization
          {
            // sanity checks
            SanityChecks(eleDofManager_, eleDofManager_uncondensed_);

            const RCP<const Epetra_Vector>  iterinciface   = ih_->cutterdis()->GetState("interface nodal iterinc");
            const RCP<const Epetra_Vector>  iterincxdomain = discretization.GetState("velpres nodal iterinc");

            // stress update
            UpdateOldDLMAndDLMRHS(iterincxdomain, iterinciface, lm, *ifacepatchlm, mystate, true);

            // create uncondensed element matrix and vector
            const int numdof_uncond = eleDofManager_uncondensed_->NumDofElemAndNode();
            Epetra_SerialDenseMatrix elemat1_uncond(numdof_uncond,numdof_uncond);
            Epetra_SerialDenseVector elevec1_uncond(numdof_uncond);

            RCP<Epetra_SerialDenseMatrix> Cud;
            RCP<Epetra_SerialDenseMatrix> Cdu;
            RCP<Epetra_SerialDenseMatrix> Cdd;
            RCP<Epetra_SerialDenseVector> rhsd;

            if (ih_->ElementIntersected(Id()))
            {
              const size_t nui = ifacepatchlm->size();
              const size_t nus = eleDofManager_uncondensed_->NumDofElemAndNode();

              fluidfluidmatrices_.Guis_uncond  = rcp(new Epetra_SerialDenseMatrix(nui, nus));
              fluidfluidmatrices_.Gsui_uncond  = rcp(new Epetra_SerialDenseMatrix(nus, nui));
              fluidfluidmatrices_.rhsui_uncond = rcp(new Epetra_SerialDenseVector(nui));
              if (monolithic_FSI)
              {
                fluidfluidmatrices_.GNudi_uncond  = rcp(new Epetra_SerialDenseMatrix(nus, nui));
                fluidfluidmatrices_.GNsdi_uncond  = rcp(new Epetra_SerialDenseMatrix(nus, nui));
                fluidfluidmatrices_.GNdidi_uncond = rcp(new Epetra_SerialDenseMatrix(nui, nui));
              }
              Cud = rcp(new Epetra_SerialDenseMatrix(lm.size(), nui));
              Cdu = rcp(new Epetra_SerialDenseMatrix(nui, lm.size()));
              Cdd = rcp(new Epetra_SerialDenseMatrix(nui, nui));
              rhsd = rcp(new Epetra_SerialDenseVector(nui));
            }

            const XFEM::AssemblyType assembly_type = XFEM::ComputeAssemblyType(
                               *eleDofManager_uncondensed_, NumNode(), NodeIds());

            // calculate element coefficient matrix and rhs
            XFLUID::callSysmat(DRT::INPUT::get<INPAR::XFEM::BoundaryIntegralType>(params, "EMBEDDED_BOUNDARY"), params, assembly_type,
                this, ih_, *eleDofManager_uncondensed_, mystate, iforcecol, elemat1_uncond, elevec1_uncond,
                mat, timealgo, dt, theta, newton, pstab, supg, cstab, false, monolithic_FSI, L2, fluidfluidmatrices_);

            // condensation
            CondenseElementStressAndStoreOldIterationStep(
                elemat1_uncond, elevec1_uncond,
                *fluidfluidmatrices_.Gsui_uncond,
                *fluidfluidmatrices_.Guis_uncond,
                *fluidfluidmatrices_.rhsui_uncond,
                fluidfluidmatrices_.GNudi_uncond,
                fluidfluidmatrices_.GNsdi_uncond,
                fluidfluidmatrices_.GNdidi_uncond,
                elemat1, elevec1,
                *Cud, *Cdu, *Cdd, *rhsd,
                *ifacepatchlm,
                true,
                monolithic_FSI
            );

            if (ih_->ElementIntersected(Id()))
            {
              params.set("Cdu",Cdu);
              params.set("Cud",Cud);
              params.set("Cdd",Cdd);
              params.set("rhsd",rhsd);
            }
           }

          params.set<double>("L2",L2);
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

      const XFEM::AssemblyType assembly_type = XFEM::ComputeAssemblyType(
              *eleDofManager_, NumNode(), NodeIds());

      // calculate element coefficient matrix and rhs
      XFLUID::callSysmatProjection(assembly_type,
              this, ih_, *eleDofManager_, mystate, elemat1, elemat2, elevec1, elevec2,
              pstab, ifaceForceContribution);

#ifdef DEBUG
      if (std::isnan(elevec1.Norm2()))   { cout << *this << endl; dserror("NaNs in elevec1 detected! Quitting..."); }
      if (std::isnan(elemat1.InfNorm())) { cout << *this << endl; dserror("NaNs in elemat1 detected! Quitting..."); }
#endif

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

      const XFEM::AssemblyType assembly_type = XFEM::ComputeAssemblyType(
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

  const double a      = M_PI/4.0;
  const double d      = M_PI/2.0;

  // get viscosity
  double  kinvisc = 0.0;
  double  dens = 0.0;
  if(material->MaterialType() == INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = dynamic_cast<const MAT::NewtonianFluid*>(material.get());
    kinvisc = actmat->Viscosity();
    dens = actmat->Density();
  }
  else
    dserror("Cannot handle material of type %d", material->MaterialType());

  const double dynvisc = kinvisc * dens;

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
        )* exp(-2.0*dynvisc*d*d*t);

    // compute analytical velocities
    u[0] = -a * ( exp(a*xint[0]) * sin(a*xint[1] + d*xint[2]) +
                  exp(a*xint[2]) * cos(a*xint[0] + d*xint[1]) ) * exp(-dynvisc*d*d*t);
    u[1] = -a * ( exp(a*xint[1]) * sin(a*xint[2] + d*xint[0]) +
                  exp(a*xint[0]) * cos(a*xint[1] + d*xint[2]) ) * exp(-dynvisc*d*d*t);
    u[2] = -a * ( exp(a*xint[2]) * sin(a*xint[0] + d*xint[1]) +
                  exp(a*xint[1]) * cos(a*xint[2] + d*xint[0]) ) * exp(-dynvisc*d*d*t);

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
    const RCP<const Epetra_Vector>  iterincxdomain,
    const RCP<const Epetra_Vector>  iterinciface,
    const std::vector<int>&         lm,
    const std::vector<int>&         lmiface,
    MyState&                        mystate,
    const bool                      interface_unknowns
    ) const
{
  const std::size_t nu = eleDofManager_uncondensed_->NumNodeDof();
  const std::size_t ns = eleDofManager_uncondensed_->NumElemDof();
  const std::size_t nui = lmiface.size();

  if (ns > 0)
  {
    if (nu != lm.size())
      dserror("mismatch in fluid velocity and pressure number of DOFs. This is a bug!");

    static const Epetra_BLAS blas;

    // add Kda . inc_velnp to feas
    // new alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas
    {
      vector<double> iterinc_velnp(lm.size());
      DRT::UTILS::ExtractMyValues(*iterincxdomain,iterinc_velnp,lm);
      DRT::DEBUGGING::NaNChecker(iterinc_velnp);

      // update old iteration residual of the stresses from velocity and pressure increments
      // DLM_info_->oldrs_(i) += DLM_info_->oldKsu_(i,j)*inc_velnp[j];
      blas.GEMV('N', ns, nu,-1.0, DLM_info_->oldKGsu_.A(), DLM_info_->oldKGsu_.LDA(), &iterinc_velnp[0], 1.0, DLM_info_->oldrs_.A());
    }

    if (nui > 0 and interface_unknowns)
    {
      vector<double> iterinc_velnp_iface(nui);
      DRT::UTILS::ExtractMyValues(*iterinciface,iterinc_velnp_iface,lmiface);
      DRT::DEBUGGING::NaNChecker(iterinc_velnp_iface);

      // update old iteration residual of the stresses from interface velocity increments
      for (std::size_t i=0;i<ns;i++)
        for (std::size_t j=0;j<nui;j++)
          DLM_info_->oldrs_(i) += DLM_info_->oldGsui_(i,j)*iterinc_velnp_iface[j];
//      blas.GEMV('N', ns, nui,-1.0, DLM_info_->oldGsui_.A(), DLM_info_->oldGsui_.LDA(), &iterinc_velnp_iface[0], 1.0, DLM_info_->oldrs_.A());
    }

    // compute element stresses
    // DLM_info_->stressdofs_(i) -= DLM_info_->oldKssinv_(i,j)*DLM_info_->oldrs_(j);
    blas.GEMV('N', ns, ns,1.0, DLM_info_->oldKssinv_.A(), DLM_info_->oldKssinv_.LDA(), DLM_info_->oldrs_.A(), 1.0, DLM_info_->stressdofs_.A());

    // increase size of element vector (old values stay and zeros are added)
    const std::size_t numdof_uncond = eleDofManager_uncondensed_->NumDofElemAndNode();
    DRT::DEBUGGING::NaNChecker(mystate.velnp);
    mystate.velnp.resize(numdof_uncond,0.0);
    mystate.veln .resize(numdof_uncond,0.0);
    mystate.velnm.resize(numdof_uncond,0.0);
    mystate.accn .resize(numdof_uncond,0.0);
    for (std::size_t i=0;i<ns;i++)
    {
      mystate.velnp[nu+i] = DLM_info_->stressdofs_(i);
    }
//    cout << "DLM_info_->stressdofs_" << endl;
//    cout << DLM_info_->stressdofs_ << endl;
  }
  DRT::DEBUGGING::NaNChecker(DLM_info_->stressdofs_);
  DRT::DEBUGGING::NaNChecker(mystate.velnp);
}

/*---------------------------------------------------------------------*
 *---------------------------------------------------------------------*/
void DRT::ELEMENTS::XFluid3::CondenseElementStressAndStoreOldIterationStep(
    const Epetra_SerialDenseMatrix& elemat1_uncond,
    const Epetra_SerialDenseVector& elevec1_uncond,
    const Epetra_SerialDenseMatrix& Gsui_uncond,
    const Epetra_SerialDenseMatrix& Guis_uncond,
    const Epetra_SerialDenseVector& rhsui_uncond,
    RCP<Epetra_SerialDenseMatrix> GNudi_uncond,
    RCP<Epetra_SerialDenseMatrix> GNsdi_uncond,
    RCP<Epetra_SerialDenseMatrix> GNdidi_uncond,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix& Cuui,
    Epetra_SerialDenseMatrix& Cuiu,
    Epetra_SerialDenseMatrix& Cuiui,
    Epetra_SerialDenseVector& rhsui,
    const std::vector<int>&   lmiface,
    const bool iface_unknowns,
    const bool monolithic_FSI
) const
{
  // for matrix vector and matrix matrix computations
  static const Epetra_BLAS blas;

  const size_t nu = eleDofManager_uncondensed_->NumNodeDof();
  const size_t ns = eleDofManager_uncondensed_->NumElemDof();
  const size_t nui = lmiface.size();

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
    Epetra_SerialDenseSolver solve_for_inverseKss;
    solve_for_inverseKss.SetMatrix(Kssinv);
    solve_for_inverseKss.Invert();


    LINALG::SerialDenseMatrix GusKssinv(nu,ns,true); // temporary Gus.Kss^{-1}

    // GusKssinv(i,j) = Gus(i,k)*Kssinv(k,j);
    blas.GEMM('N','N',nu,ns,ns,1.0,Gus.A(),Gus.LDA(),Kssinv.A(),Kssinv.LDA(),0.0,GusKssinv.A(),GusKssinv.LDA());

    // elemat1(i,j) += - GusKssinv(i,k)*KGsu(k,j);   // note that elemat1 = Cuu below
    blas.GEMM('N','N',nu,nu,ns,-1.0,GusKssinv.A(),GusKssinv.LDA(),KGsu.A(),KGsu.LDA(),1.0,elemat1.A(),elemat1.LDA());

    // elevec1(i) += - GusKssinv(i,j)*fs(j);
    blas.GEMV('N', nu, ns,-1.0, GusKssinv.A(), GusKssinv.LDA(), fs.A(), 1.0, elevec1.A());

    if (iface_unknowns)
    {
      if (nui == 0)
        dserror("think");
      LINALG::SerialDenseMatrix Gsui(ns,nui);
      LINALG::SerialDenseMatrix Guis(nui,ns);
      for (size_t i=0;i<nui;i++)
      {
        rhsui(i) = rhsui_uncond(i);
        for (size_t j=0;j<ns;j++)
        {
          Gsui(j,i) = Gsui_uncond(nu+j, i);
          if (monolithic_FSI)
            Gsui(j,i) += (*GNsdi_uncond)(nu+j, i);
          Guis(i,j) = Guis_uncond(i, nu+j);
        }
      }

      LINALG::SerialDenseMatrix GuisKssinv(nui,ns,true);

      // GuisKssinv(i,j) = Kuis(i,k)*Kssinv(k,j);
      blas.GEMM('N','N',nui,ns,ns,1.0,Guis.A(),Guis.LDA(),Kssinv.A(),Kssinv.LDA(),0.0,GuisKssinv.A(),GuisKssinv.LDA());

      for (size_t i=0;i<nu;i++)
        for (size_t j=0;j<nui;j++)
        {
          if (monolithic_FSI)
            Cuui(i,j) += (*GNudi_uncond)(i,j);
          for (size_t k=0;k<ns;k++)
            Cuui(i,j) += -GusKssinv(i,k)*Gsui(k,j);
        }
      for (size_t i=0;i<nui;i++)
        for (size_t j=0;j<nu;j++)
          for (size_t k=0;k<ns;k++)
            Cuiu(i,j) += -GuisKssinv(i,k)*KGsu(k,j);
      for (size_t i=0;i<nui;i++)
        for (size_t j=0;j<nui;j++)
        {
          if (monolithic_FSI)
            Cuiui(i,j) += (*GNdidi_uncond)(i,j);
          for (size_t k=0;k<ns;k++)
            Cuiui(i,j) += -GuisKssinv(i,k)*Gsui(k,j);
        }
      for (size_t i=0;i<nui;i++)
        for (size_t j=0;j<ns;j++)
          rhsui(i) += - GuisKssinv(i,j)*fs(j);

      //DLM_info_->oldGsui_.Update(1.0,Gsui,0.0);
      blas.COPY(DLM_info_->oldGsui_.M()*DLM_info_->oldGsui_.N(), Gsui.A(), DLM_info_->oldGsui_.A());
    }

    // store current DLM data in iteration history
    //DLM_info_->oldKssinv_.Update(1.0,Kssinv,0.0);
    blas.COPY(DLM_info_->oldKssinv_.M()*DLM_info_->oldKssinv_.N(), Kssinv.A(), DLM_info_->oldKssinv_.A());
    //DLM_info_->oldKsu_.Update(1.0,KGsu,0.0);
    blas.COPY(DLM_info_->oldKGsu_.M()*DLM_info_->oldKGsu_.N(), KGsu.A(), DLM_info_->oldKGsu_.A());
    //DLM_info_->oldrs_.Update(1.0,fa,0.0);
    blas.COPY(DLM_info_->oldrs_.M()*DLM_info_->oldrs_.N(), fs.A(), DLM_info_->oldrs_.A());
  }
  DRT::DEBUGGING::NaNChecker(DLM_info_->oldKssinv_);
  DRT::DEBUGGING::NaNChecker(DLM_info_->oldKGsu_);
  DRT::DEBUGGING::NaNChecker(DLM_info_->oldGsui_);
  DRT::DEBUGGING::NaNChecker(DLM_info_->oldrs_);
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
        XFEM::InterfaceHandleXFSI*  ih,      ///< connection to the interface handler
        const XFEM::ElementDofManager&                  dofman,  ///< dofmanager of the current element
        Epetra_SerialDenseMatrix&                       estif,   ///< element matrix to calculate
        Epetra_SerialDenseVector&                       eforce   ///< element rhs to calculate
        )
{
  // initialize arrays
  estif.Scale(0.0);
  eforce.Scale(0.0);

  const int NUMDOF = 4;

  XFLUID::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF> assembler(dofman, estif, eforce);

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

    const XFEM::ElementEnrichmentValues enrvals(
        *ele,
        ih,
        dofman,
        cellcenter_xyz, false, -1);

    const DRT::UTILS::GaussRule3D gaussrule = XFEM::getXFEMGaussrule<DISTYPE>(ele, xyze, ih->ElementIntersected(ele->Id()),cell->Shape());

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

      static XFEM::ApproxFunc<1,shpVecSize> shp;

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
        XFEM::InterfaceHandleXFSI*  ih,      ///< connection to the interface handler
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

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
