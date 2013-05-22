/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_impl.cpp

\brief Internal implementation of scalar transport boundary elements

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 */
/*----------------------------------------------------------------------*/

// general Butler-Volmer is activated if the define-flag ButlerVolmer_Shifted is off
// the shifted version of Butler-Volmer is activated if the define-flag ButlerVolmer_Shifted is on
// define-flag PERCENT: how much is the curve shifted


#include <cstdlib>
#include "scatra_ele_boundary_impl.H"
#include "scatra_ele_impl_utils.H"
#include "scatra_ele_action.H"
#include "scatra_element.H"

#include "../drt_lib/drt_globalproblem.H" // for curves and functions
#include "../drt_lib/standardtypes_cpp.H" // for EPS12 and so on
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"
#include "../drt_geometry/position_array.H"

// material headers
#include "../drt_mat/scatra_mat.H"
#include "../drt_mat/myocard.H"
#include "../drt_mat/mixfrac.H"
#include "../drt_mat/sutherland.H"
#include "../drt_mat/arrhenius_spec.H"
#include "../drt_mat/arrhenius_temp.H"
#include "../drt_mat/arrhenius_pv.H"
#include "../drt_mat/ferech_pv.H"
#include "../drt_mat/ion.H"
#include "../drt_mat/fourieriso.H"
#include "../drt_mat/thermostvenantkirchhoff.H"
#include "../drt_mat/yoghurt.H"
#include "../drt_mat/matlist.H"

#include "../drt_mat/elchmat.H"
#include "../drt_mat/diffcond.H"
#include "../drt_mat/newman.H"
#include "../drt_mat/elchphase.H"

//#define DEBUG_BATTERY


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraBoundaryImplInterface* DRT::ELEMENTS::ScaTraBoundaryImplInterface::Impl(
    const DRT::Element* ele,
    const enum INPAR::SCATRA::ScaTraType scatratype)
{
  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the computations
  const int numdofpernode = ele->NumDofPerNode(*(ele->Nodes()[0]));
  int numscal = numdofpernode;

  if (SCATRA::IsElchProblem(scatratype))
  {
    numscal -= 1;

    // get the material of the first element
    // we assume here, that the material is equal for all elements in this discretization
    // get the parent element including its material
    const DRT::ELEMENTS::TransportBoundary* transele = static_cast<const DRT::ELEMENTS::TransportBoundary*>(ele);
    Teuchos::RCP<MAT::Material> material = transele->ParentElement()->Material();
    if (material->MaterialType() == INPAR::MAT::m_elchmat)
    {
      const MAT::ElchMat* actmat = dynamic_cast<const MAT::ElchMat*>(material.get());

      if (actmat->Current())
        numscal -= DRT::UTILS::getDimension(transele->ParentElement()->Shape());
    }
  }

  switch (ele->Shape())
  {
  case DRT::Element::quad4:
  {
    return ScaTraBoundaryImpl<DRT::Element::quad4>::Instance(numdofpernode,numscal);
  }
  case DRT::Element::quad8:
  {
    return ScaTraBoundaryImpl<DRT::Element::quad8>::Instance(numdofpernode,numscal);
  }
  case DRT::Element::quad9:
  {
    return ScaTraBoundaryImpl<DRT::Element::quad9>::Instance(numdofpernode,numscal);
  }
  case DRT::Element::tri3:
  {
    return ScaTraBoundaryImpl<DRT::Element::tri3>::Instance(numdofpernode,numscal);
  }
  /*  case DRT::Element::tri6:
  {
    return ScaTraBoundaryImpl<DRT::Element::tri6>::Instance(numdofpernode,numscal);
  }*/
  case DRT::Element::line2:
  {
    return ScaTraBoundaryImpl<DRT::Element::line2>::Instance(numdofpernode,numscal);
  }
  case DRT::Element::line3:
  {
    return ScaTraBoundaryImpl<DRT::Element::line3>::Instance(numdofpernode,numscal);
  }
 /* case DRT::Element::nurbs2:    // 1D nurbs boundary element
  {
    return ScaTraBoundaryImpl<DRT::Element::nurbs2>::Instance(numdofpernode,numscal);
  } */
  case DRT::Element::nurbs3:    // 1D nurbs boundary element
  {
    return ScaTraBoundaryImpl<DRT::Element::nurbs3>::Instance(numdofpernode,numscal);
  }
 /*  case DRT::Element::nurbs4:    // 2D nurbs boundary element
  {
    return ScaTraBoundaryImpl<DRT::Element::nurbs4>::Instance(numdofpernode,numscal);
  } */
  case DRT::Element::nurbs9:    // 2D nurbs boundary element
  {
    return ScaTraBoundaryImpl<DRT::Element::nurbs9>::Instance(numdofpernode,numscal);
  }
  default:
    dserror("Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
    break;
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraBoundaryImpl<distype> * DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::Instance(
    const int numdofpernode,
    const int numscal,
    const bool create
    )
{
  static ScaTraBoundaryImpl<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraBoundaryImpl<distype>(numdofpernode,numscal);
    }
  }
  else
  { // proper destruction
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(0,0,false );
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::ScaTraBoundaryImpl
(int numdofpernode,
    int numscal)
    : numdofpernode_(numdofpernode),
    numscal_(numscal),
    isale_(false),
    is_stationary_(false),
    is_genalpha_(false),
    is_incremental_(true),
    xyze_(true),  // initialize to zero
    weights_(true),
    myknots_(nsd_),
    mypknots_(nsd_+1),
    normalfac_(1.0),
    edispnp_(true),
    diffus_(numscal_,0),
    //valence_(numscal_,0),
    shcacp_(0.0),
    xsi_(true),
    funct_(true),
    deriv_(true),
    derxy_(true),
    normal_(true),
    velint_(true),
    metrictensor_(true),
    thermpress_(0.0),
    equpot_(INPAR::ELCH::equpot_enc),
    eps_(1.0),
    tort_(1.0),
    epstort_(1.0)
    {
        return;
    }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::Evaluate(
    DRT::ELEMENTS::TransportBoundary* ele,
    Teuchos::ParameterList&           params,
    DRT::Discretization&              discretization,
    std::vector<int>&                 lm,
    Epetra_SerialDenseMatrix&         elemat1_epetra,
    Epetra_SerialDenseMatrix&         elemat2_epetra,
    Epetra_SerialDenseVector&         elevec1_epetra,
    Epetra_SerialDenseVector&         elevec2_epetra,
    Epetra_SerialDenseVector&         elevec3_epetra
)
{
  // First, do the things that are needed for all actions:

  // get the parent element including its material
  DRT::ELEMENTS::Transport* parentele = ele->ParentElement();
  Teuchos::RCP<MAT::Material> mat = parentele->Material();

  // the type of scalar transport problem has to be provided for all actions!
  const INPAR::SCATRA::ScaTraType scatratype = DRT::INPUT::get<INPAR::SCATRA::ScaTraType>(params, "scatratype");
  if (scatratype == INPAR::SCATRA::scatratype_undefined)
    dserror("Element parameter SCATRATYPE has not been set!");

  // get node coordinates (we have a nsd_+1 dimensional domain!)
  GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,nen_> >(ele,xyze_);

  // get additional state vector for ALE case: grid displacement
  isale_ = params.get<bool>("isale");
  if (isale_)
  {
    const RCP<Epetra_MultiVector> dispnp = params.get< RCP<Epetra_MultiVector> >("dispnp",Teuchos::null);
    if (dispnp==Teuchos::null) dserror("Cannot get state vector 'dispnp'");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,edispnp_,dispnp,nsd_+1);
    // add nodal displacements
    xyze_ += edispnp_;
  }
  else edispnp_.Clear();

  // Now do the nurbs specific stuff (for isogeometric elements)
  if(DRT::NURBS::IsNurbs(distype))
  {
    // for isogeometric elements --- get knotvectors for parent
    // element and boundary element, get weights
    bool zero_size = DRT::NURBS::GetKnotVectorAndWeightsForNurbsBoundary(
        ele, ele->BeleNumber(), ele->ParentElement()->Id(), discretization, mypknots_, myknots_, weights_, normalfac_);

    // if we have a zero sized element due to a interpolated point -> exit here
    if(zero_size) return(0);
  } // Nurbs specific stuff

  // check for the action parameter
  const SCATRA::BoundaryAction action = DRT::INPUT::get<SCATRA::BoundaryAction>(params,"action");
  switch (action)
  {
  case SCATRA::bd_calc_normal_vectors:
  {
    // access the global vector
    const RCP<Epetra_MultiVector> normals = params.get< RCP<Epetra_MultiVector> >("normal vectors",Teuchos::null);
    if (normals == Teuchos::null) dserror("Could not access vector 'normal vectors'");

    // determine constant outer normal to this element
    GetConstNormal(normal_,xyze_);

    // loop over the element nodes
    for (int j=0;j<nen_;j++)
    {
      const int nodegid = (ele->Nodes()[j])->Id();
      if (normals->Map().MyGID(nodegid) )
      { // OK, the node belongs to this processor

        // scaling to a unit vector is performed on the global level after
        // assembly of nodal contributions since we have no reliable information
        // about the number of boundary elements adjacent to a node
        for (int dim=0; dim<(nsd_+1); dim++)
        {
          normals->SumIntoGlobalValue(nodegid,dim,normal_(dim));
        }
      }
      //else: the node belongs to another processor; the ghosted
      //      element will contribute the right value on that proc
    }

    break;
  }
  case SCATRA::bd_calc_elch_electrode_kinetics:
  {
    if (mat->MaterialType() == INPAR::MAT::m_elchmat)
    {
      Teuchos::ParameterList& diffcondparams_ = params.sublist("DIFFCOND");
      equpot_ = DRT::INPUT::IntegralValue<INPAR::ELCH::EquPot>(diffcondparams_,"EQUPOT");
    }

    // get actual values of transported scalars
    RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

    // extract local values from the global vector
    std::vector<double> ephinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,ephinp,lm);

    // get current condition
    Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition> >("condition");
    if (cond == Teuchos::null) dserror("Cannot access condition 'ElectrodeKinetics'");

    // access parameters of the condition
    const int                 kinetics = cond->GetInt("kinetic model");
    double                    pot0 = cond->GetDouble("pot");
    const int                 curvenum = cond->GetInt("curve");
    const int                 nume = cond->GetInt("e-");
    if(nume < 0)
      dserror("The convention for electrochemical reactions at the electrodes does not allow \n"
          "a negative number of transfered electrones");

    const std::vector<int>*   stoich = cond->GetMutable<std::vector<int> >("stoich");
    if((unsigned int)numscal_ != (*stoich).size())
      dserror("Electrode kinetics: number of stoichiometry coefficients %u does not match"
              " the number of ionic species %d", (*stoich).size(), numscal_);

    // the classical implementations of kinetic electrode models does not support
    // more than one reagent or product!! There are alternative formulations
    // as e.g. Newman (2004), pp. 205, eq. 8.6 with 8.10
    {
      int reactspecies = 0;
      for(int kk=0; kk<numscal_; ++kk)
        reactspecies += abs((*stoich)[kk]);

      if(reactspecies>1 and (kinetics==INPAR::SCATRA::butler_volmer or kinetics == INPAR::SCATRA::butler_volmer_yang1997 or
          kinetics == INPAR::SCATRA::tafel or kinetics == INPAR::SCATRA::linear))
        dserror("Kinetic model Butler-Volmer / Butler-Volmer-Yang / Tafel and Linear: \n"
            "Only one educt and no product is allowed in the implemented version");
    }

    // access input parameter
    const double frt = params.get<double>("frt"); // = F/RT
    if (frt<0.0)
      dserror("A negative factor frt is not possible by definition");

    // get control parameter from parameter list
    bool iselch(true);
    if (not SCATRA::IsElchProblem(scatratype))
      iselch = false;
    const bool   is_stationary = params.get<bool>("using stationary formulation");
    const double time = params.get<double>("total time");
    double       timefac = 1.0;
    double       alphaF  = 1.0;
    double       rhsfac  = 1.0;
    // find out whether we shell use a time curve and get the factor
    // this feature can be also used for stationary "pseudo time loops"
    if (curvenum>=0)
    {
      const double curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
      // adjust potential at metal side accordingly
      pot0 *= curvefac;
    }

    if (iselch) // this is not necessary for secondary current distributions
    {
      if (mat->MaterialType() == INPAR::MAT::m_elchmat)
      {
        const MAT::ElchMat* actmat = static_cast<const MAT::ElchMat*>(mat.get());

        for (int iphase=0; iphase < actmat->NumPhase();++iphase)
        {
          const int phaseid = actmat->PhaseID(iphase);
          Teuchos::RCP<const MAT::Material> singlemat = actmat->PhaseById(phaseid);

          if(singlemat->MaterialType() == INPAR::MAT::m_elchphase)
          {
            const MAT::ElchPhase* actsinglemat = static_cast<const MAT::ElchPhase*>(singlemat.get());

            eps_ = actsinglemat->Epsilon();
            tort_ = actsinglemat->Tortuosity();
            epstort_ =eps_*tort_;

            //if(eps_ != 1.0 or tort_!=1.0)
            //  dserror("It is not clear what happens with the BV condition in case of homogenization");
          }
        }
      }
    }

    const bool calc_status = params.get<bool>("calc_status",false);
    if (!calc_status)
    {
      if (not is_stationary)
      {
        // One-step-Theta:    timefac = theta*dt
        // BDF2:              timefac = 2/3 * dt
        // generalized-alpha: timefac = (gamma*alpha_F/alpha_M) * dt
        timefac = params.get<double>("time factor");
        alphaF = params.get<double>("alpha_F");
        timefac *= alphaF;
        if (timefac < 0.0) dserror("time factor is negative.");
        // for correct scaling of rhs contribution (see below)
        rhsfac =  1.0/alphaF;
      }

      EvaluateElectrodeKinetics(
          ele,
          elemat1_epetra,
          elevec1_epetra,
          ephinp,
          timefac,
          mat,
          cond,
          nume,
          *stoich,
          kinetics,
          pot0,
          frt,
          iselch,
          scatratype
      );

      // realize correct scaling of rhs contribution for gen.alpha case
      // with dt*(gamma/alpha_M) = timefac/alpha_F
      // matrix contributions are already scaled correctly with
      // timefac=dt*(gamma*alpha_F/alpha_M)
      elevec1_epetra.Scale(rhsfac);

    }
    else
    {
      // NOTE: add integral value only for elements which are NOT ghosted!
      if(ele->Owner() == discretization.Comm().MyPID())
      {
        // get actual values of transported scalars
        RCP<const Epetra_Vector> phidtnp = discretization.GetState("timederivative");
        if (phidtnp==Teuchos::null) dserror("Cannot get state vector 'ephidtnp'");
        // extract local values from the global vector
        std::vector<double> ephidtnp(lm.size());
        DRT::UTILS::ExtractMyValues(*phidtnp,ephidtnp,lm);

        if (not is_stationary)
        {
          // One-step-Theta:    timefac = theta*dt
          // BDF2:              timefac = 2/3 * dt
          // generalized-alpha: timefac = (gamma*alpha_F/alpha_M) * dt
          timefac = params.get<double>("time factor");
          alphaF = params.get<double>("alpha_F");
          // realize correct scaling of vector with timefac/alphaF = (gamma/alpha_M) * dt
          //timefac *= alphaF;
          if (timefac < 0.0) dserror("time factor is negative.");
          // for correct scaling of rhs
          rhsfac = 1.0; //1.0/alphaF;
        }

        ElectrodeStatus(
            ele,
            params,
            cond,
            ephinp,
            ephidtnp,
            kinetics,
            *stoich,
            nume,
            pot0,
            frt,
            iselch,
            timefac);
      }
    }
    break;
  }
  case SCATRA::bd_calc_elch_adapt_DC:
  {
    // get current condition
    Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition> >("condition");
    if (cond == Teuchos::null) dserror("Cannot access condition 'ElectrodeKinetics'");

    const int    kinetics = cond->GetInt("kinetic model");

    if(kinetics == INPAR::SCATRA::nernst)
    {
      if (mat->MaterialType() == INPAR::MAT::m_elchmat)
      {
        Teuchos::ParameterList& diffcondparams_ = params.sublist("DIFFCOND");
        equpot_ = DRT::INPUT::IntegralValue<INPAR::ELCH::EquPot>(diffcondparams_,"EQUPOT");
      }

      // get actual values of transported scalars
      RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

      // extract local values from the global vector
      std::vector<double> ephinp(lm.size());
      DRT::UTILS::ExtractMyValues(*phinp,ephinp,lm);

      // access parameters of the condition
      double       pot0 = cond->GetDouble("pot");
      const int    curvenum = cond->GetInt("curve");
      const int    nume = cond->GetInt("e-");
      const double e0 = cond->GetDouble("e0");
      const double c0 = cond->GetDouble("c0");

      if(nume < 0)
        dserror("The convention for electrochemical reactions at the electrodes does not allow \n"
            "a negative number of transfered electrones");

      const std::vector<int>*   stoich = cond->GetMutable<std::vector<int> >("stoich");
      if((unsigned int)numscal_ != (*stoich).size())
        dserror("Electrode kinetics: number of stoichiometry coefficients %u does not match"
                " the number of ionic species %d", (*stoich).size(), numscal_);

      // access input parameter
      const double frt = params.get<double>("frt"); // = F/RT
      if (frt<0.0)
        dserror("A negative factor frt is not possible by definition");

      const double time = params.get<double>("total time");

      // get control parameter from parameter list
      if (not SCATRA::IsElchProblem(scatratype))
        dserror("Only available ELCH");

      if (curvenum>=0)
      {
        const double curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
        // adjust potential at metal side accordingly
        pot0 *= curvefac;
      }

      // concentration values of reactive species at element nodes
      LINALG::Matrix<nen_,1> conc(true);

      int onespecies = 0;

      for(int k=0; k<numscal_; ++k)
      {
        // only the first oxidized species O is considered for statistics
        // statistics of other species result directly from the oxidized species (current density, ...)
        // or need to be implemented (surface concentration, OCV, ...)
        if((*stoich)[k]==0)
          continue;

        if(onespecies!=0)
          dserror("More than one reacting species");

        onespecies +=1;

        double conint = 0.0;
        double equpot = 0.0;
        double area = 0.0;

        for (int inode=0; inode< nen_;++inode)
        {
          conc(inode) = ephinp[inode*numdofpernode_+k];
        }

        // integrations points and weights
        DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

        // loop over integration points
        for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
        {
          const double fac = EvalShapeFuncAndIntFac(intpoints,gpid,ele->Id());

          // elch-specific values at integration point:
          conint += funct_.Dot(conc)*fac;
          area += fac;

        }

        equpot = e0 + (log(conint/area/c0))/(frt*nume);

        for (int vi=0; vi<nen_; ++vi)
          elevec2_epetra[vi*numdofpernode_+numscal_] = pot0 - equpot;

        std::cout << "applied pot.: " << pot0 << "  conc: " << conint/area << "  Nernst pot.: " << equpot << "  new pot.: " << elevec2_epetra[0+numscal_] << std::endl;
      }
    }
    break;
  }
  case SCATRA::bd_calc_loma_therm_press:
  {
    // we dont know the parent element's lm vector; so we have to build it here
    const int nenparent = parentele->NumNode();
    std::vector<int> lmparent(nenparent);
    std::vector<int> lmparentowner;
    std::vector<int> lmparentstride;
    parentele->LocationVector(discretization,lmparent,lmparentowner,lmparentstride);

    // get velocity values at nodes
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("convective velocity field",Teuchos::null);

    // we deal with a (nsd_+1)-dimensional flow field
    Epetra_SerialDenseVector evel((nsd_+1)*nenparent);
    DRT::UTILS::ExtractMyNodeBasedValues(parentele,evel,velocity,nsd_+1);

    // get values of scalar
    RCP<const Epetra_Vector> phinp  = discretization.GetState("phinp");
    if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

    // extract local values from the global vectors for the parent(!) element
    std::vector<double> myphinp(lmparent.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lmparent);

    // define vector for normal diffusive and velocity fluxes
    std::vector<double> mynormdiffflux(lm.size());
    std::vector<double> mynormvel(lm.size());

    // determine constant outer normal to this element
    GetConstNormal(normal_,xyze_);

    // extract temperature flux vector for each node of the parent element
    LINALG::SerialDenseMatrix eflux(3,nenparent);
    DRT::Element* peleptr = (DRT::Element*) parentele;
    int k=numscal_-1;     // temperature is always last degree of freedom!!
    std::ostringstream temp;
    temp << k;
    std::string name = "flux_phi_"+temp.str();
    // try to get the pointer to the entry (and check if type is RCP<Epetra_MultiVector>)
    RCP<Epetra_MultiVector>* f = params.getPtr< RCP<Epetra_MultiVector> >(name);
    // check: field has been set and is not of type Teuchos::null
    if (f!= NULL) DRT::UTILS::ExtractMyNodeBasedValues(peleptr,eflux,*f,3);
    else          dserror("MultiVector %s has not been found!",name.c_str());

    // calculate normal diffusive and velocity flux at each node of the
    // present boundary element
    for (int i=0; i<nen_; ++i)
    {
      for(int j = 0; j<nenparent;++j)
      {
        mynormdiffflux[i] = 0.0;
        mynormvel[i]      = 0.0;
        for (int l=0; l<nsd_+1; l++)
        {
          mynormdiffflux[i] += eflux(l,j)*normal_(l);
          mynormvel[i]      += evel[i*(nsd_+1)+l]*normal_(l);
        }
      }
    }

    // calculate integral of normal diffusive and velocity flux
    // NOTE: add integral value only for elements which are NOT ghosted!
    if(ele->Owner() == discretization.Comm().MyPID())
    {
      NormDiffFluxAndVelIntegral(ele,params,mynormdiffflux,mynormvel);
    }

    break;
  }
  case SCATRA::bd_integrate_shape_functions:
  {
    // NOTE: add area value only for elements which are NOT ghosted!
    const bool addarea = (ele->Owner() == discretization.Comm().MyPID());
    IntegrateShapeFunctions(ele,params,elevec1_epetra,addarea);

    break;
  }
  case SCATRA::bd_calc_Neumann_inflow:
  {
    // get control parameters
    is_stationary_  = params.get<bool>("using stationary formulation");
    is_genalpha_    = params.get<bool>("using generalized-alpha time integration");
    is_incremental_ = params.get<bool>("incremental solver");

    // get time factor and alpha_F if required
    // one-step-Theta:    timefac = theta*dt
    // BDF2:              timefac = 2/3 * dt
    // generalized-alpha: timefac = alphaF * (gamma*/alpha_M) * dt
    double timefac = 1.0;
    double alphaF  = 1.0;
    if (not is_stationary_)
    {
      timefac = params.get<double>("time factor");
      if (is_genalpha_)
      {
        alphaF = params.get<double>("alpha_F");
        timefac *= alphaF;
      }
      if (timefac < 0.0) dserror("time factor is negative.");
    }

    // set thermodynamic pressure
    thermpress_ = 0.0;
    if (scatratype==INPAR::SCATRA::scatratype_loma)
      thermpress_ = params.get<double>("thermodynamic pressure");

    // get values of scalar
    RCP<const Epetra_Vector> phinp  = discretization.GetState("phinp");
    if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

    // extract local values from global vector
    std::vector<double> ephinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,ephinp,lm);

    // we dont know the parent element's lm vector; so we have to build it here
    const int nenparent = parentele->NumNode();
    std::vector<int> lmparent(nenparent);
    std::vector<int> lmparentowner;
    std::vector<int> lmparentstride;
    parentele->LocationVector(discretization,lmparent,lmparentowner,lmparentstride);

    // get velocity values at nodes
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("convective velocity field",Teuchos::null);

    // we deal with a (nsd_+1)-dimensional flow field
    Epetra_SerialDenseVector evel((nsd_+1)*nenparent);
    DRT::UTILS::ExtractMyNodeBasedValues(parentele,evel,velocity,nsd_+1);

    // insert velocity field into element array
    LINALG::Matrix<nsd_+1,nen_> evelnp;
    for (int i=0;i<nen_;++i)
    {
      for (int idim=0 ; idim < nsd_+1 ; idim++)
      {
        evelnp(idim,i) = evel[idim + i*(nsd_+1)];
      }
    }

    NeumannInflow(ele,
                  mat,
                  ephinp,
                  evelnp,
                  elemat1_epetra,
                  elevec1_epetra,
                  timefac,
                  alphaF);

    break;
  }
  case SCATRA::bd_calc_convective_heat_transfer:
  {
    // get control parameters
    is_stationary_  = params.get<bool>("using stationary formulation");
    is_genalpha_    = params.get<bool>("using generalized-alpha time integration");
    is_incremental_ = params.get<bool>("incremental solver");

    // get time factor and alpha_F if required
    // one-step-Theta:    timefac = theta*dt
    // BDF2:              timefac = 2/3 * dt
    // generalized-alpha: timefac = alphaF * (gamma*/alpha_M) * dt
    double timefac = 1.0;
    double alphaF  = 1.0;
    if (not is_stationary_)
    {
      timefac = params.get<double>("time factor");
      if (is_genalpha_)
      {
        alphaF = params.get<double>("alpha_F");
        timefac *= alphaF;
      }
      if (timefac < 0.0) dserror("time factor is negative.");
    }

    // get values of scalar
    RCP<const Epetra_Vector> phinp  = discretization.GetState("phinp");
    if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

    // extract local values from global vector
    std::vector<double> ephinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,ephinp,lm);

    // get condition
    Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition> >("condition");
    if (cond == Teuchos::null)
      dserror("Cannot access condition 'ThermoConvections'!");

    // get heat transfer coefficient and surrounding temperature
    const double heatranscoeff = cond->GetDouble("coeff");
    const double surtemp       = cond->GetDouble("surtemp");

    ConvectiveHeatTransfer(ele,
                           mat,
                           ephinp,
                           elemat1_epetra,
                           elevec1_epetra,
                           heatranscoeff,
                           surtemp,
                           timefac,
                           alphaF);

    break;
  }
  case SCATRA::bd_calc_weak_Dirichlet:
  {
    switch (distype)
    {
      if (numscal_>1) dserror("not yet implemented for more than one scalar\n");

      // 2D:
      case DRT::Element::line2:
      {
        if(ele->ParentElement()->Shape()==DRT::Element::quad4)
        {
          WeakDirichlet<DRT::Element::line2,DRT::Element::quad4>(ele,
                                                                 params,
                                                                 discretization,
                                                                 mat,
                                                                 elemat1_epetra,
                                                                 elevec1_epetra);
        }
        else
        {
          dserror("expected combination quad4/hex8 or line2/quad4 for surface/parent pair");
        }
        break;
      }
      // 3D:
      case DRT::Element::quad4:
      {
        if(ele->ParentElement()->Shape()==DRT::Element::hex8)
        {
          WeakDirichlet<DRT::Element::quad4,DRT::Element::hex8>(ele,
                                                                params,
                                                                discretization,
                                                                mat,
                                                                elemat1_epetra,
                                                                elevec1_epetra);
        }
        else
        {
          dserror("expected combination quad4/hex8 or line2/quad4 for surface/parent pair");
        }
        break;
      }
      default:
      {
        dserror("not implemented yet\n");
      }
    }

    break;
  }
  case SCATRA::bd_calc_surface_permeability:
  {
    // get control parameters
    is_stationary_  = params.get<bool>("using stationary formulation");
    is_genalpha_    = params.get<bool>("using generalized-alpha time integration");
    is_incremental_ = params.get<bool>("incremental solver");

    double timefac = 1.0;
    if (is_genalpha_ or not is_incremental_)
      dserror("calc_surface_permeability: chosen option not available");
    if (not is_stationary_)
    {
      timefac = params.get<double>("time factor");
      if (timefac < 0.0) dserror("time factor is negative.");
    }

    // get values of scalar
    RCP<const Epetra_Vector> phinp  = discretization.GetState("phinp");
    if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

    // extract local values from global vector
    std::vector<double> ephinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,ephinp,lm);

    // get current condition
    Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition> >("condition");
    if (cond == Teuchos::null) dserror("Cannot access condition 'SurfacePermeability'");

    const double perm = cond->GetDouble("permeability coefficient");

    EvaluateSurfacePermeability(
      ele,
      ephinp,
      elemat1_epetra,
      elevec1_epetra,
      timefac,
      perm
      );

    break;
  }
  case SCATRA::bd_calc_TG_outflow:
  {
    // implements the boundary terms for the Taylor Galerkin time integration
    // just available for the level set equation
    switch (distype)
    {
    if (numscal_>1) dserror("not yet implemented for more than one scalar\n");

    // 2D:
    case DRT::Element::line2:
    {
      dserror("combination not implemented for 2D case yet!");
      break;
    }
    // 3D:
    case DRT::Element::quad4:
    {
      if(ele->ParentElement()->Shape()==DRT::Element::hex8)
      {
        TaylorGalerkinBoundaryOutflow<DRT::Element::quad4,DRT::Element::hex8>(ele,
            params,
            discretization,
            mat,
            elemat1_epetra,
            elevec1_epetra);
      }
      else
      {
        dserror("expected combination quad4/hex8 or line2/quad4 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::quad8:
    {
      if(ele->ParentElement()->Shape()==DRT::Element::hex20)
      {
        TaylorGalerkinBoundaryOutflow<DRT::Element::quad8,DRT::Element::hex20>(ele,
            params,
            discretization,
            mat,
            elemat1_epetra,
            elevec1_epetra);
      }
      else
      {
        dserror("expected combination quad8/hex20 for surface/parent pair");
      }
      break;
    }
    default:
    {
      dserror("not implemented yet\n");
    }
    }

    break;
  }
  case SCATRA::bd_reinitialize_levelset:
  {
    // implements the boundary terms for the implicit Characteristic Galerkin time integration
    // just available for the reinitialization equation
    switch (distype)
    {
    if (numscal_>1) dserror("not yet implemented for more than one scalar\n");

    // 2D:
    case DRT::Element::line2:
    {
      dserror("combination not implemented for 2D case yet!");
      break;
    }
    // 3D:
    case DRT::Element::quad4:
    {
      if(ele->ParentElement()->Shape()==DRT::Element::hex8)
      {
        ReinitCharacteristicGalerkinBoundary<DRT::Element::quad4,DRT::Element::hex8>(ele,
            params,
            discretization,
            mat,
            elemat1_epetra,
            elevec1_epetra);
      }
      else
      {
        dserror("expected combination quad4/hex8 or line2/quad4 for surface/parent pair");
      }
      break;
    }
    case DRT::Element::quad8:
    {
      if(ele->ParentElement()->Shape()==DRT::Element::hex20)
      {
        ReinitCharacteristicGalerkinBoundary<DRT::Element::quad8,DRT::Element::hex20>(ele,
            params,
            discretization,
            mat,
            elemat1_epetra,
            elevec1_epetra);
      }
      else
      {
        dserror("expected combination quad8/hex20 surface/parent pair");
      }
      break;
    }
    default:
    {
      dserror("not implemented yet\n");
    }
    }

    break;
  }
  case SCATRA::bd_add_convective_mass_flux:
  {
    //calculate integral of convective mass/heat flux
    // NOTE: since results are added to a global vector via normal assembly
    //       it would be wrong to suppress results for a ghosted boundary!

    // get actual values of transported scalars
    RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

    // extract local values from the global vector
    std::vector<double> ephinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,ephinp,lm);

    // get velocity values at nodes
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",Teuchos::null);

    // we deal with a (nsd_+1)-dimensional flow field
    LINALG::Matrix<nsd_+1,nen_>  evel(true);
    DRT::UTILS::ExtractMyNodeBasedValues(ele,evel,velocity,nsd_+1);

    // for the moment we ignore the return values of this method
    CalcConvectiveFlux(ele,ephinp,evel,elevec1_epetra);
    //vector<double> locfluxintegral = CalcConvectiveFlux(ele,ephinp,evel,elevec1_epetra);
    //cout<<"locfluxintegral[0] = "<<locfluxintegral[0]<<endl;

    // NOTE: add value only for boundary elements which are NOT ghosted!

    // if the flux integral for the whole boundary shell be computed as well
    // by summing up locfluxintegral[k] for each element, then, of course,
    // values only for boundary elements which are NOT ghosted should be summed up
    // and added to the parameter list for transport to the outside world
    //    if(ele->Owner() == discretization.Comm().MyPID())

    break;
  }
  default:
  {
    dserror("Not acting on this boundary action. Forgot implementation?");
    break;
  }
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Surface/Line Neumann boundary condition       gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::EvaluateNeumann(
    DRT::ELEMENTS::TransportBoundary*   ele,
    Teuchos::ParameterList&             params,
    DRT::Discretization&                discretization,
    DRT::Condition&                     condition,
    std::vector<int>&                   lm,
    Epetra_SerialDenseVector&           elevec1)
{
  // get node coordinates (we have a nsd_+1 dimensional computational domain!)
  GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,nen_> >(ele,xyze_);

  // get additional state vector for ALE case: grid displacement
  isale_ = params.get<bool>("isale");
  if (isale_)
  {
    const RCP<Epetra_MultiVector> dispnp = params.get< RCP<Epetra_MultiVector> >("dispnp",Teuchos::null);
    if (dispnp==Teuchos::null) dserror("Cannot get state vector 'dispnp'");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,edispnp_,dispnp,nsd_+1);
    // add nodal displacements to point coordinates
    xyze_ += edispnp_;
  }
  else edispnp_.Clear();

  // Now do the nurbs specific stuff (for isogeometric elements)
  if(DRT::NURBS::IsNurbs(distype))
  {
    // for isogeometric elements --- get knotvectors for parent
    // element and boundary element, get weights
    bool zero_size = DRT::NURBS::GetKnotVectorAndWeightsForNurbsBoundary(
        ele, ele->BeleNumber(), ele->ParentElement()->Id(), discretization, mypknots_, myknots_, weights_, normalfac_);

    // if we have a zero sized element due to a interpolated point -> exit here
    if(zero_size) return(0);
  } // Nurbs specific stuff

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const std::vector<int>* curve  = condition.Get<std::vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // get values, switches and spatial functions from the condition
  // (assumed to be constant on element boundary)
  const std::vector<int>*    onoff = condition.Get<std::vector<int> >   ("onoff");
  const std::vector<double>* val   = condition.Get<std::vector<double> >("val"  );
  const std::vector<int>*    func  = condition.Get<std::vector<int> >   ("funct");

  // integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    double fac = EvalShapeFuncAndIntFac(intpoints,iquad,ele->Id());

    // multiply integration factor with the timecurve factor
    fac *= curvefac;

    // factor given by spatial function
    double functfac = 1.0;
    // determine global coordinates of current Gauss point
    double coordgp[3]; // we always need three coordinates for function evaluation!
    for (int i = 0; i< 3; i++)
      coordgp[i] = 0.0;

    for (int i = 0; i< nsd_; i++)
    {
      coordgp[i] = 0.0;
      for (int j = 0; j < nen_; j++)
      {
        coordgp[i] += xyze_(i,j) * funct_(j);
      }
    }

    int functnum = -1;
    const double* coordgpref = &coordgp[0]; // needed for function evaluation


    for(int dof=0;dof<numdofpernode_;dof++)
    {
      if ((*onoff)[dof]) // is this dof activated?
      {
        // factor given by spatial function
        if (func) functnum = (*func)[dof];
        {
          if (functnum>0)
          {
            // evaluate function at current gauss point (provide always 3D coordinates!)
            functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dof,coordgpref,time,NULL);
          }
          else
            functfac = 1.0;
        }

        const double val_fac_functfac = (*val)[dof]*fac*functfac;

        for (int node=0;node<nen_;++node)
        {
          elevec1[node*numdofpernode_+dof] += funct_(node)*val_fac_functfac;
        }
      } // if ((*onoff)[dof])
    }
  } //end of loop over integration points


  // **********************************************************************
  // add boundary flux contributions to the potential equation as well!
  // **********************************************************************
  const INPAR::SCATRA::ScaTraType scatratype = DRT::INPUT::get<INPAR::SCATRA::ScaTraType>(params, "scatratype");
  if (scatratype == INPAR::SCATRA::scatratype_undefined)
    dserror("Element parameter SCATRATYPE has not been set!");

  // this has to be done only for the following problem formulations:
  if ((scatratype==INPAR::SCATRA::scatratype_elch_enc_pde) or
      (scatratype==INPAR::SCATRA::scatratype_elch_enc_pde_elim))
  {
    // access the parent element's material
    DRT::ELEMENTS::Transport* parentele = ele->ParentElement();
    Teuchos::RCP<MAT::Material> material = parentele->Material();

    for (int k = 0; k < numscal_; k++)
    {
      // get valence
      double valence_k(0.0);
      if (material->MaterialType() == INPAR::MAT::m_matlist)
      {
        const MAT::MatList* actmat = static_cast<const MAT::MatList*>(material.get());

        const int matid = actmat->MatID(k);
        Teuchos::RCP<const MAT::Material> singlemat = actmat->MaterialById(matid);
        if (singlemat->MaterialType() == INPAR::MAT::m_ion)
        {
          const MAT::Ion* actsinglemat = static_cast<const MAT::Ion*>(singlemat.get());
          valence_k = actsinglemat->Valence();
        }
        else
          dserror("single material type is not 'ion'");
      }
      else
        dserror("material type is not a 'matlist' material");

      // get corresponding Neumann values, multiply with z_k and add to
      // the row of the electric potential equation
      double val(0.0);
      for (int vi=0; vi<nen_; ++vi)
      {
        val = elevec1[vi*numdofpernode_+k];
        elevec1[vi*numdofpernode_+numscal_] += valence_k*val;
      }
    } // loop over scalars
  }

  // get the parent element including its material
  DRT::ELEMENTS::Transport* parentele = ele->ParentElement();
  Teuchos::RCP<MAT::Material> mat = parentele->Material();

  if (mat->MaterialType() == INPAR::MAT::m_elchmat)
  {
    Teuchos::ParameterList& diffcondparams_ = params.sublist("DIFFCOND");
    equpot_ = DRT::INPUT::IntegralValue<INPAR::ELCH::EquPot>(diffcondparams_,"EQUPOT");
  }

  // the same procedure is also necessary for the concentrated solution theory based on div i
  if(equpot_==INPAR::ELCH::equpot_divi)
  {
    for (int k = 0; k < numscal_; k++)
    {
      // get valence
      double valence_k(0.0);
      if (mat->MaterialType() == INPAR::MAT::m_elchmat)
      {
        const MAT::ElchMat* actmat = static_cast<const MAT::ElchMat*>(mat.get());

        const int specid = actmat->SpecID(k);
        Teuchos::RCP<const MAT::Material> singlemat = actmat->SpecById(specid);
        if (singlemat->MaterialType() == INPAR::MAT::m_diffcond)
        {
          const MAT::DiffCond* actsinglemat = static_cast<const MAT::DiffCond*>(singlemat.get());
          valence_k = actsinglemat->Valence();
          if (abs(valence_k)< EPS14) dserror ("division by zero charge number");
        }
        else if  (singlemat->MaterialType() == INPAR::MAT::m_newman)
        {
          const MAT::Newman* actsinglemat = static_cast<const MAT::Newman*>(singlemat.get());
          valence_k = actsinglemat->Valence();
          if (abs(valence_k)< EPS14) dserror ("division by zero charge number");
        }
        else
          dserror("");
      }
      else
        dserror("material type is not a 'matlist' material");

      // get corresponding Neumann values, multiply with z_k and add to
      // the row of the electric potential equation
      double val(0.0);
      for (int vi=0; vi<nen_; ++vi)
      {
        val = elevec1[vi*numdofpernode_+k];
        elevec1[vi*numdofpernode_+numscal_] += valence_k*val;
      }
    } // loop over scalars
  }

  return 0;
}


/*----------------------------------------------------------------------*
 | calculate Neumann inflow boundary conditions                vg 03/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::NeumannInflow(
    const DRT::Element*                 ele,
    Teuchos::RCP<const MAT::Material>   material,
    const std::vector<double>&           ephinp,
    const LINALG::Matrix<nsd_+1,nen_>&  evelnp,
    Epetra_SerialDenseMatrix&           emat,
    Epetra_SerialDenseVector&           erhs,
    const double                        timefac,
    const double                        alphaF)
{
  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // define vector for scalar values at nodes
  LINALG::Matrix<nen_,1> phinod(true);

  // loop over all scalars
  for(int k=0;k<numdofpernode_;++k)
  {
    // compute scalar values at nodes
    for (int inode=0; inode< nen_;++inode)
    {
      phinod(inode) = ephinp[inode*numdofpernode_+k];
    }

    // loop over all integration points
    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {
      const double fac = EvalShapeFuncAndIntFac(intpoints,iquad,ele->Id(),&normal_);

      // get velocity at integration point
      velint_.Multiply(evelnp,funct_);

      // normal velocity
      const double normvel = velint_.Dot(normal_);

      if (normvel<-0.0001)
      {
        // set density to 1.0
        double dens = 1.0;

        // get density if not equal one
        if (material->MaterialType() == INPAR::MAT::m_matlist)
        {
          const MAT::MatList* actmat = static_cast<const MAT::MatList*>(material.get());

          const int matid = actmat->MatID(0);
          Teuchos::RCP<const MAT::Material> singlemat = actmat->MaterialById(matid);

          if (singlemat->MaterialType() == INPAR::MAT::m_arrhenius_temp)
          {
            const MAT::ArrheniusTemp* actsinglemat = static_cast<const MAT::ArrheniusTemp*>(singlemat.get());

            // compute temperature values at nodes (always last scalar)
            LINALG::Matrix<nen_,1> tempnod(true);
            for (int inode=0; inode< nen_;++inode)
            {
              tempnod(inode) = ephinp[(inode+1)*numdofpernode_-1];
            }

            // compute temperature
            const double temp = funct_.Dot(tempnod);

            // compute density based on temperature and thermodynamic pressure
            dens = actsinglemat->ComputeDensity(temp,thermpress_);
          }
          else dserror("type of material found in material list is not supported");
        }
        else if (material->MaterialType() == INPAR::MAT::m_scatra)
        {
          // set density
          dens = 1.0;
        }
        else if (material->MaterialType() == INPAR::MAT::m_mixfrac)
        {
          const MAT::MixFrac* actmat = static_cast<const MAT::MixFrac*>(material.get());

          // compute mixture fraction
          const double mixfrac = funct_.Dot(phinod);

          // compute density based on mixture fraction
          dens = actmat->ComputeDensity(mixfrac);
        }
        else if (material->MaterialType() == INPAR::MAT::m_sutherland)
        {
          const MAT::Sutherland* actmat = static_cast<const MAT::Sutherland*>(material.get());

          // compute temperature
          const double temp = funct_.Dot(phinod);

          // compute density based on temperature and thermodynamic pressure
          dens = actmat->ComputeDensity(temp,thermpress_);
        }
        else if (material->MaterialType() == INPAR::MAT::m_arrhenius_pv)
        {
          const MAT::ArrheniusPV* actmat = static_cast<const MAT::ArrheniusPV*>(material.get());

          // compute progress variable
          const double provar = funct_.Dot(phinod);

          // compute density
          dens = actmat->ComputeDensity(provar);
        }
        else if (material->MaterialType() == INPAR::MAT::m_ferech_pv)
        {
          const MAT::FerEchPV* actmat = static_cast<const MAT::FerEchPV*>(material.get());

          // compute progress variable
          const double provar = funct_.Dot(phinod);

          // compute density
          dens = actmat->ComputeDensity(provar);
        }
        else if (material->MaterialType() == INPAR::MAT::m_yoghurt)
        {
          const MAT::Yoghurt* actmat = static_cast<const MAT::Yoghurt*>(material.get());

          // get constant density
          dens = actmat->Density();
        }
        else dserror("Material type is not supported for Neumann inflow!");

        // integration factor for left-hand side
        const double lhsfac = dens*normvel*timefac*fac;

        // integration factor for right-hand side
        double rhsfac = 0.0;
        if (is_incremental_ and is_genalpha_)
          rhsfac = lhsfac/alphaF;
        else if (not is_incremental_ and is_genalpha_)
          rhsfac = lhsfac*(1.0-alphaF)/alphaF;
        else if (is_incremental_ and not is_genalpha_)
          rhsfac = lhsfac;

        // matrix
        for (int vi=0; vi<nen_; ++vi)
        {
          const double vlhs = lhsfac*funct_(vi);

          const int fvi = vi*numdofpernode_+k;

          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            emat(fvi,fui) -= vlhs*funct_(ui);
          }
        }

        // scalar at integration point
        const double phi = funct_.Dot(phinod);

        // rhs
        const double vrhs = rhsfac*phi;
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          erhs[fvi] += vrhs*funct_(vi);
        }
      }
    }
  }

  return;
} //ScaTraBoundaryImpl<distype>::NeumannInflow


/*----------------------------------------------------------------------*
 | calculate integral of convective flux across boundary      gjb 11/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
std::vector<double> DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::CalcConvectiveFlux(
    const DRT::Element*                 ele,
    const std::vector<double>&          ephinp,
    const LINALG::Matrix<nsd_+1,nen_>&  evelnp,
    Epetra_SerialDenseVector&           erhs
)
{
  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // define vector for scalar values at nodes
  LINALG::Matrix<nen_,1> phinod(true);

  std::vector<double> integralflux(numscal_);

  // loop over all scalars
  for(int k=0;k<numscal_;++k)
  {
    integralflux[k] = 0.0;

    // compute scalar values at nodes
    for (int inode=0; inode< nen_;++inode)
    {
      phinod(inode) = ephinp[inode*numdofpernode_+k];
    }

    // loop over all integration points
    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {
      const double fac = EvalShapeFuncAndIntFac(intpoints,iquad,ele->Id(),&normal_);

      // get velocity at integration point
      velint_.Multiply(evelnp,funct_);

      // normal velocity (note: normal_ is already a unit(!) normal)
      const double normvel = velint_.Dot(normal_);

      // scalar at integration point
      const double phi = funct_.Dot(phinod);

      const double val = phi*normvel*fac;
      integralflux[k] += val;
      // add contribution to provided vector (distribute over nodes using shape fct.)
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;
        erhs[fvi] += val*funct_(vi);
      }
    }
  }

  return integralflux;

} //ScaTraBoundaryImpl<distype>::ConvectiveFlux

/*----------------------------------------------------------------------*
 | calculate boundary cond. due to convective heat transfer    vg 10/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::ConvectiveHeatTransfer(
    const DRT::Element*                 ele,
    Teuchos::RCP<const MAT::Material>   material,
    const std::vector<double>&           ephinp,
    Epetra_SerialDenseMatrix&           emat,
    Epetra_SerialDenseVector&           erhs,
    const double                        heatranscoeff,
    const double                        surtemp,
    const double                        timefac,
    const double                        alphaF)
{
  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // define vector for scalar values at nodes
  LINALG::Matrix<nen_,1> phinod(true);

  // loop over all scalars
  for(int k=0;k<numdofpernode_;++k)
  {
    // compute scalar values at nodes
    for (int inode=0; inode< nen_;++inode)
    {
      phinod(inode) = ephinp[inode*numdofpernode_+k];
    }

    // loop over all integration points
    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {
      const double fac = EvalShapeFuncAndIntFac(intpoints,iquad,ele->Id(),&normal_);

      // get specific heat capacity at constant volume
      double shc = 0.0;
      if (material->MaterialType() == INPAR::MAT::m_th_fourier_iso)
      {
        const MAT::FourierIso* actmat = static_cast<const MAT::FourierIso*>(material.get());

        shc = actmat->Capacity();
      }
      else if (material->MaterialType() == INPAR::MAT::m_thermostvenant)
      {
        const MAT::ThermoStVenantKirchhoff* actmat = static_cast<const MAT::ThermoStVenantKirchhoff*>(material.get());

        shc = actmat->Capacity();
      }
      else dserror("Material type is not supported for convective heat transfer!");

      // integration factor for left-hand side
      const double lhsfac = heatranscoeff*timefac*fac/shc;

      // integration factor for right-hand side
      double rhsfac = 0.0;
      if (is_incremental_ and is_genalpha_)
        rhsfac = lhsfac/alphaF;
      else if (not is_incremental_ and is_genalpha_)
        rhsfac = lhsfac*(1.0-alphaF)/alphaF;
      else if (is_incremental_ and not is_genalpha_)
        rhsfac = lhsfac;

      // matrix
      for (int vi=0; vi<nen_; ++vi)
      {
        const double vlhs = lhsfac*funct_(vi);

        const int fvi = vi*numdofpernode_+k;

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          emat(fvi,fui) -= vlhs*funct_(ui);
        }
      }

      // scalar at integration point
      const double phi = funct_.Dot(phinod);

      // rhs
      const double vrhs = rhsfac*(phi-surtemp);
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        erhs[fvi] += vrhs*funct_(vi);
      }
    }
  }

  return;
} //ScaTraBoundaryImpl<distype>::ConvectiveHeatTransfer


/*----------------------------------------------------------------------*
 | evaluate shape functions and int. factor at int. point     gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::EvalShapeFuncAndIntFac(
    const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
    const int                                    iquad,      ///< id of current Gauss point
    const int                                    eleid,      ///< the element id
    LINALG::Matrix<1 + nsd_,1>*         normalvec ///< normal vector at Gauss point(optional)
)
{
  // coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (int idim=0;idim<nsd_;idim++)
  {xsi_(idim) = gpcoord[idim];}

  if(not DRT::NURBS::IsNurbs(distype))
  {
    // shape functions and their first derivatives
    DRT::UTILS::shape_function<distype>(xsi_,funct_);
    DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);
  }
  else // nurbs elements are always somewhat special...
  {
    DRT::NURBS::UTILS::nurbs_get_funct_deriv(
        funct_  ,
        deriv_  ,
        xsi_    ,
        myknots_,
        weights_,
        distype );
  }

  // the metric tensor and the area of an infinitesimal surface/line element
  // optional: get normal at integration point as well
  // Note: this is NOT yet a unit normal. Its norm corresponds to the area/length of the element
  double drs(0.0);
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>(xyze_,deriv_,metrictensor_,drs,normalvec);

  // for nurbs elements the normal vector must be scaled with a special orientation factor!!
  if(DRT::NURBS::IsNurbs(distype))
  {
    if (normalvec != NULL)
      normal_.Scale(normalfac_);
  }

  // return the integration factor
  return intpoints.IP().qwgt[iquad] * drs;
}


/*----------------------------------------------------------------------*
 | evaluate an electrode kinetics boundary condition (private) gjb 01/09|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::EvaluateElectrodeKinetics(
    const DRT::Element*        ele,
    Epetra_SerialDenseMatrix& emat,
    Epetra_SerialDenseVector& erhs,
    const std::vector<double>& ephinp,
    double  timefac,
    Teuchos::RCP<const MAT::Material> material,
    Teuchos::RCP<DRT::Condition>      cond,
    const int                 nume,
    const std::vector<int>  stoich,
    const int             kinetics,
    const double              pot0,
    const double               frt,
    const bool              iselch,
    const INPAR::SCATRA::ScaTraType scatratype
)
{
  //for pre-multiplication of i0 with 1/(F z_k)
  double faraday = INPAR::SCATRA::faraday_const;    // unit of F: C/mol or mC/mmol or muC / mumol

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // concentration values of reactive species at element nodes
  std::vector<LINALG::Matrix<nen_,1> > conreact(numscal_);

  // el. potential values at element nodes
  LINALG::Matrix<nen_,1> pot(true);
  if(iselch)
  {
    for (int inode=0; inode< nen_;++inode)
    {
      for(int kk=0; kk<numscal_; kk++)
        conreact[kk](inode) = ephinp[inode*numdofpernode_+kk];

      pot(inode) = ephinp[inode*numdofpernode_+numscal_];
    }
  }
  else
  {
    for (int inode=0; inode< nen_;++inode)
    {
      // vector conreact is not initialized
      for(int kk=0; kk<numscal_; kk++)
        conreact[kk](inode) = 1.0;
      pot(inode) = ephinp[inode*numdofpernode_];
    }
  }

  // concentration of active species at integration point
  std::vector<double> conint(numscal_,0.0);
  // el. potential at integration point
  double potint(0.0);
  // a 'working variable'
  double fac_fns_i0_funct_vi(0.0);

  // index of reactive species (starting from zero)
  // convention: n need to be positiv
  //
  // Sum_i (s_i  M_i^(z_i)) -> n e-

  // loop over all scalars
  for(int k=0; k<numscal_; ++k)
  {
    if(stoich[k]==0)
      continue;

    //(- N^(d+m)*n) = j = s_k / (nume * faraday * z_e-) * i
    //                  = s_k / (nume * faraday * (-1)) * i
    //                    |_______fns__________|
    // see, e.g. in Ehrl et al., "A computational approach for the simulation of natural convection in
    // electrochemical cells", JCP, 2012
    double fns = -1.0/faraday/nume;
    //stichometry as a consequence of the reaction convention
    fns*=stoich[k];

    // only used as an sanity check!!
    double valence_k(0.0);
    if (iselch) // this is not necessary for secondary current distributions
    {
      // get valence of the single(!) reactant
      if (material->MaterialType() == INPAR::MAT::m_matlist)
      {
        const MAT::MatList* actmat = static_cast<const MAT::MatList*>(material.get());

        const int matid = actmat->MatID(k);
        Teuchos::RCP<const MAT::Material> singlemat = actmat->MaterialById(matid);
        if (singlemat->MaterialType() == INPAR::MAT::m_ion)
        {
          const MAT::Ion* actsinglemat = static_cast<const MAT::Ion*>(singlemat.get());
          valence_k = actsinglemat->Valence();
          if (abs(valence_k)< EPS14) dserror ("division by zero charge number");
        }
        else
          dserror("single material type is not 'ion'");
      }
      else if (material->MaterialType() == INPAR::MAT::m_elchmat)
      {
        const MAT::ElchMat* actmat = static_cast<const MAT::ElchMat*>(material.get());

        const int specid = actmat->SpecID(k);
        Teuchos::RCP<const MAT::Material> singlemat = actmat->SpecById(specid);
        if (singlemat->MaterialType() == INPAR::MAT::m_diffcond)
        {
          const MAT::DiffCond* actsinglemat = static_cast<const MAT::DiffCond*>(singlemat.get());
          valence_k = actsinglemat->Valence();
          if (abs(valence_k)< EPS14) dserror ("division by zero charge number");
        }
        else if  (singlemat->MaterialType() == INPAR::MAT::m_newman)
        {
          const MAT::Newman* actsinglemat = static_cast<const MAT::Newman*>(singlemat.get());
          valence_k = actsinglemat->Valence();
          if (abs(valence_k)< EPS14) dserror ("division by zero charge number");
        }
        else
          dserror("");
      }
      else
        dserror("material type is not a 'matlist' material");
    }

 /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    const double fac = EvalShapeFuncAndIntFac(intpoints,gpid,ele->Id());

    // elch-specific values at integration point:
    // concentration is evaluated at all GP since some reaction models depend on all concentrations
    for(int kk=0;kk<numscal_;++kk)
      conint[kk] = funct_.Dot(conreact[kk]);
    potint = funct_.Dot(pot);

    // electrode potential difference (epd) at integration point
    const double epd = (pot0 - potint);

    if (iselch) // tertiary current distribution
    {
      // concentration-dependent Butler-Volmer law(s)
      switch(kinetics)
      {
      case INPAR::SCATRA::butler_volmer:
      case INPAR::SCATRA::butler_volmer_yang1997:
      {
        // read model-specific parameter
        const double alphaa = cond->GetDouble("alpha_a");
        const double alphac = cond->GetDouble("alpha_c");
        double       i0 = cond->GetDouble("i0");
        if (i0 < -EPS14) dserror("i0 is negative, \n"
                                 "a positive definition is necessary due to generalized reaction models: %f",i0);
        // add time factor
        i0*=timefac;
        const double gamma = cond->GetDouble("gamma");
        const double refcon = cond->GetDouble("refcon");
        if (refcon < EPS12) dserror("reference concentration is too small: %f",refcon);

        if(valence_k!=nume)
          dserror("Kinetic model Butler-Volmer: The number of transfered electrodes need to be  \n "
              "the same as the charge number of the reacting species %i", k);
        // opencircuit potential is assumed to be zero
        const double ocp = 0.0;
        // overpotential based on opencircuit potential
        const double eta = epd - ocp;

# if 0
        // print all parameters read from the current condition
        cout<<"kinetic model  = "<<*kinetics<<endl;
        cout<<"react. species = "<<speciesid<<endl;
        cout<<"pot0(mod.)     = "<<pot0<<endl;
        cout<<"curvenum       = "<<curvenum<<endl;
        cout<<"alpha_a        = "<<alphaa<<endl;
        cout<<"alpha_c        = "<<alphac<<endl;
        cout<<"i0(mod.)       = "<<i0<<endl;
        cout<<"gamma          = "<<gamma<<endl;
        cout<<"refcon         = "<<refcon<<endl;
        cout<<"F/RT           = "<<frt<<endl<<endl;
        cout<<"time factor    = "<<timefac<<endl;
        cout<<"alpha_F        = "<<alphaF<<endl;
#endif

#ifdef DEBUG
        // some safety checks/ user warnings
        if ((alphaa*frt*eta) > 100.0)
          cout<<"WARNING: Exp(alpha_a...) in Butler-Volmer law is near overflow!"
          <<exp(alphaa*frt*eta)<<endl;
        if (((-alphac)*frt*eta) > 100.0)
          cout<<"WARNING: Exp(alpha_c...) in Butler-Volmer law is near overflow!"
          <<exp((-alphac)*frt*eta)<<endl;
#endif
        double pow_conint_gamma_k = 0.0;
        if ((conint[k]/refcon) < EPS13)
        {
          pow_conint_gamma_k = std::pow(EPS13,gamma);
#ifdef DEBUG
          cout<<"WARNING: Rel. Conc. in Butler-Volmer formula is zero/negative: "<<(conint[k]/refcon)<<endl;
          cout<<"-> Replacement value: pow(EPS,gamma) = "<< pow_conint_gamma_k <<endl;
#endif
        }
        else
          pow_conint_gamma_k = std::pow(conint[k]/refcon,gamma);

        if (kinetics==INPAR::SCATRA::butler_volmer)
        {
          // note: gamma==0 deactivates concentration dependency in Butler-Volmer!
          const double expterm = exp(alphaa*frt*eta)-exp((-alphac)*frt*eta);

          double concterm = 0.0;
          if (conint[k] > EPS13)
            concterm = gamma*pow(conint[k],(gamma-1.0))/pow(refcon,gamma);
          else
            concterm = gamma*pow(EPS13,(gamma-1.0))/pow(refcon,gamma);

          for (int vi=0; vi<nen_; ++vi)
          {
            fac_fns_i0_funct_vi = fac*fns*i0*funct_(vi)*epstort_;
            //double fac_i0_funct_vi = fac*i0*funct_(vi);
            // ------matrix: d(R_k)/d(x) = (theta*dt*(-1)*(w_k,j_k)
            for (int ui=0; ui<nen_; ++ui)
            {
              emat(vi*numdofpernode_+k,ui*numdofpernode_+k) += -fac_fns_i0_funct_vi*concterm*funct_(ui)*expterm;
              emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) += -fac_fns_i0_funct_vi*pow_conint_gamma_k*(((-alphaa)*frt*exp(alphaa*frt*eta))+((-alphac)*frt*exp((-alphac)*frt*eta)))*funct_(ui);

              if(equpot_==INPAR::ELCH::equpot_divi)
              {
                // equation for potential is scaled with inverse of Faraday constant
                emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+k) += -fac_fns_i0_funct_vi*concterm*funct_(ui)*expterm;
                emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+numscal_) += -fac_fns_i0_funct_vi*pow_conint_gamma_k*(((-alphaa)*frt*exp(alphaa*frt*eta))+((-alphac)*frt*exp((-alphac)*frt*eta)))*funct_(ui);
              }
            }
            // -----right-hand-side: -R_k = -(theta*dt*(-1)*(w_k,j_k)
            erhs[vi*numdofpernode_+k] -= -fac_fns_i0_funct_vi*pow_conint_gamma_k*expterm;

            if(equpot_==INPAR::ELCH::equpot_divi)
              erhs[vi*numdofpernode_+numscal_] -= -fac_fns_i0_funct_vi*pow_conint_gamma_k*expterm;
          }

        } // end if(kinetics=="Butler-Volmer")
        else if (kinetics==INPAR::SCATRA::butler_volmer_yang1997)
        {
          // note: gamma==0 deactivates concentration dependency in Butler-Volmer!
          double concterm = 0.0;
          if ((conint[k]/refcon) > EPS13)
            concterm = gamma*pow(conint[k],(gamma-1.0))/pow(refcon,gamma);
          else
            concterm = gamma*pow(EPS13,(gamma-1.0))/pow(refcon,gamma);

          for (int vi=0; vi<nen_; ++vi)
          {
            fac_fns_i0_funct_vi = fac*fns*i0*funct_(vi);
            // ------matrix: d(R_k)/d(x) = (theta*dt*(-1)*(w_k,j_k)
            for (int ui=0; ui<nen_; ++ui)
            {
              emat(vi*numdofpernode_+k,ui*numdofpernode_+k) += -fac_fns_i0_funct_vi*funct_(ui)*(-(concterm*exp((-alphac)*frt*eta)));
              emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) += -fac_fns_i0_funct_vi*(((-alphaa)*frt*exp(alphaa*frt*eta))+(pow_conint_gamma_k*(-alphac)*frt*exp((-alphac)*frt*eta)))*funct_(ui);
            }
            // -----right-hand-side: -R_k = -(theta*dt*(-1)*(w_k,j_k)
            erhs[vi*numdofpernode_+k] -= -fac_fns_i0_funct_vi*(exp(alphaa*frt*eta)-(pow_conint_gamma_k*exp((-alphac)*frt*eta)));
          }
        } // if (kinetics=="Butler-Volmer-Yang1997")
        else
          dserror("You should not be here!! Two ptions: Butler-Volmer-Yang1997 and Butler-Volmer-Yang1997 ");
        break;
      }
      case INPAR::SCATRA::tafel: // Tafel law (= remove anodic term from Butler-Volmer!)
      {
        // read model-specific parameter
        const double alpha = cond->GetDouble("alpha");
        double       i0 = cond->GetDouble("i0");
        i0*=timefac;
        if (i0 < -EPS14) dserror("i0 is negative, \n"
                                 "a positive definition is necessary due to generalized reaction models: %f",i0);
        const double gamma = cond->GetDouble("gamma");
        const double refcon = cond->GetDouble("refcon");
        if (refcon < EPS12) dserror("reference concentration is too small: %f",refcon);

        if(valence_k!=nume)
          dserror("Kinetic model Butler-Volmer: The number of transfered electrodes need to be  \n "
              "the same as the charge number of the reacting species %i", k);
        // opencircuit potential is assumed to be zero
        const double ocp = 0.0;
        // overpotential based on opencircuit potential
        const double eta = epd - ocp;

        // concentration-dependent Tafel law
        double pow_conint_gamma_k(0.0);

#ifdef DEBUG
        // some safety checks/ user warnings
        if (((-alpha)*frt*eta) > 100.0)
            cout<<"WARNING: Exp(alpha_c...) in Butler-Volmer law is near overflow!"
            <<exp((-alpha)*frt*eta)<<endl;
#endif
        if ((conint[k]/refcon) < EPS13)
        {
          pow_conint_gamma_k = std::pow(EPS13,gamma);
#ifdef DEBUG
            cout<<"WARNING: Rel. Conc. in Tafel formula is zero/negative: "<<(conint[k]/refcon)<<endl;
            cout<<"-> Replacement value: pow(EPS,gamma) = "<< pow_conint_gamma_k <<endl;
#endif
        }
        else
          pow_conint_gamma_k = std::pow(conint[k]/refcon,gamma);

        // note: gamma==0 deactivates concentration dependency in Butler-Volmer!
        const double expterm = -exp((-alpha)*frt*eta);

        double concterm = 0.0;
        if (conint[k] > EPS13)
          concterm = gamma*pow(conint[k],(gamma-1.0))/pow(refcon,gamma);
        else
          concterm = gamma*pow(EPS13,(gamma-1.0))/pow(refcon,gamma);

        for (int vi=0; vi<nen_; ++vi)
        {
          fac_fns_i0_funct_vi = fac*fns*i0*funct_(vi);
          // ---------------------matrix
          for (int ui=0; ui<nen_; ++ui)
          {
            emat(vi*numdofpernode_+k,ui*numdofpernode_+k) += -fac_fns_i0_funct_vi*concterm*funct_(ui)*expterm;
            emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) += -fac_fns_i0_funct_vi*pow_conint_gamma_k*(-alpha)*frt*exp((-alpha)*frt*eta)*funct_(ui);
          }
          // ------------right-hand-side
          erhs[vi*numdofpernode_+k] -= -fac_fns_i0_funct_vi*pow_conint_gamma_k*expterm;
        }
        break;
      }
      case INPAR::SCATRA::linear: // linear law:  i_n = i_0*(alphaa*frt*(V_M - phi) + 1.0)
      {
        // read model-specific parameter
        const double alphaa = cond->GetDouble("alpha_a");
        double       i0 = cond->GetDouble("i0");
        i0*=timefac;
        if (i0 < -EPS14) dserror("i0 is negative, \n"
                                 "a positive definition is necessary due to generalized reaction models: %f",i0);
        const double gamma = cond->GetDouble("gamma");
        const double refcon = cond->GetDouble("refcon");
        if (refcon < EPS12) dserror("reference concentration is too small: %f",refcon);

        if(valence_k!=nume)
          dserror("Kinetic model Butler-Volmer: The number of transfered electrodes need to be  \n "
              "the same as the charge number of the reacting species %i", k);
        // opencircuit potential is assumed to be zero
        const double ocp = 0.0;
        // overpotential based on opencircuit potential
        const double eta = epd - ocp;

        double pow_conint_gamma_k = 0.0;
        if ((conint[k]/refcon) < EPS13)
        {
          pow_conint_gamma_k = std::pow(EPS13,gamma);
#ifdef DEBUG
          cout<<"WARNING: Rel. Conc. in Tafel formula is zero/negative: "<<(conint[k]/refcon)<<endl;
          cout<<"-> Replacement value: pow(EPS,gamma) = "<< pow_conint_gamma_k <<endl;
#endif
        }
        else
          pow_conint_gamma_k = std::pow(conint[k]/refcon,gamma);
        const double linearfunct = (alphaa*frt*eta + 1.0);
        // note: gamma==0 deactivates concentration dependency
        double concterm = 0.0;
        if (conint[k] > EPS13)
          concterm = gamma*pow(conint[k],(gamma-1.0))/pow(refcon,gamma);
        else
          dserror("Better stop here!");

        for (int vi=0; vi<nen_; ++vi)
        {
          fac_fns_i0_funct_vi = fac*fns*i0*funct_(vi);
          const int fvi = vi*numdofpernode_+k;
          // ---------------------matrix
          for (int ui=0; ui<nen_; ++ui)
          {
            emat(fvi,ui*numdofpernode_+k) += -fac_fns_i0_funct_vi*concterm*funct_(ui)*linearfunct;
            emat(fvi,ui*numdofpernode_+numscal_) += -fac_fns_i0_funct_vi*pow_conint_gamma_k*(-alphaa)*frt*funct_(ui);
          }
          // ------------right-hand-side
          erhs[fvi] -= -fac_fns_i0_funct_vi*pow_conint_gamma_k*linearfunct;
        }
        break;
      }
      case INPAR::SCATRA::butler_volmer_newman:
      {
        // "Electrochemical systems"
        // Newman and Thomas-Alyea, 2004
        // General stoichiometry: pp. 212-213, e.q. 8.26
        // consideration of a elementary step of the form:
        // Sum_i s_i M_i ->  ne-
        // n is one if charge transfer is involved, multiple electron transfers "being unlikely in
        // an elementary step

        const double k_a = cond->GetDouble("k_a");
        const double k_c = cond->GetDouble("k_c");
        const double beta = cond->GetDouble("beta");

        //reaction order of the cathodic and anodic reactants of ionic species k
        std::vector<int> q(numscal_,0);
        std::vector<int> p(numscal_,0);

        for(int ii=0; ii<numscal_; ii++)
        {
          //according to the convention: anodic reactant is positiv
          if(stoich[ii] > 0)
          {
            q[ii] = 0;
            p[ii] = stoich[ii];
          }
          //according to the convention: cathodic reactant is negative
          else
          {
            q[ii]= -stoich[ii];
            p[ii] = 0;
          }
        }

#ifdef DEBUG
        // some safety checks/ user warnings
        if (((1-beta)*frt*epd) > 100.0)
          cout<<"WARNING: Exp((1-beta)...) in Butler-Volmer law is near overflow!"
          <<exp((1-beta)*frt*epd)<<endl;
        if (((-beta)*frt*epd) > 100.0)
          cout<<"WARNING: Exp(-beta...) in Butler-Volmer law is near overflow!"
          <<exp((-beta)*frt*epd)<<endl;
#endif

        double pow_conint_p(1.0);      //product over i (c_i)^(p_i)
        double pow_conint_q(1.0);      //product over i (c_i)^(q_i)
        std::vector<double> pow_conint_p_derivative(numscal_,1.0);  //pow_conint_p derivated after conint[nspec]
        std::vector<double> pow_conint_q_derivative(numscal_,1.0); //pow_conint_q derivated after conint[nspec]

        //concentration term (product of cathodic and anodic species)
        for(int kk=0; kk<numscal_; ++kk)
        {
          if ((conint[kk]) < EPS13) // 1.0E-16)
          {
            pow_conint_p *= std::pow(EPS13,p[kk]);
            pow_conint_q *= std::pow(EPS13,q[kk]);
#ifdef DEBUG
            cout<<"WARNING: Rel. Conc. of species" <<k<<" in Butler-Volmer formula is zero/negative: "<<(conint[k])<<endl;
            cout<<"-> Replacement value: pow(1.0E-16,p[ispec]) = "<< pow(EPS13,p[k]) << " pow(1.0E-13,q[k]) = "<< pow(EPS13,q[k]) <<endl;
#endif
          }
          else
          {
            pow_conint_p *= std::pow((conint[kk]),p[kk]);
            pow_conint_q *= std::pow((conint[kk]),q[kk]);
          }
        }

        //derivation of concentration term  with respect to ionic species kk
        for(int kk=0; kk<numscal_; ++kk)
        {
          pow_conint_p_derivative[kk] = pow_conint_p*p[kk]/conint[kk];
          pow_conint_q_derivative[kk] = pow_conint_q*q[kk]/conint[kk];
        }

        // loop over reacting species; determines the line of the matrix
        const double expterma = exp((1-beta)*nume*frt*epd);
        const double exptermc = exp((-beta)*nume*frt*epd);

        for (int vi=0; vi<nen_; ++vi)
        {
          // see Wittmann, Erweiterte Reaktionsmodelle für die numerische Simulation von
          // elektrochemischen Systemen, p.20, equ. 3.4
          const double fac_fns_funct_vi = faraday*nume*fac*fns*funct_(vi);
          for (int ui=0; ui<nen_; ++ui)
          {
            //loop over the columns of the matrix, makes sure that the linearisation w.r.t the first concentration is added to the first column
            for(int kk=0; kk<numscal_; ++kk)
            {
              emat(vi*numdofpernode_+k,ui*numdofpernode_+kk) += -fac_fns_funct_vi*((k_a*expterma*pow_conint_p_derivative[kk]) - (k_c*exptermc*pow_conint_q_derivative[kk]))*funct_(ui)*timefac;
            }
            //linearisation w.r.t the potential
            emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) += -fac_fns_funct_vi*((-k_a*(1-beta)*nume*frt*expterma*pow_conint_p) - (k_c*beta*nume*frt*exptermc*pow_conint_q))*funct_(ui)*timefac;
          }
          // ------------right-hand-side
          erhs[vi*numdofpernode_+k] -= -(fac_fns_funct_vi*((k_a*expterma*pow_conint_p)-(k_c*exptermc*pow_conint_q)))*timefac;
        }
        break;
      }
      case INPAR::SCATRA::butler_volmer_bard:
      {
        // "Electrochemcial Methods Fundamentals and Applications"
        // Bard and Faulkner, 2001, pp. 94 ff; pp. 99 eq. 3.4.10
        // reaction model for a one-step, one-electron process (elementar step)
        // O + e -> R
        const double e0 = cond->GetDouble("e0");
        const double k0 = cond->GetDouble("k0");
        const double beta = cond->GetDouble("beta");
        const double c_c0 = cond->GetDouble("c_c0");
        const double c_a0 = cond->GetDouble("c_a0");

        if(nume!=1)
          dserror("electron != 1; \n "
              "this Butler-Volmer-equation (Bard/Faulkner) works for elementary steps (one electron) only!");

        // only one reactant and product are supported by the basic model
        // only stoichometry of 1
        {
          int check1 = 0;
          int check2 = 0;
          for(int kk=0; kk<numscal_;kk++)
          {
            if(abs(stoich[kk])>1)
              dserror("Stoichometry is larger than 1!! \n"
                      "This is not supported by the reaction model based on Bard");

            check1 += abs(stoich[kk]);
            check2 += stoich[kk];
          }
          if (check1>2 or check1==0)
            dserror("More than one reactant or product defined!! \n"
                    "This is not supported by the reaction model based on Bard");

          // In the moment it is not checked if two products or reactants are defined
        }

        // equilibrium potential (equpot):
        // defined in Bard, 2001, p.98, eq. 3.4.3
        const double equpot = e0 + (log(c_c0/c_a0))/(frt*nume);
        // overpotential based on equilibrium potential
        const double eta_equpot = epd - equpot;

        // negative sign: we look at electon flow
        const double i0 = k0*pow(c_c0,1-beta)*pow(c_a0,beta)*nume*faraday;

        // reactant or product not a species in the electrolyte
        // -> concentration = 1.0
        double conctermc = 1.0;
        double concterma = 1.0;
        double conctermc_der = 1.0;
        double concterma_der = 1.0;
        //species id of the anodic and cathodic reactant
        int anodic = 0;
        int cathodic = 0;
        bool checkc = 0;
        bool checka = 0;

        // concentration terms for anodic and cathodic reaction
        // only one reactant and product are supported by the basic model
        // only stoichometry of 1
        for(int kk=0; kk<numscal_;kk++)
        {
          if(stoich[kk]==1)
          {
            concterma = conint[kk]/c_a0;
            concterma_der = 1.0/c_a0;
            anodic = kk;
            checka = true;
          }
          else if(stoich[kk]==-1)
          {
            conctermc = conint[kk]/c_c0;
            conctermc_der = 1.0/c_c0;
            cathodic = kk;
            checkc = true;
          }
        }

#ifdef DEBUG
        // some safety checks/ user warnings
        if (((1-beta)*(frt*nume)*eta_equpot) > 100.0)
          cout<<"WARNING: Exp((1-beta)...) in Butler-Volmer law is near overflow!"
          <<exp((1-beta)*(frt*nume)*eta_equpot)<<endl;
        if (((-beta)*(frt*nume)*eta_equpot) > 100.0)
          cout<<"WARNING: Exp(-beta...) in Butler-Volmer law is near overflow!"
          <<exp((-beta)*(frt*nume)*eta_equpot)<<endl;
#endif

        const double expterma = exp((1-beta) * (frt*nume) * eta_equpot);
        const double exptermc = exp((-beta) * (frt*nume) * eta_equpot);

        for (int vi=0; vi<nen_; ++vi)
        {
          const double fac_i0_funct_vi = fac*fns*i0*funct_(vi);
          // ---------------------matrix
          for (int ui=0; ui<nen_; ++ui)
          {
            //derivation wrt concentration
            if(checkc == true)
              emat(vi*numdofpernode_+k,ui*numdofpernode_+cathodic) += fac_i0_funct_vi*conctermc_der*exptermc*funct_(ui)*timefac;
            if(checka == true)
              emat(vi*numdofpernode_+k,ui*numdofpernode_+anodic)   += -fac_i0_funct_vi*concterma_der*expterma*funct_(ui)*timefac;
            //derivation wrt potential
            emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_)
              += -fac_i0_funct_vi*(-concterma*(1-beta)*frt*nume*expterma-conctermc*beta*nume*frt*exptermc)*funct_(ui)*timefac;
          }
          // ------------right-hand-side
          erhs[vi*numdofpernode_+k] -= -fac_i0_funct_vi*(concterma*expterma - conctermc*exptermc)*timefac;
        }
        break;
      }
      case INPAR::SCATRA::zero:
        break;
      case INPAR::SCATRA::nernst:
      {
#if 0
        for (int vi=0; vi<nen_; ++vi)
        {
          for (int ui=0; ui<nen_; ++ui)
          {
            emat(vi*numdofpernode_+k,ui*numdofpernode_+k) += -fac*funct_(vi)/conint[k]*funct_(ui);
          }
          erhs[vi*numdofpernode_+k] -= -fac*funct_(vi)/frt/conint[k];
        }
#endif
        break;
      }
      default:
        dserror("Kinetic model not implemented");
        break;
    }
    }
    else // secondary current distribution
    {
      switch(kinetics)
      {
      case INPAR::SCATRA::butler_volmer:  // concentration-dependent Butler-Volmer law
      {
        // read model-specific parameter
        const double alphaa = cond->GetDouble("alpha_a");
        const double alphac = cond->GetDouble("alpha_c");
        double       i0 = cond->GetDouble("i0");
        i0*=timefac;
        if (i0 < -EPS14) dserror("i0 is negative, \n"
                                 "a positive definition is necessary due to generalized reaction models: %f",i0);
        const double refcon = cond->GetDouble("refcon");
        if (refcon < EPS12) dserror("reference concentration is too small: %f",refcon);
        // opencircuit potential is assumed to be zero
        const double ocp = 0.0;
        // overpotential based on opencircuit potential
        const double eta = epd - ocp;

#ifdef DEBUG
        // some safety checks/ user warnings
        if ((alphaa*eta) > 100.0)
          cout<<"WARNING: Exp(alpha_a...) in Butler-Volmer law is near overflow!"
          <<exp(alphaa*eta)<<endl;
        if (((-alphac)*eta) > 100.0)
          cout<<"WARNING: Exp(alpha_c...) in Butler-Volmer law is near overflow!"
          <<exp((-alphac)*eta)<<endl;
#endif
        // Butler-Volmer kinetics
        const double expterm = exp(alphaa*eta)-exp((-alphac)*eta);
        const double exptermderiv = (((-alphaa)*exp(alphaa*eta))+((-alphac)*exp((-alphac)*eta)));

        for (int vi=0; vi<nen_; ++vi)
        {
          const double fac_i0_funct_vi = fac*i0*funct_(vi);
          // ---------------------matrix
          for (int ui=0; ui<nen_; ++ui)
          {
            emat(vi*numdofpernode_,ui*numdofpernode_) += -fac_i0_funct_vi*exptermderiv*funct_(ui);
          }
          // ------------right-hand-side
          erhs[vi*numdofpernode_] -= -fac_i0_funct_vi*expterm;
        }
        break;
      }
      case INPAR::SCATRA::tafel: // Tafel kinetics
      {
        // read model-specific parameter
        const double alpha = cond->GetDouble("alpha");
        double       i0 = cond->GetDouble("i0");
        i0*=timefac;
        if (i0 < -EPS14) dserror("i0 is negative, \n"
                                 "a positive definition is necessary due to generalized reaction models: %f",i0);
        const double refcon = cond->GetDouble("refcon");
        if (refcon < EPS12) dserror("reference concentration is too small: %f",refcon);
        // opencircuit potential is assumed to be zero
        const double ocp = 0.0;
        // overpotential based on opencircuit potential
        const double eta = epd - ocp;

#ifdef DEBUG
        // some safety checks/ user warnings
        if ((-alpha)*eta > 100.0)
          cout<<"WARNING: Exp(alpha_c...) in Tafel law is near overflow!"
          <<exp((-alpha)*eta)<<endl;
#endif
        const double expterm = -exp((-alpha)*eta);
        const double exptermderiv = alpha*expterm; // do not forget the (-1) from differentiation of eta!

        for (int vi=0; vi<nen_; ++vi)
        {
          const double fac_i0_funct_vi = fac*i0*funct_(vi);
          // ---------------------matrix
          for (int ui=0; ui<nen_; ++ui)
          {
            emat(vi*numdofpernode_,ui*numdofpernode_) += -fac_i0_funct_vi*exptermderiv*funct_(ui);
          }
          // ------------right-hand-side
          erhs[vi*numdofpernode_] -= -fac_i0_funct_vi*expterm;
        }
        break;
      }
      case INPAR::SCATRA::linear: // linear law:  i_n = i_0*(alphaa*(V_M - phi) + 1.0))
      {
        // read model-specific parameter
        const double alphaa = cond->GetDouble("alpha_a");
        double       i0 = cond->GetDouble("i0");
        i0*=timefac;
        if (i0 < -EPS14) dserror("i0 is negative, \n"
                                 "a positive definition is necessary due to generalized reaction models: %f",i0);
        const double refcon = cond->GetDouble("refcon");
        if (refcon < EPS12) dserror("reference concentration is too small: %f",refcon);
        // opencircuit potential is assumed to be zero
        const double ocp = 0.0;
        // overpotential based on opencircuit potential
        const double eta = epd - ocp;

        for (int vi=0; vi<nen_; ++vi)
        {
          const double fac_i0_funct_vi = fac*i0*funct_(vi);
          // ---------------------matrix
          for (int ui=0; ui<nen_; ++ui)
          {
            emat(vi*numdofpernode_,ui*numdofpernode_) += -fac_i0_funct_vi*(-alphaa)*funct_(ui);
          }
          // ------------right-hand-side
          erhs[vi*numdofpernode_] -= -fac_i0_funct_vi*((alphaa*eta)+1.0);
        }
        break;
      }
      case INPAR::SCATRA::butler_volmer_yang1997:
      case INPAR::SCATRA::butler_volmer_newman:
      case INPAR::SCATRA::butler_volmer_bard:
      case INPAR::SCATRA::zero:
      {
        dserror("Kinetic model not implemented: %i",kinetics);
        break;
      }
      default:
        dserror("How did you come here? All kinetic models have been already addressed!");
        break;
      }
        //dserror("Kinetic model not implemented: %s",kinetics.c_str());
    } // if iselch
  } // end of loop over integration points gpid

  // TODO (ehrl): Is it ok to use nume instead of valence?
  if (iselch)
  {
    if ((scatratype==INPAR::SCATRA::scatratype_elch_enc_pde) or
        (scatratype==INPAR::SCATRA::scatratype_elch_enc_pde_elim))
    {
      // we have to add boundary contributions to the potential equation as well!
      // and do not forget the corresponding matrix contributions ;-)
      double val(0.0);
      for (int vi=0; vi<nen_; ++vi)
      {
        // ---------------------matrix
        for (int ui=0; ui<nen_; ++ui)
        {
          val = emat(vi*numdofpernode_+k,ui*numdofpernode_+k);
          emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+k) += nume*val;
          val = emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_);
          emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+numscal_) += nume*val;
        }
        // ------------right-hand-side
        val = erhs[vi*numdofpernode_+k];
        erhs[vi*numdofpernode_+numscal_] += nume*val;
      }
    }
  } // end if(iselch): adaptation of kinetics for eliminated scalar transport equation
  } // end for(int k=0; k<numcal;++k) loop over scalars
  return;
} // ScaTraBoundaryImpl<distype>::EvaluateElectrodeKinetics()


/*----------------------------------------------------------------------*
 | calculate electrode kinetics status information             gjb 01/09|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::ElectrodeStatus(
    const DRT::Element*        ele,
    Teuchos::ParameterList& params,
    Teuchos::RCP<DRT::Condition>  cond,
    const std::vector<double>&   ephinp,
    const std::vector<double>& ephidtnp,
    const int             kinetics,
    const std::vector<int>  stoich,
    const int                 nume,
    const double              pot0,
    const double               frt,
    const bool              iselch,
    const double           timefac
)
{
  double faraday = INPAR::SCATRA::faraday_const;    // unit of F: C/mol or mC/mmol or muC / mumol

  // get variables with their current values
  double currentintegral   = params.get<double>("currentintegral");
  double boundaryint       = params.get<double>("boundaryintegral");
  double overpotentialint  = params.get<double>("overpotentialintegral");
  double electdiffpotint   = params.get<double>("electrodedifferencepotentialintegral");
  double opencircuitpotint = params.get<double>("opencircuitpotentialintegral");
  double concentrationint  = params.get<double>("concentrationintegral");
  double currderiv         = params.get<double>("currentderiv");
  double currentresidual   = params.get<double>("currentresidual");

  double dlcapacitance = 0.0;
  if(kinetics != INPAR::SCATRA::zero)
  {
    dlcapacitance = cond->GetDouble("dlcap");
    if (dlcapacitance < -EPS12) dserror("double-layer capacitance is negative: %f",dlcapacitance);
  }

  double pot0hist(0.0);
  // access history of electrode potential
  if (dlcapacitance > EPS12)
  {
    pot0hist = cond->GetDouble("pothist");
  }

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // concentration values of reactive species at element nodes
  std::vector<LINALG::Matrix<nen_,1> > conreact(numscal_);

  // el. potential values at element nodes
  LINALG::Matrix<nen_,1> pot(true);
  LINALG::Matrix<nen_,1> potdtnp(true);
  if(iselch)
  {
    for (int inode=0; inode< nen_;++inode)
    {
      for (int kk=0; kk<numscal_; ++kk)
        conreact[kk](inode) = ephinp[inode*numdofpernode_+kk];

      pot(inode) = ephinp[inode*numdofpernode_+numscal_];
      potdtnp(inode) = ephidtnp[inode*numdofpernode_+numscal_];
    }
  }
  else
  {
    for (int inode=0; inode< nen_;++inode)
    {
      for (int kk=0; kk<numscal_; ++kk)
        conreact[kk](inode) = 1.0;

      pot(inode) = ephinp[inode*numdofpernode_];
      potdtnp(inode) = ephidtnp[inode*numdofpernode_];
    }
  }

  // concentration of active species at integration point
  std::vector<double> conint(numscal_,0.0);
  // el. potential at integration point
  double potint;
  // history term of el. potential at integration point
  double potdtnpint;

  bool statistics = false;
  // index of reactive species (starting from zero)
  for(int k=0; k<numscal_; ++k)
  {
    // only the first oxidized species O is considered for statistics
    // statistics of other species result directly from the oxidized species (current density, ...)
    // or need to be implemented (surface concentration, OCV, ...)
    if(stoich[k]>=0)
      continue;

    statistics = true;

    // loop over integration points
    for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
    {
      const double fac = EvalShapeFuncAndIntFac(intpoints,gpid,ele->Id());

      // elch-specific values at integration point:
      for (int kk=0; kk<numscal_; ++kk)
        conint[kk] = funct_.Dot(conreact[kk]);

      potint = funct_.Dot(pot);
      potdtnpint = funct_.Dot(potdtnp);

      // surface overpotential eta at integration point
      const double epd = (pot0 - potint);

      // linearization of current w.r.t applied electrode potential "pot0"
      double linea(0.0);

      // concentration-dependent Butler-Volmer law(s)
      switch(kinetics)
      {
      case INPAR::SCATRA::butler_volmer:
      case INPAR::SCATRA::butler_volmer_yang1997:
      {
        // read model-specific parameter
        const double alphaa = cond->GetDouble("alpha_a");
        const double alphac = cond->GetDouble("alpha_c");
        double       i0 = cond->GetDouble("i0");
        if (i0 < -EPS14) dserror("i0 is negative, \n"
                                 "a positive definition is necessary due to generalized reaction models: %f",i0);
        const double gamma = cond->GetDouble("gamma");
        const double refcon = cond->GetDouble("refcon");
        if (refcon < EPS12) dserror("reference concentration is too small: %f",refcon);

        // opencircuit potential is assumed to be zero
        const double ocp = 0.0;
        // overpotential based on opencircuit potential
        const double eta = epd - ocp;

        double expterm(0.0);
        if (iselch)
        {
          if (kinetics==INPAR::SCATRA::butler_volmer)
          {
            // general Butler-Volmer
            expterm = std::pow(conint[k]/refcon,gamma) * (exp(alphaa*frt*eta)-exp((-alphac)*frt*eta));
            linea = std::pow(conint[k]/refcon,gamma) * frt*((alphaa*exp(alphaa*frt*eta)) + (alphac*exp((-alphac)*frt*eta)));
          }
          if (kinetics==INPAR::SCATRA::butler_volmer_yang1997)
          {
            if (((conint[k]/refcon)<EPS13) && (gamma < 1.0))
            {// prevents NaN's in the current density evaluation
              expterm = (exp(alphaa*frt*eta)-(pow(EPS13/refcon,gamma)*exp((-alphac)*frt*eta)));
              linea = ((alphaa)*frt*exp(alphaa*frt*eta))+(pow(EPS13/refcon,gamma)*alphac*frt*exp((-alphac)*frt*eta));
            }
            else
            {
              expterm = (exp(alphaa*frt*eta)-(pow(conint[k]/refcon,gamma)*exp((-alphac)*frt*eta)));
              linea = ((alphaa)*frt*exp(alphaa*frt*eta))+(pow(conint[k]/refcon,gamma)*alphac*frt*exp((-alphac)*frt*eta));
            }
          }
       }
       else // secondary current distribution
       {
         expterm = exp(alphaa*eta)-exp((-alphac)*eta);
         linea = (alphaa*exp(alphaa*eta))+(alphac*exp((-alphac)*eta));
       }

        // scan for NaNs due to negative concentrations under exponent gamma
        if (std::isnan(expterm) or std::isnan(linea))
          dserror("NaN detected in electrode status calculation");

        // compute integrals
        overpotentialint += eta*fac;
        electdiffpotint += epd*fac;
        opencircuitpotint += ocp*fac;
        currentintegral += i0*expterm*fac; // the negative(!) normal flux density
        boundaryint += fac;
        concentrationint += conint[k]*fac;

        // tangent and rhs (= negative residual) for galvanostatic equation
        currderiv += i0*linea*timefac*fac;
        currentresidual += i0 * expterm * timefac *fac;

        if (dlcapacitance > EPS12)
        {
          // add contributions due to double-layer capacitance
          // positive due to redefinition of the exchange current density
          currderiv += fac*dlcapacitance;
          currentresidual += fac*dlcapacitance*(pot0-pot0hist-(timefac*potdtnpint));
#if 0
          cout<<"pot0-pot0hist = "<<pot0 << " - "<<pot0hist<<endl;
          cout<<"- timefac*potdtnpint = -"<<timefac<<" * "<<potdtnpint<<endl;
#endif
        }
        break;
      }
      case INPAR::SCATRA::tafel: // concentration-dependent Tafel kinetics
      {
        // read model-specific parameter
        const double alpha = cond->GetDouble("alpha");
        double       i0 = cond->GetDouble("i0");
        if (i0 < -EPS14) dserror("i0 is negative, \n"
                                 "a positive definition is necessary due to generalized reaction models: %f",i0);
        const double gamma = cond->GetDouble("gamma");
        const double refcon = cond->GetDouble("refcon");
        if (refcon < EPS12) dserror("reference concentration is too small: %f",refcon);

        // opencircuit potential is assumed to be zero
        const double ocp = 0.0;
        // overpotential based on opencircuit potential
        const double eta = epd - ocp;

        if(iselch)
        {
          const double expterm = std::pow(conint[k]/refcon,gamma) * (-exp((-alpha)*frt*eta));
          linea = std::pow(conint[k]/refcon,gamma) * frt*(alpha*exp((-alpha)*frt*eta));
          // compute integrals
          overpotentialint += eta * fac;
          electdiffpotint += epd*fac;
          opencircuitpotint += ocp*fac;
          currentintegral += i0 * expterm * fac; // the negative(!) normal flux density
          boundaryint += fac;
          concentrationint += conint[k]*fac;

          // tangent and rhs (= negative residual) for galvanostatic equation
          currderiv += i0*linea*timefac*fac;
          currentresidual += i0*expterm*timefac*fac;
        }
        else
        {
          // secondary current distribution with Tafel kinetics
          double expterm = -exp((-alpha)*eta);
          //linea = (-alphac)*exp((-alphac)*eta);

          // compute integrals
          overpotentialint += eta * fac;
          electdiffpotint += epd*fac;
          opencircuitpotint += ocp*fac;
          currentintegral += i0 * expterm * fac; // the negative(!) normal flux density
          boundaryint += fac;
        }
        break;
      }
      case INPAR::SCATRA::linear: // linear: i_n = i_0*(alphaa*frt*eta + 1.0)
      {
        // read model-specific parameter
        const double alphaa = cond->GetDouble("alpha_a");
        double       i0 = cond->GetDouble("i0");
        if (i0 < -EPS14) dserror("i0 is negative, \n"
                                 "a positive definition is necessary due to generalized reaction models: %f",i0);
        const double gamma = cond->GetDouble("gamma");
        const double refcon = cond->GetDouble("refcon");
        if (refcon < EPS12) dserror("reference concentration is too small: %f",refcon);

        // opencircuit potential is assumed to be zero
        const double ocp = 0.0;
        // overpotential based on opencircuit potential
        const double eta = epd - ocp;

        if(iselch)
        {
          // compute integrals
          overpotentialint += eta * fac;
          electdiffpotint += epd*fac;
          opencircuitpotint += ocp*fac;
          currentintegral += i0 * pow(conint[k]/refcon,gamma)*(alphaa*frt*eta + 1.0) * fac; // the negative(!) normal flux density
          boundaryint += fac;
          concentrationint += conint[k]*fac;

          // tangent and rhs (= negative residual) for galvanostatic equation
          linea = std::pow(conint[k]/refcon,gamma)*(alphaa*frt);
          currderiv += i0*linea*timefac*fac;
          currentresidual += i0*(alphaa*frt*eta + 1.0)*timefac*fac;
        }
        else
        {
          // secondary current distribution with linear kinetics
          // compute integrals
          overpotentialint += eta * fac;
          electdiffpotint += epd*fac;
          opencircuitpotint += ocp*fac;
          currentintegral += i0 *((alphaa*eta) +1.0)* fac; // the negative(!) normal flux density
          boundaryint += fac;
        }
        break;
      }
      case INPAR::SCATRA::butler_volmer_newman:
      {
        // "Electrochemical systems"
        // Newman ad Thomas-Alyea, 2004
        // General stoichiometry: pp. 212-213, e.q. 8.26
        // consideration of a elementary step of the form:
        // Sum_i s_i M_i ->  ne-
        // n is one if charge transfer is involved, multiple electron transfers "being unlikely in
        // an elementary step

        const double k_a = cond->GetDouble("k_a");
        const double k_c = cond->GetDouble("k_c");
        const double beta = cond->GetDouble("beta");

        //reaction order of the cathodic and anodic reactants of ionic species k
        std::vector<int> q(numscal_,0);
        std::vector<int> p(numscal_,0);

        for(int kk=0; kk<numscal_; kk++)
        {
          //according to the convention: anodic reactant is positiv
          if(stoich[kk] > 0)
          {
            q[kk] = 0;
            p[kk] = stoich[kk];
          }
          //according to the convention: cathodic reactant is negative
          else
          {
            q[kk]= -stoich[kk];
            p[kk] = 0;
          }
        }

        // linearization of current w.r.t applied electrode potential "pot0"
        double linea(0.0);
        double expterma(0.0);
        double exptermc(0.0);
        double expterm(0.0);
        double pow_conint_p = 1.0;      //product over i (c_i)^(p_i)
        double pow_conint_q = 1.0;      //product over i (c_i)^(q_i)

        if (iselch)
        {
          for(int kk=0; kk<numscal_; ++kk)
          {
            if ((conint[kk]) < EPS13)
            {
              pow_conint_p *= std::pow(EPS13,p[kk]);
              pow_conint_q *= std::pow(EPS13,q[kk]);
#ifdef DEBUG
              cout<<"WARNING: Rel. Conc. of species"<<kk<<" in Butler-Volmer formula is zero/negative: "<<(conint[kk])<<endl;
              cout<<"-> Replacement value: pow(EPS,p[ispec]) = "<< pow(EPS13,p[kk]) << " pow(1.0E-16,q[i]) = "<< pow(EPS13,q[kk]) <<endl;
#endif
            }
            else
            {
              pow_conint_p *= std::pow((conint[kk]),p[kk]);
              pow_conint_q *= std::pow((conint[kk]),q[kk]);
            }
          }
          expterma = exp((1-beta)*nume*frt*epd);
          exptermc = exp(-beta*nume*frt*epd);
          linea =  nume*faraday*(frt*nume*((k_a*(1-beta)*expterma*pow_conint_p)-(k_c*(-1)*beta*exptermc*pow_conint_q)));
        }
        else // secondary current distribution
          dserror("not implemented in the new function");

        // scan for NaNs due to negative concentrations under exponent gamma
        if (std::isnan(expterm) or std::isnan(linea))
          dserror("NaN detected in electrode status calculation");

        // open circuit potential (ocp): time dependent electrode surface concentrations
        // defined in Newman, 2004, p. 211, eq. 8.20
        double ocp = 1/frt/nume*log(k_c/k_a);
        for(int kk=0;kk<numscal_;++kk)
        {
          ocp +=  1/frt/nume*(q[kk]-p[kk])*log(conint[kk]);
          //safety check
          if((q[kk]-p[kk])!=-stoich[kk])
            dserror("stoichiometry factors and the factors q,p do not correlate!!");
        }
        // overpotential based on open circuit potential
        const double eta = epd - ocp;

        currentintegral += nume*faraday*((k_a*expterma*pow_conint_p)-(k_c*exptermc*pow_conint_q))*fac;

        boundaryint += fac;
        overpotentialint += eta * fac;
        electdiffpotint += epd*fac;
        opencircuitpotint += ocp*fac;
        concentrationint += conint[k]*fac;

        // tangent and rhs (= negative residual) for galvanostatic equation
        currderiv += linea*fac*timefac;
        currentresidual += nume*faraday*((k_a*expterma*pow_conint_p)-(k_c*exptermc*pow_conint_q))*timefac*fac;

        if (dlcapacitance > EPS13)
        {
          dserror("Capacity of the double layer are not tested, but implemented");
          // add contributions due to double-layer capacitance
          // positive due to redefinition of the exchange current density
          currderiv += fac*dlcapacitance;
          currentresidual += fac*dlcapacitance*(pot0-pot0hist-(timefac*potdtnpint));
#if 0
          cout<<"pot0-pot0hist = "<<pot0 << " - "<<pot0hist<<endl;
          cout<<"- timefac*potdtnpint = -"<<timefac<<" * "<<potdtnpint<<endl;
#endif
        }
        break;
      }
      case INPAR::SCATRA::butler_volmer_bard:
      {
        // "Electrochemcial Methods Fundamentals and Applications"
        // Bard and Faulkner, 2001, pp. 94 ff; pp. 99 eq. 3.4.10
        // reaction model for a one-step, one-electron process
        // O + e -> R
        const double e0 = cond->GetDouble("e0");
        const double k0 = cond->GetDouble("k0");
        const double beta = cond->GetDouble("beta");
        const double c_c0 = cond->GetDouble("c_c0");
        const double c_a0 = cond->GetDouble("c_a0");

        if(nume!=1)
          dserror("electron != 1; \n "
              "this Butler-Volmer-equation (Bard/Faulkner) works for elementary steps (one electron) only!");

        // only one reactant and product are supported by the basic model
        // only stoichometry of 1
        {
          int check1 = 0;
          int check2 = 0;
          for(int kk=0; kk<numscal_;kk++)
          {
            if(abs(stoich[kk])>1)
              dserror("Stoichometry is larger than 1!! \n"
                      "This is not supported by the reaction model based on Bard");

            check1 += abs(stoich[kk]);
            check2 += stoich[kk];
          }
          if (check1>2 or check1==0)
            dserror("More than one reactant or product defined!! \n"
                "This is not supported by the reaction model based on Bard");

          // In the moment it is not checked if two products (and no reactants) and vis versa are defined
        }

        // reactant or product not a species in the electrolyte
        // -> concentration = 1.0
        double conctermc = 1.0;
        double concterma = 1.0;

        // concentration terms for anodic and cathodic reaction
        // only one reactant and product are supported by the basic model
        // only stoichometry of 1
        for(int kk=0; kk<numscal_;kk++)
        {
          if(stoich[kk]==1)
            concterma = conint[kk]/c_a0;
          else if(stoich[kk]==-1)
            conctermc = conint[kk]/c_c0;
        }

        // equilibrium potential (equpot):
        // defined in Bard, 2001, p.98, eq. 3.4.3
        const double equpot = e0 + (log(c_c0/c_a0))/(frt*nume);
        // overpotential based on equilibrium potential
        const double eta_equpot = epd - equpot;
        // difference between equilibrium potential and open circuit potential:
        // -> equpot: depends on initial electrode surface concentration
        // -> ocp:    depends on actual electrode surface concentration

        // open circuit potential (ocp): time dependent electrode surface concentrations
        // defined in Newman, 2004, p. 211, eq. 8.20
        const double ocp = e0 + 1/frt/nume*log(conctermc/concterma);
        // overpotential based on open circuit potential
        const double eta = epd - ocp;

        const double expterma = exp((1-beta)*(frt*nume)*eta_equpot);
        const double exptermc = exp(-beta*(frt*nume)*eta_equpot);
        const double linea = concterma*(1-beta)*(frt*nume)*expterma+conctermc*beta*(frt*nume)*exptermc;

        if (not iselch)
          dserror("iselch not true; not implemented");

        // scan for NaNs due to negative concentrations under exponent gamma
        if (std::isnan(expterma) or std::isnan(exptermc) or std::isnan(linea))
          dserror("NaN detected in electrode status calculation");

        const double i0 = faraday*k0*pow(c_c0,1-beta)*pow(c_a0,beta);

        // compute integrals
        overpotentialint += eta*fac;
        electdiffpotint += epd*fac;
        opencircuitpotint += ocp*fac;
        currentintegral += i0*(concterma*expterma-conctermc*exptermc)*fac; // the negative(!) normal flux density
        boundaryint += fac;
        concentrationint += conint[k]*fac;  //concentration-output for the first species only

        // tangent and rhs (= negative residual) for galvanostatic equation
        currderiv += i0*linea*timefac*fac;
        currentresidual += i0*(concterma*expterma-conctermc*exptermc)*timefac*fac;

        if (dlcapacitance > EPS12)
        {
          dserror("Capacity of the double layer are not tested, but implemented");
          // add contributions due to double-layer capacitance
          // positive due to redefinition of the exchange current density
          currderiv += fac*dlcapacitance;
          currentresidual += fac*dlcapacitance*(pot0-pot0hist-(timefac*potdtnpint));
#if 0
          cout<<"pot0-pot0hist = "<<pot0 << " - "<<pot0hist<<endl;
          cout<<"- timefac*potdtnpint = -"<<timefac<<" * "<<potdtnpint<<endl;
#endif
        }
        break;
      } //end Butler-Volmer-Bard
      case INPAR::SCATRA::zero: // zero
      {
        // opencircuit potential is assumed to be zero
        const double ocp = 0.0;
        // overpotential based on opencircuit potential
        const double eta = epd - ocp;

        // compute integrals
        overpotentialint -= eta * fac;
        boundaryint += fac;
        concentrationint += conint[k]*fac;
        break;
      }
      case INPAR::SCATRA::nernst: // zero
      {
        // compute integrals
        overpotentialint += potint * fac;
        boundaryint += fac;
        concentrationint += conint[k]*fac;
        break;
      }
      default:
        dserror("Kinetic model not implemented");
        break;
      }
    }  // loop over integration points
    //stop loop over ionic species after one evaluation
    break;
  }  // loop over scalars

  if(statistics == false)
    dserror("There is no oxidized species O (stoich<0) defined in your input file!! \n"
            " Statistics could not be evaluated");

  // add contributions to the global values
  params.set<double>("currentintegral",currentintegral);
  params.set<double>("boundaryintegral",boundaryint);
  params.set<double>("overpotentialintegral",overpotentialint);
  params.set<double>("electrodedifferencepotentialintegral",electdiffpotint);
  params.set<double>("opencircuitpotentialintegral",opencircuitpotint);
  params.set<double>("concentrationintegral",concentrationint);
  params.set<double>("currentderiv",currderiv);
  params.set<double>("currentresidual",currentresidual);

  return;
} //ScaTraBoundaryImpl<distype>::ElectrodeStatus


/*----------------------------------------------------------------------*
 | get constant normal                                        gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::GetConstNormal(
    LINALG::Matrix<nsd_+1,1>&          normal,
    const LINALG::Matrix<nsd_+1,nen_>&  xyze
    )
{
  // determine normal to this element
  if(not DRT::NURBS::IsNurbs(distype))
  {
  switch(nsd_)
  {
  case 2:
  {
    LINALG::Matrix<3,1> dist1(true), dist2(true);
    for (int i=0; i<3; i++)
    {
      dist1(i) = xyze(i,1)-xyze(i,0);
      dist2(i) = xyze(i,2)-xyze(i,0);
    }

    normal(0) = dist1(1)*dist2(2) - dist1(2)*dist2(1);
    normal(1) = dist1(2)*dist2(0) - dist1(0)*dist2(2);
    normal(2) = dist1(0)*dist2(1) - dist1(1)*dist2(0);
  }
  break;
  case 1:
  {
    normal(0) = xyze(1,1) - xyze(1,0);
    normal(1) = (-1.0)*(xyze(0,1) - xyze(0,0));
  }
  break;
  default:
    dserror("Illegal number of space dimensions: %d",nsd_);
    break;
  } // switch(nsd)
  }
  else // NURBS case
  {
    // ToDo: this is only a temporary solution in order to have something here.
    // Current handling of node-based normal vectors not applicable in NURBS case
#if 0
    // use one integration point at element center
    DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToStabGaussRule<distype>::rule);
    // hack: ele-id = -1
    // for nurbs elements the normal vector must be scaled with a special orientation factor!!
    // this is already part of this function call
    EvalShapeFuncAndIntFac(intpoints,0,-1,&normal);
#endif
    normal(0)= 1.0;
  }

  // length of normal to this element
  const double length = normal.Norm2();
  // outward-pointing normal of length 1.0
  if (length > EPS10)
    normal.Scale(1/length);
  else
    dserror("Zero length for element normal");

  return;
} // ScaTraBoundaryImpl<distype>::


/*----------------------------------------------------------------------*
 |  Evaluate surface/interface permeability   		                |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::EvaluateSurfacePermeability(
    const DRT::Element*        ele,
    const std::vector<double>&  ephinp,
    Epetra_SerialDenseMatrix&  emat,
    Epetra_SerialDenseVector&  erhs,
    const double               timefac,
    const double               perm
    )
{
  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // define vector for scalar values at nodes
  LINALG::Matrix<nen_,1> phinod(true);

  // loop over all scalars
  for(int k=0;k<numdofpernode_;++k)
  {
    // compute scalar values at nodes
    for (int inode=0; inode< nen_;++inode)
    {
      phinod(inode) = ephinp[inode*numdofpernode_+k];
    }

    // loop over all integration points
    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {
      const double fac = EvalShapeFuncAndIntFac(intpoints,iquad,ele->Id(),&normal_);

      // integration factor for left-hand side
      const double lhsfac = timefac*fac*perm;

      // integration factor for right-hand side
      double rhsfac = 0.0;
      if (is_incremental_ and not is_genalpha_)
        rhsfac = lhsfac;
      else
        dserror("EvaluateSurfacePermeability: Requested scheme not yet implemented");

      // matrix
      for (int vi=0; vi<nen_; ++vi)
      {
        const double vlhs = lhsfac*funct_(vi);

        const int fvi = vi*numdofpernode_+k;

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          emat(fvi,fui) += vlhs*funct_(ui);
        }
      }

      // scalar at integration point
      const double phi = funct_.Dot(phinod);

      // rhs
      const double vrhs = rhsfac*phi;
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        erhs[fvi] -= vrhs*funct_(vi);
      }
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 |  Integrate shapefunctions over surface (private)           gjb 02/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::IntegrateShapeFunctions(
    const DRT::Element*        ele,
    Teuchos::ParameterList&    params,
    Epetra_SerialDenseVector&  elevec1,
    const bool                 addarea
)
{
  // access boundary area variable with its actual value
  double boundaryint = params.get<double>("boundaryint");

  bool outputall = false;
  if(params.isParameter("alldof"))
    outputall = true;

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    const double fac = EvalShapeFuncAndIntFac(intpoints,gpid,ele->Id());

    // compute integral of shape functions
    for (int node=0;node<nen_;++node)
    {
      for (int k=0; k< numscal_; k++)
      {
        elevec1[node*numdofpernode_+k] += funct_(node) * fac;
      }
      if(outputall==true)
        elevec1[node*numdofpernode_+numdofpernode_-1] += funct_(node) * fac;
    }

    if (addarea)
    {
      // area calculation
      boundaryint += fac;
    }

  } //loop over integration points

  // add contribution to the global value
  params.set<double>("boundaryint",boundaryint);

  return;

} //ScaTraBoundaryImpl<distype>::IntegrateShapeFunction


/*----------------------------------------------------------------------*
 | calculate integral of normal diffusive flux and velocity     vg 09/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::NormDiffFluxAndVelIntegral(
    const DRT::Element*             ele,
    Teuchos::ParameterList&         params,
    const std::vector<double>&       enormdiffflux,
    const std::vector<double>&       enormvel
)
{
  // get variables for integrals of normal diffusive flux and velocity
  double normdifffluxint = params.get<double>("normal diffusive flux integral");
  double normvelint      = params.get<double>("normal velocity integral");

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    const double fac = EvalShapeFuncAndIntFac(intpoints,gpid,ele->Id());

    // compute integral of normal flux
    for (int node=0;node<nen_;++node)
    {
      normdifffluxint += funct_(node) * enormdiffflux[node] * fac;
      normvelint      += funct_(node) * enormvel[node] * fac;
    }
  } // loop over integration points

  // add contributions to the global values
  params.set<double>("normal diffusive flux integral",normdifffluxint);
  params.set<double>("normal velocity integral",normvelint);

  return;

} //ScaTraBoundaryImpl<distype>::NormDiffFluxAndVelIntegral


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType bdistype,
          DRT::Element::DiscretizationType pdistype>
   void  DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::WeakDirichlet(
     DRT::ELEMENTS::TransportBoundary*  ele,
     Teuchos::ParameterList&            params,
     DRT::Discretization&               discretization,
     Teuchos::RCP<const MAT::Material>  material,
     Epetra_SerialDenseMatrix&          elemat_epetra,
     Epetra_SerialDenseVector&          elevec_epetra)
{
  //------------------------------------------------------------------------
  // control parameters for time integration
  //------------------------------------------------------------------------
  is_stationary_  = params.get<bool>("using stationary formulation");
  is_incremental_ = params.get<bool>("incremental solver");

  // get time factor and alpha_F if required
  // one-step-Theta:    timefac = theta*dt
  // BDF2:              timefac = 2/3 * dt
  // generalized-alpha: timefac = alphaF * (gamma/alpha_M) * dt
  double timefac = 1.0;
  double alphaF  = 1.0;
  if (not is_stationary_)
  {
    timefac = params.get<double>("time factor");
    if (is_genalpha_)
    {
      alphaF = params.get<double>("alpha_F");
      timefac *= alphaF;
    }
    if (timefac < 0.0) dserror("time factor is negative.");
  }

  //------------------------------------------------------------------------
  // Dirichlet boundary condition
  //------------------------------------------------------------------------
  RCP<DRT::Condition> dbc = params.get<RCP<DRT::Condition> >("condition");

  // check of total time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const std::vector<int>* curve  = (*dbc).Get<std::vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // get values and spatial functions from condition
  // (assumed to be constant on element boundary)
  const std::vector<double>* val  = (*dbc).Get<std::vector<double> >("val"  );
  const std::vector<int>*    func = (*dbc).Get<std::vector<int> >   ("funct");

  // assign boundary value multiplied by time-curve factor
  double dirichval=(*val)[0]*curvefac;

  // spatial function number
  const int funcnum = (*func)[0];

  //------------------------------------------------------------------------
  // preliminary definitions for (boundary) and parent element and
  // evaluation of nodal values of velocity and scalar based on parent
  // element nodes
  //------------------------------------------------------------------------
  // get the parent element
  DRT::ELEMENTS::Transport* pele = ele->ParentElement();

  // number of spatial dimensions regarding (boundary) element
  static const int bnsd = DRT::UTILS::DisTypeToDim<bdistype>::dim;

  // number of spatial dimensions regarding parent element
  static const int pnsd = DRT::UTILS::DisTypeToDim<pdistype>::dim;

  // number of (boundary) element nodes
  static const int bnen = DRT::UTILS::DisTypeToNumNodePerEle<bdistype>::numNodePerElement;

  // number of parent element nodes
  static const int pnen = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

  // parent element lm vector
  std::vector<int>  plm ;
  std::vector<int>  plmowner;
  std::vector<int>  plmstride;
  pele->LocationVector(discretization,plm,plmowner,plmstride);

  // get velocity values at parent element nodes
  const Teuchos::RCP<Epetra_MultiVector> velocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("convective velocity field",Teuchos::null);
  Epetra_SerialDenseVector evel(pnsd*pnen);
  DRT::UTILS::ExtractMyNodeBasedValues(pele,evel,velocity,pnsd);

  // get scalar values at parent element nodes
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

  // extract local values from global vectors for parent element
  std::vector<double> myphinp(plm.size());
  DRT::UTILS::ExtractMyValues(*phinp,myphinp,plm);

  // matrix and vector definition
  LINALG::Matrix<pnsd,pnen>       evelnp;
  std::vector<LINALG::Matrix<pnen,1> > ephinp(numscal_);

  // insert into element arrays
  for (int i=0;i<pnen;++i)
  {
    for (int k = 0; k< numscal_; ++k)
    {
      // split for each tranported scalar, insert into element arrays
      ephinp[k](i,0) = myphinp[k+(i*numdofpernode_)];
    }

    // insert velocity field into element array
    for (int idim=0 ; idim < pnsd; idim++)
    {
      evelnp(idim,i) = evel[idim + i*pnsd];
    }
  }

  //------------------------------------------------------------------------
  // preliminary definitions for integration loop
  //------------------------------------------------------------------------
  // reshape element matrices and vectors and init to zero, construct views
  elemat_epetra.Shape(pnen,pnen);
  elevec_epetra.Size (pnen);
  LINALG::Matrix<pnen,pnen> emat(elemat_epetra.A(),true);
  LINALG::Matrix<pnen,   1> erhs(elevec_epetra.A(),true);

  // (boundary) element local node coordinates
  LINALG::Matrix<pnsd,bnen>  bxyze(true);
  GEO::fillInitialPositionArray<bdistype,pnsd,LINALG::Matrix<pnsd,bnen> >(ele,bxyze);

  // parent element local node coordinates
  LINALG::Matrix<pnsd,pnen>  pxyze(true);
  GEO::fillInitialPositionArray<pdistype,pnsd,LINALG::Matrix<pnsd,pnen> >(pele,pxyze);

  // coordinates of integration points for (boundary) and parent element
  LINALG::Matrix<bnsd,   1>  bxsi(true);
  LINALG::Matrix<pnsd,   1>  pxsi(true);

  // transposed jacobian "dx/ds" and inverse of transposed jacobian "ds/dx"
  // for parent element
  LINALG::Matrix<pnsd,pnsd>  pxjm(true);
  LINALG::Matrix<pnsd,pnsd>  pxji(true);

  // metric tensor for (boundary) element
  LINALG::Matrix<bnsd,bnsd>  bmetrictensor(true);

  // (outward-pointing) unit normal vector to (boundary) element
  LINALG::Matrix<pnsd,   1>  bnormal(true);

  // velocity vector at integration point
  LINALG::Matrix<pnsd,   1>  velint;

  // gradient of scalar value at integration point
  LINALG::Matrix<pnsd,1> gradphi;

  // (boundary) element shape functions, local and global derivatives
  LINALG::Matrix<bnen,   1>  bfunct(true);
  LINALG::Matrix<bnsd,bnen>  bderiv(true);
  LINALG::Matrix<bnsd,bnen>  bderxy(true);

  // parent element shape functions, local and global derivatives
  LINALG::Matrix<pnen,   1>  pfunct(true);
  LINALG::Matrix<pnsd,pnen>  pderiv(true);
  LINALG::Matrix<pnsd,pnen>  pderxy(true);

  //------------------------------------------------------------------------
  // additional matrices and vectors for mixed-hybrid formulation
  //------------------------------------------------------------------------
  // for volume integrals
  LINALG::Matrix<pnsd*pnen,pnsd*pnen> mat_s_q(true);
  LINALG::Matrix<pnsd*pnen,     pnen> mat_s_gradphi(true);

  LINALG::Matrix<pnsd*pnen,        1> vec_s_gradphi(true);

  // for boundary integrals
  LINALG::Matrix<     pnen,pnsd*pnen> mat_w_q_o_n(true);
  LINALG::Matrix<pnsd*pnen,     pnen> mat_s_o_n_phi(true);

  LINALG::Matrix<pnsd*pnen,        1> vec_s_o_n_phi_minus_g(true);

  // inverse matrix
  LINALG::Matrix<pnsd*pnen,pnsd*pnen> inv_s_q(true);

  //------------------------------------------------------------------------
  // check whether Nitsche (default) or mixed-hybrid formulation as well as
  // preliminary definitions and computations for Nitsche stabilization term
  //------------------------------------------------------------------------
  // default is Nitsche formulation
  bool mixhyb = false;

  // stabilization parameter for Nitsche term
  const double nitsche_stab_para = (*dbc).GetDouble("TauBscaling");

  // if stabilization parameter negative: mixed-hybrid formulation
  if (nitsche_stab_para < 0.0) mixhyb = true;

  // pre-factor for adjoint-consistency term:
  // either 1.0 (adjoint-consistent, default) or -1.0 (adjoint-inconsistent)
  double gamma = 1.0;
  const string* consistency = (*dbc).Get<string>("Choice of gamma parameter");
  if      (*consistency=="adjoint-consistent") gamma = 1.0;
  else if (*consistency=="diffusive-optimal")  gamma = -1.0;
  else dserror("unknown definition for gamma parameter: %s",(*consistency).c_str());

  // use one-point Gauss rule to do calculations at element center
  DRT::UTILS::IntPointsAndWeights<bnsd> intpoints_tau(SCATRA::DisTypeToStabGaussRule<bdistype>::rule);

  // element surface area (1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double* gpcoord = (intpoints_tau.IP().qxg)[0];
  for (int idim=0;idim<bnsd;idim++)
  {
    bxsi(idim) = gpcoord[idim];
  }
  DRT::UTILS::shape_function_deriv1<bdistype>(bxsi,bderiv);
  double drs = 0.0;
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<bdistype>(bxyze,bderiv,bmetrictensor,drs,&bnormal);
  const double area = intpoints_tau.IP().qwgt[0]*drs;

  // get number of dimensions for (boundary) element (convert from int to double)
  const double dim = (double) bnsd;

  // computation of characteristic length of (boundary) element
  // (2D: square root of element area, 1D: element length)
  const double h = std::pow(area,(1.0/dim));

  //------------------------------------------------------------------------
  // preliminary computations for integration loop
  //------------------------------------------------------------------------
  // integrations points and weights for (boundary) element and parent element
  const DRT::UTILS::IntPointsAndWeights<bnsd> bintpoints(SCATRA::DisTypeToOptGaussRule<bdistype>::rule);

  const DRT::UTILS::IntPointsAndWeights<pnsd> pintpoints(SCATRA::DisTypeToOptGaussRule<pdistype>::rule);

  // transfer integration-point coordinates of (boundary) element to parent element
  Epetra_SerialDenseMatrix pqxg(pintpoints.IP().nquad,pnsd);
  {
    Epetra_SerialDenseMatrix gps(bintpoints.IP().nquad,bnsd);

    for (int iquad=0; iquad<bintpoints.IP().nquad; ++iquad)
    {
      const double* gpcoord = (bintpoints.IP().qxg)[iquad];

      for (int idim=0;idim<bnsd ;idim++)
      {
        gps(iquad,idim) = gpcoord[idim];
      }
    }
    if(pnsd==2)
    {
      DRT::UTILS::BoundaryGPToParentGP2(pqxg,gps,pdistype,bdistype,ele->BeleNumber());
    }
    else if (pnsd==3)
    {
      DRT::UTILS::BoundaryGPToParentGP3(pqxg,gps,pdistype,bdistype,ele->BeleNumber());
    }

  }

  //------------------------------------------------------------------------
  // integration loop 1: volume integrals (only for mixed-hybrid formulation)
  //------------------------------------------------------------------------
  if (mixhyb)
  {
    for (int iquad=0; iquad<pintpoints.IP().nquad; ++iquad)
    {
      // reference coordinates of integration point from (boundary) element
      const double* gpcoord = (pintpoints.IP().qxg)[iquad];
      for (int idim=0;idim<pnsd;idim++)
      {
        pxsi(idim) = gpcoord[idim];
      }

      // parent element shape functions and local derivatives
      DRT::UTILS::shape_function       <pdistype>(pxsi,pfunct);
      DRT::UTILS::shape_function_deriv1<pdistype>(pxsi,pderiv);

      // Jacobian matrix and determinant of parent element (including check)
      pxjm.MultiplyNT(pderiv,pxyze);
      const double det = pxji.Invert(pxjm);
      if (det < 1E-16) dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", pele->Id(), det);

      // compute integration factor
      const double fac = pintpoints.IP().qwgt[iquad]*det;

      // compute global derivatives
      pderxy.Multiply(pxji,pderiv);

      //--------------------------------------------------------------------
      // loop over scalars (not yet implemented for more than one scalar)
      //--------------------------------------------------------------------
      // for(int k=0;k<numdofpernode_;++k)
      int k=0;
      {
        // get viscosity
        if (material->MaterialType() == INPAR::MAT::m_scatra)
        {
          const MAT::ScatraMat* actmat = static_cast<const MAT::ScatraMat*>(material.get());

          dsassert(numdofpernode_==1,"more than 1 dof per node for SCATRA material");

          // get constant diffusivity
          diffus_[k] = actmat->Diffusivity();
        }
        else dserror("Material type is not supported");

        // gradient of current scalar value
        gradphi.Multiply(pderxy,ephinp[k]);

        // integration factor for left-hand side
        const double lhsfac = timefac*fac;

        // integration factor for right-hand side
        double rhsfac = 0.0;
        if (is_incremental_ and is_genalpha_)
          rhsfac = lhsfac/alphaF;
        else if (not is_incremental_ and is_genalpha_)
          rhsfac = lhsfac*(1.0-alphaF)/alphaF;
        else if (is_incremental_ and not is_genalpha_)
          rhsfac = lhsfac;

        //--------------------------------------------------------------------
        //  matrix and vector additions due to mixed-hybrid formulation
        //--------------------------------------------------------------------
        /*
                       /         \
                  1   |   h   h  |
              - ----- |  s , q   |
                kappa |          |
                      \          / Omega
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          //const double vlhs = lhsfac*pfunct(vi);
          const double vlhs = lhsfac*(1.0/diffus_[k])*pfunct(vi);

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            for(int i=0;i<pnsd;++i)
            {
              mat_s_q(fvi*pnsd+i,fui*pnsd+i) -= vlhs*pfunct(ui);
            }
          }
        }

        /*
                       /                  \
                      |  h         /   h\  |
                    + | s  , grad | phi  | |
                      |            \    /  |
                       \                  / Omega
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          //const double vlhs = lhsfac*diffus_[k]*pfunct(vi);
          const double vlhs = lhsfac*pfunct(vi);

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            for(int i=0;i<pnsd;++i)
            {
              mat_s_gradphi(fvi*pnsd+i,fui) += vlhs*pderxy(i,ui);
            }
          }

          //const double vrhs = rhsfac*diffus_[k]*pfunct(vi);
          const double vrhs = rhsfac*pfunct(vi);

          for(int i=0;i<pnsd;++i)
          {
            vec_s_gradphi(fvi*pnsd+i) += vrhs*gradphi(i);
          }
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // integration loop 2: boundary integrals
  //------------------------------------------------------------------------
  for (int iquad=0; iquad<bintpoints.IP().nquad; ++iquad)
  {
    // reference coordinates of integration point from (boundary) element
    const double* gpcoord = (bintpoints.IP().qxg)[iquad];
    for (int idim=0;idim<bnsd;idim++)
    {
      bxsi(idim) = gpcoord[idim];
    }

    // (boundary) element shape functions
    DRT::UTILS::shape_function       <bdistype>(bxsi,bfunct);
    DRT::UTILS::shape_function_deriv1<bdistype>(bxsi,bderiv);

    // global coordinates of current integration point from (boundary) element
    LINALG::Matrix<pnsd,1> coordgp(true);
    for (int A=0;A<bnen;++A)
    {
      for(int j=0;j<pnsd;++j)
      {
        coordgp(j)+=bxyze(j,A)*bfunct(A);
      }
    }

    // reference coordinates of integration point from parent element
    for (int idim=0;idim<pnsd;idim++)
    {
      pxsi(idim) = pqxg(iquad,idim);
    }

    // parent element shape functions and local derivatives
    DRT::UTILS::shape_function       <pdistype>(pxsi,pfunct);
    DRT::UTILS::shape_function_deriv1<pdistype>(pxsi,pderiv);

    // Jacobian matrix and determinant of parent element (including check)
    pxjm.MultiplyNT(pderiv,pxyze);
    const double det = pxji.Invert(pxjm);
    if (det < 1E-16) dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", pele->Id(), det);

    // compute measure tensor for surface element, infinitesimal area element drs
    // and (outward-pointing) unit normal vector
    DRT::UTILS::ComputeMetricTensorForBoundaryEle<bdistype>(bxyze,bderiv,bmetrictensor,drs,&bnormal);

    // for nurbs elements the normal vector must be scaled with a special orientation factor!!
    if(DRT::NURBS::IsNurbs(distype))
      bnormal.Scale(normalfac_);

    // compute integration factor
    const double fac = bintpoints.IP().qwgt[iquad]*drs;

    // compute global derivatives
    pderxy.Multiply(pxji,pderiv);

#if 1
    //--------------------------------------------------------------------
    // check whether integration-point coordinates evaluated from
    // (boundary) and parent element match
    //--------------------------------------------------------------------
    LINALG::Matrix<pnsd,1> check(true);
    LINALG::Matrix<pnsd,1> diff(true);

    for (int A=0;A<pnen;++A)
    {
      for(int j=0;j<pnsd;++j)
      {
        check(j)+=pxyze(j,A)*pfunct(A);
      }
    }

    diff=check;
    diff-=coordgp;

    const double norm=diff.Norm2();

    if (norm>1e-9)
    {
      for (int j=0;j<pnsd;++j)
      {
        printf("%12.5e %12.5e\n",check(j),coordgp(j));
      }
      dserror("Gausspoint matching error %12.5e\n",norm);
    }
#endif

    //--------------------------------------------------------------------
    // factor for Dirichlet boundary condition given by spatial function
    //--------------------------------------------------------------------
    double functfac = 1.0;
    if (funcnum > 0)
    {
      // evaluate function at current integration point (important: a 3D position vector is required)
      double coordgp3D[3];
      coordgp3D[0]=0.0;
      coordgp3D[1]=0.0;
      coordgp3D[2]=0.0;
      for (int i=0; i<pnsd;i++)
        coordgp3D[i]=coordgp(i);

      functfac = DRT::Problem::Instance()->Funct(funcnum-1).Evaluate(0,&(coordgp3D[0]),time,NULL);
    }
    else functfac = 1.0;
    dirichval *= functfac;

    //--------------------------------------------------------------------
    // loop over scalars (not yet implemented for more than one scalar)
    //--------------------------------------------------------------------
    // for(int k=0;k<numdofpernode_;++k)
    int k=0;
    {
      // get viscosity
      if (material->MaterialType() == INPAR::MAT::m_scatra)
      {
        const MAT::ScatraMat* actmat = static_cast<const MAT::ScatraMat*>(material.get());

        dsassert(numdofpernode_==1,"more than 1 dof per node for SCATRA material");

        // get constant diffusivity
        diffus_[k] = actmat->Diffusivity();
      }
      else dserror("Material type is not supported");

      // get scalar value at integration point
      const double phi = pfunct.Dot(ephinp[k]);

      // integration factor for left-hand side
      const double lhsfac = timefac*fac;

      // integration factor for right-hand side
      double rhsfac = 0.0;
      if (is_incremental_ and is_genalpha_)
        rhsfac = lhsfac/alphaF;
      else if (not is_incremental_ and is_genalpha_)
        rhsfac = lhsfac*(1.0-alphaF)/alphaF;
      else if (is_incremental_ and not is_genalpha_)
        rhsfac = lhsfac;

      if (mixhyb)
      {
        //--------------------------------------------------------------------
        //  matrix and vector additions due to mixed-hybrid formulation
        //--------------------------------------------------------------------
        /*  consistency term
                    /           \
                   |  h   h     |
                 - | w , q  o n |
                   |            |
                   \            / Gamma
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          const double vlhs = lhsfac*pfunct(vi);

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            for(int i=0;i<pnsd;++i)
            {
              mat_w_q_o_n(fvi,fui*pnsd+i) -= vlhs*pfunct(ui)*bnormal(i);
            }
          }
        }

        /*  adjoint consistency term
                    /                 \
                   |  h          h    |
                 - | s  o n , phi - g |
                   |                  |
                   \                  / Gamma
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          const double vlhs = lhsfac*pfunct(vi);

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            for(int i=0;i<pnsd;++i)
            {
              mat_s_o_n_phi(fvi*pnsd+i,fui) -= vlhs*pfunct(ui)*bnormal(i);
            }
          }

          for(int i=0;i<pnsd;++i)
          {
            vec_s_o_n_phi_minus_g(fvi*pnsd+i) -= pfunct(vi)*bnormal(i)*(rhsfac*phi - timefac*fac*dirichval);
          }
        }
      }
      else
      {
        // parameter alpha for Nitsche stabilization term
        const double alpha = nitsche_stab_para*diffus_[k]/h;

        // get velocity at integration point
        velint.Multiply(evelnp,pfunct);

        // normal velocity
        const double normvel = velint.Dot(bnormal);

        // gradient of current scalar value
        gradphi.Multiply(pderxy,ephinp[k]);

        // gradient of current scalar value in normal direction
        const double gradphi_norm = bnormal.Dot(gradphi);

        //--------------------------------------------------------------------
        //  matrix and vector additions due to Nitsche formulation
        //--------------------------------------------------------------------
        /*  consistency term
                    /                           \
                   |  h                  h      |
                 - | w , kappa * grad(phi ) o n |
                   |                            |
                   \                            / Gamma
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          const double vlhs = lhsfac*pfunct(vi)*diffus_[k];

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            for(int i=0;i<pnsd;++i)
            {
              emat(fvi,fui) -= vlhs*pderxy(i,ui)*bnormal(i);
            }
          }

          const double vrhs = rhsfac*diffus_[k];

          erhs(fvi) += vrhs*pfunct(vi)*gradphi_norm;
        }

        /*  adjoint consistency term, inflow/outflow part
              / --          --                                        \
             |  |         h  |                      h           h     |
           - |  |(a o n) w  +| gamma * kappa *grad(w ) o n , phi - g  |
             |  |            |                                        |
             \  --          --                                        / Gamma_in/out
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          // compute diffusive part
          double prefac = 0.0;
          for(int i=0;i<pnsd;++i)
          {
            prefac += gamma*diffus_[k]*pderxy(i,vi)*bnormal(i);
          }

          // add convective part in case of inflow boundary
          if (normvel<-0.0001) prefac += normvel*pfunct(vi);

          const double vlhs = lhsfac*prefac;

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            emat(fvi,fui) -= vlhs*pfunct(ui);
          }

          erhs(fvi) += prefac*(rhsfac*phi - timefac*fac*dirichval);
        }

        /*  stabilization term
                            /             \
                           |  h     h     |
                 + alpha * | w , phi - g  |
                           |              |
                           \              / Gamma
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          const double prefac = alpha*pfunct(vi);

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            emat(fvi,fui) += lhsfac*prefac*pfunct(ui);
          }

          erhs(fvi) -= prefac*(rhsfac*phi - timefac*fac*dirichval);
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // local condensation (only for mixed-hybrid formulation)
  //------------------------------------------------------------------------
  if (mixhyb)
  {
    // matrix inversion of flux-flux block
    inv_s_q = mat_s_q;

    LINALG::FixedSizeSerialDenseSolver<pnsd*pnen,pnsd*pnen> solver;

    solver.SetMatrix(inv_s_q);
    solver.Invert();

    // computation of matrix-matrix and matrix vector products, local assembly
    for (int vi=0; vi<pnen; ++vi)
    {
      for (int ui=0; ui<pnen; ++ui)
      {
        for(int rr=0; rr<pnsd*pnen; ++rr)
        {
          for(int mm=0; mm<pnsd*pnen; ++mm)
          {
            emat(vi,ui)
              -=mat_w_q_o_n(vi,rr)*inv_s_q(rr,mm)*(mat_s_gradphi(mm,ui)+mat_s_o_n_phi(mm,ui));
          }
        }
      }
    }

    for (int vi=0; vi<pnen; ++vi)
    {
      for(int rr=0; rr<pnsd*pnen; ++rr)
      {
        for(int mm=0; mm<pnsd*pnen; ++mm)
        {
          erhs(vi)-= mat_w_q_o_n(vi,rr)*inv_s_q(rr,mm)*(-vec_s_o_n_phi_minus_g(mm)-vec_s_gradphi(mm));
        }
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | calculate boundary conditions for expl TaylorGalerkin3 and           |
 | implicit Characteristic Galerkin time integration                    |
 | just for the linear transport equation                  schott 04/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType bdistype,
          DRT::Element::DiscretizationType pdistype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::TaylorGalerkinBoundaryOutflow(
    DRT::ELEMENTS::TransportBoundary*  ele,                  //!< transport element
    Teuchos::ParameterList&            params,               //!< parameter list
    DRT::Discretization&               discretization,       //!< discretization
    Teuchos::RCP<const MAT::Material>  material,             //!< material
    Epetra_SerialDenseMatrix&          elemat_epetra,        //!< ele sysmat
    Epetra_SerialDenseVector&          elevec_epetra         //!< ele rhs
 )
{
  INPAR::SCATRA::TimeIntegrationScheme timealgo = DRT::INPUT::get<INPAR::SCATRA::TimeIntegrationScheme>(params,"timealgo");


  //------------------------------------------------------------------------
  // preliminary definitions for (boundary) and parent element and
  // evaluation of nodal values of velocity and scalar based on parent
  // element nodes
  //------------------------------------------------------------------------
  // get the parent element
  DRT::ELEMENTS::Transport* pele = ele->ParentElement();

  // number of spatial dimensions regarding (boundary) element
  static const int bnsd = DRT::UTILS::DisTypeToDim<bdistype>::dim;

  // number of spatial dimensions regarding parent element
  static const int pnsd = DRT::UTILS::DisTypeToDim<pdistype>::dim;

  // number of (boundary) element nodes
  static const int bnen = DRT::UTILS::DisTypeToNumNodePerEle<bdistype>::numNodePerElement;

  // number of parent element nodes
  static const int pnen = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

  // parent element lm vector
  std::vector<int>  plm;
  std::vector<int>  plmowner;
  std::vector<int>  plmstride;
  pele->LocationVector(discretization,plm,plmowner,plmstride);

  // get velocity values at parent element nodes
  const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("convective velocity field",Teuchos::null);
  Epetra_SerialDenseVector evel(pnsd*pnen);
  DRT::UTILS::ExtractMyNodeBasedValues(pele,evel,velocity,pnsd);

  // get scalar values at parent element nodes
  RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");
  RCP<const Epetra_Vector> phin = discretization.GetState("phin");
  if (phinp==Teuchos::null) dserror("Cannot get state vector 'phin'");


  // extract local values from global vectors for parent element
  std::vector<double> myphinp(plm.size());
  DRT::UTILS::ExtractMyValues(*phinp,myphinp,plm);

  std::vector<double> myphin(plm.size());
  DRT::UTILS::ExtractMyValues(*phin,myphin,plm);

  //	  // matrix and vector definition
  LINALG::Matrix<pnsd,pnen>       evelnp;
  std::vector<LINALG::Matrix<pnen,1> > ephinp(numscal_);
  std::vector<LINALG::Matrix<pnen,1> > ephin(numscal_);

  // insert into element arrays
  for (int i=0;i<pnen;++i)
  {
    for (int k = 0; k< numscal_; ++k)
    {
      // split for each tranported scalar, insert into element arrays
      ephinp[k](i,0) = myphinp[k+(i*numdofpernode_)];
      ephin[k](i,0)  = myphin[k+(i*numdofpernode_)];
    }

    // insert velocity field into element array
    for (int idim=0 ; idim < pnsd; idim++)
    {
      evelnp(idim,i) = evel[idim + i*pnsd];
    }
  }

  //------------------------------------------------------------------------
  // preliminary definitions for integration loop
  //------------------------------------------------------------------------
  // reshape element matrices and vectors and init to zero, construct views
  elemat_epetra.Shape(pnen,pnen);
  elevec_epetra.Size (pnen);
  LINALG::Matrix<pnen,pnen> emat(elemat_epetra.A(),true);
  LINALG::Matrix<pnen,   1> erhs(elevec_epetra.A(),true);

  // (boundary) element local node coordinates
  LINALG::Matrix<pnsd,bnen>  bxyze(true);
  GEO::fillInitialPositionArray<bdistype,pnsd,LINALG::Matrix<pnsd,bnen> >(ele,bxyze);

  // parent element local node coordinates
  LINALG::Matrix<pnsd,pnen>  pxyze(true);
  GEO::fillInitialPositionArray<pdistype,pnsd,LINALG::Matrix<pnsd,pnen> >(pele,pxyze);

  // coordinates of integration points for (boundary) and parent element
  LINALG::Matrix<bnsd,   1>  bxsi(true);
  LINALG::Matrix<pnsd,   1>  pxsi(true);

  // transposed jacobian "dx/ds" and inverse of transposed jacobian "ds/dx"
  // for parent element
  LINALG::Matrix<pnsd,pnsd>  pxjm(true);
  LINALG::Matrix<pnsd,pnsd>  pxji(true);

  // metric tensor for (boundary) element
  LINALG::Matrix<bnsd,bnsd>  bmetrictensor(true);

  // (outward-pointing) unit normal vector to (boundary) element
  LINALG::Matrix<pnsd,   1>  bnormal(true);

  // velocity vector at integration point
  LINALG::Matrix<pnsd,   1>  velint;

  // gradient of scalar value at integration point
  LINALG::Matrix<pnsd,1> gradphi;

  // (boundary) element shape functions, local and global derivatives
  LINALG::Matrix<bnen,   1>  bfunct(true);
  LINALG::Matrix<bnsd,bnen>  bderiv(true);
  LINALG::Matrix<bnsd,bnen>  bderxy(true);

  // parent element shape functions, local and global derivatives
  LINALG::Matrix<pnen,   1>  pfunct(true);
  LINALG::Matrix<pnsd,pnen>  pderiv(true);
  LINALG::Matrix<pnsd,pnen>  pderxy(true);


  // use one-point Gauss rule to do calculations at element center
  DRT::UTILS::IntPointsAndWeights<bnsd> intpoints_tau(SCATRA::DisTypeToStabGaussRule<bdistype>::rule);

  // element surface area (1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double* gpcoord = (intpoints_tau.IP().qxg)[0];
  for (int idim=0;idim<bnsd;idim++)
  {
    bxsi(idim) = gpcoord[idim];
  }
  DRT::UTILS::shape_function_deriv1<bdistype>(bxsi,bderiv);
  double drs = 0.0;
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<bdistype>(bxyze,bderiv,bmetrictensor,drs,&bnormal);

  //------------------------------------------------------------------------
  // preliminary computations for integration loop
  //------------------------------------------------------------------------
  // integrations points and weights for (boundary) element and parent element
  const DRT::UTILS::IntPointsAndWeights<bnsd> bintpoints(SCATRA::DisTypeToOptGaussRule<bdistype>::rule);

  const DRT::UTILS::IntPointsAndWeights<pnsd> pintpoints(SCATRA::DisTypeToOptGaussRule<pdistype>::rule);

  // transfer integration-point coordinates of (boundary) element to parent element
  Epetra_SerialDenseMatrix pqxg(pintpoints.IP().nquad,pnsd);
  {
    Epetra_SerialDenseMatrix gps(bintpoints.IP().nquad,bnsd);

    for (int iquad=0; iquad<bintpoints.IP().nquad; ++iquad)
    {
      const double* gpcoord = (bintpoints.IP().qxg)[iquad];

      for (int idim=0;idim<bnsd ;idim++)
      {
        gps(iquad,idim) = gpcoord[idim];
      }
    }
    if(pnsd==2)
    {
      DRT::UTILS::BoundaryGPToParentGP2(pqxg,gps,pdistype,bdistype,ele->BeleNumber());
    }
    else if (pnsd==3)
    {
      DRT::UTILS::BoundaryGPToParentGP3(pqxg,gps,pdistype,bdistype,ele->BeleNumber());
    }
  }


  const double dt = params.get<double>("time_step_size");

  //------------------------------------------------------------------------
  // integration loop: boundary integrals
  //------------------------------------------------------------------------
  for (int iquad=0; iquad<bintpoints.IP().nquad; ++iquad)
  {
    // reference coordinates of integration point from (boundary) element
    const double* gpcoord = (bintpoints.IP().qxg)[iquad];
    for (int idim=0;idim<bnsd;idim++)
    {
      bxsi(idim) = gpcoord[idim];
    }

    // (boundary) element shape functions
    DRT::UTILS::shape_function       <bdistype>(bxsi,bfunct);
    DRT::UTILS::shape_function_deriv1<bdistype>(bxsi,bderiv);

    // global coordinates of current integration point from (boundary) element
    LINALG::Matrix<pnsd,1> coordgp(true);
    for (int A=0;A<bnen;++A)
    {
      for(int j=0;j<pnsd;++j)
      {
        coordgp(j)+=bxyze(j,A)*bfunct(A);
      }
    }

    // reference coordinates of integration point from parent element
    for (int idim=0;idim<pnsd;idim++)
    {
      pxsi(idim) = pqxg(iquad,idim);
    }

    // parent element shape functions and local derivatives
    DRT::UTILS::shape_function       <pdistype>(pxsi,pfunct);
    DRT::UTILS::shape_function_deriv1<pdistype>(pxsi,pderiv);

    // Jacobian matrix and determinant of parent element (including check)
    pxjm.MultiplyNT(pderiv,pxyze);
    const double det = pxji.Invert(pxjm);
    if (det < 1E-16) dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", pele->Id(), det);

    // compute measure tensor for surface element, infinitesimal area element drs
    // and (outward-pointing) unit normal vector
    DRT::UTILS::ComputeMetricTensorForBoundaryEle<bdistype>(bxyze,bderiv,bmetrictensor,drs,&bnormal);

    // for nurbs elements the normal vector must be scaled with a special orientation factor!!
    if(DRT::NURBS::IsNurbs(distype))
      bnormal.Scale(normalfac_);

    // compute integration factor
    const double fac_surface = bintpoints.IP().qwgt[iquad]*drs;

    // compute global derivatives
    pderxy.Multiply(pxji,pderiv);


    // decide if inflow or outflow



    if(timealgo == INPAR::SCATRA::timeint_tg2)
    {
      //--------------------------------------------------------------------
      // loop over scalars (not yet implemented for more than one scalar)
      //--------------------------------------------------------------------
      for(int dofindex=0;dofindex<numdofpernode_;++dofindex)
      {

        // get velocity at integration point
        velint.Multiply(evelnp,pfunct);

        // normal velocity
        const double normvel = velint.Dot(bnormal);
        //----------  --------------      |                          |
        //  mat              -1/4* d^2    | w (u*n), u*grad(D(psi) ) |
        //--------------------------      |                          |

        LINALG::Matrix<1,pnen> derxy_vel;
        derxy_vel.Clear();
        derxy_vel.MultiplyTN(velint,pderxy);

        // decide if inflow or outflow
        bool assemble_inflow_outflow = true;

        if(assemble_inflow_outflow)
        {

          for (int vi=0; vi<pnen; ++vi)
          {
            const int fvi = vi*numdofpernode_+dofindex;

            for (int ui=0; ui<pnen; ++ui)
            {
              const int fui = ui*numdofpernode_+dofindex;

              emat(fvi,fui) -= pfunct(vi)* normvel*(fac_surface*dt*dt/4.0) * derxy_vel(0,ui);
            }
          }

          //----------  --------------      |                   n+1     n      |
          //  rhs               1/4* dt^2   | w(u*n), u*grad(phi   + phi   )   |
          //--------------------------      |                                  |

          // update grad_dist_n
          LINALG::Matrix<pnsd,1> grad_dist_n(true);
          grad_dist_n.Multiply(pderxy,ephin[dofindex]);

          LINALG::Matrix<pnsd,1> grad_dist_npi(true);
          grad_dist_npi.Multiply(pderxy,ephinp[dofindex]);

          LINALG::Matrix<pnsd,1> grad_dist_sum(true);
          grad_dist_sum.Update(1.0,grad_dist_n, 1.0, grad_dist_npi);


          LINALG::Matrix<1,1> uGradDistSum;
          uGradDistSum.Clear();
          uGradDistSum.MultiplyTN(grad_dist_sum,velint);

          for (int vi=0; vi<pnen; ++vi)
          {
            const int fvi = vi*numdofpernode_+dofindex;

            erhs(fvi) += pfunct(vi)*dt*dt/4.0*fac_surface *normvel* uGradDistSum(0,0);
          }
        } // end if(assemble_inflow_outflow)

      } // loop over scalars
    } // end if tg2
    else if(timealgo == INPAR::SCATRA::timeint_tg3)
    {
      //--------------------------------------------------------------------
      // loop over scalars (not yet implemented for more than one scalar)
      //--------------------------------------------------------------------
      for(int dofindex=0;dofindex<numdofpernode_;++dofindex)
      {

        // get velocity at integration point
        velint.Multiply(evelnp,pfunct);

        // get velocity at integration point
        double phin  =  ephin[dofindex].Dot(pfunct);

        // normal velocity
        const double normvel = velint.Dot(bnormal);

        //bool outflow_point = false;
        //if(normvel > 0.0) outflow_point=true;


        LINALG::Matrix<1,pnen> derxy_vel;
        derxy_vel.Clear();
        derxy_vel.MultiplyTN(velint,pderxy);

        // decide if inflow or outflow
        //			    bool assemble_inflow_outflow = true;
        //		    if(normvel<0.0) assemble_inflow_outflow=false;


        //----------  --------------      |                          |
        //  mat              -1/6*dt^2    | w (a*n), a*grad(D(phi) ) |
        //--------------------------      |                          |


        double fac_term = fac_surface*dt*dt/6.0;

        //			    	if(outflow_point)
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+dofindex;

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+dofindex;

            emat(fvi,fui) -= fac_term * pfunct(vi)* normvel* derxy_vel(0,ui);
          }
        }


        // update grad_dist_n
        LINALG::Matrix<pnsd,1> grad_dist_n(true);
        grad_dist_n.Multiply(pderxy,ephin[dofindex]);

        LINALG::Matrix<pnsd,1> grad_dist_npi(true);
        grad_dist_npi.Multiply(pderxy,ephinp[dofindex]);

        // a*grad(phi_n)
        double a_phi_n = velint.Dot(grad_dist_n);

        // a*grad(phi_npi)
        double a_phi_npi = velint.Dot(grad_dist_npi);

        //----------  --------------      |                   n+1    n |
        //  rhs               1/6* dt^2   | w(a*n), a*grad(phi   -phi) |
        //--------------------------      |                   i        |

        fac_term = fac_surface*dt*dt/6.0 * (a_phi_npi-a_phi_n) *normvel;

        // if(outflow_point)
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+dofindex;
          erhs(fvi) += pfunct(vi) * fac_term;
        }


        //----------  --------------      |            n  |
        //  rhs                     -dt   | w(a*n), phi   |
        //--------------------------      |               |

        fac_term = dt*fac_surface *normvel* phin;

        // if(outflow_point)
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+dofindex;
          erhs(fvi) -= pfunct(vi)*fac_term;
        }


        //----------  --------------      |                    n  |
        //  mat              +1/2*dt^2    | w (a*n), a*grad(phi ) |
        //--------------------------      |                       |

        fac_term= 0.5*dt*dt*normvel*fac_surface*a_phi_n;

        // if(outflow_point)
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+dofindex;
          erhs(fvi) += pfunct(vi)*fac_term;
        }
      } //dofindex
    } // end if timeint_tg3
    else dserror("no valid timealgo for TaylorGalerkinBoundaryOutflow call");

  } // loop over integration points

  return;
}


/*----------------------------------------------------------------------*
 | calculate boundary conditions for                                    |
 | impl. Characteristic Galerkin (2nd order)                            |
 | time integration just for the reinitialization equation schott 04/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType bdistype,
          DRT::Element::DiscretizationType pdistype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::ReinitCharacteristicGalerkinBoundary(
    DRT::ELEMENTS::TransportBoundary*  ele,                  //!< transport element
    Teuchos::ParameterList&            params,               //!< parameter list
    DRT::Discretization&               discretization,       //!< discretization
    Teuchos::RCP<const MAT::Material>  material,             //!< material
    Epetra_SerialDenseMatrix&          elemat_epetra,        //!< ele sysmat
    Epetra_SerialDenseVector&          elevec_epetra         //!< ele rhs
    )
{

  //------------------------------------------------------------------------
  // preliminary definitions for (boundary) and parent element and
  // evaluation of nodal values of velocity and scalar based on parent
  // element nodes
  //------------------------------------------------------------------------
  // get the parent element
  DRT::ELEMENTS::Transport* pele = ele->ParentElement();

  // number of spatial dimensions regarding (boundary) element
  static const int bnsd = DRT::UTILS::DisTypeToDim<bdistype>::dim;

  // number of spatial dimensions regarding parent element
  static const int pnsd = DRT::UTILS::DisTypeToDim<pdistype>::dim;

  // number of (boundary) element nodes
  static const int bnen = DRT::UTILS::DisTypeToNumNodePerEle<bdistype>::numNodePerElement;

  // number of parent element nodes
  static const int pnen = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

  // parent element lm vector
  std::vector<int>  plm ;
  std::vector<int>  plmowner;
  std::vector<int>  plmstride;
  pele->LocationVector(discretization,plm,plmowner,plmstride);

  // get scalar values at parent element nodes
  RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");
  RCP<const Epetra_Vector> phin = discretization.GetState("phin");
  if (phinp==Teuchos::null) dserror("Cannot get state vector 'phin'");

  // extract local values from global vectors for parent element
  std::vector<double> myphinp(plm.size());
  DRT::UTILS::ExtractMyValues(*phinp,myphinp,plm);

  std::vector<double> myphin(plm.size());
  DRT::UTILS::ExtractMyValues(*phin,myphin,plm);

  //	  // matrix and vector definition
  //	  LINALG::Matrix<pnsd,pnen>       evelnp;
  std::vector<LINALG::Matrix<pnen,1> > ephinp(numscal_);
  std::vector<LINALG::Matrix<pnen,1> > ephin(numscal_);

  // insert into element arrays
  for (int i=0;i<pnen;++i)
  {
    for (int k = 0; k< numscal_; ++k)
    {
      // split for each tranported scalar, insert into element arrays
      ephinp[k](i,0) = myphinp[k+(i*numdofpernode_)];
      ephin[k](i,0)  = myphin[k+(i*numdofpernode_)];
    }
  }

  //------------------------------------------------------------------------
  // preliminary definitions for integration loop
  //------------------------------------------------------------------------
  // reshape element matrices and vectors and init to zero, construct views
  elemat_epetra.Shape(pnen,pnen);
  elevec_epetra.Size (pnen);
  LINALG::Matrix<pnen,pnen> emat(elemat_epetra.A(),true);
  LINALG::Matrix<pnen,   1> erhs(elevec_epetra.A(),true);

  // (boundary) element local node coordinates
  LINALG::Matrix<pnsd,bnen>  bxyze(true);
  GEO::fillInitialPositionArray<bdistype,pnsd,LINALG::Matrix<pnsd,bnen> >(ele,bxyze);

  // parent element local node coordinates
  LINALG::Matrix<pnsd,pnen>  pxyze(true);
  GEO::fillInitialPositionArray<pdistype,pnsd,LINALG::Matrix<pnsd,pnen> >(pele,pxyze);

  // coordinates of integration points for (boundary) and parent element
  LINALG::Matrix<bnsd,   1>  bxsi(true);
  LINALG::Matrix<pnsd,   1>  pxsi(true);

  // transposed jacobian "dx/ds" and inverse of transposed jacobian "ds/dx"
  // for parent element
  LINALG::Matrix<pnsd,pnsd>  pxjm(true);
  LINALG::Matrix<pnsd,pnsd>  pxji(true);

  // metric tensor for (boundary) element
  LINALG::Matrix<bnsd,bnsd>  bmetrictensor(true);

  // (outward-pointing) unit normal vector to (boundary) element
  LINALG::Matrix<pnsd,   1>  bnormal(true);

  // velocity vector at integration point
  LINALG::Matrix<pnsd,   1>  velint;

  // gradient of scalar value at integration point
  LINALG::Matrix<pnsd,1> gradphi;

  // (boundary) element shape functions, local and global derivatives
  LINALG::Matrix<bnen,   1>  bfunct(true);
  LINALG::Matrix<bnsd,bnen>  bderiv(true);
  LINALG::Matrix<bnsd,bnen>  bderxy(true);

  // parent element shape functions, local and global derivatives
  LINALG::Matrix<pnen,   1>  pfunct(true);
  LINALG::Matrix<pnsd,pnen>  pderiv(true);
  LINALG::Matrix<pnsd,pnen>  pderxy(true);


  // use one-point Gauss rule to do calculations at element center
  DRT::UTILS::IntPointsAndWeights<bnsd> intpoints_tau(SCATRA::DisTypeToStabGaussRule<bdistype>::rule);

  // element surface area (1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double* gpcoord = (intpoints_tau.IP().qxg)[0];
  for (int idim=0;idim<bnsd;idim++)
  {
    bxsi(idim) = gpcoord[idim];
  }
  DRT::UTILS::shape_function_deriv1<bdistype>(bxsi,bderiv);
  double drs = 0.0;
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<bdistype>(bxyze,bderiv,bmetrictensor,drs,&bnormal);

  //------------------------------------------------------------------------
  // preliminary computations for integration loop
  //------------------------------------------------------------------------
  // integrations points and weights for (boundary) element and parent element
  const DRT::UTILS::IntPointsAndWeights<bnsd> bintpoints(SCATRA::DisTypeToOptGaussRule<bdistype>::rule);

  const DRT::UTILS::IntPointsAndWeights<pnsd> pintpoints(SCATRA::DisTypeToOptGaussRule<pdistype>::rule);

  // transfer integration-point coordinates of (boundary) element to parent element
  Epetra_SerialDenseMatrix pqxg(pintpoints.IP().nquad,pnsd);
  {
    Epetra_SerialDenseMatrix gps(bintpoints.IP().nquad,bnsd);

    for (int iquad=0; iquad<bintpoints.IP().nquad; ++iquad)
    {
      const double* gpcoord = (bintpoints.IP().qxg)[iquad];

      for (int idim=0;idim<bnsd ;idim++)
      {
        gps(iquad,idim) = gpcoord[idim];
      }
    }
    if(pnsd==2)
    {
      DRT::UTILS::BoundaryGPToParentGP2(pqxg,gps,pdistype,bdistype,ele->BeleNumber());
    }
    else if (pnsd==3)
    {
      DRT::UTILS::BoundaryGPToParentGP3(pqxg,gps,pdistype,bdistype,ele->BeleNumber());
    }
  }


  const double reinit_pseudo_timestepsize_factor = params.get<double>("pseudotimestepsize_factor");

  const double meshsize = getEleDiameter<pdistype>(pxyze);

  const double pseudo_timestep_size = meshsize * reinit_pseudo_timestepsize_factor;

  //------------------------------------------------------------------------
  // integration loop: boundary integrals
  //------------------------------------------------------------------------
  for (int iquad=0; iquad<bintpoints.IP().nquad; ++iquad)
  {
    // reference coordinates of integration point from (boundary) element
    const double* gpcoord = (bintpoints.IP().qxg)[iquad];
    for (int idim=0;idim<bnsd;idim++)
    {
      bxsi(idim) = gpcoord[idim];
    }

    // (boundary) element shape functions
    DRT::UTILS::shape_function       <bdistype>(bxsi,bfunct);
    DRT::UTILS::shape_function_deriv1<bdistype>(bxsi,bderiv);

    // global coordinates of current integration point from (boundary) element
    LINALG::Matrix<pnsd,1> coordgp(true);
    for (int A=0;A<bnen;++A)
    {
      for(int j=0;j<pnsd;++j)
      {
        coordgp(j)+=bxyze(j,A)*bfunct(A);
      }
    }

    // reference coordinates of integration point from parent element
    for (int idim=0;idim<pnsd;idim++)
    {
      pxsi(idim) = pqxg(iquad,idim);
    }

    // parent element shape functions and local derivatives
    DRT::UTILS::shape_function       <pdistype>(pxsi,pfunct);
    DRT::UTILS::shape_function_deriv1<pdistype>(pxsi,pderiv);

    // Jacobian matrix and determinant of parent element (including check)
    pxjm.MultiplyNT(pderiv,pxyze);
    const double det = pxji.Invert(pxjm);
    if (det < 1E-16) dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", pele->Id(), det);

    // compute measure tensor for surface element, infinitesimal area element drs
    // and (outward-pointing) unit normal vector
    DRT::UTILS::ComputeMetricTensorForBoundaryEle<bdistype>(bxyze,bderiv,bmetrictensor,drs,&bnormal);

    // for nurbs elements the normal vector must be scaled with a special orientation factor!!
    if(DRT::NURBS::IsNurbs(distype))
      bnormal.Scale(normalfac_);

    // compute integration factor
    const double fac_surface = bintpoints.IP().qwgt[iquad]*drs;

    // compute global derivatives
    pderxy.Multiply(pxji,pderiv);

    //--------------------------------------------------------------------
    // loop over scalars (not yet implemented for more than one scalar)
    //--------------------------------------------------------------------
    for(int dofindex=0;dofindex<numdofpernode_;++dofindex)
    {
      //----------  --------------      |                    |
      //  mat              -1/4* dtau^2 | w, n*grad(D(psi) ) |
      //--------------------------      |                    |

      LINALG::Matrix<1,pnen> derxy_normal;
      derxy_normal.Clear();
      derxy_normal.MultiplyTN(bnormal,pderxy);

      for (int vi=0; vi<pnen; ++vi)
      {
        const int fvi = vi*numdofpernode_+dofindex;

        for (int ui=0; ui<pnen; ++ui)
        {
          const int fui = ui*numdofpernode_+dofindex;

          emat(fvi,fui) -= pfunct(vi)* (fac_surface*pseudo_timestep_size*pseudo_timestep_size/4.0) * derxy_normal(0,ui);
        }
      }

      //----------  --------------      |              m     |
      //  rhs               0.5* dtau^2 | w, n*grad(psi )    |
      //--------------------------      |                    |

      // update grad_dist_n
      LINALG::Matrix<pnsd,1> grad_dist_n(true);
      grad_dist_n.Multiply(pderxy,ephin[dofindex]);

      LINALG::Matrix<1,1> grad_dist_n_normal(true);
      grad_dist_n_normal.MultiplyTN(bnormal,grad_dist_n);

      for (int vi=0; vi<pnen; ++vi)
      {
        const int fvi = vi*numdofpernode_+dofindex;

        erhs(fvi) += pfunct(vi)*pseudo_timestep_size*pseudo_timestep_size*fac_surface/2.0 * grad_dist_n_normal(0,0);
      }


      //                    |              m+1     m  |
      //    1/4*delta_tau^2 | w, n*grad(psi   - psi ) |
      //                    |              i          |
      // update grad_dist_n
      LINALG::Matrix<pnsd,1> grad_dist_npi(true);
      grad_dist_npi.Multiply(pderxy,ephinp[dofindex]);

      LINALG::Matrix<1,1> grad_dist_npi_normal;
      grad_dist_npi_normal.Clear();
      grad_dist_npi_normal.MultiplyTN(bnormal,grad_dist_npi);

      double Grad_Dpsi_normal = grad_dist_npi_normal(0,0) - grad_dist_n_normal(0,0);


      for (int vi=0; vi<pnen; ++vi)
      {
        const int fvi = vi*numdofpernode_+dofindex;

        erhs(fvi) += pfunct(vi)*Grad_Dpsi_normal*fac_surface*pseudo_timestep_size*pseudo_timestep_size/4.0;
      }

    } // loop over scalars
  } // loop over integration points

  return;
}

