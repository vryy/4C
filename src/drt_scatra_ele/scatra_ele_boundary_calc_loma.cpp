/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_loma.cpp

\brief evaluation of ScaTra boundary elements for low Mach number problems

\level 2

<pre>
\maintainer Anh-Tu Vuong
            vuong@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15237
</pre>
 */
/*----------------------------------------------------------------------*/
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/arrhenius_pv.H"
#include "../drt_mat/arrhenius_temp.H"
#include "../drt_mat/ferech_pv.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/mixfrac.H"
#include "../drt_mat/sutherland.H"
#include "../drt_mat/tempdepwater.H"
#include "../drt_mat/thermostvenantkirchhoff.H"
#include "../drt_mat/yoghurt.H"

#include "scatra_ele.H"
#include "scatra_ele_boundary_calc_loma.H"

#include "../drt_fluid/fluid_rotsym_periodicbc.H"


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype>*
DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype>::Instance(const int numdofpernode,
    const int numscal, const std::string& disname, const ScaTraEleBoundaryCalcLoma* delete_me)
{
  static std::map<std::string, ScaTraEleBoundaryCalcLoma<distype>*> instances;

  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleBoundaryCalcLoma<distype>(numdofpernode, numscal, disname);
  }

  else
  {
    for (typename std::map<std::string, ScaTraEleBoundaryCalcLoma<distype>*>::iterator i =
             instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype>::Done()
{
  // delete singleton
  Instance(0, 0, "", this);

  return;
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype>::ScaTraEleBoundaryCalcLoma(
    const int numdofpernode, const int numscal, const std::string& disname)
    :  // constructor of base class
      my::ScaTraEleBoundaryCalc(numdofpernode, numscal, disname)
{
  return;
}


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype>::EvaluateAction(DRT::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    SCATRA::BoundaryAction action, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  // determine and evaluate action
  switch (action)
  {
    case SCATRA::bd_calc_loma_therm_press:
    {
      CalcLomaThermPress(ele, params, discretization, la);

      break;
    }

    default:
    {
      my::EvaluateAction(ele, params, discretization, action, la, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);

      break;
    }
  }  // switch (action)

  return 0;
}


/*----------------------------------------------------------------------*
 | calculate loma therm pressure                              vg 03/09  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype>::CalcLomaThermPress(DRT::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la)
{
  // get location vector associated with primary dofset
  std::vector<int>& lm = la[0].lm_;

  DRT::Element* parentele = ele->ParentElement();
  // we dont know the parent element's lm vector; so we have to build it here
  const int nenparent = parentele->NumNode();
  std::vector<int> lmparent(nenparent);
  std::vector<int> lmparentowner;
  std::vector<int> lmparentstride;
  parentele->LocationVector(discretization, lmparent, lmparentowner, lmparentstride);

  // get number of dofset associated with velocity related dofs
  const int ndsvel = params.get<int>("ndsvel");

  // get velocity values at nodes
  const Teuchos::RCP<const Epetra_Vector> convel =
      discretization.GetState(ndsvel, "convective velocity field");

  // safety check
  if (convel == Teuchos::null) dserror("Cannot get state vector convective velocity");

  // get values of velocity field from secondary dof-set
  const std::vector<int>& lmvel = la[ndsvel].lm_;
  std::vector<double> myconvel(lmvel.size());

  // extract local values of the global vectors
  DRT::UTILS::ExtractMyValues(*convel, myconvel, lmvel);

  // rotate the vector field in the case of rotationally symmetric boundary conditions
  my::rotsymmpbc_->RotateMyValuesIfNecessary(myconvel);

  // define vector for normal diffusive and velocity fluxes
  std::vector<double> mynormdiffflux(lm.size());
  std::vector<double> mynormvel(lm.size());

  // determine constant outer normal to this element
  my::GetConstNormal(my::normal_, my::xyze_);

  // extract temperature flux vector for each node of the parent element
  LINALG::SerialDenseMatrix eflux(3, nenparent);
  DRT::Element* peleptr = (DRT::Element*)parentele;
  int k = my::numscal_ - 1;  // temperature is always last degree of freedom!!
  std::ostringstream temp;
  temp << k;
  std::string name = "flux_phi_" + temp.str();
  // try to get the pointer to the entry (and check if type is Teuchos::RCP<Epetra_MultiVector>)
  Teuchos::RCP<Epetra_MultiVector>* f = params.getPtr<Teuchos::RCP<Epetra_MultiVector>>(name);
  // check: field has been set and is not of type Teuchos::null
  if (f != NULL)
    DRT::UTILS::ExtractMyNodeBasedValues(peleptr, eflux, *f, 3);
  else
    dserror("MultiVector %s has not been found!", name.c_str());

  // calculate normal diffusive and velocity flux at each node of the
  // present boundary element
  for (int i = 0; i < my::nen_; ++i)
  {
    for (int j = 0; j < nenparent; ++j)
    {
      mynormdiffflux[i] = 0.0;
      mynormvel[i] = 0.0;
      for (int l = 0; l < my::nsd_ + 1; l++)
      {
        mynormdiffflux[i] += eflux(l, j) * my::normal_(l);
        mynormvel[i] += myconvel[i * (my::nsd_ + 1) + l] * my::normal_(l);
      }
    }
  }

  // calculate integral of normal diffusive and velocity flux
  // NOTE: add integral value only for elements which are NOT ghosted!
  if (ele->Owner() == discretization.Comm().MyPID())
  {
    NormDiffFluxAndVelIntegral(ele, params, mynormdiffflux, mynormvel);
  }
}


/*----------------------------------------------------------------------*
 | calculate Neumann inflow boundary conditions              fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype>::NeumannInflow(const DRT::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& emat, Epetra_SerialDenseVector& erhs)
{
  // set thermodynamic pressure
  thermpress_ = params.get<double>("thermodynamic pressure");

  // call base class routine
  my::NeumannInflow(ele, params, discretization, la, emat, erhs);

  return;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype>::NeumannInflow


/*----------------------------------------------------------------------*
 | get density at integration point                          fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype>::GetDensity(
    Teuchos::RCP<const MAT::Material> material,
    const std::vector<LINALG::Matrix<my::nen_, 1>>& ephinp, const int k)
{
  // initialization
  double density(0.);

  // get density depending on material
  switch (material->MaterialType())
  {
    case INPAR::MAT::m_matlist:
    {
      const MAT::MatList* actmat = static_cast<const MAT::MatList*>(material.get());

      const int lastmatid = actmat->NumMat() - 1;

      if (actmat->MaterialById(lastmatid)->MaterialType() == INPAR::MAT::m_arrhenius_temp)
      {
        // compute temperature and check whether it is positive
        const double temp = my::funct_.Dot(ephinp[my::numscal_ - 1]);
        if (temp < 0.0)
          dserror(
              "Negative temperature in ScaTra Arrhenius temperature density evaluation on "
              "boundary!");

        // compute density based on temperature and thermodynamic pressure
        density = static_cast<const MAT::ArrheniusTemp*>(actmat->MaterialById(lastmatid).get())
                      ->ComputeDensity(temp, thermpress_);
      }
      else
        dserror(
            "Type of material found in material list not supported, should be Arrhenius-type "
            "temperature!");

      break;
    }

    case INPAR::MAT::m_mixfrac:
    {
      // compute density based on mixture fraction
      density = static_cast<const MAT::MixFrac*>(material.get())
                    ->ComputeDensity(my::funct_.Dot(ephinp[k]));

      break;
    }

    case INPAR::MAT::m_sutherland:
    {
      // compute temperature and check whether it is positive
      const double temp = my::funct_.Dot(ephinp[k]);
      if (temp < 0.0)
        dserror("Negative temperature in ScaTra Sutherland density evaluation on boundary!");

      // compute density based on temperature and thermodynamic pressure
      density =
          static_cast<const MAT::Sutherland*>(material.get())->ComputeDensity(temp, thermpress_);

      break;
    }

    case INPAR::MAT::m_tempdepwater:
    {
      // compute temperature and check whether it is positive
      const double temp = my::funct_.Dot(ephinp[k]);
      if (temp < 0.0)
        dserror(
            "Negative temperature in ScaTra temperature-dependent water density evaluation on "
            "boundary!");

      // compute density based on temperature
      density = static_cast<const MAT::TempDepWater*>(material.get())->ComputeDensity(temp);

      break;
    }

    case INPAR::MAT::m_arrhenius_pv:
    {
      // compute density based on progress variable
      density = static_cast<const MAT::ArrheniusPV*>(material.get())
                    ->ComputeDensity(my::funct_.Dot(ephinp[k]));

      break;
    }

    case INPAR::MAT::m_ferech_pv:
    {
      // compute density based on progress variable
      density = static_cast<const MAT::FerEchPV*>(material.get())
                    ->ComputeDensity(my::funct_.Dot(ephinp[k]));

      break;
    }

    case INPAR::MAT::m_thermostvenant:
    {
      // get constant density
      density = static_cast<const MAT::ThermoStVenantKirchhoff*>(material.get())->Density();

      break;
    }

    case INPAR::MAT::m_yoghurt:
    {
      // get constant density
      density = static_cast<const MAT::Yoghurt*>(material.get())->Density();

      break;
    }

    default:
    {
      dserror("Invalid material type!");
      break;
    }
  }

  return density;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype>::GetDensity


/*----------------------------------------------------------------------*
 | calculate integral of normal diffusive flux and velocity     vg 09/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype>::NormDiffFluxAndVelIntegral(
    const DRT::Element* ele, Teuchos::ParameterList& params,
    const std::vector<double>& enormdiffflux, const std::vector<double>& enormvel)
{
  // get variables for integrals of normal diffusive flux and velocity
  double normdifffluxint = params.get<double>("normal diffusive flux integral");
  double normvelint = params.get<double>("normal velocity integral");

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    const double fac = my::EvalShapeFuncAndIntFac(intpoints, gpid);

    // compute integral of normal flux
    for (int node = 0; node < my::nen_; ++node)
    {
      normdifffluxint += my::funct_(node) * enormdiffflux[node] * fac;
      normvelint += my::funct_(node) * enormvel[node] * fac;
    }
  }  // loop over integration points

  // add contributions to the global values
  params.set<double>("normal diffusive flux integral", normdifffluxint);
  params.set<double>("normal velocity integral", normvelint);

  return;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype>::NormDiffFluxAndVelIntegral


// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<DRT::Element::nurbs9>;
