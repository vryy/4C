/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_ls.cpp

\brief evaluations for level sets

\level 2

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*--------------------------------------------------------------------------*/
#include "scatra_ele_calc_ls.H"
#include "scatra_ele_parameter_std.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcLS<distype>* DRT::ELEMENTS::ScaTraEleCalcLS<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname,
    const ScaTraEleCalcLS* delete_me)
{
  static std::map<std::string, ScaTraEleCalcLS<distype>*> instances;

  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleCalcLS<distype>(numdofpernode, numscal, disname);
  }

  else
  {
    for (typename std::map<std::string, ScaTraEleCalcLS<distype>*>::iterator i = instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLS<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(0, 0, "", this);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcLS<distype>::ScaTraEleCalcLS(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode, numscal, disname)
{
  // safety check
  if (my::scatrapara_->RBSubGrVel()) dserror("CalcSubgrVelocityLevelSet not available anymore");

  return;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::quad4>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::nurbs27>;
