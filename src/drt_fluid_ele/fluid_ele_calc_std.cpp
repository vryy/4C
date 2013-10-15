/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_std.cpp

\brief standard routines for calculation of fluid element

<pre>
Maintainer: Ursula Rasthofer & Volker Gravemeier
            {rasthofer,vgravem}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236/-245
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_ele_calc_std.H"
#include "fluid_ele_parameter.H"
#include "../drt_inpar/inpar_fpsi.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcStd<distype>* DRT::ELEMENTS::FluidEleCalcStd<distype>::Instance( bool create, int num )
{
  if ( create )
   {
     if(static_cast<int>(MY::instances_.count(num))==0)
     {
         MY::instances_.insert(std::pair<int,std::map<int,MY* >* >(num, new std::map<int,MY* >));
         MY::instances_.at(num)->insert(std::pair<int,MY* >((int)distype,new FluidEleCalcStd<distype>(num)));
     }
     else if ( MY::instances_.count(num) > 0 and MY::instances_.at(num)->count((int)distype) == 0 )
     {
       MY::instances_.at(num)->insert(std::pair<int,MY* >((int)distype, new FluidEleCalcStd<distype>(num)));
     }

     return static_cast<DRT::ELEMENTS::FluidEleCalcStd<distype>* >(MY::instances_.at(num)->at((int)distype));
   }
   else
   {
     if ( MY::instances_.at(num)->size())
     {
       delete MY::instances_.at(num)->at((int)distype);
       MY::instances_.at(num)->erase((int)distype);

       if ( !(MY::instances_.at(num)->size()) )
         MY::instances_.erase(num);
     }

     return NULL;
   }

  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcStd<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
    Instance( false, num_ );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcStd<distype>::FluidEleCalcStd(int num)
  : DRT::ELEMENTS::FluidEleCalc<distype>::FluidEleCalc(num),
    num_(num)
{

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/

// template classes
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::hex8>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::hex20>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::hex27>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::tet4>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::tet10>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::wedge6>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::quad4>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::quad8>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::quad9>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::tri3>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::tri6>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::nurbs27>;
