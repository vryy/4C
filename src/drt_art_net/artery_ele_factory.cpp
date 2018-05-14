/*--------------------------------------------------------------------------*/
/*!
\file artery_ele_factory.cpp

\brief Factory of artery elements

<pre>
\level 3

\maintainer Johannes Kremheller
            kremheller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
</pre>
*/
/*--------------------------------------------------------------------------*/
#include "artery_ele_factory.H"
#include "artery_ele_calc_lin_exp.H"
#include "artery_ele_calc_pres_based.H"
#include "artery_ele_interface.H"

/*--------------------------------------------------------------------------*
 | (public) kremheller                                                03/18 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::ArteryEleInterface* DRT::ELEMENTS::ArtNetFactory::ProvideImpl(
  DRT::Element::DiscretizationType distype,
  INPAR::ARTDYN::ImplType problem,
  const std::string& disname
  )
{
  switch(distype)
  {
  case DRT::Element::line2:
  {
    return DefineProblemType<DRT::Element::line2>(problem,disname);

    break;
  }
  default:
    dserror("Only line2 elements available so far");
    break;
  }
  return NULL;

}


/*--------------------------------------------------------------------------*
 | (public) kremheller                                                03/18 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ArteryEleInterface* DRT::ELEMENTS::ArtNetFactory::DefineProblemType(
    INPAR::ARTDYN::ImplType problem,
    const std::string& disname)
{

  switch(problem)
  {
  case INPAR::ARTDYN::ImplType::impltype_lin_exp:
  {
    // 2 dofs per node
    return DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::Instance(2,disname);
    break;
  }
  case INPAR::ARTDYN::ImplType::impltype_pressure_based:
  {
    // 1 dof per node (only pressure)
    return DRT::ELEMENTS::ArteryEleCalcPresBased<distype>::Instance(1,disname);
    break;
  }
  default:
  {
    dserror("Defined problem type %d does not exist!!", problem);
    break;
  }
  }

  return NULL;
}
