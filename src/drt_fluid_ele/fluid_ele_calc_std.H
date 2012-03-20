/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_std.H

\brief standard routines for calculation of fluid element

<pre>
Maintainer: Volker Gravemeier & Andreas Ehrl
            {vgravem,ehrl}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15245/15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_STD_H
#define FLUID_ELE_CALC_STD_H

#include "fluid_ele_calc.H"

namespace DRT
{
namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype>
class FluidEleCalcStd : public FluidEleCalc<distype>
{
    /// private constructor, since we are a Singleton.
    FluidEleCalcStd();

    typedef FluidEleCalc<distype> my;

public:
    virtual ~FluidEleCalcStd() {}

    /// Singleton access method
    static FluidEleCalcStd<distype> * Instance( bool create=true );

    /// called upon destruction
    virtual void Done();

    virtual int CalcDissipation(Fluid3*                    ele,
                                ParameterList&             params,
                                DRT::Discretization&       discretization,
                                vector<int>&               lm,
                                RefCountPtr<MAT::Material> mat);

    /*!
      \brief finite difference check for debugging
    */

    void FDcheck(
      int                                                                 eid,
      const LINALG::Matrix<my::nsd_,my::nen_>&                            evelaf,
      const LINALG::Matrix<my::nsd_,my::nen_>&                            eveln,
      const LINALG::Matrix<my::nsd_,my::nen_>&                            fsevelaf,
      const LINALG::Matrix<my::nen_,1>&                                   epreaf,
      const LINALG::Matrix<my::nsd_,my::nen_>&                            eaccam,
      const LINALG::Matrix<my::nen_,1>&                                   escaaf,
      const LINALG::Matrix<my::nen_,1>&                                   escaam,
      const LINALG::Matrix<my::nen_,1>&                                   escadtam,
      const LINALG::Matrix<my::nsd_,my::nen_>&                            emhist,
      const LINALG::Matrix<my::nsd_,my::nen_>&                            edispnp,
      const LINALG::Matrix<my::nsd_,my::nen_>&                            egridv,
      const LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_>&  estif,
      const LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_>&  emesh,
      const LINALG::Matrix<(my::nsd_+1)*my::nen_,    1>&                  eforce,
      const double                                                        thermpressaf,
      const double                                                        thermpressam,
      const double                                                        thermpressdtaf,
      const double                                                        thermpressdtam,
      const Teuchos::RCP<const MAT::Material>                             material,
      const double                                                        timefac,
      const double&                                                       Cs,
      const double&                                                       Cs_delta_sq,
      const double&                                                       l_tau);



private:

};
}
}

#endif

