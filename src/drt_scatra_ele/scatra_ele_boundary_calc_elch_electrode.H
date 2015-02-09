/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_elch_electrode.H

\brief evaluation of ScaTra boundary elements for electrodes

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_H
#define SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_H

#include "scatra_ele_boundary_calc_elch.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    template<DRT::Element::DiscretizationType distype> class ScaTraEleBoundaryCalcElchElectrode : public ScaTraEleBoundaryCalcElch<distype>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype> my;
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype> myelch;

    public:

      //! singleton access method
      static ScaTraEleBoundaryCalcElchElectrode<distype>* Instance(
          const int numdofpernode,
          const int numscal,
          bool create=true
          );

      //! destructor
      virtual ~ScaTraEleBoundaryCalcElchElectrode(){return;};

      //! singleton destruction
      virtual void Done();

    protected:

      //! protected constructor for singletons
      ScaTraEleBoundaryCalcElchElectrode(const int numdofpernode,const int numscal);

      //! extract valence of species k from element material
      virtual const double GetValence(
          const Teuchos::RCP<const MAT::Material>&   material,   //! element material
          const int                                  k           //! species number
      ) const;
    }; // class ScaTraEleBoundaryCalcElchElectrode
  } // namespace ELEMENTS
} // namespace DRT
#endif
