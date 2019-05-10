/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch_diffcond_multiscale.H

\brief evaluation of ScaTra elements for diffusion-conduction ion-transport equations within a
multi-scale framework

\level 2

<pre>
\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef SCATRA_ELE_CALC_ELCH_DIFFCOND_MULTISCALE_H
#define SCATRA_ELE_CALC_ELCH_DIFFCOND_MULTISCALE_H

#include "scatra_ele_calc_elch_diffcond.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    class ScaTraEleDiffManagerElchDiffCondMultiScale;

    // implementation of class ScaTraEleCalcElchDiffCondMultiScale
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleCalcElchDiffCondMultiScale : public ScaTraEleCalcElchDiffCond<distype>
    {
     public:
      //! singleton access method
      static ScaTraEleCalcElchDiffCondMultiScale<distype>* Instance(const int numdofpernode,
          const int numscal, const std::string& disname,
          const ScaTraEleCalcElchDiffCondMultiScale* delete_me = NULL);

      //! singleton destruction
      virtual void Done();

     private:
      //! abbreviation
      typedef ScaTraEleCalc<distype> my;
      typedef ScaTraEleCalcElchDiffCond<distype> mydiffcond;

      //! constructor
      ScaTraEleCalcElchDiffCondMultiScale(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! destructor
      virtual ~ScaTraEleCalcElchDiffCondMultiScale() { return; };

      //! macro-scale matrix and vector contributions arising from macro-micro coupling in
      //! multi-scale simulations
      void CalcMatAndRhsMultiScale(const DRT::Element* const ele,  //!< element
          Epetra_SerialDenseMatrix& emat,                          //!< element matrix
          Epetra_SerialDenseVector& erhs,  //!< element right-hand side vector
          const int k,                     //!< species index
          const int iquad,                 //!< Gauss point index
          const double timefacfac,  //!< domain integration factor times time integration factor
          const double rhsfac       //!< domain integration factor times time integration factor for
                                    //!< right-hand side vector
      );

      //! calculate electrode state of charge and C rate
      void CalculateElectrodeSOCAndCRate(
          const DRT::Element* const& ele,             //!< the element we are dealing with
          const DRT::Discretization& discretization,  //!< discretization
          DRT::Element::LocationArray& la,            //!< location array
          Epetra_SerialDenseVector& scalars  //!< result vector for scalar integrals to be computed
          ) final;

      //! get diffusion manager
      Teuchos::RCP<ScaTraEleDiffManagerElchDiffCondMultiScale> DiffManager() const
      {
        return Teuchos::rcp_static_cast<ScaTraEleDiffManagerElchDiffCondMultiScale>(
            my::diffmanager_);
      };

      //! evaluate action
      int EvaluateAction(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const SCATRA::Action& action,
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra);

      //! compute element matrix and element right-hand side vector
      void Sysmat(DRT::Element* ele,           //!< element
          Epetra_SerialDenseMatrix& emat,      //!< element matrix
          Epetra_SerialDenseVector& erhs,      //!< element right-hand side vector
          Epetra_SerialDenseVector& subgrdiff  //!< subgrid diffusivity vector
      );
    };  // class implementation


    // implementation of class ScaTraEleDiffManagerElchDiffCondMultiScale
    class ScaTraEleDiffManagerElchDiffCondMultiScale : public ScaTraEleDiffManagerElchDiffCond
    {
     public:
      //! constructor
      ScaTraEleDiffManagerElchDiffCondMultiScale(int numscal)
          : ScaTraEleDiffManagerElchDiffCond(numscal), sigma_(0.)
      {
        return;
      };

      //! destructor
      virtual ~ScaTraEleDiffManagerElchDiffCondMultiScale() { return; };

      //! set electronic conductivity
      void SetSigma(const double sigma)
      {
        sigma_ = sigma;
        return;
      };

      //! access electronic conductivity
      double GetSigma() const { return sigma_; };

      //! Output of transport parameter (to screen)
      void OutputTransportParams(const int numscal)
      {
        // call base class routine
        ScaTraEleDiffManagerElchDiffCond::OutputTransportParams(numscal);

        // additional output
        std::cout << "electronic conductivity:   " << sigma_ << std::endl;

        return;
      };

     private:
      //! electrolyte conductivity
      double sigma_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
#endif
