/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra elements for diffusion-conduction ion-transport equations within a
multi-scale framework

\level 2

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
    template <DRT::Element::DiscretizationType distype,
        int probdim = CORE::DRT::UTILS::DisTypeToDim<distype>::dim>
    class ScaTraEleCalcElchDiffCondMultiScale : public ScaTraEleCalcElchDiffCond<distype, probdim>
    {
     public:
      //! singleton access method
      static ScaTraEleCalcElchDiffCondMultiScale<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);



     private:
      //! abbreviation
      using my = ScaTraEleCalc<distype, probdim>;
      using mydiffcond = ScaTraEleCalcElchDiffCond<distype, probdim>;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

      //! constructor
      ScaTraEleCalcElchDiffCondMultiScale(
          const int numdofpernode, const int numscal, const std::string& disname);

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

      void CalculateMeanElectrodeConcentration(const DRT::Element* const& ele,
          const DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseVector& conc) override;

      void CalculateScalars(const DRT::Element* ele, Epetra_SerialDenseVector& scalars,
          bool inverting, bool calc_grad_phi) override;

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
          : ScaTraEleDiffManagerElchDiffCond(numscal), sigma_(0.0){};

      //! set electronic conductivity
      void SetSigma(const double sigma) { sigma_ = sigma; };

      //! access electronic conductivity
      double GetSigma() const { return sigma_; };

      //! Output of transport parameter (to screen)
      void OutputTransportParams(const int numscal) override
      {
        // call base class routine
        ScaTraEleDiffManagerElchDiffCond::OutputTransportParams(numscal);

        // additional output
        std::cout << "electronic conductivity:   " << sigma_ << std::endl;
      };

     private:
      //! electrolyte conductivity
      double sigma_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
#endif
