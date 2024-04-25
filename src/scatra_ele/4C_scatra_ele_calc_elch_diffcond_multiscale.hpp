/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra elements for diffusion-conduction ion-transport equations within a
multi-scale framework

\level 2

*/
/*--------------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_CALC_ELCH_DIFFCOND_MULTISCALE_HPP
#define FOUR_C_SCATRA_ELE_CALC_ELCH_DIFFCOND_MULTISCALE_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_calc_elch_diffcond.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    class ScaTraEleDiffManagerElchDiffCondMultiScale;

    // implementation of class ScaTraEleCalcElchDiffCondMultiScale
    template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype>>
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
          CORE::LINALG::SerialDenseMatrix& emat,                   //!< element matrix
          CORE::LINALG::SerialDenseVector& erhs,  //!< element right-hand side vector
          const int k,                            //!< species index
          const int iquad,                        //!< Gauss point index
          const double timefacfac,  //!< domain integration factor times time integration factor
          const double rhsfac       //!< domain integration factor times time integration factor for
                                    //!< right-hand side vector
          ) override;

      //! calculate electrode state of charge and C rate
      void CalculateElectrodeSOCAndCRate(
          const DRT::Element* const& ele,             //!< the element we are dealing with
          const DRT::Discretization& discretization,  //!< discretization
          DRT::Element::LocationArray& la,            //!< location array
          CORE::LINALG::SerialDenseVector&
              scalars  //!< result vector for scalar integrals to be computed
          ) final;

      void CalculateMeanElectrodeConcentration(const DRT::Element* const& ele,
          const DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          CORE::LINALG::SerialDenseVector& conc) override;

      void CalculateScalars(const DRT::Element* ele, CORE::LINALG::SerialDenseVector& scalars,
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
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) override;

      //! compute element matrix and element right-hand side vector
      void Sysmat(DRT::Element* ele,                  //!< element
          CORE::LINALG::SerialDenseMatrix& emat,      //!< element matrix
          CORE::LINALG::SerialDenseVector& erhs,      //!< element right-hand side vector
          CORE::LINALG::SerialDenseVector& subgrdiff  //!< subgrid diffusivity vector
          ) override;
    };  // class implementation


    // implementation of class ScaTraEleDiffManagerElchDiffCondMultiScale
    class ScaTraEleDiffManagerElchDiffCondMultiScale : public ScaTraEleDiffManagerElchDiffCond
    {
     public:
      //! constructor
      ScaTraEleDiffManagerElchDiffCondMultiScale(int numscal)
          : ScaTraEleDiffManagerElchDiffCond(numscal){};

      //! Output of transport parameter (to screen)
      void OutputTransportParams(const int numscal) override
      {
        // call base class routine
        ScaTraEleDiffManagerElchDiffCond::OutputTransportParams(numscal);
      };
    };
  }  // namespace ELEMENTS
}  // namespace DRT
FOUR_C_NAMESPACE_CLOSE

#endif
