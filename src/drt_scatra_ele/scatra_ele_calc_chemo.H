/*----------------------------------------------------------------------*/
/*! \file
 \brief main file containing routines for calculation of scatra element with chemotactic terms

\level 2

 *----------------------------------------------------------------------*/


#ifndef SCATRA_ELE_CALC_CHEMO_H_
#define SCATRA_ELE_CALC_CHEMO_H_

#include "scatra_ele_calc.H"
#include "../drt_mat/scatra_chemotaxis_mat.H"

namespace DRT
{
  namespace ELEMENTS
  {
    template <DRT::Element::DiscretizationType distype,
        int probdim = DRT::UTILS::DisTypeToDim<distype>::dim>
    class ScaTraEleCalcChemo : public virtual ScaTraEleCalc<distype, probdim>
    {
     protected:
      //! (private) protected constructor, since we are a Singleton.
      ScaTraEleCalcChemo(const int numdofpernode, const int numscal, const std::string& disname);

     private:
      typedef ScaTraEleCalc<distype, probdim> my;
      typedef ScaTraEleInternalVariableManager<my::nsd_, my::nen_> varmanager;

     public:
      virtual ~ScaTraEleCalcChemo() {}

      //! Singleton access method
      static ScaTraEleCalcChemo<distype, probdim>* Instance(const int numdofpernode,
          const int numscal, const std::string& disname,
          const ScaTraEleCalcChemo* delete_me = NULL);

      //! called upon destruction
      virtual void Done();

     protected:
      //! calculation of chemotactic element matrix
      virtual void CalcMatChemo(Epetra_SerialDenseMatrix& emat, const int k,
          const double timefacfac, const double timetaufac, const double densnp,
          const double scatrares, const LINALG::Matrix<my::nen_, 1>& sgconv,
          const LINALG::Matrix<my::nen_, 1>& diff);

      //! calculation of chemotactic RHS
      virtual void CalcRHSChemo(Epetra_SerialDenseVector& erhs, const int k, const double rhsfac,
          const double rhstaufac, const double scatrares, const double densnp);

      //! get the material parameters
      virtual void GetMaterialParams(const DRT::Element* ele,  //!< the element we are dealing with
          std::vector<double>& densn,                          //!< density at t_(n)
          std::vector<double>& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,  //!< density at t_(n+alpha_M)
          double& visc,                 //!< fluid viscosity
          const int iquad               //!< id of current gauss point
      );

      //! Clear all chemotaxtis related class variable (i.e. set them to zero)
      void ClearChemotaxisTerms();

      //! Get and save all chemotaxtis related class variable
      virtual void GetChemotaxisCoefficients(
          const Teuchos::RCP<const MAT::Material> material  //!< pointer to current material
      );

      //! Get ID of attractant (i.e. scalar which gradient the current scalar shall follow)
      virtual int GetPartner(const std::vector<int> pair);

      //! calculation of strong residual for stabilization
      virtual void CalcStrongResidual(const int k,  //!< index of current scalar
          double& scatrares,     //!< residual of convection-diffusion-reaction eq
          const double densam,   //!< density at t_(n+am)
          const double densnp,   //!< density at t_(n+1)
          const double rea_phi,  //!< reactive contribution
          const double rhsint,   //!< rhs at integration point
          const double tau       //!< the stabilisation parameter
      );


      int numcondchemo_;                    //!< number of chemotactic conditions
      std::vector<std::vector<int>> pair_;  //!< vector containing the pairings
      std::vector<double> chemocoeff_;  //!< constants by which the chemotactic terms are multiplied

    };  // end Chemotaxis

  }  // namespace ELEMENTS

}  // namespace DRT


#endif /* SCATRA_ELE_CALC_CHEMO_H_ */
