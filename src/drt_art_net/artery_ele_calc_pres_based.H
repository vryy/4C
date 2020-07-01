/*----------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of PressureBased artery element


\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef SRC_DRT_ART_NET_ARTERY_ELE_CALC_PRES_BASED_H_
#define SRC_DRT_ART_NET_ARTERY_ELE_CALC_PRES_BASED_H_

#include "artery.H"
#include "artery_ele_action.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "artery_ele_calc.H"


namespace DRT
{
  namespace ELEMENTS
  {
    /// Internal artery implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the artery element. Additionally the method Sysmat()
      provides a clean and fast element implementation.

      <h3>Purpose</h3>

      \author kremheller
      \date 03/18
    */

    template <DRT::Element::DiscretizationType distype>
    class ArteryEleCalcPresBased : public ArteryEleCalc<distype>
    {
     private:
      typedef ArteryEleCalc<distype> my;

      /// private constructor, since we are a Singleton.
      ArteryEleCalcPresBased(const int numdofpernode, const std::string& disname);

     public:
      /// virtual destructor
      virtual ~ArteryEleCalcPresBased(){};

      //! Singleton access method
      static ArteryEleCalcPresBased<distype>* Instance(const int numdofpernode,
          const std::string& disname, const ArteryEleCalcPresBased* delete_me = NULL);

      /// called upon destruction
      virtual void Done();

      int Evaluate(Artery* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<MAT::Material> mat) override;

      int ScatraEvaluate(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat) override;

      int EvaluateService(Artery* ele, const ARTERY::Action action, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat) override;


     protected:
      /*!
        \brief calculate element matrix and rhs

        \param ele[in]            the element which is evaluated
        \param discretization[in] discretization to which element belongs
        \param la[in]             element location array
        \param sysmat[in,out]     element matrix to calculate
        \param rhs[in,out]        element rhs to calculate
        \param material[in]       artery material
        */
      void Sysmat(Artery* ele, DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          LINALG::Matrix<my::iel_, my::iel_>& sysmat, LINALG::Matrix<my::iel_, 1>& rhs,
          Teuchos::RCP<const MAT::Material> material);

      /*!
        \brief Evaluate volumetric flow inside the element (for post-processing)

        \param ele[in]            the element which is evaluated
        \param discretization[in] discretization to which element belongs
        \param la[in]             element location array
        \param flowVec[in,out]    element flow to calculate
        \param material[in]       artery material/dimesion

        \note  only checked for line2 elements
        */
      void EvaluateFlow(Artery* ele, DRT::Discretization& discretization,
          DRT::Element::LocationArray& la, Epetra_SerialDenseVector& flowVec,
          Teuchos::RCP<const MAT::Material> material);

      /*!
        \brief Calculate element length (either in current or deformed configuration)
        \param ele[in]            the element which is evaluated
        \param discretization[in] discretization to which element belongs
        \param la[in]             element location array
        \return                   element length (either in current or deformed configuration)
       */
      double CalculateEleLength(
          Artery* ele, DRT::Discretization& discretization, DRT::Element::LocationArray& la);
    };

  }  // namespace ELEMENTS

}  // namespace DRT



#endif /* SRC_DRT_ART_NET_ARTERY_ELE_CALC_PRES_BASED_H_ */
