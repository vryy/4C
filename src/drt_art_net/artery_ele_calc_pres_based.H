/*----------------------------------------------------------------------*/
/*!
\file artery_ele_calc_pres_based.H

\brief Internal implementation of PressureBased artery element

\maintainer Johannes Kremheller
            kremheller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/

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
  class ArteryEleCalcPresBased: public ArteryEleCalc<distype>
  {
  private:

    typedef ArteryEleCalc<distype> my;

    /// private constructor, since we are a Singleton.
    ArteryEleCalcPresBased(const int numdofpernode, const std::string& disname);

  public:

    /// virtual destructor
    virtual ~ArteryEleCalcPresBased() {};

    //! Singleton access method
    static ArteryEleCalcPresBased<distype> * Instance(
        const int numdofpernode,
        const std::string& disname,
        const ArteryEleCalcPresBased* delete_me=NULL);

    /// called upon destruction
    virtual void Done();

    /// Evaluate
    /*!
      The evaluate function for the general artery case.
     */
    virtual int Evaluate(Artery*                   ele,
                         Teuchos::ParameterList&   params,
                         DRT::Discretization&      discretization,
                         std::vector<int>&         lm,
                         Epetra_SerialDenseMatrix& elemat1_epetra,
                         Epetra_SerialDenseMatrix& elemat2_epetra,
                         Epetra_SerialDenseVector& elevec1_epetra,
                         Epetra_SerialDenseVector& elevec2_epetra,
                         Epetra_SerialDenseVector& elevec3_epetra,
                         Teuchos::RCP<MAT::Material> mat);

    /// ScatraEvaluate
    /*!
      Evaluate scalar transport in arteries, this is not available for the pressure-based formulation
      Instead, ScaTra-discretization is cloned and work is done by ScaTraEleCalcArtery
     */
    virtual int ScatraEvaluate(Artery*                   ele,
                         Teuchos::ParameterList&   params,
                         DRT::Discretization&      discretization,
                         std::vector<int>&         lm,
                         Epetra_SerialDenseMatrix& elemat1_epetra,
                         Epetra_SerialDenseMatrix& elemat2_epetra,
                         Epetra_SerialDenseVector& elevec1_epetra,
                         Epetra_SerialDenseVector& elevec2_epetra,
                         Epetra_SerialDenseVector& elevec3_epetra,
                         Teuchos::RCP<MAT::Material> mat);

    /// EvaluateService
    /*!
      Evaluate other stuff that might be needed
     */
    virtual int EvaluateService(Artery*                    ele,
                                const ARTERY::Action       action,
                                Teuchos::ParameterList&    params,
                                DRT::Discretization&       discretization,
                                std::vector<int>&          lm,
                                Epetra_SerialDenseMatrix&  elemat1_epetra,
                                Epetra_SerialDenseMatrix&  elemat2_epetra,
                                Epetra_SerialDenseVector&  elevec1_epetra,
                                Epetra_SerialDenseVector&  elevec2_epetra,
                                Epetra_SerialDenseVector&  elevec3_epetra,
                                Teuchos::RCP<MAT::Material> mat);


  protected:

    /*!
      \brief calculate element matrix and rhs

      \param ele              (i) the element those matrix is calculated
      \param sysmat           (o) element matrix to calculate
      \param rhs              (o) element rhs to calculate
      \param material         (i) artery material/dimesion
      */
    void Sysmat(Artery*                                  ele,
                DRT::Discretization&                     discretization,
                std::vector<int>&                        lm,
                LINALG::Matrix<my::iel_,my::iel_>&       sysmat,
                LINALG::Matrix<my::iel_,       1>&       rhs,
                Teuchos::RCP<const MAT::Material>        material
                );

  };

}

}



#endif /* SRC_DRT_ART_NET_ARTERY_ELE_CALC_PRES_BASED_H_ */
