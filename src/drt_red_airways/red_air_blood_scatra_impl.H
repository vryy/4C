/*---------------------------------------------------------------------*/
/*!

\brief Incomplete! - Purpose: Internal implementation of RedAirBloodScatra element

\maintainer Carolin Geitner

\level 3

*/
/*---------------------------------------------------------------------*/



#ifndef RED_AIR_BLOOD_SCATRA_IMPL_H_
#define RED_AIR_BLOOD_SCATRA_IMPL_H_

#include "red_airway.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../linalg/linalg_fixedsizematrix.H"


namespace DRT
{
  namespace ELEMENTS
  {
    /// Interface base class for red_air_blood_scatra_impl
    /*!
      This class exists to provide a common interface for all template
      versions of red_air_blood_scatra_impl. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of ired_air_blood_scatra_impl.
     */
    class RedAirBloodScatraImplInterface
    {
     public:
      /// Empty constructor
      RedAirBloodScatraImplInterface() {}
      /// Empty destructor
      virtual ~RedAirBloodScatraImplInterface() {}
      /// Evaluate the element
      virtual int Evaluate(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat) = 0;

      virtual void Initial(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<const MAT::Material> material) = 0;

      virtual void EvaluateTerminalBC(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1_epetra, Teuchos::RCP<MAT::Material> mat) = 0;

      virtual void CalcFlowRates(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& a_volumen,
          Epetra_SerialDenseVector& a_volumenp, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat) = 0;

      virtual void GetCoupledValues(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void SolveBloodAirTransport(RedAirBloodScatra* ele, Epetra_SerialDenseVector& dscatra,
          Epetra_SerialDenseVector& dvO2, Epetra_SerialDenseVector& scatraAcinus,
          Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      /// Internal implementation class for acinus element
      static RedAirBloodScatraImplInterface* Impl(DRT::ELEMENTS::RedAirBloodScatra* elem);
    };


    /// Internal acinus implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the acinus element. Additionally the method Sysmat()
      provides a clean and fast element implementation.

      <h3>Purpose</h3>

      \author ismail
      \date 01/09
    */

    template <DRT::Element::DiscretizationType distype>
    class RedAirBloodScatraImpl : public RedAirBloodScatraImplInterface
    {
     public:
      /// Constructor
      explicit RedAirBloodScatraImpl();

      //! number of nodes
      static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;


      /// Evaluate
      /*!
        The evaluate function for the general acinus case.
       */
      virtual int Evaluate(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat);

      /*!
        \brief calculate element matrix and rhs

        \param ele              (i) the element those matrix is calculated
        \param eqnp             (i) nodal volumetric flow rate at n+1
        \param evelnp           (i) nodal velocity at n+1
        \param eareanp          (i) nodal cross-sectional area at n+1
        \param eprenp           (i) nodal pressure at n+1
        \param estif            (o) element matrix to calculate
        \param eforce           (o) element rhs to calculate
        \param material         (i) acinus material/dimesion
        \param time             (i) current simulation time
        \param dt               (i) timestep
        */
      void Sysmat(RedAirBloodScatra* ele, Epetra_SerialDenseVector& epnp,
          Epetra_SerialDenseVector& epn, Epetra_SerialDenseVector& epnm,
          Epetra_SerialDenseMatrix& estif, Epetra_SerialDenseVector& eforce,
          Teuchos::RCP<const MAT::Material> material, Teuchos::ParameterList& params, double time,
          double dt);


      virtual void EvaluateTerminalBC(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& disctretization, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1_epetra, Teuchos::RCP<MAT::Material> mat);

      /*!
        \brief get the initial values of the degrees of freedome at the node

        \param ele              (i) the element those matrix is calculated
        \param eqnp             (i) nodal volumetric flow rate at n+1
        \param evelnp           (i) nodal velocity at n+1
        \param eareanp          (i) nodal cross-sectional area at n+1
        \param eprenp           (i) nodal pressure at n+1
        \param estif            (o) element matrix to calculate
        \param eforce           (o) element rhs to calculate
        \param material         (i) acinus material/dimesion
        \param time             (i) current simulation time
        \param dt               (i) timestep
        */
      virtual void Initial(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<const MAT::Material> material);

      /*!
       \Essential functions to compute the results of essentail matrices
      */
      virtual void CalcFlowRates(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& a_volumen_strain_np,
          Epetra_SerialDenseVector& a_volumenp, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat);

      /*!
       \Essential functions to evaluate the coupled results
      */
      virtual void GetCoupledValues(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material);

      virtual void SolveBloodAirTransport(RedAirBloodScatra* ele, Epetra_SerialDenseVector& dscatra,
          Epetra_SerialDenseVector& dvO2, Epetra_SerialDenseVector& scatraAcinus,
          Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material);

     private:
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // RED_AIR_BLOOD_SCATRA_IMPL_H_
