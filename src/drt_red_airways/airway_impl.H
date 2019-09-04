/*---------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of airway_impl element

\maintainer Carolin Geitner

\level 3

*/
/*---------------------------------------------------------------------*/



#ifndef AIRWAY_IMPL_H
#define AIRWAY_IMPL_H

#include "red_airway.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../linalg/linalg_fixedsizematrix.H"


namespace DRT
{
  namespace ELEMENTS
  {
    /// Interface base class for airway_impl
    /*!
      This class exists to provide a common interface for all template
      versions of airway_impl. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of airway_impl.
     */
    class RedAirwayImplInterface
    {
     public:
      /// Empty constructor
      RedAirwayImplInterface() {}
      /// Empty destructor
      virtual ~RedAirwayImplInterface() {}
      /// Evaluate the element
      virtual int Evaluate(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat) = 0;

      virtual void Initial(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseVector& radii_in, Epetra_SerialDenseVector& radii_out,
          Teuchos::RCP<const MAT::Material> material) = 0;

      virtual void EvaluateTerminalBC(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1_epetra, Teuchos::RCP<MAT::Material> mat) = 0;

      virtual void CalcFlowRates(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat) = 0;

      virtual void CalcElemVolume(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat) = 0;

      virtual bool GetCurveValAtCond(double& bcVal, DRT::Node* node, std::string conditionName,
          std::string optionName, std::string conditionType, double time) = 0;

      virtual void GetCoupledValues(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void GetJunctionVolumeMix(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& volumeMix_np,
          std::vector<int>& lm, Teuchos::RCP<MAT::Material> material) = 0;

      virtual void SolveScatra(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& scatra_np,
          Epetra_SerialDenseVector& volumeMix_np, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;


      virtual void SolveScatraBifurcations(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& scatra_np,
          Epetra_SerialDenseVector& volumeMix_np, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void CalcCFL(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void UpdateScatra(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void UpdateElem12Scatra(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void EvalPO2FromScatra(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void EvalNodalEssentialValues(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& nodal_surface,
          Epetra_SerialDenseVector& nodal_volume, Epetra_SerialDenseVector& nodal_flow,
          std::vector<int>& lm, Teuchos::RCP<MAT::Material> material) = 0;

      virtual double GetElementLength(RedAirway* ele) = 0;


      /// Internal implementation class for airway element
      static RedAirwayImplInterface* Impl(DRT::ELEMENTS::RedAirway* airway);
    };


    /// Internal airway implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the airway element. Additionally the method Sysmat()
      provides a clean and fast element implementation.

      <h3>Purpose</h3>

      \author ismail
      \date 01/09
    */

    template <DRT::Element::DiscretizationType distype>
    class AirwayImpl : public RedAirwayImplInterface
    {
     public:
      /// Constructor
      explicit AirwayImpl();

      //! number of nodes
      static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;


      /// Evaluate
      /*!
        The evaluate function for the general airway case.
       */
      virtual int Evaluate(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat);

      /*!
        \brief Calculate virtual trajectory xnp and state of airway (open/closed)

        \param ele              (i) the element those matrix is calculated
        \param epn              (i) nodal pressure at n
        \param dt               (i) timestep
        */
      void EvaluateCollapse(
          RedAirway* ele, Epetra_SerialDenseVector& epn, Teuchos::ParameterList& params, double dt);
      /*!
        \brief Calculate Pextn and Pextnp

        \param ele              (i) the element those matrix is calculated
        \param epn              (i) nodal pressure at n
        \param dt               (i) timestep
        */
      void ComputePext(RedAirway* ele, Teuchos::RCP<const Epetra_Vector> pn,
          Teuchos::RCP<const Epetra_Vector> pnp, Teuchos::ParameterList& params);
      /*!
        \brief calculate element matrix and rhs

        \param ele              (i) the element those matrix is calculated
        \param eqnp             (i) nodal volumetric flow rate at n+1
        \param evelnp           (i) nodal velocity at n+1
        \param eareanp          (i) nodal cross-sectional area at n+1
        \param eprenp           (i) nodal pressure at n+1
        \param estif            (o) element matrix to calculate
        \param eforce           (o) element rhs to calculate
        \param material         (i) airway material/dimesion
        \param time             (i) current simulation time
        \param dt               (i) timestep
        \param compute_awacinter(i) computing airway-acinus interdependency
        */
      void Sysmat(RedAirway* ele, Epetra_SerialDenseVector& epnp, Epetra_SerialDenseVector& epn,
          Epetra_SerialDenseVector& epnm, Epetra_SerialDenseMatrix& estif,
          Epetra_SerialDenseVector& eforce, Teuchos::RCP<const MAT::Material> material,
          Teuchos::ParameterList& params, double time, double dt, bool compute_awacinter);


      virtual void EvaluateTerminalBC(RedAirway* ele, Teuchos::ParameterList& params,
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
        \param material         (i) airway material/dimesion
        \param time             (i) current simulation time
        \param dt               (i) timestep
        */
      virtual void Initial(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseVector& radii_in, Epetra_SerialDenseVector& radii_out,
          Teuchos::RCP<const MAT::Material> material);

      /*!
       \Essential functions to compute the results of essentail matrices
      */
      virtual void CalcFlowRates(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat);

      /*!
       \Essential functions to compute the volume of elements
      */
      virtual void CalcElemVolume(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat);

      /*!
       \Essential functions to compute curve value at a Nodal Prescribed condition
      */
      virtual bool GetCurveValAtCond(double& bcVal, DRT::Node* node, std::string conditionName,
          std::string optionName, std::string conditionType, double time);

      /*!
       \Essential functions to compute the volume mixing and  flowing into a junction
      */
      virtual void GetJunctionVolumeMix(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& junctionVolumeMix_np,
          std::vector<int>& lm, Teuchos::RCP<MAT::Material> material);

      /*!
       \Essential functions to compute the volume mixing and  flowing into a junction
      */
      virtual void SolveScatra(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& scatra_np,
          Epetra_SerialDenseVector& volumeMix_np, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material);

      /*!
       \Essential functions to compute the volume mixing and  flowing into a junction
      */
      virtual void SolveScatraBifurcations(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& scatra_np,
          Epetra_SerialDenseVector& volumeMix_np, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material);
      /*!
       \Essential functions to evaluate the coupled results
      */
      virtual void GetCoupledValues(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material);

      virtual void CalcCFL(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material);

      virtual void UpdateScatra(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material);

      virtual void UpdateElem12Scatra(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material);

      virtual void EvalPO2FromScatra(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material);

      virtual void EvalNodalEssentialValues(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& nodal_surface,
          Epetra_SerialDenseVector& nodal_volume, Epetra_SerialDenseVector& nodal_flow,
          std::vector<int>& lm, Teuchos::RCP<MAT::Material> material);

      virtual double GetElementLength(RedAirway* ele);

     private:
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
