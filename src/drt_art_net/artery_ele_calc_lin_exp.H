
/*----------------------------------------------------------------------*/
/*!
\file artery_ele_calc_lin_exp.H

\brief Internal implementation of LinExp artery element

\maintainer Lena Yoshihara

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef ARTERY_LIN_EXP_H
#define ARTERY_LIN_EXP_H

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

      \author ismail
      \date 01/09
    */

    template <DRT::Element::DiscretizationType distype>
    class ArteryEleCalcLinExp : public ArteryEleCalc<distype>
    {
     private:
      typedef ArteryEleCalc<distype> my;

      /// private constructor, since we are a Singleton.
      ArteryEleCalcLinExp(const int numdofpernode, const std::string& disname);

     public:
      /// virtual destructor
      virtual ~ArteryEleCalcLinExp(){};

      //! Singleton access method
      static ArteryEleCalcLinExp<distype>* Instance(const int numdofpernode,
          const std::string& disname, const ArteryEleCalcLinExp* delete_me = NULL);

      /// called upon destruction
      virtual void Done();

      /// Evaluate
      /*!
        The evaluate function for the general artery case.
       */
      virtual int Evaluate(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat);

      virtual int ScatraEvaluate(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat);

      virtual int EvaluateService(Artery* ele, const ARTERY::Action action,
          Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
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
        \param material         (i) artery material/dimesion
        \param time             (i) current simulation time
        \param dt               (i) timestep
        */
      void Sysmat(Artery* ele, const LINALG::Matrix<my::iel_, 1>& eqnp,
          const LINALG::Matrix<my::iel_, 1>& eareanp,
          LINALG::Matrix<2 * my::iel_, 2 * my::iel_>& estif,
          LINALG::Matrix<2 * my::iel_, 1>& eforce, Teuchos::RCP<const MAT::Material> material,
          double dt);

      void ScatraSysmat(Artery* ele, const LINALG::Matrix<2 * my::iel_, 1>& escatran,
          const LINALG::Matrix<my::iel_, 1>& eqnp, const LINALG::Matrix<my::iel_, 1>& eqn,
          const LINALG::Matrix<my::iel_, 1>& eareanp, const LINALG::Matrix<my::iel_, 1>& earean,
          LINALG::Matrix<2 * my::iel_, 2 * my::iel_>& estif,
          LINALG::Matrix<2 * my::iel_, 1>& eforce, Teuchos::RCP<const MAT::Material> material,
          double dt);

      virtual bool SolveRiemann(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<const MAT::Material> mat);

      virtual void EvaluateTerminalBC(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& disctretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat);

      virtual void EvaluateScatraBC(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& disctretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat);

      virtual void CalcPostprocessingValues(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat);

      virtual void CalcScatraFromScatraFW(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material);

      virtual void EvaluateWfAndWb(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material);

      virtual void SolveScatraAnalytically(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material);

      /*!
        \brief get the initial values of the degrees of freedome at the node

        \param ele              (i) the element those matrix is calculated
        \param eqnp             (i) nodal volumetric flow rate at n+1
        \param evelnp           (i) nodal velocity at n+1
        \param eareanp          (i) nodal cross-sectional area at n+1
        \param eprenp           (i) nodal pressure at n+1
        \param estif            (o) element matrix to calculate
        \param eforce           (o) element rhs to calculate
        \param material         (i) artery material/dimesion
        \param time             (i) current simulation time
        \param dt               (i) timestep
        */
      virtual void Initial(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<const MAT::Material> material);

      /*!
       \Essential functions to compute the results of essentail matrices
      */
#if 0
    LINALG::Matrix<2,1> F    (double xi);
    LINALG::Matrix<2,1> B    (double xi);
    LINALG::Matrix<2,2> Bu   (double xi);
    LINALG::Matrix<2,2> H    (double xi);
    LINALG::Matrix<2,1> dFdxi(double xi);
    LINALG::Matrix<2,1> FLW  (double xi, double dt);
    LINALG::Matrix<2,1> BLW  (double xi, double dt);
#endif

     private:
      //! nodal volumetric flow rate at time step "n"
      LINALG::Matrix<my::iel_, 1> qn_;
      //! nodal cross-sectional area at time step "n"
      LINALG::Matrix<my::iel_, 1> an_;
      //! vector containing the initial cross-sectional area at the element nodes
      LINALG::Matrix<my::iel_, 1> area0_;
      //! vector containing the initial thickness at the element nodes
      LINALG::Matrix<my::iel_, 1> th_;
      //! vector containing the initial Youngs modulus at the element nodes
      LINALG::Matrix<my::iel_, 1> young_;
      //! vector containing the fixed external pressure
      LINALG::Matrix<my::iel_, 1> pext_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
