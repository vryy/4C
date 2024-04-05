/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of HDG transport element

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_ELE_CALC_HDG_HPP
#define FOUR_C_SCATRA_ELE_CALC_HDG_HPP


#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_shapevalues_hdg.hpp"
#include "baci_inpar_scatra.hpp"
#include "baci_scatra_ele_calc.hpp"
#include "baci_scatra_ele_hdg.hpp"
#include "baci_scatra_ele_interface.hpp"

BACI_NAMESPACE_OPEN



namespace DRT
{
  namespace ELEMENTS
  {
    //! Scatra HDG element implementation
    template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype>>
    class ScaTraEleCalcHDG : public ScaTraEleInterface
    {
     public:
      //! nen_: number of element nodes (T. Hughes: The Finite Element Method)
      static constexpr unsigned int nen_ = CORE::FE::num_nodes<distype>;

      //! number of space dimensions
      static constexpr unsigned int nsd_ = probdim;

      //! number of faces on element
      static constexpr unsigned int nfaces_ = CORE::FE::num_faces<distype>;

      /// Evaluate supporting methods of the element
      /*!
        Interface function for supporting methods of the element
       */

      //! Singleton access method
      static ScaTraEleCalcHDG<distype, probdim>* Instance(const int numdofpernode,
          const int numscal, const std::string& disname, bool create = true);

      //! evaluate service routine
      int EvaluateService(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) override;

      //! interpolates an HDG solution to the element nodes for output
      virtual int NodeBasedValues(DRT::Element* ele, DRT::Discretization& discretization,
          CORE::LINALG::SerialDenseVector& elevec1);

      //! initialize the shape functions and solver to the given element (degree is runtime
      //! parameter)
      void InitializeShapes(const DRT::Element* ele, const std::string& disname);

      //! Evaluate the element (Generic virtual interface function. Called via base pointer.)
      int Evaluate(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
          CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override;

      //! evaluate action for off-diagonal system matrix block
      int EvaluateOD(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
          CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override
      {
        dserror("Not implemented!");
        return -1;
      }

      //! Setup element evaluation
      int SetupCalc(DRT::Element* ele, DRT::Discretization& discretization) override { return 0; }

      //! projection of Dirichlet function field
      int ProjectDirichField(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          CORE::LINALG::SerialDenseVector& elevec1);

      //! update interior variables
      int UpdateInteriorVariables(DRT::ELEMENTS::ScaTraHDG* ele, Teuchos::ParameterList& params,
          CORE::LINALG::SerialDenseVector& elevec);

      //! set initial field
      int SetInitialField(const DRT::Element* ele, Teuchos::ParameterList& params,
          CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2);

      //! project field
      int ProjectField(const DRT::Element* ele, DRT::Discretization& discretization,
          Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2, DRT::Element::LocationArray& la);

      //! project material field
      virtual int ProjectMaterialField(const DRT::Element* ele) { return 0; };

      //! calc p-adaptivity
      int CalcPAdaptivity(const DRT::Element* ele, DRT::Discretization& discretization,
          Teuchos::ParameterList& params);

      //! calc error
      int CalcError(const DRT::Element* ele, Teuchos::ParameterList& params,
          CORE::LINALG::SerialDenseVector& elevec);

     protected:
      /// (private) protected constructor, since we are a Singleton.
      /// this constructor is called from a derived class
      /// -> therefore, it has to be protected instead of private
      ScaTraEleCalcHDG(const int numdofpernode, const int numscal, const std::string& disname);

      //! get the material parameters
      virtual void GetMaterialParams(DRT::Element* ele  //!< the element we are dealing with
      );

      //! get the material parameters before first timestep
      virtual void PrepareMaterialParams(DRT::Element* ele  //!< the element we are dealing with
      );

      //! evaluate material
      virtual void Materials(
          const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
          const int k,                                       //!< id of current scalar
          CORE::LINALG::SerialDenseMatrix& difftensor,       //!< diffusion tensor
          CORE::LINALG::SerialDenseVector& ivecn,            //!< reaction term at time n
          CORE::LINALG::SerialDenseVector& ivecnp,           //!< reaction term at time n+1
          CORE::LINALG::SerialDenseMatrix& ivecnpderiv       //!< reaction term derivative
      )
      {
        return;
      };

      //! evaluate material before first timestep
      virtual void PrepareMaterials(DRT::Element* ele,       //!< the element we are dealing with
          const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
          const int k,                                       //!< id of current scalar
          Teuchos::RCP<std::vector<CORE::LINALG::SerialDenseMatrix>>
              difftensor  //!< diffusion tensor
      );

      //! stores the material internal state in a vector for output and restart
      virtual void GetMaterialInternalState(const DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization)
      {
        return;
      };

      //! stores the restart information in the material internal state
      virtual void SetMaterialInternalState(const DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization)
      {
        return;
      };

      //! local data object
      Teuchos::RCP<CORE::FE::ShapeValues<distype>> shapes_;
      Teuchos::RCP<CORE::FE::ShapeValuesFace<distype>> shapesface_;

      //! extracted values from concentrations and gradients
      CORE::LINALG::SerialDenseVector interiorPhin_;
      //! extracted values from concentrations
      CORE::LINALG::SerialDenseVector interiorPhinp_;

      //! get time step
      double Dt() { return localSolver_->scatraparatimint_->Dt(); }



      //! update time dependent material
      virtual void TimeUpdateMaterial(const DRT::Element* ele  //!< the element we are dealing with
      )
      {
        return;
      };

      //! element initialization at the first time step
      void ElementInit(DRT::Element* ele);

      /*========================================================================*/
      //! @name dofs and nodes
      /*========================================================================*/

      //! number of dof per node
      const int numdofpernode_;

      //! number of transported scalars (numscal_ <= numdofpernode_)
      const int numscal_;

      //! use complete polynomial space
      bool usescompletepoly_;

      //! pointer to general scalar transport parameter class
      DRT::ELEMENTS::ScaTraEleParameterStd* scatrapara_;

     private:
      //! local solver that inverts local problem on an element and can solve with various vectors
      struct LocalSolver
      {
        static constexpr unsigned int nsd_ = ScaTraEleCalcHDG<distype, probdim>::nsd_;
        static constexpr unsigned int nfaces_ = ScaTraEleCalcHDG<distype, probdim>::nfaces_;
        int onfdofs_;

        LocalSolver(const DRT::Element* ele, CORE::FE::ShapeValues<distype>& shapeValues,
            CORE::FE::ShapeValuesFace<distype>& shapeValuesFace, bool completepoly,
            const std::string& disname, int numscal);

        //! compute the residual
        void ComputeResidual(Teuchos::ParameterList& params,
            CORE::LINALG::SerialDenseVector& elevec, CORE::LINALG::SerialDenseMatrix& elemat1,
            CORE::LINALG::SerialDenseVector& interiorPhin, CORE::LINALG::SerialDenseVector& tracen,
            CORE::LINALG::SerialDenseVector& tracenp, const DRT::ELEMENTS::ScaTraHDG* hdgele);

        //! compute Neumann boundary conditions
        void ComputeNeumannBC(DRT::Element* ele, Teuchos::ParameterList& params, int face,
            CORE::LINALG::SerialDenseVector& elevec, int indexstart);

        //! compute interior matrices
        void ComputeInteriorMatrices(DRT::ELEMENTS::ScaTraHDG* hdgele);

        //! compute interior matrices for Tet elements
        void ComputeInteriorMatricesTet(DRT::ELEMENTS::ScaTraHDG* hdgele);

        //! compute interior matrices
        void ComputeInteriorMatricesAll(DRT::ELEMENTS::ScaTraHDG* hdgele);

        //! calls local solver to compute matrices: internal and face
        void ComputeMatrices(DRT::Element* ele);

        //! compute face matrices
        void ComputeFaceMatrices(const int face, int indexstart, DRT::ELEMENTS::ScaTraHDG* hdgele);

        //! condense the local matrix (involving interior concentration gradients and
        //! concentrations) into the element matrix for the trace and similarly for the residuals
        void CondenseLocalPart(DRT::ELEMENTS::ScaTraHDG* hdgele);

        //! Compute divergence of current source (ELEMAG)
        void ComputeSource(
            const DRT::Element* ele, CORE::LINALG::SerialDenseVector& elevec1, const double time);

        //! add diffusive term to element matrix
        void AddDiffMat(
            CORE::LINALG::SerialDenseMatrix& elemat, const DRT::ELEMENTS::ScaTraHDG* hdgele);

        //! add reaction term to element matrix
        void AddReacMat(
            CORE::LINALG::SerialDenseMatrix& elemat, const DRT::ELEMENTS::ScaTraHDG* hdgele);

        //! set material parameter
        void SetMaterialParameter(DRT::ELEMENTS::ScaTraHDG* hdgele,
            CORE::LINALG::SerialDenseVector& ivecn, CORE::LINALG::SerialDenseVector& ivecnp,
            CORE::LINALG::SerialDenseMatrix& ivecnpderiv);

        //! prepare material parameter in first timestep
        void PrepareMaterialParameter(
            DRT::ELEMENTS::ScaTraHDG* hdgele, CORE::LINALG::SerialDenseMatrix& difftensor);


        // convention: we sort the entries in the matrices the following way:
        // first come the concentration, then the concentration  gradients, and finally the trace

        //! evaluated shape values
        Teuchos::RCP<CORE::FE::ShapeValues<distype>> shapes_;

        //! evaluated shape values on face
        Teuchos::RCP<CORE::FE::ShapeValuesFace<distype>> shapesface_;  /// evaluated shape values

        // Element matrices if one wants to compute them on the fly instead of storing them on the
        // element
        //      CORE::LINALG::SerialDenseMatrix  Amat;     /// concentrations - concentrations
        //      CORE::LINALG::SerialDenseMatrix  Bmat;     /// concentrations - concentrations
        //      gradients CORE::LINALG::SerialDenseMatrix  Cmat;     /// concentration - trace
        //      CORE::LINALG::SerialDenseMatrix  Dmat;     /// concentrations gradients -
        //      concentrations gradients CORE::LINALG::SerialDenseMatrix  Emat;     /// trace -
        //      concentrations gradients CORE::LINALG::SerialDenseMatrix  Gmat;     ///
        //      concentrations gradients CORE::LINALG::SerialDenseMatrix  Hmat;     /// trace -trace
        //      CORE::LINALG::SerialDenseMatrix  Mmat;     /// mass matrix (concentrations -
        //      concentrations) CORE::LINALG::SerialDenseMatrix  EmatT;    /// trace -
        //      concentrations gradients (E^T) CORE::LINALG::SerialDenseMatrix  BmatMT;   ///
        //      concentrations gradients- concentrations (-B^T) CORE::LINALG::SerialDenseMatrix
        //      Kmat;   /// condensed matrix

        // @name variables for the reaction term
        //!@{
        //! reaction term at time n
        //      CORE::LINALG::SerialDenseVector  Ivecn_;

        //! reaction term at time n+1
        // CORE::LINALG::SerialDenseVector  Ivecnp_;

        //! derivative of reaction term at time n+1
        //      CORE::LINALG::SerialDenseMatrix  Imatnpderiv_;
        //!@}

        ////      CORE::LINALG::SerialDenseMatrix  invAmat;     /// inverse of Amat
        //      CORE::LINALG::SerialDenseMatrix  invAMmat;     /// inverse of [A + (1/(dt*theta))*M]
        //
        //      // auxiliary stuff
        //      CORE::LINALG::SerialDenseMatrix  massPart;
        //      CORE::LINALG::SerialDenseMatrix  massPartW;
        //      CORE::LINALG::SerialDenseMatrix  BTAMmat;
        //      CORE::LINALG::SerialDenseMatrix  invCondmat;
        //      CORE::LINALG::SerialDenseMatrix  Xmat;

        //! pointer to general scalar transport parameter class
        DRT::ELEMENTS::ScaTraEleParameterStd* scatrapara_;

        //      Teuchos::RCP<DRT::ELEMENTS::ScaTraEleParameterBase> scatrapara_; //! pointer to
        //      parameter list
        //! pointer to time parameter list
        Teuchos::RCP<DRT::ELEMENTS::ScaTraEleParameterTimInt> scatraparatimint_;

        /*========================================================================*/
        //! @name diffusions and reaction coefficient
        /*========================================================================*/

        // diffusion tensor stored on the element, if necessary this can be changed
        //      //! diffusion tensor
        //      CORE::LINALG::SerialDenseMatrix diff_;
        //      //! inverse diffusion tensor
        //      CORE::LINALG::SerialDenseMatrix invdiff_;

        //! scalar raeaction coefficient
        //      std::vector<double> reacoeff_;
      };

      //! reads from global vectors
      void ReadGlobalVectors(
          DRT::Element* ele, DRT::Discretization& discretization, DRT::Element::LocationArray& la);

      //! local solver object
      Teuchos::RCP<LocalSolver> localSolver_;

      /*========================================================================*/
      //! @name trace and interior concentrations and gradients
      /*========================================================================*/

      //! extracted values from concentrations
      CORE::LINALG::SerialDenseVector tracen_;

      //! extracted local values (concentration gradients) at n+alpha_f
      CORE::LINALG::SerialDenseVector interiorGradPhin_;

      //! extracted values from trace solution vector at n-m
      CORE::LINALG::SerialDenseVector tracenm_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif
