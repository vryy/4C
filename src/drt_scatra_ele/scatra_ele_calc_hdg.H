/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of HDG transport element

\level 3

\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_HDG_H
#define SCATRA_ELE_CALC_HDG_H


#include "scatra_ele_interface.H"
#include "scatra_ele_hdg.H"
#include "scatra_ele_calc.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_fem_general/drt_utils_shapevalues_hdg.H"



namespace DRT
{
  namespace ELEMENTS
  {
    //! Scatra HDG element implementation
    template <DRT::Element::DiscretizationType distype,
        int probdim = DRT::UTILS::DisTypeToDim<distype>::dim>
    class ScaTraEleCalcHDG : public ScaTraEleInterface
    {
     public:
      //! nen_: number of element nodes (T. Hughes: The Finite Element Method)
      static const unsigned int nen_ =
          DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      //! number of space dimensions
      static const unsigned int nsd_ = probdim;

      //! number of faces on element
      static const unsigned int nfaces_ = DRT::UTILS::DisTypeToNumFaces<distype>::numFaces;

      /// Evaluate supporting methods of the element
      /*!
        Interface function for supporting methods of the element
       */

      virtual ~ScaTraEleCalcHDG() { return; }

      //! Singleton access method
      static ScaTraEleCalcHDG<distype, probdim>* Instance(const int numdofpernode,
          const int numscal, const std::string& disname, bool create = true);

      //! called upon destruction
      virtual void Done();

      //! evaluate service routine
      virtual int EvaluateService(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra);

      //! interpolates an HDG solution to the element nodes for output
      virtual int NodeBasedValues(DRT::Element* ele, DRT::Discretization& discretization,
          Epetra_SerialDenseVector& elevec1);

      //! initialize the shape functions and solver to the given element (degree is runtime
      //! parameter)
      void InitializeShapes(const DRT::Element* ele, const std::string& disname);

      //! Evaluate the element (Generic virtual interface function. Called via base pointer.)
      virtual int Evaluate(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
          Epetra_SerialDenseVector& elevec3);

      //! evaluate action for off-diagonal system matrix block
      virtual int EvaluateOD(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
          Epetra_SerialDenseVector& elevec3)
      {
        dserror("Not implemented!");
        return -1;
      }

      //! Setup element evaluation
      virtual int SetupCalc(DRT::Element* ele, DRT::Discretization& discretization) { return 0; }

      //! projection of Dirichlet function field
      int ProjectDirichField(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseVector& elevec1);

      //! update interior variables
      int UpdateInteriorVariables(DRT::ELEMENTS::ScaTraHDG* ele, Teuchos::ParameterList& params,
          Epetra_SerialDenseVector& elevec);

      //! set initial field
      int SetInitialField(const DRT::Element* ele, Teuchos::ParameterList& params,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2);

      //! project field
      int ProjectField(const DRT::Element* ele, DRT::Discretization& discretization,
          Teuchos::ParameterList& params, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, DRT::Element::LocationArray& la);

      //! project material field
      virtual int ProjectMaterialField(const DRT::Element* ele) { return 0; };

      //! calc p-adaptivity
      int CalcPAdaptivity(const DRT::Element* ele, DRT::Discretization& discretization,
          Teuchos::ParameterList& params);

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
          Epetra_SerialDenseMatrix& difftensor,              //!< diffusion tensor
          Epetra_SerialDenseVector& ivecn,                   //!< reaction term at time n
          Epetra_SerialDenseVector& ivecnp,                  //!< reaction term at time n+1
          Epetra_SerialDenseMatrix& ivecnpderiv              //!< reaction term derivative
      )
      {
        return;
      };

      //! evaluate material before first timestep
      virtual void PrepareMaterials(DRT::Element* ele,       //!< the element we are dealing with
          const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
          const int k,                                       //!< id of current scalar
          Teuchos::RCP<std::vector<Epetra_SerialDenseMatrix>> difftensor  //!< diffusion tensor
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
      Teuchos::RCP<DRT::UTILS::ShapeValues<distype>> shapes_;
      Teuchos::RCP<DRT::UTILS::ShapeValuesFace<distype>> shapesface_;

      //! extracted values from concentrations and gradients
      Epetra_SerialDenseVector interiorPhin_;
      //! extracted values from concentrations
      Epetra_SerialDenseVector interiorPhinp_;

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
        static const unsigned int nsd_ = ScaTraEleCalcHDG<distype, probdim>::nsd_;
        static const unsigned int nfaces_ = ScaTraEleCalcHDG<distype, probdim>::nfaces_;
        int onfdofs_;

        LocalSolver(const DRT::Element* ele, DRT::UTILS::ShapeValues<distype>& shapeValues,
            DRT::UTILS::ShapeValuesFace<distype>& shapeValuesFace, bool completepoly,
            const std::string& disname, int numscal);

        //! compute the residual
        void ComputeResidual(Teuchos::ParameterList& params, Epetra_SerialDenseVector& elevec,
            Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseVector& interiorPhin,
            Epetra_SerialDenseVector& tracen, Epetra_SerialDenseVector& tracenp,
            const DRT::ELEMENTS::ScaTraHDG* hdgele);

        //! compute Neumann boundary conditions
        void ComputeNeumannBC(DRT::Element* ele, Teuchos::ParameterList& params, int face,
            Epetra_SerialDenseVector& elevec, int indexstart);

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

        //! add diffusive term to element matrix
        void AddDiffMat(Epetra_SerialDenseMatrix& elemat, const DRT::ELEMENTS::ScaTraHDG* hdgele);

        //! add reaction term to element matrix
        void AddReacMat(Epetra_SerialDenseMatrix& elemat, const DRT::ELEMENTS::ScaTraHDG* hdgele);

        //! set material parameter
        void SetMaterialParameter(DRT::ELEMENTS::ScaTraHDG* hdgele, Epetra_SerialDenseVector& ivecn,
            Epetra_SerialDenseVector& ivecnp, Epetra_SerialDenseMatrix& ivecnpderiv);

        //! prepare material parameter in first timestep
        void PrepareMaterialParameter(
            DRT::ELEMENTS::ScaTraHDG* hdgele, Epetra_SerialDenseMatrix& difftensor);


        // convention: we sort the entries in the matrices the following way:
        // first come the concentration, then the concentration  gradients, and finally the trace

        //! evaluated shape values
        Teuchos::RCP<DRT::UTILS::ShapeValues<distype>> shapes_;

        //! evaluated shape values on face
        Teuchos::RCP<DRT::UTILS::ShapeValuesFace<distype>> shapesface_;  /// evaluated shape values

        // Element matrices if one wants to compute them on the fly instead of storing them on the
        // element
        //      Epetra_SerialDenseMatrix  Amat;     /// concentrations - concentrations
        //      Epetra_SerialDenseMatrix  Bmat;     /// concentrations - concentrations gradients
        //      Epetra_SerialDenseMatrix  Cmat;     /// concentration - trace
        //      Epetra_SerialDenseMatrix  Dmat;     /// concentrations gradients - concentrations
        //      gradients Epetra_SerialDenseMatrix  Emat;     /// trace - concentrations gradients
        //      Epetra_SerialDenseMatrix  Gmat;     /// concentrations gradients
        //      Epetra_SerialDenseMatrix  Hmat;     /// trace -trace
        //      Epetra_SerialDenseMatrix  Mmat;     /// mass matrix (concentrations -
        //      concentrations) Epetra_SerialDenseMatrix  EmatT;    /// trace - concentrations
        //      gradients (E^T) Epetra_SerialDenseMatrix  BmatMT;   /// concentrations gradients-
        //      concentrations (-B^T) Epetra_SerialDenseMatrix  Kmat;   /// condensed matrix

        // @name variables for the reaction term
        //!@{
        //! reaction term at time n
        //      Epetra_SerialDenseVector  Ivecn_;

        //! reaction term at time n+1
        // Epetra_SerialDenseVector  Ivecnp_;

        //! derivative of reaction term at time n+1
        //      Epetra_SerialDenseMatrix  Imatnpderiv_;
        //!@}

        ////      Epetra_SerialDenseMatrix  invAmat;     /// inverse of Amat
        //      Epetra_SerialDenseMatrix  invAMmat;     /// inverse of [A + (1/(dt*theta))*M]
        //
        //      // auxiliary stuff
        //      Epetra_SerialDenseMatrix  massPart;
        //      Epetra_SerialDenseMatrix  massPartW;
        //      Epetra_SerialDenseMatrix  BTAMmat;
        //      Epetra_SerialDenseMatrix  invCondmat;
        //      Epetra_SerialDenseMatrix  Xmat;

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
        //      Epetra_SerialDenseMatrix diff_;
        //      //! inverse diffusion tensor
        //      Epetra_SerialDenseMatrix invdiff_;

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
      Epetra_SerialDenseVector tracen_;

      //! extracted local values (concentration gradients) at n+alpha_f
      Epetra_SerialDenseVector interiorGradPhin_;

      //! extracted values from trace solution vector at n-m
      Epetra_SerialDenseVector tracenm_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif /* SCATRA_ELE_CALC_HDG_H_ */
