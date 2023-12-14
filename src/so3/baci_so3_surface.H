/*----------------------------------------------------------------------*/
/*! \file

\brief class for evaluation of equations on the structural surface
\level 1


*----------------------------------------------------------------------*/
#ifndef BACI_SO3_SURFACE_H
#define BACI_SO3_SURFACE_H

#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_lib_element.H"
#include "baci_lib_elementtype.H"
#include "baci_lib_node.H"
#include "baci_linalg_serialdensematrix.H"

BACI_NAMESPACE_OPEN


namespace DRT
{
  namespace ELEMENTS
  {
    class StructuralSurfaceType : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "StructuralSurfaceType"; }

      static StructuralSurfaceType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override
      {
      }

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override
      {
        CORE::LINALG::SerialDenseMatrix nullspace;
        dserror("method ComputeNullSpace not implemented!");
        return nullspace;
      }

     private:
      static StructuralSurfaceType instance_;
    };

    /*!
    \brief A general structural surface element (tri3/tri6/quad4/quad8/quad9/nurbs9)

    This element implements all common shapes. It also implements all
    common types of loads currently demanded in baci

    */
    class StructuralSurface : public DRT::FaceElement
    {
     public:
      //! @name Friends
      friend class StructuralSurfaceType;

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner: Processor owning this surface
      \param nnode: Number of nodes attached to this element
      \param nodeids: global ids of nodes attached to this element
      \param nodes: the discretizations map of nodes to build ptrs to nodes from
      \param parent: The parent fluid element of this surface
      \param lsurface: the local surface number of this surface w.r.t. the parent element
      */
      explicit StructuralSurface(int id, int owner, int nnode, const int* nodeids,
          DRT::Node** nodes, DRT::Element* parent, const int lsurface);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      explicit StructuralSurface(const StructuralSurface& old);

      /*!
      \brief Deep copy this instance of an element and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const override;

      /*!
      \brief Get shape type of element
      */
      CORE::FE::CellType Shape() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the parobject.H file.
      */
      inline int UniqueParObjectId() const override
      {
        return StructuralSurfaceType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(CORE::COMM::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Unpack(const std::vector<char>& data) override;

      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      inline int NumDofPerNode(const DRT::Node& node) const override { return numdofpernode_; }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual DRT::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      inline int NumDofPerElement() const override { return 0; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return StructuralSurfaceType::Instance(); }

      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      int NumLine() const override;

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a surface Neumann condition on the fluid3 element

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param condition (in)     : The condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If nullptr on input,

      \return 0 if successful, negative otherwise
      */
      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;

      //! Evaluate method for StructuralSurface-Elements
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elematrix1,
          CORE::LINALG::SerialDenseMatrix& elematrix2, CORE::LINALG::SerialDenseVector& elevector1,
          CORE::LINALG::SerialDenseVector& elevector2,
          CORE::LINALG::SerialDenseVector& elevector3) override;

      //! Evaluate method for StructuralSurface-Elements
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          LocationArray& la, CORE::LINALG::SerialDenseMatrix& elematrix1,
          CORE::LINALG::SerialDenseMatrix& elematrix2, CORE::LINALG::SerialDenseVector& elevector1,
          CORE::LINALG::SerialDenseVector& elevector2,
          CORE::LINALG::SerialDenseVector& elevector3) override;

      //! Evaluate trace inequality and return the maximal eigenvalue
      virtual double EstimateNitscheTraceMaxEigenvalueCombined(std::vector<double>& parent_disp);

      //! Evaluate trace inequality and return the maximal eigenvalue
      virtual double EstimateNitscheTraceMaxEigenvalueTSI(std::vector<double>& parent_disp);

      //! Return a pointer to the parent element of this boundary element
      virtual DRT::Element* ParentElement() const { return ParentMasterElement(); }

      //! Return local surface number
      int LSurfNumber() const { return FaceMasterNumber(); }

     public:  // changed to public rauch 05/2014
      //! action for surface evaluation
      enum ActionType
      {
        none,
        calc_struct_constrvol,
        calc_struct_volconstrstiff,
        calc_struct_constrarea,
        calc_init_vol,
        calc_struct_areaconstrstiff,
        calc_struct_monitarea,
        prestress_update,
        prestress_writerestart,
        prestress_readrestart,
        calc_brownian_motion,
        calc_brownian_motion_damping,
        calc_struct_centerdisp,
        calc_struct_rotation,
        calc_undo_struct_rotation,
        calc_struct_area,
        calc_ref_nodal_normals,
        calc_cur_nodal_normals,
        calc_struct_area_poro,
        calc_fluid_traction,
        calc_cur_normal_at_point,
        mark_immersed_elements,
        calc_struct_robinforcestiff
      };

     protected:
      // don't want = operator
      StructuralSurface& operator=(const StructuralSurface& old);

      // variables

      //! discretization type
      CORE::FE::CellType distype_;

      //! numdofpernode
      int numdofpernode_;

      //! gaussian integration to be used
      CORE::DRT::UTILS::GaussRule2D gaussrule_;

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates normal and detA at gaussian point

      \param detA   (out) : are at gaussian point
      \param normal (out) : The normal at gaussian point, length is detA!
      \param x      (in)  : nodal coords in either material or spatial frame
      \param deriv  (in)  : derivatives of shape functions
      */
      void SurfaceIntegration(double& detA, std::vector<double>& normal,
          const CORE::LINALG::SerialDenseMatrix& x, const CORE::LINALG::SerialDenseMatrix& deriv);

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates the normal at gaussian point

      \param normal (out) : The normal at gaussian point, length is detA!
      \param x      (in)  : nodal coords in either material or spatial frame
      \param deriv  (in)  : derivatives of shape functions
      */
      void SurfaceIntegration(std::vector<double>& normal, const CORE::LINALG::SerialDenseMatrix& x,
          const CORE::LINALG::SerialDenseMatrix& deriv);

      /*!
      \brief Linearize a Neumann boundary condition analytically

      This method computes the analytical surface derivative if necessary (for follower loads).
      */
      void analytical_DSurfaceIntegration(CORE::LINALG::SerialDenseMatrix& d_normal,
          const CORE::LINALG::SerialDenseMatrix& x, const CORE::LINALG::SerialDenseMatrix& deriv);

      /*!
      \brief Linearize a Neumann boundary condition with automatic differentiation

      This method computes the automatic (Sacado) surface derivative if necessary (for follower
      loads). This method is for development purposes only, and currently not used.
      */
      void automatic_DSurfaceIntegration(CORE::LINALG::SerialDenseMatrix& d_normal,
          const CORE::LINALG::SerialDenseMatrix& x, const CORE::LINALG::SerialDenseMatrix& deriv);

      /*!
      \brief Create matrix with material configuration

      \param x  (out)  : nodal coords in material frame
      */
      inline void MaterialConfiguration(CORE::LINALG::SerialDenseMatrix& x) const
      {
        const int numnode = NumNode();
        for (int i = 0; i < numnode; ++i)
        {
          x(i, 0) = Nodes()[i]->X()[0];
          x(i, 1) = Nodes()[i]->X()[1];
          x(i, 2) = Nodes()[i]->X()[2];
        }
        return;
      }

      /*!
      \brief Create matrix with spatial configuration

      \param x     (out)  : nodal coords in spatial frame
      \param disp  (int)  : displacements
      */
      inline void SpatialConfiguration(
          CORE::LINALG::SerialDenseMatrix& x, const std::vector<double>& disp) const
      {
        const int numnode = NumNode();
        for (int i = 0; i < numnode; ++i)
        {
          x(i, 0) = Nodes()[i]->X()[0] + disp[i * 3 + 0];
          x(i, 1) = Nodes()[i]->X()[1] + disp[i * 3 + 1];
          x(i, 2) = Nodes()[i]->X()[2] + disp[i * 3 + 2];
        }
        return;
      }

      /*!
      \brief Create matrix with spatial configuration

      \param x     (out)  : nodal coords in spatial frame
      \param xrefe (out)  : nodal coords in material frame
      \param disp  (int)  : displacements
      */
      inline void SpatialConfiguration(CORE::LINALG::SerialDenseMatrix& x,
          CORE::LINALG::SerialDenseMatrix& xrefe, const std::vector<double>& disp) const
      {
        const int numnode = NumNode();
        for (int i = 0; i < numnode; ++i)
        {
          x(i, 0) = xrefe(i, 0) + disp[i * 3 + 0];
          x(i, 1) = xrefe(i, 1) + disp[i * 3 + 1];
          x(i, 2) = xrefe(i, 2) + disp[i * 3 + 2];
        }
        return;
      }

      //! Submethod to compute the enclosed volume for volume constraint boundary condition
      double ComputeConstrVols(
          const CORE::LINALG::SerialDenseMatrix& xc,  ///< current configuration
          const int numnode                           ///< num nodes
      );


      /*!
      \brief Submethod to compute interfacial area and its first and second
             derivatives w.r.t. the displacements (this is needed for surface energy problems)
      */
      void ComputeAreaDeriv(const CORE::LINALG::SerialDenseMatrix& x,  ///< spatial configuration
          const int numnode,                                           ///< number of nodes
          const int ndof,  ///< number of degrees of freedom
          double& A,       ///< area
          const Teuchos::RCP<CORE::LINALG::SerialDenseVector>& Adiff,  ///< first derivative
          const Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& Adiff2  ///< second derivative
      );

      //! Submethod to compute constraint volume and its first and second derivatives w.r.t. the
      //! displacements
      void ComputeVolDeriv(const CORE::LINALG::SerialDenseMatrix& x,  ///< spatial configuration
          const int numnode,                                          ///< number of nodes
          const int ndof,  ///< number of degrees of freedom
          double& V,       ///< volume
          const Teuchos::RCP<CORE::LINALG::SerialDenseVector>& Vdiff,   ///< first derivative
          const Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& Vdiff2,  ///< second derivative
          const int minind = 0,  ///< minimal index to compute enclosed volume with
          const int maxind = 2   ///< maximal index to compute enclosed volume with
      );

      /// Submethod to compute normal vectors at nodes (for reference configuration) which can be
      /// assembled (i.e., added) to obtain averaged nodal normal vectors
      void BuildNormalsAtNodes(CORE::LINALG::SerialDenseVector& nodenormals,
          const std::vector<double>& mydisp, bool refconfig);

      //! Submethod to compute surface porosity
      void CalculateSurfacePorosity(
          Teuchos::ParameterList& params, DRT::Discretization& discretization, LocationArray& la);

     private:
      //! Private Constructor to be called from StructuralSurfaceType
      explicit StructuralSurface(int id, int owner);

      //! Set Discretization Type of the Surface Element
      void SetDistype();

      //! Set Gaussrule dependent on shape of the structural surface
      void SetGaussrule();


      /// @name private methods for calculation of trace inequality matrices
      //@{

      //! Templated version: parent and surface discretization type
      template <CORE::FE::CellType dt_vol, CORE::FE::CellType dt_surf>
      double EstimateNitscheTraceMaxEigenvalueCombined(std::vector<double>& parent_disp);

      //! the volume stiffness matrix
      //! unlike the "full" stiffness matrix we don't use the geometric term here
      template <CORE::FE::CellType dt_vol>
      void TraceEstimateVolMatrix(const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xrefe,
          const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xcurr,
          CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol> * 3, CORE::FE::num_nodes<dt_vol> * 3>&
              vol);

      //! the surface stiffness matrix
      template <CORE::FE::CellType dt_vol, CORE::FE::CellType dt_surf>
      void TraceEstimateSurfMatrix(
          const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xrefe,
          const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xcurr,
          CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol> * 3, CORE::FE::num_nodes<dt_vol> * 3>&
              surf);

      //! evaluate the kinematics
      template <CORE::FE::CellType dt_vol>
      void Strains(const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xrefe,
          const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xcurr,
          const CORE::LINALG::Matrix<3, 1>& xi, double& jac, CORE::LINALG::Matrix<3, 3>& defgrd,
          CORE::LINALG::Matrix<6, 1>& glstrain, CORE::LINALG::Matrix<3, 3>& rcg,
          CORE::LINALG::Matrix<6, CORE::FE::num_nodes<dt_vol> * 3>& bop,
          CORE::LINALG::Matrix<3, CORE::FE::num_nodes<dt_vol>>& N_XYZ);


      //! setup projector for removing the rigid body modes from the generalized eigenvalue problem
      template <CORE::FE::CellType dt_vol>
      void SubspaceProjector(const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xcurr,
          CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol> * CORE::FE::dim<dt_vol>,
              CORE::FE::num_nodes<dt_vol> * CORE::FE::dim<dt_vol> -
                  CORE::FE::dim<dt_vol>*(CORE::FE::dim<dt_vol> + 1) / 2>& proj);


      //! Templated version: parent and surface discretization type
      template <CORE::FE::CellType dt_vol, CORE::FE::CellType dt_surf>
      double EstimateNitscheTraceMaxEigenvalueTSI(std::vector<double>& parent_disp);

      //! the volume stiffness matrix
      //! unlike the "full" stiffness matrix we don't use the geometric term here
      template <CORE::FE::CellType dt_vol>
      void TraceEstimateVolMatrixTSI(
          const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xrefe,
          const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xcurr,
          CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, CORE::FE::num_nodes<dt_vol>>& vol);

      //! the surface stiffness matrix
      template <CORE::FE::CellType dt_vol, CORE::FE::CellType dt_surf>
      void TraceEstimateSurfMatrixTSI(
          const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xrefe,
          const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xcurr,
          CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, CORE::FE::num_nodes<dt_vol>>& surf);


      //! setup projector for removing the rigid body modes from the generalized eigenvalue problem
      template <CORE::FE::CellType dt_vol>
      void SubspaceProjectorScalar(
          CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, CORE::FE::num_nodes<dt_vol> - 1>& proj);
      //@}

    };  // class StructuralSurface : public DRT::Element


  }  // namespace ELEMENTS
}  // namespace DRT


BACI_NAMESPACE_CLOSE

#endif  // SO3_SURFACE_H
