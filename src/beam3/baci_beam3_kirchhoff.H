/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief three dimensional nonlinear Kirchhoff beam element based on a C1 curve

\level 2

*/
/*-----------------------------------------------------------------------------------------------*/

/*
3D nonlinear Kirchhoff-like beam element. It can be switched between a variant with weak enforcement
of the Kirchhoff constraint and variant with strong enforcement of the Kirchhoff constraint. This
variant with weak constraint enforcement is based on a rotation interpolation that is similar to
beam3r. As the beam curve has to be C^1-continous, it is interpolated with Hermite polynomials of
order 3. Therefore each of the two boundary nodes has 7 dofs. With the flag rotvec_ is can be
switched between to sets of degrees of freedom on the boundary node. The first set (rotvec_==true)
is
[\v{d}_1, \v{theta}_1, t_1, \v{d}_2, \v{theta}_2, t_2, \alpha_3], where \v{d}_i is the vector of
nodal positions on the boundary nodes,\v{theta}_i is a pseudo rotation vector describing the nodal
triad orientation on the boundary nodes (and therewith also the orientation of the boundary tangent
vectors), t_i is the length of the boudary tangent vectors, and alpha_3 is the scalar relative
rotation angle between reference and matrial triad at the interior node. In contrary, the second
variant (rotvec_==false) of this element has the dofs:
[\v{d}_1, \v{t}_1, \alpha_1, \v{d}_2, \v{t}_2, \alpha_2, \alpha_3], where \v{t}_i is the nodal
tangent vector (orientation and length) at the boundary nodes, and alpha_i are the scalar relative
rotation angles between reference and matrial triad at the boundary nodes and the interior node.
Besides these two boundary nodes the element has BEAM3K_COLLOCATION_POINTS-2 interior nodes which
one scalar DoF alpha_i, respectively.


Attention: Since so far linearizations are calculated with FAD, the rotation increments in the case
(rotvec_==true) are of an additive nature, which is in strong contrast to the beam3r implementation,
where the iterative rotation increments are multiplicative. Consequently, the inhomogeneous
rotational Dirichlet conditions of beam3k can be interpreted as additive increments added to the
initial values (i.e. if the initial value is zero, the Dirichlet values in the input file are the
total nodal rotation angles). This is not true for beam3r, where prescribed 3D rotation values have
no direct physical interpretation. For 2D rotations both variants are identical.
*/

#ifndef BACI_BEAM3_KIRCHHOFF_H
#define BACI_BEAM3_KIRCHHOFF_H


// Todo @grill: get rid of header in header inclusions
#include "baci_beam3_base.H"
#include "baci_discretization_fem_general_largerotations.H"
#include "baci_lib_elementtype.H"
#include "baci_lib_node.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_linalg_serialdensevector.H"
#include "baci_utils_fad.H"

#include <Epetra_Vector.h>
#include <Sacado.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

typedef Sacado::Fad::DFad<double> FAD;

#define REFERENCE_NODE \
  2  // number of reference-node for the calculation of the rotation of the other triads (local
     // numbering: 0 2 3 4 1)
// DEFAULT: REFERENCE_NODE 2   //REFERENCE_NODE=2 represents the midpoint node in the case of 3
// BEAM3K_COLLOCATION_POINTS -> standard choice!

#if defined(REFERENCE_NODE) && (REFERENCE_NODE != 2)
dserror(
    "Beam3k REFERENCE_NODE: Only the value 2 is covered by tests and has therefore been "
    "cultivated in subsequent modifications to the code; carefully check correctness of code "
    "before using other values than 2!");
// note: especially look into the code for triad interpolation based on local rotation vectors
//       (as applied in WK case); reference triad according to Crisfield/Jelenic1999 is chosen as
//       triad at the 'mid'-node for uneven number of nodes
//       => REFERENCE_NODE=2 for BEAM3K_COLLOCATION_POINTS=3
//       generalize the class for triad interpolation if other choices of the reference nodes are
//       desired here in Beam3k
#endif

#define MYGAUSSRULEBEAM3K \
  ::CORE::DRT::UTILS::GaussRule1D::line_4point  // define gauss rule; intrule_line_1point -
                                                // intrule_line_10point is implemented.
// DEFAULT: intrule_line_4point

#define BEAM3K_COLLOCATION_POINTS \
  3  // defines type of element. 2,3,4 are supported. A value of 3 or 4 means that further inner
     // nodes are introduced
// DEFAULT: BEAM3K_COLLOCATION_POINTS 3   //to interpolate the torsional degree of freedom alpha_.
// Furthermore, it specifies the number of collocation points defining the number of material triads
// used to interpolate the triad field.

#if defined(BEAM3K_COLLOCATION_POINTS) && (BEAM3K_COLLOCATION_POINTS != 3)
dserror(
    "BEAM3K_COLLOCATION_POINTS: Only the value 3 is covered by tests and has therefore been "
    "cultivated in subsequent modifications to the code; carefully check correctness of code "
    "before using other values than 3!");
#endif


// #define CONSISTENTSPINSK        //Apply variationally consistent variant of first spin vector
//  component as test function for strong Kirchhoff
//-> this interpolation enables exact energy and angular momentum balances

// forward declaration
namespace LARGEROTATIONS
{
  template <unsigned int numnodes, typename T>
  class TriadInterpolationLocalRotationVectors;
}

namespace DRT
{
  namespace ELEMENTS
  {
    class Beam3kType : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "Beam3kType"; }

      static Beam3kType& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& actnode, const double* x0, const int numdof, const int dimnsp) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static Beam3kType instance_;
    };

    /*!
    \brief 3D nonlinear Kirchhoff-like beam element that can display initially curved beams.

    */
    class Beam3k : public Beam3Base
    {
     public:
      //! @name Friends
      friend class Beam3kType;

      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id    (in): A globally unique element id
      \param etype (in): Type of element
      \param owner (in): owner processor of the element
      */
      Beam3k(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element
      */
      Beam3k(const Beam3k& old);



      /*!
      \brief Deep copy this instance of Beam3k and return pointer to the copy

      The Clone() method is used by the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed
    .
      */
      DRT::Element* Clone() const override;

      /*!
     \brief Get shape type of element
     */
      DiscretizationType Shape() const override;


      /*!
      \brief Return unique ParObject id

      Every class implementing ParObject needs a unique id defined at the
      top of parobject.H
      */
      int UniqueParObjectId() const override { return Beam3kType::Instance().UniqueParObjectId(); }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(DRT::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Unpack(const std::vector<char>& data) override;

      DRT::ElementType& ElementType() const override { return Beam3kType::Instance(); }

      //@}

      /*!
       \brief Get bool inidcating usage of rotation vectors
       */
      const bool& RotVec() const { return rotvec_; }

      /*!
      \brief get reference rotation vector i.e. theta0_
      */
      std::vector<CORE::LINALG::Matrix<3, 1>> Theta0() const { return theta0_; }

      /*!
       \brief Get (non-unit) tangent vectors at the two boundary nodes
      */
      std::vector<CORE::LINALG::Matrix<3, 1>> GetNodalTangents() const { return T_; }

      /** \brief get unit tangent vector in reference configuration at i-th node of beam element
       * (element-internal numbering)
       *
       *  \author grill
       *  \date 06/16 */
      inline void GetRefTangentAtNode(
          CORE::LINALG::Matrix<3, 1>& Tref_i, const int& i) const override
      {
        if (not((unsigned)i < Tref().size()))
          dserror("asked for tangent at node index %d, but only %d centerline nodes existing", i,
              Tref().size());
        Tref_i = Tref()[i];
      }

      /** \brief get centerline position at xi \in [-1,1] (element parameter space)
       *
       *  \author grill
       *  \date 06/16 */
      void GetPosAtXi(CORE::LINALG::Matrix<3, 1>& pos, const double& xi,
          const std::vector<double>& disp) const override;

      /** \brief get triad at xi \in [-1,1] (element parameter space)
       *
       *  \author grill
       *  \date 01/17 */
      void GetTriadAtXi(CORE::LINALG::Matrix<3, 3>& triad, const double& xi,
          const std::vector<double>& disp) const override;

      /** \brief Get scaled base vectors describing the cross-section orientation and size at given
       * parameter coordinate xi
       *
       *  note: this method is only used for visualization so far and limited to a rectangular (?)
       *        cross-section shape; the length of the base vectors indicates the size of the
       * cross-section in the direction of the base vector
       *
       *  \author meier
       *  \date 03/16 */
      void GetScaledSecondAndThirdBaseVectorAtXi(const double& xi, const std::vector<double>& disp,
          CORE::LINALG::Matrix<3, 2>& scaledbasevectors) const;

      /** \brief get generalized interpolation matrix which yields the variation of the position and
       *         orientation at xi \in [-1,1] if multiplied with the vector of primary DoF
       * variations
       *
       *  \author grill
       *  \date 01/17 */
      void GetGeneralizedInterpolationMatrixVariationsAtXi(CORE::LINALG::SerialDenseMatrix& Ivar,
          const double& xi, const std::vector<double>& disp) const override;

      /** \brief get linearization of the product of (generalized interpolation matrix for
       * variations (see above) and applied force vector) with respect to the primary DoFs of this
       * element
       *
       *  \author grill
       *  \date 01/17 */
      void GetStiffmatResultingFromGeneralizedInterpolationMatrixAtXi(
          CORE::LINALG::SerialDenseMatrix& stiffmat, const double& xi,
          const std::vector<double>& disp,
          const CORE::LINALG::SerialDenseVector& force) const override;

      /** \brief get generalized interpolation matrix which yields the increments of the position
       * and orientation at xi \in [-1,1] if multiplied with the vector of primary DoF increments
       *
       *  \author grill
       *  \date 01/17 */
      void GetGeneralizedInterpolationMatrixIncrementsAtXi(CORE::LINALG::SerialDenseMatrix& Iinc,
          const double& xi, const std::vector<double>& disp) const override;

      /** \brief get access to the reference length
       *
       *  \author grill
       *  \date 05/16 */
      inline double RefLength() const override { return length_; }

      /*!
      \brief Get jacobi factor of first Gauss point
      */
      const double& GetJacobi() const { return jacobi_[0]; }

      /** \brief get Jacobi factor ds/dxi(xi) at xi \in [-1;1]
       *
       *  \author grill
       *  \date 06/16 */
      double GetJacobiFacAtXi(const double& xi) const override;

      /** \brief Get material cross-section deformation measures, i.e. strain resultants
       *
       *  \author grill
       *  \date 04/17 */
      inline void GetMaterialStrainResultantsAtAllGPs(std::vector<double>& axial_strain_GPs,
          std::vector<double>& shear_strain_2_GPs, std::vector<double>& shear_strain_3_GPs,
          std::vector<double>& twist_GPs, std::vector<double>& curvature_2_GPs,
          std::vector<double>& curvature_3_GPs) const override
      {
        axial_strain_GPs = axial_strain_GP_;
        // note: shear deformations are zero by definition for Kirchhoff beam formulation
        shear_strain_2_GPs.clear();
        shear_strain_3_GPs.clear();

        twist_GPs = twist_GP_;
        curvature_2_GPs = curvature_2_GP_;
        curvature_3_GPs = curvature_3_GP_;
      }

      /** \brief Get material cross-section stress resultants
       *
       *  \author grill
       *  \date 04/17 */
      inline void GetMaterialStressResultantsAtAllGPs(std::vector<double>& axial_force_GPs,
          std::vector<double>& shear_force_2_GPs, std::vector<double>& shear_force_3_GPs,
          std::vector<double>& torque_GPs, std::vector<double>& bending_moment_2_GPs,
          std::vector<double>& bending_moment_3_GPs) const override
      {
        axial_force_GPs = axial_force_GP_;
        // note: shear deformations are zero by definition for Kirchhoff beam formulation
        shear_force_2_GPs.clear();
        shear_force_3_GPs.clear();

        torque_GPs = torque_GP_;
        bending_moment_2_GPs = bending_moment_2_GP_;
        bending_moment_3_GPs = bending_moment_3_GP_;
      }

      //! get internal (elastic) energy of element
      double GetInternalEnergy() const override { return Eint_; };

      //! get kinetic energy of element
      double GetKineticEnergy() const override { return Ekin_; };

      /** \brief get number of nodes used for centerline interpolation
       *
       *  \author grill
       *  \date 05/16 */
      inline int NumCenterlineNodes() const override { return 2; }

      /** \brief find out whether given node is used for centerline interpolation
       *
       *  \author grill
       *  \date 10/16 */
      inline bool IsCenterlineNode(const DRT::Node& node) const override
      {
        if (node.Id() == this->Nodes()[0]->Id() or node.Id() == this->Nodes()[1]->Id())
          return true;
        else
          return false;
      }

      /*!
      \brief Return number of lines to this element
      */
      int NumLine() const override { return 1; }


      /*!
      \brief Get vector of RCPs to the lines of this element
      */
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;


      /*!
      \brief Get number of degrees of freedom of a single node
      */
      int NumDofPerNode(const DRT::Node& node) const override
      {
        /*note: this is not necessarily the number of DOF assigned to this node by the
         *discretization finally, but only the number of DOF requested for this node by this
         *element; the discretization will finally assign the maximal number of DOF to this node
         *requested by any element connected to this node*/
        int dofpn = 0;

        if (node.Id() == this->Nodes()[0]->Id() or node.Id() == this->Nodes()[1]->Id())
          dofpn = 7;
        else
          dofpn = 1;

        return dofpn;
      }

      /*!
      \brief Get number of degrees of freedom per element not including nodal degrees of freedom
      */
      int NumDofPerElement() const override { return 0; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      //@}

      //! @name Construction


      /*!
      \brief Read input for this element

      This class implements a dummy of this method that prints a warning and
      returns false. A derived class would read one line from the input file and
      store all necessary information.

      */
      // virtual bool ReadElement();

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef) override;

      //@}


      //! @name Evaluation methods

      /*!
      \brief Evaluate an element

      An element derived from this class uses the Evaluate method to receive commands
      and parameters from some control routine in params and evaluates element matrices and
      vectors accoring to the command in params.

      \note This class implements a dummy of this method that prints a warning and
            returns false.

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param lm (in)            : location vector of this element
      \param elemat1 (out)      : matrix to be filled by element depending on commands
                                  given in params
      \param elemat2 (out)      : matrix to be filled by element depending on commands
                                  given in params
      \param elevec1 (out)      : vector to be filled by element depending on commands
                                  given in params
      \param elevec2 (out)      : vector to be filled by element depending on commands
                                  given in params
      \param elevec3 (out)      : vector to be filled by element depending on commands
                                  given in params
      \return 0 if successful, negative otherwise
      */
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override;

      /*!
      \brief Evaluate a Neumann boundary condition

      An element derived from this class uses the EvaluateNeumann method to receive commands
      and parameters from some control routine in params and evaluates a Neumann boundary condition
      given in condition

      \note This class implements a dummy of this method that prints a warning and
            returns false.

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param condition (in)     : The condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : Force vector to be filled by element

      \return 0 if successful, negative otherwise

      \author meier
      \date 01/16 */
      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1) override;

      /*!
      \brief Set rotations in the initial configuration

       Set the initial rotations based on nodal rotation (pseudo) vectors.
       The nodal rotation vectors are independent of the subsequent used
       rotational interpolation method.
      */
      void SetUpInitialRotations(const std::vector<double>& nodal_thetas);


      //! sets up from current nodal position all geometric parameters (considering current position
      //! as reference configuration)
      void SetUpReferenceGeometry(
          const std::vector<CORE::LINALG::Matrix<3, 1>>& xrefe, const bool secondinit = false);

      //! computes the artifical damping contributions for element based PTC
      void CalcStiffContributionsPTC(CORE::LINALG::SerialDenseMatrix& elemat1);

      /*!
      \brief get (non-unit) tangent vectors at the two boundary nodes in the initial configuration
      */
      std::vector<CORE::LINALG::Matrix<3, 1>> Tref() const { return Tref_; }
      //@}

      /** \brief add indices of those DOFs of a given node that are positions
       *
       *  \author grill
       *  \date 07/16 */
      inline void PositionDofIndices(
          std::vector<int>& posdofs, const DRT::Node& node) const override
      {
        if (node.Id() == this->Nodes()[0]->Id() or node.Id() == this->Nodes()[1]->Id())
        {
          posdofs.push_back(0);
          posdofs.push_back(1);
          posdofs.push_back(2);
        }
        return;
      }

      /** \brief add indices of those DOFs of a given node that are tangents (in the case of Hermite
       * interpolation)
       *
       *  \author grill
       *  \date 07/16 */
      inline void TangentDofIndices(
          std::vector<int>& tangdofs, const DRT::Node& node) const override
      {
        if ((not rotvec_) and
            (node.Id() == this->Nodes()[0]->Id() or node.Id() == this->Nodes()[1]->Id()))
        {
          tangdofs.push_back(3);
          tangdofs.push_back(4);
          tangdofs.push_back(5);
        }
        return;
      }

      /** \brief add indices of those DOFs of a given node that are rotation DOFs (non-additive
       * rotation vectors)
       *
       *  \author grill
       *  \date 07/16 */
      inline void RotationVecDofIndices(
          std::vector<int>& rotvecdofs, const DRT::Node& node) const override
      {
        if ((rotvec_) and
            (node.Id() == this->Nodes()[0]->Id() or node.Id() == this->Nodes()[1]->Id()))
        {
          rotvecdofs.push_back(3);
          rotvecdofs.push_back(4);
          rotvecdofs.push_back(5);
        }
        return;
      }

      /** \brief add indices of those DOFs of a given node that are 1D rotation DOFs
       *         (planar rotations are additive, e.g. in case of relative twist DOF of beam3k with
       * rotvec=false)
       *
       *  \author grill
       *  \date 07/16 */
      inline void Rotation1DDofIndices(
          std::vector<int>& twistdofs, const DRT::Node& node) const override
      {
        if ((not rotvec_) and
            (node.Id() == this->Nodes()[0]->Id() or node.Id() == this->Nodes()[1]->Id()))
        {
          twistdofs.push_back(6);
        }
        else if (node.Id() == this->Nodes()[2]->Id())
        {
          twistdofs.push_back(0);
        }

        return;
      }

      /** \brief add indices of those DOFs of a given node that represent norm of tangent vector
       *         (additive, e.g. in case of beam3k with rotvec=true)
       *
       *  \author grill
       *  \date 07/16 */
      inline void TangentLengthDofIndices(
          std::vector<int>& tangnormdofs, const DRT::Node& node) const override
      {
        if (rotvec_)
        {
          if (node.Id() == this->Nodes()[0]->Id() or node.Id() == this->Nodes()[1]->Id())
          {
            tangnormdofs.push_back(6);
          }
        }
        return;
      }

      /** \brief get element local indices of those Dofs that are used for centerline interpolation
       *
       *  \author grill
       *  \date 12/16 */
      inline void CenterlineDofIndicesOfElement(
          std::vector<unsigned int>& centerlinedofindices) const override
      {
        if (rotvec_)
          dserror(
              "The logic of this implementation does not apply for Beam3k with rotation vector "
              "Dofs! "
              "Be careful and find a solution to convert force/stiffness contributions to tangent "
              "Dofs "
              "automatically and apply them consistently to rotvec Dofs");

        const unsigned int nnodecl = this->NumCenterlineNodes();
        centerlinedofindices.resize(6 * nnodecl, 0);

        for (unsigned int inodecl = 0; inodecl < nnodecl; ++inodecl)
          for (unsigned int idof = 0; idof < 6; ++idof)
            centerlinedofindices[6 * inodecl + idof] = 7 * inodecl + idof;
      }

      /** \brief extract values for those Dofs relevant for centerline-interpolation from total
       * state vector
       *
       *  \author grill
       *  \date 11/16 */
      void ExtractCenterlineDofValuesFromElementStateVector(const std::vector<double>& dofvec,
          std::vector<double>& dofvec_centerline, bool add_reference_values = false) const override;

      //! computes the number of different random numbers required in each time step for generation
      //! of stochastic forces
      int HowManyRandomNumbersINeed() const override;

     private:
      // Functions

      //! @name methods for initilization of the element

      //! sets up from current nodal position all geometric parameters (considering current position
      //! as reference configuration) in case of a Weak Kirchhoff Constraint
      void SetUpReferenceGeometryWK(
          const std::vector<CORE::LINALG::Matrix<3, 1>>& xrefe, const bool secondinit);

      //! sets up from current nodal position all geometric parameters (considering current position
      //! as reference configuration) in case of a Strong Kirchhoff Constraint
      void SetUpReferenceGeometrySK(
          const std::vector<CORE::LINALG::Matrix<3, 1>>& xrefe, const bool secondinit);

      //@}

      //! @name auxiliary computation methods for non-additive, large rotation variables

      //@}

      //! @name Internal calculation methods

      /** \brief Calculate internal forces and stiffness matrix
       *
       *  \author grill
       *  \date 01/17 */
      void CalcInternalAndInertiaForcesAndStiff(Teuchos::ParameterList& params,
          std::vector<double>& disp, CORE::LINALG::SerialDenseMatrix* stiffmatrix,
          CORE::LINALG::SerialDenseMatrix* massmatrix, CORE::LINALG::SerialDenseVector* force,
          CORE::LINALG::SerialDenseVector* force_inert);

      /** \brief Calculate internal forces and stiffness matrix in case of a Weak Kirchhoff
       * Constraint
       *
       *  \author meier
       *  \date 01/16 */
      template <unsigned int nnodecl, typename T>
      void CalculateInternalForcesAndStiffWK(Teuchos::ParameterList& params,
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>&
              disp_totlag_centerline,
          const std::vector<CORE::LINALG::Matrix<3, 3, T>>& triad_mat_cp,
          CORE::LINALG::SerialDenseMatrix* stiffmatrix,
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& internal_force,
          std::vector<CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, T>>&
              v_theta_gp,
          std::vector<CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, T>>&
              lin_theta_gp,
          std::vector<CORE::LINALG::Matrix<3, 3, T>>& triad_mat_gp);

      /** \brief Calculate internal forces and stiffness matrix in case of a Strong Kirchhoff
       * Constraint
       *
       *  \author meier
       *  \date 02/16 */
      template <unsigned int nnodecl>
      void CalculateInternalForcesAndStiffSK(Teuchos::ParameterList& params,
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, FAD>&
              disp_totlag_centerline,
          const std::vector<CORE::LINALG::Matrix<3, 3, FAD>>& triad_mat_cp,
          CORE::LINALG::SerialDenseMatrix* stiffmatrix,
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, FAD>& internal_force,
          std::vector<CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, FAD>>&
              v_theta_gp,
          std::vector<CORE::LINALG::Matrix<3, 3, FAD>>& triad_mat_gp);

      /** \brief Calculate contributions to the stiffness matrix at a Gauss point analytically
       *         in case of weak Kirchhoff constraint
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void CalculateStiffmatContributionsAnalyticWK(CORE::LINALG::SerialDenseMatrix& stiffmatrix,
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, double>&
              disp_totlag_centerline,
          const LARGEROTATIONS::TriadInterpolationLocalRotationVectors<BEAM3K_COLLOCATION_POINTS,
              double>& triad_intpol,
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, double>&
              v_theta_s_bar,
          const std::vector<CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
              double>>& lin_theta_cp,
          CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_theta_bar,
          const std::vector<CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
              6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>>& lin_v_epsilon_cp,
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, double>&
              v_epsilon_bar,
          double axial_force_bar, const CORE::LINALG::Matrix<3, 1, double>& moment_resultant,
          double axial_rigidity,
          const CORE::LINALG::Matrix<3, 3, double>& constitutive_matrix_moment_material,
          const CORE::LINALG::Matrix<3, 1, double>& theta_gp,
          const CORE::LINALG::Matrix<3, 1, double>& theta_s_gp,
          const CORE::LINALG::Matrix<3, 3, double>& triad_mat_gp, double xi_gp, double jacobifac_gp,
          double GPwgt) const;

      template <unsigned int nnodecl>
      void CalculateStiffmatContributionsAnalyticWK(CORE::LINALG::SerialDenseMatrix& stiffmatrix,
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, FAD>&
              disp_totlag_centerline,
          const LARGEROTATIONS::TriadInterpolationLocalRotationVectors<BEAM3K_COLLOCATION_POINTS,
              FAD>& triad_intpol,
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, FAD>&
              v_theta_s_bar,
          const std::vector<CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, FAD>>&
              lin_theta_cp,
          CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, FAD>& lin_theta_bar,
          const std::vector<CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
              6 * nnodecl + BEAM3K_COLLOCATION_POINTS, FAD>>& lin_v_epsilon_cp,
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, FAD>&
              v_epsilon_bar,
          FAD axial_force_bar, const CORE::LINALG::Matrix<3, 1, FAD>& moment_resultant,
          FAD axial_rigidity,
          const CORE::LINALG::Matrix<3, 3, FAD>& constitutive_matrix_moment_material,
          const CORE::LINALG::Matrix<3, 1, FAD>& theta_gp,
          const CORE::LINALG::Matrix<3, 1, FAD>& theta_s_gp,
          const CORE::LINALG::Matrix<3, 3, FAD>& triad_mat_gp, double xi_gp, double jacobifac_gp,
          double GPwgt) const
      {
        // this is a dummy because in case that the pre-calculated values are of type Fad,
        // we use automatic differentiation and consequently there is no need for analytic stiffmat
      }

      /** \brief pre-compute quantities required for analytic computation of stiffness matrix
       *         in case of weak Kirchhoff constraint
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void PreComputeTermsAtCPForStiffmatContributionsAnalyticWK(
          CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_theta,
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
              6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_v_epsilon,
          const CORE::LINALG::Matrix<1, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& L,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
          const CORE::LINALG::Matrix<3, 1, double>& r_s, double abs_r_s,
          const CORE::LINALG::Matrix<4, 1, double>& Qref_conv) const;

      template <unsigned int nnodecl>
      void PreComputeTermsAtCPForStiffmatContributionsAnalyticWK(
          CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, FAD>& lin_theta,
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
              6 * nnodecl + BEAM3K_COLLOCATION_POINTS, FAD>& lin_v_epsilon,
          const CORE::LINALG::Matrix<1, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, FAD>& L,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, FAD>& N_s,
          const CORE::LINALG::Matrix<3, 1, FAD>& r_s, FAD abs_r_s,
          const CORE::LINALG::Matrix<4, 1, double>& Qref_conv) const
      {
        // this is empty because in case that the pre-calculated values are of type Fad,
        // we use automatic differentiation and consequently there is no need for analytic stiffmat
      }


      /** \brief Calculate inertia forces and mass matrix
       *
       *  \author meier
       *  \date 02/16 */
      template <unsigned int nnodecl, typename T>
      void CalculateInertiaForcesAndMassMatrix(Teuchos::ParameterList& params,
          const std::vector<CORE::LINALG::Matrix<3, 3, T>>& triad_mat_gp,
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>&
              disp_totlag_centerline,
          const std::vector<CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, T>>&
              v_theta_gp,
          const std::vector<CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, T>>&
              lin_theta_gp,
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& f_inert,
          CORE::LINALG::SerialDenseMatrix* massmatrix);

      /** \brief Calculate analytic linearization of inertia forces, i.e. mass matrix
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void CalculateMassMatrixContributionsAnalyticWK(CORE::LINALG::SerialDenseMatrix& massmatrix,
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, double>&
              disp_totlag_centerline,
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, double>&
              v_theta_bar,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>&
              lin_theta_bar,
          const CORE::LINALG::Matrix<3, 1, double>& moment_rho,
          const CORE::LINALG::Matrix<3, 1, double>& deltatheta,
          const CORE::LINALG::Matrix<3, 1, double>& angular_velocity_material,
          const CORE::LINALG::Matrix<3, 3, double>& triad_mat,
          const CORE::LINALG::Matrix<3, 3, double>& triad_mat_conv,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N,
          double mass_inertia_translational,
          const CORE::LINALG::Matrix<3, 3, double>& tensor_mass_moment_of_inertia,
          double lin_prefactor_acc, double lin_prefactor_vel, double xi_gp, double jacobifac_gp,
          double GPwgt) const;

      template <unsigned int nnodecl>
      void CalculateMassMatrixContributionsAnalyticWK(CORE::LINALG::SerialDenseMatrix& massmatrix,
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, FAD>&
              disp_totlag_centerline,
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, FAD>& v_theta_bar,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, FAD>&
              lin_theta_bar,
          const CORE::LINALG::Matrix<3, 1, FAD>& moment_rho,
          const CORE::LINALG::Matrix<3, 1, FAD>& deltatheta,
          const CORE::LINALG::Matrix<3, 1, FAD>& angular_velocity_material,
          const CORE::LINALG::Matrix<3, 3, FAD>& triad_mat,
          const CORE::LINALG::Matrix<3, 3, FAD>& triad_mat_conv,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, FAD>& N,
          double density, const CORE::LINALG::Matrix<3, 3, FAD>& tensor_mass_moment_of_inertia,
          double lin_prefactor_acc, double lin_prefactor_vel, double xi_gp, double jacobifac_gp,
          double GPwgt) const
      {
        // this is empty because in case that the pre-calculated values are of type Fad,
        // we use automatic differentiation and consequently there is no need for analytic mass
        // matrix
      }

      /** \brief evaluate contributions to element residual vector and stiffmat from point Neumann
       *         condition
       *         note: we need to evaluate this on element level because point moments need to be
       *         linearized in case of tangent-based formulation (rotvec_=false)
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void EvaluatePointNeumannEB(CORE::LINALG::SerialDenseVector& forcevec,
          CORE::LINALG::SerialDenseMatrix* stiffmat,
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, double>&
              disp_totlag,
          const CORE::LINALG::Matrix<6, 1, double>& load_vector_neumann, int node) const;

      /** \brief evaluate contributions to element residual vector from point Neumann moment
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl, typename T>
      void EvaluateResidualFromPointNeumannMoment(
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& force_ext,
          const CORE::LINALG::Matrix<3, 1, T>& moment_ext, const CORE::LINALG::Matrix<3, 1, T>& r_s,
          T abs_r_s, int node) const;

      /** \brief evaluate contributions to element stiffness matrix from point Neumann moment
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void EvaluateStiffMatrixAnalyticFromPointNeumannMoment(
          CORE::LINALG::SerialDenseMatrix& stiffmat,
          const CORE::LINALG::Matrix<3, 1, double>& moment_ext,
          const CORE::LINALG::Matrix<3, 1, double>& r_s, double abs_r_s, int node) const;

      /** \brief evaluate contributions to element residual vector and stiffmat from line Neumann
       *         condition
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void EvaluateLineNeumann(CORE::LINALG::SerialDenseVector& forcevec,
          CORE::LINALG::SerialDenseMatrix* stiffmat,
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, double>&
              disp_totlag,
          const CORE::LINALG::Matrix<6, 1, double>& load_vector_neumann,
          const std::vector<int>* function_numbers, double time) const;

      /** \brief evaluate contributions to element residual vector from line Neumann condition
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl, typename T>
      void EvaluateLineNeumannForces(
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& force_ext,
          const CORE::LINALG::Matrix<6, 1, double>& load_vector_neumann,
          const std::vector<int>* function_numbers, double time) const;

      /** \brief evaluate all contributions from brownian dynamics (thermal & viscous
       * forces/moments)
       *
       *  \author grill
       *  \date 12/16 */
      template <unsigned int nnode, unsigned int vpernode, unsigned int ndim>
      void CalcBrownianForcesAndStiff(Teuchos::ParameterList& params,
          std::vector<double>& vel,                      //!< element velocity vector
          std::vector<double>& disp,                     //!< element displacement vector
          CORE::LINALG::SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
          CORE::LINALG::SerialDenseVector* force);       //!< element internal force vector

      /** \brief evaluate all contributions from translational damping forces
       *
       *  \author grill
       *  \date 12/16 */
      template <typename T, unsigned int nnode, unsigned int vpernode, unsigned int ndim>
      void EvaluateTranslationalDamping(Teuchos::ParameterList& params,  //!< parameter list
          const CORE::LINALG::Matrix<ndim * vpernode * nnode, 1, double>& vel,
          const CORE::LINALG::Matrix<ndim * vpernode * nnode, 1, T>& disp_totlag,
          CORE::LINALG::SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
          CORE::LINALG::Matrix<ndim * vpernode * nnode + BEAM3K_COLLOCATION_POINTS, 1, T>&
              f_int);  //!< element internal force vector

      /** \brief evaluate contributions to element stiffness matrix from translational damping
       * forces
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnode, unsigned int vpernode, unsigned int ndim>
      void EvaluateAnalyticStiffmatContributionsFromTranslationalDamping(
          CORE::LINALG::SerialDenseMatrix& stiffmatrix,
          const CORE::LINALG::Matrix<ndim, ndim, double>& damping_matrix,
          const CORE::LINALG::Matrix<ndim, 1, double>& r_s,
          const CORE::LINALG::Matrix<ndim, 1, double>& vel_rel,
          const CORE::LINALG::Matrix<ndim, 1, double>& gamma,
          const CORE::LINALG::Matrix<ndim, ndim, double>& velbackgroundgrad,
          const CORE::LINALG::Matrix<1, nnode * vpernode, double>& N_i,
          const CORE::LINALG::Matrix<1, nnode * vpernode, double>& N_i_xi, double jacobifactor,
          double gp_weight) const;

      template <unsigned int nnode, unsigned int vpernode, unsigned int ndim>
      void EvaluateAnalyticStiffmatContributionsFromTranslationalDamping(
          CORE::LINALG::SerialDenseMatrix& stiffmatrix,
          const CORE::LINALG::Matrix<ndim, ndim, FAD>& damping_matrix,
          const CORE::LINALG::Matrix<ndim, 1, FAD>& r_s,
          const CORE::LINALG::Matrix<ndim, 1, FAD>& vel_rel,
          const CORE::LINALG::Matrix<ndim, 1, double>& gamma,
          const CORE::LINALG::Matrix<ndim, ndim, FAD>& velbackgroundgrad,
          const CORE::LINALG::Matrix<1, nnode * vpernode, double>& N_i,
          const CORE::LINALG::Matrix<1, nnode * vpernode, double>& N_i_xi, double jacobifactor,
          double gp_weight) const
      {
        // this is empty because in case that the pre-calculated values are of type Fad,
        // we use automatic differentiation and consequently there is no need for analytic stiffmat
      }

      /** \brief evaluate all contributions from thermal/stochastic forces
       *
       *  \author grill
       *  \date 12/16 */
      template <typename T, unsigned int nnode, unsigned int vpernode, unsigned int ndim,
          unsigned int randompergauss>
      void EvaluateStochasticForces(
          const CORE::LINALG::Matrix<ndim * vpernode * nnode, 1, T>& disp_totlag,
          CORE::LINALG::SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
          CORE::LINALG::Matrix<ndim * vpernode * nnode + BEAM3K_COLLOCATION_POINTS, 1, T>&
              f_int);  //!< element internal force vector


      /** \brief evaluate contributions to element stiffness matrix from thermal/stochastic forces
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnode, unsigned int vpernode, unsigned int ndim>
      void EvaluateAnalyticStiffmatContributionsFromStochasticForces(
          CORE::LINALG::SerialDenseMatrix& stiffmatrix,
          const CORE::LINALG::Matrix<ndim, 1, double>& r_s,
          const CORE::LINALG::Matrix<ndim, 1, double>& randnumvec,
          const CORE::LINALG::Matrix<ndim, 1, double>& gamma,
          const CORE::LINALG::Matrix<1, nnode * vpernode, double>& N_i,
          const CORE::LINALG::Matrix<1, nnode * vpernode, double>& N_i_xi, double jacobifactor,
          double gp_weight) const;

      template <unsigned int nnode, unsigned int vpernode, unsigned int ndim>
      void EvaluateAnalyticStiffmatContributionsFromStochasticForces(
          CORE::LINALG::SerialDenseMatrix& stiffmatrix,
          const CORE::LINALG::Matrix<ndim, 1, FAD>& r_s,
          const CORE::LINALG::Matrix<ndim, 1, double>& randnumvec,
          const CORE::LINALG::Matrix<ndim, 1, double>& gamma,
          const CORE::LINALG::Matrix<1, nnode * vpernode, double>& N_i,
          const CORE::LINALG::Matrix<1, nnode * vpernode, double>& N_i_xi, double jacobifactor,
          double gp_weight) const
      {
        // this is empty because in case that the pre-calculated values are of type Fad,
        // we use automatic differentiation and consequently there is no need for analytic stiffmat
      }

      /** \brief evaluate all contributions from rotational damping moment/torque
       *
       *  \author grill
       *  \date 12/16 */
      template <typename T, unsigned int nnode, unsigned int vpernode, unsigned int ndim>
      void EvaluateRotationalDamping(
          const CORE::LINALG::Matrix<ndim * vpernode * nnode + BEAM3K_COLLOCATION_POINTS, 1, T>&
              disp_totlag_centerline,
          const std::vector<CORE::LINALG::Matrix<ndim, ndim, T>>& triad_mat_cp,
          CORE::LINALG::SerialDenseMatrix* stiffmatrix,
          CORE::LINALG::Matrix<ndim * vpernode * nnode + BEAM3K_COLLOCATION_POINTS, 1, T>& f_int);

      /** \brief evaluate contributions to element stiffness matrix from rotational damping moment
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl, unsigned int vpernode, unsigned int ndim>
      void EvaluateAnalyticStiffmatContributionsFromRotationalDamping(
          CORE::LINALG::SerialDenseMatrix& stiffmatrix,
          const CORE::LINALG::Matrix<ndim * vpernode * nnodecl + BEAM3K_COLLOCATION_POINTS, 1,
              double>& disp_totlag_centerline,
          const LARGEROTATIONS::TriadInterpolationLocalRotationVectors<BEAM3K_COLLOCATION_POINTS,
              double>& triad_intpol,
          const CORE::LINALG::Matrix<3, 1, double> theta_gp,
          const CORE::LINALG::Matrix<3, 1, double>& deltatheta_gp,
          const CORE::LINALG::Matrix<3, 3, double>& triad_mat_gp,
          const CORE::LINALG::Matrix<3, 3, double>& triad_mat_conv_gp,
          const CORE::LINALG::Matrix<ndim * vpernode * nnodecl + BEAM3K_COLLOCATION_POINTS, ndim,
              double>& v_theta_par_bar,
          const std::vector<CORE::LINALG::Matrix<ndim,
              ndim * vpernode * nnodecl + BEAM3K_COLLOCATION_POINTS, double>>& lin_theta_cp,
          const CORE::LINALG::Matrix<3, 1, double> moment_viscous, double gamma_polar, double dt,
          double xi_gp, double jacobifac_GPwgt) const;

      template <unsigned int nnodecl, unsigned int vpernode, unsigned int ndim>
      void EvaluateAnalyticStiffmatContributionsFromRotationalDamping(
          CORE::LINALG::SerialDenseMatrix& stiffmatrix,
          const CORE::LINALG::Matrix<ndim * vpernode * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, FAD>&
              disp_totlag_centerline,
          const LARGEROTATIONS::TriadInterpolationLocalRotationVectors<BEAM3K_COLLOCATION_POINTS,
              FAD>& triad_intpol,
          const CORE::LINALG::Matrix<3, 1, FAD> theta_gp,
          const CORE::LINALG::Matrix<3, 1, FAD>& deltatheta_gp,
          const CORE::LINALG::Matrix<3, 3, FAD>& triad_mat_gp,
          const CORE::LINALG::Matrix<3, 3, FAD>& triad_mat_conv_gp,
          const CORE::LINALG::Matrix<ndim * vpernode * nnodecl + BEAM3K_COLLOCATION_POINTS, ndim,
              FAD>& v_theta_par_bar,
          const std::vector<CORE::LINALG::Matrix<ndim,
              ndim * vpernode * nnodecl + BEAM3K_COLLOCATION_POINTS, FAD>>& lin_theta_cp,
          const CORE::LINALG::Matrix<3, 1, FAD> moment_viscous, double gamma_polar, double dt,
          double xi_gp, double jacobifac_GPwgt) const
      {
        // this is empty because in case that the pre-calculated values are of type Fad,
        // we use automatic differentiation and consequently there is no need for analytic stiffmat
      }

      /** \brief pre-compute quantities required for linearization of rotational damping moment
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnode, unsigned int vpernode, unsigned int ndim>
      void PreComputeTermsAtCPForAnalyticStiffmatContributionsFromRotationalDamping(
          CORE::LINALG::Matrix<ndim, ndim * vpernode * nnode + BEAM3K_COLLOCATION_POINTS, double>&
              lin_theta,
          const CORE::LINALG::Matrix<1, ndim * vpernode * nnode + BEAM3K_COLLOCATION_POINTS,
              double>& L,
          const CORE::LINALG::Matrix<ndim, ndim * vpernode * nnode + BEAM3K_COLLOCATION_POINTS,
              double>& N_s,
          const CORE::LINALG::Matrix<ndim, 1, double>& r_s, double abs_r_s,
          const CORE::LINALG::Matrix<4, 1, double>& Qref_conv) const;

      template <unsigned int nnode, unsigned int vpernode, unsigned int ndim>
      void PreComputeTermsAtCPForAnalyticStiffmatContributionsFromRotationalDamping(
          CORE::LINALG::Matrix<ndim, ndim * vpernode * nnode + BEAM3K_COLLOCATION_POINTS, FAD>&
              lin_theta,
          const CORE::LINALG::Matrix<1, ndim * vpernode * nnode + BEAM3K_COLLOCATION_POINTS, FAD>&
              L,
          const CORE::LINALG::Matrix<ndim, ndim * vpernode * nnode + BEAM3K_COLLOCATION_POINTS,
              FAD>& N_s,
          const CORE::LINALG::Matrix<ndim, 1, FAD>& r_s, FAD abs_r_s,
          const CORE::LINALG::Matrix<4, 1, double>& Qref_conv) const
      {
        // this is empty because in case that the pre-calculated values are of type Fad,
        // we use automatic differentiation and consequently there is no need for analytic stiffmat
      }

      //! compute (material) strain K
      template <typename T>
      void computestrain(const CORE::LINALG::Matrix<3, 1, T>& theta,
          const CORE::LINALG::Matrix<3, 1, T>& theta_deriv, CORE::LINALG::Matrix<3, 1, T>& K) const
      {
        CORE::LINALG::Matrix<3, 3, T> Tinv = CORE::LARGEROTATIONS::Tinvmatrix(theta);

        K.Clear();
        K.MultiplyTN(Tinv, theta_deriv);
      }

      //! calculate material stress resultants M,N from material strain resultants K, epsilon
      template <typename T>
      void straintostress(const CORE::LINALG::Matrix<3, 1, T>& Omega, const T& epsilon,
          const CORE::LINALG::Matrix<3, 3, T>& Cn, const CORE::LINALG::Matrix<3, 3, T>& Cm,
          CORE::LINALG::Matrix<3, 1, T>& M, T& f_par) const;

      //! Compute the material triad in case of the strong Kirchhoff (SK)beam formulation
      template <typename T>
      void ComputeTriadSK(const T& phi, const CORE::LINALG::Matrix<3, 1, T>& r_s,
          const CORE::LINALG::Matrix<3, 3, T>& triad_ref,
          CORE::LINALG::Matrix<3, 3, T>& triad) const
      {
        CORE::LINALG::Matrix<3, 3, T> triad_bar(true);

        // Compute triad_bar via SR mapping from triad_ref onto r_s
        CORE::LARGEROTATIONS::CalculateSRTriads<T>(r_s, triad_ref, triad_bar);

        // Compute triad via relative rotation of triad_bar
        CORE::LARGEROTATIONS::RotateTriad<T>(triad_bar, phi, triad);
      }

      template <typename T1, typename T2>
      void AssembleShapefunctionsL(CORE::LINALG::Matrix<1, BEAM3K_COLLOCATION_POINTS, T1>& L_i,
          CORE::LINALG::Matrix<1, 2 * 6 + BEAM3K_COLLOCATION_POINTS, T2>& L) const;

      template <typename T1, typename T2>
      void AssembleShapefunctionsNss(CORE::LINALG::Matrix<1, 4, T1>& N_i_xi,
          CORE::LINALG::Matrix<1, 4, T1>& N_i_xixi, double jacobi, double jacobi2,
          CORE::LINALG::Matrix<3, 2 * 6 + BEAM3K_COLLOCATION_POINTS, T2>& N_ss) const;

      template <typename T1, typename T2>
      void AssembleShapefunctionsNs(CORE::LINALG::Matrix<1, 4, T1>& N_i_xi, double jacobi,
          CORE::LINALG::Matrix<3, 2 * 6 + BEAM3K_COLLOCATION_POINTS, T2>& N_s) const;

      template <typename T1, typename T2>
      void AssembleShapefunctionsN(CORE::LINALG::Matrix<1, 4, T1>& N_i,
          CORE::LINALG::Matrix<3, 2 * 6 + BEAM3K_COLLOCATION_POINTS, T2>& N) const;

      //! update absolute values for primary Dof vector based on the given displacement vector
      template <unsigned int nnodecl, typename T>
      void UpdateDispTotlag(const std::vector<double>& disp,
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag) const;

      //! update positions vectors and tangents at boundary nodes and triads at all CPs
      //  based on the given element displacement vector
      template <unsigned int nnodecl, typename T>
      void UpdateNodalVariables(
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag,
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>&
              disp_totlag_centerline,
          std::vector<CORE::LINALG::Matrix<3, 3, T>>& triad_mat_cp,
          std::vector<CORE::LINALG::Matrix<4, 1>>& Qref_new) const;

      //! extract those Dofs relevant for centerline-interpolation from total state vector
      template <unsigned int nnodecl, unsigned int vpernode, typename T>
      void ExtractCenterlineDofValuesFromElementStateVector(
          const CORE::LINALG::Matrix<3 * vpernode * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>&
              dofvec,
          CORE::LINALG::Matrix<3 * vpernode * nnodecl, 1, T>& dofvec_centerline,
          bool add_reference_values = false) const;

      //! "add" reference values to displacement state vector (multiplicative in case of rotation
      //! pseudo vectors)
      template <unsigned int nnodecl, typename T>
      void AddRefValuesDisp(
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& dofvec) const;

      //! set positions at boundary nodes
      template <unsigned int nnodecl, typename T>
      void SetPositionsAtBoundaryNodes(
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag,
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>&
              disp_totlag_centerline) const;

      //! set tangents and triads and reference triads in quaternion form at boundary nodes
      template <unsigned int nnodecl, typename T>
      void SetTangentsAndTriadsAndReferenceTriadsAtBoundaryNodes(
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag,
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>&
              disp_totlag_centerline,
          std::vector<CORE::LINALG::Matrix<3, 3, T>>& triad_mat_cp,
          std::vector<CORE::LINALG::Matrix<4, 1>>& Qref_new) const;

      //! set triads and reference triads in quaternion form at all CPs except boundary nodes
      template <unsigned int nnodecl, typename T>
      void SetTriadsAndReferenceTriadsAtRemainingCollocationPoints(
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag,
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>&
              disp_totlag_centerline,
          std::vector<CORE::LINALG::Matrix<3, 3, T>>& triad_mat_cp,
          std::vector<CORE::LINALG::Matrix<4, 1>>& Qref_new) const;

      //! set differentiation variables for automatic differentiation via FAD
      template <unsigned int nnodecl>
      void SetAutomaticDifferentiationVariables(
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, FAD>& disp_totlag) const;

      //! Pre-multiply trafo matrix if rotvec_==true: \tilde{\vec{f}_int}=\mat{T}^T*\vec{f}_int
      template <unsigned int nnodecl, typename T>
      void ApplyRotVecTrafo(
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>&
              disp_totlag_centerline,
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& f_int) const;

      //! Transform stiffness matrix in order to solve for multiplicative rotation vector increments
      template <unsigned int nnodecl, typename T>
      void TransformStiffMatrixMultipl(CORE::LINALG::SerialDenseMatrix* stiffmatrix,
          const CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag)
          const;

      //! lump mass matrix
      void Lumpmass(CORE::LINALG::SerialDenseMatrix* emass);

      template <typename T>
      void Calculate_clcurvature(CORE::LINALG::Matrix<3, 1, T>& r_s,
          CORE::LINALG::Matrix<3, 1, T>& r_ss, CORE::LINALG::Matrix<3, 1, T>& kappacl) const
      {
        // spinmatrix Sr' = r'x
        CORE::LINALG::Matrix<3, 3, T> Srs(true);
        CORE::LARGEROTATIONS::computespin(Srs, r_s);

        // cross-product r'xr''
        CORE::LINALG::Matrix<3, 1, T> Srsrss(true);
        Srsrss.Multiply(Srs, r_ss);
        T rsTrs = 0.0;

        for (int i = 0; i < 3; i++) rsTrs += r_s(i) * r_s(i);

        for (int i = 0; i < 3; i++)
        {
          kappacl(i) = Srsrss(i) / rsTrs;
        }

        return;
      }

      template <typename T>
      void computestrainSK(const T& phi_s, const CORE::LINALG::Matrix<3, 1, T>& kappacl,
          const CORE::LINALG::Matrix<3, 3, T>& triadref,
          const CORE::LINALG::Matrix<3, 3, T>& triad_mat, CORE::LINALG::Matrix<3, 1, T>& K) const
      {
        CORE::LINALG::Matrix<1, 1, T> scalar_aux(true);
        CORE::LINALG::Matrix<3, 1, T> g1(true);
        CORE::LINALG::Matrix<3, 1, T> g2(true);
        CORE::LINALG::Matrix<3, 1, T> g3(true);
        CORE::LINALG::Matrix<3, 1, T> gref1(true);
        T KR1 = 0.0;

        for (int i = 0; i < 3; i++)
        {
          g1(i) = triad_mat(i, 0);
          g2(i) = triad_mat(i, 1);
          g3(i) = triad_mat(i, 2);
          gref1(i) = triadref(i, 0);
        }

        scalar_aux.MultiplyTN(kappacl, gref1);
        KR1 = -scalar_aux(0, 0);
        scalar_aux.Clear();
        scalar_aux.MultiplyTN(g1, gref1);
        KR1 = KR1 / (1.0 + scalar_aux(0, 0));
        K(0) = KR1 + phi_s;

        scalar_aux.Clear();
        scalar_aux.MultiplyTN(kappacl, g2);
        K(1) = scalar_aux(0, 0);

        scalar_aux.Clear();
        scalar_aux.MultiplyTN(kappacl, g3);
        K(2) = scalar_aux(0, 0);

        return;
      }

      void ResizeClassVariables(const int& n)
      {
        Qrefconv_.resize(BEAM3K_COLLOCATION_POINTS);
        Qrefnew_.resize(BEAM3K_COLLOCATION_POINTS);
        K0_.resize(n);
        jacobi_.resize(n);
        jacobi2_.resize(n);
        jacobi_cp_.resize(BEAM3K_COLLOCATION_POINTS);
        jacobi2_cp_.resize(BEAM3K_COLLOCATION_POINTS);
        Qconvmass_.resize(n);
        Qnewmass_.resize(n);
        wconvmass_.resize(n);
        wnewmass_.resize(n);
        aconvmass_.resize(n);
        anewmass_.resize(n);
        amodconvmass_.resize(n);
        amodnewmass_.resize(n);
        rttconvmass_.resize(n);
        rttnewmass_.resize(n);
        rttmodconvmass_.resize(n);
        rttmodnewmass_.resize(n);
        rtconvmass_.resize(n);
        rtnewmass_.resize(n);
        rconvmass_.resize(n);
        rnewmass_.resize(n);



        return;
      }

      void SetInitialDynamicClassVariables(const int& num,
          const CORE::LINALG::Matrix<3, 3>& triad_mat, const CORE::LINALG::Matrix<3, 1>& r)
      {
        Qconvmass_[num].Clear();
        Qnewmass_[num].Clear();
        rconvmass_[num].Clear();
        rnewmass_[num].Clear();
        wconvmass_[num].Clear();
        wnewmass_[num].Clear();
        aconvmass_[num].Clear();
        anewmass_[num].Clear();
        amodconvmass_[num].Clear();
        amodnewmass_[num].Clear();
        rtconvmass_[num].Clear();
        rtnewmass_[num].Clear();
        rttconvmass_[num].Clear();
        rttnewmass_[num].Clear();
        rttmodconvmass_[num].Clear();
        rttmodnewmass_[num].Clear();

        CORE::LARGEROTATIONS::triadtoquaternion(triad_mat, Qconvmass_[num]);
        Qnewmass_[num] = Qconvmass_[num];
        rconvmass_[num] = r;
        rnewmass_[num] = r;

        return;
      }

      // don't want = operator
      Beam3k& operator=(const Beam3k& old);

      template <int dim>
      void ComputeTripleProduct(CORE::LINALG::Matrix<3, dim, FAD>& mat1,
          CORE::LINALG::Matrix<3, 1, FAD>& vec1, CORE::LINALG::Matrix<3, 1, FAD>& vec2,
          CORE::LINALG::Matrix<dim, 1, FAD>& vec_out) const
      {
        CORE::LINALG::Matrix<3, 3, FAD> auxmatrix1(true);
        CORE::LINALG::Matrix<3, 1, FAD> auxvec1(true);
        CORE::LARGEROTATIONS::computespin(auxmatrix1, vec1);
        auxvec1.Multiply(auxmatrix1, vec2);
        vec_out.MultiplyTN(mat1, auxvec1);

        return;
      }

      /** \brief compute interpolated velocity vector from element state vector
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl, unsigned int vpernode, unsigned int ndim>
      void Calc_velocity(
          const CORE::LINALG::Matrix<ndim * vpernode * nnodecl, 1, double>& velocity_dofvec,
          const CORE::LINALG::Matrix<1, vpernode * nnodecl, double>& N_i,
          CORE::LINALG::Matrix<ndim, 1, double>& velocity,
          const CORE::LINALG::Matrix<ndim, 1, double>& position, int gausspoint_index) const;

      /** \brief compute interpolated velocity vector if Fad is used
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl, unsigned int vpernode, unsigned int ndim>
      void Calc_velocity(
          const CORE::LINALG::Matrix<ndim * vpernode * nnodecl, 1, double>& velocity_dofvec,
          const CORE::LINALG::Matrix<1, vpernode * nnodecl, double>& N_i,
          CORE::LINALG::Matrix<ndim, 1, FAD>& velocity,
          const CORE::LINALG::Matrix<ndim, 1, FAD>& position, int gausspoint_index);

      /** \brief compute discrete strain variations v_thetaperp
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl, typename T>
      void Calc_v_thetaperp(
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, T>& v_thetaperp,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, T>& N_s,
          const CORE::LINALG::Matrix<3, 1, T>& r_s, T abs_r_s) const;

      /** \brief compute discrete strain variations v_thetapartheta
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl, typename T>
      void Calc_v_thetapartheta(
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, T>& v_thetapartheta,
          const CORE::LINALG::Matrix<1, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, T>& L,
          const CORE::LINALG::Matrix<3, 1, T>& r_s, T abs_r_s) const;

      /** \brief compute discrete strain increments v_lin_thetaperp
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void Calc_lin_thetaperp(
          CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_thetaperp,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
          const CORE::LINALG::Matrix<3, 1, double>& r_s, double abs_r_s) const;

      /** \brief compute discrete strain increments v_lin_thetapar
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void Calc_lin_thetapar(
          CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_thetapar,
          const CORE::LINALG::Matrix<1, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& L,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
          const CORE::LINALG::Matrix<3, 1, double>& g_1,
          const CORE::LINALG::Matrix<3, 1, double>& g_1_bar, double abs_r_s) const;

      /** \brief compute linearization of scaled tangent vector
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void Calc_lin_tangent_tilde(
          CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>&
              lin_tangent_tilde,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
          const CORE::LINALG::Matrix<3, 1, double>& g_1, double abs_r_s) const;

      /** \brief compute linearization of first arc-length derivative of scaled tangent vector
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void Calc_lin_tangent_tilde_s(
          CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>&
              lin_tangent_tilde_s,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_ss,
          const CORE::LINALG::Matrix<3, 1, double>& g_1,
          const CORE::LINALG::Matrix<3, 1, double>& g_1_s,
          const CORE::LINALG::Matrix<3, 1, double>& r_s,
          const CORE::LINALG::Matrix<3, 1, double>& r_ss, double abs_r_s) const;

      /** \brief compute linearization of first base vector
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void Calc_lin_g_1(
          CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_g_1,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
          const CORE::LINALG::Matrix<3, 1, double>& g_1, double abs_r_s) const;

      /** \brief compute linearization of first arc-length derivative of first base vector
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void Calc_lin_g_1_s(
          CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_g_1_s,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_ss,
          const CORE::LINALG::Matrix<3, 1, double>& g_1,
          const CORE::LINALG::Matrix<3, 1, double>& g_1_s,
          const CORE::LINALG::Matrix<3, 1, double>& r_s,
          const CORE::LINALG::Matrix<3, 1, double>& r_ss, double abs_r_s) const;

      /** \brief compute linearization of v_epsilon
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void Calc_lin_v_epsilon(CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
                                  6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_v_epsilon,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
          const CORE::LINALG::Matrix<3, 1, double>& g_1, double abs_r_s) const;

      /** \brief compute linearization of moment resultant
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void Calc_lin_moment_resultant(
          CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>&
              lin_moment_resultant,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_theta,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>&
              lin_theta_s,
          const CORE::LINALG::Matrix<3, 3, double>& spinmatrix_of_moment,
          const CORE::LINALG::Matrix<3, 3, double>& cm) const;

      /** \brief compute linearization of inertia moment
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void Calc_lin_moment_inertia(
          CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>&
              lin_moment_inertia,
          const CORE::LINALG::Matrix<3, 3, double>& triad_mat,
          const CORE::LINALG::Matrix<3, 3, double>& triad_mat_conv,
          const CORE::LINALG::Matrix<3, 1, double>& deltatheta,
          const CORE::LINALG::Matrix<3, 1, double>& angular_velocity_material,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_theta,
          const CORE::LINALG::Matrix<3, 3, double>& spinmatrix_of_moment,
          const CORE::LINALG::Matrix<3, 3, double>& C_rho, double lin_prefactor_acc,
          double lin_prefactor_vel) const;

      /** \brief compute linearization of moment from rotational damping
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void Calc_lin_moment_viscous(
          CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>&
              lin_moment_viscous,
          const CORE::LINALG::Matrix<3, 3, double>& triad_mat,
          const CORE::LINALG::Matrix<3, 3, double>& triad_mat_conv,
          const CORE::LINALG::Matrix<3, 1, double>& deltatheta,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_theta,
          const CORE::LINALG::Matrix<3, 3, double>& spinmatrix_of_moment, double gamma_polar,
          double dt) const;

      /** \brief compute linearization of v_theta_perp multiplied with moment vector
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void Calc_lin_v_thetaperp_moment(
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
              6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_v_thetaperp_moment,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
          const CORE::LINALG::Matrix<3, 1, double>& g_1, double abs_r_s,
          const CORE::LINALG::Matrix<3, 3, double>& spinmatrix_of_moment) const;

      /** \brief compute linearization of v_theta_perp_s multiplied with moment vector
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void Calc_lin_v_thetaperp_s_moment(
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
              6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_v_thetaperp_s_moment,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_ss,
          const CORE::LINALG::Matrix<3, 1, double>& g_1,
          const CORE::LINALG::Matrix<3, 1, double>& g_1_s,
          const CORE::LINALG::Matrix<3, 1, double>& r_s,
          const CORE::LINALG::Matrix<3, 1, double>& r_ss, double abs_r_s,
          const CORE::LINALG::Matrix<3, 3, double>& spinmatrix_of_moment) const;

      /** \brief compute linearization of v_theta_par multiplied with moment vector
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void Calc_lin_v_thetapar_moment(
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
              6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_v_thetapar_moment,
          CORE::LINALG::Matrix<1, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& L,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
          const CORE::LINALG::Matrix<3, 1, double>& g_1, double abs_r_s,
          const CORE::LINALG::Matrix<3, 1, double>& moment) const;

      /** \brief compute linearization of v_theta_par_s multiplied with moment vector
       *
       *  \author grill
       *  \date 02/17 */
      template <unsigned int nnodecl>
      void Calc_lin_v_thetapar_s_moment(
          CORE::LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
              6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_v_thetapar_s_moment,
          CORE::LINALG::Matrix<1, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& L,
          CORE::LINALG::Matrix<1, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& L_s,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
          const CORE::LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_ss,
          const CORE::LINALG::Matrix<3, 1, double>& g_1,
          const CORE::LINALG::Matrix<3, 1, double>& g_1_s,
          const CORE::LINALG::Matrix<3, 1, double>& r_s,
          const CORE::LINALG::Matrix<3, 1, double>& r_ss, double abs_r_s,
          const CORE::LINALG::Matrix<3, 1, double>& moment) const;

     private:
      // Variables

      //! bool storing whether automatic differentiation shall be used for this element evaluation
      bool useFAD_;

      //! variable saving whether element has already been initialized (then isinit_ == true)
      bool isinit_;

      //! current (non-unit) tangent vectors at the two boundary nodes
      std::vector<CORE::LINALG::Matrix<3, 1>> T_;
      // Variables needed by all local triads
      //! Matrix holding pseudo rotation vectors describing the material triads in the initial
      //! configuration at each node
      std::vector<CORE::LINALG::Matrix<3, 1>> theta0_;
      //! quaternion describing the nodal reference triads (for the case BEAM3EK_ROT ==false) of the
      //! converged configuration of the last time step
      std::vector<CORE::LINALG::Matrix<4, 1>> Qrefconv_;
      //! quaternion describing the nodal reference triads (for the case BEAM3EK_ROT ==false) of the
      //! current configuration
      std::vector<CORE::LINALG::Matrix<4, 1>> Qrefnew_;
      //! Matrix with the material curvature in the initial configuration at each gp
      std::vector<CORE::LINALG::Matrix<3, 1>> K0_;

      //! Length of the element
      double length_;
      //! Jacobi determinant of the element at the Gauss points
      std::vector<double> jacobi_;
      //! Additional jacobi factor appearing in the second derivatives (required in case of strong
      //! Kirchhoff constraint)
      std::vector<double> jacobi2_;
      //! Jacobi determinant at collocation points
      std::vector<double> jacobi_cp_;
      //! Additional jacobi factor appearing in the second derivatives (required in linearization)
      std::vector<double> jacobi2_cp_;
      //! bool indicating, if the DoFs at the element boundary are described by rotation vectors or
      //! tangent vectors plus relative angle
      bool rotvec_;
      //! bool indicating, if Kirchhoff constraint is enforced in a strong (weakkirchhoff_==false)
      //! or weak manner (weakkirchhoff_==true)
      bool weakkirchhoff_;
      //! internal energy
      double Eint_;
      //! kinetic energy
      double Ekin_;
      //! temporarily storing the rot_damp_stiffness matrix for use in the PTC scaling operator
      CORE::LINALG::SerialDenseMatrix stiff_ptc_;


      //******************************Begin: Class variables required for time
      // integration**************************************************//
      //! triads at Gauss points for exact integration in quaternion at the end of the preceeding
      //! time step (required for computation of angular velocity)
      std::vector<CORE::LINALG::Matrix<4, 1>> Qconvmass_;
      //! current triads at Gauss points for exact integration in quaternion (required for
      //! computation of angular velocity)
      std::vector<CORE::LINALG::Matrix<4, 1>> Qnewmass_;
      //! spatial angular velocity vector at Gauss points for exact integration at the end of the
      //! preceeding time step (required for computation of inertia terms)
      std::vector<CORE::LINALG::Matrix<3, 1>> wconvmass_;
      //! current spatial angular velocity vector at Gauss points for exact integration (required
      //! for computation of inertia terms)
      std::vector<CORE::LINALG::Matrix<3, 1>> wnewmass_;
      //! spatial angular acceleration vector at Gauss points for exact integration at the end of
      //! the preceeding time step (required for computation of inertia terms)
      std::vector<CORE::LINALG::Matrix<3, 1>> aconvmass_;
      //! current spatial angular acceleration vector at Gauss points for exact integration
      //! (required for computation of inertia terms)
      std::vector<CORE::LINALG::Matrix<3, 1>> anewmass_;
      //! modified spatial angular acceleration vector (according to gen-alpha time integration) at
      //! Gauss points for exact integration at the end of the preceeding time step (required for
      //! computation of inertia terms)
      std::vector<CORE::LINALG::Matrix<3, 1>> amodconvmass_;
      //! current modified spatial angular acceleration vector (according to gen-alpha time
      //! integration) at Gauss points for exact integration (required for computation of inertia
      //! terms)
      std::vector<CORE::LINALG::Matrix<3, 1>> amodnewmass_;
      //! translational acceleration vector at Gauss points for exact integration at the end of the
      //! preceeding time step (required for computation of inertia terms)
      std::vector<CORE::LINALG::Matrix<3, 1>> rttconvmass_;
      //! current translational acceleration vector at Gauss points for exact integration (required
      //! for computation of inertia terms)
      std::vector<CORE::LINALG::Matrix<3, 1>> rttnewmass_;
      //! modified translational acceleration vector at Gauss points for exact integration at the
      //! end of the preceeding time step (required for computation of inertia terms)
      std::vector<CORE::LINALG::Matrix<3, 1>> rttmodconvmass_;
      //! modified current translational acceleration vector at Gauss points for exact integration
      //! (required for computation of inertia terms)
      std::vector<CORE::LINALG::Matrix<3, 1>> rttmodnewmass_;
      //! translational velocity vector at Gauss points for exact integration at the end of the
      //! preceeding time step (required for computation of inertia terms)
      std::vector<CORE::LINALG::Matrix<3, 1>> rtconvmass_;
      //! current translational velocity vector at Gauss points for exact integration (required for
      //! computation of inertia terms)
      std::vector<CORE::LINALG::Matrix<3, 1>> rtnewmass_;
      //! translational displacement vector at Gauss points for exact integration at the end of the
      //! preceeding time step (required for computation of inertia terms)
      std::vector<CORE::LINALG::Matrix<3, 1>> rconvmass_;
      //! current translational displacement vector at Gauss points for exact integration (required
      //! for computation of inertia terms)
      std::vector<CORE::LINALG::Matrix<3, 1>> rnewmass_;
      //******************************End: Class variables required for time
      // integration**************************************************//


      //! strain resultant values at GPs
      std::vector<double> axial_strain_GP_;

      std::vector<double> twist_GP_;
      std::vector<double> curvature_2_GP_;
      std::vector<double> curvature_3_GP_;

      //! stress resultant values at GPs
      std::vector<double> axial_force_GP_;

      std::vector<double> torque_GP_;
      std::vector<double> bending_moment_2_GP_;
      std::vector<double> bending_moment_3_GP_;
    };

    // << operator
    std::ostream& operator<<(std::ostream& os, const DRT::Element& ele);

  }  // namespace ELEMENTS

}  // namespace DRT

#endif
