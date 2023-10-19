/*----------------------------------------------------------------------------*/
/*! \file
\brief wall1 2-D displacement element

\level 1


*/
/*---------------------------------------------------------------------------*/
/* macros */
#ifndef BACI_W1_H
#define BACI_W1_H

/*----------------------------------------------------------------------*/
/* headers */
#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_discretization_fem_general_utils_nurbs_shapefunctions.H"
#include "baci_inpar_structure.H"
#include "baci_lib_element.H"
#include "baci_lib_elementtype.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_node.H"
#include "baci_lib_parobjectfactory.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_linalg_serialdensevector.H"
#include "baci_so3_base.H"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    // forward declarations
    class Wall1Line;

    class Wall1Type : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "Wall1Type"; }

      static Wall1Type& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, int const numdof, int const dimnsp) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static Wall1Type instance_;
    };

    /*======================================================================*/
    /*!
     * \brief A C++ wrapper for the wall1 element
     *
     * A planar solid element for geometrically non-linear analyses.
     *
     *
     * <b>Typical quantities and their storage pattern</b>
     * - Green-Lagrange strain vector
     *   \f$E = \left[\begin{array}{c} E_{11} \\ E_{22} \\ 2E_{12} \end{array}\right]\f$
     * - Deformation gradient vector
     *   \f$\tilde{F} = \left[\begin{array}{c} F_{11} \\ F_{22} \\ F_{12} \\ F_{21}
     * \end{array}\right]\f$
     * - Constitutive matrix
     *   \f$C = \left[\begin{array}{ccc} C_{11} & C_{12} & C_{13} \\ C_{21} & C_{22} & C_{23} \\
     * C_{31} & C_{32} & C_{33} \end{array}\right]\f$
     * - Second Piola--Kirchhoff stress vector
     *   \f$S = \left[\begin{array}{c} S_{11} \\ S_{22} \\ S_{12} \end{array}\right]\f$
     * - Second Piola--Kirchhoff stress matrix
     *   \f$\bar{S} = \left[\begin{array}{cccc} S_{11} & 0 & S_{12} & 0 \\ 0 & S_{22} & 0 & S_{12}
     * \\ S_{12} & 0 & S_{22} & 0 \\ 0 & S_{12} & 0 & S_{11} \end{array}\right]\f$
     * - First Piola-Kirchhoff stress vector
     *   \f$\tilde{P} = \left[\begin{array}{c} P_{11} \\ P_{22} \\ P_{12} \\ P_{21}
     * \end{array}\right]\f$
     *
     * <b>References</b>
     * - [1] JC Simo and F Armero, Geometrically non-linear enhanced strain
     *   mixed methods and the method of incompatible modes, International
     *   Journal for Numerical Methods in Engineering, 33:1413-1449, 1992.
     * - [2] S Glaser and F Armero, On the formulation of enhanced strain finite
     *   elements in finite deformations, Engineering Computations, 14:759-791,
     *   1997.
     * - [3] P Wriggers, Nichtlineare Finite-Element-Methoden, Springer, 2001.
     * - [4] WA Wall and B Bornemann, Nichtlineare Finite-Element-Methoden,
     *   Vorlesungsskript, Lehrstuhl fuer Numerische Mechanik, SS 2007.
     */
    class Wall1 : public So_base
    {
     public:
      /// @name Friends
      //@{
      friend class Wall1Line;
      //@}

      /// @name object-wide constants
      //@{
      static const int maxnod_ = 9;  ///< maximally permitted number of nodes
      static const int numdim_ = 2;  ///< 2-dimensional/planar element
      static const int noddof_ = 2;  ///< number of nodal DOFs
      static const int numstr_ =
          4;  ///< number of (symmetric) strain/stress components --- THIS SHOULD BE 3
      static const int numnstr_ = 4;  ///< number of non-symmetric strain/stress components
      static const int neas_ = 4;     ///< number of EAS parameters
      //@}

      /// @name Constructors and destructors and related methods

      /// Standard Constructor
      Wall1(int id,  ///< A unique global id
          int owner);

      /// Copy Constructor
      ///
      /// Makes a deep copy of a Element
      Wall1(const Wall1& old);

      /// Deep copy this instance of Wall1 and return pointer to the copy
      ///
      /// The Clone() method is used from the virtual base class Element in cases
      /// where the type of the derived class is unknown and a copy-ctor is needed
      DRT::Element* Clone() const override;

      /// Get shape type of element
      DiscretizationType Shape() const override;

      /// Set discretization type of element
      virtual void SetDisType(DRT::Element::DiscretizationType shape)
      {
        distype_ = shape;
        return;
      };


      /// Return number of lines of this element
      int NumLine() const override
      {
        if (NumNode() == 4 || NumNode() == 8 || NumNode() == 9)
          return 4;
        else
          return 3;
      }

      /// Return number of surfaces of this element
      int NumSurface() const override { return 1; }

      /// Get vector of Teuchos::RCPs to the lines of this element
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      /// Get vector of Teuchos::RCPs to the surfaces of this element
      std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

      /// Return unique ParObject id
      ///
      /// every class implementing ParObject needs a unique id defined at the
      /// top of this file.
      int UniqueParObjectId() const override { return Wall1Type::Instance().UniqueParObjectId(); }

      /// Pack this class so it can be communicated
      ///
      /// \ref Pack and \ref Unpack are used to communicate this element
      void Pack(DRT::PackBuffer& data) const override;

      /// Unpack data from a char vector into this class
      ///
      /// \ref Pack and \ref Unpack are used to communicate this element
      void Unpack(const std::vector<char>& data) override;

      /// Destructor
      ~Wall1() override;

      //@}

      /// @name Acess methods
      //@{
      /*!
      \brief Does this element use EAS?
      */
      bool HaveEAS() const override { return (eastype_ != eas_vague); };

      /// Get number of degrees of freedom of a certain node
      /// (implements pure virtual DRT::Element)
      ///
      /// The element decides how many degrees of freedom its nodes must have.
      /// As this may vary along a simulation, the element can redecide the
      /// number of degrees of freedom per node along the way for each of it's nodes
      /// separately.
      int NumDofPerNode(const DRT::Node& node) const override { return Wall1::noddof_; }

      /// Get number of degrees of freedom per element
      /// (implements pure virtual DRT::Element)
      ///
      /// The element decides how many element degrees of freedom it has.
      /// It can redecide along the way of a simulation.
      ///
      /// \note Element degrees of freedom mentioned here are dofs that are visible
      ///      at the level of the total system of equations. Purely internal
      ///      element dofs that are condensed internally should NOT be considered.
      int NumDofPerElement() const override { return 0; }

      DRT::ElementType& ElementType() const override { return Wall1Type::Instance(); }

      //@}

      /// @name Input and Creation
      //@{

      /// Read input for this element
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef) override;

      //@}

      /// @name Evaluation
      //@{

      /// Evaluate an element
      ///
      /// Evaluate Wall1 element stiffness, mass, internal forces etc
      ///
      /// \return 0 if successful, negative otherwise
      int Evaluate(Teuchos::ParameterList& params,  ///< (in/out) ParameterList for communication
                                                    ///< between control routine and elements
          DRT::Discretization& discretization,  ///< A reference to the underlying discretization
          std::vector<int>& lm,                 ///< location vector of this element
          CORE::LINALG::SerialDenseMatrix& elemat1,  ///< matrix to be filled by element. If nullptr
                                                     ///< on input, the controling method does not
                                                     ///< epxect the element to fill this matrix.
          CORE::LINALG::SerialDenseMatrix& elemat2,  ///< matrix to be filled by element. If nullptr
                                                     ///< on input, the controling method does not
                                                     ///< epxect the element to fill this matrix.
          CORE::LINALG::SerialDenseVector& elevec1,  ///< vector to be filled by element. If nullptr
                                                     ///< on input, the controlling method does not
                                                     ///< epxect the element to fill this vector
          CORE::LINALG::SerialDenseVector& elevec2,  ///< vector to be filled by element. If nullptr
                                                     ///< on input, the controlling method does not
                                                     ///< epxect the element to fill this vector
          CORE::LINALG::SerialDenseVector& elevec3   ///< vector to be filled by element. If nullptr
                                                     ///< on input, the controlling method does not
                                                     ///< epxect the element to fill this vector
          ) override;

      /// Evaluate a Neumann boundary condition
      ///
      /// this method evaluates a surfaces Neumann condition on the wall element
      ///
      /// \return 0 if successful, negative otherwise
      int EvaluateNeumann(
          Teuchos::ParameterList& params,  ///< (in/out) ParameterList for communication between
                                           ///< control routine and elements
          DRT::Discretization& discretization,  ///< A reference to the underlying discretization
          DRT::Condition& condition,            ///<  The condition to be evaluated
          std::vector<int>& lm,                 ///< location vector of this element
          CORE::LINALG::SerialDenseVector&
              elevec1,  ///< vector to be filled by element. If nullptr on input
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;

      //@}

     protected:
      /// type of 2D dimension reduction
      enum DimensionalReduction
      {
        plane_none,    ///< undetermined
        plane_stress,  ///< plane stress, i.e. lateral stress is zero \f$S_{33}=S_{13}=S_{23}=0\f$
        plane_strain   ///< plane strain, i.e. lateral strain is zero \f$E_{33}=E_{13}=E_{23}=0\f$
      };

      /// type of stresses to be calculatecd in postprocessing
      ///
      /// describes in co-ordinate system the stress components are expressed
      enum StressType
      {
        w1_none,  ///< undetermined
        w1_xy,    ///< Cartesian
        w1_rs     ///< local/parametric/natural
      };

      /// EAS type
      enum EasType
      {
        eas_vague,  ///< unknown type
        eas_q1e4,   ///< Q1E4 formulation due to Simo and Armero [1]
        eas_q1et4   ///< Q1ET4 formulation due to Glaser and Armero [2]
      };

      /// container containing all kind of wall1 thisandthat
      DRT::Container data_;

      /// number of the material law
      int material_;
      /// element thickness
      double thickness_;
      // line search parameter (old step length)
      double old_step_length_;
      /// gaussian points
      CORE::DRT::UTILS::GaussRule2D gaussrule_;
      /// problem type
      DimensionalReduction wtype_;
      /// type of stress calculation
      StressType stresstype_;
      /// eas or not
      bool iseas_;
      /// EAS type
      enum EasType eastype_;
      //! struct_ale
      bool structale_;

      //! the element discretization type (shape)
      DRT::Element::DiscretizationType distype_;

      /// @name Internal calculation methods
      //@{
      /** recover elementwise stored stuff
       *
       * \author hiermeier
       * \date 05/16 */
      void w1_recover(const std::vector<int>& lm, const std::vector<double>& disp,
          const std::vector<double>& residual);

      /// evaluate the element forces and stiffness and mass
      /// \author mgit \date 03/07
      void w1_nlnstiffmass(const std::vector<int>& lm,  ///< location vector
          const std::vector<double>& disp,              ///< element displacements
          const std::vector<double>& residual,          ///< residual displacements
          const std::vector<double>& dispmat,           ///< residual displacements
          std::vector<CORE::LINALG::SerialDenseVector>&
              myknots,                                   ///< knot vector for nurbs elements
          CORE::LINALG::SerialDenseMatrix* stiffmatrix,  ///< element stiffness matrix
          CORE::LINALG::SerialDenseMatrix* massmatrix,   ///< element mass matrix
          CORE::LINALG::SerialDenseVector* force,        ///< element internal force vector
          CORE::LINALG::SerialDenseMatrix* elestress,    ///< element stresses
          CORE::LINALG::SerialDenseMatrix* elestrain,    ///< element strains
          Teuchos::RCP<const MAT::Material> material,    ///< element material
          Teuchos::ParameterList& params,                ///< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,         ///< stress output option
          const INPAR::STR::StrainType iostrain          ///< strain output option
      );

      /// evaluate the geometrically linear element forces and stiffness and mass
      void w1_linstiffmass(const std::vector<int>& lm,  ///< location vector
          const std::vector<double>& disp,              ///< element displacements
          const std::vector<double>& residual,          ///< residual displacements
          const std::vector<double>& dispmat,           ///< residual displacements
          std::vector<CORE::LINALG::SerialDenseVector>&
              myknots,                                   ///< knot vector for nurbs elements
          CORE::LINALG::SerialDenseMatrix* stiffmatrix,  ///< element stiffness matrix
          CORE::LINALG::SerialDenseMatrix* massmatrix,   ///< element mass matrix
          CORE::LINALG::SerialDenseVector* force,        ///< element internal force vector
          CORE::LINALG::SerialDenseMatrix* elestress,    ///< element stresses
          CORE::LINALG::SerialDenseMatrix* elestrain,    ///< element strains
          Teuchos::RCP<const MAT::Material> material,    ///< element material
          Teuchos::ParameterList& params,                ///< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,         ///< stress output option
          const INPAR::STR::StrainType iostrain          ///< strain output option
      );

      /// Jacobian matrix for mapping from parameter space in physical material space
      /// at point parameter space
      /// \author mgit \date 04/07
      void w1_jacobianmatrix(const CORE::LINALG::SerialDenseMatrix&
                                 xrefe,  ///< reference/material co-ordinates of element nodes
          const CORE::LINALG::SerialDenseMatrix&
              deriv,  ///< derivatives of shape functions at parameter point
          CORE::LINALG::SerialDenseMatrix& xjm,  ///< Jacobi matrix
          double* det,                           ///< determinant of Jacobi matrix
          const int iel                          ///< actual number of element nodes
      );

      /// Linear B-operator in reference configuration at point parameter space
      /// \author mgit \date 04/07
      void w1_boplin(CORE::LINALG::SerialDenseMatrix& boplin,  ///< the B-operator
          CORE::LINALG::SerialDenseMatrix&
              deriv,  ///< derivatives of shape functions at parameter point
          CORE::LINALG::SerialDenseMatrix& xjm,  ///< Jacobian at parameter point
          double& det,                           ///< Jacobian determinant at parameter point
          const int iel                          ///< number of element nodes
      );

      /// (Material) Deformation gradient \f$F\f$ and Green-Lagrange strains \f$E\f$
      /// at parameter point
      /// \author mgit \date 04/07
      void w1_defgrad(CORE::LINALG::SerialDenseVector& F,  ///< deformation gradient
          CORE::LINALG::SerialDenseVector&
              strain,  ///< GL strain \f$E^T=[E_{11} \; E_{22} \; E_{12}]\f$
          const CORE::LINALG::SerialDenseMatrix&
              xrefe,  ///< reference/material co-ordinates of element nodes
          const CORE::LINALG::SerialDenseMatrix&
              xcure,  ///< current/spatial co-ordinates of element nodes
          CORE::LINALG::SerialDenseMatrix& boplin,  ///< linear B-operator
          const int iel                             ///< number of element nodes
      );

      /// Deformation gradient, measures distortion of the mesh of the
      /// material configuration with respect to the referential configuration
      /// in structure with ale approaches (fractional step method)
      /// \author mgit \date 04/11
      void w1_defgradmat(
          CORE::LINALG::SerialDenseVector&
              F,  ///< deformation gradient (mesh distortion of material configuration)
          CORE::LINALG::SerialDenseVector&
              Fmat,  ///< deformation gradient (mesh distortion of spatial configuration)
          CORE::LINALG::SerialDenseVector& FFmatinv,  ///< product of F and Fmat(inv)
          CORE::LINALG::SerialDenseVector&
              strain,  ///< GL strain \f$E^T=[E_{11} \; E_{22} \; E_{12}]\f$
          const CORE::LINALG::SerialDenseMatrix& xrefe,  ///< reference coordinates of element nodes
          const CORE::LINALG::SerialDenseMatrix& xmat,   ///< material coordinates of element nodes
          CORE::LINALG::SerialDenseMatrix& boplin,       ///< linear B-operator
          const int iel                                  ///< number of element nodes
      );

      /// Non-linear B-operator in reference configuration
      /// \author mgit \date 04/07
      void w1_boplin_cure(CORE::LINALG::SerialDenseMatrix& b_cure,  ///< non-linear B-operator
          const CORE::LINALG::SerialDenseMatrix& boplin,            ///< linear B-operator
          const CORE::LINALG::SerialDenseVector& F,  ///< deformation gradient as Voigt-vector
          const int numeps,                          ///< number of strains
          const int nd                               ///< number of element nodes
      );

      /// Geometric stiffness constribution (total Lagrange)
      /// \author mgit \date 05/07
      void w1_kg(CORE::LINALG::SerialDenseMatrix& estif,  ///< (in/out) element stiffness matrix
          const CORE::LINALG::SerialDenseMatrix& boplin,  ///< linear B-operator
          const CORE::LINALG::SerialDenseMatrix& stress,  ///< PK2 stress vector
          const double fac,                               ///< integration factor
          const int nd,                                   ///< number of element DOFs
          const int numeps                                ///< number of strains
      );

      /// elastic and initial displacement stiffness contribution (total Lagrange)
      /// \author mgit \date 05/07
      void w1_keu(CORE::LINALG::SerialDenseMatrix& estif,  ///< (in/out) element stiffness matrix
          const CORE::LINALG::SerialDenseMatrix& b_cure,   ///< non-linear B-operator
          const CORE::LINALG::SerialDenseMatrix& C,        ///< elasticity matrix
          const double fac,                                ///< integration factor
          const int nd,                                    ///< number of element DOFs
          const int numeps                                 ///< number of strains
      );

      /// Evaluate internal element forces for large def (total Lagr)
      /// \author mgit \date 05/07
      void w1_fint(const CORE::LINALG::SerialDenseMatrix& stress,  ///< PK2 stress vector
          const CORE::LINALG::SerialDenseMatrix& b_cure,           ///< non-linear B-op
          CORE::LINALG::SerialDenseVector& intforce,               ///< internal force vector
          const double fac,                                        ///< integration factor
          const int nd                                             ///< number of element DOFs
      );

      /// lump mass matrix
      void w1_lumpmass(CORE::LINALG::SerialDenseMatrix* emass  ///< (in/out) element mass matrix
      );

      /// determine cauchy stress and store it
      void StressCauchy(const int ip,                     ///< Gauss point index
          const double& F11,                              ///< F_{11} component of def.grad.
          const double& F22,                              ///< F_{22} component of def.grad.
          const double& F12,                              ///< F_{12} component of def.grad.
          const double& F21,                              ///< F_{21} component of def.grad.
          const CORE::LINALG::SerialDenseMatrix& stress,  ///< PK2 stress matrix
          CORE::LINALG::SerialDenseMatrix* elestress      ///< element stress array
      );

      /// determine internal energy
      void Energy(Teuchos::ParameterList& params,  ///<  a list containing GEMM coefficients
          const std::vector<int>& lm,              ///< location vector
          const std::vector<double>&
              dis,  ///< element displacements \f$d_{n}^{(e)}\f$ at \f$t_{n}\f$
          CORE::LINALG::SerialDenseVector* energies,  ///< (in/out) energies
          Teuchos::RCP<const MAT::Material> material  ///< element material
      );

      //@}

      /// @name EAS tools
      //@{
      /// return the number of eas dofs
      int w1_neas() const { return neas_; };

      /// setup of constant EAS data
      void w1_eassetup(CORE::LINALG::SerialDenseMatrix& boplin0,  ///< linear B-op at origin
          CORE::LINALG::SerialDenseVector& F0,             ///< deformation gradient at origin
          CORE::LINALG::SerialDenseMatrix& xjm0,           ///< jacobian matrix at origin
          double& detJ0,                                   ///< det of Jacobian at origin
          const CORE::LINALG::SerialDenseMatrix& xrefe,    ///< material element coords
          const CORE::LINALG::SerialDenseMatrix& xcure,    ///< current element coords
          const DRT::Element::DiscretizationType& distype  ///< discretisation type
      );

      /// Get the enhanced deformation gradient and
      /// also the operators G,W0 and Z
      /// at point in parameter space
      /// \author mgit \date 01/08
      void w1_call_defgrad_enh(
          CORE::LINALG::SerialDenseMatrix& F_enh,       ///< enhanced deformation gradient vector
          const CORE::LINALG::SerialDenseMatrix xjm0,   ///< Jacobian at origin
          const CORE::LINALG::SerialDenseMatrix xjm,    ///< Jacobian at parameter point
          const double detJ0,                           ///< det of Jacobian at origin
          const double det,                             ///< det of Jacobian at parameter point
          const CORE::LINALG::SerialDenseVector F0,     ///< deformation gradient at origin
          const CORE::LINALG::SerialDenseMatrix alpha,  ///< alpha parameters (EAS params)
          const double e1,                      ///< \f$\xi\f$ co-ordinate of parameter point
          const double e2,                      ///< \f$\eta\f$ co-ordinate of parameter point
          CORE::LINALG::SerialDenseMatrix& G,   ///< G-operator
          CORE::LINALG::SerialDenseMatrix& W0,  ///< W_0-operator
          const CORE::LINALG::SerialDenseMatrix boplin0,  ///< linear B-op at origin
          CORE::LINALG::SerialDenseMatrix& Z              ///< Z-operator
      );

      /// Total deformation gradient and (enhanced) Green-Lagrange strain
      /// \author mgit \date 01/08
      void w1_call_defgrad_tot(
          const CORE::LINALG::SerialDenseMatrix& F_enh,  ///< enhanced def.grad.
          CORE::LINALG::SerialDenseMatrix& F_tot,        ///< total def.grad.
          const CORE::LINALG::SerialDenseVector& F,      ///< displ.-based def.grad.
          CORE::LINALG::SerialDenseVector& strain        ///< GL strains
      );

      /// first Piola-Kirchhoff stress vector
      /// \author mgit \author 02/08
      void w1_stress_eas(const CORE::LINALG::SerialDenseMatrix& stress,  ///< PK2 stress vector
          const CORE::LINALG::SerialDenseMatrix& F_tot,                  ///< total def.grad.
          CORE::LINALG::SerialDenseMatrix& p_stress                      ///< PK1 stress vector
      );

      /// calculate stiffness matrix kdd \f$\partial f_{int}/\partial d\f$
      /// \author mgit \date 03/08
      void w1_kdd(const CORE::LINALG::SerialDenseMatrix& boplin,  ///< linear B-op
          const CORE::LINALG::SerialDenseMatrix& W0,              ///< W-operator at origin
          const CORE::LINALG::SerialDenseMatrix& F,               ///< total def.grad.
          const CORE::LINALG::SerialDenseMatrix& C,               ///< consititutive matrix
          const CORE::LINALG::SerialDenseMatrix& stress,          ///< PK2 stress vector
          CORE::LINALG::SerialDenseMatrix& FCF,                   ///< \f$F^T \dot C \dot F\f$
          CORE::LINALG::SerialDenseMatrix& estif,                 ///< element stiff matrix
          const double fac                                        ///< integration factor
      );

      /// calculate tangential matrix kda \f$\partial f_{int}/\partial \alpha\f$
      /// \author mgit \date 03/08
      void w1_kda(const CORE::LINALG::SerialDenseMatrix& FCF,  ///< a product
          const CORE::LINALG::SerialDenseMatrix& W0,           ///< W-operator at origin
          const CORE::LINALG::SerialDenseMatrix& boplin,       ///< linear B-op
          const CORE::LINALG::SerialDenseMatrix& stress,       ///< PK2 stress
          const CORE::LINALG::SerialDenseMatrix& G,            ///< G-operator
          const CORE::LINALG::SerialDenseMatrix& Z,            ///< Z-operator
          CORE::LINALG::SerialDenseMatrix& Kda,                ///< target: kda
          const CORE::LINALG::SerialDenseMatrix& p_stress,     ///< PK1 stress
          const double fac                                     ///< integration factor
      );

      /// calculate tangential matrix kaa \f$\partial s/\partial \alpha\f$
      /// \author mgit \date 03/08
      void w1_kaa(
          const CORE::LINALG::SerialDenseMatrix& FCF,     ///< a product \f$F^T \dot C \dot F\f$
          const CORE::LINALG::SerialDenseMatrix& stress,  ///< PK2 stress
          const CORE::LINALG::SerialDenseMatrix& G,       ///< G-op
          CORE::LINALG::SerialDenseMatrix& Kaa,           ///< target: kaa
          const double fac                                ///< integration factor
      );

      /// calculate internal forces fint(displacements u) and feas
      void w1_fint_eas(const CORE::LINALG::SerialDenseMatrix& W0,  ///< W-op at origin
          const CORE::LINALG::SerialDenseMatrix& boplin,           ///< linear B-op
          const CORE::LINALG::SerialDenseMatrix& G,                ///< G-op
          const CORE::LINALG::SerialDenseMatrix& p_stress,         ///< PK1 stress
          CORE::LINALG::SerialDenseVector& intforce,               ///< internal force vector
          CORE::LINALG::SerialDenseVector& feas,  ///< internal EAS constraint condition
          const double fac                        ///< integration factor
      );

      //@}

      /// @name Material matters
      //@{

      /// Constitutive matrix \f$C\f$ and stresses
      /// \author mgit \date 05/07
      void w1_call_matgeononl(
          const CORE::LINALG::SerialDenseVector& strain,  ///< Green-Lagrange strain vector
          CORE::LINALG::SerialDenseMatrix& stress,        ///< stress matrix
          CORE::LINALG::SerialDenseMatrix& C,             ///< elasticity matrix
          const int numeps,                               ///< number of strains
          Teuchos::RCP<const MAT::Material> material,     ///< the material data
          Teuchos::ParameterList& params,                 ///< element parameter list
          int gp                                          ///< Gauss point
      );

      /// Stress and constitutive matrix mapper from 3d to 2d
      /// due to dimensional reduction #wtype_
      ///
      /// \author bborn \date 06/09
      void MaterialResponse3dPlane(CORE::LINALG::SerialDenseMatrix& stress,  ///< stress matrix
          CORE::LINALG::SerialDenseMatrix& C,                                ///< elasticity matrix
          const CORE::LINALG::SerialDenseVector& strain,  ///< Green-Lagrange strain vector
          Teuchos::ParameterList& params,                 ///< element parameter list
          int gp                                          ///< Gauss point
      );

      /// Generic 3D stress response
      ///
      /// \author bborn \date 06/09
      void MaterialResponse3d(
          CORE::LINALG::Matrix<6, 1>* stress,          ///< 3D 2nd Piola-Kirchhoff stress vector
          CORE::LINALG::Matrix<6, 6>* cmat,            ///< 3D elasticity matrix
          const CORE::LINALG::Matrix<6, 1>* glstrain,  ///< 3D Green-Lagrange strain vector
          Teuchos::ParameterList& params,              ///< element parameter list
          int gp                                       ///< Gauss point
      );

      /// Map plane Green-Lagrange strains to 3d
      void GreenLagrangePlane3d(const CORE::LINALG::SerialDenseVector&
                                    glplane,  ///< 2d version of GL strain (Voigt notation)
          CORE::LINALG::Matrix<6, 1>& gl3d    ///< 3d version of GL strain (Voigt notation)
      );

      /// Internal/strain energy
      double EnergyInternal(Teuchos::RCP<const MAT::Material> material,  ///< element material
          Teuchos::ParameterList& params,                                ///< element parameter list
          const CORE::LINALG::SerialDenseVector& Ev,  ///< Green-Lagrange strain vector
          int gp                                      ///< Gauss point
      );

      /// Kinetic Energy
      /// \return kinetic energy of element
      double EnergyKinetic(const CORE::LINALG::SerialDenseMatrix& mass,  ///< element mass matrix
          const std::vector<double>& vel  ///< element velocity vector
      );

      //@}



      /// @name Tools for generalised energy-momentum method
      //@{

      /// evaluate the element forces and stiffness and mass for GEMM
      /// \author bborn \date 08/08
      void FintStiffMassGEMM(
          Teuchos::ParameterList& params,  ///<  a list containing GEMM coefficients
          const std::vector<int>& lm,      ///< location vector
          const std::vector<double>&
              dispo,  ///< element displacements \f$d_{n+1}^{(e)}\f$ at \f$t_{n}\f$
          const std::vector<double>&
              disp,  ///< element displacements \f$d_{n}^{(e)}\f$ at \f$t_{n+1}\f$
          const std::vector<double>&
              residual,  ///< residual displacements \f$\Delta d_{n+1}^{<k>(e)}\f$
          CORE::LINALG::SerialDenseMatrix* stiffmatrix,  ///< element stiffness matrix \f$k^{(e)}\f$
          CORE::LINALG::SerialDenseMatrix* massmatrix,   ///< element mass matrix \f$m^{(e)}\f$
          CORE::LINALG::SerialDenseVector*
              force,  ///< element internal mid-force vector \f$f_{int;m}^{(e)}\f$
          CORE::LINALG::SerialDenseMatrix* elestress,  ///< element stresses
          CORE::LINALG::SerialDenseMatrix* elestrain,  ///< element strains
          Teuchos::RCP<const MAT::Material> material,  ///< element material
          const INPAR::STR::StressType iostress,       ///< stress output option
          const INPAR::STR::StrainType iostrain        ///< strain output option
      );

      /// Tangent of internal force vector \f$f_{int;m}\f$ with respect to displacements
      /// \f$d_{n+1}\f$ in pure displacement-based approach
      void TangFintByDispGEMM(const double& alphafgemm,  ///< GEMM's \f$\alpha_f\f$ coefficient
          const double& xigemm,                          ///< GEMM's \f$\xi\f$ coefficient
          const double& fac,                             ///< integration factor
          const CORE::LINALG::SerialDenseMatrix& bopm,   ///< non-lin. mid-B-op. \f$B_{m}\f$
          const CORE::LINALG::SerialDenseMatrix&
              bop,  ///< non-linear B-operator \f$B_{n+1}\f$ at \f$t_{n+1}\f$
          const CORE::LINALG::SerialDenseMatrix& C,       ///< constitutive matrix
          const CORE::LINALG::SerialDenseMatrix& boplin,  ///< linear B-operator
          const CORE::LINALG::SerialDenseMatrix& Smm,  ///< PK2 mid-stress matrix \f$\bar{S}_{m}\f$
          CORE::LINALG::SerialDenseMatrix& estif       ///< (in/out) element stiffness matrix
      );

      /// Tangent of internal force vector \f$f_{int;m}\f$ with respect to displacements
      /// \f$d_{n+1}\f$ in EAS approach
      void TangFintByDispGEMM(const double& alphafgemm,   ///< GEMM's \f$\alpha_f\f$ coefficient
          const double& xigemm,                           ///< GEMM's \f$\xi\f$ coefficient
          const double& fac,                              ///< integration factor
          const CORE::LINALG::SerialDenseMatrix& boplin,  ///< linear B-operator \f$B_L\f$
          const CORE::LINALG::SerialDenseMatrix& W0m,     ///< mid-W-operator \f$W_{0;m}\f$
          const CORE::LINALG::SerialDenseMatrix& W0,      ///< W-operator \f$W_{0;n+1}\f$
          const CORE::LINALG::SerialDenseMatrix& Fmm,     ///< mid def.grad. matrix \f$\bar{F}_m\f$
          const CORE::LINALG::SerialDenseMatrix& Fm,      ///< def.grad. matrix \f$\bar{F}_{n+1}\f$
          const CORE::LINALG::SerialDenseMatrix& C,       ///< constititive matrix \f$C\f$
          const CORE::LINALG::SerialDenseMatrix&
              Smm,  ///< 2nd Piola-Kirchhof mid-stress matrix \f$\bar{S}_m\f$
          CORE::LINALG::SerialDenseMatrix&
              FmCF,  ///< the product \f$\bar{F}_m \cdot C \cdot \bar{F}_{n+1}^T\f$
          CORE::LINALG::SerialDenseMatrix& estif  ///< (in/out) element stiffness matrix
      );

      /// Tangent of internal force vector \f$f_{int,m}\f$ with respect to enhanced strain
      /// parameters \f$\alpha_{n+1}\f$
      void TangFintByEnhGEMM(const double& alphafgemm,    ///< GEMM's \f$\alpha_f\f$ coefficient
          const double& xigemm,                           ///< GEMM's \f$\xi\f$ coefficient
          const double& fac,                              ///< integration factor
          const CORE::LINALG::SerialDenseMatrix& boplin,  ///< linear B-operator \f$B_L\f$
          const CORE::LINALG::SerialDenseMatrix& W0m,     ///< mid-W-operator \f$W_{0;m}\f$
          const CORE::LINALG::SerialDenseMatrix&
              FmCF,  ///< the product \f$\bar{F}_m \cdot C \cdot \bar{F}_{n+1}^T\f$
          const CORE::LINALG::SerialDenseMatrix&
              Smm,  ///< 2nd Piola-Kirchhof mid-stress matrix \f$\bar{S}_m\f$
          const CORE::LINALG::SerialDenseMatrix& G,  ///< G-operator \f$G_{n+1}\f$
          const CORE::LINALG::SerialDenseMatrix& Z,  ///< Z-operator \f$Z_{n+1}\f$
          const CORE::LINALG::SerialDenseMatrix&
              Pvmm,  ///< mid-mid 1st Piola-Kirchhoff stress vector
                     ///< \f$\tilde{P}_{mm} = \bar{F}_m \cdot S_m\f$
          CORE::LINALG::SerialDenseMatrix&
              kda  ///< (in/out) tangent \f$\partial f_{int;m}/\partial \alpha_{n+1}\f$
      );

      /// Tangent of enhancemed constraint vector \f$s_m\f$ with respect to displacements
      /// \f$d_{n+1}\f$
      void TangEconByDispGEMM(const double& alphafgemm,   ///< GEMM's \f$\alpha_f\f$ coefficient
          const double& xigemm,                           ///< GEMM's \f$\xi\f$ coefficient
          const double& fac,                              ///< integration factor
          const CORE::LINALG::SerialDenseMatrix& boplin,  ///< linear B-operator \f$B_L\f$
          const CORE::LINALG::SerialDenseMatrix& W0,      ///< W-operator \f$W_{0;n+1}\f$
          const CORE::LINALG::SerialDenseMatrix&
              FmCF,  ///< the product \f$\bar{F}_m \cdot C \cdot \bar{F}_{n+1}^T\f$
          const CORE::LINALG::SerialDenseMatrix&
              Smm,  ///< 2nd Piola-Kirchhof mid-stress matrix \f$\bar{S}_m\f$
          const CORE::LINALG::SerialDenseMatrix& Gm,  ///< mid-G-operator \f$G_{m}\f$
          const CORE::LINALG::SerialDenseMatrix& Z,   ///< Z-operator \f$Z_{n+1}\f$
          const CORE::LINALG::SerialDenseMatrix&
              Pvmm,  ///< mid-mid 1st Piola-Kirchhoff stress vector
                     ///< \f$\tilde{P}_{mm} = \bar{F}_m \cdot S_m\f$
          CORE::LINALG::SerialDenseMatrix&
              kad  ///< (in/out) tangent \f$\partial s_{m}/\partial d_{n+1}\f$
      );

      /// Tangent of enhancemed constraint vector \f$s_m\f$ with respect to enhanced strain
      /// parameters \f$\alpha_{n+1}\f$
      void TangEconByEnhGEMM(const double& alphafgemm,  ///< GEMM's \f$\alpha_f\f$ coefficient
          const double& xigemm,                         ///< GEMM's \f$\xi\f$ coefficient
          const double& fac,                            ///< integration factor
          const CORE::LINALG::SerialDenseMatrix&
              FmCF,  ///< the product \f$\bar{F}_m \cdot C \cdot \bar{F}_{n+1}^T\f$
          const CORE::LINALG::SerialDenseMatrix&
              Smm,  ///< 2nd Piola-Kirchhof mid-stress matrix \f$\bar{S}_m\f$
          const CORE::LINALG::SerialDenseMatrix& G,   ///< G-operator \f$G_{n+1}\f$
          const CORE::LINALG::SerialDenseMatrix& Gm,  ///< mid-G-operator \f$G_{m}\f$
          CORE::LINALG::SerialDenseMatrix&
              kaa  ///< (in/out) tangent \f$\partial s_{m}/\partial \alpha_{n+1}\f$
      );

      //@}

      /// don't want = operator
      Wall1& operator=(const Wall1& old);


      /// set number of gauss points to element shape default
      CORE::DRT::UTILS::GaussRule2D getGaussrule(int* ngp  ///< number of Gauss points
      );
    };  // class Wall1



    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================

    class Wall1LineType : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "Wall1LineType"; }

      static Wall1LineType& Instance();

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override
      {
      }

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& actnode, const double* x0, const int numdof, const int dimnsp) override
      {
        CORE::LINALG::SerialDenseMatrix nullspace;
        dserror("method ComputeNullSpace not implemented!");
        return nullspace;
      }

     private:
      static Wall1LineType instance_;
    };


    /*======================================================================*/
    /*!
     * \brief An element representing a line edge of a Wall1 element
     *
     * \note This is a pure Neumann boundary condition element. It's only
     *      purpose is to evaluate line Neumann boundary conditions that might be
     *      adjacent to a parent wall1 element. It therefore does not implement
     *      the DRT::Element::Evaluate method and does not have its own ElementRegister class.
     *
     */
    class Wall1Line : public DRT::FaceElement
    {
     public:
      /// @name Constructors and destructors and related methods
      //@{

      /// Standard Constructor
      Wall1Line(int id,        ///< A unique global id
          int owner,           ///< Processor owning this line
          int nnode,           ///< Number of nodes attached to this element
          const int* nodeids,  ///< global ids of nodes attached to this element
          DRT::Node** nodes,   ///<  the discretizations map of nodes to build ptrs to nodes from
          DRT::ELEMENTS::Wall1* parent,  ///< The parent wall element of this line
          const int lline  ///< the local line number of this line w.r.t. the parent element
      );

      /// Copy Constructor
      ///
      /// Makes a deep copy of a Element
      Wall1Line(const Wall1Line& old);


      /// Deep copy this instance of an element and return pointer to the copy
      ///
      /// The Clone() method is used from the virtual base class Element in cases
      /// where the type of the derived class is unknown and a copy-ctor is needed
      DRT::Element* Clone() const override;


      /// Get shape type of element
      DiscretizationType Shape() const override;


      /// Return unique ParObject id
      ///
      /// every class implementing ParObject needs a unique id defined at the
      /// top of the parobject.H file.
      int UniqueParObjectId() const override
      {
        return Wall1LineType::Instance().UniqueParObjectId();
      }

      /// Pack this class so it can be communicated
      ///
      /// \ref Pack and \ref Unpack are used to communicate this element
      void Pack(DRT::PackBuffer& data) const override;

      /// Unpack data from a char vector into this class
      ///
      /// \ref Pack and \ref Unpack are used to communicate this element
      void Unpack(const std::vector<char>& data) override;

      /// Destructor
      ~Wall1Line() override;

      //@}

      /// @name Access methods
      //@{

      /// Get number of degrees of freedom of a certain node
      /// (implements pure virtual DRT::Element)
      ///
      /// The element decides how many degrees of freedom its nodes must have.
      /// As this may vary along a simulation, the element can redecide the
      /// number of degrees of freedom per node along the way for each of it's nodes
      /// separately.
      int NumDofPerNode(const DRT::Node& node) const override
      {
        return ParentMasterElement()->NumDofPerNode(node);
      }

      /// Get number of degrees of freedom per element
      /// (implements pure virtual DRT::Element)
      ///
      /// The element decides how many element degrees of freedom it has.
      /// It can redecide along the way of a simulation.
      ///
      /// \note Element degrees of freedom mentioned here are dofs that are visible
      ///       at the level of the total system of equations. Purely internal
      ///       element dofs that are condensed internally should NOT be considered.
      int NumDofPerElement() const override { return 0; }

      /// Print this element
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return Wall1LineType::Instance(); }

      //@}

      /// @name Evaluate methods
      //@{

      /// Evaluate a Neumann boundary condition
      ///
      /// this method evaluates a line Neumann condition on the wall element
      ///
      /// \return 0 if successful, negative otherwise
      int EvaluateNeumann(
          Teuchos::ParameterList&
              params,  ///< (in/out) ParameterList for communication between control routine
                       ///<   and elements
          DRT::Discretization&
              discretization,         ///< (in) A reference to the underlying discretization
          DRT::Condition& condition,  ///< (in) The condition to be evaluated
          std::vector<int>& lm,       ///< (in) location vector of this element
          CORE::LINALG::SerialDenseVector&
              elevec1,  ///< (out) vector to be filled by element. If nullptr on input
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;

      /// Evaluate a boundary condition
      ///
      /// this method evaluates a line Neumann condition on the wall element
      ///
      /// \return 0 if successful, negative otherwise
      int Evaluate(Teuchos::ParameterList& params,  ///< (in/out) ParameterList for communication
                                                    ///< between control routine and elements
          DRT::Discretization&
              discretization,    ///< (in) A reference to the underlying discretization
          std::vector<int>& lm,  ///< (in) location vector of this element
          CORE::LINALG::SerialDenseMatrix& elematrix1, CORE::LINALG::SerialDenseMatrix& elematrix2,
          CORE::LINALG::SerialDenseVector& elevector1, CORE::LINALG::SerialDenseVector& elevector2,
          CORE::LINALG::SerialDenseVector& elevector3) override;

      //! Evaluate method on mutliple dofsets for wall element
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          LocationArray& la, CORE::LINALG::SerialDenseMatrix& elematrix1,
          CORE::LINALG::SerialDenseMatrix& elematrix2, CORE::LINALG::SerialDenseVector& elevector1,
          CORE::LINALG::SerialDenseVector& elevector2,
          CORE::LINALG::SerialDenseVector& elevector3) override;

      //@}

     private:
      /// digestible actions
      enum ActionType
      {
        none,  ///< undetermined
        calc_struct_areaconstrstiff,
        calc_struct_centerdisp,
        calc_struct_constrarea,
        calc_struct_area_poro
      };

      /// don't want = operator
      Wall1Line& operator=(const Wall1Line& old);

      /// Submethod to compute necessary change to stiffness matrix due to the constraint
      void ComputeAreaConstrStiff(
          CORE::LINALG::SerialDenseMatrix xscurr, CORE::LINALG::SerialDenseMatrix& elematrix);

      /// Submethod to compute first derivatives of constraint area w.r.t. the displacements
      void ComputeAreaConstrDeriv(
          CORE::LINALG::SerialDenseMatrix xscurr, CORE::LINALG::SerialDenseVector& elevector);

      /// compute infintesimal line element dr for integration along the line
      double w1_substitution(const CORE::LINALG::SerialDenseMatrix& xyze,
          const CORE::LINALG::SerialDenseMatrix& deriv, std::vector<double>* unrm, const int iel);

      /// set number of gauss points to element shape default
      CORE::DRT::UTILS::GaussRule1D getOptimalGaussrule(const DiscretizationType& distype);

    };  // class Wall1Line



  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // W1_H
