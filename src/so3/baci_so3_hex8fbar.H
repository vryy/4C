/*----------------------------------------------------------------------*/
/*! \file

\brief Solid Hex8 element with F-bar modification

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_SO3_HEX8FBAR_H
#define BACI_SO3_HEX8FBAR_H

#include "baci_so3_hex8.H"

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {

    class So_hex8fbarType : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "So_hex8fbarType"; }

      static So_hex8fbarType& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Teuchos::SerialDenseMatrix<int, double> ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_hex8fbarType instance_;

      std::string GetElementTypeString() const { return "SOLIDH8FBAR"; }
    };

    /*!
    \brief A C++ version of a 8-node hex solid element with F-Bar modification

    A structural 8-node hexahedral solid displacement element for large deformations
    and (near)-incompressibility. The F-bar technique is used to avoid volumetric locking.
    The volumetric part of the deformation gradient is only evaluated at the element
    center, which can be interpreted as a kind of selective reduced integration approach.
    The deviatoric part of the deformation gradient remains untouched.

    Refer also to the HiWi report of Stefanos Tsoukalas, 2010

    */
    class So_hex8fbar : public So_hex8
    {
     public:
      //! @name Friends
      friend class So_hex8fbarType;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      So_hex8fbar(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      So_hex8fbar(const So_hex8fbar& old);

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override
      {
        return So_hex8fbarType::Instance().UniqueParObjectId();
      }

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

      /*!
      \brief Destructor

      */
      ~So_hex8fbar() override;

      //@}

      //! @name Access methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return So_hex8fbarType::Instance(); }

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef) override;

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate So_hex8fbar element stiffness, mass, internal forces, etc.

      If NULL on input, the controlling method does not expect the element
      to fill these matrices or vectors.

      \return 0 if successful, negative otherwise
      */
      int Evaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  ///< pointer to discretization for de-assembly
          std::vector<int>& lm,                 ///< location matrix for de-assembly
          Epetra_SerialDenseMatrix& elemat1,    ///< (stiffness-)matrix to be filled by element.
          Epetra_SerialDenseMatrix& elemat2,    ///< (mass-)matrix to be filled by element.
          Epetra_SerialDenseVector& elevec1,    ///< (internal force-)vector to be filled by element
          Epetra_SerialDenseVector& elevec2,    ///< vector to be filled by element
          Epetra_SerialDenseVector& elevec3     ///< vector to be filled by element
          ) override;

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a surface Neumann condition on the solid3 element

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param condition (in)     : The condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If NULL on input,

      \return 0 if successful, negative otherwise
      */
      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseMatrix* elemat1 = nullptr) override;

      //@}

     protected:
      //! don't want = operator
      So_hex8fbar& operator=(const So_hex8fbar& old);

      // compute Jacobian mapping wrt to deformed configuration
      void UpdateJacobianMapping(
          const std::vector<double>& disp, DRT::ELEMENTS::PreStress& prestress) override;
      // compute defgrd in all gp for given disp
      void DefGradient(const std::vector<double>& disp, Epetra_SerialDenseMatrix& gpdefgrd,
          DRT::ELEMENTS::PreStress& prestress) override;

      //! Calculate nonlinear stiffness and mass matrix
      virtual void nlnstiffmass(std::vector<int>& lm,  ///< location matrix
          std::vector<double>& disp,                   ///< current displacements
          std::vector<double>* acc,                    // current accelerations
          std::vector<double>& residual,               ///< current residual displ
          CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>*
              stiffmatrix,                                             ///< element stiffness matrix
          CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* massmatrix,  ///< element mass matrix
          CORE::LINALG::Matrix<NUMDOF_SOH8, 1>* force,       ///< element internal force vector
          CORE::LINALG::Matrix<NUMDOF_SOH8, 1>* forceinert,  // element inertial force vector
          CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>*
              eleplstrain,                           ///< plastic strains at GP
          Teuchos::ParameterList& params,            ///< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,     ///< stress output option
          const INPAR::STR::StrainType iostrain,     ///< strain output option
          const INPAR::STR::StrainType ioplstrain);  ///< strain output option

      //! Update history variables at the end of time step (inelastic deformation) (braeu 07/16)
      void Update_element(std::vector<double>& disp,  // current displacements
          Teuchos::ParameterList& params,             // algorithmic parameters e.g. time
          const Teuchos::RCP<MAT::Material>& mat);    // material

      //! init the inverse of the jacobian and its determinant in the material configuration
      void InitJacobianMapping() override;

     private:
      std::string GetElementTypeString() const { return "SOLIDH8FBAR"; }
    };  // class So_hex8fbar

    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================

  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // SO3_HEX8FBAR_H
