/*----------------------------------------------------------------------*/
/*! \file
\brief 18 node solid shell with plasticity
\level 3
*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              seitz 11/14 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_PLAST_SSN_SOSH18_HPP
#define FOUR_C_SO3_PLAST_SSN_SOSH18_HPP

#include "4C_config.hpp"

#include "4C_so3_plast_ssn.hpp"
#include "4C_so3_plast_ssn_eletypes.hpp"
#include "4C_so3_sh18.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  // forward declarations
  class So_sh18Plast;
  class Discretization;

  namespace ELEMENTS
  {
    class SoSh18PlastType : public SoSh18Type
    {
     public:
      std::string Name() const override { return "SoSh18PlastType"; }

      static SoSh18PlastType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<CORE::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<CORE::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoSh18PlastType instance_;

      std::string get_element_type_string() const { return "SOLIDSH18PLAST"; }
    };  // class SoSh18PlastType

    class SoSh18Plast : public virtual So3Plast<CORE::FE::CellType::hex18>, public virtual SoSh18
    {
     public:
      //! @name Friends
      friend class SoSh18PlastType;


      //! Standard Constructor
      SoSh18Plast(int id,  //!< (i) this element's global id
          int owner        //!< elements owner
      );

      //! Copy Constructor
      //! Makes a deep copy of a Element
      SoSh18Plast(const SoSh18Plast& old);

      bool HaveEAS() const override { return (eastype_ != soh8p_easnone); };

      //! resolve "no unique final overrider"
      int NumVolume() const override { return SoSh18::NumVolume(); }
      CORE::FE::CellType Shape() const override { return CORE::FE::CellType::hex18; };
      int NumSurface() const override { return SoSh18::NumSurface(); }
      int NumLine() const override { return SoSh18::NumLine(); }
      std::vector<Teuchos::RCP<CORE::Elements::Element>> Lines() override
      {
        return SoSh18::Lines();
      }
      std::vector<Teuchos::RCP<CORE::Elements::Element>> Surfaces() override
      {
        return SoSh18::Surfaces();
      }
      int NumDofPerNode(const DRT::Node& node) const override
      {
        return SoSh18::NumDofPerNode(node);
      }
      int num_dof_per_element() const override { return SoSh18::num_dof_per_element(); }
      void VisNames(std::map<std::string, int>& names) override { return SoSh18::VisNames(names); }
      bool VisData(const std::string& name, std::vector<double>& data) override
      {
        return SoSh18::VisData(name, data);
      }
      int evaluate_neumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          CORE::Conditions::Condition& condition, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override
      {
        return SoSh18::evaluate_neumann(params, discretization, condition, lm, elevec1, elemat1);
      }

      //! Deep copy this instance of Solid3 and return pointer to the copy
      //!
      //! The Clone() method is used from the virtual base class Element in cases
      //! where the type of the derived class is unknown and a copy-ctor is needed
      CORE::Elements::Element* Clone() const override;


      //! Return unique ParObject id
      //!
      //! every class implementing ParObject needs a unique id defined at the top of
      //! this file.
      int UniqueParObjectId() const override
      {
        return SoSh18PlastType::Instance().UniqueParObjectId();
      }

      //! Pack this class so it can be communicated
      //! Pack and \ref Unpack are used to communicate this element
      void Pack(CORE::COMM::PackBuffer& data) const override;

      //! Unpack data from a char vector into this class
      //! Pack and \ref Unpack are used to communicate this element
      void Unpack(const std::vector<char>& data) override;

      //! Print this element
      void Print(std::ostream& os) const override;

      //! return elementtype
      SoSh18PlastType& ElementType() const override { return SoSh18PlastType::Instance(); }

      //! read input for this element
      bool ReadElement(const std::string& eletype, const std::string& distype,
          INPUT::LineDefinition* linedef) override;

      //! synchronize the eas variables in the two base-classes
      void SyncEAS();

      //! evaluate an element
      //! evaluate element stiffness, mass, internal forces, etc.
      //!
      //! if nullptr on input, the controlling method does not expect the element
      //!  to fill these matrices or vectors.
      //!
      //!  \return 0 if successful, negative otherwise
      int Evaluate(
          Teuchos::ParameterList&
              params,  //!< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  //!< pointer to discretization for de-assembly
          CORE::Elements::Element::LocationArray& la,  //!< location array for de-assembly
          CORE::LINALG::SerialDenseMatrix&
              elemat1_epetra,  //!< (stiffness-)matrix to be filled by element.
          CORE::LINALG::SerialDenseMatrix&
              elemat2_epetra,  //!< (mass-)matrix to be filled by element.
          CORE::LINALG::SerialDenseVector&
              elevec1_epetra,  //!< (internal force-)vector to be filled by element
          CORE::LINALG::SerialDenseVector& elevec2_epetra,  //!< vector to be filled by element
          CORE::LINALG::SerialDenseVector& elevec3_epetra   //!< vector to be filled by element
          ) override
      {
        return DRT::ELEMENTS::So3Plast<CORE::FE::CellType::hex18>::Evaluate(params, discretization,
            la, elemat1_epetra, elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra);
      }


     private:
      // don't want = operator
      SoSh18Plast& operator=(const SoSh18Plast& old) = delete;

      std::string get_element_type_string() const { return "SOLIDSH18PLAST"; }

     protected:
      //! Calculate nonlinear stiffness and mass matrix with condensed plastic matrices
      void nln_stiffmass(std::vector<double>& disp,  // current displacements
          std::vector<double>& vel,                  // current velocities
          std::vector<double>& temp,                 // current temperatures
          CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>*
              stiffmatrix,  // element stiffness matrix
          CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>*
              massmatrix,                                         // element mass matrix
          CORE::LINALG::Matrix<numdofperelement_, 1>* force,      // element internal force vector
          CORE::LINALG::Matrix<numgpt_post, numstr_>* elestress,  // stresses at GP
          CORE::LINALG::Matrix<numgpt_post, numstr_>* elestrain,  // strains at GP
          Teuchos::ParameterList& params,         // algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,  // stress output option
          const INPAR::STR::StrainType iostrain   // strain output option
          ) override;

      //! don't want sosh18 nlnstiffmass
      void nlnstiffmass(std::vector<int>& lm,  ///< location matrix
          std::vector<double>& disp,           ///< current displacements
          std::vector<double>& residual,       ///< current residual displ
          CORE::LINALG::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>*
              stiffmatrix,  ///< element stiffness matrix
          CORE::LINALG::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* massmatrix,  ///< element mass matrix
          CORE::LINALG::Matrix<NUMDOF_SOH18, 1>* force,  ///< element internal force vector
          CORE::LINALG::Matrix<NUMGPT_SOH18, MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          CORE::LINALG::Matrix<NUMGPT_SOH18, MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          Teuchos::ParameterList& params,         ///< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,  ///< stress output option
          const INPAR::STR::StrainType iostrain   ///< strain output option
          ) override
      {
        FOUR_C_THROW("don't want this");
      }
    };

  }  // namespace ELEMENTS



}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
