/*----------------------------------------------------------------------*/
/*! \file
\brief 18 node solid shell with plasticity
\level 3
*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              seitz 11/14 |
 *----------------------------------------------------------------------*/
#ifndef BACI_SO3_PLAST_SSN_SOSH18_H
#define BACI_SO3_PLAST_SSN_SOSH18_H

#include "baci_so3_plast_ssn.H"
#include "baci_so3_plast_ssn_eletypes.H"
#include "baci_so3_sh18.H"

namespace DRT
{
  // forward declarations
  class So_sh18Plast;
  class Discretization;

  namespace ELEMENTS
  {
    class So_sh18PlastType : public So_sh18Type
    {
     public:
      std::string Name() const override { return "So_sh18PlastType"; }

      static So_sh18PlastType& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_sh18PlastType instance_;

      std::string GetElementTypeString() const { return "SOLIDSH18PLAST"; }
    };  // class So_sh18PlastType

    class So_sh18Plast : public virtual So3_Plast<CORE::FE::CellType::hex18>, public virtual So_sh18
    {
     public:
      //! @name Friends
      friend class So_sh18PlastType;


      //! Standard Constructor
      So_sh18Plast(int id,  //!< (i) this element's global id
          int owner         //!< elements owner
      );

      //! Copy Constructor
      //! Makes a deep copy of a Element
      So_sh18Plast(const So_sh18Plast& old);

      bool HaveEAS() const override { return (eastype_ != soh8p_easnone); };

      //! resolve "no unique final overrider"
      int NumVolume() const override { return So_sh18::NumVolume(); }
      CORE::FE::CellType Shape() const override { return CORE::FE::CellType::hex18; };
      int NumSurface() const override { return So_sh18::NumSurface(); }
      int NumLine() const override { return So_sh18::NumLine(); }
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override { return So_sh18::Lines(); }
      std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override { return So_sh18::Surfaces(); }
      int NumDofPerNode(const DRT::Node& node) const override
      {
        return So_sh18::NumDofPerNode(node);
      }
      int NumDofPerElement() const override { return So_sh18::NumDofPerElement(); }
      void VisNames(std::map<std::string, int>& names) override { return So_sh18::VisNames(names); }
      bool VisData(const std::string& name, std::vector<double>& data) override
      {
        return So_sh18::VisData(name, data);
      }
      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override
      {
        return So_sh18::EvaluateNeumann(params, discretization, condition, lm, elevec1, elemat1);
      }

      //! Deep copy this instance of Solid3 and return pointer to the copy
      //!
      //! The Clone() method is used from the virtual base class Element in cases
      //! where the type of the derived class is unknown and a copy-ctor is needed
      DRT::Element* Clone() const override;


      //! Return unique ParObject id
      //!
      //! every class implementing ParObject needs a unique id defined at the top of
      //! this file.
      int UniqueParObjectId() const override
      {
        return So_sh18PlastType::Instance().UniqueParObjectId();
      }

      //! Pack this class so it can be communicated
      //! Pack and \ref Unpack are used to communicate this element
      void Pack(DRT::PackBuffer& data) const override;

      //! Unpack data from a char vector into this class
      //! Pack and \ref Unpack are used to communicate this element
      void Unpack(const std::vector<char>& data) override;

      //! Print this element
      void Print(std::ostream& os) const override;

      //! return elementtype
      So_sh18PlastType& ElementType() const override { return So_sh18PlastType::Instance(); }

      //! read input for this element
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef) override;

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
          DRT::Element::LocationArray& la,      //!< location array for de-assembly
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
        return DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex18>::Evaluate(params, discretization,
            la, elemat1_epetra, elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra);
      }


     private:
      // don't want = operator
      So_sh18Plast& operator=(const So_sh18Plast& old) = delete;

      std::string GetElementTypeString() const { return "SOLIDSH18PLAST"; }

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
        dserror("don't want this");
      }
    };

  }  // namespace ELEMENTS



}  // namespace DRT

#endif  // SO3_SSN_PLAST_SOSH18_H
