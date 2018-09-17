/*----------------------------------------------------------------------*/
/*!
\file so3_ssn_plast_sosh18.H
\brief 18 node solid shell with plasticity
\level 3
\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              seitz 11/14 |
 *----------------------------------------------------------------------*/
#ifndef SO3_SSN_PLAST_SOSH18_H
#define SO3_SSN_PLAST_SOSH18_H

#include "so3_ssn_plast.H"
#include "so3_ssn_plast_eletypes.H"
#include "../so_sh18.H"

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
      std::string Name() const { return "So_sh18PlastType"; }

      static So_sh18PlastType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_sh18PlastType instance_;
    };  // class So_sh18PlastType

    class So_sh18Plast : public virtual So3_Plast<DRT::Element::hex18>, public virtual So_sh18
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

      bool HaveEAS() const { return (eastype_ != soh8p_easnone); };

      //! resolve "no unique final overrider"
      virtual int NumVolume() const { return So_sh18::NumVolume(); }
      virtual DiscretizationType Shape() const { return hex18; };
      virtual int NumSurface() const { return So_sh18::NumSurface(); }
      virtual int NumLine() const { return So_sh18::NumLine(); }
      virtual std::vector<Teuchos::RCP<DRT::Element>> Lines() { return So_sh18::Lines(); }
      virtual std::vector<Teuchos::RCP<DRT::Element>> Surfaces() { return So_sh18::Surfaces(); }
      virtual std::vector<Teuchos::RCP<DRT::Element>> Volumes() { return So_sh18::Volumes(); }
      virtual int NumDofPerNode(const DRT::Node& node) const
      {
        return So_sh18::NumDofPerNode(node);
      }
      virtual int NumDofPerElement() const { return So_sh18::NumDofPerElement(); }
      virtual void VisNames(std::map<std::string, int>& names) { return So_sh18::VisNames(names); }
      virtual bool VisData(const std::string& name, std::vector<double>& data)
      {
        return So_sh18::VisData(name, data);
      }
      virtual int EvaluateNeumann(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1 = NULL)
      {
        return So_sh18::EvaluateNeumann(params, discretization, condition, lm, elevec1, elemat1);
      }

      //! Deep copy this instance of Solid3 and return pointer to the copy
      //!
      //! The Clone() method is used from the virtual base class Element in cases
      //! where the type of the derived class is unknown and a copy-ctor is needed
      DRT::Element* Clone() const;

      //! Destructor
      virtual ~So_sh18Plast();

      //! Return unique ParObject id
      //!
      //! every class implementing ParObject needs a unique id defined at the top of
      //! this file.
      virtual int UniqueParObjectId() const
      {
        return So_sh18PlastType::Instance().UniqueParObjectId();
      }

      //! Pack this class so it can be communicated
      //! Pack and \ref Unpack are used to communicate this element
      virtual void Pack(DRT::PackBuffer& data) const;

      //! Unpack data from a char vector into this class
      //! Pack and \ref Unpack are used to communicate this element
      virtual void Unpack(const std::vector<char>& data);

      //! Print this element
      void Print(std::ostream& os) const;

      //! return elementtype
      virtual So_sh18PlastType& ElementType() const { return So_sh18PlastType::Instance(); }

      //! read input for this element
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);

      //! synchronize the eas variables in the two base-classes
      void SyncEAS();

      //! evaluate an element
      //! evaluate element stiffness, mass, internal forces, etc.
      //!
      //! if NULL on input, the controlling method does not expect the element
      //!  to fill these matrices or vectors.
      //!
      //!  \return 0 if successful, negative otherwise
      int Evaluate(
          Teuchos::ParameterList&
              params,  //!< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  //!< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la,      //!< location array for de-assembly
          Epetra_SerialDenseMatrix&
              elemat1_epetra,  //!< (stiffness-)matrix to be filled by element.
          Epetra_SerialDenseMatrix& elemat2_epetra,  //!< (mass-)matrix to be filled by element.
          Epetra_SerialDenseVector&
              elevec1_epetra,  //!< (internal force-)vector to be filled by element
          Epetra_SerialDenseVector& elevec2_epetra,  //!< vector to be filled by element
          Epetra_SerialDenseVector& elevec3_epetra   //!< vector to be filled by element
      )
      {
        return DRT::ELEMENTS::So3_Plast<DRT::Element::hex18>::Evaluate(params, discretization, la,
            elemat1_epetra, elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra);
      }


     private:
      // don't want = operator
      So_sh18Plast& operator=(const So_sh18Plast& old);

     protected:
      //! Calculate nonlinear stiffness and mass matrix with condensed plastic matrices
      virtual void nln_stiffmass(std::vector<double>& disp,  // current displacements
          std::vector<double>& vel,                          // current velocities
          std::vector<double>& temp,                         // current temperatures
          LINALG::Matrix<numdofperelement_, numdofperelement_>*
              stiffmatrix,  // element stiffness matrix
          LINALG::Matrix<numdofperelement_, numdofperelement_>* massmatrix,  // element mass matrix
          LINALG::Matrix<numdofperelement_, 1>* force,      // element internal force vector
          LINALG::Matrix<numgpt_post, numstr_>* elestress,  // stresses at GP
          LINALG::Matrix<numgpt_post, numstr_>* elestrain,  // strains at GP
          Teuchos::ParameterList& params,                   // algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,            // stress output option
          const INPAR::STR::StrainType iostrain             // strain output option
      );

      //! don't want sosh18 nlnstiffmass
      virtual void nlnstiffmass(std::vector<int>& lm,               ///< location matrix
          std::vector<double>& disp,                                ///< current displacements
          std::vector<double>& residual,                            ///< current residual displ
          LINALG::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* stiffmatrix,  ///< element stiffness matrix
          LINALG::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* massmatrix,   ///< element mass matrix
          LINALG::Matrix<NUMDOF_SOH18, 1>* force,  ///< element internal force vector
          LINALG::Matrix<NUMGPT_SOH18, MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          LINALG::Matrix<NUMGPT_SOH18, MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          Teuchos::ParameterList& params,         ///< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,  ///< stress output option
          const INPAR::STR::StrainType iostrain   ///< strain output option
      )
      {
        dserror("don't want this");
      }
    };

  }  // namespace ELEMENTS



}  // namespace DRT

#endif  // SO3_SSN_PLAST_SOSH18_H
