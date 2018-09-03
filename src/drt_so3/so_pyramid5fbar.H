/*!----------------------------------------------------------------------
\file so_pyramid5fbar.H
\brief
\level 2

<pre>
\maintainer Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15271
</pre>

*----------------------------------------------------------------------*/
#ifndef SOP5FBAR_H
#define SOP5FBAR_H

#include "so_pyramid5.H"

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    // forward declarations
    class InvDesign;

    class So_pyramid5fbarType : public DRT::ElementType
    {
     public:
      std::string Name() const { return "So_pyramid5fbarType"; }

      static So_pyramid5fbarType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

      virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
          const double* x0, int numdf, int dimns);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_pyramid5fbarType instance_;
    };



    class So_pyramid5fbar : public So_pyramid5
    {
     public:
      //! @name Friends
      friend class So_pyramid5fbarType;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      So_pyramid5fbar(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      So_pyramid5fbar(const So_pyramid5fbar& old);

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const
      {
        return So_pyramid5fbarType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      virtual void Pack(DRT::PackBuffer& data) const;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      virtual void Unpack(const std::vector<char>& data);

      /*!
      \brief Destructor

      */
      virtual ~So_pyramid5fbar();

      //@}

      //! @name Access methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const { return So_pyramid5fbarType::Instance(); }

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate So_pyramid5fbar element stiffness, mass, internal forces, etc.

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
      );

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
      virtual int EvaluateNeumann(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1 = NULL);

      //@}

     protected:
      //! don't want = operator
      So_pyramid5fbar& operator=(const So_pyramid5fbar& old);

      // compute Jacobian mapping wrt to deformed configuration
      virtual void UpdateJacobianMapping(
          const std::vector<double>& disp, DRT::ELEMENTS::PreStress& prestress);
      // compute defgrd in all gp for given disp
      virtual void DefGradient(const std::vector<double>& disp, Epetra_SerialDenseMatrix& gpdefgrd,
          DRT::ELEMENTS::PreStress& prestress);

      //! Calculate nonlinear stiffness and mass matrix
      virtual void nlnstiffmass(std::vector<int>& lm,  ///< location matrix
          std::vector<double>& disp,                   ///< current displacements
          std::vector<double>* vel,                    ///< current velocities
          std::vector<double>* acc,                    ///< current accelerations
          std::vector<double>& residual,               ///< current residual displ
          std::vector<double>& dispmat,                ///< current material displacements
          LINALG::Matrix<NUMDOF_SOP5, NUMDOF_SOP5>* stiffmatrix,  ///< element stiffness matrix
          LINALG::Matrix<NUMDOF_SOP5, NUMDOF_SOP5>* massmatrix,   ///< element mass matrix
          LINALG::Matrix<NUMDOF_SOP5, 1>* force,                  ///< element internal force vector
          LINALG::Matrix<NUMDOF_SOP5, 1>* forceinert,             ///< element inertial force vector
          LINALG::Matrix<NUMDOF_SOP5, 1>* force_str,  ///< element structural force vector
          LINALG::Matrix<NUMGPT_SOP5, MAT::NUM_STRESS_3D>* elestress,    ///< stresses at GP
          LINALG::Matrix<NUMGPT_SOP5, MAT::NUM_STRESS_3D>* elestrain,    ///< strains at GP
          LINALG::Matrix<NUMGPT_SOP5, MAT::NUM_STRESS_3D>* eleplstrain,  ///< plastic strains at GP
          Teuchos::ParameterList& params,            ///< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,     ///< stress output option
          const INPAR::STR::StrainType iostrain,     ///< strain output option
          const INPAR::STR::StrainType ioplstrain);  ///< plastic strain output option

      //! init the inverse of the jacobian and its determinant in the material configuration
      virtual void InitJacobianMapping();

    };  // class So_pyramid5fbar

    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================

  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // #ifndef SOP5FBAR_H
