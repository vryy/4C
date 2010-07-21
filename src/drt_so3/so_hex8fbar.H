/*!----------------------------------------------------------------------
\file so_hex8fbarH

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET
#ifndef SOH8FBAR_H
#define SOH8FBAR_H

#include "so_hex8.H"

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
class InvDesign;

class So_hex8fbarType : public DRT::ElementType
{
public:

  std::string Name() const { return "So_hex8fbarType"; }

  static So_hex8fbarType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const string eletype,
                                             const string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual int Initialize(DRT::Discretization& dis);

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static So_hex8fbarType instance_;
};

/*!
\brief A C++ version of a 8-node hex solid element with F-Bar modification

A structural 8-node hexahedral solid displacement element for large deformations
and (near)-incompressibility. F-bar techniques are used to avoid volumetric locking.
As its discretization is fixed many data structures are evaluated just once and kept
for performance. It heavily uses Epetra objects and methods and therefore relies
on their performance.

Refer also to the HiWi report of Stefanos Tsoukalas, 2010

\author popp (popp@lnm.mw.tum,de)
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
  DRT::Element* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return So_hex8fbarType::Instance().UniqueParObjectId(); }

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Unpack(const vector<char>& data);

  /*!
  \brief Destructor

  */
  virtual ~So_hex8fbar();

  //@}

  //! @name Access methods

  /*!
  \brief Print this element
  */
  void Print(ostream& os) const;

  virtual DRT::ElementType & ElementType() const
  { return So_hex8fbarType::Instance(); }

  //@}

  //! @name Input and Creation

  /*!
  \brief Read input for this element
  */

  //virtual bool ReadElement();

  /*!
  \brief Read input for this element
  */
  bool ReadElement(const std::string& eletype,
                   const std::string& distype,
                   DRT::INPUT::LineDefinition* linedef);

  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate So_hex8fbar element stiffness, mass, internal forces, etc.

  If NULL on input, the controlling method does not expect the element
  to fill these matrices or vectors.

  \return 0 if successful, negative otherwise
  */
  int Evaluate
  (
    ParameterList& params,  ///< ParameterList for communication between control routine and elements
    DRT::Discretization&      discretization,  ///< pointer to discretization for de-assembly
    vector<int>&              lm, ///< location matrix for de-assembly
    Epetra_SerialDenseMatrix& elemat1,  ///< (stiffness-)matrix to be filled by element.
    Epetra_SerialDenseMatrix& elemat2,  ///< (mass-)matrix to be filled by element.
    Epetra_SerialDenseVector& elevec1,  ///< (internal force-)vector to be filled by element
    Epetra_SerialDenseVector& elevec2,  ///< vector to be filled by element
    Epetra_SerialDenseVector& elevec3  ///< vector to be filled by element
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
  virtual int EvaluateNeumann
  (
    ParameterList& params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    vector<int>&              lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1 = NULL
  );

  //@}

protected:

  //! don't want = operator
  So_hex8fbar& operator = (const So_hex8fbar& old);

  //! Calculate nonlinear stiffness and mass matrix
  virtual void soh8fbar_nlnstiffmass(
		vector<int>&              lm,             ///< location matrix
		vector<double>&           disp,           ///< current displacements
		vector<double>&           residual,       ///< current residual displ
		LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* stiffmatrix,    ///< element stiffness matrix
		LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* massmatrix,     ///< element mass matrix
		LINALG::Matrix<NUMDOF_SOH8,1>* force,          ///< element internal force vector
		LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestress,      ///< stresses at GP
		LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestrain,      ///< strains at GP
		ParameterList&            params,         ///< algorithmic parameters e.g. time
		const INPAR::STR::StressType   iostress,  ///< stress output option
		const INPAR::STR::StrainType   iostrain); ///< strain output option

}; // class So_hex8fbar


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================



} // namespace ELEMENTS
} // namespace DRT

#endif  // #ifndef SOH8_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
