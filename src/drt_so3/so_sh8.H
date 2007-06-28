/*!----------------------------------------------------------------------
\file so_sh8.H

\class DRT::Elements::So_sh8

\brief A C++ 8-node Solid-Shell element inherited from so_hex8

The Solid-Shell element technology is based on the work of
(1) Vu-Quoc, Tan: "Optimal solid shells for non-linear analyses
                   of multilayer composites", CMAME 2003
(2) Klinkel, Gruttmann, Wagner: "A robust non-linear solid shell element
                                 based on a mixed variational fromulation"

Refer also to the Semesterarbeit of Alexander Popp, 2006

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOH8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"
#include "so_hex8.H"

// Several parameters which are fixed for Solid-Shell8
#define NUMNOD_SOH8 8      // number of nodes
#define NODDOF_SOH8 3      // number of dofs per node
#define NUMDOF_SOH8 24     // total dofs per element
#define NUMSTR_SOH8 6      // number of strains/stresses
#define NUMGPT_SOH8 8      // total gauss points per element
#define NUMDIM_SOH8 3      // number of dimensions

// forward declarations
struct _MATERIAL;
struct _SOH8_DATA;

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{
// forward declarations
class Discretization;

namespace Elements
{
// forward declarations
class Sosh8Register;

/*!
\class So_sh8

\brief A C++ 8-node Solid-Shell element inherited from so_hex8

\author frenzel (frenzel@lnm.mw.tum.de)
*/
class So_sh8 : public So_hex8
{
public:

  //! @name Friends
  friend class Sosh8Register;
  friend class Soh8Surface;
  friend class Soh83Line;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  So_sh8(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  So_sh8(const So_sh8& old);

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
  virtual int UniqueParObjectId() const { return ParObject_So_sh8;}

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
  \brief Print this element
  */
  void Print(ostream& os) const;

  //@}

  //! @name Input and Creation

  /*!
  \brief Read input for this element
  */

  bool ReadElement();

  /*!
  \brief Destructor

  */
  virtual ~So_sh8();


  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate so_sh8 element stiffness, mass, internal forces, etc.

  \param params (in/out): ParameterList for communication between control routine
                          and elements
  \param elemat1 (out)  : (stiffness-)matrix to be filled by element. If NULL on input,
                          the controling method does not expect the element to fill
                          this matrix.
  \param elemat2 (out)  : (mass-)matrix to be filled by element. If NULL on input,
                          the controling method does not expect the element to fill
                          this matrix.
  \param elevec1 (out)  : (internal force-)vector to be filled by element. If NULL on input,
                          the controlling method does not expect the element
                          to fill this vector
  \param elevec2 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not expect the element
                          to fill this vector
  \param elevec3 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not expect the element
                          to fill this vector
  \return 0 if successful, negative otherwise
  */
  int Evaluate(ParameterList& params,
               DRT::Discretization&      discretization,
               vector<int>&              lm,
               Epetra_SerialDenseMatrix& elemat1,
               Epetra_SerialDenseMatrix& elemat2,
               Epetra_SerialDenseVector& elevec1,
               Epetra_SerialDenseVector& elevec2,
               Epetra_SerialDenseVector& elevec3);



private:

  // don't want = operator
  So_sh8& operator = (const So_sh8& old);

  void sosh8_nlnstiffmass(vector<int>&              lm,
                         vector<double>&           disp,
                         vector<double>&           residual,
                         Epetra_SerialDenseMatrix* stiffmatrix,
                         Epetra_SerialDenseMatrix* massmatrix,
                         Epetra_SerialDenseVector* force);

  void sosh8_anssetup(const int numsp,              // number of sampling points
          const int numans,             // number of ans strains
          const Epetra_SerialDenseMatrix& xrefe, // material element coords
          const Epetra_SerialDenseMatrix& xcurr, // current element coords
          Epetra_SerialDenseMatrix** deriv_sp,   // derivs eval. at all sampling points
          Epetra_SerialDenseMatrix& jac_sps,     // jac at all sampling points
          Epetra_SerialDenseMatrix& jac_cur_sps, // current jac at all sampling points
          Epetra_SerialDenseMatrix&  B_ans_loc); // modified B

  void sosh8_evaluateT(const Epetra_SerialDenseMatrix jac,  // actual jacobian
                       Epetra_SerialDenseMatrix& TinvT);    // corr. T^{-T}

}; // class So_sh8


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\class Sosh8Register

\brief A register for Solid-Shell8 element, inherited from Soh8

\author frenzel (frenzel@lnm.mw.tum.de)
*/
class Sosh8Register : public Soh8Register
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Sosh8Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Soh8Register

  */
  Sosh8Register(const DRT::Elements::Sosh8Register& old);

  /*!
  \brief Deep copy this instance of Soh8Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Elements::Sosh8Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Sosh8Register;}

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
  virtual ~Sosh8Register();

  //@}

  //! @name Acess methods

  /*!
  \brief Print this element
  */
  void Print(ostream& os) const;

private:

}; // class Sosh8Register


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================



} // namespace Elements
} // namespace DRT



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
