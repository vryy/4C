/*!----------------------------------------------------------------------
\file so_sh8.H

\class DRT::ELEMENTS::So_sh8

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

namespace ELEMENTS
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
  friend class Soh8Line;

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

  /*!
  \brief Allocate and return a DesignElementRegister class

  Every implementation of DRT::Element must also implement
  an DRT::ElementRegister. This method allocates and returns
  an instance of a matching DRT::ElementRegister.
  */
  virtual RefCountPtr<DRT::ElementRegister> ElementRegister() const;

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

  //! definition of shell-thickness direction
  enum ThicknessDirection
  {
    globx,    // global x
    globy,    // global y
    globz,    // global z
    autoj,    // find automatically by Jacobian
    autor,    // automatically set to x
    autos,    // automatically set to y
    autot,    // automatically set to z
    none,     // no rearrangement
  };

  int Initialize_numbers(DRT::Discretization& dis);

private:

  // don't want = operator
  So_sh8& operator = (const So_sh8& old);

  //! shell-thickness direction
  ThicknessDirection                             thickdir_;

  //! original (input) NodeIds, stored for eventual change to define thickness dir.
  int inp_nodeIds_[NUMNOD_SOH8];

  bool nodes_rearranged_;

  void sosh8_nlnstiffmass(vector<int>&              lm,
                         vector<double>&           disp,
                         vector<double>&           residual,
                         Epetra_SerialDenseMatrix* stiffmatrix,
                         Epetra_SerialDenseMatrix* massmatrix,
                         Epetra_SerialDenseVector* force,
                         Epetra_SerialDenseMatrix* elestress);

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

  ThicknessDirection sosh8_findthickdir();
  
  void sosh8_gmshplotlabeledelement(const int LabelIds[NUMNOD_SOH8]);
  
  void sosh8_setcylinderfiberdirection(const Epetra_SerialDenseMatrix xrefe);

 
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
  Sosh8Register(const DRT::ELEMENTS::Sosh8Register& old);

  /*!
  \brief Deep copy this instance of Soh8Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::Sosh8Register* Clone() const;

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

  /*!
  \brief Call Element init routines

  The method performs all steps needed to init the solid-shell8 elements,
  e.g. defines shell-thickness direction
  */
  int Initialize(DRT::Discretization& dis);

  //@}



  void sosh8_gmshplotdis(const DRT::Discretization& dis);

}; // class Sosh8Register


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================



} // namespace ELEMENTS
} // namespace DRT



#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
