/*!----------------------------------------------------------------------
\file so_sh8.H

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#ifndef SO_SH8_H
#define SO_SH8_H
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
#include "../drt_inpar/inpar_structure.H"


// forward declarations
struct _SOH8_DATA;

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
class Sosh8Register;

/*!
\brief A C++ 8-node Solid-Shell element inherited from so_hex8

The Solid-Shell element technology is based on the work of
(1) Vu-Quoc, Tan: "Optimal solid shells for non-linear analyses
                   of multilayer composites", CMAME 2003
(2) Klinkel, Gruttmann, Wagner: "A robust non-linear solid shell element
                                 based on a mixed variational fromulation"

Refer also to the Semesterarbeit of Alexander Popp, 2006

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
  \param owner : elements owning processor
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

  //bool ReadElement();

  /*!
  \brief Read input for this element
  */
  bool ReadElement(const std::string& eletype,
                   const std::string& distype,
                   DRT::INPUT::LineDefinition* linedef);

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
  \param discretization : pointer to discretization for de-assembly
  \param lm (in)        : location matrix for de-assembly
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
    globx,    ///< global x
    globy,    ///< global y
    globz,    ///< global z
    autoj,    ///< find automatically by Jacobian
    autor,    ///< automatically set to x
    autos,    ///< automatically set to y
    autot,    ///< automatically set to z
    enfor,    ///< read-in r-direction is rearranged to t-dir
    enfos,    ///< read-in s-direction is rearranged to t-dir
    enfot,    ///< read-in t-direction stays t-dir
    undefined,///< no clear direction identified
    none      ///< no rearrangement
  };
  
  enum ANSType
  {
    anssosh8,
    ansnone
  };

private:

  // don't want = operator
  So_sh8& operator = (const So_sh8& old);

protected:

  static const int num_sp = 8;              ///< number of ANS sampling points, here 8
  static const int num_ans = 3;             ///< number of modified ANS strains (E_rt,E_st,E_tt), here 3
  //! shell-thickness direction
  ThicknessDirection                             thickdir_;
  
  ANSType                                        anstype_;

  //! in case of changed "thin" direction this is true
  bool nodes_rearranged_;

  //! vector in thickness direction for compatibility with sosh8
  vector<double>                                   thickvec_;

  //! Compute stiffness and mass matrix
  void sosh8_nlnstiffmass(
      vector<int>&              lm,             ///< location matrix
      vector<double>&           disp,           ///< current displacements
      vector<double>&           residual,       ///< current residual displ
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* stiffmatrix,    ///< element stiffness matrix
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* massmatrix,     ///< element mass matrix
      LINALG::Matrix<NUMDOF_SOH8,1>* force,          ///< element internal force vector
      LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestress,      ///< stresses at GP
      LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestrain,      ///< strains at GP
      ParameterList&            params,         ///< algorithmic parameters e.g. time
      const INPAR::STR::StressType  iostress,   ///< stress output option
      const INPAR::STR::StrainType  iostrain);  ///< strain output option

  //! Evaluate all ANS related data at the ANS sampling points
  void sosh8_anssetup(
            const LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8>& xrefe, ///< material element coords
            const LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8>& xcurr, ///< current element coords
            vector<LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> >** deriv_sp,   ///< derivs eval. at all sampling points
            vector<LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> >& jac_sps,     ///< jac at all sampling points
            vector<LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> >& jac_cur_sps, ///< current jac at all sampling points
            LINALG::Matrix<num_ans*num_sp,NUMDOF_SOH8>& B_ans_loc); ///< modified B

  //! Evaluate transformation matrix T (parameter->material) at gp
  void sosh8_evaluateT(const LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>& jac,///< actual jacobian
                       LINALG::Matrix<NUMSTR_SOH8,NUMSTR_SOH8>& TinvT);///< T^{-T}

  //! Return true Cauchy-stress at gausspoint
  void sosh8_Cauchy(LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestress,
                    const int gp,
                    const LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>& defgrd,
                    const LINALG::Matrix<NUMSTR_SOH8,1>& glstrain,
                    const LINALG::Matrix<NUMSTR_SOH8,1>& stress);

  //! Find "thin"=thickness direction
  ThicknessDirection sosh8_findthickdir();

  //! Find parametric co-ordinate which directs in enforced thickness direction
  ThicknessDirection sosh8_enfthickdir(
    LINALG::Matrix<NUMDIM_SOH8,1>& glothickdir  ///< global direction of enforced thickness direction
    );

  //! return thickness direction
  const vector<double> GetThickvec(){return thickvec_;};

  //! Debug gmsh-plot to check thickness direction
  void sosh8_gmshplotlabeledelement(const int LabelIds[NUMNOD_SOH8]);

}; // class So_sh8


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
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

  The method performs all steps needed to init the solid-shell8 elements.
  This includes defining the thickness direction of the solid-shell by looking
  at its geometric shape in reference configuration. Thus, the loop looks a bit
  complicated because a "not-thin" So_sh8 is morphed to a So_hex8 where the
  discretization is not filled-complete anymore.
  */
  int Initialize(DRT::Discretization& dis);

  //@}


  //! debug plot for defining thickness direction
  void sosh8_gmshplotdis(const DRT::Discretization& dis);

}; // class Sosh8Register


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================



} // namespace ELEMENTS
} // namespace DRT


#endif  // #ifndef SO_SH8_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
