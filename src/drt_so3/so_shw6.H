#ifdef D_SOH8
#ifdef CCADISCRET

using namespace std;
#include "so_weg6.H"

namespace DRT
{
// forward declarations
class Discretization;

namespace Elements
{
// forward declarations
class Soshw6Register;

/*!
\class So_shw6

\brief A C++ 6-node wedge Solid-Shell element inherited from so_weg6

\author frenzel (frenzel@lnm.mw.tum.de)
*/
class So_shw6 : public So_weg6
{
public:
  //! @name Friends
  friend class Soshw6Register;
  friend class Soweg6Surface;
  friend class Soweg6Line;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  So_shw6(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  So_shw6(const So_shw6& old);

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
  virtual int UniqueParObjectId() const { return ParObject_So_shw6;}

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
  virtual ~So_shw6();

  //@}

  //! @name Acess methods


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
  virtual bool ReadElement();

  //@}

  //! @name Evaluation



private:
  
  // don't want = operator
  So_shw6& operator = (const So_shw6& old);

  virtual void soshw6_nlnstiffmass(vector<int>&              lm,
                                 vector<double>&           disp,
                                 vector<double>&           residual,
                                 Epetra_SerialDenseMatrix* stiffmatrix,
                                 Epetra_SerialDenseMatrix* massmatrix,
                                 Epetra_SerialDenseVector* force,
                                 const double              time);

  void soshw6_anssetup(
            const int numsp,              // number of sampling points
            const int numans,             // number of ans strains
            const Epetra_SerialDenseMatrix& xrefe, // material element coords
            const Epetra_SerialDenseMatrix& xcurr, // current element coords
            Epetra_SerialDenseMatrix** deriv_sp,   // derivs eval. at all sampling points
            Epetra_SerialDenseMatrix& jac_sps,     // jac at all sampling points
            Epetra_SerialDenseMatrix& jac_cur_sps, // current jac at all sampling points
            Epetra_SerialDenseMatrix& B_ans_loc); // modified B
  void soshw6_evaluateT(const Epetra_SerialDenseMatrix jac,
                        Epetra_SerialDenseMatrix& TinvT);

}; // class So_shw6

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\class Soshw6Register

\brief A register for Solid-Shell wedge6 element

\author frenzel (frenzel@lnm.mw.tum.de)
*/
class Soshw6Register : public Sow6Register
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Soshw6Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Soh8Register

  */
  Soshw6Register(const DRT::Elements::Soshw6Register& old);

  /*!
  \brief Deep copy this instance of Soh8Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Elements::Soshw6Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Soshw6Register;}

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
  virtual ~Soshw6Register();

  //@}

  //! @name Acess methods

  /*!
  \brief Print this element
  */
  void Print(ostream& os) const;

  //@}



private:


}; // class Sow6Register

} // namespace Elements
} // namespace DRT


#endif /*SO_WEG6_H_*/
#endif
