/*!----------------------------------------------------------------------
\file so_shw6.H
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

using namespace std;
#include "so_weg6.H"

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
class Soshw6Register;

/*!
\brief A C++ 6-node wedge Solid-Shell element inherited from so_weg6

This element adopts the ANS approach to improve shear and curvature-thickness
(trapezoidal) locking. However, being triangular inplane, the ANS evaluation
points of the shear strains are choosen just on two of three surfaces (See
also DRT::ELEMENTS::So_shw6::soshw6_anssetup).
Thus, the element becomes mesh dependent (where in material space lies my
parameter space), but this is in favour of less pronounced shear locking
compared to a invariant formulation (locking free just for certain mesh configurations).

See Diss. Frank Koschnick for details.

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
  \param owner : elements owning processor
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

  /*! \brief EAS technology enhanecement types of so_shw6
   * Solid-Shell Wedge6 has EAS enhancement of GL-strains to avoid locking
   * \param soshw6_poissonthick: 1 parameter for correct kinematic of solid-shell,
   *                             allows linear strain in thick-dir
   * \param soshw6_none:         no EAS
   */
  enum EASType
  {
    soshw6_easpoisthick = 1,
    soshw6_easnone = 0
  };

  //! type of EAS technology
  EASType                                             eastype_;

  //! number of EAS parameters (alphas), defined by 'EASType"
  int                                                 neas_;

  //! don't want = operator
  So_shw6& operator = (const So_shw6& old);

  /*!
  \brief Evaluate an element

  Evaluate so_hex8 element stiffness, mass, internal forces, etc.

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
  virtual int Evaluate(ParameterList& params,
               DRT::Discretization&      discretization,
               vector<int>&              lm,
               Epetra_SerialDenseMatrix& elemat1,
               Epetra_SerialDenseMatrix& elemat2,
               Epetra_SerialDenseVector& elevec1,
               Epetra_SerialDenseVector& elevec2,
               Epetra_SerialDenseVector& elevec3);

  //! Compute stiffness and mass matrix
  virtual void soshw6_nlnstiffmass(
          vector<int>&              lm,             // location matrix
          vector<double>&           disp,           // current displacements
          vector<double>&           residual,       // current residual displ
          LINALG::FixedSizeSerialDenseMatrix<NUMDOF_WEG6,NUMDOF_WEG6>* stiffmatrix, // element stiffness matrix
          LINALG::FixedSizeSerialDenseMatrix<NUMDOF_WEG6,NUMDOF_WEG6>* massmatrix,  // element mass matrix
          LINALG::FixedSizeSerialDenseMatrix<NUMDOF_WEG6,1>* force,                 // element internal force vector
          LINALG::FixedSizeSerialDenseMatrix<NUMGPT_WEG6,NUMSTR_WEG6>* elestress,   // stresses at GP
          LINALG::FixedSizeSerialDenseMatrix<NUMGPT_WEG6,NUMSTR_WEG6>* elestrain,   // strains at GP
          ParameterList&            params,         // algorithmic parameters e.g. time
          const bool                cauchy=false);         // stress output option

  static const int num_sp = 5;       ///< number of ANS sampling points
  static const int num_ans = 3;      ///< number of modified ANS strains (E_rt,E_st,E_tt)

  //! Setup ANS interpolation (shear and trapezoidal)
  void soshw6_anssetup(
          const LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,NUMDIM_WEG6>& xrefe, ///< material element coords
          const LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,NUMDIM_WEG6>& xcurr, ///< current element coords
          vector<LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMNOD_WEG6> >** deriv_sp,   ///< derivs eval. at all sampling points
          vector<LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMDIM_WEG6> >& jac_sps,     ///< jac at all sampling points
          vector<LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMDIM_WEG6> >& jac_cur_sps, ///< current jac at all sampling points
          LINALG::FixedSizeSerialDenseMatrix<num_ans*num_sp,NUMDOF_WEG6>& B_ans_loc); ///< modified B

  //! Transformation matrix parameter->material space
  void soshw6_evaluateT(const LINALG::FixedSizeSerialDenseMatrix<NUMDIM_WEG6,NUMDIM_WEG6>& jac,
                        LINALG::FixedSizeSerialDenseMatrix<NUMSTR_WEG6,NUMSTR_WEG6>& TinvT);

  //! EAS technology, init
  void soshw6_easinit();

  //! EAS technology, setup necessary data
  void soshw6_eassetup(
            vector<Epetra_SerialDenseMatrix>** M_GP,    // M-matrix evaluated at GPs
            double& detJ0,                      // det of Jacobian at origin
            LINALG::FixedSizeSerialDenseMatrix<NUMSTR_WEG6,NUMSTR_WEG6>& T0invT,   // maps M(origin) local to global
            const LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,NUMDIM_WEG6>& xrefe);    // material element coords


}; // class So_shw6

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
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

  Makes a deep copy of an WEG6Register

  */
  Soshw6Register(const DRT::ELEMENTS::Soshw6Register& old);

  /*!
  \brief Deep copy this instance of WEG6Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::Soshw6Register* Clone() const;

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

} // namespace ELEMENTS
} // namespace DRT


#endif /*SO_WEG6_H_*/
#endif
