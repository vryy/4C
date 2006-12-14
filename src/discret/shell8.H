/*!----------------------------------------------------------------------
\file shell8.H

\class DRT::Elements::Shell8

\brief A C++ wrapper for the shell8 element

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef SHELL8_H
#define SHELL8_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
using namespace std;
using namespace Teuchos;

#include "drt_element.H"
#include "drt_elementregister.H"


// forward declarations
struct _MATERIAL;
struct _S8_DATA;

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{
// forward declarations
class Discretization;
class Shell8Register;

namespace Elements
{

/*!
\class Shell8

\brief A C++ wrapper for the shell8 element

\author gee (gee@lnm.mw.tum.de)
*/
class Shell8 : public DRT::Element
{
public:

  /*!
  \brief Shell8Register is a friend of Shell8
  */
  friend class Shell8Register;

  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  \param id : A unique global id
  */
  Shell8(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  
  */
  Shell8(const Shell8& old);

  /*!
  \brief Deep copy this instance of Shell8 and return pointer to the copy
  
  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed
  
  */
  DRT::Element* Clone() const;

  /*!
  \brief Return unique ParObject id
  
  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Shell8;}

  /*!
  \brief Pack this class so it can be communicated
  
  \ref Pack and \ref Unpack are used to communicate this element
  
  */
  virtual const char* Pack(int& size) const;

  /*!
  \brief Unpack data from a char vector into this class
  
  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual bool Unpack(const char* data);

  /*!
  \brief Destructor

  */
  virtual ~Shell8();

  //@}

  // @{ \name Acess methods

  
  /*!
  \brief Get number of degrees of freedom of a certain node
         (implements pure virtual DRT::Element)
  
  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  virtual int NumDofPerNode(const DRT::Node& node) const {return 6;}

  /*!
  \brief Get number of degrees of freedom per element
         (implements pure virtual DRT::Element)
  
  The element decides how many element degrees of freedom it has.
  It can redecide along the way of a simulation.
  
  \note Element degrees of freedom mentioned here are dofs that are visible
        at the level of the total system of equations. Purely internal
        element dofs that are condensed internally should NOT be considered.
  */
  virtual int NumDofPerElement() const { return 0; }

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

  // @{ \name Other

  /*!
  \brief Read input for this element
  */
  bool ReadElement();

  /*!
  \brief Evaluate an element
  
  Evaluate shell8 element stuffness, mass, internal forces etc
  
  \param params (in/out): ParameterList for communication between control routine
                          and elements
  \param elemat1 (out)  : matrix to be filled by element. If NULL on input, 
                          the controling method does not epxect the element to fill
                          this matrix.
  \param elemat2 (out)  : matrix to be filled by element. If NULL on input, 
                          the controling method does not epxect the element to fill
                          this matrix.
  \param elevec1 (out)  : vector to be filled by element. If NULL on input, 
                          the controlling method does not epxect the element
                          to fill this vector 
  \param elevec2 (out)  : vector to be filled by element. If NULL on input, 
                          the controlling method does not epxect the element
                          to fill this vector 
  \param elevec3 (out)  : vector to be filled by element. If NULL on input, 
                          the controlling method does not epxect the element
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

  //@}


private:

  // action parameters recognized by shell8
  enum ActionType
  {
    none,
    calc_struct_linstiff,
    calc_struct_nlnstiff,
    calc_struct_internalforce,
    calc_struct_linstiffmass,
    calc_struct_nlnstiffmass,
    calc_struct_stress,
    calc_struct_eleload,
    calc_struct_fsiload
  };

  // element specific data
  enum ForceType
  {
    s8_none,
    s8_xyz,
    s8_rst,
    s8_rst_ortho
  };
  ForceType                forcetype_;     // type of force calculation
  double                   thickness_;     // shell thickness
  int                      ngp_[3];        // no. gaussian points in each spatial dimension
  int                      ngptri_;        // no. gaussian points for triangle elements
  int                      nhyb_;          // no. of internal dofs
  int                      eas_[5];        // eas parameters
  int                      ans_;           // type of ANS: 0=none 1=Q 2=T 3=QT  
  double                   sdc_;           // scaled director conditioning
  int                      material_;      // number of the material law
  DRT::Container           data_;          // container containing all kind of shell8 data
  
  // internal calculation methods

  // don't want = operator
  Shell8 operator = (const Shell8& old);
  
  void s8_nlnstiffmass(vector<int>&              lm, 
                       vector<double>&           disp, 
                       vector<double>&           residual,
                       Epetra_SerialDenseMatrix* stiffmatrix,
                       Epetra_SerialDenseMatrix* massmatrix,
                       Epetra_SerialDenseVector* force,
                       struct _MATERIAL*         material);
  
  void s8_integration_points(struct _S8_DATA& data);
  
  void s8BtDB(Epetra_SerialDenseMatrix& estif, const Epetra_SerialDenseMatrix& bop,
              const double D[][12], const int iel, const int numdf, const double weight);
  
  void s8tvkg(Epetra_SerialDenseMatrix& estif, double stress_r[], 
              const vector<double>& funct, const Epetra_SerialDenseMatrix& deriv,
              const int numdf, const int iel, const double weight,
              const double e1, const double e2);
  
  void s8anstvkg(Epetra_SerialDenseMatrix& estif, double stress_r[], 
                 const vector<double>& funct, const Epetra_SerialDenseMatrix& deriv,
                 const int numdf, const int iel, const double weight,
                 const double e1, const double e2,
                 const double frq[], const double fsq[], 
                 const vector<double> funct1q[], const vector<double> funct2q[],
                 const Epetra_SerialDenseMatrix deriv1q[], const Epetra_SerialDenseMatrix deriv2q[],
                 const int ansq, const int nsansq);

  void s8intforce(Epetra_SerialDenseVector& intforce, const double stress_r[],
                  const Epetra_SerialDenseMatrix& bop, const int iel, 
                  const int numdf, const int nstress_r, const double weight);

  void s8tmas(const vector<double>& funct, const vector<double>& thick,
              Epetra_SerialDenseMatrix& emass, const int iel, const int numdf,
              const double facv, const double facw, const double facvw);

  void s8tmat(struct _MATERIAL* material, 
              double stress[], double strain[], double** C,
              double gmkovc[][3], double gmkonc[][3], 
              double gmkovr[][3], double gmkonr[][3],
              double gkovc[][3], double gkonc[][3],
              double gkovr[][3], double gkonr[][3], 
              const double detc, const double detr, 
              const double e3, const int option, 
              const int ngauss);
  
  void s8tvma(double D[][12], double** C, double stress[], double stress_r[], 
              const double e3, const double fact, const double condfac);
  
  void s8_YpluseqAx(vector<double>& y, const Epetra_SerialDenseMatrix& A,
                    const vector<double>& x, const double factor, 
                    const bool init);
  void s8_YpluseqAx(Epetra_SerialDenseVector& y, const Epetra_SerialDenseMatrix& A,
                    const vector<double>& x, const double factor, 
                    const bool init);
  
  void s8inv3(double a[][3], double* det);
  
  void s8trans3(double a[][3]);
  
  void s8unvc(double* enorm, double vec[], const int n);
  
  void s8matmatdense(Epetra_SerialDenseMatrix& R, 
                     const Epetra_SerialDenseMatrix& A,
                     const Epetra_SerialDenseMatrix& B,
                     const int ni,
                     const int nk,
                     const int nj,
                     const int init,
                     const double factor);
  void s8matmatdense(Epetra_SerialDenseMatrix& R, 
                     const double A[][12],
                     const Epetra_SerialDenseMatrix& B,
                     const int ni,
                     const int nk,
                     const int nj,
                     const int init,
                     const double factor);
  void s8mattrnmatdense(Epetra_SerialDenseMatrix& R,
                        const Epetra_SerialDenseMatrix& A,
                        const Epetra_SerialDenseMatrix& B,
                        const int ni,
                        const int nk,
                        const int nj,
                        const int init,
                        const double factor);
  void s8mattrnvecdense(vector<double>& r,
                        const Epetra_SerialDenseMatrix& A, 
                        const double b[],
                        const int ni,
                        const int nk,
                        const int init,
                        const double factor);
  
  const double s8_localcoordsofnode(const int node, const int flag,
                                    const int numnode) const;
  
  void s8_shapefunctions(vector<double>& funct, Epetra_SerialDenseMatrix& deriv,
                               const double r, const double s, const int numnode, 
                               const int doderiv) const;
  
  void s8tvbo(const double e1, const double e2,
              Epetra_SerialDenseMatrix& bop,
              const vector<double>& funct,
              const Epetra_SerialDenseMatrix& deriv,
              const int iel,
              const int numdf,
              const double akov[][3],
              const double a3kvp[][2],
              const int nsansq);
  
  void s8tvmr(const double x[][MAXNOD_SHELL8], 
              const double a3[][MAXNOD_SHELL8], 
              double akov[][3],
              double akon[][3],
              double amkov[][3],
              double amkon[][3],
              double* det,
              const vector<double>& funct,
              const Epetra_SerialDenseMatrix& deriv,
              const int iel,
              double a3kvp[][2],
              const int flag);

  void s8tmtr(const double x[][MAXNOD_SHELL8],
              const double a3[][MAXNOD_SHELL8],
              const double e3,
              double gkov[][3],
              double gkon[][3],
              double gmkov[][3],
              double gmkon[][3],
              double* det,
              const vector<double>& funct,
              const Epetra_SerialDenseMatrix& deriv,
              const int iel,
              const double condfac,
              const int flag);

  void s8tvhe(double gmkovr[][3], double gmkovc[][3], 
              double gmkonr[][3], double gmkonc[][3],
              double gkovr[][3], double gkovc[][3],
              double* detr, double* detc,
              double amkovc[][3], double amkovr[][3],
              double akovc[][3], double akovr[][3],
              double a3kvpc[][2], double a3kvpr[][2],
              const double e3, const double condfac);

  void s8_ans_colloquationpoints(const int nsansq, const int iel, const int ans,
                                 double xr1[], double xs1[], double xr2[], double xs2[],
                                 vector<double>                  funct1q[],
                                 Epetra_SerialDenseMatrix        deriv1q[],
                                 vector<double>                  funct2q[],
                                 Epetra_SerialDenseMatrix        deriv2q[],
                                 const double xrefe[][MAXNOD_SHELL8],
                                 const double a3r[][MAXNOD_SHELL8],
                                 const double xcure[][MAXNOD_SHELL8],
                                 const double a3c[][MAXNOD_SHELL8],
                                 double       akovr1q[][3][3],
                                 double       akonr1q[][3][3],
                                 double       amkovr1q[][3][3],
                                 double       amkonr1q[][3][3],
                                 double       a3kvpr1q[][3][2],
                                 double       akovc1q[][3][3],
                                 double       akonc1q[][3][3],
                                 double       amkovc1q[][3][3],
                                 double       amkonc1q[][3][3],
                                 double       a3kvpc1q[][3][2],
                                 double       akovr2q[][3][3],
                                 double       akonr2q[][3][3],
                                 double       amkovr2q[][3][3],
                                 double       amkonr2q[][3][3],
                                 double       a3kvpr2q[][3][2],
                                 double       akovc2q[][3][3],
                                 double       akonc2q[][3][3],
                                 double       amkovc2q[][3][3],
                                 double       amkonc2q[][3][3],
                                 double       a3kvpc2q[][3][2],
                                 double* detr, double* detc);

  void s8_ans_colloquationcoords(double xqr1[], double xqs1[], 
                                 double xqr2[], double xqs2[],
                                 const int iel, const int ans);

  void s8_ansqshapefunctions(double frq[], double fsq[], const double r, 
                             const double s, const int iel, const int nsansq);

  void s8ansbbarq(Epetra_SerialDenseMatrix& bop, 
                  const double frq[], const double fsq[],
                  const vector<double> funct1q[], const vector<double> funct2q[],
                  const Epetra_SerialDenseMatrix deriv1q[], const Epetra_SerialDenseMatrix deriv2q[],
                  const double akovc1q[][3][3], const double akovc2q[][3][3],
                  const double a3kvpc1q[][3][2], const double a3kvpc2q[][3][2],
                  const int& iel, const int& numdf, const int& nsansq); 

  void s8anstvheq(double gmkovr[][3], double gmkovc[][3], double gmkonr[][3], double gmkonc[][3],
                  double gkovr[][3], double gkovc[][3], double amkovc[][3], double amkovr[][3],
                  double akovc[][3], double akovr[][3], double a3kvpc[][2], double a3kvpr[][2],
                  double* detr, double* detc,
                  double amkovr1q[][3][3], double amkovc1q[][3][3],
                  double akovr1q[][3][3], double akovc1q[][3][3],
                  double a3kvpr1q[][3][2], double a3kvpc1q[][3][2],
                  double amkovr2q[][3][3], double amkovc2q[][3][3],
                  double akovr2q[][3][3], double akovc2q[][3][3],
                  double a3kvpr2q[][3][2], double a3kvpc2q[][3][2],
                  double frq[], double fsq[], const double e3, 
                  const int nansq, const int iel, const double condfac);

  void s8eas(const int nhyb, const double e1, const double e2, const int iel,
             const int* eas, Epetra_SerialDenseMatrix& P);

  void s8transeas(Epetra_SerialDenseMatrix& P, Epetra_SerialDenseMatrix& transP,
                  Epetra_SerialDenseMatrix& T, const double akovr[][3],
                  const double akonr0[][3], const double detr, const double detr0,
                  const int nhyb);

  void s8vthv(double gmkovc[][3], double gmkonc[][3], const vector<double>& epsh,
              double* detc, const double e3, const double condfac);

}; // class Shell8


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\class Shell8Register

\brief A register for shell8 element

\author gee (gee@lnm.mw.tum.de)
*/
class Shell8Register : public DRT::ElementRegister
{
public:

  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  */
  Shell8Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Element
  
  */
  Shell8Register(const DRT::Elements::Shell8Register& old);

  /*!
  \brief Deep copy this instance of Shell8Register and return pointer to the copy
  
  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed
  
  */
  DRT::Elements::Shell8Register* Clone() const;

  /*!
  \brief Return unique ParObject id
  
  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Shell8Register;}

  /*!
  \brief Pack this class so it can be communicated
  
  \ref Pack and \ref Unpack are used to communicate this element
  
  */
  virtual const char* Pack(int& size) const;

  /*!
  \brief Unpack data from a char vector into this class
  
  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual bool Unpack(const char* data);

  /*!
  \brief Destructor

  */
  virtual ~Shell8Register();

  //@}

  // @{ \name Acess methods

  /*!
  \brief Print this element
  */
  void Print(ostream& os) const;

  //@}

  // @{ \name Construction

  /*!
  \brief Call Element init routines
  
  The method performs all steps needed to init the element routines
  */
  int Initialize(DRT::Discretization& dis);

  //@}

private:





}; // class DesignElementRegister









} // namespace Elements
} // namespace DRT











#endif  // #ifndef SHELL8_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SHELL8
