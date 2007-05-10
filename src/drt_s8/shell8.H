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

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"


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

namespace Elements
{
// forward declarations
class Shell8Line;
class Shell8Register;

/*!
\class Shell8

\brief A C++ wrapper for the shell8 element

\author gee (gee@lnm.mw.tum.de)
*/
class Shell8 : public DRT::Element
{
public:

  //! @name Friends
  friend class Shell8Register;
  friend class Shell8Line;

  //@}
  //! @name Constructors and destructors and related methods

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
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const;

  /*!
  \brief Return number of lines of this element
  */
  virtual int NumLine() const
  {if (NumNode()==4 || NumNode()==9) return 4; else return 3;}

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const {return 1;}

  /*!
  \brief Get vector of ptrs to the lines of this element

  */
  virtual DRT::Element** Lines();

  /*!
  \brief Get vector of ptrs to the surfaces of this element

  */
  virtual DRT::Element** Surfaces();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Shell8;}

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
  virtual ~Shell8();

  //@}

  //! @name Acess methods


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

  //! @name Input and Creation

  /*!
  \brief Read input for this element
  */
  bool ReadElement();

  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate shell8 element stiffness, mass, internal forces etc

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


  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a surfaces Neumann condition on the shell element

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1);


  //@}

  //! @name Other

  /*!
  \brief Get ptr to nodal thicknesses
  */
  vector<double>* GetThickness() { return data_.Get<vector<double> >("thick");}

  /*!
  \brief Get reference to directors
  */
  Epetra_SerialDenseMatrix* GetDirectors() { return data_.Get<Epetra_SerialDenseMatrix>("a3ref"); }

  //@}


private:

  //! action parameters recognized by shell8
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
    calc_struct_fsiload,
    calc_struct_update_istep
  };

  //! type of forces to be calculated in postprocessing
  enum ForceType
  {
    s8_none,
    s8_xyz,
    s8_rst,
    s8_rst_ortho
  };
  //! type of force calculation
  ForceType                                       forcetype_;
  //! shell thickness
  double                                          thickness_;
  //! no. gaussian points in each spatial dimension
  int                                             ngp_[3];
  //! no. gaussian points for triangle elements
  int                                             ngptri_;
  //! no. of internal eas dofs
  int                                             nhyb_;
  //! eas parameters
  int                                             eas_[5];
  //! type of ANS: 0=none 1=Q 2=T 3=QT
  int                                             ans_;
  //! scaled director conditioning
  double                                          sdc_;
  //! number of the material law
  int                                             material_;
  //! container containing all kind of shell8 thisandthat
  DRT::Container                                  data_;
  //! vector of surfaces of this element (length 1)
  vector<DRT::Element*>                           surfaces_;
  //! vector of lines adjacent to this element (shell8: 3 or 4 lines)
  vector<RefCountPtr<DRT::Elements::Shell8Line> > lines_;
  //! vector of lines adjacent to this element (shell8: 3 or 4 lines)
  vector<DRT::Elements::Shell8Line*>              lineptrs_;

  // internal calculation methods

  // don't want = operator
  Shell8& operator = (const Shell8& old);

  void s8_nlnstiffmass(vector<int>&              lm,
                       vector<double>&           disp,
                       vector<double>&           residual,
                       Epetra_SerialDenseMatrix* stiffmatrix,
                       Epetra_SerialDenseMatrix* massmatrix,
                       Epetra_SerialDenseVector* force,
                       struct _MATERIAL*         material);

  void s8stress(struct _MATERIAL* material, vector<double>& mydisp);

  void s8_surfaceload(ParameterList& params,Epetra_SerialDenseVector& elevec1,
                    const vector<double>& disp);

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

  void s8_jaco(const vector<double>& funct,
               const Epetra_SerialDenseMatrix& deriv,
               const double x[][MAXNOD_SHELL8],
               double xjm[][3],
               const vector<double>& hte,
               const double a3ref[][MAXNOD_SHELL8],
               const double e3,
               const int iel,
               double* det,
               double* deta);

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

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Shell8Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Shell8Register

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
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Unpack(const vector<char>& data);

  /*!
  \brief Destructor

  */
  virtual ~Shell8Register();

  //@}

  //! @name Acess methods

  /*!
  \brief Print this element
  */
  void Print(ostream& os) const;

  //@}

  //! @name Construction

  /*!
  \brief Call Element init routines

  The method performs all steps needed to init the shell8 elements. it e.g. calculates the
  field of averaged directors at nodes.
  */
  int Initialize(DRT::Discretization& dis);

  //@}

private:


}; // class Shell8Register


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\class Shell8Line

\brief An element representing a line edge of a shell8 element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate line Neumann boundary conditions that might be
      adjacent to a parent shell8 element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

\author gee (gee@lnm.mw.tum.de)
*/
class Shell8Line : public DRT::Element
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner: Processor owning this line
  \param nnode: Number of nodes attached to this element
  \param nodeids: global ids of nodes attached to this element
  \param nodes: the discretizations map of nodes to build ptrs to nodes from
  \param parent: The parent shell element of this line
  \param lline: the local line number of this line w.r.t. the parent element
  */
  Shell8Line(int id, int owner, int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::Elements::Shell8* parent, const int lline);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Shell8Line(const Shell8Line& old);

  /*!
  \brief Deep copy this instance of an element and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of the drt_parobject.H file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Shell8Line;}

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
  virtual ~Shell8Line();

  //@}

  //! @name Acess methods

  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const;

  /*!
  \brief Get number of degrees of freedom of a certain node
         (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  virtual int NumDofPerNode(const DRT::Node& node) const
  {return 6;}

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
  \brief Allocate and return a ElementRegister class

  Every implementation of DRT::Element must also implement
  an DRT::ElementRegister. This method allocates and returns
  an instance of a matching DRT::ElementRegister.
  */
  virtual RefCountPtr<DRT::ElementRegister> ElementRegister() const { return null; }

  //@}

  //! @name Evaluate methods

  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a line Neumann condition on the shell element

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1);

  //@}

private:

  // don't want = operator
  Shell8Line& operator = (const Shell8Line& old);

  //! The parent element of this line
  DRT::Elements::Shell8* parent_;
  //! The local line number of this line w.r.t to the parent_ element
  int                    lline_;

}; // class Shell8Line





} // namespace Elements
} // namespace DRT











#endif  // #ifndef SHELL8_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SHELL8
