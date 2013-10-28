/*!----------------------------------------------------------------------
\file beam3eb.H

\brief three dimensional nonlinear torsionless rod based on a C1 curve

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15301
</pre>
3D nonlinear Euler-Bernoulli-like beam element based on chapter 5 of the diploma thesis "Development of a finite element for
nonlinear beams based on the formulas of Frenet-Serret" by Christoph Meier. The current formulation is only able to display axial
tension and bending curvature based on the curve describing the centerline of an initially (i.e. stress free) straight beam.
There is no shear deformation and no torsion. The expansion of the model to a full Euler Bernoulli beam (inclusive torsion) is
possible and implemented in beam3ebtor.
To be able to use this element correctly so far structural dynamic parameters need to be set to:

	LOADLIN	     Yes

since due to this beam formulation external point loads are being linearized and have an effect on the stiffness matrix. For this reason a special
type of Neumann conditions, namely DESIGN POINT MOMENT EB CONDITIONS, are needed for this element.

As the beam curve has to be C1 it is interpolated with hermitien polynomials of order 3. Therefore each node has 6 dofs: the position
vector of the node (3 dofs) and the tangent vector to the curve at the node (3 dofs). If Dirichlet BC are applied one has to make sure that
the last three flags refer to the tangent at the node (rotational defree of freedom). The flag of the tangent in beam centerline direction must
be set to 0, otherwise the axial tension at the boundary would be prescribed.

ARBITRARY PRECISSION: This element offers the method HighPrecissionCalc() which has to be uncommented in Beam3eb::Evaluate. With this method the whole
FEM calculation can be done (for simple problems) in an element intern routine which calculates with a prescribed arbitrary precision. This can
be used for debugging reasons to identify the sources of ill-conditioning and rounding errors. All intern methods with the ending ...prec belong to this
high precission calculation. The following facts have to be considered when doing a high precision calculation:

- for high precision calculations the cln libraries (cln/float.h, cln/io.h, cln/float_io.h, cln/real.h) have to be included which are standard linux
  libraries for f16 machines (and higher). An instruction manual can be found at: "http://www.ginac.de/CLN/cln.html#SEC_Top". The libraries can be
  found at "/usr/include/cln" at your f16 machine.

- The cln libraries are not compatible with the Sacado library which offers the FAD (Forward Automatic Differentiation) data type. Therefore
  automatic differentiation can not be done with arbitrary precision!

- Please make sure that the line -D CMAKE_EXE_LINKER_FLAGS:STRING="-lhdf5_hl " is replaced by -D CMAKE_EXE_LINKER_FLAGS:STRING="-lhdf5_hl -lcln"
  in your do-configure file.

- There may be some standard methods which are not compatible with arbitrary precision. The only method for which this has been recognized so far
  is the method std::fill. Therefore the following code block in linalg_fixedsizematrix.H has to be adapted:

      //this line is needed for arbitrary precision (cln) data types
      std::memset(mat,0,i*j*sizeof(value_type));
      //this line is needed for FAD data type
      //std::fill(mat,mat+i*j,0);

*----------------------------------------------------------------------*/



//header file only included if not yet included!
#ifndef BEAM3EB_H
#define BEAM3EB_H


#include "Teuchos_RCP.hpp"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "../drt_lib/drt_elementtype.H"

#include "Sacado.hpp"
typedef Sacado::Fad::DFad<double> FAD;

//#define PRECISION ;           // Flag to switch on arbitrary precision calculation. To use this functionality two
                                // slight modification in the do-configure file and in linalg_fixesizematrix.H are
                                // necessary. Read the text above for more information
#ifdef PRECISION
#include <cln/float.h>
#include <cln/io.h>
#include <cln/float_io.h>
#include <cln/real.h>
using namespace cln;
#endif

//#define SIMPLECALC ;          // simplified residuum and stiffness calculation for the small tension case
#define ScaleFactorLine 1.0;    // Scale factor for each line of the equation system for better conditioning
#define ScaleFactorColumn 1.0;  // Scale factor for each column of the stiffness matrix for better conditioning
#define PreConditioning false;  // Flag to switch preconditioning of element stiffness matrix on or off
#define NODALDOFS 2             // With this flag it can be switched between third order (2 nodal dofs) and fifth order
                                // hermitian polynomials (3 nodal dofs). So the valid values of NODALDOFS are 2 and 3
//#define ANS                     // Flag, to apply Assumed Natural Strain approach for axial tension epsilon
//#define ANSVALUES 3             //Decide, how many ANS points are applied (2 or 3 possible)

#define mygaussruleeb intrule_line_6point //define gauss rule; intrule_line_1point - intrule_line_10point is implemented

//! \brief DRT: namespace of the ccarat discretization module
namespace DRT
{

namespace ELEMENTS
{

class Beam3ebType : public DRT::ElementType
{
public:

  std::string Name() const { return "Beam3ebType"; }

  static Beam3ebType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual int Initialize(DRT::Discretization& dis);

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static Beam3ebType instance_;
};

/*!
\brief 3D nonlinear Euler-Bernoulli-like beam element based on chapter 5 of the diploma thesis "Development of a finite element for
nonlinear beams based on the formulas of Frenet-Serret" by Christoph Meier. The current formulation is only able to display axial
tension and bending curvature based on the curve describing the centerline of an initially (i.e. stress free) straight beam.
There is no shear deformation and no torsion. The expansion of the model to a full Euler Bernoulli beam (inclusive torsion) is possible.

author meier (meier@lnm.mw.tum.de)
*/
class Beam3eb : public DRT::Element
{
public:

  //! @name Friends
  friend class Beam3ebType;


  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  */
  Beam3eb(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  Beam3eb(const Beam3eb& old);

  /*!
  \brief Destructor
  */
  virtual ~Beam3eb();

  /*!
  \brief Deep copy this instance of Beam3eb and return pointer to the copy

  The Clone() method is used by the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed
.
  */
  virtual DRT::Element* Clone() const;

   /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const;


  /*!
  \brief Return unique ParObject id

  Every class implementing ParObject needs a unique id defined at the
  top of drt_parobject.H
  */
  virtual int UniqueParObjectId() const { return Beam3ebType::Instance().UniqueParObjectId(); }

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

  virtual DRT::ElementType & ElementType() const
  { return Beam3ebType::Instance(); }

  //@}

  /*!
  \brief Return number of lines to this element
  */
  virtual int NumLine() const {return 1;}


  /*!
  \brief Get vector of RCPs to the lines of this element
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();


  /*!
  \brief Get number of degrees of freedom of a single node
  */
  virtual int NumDofPerNode(const DRT::Node& node) const
  {
    /*note: this is not necessarily the number of DOF assigned to this node by the discretization finally, but only
    *the number of DOF requested for this node by this element; the discretization will finally assign the maximal
    *number of DOF to this node requested by any element connected to this node*/
      return 3*NODALDOFS;
  }

  /*!
  \brief Get number of degrees of freedom per element not including nodal degrees of freedom
  */
  virtual int NumDofPerElement() const { return 0;}

  /*!
  \brief Print this element
  */
  virtual void Print(std::ostream& os) const;

  /*!
  \brief get reference triad i.e. Tref_
  */
  std::vector<LINALG::Matrix<3,1> > Tref() const;

  /*!
  \brief get jacobi factor jacobi_
  */
  double jacobi() const;

  /*!
  \brief get moment of intertia of area Iyy_
  */
  double Iyy();

  //@}

  //! @name Construction


  /*!
  \brief Read input for this element

  This class implements a dummy of this method that prints a warning and
  returns false. A derived class would read one line from the input file and
  store all necessary information.

  */
  //virtual bool ReadElement();

  /*!
  \brief Read input for this element
  */
  bool ReadElement(const std::string& eletype,
                   const std::string& distype,
                   DRT::INPUT::LineDefinition* linedef);

  //@}

  //! @name Evaluation methods


  /*!
  \brief Evaluate an element

  An element derived from this class uses the Evaluate method to receive commands
  and parameters from some control routine in params and evaluates element matrices and
  vectors accoring to the command in params.

  \note This class implements a dummy of this method that prints a warning and
        returns false.

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elemat1 (out)      : matrix to be filled by element depending on commands
                              given in params
  \param elemat2 (out)      : matrix to be filled by element depending on commands
                              given in params
  \param elevec1 (out)      : vector to be filled by element depending on commands
                              given in params
  \param elevec2 (out)      : vector to be filled by element depending on commands
                              given in params
  \param elevec3 (out)      : vector to be filled by element depending on commands
                              given in params
  \return 0 if successful, negative otherwise
  */
  virtual int Evaluate(Teuchos::ParameterList&   params,
                       DRT::Discretization&      discretization,
                       std::vector<int>&         lm,
                       Epetra_SerialDenseMatrix& elemat1,
                       Epetra_SerialDenseMatrix& elemat2,
                       Epetra_SerialDenseVector& elevec1,
                       Epetra_SerialDenseVector& elevec2,
                       Epetra_SerialDenseVector& elevec3);

  /*!
  \brief Evaluate a Neumann boundary condition

  An element derived from this class uses the EvaluateNeumann method to receive commands
  and parameters from some control routine in params and evaluates a Neumann boundary condition
  given in condition

  \note This class implements a dummy of this method that prints a warning and
        returns false.

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : Force vector to be filled by element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(Teuchos::ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);


  //@}


private:

  //! action parameters recognized by Beam3eb
  enum ActionType
  {
    calc_none,
    calc_struct_linstiff,
    calc_struct_nlnstiff,
    calc_struct_internalforce,
    calc_struct_linstiffmass,
    calc_struct_nlnstiffmass,
    calc_struct_nlnstifflmass,  //!< internal force, its stiffness and lumped mass matrix
    calc_struct_stress,
    calc_struct_eleload,
    calc_struct_fsiload,
    calc_struct_update_istep,
    calc_struct_reset_istep,
    calc_struct_ptcstiff //calculate artificial stiffness due to PTC solution strategy
  };


  //! variable saving whether element has already been initialized (then isinit_ == true)
  bool isinit_;

  //! cross section area
  double crosssec_;

  //!moment of inertia of area around y-axis(orthognoal to rod axis in yz-plane)
  double Iyy_;

  //!moment of inertia of area around z-axis(forming cartesian system with rod- and y-axis)
  double Izz_;

  //!polar moment of inertia of area
  double Irr_;

  //! Vector holding value of Jacobi determinant jacobi for each Gauss point for underintegration
  double jacobi_;

  //! bool recognizing first element call
  bool firstcall_;

  //! Matrix holding tangents at each node in reference configuration
  std::vector<LINALG::Matrix<3,1> > Tref_;

#if NODALDOFS == 3
  //! Matrix holding the derivatives of the tangents at each node in reference configuration
  std::vector<LINALG::Matrix<3,1> > Kref_;
#endif




  // some quantities nessecary for high precission calculation
  #ifdef PRECISION
  cl_F jacobiprec_;
  cl_F crosssecprec_;
  cl_F Izzprec_;
  vector<cl_F> Xrefprec_;
  cl_F Eprec_;
  cl_F balkenradiusprec_;
  LINALG::TMatrix<cl_F,3,1> Trefprec_;
  LINALG::TMatrix<cl_F,12,12> stifftest_;
  LINALG::TMatrix<cl_F,12,1> restest_;
  #endif

  //! @name methods for initilization of the element

  //!sets up from current nodal position all geometric parameters (considering current position as reference configuration)
  void SetUpReferenceGeometry(const std::vector<double>& xrefe,const bool secondinit = false);

  //@}

  //! @name auxiliary computation methods for non-additive, large rotation variables




  //@}

  //! @name Internal calculation methods

  //!calculation of nonlinear stiffness and mass matrix
  void eb_nlnstiffmass( Teuchos::ParameterList&   params,
                        std::vector<double>&      vel,
                        std::vector<double>&      disp,
                        Epetra_SerialDenseMatrix* stiffmatrix,
                        Epetra_SerialDenseMatrix* massmatrix,
                        Epetra_SerialDenseVector* force);

  #ifdef PRECISION
  //!calculation of nonlinear stiffness matrix with high precision
  void eb_nlnstiffmassprec( LINALG::TMatrix<cl_F,12,1>* displocal,
                            LINALG::TMatrix<cl_F,12,12>* stifflocal,
                            LINALG::TMatrix<cl_F,12,1>* reslocal,
                            LINALG::TMatrix<cl_F,6,1>* xreflocal);

  //!calculation of moment Neumann BC with high precision
  void EvaluateNeumannPrec( LINALG::TMatrix<cl_F,3,1> tangentdisp,
                            LINALG::TMatrix<cl_F,3,1> mextvec,
                            LINALG::TMatrix<cl_F,3,1>* fextm,
                            LINALG::TMatrix<cl_F,3,3>* stiffextm);

  //!calculation of the entire problem with high precision
  void HighPrecissionCalc();
  #endif

  //!check of stiffness matrix via FAD
  void FADCheckStiffMatrix(std::vector<double>& disp,
                           Epetra_SerialDenseMatrix* stiffmatrix,
                           Epetra_SerialDenseVector* force);

  //!check of stiffness terms stemming from Neumann Boundary Conditions via FAD
  void FADCheckNeumann(Teuchos::ParameterList& params,
                       DRT::Discretization& discretization,
                       DRT::Condition& condition,
                       std::vector<int>& lm,
                       Epetra_SerialDenseVector& elevec1,
                       Epetra_SerialDenseMatrix* elemat1);

  //! lump mass matrix
  void lumpmass(Epetra_SerialDenseMatrix* emass);

/*  float_precision power(float_precision zahl, float_precision exponent)
  {
    if (zahl<float_precision(0.0))
      return pow(-zahl, exponent);
    else
      return pow(zahl, exponent);
  }*/


  // don't want = operator
  Beam3eb& operator = (const Beam3eb& old);


}; // class Beam3eb



// << operator
std::ostream& operator << (std::ostream& os, const DRT::Element& ele);


} // namespace Elements
} // namespace DRT


#endif  // #ifndef BEAM3EB_H
