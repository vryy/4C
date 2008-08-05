/*!----------------------------------------------------------------------
\file so_ptet.H

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#ifdef D_SOLID3
#ifdef CCADISCRET
#ifndef SO_PTET4_H
#define SO_PTET4_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/linalg_serialdensematrix.H"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"

#include "../drt_mat/material.H"

#include "so_surface.H"
#include "so_line.H"

#define NUMNOD_PTET       4    ///< number of nodes
#define NODDOF_PTET       3    ///< number of dofs per node
#define NUMDOF_PTET      12    ///< total dofs per element
#define NUMSTR_PTET       6    ///< number of strains/stresses
#define NUMGPT_PTET       1    ///< total gauss points per element
#define NUMDIM_PTET       3    ///< number of dimensions
#define NUMCOORD_PTET     4    ///< number of shape function cooordinates (ksi1-ksi4)
#define ALPHA_PTET      0.1    ///< stabilization parameter

// forward declarations
struct _MATERIAL;

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
class StructuralVolume;
class Ptet;

//---------------------------------------------------------------------
/*!
\brief A register for ptet element

\author gee (gee@lnm.mw.tum.de)
*/

class PtetRegister : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  PtetRegister(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Soh8Register

  */
  PtetRegister(const DRT::ELEMENTS::PtetRegister& old);

  /*!
  \brief Deep copy this instance of Soh8Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::PtetRegister* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_PtetRegister;}

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
  virtual ~PtetRegister();

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
  
  /*!
  \brief Call Element pre-evaluation

  */
  virtual void PreEvaluate(DRT::Discretization& dis,
                           Teuchos::ParameterList& p,
                           RCP<LINALG::SparseOperator> systemmatrix1,
                           RCP<LINALG::SparseOperator> systemmatrix2,
                           RCP<Epetra_Vector>          systemvector1,
                           RCP<Epetra_Vector>          systemvector2,
                           RCP<Epetra_Vector>          systemvector3);

  //@}

  //! allow Ptet element to access the nodal data
  friend class DRT::ELEMENTS::Ptet;

private:

  // map of row nodes adjacent to Ptet elements
  map<int,DRT::Node*>                    noderids_;
  // map of column Ptet elements
  map<int,DRT::ELEMENTS::Ptet*>          elecids_;
  // map of nodal stresses
  map<int,vector<double> >               nodestress_;
  // map of nodal stresses
  map<int,vector<double> >               nodestrain_;

#if 1 // efficiency at the price of memory

  // vector of elements adjacent to each row node
  map<int,vector<DRT::ELEMENTS::Ptet*> > adjele_;
  // map of nodes forming a patch around a node
  map<int,map<int,DRT::Node*> >          adjnode_;
  // map of location vectors for patch around a node
  map<int,vector<int> >                  adjlm_;

#endif

  void NodalIntegration(Epetra_SerialDenseMatrix*       stiff,
                        Epetra_SerialDenseVector*       force,
                        map<int,DRT::Node*>&            nodepatch,
                        vector<DRT::ELEMENTS::Ptet*>&   adjele,
                        vector<double>*                 nodalstress,
                        vector<double>*                 nodalstrain,
                        bool                            cauchy,
                        bool                            ea);


  void SelectMaterial(RCP<MAT::Material> mat,
                      Epetra_SerialDenseVector& stress,
                      Epetra_SerialDenseMatrix& cmat,
                      double& density,
                      const Epetra_SerialDenseVector& glstrain,
                      const Epetra_SerialDenseMatrix& defgrd,
                      int gp);

  // compute determinant of a 3x3
  inline static double Det(Epetra_SerialDenseMatrix& A)
  {
   double detA = A(0,0)*A(1,1)*A(2,2) +
                 A(0,1)*A(1,2)*A(2,0) +
                 A(0,2)*A(1,0)*A(2,1) -
                 A(0,2)*A(1,1)*A(2,0) -
                 A(0,0)*A(1,2)*A(2,1) -
                 A(0,1)*A(1,0)*A(2,2);
    return detA;
  }

   // compute deviatoric stresses and tangent
   static void DevStressTangent(Epetra_SerialDenseVector& Sdev,
                                Epetra_SerialDenseMatrix& CCdev,
                                Epetra_SerialDenseMatrix& CC,
                                const Epetra_SerialDenseVector& S,
                                const Epetra_SerialDenseMatrix& C);

}; // class PtetRegister


//=======================================================================
//=======================================================================

/*!
\brief A nodal-averaged strain 4-noded tet element according to Mike Puso

*/

class Ptet : public DRT::Element
{
public:

  friend class PtetRegister;


  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  Ptet(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Ptet(const Ptet& old);

  /*!
  \brief Deep copy this instance of Solid3 and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  inline DRT::Element* Clone() const
  { return new DRT::ELEMENTS::Ptet(*this); }

  /*!
  \brief Get shape type of element
  */
  inline virtual DiscretizationType Shape() const { return tet4; }

  /*!
  \brief Return number of volumes of this element
  */
  virtual int NumVolume() const {return 1;}

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const {return 4;}

  /*!
  \brief Return number of lines of this element
  */
  virtual int NumLine() const {return 6;}

  /*!
  \brief Get vector of RefCountPtrs to the lines of this element

  */
  virtual vector<RCP<DRT::Element> > Lines();

  /*!
  \brief Get vector of RefCountPtrs to the surfaces of this element

  */
  virtual vector<RCP<DRT::Element> > Surfaces();

  /*!
  \brief Get vector of RefCountPtrs to the volumes of this element

  */
  virtual vector<RCP<DRT::Element> > Volumes();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Ptet;}

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
  virtual ~Ptet();

  /*!
  \brief Get number of degrees of freedom of a certain node
         (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  virtual inline int NumDofPerNode(const DRT::Node& node) const {return 3;}

  /*!
  \brief Get number of degrees of freedom per element
         (implements pure virtual DRT::Element)

  The element decides how many element degrees of freedom it has.
  It can redecide along the way of a simulation.

  \note Element degrees of freedom mentioned here are dofs that are visible
        at the level of the total system of equations. Purely internal
        element dofs that are condensed internally should NOT be considered.
  */
  virtual inline int NumDofPerElement() const { return 0; }

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
  virtual inline RCP<DRT::ElementRegister> ElementRegister() const
  { return rcp(new DRT::ELEMENTS::PtetRegister(Type())); }

  /*!
  \brief Read input for this element
  */
  bool ReadElement();

  /*!
  \brief Evaluate an element

  Evaluate so_tet4 element stiffness, mass, internal forces, etc.

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
  virtual int EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1);


protected:

  //! action parameters recognized by this element
  enum ActionType
  {
    none,
    calc_struct_linstiff,
    calc_struct_nlnstiff,
    calc_struct_internalforce,
    calc_struct_linstiffmass,
    calc_struct_nlnstiffmass,
    calc_struct_nlnstifflmass,
    calc_struct_stress,
    postprocess_stress,
    calc_struct_eleload,
    calc_struct_fsiload,
    calc_struct_update_istep,
    calc_struct_update_imrlike,
    calc_struct_reset_istep,  //!< reset elementwise internal variables
                              //!< during iteration to last converged state
    calc_struct_B_and_F
  };

  enum StressType
  {
    so_tet4_stress_none,
    so_tet4_stress_gpxyz,
    so_tet4_stress_gprst,
    so_tet4_stress_gp123,
    so_tet4_stress_ndxyz,
    so_tet4_stress_ndrst,
    so_tet4_stress_nd123
  };

  //! type of stress evaluation
  StressType                                               stresstype_;
  //! number of the material law
  int                                                      material_;
  //! vector of history variables for each Gauss point
  vector<RCP<MAT::Material> >                              mat_;
  //! volume of element
  double                                                   V_;
  //! derivs of shape functions wrt to material coords
  LINALG::SerialDenseMatrix                                nxyz_;
  //! element's def-gradient
  bool                                                     FisNew_;
  Epetra_SerialDenseMatrix                                 F_;
  //! ptr to this element's register class
  DRT::ELEMENTS::PtetRegister*                             myregister_;

  // internal calculation methods

  // don't want = operator
  Ptet& operator = (const Ptet& old);
  
  // init the inverse of the jacobian and its determinant
  // in the material configuration
  virtual void InitElement(DRT::ELEMENTS::PtetRegister* myregister);

  // return volume of element
  virtual inline double Volume() { return V_; }

  //! Shape functions
  inline void ShapeFunction(LINALG::SerialDenseVector& funct, 
                            const double& e1,
                            const double& e2,
                            const double& e3,
                            const double& e4)
  {
    // shape function is N_i = xsi_i
    funct(0) = e1;
    funct(1) = e2;
    funct(2) = e3;
    funct(3) = e4;
    return;
  }
  
  //! Shape function derivatives
  inline void ShapeFunctionDerivatives(LINALG::SerialDenseMatrix& deriv)
  {
    // Ni,j = 1.0 for i==j, otherwise 0.0
    const int size = deriv.M()*deriv.N();
    for (int i=0; i<size; ++i) deriv.A()[i] = 0.0;
    deriv(0,0) = 1.0;
    deriv(1,1) = 1.0;
    deriv(2,2) = 1.0;
    deriv(3,3) = 1.0;
    return;
  }

  //! standards displ. tet4 calc routine
  virtual void ptetnlnstiffmass(vector<int>&              lm,
                                vector<double>&           disp,
                                vector<double>&           residual,
                                Epetra_SerialDenseMatrix* stiffmatrix,
                                Epetra_SerialDenseMatrix* massmatrix,
                                Epetra_SerialDenseVector* force);

  //! lump mass matrix (bborn 07/08)
  void ptetlumpmass(Epetra_SerialDenseMatrix* emass);

  void DeformationGradient(vector<double>& disp);


  void so_ptet_expol(Epetra_SerialDenseMatrix& stresses,
                     Epetra_SerialDenseMatrix& nodalstresses);


  
  void SelectMaterial(Epetra_SerialDenseVector& stress,
                      Epetra_SerialDenseMatrix& cmat,
                      double& density,
                      const Epetra_SerialDenseVector& glstrain,
                      const Epetra_SerialDenseMatrix& defgrd,
                      int gp);

}; // class Ptet






} // namespace ELEMENTS
} // namespace DRT

#endif  // #ifndef SO_PTET_H_
#endif  // #ifdef  CCADISCRET
#endif  // #ifdef  D_SOLID3


