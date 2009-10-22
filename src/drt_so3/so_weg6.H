/*!----------------------------------------------------------------------
\file so_weg6.H
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
#ifndef SO_WEG6_H_
#define SO_WEG6_H_

#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"

// prestressing object
#if defined(PRESTRESS) || defined(POSTSTRESS)
#include "prestress.H"
#endif

#include "../drt_mat/material.H"
#include "so_surface.H"
#include "so_line.H"
#include "../drt_inpar/inpar_structure.H"

#include "../drt_lib/linalg_fixedsizematrix.H"

// Several parameters which are fixed for Solid Wedge6
const int NUMNOD_WEG6 = 6;      ///< number of nodes
const int NODDOF_WEG6 = 3;      ///< number of dofs per node
const int NUMDOF_WEG6 = 18;     ///< total dofs per element
const int NUMSTR_WEG6 = 6;      ///< number of strains/stresses
const int NUMGPT_WEG6 = 6;      ///< total gauss points per element
const int NUMDIM_WEG6 = 3;      ///< number of dimensions

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
class Sow6Register;
class InvDesign;

/*!
\brief A C++ version of the 6-node wedge solid element

\author frenzel (frenzel@lnm.mw.tum.de)
*/
class So_weg6 : public DRT::Element
{
public:
  //! @name Friends
  friend class Sow6Register;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner : elements owning processor
  */
  So_weg6(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  So_weg6(const So_weg6& old);

  /*!
  \brief Deep copy this instance of Solid3 and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const;

  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const;

  /*!
  \brief Return number of volumes of this element
  */
  virtual inline int NumVolume() const {return 1;}

  /*!
  \brief Return number of surfaces of this element
  */
  virtual inline int NumSurface() const {return 5;}

  /*!
  \brief Return number of lines of this element
  */
  virtual inline int NumLine() const {return 9;}

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
  virtual inline int UniqueParObjectId() const { return ParObject_So_weg6;}

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
  virtual ~So_weg6();

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
  virtual int NumDofPerNode(const DRT::Node& node) const {return 3;}

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
  /*!
  \brief Query names of element data to be visualized using BINIO

  The element fills the provided map with key names of
  visualization data the element wants to visualize AT THE CENTER
  of the element geometry. The values is supposed to be dimension of the
  data to be visualized. It can either be 1 (scalar), 3 (vector), 6 (sym. tensor)
  or 9 (nonsym. tensor)

  Example:
  \code
    // Name of data is 'Owner', dimension is 1 (scalar value)
    names.insert(pair<string,int>("Owner",1));
    // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
    names.insert(pair<string,int>("StressesXYZ",6));
  \endcode

  \param names (out): On return, the derived class has filled names with
                      key names of data it wants to visualize and with int dimensions
                      of that data.
  */
  virtual void VisNames(map<string,int>& names);

  /*!
  \brief Query data to be visualized using BINIO of a given name

  The method is supposed to call this base method to visualize the owner of
  the element.
  If the derived method recognizes a supported data name, it shall fill it
  with corresponding data.
  If it does NOT recognizes the name, it shall do nothing.

  \warning The method must not change size of data

  \param name (in):   Name of data that is currently processed for visualization
  \param data (out):  data to be filled by element if element recognizes the name
  */
  virtual bool VisData(const string& name, vector<double>& data);

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
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);


  //@}

protected:

  //! action parameters recognized by so_hex8
  enum ActionType
  {
    none,
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
    calc_struct_update_imrlike,
    calc_struct_reset_istep,  //!< reset elementwise internal variables
                              //!< during iteration to last converged state
    postprocess_stress,
    prestress_update,
    inversedesign_update
  };

  //! kinematic types of so_weg6
  enum KinematicType
  {
    sow6_geolin,
    sow6_totlag,
    sow6_updlag
  };
  //! kinematic type
  KinematicType                                     kintype_;
  //! container containing all kind of solid3 thisandthat
  DRT::Container                                    data_;
  //! vector of inverses of the jacobian in material frame
  vector<LINALG::Matrix<NUMDIM_WEG6,NUMDIM_WEG6> >  invJ_;
  //! determinant of Jacobian in material frame
  vector<double>                                    detJ_;

#if defined(PRESTRESS) || defined(POSTSTRESS)
  /// Prestressing object
  RCP<DRT::ELEMENTS::PreStress>                    prestress_;

  // compute Jacobian mapping wrt to deformed configuration
  void UpdateJacobianMapping(const vector<double>& disp,
                             DRT::ELEMENTS::PreStress& prestress);

  // compute defgrd in all gp for given disp
  void DefGradient(const vector<double>& disp,
                   Epetra_SerialDenseMatrix& gpdefgrd,
                   DRT::ELEMENTS::PreStress& prestress);
#endif

  friend class DRT::ELEMENTS::InvDesign;
#if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
  RCP<DRT::ELEMENTS::InvDesign>                    invdesign_;
#endif



  // internal calculation methods

  // don't want = operator
  So_weg6& operator = (const So_weg6& old);

  //! init the inverse of the jacobian and its determinant in the material configuration
  virtual void InitJacobianMapping();

  //! Calculate nonlinear stiffness and mass matrix
  virtual void sow6_nlnstiffmass(vector<int>&              lm,
                                 vector<double>&           disp,
                                 vector<double>&           residual,
                                 LINALG::Matrix<NUMDOF_WEG6,NUMDOF_WEG6>* stiffmatrix,
                                 LINALG::Matrix<NUMDOF_WEG6,NUMDOF_WEG6>* massmatrix,
                                 LINALG::Matrix<NUMDOF_WEG6,1>* force,
                                 LINALG::Matrix<NUMGPT_WEG6,NUMSTR_WEG6>* elestress,
                                 LINALG::Matrix<NUMGPT_WEG6,NUMSTR_WEG6>* elestrain,
                                 ParameterList&            params,         // algorithmic parameters e.g. time
                                 const INPAR::STR::StressType             iostress,
                                 const INPAR::STR::StrainType             iostrain);

  //! material selection for sow6
  void sow6_mat_sel(Epetra_SerialDenseVector* stress,
                    Epetra_SerialDenseMatrix* cmat,
                    double* density,
                    const Epetra_SerialDenseVector* glstrain,
                    Epetra_SerialDenseMatrix* defgrd,
                    const int gp,
                    ParameterList&            params);         // algorithmic parameters e.g. time

  //! material selection for sow6
  void sow6_mat_sel(
                    LINALG::Matrix<6,1>* stress,
                    LINALG::Matrix<6,6>* cmat,
                    double* density,
                    LINALG::Matrix<6,1>* glstrain,
                    LINALG::Matrix<3,3>* defgrd,
                    const int gp,
                    ParameterList&  params);         // algorithmic parameters e.g. time

  //! Evaluate Wedge6 Shapefcts to keep them static
  const vector<LINALG::Matrix<NUMNOD_WEG6,1> > sow6_shapefcts();
  //! Evaluate Wedge6 Derivs to keep them static
  const vector<LINALG::Matrix<NUMDIM_WEG6,NUMNOD_WEG6> > sow6_derivs();
  //! Evaluate Wedge6 Weights to keep them static
  const vector<double> sow6_weights();


  //! calculate static shape functions and derivatives for sow6
  void sow6_shapederiv(LINALG::Matrix<NUMNOD_WEG6,NUMGPT_WEG6>** shapefct,
                       LINALG::Matrix<NUMDOF_WEG6,NUMNOD_WEG6>** deriv,
                       LINALG::Matrix<NUMGPT_WEG6,1>** weights);

  //! lump mass matrix (bborn 07/08)
  void sow6_lumpmass(LINALG::Matrix<NUMDOF_WEG6,NUMDOF_WEG6>* emass);

  /// Extrapolate Gauss-point values (e.g. stresses) to nodes and store them in elevecs
  void soweg6_expol
  (
      LINALG::Matrix<NUMGPT_WEG6,NUMSTR_WEG6>& stresses, ///< gp stresses
      LINALG::Matrix<NUMDOF_WEG6,1>& elevec1, ///< nodal normal stresses
      LINALG::Matrix<NUMDOF_WEG6,1>& elevec2  ///< nodal shear stresses
  );


}; // class So_weg6

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief A register for wedge6 element

\author frenzel (frenzel@lnm.mw.tum.de)
*/
class Sow6Register : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Sow6Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Soh8Register

  */
  Sow6Register(const DRT::ELEMENTS::Sow6Register& old);

  /*!
  \brief Deep copy this instance of Soh8Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::Sow6Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Sow6Register;}

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
  virtual ~Sow6Register();

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

  The method performs all steps needed to init the wedge elements,
  e.g. defines shell-thickness direction
  */
  virtual int Initialize(DRT::Discretization& dis);

  //@}


private:


}; // class Sow6Register
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================




} // namespace ELEMENTS
} // namespace DRT


#endif /*SO_WEG6_H_*/
#endif
#endif
