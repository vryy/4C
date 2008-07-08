/*!----------------------------------------------------------------------
\file so_disp.H
\brief pure displacement based element for higher order shape functions (quadratic)

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET
#ifndef SO_DISP_H
#define SO_DISP_H

#include "so_integrator.H"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_integration.H"

#include "../drt_mat/material.H"

// Several parameters which are fixed for Solid Wedge6
const int NUMDIM_DISP = 3;      // number of dimensions
const int NODDOF_DISP = 3;      // number of dofs per node
const int NUMSTR_DISP = 6;      // number of strains/stresses

// forward declarations
struct _MATERIAL;


namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
class SoDispVolume;
class SoDispRegister;

/*!
\brief A general quadratic displacement based solid element

\author frenzel (frenzel@lnm.mw.tum.de)
*/
class SoDisp : public DRT::Element
{
public:
  //! @name Friends
  friend class SoDispRegister;
  friend class Integrator_SoDisp;
  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  SoDisp(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  SoDisp(const SoDisp& old);

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
  virtual int NumVolume() const {return 1;}

  /*!
  \brief Return number of lines of this element
  */
  virtual int NumLine() const {return DRT::UTILS::getNumberOfElementLines(Shape());}

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const {return DRT::UTILS::getNumberOfElementSurfaces(Shape());}

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
  virtual int UniqueParObjectId() const { return ParObject_SoDisp;}

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
  virtual ~SoDisp();

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
  virtual bool ReadElement();

  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate so_hex8 element stiffness, mass, internal forces, etc.

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
    calc_init_vol
  };

  //! kinematic types of so_disp
  enum KinematicType
  {
    sodisp_geolin,
    sodisp_totlag,
    sodisp_updlag
  };
  //! stress evaluation types of so_disp
  enum StressType
  {
    sodisp_stress_none,
    sodisp_stress_gpxyz,
    sodisp_stress_gprst,
    sodisp_stress_gp123,
    sodisp_stress_ndxyz,
    sodisp_stress_ndrst,
    sodisp_stress_nd123
  };
  //! kinematic type
  KinematicType                                   kintype_;
  //! type of stress evaluation
  StressType                                      stresstype_;
  //! container containing all kind of solid3 thisandthat
  DRT::Container                                  data_;
  //! Gaussrule
  DRT::UTILS::GaussRule3D                         gaussrule_;

  int numnod_disp_;      // number of nodes
  int numdof_disp_;     // total dofs per element
  int numgpt_disp_;      // total gauss points per element


  // internal calculation methods

  // don't want = operator
  SoDisp& operator = (const SoDisp& old);

  virtual void sodisp_nlnstiffmass(vector<int>&              lm,
                                 vector<double>&           disp,
                                 vector<double>&           residual,
                                 Epetra_SerialDenseMatrix* stiffmatrix,
                                 Epetra_SerialDenseMatrix* massmatrix,
                                 Epetra_SerialDenseVector* force,
                                 ParameterList&            params);         // algorithmic parameters e.g. time

  void sodisp_mat_sel(Epetra_SerialDenseVector* stress,
                    Epetra_SerialDenseMatrix* cmat,
                    double* density,
                    const Epetra_SerialDenseVector* glstrain,
                    ParameterList&            params);         // algorithmic parameters e.g. time

  /*!
  \brief Create a integrator for the SoDisp element

  this constructor provides Object Oriented So_integrator functionality by adapting the strictly procedural
  DRT::UTILS::IntegrationPoints3D struct.
  */
  class Integrator_SoDisp : public DRT::ELEMENTS::So_integrator
  {
  	public:
  		Integrator_SoDisp(DRT::ELEMENTS::SoDisp& this_element);
  };

}; // class SoDisp

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief A register for solid3 element

\author frenzel (frenzel@lnm.mw.tum.de)
*/
class SoDispRegister : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  SoDispRegister(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Soh8Register

  */
  SoDispRegister(const DRT::ELEMENTS::SoDispRegister& old);

  /*!
  \brief Deep copy this instance of Soh8Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::SoDispRegister* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_SoDispRegister;}

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
  virtual ~SoDispRegister();

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

  The method performs all steps needed to init the element. it e.g. calculates the
  field of averaged directors at nodes.
  */
  int Initialize(DRT::Discretization& dis);


private:


}; // class SoDispRegister


} // namespace ELEMENTS
} // namespace DRT


#endif /*SO_DISP_H*/
#endif
#endif
