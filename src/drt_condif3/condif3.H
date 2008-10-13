/*!
\file condif3.H

\brief A C++ wrapper for the fluid3 element

This file contains the element specific service routines like
Pack, Unpack, NumDofPerNode etc.


<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
#ifdef D_FLUID3
#ifdef CCADISCRET
#ifndef DRT_CONDIF3_H
#define DRT_CONDIF3_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_integration.H"

// forward declarations
struct _MATERIAL;

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
class Condif3Surface;
class Condif3Register;

/*!
\brief A C++ wrapper for the condif3 element

\author g.bau (bauer@lnm.mw.tum.de)
*/
class Condif3 : public DRT::Element
{
public:

  //! @name Friends
  friend class Condif3Register;
  friend class Condif3Surface;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor
  */
  Condif3(
    int id,    ///< A unique global id
    int owner  ///< processor id who owns a certain instance of this class
    );

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Condif3(const Condif3& old);

  /*!
  \brief Deep copy this instance of Condif3 and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const;

  /// Set element material
  /*!
    Material numbers are read from the input file. The element stores
    a corresponding material object. These material objects can be
    anything from very simple (just a little calculation) to highly
    sophisticated with history data. The material is packed and
    unpacked along with its element.

    \param matnum : material number from input file

    \note reimplementation of this method, due to initialising
          numdofpernode_, since the material is known now.
   */
  virtual void SetMaterial(int matnum);

  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const;

  /*!
  \brief Return number of lines of this element
  */
  virtual int NumLine() const
  {if (NumNode()==8 || NumNode()==20 || NumNode()==27) return 12;
   else if (NumNode()==4 || NumNode()==10) return 6;
   else {
    dserror("Could not determine number of lines");
    return -1;
    }
  }

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const
  {if (NumNode()==8 || NumNode()==20 || NumNode()==27) return 6;
    else if (NumNode()==4 || NumNode()==10) return 4;
   else if (NumNode()==6 || NumNode()==15 || NumNode()==5) return 5;
   else {
   	dserror("Could not determine number of surfaces");
   	return -1;
	}
  }

  /*!
  \brief Return number of volumes of this element (always 1)
  */
  virtual int NumVolume() const {return 1;}

  /*!
  \brief Get vector of RefCountPtrs to the lines of this element
  */
  virtual vector<RCP<DRT::Element> > Lines();

  /*!
  \brief Get vector of RefCountPtrs to the surfaces of this element
  */
  virtual vector<RCP<DRT::Element> > Surfaces();

  /*!
  \brief Get vector of RefCountPtrs to the volume of this element
  */
  virtual vector<RCP<DRT::Element> > Volumes();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Condif3;}

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
  virtual ~Condif3();

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
  virtual int NumDofPerNode(const DRT::Node& node) const
  {
    dsassert(numdofpernode_!=-1,"NumDofPerNode is -1");
    return numdofpernode_;
  }

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
  \brief Evaluate an element, i.e. call the implementation to evaluate condif3
  element stiffness, mass, internal forces or evaluate errors, statistics or updates
  etc. directly.

  Following implementations of the element are allowed:
  <pre>
    o Evaluation of condif system matrix and residual for the One-Step-Theta or BDF2
      method

    o Evaluation of condif system matrix and residual for the stationary fluid solver
  </pre>

  \param params (in/out): ParameterList for communication between control routine
                          and elements
  \param discretization (in): A reference to the underlying discretization
  \param lm (in):         location vector of this element
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

  //@}

private:

  //! action parameters recognized by condif3
  enum ActionType
  {
    none,
    calc_condif_systemmat_and_residual,
    calc_condif_flux,
    initialize_one_step_theta,
    calc_temp_and_dens
  };

  //! flux parameters recognized by condif3
  enum FluxType
  {
    noflux,
    totalflux,
    diffusiveflux
  };

  //! Gaussrule
  DRT::UTILS::GaussRule3D   gaussrule_;
  //! container containing all kind of condif3 thisandthat
  DRT::Container            data_;
  //! number of dofs per node (for systems of transport equations)
  //! (storage neccessary because we dont know the material in the post filters anymore)
  int                       numdofpernode_;

  // internal calculation methods

  //! don't want = operator
  Condif3& operator = (const Condif3& old);

  //! convert action string to condif3 action type
  ActionType convertStringToActionType(const string& action) const;

  //! calculate flux vector for actual solution
  Epetra_SerialDenseMatrix CalculateFlux(
      vector<double>&           phinp,
      struct _MATERIAL*         material,
      bool                      temperature,
      Epetra_SerialDenseVector& evel,
      Condif3::FluxType         fluxtype,
      const int&                dofindex
      );

  //! calculate integrals of temperature, density and domain
 void CalculateTempAndDens(
      ParameterList&            params,
      vector<double>&           ephinp,
      vector<double>&           edensnp
      );

}; // class Condif3


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief A register for condif3 element

\author gee (gee@lnm.mw.tum.de)
*/
class Condif3Register : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Condif3Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Condif3Register

  */
  Condif3Register(const DRT::ELEMENTS::Condif3Register& old);

  /*!
  \brief Deep copy this instance of Condif3Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::Condif3Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Condif3Register;}

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
  virtual ~Condif3Register();

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

  The method performs all steps needed to init the condif3 elements. it e.g. calculates the
  field of averaged directors at nodes.
  */
  //int Initialize(DRT::Discretization& dis);

  //@}

private:


}; // class Condif3Register


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief An element representing a surface of a condif3 element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate surface Neumann boundary conditions that might be
      adjacent to a parent condif3 element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

\author gee (gee@lnm.mw.tum.de)
*/
class Condif3Surface : public DRT::Element
{
public:

  friend class Condif3Line;

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner: Processor owning this surface
  \param nnode: Number of nodes attached to this element
  \param nodeids: global ids of nodes attached to this element
  \param nodes: the discretizations map of nodes to build ptrs to nodes from
  \param parent: The parent fluid element of this surface
  \param lsurface: the local surface number of this surface w.r.t. the parent element
  */
  Condif3Surface(int id, int owner, int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::ELEMENTS::Condif3* parent, const int lsurface);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Condif3Surface(const Condif3Surface& old);

  /*!
  \brief Deep copy this instance of an element and return pointer to the copy

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
  {if (NumNode()==4 || NumNode()==8 || NumNode()==9) return 4;
  else if (NumNode()==3 || NumNode()==6) return 3;
  else
  {
    dserror("Could not determine number of lines");
    return -1;
  }
  }

  /*!
  \brief Get vector of RefCountPtrs to the lines of this element

  */
  virtual vector<RCP<DRT::Element> > Lines();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of the drt_parobject.H file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Condif3Surface;}

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
  virtual ~Condif3Surface();

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
  virtual int NumDofPerNode(const DRT::Node& node) const
  {return parent_->numdofpernode_;}

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

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate condif3 element stiffness, mass, internal forces etc

  \param params (in/out): ParameterList for communication between control routine
                          and elements
  \param discretization (in): A reference to the underlying discretization
  \param lm (in):         location vector of this element
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

  //! @name Evaluate methods

  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a surface Neumann condition on the condif3 element

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

  //! action parameters recognized by Condif3Surface
  enum ActionType
  {
    none,
    calc_condif_flux,
    calc_elch_electrode_kinetics
  };

  // don't want = operator
  Condif3Surface& operator = (const Condif3Surface& old);

  //! evaluate an electrode kinetics boundary condition
  void EvaluateElectrodeKinetics(
     Epetra_SerialDenseMatrix& emat,  ///< element-matrix
     Epetra_SerialDenseVector& erhs,  ///< element-rhs
     const vector<double>&   ephinp,  ///< actual conc. and pot. values
     struct _MATERIAL*     material,  ///< the material
     const double&             sign,  ///< sign due to electrode type
     const int&          reactantid,  ///< MATERIAL id of reactive species
     const std::string*    kinetics,  ///< desired electrode kinetics model
     const double&             pot0,  ///< actual electrode potential on metal side
     const double&           alphaa,  ///< anodic factor in Butler-Volmer model
     const double&           alphac,  ///< cathodic factor in Butler-Volmer model
     const double&               i0,  ///< exchange current density
     const double&              frt   ///< factor F/RT
     );

 /*!
  brief integral of normal flux over boundary surface

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param xyze (in)          : node coordinates
  \param enormflux (in)     : normal fluxes at element nodes
  */
  void NormalFluxIntegral(ParameterList&                  params,
                          DRT::Discretization&            discretization,
                          const vector<int>&              lm,
                          const Epetra_SerialDenseMatrix  xyze,
                          const vector<double>&           enormflux);

  //! The parent element of this surface
  DRT::ELEMENTS::Condif3* parent_;
  //! The local surface number of this surface w.r.t to the parent_ element
  int                    lsurface_;
  //! number of dofs per node (for systems of transport equations)
  //! (storage neccessary because we dont know the material in the post filters anymore)
  int                       numdofpernode_;


}; // class Condif3Surface

// put class Condif3Line in here, if needed

} // namespace ELEMENTS
} // namespace DRT



#endif  // #ifndef CONDIF3_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
