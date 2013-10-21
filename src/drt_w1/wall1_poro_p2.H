/*----------------------------------------------------------------------*/
/*!
 \file wall1_poro_p2.H

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#ifndef WALL1_PORO_P2_H_
#define WALL1_PORO_P2_H_

#include "wall1_poro.H"

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{

  /*!
  \brief A C++ version of a 2 dimensional solid element with modifications for porous media

  */
  template<DRT::Element::DiscretizationType distype>
  class Wall1_PoroP2 : public Wall1_Poro<distype>
  {

    typedef DRT::ELEMENTS::Wall1_Poro<distype> my;

  public:

    //@}
    //! @name Constructors and destructors and related methods

    /*!
    \brief Standard Constructor

    \param id : A unique global id
    \param owner : elements owner
    */
    Wall1_PoroP2(int id, int owner);

    /*!
    \brief Copy Constructor

    Makes a deep copy of a Element

    */
    Wall1_PoroP2(const Wall1_PoroP2& old);

    /*!
    \brief Destructor

    */
    virtual ~Wall1_PoroP2(){};

    //@}

    //! number of dofs per node
    static const int noddof_ = my::noddof_;

    //! total dofs per element
    static const int numdof_ = noddof_ * my::numnod_;

    //! @name Acess methods

    /*!
    \brief Deep copy this instance of Solid3 and return pointer to the copy

    The Clone() method is used from the virtual base class Element in cases
    where the type of the derived class is unknown and a copy-ctor is needed

    */
    DRT::Element* Clone() const ;

    /*!
    \brief Return unique ParObject id

    every class implementing ParObject needs a unique id defined at the
    top of this file.
    */
    virtual int UniqueParObjectId() const;

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

    //! @name Access methods

    /*!
    \brief Get number of degrees of freedom of a certain node
           (implements pure virtual DRT::Element)

    The element decides how many degrees of freedom its nodes must have.
    As this may vary along a simulation, the element can redecide the
    number of degrees of freedom per node along the way for each of it's nodes
    separately.
    */
    virtual int NumDofPerNode(const DRT::Node& node) const {return 2;}

    /*!
    \brief Get number of degrees of freedom of a certain node in case of multiple
           dofsets (implements pure virtual DRT::Element)

    The element decides how many degrees of freedom its nodes must have.
    As this may vary along a simulation, the element can redecide the
    number of degrees of freedom per node along the way for each of it's nodes
    separately.
    */
    virtual int NumDofPerNode(const unsigned nds, const DRT::Node& node, const std::string disname = std::string("dummy")) const
    {
      /* here we assume the dofset numbers to be assigned as follows
       *     structure = dofset 0
       *     fluid     = dofset 1
       *     scatra    = dofset 2
       * this is kind of implicit. One might think of a better solution ...
      */
      if (nds==0)
        // structure, 2D: 2+1 Dofs per node
        return NumDofPerNode(node);
      else if (nds==1)
        // in case of porous media: nds=1 (second discretization) fluid: 3Dof/Node
        return 3;
      else if (nds==2)
        // in case of scalar transport in porous media: nds=2 (third discretization) scatra: 1Dof/Node
        return 1;
      dserror("invalid number of dof sets");
      return -1;
    };

    /*!
    \brief Print this element
    */
    void Print(std::ostream& os) const  ;

    virtual DRT::ElementType & ElementType() const;

    //@}

    //! initialize the inverse of the jacobian and its determinant in the material configuration
    virtual void InitElement();

    //! return number of degrees of freedom
    virtual int NumDof(){return 2;};

    //! return whether this element has more degrees of freedom than standard solid elements
    virtual bool HasExtraDof(){return false;}

    //! evaluate constitutive relation for porosity and compute derivatives
   virtual void ConsitutiveDerivatives(Teuchos::ParameterList& params,       ///< (i) parameter list
                                                  double     press,   ///< (i) fluid pressure at gauss point
                                                  double     J,       ///< (i) Jacobian determinant at gauss point
                                                  double     porosity,///< (i) porosity at gauss point
                                                  double*    dW_dp,   ///< (o) derivative of potential w.r.t. pressure
                                                  double*    dW_dphi, ///< (o) derivative of potential w.r.t. porosity
                                                  double*    dW_dJ,   ///< (o) derivative of potential w.r.t. jacobian
                                                  double*    W        ///< (o) inner potential
                                                  );

  protected:

    //! don't want = operator
    Wall1_PoroP2& operator = (const Wall1_PoroP2& old);

    /*!
    \brief Evaluate an element

    Evaluate So3_poro element stiffness, mass, internal forces, etc.
    Templated evaluate routine of element matrixes

    If NULL on input, the controlling method does not expect the element
    to fill these matrices or vectors.

    \return 0 if successful, negative otherwise
    */
    virtual int MyEvaluate
    (
      Teuchos::ParameterList& params,  ///< ParameterList for communication between control routine and elements
      DRT::Discretization&      discretization,  ///< pointer to discretization for de-assembly
      DRT::Element::LocationArray& la, ///< location array for de-assembly
      Epetra_SerialDenseMatrix& elemat1,  ///< (stiffness-)matrix to be filled by element.
      Epetra_SerialDenseMatrix& elemat2,  ///< (mass-)matrix to be filled by element.
      Epetra_SerialDenseVector& elevec1,  ///< (internal force-)vector to be filled by element
      Epetra_SerialDenseVector& elevec2,  ///< vector to be filled by element
      Epetra_SerialDenseVector& elevec3  ///< vector to be filled by element
    );

    //! Calculate nonlinear stiffness and internal force for poroelasticity problems
    void nlnstiff_poroelast(
        std::vector<int>& lm, ///< location matrix
        LINALG::Matrix<my::numdim_, my::numnod_>&               disp,         // current displacements
        LINALG::Matrix<my::numdim_, my::numnod_>&               vel,          // current velocities
        LINALG::Matrix<my::numnod_, 1>*               porosity,
        //    std::vector<double>&      residual,       ///< current residual displ
        LINALG::Matrix<my::numdim_,my::numnod_> & evelnp, //< fluid velocity of element
        LINALG::Matrix<my::numnod_,1> & epreaf, //< fluid pressure of element
        LINALG::Matrix<numdof_,numdof_>* stiffmatrix, ///< element stiffness matrix
        LINALG::Matrix<numdof_,numdof_>* reamatrix, // element reactive matrix
        LINALG::Matrix<numdof_,1>* force, ///< element internal force vector
        Teuchos::ParameterList& params ///< algorithmic parameters e.g. time
        );

    //! Calculate coupling terms in nonlinear stiffness and internal force for poroelasticity problems
    void coupling_poroelast(
        std::vector<int>& lm, ///< location matrix
        LINALG::Matrix<my::numdim_, my::numnod_>&               disp,         // current displacements
        LINALG::Matrix<my::numdim_, my::numnod_>&               vel,          // current velocities
        LINALG::Matrix<my::numnod_, 1>*               porosity,
        LINALG::Matrix<my::numdim_, my::numnod_> & evelnp, //< fluid velocity of element
        LINALG::Matrix<my::numnod_,1> & epreaf, //< fluid pressure of element
        LINALG::Matrix<numdof_,(my::numdim_+1)*my::numnod_>* stiffmatrix, ///< element stiffness matrix
        LINALG::Matrix<numdof_,(my::numdim_+1)*my::numnod_>* reamatrix, //< element reactive matrix
        LINALG::Matrix<numdof_,1>* force, ///< element internal force vector
        Teuchos::ParameterList& params); ///< algorithmic parameters e.g. time

    //! Calculate coupling terms in nonlinear stiffness and internal force for poroelasticity problems
    void coupling_poroscatra(
        std::vector<int>& lm, ///< location matrix
        LINALG::Matrix<my::numdim_, my::numnod_>&               disp,         // current displacements
        LINALG::Matrix<my::numdim_, my::numnod_>&               vel,          // current velocities
        LINALG::Matrix<my::numnod_, 1>*               porosity,
        LINALG::Matrix<my::numdim_, my::numnod_> & evelnp, //< fluid velocity of element
        LINALG::Matrix<my::numnod_,1> & epreaf, //< fluid pressure of element
        LINALG::Matrix<numdof_,my::numnod_>* stiffmatrix, ///< element stiffness matrix
        LINALG::Matrix<numdof_,my::numnod_>* reamatrix, //< element reactive matrix
        LINALG::Matrix<numdof_,1>* force, ///< element internal force vector
        Teuchos::ParameterList& params); ///< algorithmic parameters e.g. time

    void GaussPointLoopPoroScatraOD(
        Teuchos::ParameterList&                 params,
        const LINALG::Matrix<my::numdim_,my::numnod_>&  xrefe,
        const LINALG::Matrix<my::numdim_,my::numnod_>&  xcurr,
        const LINALG::Matrix<my::numdim_,my::numnod_>&  nodaldisp,
        const LINALG::Matrix<my::numdim_,my::numnod_>&  nodalvel,
        const LINALG::Matrix<my::numdim_,my::numnod_>&  evelnp,
        const LINALG::Matrix<my::numnod_,1> &       epreaf,
        const LINALG::Matrix<my::numnod_, 1>*       porosity_dof,
        LINALG::Matrix<numdof_, my::numnod_>& stiffmatrix
                                            );

    //compute porosity at gausspoint and linearization of porosity w.r.t. structural displacements
    virtual void ComputePorosityAndLinearization(Teuchos::ParameterList&                params,
                                                 const double&                          press,
                                                 const double&                          J,
                                                 const int&                             gp,
                                                 const LINALG::Matrix<my::numnod_,1>&   shapfct,
                                                 const LINALG::Matrix<my::numnod_,1>*   myporosity,
                                                 const LINALG::Matrix<1,my::numdof_>&   dJ_dus,
                                                 double &                               porosity,
                                                 LINALG::Matrix<1,my::numdof_>&         dphi_dus);

    //compute porosity at gausspoint and linearization of porosity w.r.t. fluid pressure
    virtual void ComputePorosityAndLinearizationOD(Teuchos::ParameterList&              params,
                                                 const double&                          press,
                                                 const double&                          J,
                                                 const int&                             gp,
                                                 const LINALG::Matrix<my::numnod_,1>&    shapfct,
                                                 const LINALG::Matrix<my::numnod_,1>*    myporosity,
                                                 double &                               porosity,
                                                 double &                               dphi_dp);

  };
}
}


#endif /* WALL1_PORO_P2_H_ */
