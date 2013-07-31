/*----------------------------------------------------------------------*/
/*!
 \file wall1_poro.H

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#ifndef WALL1_PORO_H_
#define WALL1_PORO_H_

#include "../drt_inpar/inpar_structure.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"

#include "wall1.H"
#include "wall1_poro_eletypes.H"

#include "../drt_so3/so_poro_interface.H"

namespace MAT
{
  class StructPoro;
  class FluidPoro;
}

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
  class Wall1_Poro : public Wall1, public So_Poro_Interface
  {

  public:

    //@}
    //! @name Constructors and destructors and related methods

    /*!
    \brief Standard Constructor

    \param id : A unique global id
    \param owner : elements owner
    */
    Wall1_Poro(int id, int owner);

    /*!
    \brief Copy Constructor

    Makes a deep copy of a Element

    */
    Wall1_Poro(const Wall1_Poro& old);

    /*!
    \brief Destructor

    */
    virtual ~Wall1_Poro(){};

    //@}

    //! number of element nodes (
    static const int numnod_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    //! number of strains per node
    static const int numstr_ = 3;

    //! number or second derivatives
    static const int numderiv2_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;;

    //! number of degrees of freedom of element
    static const int numdof_=numnod_*noddof_;

    //!total gauss points per element
    int numgpt_;

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
    virtual int UniqueParObjectId() const
    {
      switch(distype)
      {
      case DRT::Element::quad4:
        return DRT::ELEMENTS::WallQuad4PoroType::Instance().UniqueParObjectId();
        break;
      case DRT::Element::quad9:
        return DRT::ELEMENTS::WallQuad9PoroType::Instance().UniqueParObjectId();
        break;
      default:
        dserror("unknown element type");
        break;
      }
      return -1;
    };

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
    \brief Get number of degrees of freedom of a certain node in case of multiple
           dofsets (implements pure virtual DRT::Element)

    The element decides how many degrees of freedom its nodes must have.
    As this may vary along a simulation, the element can redecide the
    number of degrees of freedom per node along the way for each of it's nodes
    separately.
    */
    virtual int NumDofPerNode(const unsigned nds, const DRT::Node& node) const
    {
      /* here we assume the dofset numbers to be assigned as follows
       *     structure = dofset 0
       *     fluid     = dofset 1
       *     scatra    = dofset 2
       * this is kind of implicit. One might think of a better solution ...
      */
      if (nds==0)
        // structure, 2D: 2 Dofs per node
        return Wall1::NumDofPerNode(node);
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

    virtual DRT::ElementType & ElementType() const
    {
      switch(distype)
      {
      case DRT::Element::quad4:
        return DRT::ELEMENTS::WallQuad4PoroType::Instance();
        break;
      case DRT::Element::quad9:
        return DRT::ELEMENTS::WallQuad9PoroType::Instance();
        break;
      default:
        dserror("unknown element type");
        break;
      }
      return DRT::ELEMENTS::WallQuad4PoroType::Instance();
    };

    //@}

    //! @name Evaluation

    void PreEvaluate
      (
        Teuchos::ParameterList& params,  ///< ParameterList for communication between control routine and elements
        DRT::Discretization&      discretization,  ///< pointer to discretization for de-assembly
        DRT::Element::LocationArray& la ///< location array for de-assembly
      );

    /*!
    \brief Evaluate an element

    Evaluate So_tet4fbar element stiffness, mass, internal forces, etc.

    If NULL on input, the controlling method does not expect the element
    to fill these matrices or vectors.

    \return 0 if successful, negative otherwise
    */
    int Evaluate
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


    //! initialize the inverse of the jacobian and its determinant in the material configuration
    virtual void InitElement();

    //@}

    //! @name Input and Creation

    /*!
    \brief Read input for this element
    */
    bool ReadElement(const std::string& eletype,
                     const std::string& eledistype,
                     DRT::INPUT::LineDefinition* linedef) ;

    //@}

    virtual void ComputePorosity( Teuchos::ParameterList& params,
                               double press,
                               double J,
                               int gp,
                               double& porosity,
                               double* dphi_dp,
                               double* dphi_dJ,
                               double* dphi_dJdp,
                               double* dphi_dJJ,
                               double* dphi_dpp,
                               bool save);

    virtual void ComputeSurfPorosity( Teuchos::ParameterList& params,
                           double press,
                           double J,
                           int surfnum,
                           int gp,
                           double& porosity,
                           double* dphi_dp,
                           double* dphi_dJ,
                           double* dphi_dJdp,
                           double* dphi_dJJ,
                           double* dphi_dpp,
                           bool save) ;

    virtual double RefPorosityTimeDeriv();

    virtual int NumDof(){return 2;};

    virtual bool HasExtraDof(){return false;};

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
      names.insert(std::pair<string,int>("Owner",1));
      // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
      names.insert(std::pair<string,int>("StressesXYZ",6));
    \endcode

    \param names (out): On return, the derived class has filled names with
                        key names of data it wants to visualize and with int dimensions
                        of that data.
    */
    virtual void VisNames(std::map<std::string,int>& names);

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
    virtual bool VisData(const std::string& name, std::vector<double>& data);

  protected:

    //! don't want = operator
    Wall1_Poro& operator = (const Wall1_Poro& old);

    /*!
    \brief Evaluate an element

    Evaluate So3_poro element stiffness, mass, internal forces, etc.
    Templated evaluate routine of element matrixes

    If NULL on input, the controlling method does not expect the element
    to fill these matrices or vectors.

    \return 0 if successful, negative otherwise
    */
    int MyEvaluate
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

    //! action parameters recognized by poro element
    enum ActionType
    {
      none,
      calc_struct_nlnstiff,
      calc_struct_internalforce,
      calc_struct_nlnstiffmass,
      calc_struct_stress,
      calc_struct_multidofsetcoupling,//!< structure-fluid coupling: internal force, stiffness for poroelasticity (structural part)
    };

    //! vector of inverses of the jacobian in material frame
    std::vector<LINALG::Matrix<numdim_,numdim_> >                      invJ_;
    //! determinant of Jacobian in material frame
    std::vector<double>                                                        detJ_;
    //! container containing all kind of solid3 thisandthat
    DRT::Container                                                        data_;
    //! vector of coordinates of current integration point in reference coordinates
    std::vector<LINALG::Matrix<numdim_,1> > xsi_;

    //! Calculate nonlinear stiffness and internal force for poroelasticity problems
    virtual void nlnstiff_poroelast(
        std::vector<int>& lm, ///< location matrix
        LINALG::Matrix<numdim_, numnod_>&               disp,         // current displacements
        LINALG::Matrix<numdim_, numnod_>&               vel,          // current velocities
        LINALG::Matrix<numdim_,numnod_> & evelnp, //< fluid velocity of element
        LINALG::Matrix<numnod_,1> & epreaf, //< fluid pressure of element
        LINALG::Matrix<numdof_,numdof_>* stiffmatrix, ///< element stiffness matrix
        LINALG::Matrix<numdof_,numdof_>* reamatrix, // element reactive matrix
        LINALG::Matrix<numdof_,1>* force, ///< element internal force vector
        Teuchos::ParameterList& params ///< algorithmic parameters e.g. time
        );

    //! Calculate coupling terms in nonlinear stiffness and internal force for poroelasticity problems
    virtual void coupling_poroelast(
        std::vector<int>& lm, ///< location matrix
        LINALG::Matrix<numdim_, numnod_>&               disp,         // current displacements
        LINALG::Matrix<numdim_, numnod_>&               vel,          // current velocities
        LINALG::Matrix<numdim_,numnod_> & evelnp, //< fluid velocity of element
        LINALG::Matrix<numnod_,1> & epreaf, //< fluid pressure of element
        LINALG::Matrix<numdof_,(numdim_+1)*numnod_>* stiffmatrix, ///< element stiffness matrix
        LINALG::Matrix<numdof_,(numdim_+1)*numnod_>* reamatrix, //< element reactive matrix
        LINALG::Matrix<numdof_,1>* force, ///< element internal force vector
        Teuchos::ParameterList& params); ///< algorithmic parameters e.g. time

    //! Calculate coupling stress for poroelasticity problems
    virtual void couplstress_poroelast(
        LINALG::Matrix<numdim_, numnod_>&               disp,         // current displacements
        LINALG::Matrix<numdim_,numnod_> & evelnp, //< fluid velocity of element
        LINALG::Matrix<numnod_,1> & epreaf, //< fluid pressure of element
        Epetra_SerialDenseMatrix*      elestress, ///< stresses at GP
        Epetra_SerialDenseMatrix* elestrain, ///< strains at GP
        Teuchos::ParameterList& params, ///< algorithmic parameters e.g. time
        const INPAR::STR::StressType iostress ///< stress output option
        );

    void GaussPointLoop(   Teuchos::ParameterList& params,
                                        const LINALG::Matrix<numdim_,numnod_>&  xrefe,
                                        const LINALG::Matrix<numdim_,numnod_>&  xcurr,
                                        const LINALG::Matrix<numdim_,numnod_>&  nodaldisp,
                                        const LINALG::Matrix<numdim_,numnod_>&  nodalvel,
                                        const LINALG::Matrix<numdim_,numnod_> & evelnp,
                                        const LINALG::Matrix<numnod_,1> &       epreaf,
                                        const LINALG::Matrix<numnod_, 1>*       porosity_dof,
                                        LINALG::Matrix<numdof_,numdof_>& erea_v,
                                        LINALG::Matrix<numdof_, numdof_>*  stiffmatrix,
                                        LINALG::Matrix<numdof_,numdof_>* reamatrix,
                                        LINALG::Matrix<numdof_,1>*              force
                                            );

    void GaussPointLoopOD(
        Teuchos::ParameterList&                 params,
        const LINALG::Matrix<numdim_,numnod_>&  xrefe,
        const LINALG::Matrix<numdim_,numnod_>&  xcurr,
        const LINALG::Matrix<numdim_,numnod_>&  nodaldisp,
        const LINALG::Matrix<numdim_,numnod_>&  nodalvel,
        const LINALG::Matrix<numdim_,numnod_>&  evelnp,
        const LINALG::Matrix<numnod_,1> &       epreaf,
        LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>& stiffmatrix
                                            );

    //compute porosity at gausspoint and linearization of porosity w.r.t. structural displacements
    virtual void ComputePorosityAndLinearization(Teuchos::ParameterList&                params,
                                                 const double&                          press,
                                                 const double&                          J,
                                                 const int&                             gp,
                                                 const LINALG::Matrix<numnod_,1>&       shapfct,
                                                 const LINALG::Matrix<numnod_,1>*       myporosity,
                                                 const LINALG::Matrix<1,numdof_>&       dJ_dus,
                                                 double &                               porosity,
                                                 LINALG::Matrix<1,numdof_>&             dphi_dus);

    //compute porosity at gausspoint and linearization of porosity w.r.t. fluid pressure
    virtual void ComputePorosityAndLinearizationOD(Teuchos::ParameterList&        params,
                                                 const double&                    press,
                                                 const double&                    J,
                                                 const int&                       gp,
                                                 const LINALG::Matrix<numnod_,1>& shapfct,
                                                 const LINALG::Matrix<numnod_,1>* myporosity,
                                                 double &                         porosity,
                                                 double &                         dphi_dp);

    void FillMatrixAndVectors(
        const int &                               gp,
        const LINALG::Matrix<numnod_,1>&          shapefct,
        const LINALG::Matrix<numdim_,numnod_>&    N_XYZ,
        const double&                             J,
        const double&                             press,
        const double&                             porosity,
        const LINALG::Matrix<numdim_,1>&          velint,
        const LINALG::Matrix<numdim_,1>&          fvelint,
        const LINALG::Matrix<numdim_,numdim_>&    fvelder,
        const LINALG::Matrix<numdim_,numdim_>&    defgrd_inv,
        const LINALG::Matrix<numstr_,numdof_>&    bop,
        const LINALG::Matrix<numdim_,numdim_>&    C_inv,
        const LINALG::Matrix<numdim_,1>&          Finvgradp,
        const LINALG::Matrix<1,numdof_>&          dphi_dus,
        const LINALG::Matrix<1,numdof_>&          dJ_dus,
        const LINALG::Matrix<numstr_,numdof_>&    dCinv_dus,
        const LINALG::Matrix<numdim_,numdof_>&    dFinvdus_gradp,
        LINALG::Matrix<numdof_,numdof_>&          erea_v,
        LINALG::Matrix<numdof_, numdof_>*         stiffmatrix,
        LINALG::Matrix<numdof_,1>*                force,
        LINALG::Matrix<numstr_,1>&                fstress);

    void FillMatrixAndVectorsBrinkman(
        const int &                                     gp,
        const double&                                   J,
        const double&                                   porosity,
        const LINALG::Matrix<numdim_,numdim_>&          fvelder,
        const LINALG::Matrix<numdim_,numdim_>&          defgrd_inv,
        const LINALG::Matrix<numstr_,numdof_>&          bop,
        const LINALG::Matrix<numdim_,numdim_>&          C_inv,
        const LINALG::Matrix<1,numdof_>&                dphi_dus,
        const LINALG::Matrix<1,numdof_>&                dJ_dus,
        const LINALG::Matrix<numstr_,numdof_>&          dCinv_dus,
        const LINALG::Matrix<numdim_*numdim_,numdof_>&  dFinvTdus,
        LINALG::Matrix<numdof_, numdof_>*               stiffmatrix,
        LINALG::Matrix<numdof_,1>*                      force,
        LINALG::Matrix<numstr_,1>&                      fstress);

    void FillMatrixAndVectorsOD(
        const int &                             gp,
        const LINALG::Matrix<numnod_,1>&        shapefct,
        const LINALG::Matrix<numdim_,numnod_>&  N_XYZ,
        const double&                           J,
        const double&                           porosity,
        const double&                           dphi_dp,
        const LINALG::Matrix<numdim_,1>&        velint,
        const LINALG::Matrix<numdim_,1>&        fvelint,
        const LINALG::Matrix<numdim_,numdim_>&  defgrd_inv,
        const LINALG::Matrix<numdim_,1>&        Gradp,
        const LINALG::Matrix<numstr_,numdof_>&  bop,
        const LINALG::Matrix<numdim_,numdim_>&  C_inv,
        LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>& ecoupl);

    void FillMatrixAndVectorsBrinkmanOD(
        const int &                             gp,
        const LINALG::Matrix<numnod_,1>&        shapefct,
        const LINALG::Matrix<numdim_,numnod_>&  N_XYZ,
        const double&                           J,
        const double&                           porosity,
        const double&                           dphi_dp,
        const LINALG::Matrix<numdim_,numdim_>&  fvelder,
        const LINALG::Matrix<numdim_,numdim_>&  defgrd_inv,
        const LINALG::Matrix<numstr_,numdof_>&  bop,
        const LINALG::Matrix<numdim_,numdim_>&  C_inv,
        LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>& ecoupl);

    void ComputeBOperator(
        LINALG::Matrix<numstr_,numdof_>& bop,
        const LINALG::Matrix<numdim_,numdim_>& defgrd,
        const LINALG::Matrix<numdim_,numnod_>& N_XYZ);

    void ComputeShapeFunctionsAndDerivatives(
        const int &                       gp,
        LINALG::Matrix<numnod_,1>&        shapefct,
        LINALG::Matrix<numdim_,numnod_>&  deriv ,
        LINALG::Matrix<numdim_,numnod_>&  N_XYZ
        );

    double ComputeJacobianDeterminant(
       const int & gp,
       const LINALG::Matrix<numdim_,numnod_>&   xcurr,
       const   LINALG::Matrix<numdim_,numnod_>& deriv
       );

    void ComputeLinearizationOfJacobian(
       LINALG::Matrix<1,numdof_>& dJ_dus,
       const double& J,
       const LINALG::Matrix<numdim_,numnod_>& N_XYZ,
       const LINALG::Matrix<numdim_,numdim_>& defgrd_inv
       );

    void ComputeAuxiliaryValues(
        const LINALG::Matrix<numdim_,numnod_>&    N_XYZ,
        const LINALG::Matrix<numdim_,numdim_>&    defgrd_inv,
        const LINALG::Matrix<numdim_,numdim_>&    C_inv,
        const LINALG::Matrix<numdim_,1>&          Gradp,
        LINALG::Matrix<numdim_*numdim_,numdof_>&  dFinvTdus,
        LINALG::Matrix<numdim_,1>&                Finvgradp,
        LINALG::Matrix<numdim_,numdof_>&          dFinvdus_gradp,
        LINALG::Matrix<numstr_,numdof_>&          dCinv_dus);

    //! push forward of material stresses to the current, spatial configuration
    void PK2toCauchy(
      LINALG::Matrix<Wall1::numstr_,1>& stress,
      LINALG::Matrix<numdim_,numdim_>& defgrd,
      LINALG::Matrix<numdim_,numdim_>& cauchystress
      );

    void ExtractValuesFromGlobalVector(
                                        const DRT::Discretization&          discretization, ///< discretization
                                        const int&                          dofset,         ///< number of dofset
                                        const std::vector<int>&             lm,             ///< location vector
                                        LINALG::Matrix<numdim_,numnod_> *   matrixtofill,   ///< vector field
                                        LINALG::Matrix<numnod_,1> *         vectortofill,   ///< scalar field
                                        const std::string                   state           ///< state of the global vector
                                        );

    void GetMaterials();

    //!Gauss integration rule
    DRT::UTILS::GaussIntegration intpoints_;

    //!flag indicating higher order element
    bool ishigherorder_;

    //! flag indicating initialization of element
    bool init_;

    //! flag for scatra coupling
    bool scatracoupling_;

    //!corresponding fluid material
    Teuchos::RCP< MAT::FluidPoro > fluidmat_;

    //!own poro structure material
    Teuchos::RCP< MAT::StructPoro > structmat_;

  }; // class So3_Poro

}
}
#endif /* WALL1_PORO_H_ */
