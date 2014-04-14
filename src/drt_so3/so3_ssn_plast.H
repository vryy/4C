/*----------------------------------------------------------------------*/
/*!
\file so3_ssn_plast.H

<pre>
   Maintainer: Alexander Seitz
               seitz@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15271
</pre>
*/


/*----------------------------------------------------------------------*
 | definitions                                              seitz 07/13 |
 *----------------------------------------------------------------------*/
#ifndef SO3_SSN_PLAST_H
#define SO3_SSN_PLAST_H

/*----------------------------------------------------------------------*
 | headers                                                  seitz 07/13 |
 *----------------------------------------------------------------------*/
#include "so3_ssn_plast_eletypes.H"

#include "../drt_tsi/kinematictype.H"
#include "../drt_thermo/thermo_ele_impl_utils.H"

// include header because number of Gauss points are determined
// dependent on distype
#include "../drt_fem_general/drt_utils_gausspoints.H"

#define AS_CONVERGENCE_TOL 1.e-12

/*----------------------------------------------------------------------*
 |                                                          seitz 07/13 |
 *----------------------------------------------------------------------*/
namespace DRT
{
class Discretization;

namespace ELEMENTS
{
  //! A wrapper for 3-dimensional structural elements
  class So3_Base
  {

  public:

    //! @name Constructors and destructors and related methods

    //! Destructor
    virtual ~So3_Base(){};

    //@}



    //! kinematic type
    GenKinematicType kintype_;

    // get the kinematic type from the corresponding cloned element
    // function is called in clone strategy in tsi_utils.cpp
    GenKinematicType GetKinematicType()
    {
      return kintype_;
    }
  };  // class So3_Base


  //! A C++ version of a 3 dimensional solid element with modifications von Mises
  //! plasticity using a semi-smooth Newton method
  //!
  //! A structural 3 dimensional solid displacement element for large deformations
  //! and with small and large strains
  template<class so3_ele, DRT::Element::DiscretizationType distype>
  class So3_Plast : public so3_ele,
                     public So3_Base
  {
    //! @name Friends
    friend class So_hex8PlastType;
    friend class So_hex8fbarPlastType;
    friend class So_tet4PlastType;
    friend class So_hex27PlastType;

  public:

    //@}
    //! @name Constructors and destructors and related methods


    //! Standard Constructor
    So3_Plast(
      int id,  //!< (i) this element's global id
      int owner  //!< elements owner
      );

    //! Copy Constructor
    //! Makes a deep copy of a Element
    So3_Plast(const So3_Plast& old);

    //! Destructor
    virtual ~So3_Plast(){};

    //@}

    //! number of element nodes
    static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
    //! number of space dimensions
    static const int nsd_ = 3;
    //! number of dofs per node
    static const int numdofpernode_ = 3;
    //! total dofs per element
    static const int numdofperelement_ = numdofpernode_ * nen_;
    //! number of strains/stresses
    static const int numstr_ = 6;
    //! number of Gauss points per element (value is added in so3_thermo.cpp)
    int numgpt_;
    //! static const is required for fixedsizematrices
    //! TODO maybe more beauty is possible
    static const int numgpt_post = THR::DisTypeToSTRNumGaussPoints<distype>::nquad;


    //! @name Acess methods

    //! Deep copy this instance of Solid3 and return pointer to the copy
    //!
    //! The Clone() method is used from the virtual base class Element in cases
    //! where the type of the derived class is unknown and a copy-ctor is needed
    DRT::Element* Clone() const;

    //! Return unique ParObject id
    //!
    //! every class implementing ParObject needs a unique id defined at the top of
    //! this file.
    virtual int UniqueParObjectId() const;


    //! Pack this class so it can be communicated
    //! Pack and \ref Unpack are used to communicate this element
    virtual void Pack(DRT::PackBuffer& data) const;

    //! Unpack data from a char vector into this class
    //! Pack and \ref Unpack are used to communicate this element
    virtual void Unpack(const std::vector<char>& data);

    //! Get number of degrees of freedom of a certain node in case of multiple
    //! dofsets (implements pure virtual DRT::Element)
    //!
    //! The element decides how many degrees of freedom its nodes must have.
    virtual int NumDofPerNode() const
    {
        return nsd_;
    };

    //@}

    //! @name Access methods

    //! Print this element
    void Print(std::ostream& os) const  ;

    //! return elementtype
    virtual DRT::ElementType& ElementType() const;

    //@}

    //! @name Input and Creation

    //! Query names of element data to be visualized using BINIO
    //!
    //! The element fills the provided map with key names of
    //! visualization data the element wants to visualize AT THE CENTER
    //! of the element geometry. The values is supposed to be dimension of the
    //! data to be visualized. It can either be 1 (scalar), 3 (vector), 6 (sym. tensor)
    //! or 9 (nonsym. tensor)
    //!
    //! Example:
    //! \code
    //!   // Name of data is 'Owner', dimension is 1 (scalar value)
    //!   names.insert(std::pair<string,int>("Owner",1));
    //!   // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
    //!   names.insert(std::pair<string,int>("StressesXYZ",6));
    //! \endcode
    //!
    //!  names (out): On return, the derived class has filled names with key
    //!               names of data it wants to visualize and with int dimensions
    //!               of that data.
    virtual void VisNames(std::map<std::string,int>& names);

    //!  Query data to be visualized using BINIO of a given name
    //!
    //! The method is supposed to call this base method to visualize the owner of
    //! the element.
    //! If the derived method recognizes a supported data name, it shall fill it
    //! with corresponding data.
    //! If it does NOT recognizes the name, it shall do nothing.
    //!
    //! warning: the method must not change size of data
    //!
    //!  name (in):   Name of data that is currently processed for visualization
    //! \param data (out):  data to be filled by element if element recognizes the name
    virtual bool VisData(const std::string& name, std::vector<double>& data);

    //! read input for this element
    bool ReadElement(
      const std::string& eletype,  //!< so3plast(fbar)
      const std::string& eledistype,  //!< hex8,tet4,...
      DRT::INPUT::LineDefinition* linedef  //!< what parameters have to be read
      );

    //@}

    //! @name Evaluation

    //! evaluate an element
    //! evaluate element stiffness, mass, internal forces, etc.
    //!
    //! if NULL on input, the controlling method does not expect the element
    //!  to fill these matrices or vectors.
    //!
    //!  \return 0 if successful, negative otherwise
    int Evaluate(
        Teuchos::ParameterList& params,//!< ParameterList for communication between control routine and elements
        DRT::Discretization& discretization, //!< pointer to discretization for de-assembly
        DRT::Element::LocationArray& la,  //!< location array for de-assembly
        Epetra_SerialDenseMatrix& elemat1_epetra, //!< (stiffness-)matrix to be filled by element.
        Epetra_SerialDenseMatrix& elemat2_epetra,//!< (mass-)matrix to be filled by element.
        Epetra_SerialDenseVector& elevec1_epetra,//!< (internal force-)vector to be filled by element
        Epetra_SerialDenseVector& elevec2_epetra, //!< vector to be filled by element
        Epetra_SerialDenseVector& elevec3_epetra//!< vector to be filled by element
      );

    //@}

    //! init the inverse of the jacobian and its determinant in the material
    //! configuration
    virtual void InitJacobianMapping();

    // get parameter list from ssn_plast_manager
    virtual void ReadParameterList(Teuchos::RCP<Teuchos::ParameterList> plparams);

    //@}

  protected:

    //! don't want = operator
    So3_Plast& operator = (const So3_Plast& old);

    //! action parameters recognized by so_hex8
    enum ActionType
    {
      none,
      calc_struct_linstiff,  //!< geometrical linear analysis: stiffness matrix
      calc_struct_nlnstiff,  //!< stiffness matrix
      calc_struct_internalforce,  //!< internal force
      calc_struct_linstiffmass,  //!< geometrical linear analysis: internal force,
                                 //!< its stiffness and mass matrix
      calc_struct_nlnstiffmass,  //!< internal force, its stiffness and mass matrix
      calc_struct_nlnstifflmass,  //!< internal force, its stiffness and lumped mass matrix
      calc_struct_stress,  //!< calculate stresses
      calc_struct_update_istep,  //!< update all at element level
      calc_struct_reset_istep,  //!< reset elementwise internal variables
                                //!< during iteration to last converged state
                                //!< needed for predictor TangDis
      postprocess_stress  //!< postprocess stresses
    };

    /*!
     * \brief EAS technology enhancement types of so_hex8
     *
     * Solid Hex8 has EAS enhancement of GL-strains to avoid locking.
     * New enum for so3_plast. Currently only for hex8
     */
    enum EASType  // with meaningfull value for matrix size info
    {
      soh8p_easnone = 0,  //!< no EAS i.e. displacement based with tremendous locking
      soh8p_easmild = 9,  //!< 9 parameters consisting of modes to alleviate
                         //!< shear locking (bending) and main incompressibility modes
                         //!< (for Solid Hex8)
      soh8p_easfull = 21,  //!< 21 parameters to prevent almost all locking modes.
                          //!< Equivalent to all 30 parameters to fully complete element
                          //!< with quadratic modes (see Andelfinger 1993 for details) and therefore
                          //!< also suitable for distorted elements.
                          //!< (for Solid Hex8)
    };

    //! vector of inverses of the jacobian in material frame
    std::vector<LINALG::Matrix<nsd_,nsd_> > invJ_;
    //! determinant of Jacobian in material frame
    std::vector<double> detJ_;
    //! container containing all kind of solid3 thisandthat
    DRT::Container data_;
    //! vector of coordinates of current integration point in reference coordinates
    std::vector<LINALG::Matrix<nsd_,1> > xsi_;

    DRT::UTILS::GaussIntegration intpoints_;

    //! @name plasticity related stuff


    //! Calculate nonlinear stiffness and mass matrix with condensed plastic matrices
    virtual void nln_stiffmass(
        std::vector<int>&              lm,             // location matrix
        std::vector<double>&           disp,           // current displacements
        std::vector<double>&           residual,       // current residual displ
        LINALG::Matrix<numdofperelement_,numdofperelement_>* stiffmatrix, // element stiffness matrix
        LINALG::Matrix<numdofperelement_,numdofperelement_>* massmatrix,  // element mass matrix
        LINALG::Matrix<numdofperelement_,1>* force,                 // element internal force vector
        LINALG::Matrix<numgpt_post,numstr_>* elestress,   // stresses at GP
        LINALG::Matrix<numgpt_post,numstr_>* elestrain,   // strains at GP
        Teuchos::ParameterList&        params,         // algorithmic parameters e.g. time
        const INPAR::STR::StressType   iostress,  // stress output option
        const INPAR::STR::StrainType   iostrain,  // strain output option
        const int MyPID  // processor id
        );

    //! Calculate nonlinear stiffness and mass matrix with condensed plastic matrices
    virtual void nln_stiffmass_hill(
        std::vector<int>&              lm,             // location matrix
        std::vector<double>&           disp,           // current displacements
        std::vector<double>&           residual,       // current residual displ
        LINALG::Matrix<numdofperelement_,numdofperelement_>* stiffmatrix, // element stiffness matrix
        LINALG::Matrix<numdofperelement_,numdofperelement_>* massmatrix,  // element mass matrix
        LINALG::Matrix<numdofperelement_,1>* force,                 // element internal force vector
        LINALG::Matrix<numgpt_post,numstr_>* elestress,   // stresses at GP
        LINALG::Matrix<numgpt_post,numstr_>* elestrain,   // strains at GP
        Teuchos::ParameterList&        params,         // algorithmic parameters e.g. time
        const INPAR::STR::StressType   iostress,  // stress output option
        const INPAR::STR::StrainType   iostrain,  // strain output option
        const int MyPID  // processor id
        );

    //! Calculate linear stiffness and mass matrix
    virtual void lin_stiffmass(
        std::vector<int>&              lm,             // location matrix
        std::vector<double>&           disp,           // current displacements
        std::vector<double>&           residual,       // current residual displ
        LINALG::Matrix<numdofperelement_,numdofperelement_>* stiffmatrix, // element stiffness matrix
        LINALG::Matrix<numdofperelement_,numdofperelement_>* massmatrix,  // element mass matrix
        LINALG::Matrix<numdofperelement_,1>* force,                 // element internal force vector
        LINALG::Matrix<numgpt_post,numstr_>* elestress,   // stresses at GP
        LINALG::Matrix<numgpt_post,numstr_>* elestrain,   // strains at GP
        Teuchos::ParameterList&        params,         // algorithmic parameters e.g. time
        const INPAR::STR::StressType   iostress,  // stress output option
        const INPAR::STR::StrainType   iostrain,  // strain output option
        const int MyPID  // processor id
        );

    //! calculate the exponential of a 3x3 matrix (symmetric or non-symmetric)
    void MatrixExponential3x3( LINALG::Matrix<3,3>& MatrixInOut );

    //! calculate the derivative of the exponential of a symmetric 3x3 matrix
    void MatrixExponentialDerivativeSym3x3(const LINALG::Matrix<3,3> MatrixIn, LINALG::Matrix<6,6>& MatrixExpDeriv);

    //! calculate the derivative of the exponential of a non-symmetric 3x3 matrix
    void MatrixExponentialDerivative3x3(const LINALG::Matrix<3,3> MatrixIn, LINALG::Matrix<9,9>& MatrixExpDeriv);

    // index notation:
    // Dexp_ijkl += 0.5*fac* (A_ik*B_lj + B_ik*A_lj + A_jk*B_li + B_jk*A_li)
    // A and B symmetric!
    void AddToSymMatrixExponentialDeriv(const double fac,
        const LINALG::Matrix<3,3> A,const LINALG::Matrix<3,3> B, LINALG::Matrix<6,6>& Dexp);

    // index notation:
    // Dexp_ijkl += fac*( A_ik * B_lj )
    // A and B non-symmetric!
    void AddToMatrixExponentialDeriv(const double fac,
        const LINALG::Matrix<3,3> A,const LINALG::Matrix<3,3> B, LINALG::Matrix<9,9>& Dexp);

    // Add plastic increment of converged state to plastic history for nonlinear kinematics
    void UpdatePlasticDeformation_nln();

    // Add plastic increment of converged state to plastic history for nonlinear kinematics and Hill plasticity
    void UpdatePlasticDeformationHill_nln();

    // Add plastic increment of converged state to plastic history for nonlinear kinematics
    void UpdatePlasticDeformation_lin();

    //@}

    //! calculate nonlinear B-operator
    void CalculateBop(
      LINALG::Matrix<numstr_,numdofperelement_>* bop,
      LINALG::Matrix<nsd_,nsd_>* defgrd,
      LINALG::Matrix<nsd_,nen_>* N_XYZ
      );

    //! calculate linear B-operator
    void CalculateBoplin(
      LINALG::Matrix<numstr_,numdofperelement_>* boplin,
      LINALG::Matrix<nsd_,nen_>* N_XYZ
      );

    // return the accumulated plastic strain at given gp
    double AccumulatedStrain(int gp) {return sqrt(2./3.)*last_alpha_isotropic_->at(gp)(0);}

    bool HaveHillPlasticity() ;

    /// Extrapolate Gauss-point values (e.g. stresses) to nodes and store results in elevectors
    /// todo: Unfortunately, there is no universal extrapolation function in all base elements,
    /// i.e. the hex8 elements would call soh8_expol, hex27 would call soh27_expol...
    /// For now, we copy the extrapolation function here and allow the use only for hex8 elements
    void soh8_expol(
        LINALG::Matrix<numgpt_post,numstr_>& stresses,  ///< gp stresses
        Epetra_MultiVector& expolstresses);  ///< nodal stresses

    // EAS element techonolgy *******************************************
    //! Initialize data for EAS (once)
    void EasInit();

    //! setup EAS for each evaluation
    void EasSetup(
        std::vector<Epetra_SerialDenseMatrix>** M_GP,    // M-matrix evaluated at GPs
        double& detJ0,                      // det of Jacobian at origin
        LINALG::Matrix<numstr_,numstr_>& T0invT,   // maps M(origin) local to global
        const LINALG::Matrix<nen_,nsd_>& xrefe)    // material element coords
    ;

    /*! \brief Calculate the deformation gradient that is consistent
     *         with modified (e.g. EAS) GL strain tensor.
     *         Expensive (two polar decomposition), but required, if
     *         the material evaluation is based on the deformation
     *         gradient rather than the GL strain tensor (e.g. plasticity).
     *
     * \param defgrd_disp  (in)  : displacement-based deformation gradient
     * \param glstrain_mod (in)  : modified GL strain tensor (strain-like Voigt notation)
     * \param defgrd_mod   (out) : consistent modified deformation gradient
     */
    void CalcConsistentDefgrd(LINALG::Matrix<3,3> defgrd_disp,
                              LINALG::Matrix<6,1> glstrain_mod,
                              LINALG::Matrix<3,3>& defgrd_mod);

    //! @name plasticity related stuff

    //@}

  private:
    DRT::Node** Nodes();

    Teuchos::RCP<MAT::Material> Material() const;

    int Id() const;

    // algoirthmic parameters
    double stab_s_;
    double cpl_;

    static const int VOIGT3X3SYM_[3][3];
    static const int VOIGT3X3NONSYM_[3][3];

    // plastic history variables ******************************************
    // inverse plastic deformation gradient for each Gauss point at last converged state
    Teuchos::RCP<std::vector<LINALG::Matrix<3,3> > > last_plastic_defgrd_inverse_;

    // accumulated plastic strain for each Gauss point at last converged state
    Teuchos::RCP<std::vector<LINALG::Matrix<1,1> > > last_alpha_isotropic_;

    // accumulated plastic strain for each Gauss point at last converged state
    Teuchos::RCP<std::vector<LINALG::Matrix<3,3> > > last_alpha_kinematic_;

    // classification, if the Gauss point is currently in the active (true) or inactive (false) set
    Teuchos::RCP<std::vector<bool> > activity_state_;

    // Hill plasticity ********************************************
    // here we have 5 scalar plastic variables at each time step
    // (symmetric, traceless 3x3 matrix)
    // Kbb^-1 at each Gauss point for recovery of inner variables
    Teuchos::RCP<std::vector<LINALG::Matrix<5,5> > > KbbInv_;

    // Kbd at each Gauss point for recovery of inner variables
    Teuchos::RCP<std::vector<LINALG::Matrix<5,numdofperelement_> > > Kbd_;

    // f_b at each Gauss point for recovery of inner variables
    Teuchos::RCP<std::vector<LINALG::Matrix<5,1> > > fbeta_;

    // plastic increment at each Gauss point at last Newton iteration
    // for a slightly simpler implementation, we use Delta alpha^k=-Delta D^p instead of Delta D^p
    Teuchos::RCP<std::vector<LINALG::Matrix<5,1> > > DalphaK_last_iter_;

    // Hill plasticity ********************************************
    // here we have 8 scalar plastic variables at each time step
    // (non-symmetric, traceless 3x3 matrix)
    // Kbb^-1 at each Gauss point for recovery of inner variables
    Teuchos::RCP<std::vector<LINALG::Matrix<8,8> > > KbbInvHill_;

    // Kbd at each Gauss point for recovery of inner variables
    Teuchos::RCP<std::vector<LINALG::Matrix<8,numdofperelement_> > > KbdHill_;

    // f_b at each Gauss point for recovery of inner variables
    Teuchos::RCP<std::vector<LINALG::Matrix<8,1> > > fbetaHill_;

    // plastic increment at each Gauss point at last Newton iteration
    Teuchos::RCP<std::vector<LINALG::Matrix<8,1> > > mDLp_last_iter_;

    // accumulated plastic strain increment for each Gauss point
    Teuchos::RCP<std::vector<double> > deltaAlphaI_;

    // EAS element technology ************************************
    Teuchos::RCP<Epetra_SerialDenseMatrix> KaaInv_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> Kad_;
    Teuchos::RCP<Epetra_SerialDenseVector> feas_;
    Teuchos::RCP<std::vector<Epetra_SerialDenseMatrix> > Kba_;
    Teuchos::RCP<Epetra_SerialDenseVector> alpha_eas_;
    Teuchos::RCP<Epetra_SerialDenseVector> alpha_eas_last_timestep_;
    Teuchos::RCP<Epetra_SerialDenseVector> alpha_eas_delta_over_last_timestep_;
    EASType                  eastype_;
    int                      neas_;

  };  // class So3_Plast

}  // namespace ELEMENTS
}  // namespace DRT


/*----------------------------------------------------------------------*/
#endif  // SO3_SSN_PLAST_H
