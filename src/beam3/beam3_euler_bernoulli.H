/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief three dimensional nonlinear torsionless rod based on a C1 curve

\level 2

*/
/*-----------------------------------------------------------------------------------------------*/

/* 3D nonlinear Euler-Bernoulli-like beam element based on chapter 5 of the diploma thesis
 * "Development of a finite element for nonlinear beams based on the formulas of Frenet-Serret" by
 * Christoph Meier. The current formulation is only able to display axial tension and bending
 * curvature based on the curve describing the centerline of an initially (i.e. stress free)
 * straight beam. There is no shear deformation and no torsion. Refer to 'beam3k' for a full
 * Kirchhoff type beam element. To be able to use this element correctly so far structural dynamic
 * parameters need to be set to:
 *
 * LOADLIN     Yes
 *
 * since due to this beam formulation external point loads are being linearized and have an effect
 * on the stiffness matrix. For this reason a special type of Neumann conditions, namely DESIGN
 * POINT MOMENT EB CONDITIONS, are needed for this element.
 *
 * As the beam curve has to be C1 it is interpolated with hermitien polynomials of order 3.
 * Therefore each node has 6 dofs: the position vector of the node (3 dofs) and the tangent vector
 * to the curve at the node (3 dofs). If Dirichlet BC are applied one has to make sure that the last
 * three flags refer to the tangent at the node (rotational degree of freedom). The flag of the
 * tangent in beam centerline direction must be set to 0, otherwise the axial tension at the
 * boundary would be prescribed.
 */

// header file only included if not yet included!
#ifndef BEAM3_EULER_BERNOULLI_H
#define BEAM3_EULER_BERNOULLI_H


#include "lib_node.H"

#include "beam3_base.H"

// #define SIMPLECALC ;          // simplified residuum and stiffness calculation for the small
//  tension case Default: Off
#define NODALDOFS \
  2  // With this flag it can be switched between third order (2 nodal dofs) and fifth order
     // hermitian polynomials (3 nodal dofs). So the valid values of NODALDOFS are 2 and 3
// Default: 2

#define ANS_BEAM3EB  // Flag, to apply Assumed Natural Strain approach for axial tension epsilon
// Default: On
#define ANSVALUES 3  // Decide, how many ANS points are applied (2 or 3 possible)
// Default: 3

#define mygaussruleeb \
  GaussRule1D::line_4point  // define gauss rule; intrule_line_1point - intrule_line_10point
                            // is implemented
// Default: 4

// #define BEAM3EBAUTOMATICDIFF      // is Sacado used or not?
//  Default: off

// #define CONSISTENTANSBEAM3EB   //Decide wether the variational correct or the simplified ANS
//  approach should be applied Default: off

// #define INEXTENSIBLE 1.0        //apply inextensibility constraint. only possible in combination
//  with FAD and ANS_BEAM3EB and ANSVALUES=3 as well as NODALDOFS=2
//  Default: off

// #define SWITCHINEXTENSIBLEON
//  Default: off

// #define ORTHOPRESSURE 3.1415926535897938
//  Default: off

// #define BEAM3EBROTPTC          //switch from standard PTC scheme to angular velocity based PTC
// for
//  tangential DoFs Default: off #define BEAM3EBCONSTSTOCHFORCE       //Flag in order to hold
//  stochastic forces constant over the element length Default: off                  //(needed in
//  order to compare with Reissner beams)

// #define BEAM3EBDISCRETELINENEUMANN 0.0 //Flag in order to interpret line Neumann condition as
//  discrete force at given element parameter position DISCRETELINENEUMANN Default: off //for all
//  elements contained in the corresponding design line!

#if (defined(CONSISTENTANSBEAM3EB) or defined(INEXTENSIBLE) or defined(ORTHOPRESSURE)) and \
    !defined(BEAM3EBAUTOMATICDIFF)
dserror(
    "CONSISTENTANSBEAM3EB, INEXTENSIBLE and ORTHOPRESSURE only work in combination with "
    "BEAM3EBAUTOMATICDIFF!");
#endif

#if defined(INEXTENSIBLE) and !defined(ANS_BEAM3EB)
dserror(Flag INEXTENSIBLE only possible in combination with flag ANS_BEAM3EB !);
#endif

namespace DRT::ELEMENTS
{
  class Beam3ebType : public DRT::ElementType
  {
   public:
    std::string Name() const override { return "Beam3ebType"; }

    static Beam3ebType& Instance();

    DRT::ParObject* Create(const std::vector<char>& data) override;

    Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
        const int id, const int owner) override;

    Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

    int Initialize(DRT::Discretization& dis) override;

    void NodalBlockInformation(
        DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

    Teuchos::SerialDenseMatrix<int, double> ComputeNullSpace(
        DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

    void SetupElementDefinition(
        std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
        override;

   private:
    static Beam3ebType instance_;
  };

  /*!
  \brief 3D nonlinear Euler-Bernoulli-like beam element based on chapter 5 of the diploma thesis
  "Development of a finite element for nonlinear beams based on the formulas of Frenet-Serret" by
  Christoph Meier. The current formulation is only able to display axial tension and bending
  curvature based on the curve describing the centerline of an initially (i.e. stress free)
  straight beam. There is no shear deformation and no torsion. Refer to 'beam3k' for a full
  Kirchhoff type element.

  */
  class Beam3eb : public Beam3Base
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
    virtual ~Beam3eb() = default;

    /*!
    \brief Deep copy this instance of Beam3eb and return pointer to the copy

    The Clone() method is used by the virtual base class Element in cases
    where the type of the derived class is unknown and a copy-ctor is needed
  .
    */
    DRT::Element* Clone() const override;

    //! \brief Get shape type of element
    DiscretizationType Shape() const override;

    /*!
    \brief Return unique ParObject id

    Every class implementing ParObject needs a unique id defined at the
    top of parobject.H
    */
    int UniqueParObjectId() const override { return Beam3ebType::Instance().UniqueParObjectId(); }

    /*!
    \brief Pack this class so it can be communicated

    \ref Pack and \ref Unpack are used to communicate this element

    */
    void Pack(DRT::PackBuffer& data) const override;

    /*!
    \brief Unpack data from a char vector into this class

    \ref Pack and \ref Unpack are used to communicate this element

    */
    void Unpack(const std::vector<char>& data) override;

    DRT::ElementType& ElementType() const override { return Beam3ebType::Instance(); }

    //@}

    //! \brief Return number of lines to this element
    int NumLine() const override { return 1; }

    //! \brief get number of nodes used for centerline interpolation
    inline int NumCenterlineNodes() const override { return 2; }

    //! \brief find out whether given node is used for centerline interpolation
    inline bool IsCenterlineNode(const DRT::Node& node) const override { return true; }

    //! \brief Get jacobian this element
    const double& GetJacobi() const { return jacobi_; }

    //! \brief get access to the reference length
    inline double RefLength() const override { return 2.0 * jacobi_; }

    /** \brief get unit tangent vector in reference configuration at i-th node of beam element
     * (element-internal numbering)
     */
    inline void GetRefTangentAtNode(LINALG::Matrix<3, 1>& Tref_i, const int& i) const override
    {
      Tref_i = Tref()[i];
    }

    //! \brief Get maximal bending curvature occuring in this element
    const double& GetKappaMax() const { return kappa_max_; }

    //! \brief Get material cross-section deformation measures, i.e. strain resultants
    inline void GetMaterialStrainResultantsAtAllGPs(std::vector<double>& axial_strain_GPs,
        std::vector<double>& shear_strain_2_GPs, std::vector<double>& shear_strain_3_GPs,
        std::vector<double>& twist_GPs, std::vector<double>& curvature_2_GPs,
        std::vector<double>& curvature_3_GPs) const override
    {
      axial_strain_GPs = axial_strain_GP_;
      // shear deformations are zero by definition for Kirchhoff beam formulation
      shear_strain_2_GPs.clear();
      shear_strain_3_GPs.clear();

      // twist deformation cannot be represented by this reduced Kirchhoff beam formulation
      twist_GPs.clear();
      // only one curvature component due to isotropic formulation of this reduced Kirchhoff beam
      curvature_2_GPs = curvature_GP_;
      curvature_3_GPs.clear();
    }

    //! \brief Get material cross-section stress resultants
    inline void GetMaterialStressResultantsAtAllGPs(std::vector<double>& axial_force_GPs,
        std::vector<double>& shear_force_2_GPs, std::vector<double>& shear_force_3_GPs,
        std::vector<double>& torque_GPs, std::vector<double>& bending_moment_2_GPs,
        std::vector<double>& bending_moment_3_GPs) const override
    {
      axial_force_GPs = axial_force_GP_;
      // note: shear deformations are zero by definition for Kirchhoff beam formulation
      shear_force_2_GPs.clear();
      shear_force_3_GPs.clear();

      // torsion moment cannot be represented by this reduced Kirchhoff beam formulation
      torque_GPs.clear();
      // only one bending moment component due to isotropic formulation of this reduced Kirchhoff
      // beam
      bending_moment_2_GPs = bending_moment_GP_;
      bending_moment_3_GPs.clear();
    }

    //! \brief get centerline position at xi \in [-1,1] (element parameter space)
    void GetPosAtXi(LINALG::Matrix<3, 1>& pos, const double& xi,
        const std::vector<double>& disp) const override;

    //! \brief get triad at xi \in [-1,1] (element parameter space)
    void GetTriadAtXi(LINALG::Matrix<3, 3>& triad, const double& xi,
        const std::vector<double>& disp) const override;

    //! \brief Get axial strain at xi for given nodal displacements
    double GetAxialStrain(double& xi, const LINALG::Matrix<12, 1>& disp_totlag) const;

    //! get internal (elastic) energy of element
    double GetInternalEnergy() const override { return Eint_; };

    //! get kinetic energy of element
    double GetKineticEnergy() const override { return Ekin_; };

    //! \brief Get vector of Teuchos::RCPs to the lines of this element
    std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

    //! \brief Get number of degrees of freedom of a single node
    int NumDofPerNode(const DRT::Node& node) const override
    {
/*note: this is not necessarily the number of DOF assigned to this node by the discretization
 *finally, but only the number of DOF requested for this node by this element; the discretization
 *will finally assign the maximal number of DOF to this node requested by any element connected to
 *this node*/
#ifndef INEXTENSIBLE
      return 3 * NODALDOFS;
#else
      if (node.Id() == this->Nodes()[0]->Id() or node.Id() == this->Nodes()[1]->Id())
      {
        return 7;
      }
      else
      {
        return 1;
      }
#endif
    }

    //! \brief Get number of degrees of freedom per element not including nodal degrees of freedom
    int NumDofPerElement() const override { return 0; }

    //! \brief Print this element
    void Print(std::ostream& os) const override;

    //! \brief get reference triad i.e. Tref_
    std::vector<LINALG::Matrix<3, 1>> Tref() const;

    //! \brief get jacobi factor jacobi_
    double jacobi() const;

    //! \brief get Jacobi factor ds/dxi(xi) at xi \in [-1;1]
    inline double GetJacobiFacAtXi(const double& xi) const override
    {
      /* beam3eb elements are restricted to initially straight beams, for which the jacobi factor
       * is constant along the entire element length */
      return jacobi();
    }

    //! @name Construction

    /*!
    \brief Read input for this element

    This class implements a dummy of this method that prints a warning and
    returns false. A derived class would read one line from the input file and
    store all necessary information.

    */
    // virtual bool ReadElement();

    //! \brief Read input for this element
    bool ReadElement(const std::string& eletype, const std::string& distype,
        DRT::INPUT::LineDefinition* linedef) override;

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
    int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
        Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
        Epetra_SerialDenseVector& elevec3) override;

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
    int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        DRT::Condition& condition, std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
        Epetra_SerialDenseMatrix* elemat1 = NULL) override;

    /*!
    \brief Evaluate PTC addition to stiffness for free Brownian motion

    An element derived from this class uses the Evaluate method to receive commands
    and parameters from some control routine in params and evaluates a  statistical Neumann
    boundary condition as used in the problem type STATISTICAL MECHANICS

    \param params (in/out)       : ParameterList for communication between control routine and
    elements \param vector<double> mydisp : current nodal displacement \param elemat1 (out) :
    artificial damping matrix to be filled by element

    \return 0 if successful, negative otherwise
    */
    template <int nnode>
    void EvaluatePTC(Teuchos::ParameterList& params, Epetra_SerialDenseMatrix& elemat1);

    //@}


    //@}
    //! sets up from current nodal position all geometric parameters (considering current position
    //! as reference configuration)
    void SetUpReferenceGeometry(const std::vector<double>& xrefe, const bool secondinit = false);

    /* computes the number of different random numbers required in each time step for generation
     of stochastic forces                                                                    */
    int HowManyRandomNumbersINeed() const override;

    //! \brief add indices of those DOFs of a given node that are positions
    inline void PositionDofIndices(std::vector<int>& posdofs, const DRT::Node& node) const override
    {
      posdofs.push_back(0);
      posdofs.push_back(1);
      posdofs.push_back(2);
    }

    /** \brief add indices of those DOFs of a given node that are tangents (in the case of Hermite
     * interpolation)
     */
    inline void TangentDofIndices(std::vector<int>& tangdofs, const DRT::Node& node) const override
    {
      tangdofs.push_back(3);
      tangdofs.push_back(4);
      tangdofs.push_back(5);
    }

    /** \brief add indices of those DOFs of a given node that are rotation DOFs (non-additive
     * rotation vectors)
     */
    inline void RotationVecDofIndices(
        std::vector<int>& rotvecdofs, const DRT::Node& node) const override
    {
    }

    /** \brief add indices of those DOFs of a given node that are 1D rotation DOFs
     *         (planar rotations are additive, e.g. in case of relative twist DOF of beam3k with
     *          rotvec=false)
     */
    inline void Rotation1DDofIndices(
        std::vector<int>& twistdofs, const DRT::Node& node) const override
    {
    }

    /** \brief add indices of those DOFs of a given node that represent norm of tangent vector
     *         (additive, e.g. in case of beam3k with rotvec=true)
     */
    inline void TangentLengthDofIndices(
        std::vector<int>& tangnormdofs, const DRT::Node& node) const override
    {
    }

    //! \brief get element local indices of those Dofs that are used for centerline interpolation
    inline void CenterlineDofIndicesOfElement(
        std::vector<unsigned int>& centerlinedofindices) const override
    {
      centerlinedofindices.resize(12, 0);

      for (unsigned int idof = 0; idof < 12; ++idof) centerlinedofindices[idof] = idof;
    }

    /** \brief extract values for those Dofs relevant for centerline-interpolation from total
     * state vector
     */
    inline void ExtractCenterlineDofValuesFromElementStateVector(const std::vector<double>& dofvec,
        std::vector<double>& dofvec_centerline, bool add_reference_values = false) const override
    {
      if (dofvec.size() != 12)
        dserror(
            "size mismatch: expected 12 values for element state vector and got %d", dofvec.size());

      dofvec_centerline.resize(12, 0.0);
      std::copy(dofvec.begin(), dofvec.end(), dofvec_centerline.begin());

      if (add_reference_values)
      {
        for (unsigned int dim = 0; dim < 3; ++dim)
          for (unsigned int node = 0; node < 2; ++node)
          {
            dofvec_centerline[6 * node + dim] += Nodes()[node]->X()[dim];

            // Hermite interpolation: add reference values for tangent DOFs as well
            dofvec_centerline[6 * node + 3 + dim] += Tref_[node](dim);
          }
      }
    }

   private:
    //! vector storing the internal force vector
    Epetra_SerialDenseVector internalforces_;

    //! variable saving whether element has already been initialized (then isinit_ == true)
    bool isinit_;

    //! Vector holding value of Jacobi determinant jacobi for each Gauss point for
    //! underintegration
    double jacobi_;

    //! bool recognizing first element call
    bool firstcall_;

    //! kinetic energy
    double Ekin_;

    //! internal energy
    double Eint_;
    //! internal energy stemming from axial tension
    double Eint_axial_;

    //! angular momentum of the element
    LINALG::Matrix<3, 1> L_;
    //! linear momentum of the element
    LINALG::Matrix<3, 1> P_;
    //! nodal tangents of last time step (necessary for PTC scheme)
    LINALG::Matrix<3, 2> t0_;
    //! nodal tangents of current time step (necessary for PTC scheme)
    LINALG::Matrix<3, 2> t_;
    //! norm of maximal curvature occuring in this element
    double kappa_max_;
    //! norm of maximal axial tension occuring in this element
    double epsilon_max_;

    //! strain resultant values at GPs
    std::vector<double> axial_strain_GP_;
    std::vector<double> curvature_GP_;

    //! stress resultant values at GPs
    std::vector<double> axial_force_GP_;
    std::vector<double> bending_moment_GP_;

#if NODALDOFS == 3
    //! Matrix holding the derivatives of the tangents at each node in reference configuration
    std::vector<LINALG::Matrix<3, 1>> Kref_;
#endif

    //! @name Internal calculation methods

    //! calculation of nonlinear stiffness and mass matrix
    void CalcInternalAndInertiaForcesAndStiff(Teuchos::ParameterList& params,
        std::vector<double>& vel, std::vector<double>& disp, Epetra_SerialDenseMatrix* stiffmatrix,
        Epetra_SerialDenseMatrix* massmatrix, Epetra_SerialDenseVector* force);

    template <unsigned int nnode, unsigned int dofpn>
    void UpdateDispTotlag(
        const std::vector<double>& disp, LINALG::Matrix<dofpn * nnode, 1>& disp_totlag) const;

    //! lump mass matrix
    void lumpmass(Epetra_SerialDenseMatrix* emass);

    //! computes translational damping forces and stiffness
    template <unsigned int nnode, unsigned int vpernode, int ndim>
    void EvaluateTranslationalDamping(Teuchos::ParameterList& params,  //!< parameter list
        const LINALG::Matrix<ndim * vpernode * nnode, 1>& vel,
        const LINALG::Matrix<ndim * vpernode * nnode, 1>& disp_totlag,
        Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
        Epetra_SerialDenseVector* force);       //!< element internal force vector


    //! computes stochastic forces and resulting stiffness
    template <unsigned int nnode, unsigned int vpernode, unsigned int ndim,
        unsigned int randompergauss>
    void EvaluateStochasticForces(Teuchos::ParameterList& params,  //!< parameter list
        const LINALG::Matrix<ndim * vpernode * nnode, 1>& disp_totlag,
        Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
        Epetra_SerialDenseVector* force);       //!< element internal force vector


    /*! \brief Assemble stochastic and viscous forces and respective stiffness according to
     *         fluctuation dissipation
     */
    template <unsigned int nnode, unsigned int vpernode, unsigned int ndim>
    void CalcBrownianForcesAndStiff(Teuchos::ParameterList& params,
        std::vector<double>& vel,               //!< element velocity vector
        std::vector<double>& disp,              //!< element displacement vector
        Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
        Epetra_SerialDenseVector* force);       //!< element internal force vector

    // don't want = operator
    Beam3eb& operator=(const Beam3eb& old) = delete;
  };

  // << operator
  std::ostream& operator<<(std::ostream& os, const DRT::Element& ele);

}  // namespace DRT::ELEMENTS

#endif
