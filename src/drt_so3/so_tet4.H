/*----------------------------------------------------------------------*/
/*! \file
\brief Solid Tet4 Element
\level 3
*----------------------------------------------------------------------*/
#ifndef SO_TET4_H
#define SO_TET4_H


#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"

#include "../drt_mat/material.H"

#include "../drt_inpar/inpar_structure.H"

#include "so_base.H"


#define NUMNOD_SOTET4 4       ///< number of nodes
#define NODDOF_SOTET4 3       ///< number of dofs per node
#define NUMDOF_SOTET4 12      ///< total dofs per element
#define NUMGPT_SOTET4 1       ///< total gauss points per element  /****/
#define NUMDIM_SOTET4 3       ///< number of dimensions/****/
#define NUMCOORD_SOTET4 4     ///< number of shape function cooordinates (ksi1-ksi4)
#define NUMNOD_SOTET4_FACE 3  ///< number of nodes on a TET4 face (which is a TRI3)
#define NUMGPT_SOTET4_FACE 1  ///< number of GP on a TET4 face (which is a TRI3)

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    // forward declarations
    class InvDesign;
    class PreStress;

    class So_tet4Type : public DRT::ElementType
    {
     public:
      std::string Name() const { return "So_tet4Type"; }

      static So_tet4Type& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

      virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
          const double* x0, int numdf, int dimns);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_tet4Type instance_;
    };

    /*!
    \brief A C++ version of the 4-node tet solid element

    */

    class So_tet4 : public So_base
    {
     public:
      //! @name Friends
      friend class So_tet4Type;


      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owning processor
      */
      So_tet4(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      So_tet4(const So_tet4& old);

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
      virtual int NumVolume() const { return 1; }

      /*!
      \brief Return number of surfaces of this element
      */
      virtual int NumSurface() const { return 4; }

      /*!
      \brief Return number of lines of this element
      */
      virtual int NumLine() const { return 6; }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element

      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Lines();

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element

      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Surfaces();

      /*!
      \brief Get vector of Teuchos::RCPs to the volumes of this element

      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Volumes();

      virtual std::vector<double> ElementCenterRefeCoords();

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const { return So_tet4Type::Instance().UniqueParObjectId(); }

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

      /*!
      \brief Destructor

      */
      virtual ~So_tet4();
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
      virtual int NumDofPerNode(const DRT::Node& node) const { return 3; }

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
      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const { return So_tet4Type::Instance(); }

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
        names.insert(std::pair<std::string,int>("Owner",1));
        // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
        names.insert(std::pair<std::string,int>("StressesXYZ",6));
      \endcode

      \param names (out): On return, the derived class has filled names with
                          key names of data it wants to visualize and with int dimensions
                          of that data.
      */
      virtual void VisNames(std::map<std::string, int>& names);

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

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);


      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate so_tet4 element stiffness, mass, internal forces, etc.

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
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);


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
      virtual int EvaluateNeumann(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1 = NULL);

      /*!
      \brief Return value how expensive it is to evaluate this element

      \param double (out): cost to evaluate this element
      */
      virtual double EvaluationCost()
      {
        if (Material()->MaterialType() == INPAR::MAT::m_struct_multiscale)
          return 25000.0;
        else
          return 10.0;
      }

      /**
       * \brief Evaluate Cauchy stress at given point in parameter space and calculate
       * linearizations
       *
       * \param xi[in]          position in parameter space xi
       * \param disp[in]        vector of displacements
       * \param n[in]           vector n (\f[\bold{n}\f])
       * \param t[in]           vector t (\f[\bold{t}\f])
       * \param sigma_nt[out]   cauchy stress tensor contracted using the vectors n and t
                                (\f[ \bold{\sigma} \cdot \bold{n} \cdot \bold{t} \f])
       * \param DsntDd[out]     derivative of sigma_nt w.r.t. displacements (\f[ \frac{ \mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{d}}
                                \f])
       * \param D2sntDd2[out]   second derivative of sigma_nt w.r.t. displacements (\f[ \frac{
                                \mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}
                                {\mathrm{d} \bold{d}^2} \f])
       * \param D2sntDdDn[out]  second derivative of sigma_nt w.r.t. displacements and vector n (\f[
                                \frac{\mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}
                                {\mathrm{d} \bold{d} \mathrm{d} \bold{n} } \f])
       * \param D2sntDdDt[out]  second derivative of sigma_nt w.r.t. displacements and vector t (\f[
                                \frac{\mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}
                                {\mathrm{d} \bold{d} \mathrm{d} \bold{t} } \f])
       * \param D2sntDdDxi[out] second derivative of sigma_nt w.r.t. displacements and local
                                parameter coordinate xi (\f[ \frac{\mathrm{d}^2 \bold{\sigma} \cdot
                                \bold{n} \cdot \bold{t}} {\mathrm{d} \bold{d} \mathrm{d} \bold{\xi}
                                } \f])
       * \param DsntDn[out]     derivative of sigma_nt w.r.t. vector n (\f[ \frac{ \mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{n}}
                                \f])
       * \param DsntDt[out]     derivative of sigma_nt w.r.t. vector n (\f[ \frac{ \mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{t}}
                                \f])
       * \param DsntDxi[out]    derivative of sigma_nt w.r.t. local parameter coordinate xi (\f[
                                \frac{ \mathrm{d} \bold{\sigma} \cdot \bold{n} \cdot
                                \bold{t}}{\mathrm{d} \bold{\xi}} \f])
       * \param temp[in]        temperature
       * \param DsntDT[out]     derivative of sigma_nt w.r.t. temperature (\f[ \frac{ \mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} T} \f])
       * \param D2sntDdDT[out]  second derivative of sigma_nt w.r.t. displacements and temperature
                                (\f[ \frac{\mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}
                                {\mathrm{d} \bold{d} \mathrm{d} T } \f])
       *
       * \note At the moment this method is only used for the nitsche contact formulation
       */
      void GetCauchyAtXi(const LINALG::Matrix<3, 1>& xi, const std::vector<double>& disp,
          const LINALG::Matrix<3, 1>& n, const LINALG::Matrix<3, 1>& t, double& sigma_nt,
          Epetra_SerialDenseMatrix* DsntDd, Epetra_SerialDenseMatrix* D2sntDd2,
          Epetra_SerialDenseMatrix* D2sntDdDn, Epetra_SerialDenseMatrix* D2sntDdDt,
          Epetra_SerialDenseMatrix* D2sntDdDxi, LINALG::Matrix<3, 1>* DsntDn,
          LINALG::Matrix<3, 1>* DsntDt, LINALG::Matrix<3, 1>* DsntDxi,
          const std::vector<double>* temp = NULL, Epetra_SerialDenseMatrix* DsntDT = NULL,
          Epetra_SerialDenseMatrix* D2sntDdDT = NULL) override;
      //@}

     protected:
      //! action parameters recognized by so_tet4
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
        struct_calc_store_istep,
        struct_calc_recover_istep,
        calc_struct_update_istep,
        calc_struct_reset_istep,     //!< reset elementwise internal variables
                                     //!< during iteration to last converged state
        calc_struct_reset_all,       //!< reset elementwise internal variables
                                     //!< to state in the beginning of the computation
        calc_struct_errornorms,      //!< compute error norms (L2,H1,energy)
        calc_global_gpstresses_map,  //! basically calc_struct_stress but with assembly of global
                                     //! gpstresses map
        prestress_update,
        calc_struct_energy,
        calc_struct_output_E,
        inversedesign_update,
        inversedesign_switch,
        multi_calc_dens,
        multi_readrestart
      };

      //! number of the material law
      // int                                                                      material_;

      //! container containing all kind of solid3 thisandthat
      // DRT::Container                                                           data_;
      //! vector of history variables for each Gauss point
      // vector<Teuchos::RCP<MAT::Material> >                                              mat_;
      //! volume of the element
      double V_;
      //! vector of partial derivatives in material frame
      // vector<LINALG::Matrix<NUMNOD_SOTET4,NUMDIM_SOTET4> > nxyz_;
      LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> nxyz_;

      /// prestressing switch & time
      INPAR::STR::PreStress pstype_;
      double pstime_;
      double time_;
      /// Prestressing object
      Teuchos::RCP<DRT::ELEMENTS::PreStress> prestress_;
      /// compute Jacobian mapping wrt to deformed configuration
      void UpdateJacobianMapping(
          const std::vector<double>& disp, DRT::ELEMENTS::PreStress& prestress);
      /// compute defgrd ypein all gp for given disp
      void DefGradient(const std::vector<double>& disp, Epetra_SerialDenseMatrix& gpdefgrd,
          DRT::ELEMENTS::PreStress& prestress);

      /*!
       * \brief Compute the deformation gradient
       *
       * \param defgrd Deformation gradient
       * \param xdisp Displacement vectir for each node (3x4)
       * \param gp Gauss point
       */
      void ComputeDeformationGradient(LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4>& defgrd,
          const LINALG::Matrix<NUMDIM_SOTET4, NUMNOD_SOTET4>& xdisp, const int gp) const;

      friend class DRT::ELEMENTS::InvDesign;
      Teuchos::RCP<DRT::ELEMENTS::InvDesign> invdesign_;

      // internal calculation methods

      // don't want = operator
      So_tet4& operator=(const So_tet4& old);

      //! init the inverse of the jacobian and its determinant in the material configuration
      virtual void InitJacobianMapping();

      //! Calculate nonlinear stiffness and mass matrix
      virtual void nlnstiffmass(std::vector<int>& lm,  // location matrix
          std::vector<double>& disp,                   // current displacements
          std::vector<double>* vel,                    // current velocities
          std::vector<double>* acc,                    // current accelerations
          std::vector<double>& residual,               // current residual displ
          std::vector<double>& dispmat,                // current material displacements
          LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4>* stiffmatrix,  // element stiffness matrix
          LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4>* massmatrix,   // element mass matrix
          LINALG::Matrix<NUMDOF_SOTET4, 1>* force,       // element internal force vector
          LINALG::Matrix<NUMDOF_SOTET4, 1>* forceinert,  // element inertial force vector
          LINALG::Matrix<NUMDOF_SOTET4, 1>* force_str,   // element structural force vector
          LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D>* elestress,    // stresses at GP
          LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D>* elestrain,    // strains at GP
          LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D>* eleplstrain,  // plastic strains at GP
          Teuchos::ParameterList& params,          // algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,   // stress output option
          const INPAR::STR::StrainType iostrain,   // strain output option
          const INPAR::STR::StrainType ioplstrain  // plastic strain output option
      );

      //! lump mass matrix (bborn 07/08)
      void so_tet4_lumpmass(
          LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4>* emass);  //!< element mass matrix


      /// Extrapolate Gauss-point values (e.g. stresses) to nodes and store them in elevecs
      void so_tet4_expol(
          LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D>& stresses,  ///< gp stresses
          Epetra_MultiVector& expolstresses                             ///< nodal stresses
      );

      //! remodeling for fibers at the end of time step (st 01/10)
      void so_tet4_remodel(std::vector<int>& lm,  // location matrix
          std::vector<double>& disp,              // current displacements
          Teuchos::ParameterList& params,         // algorithmic parameters e.g. time
          Teuchos::RCP<MAT::Material> mat);       // material

      //! Evaluate Tet4 Shapefcts at 1 gausspoint to keep them static
      const std::vector<LINALG::Matrix<NUMNOD_SOTET4, 1>> so_tet4_1gp_shapefcts();
      //! Evaluate Tet4 Derivs at 1 gausspoint to keep them static
      const std::vector<LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>> so_tet4_1gp_derivs();
      //! Evaluate Tet4 Weights at 1 gausspoint to keep them static
      const std::vector<double> so_tet4_1gp_weights();

      //! Evaluate Tet4 Shapefcts at 4 gausspoints to keep them static
      const std::vector<LINALG::Matrix<NUMNOD_SOTET4, 1>> so_tet4_4gp_shapefcts();
      //! Evaluate Tet4 Derivs at 4 gausspoints to keep them static
      const std::vector<LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>> so_tet4_4gp_derivs();
      //! Evaluate Tet4 Weights at 4 gausspoints to keep them static
      const std::vector<double> so_tet4_4gp_weights();

      //! @name Multi-scale related stuff

      /// Determine a homogenized material density for multi-scale analyses by averaging over the
      /// initial volume
      void sotet4_homog(Teuchos::ParameterList& params);

      /// Read restart on the microscale
      void sotet4_read_restart_multi();

      //@}

      /*!
       * Executes the post setup call for all materials. This method will be called once per element
       * at the first Evaluate call.
       *
       * @param params ParameterList to be passed to the materials
       */
      void MaterialPostSetup(Teuchos::ParameterList& params) override;
    };  // class So_tet4


  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // #ifndef SO_TET4_H_
