/*----------------------------------------------------------------------*/
/*! \file
\brief Solid Tet10 Element
\level 1
*----------------------------------------------------------------------*/
#ifndef SO_TET10_H
#define SO_TET10_H


#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_parobjectfactory.H"

#include "../drt_mat/material.H"
#include "../drt_inpar/inpar_structure.H"

#include "so_base.H"

#include "../drt_fem_general/drt_utils_integration.H"

// gee: note that these are also defined in so_integrator.H
#define NUMNOD_SOTET10 10       ///< number of nodes
#define NODDOF_SOTET10 3        ///< number of dofs per node
#define NUMDOF_SOTET10 30       ///< total dofs per element
#define NUMGPT_SOTET10 4        ///< number gauss points per element for stiffness integration
#define NUMGPT_MASS_SOTET10 11  ///< number gauss points per element for mass integration
#define NUMDIM_SOTET10 3        ///< number of dimensions
#define NUMCOORD_SOTET10 4      ///< number of shape function coordinates (ksi1-ksi4)
#define NUMNOD_SOTET10_FACE 6   ///< number of nodes on a TET10 face (which is a TRI6)
#define NUMGPT_SOTET10_FACE 3   ///< number of GP    on a TET10 face (which is a TRI6)

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    // forward declarations
    class InvDesign;
    class PreStress;

    class So_tet10Type : public DRT::ElementType
    {
     public:
      std::string Name() const { return "So_tet10Type"; }

      static So_tet10Type& Instance();

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
      static So_tet10Type instance_;
    };

    /*!
    \brief A C++ version of the 10-node tet solid element

    a structural 10-node tetrahedral solid displacement element
    \author biehler
    */

    class So_tet10 : public So_base
    {
     public:
      //! @name Friends
      friend class So_tet10Type;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owning processor
      */
      So_tet10(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      So_tet10(const So_tet10& old);

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
      virtual int UniqueParObjectId() const { return So_tet10Type::Instance().UniqueParObjectId(); }

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
      virtual ~So_tet10();
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

      virtual DRT::ElementType& ElementType() const { return So_tet10Type::Instance(); }

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
      \endcode*/

      virtual void VisNames(std::map<std::string, int>&
              names);  ///< to be filled with key names of data to visualize and with int dimensions

      /*!
      \brief Query data to be visualized using BINIO of a given name

      The method is supposed to call this base method to visualize the owner of
      the element.
      If the derived method recognizes a supported data name, it shall fill it
      with corresponding data.
      If it does NOT recognizes the name, it shall do nothing.

      \warning The method must not change size of data

      */
      virtual bool VisData(
          const std::string& name,  ///< Name of data that is currently processed for visualization
          std::vector<double>&
              data  ///< d ata to be filled by element if element recognizes the name
      );
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

      Evaluate so_tet10 element stiffness, mass, internal forces, etc.

      \return 0 if successful, negative otherwise
      */
      int Evaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  ///< pointer to discretization for de-assembly
          std::vector<int>& lm,                 ///< location matrix for de-assembly
          Epetra_SerialDenseMatrix& elemat1,    ///< (stiffness-)matrix to be filled by element.
          Epetra_SerialDenseMatrix& elemat2,    ///< (mass-)matrix to be filled by element.
          Epetra_SerialDenseVector& elevec1,    ///< (internal force-)vector to be filled by element
          Epetra_SerialDenseVector& elevec2,    ///< vector to be filled by element
          Epetra_SerialDenseVector& elevec3     ///< vector to be filled by element
      );


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

      //@}

     protected:
      //! action parameters recognized by so_tet10
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
        calc_struct_reset_istep,  //!< reset elementwise internal variables
                                  //!< during iteration to last converged state
        calc_struct_reset_all,    //!< reset elementwise internal variables
                                  //!< to state in the beginning of the computation
        calc_struct_energy,       // compute energy
        calc_struct_errornorms,   //!< compute error norms (L2,H1,energy)
        prestress_update,
        postprocess_stress,
        calc_global_gpstresses_map
      };


      //! container containing all kind of solid3 thisandthat
      DRT::Container data_;
      //! vector of inverses of the jacobian in material frame
      std::vector<LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10>> invJ_;
      // determinant of Jacobian in material frame
      std::vector<double> detJ_;
      // stuff consistent mass matrix
      //! vector of inverses of the jacobian in material frame
      std::vector<LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10>> invJ_mass_;
      //! determinant of Jacobian in material frame
      std::vector<double> detJ_mass_;

      //! vector of partial derivatives in material frame
      // vector<LINALG::Matrix<NUMNOD_SOTET4,NUMDIM_SOTET4> > nxyz_;
      LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> nxyz_;

      /// prestressing switch & time
      INPAR::STR::PreStress pstype_;
      double pstime_;
      double time_;
      /// Prestressing object
      Teuchos::RCP<DRT::ELEMENTS::PreStress> prestress_;
      /// compute Jacobian mapping wrt to deformed configuration
      void UpdateJacobianMapping(
          const std::vector<double>& disp, DRT::ELEMENTS::PreStress& prestress);

      //! Update history variables at the end of time step (fiber direction, inelastic deformation)
      //! (gebauer 07/19)
      void Update_element(std::vector<double>& disp,  // current displacements
          Teuchos::ParameterList& params,             // algorithmic parameters e.g. time
          Teuchos::RCP<MAT::Material> mat);           // material

      /// compute defgrd in all gp for given disp
      void DefGradient(const std::vector<double>& disp, Epetra_SerialDenseMatrix& gpdefgrd,
          DRT::ELEMENTS::PreStress& prestress);

      friend class DRT::ELEMENTS::InvDesign;
      Teuchos::RCP<DRT::ELEMENTS::InvDesign> invdesign_;


      // internal calculation methods

      // don't want = operator
      So_tet10& operator=(const So_tet10& old);

      //! init the inverse of the jacobian and its determinant in the material configuration
      virtual void InitJacobianMapping();

      //! Calculate nonlinear stiffness and mass matrix
      virtual void so_tet10_nlnstiffmass(std::vector<int>& lm,  ///< location matrix
          std::vector<double>& disp,                            ///< current displacements
          std::vector<double>* vel,                             ///< current velocities
          std::vector<double>* acc,                             ///< current accelerations
          std::vector<double>& residual,                        ///< current residual displ
          std::vector<double>& dispmat,                         ///< current material displacements
          LINALG::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10>*
              stiffmatrix,                                             ///< element stiffness matrix
          LINALG::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10>* massmatrix,  ///< element mass matrix
          LINALG::Matrix<NUMDOF_SOTET10, 1>* force,       ///< element internal force vector
          LINALG::Matrix<NUMDOF_SOTET10, 1>* forceinert,  ///< element inertial force vector
          LINALG::Matrix<NUMDOF_SOTET10, 1>* force_str,   ///< element structural force vector
          LINALG::Matrix<NUMGPT_SOTET10, MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          LINALG::Matrix<NUMGPT_SOTET10, MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          Teuchos::ParameterList& params,          ///< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,   ///< stress output option
          const INPAR::STR::StrainType iostrain);  ///< strain output option

      //! lump mass matrix (bborn 07/08)
      void so_tet10_lumpmass(LINALG::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10>* emass);

      /// Extrapolate Gauss-point values (e.g. stresses) to nodes and store them in elevecs
      void so_tet10_expol(
          LINALG::Matrix<NUMGPT_SOTET10, MAT::NUM_STRESS_3D>& stresses,  ///< gp stresses
          Epetra_MultiVector& expolstresses                              ///< nodal stresses
      );

      // These functions are basically copied from the So_integrator,
      // along with the rather inconsitent matrix sizes.

      //! Evaluate Tet10 Shapefcts at 4 gausspoints to keep them static
      const std::vector<LINALG::Matrix<NUMNOD_SOTET10, 1>> so_tet10_4gp_shapefcts();
      //! Evaluate Tet10 Derivs at 4 gausspoints to keep them static
      const std::vector<LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>>& so_tet10_4gp_derivs();

      /*!
       * \brief Evaluate the first derivative of the shape functions at the Gauss point with the
       * 4 Gauss-point integration rule
       *
       * \param derivs First derivatives of the shape functions
       * \param gp Gauss point
       */
      template <DRT::UTILS::GaussRule3D intrule>
      void so_tet10_derivs(LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>& derivs, int gp) const;

      //! Evaluate Tet10 Weights at 4 gausspoints to keep them static
      const std::vector<double>& so_tet10_4gp_weights();

      //! Evaluate Tet10 Shapefcts at 10 gausspoints to keep them static
      const std::vector<LINALG::Matrix<NUMNOD_SOTET10, 1>>& so_tet10_11gp_shapefcts();
      //! Evaluate Tet10 Derivs at 10 gausspoints to keep them static
      const std::vector<LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>>& so_tet10_11gp_derivs();
      //! Evaluate Tet10 Weights at 10 gausspoints to keep them static
      const std::vector<double>& so_tet10_11gp_weights();

      /*!
       * Executes the post setup call for all materials. This method will be called once per element
       * at the first Evaluate call.
       *
       * @param params ParameterList to be passed to the materials
       */
      void MaterialPostSetup(Teuchos::ParameterList& params) override;
    };  // class So_tet10



  }  // namespace ELEMENTS
}  // namespace DRT


#endif  // #ifndef SO_TET10_H_
