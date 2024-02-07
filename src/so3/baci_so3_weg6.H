/*----------------------------------------------------------------------*/
/*! \file
\brief Solid Wedge6 Element
\level 1


*----------------------------------------------------------------------*/
#ifndef BACI_SO3_WEG6_H
#define BACI_SO3_WEG6_H

#include "baci_config.H"

#include "baci_inpar_structure.H"
#include "baci_lib_element.H"
#include "baci_lib_elementtype.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_mat_material.H"
#include "baci_so3_base.H"

#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

// Several parameters which are fixed for Solid Wedge6
const int NUMNOD_WEG6 = 6;   ///< number of nodes
const int NODDOF_WEG6 = 3;   ///< number of dofs per node
const int NUMDOF_WEG6 = 18;  ///< total dofs per element
const int NUMGPT_WEG6 = 6;   ///< total gauss points per element
const int NUMDIM_WEG6 = 3;   ///< number of dimensions

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    // forward declarations
    class PreStress;

    class So_weg6Type : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "So_weg6Type"; }

      static So_weg6Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_weg6Type instance_;

      std::string GetElementTypeString() const { return "SOLIDW6"; }
    };

    /*!
    \brief A C++ version of the 6-node wedge solid element

    */
    class So_weg6 : public So_base
    {
     public:
      //! @name Friends
      friend class So_weg6Type;

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
      DRT::Element* Clone() const override;

      /*!
      \brief Get shape type of element
      */
      CORE::FE::CellType Shape() const override;

      /*!
      \brief Return number of volumes of this element
      */
      inline int NumVolume() const override { return 1; }

      /*!
      \brief Return number of surfaces of this element
      */
      inline int NumSurface() const override { return 5; }

      /*!
      \brief Return number of lines of this element
      */
      inline int NumLine() const override { return 9; }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element

      */
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element

      */
      std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

      /*!
        \brief Get coordinates of element center

        */
      virtual std::vector<double> ElementCenterRefeCoords();

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      inline int UniqueParObjectId() const override
      {
        return So_weg6Type::Instance().UniqueParObjectId();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(CORE::COMM::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Unpack(const std::vector<char>& data) override;


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
      int NumDofPerNode(const DRT::Node& node) const override { return 3; }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual DRT::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      int NumDofPerElement() const override { return 0; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return So_weg6Type::Instance(); }

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
      void VisNames(std::map<std::string, int>& names) override;

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
      bool VisData(const std::string& name, std::vector<double>& data) override;

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          INPUT::LineDefinition* linedef) override;

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate so_hex8 element stiffness, mass, internal forces, etc.

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param discretization : pointer to discretization for de-assembly
      \param lm (in)        : location matrix for de-assembly
      \param elemat1 (out)  : (stiffness-)matrix to be filled by element. If nullptr on input,
                              the controling method does not expect the element to fill
                              this matrix.
      \param elemat2 (out)  : (mass-)matrix to be filled by element. If nullptr on input,
                              the controling method does not expect the element to fill
                              this matrix.
      \param elevec1 (out)  : (internal force-)vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override;


      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a surface Neumann condition on the solid3 element

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param condition (in)     : The condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If nullptr on input,

      \return 0 if successful, negative otherwise
      */
      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;

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
        calc_struct_reset_istep,  //!< reset elementwise internal variables
                                  //!< during iteration to last converged state
        calc_struct_reset_all,    //!< reset elementwise internal variables
                                  //!< to state in the beginning of the computation
        calc_struct_energy,
        prestress_update,
        calc_global_gpstresses_map,  //! basically calc_struct_stress but with assembly of global
                                     //! gpstresses map
        calc_recover

      };

      //! vector of inverses of the jacobian in material frame
      std::vector<CORE::LINALG::Matrix<NUMDIM_WEG6, NUMDIM_WEG6>> invJ_;
      //! determinant of Jacobian in material frame
      std::vector<double> detJ_;

      /// prestressing switch & time
      INPAR::STR::PreStress pstype_;
      double pstime_;
      double time_;
      /// Prestressing object
      Teuchos::RCP<DRT::ELEMENTS::PreStress> prestress_;
      // compute Jacobian mapping wrt to deformed configuration
      void UpdateJacobianMapping(
          const std::vector<double>& disp, DRT::ELEMENTS::PreStress& prestress);
      // compute defgrd in all gp for given disp
      void DefGradient(const std::vector<double>& disp, CORE::LINALG::SerialDenseMatrix& gpdefgrd,
          DRT::ELEMENTS::PreStress& prestress);


      // internal calculation methods

      // don't want = operator
      So_weg6& operator=(const So_weg6& old);

      //! init the inverse of the jacobian and its determinant in the material configuration
      virtual void InitJacobianMapping();

      //! Calculate nonlinear stiffness and mass matrix
      virtual void sow6_nlnstiffmass(std::vector<int>& lm,  ///< location matrix
          std::vector<double>& disp,                        ///< current displacements
          std::vector<double>* vel,                         ///< current velocities
          std::vector<double>* acc,                         ///< current accelerations
          std::vector<double>& residual,                    ///< current residual displ
          std::vector<double>& dispmat,                     ///< current material displacements
          CORE::LINALG::Matrix<NUMDOF_WEG6, NUMDOF_WEG6>*
              stiffmatrix,                                             ///< element stiffness matrix
          CORE::LINALG::Matrix<NUMDOF_WEG6, NUMDOF_WEG6>* massmatrix,  ///< element mass matrix
          CORE::LINALG::Matrix<NUMDOF_WEG6, 1>* force,       ///< element internal force vector
          CORE::LINALG::Matrix<NUMDOF_WEG6, 1>* forceinert,  ///< element inertial force vector
          CORE::LINALG::Matrix<NUMDOF_WEG6, 1>* force_str,   ///< element structural force vector
          CORE::LINALG::Matrix<NUMGPT_WEG6, MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          CORE::LINALG::Matrix<NUMGPT_WEG6, MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          Teuchos::ParameterList& params,          ///< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,   ///< stress output option
          const INPAR::STR::StrainType iostrain);  ///< strain output option

      //! remodeling for fibers at the end of time step (st 01/10)
      void sow6_remodel(std::vector<int>& lm,       // location matrix
          std::vector<double>& disp,                // current displacements
          Teuchos::ParameterList& params,           // algorithmic parameters e.g. time
          const Teuchos::RCP<MAT::Material>& mat);  // material

      //! Evaluate Wedge6 Shapefcts to keep them static
      std::vector<CORE::LINALG::Matrix<NUMNOD_WEG6, 1>> sow6_shapefcts();
      //! Evaluate Wedge6 Derivs to keep them static
      std::vector<CORE::LINALG::Matrix<NUMDIM_WEG6, NUMNOD_WEG6>> sow6_derivs();
      //! Evaluate Wedge6 Weights to keep them static
      std::vector<double> sow6_weights();


      //! calculate static shape functions and derivatives for sow6
      void sow6_shapederiv(CORE::LINALG::Matrix<NUMNOD_WEG6, NUMGPT_WEG6>** shapefct,
          CORE::LINALG::Matrix<NUMDOF_WEG6, NUMNOD_WEG6>** deriv,
          CORE::LINALG::Matrix<NUMGPT_WEG6, 1>** weights);

      //! lump mass matrix (bborn 07/08)
      void sow6_lumpmass(CORE::LINALG::Matrix<NUMDOF_WEG6, NUMDOF_WEG6>* emass);

     private:
      std::string GetElementTypeString() const { return "SOLIDW6"; }
    };  // class So_weg6


  }  // namespace ELEMENTS
}  // namespace DRT


BACI_NAMESPACE_CLOSE

#endif  // SO3_WEG6_H
