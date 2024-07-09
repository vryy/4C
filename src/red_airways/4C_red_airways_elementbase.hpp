/*---------------------------------------------------------------------*/
/*! \file

\brief Header for implementations of the various reduced elements


\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_RED_AIRWAYS_ELEMENTBASE_HPP
#define FOUR_C_RED_AIRWAYS_ELEMENTBASE_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_elementtype.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_red_airways_elem_params.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{

  namespace ELEMENTS
  {
    // forward declarations
    class RedAirwayImplInterface;
    template <Core::FE::CellType distype>
    class AirwayImpl;

    class RedAirwayType : public Core::Elements::ElementType
    {
     public:
      std::string name() const override { return "RedAirwayType"; }

      static RedAirwayType& instance();

      Core::Communication::ParObject* create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override
      {
      }

      Core::LinAlg::SerialDenseMatrix compute_null_space(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override
      {
        Core::LinAlg::SerialDenseMatrix nullspace;
        FOUR_C_THROW("method ComputeNullSpace not implemented");
        return nullspace;
      }

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static RedAirwayType instance_;
    };

    /*!
    \brief A C++ wrapper for the airway element

    */
    class RedAirway : public Core::Elements::Element
    {
     public:
      //! @name Friends
      friend class RedAirwayImplInterface;
      friend class AirwayImpl<Core::FE::CellType::line2>;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      */
      RedAirway(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      RedAirway(const RedAirway& old);

      /*!
      \brief Deep copy this instance of RedAirway and return pointer to the copy

      The clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      Core::Elements::Element* clone() const override;

      /*!
      \brief Get shape type of element
      */
      Core::FE::CellType shape() const override;

      /*!
      \brief Return number of lines of this element
      */
      int num_line() const override
      {
        if (num_node() == 2)
          return 1;
        else
        {
          FOUR_C_THROW("Could not determine number of lines");
          return -1;
        }
      }

      /*!
      \brief Return number of surfaces of this element (always 1)
      */
      int num_surface() const override { return -1; }

      /*!
      \brief Return number of volumes of this element (always 1)
      */
      int num_volume() const override { return -1; }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      std::vector<Teuchos::RCP<Core::Elements::Element>> lines() override;


      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int unique_par_object_id() const override
      {
        return RedAirwayType::instance().unique_par_object_id();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref pack and \ref unpack are used to communicate this element

      */
      void pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref pack and \ref unpack are used to communicate this element

      */
      void unpack(const std::vector<char>& data) override;


      //@}

      //! @name Acess methods


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual Core::Elements::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      int num_dof_per_node(const Core::Nodes::Node& node) const override { return 1; }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual Core::Elements::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      int num_dof_per_element() const override { return 0; }

      /*!
      \brief Print this element
      */
      void print(std::ostream& os) const override;

      RedAirwayType& element_type() const override { return RedAirwayType::instance(); }

      //@}

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
      void vis_names(std::map<std::string, int>& names) override{};

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
      bool vis_data(const std::string& name, std::vector<double>& data) override;


      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool read_element(const std::string& eletype, const std::string& distype,
          Input::LineDefinition* linedef) override;

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate airway element stiffness, mass, internal forces etc

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param elemat1 (out)  : matrix to be filled by element. If nullptr on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elemat2 (out)  : matrix to be filled by element. If nullptr on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elevec1 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      int evaluate(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;

      /*!
      \brief Evaluate a Neumann boundary condition

      An element derived from this class uses the evaluate_neumann method to receive commands
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
      int evaluate_neumann(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseMatrix* elemat1 = nullptr) override;

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a line Neumann condition on the airway element

      \return 0 if successful, negative otherwise
      */
      virtual int evaluate_dirichlet(Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1);


      //@}

      //! @name Other
      std::string type() { return elem_type_; }

      //! @name Other
      std::string resistance() { return resistance_; }

      //! @name Other
      std::string elem_solving_type() { return elemsolving_type_; }

      //@}

      //! get cross-sectional area
      //  double         getdata(){ return A_;}

      //! get youngs modulus of the wall
      //  double         getEw(){ return Ew_;}

      //! get youngs modulus of the air
      //  double         getEa(){ return Ea_;}

      //! get wall thickness
      //  double         gettw(){ return tw_;}

      //! set qin
      //  void           setqin (double Qin ){ qin_ = Qin;}

      //! set qout
      //  void           setqout(double Qout){ qout_ = Qout;}

      /*!
       * \brief Get fixed airway parameters of the RedAirway element
       */
      [[nodiscard]] const ReducedLung::AirwayParams& get_airway_params() const;

     private:
      //! action parameters recognized by airway
      enum ActionType
      {
        none,
        calc_sys_matrix_rhs,
        calc_sys_matrix_rhs_iad,
        get_initial_state,
        set_bc,
        calc_flow_rates,
        get_coupled_values,
        get_junction_volume_mix,
        solve_scatra,
        solve_junction_scatra,
        calc_cfl,
        solve_blood_air_transport,
        update_scatra,
        update_elem12_scatra,
        eval_nodal_ess_vals,
        eval_PO2_from_concentration,
        calc_elem_volumes
      };



      //! Element Type
      std::string elem_type_;

      //! Resistance Type
      std::string resistance_;

      //! Solver Type
      std::string elemsolving_type_;

      //! Airway-specific parameters
      Discret::ReducedLung::AirwayParams airway_params_;


      // internal calculation methods

      // don't want = operator
      RedAirway& operator=(const RedAirway& old);


      /// set number of gauss points to element shape default
      Core::FE::GaussRule1D get_optimal_gaussrule(const Core::FE::CellType& distype);

      /*!
       * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are
       * necessary \return boolean indicating higher order status
       */
      bool is_higher_order_element(const Core::FE::CellType distype  ///< discretization type
      ) const;


    };  // class RedAirway


    class RedAcinusType : public Core::Elements::ElementType
    {
     public:
      std::string name() const override { return "RedAcinusType"; }

      static RedAcinusType& instance();

      Core::Communication::ParObject* create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override
      {
      }

      Core::LinAlg::SerialDenseMatrix compute_null_space(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override
      {
        Core::LinAlg::SerialDenseMatrix nullspace;
        FOUR_C_THROW("method ComputeNullSpace not implemented");
        return nullspace;
      }

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static RedAcinusType instance_;
    };

    /*!
    \brief A C++ wrapper for the acinus element

    */
    class RedAcinus : public Core::Elements::Element
    {
     public:
      //! @name Friends
      friend class RedAirwayImplInterface;
      friend class AirwayImpl<Core::FE::CellType::line2>;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      */
      RedAcinus(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      RedAcinus(const RedAcinus& old);

      /*!
      \brief Deep copy this instance of Redacinus and return pointer to the copy

      The clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      Core::Elements::Element* clone() const override;

      /*!
      \brief Get shape type of element
      */
      Core::FE::CellType shape() const override;

      /*!
      \brief Return number of lines of this element
      */
      int num_line() const override
      {
        if (num_node() == 2)
          return 1;
        else
        {
          FOUR_C_THROW("Could not determine number of lines");
          return -1;
        }
      }

      /*!
      \brief Return number of surfaces of this element (always 1)
      */
      int num_surface() const override { return -1; }

      /*!
      \brief Return number of volumes of this element (always 1)
      */
      int num_volume() const override { return -1; }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      std::vector<Teuchos::RCP<Core::Elements::Element>> lines() override;

      /*!
      \brief Return center coordinates of element
      */
      virtual std::vector<double> element_center_refe_coords();

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int unique_par_object_id() const override
      {
        return RedAcinusType::instance().unique_par_object_id();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref pack and \ref unpack are used to communicate this element

      */
      void pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref pack and \ref unpack are used to communicate this element

      */
      void unpack(const std::vector<char>& data) override;


      //@}

      //! @name Acess methods


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual Core::Elements::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      int num_dof_per_node(const Core::Nodes::Node& node) const override { return 1; }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual Core::Elements::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      int num_dof_per_element() const override { return 0; }

      /*!
      \brief Print this element
      */
      void print(std::ostream& os) const override;

      RedAcinusType& element_type() const override { return RedAcinusType::instance(); }

      //@}

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
      void vis_names(std::map<std::string, int>& names) override;

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
      bool vis_data(const std::string& name, std::vector<double>& data) override;


      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool read_element(const std::string& eletype, const std::string& distype,
          Input::LineDefinition* linedef) override;

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate acinus element stiffness, mass, internal forces etc

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param elemat1 (out)  : matrix to be filled by element. If nullptr on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elemat2 (out)  : matrix to be filled by element. If nullptr on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elevec1 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      int evaluate(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;

      /*!
      \brief Evaluate a Neumann boundary condition

      An element derived from this class uses the evaluate_neumann method to receive commands
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
      int evaluate_neumann(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseMatrix* elemat1 = nullptr) override;

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a line Neumann condition on the acinus element

      \return 0 if successful, negative otherwise
      */
      virtual int evaluate_dirichlet(Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1);


      //@}


      //! @name Other
      std::string type() { return elem_type_; }

      //! @name Other
      std::string resistance() { return resistance_; }

      //@}

      //! get cross-sectional area
      //  double         getdata(){ return A_;}

      //! get youngs modulus of the wall
      //  double         getEw(){ return Ew_;}

      //! get youngs modulus of the air
      //  double         getEa(){ return Ea_;}

      //! get wall thickness
      //  double         gettw(){ return tw_;}

      //! set qin
      //  void           setqin (double Qin ){ qin_ = Qin;}

      //! set qout
      //  void           setqout(double Qout){ qout_ = Qout;}

      /*!
       * \brief Update relaxed reference volume of the RedAcinus element.
       *
       * Prestressed volume states in the input file may differ from the actual V0 used in many
       * formulas.
       */
      void update_relaxed_volume(double newVol);

      /*!
       * \brief Get fixed acinus parameters of the RedAcinus element
       */
      [[nodiscard]] const ReducedLung::AcinusParams& get_acinus_params() const;

     private:
      //! action parameters recognized by acinus
      enum ActionType
      {
        none,
        calc_sys_matrix_rhs,
        calc_sys_matrix_rhs_iad,
        get_initial_state,
        set_bc,
        calc_flow_rates,
        get_coupled_values,
        get_junction_volume_mix,
        solve_scatra,
        solve_junction_scatra,
        calc_cfl,
        solve_blood_air_transport,
        update_scatra,
        update_elem12_scatra,
        eval_nodal_ess_vals,
        eval_PO2_from_concentration,
        calc_elem_volumes
      };


      //! Element Type
      std::string elem_type_;

      //! Resistance Type
      std::string resistance_;

      //! Acinus-specific parameters
      Discret::ReducedLung::AcinusParams acinus_params_;


      // internal calculation methods

      // don't want = operator
      RedAcinus& operator=(const RedAcinus& old);


      /// set number of gauss points to element shape default
      Core::FE::GaussRule1D get_optimal_gaussrule(const Core::FE::CellType& distype);

      /*!
       * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are
       * necessary \return boolean indicating higher order status
       */
      bool is_higher_order_element(const Core::FE::CellType distype  ///< discretization type
      ) const;


    };  // class RedAcinus

    /*!
    \brief A C++ wrapper for the inter acinar dependency element

    */

    class RedInterAcinarDepType : public Core::Elements::ElementType
    {
     public:
      std::string name() const override { return "RedInterAcinarDepType"; }

      static RedInterAcinarDepType& instance();

      Core::Communication::ParObject* create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override
      {
      }

      Core::LinAlg::SerialDenseMatrix compute_null_space(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override
      {
        Core::LinAlg::SerialDenseMatrix nullspace;
        FOUR_C_THROW("method ComputeNullSpace not implemented");
        return nullspace;
      }

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static RedInterAcinarDepType instance_;
    };

    /*!
    \brief A C++ wrapper for the inter acinar dependency element

    */
    class RedInterAcinarDep : public Core::Elements::Element
    {
     public:
      //! @name Friends
      friend class RedAirwayImplInterface;
      friend class AirwayImpl<Core::FE::CellType::line2>;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      */
      RedInterAcinarDep(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      RedInterAcinarDep(const RedInterAcinarDep& old);

      /*!
      \brief Deep copy this instance of Redacinus and return pointer to the copy

      The clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      Core::Elements::Element* clone() const override;

      /*!
      \brief Get shape type of element
      */
      Core::FE::CellType shape() const override;

      /*!
      \brief Return number of lines of this element
      */
      int num_line() const override
      {
        if (num_node() == 2)
          return 1;
        else
        {
          FOUR_C_THROW("Could not determine number of lines");
          return -1;
        }
      }

      /*!
      \brief Return number of surfaces of this element (always 1)
      */
      int num_surface() const override { return -1; }

      /*!
      \brief Return number of volumes of this element (always 1)
      */
      int num_volume() const override { return -1; }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      std::vector<Teuchos::RCP<Core::Elements::Element>> lines() override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int unique_par_object_id() const override
      {
        return RedInterAcinarDepType::instance().unique_par_object_id();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref pack and \ref unpack are used to communicate this element

      */
      void pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref pack and \ref unpack are used to communicate this element

      */
      void unpack(const std::vector<char>& data) override;


      //@}

      //! @name Acess methods


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual Core::Elements::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      int num_dof_per_node(const Core::Nodes::Node& node) const override { return 1; }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual Core::Elements::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      int num_dof_per_element() const override { return 0; }

      /*!
      \brief Print this element
      */
      void print(std::ostream& os) const override;

      RedInterAcinarDepType& element_type() const override
      {
        return RedInterAcinarDepType::instance();
      }

      //@}

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
      void vis_names(std::map<std::string, int>& names) override;

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
      bool vis_data(const std::string& name, std::vector<double>& data) override;


      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool read_element(const std::string& eletype, const std::string& distype,
          Input::LineDefinition* linedef) override;

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate inter acinar dependency element stiffness, mass, internal forces etc

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param elemat1 (out)  : matrix to be filled by element. If nullptr on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elemat2 (out)  : matrix to be filled by element. If nullptr on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elevec1 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      int evaluate(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;

      /*!
      \brief Evaluate a Neumann boundary condition

      An element derived from this class uses the evaluate_neumann method to receive commands
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
      int evaluate_neumann(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseMatrix* elemat1 = nullptr) override;

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a line Neumann condition on the acinus element

      \return 0 if successful, negative otherwise
      */
      virtual int evaluate_dirichlet(Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1);


      //@}


      //! get cross-sectional area
      //  double         getdata(){ return A_;}

      //! get youngs modulus of the wall
      //  double         getEw(){ return Ew_;}

      //! get youngs modulus of the air
      //  double         getEa(){ return Ea_;}

      //! get wall thickness
      //  double         gettw(){ return tw_;}

      //! set qin
      //  void           setqin (double Qin ){ qin_ = Qin;}

      //! set qout
      //  void           setqout(double Qout){ qout_ = Qout;}

      //! get element parameters
      void get_params(std::string name, double& var);

      //! get element parameters
      void get_params(std::string name, int& var);

     private:
      //! action parameters recognized by inter acinar dependency
      enum ActionType
      {
        none,
        calc_sys_matrix_rhs,
        calc_sys_matrix_rhs_iad,
        get_initial_state,
        set_bc,
        calc_flow_rates,
        get_coupled_values,
        get_junction_volume_mix,
        solve_scatra,
        calc_cfl,
        solve_blood_air_transport,
        update_scatra,
        update_elem12_scatra,
        eval_nodal_ess_vals,
        eval_PO2_from_concentration,
        calc_elem_volumes
      };

      // data
      std::map<std::string, double> elem_params_;

      // element tree generation
      int generation_;

      // internal calculation methods

      // don't want = operator
      RedInterAcinarDep& operator=(const RedInterAcinarDep& old);


      /// set number of gauss points to element shape default
      Core::FE::GaussRule1D get_optimal_gaussrule(const Core::FE::CellType& distype);

      /*!
       * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are
       * necessary \return boolean indicating higher order status
       */
      bool is_higher_order_element(const Core::FE::CellType distype  ///< discretization type
      ) const;


    };  // class RedInterAcinarDep

    /*!
    \brief A C++ wrapper for the Air/Blood diffusion element

    */

    class RedAirBloodScatraType : public Core::Elements::ElementType
    {
     public:
      std::string name() const override { return "RedAirBloodScatraType"; }

      static RedAirBloodScatraType& instance();

      Core::Communication::ParObject* create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override
      {
      }

      Core::LinAlg::SerialDenseMatrix compute_null_space(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override
      {
        Core::LinAlg::SerialDenseMatrix nullspace;
        FOUR_C_THROW("method ComputeNullSpace not implemented");
        return nullspace;
      }

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static RedAirBloodScatraType instance_;
    };  // class RedAirBloodScatraType

    class RedAirBloodScatra : public Core::Elements::Element
    {
     public:
      //! @name Friends
      friend class RedAirwayImplInterface;
      friend class AirwayImpl<Core::FE::CellType::line2>;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      */
      RedAirBloodScatra(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      RedAirBloodScatra(const RedAirBloodScatra& old);

      /*!
      \brief Deep copy this instance of Redacinus and return pointer to the copy

      The clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      Core::Elements::Element* clone() const override;

      /*!
      \brief Get shape type of element
      */
      Core::FE::CellType shape() const override;

      /*!
      \brief Return number of lines of this element
      */
      int num_line() const override
      {
        if (num_node() == 2)
          return 1;
        else
        {
          FOUR_C_THROW("Could not determine number of lines");
          return -1;
        }
      }

      /*!
      \brief Return number of surfaces of this element (always 1)
      */
      int num_surface() const override { return -1; }

      /*!
      \brief Return number of volumes of this element (always 1)
      */
      int num_volume() const override { return -1; }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      std::vector<Teuchos::RCP<Core::Elements::Element>> lines() override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int unique_par_object_id() const override
      {
        return RedAirBloodScatraType::instance().unique_par_object_id();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref pack and \ref unpack are used to communicate this element

      */
      void pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref pack and \ref unpack are used to communicate this element

      */
      void unpack(const std::vector<char>& data) override;


      //@}

      //! @name Acess methods


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual Core::Elements::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      int num_dof_per_node(const Core::Nodes::Node& node) const override { return 1; }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual Core::Elements::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      int num_dof_per_element() const override { return 0; }

      /*!
      \brief Print this element
      */
      void print(std::ostream& os) const override;

      RedAirBloodScatraType& element_type() const override
      {
        return RedAirBloodScatraType::instance();
      }

      //@}

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
      void vis_names(std::map<std::string, int>& names) override;

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
      bool vis_data(const std::string& name, std::vector<double>& data) override;


      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool read_element(const std::string& eletype, const std::string& distype,
          Input::LineDefinition* linedef) override;

      //@}

      //@}


      //! get cross-sectional area
      //  double         getdata(){ return A_;}

      //! get youngs modulus of the wall
      //  double         getEw(){ return Ew_;}

      //! get youngs modulus of the air
      //  double         getEa(){ return Ea_;}

      //! get wall thickness
      //  double         gettw(){ return tw_;}

      //! set qin
      //  void           setqin (double Qin ){ qin_ = Qin;}

      //! set qout
      //  void           setqout(double Qout){ qout_ = Qout;}

      //! get element parameters
      void get_params(std::string name, double& var);

      //! get element parameters
      void get_params(std::string name, int& var);

     private:
      //! action parameters recognized by inter acinar dependency
      enum ActionType
      {
        none,
        calc_sys_matrix_rhs,
        calc_sys_matrix_rhs_iad,
        get_initial_state,
        set_bc,
        calc_flow_rates,
        get_coupled_values,
        get_junction_volume_mix,
        solve_scatra,
        calc_cfl,
        solve_junction_scatra,
        solve_blood_air_transport,
        update_scatra,
        update_elem12_scatra,
        eval_nodal_ess_vals,
        eval_PO2_from_concentration,
        calc_elem_volumes
      };

      // data
      std::map<std::string, double> elem_params_;

      // element tree generation
      int generation_;

      // internal calculation methods

      // don't want = operator
      RedAirBloodScatra& operator=(const RedAirBloodScatra& old);
    };  // class RedAirBloodScatra


    /*!
    \brief A C++ wrapper for the Air/Blood diffusion element

    */
    class RedAirBloodScatraLine3Type : public Core::Elements::ElementType
    {
     public:
      std::string name() const override { return "RedAirBloodScatraLine3Type"; }

      static RedAirBloodScatraLine3Type& instance();

      Core::Communication::ParObject* create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override
      {
      }

      Core::LinAlg::SerialDenseMatrix compute_null_space(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override
      {
        Core::LinAlg::SerialDenseMatrix nullspace;
        FOUR_C_THROW("method ComputeNullSpace not implemented");
        return nullspace;
      }

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static RedAirBloodScatraLine3Type instance_;
    };  // class RedAirBloodScatraLine3Type


    /*!
    \brief A C++ wrapper for the inter acinar dependency element

    */
    class RedAirBloodScatraLine3 : public Core::Elements::Element
    {
     public:
      //! @name Friends
      friend class RedAirwayImplInterface;
      friend class AirwayImpl<Core::FE::CellType::line3>;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      */
      RedAirBloodScatraLine3(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      RedAirBloodScatraLine3(const RedAirBloodScatraLine3& old);

      /*!
      \brief Deep copy this instance of Redacinus and return pointer to the copy

      The clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      Core::Elements::Element* clone() const override;

      /*!
      \brief Get shape type of element
      */
      Core::FE::CellType shape() const override;

      /*!
      \brief Return number of lines of this element
      */
      int num_line() const override
      {
        if (num_node() == 3)
        {
          return 2;
        }
        else
        {
          FOUR_C_THROW("Could not determine number of lines");
          return -1;
        }
      }

      /*!
      \brief Return number of surfaces of this element (always 1)
      */
      int num_surface() const override { return 1; }

      /*!
      \brief Return number of volumes of this element (always 1)
      */
      int num_volume() const override { return -1; }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      std::vector<Teuchos::RCP<Core::Elements::Element>> lines() override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int unique_par_object_id() const override
      {
        return RedAirBloodScatraLine3Type::instance().unique_par_object_id();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref pack and \ref unpack are used to communicate this element

      */
      void pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref pack and \ref unpack are used to communicate this element

      */
      void unpack(const std::vector<char>& data) override;


      //@}

      //! @name Acess methods


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual Core::Elements::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      int num_dof_per_node(const Core::Nodes::Node& node) const override { return 1; }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual Core::Elements::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      int num_dof_per_element() const override { return 0; }

      /*!
      \brief Print this element
      */
      void print(std::ostream& os) const override;

      RedAirBloodScatraLine3Type& element_type() const override
      {
        return RedAirBloodScatraLine3Type::instance();
      }

      //@}

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
      void vis_names(std::map<std::string, int>& names) override;

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
      bool vis_data(const std::string& name, std::vector<double>& data) override;


      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool read_element(const std::string& eletype, const std::string& distype,
          Input::LineDefinition* linedef) override;

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate inter acinar dependency element stiffness, mass, internal forces etc

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param elemat1 (out)  : matrix to be filled by element. If nullptr on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elemat2 (out)  : matrix to be filled by element. If nullptr on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elevec1 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      int evaluate(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;

      /*!
      \brief Evaluate a Neumann boundary condition

      An element derived from this class uses the evaluate_neumann method to receive commands
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
      int evaluate_neumann(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseMatrix* elemat1 = nullptr) override;

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a line Neumann condition on the acinus element

      \return 0 if successful, negative otherwise
      */
      virtual int evaluate_dirichlet(Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1);


      //@}


      //! get cross-sectional area
      //  double         getdata(){ return A_;}

      //! get youngs modulus of the wall
      //  double         getEw(){ return Ew_;}

      //! get youngs modulus of the air
      //  double         getEa(){ return Ea_;}

      //! get wall thickness
      //  double         gettw(){ return tw_;}

      //! set qin
      //  void           setqin (double Qin ){ qin_ = Qin;}

      //! set qout
      //  void           setqout(double Qout){ qout_ = Qout;}

      //! get element parameters
      void get_params(std::string name, double& var);

      //! get element parameters
      void get_params(std::string name, int& var);

     private:
      //! action parameters recognized by inter acinar dependency
      enum ActionType
      {
        none,
        calc_sys_matrix_rhs,
        calc_sys_matrix_rhs_iad,
        get_initial_state,
        set_bc,
        calc_flow_rates,
        get_coupled_values,
        get_junction_volume_mix,
        solve_scatra,
        calc_cfl,
        solve_junction_scatra,
        solve_blood_air_transport,
        update_scatra,
        update_elem12_scatra,
        eval_nodal_ess_vals,
        eval_PO2_from_concentration,
        calc_elem_volumes
      };

      // data
      std::map<std::string, double> elem_params_;

      // element tree generation
      int generation_;

      // internal calculation methods

      // don't want = operator
      RedAirBloodScatraLine3& operator=(const RedAirBloodScatraLine3& old);


      /// set number of gauss points to element shape default
      Core::FE::GaussRule1D get_optimal_gaussrule(const Core::FE::CellType& distype);

      /*!
       * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are
       * necessary \return boolean indicating higher order status
       */
      bool is_higher_order_element(const Core::FE::CellType distype  ///< discretization type
      ) const;

    };  // class RedAirBloodScatraLine3


  }  // namespace ELEMENTS
}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
