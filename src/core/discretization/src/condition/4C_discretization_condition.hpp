/*---------------------------------------------------------------------*/
/*! \file

\brief A condition of any kind

\level 1


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_DISCRETIZATION_CONDITION_HPP
#define FOUR_C_DISCRETIZATION_CONDITION_HPP


#include "4C_config.hpp"

#include "4C_io_input_parameter_container.hpp"
#include "4C_legacy_enum_definitions_conditions.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

#include <tuple>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  // forward declarations
  class Element;
  class Discretization;
}  // namespace DRT

namespace CORE::Conditions
{

  /*!
   * A condition is mainly used to realize boundary conditions. Parameters for the condition
   * are stored in a InputParameterContainer.
   * The condition can additionally store a discretization of the condition which is
   * driven by the Discretization class that is evaluating this condition.
   * The Discretization class is therefore a friend of the Condition and has access to
   * the protected methods dealing with the discretization of this condition.
   */
  class Condition
  {
   public:
    //! @name Enums and Friends

    /*!
    \brief Discretization is a friend of the condition to have access
           to the private methods that would otherwise have to be public.

    */
    friend class DRT::Discretization;

    //@}

    //! @name Constructors and destructors

    /*!
    \brief Standard Constructor

    The way a condition is treated later on depends on the type of the
    condition. E.g. Dirichlet conditions are treated differently from
    Neumann conditions. How they are treated is not described here but in
    DRT::Discretization.

    \note In case you might wonder where this condition class actually stores
          data necessary for the condition: This class implements IO::InputParameterContainer.

    \param id (in): a unique id for this condition
    \param type (in): type of the condition
    \param buildgeometry (in): flag indicating whether explicit condition geometry
                               (elements) have to be build
    \param gtype (in): type of geometric entity this condition lives on
    */
    Condition(const int id, const CORE::Conditions::ConditionType type, const bool buildgeometry,
        const CORE::Conditions::GeometryType gtype);

    /*!
    \brief Default constructor with type condition_none
    */
    Condition() = default;


    //@}

    //! @name Query methods

    /*!
    \brief Return condition id
    */
    [[nodiscard]] inline int Id() const { return id_; }

    /*!
    \brief Return vector of my global node ids
    */
    [[nodiscard]] const std::vector<int>* GetNodes() const { return &nodes_; }

    /*!
    \brief Set vector of my global node ids
    */
    void SetNodes(const std::vector<int>& nodes) { nodes_ = nodes; }

    /*!
      \brief Return if a node gid is contained in this condition
     */
    [[nodiscard]] bool ContainsNode(int ngid) const
    {
      const std::vector<int>* n = GetNodes();
      // CORE::Conditions::Condition nodes are ordered by design! So we can perform a binary
      // search here.
      return std::binary_search(n->begin(), n->end(), ngid);
    }

    /*!
    \brief Return flag indicating whether this condition needs to build a geometry
           description

    Some boundary conditions such as e.g. Neumann BCs need a geometry description
    to perform an integration on the boundary. Some BCs such as Dirichlet BCs
    don't need such a geometry description as it is sufficient to have access to
    the nodes only.<br>
    In case the condition needs to build elements describing the geometry of the
    condition the returned flag is true, otherwise its false;

    */
    [[nodiscard]] inline bool GeometryDescription() const { return buildgeometry_; }

    /*!
    \brief Return type of geometry this condition lives on

    The type of geometry this condition lives on determines what type of
    geometry description is build for this condition iff GeometryDescription()==true

    */
    [[nodiscard]] inline CORE::Conditions::GeometryType GType() const { return gtype_; }

    /*!
    \brief Print this Condition
    */
    void Print(std::ostream& os) const;

    /*!
    \brief Return type of condition
    */
    [[nodiscard]] inline CORE::Conditions::ConditionType Type() const { return type_; }

    /*!
    \brief Get a reference to the geometry description of the condition

    */
    std::map<int, Teuchos::RCP<DRT::Element>>& Geometry() { return *geometry_; }

    [[nodiscard]] const std::map<int, Teuchos::RCP<DRT::Element>>& Geometry() const
    {
      return *geometry_;
    }

    //! Access the container that stores the input parameters.
    [[nodiscard]] const IO::InputParameterContainer& parameters() const { return container_; }

    //! Access the container that stores the input parameters.
    IO::InputParameterContainer& parameters() { return container_; }

    /*!
    \brief Adjust IDs of associated elements in order to obtain global
    unique IDs within one condition type
    */
    void AdjustId(const int shift);

    /**
     * Create a copy of this object but do not copy the geometry.
     */
    [[nodiscard]] Teuchos::RCP<CORE::Conditions::Condition> copy_without_geometry() const;

    //! Comparison operator.
    friend bool operator<(const Condition& lhs, const Condition& rhs);

    //@}

   private:
    //! @name Construction methods
    /*!
    \brief Add a geometry description to the condition

    A geometry description can be added to the condition.
    In case the condition refers to lines, surfaces or volumes, a
    geometry description might be needed to properly evaluate the condition
    (e.g. in the case of Neumann conditions).
    Such a geometry description is build in \ref DRT::Discretization::BoundaryConditionsGeometry
    and then added to this Condition.
    The geometry description consists of elements that are capable to
    perform the necessary operations on the condition (e.g. integrate a Neumann BC
    along a line). The matching nodes are taken from the
    underlying discretization itself. Also, it is actually the Discretization class
    that drives this process, so do not add elements yourself to the condition, let
    the Discretization do it for you.

    \param geom (in): Map of elements describing the geometry.
                      A deep copy of the map is made and stored.
                      Normally though, these elements are a line, surface or
                      volume elements produced by and shared with the discretization.
                      Do not mess with their Teuchos::RCP!

    */
    void AddGeometry(Teuchos::RCP<std::map<int, Teuchos::RCP<DRT::Element>>> geom)
    {
      geometry_ = geom;
    }

    /*!
    \brief Delete a geometry description of the condition
    */
    void ClearGeometry() { geometry_ = Teuchos::null; }

    //@}

    Condition(const Condition& old) = default;

    Condition& operator=(const Condition& old) = default;

    //! Unique id of this condition, no second condition of the same type with same id may exist
    int id_{-1};

    //! global node ids
    std::vector<int> nodes_{};

    //! flag indicating whether this condition builds a geometry description or not
    bool buildgeometry_{};

    //! Type of this condition
    CORE::Conditions::ConditionType type_{};

    //! Type of geometry the condition lives on
    CORE::Conditions::GeometryType gtype_{};

    //! Geometry description of this condition
    Teuchos::RCP<std::map<int, Teuchos::RCP<DRT::Element>>> geometry_{};

    IO::InputParameterContainer container_;
  };  // class Condition

  inline bool operator<(const Condition& lhs, const Condition& rhs)
  {
    return std::tie(lhs.type_, lhs.id_) < std::tie(rhs.type_, rhs.id_);
  }

}  // namespace CORE::Conditions


//! << operator
std::ostream& operator<<(std::ostream& os, const CORE::Conditions::Condition& cond);


FOUR_C_NAMESPACE_CLOSE

#endif
