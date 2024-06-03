/*---------------------------------------------------------------------*/
/*! \file

\brief Main class to manage a discretization

\level 0

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_LIB_DISCRET_HPP
#define FOUR_C_LIB_DISCRET_HPP

#include "4C_config.hpp"

#include "4C_discretization_dofset_interface.hpp"
#include "4C_discretization_fem_general_element.hpp"
#include "4C_discretization_fem_general_node.hpp"
#include "4C_discretization_fem_general_shape_function_type.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_std_cxx20_ranges.hpp"

#include <Epetra_Comm.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_IntVector.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <functional>
#include <vector>

// forward declarations
namespace Teuchos
{
  class ParameterList;

  template <class T>
  class RCP;
}  // namespace Teuchos

FOUR_C_NAMESPACE_OPEN

class PostProblem;

namespace CORE::COMM
{
  class Communicators;
}  // namespace CORE::COMM

namespace CORE::LINALG
{
  class SparseOperator;
  class MapExtractor;
  class SerialDenseVector;
}  // namespace CORE::LINALG

namespace CORE::FE
{
  class AssembleStrategy;
}

namespace CORE::Dofsets
{
  class DofSetProxy;
}

namespace CORE::UTILS
{
  class FunctionManager;
}

namespace CORE::IO
{
  class DiscretizationWriter;
}  // namespace CORE::IO

namespace GLOBAL
{
  class Problem;
}

namespace DRT
{
  /*!
  \brief A class to manage a discretization in parallel

  A discretization describes a discretized physical domain including any
  boundary conditions. This is a huge concept that includes four different parts:

  - a mesh consisting of elements and nodes
  - a number of sets of degrees of freedom (dofsets)
  - conditions that point to nodes (nodal clouds) but really mean dofs
  - state vectors that belong to one of the dofsets

  <h3>Parallelization</h3>

  A huge point with the discretization is its parallel distribution. In normal
  circumstances each processor will hole a slice of the distribution. There will
  be some overlapping elements. This is known as ghosting as all elements belong
  to exactly one processor and may appear as ghosts (guests) on other processors
  as well.

  The parallel distribution is handled via Epetra_Map objects. There are
  row-maps that are unique without any overlap. There are column-maps that
  extend the row-maps with the overlap information. Elements, nodes and dofs are
  managed by using maps. That is all three items have a globally unique id,
  oftentimes called gid. Global ids are arbitrary numbers. Apart from the fact
  that these numbers are unique there comes no guarantee whatsoever with these
  numbers. You must never use the numerical value of a gid.

  In addition to the gids there are local ids, called lids. These are always
  consecutive numbers starting from 0. Local ids just count local items. See
  Epetra_Map for details.

  <h3>Initialization</h3>

  The initialization of a discretization is a major effort since the parallel
  distribution needs to be established. Normally the discretization is read from
  a dat-file with the help of DRT::DatFileReader and related classes. The
  discretization class comes with a bunch of helper method for its setup phase.

  <h3>Filled-State</h3>

  A discretization can be in its filled or non-filled state. After construction
  a discretization starts in its non-filled state. In the non-filled state there
  are no lids calculated yet. The discretization is not fit for any calculation.

  A call for fill_complete() fills the discretization, that is brings it to its
  filled state. Now all lids are available. See below for details.

  <h3>DofSets</h3>

  A discretization can have multiple dof sets. At construction time the first
  DofSet is created. Further dof sets can be added later on.

  The first DofSet is always special as it describes the set of unknowns of this
  discretization. These are the values that are to be calculated later
  on. Additional dof sets can describe additional sets of variables that
  supplement the calculation. In normal single field calculations only one
  DofSet is needed.

  For multifield calculations with volume coupled fully overlapping matching
  discretizations that share the same parallel distribution (that is different
  fields that are discretized with the same mesh -- a rather special case that
  turns out to be quite common) there is the DofSetProxy class that can be used
  to introduce the DofSet of the respective other field.

  <h3>Misc</h3>

  The \ref DRT::Discretization class supports the ostream& operator <<

  */
  class Discretization
  {
   public:
    /*!
    \brief Standard Constructor

    \param name: name of this discretization
    \param comm: Epetra comm object associated with this discretization
    \param n_dim: number of space dimensions of this discretization
    */
    Discretization(const std::string& name, Teuchos::RCP<Epetra_Comm> comm, unsigned int n_dim);

    /**
     * Virtual destructor.
     */
    virtual ~Discretization() = default;

    /**
     * The discretization is a heavy object that should not be copied (accidentally).
     */
    Discretization(const DRT::Discretization&) = delete;

    /**
     * The discretization is a heavy object that should not be copied (accidentally).
     */
    Discretization operator=(const Discretization&) = delete;

    /**
     * Moving a discretization is not supported.
     */
    Discretization(DRT::Discretization&&) = delete;

    /**
     * Moving a discretization is not supported.
     */
    Discretization operator=(Discretization&&) = delete;

    //! @name Query methods

    /*!
    \brief Get communicator associated with this class
    */
    [[nodiscard]] virtual const Epetra_Comm& Comm() const { return *comm_; }

    /*!
    \brief Get output writer for this discretization

    \warning This routine does not verify if a valid CORE::IO::DiscretizationWriter has
    been set. If not, this will cause a segmentation fault.
    */
    [[nodiscard]] virtual Teuchos::RCP<CORE::IO::DiscretizationWriter> Writer() const
    {
      return writer_;
    }

    /*!
    \brief Get flag indicating whether fill_complete() has been called
    */
    [[nodiscard]] virtual bool Filled() const { return filled_; }

    /*!
    \brief Get name of this discretization
    */
    [[nodiscard]] virtual const std::string& Name() const { return name_; }

    /*!
    \brief Get flag indicating whether degrees of freedom where assigned

    Degrees of freedom need to be assigned using assign_degrees_of_freedom()
    before any calculations using this discretization can be made
    */
    [[nodiscard]] virtual bool HaveDofs() const { return havedof_; }

    [[nodiscard]] virtual int NumDofSets() const { return static_cast<int>(dofsets_.size()); }
    /// @name Dof query methods for single dof set discretizations

    /*!
    \brief Get number of dofs for given element.

    Ask the current DofSet for the number of dofs of this element.
    There is a variable number of DofSets and one is currently selected.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param element (in)      : the element those number of dofs are requested
    */
    virtual int NumDof(const CORE::Elements::Element* element) const
    {
      FOUR_C_ASSERT(
          dofsets_.size() == 1, "Discretization %s expects just one dof set!", name_.c_str());
      return NumDof(0, element);
    }

    /*!
    \brief Get the gid of a dof for given element.

    Ask the current DofSet for the gid of the dof of this element.
    There is a variable number of DofSets and one is currently selected.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param element (in)      : the element
    \param dof (in)          : the element local dof number
    */
    virtual int Dof(const CORE::Elements::Element* element, const int dof) const
    {
      FOUR_C_ASSERT(
          dofsets_.size() == 1, "Discretization %s expects just one dof set!", name_.c_str());
      return Dof(0, element, dof);
    }

    /*!
    \brief Get the gid of all dofs of a element.

    Ask the current DofSet for the gids of the dofs of this element. The
    required vector is created and filled on the fly. So better keep it
    if you need more that one dof gid.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param element (in)      : the element
    */
    virtual std::vector<int> Dof(const CORE::Elements::Element* element) const
    {
      FOUR_C_ASSERT(
          dofsets_.size() == 1, "Discretization %s expects just one dof set!", name_.c_str());
      return Dof(0, element);
    }

    /*!
    \brief Get the gid of all dofs of a node.

    Ask the current DofSet for the gids of the dofs of this node. The
    required vector is created and filled on the fly. So better keep it
    if you need more than one dof gid.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param element (in)   : the element operated on
    \param node (in)      : the node
    \param lm (in/out)    : lm vector the dofs are appended to
    */
    virtual void Dof(const CORE::Elements::Element* element, const CORE::Nodes::Node* node,
        std::vector<int>& lm) const
    {
      FOUR_C_ASSERT(
          dofsets_.size() == 1, "Discretization %s expects just one dof set!", name_.c_str());
      Dof(0, element, node, lm);
    }

    /*!
    \brief Get the gid of all dofs of a element.

    Ask the current DofSet for the gids of the dofs of this element. The
    required vector is created and filled on the fly. So better keep it
    if you need more that one dof gid.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param element (in)      : the element
    \param lm (in/out)    : lm vector the dofs are appended to
    */
    virtual void Dof(const CORE::Elements::Element* element, std::vector<int>& lm) const
    {
      FOUR_C_ASSERT(
          dofsets_.size() == 1, "Discretization %s expects just one dof set!", name_.c_str());
      Dof(0, element, lm);
    }

    /*!
    \brief Get the gid of all dofs of a node.

    Ask the current DofSet for the gids of the dofs of this node. The
    required vector is created and filled on the fly. So better keep it
    if you need more than one dof gid.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param node (in)      : the node
    \param startindex (in): first index of vector at which will be written to end
    \param lm (in/out)    : lm vector the dofs are appended to
    */
    virtual void Dof(
        const CORE::Nodes::Node* node, const unsigned startindex, std::vector<int>& lm) const
    {
      FOUR_C_ASSERT(
          dofsets_.size() == 1, "Discretization %s expects just one dof set!", name_.c_str());
      Dof(0, node, startindex, lm);
    }

    /*!
    \brief Replace the dofset associated with the discretisation by a new dofset.
                   Sets havedof_ to false.

    @param[in] newdofset New DofSet to be used in this discretization
    @param[in] replaceinstatdofsets Replace also in static dofsets?
    */
    virtual void ReplaceDofSet(
        Teuchos::RCP<CORE::Dofsets::DofSetInterface> newdofset, bool replaceinstatdofsets = false);

    /*!
    \brief Get master to slave coupling in case of periodic boundary conditions
    */
    std::map<int, std::vector<int>>* get_all_pbc_coupled_col_nodes();

    /*!
    \brief Get slave to master connectivity in case of periodic boundary conditions
    */
    virtual Teuchos::RCP<std::map<int, int>> get_pbc_slave_to_master_node_connectivity();

    //@}

    /// @name Dof query methods for multi dof set discretizations

    /*!
    \brief Get number of dofs for given node at the default dofset (0).

    Ask the current DofSet for the number of dofs of this node.
    There is a variable number of DofSets and one is currently selected.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param node (in)      : the node those number of dofs are requested
    */
    int NumDof(const CORE::Nodes::Node* node) const
    {
      FOUR_C_ASSERT(
          dofsets_.size() == 1, "Discretization %s expects just one dof set!", name_.c_str());
      return NumDof(0, node);
    }

    /*!
    \brief Get number of dofs for given node.

    Ask the current DofSet for the number of dofs of this node.
    There is a variable number of DofSets and one is currently selected.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param nds (in)       : number of dofset
    \param node (in)      : the node those number of dofs are requested
    */
    virtual int NumDof(unsigned nds, const CORE::Nodes::Node* node) const
    {
      FOUR_C_ASSERT(
          nds < dofsets_.size(), "undefined dof set found in discretization %s!", name_.c_str());
      FOUR_C_ASSERT(havedof_, "no dofs assigned in discretization %s!", name_.c_str());
      return dofsets_[nds]->NumDof(node);
    }

    /*!
    \brief Get number of dofs for given element.

    Ask the current DofSet for the number of dofs of this element.
    There is a variable number of DofSets and one is currently selected.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param nds (in)          : number of dofset
    \param element (in)      : the element those number of dofs are requested
    */
    virtual int NumDof(unsigned nds, const CORE::Elements::Element* element) const
    {
      FOUR_C_ASSERT(
          nds < dofsets_.size(), "undefined dof set found in discretization %s!", name_.c_str());
      FOUR_C_ASSERT(havedof_, "no dofs assigned in discretization %s!", name_.c_str());
      return dofsets_[nds]->NumDof(element);
    }

    /** \brief Get number of standard (w/o enrichment) dofs for given node.
     *
     *  In the default case, we can use the NumDof() routine. This method
     *  gives you only for enriched nodes a different result!
     *  ( see XFEM::DiscretizationXFEM )
     *
     *  \param nds (in)       : number of dofset
     *  \param node (in)      : the node those number of DoF's are requested
     *
     *  \author hiermeier \date 10/16 */
    virtual int NumStandardDof(const unsigned& nds, const CORE::Nodes::Node* node) const
    {
      return NumDof(nds, node);
    }

    /*!
    \brief Get the gid of a dof for given node at the default dofset (0).

    Ask the current DofSet for the gid of the dof of this node.
    There is a variable number of DofSets and one is currently selected.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param node (in)      : the node
    \param dof (in)       : the node local dof number
    */
    int Dof(const CORE::Nodes::Node* node, const int ldof) const
    {
      FOUR_C_ASSERT(
          NumDofSets() == 1, "Discretization %s expects just one dof set!", name_.c_str());
      return Dof(0, node, ldof);
    }

    /*!
    \brief Get the gid of a dof for given node.

    Ask the current DofSet for the gid of the dof of this node.
    There is a variable number of DofSets and one is currently selected.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param nds (in)       : number of dofset
    \param node (in)      : the node
    \param dof (in)       : the node local dof number
    */
    virtual int Dof(unsigned nds, const CORE::Nodes::Node* node, const int dof) const
    {
      FOUR_C_ASSERT(
          nds < dofsets_.size(), "undefined dof set found in discretization %s!", name_.c_str());
      FOUR_C_ASSERT(havedof_, "no dofs assigned in discretization %s!", name_.c_str());
      return dofsets_[nds]->Dof(node, dof);
    }

    /*!
    \brief Get the gid of a dof for given element.

    Ask the current DofSet for the gid of the dof of this element.
    There is a variable number of DofSets and one is currently selected.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param nds (in)          : number of dofset
    \param element (in)      : the element
    \param dof (in)          : the element local dof number
    */
    virtual int Dof(unsigned nds, const CORE::Elements::Element* element, const int dof) const
    {
      FOUR_C_ASSERT(
          nds < dofsets_.size(), "undefined dof set found in discretization %s!", name_.c_str());
      FOUR_C_ASSERT(havedof_, "no dofs assigned in discretization %s!", name_.c_str());
      return dofsets_[nds]->Dof(element, dof);
    }

    /*!
    \brief Get the gid of all dofs of a node at the default dofset (0).

    Ask the current DofSet for the gids of the dofs of this node. The
    required vector is created and filled on the fly. So better keep it
    if you need more than one dof gid.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param node (in)      : the node
    */
    std::vector<int> Dof(const CORE::Nodes::Node* node) const
    {
      FOUR_C_ASSERT(
          NumDofSets() == 1, "Discretization %s expects just one dof set!", name_.c_str());
      return Dof(0, node);
    }

    /*!
    \brief Get the gid of all dofs of a node.

    Ask the current DofSet for the gids of the dofs of this node. The
    required vector is created and filled on the fly. So better keep it
    if you need more than one dof gid.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param nds (in)       : number of dofset
    \param node (in)      : the node
    */
    virtual std::vector<int> Dof(unsigned nds, const CORE::Nodes::Node* node) const
    {
      FOUR_C_ASSERT(
          nds < dofsets_.size(), "undefined dof set found in discretization %s!", name_.c_str());
      FOUR_C_ASSERT(havedof_, "no dofs assigned in discretization %s!", name_.c_str());
      return dofsets_[nds]->Dof(node);
    }

    /*!
    \brief Get the gid of all dofs of a node.

    Ask the current DofSet for the gids of the dofs of this node. The
    required vector is created and filled on the fly. So better keep it
    if you need more than one dof gid.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))

    Additional input nodal dof set: If the node contains more than one set of dofs, which can be
    evaluated, the number of the set needs to be given. Currently only the case for XFEM.

    \param dof         (out): vector of dof gids (to be filled)
    \param node        (in) : the node
    \param nds         (in) : number of dofset
    \param nodaldofset (in) : number of nodal dofset
    \param element     (in) : the element (optionally)
    */
    virtual void Dof(std::vector<int>& dof, const CORE::Nodes::Node* node, unsigned nds,
        unsigned nodaldofset, const CORE::Elements::Element* element = nullptr) const
    {
      FOUR_C_ASSERT(
          nds < dofsets_.size(), "undefined dof set found in discretization %s!", name_.c_str());
      FOUR_C_ASSERT(havedof_, "no dofs assigned in discretization %s!", name_.c_str());
      return dofsets_[nds]->Dof(dof, node, nodaldofset);
    }

    /*!
    \brief Get the GID of all dofs of a element.

    Ask the current DofSet for the gids of the dofs of this element. The
    required vector is created and filled on the fly. So better keep it
    if you need more that one dof gid.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param nds     (in) : number of dofset
    \param element (in) : the element
    */
    virtual std::vector<int> Dof(unsigned nds, const CORE::Elements::Element* element) const
    {
      FOUR_C_ASSERT(
          nds < dofsets_.size(), "undefined dof set found in discretization %s!", name_.c_str());
      FOUR_C_ASSERT(havedof_, "no dofs assigned in discretization %s!", name_.c_str());
      return dofsets_[nds]->Dof(element);
    }

    /*!
    \brief Get the gid of all dofs of a node at the default dofset (0).

    Ask the current DofSet for the gids of the dofs of this node. The
    required vector is created and filled on the fly. So better keep it
    if you need more than one dof gid.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param node (in)      : the node
    \param lm (in/out)    : lm vector the dofs are appended to
    */
    void Dof(const CORE::Nodes::Node* node, std::vector<int>& lm) const
    {
      FOUR_C_ASSERT(
          NumDofSets() == 1, "Discretization %s expects just one dof set!", name_.c_str());
      Dof((unsigned)0, node, lm);
    }

    /*!
    \brief Get the gid of all dofs of a node.

    Ask the current DofSet for the gids of the dofs of this node. The
    required vector is created and filled on the fly. So better keep it
    if you need more than one dof gid.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param nds (in)       : number of dofset
    \param node (in)      : the node
    \param lm (in/out)    : lm vector the dofs are appended to
    */
    virtual void Dof(unsigned nds, const CORE::Nodes::Node* node, std::vector<int>& lm) const
    {
      FOUR_C_ASSERT(
          nds < dofsets_.size(), "undefined dof set found in discretization %s!", name_.c_str());
      FOUR_C_ASSERT(havedof_, "no dofs assigned in discretization %s!", name_.c_str());
      dofsets_[nds]->Dof(node, lm);
    }

    /*!
    \brief Get the gid of all dofs of a node.

    Ask the current DofSet for the gids of the dofs of this node. The
    required vector is created and filled on the fly. So better keep it
    if you need more than one dof gid.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param nds (in)       : number of dofset
    \param element (in)   : the element operated on
    \param node (in)      : the node
    \param lm (in/out)    : lm vector the dofs are appended to
    */
    virtual void Dof(unsigned nds, const CORE::Elements::Element* element,
        const CORE::Nodes::Node* node, std::vector<int>& lm) const
    {
      FOUR_C_ASSERT(
          nds < dofsets_.size(), "undefined dof set found in discretization %s!", name_.c_str());
      FOUR_C_ASSERT(havedof_, "no dofs assigned in discretization %s!", name_.c_str());
      dofsets_[nds]->Dof(element, node, lm);
    }

    /*!
    \brief Get the gid of all dofs of a element.

    Ask the current DofSet for the gids of the dofs of this element. The
    required vector is created and filled on the fly. So better keep it
    if you need more that one dof gid.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param nds (in)       : number of dofset
    \param element (in)      : the element
    \param lm (in/out)    : lm vector the dofs are appended to
    */
    virtual void Dof(
        unsigned nds, const CORE::Elements::Element* element, std::vector<int>& lm) const
    {
      FOUR_C_ASSERT(
          nds < dofsets_.size(), "undefined dof set found in discretization %s!", name_.c_str());
      FOUR_C_ASSERT(havedof_, "no dofs assigned in discretization %s!", name_.c_str());
      dofsets_[nds]->Dof(element, lm);
    }

    /*!
    \brief Get the gid of all dofs of a node.

    Ask the current DofSet for the gids of the dofs of this node. The
    required vector is created and filled on the fly. So better keep it
    if you need more than one dof gid.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))
    \param nds (in)       : number of dofset
    \param node (in)      : the node
    \param startindex (in): first index of vector at which will be written to end
    \param lm (in/out)    : lm vector the dofs are appended to
    */
    virtual void Dof(unsigned nds, const CORE::Nodes::Node* node, const unsigned startindex,
        std::vector<int>& lm) const
    {
      FOUR_C_ASSERT(
          nds < dofsets_.size(), "undefined dof set found in discretization %s!", name_.c_str());
      FOUR_C_ASSERT(havedof_, "no dofs assigned in discretization %s!", name_.c_str());
      dofsets_[nds]->Dof(node, startindex, lm);
    }

    /*!
    \brief Replace the dofset associated with the discretisation by a new dofset.
                   Sets havedof_ to false.
    */
    virtual void ReplaceDofSet(unsigned nds, Teuchos::RCP<CORE::Dofsets::DofSetInterface> newdofset,
        bool replaceinstatdofsets = false);

    /*!
    \brief Add a new dofset to the discretisation.

    Sets havedof_ to false only if the new dofset is not properly filled yet.
    */
    virtual int AddDofSet(Teuchos::RCP<CORE::Dofsets::DofSetInterface> newdofset);

    /*!
    \brief Get proxy to dof set.
    */
    virtual Teuchos::RCP<CORE::Dofsets::DofSetInterface> GetDofSetProxy(int nds = 0);

    /*!
    \brief Get degree of freedom row map (Filled()==true prerequisite)

    Return ptr to degree of freedom row distribution map of this discretization.
    If it does not exist yet, build it.

    - Filled()==true prerequisite
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))

    */
    [[nodiscard]] virtual const Epetra_Map* dof_row_map(unsigned nds = 0) const;

    /*!
    \brief Get degree of freedom column map (Filled()==true prerequisite)

    Return ptr to degree of freedom column distribution map of this discretization.
    If it does not exist yet, build it.

    - Filled()==true prerequisite
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))

    */
    virtual const Epetra_Map* DofColMap(unsigned nds = 0) const;

    //@}

    /*!
    \brief Print this discretization to os (Filled()==true NOT prerequisite)
           (ostream << also supported)

    \note This is a collective call
    */
    virtual void Print(std::ostream& os) const;


    /*!
    \brief Get map associated with the distribution of the ownership of nodes
           (Filled()==true prerequisite)

    This map includes all nodes stored on this proc and also owned by this proc.
    This map is non-ambiguous, meaning that it is a non-overlapping map.

    \return nullptr if Filled() is false. A call to fill_complete() is a prerequisite.
    */
    [[nodiscard]] virtual const Epetra_Map* NodeRowMap() const;

    /*!
    \brief Get map associated with the distribution of nodes including ghosted nodes
           (Filled()==true prerequisite)

    This map includes all nodes stored on this proc including any ghosted nodes
    This map is ambiguous, meaning that it is an overlapping map

    \return nullptr if Filled() is false. A call to fill_complete() is a prerequisite.
    */
    virtual const Epetra_Map* NodeColMap() const;
    /*!
    \brief Get map associated with the distribution of the ownership of elements
           (Filled()==true prerequisite)

    This map includes all elements stored on this proc and also owned by this proc.
    This map is non-ambiguous, meaning that it is a non-overlapping map.

    \return nullptr if Filled() is false. A call to fill_complete() is a prerequisite.
    */
    [[nodiscard]] virtual const Epetra_Map* ElementRowMap() const;
    /*!
    \brief Get map associated with the distribution of elements including ghosted elements
           (Filled()==true prerequisite)

    This map includes all elements stored on this proc including any ghosted elements
    This map is ambiguous, meaning that it is an overlapping map

    \return nullptr if Filled() is false. A call to fill_complete() is a prerequisite.
    */
    virtual const Epetra_Map* ElementColMap() const;

    /*!
    \brief Get global number of elements (true number of total elements)
           (Filled()==true prerequisite)

    This is a collective call
    */
    virtual int NumGlobalElements() const;

    /*!
    \brief Get processor local number of elements owned by this processor
           (Filled()==true prerequisite)
    */
    [[nodiscard]] virtual int NumMyRowElements() const;

    /*!
    \brief Get processor local number of elements including ghost elements
           (Filled()==true NOT prerequisite)
    */
    virtual int NumMyColElements() const;

    /*!
    \brief Get global number of nodes (true number of total nodes without ghosting)
           (Filled()==true prerequisite)
    */
    virtual int NumGlobalNodes() const;

    /*!
    \brief Get processor local number of nodes owned by this processor
           (Filled()==true prerequisite)
    */
    [[nodiscard]] virtual int NumMyRowNodes() const;

    /*!
    \brief Get processor local number of nodes including ghost nodes
           (Filled()==true NOT prerequisite)
    */
    virtual int NumMyColNodes() const;

    /*!
    \brief Query whether an element with global id @p gid is stored (i.e. owned or ghosted) on this
    proc.

    \note: this query does not tell, if the element is owned by this proc. Use
    CORE::Elements::Element.Owner() for this.

    */
    virtual bool HaveGlobalElement(int gid) const;

    /*!
    \brief Get the element with global id gid (Filled()==true NOT prerequisite)

    Returns the element with global row id gid if element is on this proc.
    Will return row or column element, ghosted or not.
    This is an individual call. Will test for existence of element in
    DEBUG version and throw an error message if not. Will crash in non-DEBUG
    version if element does not exist on calling processor.

    \return Address of element if element is owned by calling proc, returns nullptr
            otherwise
    */
    [[nodiscard]] virtual CORE::Elements::Element* gElement(int gid) const;

    /*!
    \brief Get the element with local row id lid (Filled()==true prerequisite)

    Returns the element with local row index lid.
    Will not return any ghosted element.
    This is an individual call and Filled()=true is a prerequisite

    \return Address of element if element is owned by calling proc
    */
    [[nodiscard]] virtual CORE::Elements::Element* lRowElement(int lid) const
    {
      FOUR_C_ASSERT(Filled(), "Discretization %s not Filled()!", name_.c_str());
      return elerowptr_[lid];
    }

    /*!
    \brief Get the element with local column id lid (Filled()==true prerequisite)

    Returns the element with local column index lid.
    Will also return any ghosted element.
    This is an individual call and Filled()=true is a prerequisite

    \return Address of element if element is stored by calling proc
    */
    virtual CORE::Elements::Element* lColElement(int lid) const
    {
      FOUR_C_ASSERT(Filled(), "Discretization %s not Filled()!", name_.c_str());
      return elecolptr_[lid];
    }

    /**
     * This function is useful for range based for-loops over all row elements.
     *
     * \code
     *      for (CORE::Elements::Element* actele : MyRowElementRange()) {}
     * \endcode
     *
     * \return A range of all local row elements.
     */
    [[nodiscard]] auto MyRowElementRange() const
    {
      FOUR_C_ASSERT(Filled(), "Discretization %s not Filled()!", name_.c_str());
      return std_20::ranges::views::all(elerowptr_);
    }

    /**
     * This function is equivalent to MyRowElementRange(), but applied to column elements
     */
    [[nodiscard]] auto MyColElementRange() const
    {
      FOUR_C_ASSERT(Filled(), "Discretization %s not Filled()!", name_.c_str());
      return std_20::ranges::views::all(elecolptr_);
    }

    /*!
    \brief Query whether a node with global id @p gid is stored (i.e. owned or ghosted) on this proc
           (Filled()==true NOT prerequisite)

    \note: this query does not tell, if the node is owned by this proc. Use
    CORE::Nodes::Node.Owner() for this.

    */
    [[nodiscard]] virtual bool HaveGlobalNode(int gid) const;

    /*!
    \brief Get the node with global row id gid (Filled()==true NOT prerequisite)

    Returns the node with global row id gid if node is on this proc.
    Will return row or column node, ghosted or  not.
    This is an individual call

    \return Address of node if node is stored on calling proc
    */
    [[nodiscard]] virtual CORE::Nodes::Node* gNode(int gid) const;

    /*!
    \brief Get the node with local row id lid (Filled()==true prerequisite)

    Returns the node with local row index lid.
    Will not return any ghosted node.
    This is an individual call and Filled()=true is a prerequisite

    \return Address of node if node is owned and stored by calling proc
    */
    [[nodiscard]] virtual CORE::Nodes::Node* lRowNode(int lid) const
    {
      FOUR_C_ASSERT(Filled(), "Discretization %s not Filled()!", name_.c_str());
      return noderowptr_[lid];
    }

    /*!
    \brief Get the node with local column id lid (Filled()==true prerequisite)

    Returns the node with local column index lid.
    Will return any node stored on this proc.
    This is an individual call and Filled()=true is a prerequisite

    \return Address of node if node is stored by calling proc
    */
    virtual CORE::Nodes::Node* lColNode(int lid) const
    {
      FOUR_C_ASSERT(Filled(), "Discretization %s not Filled()!", name_.c_str());
      return nodecolptr_[lid];
    }


    [[nodiscard]] auto MyRowNodeRange() const
    {
      FOUR_C_ASSERT(Filled(), "Discretization %s not Filled()!", name_.c_str());
      return std_20::ranges::views::all(noderowptr_);
    }


    [[nodiscard]] auto MyColNodeRange() const
    {
      FOUR_C_ASSERT(Filled(), "Discretization %s not Filled()!", name_.c_str());
      return std_20::ranges::views::all(nodecolptr_);
    }

    unsigned int n_dim() const { return n_dim_; }

    //@}

    //! @name Construction methods

    /*!
    \brief Set an Epetra_Comm object which can be either
           Epetra_SerialComm or Epetra_MpiComm (Filled()==true NOT prerequisite)

    \note Sets Filled()=false
    */
    virtual void SetComm(Teuchos::RCP<Epetra_Comm> comm)
    {
      filled_ = false;
      comm_ = comm;
    }

    /*!
    \brief Set a DiscretizationWriter
    */
    virtual void SetWriter(Teuchos::RCP<CORE::IO::DiscretizationWriter> writer)
    {
      writer_ = writer;
    }

    /*!
    \brief Add an element to the discretization (Filled()==true NOT prerequisite)

    The discretization takes ownership of the added element.
    Note that if an element with the same Id() exists, it will be
    deleted and replaced by the new one;
    note furthermore that in this method Reset() is called only on the processor where
    an element has been added actually; however, as such a modification may affect the
    discretization as a whole it may be required to call Reset() on each processor
    subsequently; to do so, please, call CheckFilledGlobally() if required

    \note Sets Filled()=false
    */
    virtual void add_element(Teuchos::RCP<CORE::Elements::Element> ele);

    /*!
    \brief Add a node to the discretization  (Filled()==true NOT prerequisite)

    The discretization takes ownership of the added node.
    Note that if a node with the same Id() exists, it will be
    deleted and replaced by the new one;
    note furthermore that in this method Reset() is called only on the processor where
    a node has been added actually; however, as such a modification may affect the
    discretization as a whole it may be required to call Reset() on each processor
    subsequently; to do so, please, call CheckFilledGlobally() if required

    \note Sets Filled()=false
    */
    virtual void AddNode(Teuchos::RCP<CORE::Nodes::Node> node);

    /*!
    \brief Delete an node from the discretization (Filled()==true NOT prerequisite)

    Delete an node from the discretization. Node can either be ghosted or not.
    If node on calling processor is not found nothing is done and false is returned.
    Note that this is not a fatal error and no message will be posted;
    note furthermore that in this method Reset() is called only on the processor where
    an node has been deleted actually; however, as such a modification may affect the
    discretization as a whole it may be required to call Reset() on each processor
    subsequently; to do so, please, call CheckFilledGlobally() if required

    \return Returns true upon successful deletion or node, returns false if node
            was not found on calling proc.

    \note Sets Filled()=false and calls Reset() upon discretization.
    */
    virtual bool DeleteNode(Teuchos::RCP<CORE::Nodes::Node> node);

    /*!
    \brief Delete an node with global id gid from the discretization
           (Filled()==true NOT prerequisite)

    Delete an node from the discretization. Node can either be ghosted or not.
    If node on calling processor is not found nothing is done and false is returned.
    Note that this is not a fatal error and no message will be posted;
    note furthermore that in this method Reset() is called only on the processor where
    an node has been deleted actually; however, as such a modification may affect the
    discretization as a whole it may be required to call Reset() on each processor
    subsequently; to do so, please, call CheckFilledGlobally() if required

    \return Returns true upon successful deletion or node, returns false if node
            was not found on calling proc.

    \note Sets Filled()=false and calls Reset() upon discretization.
    */
    virtual bool DeleteNode(const int gid);

    /*!
    \brief Removes all nodes from discretization (Filled()==true NOT prerequisite)

    \return Returns true upon successful deletion of all nodes

    \note CheckFilledGlobally() is called to make sure all processors are informed.
    */
    virtual bool DeleteNodes();

    /*!
    \brief Removes all elements from discretization (Filled()==true NOT prerequisite)

    \return Returns true upon successful deletion of all elements

    \note CheckFilledGlobally() is called to make sure all processors are informed.
    */
    virtual bool DeleteElements();

    /*!
    \brief Delete an element from the discretization (Filled()==true NOT prerequisite)

    Delete an element from the discretization. Element can either be ghosted or not.
    If element on calling processor is not found nothing is done and false is returned.
    Note that this is not a fatal error and no message will be posted;
    note furthermore that in this method Reset() is called only on the processor where
    an element has been deleted actually; however, as such a modification may affect the
    discretization as a whole it may be required to call Reset() on each processor
    subsequently; to do so, please, call CheckFilledGlobally() if required

    \return Returns true upon successful deletion or element, returns false if element
            was not found on calling proc.

    \note Sets Filled()=false and calls Reset() upon discretization.
    */
    virtual bool DeleteElement(Teuchos::RCP<CORE::Elements::Element> ele);

    /*!
    \brief Delete an element with global id gid from the discretization
           (Filled()==true NOT prerequisite)

    Delete an element from the discretization. Element can either be ghosted or not.
    If element on calling processor is not found nothing is done and false is returned.
    Note that this is not a fatal error and no message will be posted;
    note furthermore that in this method Reset() is called only on the processor where
    an element has been deleted actually; however, as such a modification may affect the
    discretization as a whole it may be required to call Reset() on each processor
    subsequently; to do so, please, call CheckFilledGlobally() if required

    \return Returns true upon successful deletion or element, returns false if element
            was not found on calling proc.

    \note Sets Filled()=false and calls Reset() upon discretization.
    */
    virtual bool DeleteElement(const int gid);

    /*!
    \brief Removes all nodes and elements from discretization (Filled()==true NOT prerequisite)

    \return Returns true upon successful deletion of all nodes and elements

    \note CheckFilledGlobally() is called to make sure all processors are informed.
    */
    virtual bool ClearDiscret();

    /*!
    \brief Complete construction of a discretization  (Filled()==true NOT prerequisite)

    After adding or deleting nodes or elements or redistributing them in parallel,
    or adding/deleting boundary conditions, this method has to be called to (re)construct
    pointer topologies.<br>
    It builds in this order:<br>
    - row map of nodes
    - column map of nodes
    - row map of elements
    - column map of elements
    - pointers from elements to nodes
    - pointers from nodes to elements
    - assigns degrees of freedoms
    - map of element register classes
    - calls all element register initialize methods
    - build geometries of all Dirichlet and Neumann boundary conditions

    \param assigndegreesoffreedom (in) : if true, resets existing dofsets and performs
                                         assigning of degrees of freedoms to nodes and
                                         elements.
    \param initelements (in) : if true, build element register classes and call Initialize()
                               on each type of finite element present
    \param doboundaryconditions (in) : if true, build geometry of boundary conditions
                                       present.

    \note In order to receive a fully functional discretization, this method must be called
          with all parameters set to true (the default). The parameters though can be
          used to turn off specific tasks to allow for more flexibility in the
          construction of a discretization, where it is known that this method will
          be called more than once.

    \note Sets Filled()=true
    */
    virtual int fill_complete(bool assigndegreesoffreedom = true, bool initelements = true,
        bool doboundaryconditions = true);

    /*!
    \brief Synchronize filled_ flag on all processors

    Modifying the discretization on some processors only (e.g. using add_element, AddNode,
    DeleteElement) may amount to filled_ == false on some processors only, but not on each
    processor; as one may consider the discretization as incomplete already if filled_ == false on
    at least one processor, it is recommended to set filled_ == false on each processor as soon as
    it occurs on at least one processor and to call Reset() subsequently on each processor; to do so
    just call CheckFilledGlobally
    */
    virtual void CheckFilledGlobally();


    //@}

    //! @name Boundary condition construction methods

    /*!
    \brief Set a condition with a certain name (Filled()==false on exit)

    Store a condition with a certain name in the discretization. The name need not
    be unique, meaning multiple conditions with the same name can be stored.
    Conditions can then be accessed with the GetCondition methods.

    \note Conditions attached to the discretization have to be
          completely redundant meaning that nodal cloud in the
          condition is the same on each processor and spans all
          nodes that hold this condition.<br>
          Also, setting a condition to the discretization sets the Filled() flag
          to false.

    \note

    \param name : Name of condition
    \param cond : The Condition class

    \warning If a condition with the exact same name already exists, it will
             NOT be overwritten but stored twice in the discretization

    */
    virtual void SetCondition(
        const std::string& name, Teuchos::RCP<CORE::Conditions::Condition> cond);

    /*!
    \brief Replace a condition with a certain name (Filled()==false on exit)

    First we look up if there are already conditions with the given name, and if
    yes we delete them. Afterwards we set the new conditions with the given name.

    \param name : name of the condition
    \param cond : vector of new conditions

    \author hiermeier \date 11/16 */
    void ReplaceConditions(const std::string& name,
        const std::vector<Teuchos::RCP<CORE::Conditions::Condition>>& conds);

    /*!
    \brief Get all conditions with a certain name

    Get all conditions with a certain name. A vector of ptrs to all conditions
    with name name is returned in out. The number of conditions found with name
    name is out.size(). out.size() is 0 if no condition with that name is found.

    \note Conditions attached to the discretization have to be
          completely redundant meaning that nodal cloud in the
          condition is the same on each processor and spans all
          nodes that hold this condition

    \param name (in): Name of condition
    \param out  (out): vector of pointers to all conditions with that name

    \return Returns out.size()=0 if condition with that name does not exist
    */
    virtual void GetCondition(
        const std::string& name, std::vector<CORE::Conditions::Condition*>& out) const;

    virtual void GetCondition(
        const std::string& name, std::vector<Teuchos::RCP<CORE::Conditions::Condition>>& out) const;

    /*! \brief Get a condition with a certain name

    Returns the first condition with name name found in the multimap.
    If multiple conditions with the same name exist, the first condition is
    returned and behaviour is therefore non-deterministic. This method should
    therefore only be used in cases where the user is sure that name is unique.

    \note Conditions attached to the discretization have to be
          completely redundant meaning that nodal cloud in the
          condition is the same on each processor and spans all
          nodes that hold this condition

    \param name (in): Name of condition

    \return Returns nullptr if condition with that name does not exist
    */
    [[nodiscard]] virtual CORE::Conditions::Condition* GetCondition(const std::string& name) const;

    /// return all condition names defined in this discretization
    virtual void GetConditionNames(std::vector<std::string>& names) const;

    /// return all conditions defined on this discretization
    virtual std::multimap<std::string, Teuchos::RCP<CORE::Conditions::Condition>>&
    GetAllConditions()
    {
      return condition_;
    }

    //@}

    //! @name Parallel (re)distribution

    /*!
    \brief Redistribute the discretization according to provided maps
           (Filled()==NOT true prerequisite)

    Steps taken in this method are:<br>
    - build element maps (row and column)
    - do export of row/col nodes and row/col elements
    - call Fillcomplete(assigndegreesoffreedom,initelements,doboundaryconditions)

    \param noderowmap (in): new node map the discretization shall have on exit
    \param nodecolmap (in): new column map the discretization shall have on exit

    \param assigndegreesoffreedom (in) : if true, resets existing dofsets and performs
                                         assigning of degrees of freedoms to nodes and
                                         elements.
    \param initelements           (in) : if true, build element register classes and
                                         call Initialize() on each type of finite element
                                         present
    \param doboundaryconditions   (in) : if true, build geometry of boundary conditions
                                         present.

    \note Filled()==true is a prerequisite, Filled()==true on exit
    */
    virtual void Redistribute(const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap,
        bool assigndegreesoffreedom = true, bool initelements = true,
        bool doboundaryconditions = true, bool killdofs = true, bool killcond = true);

    /*!
    \brief Redistribute the discretization according to provided maps
           (Filled()==NOT true prerequisite)

    \param noderowmap (in): new node map the discretization shall have on exit
    \param nodecolmap (in): new column map the discretization shall have on exit
    \param elerowmap  (in): new element row map the discretization shall have on exit
    \param elecolmap  (in): new element column map the discretization shall have on exit

    \param assigndegreesoffreedom (in) : if true, resets existing dofsets and performs
                                         assigning of degrees of freedoms to nodes and
                                         elements.
    \param initelements           (in) : if true, build element register classes and
                                         call Initialize() on each type of finite element
                                         present
    \param doboundaryconditions   (in) : if true, build geometry of boundary conditions
                                         present.

    \note Filled()==true is a prerequisite, Filled()==true on exit
    */
    virtual void Redistribute(const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap,
        const Epetra_Map& elerowmap, const Epetra_Map& elecolmap,
        bool assigndegreesoffreedom = true, bool initelements = true,
        bool doboundaryconditions = true, bool killdofs = true, bool killcond = true);

    /*!
      \brief Ghost elements on processors according to provided element column map

      Steps taken in this method are:<br>
      - build element maps (row and column)
      - do export of row/col nodes and row/col elements
      - call Fillcomplete(assigndegreesoffreedom,initelements,doboundaryconditions)

      \param elecolmap  (in): new element column map the discretization shall have on exit
      \param checkghosting (in): additional check can be performed

      */
    virtual void ExtendedGhosting(const Epetra_Map& elecolmap, bool assigndegreesoffreedom = true,
        bool initelements = true, bool doboundaryconditions = true, bool checkghosting = true);

    // Setup ghosting if you have a proper distribution of rownodes and elements
    // only use this if you are sure that the parallel distribution is ok!
    virtual void SetupGhosting(
        bool assigndegreesoffreedom, bool initelements, bool doboundaryconditions);

    //-----------------------------------------------------------------------------------
    // this is temporary to hide away this method SetupGhosting and only give access to the
    // class PostProblem. It is supposed to go away and PostProblem currently is the only place
    // where it is used.
    // Naming "SetupGhosting" is misleading, implementation duplicates functionality that
    // can be build using some of the other methods
    // mwgee 4/12
    friend class PostProblem;

   public:
    //-----------------------------------------------------------------------------------
    /*!
    \brief Build element row and column map from nodal row and column maps
           (Filled()==true NOT prerequisite)

    Create a unique element map elerowmap assigning each element an owner.
    Create an overlapping element map elecolmap representing a one layer
    overlap of element.
    Maps are created such that they match the nodal row and column maps provided.

    \param noderowmap (in): unique nodal row map of some distribution
    \param nodecolmap (in): overlapping nodal column map
    \param ExtendedGhosting (in): allow for ghosting of elements on processors, which don't own any
    of the elements' nodes (default: false) // rauch 09/13

    @return unique element row map and overlapping element column map

    \note The provided noderowmap and nodecolmap do not need to match the
    distribution of nodes in this discretization class. Also, the output
    element maps do not match the distribution of elements in this class.
    Total numbers of nodes and elements have to match nodes and elements in
    this class. The status of this->Filled() is not changed by this method.
    Neither nodes nor elements are actually redistributed here, this method
    only builds maps!

    */
    virtual std::pair<Teuchos::RCP<Epetra_Map>, Teuchos::RCP<Epetra_Map>> build_element_row_column(
        const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap) const;

    /*!
    \brief Export the nodes to a different parallel layout
           (Filled()==true NOT prerequisite)

    The discretization has a parallel layout of nodes reflected in NodeRowMap().
    This method communicates the nodes such that after
    the export the nodes are stored as provided in newmap.<br>
    There are some important aspects to this method:<br>
    - Filled()=false on exit. This means fill_complete() needs to be called again.<br>
    - This is a dull export meaning that there is no notion of whether the
      exported node distribution still matches the element distribution. A call
      to fill_complete() might therefore not be possible until elements are
      distributed accordingly.<br>
    - newmap has to be a non-overlapping map (will be tested) as
      the ownership of an exported node changes to the receiving proc<br>
    - All ghosted nodes on all processors will be destroyed.
      Ghosting has to be recreated afterwards<br>

    \param newmap (in): new nodal row map the discretization should use

    \note Sets Filled()=false and deletes noderowmap_ and nodecolmap_
    */
    virtual void ExportRowNodes(
        const Epetra_Map& newmap, bool killdofs = true, bool killcond = true);


    /*!
    \brief Export overlap of nodes
          (Filled()==true NOT prerequisite)

    The discretization has a parallel layout of nodes reflected in NodeRowMap().
    This method communicates the nodes such that after
    the export the nodes are stored as provided in newmap.<br>
    There are some important aspects to this method:<br>
    - All existing ghosted nodes on all processors will be destroyed
      before the communication.<br>
    - newmap must contain all nodes of noderowmap_ (will be tested) because
      otherwise a node is shipped to a different proc and deleted from the
      originating proc. It then merely exists as a ghost node on the receiving proc.
      (which is a state not tolerated)<br>
    - Filled()=false on exit. This means fill_complete() needs to be called again.<br>
    - This is a dull export meaning that there is no notion of whether the
      exported node distribution still matches the element distribution. A call
      to fill_complete() might therefore not be possible until elements are
      distributed accordingly.<br>
    - newmap should be an overlapping map.<br>
    - The ownership of an exported node does not change on the receiving proc.<br>
      The received node becomes a ghost node on the receiving proc.

    \note Sets Filled()=false and deletes noderowmap_ and nodecolmap_
    */
    virtual void ExportColumnNodes(
        const Epetra_Map& newmap, bool killdofs = true, bool killcond = true);

    /*!
    \brief Export the elements from proc 0 to a different parallel row layout
          (Filled()==true NOT prerequisite)

    This routine is mainly useful for input read on proc 0

    \param target (in): desired distirbution of elements
    \param gidlist (in): list of element gids to be distributed

    */
    void proc_zero_distribute_elements_to_all(Epetra_Map& target, std::vector<int>& gidlist);

    /*!
    \brief Export the nodes from proc 0 to a different parallel row layout
          (Filled()==true NOT prerequisite)

    This routine is mainly useful for input read on proc 0

    \param target (in): desired distirbution of elements

    */
    void proc_zero_distribute_nodes_to_all(Epetra_Map& target);

    /*!
    \brief Export the elements to a different parallel row layout
          (Filled()==true NOT prerequisite)

    The discretization has a parallel layout of elements reflected in ElementRowMap().
    This method communicates the elements in this discretization such that after
    the export the elements are stored as provided in newmap.<br>
    There are some important aspects to this method:<br>
    - Filled()=false on exit. This means fill_complete() needs to be called again.<br>
    - This is a dull export meaning that there is no notion of whether the
      exported element distribution still matches the node distribution. A call
      to fill_complete() might therefore not be possible until nodes are
      distributed accordingly.<br>
    - newmap has to be a non-overlapping map (will be tested) as
      the ownership of an exported element changes to the receiving proc<br>
    - All ghosted elements on all processors will be destroyed.
      Ghosting has to be recreated afterwards<br>

    \note Sets Filled()=false and deletes elerowmap_ and elecolmap_
    */
    virtual void ExportRowElements(
        const Epetra_Map& newmap, bool killdofs = true, bool killcond = true);


    /*!
    \brief Export overlap of elements
          (Filled()==true NOT prerequisite)

    The discretization has a parallel layout of elements reflected in ElementRowMap().
    This method communicates the elements such that after
    the export the elements are stored as provided in newmap.<br>
    There are some important aspects to this method:<br>
    - All existing ghosted elements on all processors will be destroyed
      before the communication.<br>
    - newmap must contain all elements of elerowmap_ (will be tested) because
      otherwise an element is shipped to a different proc and deleted from the
      originating proc. It then merely exists as a ghost element on the receiving proc.<br>
    - Filled()=false on exit. This means fill_complete() needs to be called again.<br>
    - This is a dull export meaning that there is no notion of whether the
      exported element distribution still matches the node distribution. A call
      to fill_complete() might therefore not be possible until elements are
      distributed accordingly.<br>
    - newmap should be an overlapping map<br>
    - The ownership of an exported element does not change on the receiving proc.<br>
      The received element becomes a ghost element on the receiving processor.

    \note Sets Filled()=false and deletes elerowmap_ and elecolmap_
    */
    virtual void export_column_elements(
        const Epetra_Map& newmap, bool killdofs = true, bool killcond = true);

    /*!
    \brief Build nodal graph of discretization (Filled()==true prerequisite)

    Build a nodal graph of the discretization in parallel.<br>
    The graph has a row map of NodeRowMap().<br>
    The graph is build from elements stored on each proc, nodes are
    not referenced.<br>
    If a proc stores the appropiate ghosted elements the resulting graph
    will be the complete graph of the distributed discretization.<br>
    If procs do not store appropiate ghosted elements, the resulting
    graph is decoupled or partially decoupled among procs.<br>
    This might also lead to an unsymmetric graph.

    \note Filled()=true is a prerequisite

    \return Graph of discretization distributed across processors according to
            the discretization distribution
    */
    virtual Teuchos::RCP<Epetra_CrsGraph> BuildNodeGraph() const;

    /*!
    \brief Build nodal coordinate vector of this discretization based on a nodal rowmap

    Per default the node map of the underlying discretization is used to build the nodal
    coordinate vector. The user can specify a node rowmap to calculate the coordinate vector on.
    The map has to be a submap of the overall full node rowmap of this discretization.

    \pre discretization does have to be fill_complete().

    @param[in] noderowmap Map representing the nodal distribution of the discretization (can be a
               submap of the full node rowmap returned from NodeRowMap())
    @return Vector containing the coordinates of all nodes which are present in the given noderowmap
     */
    Teuchos::RCP<Epetra_MultiVector> build_node_coordinates(
        Teuchos::RCP<const Epetra_Map> noderowmap = Teuchos::null) const;

    //@}

    //! @name Evaluate methods

    /*!
    \brief Set a reference to a data vector at the default dofset (0)

    Using this method, a reference to a vector can
    be supplied to the discretization. The elements can access
    this vector by using the name of that vector.
    The method expects state to be either of dof row map or of
    dof column map.
    If the vector is supplied in DofColMap() a reference to it will be stored.
    If the vector is NOT supplied in DofColMap(), but in dof_row_map(),
     a vector with column map is allocated and the supplied vector is exported to it.
    Everything is stored/referenced using Teuchos::RCP.

    \param name (in): Name of data
    \param state (in): vector of some data

    \note This class will not take ownership or in any way modify the solution vector.
    */
    void set_state(const std::string& name, Teuchos::RCP<const Epetra_Vector> state)
    {
      set_state(0, name, state);
    }

    /*!
    \brief Set a reference to a data vector

    Using this method, a reference to a vector can
    be supplied to the discretization. The elements can access
    this vector by using the name of that vector.
    The method expects state to be either of dof row map or of
    dof column map.
    If the vector is supplied in DofColMap() a reference to it will be stored.
    If the vector is NOT supplied in DofColMap(), but in dof_row_map(),
     a vector with column map is allocated and the supplied vector is exported to it.
    Everything is stored/referenced using Teuchos::RCP.

    \param nds (in): number of dofset
    \param name (in): Name of data
    \param state (in): vector of some data

    \note This class will not take ownership or in any way modify the solution vector.
    */
    virtual void set_state(
        unsigned nds, const std::string& name, Teuchos::RCP<const Epetra_Vector> state);

    /*!
    \brief Get a reference to a data vector at the default dofset (0)

    Providing a name of a solution state, get a reference to the solution vector.
    If a vector under the provided name does not exist, the method will throw
    a lethal error message.

    \param name (in): Name of solution state

    \return Reference to solution state
    */
    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> GetState(const std::string& name) const
    {
      return GetState(0, name);
    }

    /*!
    \brief Get a reference to a data vector

    Providing a name of a solution state, get a reference to the solution vector.
    If a vector under the provided name does not exist, the method will throw
    a lethal error message.

    \param nds (in): number of dofset
    \param name (in): Name of solution state

    \return Reference to solution state
    */
    [[nodiscard]] virtual Teuchos::RCP<const Epetra_Vector> GetState(
        unsigned nds, const std::string& name) const
    {
      FOUR_C_ASSERT(
          nds < dofsets_.size(), "undefined dof set found in discretization %s!", name_.c_str());
      if (state_.size() <= nds) return Teuchos::null;

      auto state_iterator = state_[nds].find(name);
      FOUR_C_THROW_UNLESS(state_iterator != state_[nds].end(),
          "Cannot find state %s in discretization %s", name.data(), name_.c_str());
      return state_iterator->second;
    }

    /*!
      \brief Tell whether a state vector has been set
      \param nds (in): number of dofset
      \param name (in): Name of solution states
     */
    virtual bool HasState(unsigned nds, const std::string& name) const
    {
      FOUR_C_ASSERT(
          nds < dofsets_.size(), "undefined dof set found in discretization %s!", name_.c_str());
      if (state_.size() <= nds) return false;

      return state_[nds].count(name) == 1;
    }

    /*!
    \brief Clear solution state references

    The method deletes references to any solution data
    */
    virtual void ClearState(bool clearalldofsets = false)
    {
      // clear all states
      if (clearalldofsets) state_.clear();
      // clear states that belong to own dofset only
      else if (!state_.empty())
        state_[0].clear();
    }

    virtual void RedistributeState(unsigned nds, const std::string& name);


    /*!
      \brief Tell whether a state vector has been set
      \param name (in): Name of solution state
     */
    virtual bool HasState(const std::string& name) const { return HasState(0, name); }


    /*!
    \brief Export a node based Epetra_MultiVector to column layout
           and store it in the given ParameterList

    Using this method, the node based vector is exported to column
    layout and stored in the given ParameterList. The elements can access
    this vector by using the name of that vector in the corresponding
    ParameterList.

    \param nds (in): ParameterList in which vector is stored
    \param name (in): Name of data
    \param state (in): vector of some data

    */
    virtual void add_multi_vector_to_parameter_list(Teuchos::ParameterList& p,
        const std::string name, Teuchos::RCP<const Epetra_MultiVector> vec);

    /*!
    \brief Call elements to evaluate

    Call element routines to perform integration and return element contributions to
    system vectors and matrices. Type of action taken by the elements is
    controlled by the params parameter.<br>
    Parameters that control element behavior are:<br>
    \code
    params.set("action","<element_action>"); // <element_action> somthing that elements understand
    \endcode
    Other parameters eventually recognized by the elements:<br>
    \code
    params.set("total time",1.23);     // current total time of simulation
    params.set("delta time",0.01);     // time increment
    \endcode


    \param params (in): Parameter list past to the elements containing
                        commands and parameters for the elements and
                        containing assembly intructions
    \param systemmatrix1 (out)   : Sparse matrix that may be filled by
                                   assembly of element contributions.
                                   May be Teuchos::null on entry.
                                   Matrix must be systemmatrix1->Filled()==false on input.
    \param systemmatrix2 (out):    Sparse matrix that may be filled by
                                   assembly of element contributions.
                                   May be Teuchos::null on entry.
                                   Matrix must be systemmatrix2->Filled()==false on input.
    \param systemvector1 (out):    Distributed vector that may be filled by
                                   aasembly of element contributions.
                                   May be Teuchos::null on entry.
                                   Vector will NOT be initialized to zero by
                                   the underlying assembly methods that add element
                                   contributions.
    \param systemvector2 (out):    Distributed vector that may be filled by
                                   aasembly of element contributions.
                                   May be Teuchos::null on entry.
                                   Vector will NOT be initialized to zero by
                                   the underlying assembly methods that add element
                                   contributions.
    \param systemvector3 (out):    Distributed vector that may be filled by
                                   aasembly of element contributions.
                                   May be Teuchos::null on entry.
                                   Vector will NOT be initialized to zero by
                                   the underlying assembly methods that add element
                                   contributions.
    */
    virtual void Evaluate(Teuchos::ParameterList& params,
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1,
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2,
        Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
        Teuchos::RCP<Epetra_Vector> systemvector3);

    /// Call elements to evaluate
    virtual void Evaluate(Teuchos::ParameterList& params, CORE::FE::AssembleStrategy& strategy);

    /**
     * Loop over all elements of the discretization and perform the given @p element_action. In
     * contrast to the other overloads of Evaluate(), this function allows to perform any local
     * action on an Element that can be encoded within the passed function object @p element_action.
     * This is very useful for one-off actions, that one does not want to implement inside the
     * actual Element's Evaluate call.
     */
    virtual void Evaluate(Teuchos::ParameterList& params, CORE::FE::AssembleStrategy& strategy,
        const std::function<void(CORE::Elements::Element&, CORE::Elements::Element::LocationArray&,
            CORE::LINALG::SerialDenseMatrix&, CORE::LINALG::SerialDenseMatrix&,
            CORE::LINALG::SerialDenseVector&, CORE::LINALG::SerialDenseVector&,
            CORE::LINALG::SerialDenseVector&)>& element_action);

    /// Call elements to evaluate
    /*!
      Abbreviated Evaluate() call that always assembles one matrix and
      one vector. No need to set assemble instructions to the
      Teuchos::ParameterList.

      \param params (in): Parameter list past to the elements containing
                          commands and parameters for the elements and
                          containing assembly intructions
      \param systemmatrix (out) : Sparse matrix that may be filled by
                                  assembly of element contributions.
                                  May not be Teuchos::null.
                                  Matrix must be systemmatrix->Filled()==false on input.
      \param systemvector (out) : Distributed vector that may be filled by
                                  assembly of element contributions.
                                  May not be Teuchos::null.
                                  Vector will NOT be initialized to zero by
                                  the underlying assembly methods that add element
                                  contributions.
     */
    void Evaluate(Teuchos::ParameterList& params,
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<Epetra_Vector> systemvector);


    /// Call elements to evaluate
    /*!
    Call element routines to perform tasks defined in the parameter list and
    return element information via the parameter list. The element loop
    is over all col elements, hence, Filled() status is tested.
    HaveDofs() is NOT tested and the lm vector is empty. No assembly is performed.
    All states are cleared to avoid errors.

    \param params (in/out): Parameter list past to the elements containing
                            commands and parameters for the elements and
                            containing assembly instructions
   */
    void Evaluate(Teuchos::ParameterList& params);

    virtual void Evaluate(const std::function<void(CORE::Elements::Element&)>& element_action);

    /** \brief Evaluate Neumann boundary conditions
     *
     *  Loop all Neumann conditions attached to the discretization and evaluate them.
     *  This method considers all conditions in condition_ with the names
     *  "PointNeumann", LineNeumann", "SurfaceNeumann" and "VolumeNeumann".
     *  It takes a current time from the parameter list params named "total time"
     *  and evaluates the appropiate time curves at that time for each
     *  Neumann condition separately. If "total time" is not included
     *  in the parameters, no time curves are used.
     *  Parameters recognized by this method:
     *  \code
     *  params.set("total time",acttime); // current total time
     *  \endcode
     *
     *  \param params        (in): List of parameters
     *  \param systemvector (out): Vector to assemble Neumann BCs to.
     *                             The vector is NOT initialized to zero by
     *                             this method.
     *  \param systemvector (out): Optional matrix to assemble linearization
     *                             of Neumann BCs to. */
    void evaluate_neumann(Teuchos::ParameterList& params, Teuchos::RCP<Epetra_Vector> systemvector,
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix = Teuchos::null)
    {
      if (systemmatrix.is_null())
        evaluate_neumann(params, *systemvector);
      else
        evaluate_neumann(params, *systemvector, systemmatrix.get());
    }

    /*!
    \brief Evaluate Neumann boundary conditions

    Loop all Neumann conditions attached to the discretization and evaluate them.
    This method considers all conditions in condition_ with the names
    "PointNeumann", LineNeumann", "SurfaceNeumann" and "VolumeNeumann".
    It takes a current time from the parameter list params named "total time"
    and evaluates the appropiate time curves at that time for each
    Neumann condition separately. If "total time" is not included
    in the parameters, no time curves are used.
    Parameters recognized by this method:
    \code
      params.set("total time",acttime); // current total time
    \endcode

    \param params (in): List of parameters
    \param systemvector (out): Vector to assemble Neumann BCs to.
                               The vector is NOT initialized to zero by this method.
    */
    virtual void evaluate_neumann(Teuchos::ParameterList& params, Epetra_Vector& systemvector,
        CORE::LINALG::SparseOperator* systemmatrix = nullptr);

    /*!
    \brief Evaluate Dirichlet boundary conditions

    Loop all Dirichlet conditions attached to the discretization and evaluate them.
    This method considers all conditions in condition_ with the names
    "PointDirichlet", "LineDirichlet", "SurfaceDirichlet" and "VolumeDirichlet".
    It takes a current time from the parameter list params named "total time"
    and evaluates the appropiate time curves at that time for each
    Dirichlet condition separately. If "total time" is not included
    in the parameters, no time curves are used.

    \note Opposed to the other 'Evaluate' method does this one NOT assembly but
          OVERWRITE values in the output vector systemvector. For this reason,
          dirichlet boundary conditions are evaluated in the following order:
          First "VolumeDirichlet", then "SurfaceDirichlet", then "LineDirichlet"
          and finally "PointDirichlet". This way, the lower entity dirichlet BCs override
          the higher ones and a point Dirichlet BCs has priority over other dirichlet
          BCs in the input file.

    Parameters recognized by this method:
    \code
    params.set("total time",acttime); // current total time
    \endcode

    \param params (in): List of parameters
    \param systemvector (out): Vector holding prescribed Dirichlet values
    \param systemvectord (out): Vector holding 1st time derivative of prescribed Dirichlet values
    \param systemvectordd (out): Vector holding 2nd time derivative prescribed Dirichlet values
    \param toggle (out): Vector containing 1.0 for each Dirichlet dof and 0 for everything else
    \param dbcmapextractor (out): Map extractor containing maps for the DOFs subjected to
                                  Dirichlet boundary conditions and the remaining/free DOFs
    */
    virtual void evaluate_dirichlet(Teuchos::ParameterList& params,
        Teuchos::RCP<Epetra_Vector> systemvector, Teuchos::RCP<Epetra_Vector> systemvectord,
        Teuchos::RCP<Epetra_Vector> systemvectordd, Teuchos::RCP<Epetra_IntVector> toggle,
        Teuchos::RCP<CORE::LINALG::MapExtractor> dbcmapextractor = Teuchos::null) const;

    /// Evaluate a specific condition using assemble strategy
    virtual void evaluate_condition(Teuchos::ParameterList& params,
        CORE::FE::AssembleStrategy& strategy, const std::string& condstring, const int condid = -1);

    /** \brief Evaluate a specified condition
     *
     *  Loop all conditions attached to the discretization and evaluate them.
     *  This method considers all conditions in condition_ with the names
     *  matching the user-provided string condstring.
     *  Calls more general evaluate_condition method, see below.
     *
     *  \param params        (in): List of parameters for use at element level
     *  \param systemvector (out): Vector to assemble BCs to.(NOT initialized to zero
     *                             by this method)
     *  \param condstring    (in): Name of condition to be evaluated
     *  \param condid        (in): condition ID */
    void evaluate_condition(Teuchos::ParameterList& params,
        Teuchos::RCP<Epetra_Vector> systemvector, const std::string& condstring,
        const int condid = -1)
    {
      evaluate_condition(params, Teuchos::null, Teuchos::null, systemvector, Teuchos::null,
          Teuchos::null, condstring, condid);
    }

    /** \brief Evaluate a specified condition
     *
     *  Loop all conditions attached to the discretization and evaluate them.
     *  This method considers all conditions in condition_ with the names
     *  matching the user-provided string condstring.
     *  Calls more general evaluate_condition method, see below.
     *
     *  \param params     (in): List of parameters for use at element level
     *  \param condstring (in): Name of condition to be evaluated
     *  \param condid     (in): condition ID */
    void evaluate_condition(
        Teuchos::ParameterList& params, const std::string& condstring, const int condid = -1)
    {
      evaluate_condition(params, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null,
          Teuchos::null, condstring, condid);
    }

    /*!
    \brief Evaluate a specific condition

    Loop all conditions attached to the discretization and evaluate them.
    This method considers all conditions in condition_ with the names
    matching the user-provided string condstring.
    It takes a current time from the parameter list params named "total time"
    and evaluates the appropriate time curves at that time for each
    condition separately. If "total time" is not included
    in the parameters, no time curves are used.

      \param params (in):        List of parameters for use at element level
      \param systemmatrix1 (out): Sparse matrix that may be changed by
                                 assembly of boundary element contributions.
                                 May not be Teuchos::null.
                                 Matrix must be systemmatrix->Filled()==false on input.
      \param systemmatrix2 (out): Sparse matrix that may be changed by
                                 assembly of boundary element contributions.
                                 May not be Teuchos::null.
                                 Matrix must be systemmatrix->Filled()==false on input.
      \param systemvector1 (out):Vector to assemble BCs to.
                                 The vector is NOT initialized to zero by this method.
      \param systemvector2 (out):Vector to assemble BCs to.
                                 The vector is NOT initialized to zero by this method.
      \param systemvector3 (out):Vector to assemble BCs to.
                                 The vector is NOT initialized to zero by this method.
      \param condstring (in):    Name of condition to be evaluated
      \param condid (in):        Condition ID
      */
    virtual void evaluate_condition(Teuchos::ParameterList& params,
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1,
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2,
        Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
        Teuchos::RCP<Epetra_Vector> systemvector3, const std::string& condstring,
        const int condid = -1);

    /*!
    \brief Compute the nullspace of the discretization

    This method looks in the solver parameters whether algebraic multigrid (AMG)
    is used as preconditioner. AMG desires the nullspace of the
    system of equations which is then computed here if it does not already exist
    in the parameter list.

    \param solveparams (in): List of parameters
    \param recompute (in)  : force method to recompute the nullspace
    */
    virtual void compute_null_space_if_necessary(
        Teuchos::ParameterList& solveparams, bool recompute = false);

    /*!
     * \brief Assemble scalar quantities across elements
     *
     * Every element is only called \b once by its owning processor. (We are
     * parsing the element row map.) At this call the element return its value(s),
     * i.e. its scalar(s), and contributes to the global
     * value of the respective scalar quantity(ies).
     *
     * Example: strain energy in structures
     *
     * \author bborn \date 08/08
     */
    virtual void EvaluateScalars(
        Teuchos::ParameterList& params, Teuchos::RCP<CORE::LINALG::SerialDenseVector> scalars);

    /*!
     * \brief Assemble scalar quantities across conditioned elements
     *
     * Every element is only called \b once by its owning processor. (We are
     * checking this inside the function.) At this call the element return its value(s),
     * i.e. its scalar(s), and contributes to the global
     * value of the respective scalar quantity(ies).
     *
     * Example: strain energy in structures
     *
     * \author fang \date 02/15
     */
    void EvaluateScalars(Teuchos::ParameterList& params,  //! (in) parameter list
        Teuchos::RCP<CORE::LINALG::SerialDenseVector>
            scalars,                    //! (out) result vector for scalar quantities to be computed
        const std::string& condstring,  //! (in) name of condition to be evaluated
        const int condid = -1           //! (in) condition ID (optional)
    );

    /*!
     * \brief Assemble scalar quantities for each element separately
     *
     * Every element is only called \b once by its owning processor. (We are
     * parsing the element row map.) At this call the element return its value,
     * i.e. its scalar. The scalar is put in a Epetra_Vector called "element scalar"
     * iun the parameter list
     *
     * Example: strain energy in structures
     *
     * \author gee \date 05/11
     */
    void EvaluateScalars(Teuchos::ParameterList& params, /*!< parameters */
        Teuchos::RCP<Epetra_MultiVector> scalars /*!< output element-wise scalar quantities */
    );

    /*!
    \brief Evaluate a specified initial field (scalar or vector field)

    Loop all intial field conditions attached to the discretization and
    evaluate them if their names match the user-provided string fieldstring.
    Information on which local DOFs ids are addressed by the condition
    MUST be pre-defined and is represented by the locids vector.
    As an example, if we provide an initial velocity for a 3D structural
    dynamics simulation, locids must contain the local DOF ids {0,1,2}.
    Another example would be prescribing an initial pressure in a 3D
    fluid dynamics simulation, where locids would have to contain only
    the local pressure DOF id, namely {3}.

    */
    void evaluate_initial_field(const CORE::UTILS::FunctionManager& function_manager,
        const std::string& fieldstring, Teuchos::RCP<Epetra_Vector> fieldvector,
        const std::vector<int>& locids) const;

    //@}

    //! @name IO methods

    /*!
      \brief Pack local elements (row map) into buffer.

      Call Pack on all local (row map) elements and put the results into
      a common vector. This is used to output the discretization.

      \note Filled()=true is a prerequisite
     */
    Teuchos::RCP<std::vector<char>> PackMyElements() const;

    /*!
      \brief Pack local nodes (row map) into buffer.

      Call Pack on all local (row map) nodes and put the results into
      a common vector. This is used to output the discretization.

      \note Filled()=true is a prerequisite
     */
    Teuchos::RCP<std::vector<char>> PackMyNodes() const;

    /*!
      \brief Unpack element buffer and create local elements.

      Interprets the argument as packed elements and unpacks them on the
      local processor. Takes the ownership of the unpacked elements.

      \param e (in): buffer of packed elements

      \note Sets Filled()=false
     */
    void UnPackMyElements(Teuchos::RCP<std::vector<char>> e);

    /*!
      \brief Unpack nodal buffer and create local nodes.

      Interprets the argument as packed nodes and unpacks them on the
      local processor. Takes the ownership of the unpacked nodes.

      \param e (in): buffer of packed nodes

      \note Sets Filled()=false
     */
    void UnPackMyNodes(Teuchos::RCP<std::vector<char>> e);

    //@}

    /*!
    \brief Assign degrees of freedom to discretization (Filled()==true prerequisite)

    Assign nodes and elements their no. of degrees of freedom as aquired
    by Element::NumDofPerNode and Element::num_dof_per_element.
    Number degrees of freedom (dofs) ascending according to global node numbers
    followed by global element numbers

    \param start (int): first dof number to assign

    \return last dof assigned + 1
    */
    virtual int assign_degrees_of_freedom(int start);

   private:
    /*!
    \brief Build noderowmap_ (Filled()==true NOT prerequisite)

    Build the parallel layout of nodes in this
    discretization and store it as an Epetra_Map in noderowmap_
    noderowmap_ is unique.
    It considers nodes owned by a proc only.

    \note This is a collective call
    */
    virtual void build_node_row_map();

    /*!
    \brief Build nodecolmap_ (Filled()==true NOT prerequisite)

    Build the parallel layout of nodes in this
    discretization and store it as an Epetra_Map in nodecolmap_
    nodecolmap_ is potentially but not necessarily overlapping.
    It considers nodes owned by a proc and its ghosted nodes

    \note This is a collective call
    */
    virtual void build_node_col_map();

    /*!
    \brief Build elerowmap_ (Filled()==true NOT prerequisite)

    Build the parallel layout of elements in this
    discretization and store it as an Epetra_Map in elerowmap_
    elerowmap_ is unique.
    It considers elements owned by a proc only

    \note This is a collective call

    */
    virtual void build_element_row_map();

    /*!
    \brief Build elecolmap_ (Filled()==true NOT prerequisite)

    Build the potentially overlapping parallel layout of elements in this
    discretization and store it as an Epetra_Map in elecolmap_
    elecolmap_ includes ghosted elements and is potentially overlapping.

    \note This is a collective call

    */
    virtual void build_element_col_map();

    /*!
    \brief Build pointers elements -> Nodes (Filled()==true NOT prerequisite)
    */
    virtual void build_element_to_node_pointers();

    /*!
    \brief Build pointers Node -> Element (Filled()==true NOT prerequisite)

    \note This is a collective call
    */
    virtual void build_node_to_element_pointers();

    /*!
    \brief Build pointers Element -> Element (Filled()==true NOT prerequisite)

    \note This is a collective call
    */
    virtual void build_element_to_element_pointers();

    /*!
    \brief Build the geometry of surfaces belonging to the structure-fluid
    volume coupling condition -> this is special since an associated volume
    conditions also needs to be considered
    */
    void find_associated_ele_i_ds(Teuchos::RCP<CORE::Conditions::Condition> cond,
        std::set<int>& VolEleIDs, const std::string& name);

   protected:
    /*!
    \brief Build the geometry of lines for a certain line condition

    */
    virtual bool build_linesin_condition(
        const std::string& name, Teuchos::RCP<CORE::Conditions::Condition> cond);

    /*!
    \brief Build the geometry of surfaces for a certain surface condition.

    \note On exit, the parallel distribution of the newly created surface geometry
    matches the distribution of the originating volume. Each surface is owned by the
    same proc as the underlying volume element (rauch 10/16).

    \author Michael Gee
    \date 01/07

    \version rework by Andreas Rauch ( rauch 10/16 )       */
    virtual bool build_surfacesin_condition(
        const std::string& name, Teuchos::RCP<CORE::Conditions::Condition> cond);

    /*!
    \brief Build the geometry of volumes for a certain volume condition

    */
    virtual bool build_volumesin_condition(
        const std::string& name, Teuchos::RCP<CORE::Conditions::Condition> cond);

    /*!
    \brief Reset all maps and set Filled()=false (Filled()==true NOT prerequisite)

    Resets all maps and sets flags filled_ and havedof_ to false.

    \param killdofs (in): if true reset existing dofsets in discretization

    \note This is a collective call
    */
    // set protected to be accessible from derived class DRT::MESHFREE::MeshfreeDiscretization
    // (nis) Jan12
    virtual void reset(bool killdofs, bool killcond);
    virtual void reset() { this->reset(true, true); }

    /*!
    \brief Initialize element routines

    \note initialize_elements might be called more then once!

    */
    virtual void initialize_elements();

    /*!
    \brief Build the geometry for boundary conditions

    */
    virtual void boundary_conditions_geometry();

    /*!
     *  A helper function for build_surfacesin_condition,
     *  build_linesin_condition, BuildInternalFaces, etc.
     *
     *  A helper method for build_linesin_condition and
     *  build_surfacesin_condition, below.
     *  Gets a map (vector_of_nodes)->Element that maps
     *
     *  (A map with globally unique ids.)
     *
     * \note The point here is to make sure the element gid are the same on any
     *  parallel distribution of the elements. Thus we allreduce thing to
     *  processor 0 and sort the element descriptions (vectors of nodal ids)
     * there.
     *
     *  \warning This routine has not been optimized for efficiency. I don't think that is
     *  needed.
     *
     *  \param comm (i) communicator
     *  \param elementmap (i) map (vector_of_nodes_ids)->(element) that maps
     *  the nodes of an element to the element itself.
     *  \param finalelements (o) map (global_id)->(element) that can be
     *  added to a condition.
     *
     *  \author h.kue
     *  \date 09/07    */
    virtual void assign_global_i_ds(const Epetra_Comm& comm,
        const std::map<std::vector<int>, Teuchos::RCP<CORE::Elements::Element>>& elementmap,
        std::map<int, Teuchos::RCP<CORE::Elements::Element>>& finalgeometry);


    //! Name of this discretization
    std::string name_;

    //! Epetra_comm
    Teuchos::RCP<Epetra_Comm> comm_;

    //! DiscretizationWriter
    Teuchos::RCP<CORE::IO::DiscretizationWriter> writer_;

    //! Flag indicating whether fill_complete() has been called
    bool filled_;

    //! Flag indicating whether degrees of freedom where assigned
    bool havedof_;


    //! @name Elements
    //! @{

    //! Unique distribution of element ownerships
    Teuchos::RCP<Epetra_Map> elerowmap_;

    //! Distribution of elements including ghost elements
    Teuchos::RCP<Epetra_Map> elecolmap_;

    //! Vector of pointers to row elements for faster access
    std::vector<CORE::Elements::Element*> elerowptr_;

    //! Vector of pointers to column elements for faster access
    std::vector<CORE::Elements::Element*> elecolptr_;

    //! Map of elements
    std::map<int, Teuchos::RCP<CORE::Elements::Element>> element_;

    //! @}

    //! @name Nodes
    //! @{

    //! Unique distribution of nodal ownerships
    Teuchos::RCP<Epetra_Map> noderowmap_;

    //! Distribution of nodes including ghost nodes
    Teuchos::RCP<Epetra_Map> nodecolmap_;

    //! Vector of pointers to row nodes for faster access
    std::vector<CORE::Nodes::Node*> noderowptr_;

    //! Vector of pointers to column nodes for faster access
    std::vector<CORE::Nodes::Node*> nodecolptr_;

    //! Map from nodal Gid to node pointers
    std::map<int, Teuchos::RCP<CORE::Nodes::Node>> node_;

    //! @}

    //! Map of references to solution states
    std::vector<std::map<std::string, Teuchos::RCP<const Epetra_Vector>>> state_;

    ///< Map of import objects for states
    std::vector<Teuchos::RCP<Epetra_Import>> stateimporter_;

    ///< Some conditions e.g. boundary conditions
    std::multimap<std::string, Teuchos::RCP<CORE::Conditions::Condition>> condition_;

    //! Vector of DofSets
    std::vector<Teuchos::RCP<CORE::Dofsets::DofSetInterface>> dofsets_;

    //! number of space dimension
    const unsigned int n_dim_;
  };  // class Discretization
}  // namespace DRT

/// << operator
std::ostream& operator<<(std::ostream& os, const DRT::Discretization& dis);

FOUR_C_NAMESPACE_CLOSE

#endif
