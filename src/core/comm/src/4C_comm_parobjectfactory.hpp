/*----------------------------------------------------------------------*/
/*! \file

\brief Central type object management.

\level 0

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_COMM_PAROBJECTFACTORY_HPP
#define FOUR_C_COMM_PAROBJECTFACTORY_HPP

#include "4C_config.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <map>
#include <set>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class SparseOperator;
}

namespace Input
{
  class LineDefinition;
}

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Elements
{
  class Element;
  class ElementType;
}  // namespace Core::Elements

namespace Core::Communication
{
  class ParObject;

  /// abstract ParObject type information
  /*!
    There needs to be one ParObjectType subclass for every ParObject. The
    ParObjectType subclass must be a singleton, that is there must be just one
    instance and it must be created at startup.

    ParObjectType is responsible to register with ParObjectFactory and provide
    means to create its ParObject.

    \author u.kue
    \date 06/10
   */
  class ParObjectType
  {
    friend class ParObjectFactory;

   public:
    ParObjectType();

    virtual ~ParObjectType() = default;
    /// Unique ParObject id. Automatically determined.
    int unique_par_object_id();

    /// Create ParObject from packed data
    virtual ParObject* create(const std::vector<char>& data) { return nullptr; }

    /// internal name of this ParObjectType.
    virtual std::string name() const = 0;

    /// test equality by comparing pointers as there are just singletons allowed
    bool operator==(const ParObjectType& other) const { return this == &other; }

    /// test inequality by comparing pointers as there are just singletons allowed
    bool operator!=(const ParObjectType& other) const { return this != &other; }

   private:
    int objectid_;
  };


  /// Singleton ParObject factory class
  /*!
    Central ParObject factory. A singleton class. Each ParObjectType registers
    itself here. The factory is able to loop all ParObjectType objects and call
    the appropriate creation methods.

    \author u.kue
    \date 06/10
   */
  class ParObjectFactory
  {
    friend class ParObjectType;

   public:
    static ParObjectFactory& instance();

    /// Virtual destructor.
    virtual ~ParObjectFactory() = default;

    /// create a parobject from its data stream
    ParObject* create(const std::vector<char>& data);

    /// create an element from its name (and dis type if needed)
    Teuchos::RCP<Core::Elements::Element> create(
        const std::string eletype, const std::string eledistype, const int id, const int owner);

    /// initialize all element types
    void initialize_elements(Core::FE::Discretization& dis);

    /// preevaluate elements (via element types)
    void pre_evaluate(Core::FE::Discretization& dis, Teuchos::ParameterList& p,
        Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix1,
        Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix2,
        Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
        Teuchos::RCP<Epetra_Vector> systemvector3);

    /// setup definition of element input file lines
    void setup_element_definition(
        std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions);

   private:
    void do_register(ParObjectType* object_type);

    void finalize_registration();

    ParObjectFactory() = default;

    /// list of parobject types

    /// Id to type object map. The central beast.
    std::map<int, ParObjectType*> type_map_;

    /// element name cache
    std::map<std::string, Core::Elements::ElementType*> element_cache_;

    /// preregistered types
    std::vector<ParObjectType*> types_;

    /// element types that are actually used
    std::map<Core::FE::Discretization*, std::set<Core::Elements::ElementType*>> active_elements_;

    // no copying

    ParObjectFactory(const ParObjectFactory&);
    ParObjectFactory& operator=(const ParObjectFactory&);
  };

}  // namespace Core::Communication

FOUR_C_NAMESPACE_CLOSE

#endif
