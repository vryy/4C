/*----------------------------------------------------------------------*/
/*!
\file drt_parobjectfactory.H

\brief Central type object management.

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef DRT_PAROBJECTFACTORY_H
#define DRT_PAROBJECTFACTORY_H

#include <vector>
#include <map>
#include <set>
#include <string>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

#include "drt_singletondestruction.H"

namespace LINALG
{
  class SparseOperator;
}

namespace DRT
{
  namespace INPUT
  {
    class LineDefinition;
  }

  class ParObject;
  class Element;
  class ElementType;
  class Discretization;


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

    virtual ~ParObjectType() {}

    /// Unique ParObject id. Automatically determined.
    int UniqueParObjectId();

    /// Create ParObject from packed data
    virtual DRT::ParObject* Create(const std::vector<char>& data) { return NULL; }

    /// internal name of this ParObjectType.
    virtual std::string Name() const = 0;

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
  class ParObjectFactory : public SingletonDestruction
  {
    friend class ParObjectType;

   public:
    static ParObjectFactory& Instance();

    /// create a parobject from its data stream
    ParObject* Create(const std::vector<char>& data);

    /// create an element from its name (and dis type if needed)
    Teuchos::RCP<Element> Create(
        const std::string eletype, const std::string eledistype, const int id, const int owner);

    /// initialize all element types
    void InitializeElements(DRT::Discretization& dis);

    /// preevaluate elements (via element types)
    void PreEvaluate(DRT::Discretization& dis, Teuchos::ParameterList& p,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
        Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
        Teuchos::RCP<Epetra_Vector> systemvector3);

    /// setup definition of element input file lines
    void SetupElementDefinition(
        std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

    //! called upon destruction (free memory of parobject factory instance)
    virtual void Done();

   private:
    void Register(ParObjectType* object_type);

    void FinalizeRegistration();

    /// singleton.
    static ParObjectFactory* instance_;

    ParObjectFactory();

    /// list of parobject types

    /// Id to type object map. The central beast.
    std::map<int, ParObjectType*> type_map_;

    /// element name cache
    std::map<std::string, ElementType*> element_cache_;

    /// preregistered types
    std::vector<ParObjectType*> types_;

    /// element types that are actually used
    std::map<DRT::Discretization*, std::set<ElementType*>> active_elements_;

    // no copying

    ParObjectFactory(const ParObjectFactory&);
    ParObjectFactory& operator=(const ParObjectFactory&);
  };

}  // namespace DRT

#endif
