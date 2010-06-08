#ifndef DRT_PAROBJECTFACTORY_H
#define DRT_PAROBJECTFACTORY_H

#include <vector>
#include <map>
#include <string>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

namespace LINALG
{
  class SparseOperator;
}

namespace DRT
{

class ParObject;
class Element;
class Discretization;


/// abstract ParObject type information
class ParObjectType
{
  friend class ParObjectFactory;
public:

  ParObjectType();

  virtual ~ParObjectType() {}

  int UniqueParObjectId();

  virtual DRT::ParObject* Create( const std::vector<char> & data )
  { return NULL; }

  virtual Teuchos::RCP<Element> Create( const std::string eletype,
                                        const std::string eledistype,
                                        const int id,
                                        const int owner )
  { return Teuchos::null; }

  virtual std::string Name() const;

private:

  int objectid_;
};


/// Singleton ParObject factory class
class ParObjectFactory
{
  friend class ParObjectType;
public:

  static ParObjectFactory& Instance();

  ParObject* Create( const std::vector<char> & data );

  Teuchos::RCP<Element> Create( const std::string eletype,
                                const std::string eledistype,
                                const int id,
                                const int owner );

  void InitializeElements( DRT::Discretization & dis );

  void PreEvaluate(DRT::Discretization& dis,
                   Teuchos::ParameterList& p,
                   Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                   Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                   Teuchos::RCP<Epetra_Vector>          systemvector1,
                   Teuchos::RCP<Epetra_Vector>          systemvector2,
                   Teuchos::RCP<Epetra_Vector>          systemvector3);

private:

  void Register( ParObjectType* object_type );

  void FinalizeRegistration();

  /// singleton.
  static ParObjectFactory * instance_;

  ParObjectFactory();

  /// list of parobject types

  std::map<int, ParObjectType*> type_map_;

  std::map<std::string, ParObjectType*> element_cache_;

  std::vector<ParObjectType*> types_;

  // no copying

  ParObjectFactory( const ParObjectFactory& );
  ParObjectFactory & operator=( const ParObjectFactory& );
};

}

#endif
