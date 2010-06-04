#ifndef DRT_PAROBJECTFACTORY_H
#define DRT_PAROBJECTFACTORY_H

#include <vector>
#include <map>
#include <string>

#include <Teuchos_RCP.hpp>

namespace DRT
{

class ParObject;
class Element;


/// abstract ParObject type information
class ParObjectType
{
  friend class ParObjectFactory;
public:

  ParObjectType();

  virtual ~ParObjectType() {}

  int UniqueParObjectId();

  virtual DRT::ParObject* Create( const vector<char> & data )
  { return NULL; }

  virtual Teuchos::RCP<Element> Create( const string eletype,
                                        const string eledistype,
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

  ParObject* Create( const vector<char> & data );

  Teuchos::RCP<Element> Create( const string eletype,
                                const string eledistype,
                                const int id,
                                const int owner );

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
