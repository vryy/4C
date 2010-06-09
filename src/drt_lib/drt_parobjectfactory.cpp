
#include "drt_dserror.H"
#include "drt_element.H"
#include "drt_parobject.H"
#include "drt_parobjectfactory.H"
#include "drt_elementtype.H"


DRT::ParObjectType::ParObjectType()
  : objectid_( -1 )
{
  // Register object with incomplete type information. The derived class is
  // not yet known at this point. Thus ParObjectFactory needs to do a lazy
  // registration.
  DRT::ParObjectFactory::Instance().Register( this );
}


int DRT::ParObjectType::UniqueParObjectId()
{
  if ( objectid_==0 )
  {
    ParObjectFactory::Instance().FinalizeRegistration();
  }
  return objectid_;
}


std::string DRT::ParObjectType::Name() const
{
#if 0
  const type_info & ti = typeid( *this );
  return ti.name();
#else
  dserror( "ParObjectType without name" );
  return "";
#endif
}


DRT::ParObjectFactory * DRT::ParObjectFactory::instance_;


DRT::ParObjectFactory::ParObjectFactory()
{

}


DRT::ParObjectFactory& DRT::ParObjectFactory::Instance()
{
  if ( instance_ == NULL )
  {
    // Create on demand. This is required since the instance will be accessed
    // by ParObjectType constructors. ParObjectType are singletons as
    // well. The singleton creation order is undefined.
    instance_ = new ParObjectFactory;
  }
  return *instance_;
}


DRT::ParObject* DRT::ParObjectFactory::Create( const vector<char> & data )
{
  FinalizeRegistration();

  // mv ptr behind the size record
  const int* ptr = (const int*)(&data[0]);
  // get the type
  const int type = *ptr;

  std::map<int, ParObjectType*>::iterator i = type_map_.find( type );
  if ( i==type_map_.end() )
  {
    dserror( "object id %d undefined", type );
  }

  DRT::ParObject* o = i->second->Create( data );

  if ( o==NULL )
  {
    dserror( "failed to create object of type %d", type );
  }

  return o;
}


Teuchos::RCP<DRT::Element> DRT::ParObjectFactory::Create( const string eletype,
                                                          const string eledistype,
                                                          const int id,
                                                          const int owner )
{
  FinalizeRegistration();

  std::map<std::string, ParObjectType*>::iterator c = element_cache_.find( eletype );
  if ( c!=element_cache_.end() )
  {
    return c->second->Create( eletype, eledistype, id, owner );
  }

  for ( std::map<int, ParObjectType*>::iterator i = type_map_.begin();
        i != type_map_.end();
        ++i )
  {
    Teuchos::RCP<DRT::Element> ele = i->second->Create( eletype, eledistype, id, owner );
    if ( ele!=Teuchos::null )
    {
      element_cache_[eletype] = i->second;
      return ele;
    }
  }

  dserror( "Unknown type '%s' of finite element", eletype.c_str() );
  return Teuchos::null;
}


void DRT::ParObjectFactory::Register( ParObjectType* object_type )
{
  // Lazy registration. Just remember the pointer here.
  types_.push_back( object_type );
}


void DRT::ParObjectFactory::FinalizeRegistration()
{
  // This is called during program execution. All types are fully
  // constructed. On first call we need to create object ids.

  if ( type_map_.size()==0 and types_.size()>0 )
  {
    for ( std::vector<ParObjectType*>::iterator iter=types_.begin();
          iter!=types_.end();
          ++iter )
    {
      ParObjectType* object_type = *iter;

      std::string name = object_type->Name();
      const unsigned char* str = reinterpret_cast<const unsigned char*>( name.c_str() );

      // simple hash
      // see http://www.cse.yorku.ca/~oz/hash.html
      int hash = 5381;
      for ( int c = 0; ( c = *str ); ++str )
      {
        hash = ((hash << 5) + hash) ^ c; /* hash * 33 ^ c */
      }

      // assume no collisions for now

      std::map<int, ParObjectType*>::iterator i = type_map_.find( hash );
      if ( i!=type_map_.end() )
      {
        dserror( "object (%s,%d) already defined: (%s,%d)",
                 name.c_str(), hash, i->second->Name().c_str(), i->first );
      }

      if ( hash==0 )
      {
        dserror( "illegal hash value" );
      }

      //std::cout << "register type object: '" << name << "': " << hash << "\n";

      type_map_[hash] = object_type;
      object_type->objectid_ = hash;
    }
  }
}


void DRT::ParObjectFactory::InitializeElements( DRT::Discretization & dis )
{
  FinalizeRegistration();

  // This is element specific code. Thus we need a down cast.

  for ( std::map<int, ParObjectType*>::iterator i=type_map_.begin();
        i!=type_map_.end();
        ++i )
  {
    ParObjectType * pot = i->second;
    ElementType * eot = dynamic_cast<ElementType*>( pot );
    if ( eot!=NULL )
    {
      int err = eot->Initialize( dis );
      if (err) dserror("Element Initialize returned err=%d",err);
    }
  }
}


void DRT::ParObjectFactory::PreEvaluate(DRT::Discretization& dis,
                                        Teuchos::ParameterList& p,
                                        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                                        Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                                        Teuchos::RCP<Epetra_Vector> systemvector1,
                                        Teuchos::RCP<Epetra_Vector> systemvector2,
                                        Teuchos::RCP<Epetra_Vector> systemvector3)
{
  FinalizeRegistration();

  // This is element specific code. Thus we need a down cast.

  for ( std::map<int, ParObjectType*>::iterator i=type_map_.begin();
        i!=type_map_.end();
        ++i )
  {
    ParObjectType * pot = i->second;
    ElementType * eot = dynamic_cast<ElementType*>( pot );
    if ( eot!=NULL )
    {
      eot->PreEvaluate( dis, p, systemmatrix1, systemmatrix2,
                        systemvector1, systemvector2, systemvector3 );
    }
  }
}
