/*----------------------------------------------------------------------*/
/*! \file

\brief Central type object management.

\level 0

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#include <Epetra_Comm.h>

#include "drt_dserror.H"
#include "drt_discret.H"
#include "drt_element.H"
#include "drt_parobject.H"
#include "drt_parobjectfactory.H"
#include "drt_elementtype.H"
#include "../linalg/linalg_utils_densematrix_communication.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
namespace DRT
{
  namespace
  {
    class ParObjectPreRegister
    {
     public:
      static ParObjectPreRegister* Instance()
      {
        if (instance_ == NULL)
        {
          instance_ = new ParObjectPreRegister;
        }
        return instance_;
      }

      void Register(ParObjectType* parobjecttype) { types_.push_back(parobjecttype); }

      static void Finalize()
      {
        if (instance_ != NULL)
        {
          std::vector<ParObjectType*>& types = instance_->types_;
          for (std::vector<ParObjectType*>::iterator i = types.begin(); i != types.end(); ++i)
          {
            ParObjectType* parobjecttype = *i;
            // DRT::ParObjectFactory::Instance().Register( parobjecttype );
            parobjecttype->UniqueParObjectId();
          }
          delete instance_;
          instance_ = NULL;
        }
      }

     private:
      static ParObjectPreRegister* instance_;

      /// preregistered types
      std::vector<ParObjectType*> types_;
    };

    ParObjectPreRegister* ParObjectPreRegister::instance_;
  }  // namespace
}  // namespace DRT

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::ParObjectType::ParObjectType() : objectid_(0)
{
  ParObjectPreRegister::Instance()->Register(this);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int DRT::ParObjectType::UniqueParObjectId()
{
  if (objectid_ == 0)
  {
    DRT::ParObjectFactory::Instance().Register(this);
    // ParObjectFactory::Instance().FinalizeRegistration();
  }
  return objectid_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::ParObjectFactory* DRT::ParObjectFactory::instance_;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::ParObjectFactory::ParObjectFactory() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::ParObjectFactory& DRT::ParObjectFactory::Instance()
{
  if (instance_ == NULL)
  {
    // Create on demand. This is required since the instance will be accessed
    // by ParObjectType constructors. ParObjectType are singletons as
    // well. The singleton creation order is undefined.
    instance_ = new ParObjectFactory;
  }
  return *instance_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ParObjectFactory::Done()
{
  if (instance_ != NULL) delete instance_;
  instance_ = NULL;
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::ParObject* DRT::ParObjectFactory::Create(const std::vector<char>& data)
{
  FinalizeRegistration();

  // mv ptr behind the size record
  const int* ptr = (const int*)(&data[0]);
  // get the type
  const int type = *ptr;

  std::map<int, ParObjectType*>::iterator i = type_map_.find(type);
  if (i == type_map_.end())
  {
    dserror("object id %d undefined. Have you extended DRT::ParObjectList()?", type);
  }

  DRT::ParObject* o = i->second->Create(data);

  if (o == NULL)
  {
    dserror("failed to create object of type %d", type);
  }

  return o;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ParObjectFactory::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  FinalizeRegistration();
  std::map<std::string, ElementType*>::iterator c = element_cache_.find(eletype);
  if (c != element_cache_.end())
  {
    return c->second->Create(eletype, eledistype, id, owner);
  }

  // This is element specific code. Thus we need a down cast.

  for (std::map<int, ParObjectType*>::iterator i = type_map_.begin(); i != type_map_.end(); ++i)
  {
    ParObjectType* pot = i->second;
    ElementType* eot = dynamic_cast<ElementType*>(pot);
    if (eot != NULL)
    {
      Teuchos::RCP<DRT::Element> ele = eot->Create(eletype, eledistype, id, owner);
      if (ele != Teuchos::null)
      {
        element_cache_[eletype] = eot;
        return ele;
      }
    }
  }

  dserror("Unknown type '%s' of finite element", eletype.c_str());
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ParObjectFactory::Register(ParObjectType* object_type)
{
  std::string name = object_type->Name();
  const unsigned char* str = reinterpret_cast<const unsigned char*>(name.c_str());

  // simple hash
  // see http://www.cse.yorku.ca/~oz/hash.html
  int hash = 5381;
  for (int c = 0; (c = *str); ++str)
  {
    hash = ((hash << 5) + hash) ^ c; /* hash * 33 ^ c */
  }

  // assume no collisions for now

  std::map<int, ParObjectType*>::iterator i = type_map_.find(hash);
  if (i != type_map_.end())
  {
    dserror("object (%s,%d) already defined: (%s,%d)", name.c_str(), hash,
        i->second->Name().c_str(), i->first);
  }

  if (hash == 0)
  {
    dserror("illegal hash value");
  }

  // std::cout << "register type object: '" << name << "': " << hash << "\n";

  type_map_[hash] = object_type;
  object_type->objectid_ = hash;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ParObjectFactory::FinalizeRegistration() { ParObjectPreRegister::Finalize(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ParObjectFactory::InitializeElements(DRT::Discretization& dis)
{
  FinalizeRegistration();

  // find participating element types such that only those element types are initialized

  std::set<int> ids;
  int numelements = dis.NumMyColElements();
  for (int i = 0; i < numelements; ++i)
  {
    ids.insert(dis.lColElement(i)->ElementType().UniqueParObjectId());
  }

  std::vector<int> localtypeids;
  std::vector<int> globaltypeids;

  localtypeids.reserve(ids.size());
  localtypeids.assign(ids.begin(), ids.end());

  LINALG::AllreduceVector(localtypeids, globaltypeids, dis.Comm());

  std::set<ElementType*>& ae = active_elements_[&dis];

  // This is element specific code. Thus we need a down cast.

  for (std::vector<int>::iterator i = globaltypeids.begin(); i != globaltypeids.end(); ++i)
  {
    ParObjectType* pot = type_map_[*i];
    ElementType* eot = dynamic_cast<ElementType*>(pot);
    if (eot != NULL)
    {
      ae.insert(eot);
      int err = eot->Initialize(dis);
      if (err) dserror("Element Initialize returned err=%d", err);
    }
    else
    {
      dserror("illegal element type id %d", *i);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ParObjectFactory::PreEvaluate(DRT::Discretization& dis, Teuchos::ParameterList& p,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix2, Teuchos::RCP<Epetra_Vector> systemvector1,
    Teuchos::RCP<Epetra_Vector> systemvector2, Teuchos::RCP<Epetra_Vector> systemvector3)
{
  FinalizeRegistration();

  std::set<ElementType*>& ae = active_elements_[&dis];

  for (std::set<ElementType*>::iterator i = ae.begin(); i != ae.end(); ++i)
  {
    (*i)->PreEvaluate(
        dis, p, systemmatrix1, systemmatrix2, systemvector1, systemvector2, systemvector3);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ParObjectFactory::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  FinalizeRegistration();

  // Here we want to visit all elements known to the factory.
  //
  // This is element specific code. Thus we need a down cast.

  for (std::map<int, ParObjectType*>::iterator i = type_map_.begin(); i != type_map_.end(); ++i)
  {
    ParObjectType* pot = i->second;
    ElementType* eot = dynamic_cast<ElementType*>(pot);
    if (eot != NULL)
    {
      eot->SetupElementDefinition(definitions);
    }
  }
}
