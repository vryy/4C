/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu factory class for BACI
\level 2
\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/

#ifndef MUELU_BACIFACTORYFACTORY_DECL_HPP
#define MUELU_BACIFACTORYFACTORY_DECL_HPP

#ifdef HAVE_MueLu

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_FactoryFactory.hpp"
#include "MueLu_FactoryBase.hpp"
#include "MueLu_FactoryManagerBase_fwd.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"

#include "MueLu_NodeDefinition.hpp"
#include "MueLu_ContactAFilterFactory_decl.hpp"
#include "MueLu_ContactSPAggregationFactory_decl.hpp"
#include "MueLu_SelectiveSaPFactory_decl.hpp"
#include "MueLu_IterationAFactory_decl.hpp"
#include "MueLu_MyTrilinosSmoother_decl.hpp"

namespace MueLu
{
  /*! class BaciFactoryFactory
    @brief FactoryFactory extension for special MueLu classes in BACI

   */
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal,
      class Node = KokkosSerialNode>
  class BaciFactoryFactory : public FactoryFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>
  {
#undef MUELU_BACIFACTORYFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

    typedef std::map<std::string, RCP<const FactoryBase>> FactoryMap;  // TODO: remove
    typedef std::map<std::string, RCP<FactoryManagerBase>> FactoryManagerMap;

   public:
    /*!
      @brief Parameter list parsing (factory design)

      Routine extends BuildFactory method from MueLu::FactoryFactory class for Baci specific MueLu
      factories

      @param[in] param: TODO
      @param[in] factoryMapIn: TODO
      @param[in] factoryManagersIn: TODO
     */
    virtual RCP<const FactoryBase> BuildFactory(const Teuchos::ParameterEntry& param,
        const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const
    {
      // Find factory
      std::string factoryName;
      Teuchos::ParameterList paramList;
      if (!param.isList())
      {
        factoryName = Teuchos::getValue<std::string>(param);
      }
      else
      {
        paramList = Teuchos::getValue<Teuchos::ParameterList>(param);
        factoryName = paramList.get<std::string>("factory");
      }

      // TODO: add new factories...
      if (factoryName == "ContactAFilterFactory")
        return myBuild<MueLu::ContactAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(
            paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "ContactSPAggregationFactory")
        return myBuild<
            MueLu::ContactSPAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(
            paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "SelectiveSaPFactory")
        return myBuild<MueLu::SelectiveSaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(
            paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "IterationAFactory")
        return myBuild<MueLu::IterationAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(
            paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "FilteredTrilinosSmoother")
        return myBuildSmoother(paramList, factoryMapIn, factoryManagersIn);

      if (factoryMapIn.find(factoryName) != factoryMapIn.end())
      {
        // factory handled here
        return factoryMapIn.find(factoryName)->second;
      }

      // call BuildFactory from basic FactoryFactory class
      return FactoryFactory::BuildFactory(param, factoryMapIn, factoryManagersIn);
    }

    template <class T>
    RCP<T> myBuild(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn,
        const FactoryManagerMap& factoryManagersIn) const
    {
      RCP<T> factory = rcp(new T());

      ParameterList paramListWithFactories;

      // Read the RCP<Factory> parameters of the class T
      RCP<const ParameterList> validParamList = factory->GetValidParameterList();
      for (ParameterList::ConstIterator param = validParamList->begin();
           param != validParamList->end(); ++param)
      {
        const std::string& pName = validParamList->name(param);

        if (!paramList.isParameter(pName))
        {
          // Ignore unknown parameters
          continue;
        }

        if (validParamList->isType<RCP<const FactoryBase>>(pName))
        {
          // Generate or get factory described by param
          RCP<const FactoryBase> generatingFact =
              BuildFactory(paramList.getEntry(pName), factoryMapIn, factoryManagersIn);
          paramListWithFactories.set(pName, generatingFact);
        }
        else if (validParamList->isType<RCP<const ParameterList>>(pName))
        {
          if (pName == "ParameterList")
          {
            // NOTE: we cannot use
            //     subList = sublist(rcpFromRef(paramList), pName)
            // here as that would result in sublist also being a reference to a temporary object.
            // The resulting dereferencing in the corresponding factory would then segfault
            RCP<const ParameterList> subList =
                Teuchos::sublist(rcp(new ParameterList(paramList)), pName);
            paramListWithFactories.set(pName, subList);
          }
        }
        else
        {
          paramListWithFactories.setEntry(pName, paramList.getEntry(pName));
        }
      }

      // Configure the factory
      factory->SetParameterList(paramListWithFactories);

      return factory;
    }

    //! myBuildSmoother (experimental, not really working yet)
    RCP<FactoryBase> myBuildSmoother(const Teuchos::ParameterList& paramList,
        const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const
    {
      TEUCHOS_TEST_FOR_EXCEPTION(
          paramList.get<std::string>("factory") != "FilteredTrilinosSmoother",
          Exceptions::RuntimeError, "");

      // Is it true? TEUCHOS_TEST_FOR_EXCEPTION(!paramList.isParameter("type"),
      // Exceptions::RuntimeError, "TrilinosSmoother: parameter 'type' is mandatory"); type="" is
      // default in TrilinosSmoother, but what happen then?

      std::string type = "";
      if (paramList.isParameter("type")) type = paramList.get<std::string>("type");
      int overlap = 0;
      if (paramList.isParameter("overlap")) overlap = paramList.get<int>("overlap");
      // std::string verbose;         if(paramList.isParameter("verbose"))       verbose     =
      // paramList.get<std::string>("verbose");
      Teuchos::ParameterList params;
      if (paramList.isParameter("ParameterList"))
        params = paramList.get<Teuchos::ParameterList>("ParameterList");
      std::string mapName = "";
      if (paramList.isParameter("map: name")) mapName = paramList.get<std::string>("map: name");
      std::string mapFactName = "";
      if (paramList.isParameter("map: factory"))
        mapFactName = paramList.get<std::string>("map: factory");

      Teuchos::RCP<const FactoryBase> mapFact = Teuchos::null;
      // check whether user has provided a specific name for the MapFactory
      if (mapFactName == "NoFactory")
      {
        mapFact = MueLu::NoFactory::getRCP();
      }
      else if (mapFactName != "null")
      {
        mapFact = factoryMapIn.at(mapFactName);
      }

      std::string AfactName = "";
      if (paramList.isParameter("A: factory")) AfactName = paramList.get<std::string>("A: factory");
      Teuchos::RCP<const FactoryBase> AFact = Teuchos::null;
      if (AfactName == "NoFactory")
      {
        AFact = MueLu::NoFactory::getRCP();
      }
      else if (AfactName != "null")
      {
        std::cout << "AfactName = " << AfactName << std::endl;
        // if (factoryMapIn.find(AfactName) != std::map::end) {
        AFact = factoryMapIn.at(AfactName);
        std::cout << "found factory " << AFact << std::endl;
        //}
      }

      Teuchos::RCP<SmootherPrototype> smooProto =
          Teuchos::rcp(new MueLu::MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
              mapName, mapFact, type, params, overlap, AFact));
      return rcp(new SmootherFactory(smooProto));
    }

  };  // end class
}  // namespace MueLu

#define MUELU_BACIFACTORYFACTORY_SHORT
#endif /* HAVE_MueLu */
#endif /* MUELU_BACIFACTORYFACTORY_DECL_HPP */
