/*----------------------------------------------------------------------*/
/*!
 \file poro_utils_clonestrategy.H

 \brief strategy to clone porofluid from porous solid

\level 2

\maintainer Ager Christoph
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
 *----------------------------------------------------------------------*/


#ifndef PORO_CLONESTRATEGY_H_
#define PORO_CLONESTRATEGY_H_

#include "../drt_ssi/ssi_clonestrategy.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace POROELAST
{
  /// POROELAST::UTILS: Random stuff that might be helpful when dealing with poroelasticity problems
  namespace UTILS
  {
    //! \brief implementation of special clone strategy for automatic generation
    //!        of fluid discretization from a given structure discretization
    class PoroelastCloneStrategy {
    public:

      /// constructor
      explicit PoroelastCloneStrategy() {}
      /// destructor
      virtual ~PoroelastCloneStrategy() {}

      /// returns conditions names to be copied (source and target name)
      std::map<std::string, std::string> ConditionsToCopy();

    protected:

      /// determine element type string and whether element is copied or not
      bool DetermineEleType(
          DRT::Element*             actele,
          const bool                ismyele,
          std::vector<std::string>& eletype);
      /// set element-specific data (material etc.)
      virtual void SetElementData(
          Teuchos::RCP<DRT::Element> newele,
          DRT::Element*              oldele,
          const int                  matid,
          const bool                 isnurbs);

      /// check for correct material
      void CheckMaterialType(const int matid);

    }; // class PoroelastCloneStrategy


    //! \brief implementation of special clone strategy for automatic generation
    //!        of scatra discretization from a given structure discretization for porous media
    class PoroScatraCloneStrategy: public SSI::ScatraStructureCloneStrategy {
    public:

      typedef SSI::ScatraStructureCloneStrategy my;

      /// constructor
      explicit PoroScatraCloneStrategy() {}
      /// destructor
      virtual ~PoroScatraCloneStrategy() {}

      //! return SCATRA::ImplType of the element
      INPAR::SCATRA::ImplType GetImplType(
          DRT::Element*                ele        //! element whose SCATRA::ImplType shall be determined
          );

    protected:

      /// determine element type string and whether element is copied or not
      bool DetermineEleType(
          DRT::Element* actele,
          const bool ismyele,
          std::vector<std::string>& eletype);

      //! provide cloned element with element specific data (material etc.)
      void SetElementData(
          Teuchos::RCP<DRT::Element>   newele,   //! current cloned element on target discretization
          DRT::Element*                oldele,   //! current element on source discretization
          const int                    matid,    //! material of cloned element
          const bool                   isnurbs   //! nurbs flag
          );

      /// returns conditions names to be copied (source and target name)
      std::map<std::string,std::string> ConditionsToCopy();
    }; // class PoroScatraCloneStrategy


    //! \brief implementation of special clone strategy for automatic generation
    //!        of immersed fluid discretization from a given structure discretization for porous media
    class PoroelastImmersedCloneStrategy: public POROELAST::UTILS::PoroelastCloneStrategy {
    public:

      /// constructor
      explicit PoroelastImmersedCloneStrategy() {}
      /// destructor
      virtual ~PoroelastImmersedCloneStrategy() {}

    protected:

      /// determine element type string and whether element is copied or not
      virtual bool DetermineEleType(
          DRT::Element* actele,
          const bool ismyele,
          std::vector<std::string>& eletype);

      /// set element-specific data (material etc.)
      virtual void SetElementData(
          Teuchos::RCP<DRT::Element> newele,
          DRT::Element*              oldele,
          const int                  matid,
          const bool                 isnurbs);

    }; // class PoroElastImmersedCloneStrategy

  } // namespace UTILS
} // namespace POROELAST

#endif /* PORO_CLONESTRATEGY_H_ */
