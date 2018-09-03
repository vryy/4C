/*!----------------------------------------------------------------------
\file fluid_ele_immersed_base.H

\brief base class for all immersed fluid elements

\level 1

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 - 15240

*----------------------------------------------------------------------*/
#include "fluid_ele.H"

#ifndef FLUID_ELE_IMMERSED_BASE_H
#define FLUID_ELE_IMMERSED_BASE_H

namespace DRT
{
  namespace ELEMENTS
  {
    class FluidTypeImmersedBase : public FluidType
    {
     public:
      /*!
      \brief Decide which element type should be created : FLUIDIMMERSED or FLUIDPOROIMMERSED
      */
      virtual DRT::ParObject* Create(const std::vector<char>& data) = 0;

      /*!
      \brief Decide which element type should be created : FLUIDIMMERSED or FLUIDPOROIMMERSED

      */
      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      /*!
      \brief Decide which element type should be created : FLUIDIMMERSED or FLUIDPOROIMMERSED
      */
      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner) = 0;

      /*!
      \brief Setup the definition line for this element
      */
      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>&
              definitions) = 0;


    };  // class FluidTypeImmersed


    class FluidImmersedBase : public virtual Fluid
    {
     public:
      //@}
      //! @name constructors and destructors and related methods

      /*!
      \brief standard constructor
      */
      FluidImmersedBase(int id,  ///< A unique global id
          int owner              ///< ???
      );

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      FluidImmersedBase(const FluidImmersedBase& old);

      /*!
      \brief Deep copy this instance of fluid and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      virtual DRT::Element* Clone() const
      {
        dserror("not implemented in base class");
        return NULL;
      };

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const
      {
        dserror("not implemented in base class");
        return -1234;
      }

      /*!
      \brief Each element whose nodes are all covered by an immersed dis. are set IsImmersed
      */
      virtual void SetIsImmersed(int isimmersed) { dserror("not implemented in base class"); };

      /*!
      \brief Each element which has nodes covered by the immersed dis. but at least one node that is
      not covered is set BoundaryIsImmersed

       An element may also be set BoundaryIsImmersed, if an integration point of the structural
      surface lies in this element.

      */
      virtual void SetBoundaryIsImmersed(int IsBoundaryImmersed)
      {
        dserror("not implemented in base class");
      };

      /*!
      \brief Each element, which is either Immersed or BoundaryImmersed is also set
      HasProjectedDirichlet
      */
      virtual void SetHasProjectedDirichlet(int has_projected_dirichletvalues)
      {
        dserror("not implemented in base class");
      };

      /*!
      \brief set if divergence needs to be projected to an integration point
      */
      virtual void SetIntPointHasProjectedDivergence(int gp, int intpoint_has_projected_divergence)
      {
        dserror("not implemented in base class");
      };

      /*!
      \brief store the projected divergence
      */
      virtual void StoreProjectedIntPointDivergence(int gp, double projected_intpoint_divergence)
      {
        dserror("not implemented in base class");
      };

      /*!
      \brief returns true if element was set IsImmersed
      */
      virtual int IsImmersed()
      {
        dserror("not implemented in base class");
        return -1234;
      };

      /*!
      \brief returns true if element was set IsBundaryImmersed
      */
      virtual int IsBoundaryImmersed()
      {
        dserror("not implemented in base class");
        return -1234;
      };

      /*!
      \brief returns true if element needs to get projected Dirichlet values
      */
      virtual int HasProjectedDirichlet()
      {
        dserror("not implemented in base class");
        return -1234;
      };

      /*!
      \brief returns true if element needs to get projected divergence at integration point

      */
      virtual int IntPointHasProjectedDivergence(int gp)
      {
        dserror("not implemented in base class");
        return -1234;
      };

      /*!
      \brief returns projected divergence at integration point
      */
      virtual double ProjectedIntPointDivergence(int gp)
      {
        dserror("not implemented in base class");
        return -1234.0;
      };

      /*!
      \brief returns rcp to vector containing gps with projected divergence
      */
      virtual Teuchos::RCP<std::vector<int>> GetRCPIntPointHasProjectedDivergence()
      {
        dserror("not implemented in base class");
        return Teuchos::null;
      };

      /*!
      \brief returns rcp to vector containing projected divergence values
      */
      virtual Teuchos::RCP<std::vector<double>> GetRCPProjectedIntPointDivergence()
      {
        dserror("not implemented in base class");
        return Teuchos::null;
      };

      /*!
      \brief related to output of IsImmersed information
      */
      virtual void VisIsImmersed(std::map<std::string, int>& names)
      {
        dserror("not implemented in base class");
        return;
      }

      /*!
      \brief related to output of IsBoundaryImmersed information
      */
      virtual void VisIsBoundaryImmersed(std::map<std::string, int>& names)
      {
        dserror("not implemented in base class");
        return;
      }

      /*!
      \brief construct rcp to vector for divergence projection handling
      */
      virtual void ConstructElementRCP(int size)
      {
        dserror("not implemented in base class");
        return;
      };

      /*!
      \brief Clean up the element rcp
      */
      virtual void DestroyElementRCP()
      {
        dserror("not implemented in base class");
        return;
      };

      /*!
      \brief Query data to be visualized using BINIO of a given name

      This method is to be overloaded by a derived method.
      The derived method is supposed to call this base method to visualize the owner of
      the element.
      If the derived method recognizes a supported data name, it shall fill it
      with corresponding data.
      If it does NOT recognizes the name, it shall do nothing.

      \warning The method must not change size of variable data

      \param name (in):   Name of data that is currently processed for visualization
      \param data (out):  data to be filled by element if it recognizes the name
      */
      virtual bool VisData(const std::string& name, std::vector<double>& data)
      {
        dserror("not implemented in base class");
        return false;
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      virtual void Pack(DRT::PackBuffer& data) const { dserror("not implemented in base class"); };

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element
      */
      virtual void Unpack(const std::vector<char>& data)
      {
        dserror("not implemented in base class");
      };

      /*!
      \brief Destructor
      */
      ~FluidImmersedBase() { return; }


    };  // class FluidImmersedBase
  }     // namespace ELEMENTS
}  // namespace DRT

#endif
