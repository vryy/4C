/*!----------------------------------------------------------------------
\file fluid_ele_immersed_base.H

\brief base class for all immersed fluid elements

<pre>
Maintainer:  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 - 15240
</pre>
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

  virtual DRT::ParObject* Create( const std::vector<char> & data) = 0;

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner ) = 0;

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions ) = 0;


private:

  //static FluidTypeImmersed instance_;

}; //class FluidTypeImmersed


class FluidImmersedBase : public virtual Fluid
{
public:

  //@}
  //! @name constructors and destructors and related methods

  /*!
  \brief standard constructor
  */
  FluidImmersedBase(
    int id,    ///< A unique global id
    int owner  ///< ???
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
  virtual DRT::Element* Clone() const {dserror("not implemented in base class"); return NULL;};

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { dserror("not implemented in base class"); return -1234; }

  virtual void SetIsImmersed(int isimmersed){dserror("not implemented in base class");};

  virtual void SetIsImmersedBoundary(int isimmersedboundary){dserror("not implemented in base class");};

  virtual void SetHasProjectedDirichlet(int has_projected_dirichletvalues){dserror("not implemented in base class");};

  virtual int IsImmersed(){dserror("not implemented in base class"); return -1234;};

  virtual int IsImmersedBoundary(){dserror("not implemented in base class"); return -1234;};

  virtual int HasProjectedDirichlet(){dserror("not implemented in base class"); return -1234;};

  virtual void VisIsImmersed(std::map<std::string,int>& names)
  {
    dserror("not implemented in base class");
    return; }

  virtual void VisIsImmersedBoundary(std::map<std::string,int>& names)
  {
    dserror("not implemented in base class");
    return; }

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
  virtual void Pack(DRT::PackBuffer& data) const {dserror("not implemented in base class");};

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element
  */
  virtual void Unpack(const std::vector<char>& data){dserror("not implemented in base class");};

  /*!
  \brief Destructor

  */
  ~FluidImmersedBase()
  {
    return;
  }


}; // class FluidImmersedBase
}
}

#endif
