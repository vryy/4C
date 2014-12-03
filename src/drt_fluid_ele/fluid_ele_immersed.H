/*!----------------------------------------------------------------------
\file fluid_ele_immersed.H

\brief base class for all immersed algorithms

<pre>
Maintainers: Andreas Rauch & Anh-Tu Vuong
             {rauch,vuong}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289--15240 / 15264
</pre>
*----------------------------------------------------------------------*/
#include "fluid_ele.H"

#ifndef FLUID_ELE_IMMERSED_H
#define FLUID_ELE_IMMERSED_H

namespace DRT
{

namespace ELEMENTS
{

class FluidTypeImmersed : public FluidType
{
public:

  std::string Name() const { return "FluidTypeImmersed"; }

  static FluidType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data);

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );


private:

  static FluidTypeImmersed instance_;

}; //class FluidTypeImmersed

class FluidImmersed : public Fluid
{
public:

  //@}
  //! @name constructors and destructors and related methods

  /*!
  \brief standard constructor
  */
  FluidImmersed(
    int id,    ///< A unique global id
    int owner  ///< ???
    );

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  FluidImmersed(const FluidImmersed& old);

  /*!
  \brief Deep copy this instance of fluid and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  virtual DRT::Element* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return FluidTypeImmersed::Instance().UniqueParObjectId(); }

  virtual void SetIsImmersed(int isimmersed){is_immersed_ = isimmersed;};

  virtual void SetIsImmersedBoundary(int isimmersedboundary){is_immersed_bdry_ = isimmersedboundary;};

  virtual void SetHasProjectedDirichlet(int has_projected_dirichletvalues){has_projected_dirichletvalues_ = has_projected_dirichletvalues;};

  virtual int IsImmersed(){return is_immersed_;};

  virtual int IsImmersedBoundary(){return is_immersed_bdry_;};

  virtual int HasProjectedDirichlet(){return has_projected_dirichletvalues_;};

  virtual void VisIsImmersed(std::map<std::string,int>& names)
  {
    names.insert(std::pair<std::string,int>("IsImmersed",1));
    return; }

  virtual void VisIsImmersedBoundary(std::map<std::string,int>& names)
  {
    names.insert(std::pair<std::string,int>("IsImmersedBoundary",1));
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
    if (name=="Owner")
    {
      if ((int)data.size()<1) dserror("Size mismatch");
      data[0] = Owner();
      return true;
    }
      if (name=="IsImmersed")
      {
        if ((int)data.size()<1) dserror("Size mismatch");
        data[0] = IsImmersed();
        return true;
      }
      if (name=="IsImmersedBoundary")
          {
            if ((int)data.size()<1) dserror("Size mismatch");
            data[0] = IsImmersedBoundary();
            return true;
          }
    if (name=="EleGId")
    {
      if ((int)data.size()<1) dserror("Size mismatch");
      data[0] = Id();
      return true;
    }
    return false;
  }

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element
  */
  virtual void Unpack(const std::vector<char>& data);

  /*!
  \brief Destructor

  */
  ~FluidImmersed()
  {
    return;
  }

private:

  int is_immersed_;
  int is_immersed_bdry_;
  int has_projected_dirichletvalues_;
  std::vector<int> matchednodes_;

};
}
}

#endif
