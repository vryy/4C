/*!----------------------------------------------------------------------
\file fluid_ele_immersed.H

\brief specialized immersed element used in immersed fsi

\level 1

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 - 15240

*----------------------------------------------------------------------*/
#include "fluid_ele_immersed_base.H"

#ifndef FLUID_ELE_IMMERSED_H
#define FLUID_ELE_IMMERSED_H

namespace DRT
{

namespace ELEMENTS
{

class FluidTypeImmersed : public FluidTypeImmersedBase
{
public:

  std::string Name() const { return "FluidTypeImmersed"; }

  static FluidTypeImmersed & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data);

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );


private:

  static FluidTypeImmersed instance_;

}; //class FluidTypeImmersed

class FluidImmersed : public FluidImmersedBase
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

  /*!
  \brief Each element whose nodes are all covered by an immersed dis. are set IsImmersed
  */
  virtual void SetIsImmersed(int isimmersed){is_immersed_ = isimmersed;};

  /*!
  \brief Each element which has nodes covered by the immersed dis. but at least one node that is not covered is set BoundaryIsImmersed

   An element may also be set BoundaryIsImmersed, if an integration point of the structural surface lies in this element.

  */
  virtual void SetBoundaryIsImmersed(int IsBoundaryImmersed){is_immersed_bdry_ = IsBoundaryImmersed;};

  /*!
  \brief Each element, which is either Immersed or BoundaryImmersed is also set HasProjectedDirichlet
  */
  virtual void SetHasProjectedDirichlet(int has_projected_dirichletvalues){has_projected_dirichletvalues_ = has_projected_dirichletvalues;};

  /*!
  \brief set if divergence needs to be projected to an integration point
  */
  virtual void SetIntPointHasProjectedDivergence(int gp, int intpoint_has_projected_divergence){intpoint_has_projected_divergence_->at(gp) = intpoint_has_projected_divergence;};

  /*!
  \brief store the projected divergence
  */
  virtual void StoreProjectedIntPointDivergence(int gp, double projected_intpoint_divergence){stored_projected_intpoint_divergence_->at(gp) = projected_intpoint_divergence;};

  /*!
  \brief returns true if element was set IsImmersed
  */
  virtual int IsImmersed(){return is_immersed_;};

  /*!
  \brief returns true if element was set IsBundaryImmersed
  */
  virtual int IsBoundaryImmersed(){return is_immersed_bdry_;};

  /*!
  \brief returns true if element needs to get projected Dirichlet values
  */
  virtual int HasProjectedDirichlet(){return has_projected_dirichletvalues_;};

  /*!
  \brief returns true if element needs to get projected divergence at integration point

  */
  virtual int IntPointHasProjectedDivergence(int gp){return intpoint_has_projected_divergence_->at(gp);};

  /*!
  \brief returns projected divergence at integration point
  */
  virtual double ProjectedIntPointDivergence(int gp){return stored_projected_intpoint_divergence_->at(gp);};

  /*!
  \brief returns rcp to vector containing gps with projected divergence
  */
  virtual Teuchos::RCP<std::vector<int> > GetRCPIntPointHasProjectedDivergence(){return intpoint_has_projected_divergence_;};

  /*!
  \brief returns rcp to vector containing projected divergence values
  */
  virtual Teuchos::RCP<std::vector<double> > GetRCPProjectedIntPointDivergence(){return stored_projected_intpoint_divergence_;};

  /*!
  \brief construct rcp to vector for divergence projection handling
  */
  virtual void ConstructElementRCP(int size){
    if(intpoint_has_projected_divergence_ == Teuchos::null)
      intpoint_has_projected_divergence_   = Teuchos::rcp(new std::vector<int>(size,0));
    else
      intpoint_has_projected_divergence_->resize(size,0);

    if(stored_projected_intpoint_divergence_ == Teuchos::null)
      stored_projected_intpoint_divergence_= Teuchos::rcp(new std::vector<double>(size,0.0));
    else
      stored_projected_intpoint_divergence_->resize(size,0.0);
   };

  /*!
  \brief Clean up the element rcp
  */
  virtual void DestroyElementRCP(){
   intpoint_has_projected_divergence_->clear();
   stored_projected_intpoint_divergence_->clear();

   if(intpoint_has_projected_divergence_->size() > 0)
     dserror("intpoint_has_projected_divergence_ not cleared properly");

   if(stored_projected_intpoint_divergence_->size() > 0)
     dserror("stored_projected_intpoint_divergence_ not cleared properly");
  };

  /*!
  \brief related to output of IsImmersed information
  */
  virtual void VisIsImmersed(std::map<std::string,int>& names)
  {
    names.insert(std::pair<std::string,int>("IsImmersed",1));
    return; }

  /*!
  \brief related to output of IsBoundaryImmersed information
  */
  virtual void VisIsBoundaryImmersed(std::map<std::string,int>& names)
  {
    names.insert(std::pair<std::string,int>("IsBoundaryImmersed",1));
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
    if (name=="IsBoundaryImmersed")
    {
      if ((int)data.size()<1) dserror("Size mismatch");
      data[0] = IsBoundaryImmersed();
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

  int is_immersed_; //!< true if all nodes of this element are covered by immersed discretization
  int is_immersed_bdry_; //!< true if boundary intersects this element
  int has_projected_dirichletvalues_; //!< true if Dirichlet values need to be projected to this element
  Teuchos::RCP<std::vector<int> >    intpoint_has_projected_divergence_; //!< 1000 max number of gp
  Teuchos::RCP<std::vector<double> > stored_projected_intpoint_divergence_; //!< 1000 max number of gp

};
}
}

#endif
