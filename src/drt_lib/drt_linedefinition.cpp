/*----------------------------------------------------------------------*/
/*!
\file drt_linedefinition.cpp

\brief Definition of one line of an input file.

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/



#include "drt_linedefinition.H"
#include "drt_dserror.H"
#include <iterator>


namespace DRT
{
namespace INPUT
{

  /// line component to describe a single string
  class TagComponent : public LineComponent
  {
  public:
    TagComponent(std::string name) : name_(name) {}

    virtual ~TagComponent() {}

    virtual LineComponent* Clone() { return new TagComponent(name_); }

    virtual void Print(std::ostream& stream) { stream << name_; }

    virtual bool Read(LineDefinition& definition, std::istream& stream);

    virtual bool Read(LineDefinition& definition, std::string& name, std::istream& stream) ;

    virtual bool Read(std::istream& stream);

    virtual bool IsNamed(std::string name) { return name==name_; }

  private:
    std::string name_;
  };

  /// line component to describe a string followed by some value
  template <class type>
  class NamedComponent : public LineComponent
  {
  public:
    NamedComponent(std::string name, type value)
      : name_(name), value_(value) {}

    virtual ~NamedComponent() {}

    virtual LineComponent* Clone() { return new NamedComponent(name_,value_); }

    virtual void Print(std::ostream& stream) { stream << name_ << " " << value_; }

    virtual bool Read(LineDefinition& definition, std::istream& stream);

    virtual bool Read(LineDefinition& definition, std::string& name, std::istream& stream) ;

    virtual bool Read(std::istream& stream);

    virtual bool IsNamed(std::string name) { return name==name_; }

    type Value() const { return value_; }

  protected:
    std::string name_;
    type value_;
  };

  /// line component to describe some value without a (visible) string
  template <class type>
  class UnnamedComponent : public NamedComponent<type>
  {
  public:
    UnnamedComponent(std::string name, type value)
      : NamedComponent<type>(name,value) {}

    virtual ~UnnamedComponent() {}

    virtual LineComponent* Clone() { return new UnnamedComponent(this->name_,this->value_); }

    virtual void Print(std::ostream& stream) { stream << this->value_; }

    virtual bool Read(LineDefinition& definition, std::istream& stream)
    {
      return NamedComponent<type>::Read(stream);
    }
  };

//  class SpecialFunctionComponent : public LineComponent
//  {
//    SpecialFunctionComponent(std::string name, int length) : name_(name), values_(length) {}
//
//    SpecialFunctionComponent(std::string name, const std::vector<type>& values)
//      : name_(name), values_(values.begin(),values.end()) {}
//
//    virtual ~SpecialFunctionComponent() {}
//
//    virtual LineComponent* Clone() { return new NamedVectorComponent<type>(name_,values_); }
//
//    virtual void Print(std::ostream& stream)
//    {
//      stream << name_ << " ";
//      std::copy(values_.begin(), values_.end(), std::ostream_iterator<type>(stream, " "));
//    }
//
//    virtual bool Read(LineDefinition& definition, std::istream& stream);
//
//    virtual bool Read(LineDefinition& definition, std::string& name, std::istream& stream) ;
//
//    virtual bool Read(std::istream& stream);
//
//    virtual bool IsNamed(std::string name) { return name==name_; }
//
//    std::pair<std::vector<double>, std::vector<std::string> > Value() const { return std::make_pair(slices_,expressions_); }
//
//  protected:
//    std::string name_;
//    std::vector<double> slices_;
//    std::vector<std::string> expressions_;
//  };


  /// line component to describe a string followed by a vector of values
  template <class type>
  class NamedVectorComponent : public LineComponent
  {
  public:
    NamedVectorComponent(std::string name, int length) : name_(name), values_(length) {}

    NamedVectorComponent(std::string name, const std::vector<type>& values)
      : name_(name), values_(values.begin(),values.end()) {}

    virtual ~NamedVectorComponent() {}

    virtual LineComponent* Clone() { return new NamedVectorComponent<type>(name_,values_); }

    virtual void Print(std::ostream& stream)
    {
      stream << name_ << " ";
      std::copy(values_.begin(), values_.end(), std::ostream_iterator<type>(stream, " "));
    }

    virtual bool Read(LineDefinition& definition, std::istream& stream);

    virtual bool Read(LineDefinition& definition, std::string& name, std::istream& stream) ;

    virtual bool Read(std::istream& stream);

    virtual bool IsNamed(std::string name) { return name==name_; }

    std::vector<type> Value() const { return values_; }

  protected:
    std::string name_;
    std::vector<type> values_;
  };

  /// line component to describe a string followed by a vector of values
  /// specialization for pair
  template <typename T1,typename T2>
  class NamedVectorComponent< std::pair<T1,T2> > : public LineComponent
  {
  public:
    NamedVectorComponent(std::string name, int length) : name_(name), values_(length) {}

    NamedVectorComponent(std::string name, const std::vector< std::pair<T1,T2> >& values)
      : name_(name), values_(values.begin(),values.end()) {}

    virtual ~NamedVectorComponent() {}

    virtual LineComponent* Clone() { return new NamedVectorComponent< std::pair<T1,T2> >(name_,values_); }

    virtual void Print(std::ostream& stream)
    {
      stream << name_ << " ";
      for(typename std::vector< std::pair<T1,T2> >::iterator it = values_.begin();it!=values_.end();it++)
      {
        stream << it->first << " ";
        stream << it->second << " ";
      };
    }

    virtual bool Read(LineDefinition& definition, std::istream& stream);

    virtual bool Read(LineDefinition& definition, std::string& name, std::istream& stream) ;

    virtual bool Read(std::istream& stream);

    virtual bool IsNamed(std::string name) { return name==name_; }

    std::vector< std::pair<T1,T2> > Value() const { return values_; }

  protected:
    std::string name_;
    std::vector< std::pair<T1,T2> > values_;
  };

  /// line component to describe a vector of values without a (visible) string
  template <class type>
  class UnnamedVectorComponent : public NamedVectorComponent<type>
  {
  public:
    UnnamedVectorComponent(std::string name, int length)
      : NamedVectorComponent<type>(name, length) {}

    UnnamedVectorComponent(std::string name, const std::vector<type>& values)
      : NamedVectorComponent<type>(name, values) {}

    virtual ~UnnamedVectorComponent() {}

    virtual LineComponent* Clone() { return new UnnamedVectorComponent<type>(this->name_,this->values_); }

    virtual void Print(std::ostream& stream)
    {
      std::copy(this->values_.begin(), this->values_.end(), std::ostream_iterator<type>(stream, " "));
    }

    virtual bool Read(LineDefinition& definition, std::istream& stream)
    {
      return NamedVectorComponent<type>::Read(stream);
    }
  };

  /// line component to describe a string followed by a vector of values with
  /// arbitrary length
  template <class type>
  class NamedVariableVectorComponent : public NamedVectorComponent<type>
  {
  public:
    NamedVariableVectorComponent(std::string name, std::string lengthdef)
      : NamedVectorComponent<type>(name,0), lengthdef_(lengthdef) {}

    NamedVariableVectorComponent(std::string name, const std::vector<type>& values, std::string lengthdef)
      : NamedVectorComponent<type>(name,values), lengthdef_(lengthdef) {}

    virtual ~NamedVariableVectorComponent() {}

    virtual LineComponent* Clone() { return new NamedVariableVectorComponent<type>(this->name_,this->values_,lengthdef_); }

    virtual void Print(std::ostream& stream)
    {
      stream << this->name_ << " [...]";
      //std::copy(values_.begin(), values_.end(), std::ostream_iterator<type>(stream, " "));
    }

    virtual bool Read(LineDefinition& definition, std::string& name, std::istream& stream) ;

  private:
    std::string lengthdef_;
  };
}
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::INPUT::TagComponent::Read(DRT::INPUT::LineDefinition& definition, std::istream& stream)
{
  std::string tag;
  stream >> tag;
  return Read(definition,tag,stream);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::INPUT::TagComponent::Read(DRT::INPUT::LineDefinition& definition, std::string& name, std::istream& stream)
{
  return name==name_ and !(stream.fail());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::INPUT::TagComponent::Read(std::istream& stream)
{
  std::string tag;
  stream >> tag;
  return tag==name_ and !(stream.fail());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <class type>
bool DRT::INPUT::NamedComponent<type>::Read(DRT::INPUT::LineDefinition& definition, std::istream& stream)
{
  std::string name;
  stream >> name;

 return Read(definition,name,stream);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <class type>
bool DRT::INPUT::NamedComponent<type>::Read(DRT::INPUT::LineDefinition& definition, std::string& name, std::istream& stream)
{
  // here we do not require any whitespace
  // This is needed for "FUNCT1" and such nonsense

  if (name.length() > name_.length())
  {
    if (name.substr(0,name_.length()) != name_)
      return false;
    std::stringstream vstream(name.substr(name_.length()));
    return Read(vstream);
  }
  else
  {
    if (name != name_)
      return false;
    return Read(stream);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <class type>
bool DRT::INPUT::NamedComponent<type>::Read(std::istream& stream)
{
  stream >> value_;
  return !(stream.fail());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <class type>
bool DRT::INPUT::NamedVectorComponent<type>::Read(DRT::INPUT::LineDefinition& definition, std::istream& stream)
{
  std::string name;
  stream >> name;

  return Read(definition,name,stream);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <class type>
bool DRT::INPUT::NamedVectorComponent<type>::Read(DRT::INPUT::LineDefinition& definition, std::string& name, std::istream& stream)
{
  // here we require whitespaces between name and values

  if (name != name_)
    return false;

  return Read(stream);
}

namespace{
template <typename T, typename U>
struct types_are_equal
{
  static const bool value = false;
};

template <typename T>
struct types_are_equal<T,T>
{
  static const bool value = true;
};
}
//  && types_are_equal<type,double>::value

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <class type>
bool DRT::INPUT::NamedVectorComponent<type>::Read(std::istream& stream)
{
  for (unsigned i=0; i<values_.size(); ++i)
  {
    stream >> values_[i];
  }

  return !(stream.fail());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<typename T1,typename T2>
bool DRT::INPUT::NamedVectorComponent< std::pair<T1,T2> >::Read(DRT::INPUT::LineDefinition& definition, std::istream& stream)
{
  std::string name;
  stream >> name;

  return Read(definition,name,stream);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<typename T1,typename T2>
bool DRT::INPUT::NamedVectorComponent< std::pair<T1,T2> >::Read(
    DRT::INPUT::LineDefinition& definition,
    std::string& name,
    std::istream& stream)
{
  // here we require whitespaces between name and values

  if (name != name_)
    return false;

  return Read(stream);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<typename T1,typename T2>
bool DRT::INPUT::NamedVectorComponent< std::pair<T1,T2> >::Read(std::istream& stream)
{
  for (unsigned i=0; i<values_.size(); ++i)
  {
    stream >> values_[i].first;
    // here we require whitespaces between name and values
    stream >> values_[i].second;
  }

  return !(stream.fail());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <class type>
bool DRT::INPUT::NamedVariableVectorComponent<type>::Read(
    DRT::INPUT::LineDefinition& definition,
    std::string& name,
    std::istream& stream)
{
  // Find expected vector on line. It has to be read already!

  int length=0;
  definition.ExtractInt(lengthdef_,length);

  // for multifunction variables the string vector has NUMPOINTS-1 component
  if (name.compare("DESCRIPTION") == 0)
  {
    length = length - 1;
  }

  this->values_.resize(length);

  return NamedVectorComponent<type>::Read(definition,name,stream);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition::LineDefinition()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition::~LineDefinition()
{
  Clear();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition::LineDefinition(const DRT::INPUT::LineDefinition& other)
{
  for (unsigned i=0; i<other.components_.size(); ++i)
  {
    components_.push_back(other.components_[i]->Clone());
  }
  for (std::map<std::string,LineComponent*>::const_iterator i=other.optionaltail_.begin();
       i!=other.optionaltail_.end();
       ++i)
  {
    optionaltail_[i->first] = i->second->Clone();
  }
  std::copy(other.readtailcomponents_.begin(),
            other.readtailcomponents_.end(),
            std::inserter(readtailcomponents_,
                          readtailcomponents_.begin()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::operator=(const DRT::INPUT::LineDefinition& other)
{
  Clear();
  for (unsigned i=0; i<other.components_.size(); ++i)
  {
    components_.push_back(other.components_[i]->Clone());
  }
  for (std::map<std::string,LineComponent*>::const_iterator i=other.optionaltail_.begin();
       i!=other.optionaltail_.end();
       ++i)
  {
    optionaltail_[i->first] = i->second->Clone();
  }
  std::copy(other.readtailcomponents_.begin(),
            other.readtailcomponents_.end(),
            std::inserter(readtailcomponents_,
                          readtailcomponents_.begin()));
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::INPUT::LineDefinition> DRT::INPUT::LineDefinition::Clone()
{
  return Teuchos::rcp(new LineDefinition(*this));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::LineDefinition::Clear()
{
  for (unsigned i=0; i<components_.size(); ++i)
  {
    delete components_[i];
  }
  components_.resize(0);
  for (std::map<std::string,LineComponent*>::iterator i=optionaltail_.begin();
       i!=optionaltail_.end();
       ++i)
  {
    delete i->second;
  }
  optionaltail_.clear();
  readtailcomponents_.clear();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddTag(std::string name)
{
  components_.push_back(new TagComponent(name));
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddString(std::string name)
{
  components_.push_back(new UnnamedComponent<std::string>(name,"''"));
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddInt(std::string name)
{
  components_.push_back(new UnnamedComponent<int>(name,0));
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddIntVector(std::string name, int length)
{
  components_.push_back(new UnnamedVectorComponent<int>(name,length));
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddDoubleVector(std::string name, int length)
{
  components_.push_back(new UnnamedVectorComponent<double>(name,length));
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddNamedString(std::string name)
{
  components_.push_back(new NamedComponent<std::string>(name,"''"));
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddNamedInt(std::string name)
{
  components_.push_back(new NamedComponent<int>(name,0));
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddNamedIntVector(std::string name, int length)
{
  components_.push_back(new NamedVectorComponent<int>(name,length));
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddNamedDouble(std::string name)
{
  components_.push_back(new NamedComponent<double>(name,0.0));
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddNamedDoubleVector(std::string name, int length)
{
  components_.push_back(new NamedVectorComponent<double>(name,length));
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddNamedDoubleVector(std::string name, std::string lengthdef)
{
  components_.push_back(new NamedVariableVectorComponent<double>(name,lengthdef));
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddOptionalTag(std::string name)
{
  if (optionaltail_.find(name)!=optionaltail_.end())
    dserror("optional component '%s' already defined",name.c_str());
  optionaltail_[name] = new TagComponent(name);
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddOptionalNamedString(std::string name)
{
  if (optionaltail_.find(name)!=optionaltail_.end())
    dserror("optional component '%s' already defined",name.c_str());
  optionaltail_[name] = new NamedComponent<std::string>(name,"''");
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddOptionalNamedInt(std::string name)
{
  if (optionaltail_.find(name)!=optionaltail_.end())
    dserror("optional component '%s' already defined",name.c_str());
  optionaltail_[name] = new NamedComponent<int>(name,0);
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddOptionalNamedIntVector(std::string name, int length)
{
  if (optionaltail_.find(name)!=optionaltail_.end())
    dserror("optional component '%s' already defined",name.c_str());
  optionaltail_[name] = new NamedVectorComponent<int>(name,length);
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddOptionalNamedDouble(std::string name)
{
  if (optionaltail_.find(name)!=optionaltail_.end())
    dserror("optional component '%s' already defined",name.c_str());
  optionaltail_[name] = new NamedComponent<double>(name,0);
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddOptionalNamedDoubleVector(std::string name, int length)
{
  if (optionaltail_.find(name)!=optionaltail_.end())
    dserror("optional component '%s' already defined",name.c_str());
  optionaltail_[name] = new NamedVectorComponent<double>(name,length);
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddOptionalNamedDoubleVector(std::string name, std::string lengthdef)
{
  if (optionaltail_.find(name)!=optionaltail_.end())
    dserror("optional component '%s' already defined",name.c_str());
  optionaltail_[name] = new NamedVariableVectorComponent<double>(name,lengthdef);
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddOptionalNamedStringVector(std::string name, int length)
{
  if (optionaltail_.find(name)!=optionaltail_.end())
    dserror("optional component '%s' already defined",name.c_str());
  optionaltail_[name] = new NamedVectorComponent<std::string>(name,length);
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddOptionalNamedStringVector(std::string name, std::string lengthdef)
{
  if (optionaltail_.find(name)!=optionaltail_.end())
    dserror("optional component '%s' already defined",name.c_str());
  optionaltail_[name] = new NamedVariableVectorComponent<std::string>(name,lengthdef);
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition& DRT::INPUT::LineDefinition::AddOptionalNamedPairOfStringAndDoubleVector(
    std::string name, std::string lengthdef)
{
  if (optionaltail_.find(name)!=optionaltail_.end())
    dserror("optional component '%s' already defined",name.c_str());
  optionaltail_[name] = new NamedVariableVectorComponent< std::pair<std::string,double> >(name,lengthdef);
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::LineDefinition::Print(std::ostream& stream)
{
  for (unsigned i=0; i<components_.size(); ++i)
  {
    components_[i]->Print(stream);
    stream << ' ';
  }
  if (optionaltail_.size()>0)
  {
    stream << "[ ";
    for (std::map<std::string,LineComponent*>::iterator i=optionaltail_.begin();
         i!=optionaltail_.end();
         ++i)
    {
      i->second->Print(stream);
      stream << ' ';
    }
    stream << "] ";
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::INPUT::LineDefinition::Read(std::istream& stream)
{
  return Read(stream, 0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::INPUT::LineDefinition::Read(std::istream& stream, const std::string* skipname)
{
  readtailcomponents_.clear();
  for (unsigned i=0; i<components_.size(); ++i)
  {
    if ( !skipname || !components_[i]->IsNamed(*skipname))
    {
      if (not components_[i]->Read(*this,stream))
      {
        return false;
      }
    }
  }

  // we expect as much optional components as are defined (or less)
  for (unsigned a=0; a<optionaltail_.size(); ++a)
  {
    std::string name;
    stream >> name;
    if (not stream)
      break;
    if (readtailcomponents_.find(name)!=readtailcomponents_.end())
    {
      // duplicated optional component
      return false;
    }
    std::map<std::string,LineComponent*>::iterator i = optionaltail_.find(name);
    if (i==optionaltail_.end())
      return false;
    if (not i->second->Read(*this,name,stream))
      return false;
    readtailcomponents_.insert(name);
  }
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::INPUT::LineDefinition::HaveNamed(std::string name)
{
  return FindNamed(name)!=NULL;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::LineDefinition::ExtractString(std::string name, std::string& value)
{
  NamedComponent<std::string>* c = dynamic_cast<NamedComponent<std::string>*>(FindNamed(name));
  if (c!=NULL)
  {
    value = c->Value();
    return;
  }
  dserror("string '%s' not found", name.c_str());
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::INPUT::LineDefinition::FindString(std::string name)
{
  NamedComponent<std::string>* c = dynamic_cast<NamedComponent<std::string>*>(FindNamed(name));
  if (c!=NULL)
  {
    return true;
  }
  else
  {
    return false;
  }
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::LineDefinition::ExtractInt(std::string name, int& value)
{
  NamedComponent<int>* c = dynamic_cast<NamedComponent<int>*>(FindNamed(name));
  if (c!=NULL)
    value = c->Value();
  {
    return ;
  }
  dserror("int '%s' not found", name.c_str());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::LineDefinition::ExtractIntVector(std::string name, std::vector<int>& v)
{
  NamedVectorComponent<int>* c = dynamic_cast<NamedVectorComponent<int>*>(FindNamed(name));
  if (c!=NULL)
  {
    v = c->Value();
    return;
  }
  dserror("int vector '%s' not found", name.c_str());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::LineDefinition::ExtractDouble(std::string name, double& value)
{
  NamedComponent<double>* c = dynamic_cast<NamedComponent<double>*>(FindNamed(name));
  if (c!=NULL)
  {
    value = c->Value();
    return;
  }
  dserror("double '%s' not found", name.c_str());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::LineDefinition::ExtractDoubleVector(std::string name, std::vector<double>& v)
{
  NamedVectorComponent<double>* c = dynamic_cast<NamedVectorComponent<double>*>(FindNamed(name));
  if (c!=NULL)
  {
    v = c->Value();
    return;
  }
  dserror("double vector '%s' not found", name.c_str());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::LineDefinition::ExtractStringVector(std::string name, std::vector<std::string>& v)
{
  NamedVectorComponent<std::string>* c = dynamic_cast<NamedVectorComponent<std::string>*>(FindNamed(name));
  if (c!=NULL)
  {
    v = c->Value();
    return;
  }
  else
  {
    NamedComponent<std::string>* c1 = dynamic_cast<NamedComponent<std::string>*>(FindNamed(name));
    if (c1!=NULL)
    {
      v.push_back(c1->Value());
      return;
    }
  }
  dserror("string vector '%s' not found", name.c_str());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::LineDefinition::ExtractPairOfStringAndDoubleVector(std::string name, std::vector< std::pair<std::string,double> >& v)
{
  NamedVectorComponent< std::pair<std::string,double> >* c = dynamic_cast<NamedVectorComponent< std::pair<std::string,double> >*>(FindNamed(name));
  if (c!=NULL)
  {
    v = c->Value();
    return;
  }
  dserror("double vector '%s' not found", name.c_str());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineComponent* DRT::INPUT::LineDefinition::FindNamed(std::string name)
{
  if (readtailcomponents_.find(name)!=readtailcomponents_.end())
  {
    std::map<std::string,LineComponent*>::iterator i = optionaltail_.find(name);
    if (i!=optionaltail_.end())
      return i->second;
  }
  else
  {
    for (unsigned i=0; i<components_.size(); ++i)
    {
      if (components_[i]->IsNamed(name))
      {
        return components_[i];
      }
    }
  }
  return NULL;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::Lines::Print(std::ostream& stream)
{
  std::string blue2light = "";
  std::string bluelight = "";
  std::string redlight = "";
  std::string yellowlight = "";
  std::string greenlight = "";
  std::string magentalight = "";
  std::string endcolor = "";

  unsigned l = sectionname_.length();
  stream << redlight << "--";
  for (int i=0; i<std::max<int>(65-l,0); ++i) stream << '-';
  stream << greenlight << sectionname_ << endcolor << '\n';

  for (unsigned i=0; i<definitions_.size(); ++i)
  {
    stream << "// ";
    definitions_[i].Print(stream);
    stream << '\n';
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition> >
DRT::INPUT::Lines::Read(DatFileReader& reader, int suffix)
{
  std::ostringstream name;
  name << "--" << sectionname_;
  if (suffix>0)
    name << suffix;

  std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition> > lines;

  std::vector<const char*> section = reader.Section(name.str());
  if (section.size()>0)
  {
    for (std::vector<const char*>::iterator i=section.begin();
         i!=section.end();
         ++i)
    {
      Teuchos::RCP<DRT::INPUT::LineDefinition> line = Read(*i);
      if (line==Teuchos::null)
      {
        dserror("read failed in section '%s': line '%s'",
                name.str().c_str(), *i);
      }
      lines.push_back(line);
    }
  }

  return lines;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::INPUT::LineDefinition>
DRT::INPUT::Lines::Read(const char* line)
{
  for (unsigned i=0; i<definitions_.size(); ++i)
  {
    Teuchos::RCP<std::stringstream> l = Teuchos::rcp(new std::stringstream(line));
    if (definitions_[i].Read(*l))
    {
      return definitions_[i].Clone();
    }
  }
  return Teuchos::null;
}
