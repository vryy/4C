
#ifdef CCADISCRET

#include "drt_linedefinition.H"
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

    virtual LineComponent* Clone() { return new TagComponent(name_); }

    virtual void Print(std::ostream& stream) { stream << name_; }

    virtual bool Read(LineDefinition& definition, std::istream& stream);

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

    virtual LineComponent* Clone() { return new NamedComponent(name_,value_); }

    virtual void Print(std::ostream& stream) { stream << name_ << " " << value_; }

    virtual bool Read(LineDefinition& definition, std::istream& stream);

    virtual bool IsNamed(std::string name) { return name==name_; }

    type Value() const { return value_; }

  private:
    std::string name_;
    type value_;
  };

  /// line component to describe a string followed by a vector of values
  template <class type>
  class NamedVectorComponent : public LineComponent
  {
  public:
    NamedVectorComponent(std::string name, int length) : name_(name), values_(length) {}

    NamedVectorComponent(std::string name, const std::vector<type>& values)
      : name_(name), values_(values.begin(),values.end()) {}

    virtual LineComponent* Clone() { return new NamedVectorComponent<type>(name_,values_); }

    virtual void Print(std::ostream& stream)
    {
      stream << name_ << " ";
      std::copy(values_.begin(), values_.end(), std::ostream_iterator<type>(stream, " "));
    }

    virtual bool Read(LineDefinition& definition, std::istream& stream);

    virtual bool IsNamed(std::string name) { return name==name_; }

    std::vector<type> Value() const { return values_; }

  protected:
    std::string name_;
    std::vector<type> values_;
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

    virtual LineComponent* Clone() { return new NamedVariableVectorComponent<type>(this->name_,this->values_,lengthdef_); }

    virtual void Print(std::ostream& stream)
    {
      stream << this->name_ << " [...]";
      //std::copy(values_.begin(), values_.end(), std::ostream_iterator<type>(stream, " "));
    }

    virtual bool Read(LineDefinition& definition, std::istream& stream);

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
  return tag==name_ and stream;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <class type>
bool DRT::INPUT::NamedComponent<type>::Read(DRT::INPUT::LineDefinition& definition, std::istream& stream)
{
  std::string name;
  stream >> name;

  // here we do not require any whitespace
  // This is needed for "FUNCT1" and such nonsense

  if (name.length() > name_.length())
  {
    if (name.substr(0,name_.length()) != name_)
      return false;
    std::stringstream vstream(name.substr(name_.length()));
    vstream >> value_;
  }
  else
  {
    if (name != name_)
      return false;
    stream >> value_;
  }

  return stream;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <class type>
bool DRT::INPUT::NamedVectorComponent<type>::Read(DRT::INPUT::LineDefinition& definition, std::istream& stream)
{
  std::string name;
  stream >> name;

  // here we require whitespaces between name and values

  if (name != name_)
    return false;

  for (unsigned i=0; i<values_.size(); ++i)
  {
    stream >> values_[i];
  }

  return stream;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <class type>
bool DRT::INPUT::NamedVariableVectorComponent<type>::Read(DRT::INPUT::LineDefinition& definition, std::istream& stream)
{
  // Find expexted vector on line. It has to be read already!

  int length;
  definition.ExtractNamedInt(lengthdef_,length);
  this->values_.resize(length);

  return NamedVectorComponent<type>::Read(definition,stream);
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
  return *this;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::LineDefinition::Print(std::ostream& stream)
{
  stream << "// ";
  for (unsigned i=0; i<components_.size(); ++i)
  {
    components_[i]->Print(stream);
    stream << ' ';
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::INPUT::LineDefinition::Read(std::istream& stream)
{
  for (unsigned i=0; i<components_.size(); ++i)
  {
    if (not components_[i]->Read(*this,stream))
    {
      return false;
    }
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
void DRT::INPUT::LineDefinition::ExtractNamedString(std::string name, std::string& value)
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
void DRT::INPUT::LineDefinition::ExtractNamedInt(std::string name, int& value)
{
  NamedComponent<int>* c = dynamic_cast<NamedComponent<int>*>(FindNamed(name));
  if (c!=NULL)
  {
    value = c->Value();
    return ;
  }
  dserror("int '%s' not found", name.c_str());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::LineDefinition::ExtractNamedDouble(std::string name, double& value)
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
void DRT::INPUT::LineDefinition::ExtractNamedDoubleVector(std::string name, std::vector<double>& v)
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
DRT::INPUT::LineComponent* DRT::INPUT::LineDefinition::FindNamed(std::string name)
{
  for (unsigned i=0; i<components_.size(); ++i)
  {
    if (components_[i]->IsNamed(name))
    {
      return components_[i];
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
    definitions_[i].Print(stream);
    stream << '\n';
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition> >
DRT::INPUT::Lines::Read(const DatFileReader& reader, int suffix)
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


#endif
