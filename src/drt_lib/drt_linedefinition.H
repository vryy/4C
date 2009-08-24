#ifndef DRT_LINEDEFINITION_H
#define DRT_LINEDEFINITION_H

#include <iostream>
#include <string>
#include <vector>

#include "drt_inputreader.H"

namespace DRT
{
namespace INPUT
{

  class LineDefinition;
  class Lines;

  class LineComponent
  {
  public:
    virtual ~LineComponent() {}

    virtual LineComponent* Clone() = 0;

    virtual void Print(std::ostream& stream) = 0;

    virtual bool Read(LineDefinition& definition, std::istream& stream) = 0;

    virtual bool IsNamed(std::string name) { return false; }
  };


  class LineDefinition
  {
  public:
    LineDefinition();
    ~LineDefinition();

    LineDefinition(const LineDefinition& other);

    LineDefinition& operator=(const LineDefinition& other);

    Teuchos::RCP<DRT::INPUT::LineDefinition> Clone();

    LineDefinition& AddTag(std::string name);
    LineDefinition& AddNamedString(std::string name);
    LineDefinition& AddNamedInt(std::string name);
    LineDefinition& AddNamedDouble(std::string name);
    LineDefinition& AddNamedDoubleVector(std::string name, int length);
    LineDefinition& AddNamedDoubleVector(std::string name, std::string lengthdef);

    void Print(std::ostream& stream);

    bool Read(std::istream& stream);

    bool HaveNamed(std::string name);

    void ExtractNamedString(std::string name, std::string& value);
    void ExtractNamedInt(std::string name, int& value);
    void ExtractNamedDouble(std::string name, double& value);
    void ExtractNamedDoubleVector(std::string name, std::vector<double>& v);

  private:

    LineComponent* FindNamed(std::string name);

    void Clear();

    std::vector<LineComponent*> components_;
  };


  class Lines
  {
  public:

    explicit Lines(std::string sectionname) : sectionname_(sectionname) {}

    void Add(const LineDefinition& ld) { definitions_.push_back(ld); }

    void Print(std::ostream& stream);

    std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition> >
    Read(const DatFileReader& reader, int suffix);

  private:

    Teuchos::RCP<DRT::INPUT::LineDefinition> Read(const char* line);

    std::string sectionname_;
    std::vector<LineDefinition> definitions_;
  };

}
}

#endif
