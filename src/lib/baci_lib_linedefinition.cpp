/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of one line of an input file.

\level 0

*/
/*----------------------------------------------------------------------*/


#include "baci_lib_linedefinition.H"

#include "baci_utils_demangle.H"
#include "baci_utils_exceptions.H"

#include <functional>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <utility>
#include <variant>

namespace
{
  using namespace DRT::INPUT;

  template <class T>
  std::string StringFromDataType()
  {
    if constexpr (std::is_same_v<T, int>)
      return "integer";
    else if constexpr (std::is_same_v<T, double>)
      return "double";
    else if constexpr (std::is_same_v<T, std::string>)
      return "string";
    else
    {
      // Trick to use static_assert in constexpr: make the expression dependent on template
      // parameter. In this case, the assert always fails if the type is not included in the cases
      // above.
      static_assert(!sizeof(T), "No match for given type.");
    }
  }

  /**
   * Typedef which allows to have no value in ComponentValue. This is useful to represent a tag
   * without a value.
   */
  using Empty = std::monostate;

  /**
   * Various types that can be stored in a concrete LineComponent.
   */
  using ComponentValue = std::variant<Empty, int, double, std::string, std::vector<int>,
      std::vector<double>, std::vector<std::string>, std::vector<std::pair<std::string, double>>>;


  //! Decide on behavior of the component name when reading or printing.
  enum class Behavior
  {
    // The name is printed to and read from lines.
    read_print_name,
    // The name is only used to name a component internally and not part of the input line.
    ignore_name,
  };


  /**
   * Visitor to read a ComponentValue from an input stream.
   */
  struct ReadValue
  {
    bool operator()(Empty& /*unused*/) { return true; }

    template <typename T>
    bool operator()(T& value)
    {
      in >> value;

      CheckStreamForUnparsedCharacters(StringFromDataType<T>());
      return !in.fail();
    }

    template <typename T>
    bool operator()(std::vector<T>& values)
    {
      for (auto& v : values)
      {
        in >> v;
      }

      CheckStreamForUnparsedCharacters(StringFromDataType<T>());
      return !in.fail();
    }

    template <typename T1, typename T2>
    bool operator()(std::vector<std::pair<T1, T2>>& values)
    {
      for (auto& [a, b] : values)
      {
        in >> a;
        CheckStreamForUnparsedCharacters(StringFromDataType<T1>());

        in >> b;
        CheckStreamForUnparsedCharacters(StringFromDataType<T2>());
      }

      return !in.fail();
    }

    std::istream& in;

    //! Store the name to throw better errors.
    const std::string& name;

   private:
    void CheckStreamForUnparsedCharacters(const std::string& correctDataType)
    {
      if (!(in.eof()))  // otherwise peek writes a failbit
      {
        if ((in.peek() != ' ') && (in.peek() != '/'))  // exclude whitespaces and comments
        {
          dserror(
              "Data type of component %s doesn't match the specified value. Needs to be of type %s",
              name.c_str(), correctDataType.c_str());
        }
      }
    }
  };


  /**
   * Visitor to print a ComponentValue.
   *
   * The call operator is overloaded on the type in ComponentValue since this type determines
   * whether and how printing happens.
   */
  struct PrintComponent
  {
    void operator()(const Empty& /*unused*/)
    {
      if (behavior == Behavior::read_print_name) out << name;
    }

    template <typename T>
    void operator()(const T& value)
    {
      if (behavior == Behavior::read_print_name) out << name << " ";
      out << value;
    }

    template <typename T>
    void operator()(const std::vector<T>& values)
    {
      if (behavior == Behavior::read_print_name) out << name << " ";
      std::copy(values.begin(), values.end(), std::ostream_iterator<T>(out, " "));
    }

    template <typename T1, typename T2>
    void operator()(const std::vector<std::pair<T1, T2>>& values)
    {
      if (behavior == Behavior::read_print_name) out << name << " ";
      for (const auto& [a, b] : values)
      {
        out << a << " " << b << " ";
      };
    }

    std::ostream& out;
    const Behavior& behavior;
    const std::string& name;
  };


  /**
   * A type-erased container for components of a LineDefinition.
   */
  class LineComponent
  {
   private:
    //! The type-erased interface.
    class LineComponentConcept
    {
     public:
      virtual ~LineComponentConcept() = default;

      [[nodiscard]] virtual std::unique_ptr<LineComponentConcept> Clone() const = 0;

      virtual void Print(std::ostream& stream) const = 0;

      virtual bool ReadRequired(LineDefinition& definition, std::istream& stream) = 0;

      virtual bool ReadOptional(
          LineDefinition& definition, std::string& name, std::istream& stream) = 0;

      [[nodiscard]] virtual bool IsNamed(const std::string& name) const = 0;

      [[nodiscard]] virtual const ComponentValue& Value() const = 0;
    };

    //! The wrapper for a concrete type T compatible with the interface.
    //! Forward all calls to the concrete type T.
    template <typename T>
    class LineComponentModel : public LineComponentConcept
    {
     public:
      template <typename T2,
          // prevent matching a copy or move constructor
          typename =
              std::enable_if_t<!std::is_same_v<std::decay_t<T2>, LineComponentModel<T>>, void>>
      LineComponentModel(T2&& in) : component_(std::forward<T2>(in))
      {
      }

      [[nodiscard]] std::unique_ptr<LineComponentConcept> Clone() const override
      {
        return std::make_unique<LineComponentModel<T>>(*this);
      }

      void Print(std::ostream& stream) const override { component_.Print(stream); }

      bool ReadRequired(LineDefinition& definition, std::istream& stream) override
      {
        return component_.ReadRequired(definition, stream);
      }

      bool ReadOptional(
          LineDefinition& definition, std::string& name, std::istream& stream) override
      {
        return component_.ReadOptional(definition, name, stream);
      }

      [[nodiscard]] bool IsNamed(const std::string& name) const override
      {
        return component_.IsNamed(name);
      }

      [[nodiscard]] const ComponentValue& Value() const override { return component_.Value(); }

     private:
      T component_;
    };

   public:
    /**
     * Templated constructor creating a bridge from concrete type T to type-erase interface
     * LinearComponentConcept.
     */
    template <typename T,
        // prevent matching a copy or move constructor
        typename = std::enable_if_t<!std::is_same_v<std::decay_t<T>, LineComponent>, void>>
    LineComponent(T&& in) : pimpl_(std::make_unique<LineComponentModel<T>>(std::forward<T>(in)))
    {
    }

    //! Copy constructor.
    LineComponent(const LineComponent& other) : pimpl_(other.pimpl_->Clone()) {}

    //! Copy assignment.
    LineComponent& operator=(const LineComponent& other)
    {
      this->pimpl_ = other.pimpl_->Clone();
      return *this;
    }

    //! Move constructor.
    LineComponent(LineComponent&& other) = default;

    //! Move assignment.
    LineComponent& operator=(LineComponent&& other) = default;

    /// print to a dat file comment
    void Print(std::ostream& stream) const { pimpl_->Print(stream); }

    /// Try to read component from input line
    /// This function is called for required components.
    bool ReadRequired(LineDefinition& definition, std::istream& stream)
    {
      return pimpl_->ReadRequired(definition, stream);
    }

    /// Try to read component from input line
    /// This function is called for optional components.
    bool ReadOptional(LineDefinition& definition, std::string& name, std::istream& stream)
    {
      return pimpl_->ReadOptional(definition, name, stream);
    }

    /// tell if the component has the specified name tag
    [[nodiscard]] bool IsNamed(const std::string& name) const { return pimpl_->IsNamed(name); }

    [[nodiscard]] const ComponentValue& Value() const { return pimpl_->Value(); }

   private:
    std::unique_ptr<LineComponentConcept> pimpl_;
  };


  /**
   * Component that works for various primitive types and STL containers.
   */
  template <typename T>
  class GenericComponent
  {
   public:
    GenericComponent(std::string name, T value, Behavior behavior = Behavior::read_print_name)
        : name_(std::move(name)), value_(std::move(value)), behavior_(behavior)
    {
    }

    void Print(std::ostream& stream) const
    {
      std::visit(PrintComponent{stream, behavior_, name_}, value_);
    }

    bool ReadOptional(LineDefinition& definition, std::string& name, std::istream& stream)
    {
      if (name != name_) return false;
      return std::visit(ReadValue{stream, name_}, value_);
    }

    bool ReadRequired(LineDefinition& definition, std::istream& stream)
    {
      if (behavior_ == Behavior::read_print_name)
      {
        std::string name;
        stream >> name;
        if (name != name_) return false;
      }

      return std::visit(ReadValue{stream, name_}, value_);
    }

    [[nodiscard]] bool IsNamed(const std::string& name) const { return name == name_; }

    [[nodiscard]] const ComponentValue& Value() const { return value_; }

   private:
    //! Name of the component.
    std::string name_;

    //! The value stored as a variant type. This might look strange at first, given that we know the
    //! type T. However, this design enables us to query values generically through the type-erased
    //! interface LineComponent.
    ComponentValue value_{T{}};

    //! Store how the component behaves when reading or writing.
    Behavior behavior_{Behavior::read_print_name};
  };

  /// Special LineComponent to describe a string followed by a vector of values with
  /// arbitrary length. The length is defined by a function object.
  template <class T>
  class NamedVariableVectorComponent
  {
   public:
    NamedVariableVectorComponent(std::string name,
        LineDefinition::Builder::LengthDefinition lengthdef,
        Behavior behavior = Behavior::read_print_name)
        : name_(std::move(name)),
          lengthdef_(std::move(lengthdef)),
          behavior_(behavior),
          value_(std::vector<T>{})
    {
    }

    void Print(std::ostream& stream) const
    {
      if (behavior_ == Behavior::read_print_name) stream << this->name_ << " ";
      stream << "[...]";
    }

    bool ReadRequired(LineDefinition& definition, std::istream& stream)
    {
      if (behavior_ == Behavior::read_print_name)
      {
        std::string name;
        stream >> name;
        if (name != name_) return false;
      }

      return ReadOptional(definition, name_, stream);
    }

    bool ReadOptional(LineDefinition& definition, std::string& name, std::istream& stream)
    {
      // Find expected vector on line. It has to be read already!
      const std::size_t length = lengthdef_(definition);

      dsassert(std::holds_alternative<std::vector<T>>(value_), "Internal error.");
      std::get<std::vector<T>>(value_).resize(length);

      return std::visit(ReadValue{stream, name_}, value_);
    }

    [[nodiscard]] bool IsNamed(const std::string& name) const { return name == name_; }

    [[nodiscard]] const ComponentValue& Value() const { return value_; }

   private:
    std::string name_;
    LineDefinition::Builder::LengthDefinition lengthdef_;
    Behavior behavior_;
    ComponentValue value_;
  };

}  // namespace


namespace DRT::INPUT
{
  namespace INTERNAL
  {
    /**
     * Internal data used in the implementation. This type has value semantics.
     */
    class LineDefinitionImplementation
    {
     public:
      [[nodiscard]] const LineComponent* FindNamed(const std::string& name) const
      {
        if (readtailcomponents_.find(name) != readtailcomponents_.end())
        {
          const auto i = optionaltail_.find(name);
          if (i != optionaltail_.end()) return &i->second;
        }
        else
        {
          for (const auto& component : components_)
          {
            if (component.IsNamed(name))
            {
              return &component;
            }
          }
        }
        return nullptr;
      }

      //! Return if a component of given name and type T exists.
      template <typename T>
      bool HasNamed(const std::string& name)
      {
        const auto* c = FindNamed(name);
        return c && std::holds_alternative<T>(c->Value());
      }

      template <typename T>
      void TryExtract(const std::string& name, T& dst)
      {
        if (const auto* c = FindNamed(name))
        {
          if (const auto* p_value = std::get_if<T>(&c->Value()))
          {
            dst = *p_value;
          }
          else
          {
            const std::string actual = std::visit([](const auto& v)
                { return CORE::UTILS::TryDemangle(typeid(std::decay_t<decltype(v)>).name()); },
                c->Value());

            const std::string tried = CORE::UTILS::TryDemangle(typeid(T).name());

            dserror("Line component '%s' has type '%s' but queried type is '%s'.", name.c_str(),
                actual.c_str(), tried.c_str());
          }
        }
        else
          dserror("Component names '%s' not found.", name.c_str());
      }

      /// Gather all added required components.
      std::vector<LineComponent> components_;

      /// Gather all added optional components.
      std::unordered_map<std::string, LineComponent> optionaltail_;

      /// Store which optional components have been read.
      std::set<std::string> readtailcomponents_;
    };
  }  // namespace INTERNAL


  LineDefinition::LineDefinition()
      : pimpl_(std::make_unique<INTERNAL::LineDefinitionImplementation>())
  {
  }


  // The PIMPL idiom forces us to default a few special members in the implementation file.
  LineDefinition::~LineDefinition() = default;

  LineDefinition::LineDefinition(LineDefinition&&) noexcept = default;

  LineDefinition& LineDefinition::operator=(LineDefinition&&) noexcept = default;

  LineDefinition::LineDefinition(const LineDefinition& other)
      : pimpl_(std::make_unique<INTERNAL::LineDefinitionImplementation>(*other.pimpl_))
  {
  }

  LineDefinition& LineDefinition::operator=(const LineDefinition& other)
  {
    pimpl_ = std::make_unique<INTERNAL::LineDefinitionImplementation>(*other.pimpl_);
    return *this;
  }

  LineDefinition::LineDefinition(std::unique_ptr<INTERNAL::LineDefinitionImplementation>&& pimpl)
      : pimpl_(std::move(pimpl))
  {
  }



  LineDefinition::Builder::Builder()
      : pimpl_(std::make_unique<INTERNAL::LineDefinitionImplementation>())
  {
  }

  // The PIMPL idiom forces us to default this in the implementation file.
  LineDefinition::Builder::~Builder() = default;

  LineDefinition::Builder::Builder(const LineDefinition& line_definition)
      : pimpl_(std::make_unique<INTERNAL::LineDefinitionImplementation>(*line_definition.pimpl_))
  {
  }

  LineDefinition LineDefinition::Builder::Build() &&
  {
    // Steal the internal data since this operation is performed on an rvalue.
    return LineDefinition(std::move(pimpl_));
  }

  LineDefinition LineDefinition::Builder::Build() const&
  {
    // Make a copy of the implementation details
    return LineDefinition(std::make_unique<INTERNAL::LineDefinitionImplementation>(*pimpl_));
  }


  LineDefinition::Builder& LineDefinition::Builder::AddTag(std::string name)
  {
    pimpl_->components_.emplace_back(GenericComponent<Empty>{std::move(name), Empty()});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddString(std::string name)
  {
    pimpl_->components_.emplace_back(
        GenericComponent<std::string>{std::move(name), "''", Behavior::ignore_name});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddInt(std::string name)
  {
    pimpl_->components_.emplace_back(
        GenericComponent<int>{std::move(name), 0, Behavior::ignore_name});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddIntVector(std::string name, int length)
  {
    pimpl_->components_.emplace_back(
        GenericComponent(std::move(name), std::vector<int>(length), Behavior::ignore_name));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddDoubleVector(std::string name, int length)
  {
    pimpl_->components_.emplace_back(
        GenericComponent(std::move(name), std::vector<double>(length), Behavior::ignore_name));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddNamedString(std::string name)
  {
    pimpl_->components_.emplace_back(GenericComponent<std::string>{std::move(name), "''"});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddNamedInt(std::string name)
  {
    pimpl_->components_.emplace_back(GenericComponent<int>{std::move(name), 0});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddNamedIntVector(std::string name, int length)
  {
    pimpl_->components_.emplace_back(GenericComponent(std::move(name), std::vector<int>(length)));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddNamedDouble(std::string name)
  {
    pimpl_->components_.emplace_back(GenericComponent<double>{std::move(name), 0.0});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddNamedDoubleVector(
      std::string name, int length)
  {
    pimpl_->components_.emplace_back(
        GenericComponent(std::move(name), std::vector<double>(length)));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddNamedDoubleVector(
      std::string name, LengthDefinition length_definition)
  {
    pimpl_->components_.emplace_back(
        NamedVariableVectorComponent<double>(std::move(name), std::move(length_definition)));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalTag(const std::string& name)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      dserror("optional component '%s' already defined", name.c_str());
    pimpl_->optionaltail_.emplace(name, GenericComponent<Empty>{name, Empty()});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalNamedString(const std::string& name)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      dserror("optional component '%s' already defined", name.c_str());
    pimpl_->optionaltail_.emplace(name, GenericComponent<std::string>{name, "''"});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalNamedInt(const std::string& name)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      dserror("optional component '%s' already defined", name.c_str());
    pimpl_->optionaltail_.emplace(name, GenericComponent<int>{name, 0});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalNamedIntVector(
      const std::string& name, int length)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      dserror("optional component '%s' already defined", name.c_str());
    pimpl_->optionaltail_.emplace(name, GenericComponent(name, std::vector<int>(length)));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalNamedDouble(const std::string& name)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      dserror("optional component '%s' already defined", name.c_str());
    pimpl_->optionaltail_.emplace(name, GenericComponent<double>{name, 0.});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalNamedDoubleVector(
      const std::string& name, int length)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      dserror("optional component '%s' already defined", name.c_str());
    pimpl_->optionaltail_.emplace(name, GenericComponent(name, std::vector<double>(length)));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalNamedDoubleVector(
      const std::string& name, LengthDefinition lengthdef)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      dserror("optional component '%s' already defined", name.c_str());
    pimpl_->optionaltail_.emplace(
        name, NamedVariableVectorComponent<double>(name, std::move(lengthdef)));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalNamedStringVector(
      const std::string& name, int length)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      dserror("optional component '%s' already defined", name.c_str());
    pimpl_->optionaltail_.emplace(
        name, GenericComponent(name, std::vector<std::string>(length, "''")));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalNamedStringVector(
      const std::string& name, LengthDefinition lengthdef)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      dserror("optional component '%s' already defined", name.c_str());
    pimpl_->optionaltail_.emplace(
        name, NamedVariableVectorComponent<std::string>(name, std::move(lengthdef)));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalNamedPairOfStringAndDoubleVector(
      const std::string& name, LengthDefinition lengthdef)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      dserror("optional component '%s' already defined", name.c_str());
    pimpl_->optionaltail_.emplace(name,
        NamedVariableVectorComponent<std::pair<std::string, double>>(name, std::move(lengthdef)));
    return *this;
  }



  void LineDefinition::Print(std::ostream& stream) const
  {
    for (const auto& component : pimpl_->components_)
    {
      component.Print(stream);
      stream << ' ';
    }
    if (!pimpl_->optionaltail_.empty())
    {
      stream << "[ ";
      for (const auto& [_, component] : pimpl_->optionaltail_)
      {
        component.Print(stream);
        stream << ' ';
      }
      stream << "] ";
    }
  }



  bool LineDefinition::Read(std::istream& stream) { return Read(stream, nullptr); }



  bool LineDefinition::Read(std::istream& stream, const std::string* skipname)
  {
    pimpl_->readtailcomponents_.clear();
    for (auto& component : pimpl_->components_)
    {
      if (!skipname || !component.IsNamed(*skipname))
      {
        if (not component.ReadRequired(*this, stream))
        {
          return false;
        }
      }
    }

    // we expect as many optional components as are defined (or less)
    for (unsigned a = 0; a < pimpl_->optionaltail_.size(); ++a)
    {
      std::string name;
      stream >> name;
      if (not stream) break;
      if (pimpl_->readtailcomponents_.find(name) != pimpl_->readtailcomponents_.end())
      {
        // duplicated optional component
        return false;
      }
      auto i = pimpl_->optionaltail_.find(name);
      if (i == pimpl_->optionaltail_.end()) return false;
      if (not i->second.ReadOptional(*this, name, stream)) return false;
      pimpl_->readtailcomponents_.insert(name);
    }

    // check if any other unused strings except from comments and whitespaces are given
    std::string superfluousstring;
    stream >> superfluousstring;  // stream strips whitespaces

    // Check that remaining string is either empty or is a comment, i.e., starts with "//"
    return superfluousstring.empty() || (superfluousstring.rfind("//", 0) == 0);
  }



  bool LineDefinition::HaveNamed(const std::string& name) const
  {
    return pimpl_->FindNamed(name) != nullptr;
  }



  void LineDefinition::ExtractString(const std::string& name, std::string& value) const
  {
    pimpl_->TryExtract(name, value);
  }



  bool LineDefinition::HasString(const std::string& name) const
  {
    return pimpl_->HasNamed<std::string>(name);
  }



  void LineDefinition::ExtractInt(const std::string& name, int& value) const
  {
    pimpl_->TryExtract(name, value);
  }



  void LineDefinition::ExtractIntVector(const std::string& name, std::vector<int>& v) const
  {
    pimpl_->TryExtract(name, v);
  }



  void LineDefinition::ExtractDouble(const std::string& name, double& value) const
  {
    pimpl_->TryExtract(name, value);
  }



  void LineDefinition::ExtractDoubleVector(const std::string& name, std::vector<double>& v) const
  {
    pimpl_->TryExtract(name, v);
  }



  void LineDefinition::ExtractStringVector(
      const std::string& name, std::vector<std::string>& v) const
  {
    // special case: search for a string vector first but fall back to a string value
    pimpl_->TryExtract(name, v);
  }



  void LineDefinition::ExtractPairOfStringAndDoubleVector(
      const std::string& name, std::vector<std::pair<std::string, double>>& v) const
  {
    pimpl_->TryExtract(name, v);
  }



  void Lines::Print(std::ostream& stream) const
  {
    constexpr std::size_t max_line_width = 65ul;
    stream << "--";
    stream << std::string(std::max(max_line_width - sectionname_.length(), 0ul), '-');
    stream << sectionname_ << '\n';

    for (const auto& definition : definitions_)
    {
      stream << "// ";
      definition.Print(stream);
      stream << '\n';
    }
  }



  std::vector<LineDefinition> Lines::Read(DatFileReader& reader, int suffix)
  {
    std::ostringstream name;
    name << "--" << sectionname_;
    if (suffix > 0) name << suffix;

    std::vector<LineDefinition> lines;

    std::vector<const char*> section = reader.Section(name.str());
    // these are the lines of the section stored in [0],[1],...
    if (!section.empty())
    {
      for (const auto& input_line : section)
      {
        const LineDefinition& filled_line_definition = std::invoke(
            [&]()
            {
              for (auto& definition : definitions_)
              {
                std::stringstream l{input_line};
                if (definition.Read(l))
                {
                  return definition;
                }
              }
              {
                std::stringstream out;
                out << "Read failed in section " << std::quoted(name.str()) << ": line "
                    << std::quoted(input_line) << "\n";
                out << "Valid lines are:\n\n";
                std::for_each(definitions_.begin(), definitions_.end(),
                    [&](const LineDefinition& def)
                    {
                      def.Print(out);
                      out << "\n";
                    });
                dserror(out.str().c_str());
              }
            });

        lines.push_back(filled_line_definition);
      }
    }

    return lines;
  }


  LengthFromIntNamed::LengthFromIntNamed(std::string definition_name)
      : definition_name_(std::move(definition_name))
  {
  }


  std::size_t LengthFromIntNamed::operator()(const LineDefinition& already_read_line)
  {
    int length = 0;
    already_read_line.ExtractInt(definition_name_, length);
    return length;
  }
}  // namespace DRT::INPUT
