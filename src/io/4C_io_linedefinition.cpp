/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of one line of an input file.

\level 0

*/
/*----------------------------------------------------------------------*/


#include "4C_io_linedefinition.hpp"

#include "4C_inpar_container.hpp"
#include "4C_utils_demangle.hpp"
#include "4C_utils_exceptions.hpp"

#include <functional>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <utility>
#include <variant>

FOUR_C_NAMESPACE_OPEN

// Note on the namespace: while implementation details can often be put into an anonymous namespace,
// this is trickier when these details are part of a forward declared implementation class (in this
// case LineDefinitionImplementation). Compiling this file multiple times would then lead to
// different implementations of that class because the anonymous namespace is always unique.
// Compilers do not always warn about this and this is not guaranteed to be a problem (after all,
// why would we compile the class twice?). Here we encountered the warning in unity builds, so this
// is not an anonymous namespace but the INTERNAL one.
namespace INPUT::INTERNAL
{
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
          FOUR_C_THROW(
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
  class LineDefinitionComponent
  {
   private:
    //! The type-erased interface.
    class LineDefinitionComponentConcept
    {
     public:
      virtual ~LineDefinitionComponentConcept() = default;

      [[nodiscard]] virtual std::unique_ptr<LineDefinitionComponentConcept> Clone() const = 0;

      virtual void Print(std::ostream& stream) const = 0;

      virtual bool ReadRequired(
          INPAR::InputParameterContainer& container, std::istream& stream) = 0;

      virtual bool ReadOptional(
          INPAR::InputParameterContainer& container, std::string& name, std::istream& stream) = 0;

      [[nodiscard]] virtual bool IsNamed(const std::string& name) const = 0;
    };

    //! The wrapper for a concrete type T compatible with the interface.
    //! Forward all calls to the concrete type T.
    template <typename T>
    class LineDefinitionComponentModel : public LineDefinitionComponentConcept
    {
     public:
      template <typename T2,
          // prevent matching a copy or move constructor
          typename = std::enable_if_t<
              !std::is_same_v<std::decay_t<T2>, LineDefinitionComponentModel<T>>, void>>
      LineDefinitionComponentModel(T2&& in) : component_(std::forward<T2>(in))
      {
      }

      [[nodiscard]] std::unique_ptr<LineDefinitionComponentConcept> Clone() const override
      {
        return std::make_unique<LineDefinitionComponentModel<T>>(*this);
      }

      void Print(std::ostream& stream) const override { component_.Print(stream); }

      bool ReadRequired(INPAR::InputParameterContainer& container, std::istream& stream) override
      {
        return component_.ReadRequired(container, stream);
      }

      bool ReadOptional(INPAR::InputParameterContainer& container, std::string& name,
          std::istream& stream) override
      {
        return component_.ReadOptional(container, name, stream);
      }

      [[nodiscard]] bool IsNamed(const std::string& name) const override
      {
        return component_.IsNamed(name);
      }

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
        typename =
            std::enable_if_t<!std::is_same_v<std::decay_t<T>, LineDefinitionComponent>, void>>
    LineDefinitionComponent(T&& in)
        : pimpl_(std::make_unique<LineDefinitionComponentModel<T>>(std::forward<T>(in)))
    {
    }

    //! Copy constructor.
    LineDefinitionComponent(const LineDefinitionComponent& other) : pimpl_(other.pimpl_->Clone()) {}

    //! Copy assignment.
    LineDefinitionComponent& operator=(const LineDefinitionComponent& other)
    {
      this->pimpl_ = other.pimpl_->Clone();
      return *this;
    }

    //! Move constructor.
    LineDefinitionComponent(LineDefinitionComponent&& other) = default;

    //! Move assignment.
    LineDefinitionComponent& operator=(LineDefinitionComponent&& other) = default;

    /// print to a dat file comment
    void Print(std::ostream& stream) const { pimpl_->Print(stream); }

    /// Try to read component from input line
    /// This function is called for required components.
    bool ReadRequired(INPAR::InputParameterContainer& container, std::istream& stream)
    {
      return pimpl_->ReadRequired(container, stream);
    }

    /// Try to read component from input line
    /// This function is called for optional components.
    bool ReadOptional(
        INPAR::InputParameterContainer& container, std::string& name, std::istream& stream)
    {
      return pimpl_->ReadOptional(container, name, stream);
    }

    /// tell if the component has the specified name tag
    [[nodiscard]] bool IsNamed(const std::string& name) const { return pimpl_->IsNamed(name); }

   private:
    std::unique_ptr<LineDefinitionComponentConcept> pimpl_;
  };


  /**
   * Component that works for various primitive types and STL containers.
   */
  template <typename T>
  class GenericComponent
  {
   public:
    GenericComponent(std::string name, T value, Behavior behavior = Behavior::read_print_name,
        std::function<void(INPAR::InputParameterContainer&, T& value)> value_prepare = {})
        : name_(std::move(name)),
          behavior_(behavior),
          value_prepare_(std::move(value_prepare)),
          default_value_(std::move(value))
    {
    }

    void Print(std::ostream& stream) const
    {
      if (value_prepare_)
      {
        // Assume that the custom parsing logic is not printable. Print a placeholder instead.
        if (behavior_ == Behavior::read_print_name) stream << this->name_ << " ";
        stream << "[...]";
      }
      else
        PrintComponent{stream, behavior_, name_}(default_value_);
    }

    bool ReadOptional(
        INPAR::InputParameterContainer& container, std::string& name, std::istream& stream)
    {
      if (name != name_) return false;

      T value = default_value_;

      if (value_prepare_) value_prepare_(container, value);

      bool read_succeeded = ReadValue{stream, name_}(value);

      if (read_succeeded) container.Add(name, value);

      return read_succeeded;
    }

    bool ReadRequired(INPAR::InputParameterContainer& container, std::istream& stream)
    {
      if (behavior_ == Behavior::read_print_name)
      {
        std::string name;
        stream >> name;
        if (name != name_) return false;
      }
      return ReadOptional(container, name_, stream);
    }

    [[nodiscard]] bool IsNamed(const std::string& name) const { return name == name_; }

   private:
    //! Name of the component.
    std::string name_;

    //! Store how the component behaves when reading or writing.
    Behavior behavior_{Behavior::read_print_name};

    //! An optional function that determines how to parse the value. This is useful to
    //! do prelimiary work, like querying already read data.
    std::function<void(INPAR::InputParameterContainer&, T&)> value_prepare_;

    T default_value_;
  };
}  // namespace INPUT::INTERNAL


namespace INPUT
{
  namespace INTERNAL
  {
    /**
     * Internal data used in the implementation. This type has value semantics.
     */
    class LineDefinitionImplementation
    {
     public:
      [[nodiscard]] const LineDefinitionComponent* FindNamed(const std::string& name) const
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
        return container_.GetIf<T>(name) != nullptr;
      }

      template <typename T>
      void TryExtract(const std::string& name, T& dst)
      {
        dst = *container_.Get<T>(name);
      }

      /// Gather all added required components.
      std::vector<LineDefinitionComponent> components_;

      /// Gather all added optional components.
      std::unordered_map<std::string, LineDefinitionComponent> optionaltail_;

      /// Store which optional components have been read.
      std::set<std::string> readtailcomponents_;

      /// Store the read data.
      INPAR::InputParameterContainer container_;
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

  LineDefinition::Builder::Builder(LineDefinition::Builder&&) noexcept = default;

  LineDefinition::Builder& LineDefinition::Builder::operator=(
      LineDefinition::Builder&&) noexcept = default;

  LineDefinition::Builder::Builder(const LineDefinition::Builder& other)
      : pimpl_(std::make_unique<INTERNAL::LineDefinitionImplementation>(*other.pimpl_))
  {
  }

  LineDefinition::Builder& LineDefinition::Builder::operator=(const LineDefinition::Builder& other)
  {
    pimpl_ = std::make_unique<INTERNAL::LineDefinitionImplementation>(*other.pimpl_);
    return *this;
  }

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
    pimpl_->components_.emplace_back(
        INTERNAL::GenericComponent<INTERNAL::Empty>{std::move(name), INTERNAL::Empty()});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddString(std::string name)
  {
    pimpl_->components_.emplace_back(INTERNAL::GenericComponent<std::string>{
        std::move(name), "''", INTERNAL::Behavior::ignore_name});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddInt(std::string name)
  {
    pimpl_->components_.emplace_back(
        INTERNAL::GenericComponent<int>{std::move(name), 0, INTERNAL::Behavior::ignore_name});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddIntVector(std::string name, int length)
  {
    pimpl_->components_.emplace_back(INTERNAL::GenericComponent(
        std::move(name), std::vector<int>(length), INTERNAL::Behavior::ignore_name));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddDoubleVector(std::string name, int length)
  {
    pimpl_->components_.emplace_back(INTERNAL::GenericComponent(
        std::move(name), std::vector<double>(length), INTERNAL::Behavior::ignore_name));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddNamedString(std::string name)
  {
    pimpl_->components_.emplace_back(
        INTERNAL::GenericComponent<std::string>{std::move(name), "''"});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddNamedInt(std::string name)
  {
    pimpl_->components_.emplace_back(INTERNAL::GenericComponent<int>{std::move(name), 0});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddNamedIntVector(std::string name, int length)
  {
    pimpl_->components_.emplace_back(
        INTERNAL::GenericComponent(std::move(name), std::vector<int>(length)));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddNamedDouble(std::string name)
  {
    pimpl_->components_.emplace_back(INTERNAL::GenericComponent<double>{std::move(name), 0.0});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddNamedDoubleVector(
      std::string name, int length)
  {
    pimpl_->components_.emplace_back(
        INTERNAL::GenericComponent(std::move(name), std::vector<double>(length)));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddNamedDoubleVector(
      std::string name, LengthDefinition length_definition)
  {
    pimpl_->components_.emplace_back(INTERNAL::GenericComponent<std::vector<double>>(name,
        std::vector<double>{}, INTERNAL::Behavior::read_print_name,
        [lengthdef = std::move(length_definition)](
            INPAR::InputParameterContainer& linedef, std::vector<double>& values)
        {
          // Find expected vector on line. It has to be read already!
          const std::size_t length = lengthdef(linedef);
          values.resize(length);
        }));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalTag(const std::string& name)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      FOUR_C_THROW("optional component '%s' already defined", name.c_str());
    pimpl_->optionaltail_.emplace(
        name, INTERNAL::GenericComponent<INTERNAL::Empty>{name, INTERNAL::Empty()});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalNamedString(const std::string& name)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      FOUR_C_THROW("optional component '%s' already defined", name.c_str());
    pimpl_->optionaltail_.emplace(name, INTERNAL::GenericComponent<std::string>{name, "''"});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalNamedInt(const std::string& name)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      FOUR_C_THROW("optional component '%s' already defined", name.c_str());
    pimpl_->optionaltail_.emplace(name, INTERNAL::GenericComponent<int>{name, 0});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalNamedIntVector(
      const std::string& name, int length)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      FOUR_C_THROW("optional component '%s' already defined", name.c_str());
    pimpl_->optionaltail_.emplace(name, INTERNAL::GenericComponent(name, std::vector<int>(length)));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalNamedDouble(const std::string& name)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      FOUR_C_THROW("optional component '%s' already defined", name.c_str());
    pimpl_->optionaltail_.emplace(name, INTERNAL::GenericComponent<double>{name, 0.});
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalNamedDoubleVector(
      const std::string& name, int length)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      FOUR_C_THROW("optional component '%s' already defined", name.c_str());
    pimpl_->optionaltail_.emplace(
        name, INTERNAL::GenericComponent(name, std::vector<double>(length)));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalNamedDoubleVector(
      const std::string& name, LengthDefinition lengthdef)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      FOUR_C_THROW("optional component '%s' already defined", name.c_str());

    pimpl_->optionaltail_.emplace(
        name, INTERNAL::GenericComponent<std::vector<double>>(name, std::vector<double>{},
                  INTERNAL::Behavior::read_print_name,
                  [lengthdef = std::move(lengthdef)](
                      INPAR::InputParameterContainer& linedef, std::vector<double>& values)
                  {
                    // Find expected vector on line. It has to be read already!
                    const std::size_t length = lengthdef(linedef);
                    values.resize(length);
                  }));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalNamedStringVector(
      const std::string& name, int length)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      FOUR_C_THROW("optional component '%s' already defined", name.c_str());
    pimpl_->optionaltail_.emplace(
        name, INTERNAL::GenericComponent(name, std::vector<std::string>(length, "''")));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalNamedStringVector(
      const std::string& name, LengthDefinition lengthdef)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      FOUR_C_THROW("optional component '%s' already defined", name.c_str());
    using T = std::vector<std::string>;

    pimpl_->optionaltail_.emplace(name,
        INTERNAL::GenericComponent<T>(name, T{}, INTERNAL::Behavior::read_print_name,
            [lengthdef = std::move(lengthdef)](INPAR::InputParameterContainer& linedef, T& values)
            {
              // Find expected vector on line. It has to be read already!
              const std::size_t length = lengthdef(linedef);
              values.resize(length);
            }));
    return *this;
  }



  LineDefinition::Builder& LineDefinition::Builder::AddOptionalNamedPairOfStringAndDoubleVector(
      const std::string& name, LengthDefinition lengthdef)
  {
    if (pimpl_->optionaltail_.find(name) != pimpl_->optionaltail_.end())
      FOUR_C_THROW("optional component '%s' already defined", name.c_str());
    using T = std::vector<std::pair<std::string, double>>;

    pimpl_->optionaltail_.emplace(name,
        INTERNAL::GenericComponent<T>(name, T{}, INTERNAL::Behavior::read_print_name,
            [lengthdef = std::move(lengthdef)](INPAR::InputParameterContainer& linedef, T& values)
            {
              // Find expected vector on line. It has to be read already!
              const std::size_t length = lengthdef(linedef);
              values.resize(length);
            }));
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



  std::optional<INPAR::InputParameterContainer> LineDefinition::Read(std::istream& stream)
  {
    return Read(stream, nullptr);
  }



  std::optional<INPAR::InputParameterContainer> LineDefinition::Read(
      std::istream& stream, const std::string* skipname)
  {
    pimpl_->readtailcomponents_.clear();
    for (auto& component : pimpl_->components_)
    {
      if (!skipname || !component.IsNamed(*skipname))
      {
        if (not component.ReadRequired(pimpl_->container_, stream))
        {
          return std::nullopt;
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
        return std::nullopt;
      }
      auto i = pimpl_->optionaltail_.find(name);
      if (i == pimpl_->optionaltail_.end()) return std::nullopt;
      if (not i->second.ReadOptional(pimpl_->container_, name, stream)) return std::nullopt;
      pimpl_->readtailcomponents_.insert(name);
    }

    // check if any other unused strings except from comments and whitespaces are given
    std::string superfluousstring;
    stream >> superfluousstring;  // stream strips whitespaces


    // Check that remaining string is either empty or is a comment, i.e., starts with "//"
    if (!superfluousstring.empty() && (superfluousstring.rfind("//", 0) != 0))
      return std::nullopt;
    else
      return pimpl_->container_;
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
    name << sectionname_;
    if (suffix > 0) name << suffix;

    auto [parsed_lines, unparsed_lines] = ReadMatchingLines(reader, name.str(), definitions_);

    // In this function, encountering unparsed lines is an error, so construct a nice message.
    if (unparsed_lines.size() > 0)
    {
      std::stringstream out;
      out << "Read failed in section " << std::quoted(name.str()) << '\n';
      for (const auto& unparsed : unparsed_lines)
      {
        out << "  " << std::quoted(unparsed) << '\n';
      }
      out << "Valid lines are:\n";
      std::for_each(definitions_.begin(), definitions_.end(),
          [&](const LineDefinition& def)
          {
            def.Print(out);
            out << '\n';
          });
      FOUR_C_THROW(out.str().c_str());
    }

    return parsed_lines;
  }


  std::pair<std::vector<LineDefinition>, std::vector<std::string>> ReadMatchingLines(
      DatFileReader& reader, const std::string& section,
      const std::vector<LineDefinition>& possible_lines)
  {
    const std::vector<const char*> lines_in_section = reader.Section("--" + section);

    std::vector<std::string> unparsed_lines;
    std::vector<LineDefinition> parsed_lines;

    const auto process_line = [&](const std::string& input_line)
    {
      for (const auto& definition : possible_lines)
      {
        std::stringstream l{input_line};

        // Make a copy that potentially gets filled by the Read.
        auto parsed_definition = definition;
        if (parsed_definition.Read(l))
        {
          parsed_lines.emplace_back(std::move(parsed_definition));
          return;
        }
      }
      unparsed_lines.emplace_back(input_line);
    };

    for (const auto& input_line : lines_in_section)
    {
      process_line(input_line);
    }

    FOUR_C_ASSERT(
        unparsed_lines.size() + parsed_lines.size() == lines_in_section.size(), "Internal error.");

    return {parsed_lines, unparsed_lines};
  }


  LengthFromIntNamed::LengthFromIntNamed(std::string definition_name)
      : definition_name_(std::move(definition_name))
  {
  }


  std::size_t LengthFromIntNamed::operator()(
      const INPAR::InputParameterContainer& already_read_line)
  {
    int length = *already_read_line.Get<int>(definition_name_);
    return length;
  }
}  // namespace INPUT

FOUR_C_NAMESPACE_CLOSE
