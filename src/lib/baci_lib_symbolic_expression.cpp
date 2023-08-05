/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluating of arbitrary symbolic expressions, which contain literals
       ('1.0', 'pi', etc) operations ('+', '-','sin', etc.).
       The parsed expression is organised in a syntax tree whose nodes can
       either hold an operation or a literal.

\level 0

*/
/*---------------------------------------------------------------------*/


#include "baci_lib_symbolic_expression.H"

#include "baci_utils_exceptions.H"

#include <Sacado.hpp>

#include <cmath>
#include <cstring>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>

namespace
{
  const std::array<std::string, 7> valid_operators = {"+", "-", "*", "/", "^", ".", ","};

  const std::array<std::string, 16> valid_functions = {"acos", "asin", "atan", "cos", "sin", "tan",
      "cosh", "sinh", "tanh", "exp", "log", "log10", "sqrt", "heaviside", "fabs", "atan2"};

  const std::array<std::string, 2> valid_brackets = {"(", ")"};

  const std::array<std::string, 6> reserved_words = {"pi", "e", "E", "\t", "\n", " "};

  template <std::size_t n>
  [[nodiscard]] bool Contains(const std::array<std::string, n>& array, const std::string& element)
  {
    return std::find(array.begin(), array.end(), element) != array.end();
  }

}  // namespace



namespace DRT::UTILS::SYMBOLICEXPRESSIONDETAILS
{
  /*----------------------------------------------------------------------*/
  /*!
  \brief Syntax tree node holding binary, unary operator or literals
  */
  template <class T>
  class SyntaxTreeNode
  {
   public:
    using NodePtr = std::unique_ptr<SyntaxTreeNode<T>>;

    //! destructor
    ~SyntaxTreeNode() = default;

    //! copy constructor
    SyntaxTreeNode(const SyntaxTreeNode& other)
        : type_{other.type_},
          variable_{other.variable_},
          function_{other.function_},
          lhs_{std::make_unique<SyntaxTreeNode>(*other.lhs_)},
          rhs_{std::make_unique<SyntaxTreeNode>(*other.rhs_)} {};

    //! copy assignment operator
    SyntaxTreeNode& operator=(const SyntaxTreeNode& other)
    {
      if (&other == this) return *this;
      type_ = other.type;
      variable_ = other.variable_;
      function_ = other.function_;
      lhs_ = std::make_unique<SyntaxTreeNode>(*other.lhs_);
      rhs_ = std::make_unique<SyntaxTreeNode>(*other.rhs_);
      return *this;
    };

    //! move constructor
    SyntaxTreeNode(SyntaxTreeNode&& other) noexcept = default;

    //! move assignment operator
    SyntaxTreeNode& operator=(SyntaxTreeNode&& other) noexcept = default;

    enum NodeType
    {
      lt_variable,  // independent variable: 't', 'x', ...
      lt_number,    // number (literal)
      lt_function,
      lt_operator
    };

    //! Construct a SyntaxTreeNode with content @p type and a @p lhs and  @p rhs operand.
    SyntaxTreeNode(NodeType type, NodePtr lhs, NodePtr rhs)
        : type_(type), lhs_(std::move(lhs)), rhs_(std::move(rhs))
    {
    }

    //! type of the node, ie operator, literal, ...
    NodeType type_;

    //! particularise the node content: a literal number is double, an operatir has a character
    union
    {
      double number;  // holds 1.0, 7.e-9, etc
      char op;        // hold '+', '*', etc
    } v_;

    //! input string expressing the variable, holds 't', 'x', etc
    std::string variable_;
    //! input string expressing the function
    std::string function_;

    NodePtr lhs_;  // left hand side node
    NodePtr rhs_;  // right hand side node
  };


  /*----------------------------------------------------------------------*/
  /*!
  \brief Class holds auxiliar variables for Lexan method which steps through
         the string destilling the function tokens
  */
  class Lexer
  {
   public:
    //! constructor
    Lexer(std::string funct) : funct_(std::move(funct)), pos_(0) {}

    //! delivers funct_ character at position pos_++
    int GetNext();

    //! identifies a token (value and kind) in the funct_ string
    void Lexan();

    //! type of identifiable tokens in string funct_
    enum TokenType
    {
      tok_none,
      tok_done,
      tok_name,  // operator name, e.g. 'sin'
      tok_int,   // integer number
      tok_real,  // reals number
      tok_add,   // addition '+'
      tok_sub,   // subtraction and negation '-'
      tok_mul,   // multiplication '*'
      tok_div,   // division '/'
      tok_pow,   // power '^'
      tok_lpar,  // left parenthesis '('
      tok_rpar,  // right parenthesis ')'
      tok_comma  // comma ',' (used to separate function arguments)
    };

    std::string funct_;  // function description as string, i.e. "t^2", "sin(t)", etc
    unsigned pos_;       // current position in string funct_
    TokenType tok_;      // current token of string funct_
    char* str_;          // pointer to current character in funct_
    int integer_;        // translated integer number or length of operator word
    double real_;        // translated real number
  };

  //! helper structs
  template <typename T>
  struct IsFAD : public std::false_type
  {
  };

  template <typename T>
  struct IsFAD<Sacado::Fad::DFad<T>> : public std::true_type
  {
  };


  /*----------------------------------------------------------------------*/
  /*!
  \brief Parser
  */
  template <class T>
  class Parser
  {
   public:
    using NodePtr = typename SyntaxTreeNode<T>::NodePtr;

    //! constructor
    Parser(std::string funct);

    //! destructor
    ~Parser() = default;

    //! copy constuctor
    Parser(const Parser& other)
        : symbolicexpression_{other.symbolicexpression_},
          expr_{std::make_unique<SyntaxTreeNode<T>>(*other.expr_)},
          parsed_variable_constant_names_{other.parsed_variable_constant_names_} {};

    //! copy assignment operator
    Parser& operator=(const Parser& other)
    {
      if (&other == this) return *this;
      symbolicexpression_ = other.symbolicexpression_;
      expr_ = std::make_unique<NodePtr>(*other.expr_);
      parsed_variable_constant_names_ = other.parsed_variable_constant_names_;
      return *this;
    };

    //! move constructor
    Parser(Parser&& other) noexcept = default;

    //! move assignment operator
    Parser& operator=(Parser&& other) noexcept = default;

    /*!
     * @brief evaluates the parsed expression for a given set of variables
     *
     * @param[in] variable_values A map containing all variables (variablename, value) necessary to
     * evaluate the parsed expression
     * @return Value of the parsed expression
     */
    T EvaluateExpression(const std::map<std::string, T>& variable_values) const;

    /*!
     * @brief evaluates the derivative of the parsed expression with respect to a given set of
     * variables
     *
     * @param[in] variable_values A map containing all variables (variablename, value) necessary to
     * evaluate the parsed expression. Since the derivative of the parsed expression is evaluated,
     * only  Sacado::Fad::DFad<T> types are allowed
     * @param[in] constants A map containing all constants (constantname, value) necessary
     * to evaluate the parsed expression
     * @return  Derivative of the parsed expression with respect to the variables
     */
    template <typename T2, std::enable_if_t<IsFAD<T2>::value, void*> = nullptr>
    T EvaluateDerivative(const std::map<std::string, T2>& variable_values,
        const std::map<std::string, double>& constants = {}) const;

    //! Check if a variable with name 'varname' exists
    [[nodiscard]] bool IsVariable(const std::string& varname) const;

   private:
    NodePtr ParsePrimary(Lexer& lexer);
    NodePtr ParsePow(Lexer& lexer);
    NodePtr ParseTerm(Lexer& lexer);
    NodePtr ParseExpr(Lexer& lexer);
    NodePtr Parse(Lexer& lexer);

    //! given symbolic expression
    std::string symbolicexpression_;

    //! evaluates the parsed expression
    T Evaluate(const std::map<std::string, T>& variable_values,
        const std::map<std::string, double>& constants = {}) const;

    //! recursively extract corresponding number out of a syntax tree node
    T Interpret(const SyntaxTreeNode<T>& node) const;

    //! syntax tree root
    NodePtr expr_;

    //! set of all parsed variables
    std::set<std::string> parsed_variable_constant_names_;

    //! necessary variable values for evaluating the parsed expression
    mutable const std::map<std::string, T>* variable_values_{nullptr};

    //! necessary constant values for evaluating the parsed expression
    mutable const std::map<std::string, double>* constant_values_{nullptr};
  };


  /*======================================================================*/
  /* Lexer methods */

  /*----------------------------------------------------------------------*/
  /*!
  \brief method used to step through std::string funct_
         delivers its character at position pos_++
  */
  int Lexer::GetNext()
  {
    if (pos_ < funct_.length())
    {
      return funct_[pos_++];
    }
    else
    {
      return EOF;
    }
  }

  /*----------------------------------------------------------------------*/
  /*!
  \brief Identify current token
         type: tok_,
         value: integer_, real_,
         operator name: str_
  */
  void Lexer::Lexan()
  {
    for (;;)
    {
      int t = GetNext();
      if ((t == ' ') || (t == '\t'))
      {
        /* ignore whitespaces */
        /* this should never happen because we cannot read strings with
         * whitespaces from .dat files. :( */
      }
      else if (t == '\n')
      {
        dserror("newline in function definition");
      }
      else if (t == EOF)
      {
        tok_ = Lexer::Lexer::tok_done;
        return;
      }
      else
      {
        if (isdigit(t))
        {
          str_ = &(funct_[pos_ - 1]);
          while (isdigit(t))
          {
            t = GetNext();
          }
          if ((t != '.') && (t != 'E') && (t != 'e'))
          {
            if (t != EOF)
            {
              pos_--;
            }
            integer_ = atoi(str_);
            tok_ = Lexer::tok_int;
            return;
          }
          if (t == '.')
          {
            t = GetNext();
            if (isdigit(t))
            {
              while (isdigit(t))
              {
                t = GetNext();
              }
            }
            else
            {
              dserror("no digits after point at pos %d", pos_);
            }
          }
          if ((t == 'E') || (t == 'e'))
          {
            t = GetNext();
            if ((t == '-') || (t == '+'))
            {
              t = GetNext();
            }
            if (isdigit(t))
            {
              while (isdigit(t))
              {
                t = GetNext();
              }
            }
            else
            {
              dserror("no digits after exponent at pos %d", pos_);
            }
          }
          if (t != EOF)
          {
            pos_--;
          }
          real_ = strtod(str_, nullptr);
          tok_ = Lexer::tok_real;
          return;
        }
        else if (isalpha(t) || (t == '_'))
        {
          str_ = &(funct_[pos_ - 1]);
          while (isalnum(t) || (t == '_'))
          {
            t = GetNext();
          }
          if (t != EOF)
          {
            pos_--;
          }
          tok_ = Lexer::tok_name;
          integer_ = &(funct_[pos_]) - str_;  // length of operator name, e.g. 'sin' has '3'
          return;
        }
        else if (t == '+')
        {
          tok_ = Lexer::tok_add;
          return;
        }
        else if (t == '-')
        {
          tok_ = Lexer::tok_sub;
          return;
        }
        else if (t == '*')
        {
          tok_ = Lexer::tok_mul;
          return;
        }
        else if (t == '/')
        {
          tok_ = Lexer::tok_div;
          return;
        }
        else if (t == '^')
        {
          tok_ = Lexer::tok_pow;
          return;
        }
        else if (t == '(')
        {
          tok_ = Lexer::tok_lpar;
          return;
        }
        else if (t == ')')
        {
          tok_ = Lexer::tok_rpar;
          return;
        }
        else if (t == ',')
        {
          tok_ = Lexer::tok_comma;
          return;
        }
        else
        {
          if (t >= 32)
            dserror("unexpected char '%c' at pos %d", t, pos_);
          else
            dserror("unexpected char '%d' at pos %d", t, pos_);
          tok_ = Lexer::tok_none;
          return;
        }
      }
    }
  }

  /*===================================== da=================================*/
  /* Parser methods */

  /*----------------------------------------------------------------------*/
  /*!
  \brief Constructor of parser object
  */
  template <class T>
  Parser<T>::Parser(std::string funct)
  {
    //! set symbolic expression
    symbolicexpression_ = funct;

    //! create Lexer which stores all token of funct
    Lexer lexer{std::move(funct)};

    //! retrieve first token of funct
    lexer.Lexan();

    //! create syntax tree equivalent to funct
    expr_ = Parse(lexer);
  }

  /*----------------------------------------------------------------------*/
  /*!
  \brief Check if input is a variable name
  */
  template <class T>
  bool Parser<T>::IsVariable(const std::string& varname) const
  {
    //! check if variable exists
    return (parsed_variable_constant_names_.count(varname) != 0);
  }


  template <class T>
  T Parser<T>::EvaluateExpression(const std::map<std::string, T>& variable_values) const
  {
    return Evaluate(variable_values);
  }


  template <class T>
  template <typename T2, std::enable_if_t<IsFAD<T2>::value, void*>>
  T Parser<T>::EvaluateDerivative(const std::map<std::string, T2>& variable_values,
      const std::map<std::string, double>& constants) const
  {
    return Evaluate(variable_values, constants);
  }


  template <class T>
  T Parser<T>::Evaluate(const std::map<std::string, T>& variable_values,
      const std::map<std::string, double>& constants) const
  {
#ifdef DEBUG
    const bool all_required_variables_passed =
        std::all_of(parsed_variable_constant_names_.begin(), parsed_variable_constant_names_.end(),
            [&](const auto& var_name)
            { return (variable_values.count(var_name) + constants.count(var_name)) == 1; });

    if (!all_required_variables_passed)
    {
      std::string evaluate_variable_names =
          std::accumulate(variable_values.begin(), variable_values.end(), std::string(),
              [](const std::string& acc, const auto& v)
              { return acc.empty() ? v.first : acc + ", " + v.first; });

      std::string evaluate_constant_names =
          std::accumulate(constants.begin(), constants.end(), std::string(),
              [](const std::string& acc, const auto& v)
              { return acc.empty() ? v.first : acc + ", " + v.first; });

      dserror(
          "Some variables that this parser encountered in the expression are not passed to "
          "the Evaluate function.\n\n"
          "Expression:  %s \n"
          "Variables passed to Evaluate: %s \n"
          "Constants passed to Evaluate: %s",
          symbolicexpression_.c_str(), evaluate_variable_names.c_str(),
          evaluate_constant_names.c_str());
    }
#endif

    //! safety check if variable_values_ and constant_values_ are nullptr
    dsassert(variable_values_ == nullptr, "Internal error");
    dsassert(constant_values_ == nullptr, "Internal error");

    //! set variable values
    variable_values_ = &variable_values;
    //! set constant values if map of constants is not empty
    if (constants.size() != 0) constant_values_ = &constants;

    //! check if function has been parsed
    dsassert(expr_ != nullptr, "Internal error");

    //! evaluate syntax tree of function depending on set variables
    auto result = this->Interpret(*expr_);

    //! set variable_values_ and constant_values_ as nullptr again for safety reasons
    variable_values_ = nullptr;
    constant_values_ = nullptr;

    return result;
  }

  /*----------------------------------------------------------------------*/
  /*!
  \brief Parse primary entities, i.e. literals and unary operators,
         such as numbers, parentheses, independent variables, operator names
  */
  template <class T>
  auto Parser<T>::ParsePrimary(Lexer& lexer) -> NodePtr
  {
    NodePtr lhs = nullptr;

    switch (lexer.tok_)
    {
      case Lexer::tok_lpar:
        lexer.Lexan();
        lhs = ParseExpr(lexer);
        if (lexer.tok_ != Lexer::tok_rpar) dserror("')' expected");
        lexer.Lexan();
        break;
      case Lexer::tok_int:
        lhs = std::make_unique<SyntaxTreeNode<T>>(SyntaxTreeNode<T>::lt_number, nullptr, nullptr);
        lhs->v_.number = lexer.integer_;
        lexer.Lexan();
        break;
      case Lexer::tok_real:
        lhs = std::make_unique<SyntaxTreeNode<T>>(SyntaxTreeNode<T>::lt_number, nullptr, nullptr);
        lhs->v_.number = lexer.real_;
        lexer.Lexan();
        break;
      case Lexer::tok_sub:
      {
        NodePtr rhs;
        lexer.Lexan();
        /*rhs = parse_primary();*/
        rhs = ParsePow(lexer);
        if (rhs->type_ == SyntaxTreeNode<T>::lt_number)
        {
          rhs->v_.number *= -1;
          lhs = std::move(rhs);
        }
        else
        {
          lhs = std::make_unique<SyntaxTreeNode<T>>(SyntaxTreeNode<T>::lt_number, nullptr, nullptr);
          lhs->v_.number = -1;
          lhs = std::make_unique<SyntaxTreeNode<T>>(
              SyntaxTreeNode<T>::lt_operator, std::move(lhs), std::move(rhs));
          lhs->v_.op = '*';
        }
        break;
      }
      case Lexer::tok_name:
      {
        // get substring starting from str_ with length of lexer.integer_
        std::string name(lexer.str_, lexer.integer_);
        if ((lexer.integer_ == 2) && (std::strncmp("pi", lexer.str_, lexer.integer_) == 0))
        {
          lhs = std::make_unique<SyntaxTreeNode<T>>(SyntaxTreeNode<T>::lt_number, nullptr, nullptr);
          lhs->v_.number = M_PI;
          lexer.Lexan();
          break;
        }
        else
        {
          if (name == "acos" or name == "asin" or name == "atan" or name == "cos" or
              name == "sin" or name == "tan" or name == "cosh" or name == "sinh" or
              name == "tanh" or name == "exp" or name == "log" or name == "log10" or
              name == "sqrt" or name == "ceil" or name == "heaviside" or name == "fabs" or
              name == "floor")
          {
            lhs = std::make_unique<SyntaxTreeNode<T>>(
                SyntaxTreeNode<T>::lt_function, nullptr, nullptr);
            lhs->function_ = name;
            lexer.Lexan();
            if (lexer.tok_ != Lexer::tok_lpar)
              dserror("'(' expected after function name '%s'", name.c_str());
            lexer.Lexan();
            lhs->lhs_ = ParseExpr(lexer);
            if (lexer.tok_ != Lexer::tok_rpar) dserror("')' expected");
            lexer.Lexan();
            break;
          }
          else if (name == "atan2")
          {
            lhs = std::make_unique<SyntaxTreeNode<T>>(
                SyntaxTreeNode<T>::lt_function, nullptr, nullptr);
            lhs->function_ = name;
            lexer.Lexan();
            if (lexer.tok_ != Lexer::tok_lpar)
              dserror("'(' expected after function name '%s'", name.c_str());
            lexer.Lexan();
            lhs->lhs_ = ParseExpr(lexer);
            if (lexer.tok_ != Lexer::tok_comma) dserror("',' expected");
            lexer.Lexan();
            lhs->rhs_ = ParseExpr(lexer);
            if (lexer.tok_ != Lexer::tok_rpar)
              dserror("')' expected for function name '%s'", name.c_str());
            lexer.Lexan();
            break;
          }
          else
          {
            lhs = std::make_unique<SyntaxTreeNode<T>>(
                SyntaxTreeNode<T>::lt_variable, nullptr, nullptr);
            lhs->variable_ = name;
            parsed_variable_constant_names_.insert(name);
            lexer.Lexan();
          }
        }
        break;
      }
      default:
        dserror("unexpected token %d", lexer.tok_);
        break;
    }

    return lhs;
  }


  /*----------------------------------------------------------------------*/
  /*!
  \brief Parse entities connected by power: a^b
  */
  template <class T>
  auto Parser<T>::ParsePow(Lexer& lexer) -> NodePtr
  {
    NodePtr lhs;
    NodePtr rhs;

    lhs = ParsePrimary(lexer);
    for (;;)
    {
      if (lexer.tok_ == Lexer::tok_pow)
      {
        lexer.Lexan();
        rhs = ParsePrimary(lexer);
        lhs = std::make_unique<SyntaxTreeNode<T>>(
            SyntaxTreeNode<T>::lt_operator, std::move(lhs), std::move(rhs));
        lhs->v_.op = '^';
      }
      else
      {
        break;
      }
    }

    return lhs;
  }


  /*----------------------------------------------------------------------*/
  /*!
  \brief Parse entities connected by multiplication or division: a*b, a/b
  */
  template <class T>
  auto Parser<T>::ParseTerm(Lexer& lexer) -> NodePtr
  {
    NodePtr lhs;
    NodePtr rhs;

    lhs = ParsePow(lexer);
    for (;;)
    {
      if (lexer.tok_ == Lexer::tok_mul)
      {
        lexer.Lexan();
        rhs = ParsePow(lexer);
        lhs = std::make_unique<SyntaxTreeNode<T>>(
            SyntaxTreeNode<T>::lt_operator, std::move(lhs), std::move(rhs));
        lhs->v_.op = '*';
      }
      else if (lexer.tok_ == Lexer::tok_div)
      {
        lexer.Lexan();
        rhs = ParsePow(lexer);
        lhs = std::make_unique<SyntaxTreeNode<T>>(
            SyntaxTreeNode<T>::lt_operator, std::move(lhs), std::move(rhs));
        lhs->v_.op = '/';
      }
      else
      {
        break;
      }
    }

    return lhs;
  }

  /*----------------------------------------------------------------------*/
  /*!
  \brief Parse entity
  */
  template <class T>
  auto Parser<T>::Parse(Lexer& lexer) -> NodePtr
  {
    NodePtr lhs;

    lhs = ParseExpr(lexer);

    // check for invalid tokens at the beginning of a parse entities
    if (lexer.tok_ == Lexer::tok_comma or lexer.tok_ == Lexer::tok_rpar)
    {
      dserror(
          "unexpected token %d. Invalid syntax. Missing brackets or comma instead of decimal "
          "point?",
          lexer.tok_);
    }

    return lhs;
  }

  /*----------------------------------------------------------------------*/
  /*!
  \brief Parse entities connected by addition or subtraction: a+b, a-b
  */
  template <class T>
  auto Parser<T>::ParseExpr(Lexer& lexer) -> NodePtr
  {
    typename SyntaxTreeNode<T>::NodePtr lhs;
    typename SyntaxTreeNode<T>::NodePtr rhs;

    lhs = ParseTerm(lexer);
    for (;;)
    {
      if (lexer.tok_ == Lexer::tok_add)
      {
        lexer.Lexan();
        rhs = ParseTerm(lexer);
        lhs = std::make_unique<SyntaxTreeNode<T>>(
            SyntaxTreeNode<T>::lt_operator, std::move(lhs), std::move(rhs));
        lhs->v_.op = '+';
      }
      else if (lexer.tok_ == Lexer::tok_sub)
      {
        lexer.Lexan();
        rhs = ParseTerm(lexer);
        lhs = std::make_unique<SyntaxTreeNode<T>>(
            SyntaxTreeNode<T>::lt_operator, std::move(lhs), std::move(rhs));
        lhs->v_.op = '-';
      }
      else
      {
        break;
      }
    }

    return lhs;
  }


  /*----------------------------------------------------------------------*/
  /*!
  \brief Recursively extract corresponding number out of a syntax tree node
  */
  template <class T>
  T Parser<T>::Interpret(const SyntaxTreeNode<T>& node) const
  {
    T res = 0;  // the result

    switch (node.type_)
    {
      // literal numbers: leaf of syntax tree node
      case SyntaxTreeNode<T>::lt_number:
        res = node.v_.number;
        break;
      // binary operators: bifurcating branch of syntax tree node
      case SyntaxTreeNode<T>::lt_operator:
      {
        T lhs;
        T rhs;

        // recursively visit branches and obtain sub-results
        lhs = Interpret(*node.lhs_);
        rhs = Interpret(*node.rhs_);

        // evaluate the node operator
        switch (node.v_.op)
        {
          case '+':
            res = lhs + rhs;
            break;
          case '-':
            res = lhs - rhs;
            break;
          case '*':
            res = lhs * rhs;
            break;
          case '/':
            /* check for rhs==0.0? */
            res = lhs / rhs;
            break;
          case '^':
            res = std::pow(lhs, rhs);
            break;
          default:
            dserror("unsupported operator '%c'", node.v_.op);
        }
        break;
      }
      // independent variables: as set by user
      case SyntaxTreeNode<T>::lt_variable:
      {
        if (parsed_variable_constant_names_.count(node.variable_) != 0)
        {
          if (constant_values_ == nullptr)
          {
            if (variable_values_->find(node.variable_) == variable_values_->end())
            {
              dserror("variable or constant '%s' not given as input in Evaluate()",
                  node.variable_.c_str());
            }
            else
            {
              res = variable_values_->at(std::string(node.variable_));
            }
          }
          else
          {
            if ((variable_values_->find(node.variable_) == variable_values_->end()) &&
                (constant_values_->find(node.variable_) == constant_values_->end()))
            {
              dserror("variable or constant '%s' not given as input in EvaluateDeriv()",
                  node.variable_.c_str());
            }
            else
            {
              if (variable_values_->find(node.variable_) != variable_values_->end())
              {
                res = variable_values_->at(node.variable_);
              }
              else if (constant_values_->find(node.variable_) != constant_values_->end())
              {
                res = constant_values_->at(node.variable_);
              }
              else
                dserror("Something went really wrong!");
            }
          }
        }
        else
          dserror("unknown variable '%s'", node.variable_.c_str());
        break;
      }
      // unary operators
      case SyntaxTreeNode<T>::lt_function:
      {
        T arg;
        arg = Interpret(*node.lhs_);
        if (node.function_ == "acos")
          res = acos(arg);
        else if (node.function_ == "asin")
          res = asin(arg);
        else if (node.function_ == "atan")
          res = atan(arg);
        else if (node.function_ == "cos")
          res = cos(arg);
        else if (node.function_ == "sin")
          res = sin(arg);
        else if (node.function_ == "tan")
          res = tan(arg);
        else if (node.function_ == "cosh")
          res = cosh(arg);
        else if (node.function_ == "sinh")
          res = sinh(arg);
        else if (node.function_ == "tanh")
          res = tanh(arg);
        else if (node.function_ == "exp")
          res = exp(arg);
        else if (node.function_ == "log")
          res = log(arg);
        else if (node.function_ == "log10")
          res = log10(arg);
        else if (node.function_ == "sqrt")
          res = sqrt(arg);
        else if (node.function_ == "atan2")
        {
          T arg2;
          // recursively visit branches and obtain sub-results
          arg2 = Interpret(*node.rhs_);
          res = atan2(arg, arg2);
        }
        else if (node.function_ == "fabs")
          res = fabs(arg);
        else if (node.function_ == "heaviside")
        {
          if (arg > 0)
          {
            res = 1.0;
          }
          else
          {
            res = 0.0;
          }
        }
        else
          dserror("unknown function_ '%s'", node.function_.c_str());
        break;
      }
      default:
        dserror("unknown syntax tree node type");
        break;
    }

    return res;
  }

}  // namespace DRT::UTILS::SYMBOLICEXPRESSIONDETAILS

template <typename T>
DRT::UTILS::SymbolicExpression<T>::SymbolicExpression(const std::string& expression)
    : parser_for_value_(
          std::make_unique<DRT::UTILS::SYMBOLICEXPRESSIONDETAILS::Parser<ValueType>>(expression)),
      parser_for_firstderivative_(
          std::make_unique<DRT::UTILS::SYMBOLICEXPRESSIONDETAILS::Parser<FirstDerivativeType>>(
              expression)),
      parser_for_secondderivative_(
          std::make_unique<DRT::UTILS::SYMBOLICEXPRESSIONDETAILS::Parser<SecondDerivativeType>>(
              expression))
{
}



template <typename T>
auto DRT::UTILS::SymbolicExpression<T>::Value(
    const std::map<std::string, ValueType>& variable_values) const -> ValueType
{
  return parser_for_value_->EvaluateExpression(variable_values);
}


template <typename T>
auto DRT::UTILS::SymbolicExpression<T>::FirstDerivative(
    std::map<std::string, FirstDerivativeType> variable_values,
    const std::map<std::string, ValueType>& constant_values) const -> FirstDerivativeType
{
  return parser_for_firstderivative_->EvaluateDerivative(variable_values, constant_values);
}


template <typename T>
auto DRT::UTILS::SymbolicExpression<T>::SecondDerivative(
    const std::map<std::string, SecondDerivativeType>& variable_values,
    const std::map<std::string, ValueType>& constant_values) const -> SecondDerivativeType
{
  return parser_for_secondderivative_->EvaluateDerivative(variable_values, constant_values);
}


template <typename Number>
DRT::UTILS::SymbolicExpression<Number>::SymbolicExpression(
    const DRT::UTILS::SymbolicExpression<Number>& other)
    : parser_for_value_{std::make_unique<DRT::UTILS::SYMBOLICEXPRESSIONDETAILS::Parser<ValueType>>(
          *other.parser_for_value_)},
      parser_for_firstderivative_{
          std::make_unique<DRT::UTILS::SYMBOLICEXPRESSIONDETAILS::Parser<FirstDerivativeType>>(
              *other.parser_for_firstderivative_)},
      parser_for_secondderivative_{
          std::make_unique<DRT::UTILS::SYMBOLICEXPRESSIONDETAILS::Parser<SecondDerivativeType>>(
              *other.parser_for_secondderivative_)}
{
}


template <typename Number>
DRT::UTILS::SymbolicExpression<Number>& DRT::UTILS::SymbolicExpression<Number>::operator=(
    const DRT::UTILS::SymbolicExpression<Number>& other)
{
  parser_for_value_ = std::make_unique<DRT::UTILS::SYMBOLICEXPRESSIONDETAILS::Parser<ValueType>>(
      *other.parser_for_value_);
  parser_for_firstderivative_ =
      std::make_unique<DRT::UTILS::SYMBOLICEXPRESSIONDETAILS::Parser<FirstDerivativeType>>(
          *other.parser_for_firstderivative_);
  parser_for_secondderivative_ =
      std::make_unique<DRT::UTILS::SYMBOLICEXPRESSIONDETAILS::Parser<SecondDerivativeType>>(
          *other.parser_for_secondderivative_);
  return *this;
}



template <typename Number>
DRT::UTILS::SymbolicExpression<Number>::~SymbolicExpression() = default;

// explicit instantiations
template class DRT::UTILS::SymbolicExpression<double>;
