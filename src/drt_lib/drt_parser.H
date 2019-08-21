/*----------------------------------------------------------------------*/
/*! \file

\brief Parser for mathematical expressions, which contain literals
       ('1.0', 'pi', etc) operations ('+', '-', 'sin', etc.).
       The parsed expression is organised in a syntax tree whose nodes
       can either hold an operation or a literal.

\level 0

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#ifndef DRT_PARSER_H
#define DRT_PARSER_H

#include <iostream>
#include <map>
#include <cmath>
#include <string>
#include <cstring>
#include <cstdlib>

#include "drt_dserror.H"

namespace DRT
{
  namespace PARSER
  {
    /*----------------------------------------------------------------------*/
    /*!
    \brief Syntax tree node holding binary, unary operator or literals
    \author u.kue
    \date 10/07
    */
    template <class T>
    class SyntaxTreeNode
    {
     public:
      enum NodeType
      {
        lt_none,
        lt_variable,  // independent variabel: 't', 'x', ...
        lt_number,    // number (literal)
        lt_function,
        lt_operator
      };

      //! constructor
      SyntaxTreeNode(NodeType type, SyntaxTreeNode* lhs = NULL, SyntaxTreeNode* rhs = NULL)
          : type_(type), lhs_(lhs), rhs_(rhs)
      {
      }

      //! destructor
      ~SyntaxTreeNode()
      {
        delete lhs_;
        delete rhs_;
      }

      //! type of the node, ie operator, literal, ...
      NodeType type_;

      //! particularise the node content: a literal number is double, an operatir has a character
      union {
        double number;  // holds 1.0, 7.e-9, etc
        char op;        // hold '+', '*', etc
      } v_;

      //! input string expressing the variable, holds 't', 'x', etc
      std::string variable_;
      //! input string expressing the function
      std::string function_;

      SyntaxTreeNode<T>* lhs_;  // left hand side node
      SyntaxTreeNode<T>* rhs_;  // right hand side node
    };


    /*----------------------------------------------------------------------*/
    /*!
    \brief Class holds auxiliar variables for Lexan method which steps through
           the string destilling the function tokens
    \author u.kue
    \date 10/07
    */
    class Lexer
    {
     protected:
      //! constructor
      Lexer(std::string funct) : funct_(funct), pos_(0) {}

      //! delivers funct_ character at position pos_++
      int GetNext();

      //! identifies a token (value and kind) in the funct_ string
      void Lexan();

      //! helper methods
      bool IsOperator(const std::string& s) const;
      bool IsFunction(const std::string& s) const;
      bool IsBracket(const std::string& s) const;
      bool IsReservedWord(const std::string& s) const;

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
        tok_mod,   // modulus ??? (NOT AVAILABLE)
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

      static const std::string operator_list_[];     // list of all valid operators
      static const std::size_t operator_list_size_;  // number of all valid operators

      static const std::string function_list_[];     // list of all valid functions
      static const std::size_t function_list_size_;  // number of all valid functions

      static const std::string bracket_list_[];     // list of all valid brackets
      static const std::size_t bracket_list_size_;  // number of all valid brackets

      static const std::string reserved_words_list_[];     // list of all special words
      static const std::size_t reserved_words_list_size_;  // number of all special words
    };

    /*----------------------------------------------------------------------*/
    /*!
    \brief Parser
    \author u.kue
    \date 10/07
    */
    template <class T>
    class Parser : public Lexer
    {
     public:
      //! constructor
      Parser(std::string funct);
      //! destructor
      ~Parser();

      //! parse the function
      void ParseFunction();

      //! evaluate function depending on set variables
      T Evaluate();

      //! Add variable to parser object
      void AddVariable(const std::string& varname, const T& value);
      //! Set value of an existing variable to parser object
      void SetValue(const std::string& varname, const T& value);
      //! Check if a variable with name 'varname' exists
      bool IsVariable(const std::string& varname) const;
      //! Check if variable name is valid
      bool IsValidVariableName(const std::string& varname) const;

     private:
      SyntaxTreeNode<T>* ParsePrimary();
      SyntaxTreeNode<T>* ParsePow();
      SyntaxTreeNode<T>* ParseTerm();
      SyntaxTreeNode<T>* ParseExpr();
      SyntaxTreeNode<T>* Parse();

      //! recursively extract corresponding number out of a syntax tree node
      T Interpret(SyntaxTreeNode<T>* node);

      //! syntax tree root
      SyntaxTreeNode<T>* expr_;

      //! map of all independent variables
      std::map<std::string, T> variable_list_;
    };


    /*======================================================================*/
    /* Parser methods */

    /*----------------------------------------------------------------------*/
    /*!
    \brief Constructor of parser object
    \author u.kue
    \date 10/07
    */
    template <class T>
    Parser<T>::Parser(std::string funct) : Lexer(funct), expr_(NULL)
    {
    }


    /*----------------------------------------------------------------------*/
    /*!
    \brief Destructor of parser object
    \author u.kue
    \date 10/07
    */
    template <class T>
    Parser<T>::~Parser()
    {
      delete expr_;
    }

    /*----------------------------------------------------------------------*/
    /*!
    \brief Parse the function, i.e. create the syntax tree
    \author vuong
    \date 08/16
    */
    template <class T>
    void Parser<T>::ParseFunction()
    {
      Lexan();          // retrieve first token of funct_
                        // it is stored in (base) Lexer variables of Parser object
      expr_ = Parse();  // create syntax tree equivalent to funct_
    }


    /*----------------------------------------------------------------------*/
    /*!
    \brief Add a variable definition to the parser object
    \author vuong
    \date 08/16
    */
    template <class T>
    void Parser<T>::AddVariable(const std::string& varname, const T& value)
    {
      // check if variable is valid
      if (IsValidVariableName(varname))
        // set default value
        variable_list_[varname] = value;
      else
        dserror("Variable name '%s' already defined or invalid!", varname.c_str());
    }

    /*----------------------------------------------------------------------*/
    /*!
    \brief Set value of a variable in the parser object
    \author vuong
    \date 08/16
    */
    template <class T>
    void Parser<T>::SetValue(const std::string& varname, const T& value)
    {
      // check if variable exists
      if (IsVariable(varname))
        // set value
        variable_list_[varname] = value;
      else
        dserror("Variable '%s' not defined or invalid!", varname.c_str());
    }

    /*----------------------------------------------------------------------*/
    /*!
    \brief Check if input is a variable name
    \author vuong
    \date 08/16
    */
    template <class T>
    bool Parser<T>::IsVariable(const std::string& varname) const
    {
      // check if variable exists
      return (variable_list_.find(varname) != variable_list_.end());
    }

    /*----------------------------------------------------------------------*/
    /*!
    \brief Check if input is a valid variable name
    \author vuong
    \date 08/16
    */
    template <class T>
    bool Parser<T>::IsValidVariableName(const std::string& varname) const
    {
      // check if variable is valid
      return (not varname.empty() and      // there must be a name
              not isdigit(varname[0]) and  // the first character must not be a digit
              not iscntrl(varname[0]) and  // the first character must not be a control character
              not IsVariable(varname) and  // it must not be a variable defined yet
              not IsOperator(varname) and  // it must not be a operator
              not IsFunction(varname) and  // it must not be a function name
              not IsBracket(varname) and   // it must not be a bracket
              not IsReservedWord(varname)  // it must not be a reserved 'word'
      );
    }

    /*----------------------------------------------------------------------*/
    /*!
    \brief global interpreter call.
    \author vuong
    \date 08/16
    */
    template <class T>
    T Parser<T>::Evaluate()
    {
      if (expr_ == NULL) dserror("Function has not been parsed! Call ParseFunction()!");
      // evaluate syntax tree of function depending on set variables
      return Interpret(expr_);
    }

    /*----------------------------------------------------------------------*/
    /*!
    \brief Parse primary entities, i.e. literals and unary operators,
           such as numbers, parentheses, independent variables, operator names
    \author u.kue
    \date 10/07
    */
    template <class T>
    SyntaxTreeNode<T>* Parser<T>::ParsePrimary()
    {
      SyntaxTreeNode<T>* lhs = NULL;

      switch (tok_)
      {
        case tok_lpar:
          Lexan();
          lhs = ParseExpr();
          if (tok_ != tok_rpar) dserror("')' expected");
          Lexan();
          break;
        case tok_int:
          lhs = new SyntaxTreeNode<T>(SyntaxTreeNode<T>::lt_number, NULL, NULL);
          lhs->v_.number = integer_;
          Lexan();
          break;
        case tok_real:
          lhs = new SyntaxTreeNode<T>(SyntaxTreeNode<T>::lt_number, NULL, NULL);
          lhs->v_.number = real_;
          Lexan();
          break;
        case tok_sub:
        {
          SyntaxTreeNode<T>* rhs;
          Lexan();
          /*rhs = parse_primary();*/
          rhs = ParsePow();
          if (rhs->type_ == SyntaxTreeNode<T>::lt_number)
          {
            rhs->v_.number *= -1;
            lhs = rhs;
          }
          else
          {
            lhs = new SyntaxTreeNode<T>(SyntaxTreeNode<T>::lt_number, NULL, NULL);
            lhs->v_.number = -1;
            lhs = new SyntaxTreeNode<T>(SyntaxTreeNode<T>::lt_operator, lhs, rhs);
            lhs->v_.op = '*';
          }
          break;
        }
        case tok_name:
        {
          // get substring starting from str_ with length of integer_
          std::string name(str_, integer_);
          if (IsVariable(name))
          {
            lhs = new SyntaxTreeNode<T>(SyntaxTreeNode<T>::lt_variable, NULL, NULL);
            lhs->variable_ = name;
            Lexan();
          }
          else
          {
            if ((integer_ == 2) && (std::strncmp("pi", str_, integer_) == 0))
            {
              lhs = new SyntaxTreeNode<T>(SyntaxTreeNode<T>::lt_number, NULL, NULL);
              lhs->v_.number = M_PI;
              Lexan();
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
                lhs = new SyntaxTreeNode<T>(SyntaxTreeNode<T>::lt_function, NULL, NULL);
                lhs->function_ = name;
                Lexan();
                if (tok_ != tok_lpar)
                  dserror("'(' expected after function name '%s'", name.c_str());
                Lexan();
                lhs->lhs_ = ParseExpr();
                if (tok_ != tok_rpar) dserror("')' expected");
                Lexan();
                break;
              }
              else if (name == "atan2")
              {
                lhs = new SyntaxTreeNode<T>(SyntaxTreeNode<T>::lt_function, NULL, NULL);
                lhs->function_ = name;
                Lexan();
                if (tok_ != tok_lpar)
                  dserror("'(' expected after function name '%s'", name.c_str());
                Lexan();
                lhs->lhs_ = ParseExpr();
                if (tok_ != tok_comma) dserror("',' expected");
                Lexan();
                lhs->rhs_ = ParseExpr();
                if (tok_ != tok_rpar) dserror("')' expected for function name '%s'", name.c_str());
                Lexan();
                break;
              }
              else
              {
                str_[integer_] = '\0';
                dserror("unknown function '%s'", str_);
              }
            }
          }
          break;
        }
        default:
          dserror("unexpected token %d", tok_);
          break;
      }

      return lhs;
    }


    /*----------------------------------------------------------------------*/
    /*!
    \brief Parse entities connected by power: a^b
    \author u.kue
    \date 10/07
    */
    template <class T>
    SyntaxTreeNode<T>* Parser<T>::ParsePow()
    {
      SyntaxTreeNode<T>* lhs;
      SyntaxTreeNode<T>* rhs;

      lhs = ParsePrimary();
      for (;;)
      {
        if (tok_ == tok_pow)
        {
          Lexan();
          rhs = ParsePrimary();
          lhs = new SyntaxTreeNode<T>(SyntaxTreeNode<T>::lt_operator, lhs, rhs);
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
    \author u.kue
    \date 10/07
    */
    template <class T>
    SyntaxTreeNode<T>* Parser<T>::ParseTerm()
    {
      SyntaxTreeNode<T>* lhs;
      SyntaxTreeNode<T>* rhs;

      lhs = ParsePow();
      for (;;)
      {
        if (tok_ == tok_mul)
        {
          Lexan();
          rhs = ParsePow();
          lhs = new SyntaxTreeNode<T>(SyntaxTreeNode<T>::lt_operator, lhs, rhs);
          lhs->v_.op = '*';
        }
        else if (tok_ == tok_div)
        {
          Lexan();
          rhs = ParsePow();
          lhs = new SyntaxTreeNode<T>(SyntaxTreeNode<T>::lt_operator, lhs, rhs);
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
    \author u.kue
    \date 10/07
    */
    template <class T>
    SyntaxTreeNode<T>* Parser<T>::Parse()
    {
      SyntaxTreeNode<T>* lhs;

      lhs = ParseExpr();

      // check for invalid tokens at the beginning of a parse entities
      if (tok_ == tok_comma or tok_ == tok_rpar)
        dserror(
            "unexpected token %d. Invalid syntax. Missing brackets or comma instead of decimal "
            "point?",
            tok_);

      return lhs;
    }

    /*----------------------------------------------------------------------*/
    /*!
    \brief Parse entities connected by addition or subtraction: a+b, a-b
    \author u.kue
    \date 10/07
    */
    template <class T>
    SyntaxTreeNode<T>* Parser<T>::ParseExpr()
    {
      SyntaxTreeNode<T>* lhs;
      SyntaxTreeNode<T>* rhs;

      lhs = ParseTerm();
      for (;;)
      {
        if (tok_ == tok_add)
        {
          Lexan();
          rhs = ParseTerm();
          lhs = new SyntaxTreeNode<T>(SyntaxTreeNode<T>::lt_operator, lhs, rhs);
          lhs->v_.op = '+';
        }
        else if (tok_ == tok_sub)
        {
          Lexan();
          rhs = ParseTerm();
          lhs = new SyntaxTreeNode<T>(SyntaxTreeNode<T>::lt_operator, lhs, rhs);
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
    \author u.kue
    \date 10/07
    */
    template <class T>
    T Parser<T>::Interpret(SyntaxTreeNode<T>* node)
    {
      T res = 0;  // the result

      switch (node->type_)
      {
        // literal numbers: leaf of syntax tree node
        case SyntaxTreeNode<T>::lt_number:
          res = node->v_.number;
          break;
        // binary operators: bifurcating branch of syntax tree node
        case SyntaxTreeNode<T>::lt_operator:
        {
          T lhs;
          T rhs;

          // recursively visit branches and obtain sub-results
          lhs = Interpret(node->lhs_);
          rhs = Interpret(node->rhs_);

          // evaluate the node operator
          switch (node->v_.op)
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
              dserror("unsupported operator '%c'", node->v_.op);
          }
          break;
        }
        // independent variables: as set by user
        case SyntaxTreeNode<T>::lt_variable:
        {
          typename std::map<std::string, T>::iterator iter = variable_list_.find(node->variable_);
          if (iter != variable_list_.end())
            res = iter->second;
          else
            dserror("unknown variable '%s'", node->variable_.c_str());
          break;
        }
        // unary operators
        case SyntaxTreeNode<T>::lt_function:
        {
          T arg;
          arg = Interpret(node->lhs_);
          if (node->function_ == "acos")
            res = acos(arg);
          else if (node->function_ == "asin")
            res = asin(arg);
          else if (node->function_ == "atan")
            res = atan(arg);
          else if (node->function_ == "cos")
            res = cos(arg);
          else if (node->function_ == "sin")
            res = sin(arg);
          else if (node->function_ == "tan")
            res = tan(arg);
          else if (node->function_ == "cosh")
            res = cosh(arg);
          else if (node->function_ == "sinh")
            res = sinh(arg);
          else if (node->function_ == "tanh")
            res = tanh(arg);
          else if (node->function_ == "exp")
            res = exp(arg);
          else if (node->function_ == "log")
            res = log(arg);
          else if (node->function_ == "log10")
            res = log10(arg);
          else if (node->function_ == "sqrt")
            res = sqrt(arg);
          else if (node->function_ == "atan2")
          {
            T arg2;
            // recursively visit branches and obtain sub-results
            arg2 = Interpret(node->rhs_);
            res = atan2(arg, arg2);
          }
          // else if (node->function_=="ceil") res = ceil(arg);
          else if (node->function_ == "fabs")
            res = fabs(arg);
          // else if (node->function_=="floor") res = floor(arg);
          else if (node->function_ == "heaviside")
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
            dserror("unknown function_ '%s'", node->function_.c_str());
          break;
        }
        default:
          dserror("unknown syntax tree node type");
          break;
      }

      return res;
    }

  }  // namespace PARSER
}  // end namespace DRT

#endif  // DRT_PARSER_H
