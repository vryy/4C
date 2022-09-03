/*----------------------------------------------------------------------*/
/*! \file

\brief Parser for mathematical expressions, which contain literals
       ('1.0', 'pi', etc) and operations ('+', '-', 'sin', etc.)
       is a templated class. Thus its methods are defined in
       drt_parser.H (otherwise binding issues).
       A few non-templated methods of the Lexer base class
       are declared here.

\level 0
*/
/*---------------------------------------------------------------------*/


#include "drt_parser.H"
#include <algorithm>

namespace
{
  const std::array<std::string, 7> valid_operators = {"+", "-", "*", "/", "^", ".", ","};

  const std::array<std::string, 18> valid_functions = {"acos", "asin", "atan", "cos", "sin", "tan",
      "cosh", "sinh", "tanh", "exp", "log", "log10", "sqrt", "ceil", "heaviside", "fabs", "floor",
      "atan2"};

  const std::array<std::string, 2> valid_brackets = {"(", ")"};

  const std::array<std::string, 6> reserved_words = {"pi", "e", "E", "\t", "\n", " "};

  template <std::size_t n>
  [[nodiscard]] bool Contains(const std::array<std::string, n>& array, const std::string& element)
  {
    return std::find(array.begin(), array.end(), element) != array.end();
  }

}  // namespace

/*======================================================================*/
/* Lexer methods */

/*----------------------------------------------------------------------*/
/*!
\brief method used to step through std::string funct_
       delivers its character at position pos_++
\author u.kue
\date 10/07
*/
int DRT::PARSER::Lexer::GetNext()
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
\author u.kue
\date 10/07
*/
void DRT::PARSER::Lexer::Lexan()
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
      tok_ = tok_done;
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
          tok_ = tok_int;
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
        real_ = strtod(str_, NULL);
        tok_ = tok_real;
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
        tok_ = tok_name;
        integer_ = &(funct_[pos_]) - str_;  // length of operator name, e.g. 'sin' has '3'
        return;
      }
      else if (t == '+')
      {
        tok_ = tok_add;
        return;
      }
      else if (t == '-')
      {
        tok_ = tok_sub;
        return;
      }
      else if (t == '*')
      {
        tok_ = tok_mul;
        return;
      }
      else if (t == '/')
      {
        tok_ = tok_div;
        return;
      }
      else if (t == '^')
      {
        tok_ = tok_pow;
        return;
      }
      else if (t == '(')
      {
        tok_ = tok_lpar;
        return;
      }
      else if (t == ')')
      {
        tok_ = tok_rpar;
        return;
      }
      else if (t == ',')
      {
        tok_ = tok_comma;
        return;
      }
      else
      {
        if (t >= 32)
          dserror("unexpected char '%c' at pos %d", t, pos_);
        else
          dserror("unexpected char '%d' at pos %d", t, pos_);
        tok_ = tok_none;
        return;
      }
    }
  }
}

/*----------------------------------------------------------------------*/
/*!
\brief check if given string is an operator
\author vuong
\date 08/16
*/
bool DRT::PARSER::Lexer::IsOperator(const std::string& s) const
{
  return Contains(valid_operators, s);
}

/*----------------------------------------------------------------------*/
/*!
\brief check if given string is a function name
\author vuong
\date 08/16
*/
bool DRT::PARSER::Lexer::IsFunction(const std::string& s) const
{
  return Contains(valid_functions, s);
}

/*----------------------------------------------------------------------*/
/*!
\brief check if given string is a bracket
\author vuong
\date 08/16
*/
bool DRT::PARSER::Lexer::IsBracket(const std::string& s) const
{
  return Contains(valid_brackets, s);
}

/*----------------------------------------------------------------------*/
/*!
\brief check if given string is a reserved word
\author vuong
\date 08/16
*/
bool DRT::PARSER::Lexer::IsReservedWord(const std::string& s) const
{
  return Contains(reserved_words, s);
}
