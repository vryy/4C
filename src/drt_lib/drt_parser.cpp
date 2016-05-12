/*----------------------------------------------------------------------*/
/*!
\file drt_parser.cpp

\brief Implementation Parser for mathematical expressions, which contain literals
       ('1.0', 'pi', etc) and operations ('+', '-', 'sin', etc.)
       is a templated class. Thus its methods are defined in
       drt_parser.H (otherwise binding issues).
       A few non-templated methods of the Lexer base class
       are declared here.

<pre>
\brief Implementation
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

\author u.kue
\date 10/07
*/


#include "drt_parser.H"

/*======================================================================*/
/* Lexer methods */

/*----------------------------------------------------------------------*/
/*!
\brief Implementation method used to step through std::string funct_
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
\brief Implementation Identify current token
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
  str_ = &(funct_[pos_-1]);
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
  str_ = &(funct_[pos_-1]);
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



