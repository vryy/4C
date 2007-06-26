/*!
\file
\brief A top down parser for function evaluation.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

This is a parser for simple expressions only. Only numbers,
parenthesis, ordinary operators ( + - * / ^ ) and some unary functions
are allowed. Additionally the variables x, y and z as well as t are
substituted for their respective values.

The expression string is parsed by pss_parse and an internal syntax
tree is returned. This can be used by pss_evaluate to evaluate the
expression for particular values.

There is a lot of internal stuff here that is only used within this
file.
*/

#include "pss_parser.h"

/*----------------------------------------------------------------------*/
/*!
  \brief The mathematical functions supported by the parser.

  Right now only one argument functions are supported.

  \author u.kue
  \date 12/05
*/
/*----------------------------------------------------------------------*/
typedef DOUBLE (*mathfunc)(DOUBLE);
typedef struct _FUNCTAB
{
  char* name;
  mathfunc func;
} FUNCTAB;

FUNCTAB function[] =
{
  { "acos", acos },
  { "asin", asin },
  { "atan", atan },
  /*atan2*/

  { "cos", cos },
  { "sin", sin },
  { "tan", tan },

  { "cosh", cosh },
  { "sinh", sinh },
  { "tanh", tanh },

  { "exp", exp },
  { "log", log },
  { "log10", log10 },

  { "sqrt", sqrt },

  { "ceil", ceil },
  { "fabs", fabs },
  { "floor", floor },

  { NULL, NULL }
};


/*----------------------------------------------------------------------*/
/*!
  \brief The types of nodes in the syntax tree.

  \author u.kue
  \date 12/05
*/
/*----------------------------------------------------------------------*/
typedef enum _ST_NODE_TYPE
{
  lt_none,
  lt_variable,
  lt_number,
  lt_function,
  lt_operator
} ST_NODE_TYPE;


/*----------------------------------------------------------------------*/
/*!
  \brief The syntax tree structure.

  The whole expression is represented by a syntax tree.

  \author u.kue
  \date 12/05
*/
/*----------------------------------------------------------------------*/
typedef struct _ST_NODE
{
  ST_NODE_TYPE type;
  union
  {
    DOUBLE number;
    CHAR variable;
    mathfunc function;
    CHAR operator;
  } v;

  struct _ST_NODE* lhs;
  struct _ST_NODE* rhs;
} ST_NODE;


static ST_NODE* st_node_new(ST_NODE_TYPE type, ST_NODE* lhs, ST_NODE* rhs)
{
  ST_NODE* node;

#ifdef DEBUG
  dstrc_enter("st_node_new");
#endif

  node = CCACALLOC(1,sizeof(ST_NODE));
  node->type = type;
  node->lhs = lhs;
  node->rhs = rhs;

#ifdef DEBUG
  dstrc_exit();
#endif

  return node;
}


static void st_node_delete(ST_NODE* node)
{
#ifdef DEBUG
  dstrc_enter("st_node_delete");
#endif

  if (node->lhs)
    st_node_delete(node->lhs);
  if (node->rhs)
    st_node_delete(node->rhs);
  CCAFREE(node);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief The types of tokens recognized by the lexer.

  \author u.kue
  \date 12/05
*/
/*----------------------------------------------------------------------*/
typedef enum _ST_TOKEN_TYPE
{
  tok_none,
  tok_done,
  tok_name,
  tok_int,
  tok_real,
  tok_add,
  tok_sub,
  tok_mul,
  tok_div,
  tok_mod,
  tok_pow,
  tok_lpar,
  tok_rpar
} ST_TOKEN_TYPE;


/*----------------------------------------------------------------------*/
/*!
  \brief state variables for the lexer and the parser.

  \author u.kue
  \date 12/05
*/
/*----------------------------------------------------------------------*/
static CHAR* lex_funct;
static INT lex_length;
static INT lex_pos;
static ST_TOKEN_TYPE lex_tok;
static CHAR* lex_string;
static INT lex_int;
static DOUBLE lex_real;

static DOUBLE parse_x;
static DOUBLE parse_y;
static DOUBLE parse_z;
static DOUBLE parse_t;


/*----------------------------------------------------------------------*/
/*!
  \brief there is a simple lexer that always gets the next token from
  the expression string.

  \author u.kue
  \date 12/05
*/
/*----------------------------------------------------------------------*/
static void lex_init(CHAR* funct)
{
#ifdef DEBUG
  dstrc_enter("lex_init");
#endif

  lex_funct = funct;
  lex_length = strlen(funct);
  lex_pos = 0;

#ifdef DEBUG
  dstrc_exit();
#endif
}


static int lex_getnext()
{
  if (lex_pos < lex_length)
  {
    /* Increment the counter and return the char at the old position. */
    return lex_funct[lex_pos++];
  }
  return EOF;
}


/*----------------------------------------------------------------------*/
/*!
  \brief The lexer's main function.

  Read the expression string char by char and set the static variables.

  \author u.kue
  \date 12/05
*/
/*----------------------------------------------------------------------*/
static void lexan()
{
  INT t;

#ifdef DEBUG
  dstrc_enter("lexan");
#endif

  for (;;)
  {
    t = lex_getnext();
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
      lex_tok = tok_done;
      goto end;
    }
    else
    {
      if (isdigit(t))
      {
	lex_string = &(lex_funct[lex_pos-1]);
	while (isdigit(t))
	{
	  t = lex_getnext();
	}
	if ((t != '.') && (t != 'E') && (t != 'e'))
	{
	  if (t != EOF)
	  {
	    lex_pos--;
	  }
	  lex_int = atoi(lex_string);
	  lex_tok = tok_int;
	  goto end;
	}
	if (t == '.')
	{
	  t = lex_getnext();
	  if (isdigit(t))
	  {
	    while (isdigit(t))
	    {
	      t = lex_getnext();
	    }
	  }
	  else
	  {
	    dserror("no digits after point at pos %d", lex_pos);
	  }
	}
	if ((t == 'E') || (t == 'e'))
	{
	  t = lex_getnext();
	  if ((t == '-') || (t == '+'))
	  {
	    t = lex_getnext();
	  }
	  if (isdigit(t))
	  {
	    while (isdigit(t))
	    {
	      t = lex_getnext();
	    }
	  }
	  else
	  {
	    dserror("no digits after exponent at pos %d", lex_pos);
	  }
	}
	if (t != EOF)
	{
	  lex_pos--;
	}
	lex_real = strtod(lex_string, NULL);
	lex_tok = tok_real;
	goto end;
      }
      else if (isalpha(t) || (t == '_'))
      {
	lex_string = &(lex_funct[lex_pos-1]);
	while (isalnum(t) || (t == '_'))
	{
	  t = lex_getnext();
	}
	if (t != EOF)
	{
	  lex_pos--;
	}
	lex_tok = tok_name;
	lex_int = &(lex_funct[lex_pos]) - lex_string;
	goto end;
      }
      else if (t == '+')
      {
	lex_tok = tok_add;
	goto end;
      }
      else if (t == '-')
      {
	lex_tok = tok_sub;
	goto end;
      }
      else if (t == '*')
      {
	lex_tok = tok_mul;
	goto end;
      }
      else if (t == '/')
      {
	lex_tok = tok_div;
	goto end;
      }
      else if (t == '^')
      {
	lex_tok = tok_pow;
	goto end;
      }
      else if (t == '(')
      {
	lex_tok = tok_lpar;
	goto end;
      }
      else if (t == ')')
      {
	lex_tok = tok_rpar;
	goto end;
      }
      else
      {
	if (t >= 32)
	  dserror("unexpected char '%c' at pos %d", t, lex_pos);
	else
	  dserror("unexpected char '%d' at pos %d", t, lex_pos);
	lex_tok = tok_none;
	goto end;
      }
    }
  }

end:

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*----------------------------------------------------------------------*/
/*!
  \brief The parser function.

  Here we build the syntax tree out of the tokens the lexer gives us.

  In order to get the operator execution order right there are
  different functions for each operator level.

  \author u.kue
  \date 12/05
*/
/*----------------------------------------------------------------------*/
static ST_NODE* parse_expr();
static ST_NODE* parse_pow();

static ST_NODE* parse_primary()
{
  ST_NODE* lhs;
  INT i;

#ifdef DEBUG
  dstrc_enter("parse_primary");
#endif

  switch (lex_tok)
  {
  case tok_lpar:
    lexan();
    lhs = parse_expr();
    if (lex_tok!=tok_rpar)
      dserror("')' expected");
    lexan();
    break;
  case tok_int:
    lhs = st_node_new(lt_number,NULL,NULL);
    lhs->v.number = lex_int;
    lexan();
    break;
  case tok_real:
    lhs = st_node_new(lt_number,NULL,NULL);
    lhs->v.number = lex_real;
    lexan();
    break;
  case tok_sub:
  {
    ST_NODE* rhs;
    lexan();
    /*rhs = parse_primary();*/
    rhs = parse_pow();
    lhs = st_node_new(lt_number,NULL,NULL);
    lhs->v.number = -1;
    lhs = st_node_new(lt_operator,lhs,rhs);
    lhs->v.operator = '*';
    break;
  }
  case tok_name:
    lhs = NULL;
    if (lex_int==1)
    {
      switch (lex_string[0])
      {
      case 'x':
      case 'y':
      case 'z':
      case 't':
	lhs = st_node_new(lt_variable,NULL,NULL);
	lhs->v.variable = lex_string[0];
	lexan();
	break;
  default:
	break;
      }
    }
    if (lhs==NULL)
    {
      if ((lex_int == 2) &&
	  (strncmp("pi", lex_string, lex_int)==0))
      {
	lhs = st_node_new(lt_number,NULL,NULL);
	lhs->v.number = PI;
	lexan();
	break;
      }
      else
      {
	for (i=0; function[i].name!=NULL; ++i)
	{
	  if ((strlen(function[i].name) >= lex_int) &&
	      (strncmp(function[i].name, lex_string, lex_int)==0))
	  {
	    lhs = st_node_new(lt_function,NULL,NULL);
	    lhs->v.function = function[i].func;
	    lexan();
	    if (lex_tok!=tok_lpar)
	      dserror("'(' expected after function name '%s'", function[i].name);
	    lexan();
	    lhs->lhs = parse_expr();
	    if (lex_tok!=tok_rpar)
	      dserror("')' expected");
	    lexan();
	    break;
	  }
	}
	if (lhs==NULL)
	{
	  lex_string[lex_int] = '\0';
	  dserror("unknown function '%s'", lex_string);
	}
      }
    }
    break;
  default:
    dserror("unexpected token %d", lex_tok);
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return lhs;
}

static ST_NODE* parse_pow()
{
  ST_NODE* lhs;
  ST_NODE* rhs;

#ifdef DEBUG
  dstrc_enter("parse_pow");
#endif

  lhs = parse_primary();
  for (;;)
  {
    if (lex_tok==tok_pow)
    {
      lexan();
      rhs = parse_primary();
      lhs = st_node_new(lt_operator,lhs,rhs);
      lhs->v.operator = '^';
    }
    else
    {
      break;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return lhs;
}

static ST_NODE* parse_term()
{
  ST_NODE* lhs;
  ST_NODE* rhs;

#ifdef DEBUG
  dstrc_enter("parse_term");
#endif

  lhs = parse_pow();
  for (;;)
  {
    if (lex_tok==tok_mul)
    {
      lexan();
      rhs = parse_pow();
      lhs = st_node_new(lt_operator,lhs,rhs);
      lhs->v.operator = '*';
    }
    else if (lex_tok==tok_div)
    {
      lexan();
      rhs = parse_pow();
      lhs = st_node_new(lt_operator,lhs,rhs);
      lhs->v.operator = '/';
    }
    else
    {
      break;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return lhs;
}

static ST_NODE* parse_expr()
{
  ST_NODE* lhs;
  ST_NODE* rhs;

#ifdef DEBUG
  dstrc_enter("parse_expr");
#endif

  lhs = parse_term();
  for (;;)
  {
    if (lex_tok==tok_add)
    {
      lexan();
      rhs = parse_term();
      lhs = st_node_new(lt_operator,lhs,rhs);
      lhs->v.operator = '+';
    }
    else if (lex_tok==tok_sub)
    {
      lexan();
      rhs = parse_term();
      lhs = st_node_new(lt_operator,lhs,rhs);
      lhs->v.operator = '-';
    }
    else
    {
      break;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return lhs;
}


/*----------------------------------------------------------------------*/
/*!
  \brief global parser function

  Enter your expression string and get an (opaque) syntax tree. Your
  are not supposed to interpret the syntax tree yourself. Just hand it
  into pss_evaluate!

  \author u.kue
  \date 12/05
*/
/*----------------------------------------------------------------------*/
struct _ST_NODE* pss_parse(CHAR* funct)
{
  ST_NODE* head;

#ifdef DEBUG
  dstrc_enter("pss_parse");
#endif

  lex_init(funct);
  lexan();
  head = parse_expr();

#ifdef DEBUG
  dstrc_exit();
#endif

  return head;
}


/*----------------------------------------------------------------------*/
/*!
  \brief remove the syntax tree

  \author u.kue
  \date 12/05
*/
/*----------------------------------------------------------------------*/
void pss_parse_cleanup(struct _ST_NODE* head)
{
#ifdef DEBUG
  dstrc_enter("pss_parse_cleanup");
#endif

  st_node_delete(head);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief the internal interpreter.

  We decent the syntax tree recursively and evaluate the expression on
  the way.

  \author u.kue
  \date 12/05
*/
/*----------------------------------------------------------------------*/
static DOUBLE pss_interpret(ST_NODE* node)
{
  DOUBLE res;

#ifdef DEBUG
  dstrc_enter("pss_interpret");
#endif

  switch (node->type)
  {
  case lt_number:
    res = node->v.number;
    break;
  case lt_operator:
  {
    DOUBLE lhs;
    DOUBLE rhs;

    lhs = pss_interpret(node->lhs);
    rhs = pss_interpret(node->rhs);
    switch (node->v.operator)
    {
    case '+':
      res = lhs+rhs;
      break;
    case '-':
      res = lhs-rhs;
      break;
    case '*':
      res = lhs*rhs;
      break;
    case '/':
      /* check for rhs==0.0? */
      res = lhs/rhs;
      break;
    case '^':
      res = pow(lhs,rhs);
      break;
    default:
      dserror("unsupported operator '%c'", node->v.operator);
    }
    break;
  }
  case lt_variable:
    switch (node->v.variable)
    {
    case 'x':
      res = parse_x;
      break;
    case 'y':
      res = parse_y;
      break;
    case 'z':
      res = parse_z;
      break;
    case 't':
      res = parse_t;
      break;
    default:
      dserror("unknown variable '%c'", node->v.variable);
    }
    break;
  case lt_function:
  {
    DOUBLE arg;
    arg = pss_interpret(node->lhs);
    res = node->v.function(arg);
    break;
  }
  default:
    dserror("unknown syntax tree node type");
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return res;
}


/*----------------------------------------------------------------------*/
/*!
  \brief global interpreter call for space functions.

  \author u.kue
  \date 12/05
*/
/*----------------------------------------------------------------------*/
DOUBLE pss_evaluate_funct(struct _ST_NODE* funct, DOUBLE x, DOUBLE y, DOUBLE z)
{
  DOUBLE res;

#ifdef DEBUG
  dstrc_enter("pss_evaluate");
#endif

  parse_x = x;
  parse_y = y;
  parse_z = z;
  parse_t = 0;
  res = pss_interpret(funct);

#ifdef DEBUG
  dstrc_exit();
#endif

  return res;
}

/*----------------------------------------------------------------------*/
/*!
  \brief global interpreter call for time curves.

  \author u.kue
  \date 12/05
*/
/*----------------------------------------------------------------------*/
DOUBLE pss_evaluate_curve(struct _ST_NODE* funct, DOUBLE t)
{
  DOUBLE res;

#ifdef DEBUG
  dstrc_enter("pss_evaluate_curve");
#endif

  parse_x = 0;
  parse_y = 0;
  parse_z = 0;
  parse_t = t;
  res = pss_interpret(funct);

#ifdef DEBUG
  dstrc_exit();
#endif

  return res;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void pss_parser_print(FILE* out, struct _ST_NODE* node)
{
  switch (node->type)
  {
  case lt_number:
    fprintf(out, "%g", node->v.number);
    break;
  case lt_operator:
  {
    fprintf(out,"(");
    pss_parser_print(out,node->lhs);
    fprintf(out,")%c(",node->v.operator);
    pss_parser_print(out,node->rhs);
    fprintf(out,")");
    break;
  }
  case lt_variable:
    fprintf(out,"%c",node->v.variable);
    break;
  case lt_function:
  {
    int i;
    for (i=0; function[i].name!=NULL; ++i)
    {
      if (function[i].func==node->v.function)
      {
        fprintf(out,function[i].name);
        break;
      }
    }
    if (function[i].name==NULL)
      dserror("unknown function");
    fprintf(out,"(");
    pss_parser_print(out,node->lhs);
    fprintf(out,")");
    break;
  }
  default:
    dserror("unknown syntax tree node type");
  }
}
