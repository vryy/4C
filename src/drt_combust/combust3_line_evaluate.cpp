/*!----------------------------------------------------------------------
\file combust3_line_evaluate.cpp
\brief

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

*----------------------------------------------------------------------*/

#include "combust3.H"


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Combust3Line::Evaluate(
        Teuchos::ParameterList&   params,
        DRT::Discretization&      discretization,
        std::vector<int>&         lm,
        Epetra_SerialDenseMatrix& elemat1,
        Epetra_SerialDenseMatrix& elemat2,
        Epetra_SerialDenseVector& elevec1,
        Epetra_SerialDenseVector& elevec2,
        Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Combust3Line::ActionType act = Combust3Line::none;
  std::string action = params.get<std::string>("action","none");
  if (action == "none") dserror("No action supplied");
  else dserror("Unknown type of action for Combust3_Line");

  switch(act)
  {
  default:
    dserror("Unknown type of action for Combust3_Line");
  } // end of switch(act)

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Combust3Line::EvaluateNeumann(
        Teuchos::ParameterList&   params,
        DRT::Discretization&      discretization,
        DRT::Condition&           condition,
        std::vector<int>&         lm,
        Epetra_SerialDenseVector& elevec1,
        Epetra_SerialDenseMatrix* elemat1)
{
  std::cout << "/!\\ warning === Neumann boundary conditions in XFEM problems are not implemented in a general way!" << std::endl;
  dserror("Neumann condition on line not implemented");
  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Combust3Line::LocationVector(
    const Discretization&   dis,
    LocationArray&          la,
    bool                    doDirichlet,
    const std::string&      condstring,
    Teuchos::ParameterList& params
    ) const
{
  DRT::ELEMENTS::Combust3Line::ActionType act = Combust3Line::none;
  std::string action = params.get<std::string>("action","none");
  if (action == "none") dserror("No action supplied");

  switch(act)
  {
  default:
    DRT::Element::LocationVector(dis,la,doDirichlet);
    break;
  }
  return;
}
