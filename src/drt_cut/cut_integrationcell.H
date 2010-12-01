#ifndef CUT_INTEGRATIONCELL_H
#define CUT_INTEGRATIONCELL_H

#include <vector>
#include <set>

#include "cut_element.H"
#include "cut_boundingbox.H"

namespace GEO
{
namespace CUT
{

/// Base of (non-tetgen) integration cells
/*!
  Die sich selbst in Unterzellen zerlegende Zelle. Das Kriterium für die
  Zerlegung ist jeweils verschieden. Es gibt die

  - Hauptzellen, die dem Element gleichen und bei der die Zerlegung ihren
  Anfang nimmt

  - leere Zellen, die lediglich Raum im Element füllen und sich selbst nicht
  weiter unterteilen

  - geschnittene Zellen, die den Schnitt enthalten und aus diesem Grund weiter
  unterteilt werden müssen

  Die Unterteilung erfolgt jeweils anhand von Topologie-Informationen, die
  Geometrie braucht hier nicht mehr ausgewertet werden, da die Elemente selbst
  bereits geschnitten sind. Die Aufteilung findet in Elementkoordinaten statt.
 */
class IntegrationCell
{
public:

  virtual ~IntegrationCell() {}
};


/// empty cells that have no cut nor any subcells
class EmptyIntegrationCell : public IntegrationCell
{
public:
};

class Hex8EmptyIntegrationCell : public EmptyIntegrationCell
{
public:
};

class Tet4EmptyIntegrationCell : public EmptyIntegrationCell
{
public:
};

class Wedge6EmptyIntegrationCell : public EmptyIntegrationCell
{
public:
};

class Pyramid5EmptyIntegrationCell : public EmptyIntegrationCell
{
public:
};


/// contains the cut
class CutIntegrationCell : public IntegrationCell
{
public:
};

class Hex8CutIntegrationCell : public CutIntegrationCell
{
public:
};

class Tet4CutIntegrationCell : public CutIntegrationCell
{
public:
};

class Wedge6CutIntegrationCell : public CutIntegrationCell
{
public:
};

class Pyramid5CutIntegrationCell : public CutIntegrationCell
{
public:
};


/// matches the whole element
class MainIntegrationCell : public IntegrationCell
{
public:

  virtual Element * element() = 0;

protected:

  std::vector<Teuchos::RCP<EmptyIntegrationCell> > emptycells_;
};

class Hex8MainIntegrationCell : public MainIntegrationCell
{
public:

  Hex8MainIntegrationCell( ConcreteElement<DRT::Element::hex8> * e );

  ConcreteElement<DRT::Element::hex8> * element() { return element_; }

private:
  ConcreteElement<DRT::Element::hex8> * element_;
  Teuchos::RCP<Hex8CutIntegrationCell> cut_cell_;
};

class Tet4MainIntegrationCell : public MainIntegrationCell
{
public:

  Tet4MainIntegrationCell( ConcreteElement<DRT::Element::tet4> * e );

  ConcreteElement<DRT::Element::tet4> * element() { return element_; }

private:
  ConcreteElement<DRT::Element::tet4> * element_;
};

class Wedge6MainIntegrationCell : public MainIntegrationCell
{
public:

  Wedge6MainIntegrationCell( ConcreteElement<DRT::Element::wedge6> * e );

  ConcreteElement<DRT::Element::wedge6> * element() { return element_; }

private:
  ConcreteElement<DRT::Element::wedge6> * element_;
};

class Pyramid5MainIntegrationCell : public MainIntegrationCell
{
public:

  Pyramid5MainIntegrationCell( ConcreteElement<DRT::Element::pyramid5> * e );

  ConcreteElement<DRT::Element::pyramid5> * element() { return element_; }

private:
  ConcreteElement<DRT::Element::pyramid5> * element_;
};


}
}

#endif
