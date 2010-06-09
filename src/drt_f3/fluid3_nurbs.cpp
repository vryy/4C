/*!----------------------------------------------------------------------
\file fluid3_nurbs.cpp

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "fluid3_nurbs.H"
#include "../drt_lib/drt_utils.H"

DRT::ELEMENTS::NURBS::Fluid3NurbsType DRT::ELEMENTS::NURBS::Fluid3NurbsType::instance_;


DRT::ParObject* DRT::ELEMENTS::NURBS::Fluid3NurbsType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::NURBS::Fluid3Nurbs* object = new DRT::ELEMENTS::NURBS::Fluid3Nurbs(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NURBS::Fluid3NurbsType::Create( const string eletype,
                                                                          const string eledistype,
                                                                          const int id,
                                                                          const int owner )
{
  if ( eletype=="FLUID3" )
  {
    if ( eledistype=="NURBS8" || eledistype=="NURBS27")
    {
      return rcp(new DRT::ELEMENTS::NURBS::Fluid3Nurbs(id,owner));
    }
  }
  return Teuchos::null;
}


void DRT::ELEMENTS::NURBS::Fluid3NurbsType::NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = dwele->NumDofPerNode(*(dwele->Nodes()[0]));
  dimns = numdf;
  nv = numdf-1;
  np = 1;
}


void DRT::ELEMENTS::NURBS::Fluid3NurbsType::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeFluid3DNullSpace( dis, ns, x0, numdf, dimns );
}


DRT::ELEMENTS::NURBS::Fluid3NurbsSurfaceType DRT::ELEMENTS::NURBS::Fluid3NurbsSurfaceType::instance_;


/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 05/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Fluid3Nurbs::Fluid3Nurbs(int id, int owner) :
DRT::ELEMENTS::Fluid3::Fluid3(id,owner)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 05/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Fluid3Nurbs::Fluid3Nurbs
(const DRT::ELEMENTS::NURBS::Fluid3Nurbs& old) :
DRT::ELEMENTS::Fluid3::Fluid3(old)
{
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 05/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Fluid3Nurbs::~Fluid3Nurbs()
{
  return;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          gammi 02/09 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::NURBS::Fluid3Nurbs::Clone() const
{
  DRT::ELEMENTS::NURBS::Fluid3Nurbs* newelement
    =
    new DRT::ELEMENTS::NURBS::Fluid3Nurbs(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          gammi 05/08 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::NURBS::Fluid3Nurbs::Shape() const
{
  switch (NumNode())
  {
  case   8: return nurbs8;
  case  27: return nurbs27;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }

  return dis_none;
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                         gammi 02/09 |
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::NURBS::Fluid3Nurbs::Surfaces()
{

  return
    DRT::UTILS::ElementBoundaryFactory
    <Fluid3NurbsSurface,Fluid3Nurbs>
    (
      DRT::UTILS::buildSurfaces,
      this
    );
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 02/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Fluid3NurbsSurface::Fluid3NurbsSurface(
  int                    id      ,
  int                    owner   ,
  int                    nnode   ,
  const int*             nodeids ,
  DRT::Node**            nodes   ,
  DRT::ELEMENTS::Fluid3* parent  ,
  const int              lsurface
  )
:
DRT::ELEMENTS::Fluid3Boundary(id,owner,nnode,nodeids,nodes,parent,lsurface)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 02/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Fluid3NurbsSurface::Fluid3NurbsSurface
(const DRT::ELEMENTS::NURBS::Fluid3NurbsSurface& old) :
DRT::ELEMENTS::Fluid3Boundary::Fluid3Boundary(old)
{
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 02/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Fluid3NurbsSurface::~Fluid3NurbsSurface()
{
  return;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          gammi 02/09 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::NURBS::Fluid3NurbsSurface::Clone() const
{
  DRT::ELEMENTS::NURBS::Fluid3NurbsSurface* newelement
    =
    new DRT::ELEMENTS::NURBS::Fluid3NurbsSurface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          gammi 02/09 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::NURBS::Fluid3NurbsSurface::Shape() const
{
  switch (NumNode())
  {
  case   4: return nurbs4;
  case   9: return nurbs9;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }

  return dis_none;
}

#endif

#endif //D_FLUID3
