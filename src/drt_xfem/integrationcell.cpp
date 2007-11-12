/*!
\file integrationcell.cpp

\brief integration cell

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/

#ifdef CCADISCRET

#include "integrationcell.H"
#include <string>
#include <sstream>

using namespace std;
using namespace XFEM;

#define MFOREACH(TYPE,VAL,VALS) for( TYPE::iterator VAL = VALS.begin(); VAL != VALS.end(); ++VAL )
#define MCONST_FOREACH(TYPE,VAL,VALS) for( TYPE::const_iterator VAL = VALS.begin(); VAL != VALS.end(); ++VAL )
#define MPFOREACH(TYPE,VAL,VALS) for( TYPE::const_iterator VAL = VALS->begin(); VAL != VALS->end(); ++VAL )

//
//  ctor
//
IntCell::IntCell(const vector< vector<double> > domainCoordinates) :
	domainCoordinates_(domainCoordinates)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor                                                mwgee 11/06|
 *----------------------------------------------------------------------*/
IntCell::IntCell(const IntCell& old) : 
	domainCoordinates_(old.domainCoordinates_)
{
    return;   
}
 
/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
IntCell::~IntCell()
{
  return;
}

////
////  get coordinates
////
//vector< vector<double> >  IntCell::GetDomainCoord() const
//{
//    return domainCoordinates_;   
//}

//
// virtual Print method
//
std::string IntCell::Print() const
{
  return "";
}



//
//  ctor
//
DomainIntCell::DomainIntCell(const vector< vector<double> > coordinates) :
IntCell(coordinates)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor                                                mwgee 11/06|
 *----------------------------------------------------------------------*/
DomainIntCell::DomainIntCell(const DomainIntCell& old) : 
IntCell(old)
{
    return;   
}
     
//
//  get coordinates
//
vector< vector<double> >  DomainIntCell::GetDomainCoord() const
{
    return domainCoordinates_;   
}

string DomainIntCell::Print() const
{
    stringstream s;
    s << "DomainIntCell" << endl;
    MCONST_FOREACH(vector< vector<double> >, coordinate, domainCoordinates_)
    {
        s << "[";
        MPFOREACH(vector<double>, val, coordinate)
        {
            s << *val << " ";
        };
        s << "]" << endl;
    };
    return s.str();
}


//
//  ctor
//
BoundaryIntCell::BoundaryIntCell(const vector< vector<double> > coordinates) :
IntCell(coordinates)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor                                                mwgee 11/06|
 *----------------------------------------------------------------------*/
BoundaryIntCell::BoundaryIntCell(const BoundaryIntCell& old) : 
IntCell(old)
{
    return;   
}
     
//
//  get coordinates
//
vector< vector<double> >  BoundaryIntCell::GetDomainCoord() const
{
    return domainCoordinates_;   
}

//
//  get coordinates
//
vector< vector<double> >  BoundaryIntCell::GetBoundaryCoord() const
{
    return boundaryCoordinates_;   
}


string BoundaryIntCell::Print() const
{
    stringstream s;
    s << "BoundaryIntCell" << endl;
    MCONST_FOREACH(vector< vector<double> >, coordinate, domainCoordinates_)
    {
        s << "[";
        MPFOREACH(vector<double>, val, coordinate)
        {
            s << *val << " ";
        };
        s << "]" << endl;
    };
    return s.str();
}

#endif  // #ifdef CCADISCRET
