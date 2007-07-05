/*!
\file integrationcell.cpp

\brief integration cell

<pre>
Maintainer: 
</pre>
*/

#ifdef D_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "integrationcell.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 |  id             (in)  this integrationcell's global id               |
 *----------------------------------------------------------------------*/
Integrationcell::Integrationcell(	int id, 
									std::vector< std::vector<double> > coordinates)
{
	id_ = id;
    coordinates_ = coordinates;
    return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 |  id             (in)  this integrationcell's global id               |
 *----------------------------------------------------------------------*/
Integrationcell::Integrationcell(const Integrationcell& old) 
{
    id_ = old.id_;
    coordinates_ = old.coordinates_;
    return;   
}


     
/*----------------------------------------------------------------------*
 |  get coordinates			                                 mwgee 11/06|  
 *----------------------------------------------------------------------*/
std::vector< std::vector<double> >  Integrationcell::GetCoord() 
{
    return coordinates_;   
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_XFEM


