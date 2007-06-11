/*!----------------------------------------------------------------------
\file fluid3_surface.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fluid3_xfem.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"

extern "C"
{
#include "../headers/standardtypes.h"
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 01/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3Surface::XFluid3Surface(int id,
                                              int owner,
                                              int nnode,
                                              const int* nodeids,
                                              DRT::Node** nodes,
                                              DRT::Elements::XFluid3* parent,
                                              const int lsurface) :
DRT::Element(id,element_xfluid3surface,owner),
parent_(parent),
lsurface_(lsurface)
{
  lines_.resize(0);
  lineptrs_.resize(0);	
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3Surface::XFluid3Surface(const DRT::Elements::XFluid3Surface& old) :
DRT::Element(old),
parent_(old.parent_),
lsurface_(old.lsurface_),
lines_(old.lines_),
lineptrs_(old.lineptrs_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::XFluid3Surface::Clone() const
{
  DRT::Elements::XFluid3Surface* newelement = new DRT::Elements::XFluid3Surface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::XFluid3Surface::Shape() const
{
  switch (NumNode())
  {
		case 3: return tri3;
  		case 4: return quad4;
  		case 6: return tri6;
  		case 8: return quad8;
  		case 9: return quad9;
  		default:
    		dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Surface::Pack(vector<char>& data) const
{
  data.resize(0);
  dserror("this Fluid3Surface element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Surface::Unpack(const vector<char>& data)
{
  dserror("this XFluid3Surface element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3Surface::~XFluid3Surface()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Surface::Print(ostream& os) const
{
  os << "XFluid3Surface ";
  Element::Print(os);
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             gammi 04/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::XFluid3Surface::Lines()
{

	const int nline   = NumLine();
  	const int numnode = NumNode();
  	lines_.resize(nline);
  	lineptrs_.resize(nline);
  	int nodeids[100];
  	DRT::Node* nodes[100];

  	/* triangle elements */
  	if (nline==3)
  	{
    /* mind: The displayed element is not the reference element!
     *
     *               2
     *                X
     *                |\
     *                | \       edge1
     *   edge2       6o  o5    (line1)
     *  (line2)       |   \
     *                |    \
     *                X--o--X
     *               0   4   1
     *
     *             edge0
     *            (line0)
     *
     *
     *      X nodes for tri3
     *      o addotional nodes for tri6
     *                                                                       */
    /* linear triangles*/
		if (numnode==3)
    	{
      	/* first edge */
      	// set node id's
      	nodeids[0] = NodeIds()[0];
      	nodeids[1] = NodeIds()[1];
      	// get nodes
      	nodes[0] = Nodes()[0];
      	nodes[1] = Nodes()[1];
      	// create the line and get the line pointer
      	lines_[0] =
        		rcp(new DRT::Elements::XFluid3Line(0,Owner(),2,nodeids,nodes,this,NULL,0));
      	lineptrs_[0] = lines_[0].get();

      	/* second edge */
      	// set node id's
     	 	nodeids[0] = NodeIds()[1];
      	nodeids[1] = NodeIds()[2];
      	// get nodes
      	nodes[0] = Nodes()[1];
      	nodes[1] = Nodes()[2];
      	// create the line and get the line pointer
      	lines_[1] =
        		rcp(new DRT::Elements::XFluid3Line(1,Owner(),2,nodeids,nodes,this,NULL,1));
      	lineptrs_[1] = lines_[1].get();

	      /* third edge */
	      // set node id's
	      nodeids[0] = NodeIds()[2];
	      nodeids[1] = NodeIds()[0];
	      // get nodes
	      nodes[0] = Nodes()[2];
	      nodes[1] = Nodes()[0];
	      // create the line and get the line pointer
	      lines_[2] =
	        rcp(new DRT::Elements::XFluid3Line(2,Owner(),2,nodeids,nodes,this,NULL,2));
	      lineptrs_[2] = lines_[2].get();
    	}
    	/* quadratic triangles*/
    	else if (numnode==6)
    	{
      	dserror("TRI6 lines not implemented.");
    	}
  	}
  	/* quad elements*/
  	else if (nline==4)
  	{
    /* mind: The displayed element is not the reference element!
     *
     *                  
     *                  edge2	  
     *		       (line2)
     *           
     *              3     6     2
     *               X----o----X
     *               |         |
     *               |         |   
     *   edge3      7o    O    o5      edge1	  
     *  (line3)      |    8    |      (line1)
     *               |         |
     *               X----o----X
     *              0     4     1
     *
     *                  edge0
     *                 (line0)
     *
     *
     *      X nodes for quad
     *      o addotional nodes for quad8
     *      O addotional nodes for full quadratic quad9
     *      
     *                                                                 */
    	if (numnode==4)
    	{
	      /* first edge */
	      // set node id's
	      nodeids[0] = NodeIds()[0];
	      nodeids[1] = NodeIds()[1];
	      // get nodes
	      nodes[0] = Nodes()[0];
	      nodes[1] = Nodes()[1];
	      // create the line and get the line pointer
	      lines_[0] =
	        rcp(new DRT::Elements::XFluid3Line(0,Owner(),2,nodeids,nodes,this,NULL,0));
	      lineptrs_[0] = lines_[0].get();
	
	      /* second edge */
	      // set node id's
	      nodeids[0] = NodeIds()[1];
	      nodeids[1] = NodeIds()[2];
	      // get nodes
	      nodes[0] = Nodes()[1];
	      nodes[1] = Nodes()[2];
	      // create the line and get the line pointer
	      lines_[1] =
	        rcp(new DRT::Elements::XFluid3Line(1,Owner(),2,nodeids,nodes,this,NULL,1));
	      lineptrs_[1] = lines_[1].get();
	
	      /* third edge */
	      // set node id's
	      nodeids[0] = NodeIds()[2];
	      nodeids[1] = NodeIds()[3];
	      // get nodes
	      nodes[0] = Nodes()[2];
	      nodes[1] = Nodes()[3];
	      // create the line and get the line pointer
	      lines_[2] =
	        rcp(new DRT::Elements::XFluid3Line(2,Owner(),2,nodeids,nodes,this,NULL,2));
	      lineptrs_[2] = lines_[2].get();

	      /* fourth edge */
	      // set node id's
	      nodeids[0] = NodeIds()[3];
	      nodeids[1] = NodeIds()[0];
	      // get nodes
	      nodes[0] = Nodes()[3];
	      nodes[1] = Nodes()[0];
	      // create the line and get the line pointer
	      lines_[3] =
	        rcp(new DRT::Elements::XFluid3Line(3,Owner(),2,nodeids,nodes,this,NULL,3));
	      lineptrs_[3] = lines_[3].get();
	      
    	}
    	else if (numnode==8)
    	{
      	dserror("quad8 lines not implemented.");
    	}
    	else if (numnode==9)
    	{
      	/* first edge */
	      // set node id's
	      nodeids[0] = NodeIds()[0];
	      nodeids[1] = NodeIds()[4];
	      nodeids[2] = NodeIds()[1];
	      // get nodes
	      nodes[0] = Nodes()[0];
	      nodes[1] = Nodes()[4];
	      nodes[2] = Nodes()[1];
	      // create the line and get the line pointer
	      lines_[0] =
	        rcp(new DRT::Elements::XFluid3Line(0,Owner(),3,nodeids,nodes,this,NULL,0));
	      lineptrs_[0] = lines_[0].get();
	
	      /* second edge */
	      // set node id's
	      nodeids[0] = NodeIds()[1];
	      nodeids[1] = NodeIds()[5];
	      nodeids[2] = NodeIds()[2];
	      // get nodes
	      nodes[0] = Nodes()[1];
	      nodes[1] = Nodes()[5];
	      nodes[2] = Nodes()[2];
	      // create the line and get the line pointer
	      lines_[1] =
	        rcp(new DRT::Elements::XFluid3Line(1,Owner(),3,nodeids,nodes,this,NULL,1));
	      lineptrs_[1] = lines_[1].get();
	
	      /* third edge */
	      // set node id's
	      nodeids[0] = NodeIds()[2];
	      nodeids[1] = NodeIds()[6];
	      nodeids[2] = NodeIds()[3];
	      // get nodes
	      nodes[0] = Nodes()[2];
	      nodes[1] = Nodes()[6];
	      nodes[2] = Nodes()[3];
	      // create the line and get the line pointer
	      lines_[2] =
	        rcp(new DRT::Elements::XFluid3Line(2,Owner(),3,nodeids,nodes,this,NULL,2));
	      lineptrs_[2] = lines_[2].get();

	      /* fourth edge */
	      // set node id's
	      nodeids[0] = NodeIds()[3];
	      nodeids[1] = NodeIds()[7];
	      nodeids[2] = NodeIds()[0];
	      // get nodes
	      nodes[0] = Nodes()[3];
	      nodes[1] = Nodes()[7];
	      nodes[2] = Nodes()[0];
	      // create the line and get the line pointer
	      lines_[3] =
	        rcp(new DRT::Elements::XFluid3Line(3,Owner(),3,nodeids,nodes,this,NULL,3));
	      lineptrs_[3] = lines_[3].get();
    	}
    	else dserror("Number of nodes not supported");
  	}
  	else dserror("Number of lines not supported");

  	return (DRT::Element**)(&(lineptrs_[0]));
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
