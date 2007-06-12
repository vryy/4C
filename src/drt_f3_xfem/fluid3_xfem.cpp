/*!----------------------------------------------------------------------
\file fluid3_xfem.cpp
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
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3::XFluid3(int id, int owner) :
DRT::Element(id,element_xfluid3,owner),
material_(0),
is_ale_(false),
data_()
{
    gaussrule_ = hex_27point;
    surfaces_.resize(0);
    surfaceptrs_.resize(0);
    lines_.resize(0);
    lineptrs_.resize(0);
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3::XFluid3(const DRT::Elements::XFluid3& old) :
DRT::Element(old),
material_(old.material_),
is_ale_(old.is_ale_),
data_(old.data_),
surfaces_(old.surfaces_),
surfaceptrs_(old.surfaceptrs_),
lines_(old.lines_),
lineptrs_(old.lineptrs_)
{
    gaussrule_ = old.gaussrule_;
    return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Fluid3 and return pointer to it (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::XFluid3::Clone() const
{
    DRT::Elements::XFluid3* newelement = new DRT::Elements::XFluid3(*this);
    return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::XFluid3::Shape() const
{
    switch (NumNode())
    {
    case  4: return tet4;
    case  8: return hex8;
    case 10: return tet10;
    case 20: return hex20;
    case 27: return hex27;
    default:
        dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3::Pack(vector<char>& data) const
{
    data.resize(0);

    // pack type of this instance of ParObject
    int type = UniqueParObjectId();
    AddtoPack(data,type);
    // add base class Element
    vector<char> basedata(0);
    Element::Pack(basedata);
    AddtoPack(data,basedata);
    // ngp_
    AddtoPack(data,gaussrule_);
    // material_
    AddtoPack(data,material_);
    // is_ale_
    AddtoPack(data,is_ale_);
    // data_
    vector<char> tmp(0);
    data_.Pack(tmp);
    AddtoPack(data,tmp);

    return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3::Unpack(const vector<char>& data)
{
    int position = 0;
    // extract type
    int type = 0;
    ExtractfromPack(position,data,type);
    if (type != UniqueParObjectId()) dserror("wrong instance type data");
    // extract base class Element
    vector<char> basedata(0);
    ExtractfromPack(position,data,basedata);
    Element::Unpack(basedata);
  
    ExtractfromPack(position,data,gaussrule_);
    ExtractfromPack(position,data,material_);
    ExtractfromPack(position,data,is_ale_);
    // data_
    vector<char> tmp(0);
    ExtractfromPack(position,data,tmp);
    data_.Unpack(tmp);

    if (position != (int)data.size())
        dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
    return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3::~XFluid3()
{
    return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3::Print(ostream& os) const
{
    os << "XFluid3 ";
    Element::Print(os);
    cout << endl;
    cout << data_;
    return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Fluid3Register (public)              mwgee 12/06|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::Elements::XFluid3::ElementRegister() const
{
    return rcp(new DRT::Elements::XFluid3Register(Type()));
}


/*----------------------------------------------------------------------*
 |  get vector of lines (length 1) (public)                g.bau 03/07|
 *----------------------------------------------------------------------*/ 
DRT::Element** DRT::Elements::XFluid3::Lines()
{
	
	const int numnode = NumNode();
	const int nline = NumLine();
	lines_.resize(nline);
	lineptrs_.resize(nline);
	int nodeids[100];
	DRT::Node* nodes[100];
	
	if(numnode == 8)
	{
	  	nodeids[0] = NodeIds()[0];
	  	nodeids[1] = NodeIds()[1];
	  	nodes[0] = Nodes()[0];
	  	nodes[1] = Nodes()[1];
	  	lines_[0] =
	    	rcp(new DRT::Elements::XFluid3Line(0,Owner(),2,nodeids,nodes,NULL,this,0));
	  	lineptrs_[0] = lines_[0].get();
	
	  	nodeids[0] = NodeIds()[1];
	  	nodeids[1] = NodeIds()[2];
	  	nodes[0] = Nodes()[1];
	  	nodes[1] = Nodes()[2];
	  	lines_[1] =
	    	rcp(new DRT::Elements::XFluid3Line(1,Owner(),2,nodeids,nodes,NULL,this,1));
	  	lineptrs_[1] = lines_[1].get();
	
	  	nodeids[0] = NodeIds()[2];
	  	nodeids[1] = NodeIds()[3];
	  	nodes[0] = Nodes()[2];
	  	nodes[1] = Nodes()[3];
	  	lines_[2] =
	    	rcp(new DRT::Elements::XFluid3Line(2,Owner(),2,nodeids,nodes,NULL,this,2));
	  	lineptrs_[2] = lines_[2].get();
	
	  	nodeids[0] = NodeIds()[3];
	  	nodeids[1] = NodeIds()[0];
	  	nodes[0] = Nodes()[3];
	  	nodes[1] = Nodes()[0];
	  	lines_[3] =
	    	rcp(new DRT::Elements::XFluid3Line(3,Owner(),2,nodeids,nodes,NULL,this,3));
	  	lineptrs_[3] = lines_[3].get();
	
	  	nodeids[0] = NodeIds()[0];
	  	nodeids[1] = NodeIds()[4];
	  	nodes[0] = Nodes()[0];
	  	nodes[1] = Nodes()[4];
	  	lines_[4] =
	    	rcp(new DRT::Elements::XFluid3Line(4,Owner(),2,nodeids,nodes,NULL,this,4));
	  	lineptrs_[4] = lines_[4].get();
	
	  	nodeids[0] = NodeIds()[1];
	  	nodeids[1] = NodeIds()[5];
	  	nodes[0] = Nodes()[1];
	  	nodes[1] = Nodes()[5];
	  	lines_[5] =
	    	rcp(new DRT::Elements::XFluid3Line(5,Owner(),2,nodeids,nodes,NULL,this,5));
	  	lineptrs_[5] = lines_[5].get();
	
	  	nodeids[0] = NodeIds()[2];
	  	nodeids[1] = NodeIds()[6];
	  	nodes[0] = Nodes()[2];
	  	nodes[1] = Nodes()[6];
	  	lines_[6] =
	    	rcp(new DRT::Elements::XFluid3Line(6,Owner(),2,nodeids,nodes,NULL,this,6));
	  	lineptrs_[6] = lines_[6].get();
	
	  	nodeids[0] = NodeIds()[3];
	  	nodeids[1] = NodeIds()[7];
	  	nodes[0] = Nodes()[3];
	  	nodes[1] = Nodes()[7];
	  	lines_[7] =
	    	rcp(new DRT::Elements::XFluid3Line(7,Owner(),2,nodeids,nodes,NULL,this,7));
	  	lineptrs_[7] = lines_[7].get();
	
	  	nodeids[0] = NodeIds()[4];
	  	nodeids[1] = NodeIds()[5];
	  	nodes[0] = Nodes()[4];
	  	nodes[1] = Nodes()[5];
	  	lines_[8] =
	    	rcp(new DRT::Elements::XFluid3Line(8,Owner(),2,nodeids,nodes,NULL,this,8));
	  	lineptrs_[8] = lines_[8].get();
	
	  	nodeids[0] = NodeIds()[5];
	  	nodeids[1] = NodeIds()[6];
	  	nodes[0] = Nodes()[5];
	  	nodes[1] = Nodes()[6];
	  	lines_[9] =
	    	rcp(new DRT::Elements::XFluid3Line(9,Owner(),2,nodeids,nodes,NULL,this,9));
	  	lineptrs_[9] = lines_[9].get();
	
	  	nodeids[0] = NodeIds()[6];
	  	nodeids[1] = NodeIds()[7];
	  	nodes[0] = Nodes()[6];
	  	nodes[1] = Nodes()[7];
	  	lines_[10] =
	    	rcp(new DRT::Elements::XFluid3Line(10,Owner(),2,nodeids,nodes,NULL,this,10));
	  	lineptrs_[10] = lines_[10].get();
	
	  	nodeids[0] = NodeIds()[7];
	  	nodeids[1] = NodeIds()[4];
	  	nodes[0] = Nodes()[7];
	  	nodes[1] = Nodes()[4];
	  	lines_[11] =
	    	rcp(new DRT::Elements::XFluid3Line(11,Owner(),2,nodeids,nodes,NULL,this,11));
	  	lineptrs_[11] = lines_[11].get();
	}
	else if(numnode == 27)
	{
		// CHECK NODE NUMBERING
		nodeids[0] = NodeIds()[0];
	  	nodeids[1] = NodeIds()[1];
	  	nodeids[2] = NodeIds()[2];
	  	nodes[0] = Nodes()[0];
	  	nodes[1] = Nodes()[1];
	  	nodes[2] = Nodes()[1];
	  	lines_[0] =
	    	rcp(new DRT::Elements::XFluid3Line(0,Owner(),3,nodeids,nodes,NULL,this,0));
	  	lineptrs_[0] = lines_[0].get();
	
	  	nodeids[0] = NodeIds()[1];
	  	nodeids[1] = NodeIds()[2];
	  	nodeids[2] = NodeIds()[2];	  	
	  	nodes[0] = Nodes()[1];
	  	nodes[1] = Nodes()[2];
	  	nodes[2] = Nodes()[1];
	  	lines_[1] =
	    	rcp(new DRT::Elements::XFluid3Line(1,Owner(),3,nodeids,nodes,NULL,this,1));
	  	lineptrs_[1] = lines_[1].get();
	
	  	nodeids[0] = NodeIds()[2];
	  	nodeids[1] = NodeIds()[3];
	  	nodeids[2] = NodeIds()[2];	  	
	  	nodes[0] = Nodes()[2];
	  	nodes[1] = Nodes()[3];
	  	nodes[2] = Nodes()[1];
	  	lines_[2] =
	    	rcp(new DRT::Elements::XFluid3Line(2,Owner(),3,nodeids,nodes,NULL,this,2));
	  	lineptrs_[2] = lines_[2].get();
	
	  	nodeids[0] = NodeIds()[3];
	  	nodeids[1] = NodeIds()[0];
	  	nodeids[2] = NodeIds()[2];	  	
	  	nodes[0] = Nodes()[3];
	  	nodes[1] = Nodes()[0];
	  	nodes[2] = Nodes()[1];
	  	lines_[3] =
	    	rcp(new DRT::Elements::XFluid3Line(3,Owner(),3,nodeids,nodes,NULL,this,3));
	  	lineptrs_[3] = lines_[3].get();
	
	  	nodeids[0] = NodeIds()[0];
	  	nodeids[1] = NodeIds()[4];
	  	nodeids[2] = NodeIds()[2];	  	
	  	nodes[0] = Nodes()[0];
	  	nodes[1] = Nodes()[4];
	  	nodes[2] = Nodes()[1];
	  	lines_[4] =
	    	rcp(new DRT::Elements::XFluid3Line(4,Owner(),3,nodeids,nodes,NULL,this,4));
	  	lineptrs_[4] = lines_[4].get();
	
	  	nodeids[0] = NodeIds()[1];
	  	nodeids[1] = NodeIds()[5];
	  	nodeids[2] = NodeIds()[2];	  	
	  	nodes[0] = Nodes()[1];
	  	nodes[1] = Nodes()[5];
	  	nodes[2] = Nodes()[1];
	  	lines_[5] =
	    	rcp(new DRT::Elements::XFluid3Line(5,Owner(),3,nodeids,nodes,NULL,this,5));
	  	lineptrs_[5] = lines_[5].get();
	
	  	nodeids[0] = NodeIds()[2];
	  	nodeids[1] = NodeIds()[6];
	  	nodeids[2] = NodeIds()[2];	  	
	  	nodes[0] = Nodes()[2];
	  	nodes[1] = Nodes()[6];
	  	nodes[2] = Nodes()[1];
	  	lines_[6] =
	    	rcp(new DRT::Elements::XFluid3Line(6,Owner(),3,nodeids,nodes,NULL,this,6));
	  	lineptrs_[6] = lines_[6].get();
	
	  	nodeids[0] = NodeIds()[3];
	  	nodeids[1] = NodeIds()[7];
	  	nodeids[2] = NodeIds()[2];	  	
	  	nodes[0] = Nodes()[3];
	  	nodes[1] = Nodes()[7];
	  	nodes[2] = Nodes()[1];
	  	lines_[7] =
	    	rcp(new DRT::Elements::XFluid3Line(7,Owner(),3,nodeids,nodes,NULL,this,7));
	  	lineptrs_[7] = lines_[7].get();
	
	  	nodeids[0] = NodeIds()[4];
	  	nodeids[1] = NodeIds()[5];
	  	nodeids[2] = NodeIds()[2];	  	
	  	nodes[0] = Nodes()[4];
	  	nodes[1] = Nodes()[5];
	  	nodes[2] = Nodes()[1];
	  	lines_[8] =
	    	rcp(new DRT::Elements::XFluid3Line(8,Owner(),3,nodeids,nodes,NULL,this,8));
	  	lineptrs_[8] = lines_[8].get();
	
	  	nodeids[0] = NodeIds()[5];
	  	nodeids[1] = NodeIds()[6];
	  	nodeids[2] = NodeIds()[2];	  	
	  	nodes[0] = Nodes()[5];
	  	nodes[1] = Nodes()[6];
	  	nodes[2] = Nodes()[1];
	  	lines_[9] =
	    	rcp(new DRT::Elements::XFluid3Line(9,Owner(),3,nodeids,nodes,NULL,this,9));
	  	lineptrs_[9] = lines_[9].get();
	
	  	nodeids[0] = NodeIds()[6];
	  	nodeids[1] = NodeIds()[7];
	  	nodeids[2] = NodeIds()[2];	  	
	  	nodes[0] = Nodes()[6];
	  	nodes[1] = Nodes()[7];
	  	nodes[2] = Nodes()[1];
	  	lines_[10] =
	    	rcp(new DRT::Elements::XFluid3Line(10,Owner(),3,nodeids,nodes,NULL,this,10));
	  	lineptrs_[10] = lines_[10].get();
	
	  	nodeids[0] = NodeIds()[7];
	  	nodeids[1] = NodeIds()[4];
	  	nodeids[2] = NodeIds()[2];	  	
	  	nodes[0] = Nodes()[7];
	  	nodes[1] = Nodes()[4];
	  	nodes[2] = Nodes()[1];
	  	lines_[11] =
	    	rcp(new DRT::Elements::XFluid3Line(11,Owner(),3,nodeids,nodes,NULL,this,11));
	  	lineptrs_[11] = lines_[11].get();
	
	}
	else dserror("Number of nodes not supported");
  	return (DRT::Element**)(&(lineptrs_[0]));
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          g.bau 03/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::XFluid3::Surfaces()
{

    const int nsurf = NumSurface();
  	const int numnode = NumNode();
  	surfaces_.resize(nsurf);
  	surfaceptrs_.resize(nsurf);
  	int nodeids[100];
  	DRT::Node* nodes[100];

  	if (nsurf==4)
 	{
  		if (numnode==4)
    	{
	      nodeids[0] = NodeIds()[0];
	      nodeids[1] = NodeIds()[1];
	      nodeids[2] = NodeIds()[2];
	      nodes[0] = Nodes()[0];
	      nodes[1] = Nodes()[1];
	      nodes[2] = Nodes()[2];
	      surfaces_[0] =
	        rcp(new DRT::Elements::XFluid3Surface(0,Owner(),3,nodeids,nodes,this,0));
	      surfaceptrs_[0] = surfaces_[0].get();
	
	      nodeids[0] = NodeIds()[0];
	      nodeids[1] = NodeIds()[1];
	      nodeids[2] = NodeIds()[3];
	      nodes[0] = Nodes()[0];
	      nodes[1] = Nodes()[1];
	      nodes[2] = Nodes()[3];
	      surfaces_[1] =
	        rcp(new DRT::Elements::XFluid3Surface(1,Owner(),3,nodeids,nodes,this,1));
	      surfaceptrs_[1] = surfaces_[1].get();
	
	      nodeids[0] = NodeIds()[2];
	      nodeids[1] = NodeIds()[0];
	      nodeids[2] = NodeIds()[3];
	      nodes[0] = Nodes()[2];
	      nodes[1] = Nodes()[0];
	      nodes[2] = Nodes()[3];
	      surfaces_[2] =
	        rcp(new DRT::Elements::XFluid3Surface(2,Owner(),3,nodeids,nodes,this,2));
	      surfaceptrs_[2] = surfaces_[2].get();
	
	      nodeids[0] = NodeIds()[1];
	      nodeids[1] = NodeIds()[2];
	      nodeids[2] = NodeIds()[3];
	      nodes[0] = Nodes()[1];
	      nodes[1] = Nodes()[2];
	      nodes[2] = Nodes()[3];
	      surfaces_[3] =
	        rcp(new DRT::Elements::XFluid3Surface(3,Owner(),3,nodeids,nodes,this,3));
	      surfaceptrs_[3] = surfaces_[3].get();
    	}
    	else if (numnode==10)
    	{
      	dserror("TET10 surfaces not implemented.");
    	}
    	else dserror("Number of nodes not supported");
  	}
	else if (nsurf==6)
  	{
    	if (numnode==8)
    	{
	      nodeids[0] = NodeIds()[0];
	      nodeids[1] = NodeIds()[3];
	      nodeids[2] = NodeIds()[2];
	      nodeids[3] = NodeIds()[1];
	      nodes[0] = Nodes()[0];
	      nodes[1] = Nodes()[3];
	      nodes[2] = Nodes()[2];
	      nodes[3] = Nodes()[1];
	      surfaces_[0] =
	        rcp(new DRT::Elements::XFluid3Surface(0,Owner(),4,nodeids,nodes,this,0));
	      surfaceptrs_[0] = surfaces_[0].get();
	
	      nodeids[0] = NodeIds()[0];
	      nodeids[1] = NodeIds()[1];
	      nodeids[2] = NodeIds()[5];
	      nodeids[3] = NodeIds()[4];
	      nodes[0] = Nodes()[0];
	      nodes[1] = Nodes()[1];
	      nodes[2] = Nodes()[5];
	      nodes[3] = Nodes()[4];
	      surfaces_[1] =
	        rcp(new DRT::Elements::XFluid3Surface(1,Owner(),4,nodeids,nodes,this,1));
	      surfaceptrs_[1] = surfaces_[1].get();
	
	      nodeids[0] = NodeIds()[0];
	      nodeids[1] = NodeIds()[4];
	      nodeids[2] = NodeIds()[7];
	      nodeids[3] = NodeIds()[3];
	      nodes[0] = Nodes()[0];
	      nodes[1] = Nodes()[4];
	      nodes[2] = Nodes()[7];
	      nodes[3] = Nodes()[3];
	      surfaces_[2] =
	        rcp(new DRT::Elements::XFluid3Surface(2,Owner(),4,nodeids,nodes,this,2));
	      surfaceptrs_[2] = surfaces_[2].get();
	
	      nodeids[0] = NodeIds()[2];
	      nodeids[1] = NodeIds()[3];
	      nodeids[2] = NodeIds()[7];
	      nodeids[3] = NodeIds()[6];
	      nodes[0] = Nodes()[2];
	      nodes[1] = Nodes()[3];
	      nodes[2] = Nodes()[7];
	      nodes[3] = Nodes()[6];
	      surfaces_[3] =
	        rcp(new DRT::Elements::XFluid3Surface(3,Owner(),4,nodeids,nodes,this,3));
	      surfaceptrs_[3] = surfaces_[3].get();
	
	      nodeids[0] = NodeIds()[1];
	      nodeids[1] = NodeIds()[2];
	      nodeids[2] = NodeIds()[6];
	      nodeids[3] = NodeIds()[5];
	      nodes[0] = Nodes()[1];
	      nodes[1] = Nodes()[2];
	      nodes[2] = Nodes()[6];
	      nodes[3] = Nodes()[5];
	      surfaces_[4] =
	        rcp(new DRT::Elements::XFluid3Surface(4,Owner(),4,nodeids,nodes,this,4));
	      surfaceptrs_[4] = surfaces_[4].get();
	
	      nodeids[0] = NodeIds()[4];
	      nodeids[1] = NodeIds()[5];
	      nodeids[2] = NodeIds()[6];
	      nodeids[3] = NodeIds()[7];
	      nodes[0] = Nodes()[4];
	      nodes[1] = Nodes()[5];
	      nodes[2] = Nodes()[6];
	      nodes[3] = Nodes()[7];
	      surfaces_[5] =
	        rcp(new DRT::Elements::XFluid3Surface(5,Owner(),4,nodeids,nodes,this,5));
	      surfaceptrs_[5] = surfaces_[5].get();
    	}
    	else if (numnode==20)
    	{
      	dserror("hex20 surfaces not implemented.");
    	}
    	else if (numnode==27)
    	{
    		/*
			In this routine the shape functions and their natural first and second
			derivatives with respect to r/s/t are evaluated for H E X A H E D E R

   		Numbering of the nodes:

                           ^ t
                           |
                           |
                           |
                    8      |  15        7
                    o---------o---------o
                   /|                  /|
                  / |                 / |
                 /  |                /  |
              16o   |     o       14o   |
               /    o20       o    /    o19
              /     |             /     |
             /      |  13      6 /  |
          5 o---------o---------o   |
            |   o   |     o     |   o   |  ---------->
            |       o---------o-|-------o           s
            |      / 4       11 |      /3
            |     /             |     /
          17o    /    o         o18  /
            | 12o         o     |   o10
            |  /                |  /
            | /                 | /
            |/                  |/
            o---------o---------o
        1   /     9         2
           /
          /
         /
        r

   GiD:

                           ^ t
                           |
                           |
                           |
                    8      |  19        7
                    o---------o---------o
                   /|                  /|
                  / |                 / |
                 /  |                /  |
              20o   |   26o       18o   |
               /    o16     24o    /    o15
              /     |             /     |
             /      |  17      6 /  |
          5 o---------o---------o   23  |
            |   o   |   27o     |   o   |  ---------->
            |  25   o---------o-|-------o           s
            |      / 4       11 |      /3
            |     /             |     /
          13o    /  22o         o14  /
            | 12o         o     |   o10
            |  /         21     |  /
            | /                 | /
            |/                  |/
            o---------o---------o
        1   /     9         2
           /
          /
         /
        r



   	PROBLEM: GID has a different numbering of the element nodes than this one.
         So either the shape functions for hex20 and hex27 (see drawing)
			has to be adapted or during the input phase the numbering has to
        	be adapted to the shape functions.
        	This is all in progress and should be done for fluid3 and
        	brick1 the same way!!!!

   	There are no HEX27 Elements in brick1 so we just go ahead here and
   	use the GiD numbering for HEX27.
   	NUMBERING IS NOT YET CORRECTED 
   	*/
    		
			nodeids[0] = NodeIds()[0];
	      nodeids[1] = NodeIds()[3];
	      nodeids[2] = NodeIds()[2];
	      nodeids[3] = NodeIds()[1];
	      nodeids[4] = NodeIds()[11];
	      nodeids[5] = NodeIds()[10];
	      nodeids[6] = NodeIds()[9];
	      nodeids[7] = NodeIds()[8];
	      nodeids[8] = NodeIds()[20];
	      
	      nodes[0] = Nodes()[0];
	      nodes[1] = Nodes()[3];
	      nodes[2] = Nodes()[2];
	      nodes[3] = Nodes()[1];
	      nodes[4] = Nodes()[11];
	      nodes[5] = Nodes()[10];
	      nodes[6] = Nodes()[9];
	      nodes[7] = Nodes()[8];
	      nodes[8] = Nodes()[20];
	      surfaces_[0] =
	        rcp(new DRT::Elements::XFluid3Surface(0,Owner(),9,nodeids,nodes,this,0));
	      surfaceptrs_[0] = surfaces_[0].get();
	
	      nodeids[0] = NodeIds()[0];
	      nodeids[1] = NodeIds()[1];
	      nodeids[2] = NodeIds()[5];
	      nodeids[3] = NodeIds()[4];
	      nodeids[4] = NodeIds()[0];
	      nodeids[5] = NodeIds()[3];
	      nodeids[6] = NodeIds()[2];
	      nodeids[7] = NodeIds()[1];
	      nodeids[8] = NodeIds()[0];
	      nodes[0] = Nodes()[0];
	      nodes[1] = Nodes()[1];
	      nodes[2] = Nodes()[5];
	      nodes[3] = Nodes()[4];
	      nodes[4] = Nodes()[0];
	      nodes[5] = Nodes()[3];
	      nodes[6] = Nodes()[2];
	      nodes[7] = Nodes()[1];
	      nodes[8] = Nodes()[0];
	      surfaces_[1] =
	        rcp(new DRT::Elements::XFluid3Surface(1,Owner(),9,nodeids,nodes,this,1));
	      surfaceptrs_[1] = surfaces_[1].get();
	
	      nodeids[0] = NodeIds()[0];
	      nodeids[1] = NodeIds()[4];
	      nodeids[2] = NodeIds()[7];
	      nodeids[3] = NodeIds()[3];
	      nodeids[4] = NodeIds()[0];
	      nodeids[5] = NodeIds()[3];
	      nodeids[6] = NodeIds()[2];
	      nodeids[7] = NodeIds()[1];
	      nodeids[8] = NodeIds()[0];
	      nodes[0] = Nodes()[0];
	      nodes[1] = Nodes()[4];
	      nodes[2] = Nodes()[7];
	      nodes[3] = Nodes()[3];
	      nodes[4] = Nodes()[0];
	      nodes[5] = Nodes()[3];
	      nodes[6] = Nodes()[2];
	      nodes[7] = Nodes()[1];
	      nodes[8] = Nodes()[0];
	      surfaces_[2] =
	        rcp(new DRT::Elements::XFluid3Surface(2,Owner(),9,nodeids,nodes,this,2));
	      surfaceptrs_[2] = surfaces_[2].get();
	
	      nodeids[0] = NodeIds()[2];
	      nodeids[1] = NodeIds()[3];
	      nodeids[2] = NodeIds()[7];
	      nodeids[3] = NodeIds()[6];
	      nodes[0] = Nodes()[2];
	      nodes[1] = Nodes()[3];
	      nodes[2] = Nodes()[7];
	      nodes[3] = Nodes()[6];
	      nodes[4] = Nodes()[0];
	      nodes[5] = Nodes()[3];
	      nodes[6] = Nodes()[2];
	      nodes[7] = Nodes()[1];
	      nodes[8] = Nodes()[0];
	      surfaces_[3] =
	        rcp(new DRT::Elements::XFluid3Surface(3,Owner(),9,nodeids,nodes,this,3));
	      surfaceptrs_[3] = surfaces_[3].get();
	
	      nodeids[0] = NodeIds()[1];
	      nodeids[1] = NodeIds()[2];
	      nodeids[2] = NodeIds()[6];
	      nodeids[3] = NodeIds()[5];
	      nodeids[4] = NodeIds()[0];
	      nodeids[5] = NodeIds()[3];
	      nodeids[6] = NodeIds()[2];
	      nodeids[7] = NodeIds()[1];
	      nodeids[8] = NodeIds()[0];
	      nodes[0] = Nodes()[1];
	      nodes[1] = Nodes()[2];
	      nodes[2] = Nodes()[6];
	      nodes[3] = Nodes()[5];
	      nodes[4] = Nodes()[0];
	      nodes[5] = Nodes()[3];
	      nodes[6] = Nodes()[2];
	      nodes[7] = Nodes()[1];
	      nodes[8] = Nodes()[0];
	      surfaces_[4] =
	        rcp(new DRT::Elements::XFluid3Surface(4,Owner(),9,nodeids,nodes,this,4));
	      surfaceptrs_[4] = surfaces_[4].get();
	
	      nodeids[0] = NodeIds()[4];
	      nodeids[1] = NodeIds()[5];
	      nodeids[2] = NodeIds()[6];
	      nodeids[3] = NodeIds()[7];
	      nodeids[4] = NodeIds()[0];
	      nodeids[5] = NodeIds()[3];
	      nodeids[6] = NodeIds()[2];
	      nodeids[7] = NodeIds()[1];
	      nodeids[8] = NodeIds()[0];
	      nodes[0] = Nodes()[4];
	      nodes[1] = Nodes()[5];
	      nodes[2] = Nodes()[6];
	      nodes[3] = Nodes()[7];
	      nodes[4] = Nodes()[0];
	      nodes[5] = Nodes()[3];
	      nodes[6] = Nodes()[2];
	      nodes[7] = Nodes()[1];
	      nodes[8] = Nodes()[0];
	      surfaces_[5] =
	        rcp(new DRT::Elements::XFluid3Surface(5,Owner(),9,nodeids,nodes,this,5));
	      surfaceptrs_[5] = surfaces_[5].get();     
     	}
     	else dserror("Number of nodes not supported");
  	}
  	else dserror("Number of lines not supported");

  	return (DRT::Element**)(&(surfaceptrs_[0]));
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                g.bau 03/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::XFluid3::Volumes()
{
    volume_.resize(1);
    volume_[0] = this;     //points to XFluid3 element itself
    return &volume_[0];
}
 
   


  

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3Register::XFluid3Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
    return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3Register::XFluid3Register(
                               const DRT::Elements::XFluid3Register& old) :
ElementRegister(old)
{
    return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3Register* DRT::Elements::XFluid3Register::Clone() const
{
  return new DRT::Elements::XFluid3Register(*this);
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Register::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class ElementRegister
  vector<char> basedata(0);
  ElementRegister::Pack(basedata);
  AddtoPack(data,basedata);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Register::Unpack(const vector<char>& data)
{
    int position = 0;
    // extract type
    int type = 0;
    ExtractfromPack(position,data,type);
    if (type != UniqueParObjectId()) dserror("wrong instance type data");
    // base class ElementRegister
    vector<char> basedata(0);
    ExtractfromPack(position,data,basedata);
    ElementRegister::Unpack(basedata);

    if (position != (int)data.size())
        dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
    return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3Register::~XFluid3Register()
{
    return;
}


/*----------------------------------------------------------------------*
 |  print (public)                                           mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Register::Print(ostream& os) const
{
    os << "XFluid3Register ";
    ElementRegister::Print(os);
    return;
}




#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3_XFEM
