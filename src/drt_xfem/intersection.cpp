/*!----------------------------------------------------------------------
\file intersection.cpp

\brief  collection of intersection tools for the computation of the
        intersection of two arbitrary discretizations
        
    The class intersection handles the intersection computation of 
    Cartesian, linear and quadratic discretization. The discretiazation,
    which is intersected is refered to as xfem discretization and the
    discretization acting as a cutter is called cutter discretization.
    The intersection algorithm returns a list of quadratic integration cells
    for each intersected xfem element. 
     
    The methods are categorized as follows for clearity:
    MAIN    public method which has to be called 
            from outside to perform the intersection computation
    GM      general methods
    ICS     intersection candidate search
    CLI     contruction of the linearized interface
    CDT     contrained delaunay tetrahedralization
    RCI     recovery of the curved interface
    DB      debug methods
    
<pre>
Maintainer: Ursula Mayer
</pre>
 
*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "intersection.H"
#include "../drt_xfem/intersection_service.H"
#include "../drt_xfem/intersection_math.H"
#include "../drt_xfem/integrationcell.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_element.H"

#ifdef PARALLEL
#include "../drt_lib/drt_exporter.H"
#include <mpi.h>
#endif 

using namespace XFEM;
using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 |  MAIN:   computes the interface between the xfem          u.may 06/07|
 |          discretization and the cutter discretization.               |
 |          It returns a list of intersected xfem elements              |
 |          and their integrations cells.                               |
 *----------------------------------------------------------------------*/  
void Intersection::computeIntersection( const RCP<DRT::Discretization>  xfemdis,
                                        const RCP<DRT::Discretization>  cutterdis,  
                                        map< int, DomainIntCells >&   			domainintcells,
                                        map< int, BoundaryIntCells >&   		boundaryintcells,
                                        map< int, vector< DRT::Element* > >&    cutterElementMap,  ///< int is the xfem element global id
                                        map< int, RCP<DRT::Node> >&     cutterNodeMap      ///< int is the xfem element global id
                                        )
{
    
	
    bool xfemIntersection; 
    vector< DRT::Condition * >      		xfemConditions;
    vector< DRT::Element* > 				cutterElements;
    
    cutterElementMap.clear();
    cutterNodeMap.clear();

    countMissedPoints_ = 0;
    
    const double t_start = ds_cputime();
    //initializeIntersection(xfemdis);
    
#ifdef PARALLEL
	const int cmyrank = cutterdis->Comm().MyPID();
	const int xnumproc = xfemdis->Comm().NumProc();
	const int cnumproc = cutterdis->Comm().NumProc();
	map< int, set<int> > xfemCutterIdMap;
	vector< int > conditionEleCount;
	
	if(cnumproc != xnumproc)
		dserror("the number of processors for xfem and cutter discretizations have to equal each other");		
#endif

	//debugXFEMConditions(cutterdis);
    // obtain vector of pointers to all xfem conditions of all cutter discretizations
    cutterdis->GetCondition ("XFEMCoupling", xfemConditions);
        
    if(xfemConditions.size()==0)
        cout << "number of fsi xfem conditions = 0" << endl;
      
       
#ifdef PARALLEL
		adjustCutterElementNumbering(cutterdis, xfemConditions, conditionEleCount);
#endif         
   
    //  k < xfemdis->NumMyColElements()
    for(int k = 0; k < xfemdis->NumMyColElements(); ++k)
    {
        xfemIntersection = false;
        DRT::Element* xfemElement = xfemdis->lColElement(k);
        initializeXFEM(k, xfemElement);  
        
        const Epetra_SerialDenseMatrix xfemXAABB = computeFastXAABB(xfemElement);
            
        startPointList();
        
        //printf("size of xfem condition = %d\n", c);
        int condCounter = -1;  // was i
        for(vector<DRT::Condition*>::const_iterator conditer = xfemConditions.begin(); conditer!=xfemConditions.end(); ++conditer)
        {
            condCounter++;
            DRT::Condition* xfemCondition = *conditer;
            const map<int, RCP<DRT::Element > > geometryMap = xfemCondition->Geometry();
            map<int, RCP<DRT::Element > >::const_iterator iterGeo;   
            //if(geometryMap.size()==0)   printf("geometry does not obtain elements\n");
          	//printf("size of %d.geometry map = %d\n",i, geometryMap.size());
         
            for(iterGeo = geometryMap.begin(); iterGeo != geometryMap.end(); ++iterGeo )
            { 
                DRT::Element*  cutterElement = iterGeo->second.get();
                if(cutterElement == NULL) dserror("geometry does not obtain elements");
            	
                const Epetra_SerialDenseMatrix    cutterXAABB = computeFastXAABB(cutterElement);   
                  
                const bool intersected = intersectionOfXAABB(cutterXAABB, xfemXAABB);   
              
                if(intersected) 
                {
                	if(cutterElementMap.find(xfemElement->LID()) != cutterElementMap.end())
                	{
                    	cutterElementMap.find(xfemElement->LID())->second.push_back(cutterElement); 
                	}
                   	else
                   	{
                   		vector<DRT::Element *> cutterVector;
                   		cutterVector.push_back(cutterElement);
                   		cutterElementMap.insert(make_pair(xfemElement->LID(), cutterVector));
                   	}
#ifdef PARALLEL         
					int addToCutterId = 0;
					if(condCounter > 0) addToCutterId = conditionEleCount[condCounter];

                   	if(xfemCutterIdMap.find(xfemElement->LID()) != xfemCutterIdMap.end())
                   	{
                   		xfemCutterIdMap.find(xfemElement->LID())->second.insert(cutterElement->Id() + addToCutterId);  
                   	}
                   	else
                   	{
                   		set<int> cutterIds;
                   		cutterIds.insert(cutterElement->Id() + addToCutterId);
                   		xfemCutterIdMap.insert(make_pair(xfemElement->LID(), cutterIds));
                   	}
#endif                   	
                } 
            }// for-loop over all geometryMap
		}// for-loop over all xfemConditions        

#ifdef PARALLEL
	} // for loop over all xfem elements} 

	getCutterElementsInParallel(xfemConditions,  conditionEleCount, cutterElementMap,
	 							cutterNodeMap,  xfemCutterIdMap,  xfemdis, cutterdis);
	 										   
    for(int k = 0; k < xfemdis->NumMyColElements(); k++)
    {
        xfemIntersection = false;
        DRT::Element* xfemElement = xfemdis->lColElement(k);
        initializeXFEM(k, xfemElement);    
        startPointList();

#endif         
        
        DRT::Element** xfemElementSurfaces = xfemElement->Surfaces();
        DRT::Element** xfemElementLines = xfemElement->Lines();
        
        if(cutterElementMap.find(xfemElement->LID()) != cutterElementMap.end())
        {
            cutterElements = cutterElementMap.find(k)->second;
            //debugIntersection(xfemElement, cutterElements);
        }
        else
            cutterElements.resize(0);
          
        for(vector< DRT::Element* >::iterator i = cutterElements.begin(); i != cutterElements.end(); ++i )
        {
            DRT::Element* cutterElement = (*i);
            if(cutterElement == NULL) dserror("cutter element is null\n");
            DRT::Element*const* cutterElementLines = cutterElement->Lines();
            DRT::Node*const* cutterElementNodes = cutterElement->Nodes();
            
            numInternalPoints_= 0; 
            numBoundaryPoints_ = 0;
            vector< InterfacePoint >  interfacePoints; 
           
//            initializeCutter(cutterElement); 

            // collect internal points
            for(int m=0; m<cutterElement->NumLine() ; m++)                    
                collectInternalPoints( xfemElement, cutterElement, cutterElementNodes[m],
                                        interfacePoints, k, m);
  
            // collect intersection points                                   
            for(int m=0; m<xfemElement->NumLine() ; m++) 
                if(collectIntersectionPoints(   cutterElement, xfemElementLines[m],
                                                interfacePoints, 0, m, false, xfemIntersection))                                         
                    storeIntersectedCutterElement(cutterElement); 
            
            for(int m=0; m<cutterElement->NumLine() ; m++)                                              
                for(int p=0; p<xfemElement->NumSurface() ; p++) 
                    if(collectIntersectionPoints(   xfemElementSurfaces[p], 
                                                    cutterElementLines[m], interfacePoints,
                                                    p, m, true, xfemIntersection)) 
                        storeIntersectedCutterElement(cutterElement); 
            
            
            // order interface points           
            if( interfacePoints.size() > 0)
            {
#ifdef QHULL
                computeConvexHull( xfemElement, cutterElement, interfacePoints);
#else
                dserror("Set QHULL flag to use XFEM intersections!!!");
#endif
            }    
           
            //interfacePoints.clear();  
        }// for-loop over all cutter elements
        
        if(xfemIntersection)
        {         
            //debugTetgenDataStructure(xfemElement);
            computeCDT(xfemElement, cutterElements[0], domainintcells, boundaryintcells);
        }
        
    }// for-loop over all  actdis->NumMyColElements()
   
    //debugDomainIntCells(domainintcells,2);
    const double t_end = ds_cputime()-t_start;
    cout << endl;
    if(countMissedPoints_ > 0)
    	cout << "Number of missed points during the recovery copy = " << countMissedPoints_ << endl;
    	
    cout << "Intersection computed sucessfully in " << t_end  <<  " secs";
#ifdef PARALLEL
    flush(cout);
	cout << " rank = " << cmyrank;
    flush(cout);
#endif
    cout << endl;
    cout << endl;
  
}

void getCutterElementPlusNodes(
        map< int, vector< DRT::Element* > >& cutterElementMap,
        map< int, RCP<DRT::Node> >&  cutterNodeMap
        )
{
    
}


/*----------------------------------------------------------------------*
 |  INIT:   initializes the intersection computation         u.may 12/07|
 |                                                                      |
 *----------------------------------------------------------------------*/
//void Intersection::initializeIntersection(
//    const RCP<DRT::Discretization>  xfemdis)
//{
//    xfemOldDistype_ = xfemdis->lColElement(0)->Shape();
//}



/*----------------------------------------------------------------------*
 |  INIT:   initializes the private members of the           u.may 08/07|
 |          current xfem element                                        |
 *----------------------------------------------------------------------*/
void Intersection::initializeXFEM(
    const int             xfemId,
    const DRT::Element*   xfemElement)
{ 
    //xfemDistype_ = xfemElement->Shape();
    const DRT::Element::DiscretizationType xfemDistype = xfemElement->Shape();
    
    // TODO: check for paralleilzation
    //if(xfemOldDistype_ != xfemDistype_ || xfemId == 0) 
    //{   
        numXFEMSurfaces_ = xfemElement->NumSurface(); 
//        numXFEMLines_ = xfemElement->NumLine(); 
       
        numXFEMCornerNodes_  = getNumberOfElementCornerNodes(xfemDistype);
        
        eleLinesSurfaces_ = getEleNodeNumbering_lines_surfaces(xfemDistype);
        eleNodesSurfaces_ = getEleNodeNumbering_nodes_surfaces(xfemDistype);
        eleNumberingSurfaces_ = getEleNodeNumberingSurfaces(xfemDistype);
        eleRefCoordinates_ = getEleNodeNumbering_nodes_reference(xfemDistype);
    //}
    
    //xfemOldDistype_ = xfemElement->Shape();
    
    pointList_.clear();
    segmentList_.clear();                 
    surfacePointList_.clear();                
    triangleList_.clear();
    
    segmentList_.resize(numXFEMSurfaces_);
    surfacePointList_.resize(numXFEMSurfaces_);    
    
    intersectingCutterElements_.clear();
    faceMarker_.clear();
}



/*----------------------------------------------------------------------*
 |  INIT:   initializes the private members of the           u.may 08/07|
 |          current cutter element                                      |
 *----------------------------------------------------------------------*/
//void Intersection::initializeCutter(  
//    const DRT::Element* cutterElement)
//{  
//    cutterDistype_ = cutterElement->Shape();
//}


#ifdef PARALLEL
/*----------------------------------------------------------------------*
 |  PARALLEL:   collects all cutter elements on              u.may 10/07|
 |              each processor   										|
 *----------------------------------------------------------------------*/
void Intersection::adjustCutterElementNumbering( 
	const RCP<DRT::Discretization>&  	cutterdis,  
	vector< DRT::Condition * >&   				xfemConditions,
	vector<int>& 								conditionEleCount) const
{
	vector<int> countSend((int) xfemConditions.size());
	DRT::Exporter exporter(cutterdis->Comm());
			
	for(unsigned int i=0; i<xfemConditions.size(); i++)	
   	    countSend[i] = (int) xfemConditions[i]->Geometry().size();

	exporter.Allreduce(countSend, conditionEleCount, MPI_SUM );
	
	for(unsigned int i=1; i<xfemConditions.size(); i++)	
	   conditionEleCount[i] += conditionEleCount[i-1];
}



/*----------------------------------------------------------------------*
 |  PARALLEL:   pack cutter elements and their nodes         u.may 10/07|
 |              on each processor                                       |
 *----------------------------------------------------------------------*/
void Intersection::packData(
    const RCP<DRT::Discretization>      cutterdis,
    vector<int>&                                conditionSend, 
    vector<int>&                                lengthSend, 
    int&                                        nodeSetSizeSend,  
    vector<int>&                                nodeVectorSend,
    vector<char>&                               cutterDataSend ) const
{
    
    set<int> nodeSet;
    vector< DRT::Condition * > xfemConditions;
    // obtain vector of pointers to all xfem conditions of all cutter discretizations
    cutterdis->GetCondition ("XFEMCoupling", xfemConditions);
    
    lengthSend[0] = 0;
    // pack data on all processors
    for(unsigned int i=0; i<xfemConditions.size(); i++)
    {
        map< int, RCP<DRT::Element > >  geometryMap = xfemConditions[i]->Geometry();
        map< int, RCP<DRT::Element > >::iterator iterGeo;   
        //if(geometryMap.size()==0)   printf("geometry does not obtain elements\n");
        //cout << "size of  " << i << " geometry map = " << geometryMap.size() << " rank = " << cmyrank << endl;
        
        conditionSend.push_back((int) geometryMap.size());
        if(i>0) conditionSend[i] += conditionSend[i-1];
        // pack all cutter elements and create nodeSet
        for(iterGeo = geometryMap.begin(); iterGeo != geometryMap.end(); iterGeo++ )
        {   
            for(int inode = 0; inode < iterGeo->second.get()->NumNode(); inode++)
            {
                int nodeId = iterGeo->second.get()->Nodes()[inode]->Id();
                if(nodeSet.find(nodeId) == nodeSet.end())
                    nodeSet.insert(nodeId);     
            }
            vector<char> data;
            iterGeo->second.get()->Pack(data);
            DRT::ParObject::AddtoPack(cutterDataSend,data);
        } // for loop over geometry
    }// for loop over xfemConditions
    lengthSend[0] = cutterDataSend.size();
    nodeSetSizeSend = (int) nodeSet.size();
    
    
    // pack nodes
    if(lengthSend[0]>0)
    {
        set<int>::iterator nodeIter;   
        for(nodeIter = nodeSet.begin(); nodeIter != nodeSet.end(); nodeIter++)
        {   
            vector<char> data;
            cutterdis->gNode(*nodeIter)->Pack(data);
            DRT::ParObject::AddtoPack(cutterDataSend,data);
            nodeVectorSend.push_back(*nodeIter);
        }
        lengthSend[1] = cutterDataSend.size()-lengthSend[0];
    }
    else
        lengthSend[1] = 0;
}           



/*----------------------------------------------------------------------*
 |  PARALLEL:   unpacks the nodal data                       u.may 10/07|
 *----------------------------------------------------------------------*/
void Intersection::unpackNodes(
    int                             index, 
    vector<char>&                   cutterDataRecv, 
    vector<int>&                    nodeVectorRecv, 
    map< int, RCP<DRT::Node> >&     nodeMap ) const
{       
    int count = 0;
       
    while (index < (int) cutterDataRecv.size())
    {
        vector<char> data;
        DRT::ParObject::ExtractfromPack(index, cutterDataRecv, data);

        // allocate a node. Fill it with info from extracted element data
        DRT::ParObject* o = DRT::UTILS::Factory(data);

        // cast ParObject to node
        DRT::Node* actNode = dynamic_cast< DRT::Node* >(o);
        RCP<DRT::Node> nodeRCP = rcp(actNode);
    
        // store nodes  
        nodeMap.insert(make_pair(nodeVectorRecv[count++], nodeRCP));
    }
}



/*----------------------------------------------------------------------*
 |  PARALLEL:   collects all cutter elements on              u.may 10/07|
 |              each processor   										|
 *----------------------------------------------------------------------*/
void Intersection::getCutterElementsInParallel(	
	const vector< DRT::Condition * >&       xfemConditions,
	vector<int>&							conditionEleCount, 
	map< int, vector< DRT::Element* > >&   	cutterElementMap,
	map< int, RCP<DRT::Node> >&  	cutterNodeMap,
	map<int, set<int> >& 					xfemCutterIdMap,
	const RCP<DRT::Discretization>& xfemdis,
   	const RCP<DRT::Discretization>& cutterdis  )
{
	const int cmyrank = cutterdis->Comm().MyPID();
	const int cnumproc = cutterdis->Comm().NumProc();
	
	set<int> cutterNodeIdSet;
	set<int> cutterIdSet;
	
	vector<int> conditionSend;
	vector<int> conditionRecv;
	
	int length;
	vector<char> cutterDataSend;
	
	vector<int> lengthSend(2,0);
	vector<int> lengthRecv(2,0);

	int nodeSetSizeSend;
	vector<int> nodeSetSizeRecv(1,0);
	
	vector<int> nodeVectorSend;
	vector<int> nodeVectorRecv; 
	
	MPI_Request req, req1, req2, req3, req4;
	DRT::Exporter exporter(cutterdis->Comm());

	int dest = cmyrank+1;
	if(cmyrank == (cnumproc-1)) dest = 0;

	int source = cmyrank-1;
	if(cmyrank == 0) 	source = cnumproc-1;


    packData(cutterdis, conditionSend, lengthSend, nodeSetSizeSend, nodeVectorSend, cutterDataSend);

            
	// cnumproc - 1
	for(int num = 0; num < cnumproc - 1; num++)
	{		
		// send length of the datablock to be recieved			
		exporter.ISend(cmyrank, dest, &(lengthSend[0]) , 2, 0, req);
		
		if(lengthSend[0]>0)
		{		
			// send size of node Id set to be recieved	
			exporter.ISend(cmyrank, dest, &nodeSetSizeSend , 1, 1, req1);
			
			// send consitions size
			exporter.ISend(cmyrank, dest, &(conditionSend[0]) , (int) xfemConditions.size(), 2, req2);
		
			// send node Id set	
			exporter.ISend(cmyrank, dest, &(nodeVectorSend[0]) , nodeSetSizeSend, 3, req3);
		}
		
		// recieve length of incoming data
		length = (int) lengthRecv.size();
		exporter.Receive(source, 0, lengthRecv, length);
        exporter.Wait(req);   
        
        if(lengthRecv[0]>0)
		{
        	length = 1;
        	exporter.Receive(source, 1, nodeSetSizeRecv, length);
        	exporter.Wait(req1);  
        	
        	length = (int) xfemConditions.size();
        	exporter.Receive(source, 2, conditionRecv, length);
        	exporter.Wait(req2);  
			
        	exporter.Receive(source, 3, nodeVectorRecv, nodeSetSizeRecv[0]);
        	exporter.Wait(req3);   
		}
		else
		{
			nodeSetSizeRecv[0] = 0;
			if(!nodeVectorRecv.empty())
				nodeVectorRecv.clear();
			
			if(!conditionRecv.empty())
				conditionRecv.clear();	
		}
		
		
		if(lengthSend[0] > 0)
		{
			length = lengthSend[0] + lengthSend[1];
			exporter.ISend(cmyrank, dest, &(cutterDataSend[0]), length, 4, req4);
		}
		
        length = lengthRecv[0] + lengthRecv[1];
        vector<char> cutterDataRecv(length);
        
     
       	if(lengthRecv[0] > 0)
       	{
			int tag = 4;
        	exporter.ReceiveAny(source, tag, cutterDataRecv, length);
        	exporter.Wait(req4);  
            
            cutterDataSend = cutterDataRecv;
            
			int count = 0;
			int index = lengthRecv[0];
			map< int, RCP<DRT::Node> > nodeMap;
			
            unpackNodes(index, cutterDataRecv, nodeVectorRecv, nodeMap);
            
			index = 0;
			count = 0;

			while (index < lengthRecv[0])
			{				
  				// extract cutter element data form recieved data
  				vector<char> data;
  				DRT::ParObject::ExtractfromPack(index,cutterDataRecv,data);

  				// allocate an "empty cutter element". Fill it with info from
  				// extracted element data
  				DRT::ParObject* o = DRT::UTILS::Factory(data);

  				// cast ParObject to cutter element 
  				DRT::Element* actCutter = dynamic_cast<DRT::Element*>(o);
  				
  				// create nodal data for a particular element
				actCutter->BuildNodalPointers(nodeMap);
  				
  				int cutterIdAdd = 0;
  				for(unsigned int xf = 1; xf < xfemConditions.size(); xf++)
					if((count >= conditionRecv[xf-1]) && (count < conditionRecv[xf]) )
					{
						cutterIdAdd =  conditionEleCount[xf];
						break;
					}
				
                const int actCutterId = actCutter->Id() + cutterIdAdd;
  				count++;
  				
  				const Epetra_SerialDenseMatrix    cutterXAABB = computeFastXAABB(actCutter);
                
  				for(int k = 0; k < xfemdis->NumMyColElements(); k++)
    			{      
        			const DRT::Element* xfemElement = xfemdis->lColElement(k);
       				initializeXFEM(k, xfemElement);
  					const Epetra_SerialDenseMatrix    xfemXAABB = computeFastXAABB(xfemElement);
            		const bool intersected = intersectionOfXAABB(cutterXAABB, xfemXAABB);
            		//debugXAABBIntersection( cutterXAABB, xfemXAABB, actCutter, xfemElement, i, k);
                     
            		if(intersected)
        			{
        				if(cutterElementMap.find(xfemElement->LID()) != cutterElementMap.end())
        				{
                            set<int> currentSet =  xfemCutterIdMap.find(xfemElement->LID())->second; 
        					if(currentSet.find(actCutterId) == currentSet.end()) 
        					{
            					cutterElementMap.find(xfemElement->LID())->second.push_back(actCutter); 
            					xfemCutterIdMap.find(xfemElement->LID())->second.insert(actCutterId);   
        					}
        				}
           				else
           				{
           					vector<DRT::Element*> cutterVector;
           					cutterVector.push_back(actCutter);
           					cutterElementMap.insert(make_pair(xfemElement->LID(), cutterVector));
           					
           					set<int> cutterIds;
                   			cutterIds.insert(actCutterId);
                   			xfemCutterIdMap.insert(make_pair(xfemElement->LID(), cutterIds));
           				}
           				
           				if(cutterIdSet.find(actCutterId) == cutterIdSet.end())
           				{
           					cutterIdSet.insert(actCutterId);
           				
           					for(int inode = 0; inode < actCutter->NumNode(); inode++ )
           					{
           						const int nodeId =  actCutter->Nodes()[inode]->Id();
           						if(cutterNodeIdSet.find(nodeId) == cutterNodeIdSet.end())
           						{
           							cutterNodeIdSet.insert(nodeId);
           							cutterNodeMap.insert(make_pair(nodeId , nodeMap.find(nodeId)->second));	
           						}	
           					} // for loop over nodes
                            actCutter->BuildNodalPointers(cutterNodeMap);
           				} 
       	 			} // if intersected			
       	 		} // for - loop over all xfem elements         
			}
       	} // if recieve[0] > 0
       	else
            cutterDataSend.clear();
       
       	lengthSend = lengthRecv;
       	nodeSetSizeSend = nodeSetSizeRecv[0];
       	nodeVectorSend = nodeVectorRecv;
       	conditionSend = conditionRecv;
	}   // loop over all procs 
}

#endif




/*----------------------------------------------------------------------*
 |  CLI:    collects points that belong to the interface     u.may 06/07|
 |          and lie within an xfem element                              |
 *----------------------------------------------------------------------*/
bool Intersection::collectInternalPoints(   
    const DRT::Element*             xfemElement,
    DRT::Element*                   cutterElement,
    const DRT::Node*                node,
    std::vector< InterfacePoint >&  interfacePoints,
    const int                       elemId,
    const int                       nodeId)
{
    Epetra_SerialDenseVector xsi(3);
    Epetra_SerialDenseVector x(3);
       
    x[0] = node->X()[0];
    x[1] = node->X()[1];
    x[2] = node->X()[2];
    
    currentToElementCoordinates(xfemElement, x, xsi);
    const bool nodeWithinElement = checkPositionWithinElementParameterSpace(xsi, xfemElement->Shape());
    // debugNodeWithinElement(xfemElement,node, xsi, elemId ,nodeId, nodeWithinElement);  
    
    if(nodeWithinElement)
    {   
        InterfacePoint ip;
        //debugNodeWithinElement(xfemElement,node,xsi,elemId ,nodeId, nodeWithinElement);  
          
        numInternalPoints_++;
        
        // check if node lies on the boundary of the xfem element
        if(checkIfOnBoundary(xfemElement->Shape(), xsi, ip))  numBoundaryPoints_++;
                                   
        // intersection coordinates in the surface 
        // element element coordinate system
        getNodeCoordinates(nodeId, ip.coord, cutterElement->Shape());
                          
        interfacePoints.push_back(ip);
        
        storeIntersectedCutterElement(cutterElement); 
    }
      
    return nodeWithinElement;
}



/*----------------------------------------------------------------------*
 |  CLI:    checks if a node that lies within an element     u.may 06/07|
 |          lies on one of its surfaces or nodes                        |                                           
 *----------------------------------------------------------------------*/
bool Intersection::checkIfOnBoundary(
        const DRT::Element::DiscretizationType  xfemDistype,
        const Epetra_SerialDenseVector&         xsi,    
        InterfacePoint&                         ip
        ) const
{
    bool onSurface = false;
    const int count = getSurfaces(xsi, ip.surfaces, xfemDistype);
        
    // point lies on one surface
    if(count == 1)
    {
        onSurface = true;
        ip.nsurf = count;
        ip.pType = SURFACE;
    }
    // point lies on line, which has two neighbouring surfaces
    else if(count == 2)
    {
        onSurface = true;
        ip.nsurf = count;
        ip.pType = LINE;
    }
    // point lies on a node, which has three neighbouring surfaces
    else if(count == 3)
    {
        onSurface = true;
        ip.nsurf = count;
        ip.pType = NODE;
    }
    else
    {
        onSurface = false;
        ip.nsurf = 0;
        ip.pType = INTERNAL;
    }
    
    return onSurface;
}



/*----------------------------------------------------------------------*
 |  CLI:    collects all intersection points of a line and   u.may 06/07|
 |          and a surface                                               |
 *----------------------------------------------------------------------*/  
bool Intersection::collectIntersectionPoints(  
    const DRT::Element*             surfaceElement,
    const DRT::Element*             lineElement,
    std::vector<InterfacePoint>&    interfacePoints, 
    const int                       surfaceId,
    const int                       lineId,
    const bool                      lines,
    bool&                           xfemIntersection
    ) const
{
    Epetra_SerialDenseVector xsi(3);
    Epetra_SerialDenseVector upLimit(3);
    Epetra_SerialDenseVector loLimit(3);
  
    for(int i = 0; i < 3; i++)
    {
       upLimit[i]  =  1.0;  
       loLimit[i]  = -1.0;
       // extend for triangle
    }
    
    const bool intersected = computeCurveSurfaceIntersection(surfaceElement, lineElement, xsi, upLimit, loLimit);
                                        
    if(intersected) 
        addIntersectionPoint( surfaceElement, lineElement, xsi, upLimit, loLimit, 
                              interfacePoints, surfaceId, lineId, lines);
      
    
    // in this case a node of this line lies on the facet of the xfem element
    // but there is no intersection within the element                                          
    if(!((int) interfacePoints.size() == numBoundaryPoints_)) 
        xfemIntersection = true;
      
    return intersected;
}



/*----------------------------------------------------------------------*
 |  CLI:    computes the intersection between a              u.may 06/07|
 |          curve and a surface                    (CSI)                |
 *----------------------------------------------------------------------*/
bool Intersection::computeCurveSurfaceIntersection( 
    const DRT::Element*               surfaceElement,
    const DRT::Element*               lineElement,
    Epetra_SerialDenseVector&         xsi,
    const Epetra_SerialDenseVector&   upLimit,
    const Epetra_SerialDenseVector&   loLimit) const
{
    bool intersection = true;
    int iter = 0;
    const int maxiter = 30;
    double residual = 1.0;
    Epetra_SerialDenseMatrix A(3,3);
    Epetra_SerialDenseVector b(3);
    Epetra_SerialDenseVector dx(3);
    dx = xsi;
 
    updateRHSForCSI( b, xsi, surfaceElement, lineElement);
              
    while(residual > TOL14)
    {   
        updateAForCSI( A, xsi, surfaceElement, lineElement);
        
        if(!gaussElimination(A, b, dx, true, 3, 1))
        {    
            if(computeSingularCSI(xsi, lineElement, surfaceElement));
            {
                intersection = false;
                break;
            } 
            //for(int i = 0; i < 3; i++)  dx[i] = 0.0;
            dx.Scale(0.0);
            iter++;
            printf("SINGULAR\n");
        }
        
        //xsi = addTwoVectors(xsi,dx);
        xsi += dx;
        updateRHSForCSI( b, xsi, surfaceElement, lineElement);
        residual = b.Norm2(); 
        iter++;
        
        if(iter >= maxiter )
        {   
            intersection = false;
            break;
        }      
    } 
    
    if(intersection)
    {
        if( (xsi[0] > (upLimit[0]+TOL7)) || (xsi[1] > (upLimit[1]+TOL7)) || (xsi[2] > (upLimit[2]+TOL7))  || 
            (xsi[0] < (loLimit[0]-TOL7)) || (xsi[1] < (loLimit[1]-TOL7)) || (xsi[2] < (loLimit[2]-TOL7))) 
                intersection = false;
    }
           
    return intersection;
}



/*----------------------------------------------------------------------*
 |  CLI:    solves a singular linear system with the hel     u.may 12/07| 
 |          of a Singular Value Decomposition (SVD)                     |
 *----------------------------------------------------------------------*/
bool Intersection::computeSingularCSI(
    Epetra_SerialDenseVector&   xsi,
    const DRT::Element*         lineElement,
    const DRT::Element*         surfaceElement) const
{
    bool singular = false;
    int iter = 0;
    const int maxiter = 5;
    double residual = 1.0;
    Epetra_SerialDenseMatrix A(3,3);
    Epetra_SerialDenseVector b(3);
    Epetra_SerialDenseVector dx(3);
 
    updateRHSForCSI( b, xsi, surfaceElement, lineElement);
              
    while(residual > TOL14)
    {   
        updateAForCSI( A, xsi, surfaceElement, lineElement);
        
        if(solveLinearSystemWithSVD(A, b, dx, 3))
        {
            singular = false;
            //xsi = addTwoVectors(xsi,dx);
            xsi += dx;
            break;
        }
        
        //xsi = addTwoVectors(xsi,dx); 
        xsi += dx;
        updateRHSForCSI( b, xsi, surfaceElement, lineElement);
        residual = b.Norm2(); 
        iter++;
        
        if(iter >= maxiter )
        {
            singular = true;
            break;    
        }
    } 
    return singular;            
}



/*----------------------------------------------------------------------*
 |  CLI:    updates the systemmatrix for the computation     u.may 06/07| 
 |          of a curve-surface intersection (CSI)                       |
 *----------------------------------------------------------------------*/
void Intersection::updateAForCSI(  	Epetra_SerialDenseMatrix&         A,
                                    const Epetra_SerialDenseVector&   xsi,
                                    const DRT::Element*               surfaceElement,
                                    const DRT::Element*               lineElement) const
{	
	const int numNodesSurface = surfaceElement->NumNode();
   	const int numNodesLine = lineElement->NumNode();
   	
   	A.Scale(0.0);
 
    const BlitzMat surfaceDeriv1(shape_function_2D_deriv1(xsi[0], xsi[1], surfaceElement->Shape()));
    const DRT::Node*const* surfaceElementNodes = surfaceElement->Nodes();
    for(int inode=0; inode<numNodesSurface; inode++)
    {
        const double* x = surfaceElementNodes[inode] ->X();
        for(int isd=0; isd<3; isd++)
		{
			A[isd][0] += x[isd] * surfaceDeriv1(0,inode);
			A[isd][1] += x[isd] * surfaceDeriv1(1,inode);
		}
    }	
   
    const BlitzMat lineDeriv1(shape_function_1D_deriv1(xsi[2], lineElement->Shape()));
    const DRT::Node*const* lineElementNodes = lineElement->Nodes();
    for(int inode=0; inode<numNodesLine; inode++)
    {   
        const double* x = lineElementNodes[inode]->X();
        for(int isd=0; isd<3; isd++)
        {
            A[isd][2] -= x[isd] * lineDeriv1(0,inode);
        }   
    }
}

   

/*----------------------------------------------------------------------*
 |  CLI:    updates the right-hand-side for the              u.may 06/07|
 |          computation of                                              |
 |          a curve-surface intersection (CSI)                          |
 *----------------------------------------------------------------------*/
void Intersection::updateRHSForCSI( 
    Epetra_SerialDenseVector&         b,
    const Epetra_SerialDenseVector&   xsi,
    const DRT::Element*               surfaceElement,
    const DRT::Element*               lineElement) const
{
    const int numNodesSurface = surfaceElement->NumNode();
    const int numNodesLine = lineElement->NumNode();
 		
    b.Scale(0.0);
    
    const BlitzVec surfaceFunct(shape_function_2D( xsi[0], xsi[1], surfaceElement->Shape()));
    const DRT::Node*const* surfaceElementNodes = surfaceElement->Nodes();
    for(int i=0; i<numNodesSurface; i++)
    {
        const double* x = surfaceElementNodes[i]->X();
   	    for(int dim=0; dim<3; dim++)
			b[dim] -= x[dim] * surfaceFunct(i);
    }

    const BlitzVec lineFunct(shape_function_1D( xsi[2], lineElement->Shape()));
    const DRT::Node*const* lineElementNodes = lineElement->Nodes();
    for(int i=0; i<numNodesLine; i++)
    {
       const double* x = lineElementNodes[i]->X();
	   for(int dim=0; dim<3; dim++)   
			b[dim] += x[dim] * lineFunct(i);
    }
}



/*----------------------------------------------------------------------*
 |  CLI:    computes a new starting point for the            u.may 06/07| 
 |          Newton-method in order to find all intersection points      | 
 |          of a curve-surface intersection                             |
 *----------------------------------------------------------------------*/  
int Intersection::computeNewStartingPoint(
    const DRT::Element*                surfaceElement,
    const DRT::Element*                lineElement,
    const int                          surfaceId,
    const int                          lineId,
    const Epetra_SerialDenseVector&    xsiOld,
    const Epetra_SerialDenseVector&    upLimit,
    const Epetra_SerialDenseVector&    loLimit,
    std::vector<InterfacePoint>&       interfacePoints,
    const bool                         lines
    ) const
{	
    bool interval = true;
	int numInterfacePoints = 0;
	Epetra_SerialDenseVector xsi(3);
	
    //printf("xsi = %f   %f   %f\n", fabs(xsi[0]), fabs(xsi[1]), fabs(xsi[2]) ); 
    //printf("lolimit = %f   %f   %f\n", fabs(loLimit[0]), fabs(loLimit[1]), fabs(loLimit[2]) ); 
    //printf("uplimit = %f   %f   %f\n", fabs(upLimit[0]), fabs(upLimit[1]), fabs(upLimit[2]) ); 
    
    if(comparePoints(upLimit, loLimit))
        interval = false;
        
	for(int i = 0; i < 3; i++)
		xsi[i] = (double) (( upLimit[i] + loLimit[i] )/2.0);
         
	bool intersected = computeCurveSurfaceIntersection(surfaceElement, lineElement, xsi, upLimit, loLimit);
   
    
    if( comparePoints(xsi, xsiOld))     intersected = false;
       							 	
	if(intersected && interval)		
   		numInterfacePoints = addIntersectionPoint(	surfaceElement, lineElement,xsi, upLimit, loLimit, 
   													interfacePoints, surfaceId, lineId, lines);
   	
   	//printf("number of intersection points = %d\n", numInterfacePoints );
   	return numInterfacePoints;    						
}



/*----------------------------------------------------------------------*
 |  CLI:    adds an intersection point to the 			     u.may 07/07|
 |          list of interface points                       			    |
 *----------------------------------------------------------------------*/  
int Intersection::addIntersectionPoint(
    const DRT::Element*             surfaceElement,
    const DRT::Element*             lineElement,
    const Epetra_SerialDenseVector&	xsi,
    const Epetra_SerialDenseVector& upLimit,
    const Epetra_SerialDenseVector& loLimit,
    std::vector<InterfacePoint>& 	interfacePoints,
    const int                       surfaceId,
    const int                       lineId,							  
    const bool 					    lines
    ) const
{

	int numInterfacePoints = 0;
	
 	InterfacePoint ip;
    if(lines)
    {   
        ip.nsurf = 1;
        ip.surfaces[0] = surfaceId;
         
        getLineCoordinates(lineId, xsi[2], ip.coord, surfaceElement->Shape());
    }
    else
    {
        ip.nsurf = 2;
        ip.surfaces[0] = eleLinesSurfaces_[lineId][0];
        ip.surfaces[1] = eleLinesSurfaces_[lineId][1];
        ip.coord[0] = xsi[0]; 
        ip.coord[1] = xsi[1]; 
    }
    
    ip.coord[2] = 0.0; 
    ip.pType = INTERSECTION;
    
    vector<InterfacePoint>::iterator it;
    bool alreadyInList = false;
    for(it = interfacePoints.begin(); it != interfacePoints.end(); it++ )  
        if(comparePoints(ip.coord, it->coord,3))   
        {   
            //printf("alreadyinlist = true\n");
            alreadyInList = true;
            break;
            
        }
      
    if(!alreadyInList)      
    {
        vector< Epetra_SerialDenseVector >  upperLimits(8, Epetra_SerialDenseVector(3));
        vector< Epetra_SerialDenseVector >  lowerLimits(8, Epetra_SerialDenseVector(3));
        createNewLimits(xsi, upLimit, loLimit, upperLimits, lowerLimits);
        
        interfacePoints.push_back(ip);  
        numInterfacePoints++;
        
        // recursive call
        for(int i = 0; i < 8; i++)
            numInterfacePoints += computeNewStartingPoint(  
                                        surfaceElement, lineElement, surfaceId, lineId, xsi,
                                        upperLimits[i], lowerLimits[i], interfacePoints, lines);
        
    }
	return numInterfacePoints;
}


/*----------------------------------------------------------------------*
 |  CLI:    create new ranges for the recursive              u.may 07/07|
 |          computation of all intersection points                      |
 *----------------------------------------------------------------------*/  
void Intersection::createNewLimits(
    const Epetra_SerialDenseVector&         xsi, 
    const Epetra_SerialDenseVector&         upLimit,
    const Epetra_SerialDenseVector&         loLimit,
    vector< Epetra_SerialDenseVector >&     upperLimits, 
    vector< Epetra_SerialDenseVector >&     lowerLimits) const
{
    

/*          Surface:                                Line:
 *        (-1, 1)               (1,1)
 *          0_____________________1
 *          |          s          | 
 *          |         /\          |
 *          |          |          |                 4 ___________x__________ 5
 *          |          |          |              ( -1 )                    ( 1 )
 *          |          x ----> r  |
 *          |                     |
 *          |                     |
 *          |                     |
 *          2_____________________3
 *        (-1,-1)                (1,-1)
 */     

    // upper left corner of surface with lower part of line
    upperLimits[0][0] = xsi[0];         lowerLimits[0][0] = loLimit[0]; 
    upperLimits[0][1] = upLimit[1];     lowerLimits[0][1] = xsi[1];
    upperLimits[0][2] = xsi[2];         lowerLimits[0][2] = loLimit[2];
    
    // upper left corner of surface with upper part of line
    upperLimits[1][0] = xsi[0];         lowerLimits[1][0] = loLimit[0]; 
    upperLimits[1][1] = upLimit[1];     lowerLimits[1][1] = xsi[1];
    upperLimits[1][2] = upLimit[2];     lowerLimits[1][2] = xsi[2];
    
    // upper right corner of surface with lower part of line
    upperLimits[2][0] = upLimit[0];     lowerLimits[2][0] = xsi[0]; 
    upperLimits[2][1] = upLimit[1];     lowerLimits[2][1] = xsi[1];
    upperLimits[2][2] = xsi[2];         lowerLimits[2][2] = loLimit[2];
    
    // upper right corner of surface with upper part of line
    upperLimits[3][0] = upLimit[0];     lowerLimits[3][0] = xsi[0]; 
    upperLimits[3][1] = upLimit[1];     lowerLimits[3][1] = xsi[1];
    upperLimits[3][2] = upLimit[2];     lowerLimits[3][2] = xsi[2];

    // lower right corner of surface with lower part of line
    upperLimits[4][0] = upLimit[0];     lowerLimits[4][0] = xsi[0]; 
    upperLimits[4][1] = xsi[1];         lowerLimits[4][1] = loLimit[1];
    upperLimits[4][2] = xsi[2];         lowerLimits[4][2] = loLimit[2];
    
    // lower right corner of surface with upper part of line
    upperLimits[5][0] = upLimit[0];     lowerLimits[5][0] = xsi[0]; 
    upperLimits[5][1] = xsi[1];         lowerLimits[5][1] = loLimit[1];
    upperLimits[5][2] = upLimit[2];     lowerLimits[5][2] = xsi[2];
    
    // lower left corner of surface with lower part of line
    upperLimits[6][0] = xsi[0];         lowerLimits[6][0] = loLimit[0];
    upperLimits[6][1] = xsi[1];         lowerLimits[6][1] = loLimit[1];
    upperLimits[6][2] = xsi[2];         lowerLimits[6][2] = loLimit[2];
    
    // lower left corner of surface with upper part of line
    upperLimits[7][0] = xsi[0];         lowerLimits[7][0] = loLimit[0];
    upperLimits[7][1] = xsi[1];         lowerLimits[7][1] = loLimit[1];
    upperLimits[7][2] = upLimit[2];     lowerLimits[7][2] = xsi[2];
}



/*----------------------------------------------------------------------*
 |  ICS:    computes the convex hull of a set of             u.may 06/07|
 |          interface points and stores resulting points,               |
 |          segments and triangles for the use with Tetgen (CDT)        |
 *----------------------------------------------------------------------*/
#ifdef QHULL
void Intersection::computeConvexHull(   
        const DRT::Element*     xfemElement,
        const DRT::Element*     surfaceElement,
        vector<InterfacePoint>& interfacePoints)
{
    
    vector<int>                 positions;
    vector<double>              searchPoint(3,0);
    vector<double>              vertex(3,0);
    vector< vector<double> >    vertices;  
    InterfacePoint              midpoint;
    
           
    if(interfacePoints.size() > 2)  
    {
        midpoint = computeMidpoint(interfacePoints);
        // transform it into current coordinates
        Epetra_SerialDenseVector    curCoord(3);
        for(int j = 0; j < 2; j++)      curCoord[j]  = midpoint.coord[j];            
        elementToCurrentCoordinates(surfaceElement, curCoord);    
        currentToElementCoordinates(xfemElement, curCoord);    
        for(int j = 0; j < 3; j++)      midpoint.coord[j] = curCoord[j]; 
     
        // store coordinates in 
        // points has numInterfacePoints*dim-dimensional components
        // points[0] is the first coordinate of the first point
        // points[1] is the second coordinate of the first point
        // points[dim] is the first coordinate of the second point                     
        coordT* coordinates = (coordT *)malloc((2*interfacePoints.size())*sizeof(coordT));
        int fill = 0;
        for(vector<InterfacePoint>::iterator ipoint = interfacePoints.begin(); ipoint != interfacePoints.end(); ++ipoint)
        {
            for(int j = 0; j < 2; j++)
            {
                coordinates[fill++] = ipoint->coord[j]; 
               // printf("coord = %f\t", ipoint->coord[j]);
            }
            // printf("\n");
            // transform interface points into current coordinates
            Epetra_SerialDenseVector    curCoord(3);
            for(int j = 0; j < 2; j++)      
                curCoord[j]  = ipoint->coord[j];   
                          
            elementToCurrentCoordinates(surfaceElement, curCoord);  
            currentToElementCoordinates(xfemElement, curCoord);
                     
            for(int j = 0; j < 3; j++)         
                ipoint->coord[j] = curCoord[j];
                
        }     
      
        // compute convex hull - exitcode = 0 no error
        if (qh_new_qhull(2, interfacePoints.size(), coordinates, false, "qhull ", NULL, stderr)!=0) 
            dserror(" error in the computation of the convex hull (qhull error)"); 
                            
        if(((int) interfacePoints.size()) != qh num_vertices) 
            dserror("resulting surface is concave - convex hull does not include all points");  
       
        // copy vertices out of the facet list
        facetT* facet = qh facet_list;
        for(int i = 0; i< qh num_facets; i++)
        {
            for(int j = 0; j < 2; j++)
            {
                double* point  = SETelemt_(facet->vertices, j, vertexT)->point;
                for(int k = 0; k < 2; k++)  
                    vertex[k] = point[k];
                                
                Epetra_SerialDenseVector    curCoord(3);
                for(int m = 0; m < 2; m++)      curCoord[m]  = vertex[m];           
                // surface element coordinates to current coordinates       
                elementToCurrentCoordinates(surfaceElement, curCoord); 
                // current coordinates to xfem element element coordinates
                currentToElementCoordinates(xfemElement, curCoord);    
                for(int m = 0; m < 3; m++)      vertex[m] = curCoord[m];
                                                
                vertices.push_back(vertex);
            }                
            facet = facet->next;
        }
        
        // free memory and clear vector of interface points
        qh_freeqhull(!qh_ALL);
        int curlong, totlong;           // memory remaining after qh_memfreeshort 
        qh_memfreeshort (&curlong, &totlong);
        if (curlong || totlong) 
            printf("qhull internal warning (main): did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);

        free(coordinates);
                   
    }
    else if(interfacePoints.size() <= 2 && interfacePoints.size() > 0)
    {       
        for(vector<InterfacePoint>::iterator ipoint = interfacePoints.begin(); ipoint != interfacePoints.end(); ++ipoint)
        {
            // transform interface points into current coordinates
            Epetra_SerialDenseVector    curCoord(3);
            for(int j = 0; j < 2; j++)         
                curCoord[j]  = ipoint->coord[j];   
            
            // surface element coordinates to current coordinates       
            elementToCurrentCoordinates(surfaceElement, curCoord); 
            // current coordinates to xfem element element coordinates
            currentToElementCoordinates(xfemElement, curCoord);   
              
            for(int j = 0; j < 3; j++)      
            {   
                ipoint->coord[j] = curCoord[j];
                vertex[j] = curCoord[j];
            }
            vertices.push_back(vertex);
        }              
    }  
    else
        dserror("collection of interface points is empty");

    storePoint(vertices[0], interfacePoints, positions);
    vertices.erase(vertices.begin());
    
    if(interfacePoints.size() > 1)
    {
        // store points, segments and triangles for the computation of the
        // Constrained Delaunay Tetrahedralization with Tetgen  
        searchPoint = vertices[0];   
        storePoint(vertices[0], interfacePoints, positions );
        vertices.erase(vertices.begin());
    }
    
 
    while(vertices.size()>2)
    {                    
        findNextSegment(vertices, searchPoint);
        storePoint(searchPoint, interfacePoints, positions);      
    } 
    vertices.clear();
   
   
    storeSurfacePoints(interfacePoints);
        
    // cutter element lies on the surface of an xfem element
    if(numInternalPoints_ == numBoundaryPoints_ && numInternalPoints_ != 0)
    {
        if(numBoundaryPoints_ > 1)              
            storeSegments(positions);   
              
    }
    else
    {
        if(interfacePoints.size() > 1)
            storeSegments( positions );
       
        if(interfacePoints.size() > 2)
        {
            pointList_.push_back(midpoint);
            storeTriangles(positions);
        }
    }
    interfacePoints.clear();
    
}
#endif //QHULL
 


/*----------------------------------------------------------------------*
 |  ICS:    finds the next facet of a convex hull            u.may 06/07|
 |          and returns the point different form the searchpoint        |
 *----------------------------------------------------------------------*/  
void Intersection::findNextSegment(   
    vector< vector<double> >&   vertices, 
    vector<double>&             searchPoint) const
{     
    vector< vector<double> >::iterator it;
    bool pointfound = false;
    
    if(vertices.size()==0 || searchPoint.size()==0)
        dserror("one or both vectors are empty");   
    
    for(it = vertices.begin(); it != vertices.end(); it=it+2 )
    {      
        if(comparePoints(searchPoint, *it))
        {
            pointfound = true;
            searchPoint = *(it+1);              
            vertices.erase(it);
            vertices.erase(it); // remove it+ 1
            break; 
        }
               
        if(comparePoints(searchPoint, *(it+1)))
        {
            pointfound = true;
            searchPoint = *(it);                
            vertices.erase(it);
            vertices.erase(it); // remove it+ 1
            break;
        }
    }        
    if(!pointfound) dserror("no point found");   
}




/*----------------------------------------------------------------------*
 |  CDT:    computes the Constrained Delaunay                u.may 06/07|
 |          Tetrahedralization in 3D with help of Tetgen library        |
 |  for an intersected xfem element in element configuration          |
 |  TetGen provides the method                                          |
 |  void "tetrahedralize(char *switches, tetgenio* in, tetgenio* out)"  |
 |  as an interface for its use within other codes.                     |
 |  The char switches passes all the command line switches to TetGen.   |
 |  The most important command line switches include:                   |
 |      - d     detects intersections of PLC facets                     |
 |      - p     tetrahedralizes a PLC                                   |
 |      - q     quality mesh generation                                 |
 |      - nn    writes a list of boundary faces and their adjacent      |
 |              tetrahedra to the output tetgenio data structure        |
 |      - o2    resulting tetrahedra still have linear shape but        |
 |              have a 2nd order node distribution                      |
 |      - A     assigns region attributes                               |
 |      - Q     no terminal output except errors                        |
 |      - T     sets a tolerance                                        |
 |      - V     verbose: detailed information more terminal output      |
 |      - Y     prohibits Steiner point insertion on boundaries         |
 |              very help full for later visualization                  |
 |  The data structure tetgenio* in provides Tetgen with the input      |
 |  PLC and has to filled accordingly. tetgenio* out delivers the       |
 |  output information such as the resulting tetrahedral mesh.          |
 |  These two pointers must NOT be null at any time.                    |
 |  For further information please consult the TetGen manual            |                                         
 *----------------------------------------------------------------------*/  
void Intersection::computeCDT(  
        const DRT::Element*           			element,
        const DRT::Element*            			cutterElement,
        map< int, DomainIntCells >&	            domainintcells,
        map< int, BoundaryIntCells >&           boundaryintcells)
{
    int dim = 3;
    int nsegments = 0; 
    int nsurfPoints = 0;
    tetgenio in;
    tetgenio out;
    char switches[] = "pnno2QY";    //o2 Y
    tetgenio::facet *f;
    tetgenio::polygon *p;
    

    // allocate pointlist
    in.numberofpoints = pointList_.size();
    in.pointlist = new REAL[in.numberofpoints * dim];
       
    // fill point list
    int fill = 0;
    for(int i = 0; i <  in.numberofpoints; i++)
        for(int j = 0; j < dim; j++)  
            in.pointlist[fill++] = (REAL) pointList_[i].coord[j]; 
 
 
    in.pointmarkerlist = new int[in.numberofpoints];   
    for(int i = 0; i < numXFEMCornerNodes_; i++)
        in.pointmarkerlist[i] = 3;    // 3 : point not lying on the xfem element
        
    for(int i = numXFEMCornerNodes_; i < in.numberofpoints; i++)
        in.pointmarkerlist[i] = 2;    // 2 : point not lying on the xfem element

   
    if(triangleList_.size()>0)      in.numberoffacets = numXFEMSurfaces_ + triangleList_.size(); 
    else                            in.numberoffacets = numXFEMSurfaces_;   
      
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    
    
    // loop over all xfem element surfaces
    for(int i = 0; i < numXFEMSurfaces_; i++)
    {
        f = &in.facetlist[i];
        if(segmentList_[i].size() > 0)          nsegments = (int) (segmentList_[i].size()/2);
        else                                    nsegments = 0;
        if(surfacePointList_[i].size() > 0)     nsurfPoints = surfacePointList_[i].size();
        else                                    nsurfPoints = 0;
        f->numberofpolygons = 1 + nsegments + nsurfPoints; 
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];
        p->numberofvertices = 4;
        p->vertexlist = new int[p->numberofvertices];
        for(int j = 0; j < 4; j ++)
            p->vertexlist[j] = eleNumberingSurfaces_[i][j];
           
      
        int count = 0;
        for(int j = 1; j < 1 + nsegments; j ++)
        {
            if(segmentList_[i].size() > 0)
            {             
                p = &f->polygonlist[j];
                p->numberofvertices = 2;
                p->vertexlist = new int[p->numberofvertices];
            
                for(int k = 0; k < 2; k ++)
                {
                   p->vertexlist[k] = segmentList_[i][count]; 
                   in.pointmarkerlist[segmentList_[i][count]] = 3;  // 3: point lying on the xfem boundary
                   count++; 
                }  
            }
        } 
        
        count = 0;
        for(int j = 1 + nsegments; j < f->numberofpolygons; j++)
        {
            if(surfacePointList_[i].size() > 0)
            {             
                p = &f->polygonlist[j];
                p->numberofvertices = 1;
                p->vertexlist = new int[p->numberofvertices];
            
                p->vertexlist[0] = surfacePointList_[i][count];  
                in.pointmarkerlist[surfacePointList_[i][count]] = 3;  // 3: point lying on the xfem boundary
                count++;  
            }
        }    
    }
    
    // store triangles
    for(int i = numXFEMSurfaces_; i < in.numberoffacets; i++)
    {
        f = &in.facetlist[i];
        f->numberofpolygons = 1; 
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];
        p->numberofvertices = 3;
        p->vertexlist = new int[p->numberofvertices];
        for(int j = 0; j < 3; j ++)
            p->vertexlist[j] = triangleList_[i - element->NumSurface()][j];    
    }
      
        
    // set facetmarkers
    for(int i = 0; i < in.numberoffacets; i ++)
        in.facetmarkerlist[i] = faceMarker_[i] + facetMarkerOffset_;   
        
        
	// specify regions
/*	bool regions = true;	
	if(regions)
	{
        //double regionCoordinates[6];
		//computeRegionCoordinates(element,cutterElement,regionCoordinates);
	    fill = 0;
	    //int read = 0;
	    in.numberofregions = 2;
	    in.regionlist = new REAL[in.numberofregions*5];
	  	//for(int i = 0; i < in.numberofregions; i++)
	    //{
	    	// store coordinates 
	    	//for(int j = 0; j < 3; j++)
	    	//	in.regionlist[fill++] = regionCoordinates[read++];
            in.regionlist[fill++] = -0.9;
            in.regionlist[fill++] = -0.9;
            in.regionlist[fill++] = 0.9;
	    		
	    	// store regional attribute (switch A)   i=0 cutter i=1 fluid
	    	in.regionlist[fill++] = 1;
	    	
	    	// store volume constraint (switch a)   
	    	in.regionlist[fill++] = 0.0;
            
            in.regionlist[fill++] = -1.0;
            in.regionlist[fill++] = -1.0;
            in.regionlist[fill++] = -1.0;
                
            // store regional attribute (switch A)   i=0 cutter i=1 fluid
            in.regionlist[fill++] = 0;
            
            // store volume constraint (switch a)   
            in.regionlist[fill++] = 0.0;
		//} 
	}
	*/
    
    //in.save_nodes("tetin");   
    //in.save_poly("tetin");   
    //  Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
    //  do quality mesh generation (q) with a specified quality bound
    //  (1.414), and apply a maximum volume constraint (a0.1)
    tetrahedralize(switches, &in, &out); 
 
    
    //Debug
//    vector<int> elementIds;
//    for(int i = 388; i<389; i++)
//        elementIds.push_back(i);
    
    //debugTetgenOutput(in, out, element, elementIds);
    //printTetViewOutputPLC( element, element->Id(), in);
    
    // recover curved interface for higher order meshes
    const bool curvedInterface = true;
    if(curvedInterface)
        recoverCurvedInterface(element, boundaryintcells, out);
 
    //if(element->Id()==388)
    //	debugFaceMarker(element->Id(), out);
    	
    //printTetViewOutput(element->Id(), out);
    
    // store domain integration cells
    addCellsToDomainIntCellsMap(element, domainintcells, out);
}



/*----------------------------------------------------------------------*
 |  CDT:    fills the point list with the corner points      u.may 06/07|
 |          in element coordinates of the xfem element                  |
 *----------------------------------------------------------------------*/  
void Intersection::startPointList(
    )
{
    InterfacePoint ip;
        
    pointList_.clear();
    	
    for(int i = 0; i < numXFEMCornerNodes_; i++)
    {
        ip.nsurf = 3;
        
        // change for other element types
        for(int j = 0; j < 3; j++) 
        { 
           ip.coord[j]      =   eleRefCoordinates_[i][j]; 
           ip.surfaces[j]   =   eleNodesSurfaces_[i][j]; 
           ip.pType         =   NODE;
        }
        pointList_.push_back(ip);           
    }
    
    for(int i = 0; i < numXFEMSurfaces_; i++)
        faceMarker_.push_back(-1);   
         
}



/*----------------------------------------------------------------------*
 |  CDT:    stores a point within a list of points           u.may 06/07|
 |          which is to be copy to the tetgen data structure            |
 |          for the computation of the                                  |
 |          Constrained Delauney Triangulation                          |
 *----------------------------------------------------------------------*/  
void Intersection::storePoint(  
    const vector<double>&       point, 
    vector<InterfacePoint>&     interfacePoints, 
    vector<int>&                positions)
{
    bool alreadyInList = false; 

    for(vector<InterfacePoint>::const_iterator ipoint = interfacePoints.begin(); ipoint != interfacePoints.end(); ++ipoint )
    {   
        if(comparePoints(point, ipoint->coord))
        {         
            alreadyInList = false;
            int count = 0;
            for(vector<InterfacePoint>::const_iterator it = pointList_.begin(); it != pointList_.end(); ++it )  
            {
                if(comparePoints(point, it->coord))   
                {   
                    alreadyInList = true;
                    break;
                }
                count++;
            }  
            
            if(!alreadyInList) 
            {               
                pointList_.push_back(*ipoint);
                positions.push_back(pointList_.size()-1);
            }
            else
            {
                positions.push_back(count); 
            }
            break;
        } 
    }
}



/*----------------------------------------------------------------------*
 |  CDT:    computes the midpoint of a collection of         u.may 06/07| 
 |          InterfacePoints                                             |
 *----------------------------------------------------------------------*/  
InterfacePoint Intersection::computeMidpoint(
    const vector<InterfacePoint>& interfacePoints
    ) const
{
    
     int n = interfacePoints.size();
     InterfacePoint ip;
     
     ip.nsurf = 0;
     
     for(int i = 0; i < 3 ; i++)
        ip.coord[i] = 0.0;
     
     for(int i = 0; i < n ; i++)
        for(int j = 0; j < 3 ; j++)
         ip.coord[j] += interfacePoints[i].coord[j];
     
     for(int i = 0; i < 3 ; i++)
       ip.coord[i] = ip.coord[i]/(double)n;
     
     return ip;
}



/*----------------------------------------------------------------------*
 |  CDT:    stores a single point lying on a surface of an   u.may 06/07|
 |          xfem element, if no segments are lying in                   |
 |          that surface                                                |
 *----------------------------------------------------------------------*/  
void Intersection::storeSurfacePoints(  
    vector<InterfacePoint>&     interfacePoints)
{   

    for(unsigned int i = 0; i < interfacePoints.size(); i++)
    {
        bool singlePoint = true;
        if(interfacePoints[i].pType == SURFACE || interfacePoints[i].pType == LINE)
        {    
            for(unsigned int j = 0; j < interfacePoints.size(); j++)
            {
                if((interfacePoints[j].pType != INTERNAL) && (i != j))
                {
                    for(int k = 0; k < interfacePoints[i].nsurf; k++)
                    {
                        for(int l = 0; l < interfacePoints[j].nsurf; l++)
                        {
                            const int surf1 = interfacePoints[i].surfaces[k];
                            const int surf2 = interfacePoints[j].surfaces[l];
                        
                            if(surf1 == surf2)
                            {
                                singlePoint = false;
                                break;
                            }
                        }
                        if(!singlePoint)
                            break;  
                    }  
                }
                if(!singlePoint)
                        break;  
            }
        }
        else 
            singlePoint = false;
        
        
        if(singlePoint)
        {
            for(unsigned int jj = numXFEMCornerNodes_; jj < pointList_.size(); jj++ )
            {
                if(comparePoints(interfacePoints[i].coord, pointList_[jj].coord, 3)) 
                {
                    bool alreadyInList = false;
                    for(int kk = 0; kk < numXFEMSurfaces_; kk++)
                        for(unsigned int ll = 0; ll < surfacePointList_[kk].size(); ll++)
                            if(surfacePointList_[kk][ll] == (int) jj)
                            {
                                alreadyInList = true;
                                break;
                            }
                    
                    if(!alreadyInList)
                        surfacePointList_[interfacePoints[i].surfaces[0]].push_back(jj);   
                        
                    break;
                }
            }
        }
    } 
}



/*----------------------------------------------------------------------*
 |  CDT:    stores a segment within a list of segments       u.may 06/07|
 |          which is to be copied to the tetgen data structure          |
 |          for the computation of the                                  |
 |          Constrained Delauney Triangulation                          |
 *----------------------------------------------------------------------*/  
void Intersection::storeSegments(   
    const vector<int>&              positions)
{
    
    for(unsigned int i = 0; i < positions.size(); i++ )
    {
        const int pos1 = positions[i];
        int pos2 = 0;
        if(pos1 ==  positions[positions.size()-1])
            pos2 = positions[0]; 
        else
            pos2 = positions[i+1];
     
     
        for(int j = 0; j < pointList_[pos1].nsurf; j++ )  
            for(int k = 0; k < pointList_[pos2].nsurf; k++ ) 
            {
                const int surf1 = pointList_[pos1].surfaces[j];
                const int surf2 = pointList_[pos2].surfaces[k];
             
                if( (surf1 == surf2) ) 
                { 
                    bool alreadyInList = false;
                    
                    for(unsigned int is = 0 ; is < segmentList_[surf1].size() ; is = is + 2)
                    {
                        if( (segmentList_[surf1][is] == pos1  &&  segmentList_[surf1][is+1] == pos2)  ||
                            (segmentList_[surf1][is] == pos2  &&  segmentList_[surf1][is+1] == pos1) )
                        {
                            alreadyInList = true;
                            break;
                        }
                    }
                    
                    if(!alreadyInList)
                    {
                        segmentList_[surf1].push_back(pos1);
                        segmentList_[surf1].push_back(pos2);
                    }
                }
            }
    }
}
    

/*----------------------------------------------------------------------*
 |  CDT:    stores a triangle within a list of trianles      u.may 06/07|
 |          which is to be copy to the tetgen data structure            |
 |          for the computation of the                                  |
 |          Constrained Delauney Triangulation                          |
 *----------------------------------------------------------------------*/  
void Intersection::storeTriangles(  
    const vector<int>               positions)
{
    vector<int> triangle(3,0);
    
    for(unsigned int i = 0; i < positions.size()-1; i++ )
    {
        triangle[0] = positions[i];
        triangle[1] = positions[i+1];
        triangle[2] = pointList_.size()-1;
        
        triangleList_.push_back(triangle);
        faceMarker_.push_back(intersectingCutterElements_.size()-1);
    }
    
    triangle[0] = positions[positions.size()-1];
    triangle[1] = positions[0];
    triangle[2] = pointList_.size()-1;
        
    triangleList_.push_back(triangle);
    faceMarker_.push_back(intersectingCutterElements_.size()-1);
}


/*----------------------------------------------------------------------*
 |  RCI:    stores a pointer to each intersecting            u.may 08/07|
 |          cutter element  used for the recovery of curved interface   |     
 *----------------------------------------------------------------------*/  
void Intersection::storeIntersectedCutterElement(
        DRT::Element* surfaceElement)
{
    bool alreadyInList = false;
  
    vector<DRT::Element*>::const_iterator eleiter; 
    for(eleiter = intersectingCutterElements_.begin(); eleiter != intersectingCutterElements_.end(); ++eleiter)
        if(*eleiter == surfaceElement)
        {
            alreadyInList = true;
            break; 
        }
      
    if(!alreadyInList)
        intersectingCutterElements_.push_back(surfaceElement);  
}   



/*----------------------------------------------------------------------*
 |  RCI:    recovers the curved interface after the          u.may 08/07|
 |          Constrained Delaunay Tetrahedralization                     |
 *----------------------------------------------------------------------*/  
void Intersection::recoverCurvedInterface(
        const DRT::Element*             xfemElement,
        map< int, BoundaryIntCells >&   boundaryintcells,
        tetgenio&                       out
        )
{
    vector<int>                                     order(3,0);
    vector<int>                                     tetraCornerIndices(4,0); 
    vector< Epetra_SerialDenseVector >             	tetraCornerNodes(4, Epetra_SerialDenseVector(3));
    BoundaryIntCells                     			listBoundaryICPerElement;
    
    
    // list of point markers , if already visited = 1 , if not = 0
    vector<int> visitedPointIndexList(out.numberofpoints);      
    for(int i = 0; i<out.numberofpoints; ++i)
        visitedPointIndexList[i] = 0;
        
    // lifts all corner points into the curved interface
    liftAllSteinerPoints(xfemElement, out);
      
    for(int i=0; i<out.numberoftrifaces; i++)
    {
        // run over all faces not lying in on of the xfem element planes
        const int faceMarker = out.trifacemarkerlist[i] - facetMarkerOffset_;
        vector<vector<double> > domainCoord(6, std::vector<double>(3,0.0));
        vector<vector<double> > boundaryCoord(6, std::vector<double>(3,0.0));
                	
        if(faceMarker > -1)
        {        
            const int tetIndex = out.adjtetlist[i*2];
            //printf("tetIndex = %d\n", tetIndex); 

            getTetrahedronInformation(tetIndex, i, tetraCornerIndices, order, out );
            getTetrahedronNodes(tetraCornerNodes, tetraCornerIndices, xfemElement, out);
            
            // run over each triface
            for(int index1 = 0; index1 < 3 ;index1++)
            {                   
                int index2 = index1+1;
                if(index2 > 2) index2 = 0;
                
                //printf("edgeIndex1 = %d\n", out.tetrahedronlist[tetIndex*out.numberofcorners+order[index1]]);
                //printf("edgeIndex2 = %d\n", out.tetrahedronlist[tetIndex*out.numberofcorners+order[index2]]);
                
                const int localHigherOrderIndex = getHigherOrderIndex(order[index1], order[index2], DRT::Element::tet10); 
                //printf("localHOindex = %d\n", localHigherOrderIndex);
                const int globalHigherOrderIndex = out.tetrahedronlist[tetIndex*out.numberofcorners+localHigherOrderIndex];
                //printf("globalHOindex = %d\n", globalHigherOrderIndex);
                if(visitedPointIndexList[globalHigherOrderIndex]== 0)
                {     
                    visitedPointIndexList[globalHigherOrderIndex] = 1; 
       
                    computeHigherOrderPoint(    index1, index2, i, faceMarker, globalHigherOrderIndex, 
                                                tetraCornerIndices, tetraCornerNodes, xfemElement, out);                     
                }
                
                // store boundary integration cells
                addCellsToBoundaryIntCellsMap(	i, index1, globalHigherOrderIndex, faceMarker, 
                						domainCoord, boundaryCoord, xfemElement, out);
            }
            
            /*for(int ii=0; ii < 6; ii++)
            	for(int jj=0; jj < 3; jj++)
            		printf("boundary = %f\n",boundaryCoord[ii][jj]);
            	*/
            const int ele_gid = 1; //TODO: put in the correct number
            listBoundaryICPerElement.push_back(BoundaryIntCell(DRT::Element::tri6, ele_gid, domainCoord, boundaryCoord));        
         
        }
        
    }
    
    boundaryintcells.insert(make_pair(xfemElement->Id(),listBoundaryICPerElement));
    
    intersectingCutterElements_.clear();
}



/*----------------------------------------------------------------------*
 |  RCI:    checks if all tetrahedra corner points are lying u.may 09/07|
 |          in a surface element                                        |   
 |          if not corner points is recovered on the surface element    | 
 *----------------------------------------------------------------------*/
void Intersection::liftAllSteinerPoints(
        const DRT::Element*                             xfemElement,
        tetgenio&                                       out)
{
    Epetra_SerialDenseVector edgePoint(3);
    Epetra_SerialDenseVector oppositePoint(3);
    vector< vector <int> > adjacentFacesList;
    vector< vector <int> > adjacentFacemarkerList;
    
    locateSteinerPoints(adjacentFacesList, adjacentFacemarkerList, out);
    
    if(!adjacentFacesList.empty())
    {
        // run over all Steiner points
        for(unsigned int i=0; i<adjacentFacesList.size(); i++)
        {
            int lineIndex = -1, cutterIndex = -1;
            const int  caseSteiner = decideSteinerCase(  i, lineIndex, cutterIndex, adjacentFacesList, adjacentFacemarkerList,
                                                   edgePoint, oppositePoint, xfemElement, out);
            switch(caseSteiner)
            {
                case 1:
                {
                    liftSteinerPointOnSurface(i, adjacentFacesList, adjacentFacemarkerList, xfemElement, out);
                    break;
                }
                case 2:
                {
                    liftSteinerPointOnEdge( i, lineIndex, cutterIndex, edgePoint, oppositePoint, 
                                        adjacentFacesList, xfemElement, out);
                    break;
                }
                case 3:
                {
                    liftSteinerPointOnBoundary( i, adjacentFacesList, adjacentFacemarkerList, xfemElement, out);
                    break;
                }
         
                default:
                    dserror("case of lifting Steiner point does not exist");
            }
        }
    }
}
            


/*----------------------------------------------------------------------*
 |  RCI:    stores for each Steiner point its adjacent       u.may 11/07|
 |          faces and face markers                                      |   
 *----------------------------------------------------------------------*/
void Intersection::locateSteinerPoints(
    vector< vector <int> >&     adjacentFacesList,
    vector< vector <int> >&     adjacentFacemarkerList,
    tetgenio&                   out) const
{

    for(int i = 0; i < out.numberoftrifaces; i++)
    {   
        if( (out.trifacemarkerlist[i]-facetMarkerOffset_) > -1)
        {  
            for (int j = 0; j < 3; j++)
            {
                const int pointIndex = out.trifacelist[i*3+j];
                
                // check if point is a Steiner point
                if(out.pointmarkerlist[pointIndex] != 2 && out.pointmarkerlist[pointIndex] != 3 )
                {
                    bool alreadyInList = false;
                    const vector<int> pointIndices = getPointIndices(out, i, j);
                    
                    for(unsigned int k = 0; k < adjacentFacesList.size(); k++)
                        if(adjacentFacesList[k][0] == pointIndex)
                        {
                            alreadyInList = true;
                            
                            adjacentFacesList[k].push_back(pointIndices[0]);
                            adjacentFacesList[k].push_back(pointIndices[1]);
                            adjacentFacemarkerList[k].push_back( out.trifacemarkerlist[i]-facetMarkerOffset_ );
                            break;
                        }
                    
                    
                    if(!alreadyInList)
                    {
                        vector<int> adjacentFaces(3);
                        adjacentFaces[0] = pointIndex;
                        adjacentFaces[1] = pointIndices[0];
                        adjacentFaces[2] = pointIndices[1];
                    
                        adjacentFacesList.push_back(adjacentFaces);
                        
                        vector<int> adjacentFacemarkers(1);
                        adjacentFacemarkers[0] = out.trifacemarkerlist[i]-facetMarkerOffset_ ;
                        adjacentFacemarkerList.push_back(adjacentFacemarkers);
                    }
                }
            }       
        }
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    checks if the Steiner points lies within         u.may 11/07|
 |          the cutter element or on one of its edges                   |   
 *----------------------------------------------------------------------*/
int Intersection::decideSteinerCase( 
        const int                       steinerIndex,
        int&                            lineIndex, 
        int&                            cutterIndex,  
        const vector< vector <int> >&   adjacentFacesList,
        const vector< vector <int> >&   adjacentFacemarkerList,
        Epetra_SerialDenseVector&       edgePoint,
        Epetra_SerialDenseVector&       oppositePoint,
        const DRT::Element*             xfemElement, 
        const tetgenio&                 out) const
{

    bool                        normalSteiner = true;
    int                         pointIndex = adjacentFacesList[steinerIndex][0];
    Epetra_SerialDenseVector    x(3);
    Epetra_SerialDenseVector    xsi(3);
        
    for(int k=0; k<3; k++)
        x[k]   = out.pointlist[pointIndex*3 + k];
   
    InterfacePoint emptyIp;
    currentToElementCoordinates(xfemElement, x, xsi);
    if(checkIfOnBoundary(xfemElement-> Shape(), xsi, emptyIp))
        out.pointmarkerlist[pointIndex] = 3;    // on xfem boundary  
    else
        out.pointmarkerlist[pointIndex] = 2;    // not on xfem boundary
    
   
    for(unsigned int j = 0; j < adjacentFacemarkerList[steinerIndex].size(); j++ )
    {
        for(unsigned int k = 0; k < adjacentFacemarkerList[steinerIndex].size(); k++ )
        {
            if(adjacentFacemarkerList[steinerIndex][j] != adjacentFacemarkerList[steinerIndex][k])
            {
                //printf("a = %d and b = %d\n", adjacentFacemarkerList[steinerIndex][j], adjacentFacemarkerList[steinerIndex][k]);
                if(findCommonFaceEdge(  j, k, adjacentFacesList[steinerIndex], edgePoint, oppositePoint, out))
                {
                    if(!findCommonCutterLine(  adjacentFacemarkerList[steinerIndex][j], adjacentFacemarkerList[steinerIndex][k],
                                               lineIndex, cutterIndex))
                                               dserror("no common line element found\n");
                                                
                    normalSteiner = false;
                }
            }
            if(!normalSteiner)      break;
        }
        if(!normalSteiner)      break;
    }
    
    int caseSteiner = 0;
    
    if(normalSteiner)                                   caseSteiner = 1;
    else                                                caseSteiner = 2;
    if(out.pointmarkerlist[pointIndex] == 3)            caseSteiner = 3;
                                                                 
   
    
    return caseSteiner;
}



/*----------------------------------------------------------------------*
 |  RCI:    lifts Steiner points lying                       u.may 11/07|
 |          within a cutter element                                     |   
 *----------------------------------------------------------------------*/
void Intersection::liftSteinerPointOnSurface(
        const int                       steinerIndex,
        const vector< vector <int> >&   adjacentFacesList,
        const vector< vector <int> >&   adjacentFacemarkerList,
        const DRT::Element*             xfemElement, 
        tetgenio&                       out
        )
{
    // get Steiner point coordinates
    Epetra_SerialDenseVector    Steinerpoint(3);
    for(int j=0; j<3; ++j)
    {
        Steinerpoint[j] = out.pointlist[adjacentFacesList[steinerIndex][0]*3 + j];
    }
    elementToCurrentCoordinates(xfemElement, Steinerpoint);   
    
    Epetra_SerialDenseVector    averageNormal(3);
    averageNormal.Scale(0.0);
    vector<Epetra_SerialDenseVector> normals; 
    
    const int length = (int) ( ( (double) (adjacentFacesList[steinerIndex].size()-1))*0.5 );
    
    for(int j = 0; j < length; ++j)
    {
        const int pointIndex1 = adjacentFacesList[steinerIndex][1 + 2*j];
        const int pointIndex2 = adjacentFacesList[steinerIndex][1 + 2*j + 1];
        Epetra_SerialDenseVector    p1(3);
        Epetra_SerialDenseVector    p2(3);
     
        for(int k = 0; k < 3; ++k)
        {
            p1[k] =  out.pointlist[pointIndex1*3 + k];
            p2[k] =  out.pointlist[pointIndex2*3 + k];
        }
    
        elementToCurrentCoordinates(xfemElement, p1);   
        elementToCurrentCoordinates(xfemElement, p2);   
        const Epetra_SerialDenseVector n1 = subtractsTwoVectors(p1, Steinerpoint);
        const Epetra_SerialDenseVector n2 = subtractsTwoVectors(p2, Steinerpoint);
    
        Epetra_SerialDenseVector normal = computeCrossProduct( n1, n2);
        normalizeVector(normal);
    
//        for(int k=0; k<3; k++)
//            averageNormal[k] += normal[k];
        averageNormal += normal;
            
        normals.push_back(normal);
    }

    // compute average normal
//    for(int j=0; j<3; j++)
//        averageNormal[j] = averageNormal[j]/( (double)length);
    averageNormal.Scale(1.0 / ((double)length));
        
    vector<Epetra_SerialDenseVector> plane(4, Epetra_SerialDenseVector(3));                      
    plane[0] = addTwoVectors(Steinerpoint, averageNormal);               
    plane[1] = subtractsTwoVectors(Steinerpoint, averageNormal);
 
 
    const int faceMarker = adjacentFacemarkerList[steinerIndex][0];
    Epetra_SerialDenseVector xsi(3);
    
    const bool intersected = computeRecoveryNormal( xsi, plane, intersectingCutterElements_[faceMarker],false);
    if(intersected)
    {
        storeHigherOrderNode(   true, adjacentFacesList[steinerIndex][0], -1, xsi,
                                intersectingCutterElements_[faceMarker], xfemElement, out);
    }
    else
    {
        // loop over all individual normals
        bool intersected = false;
        vector<Epetra_SerialDenseVector>::const_iterator normalptr;
        for(normalptr = normals.begin(); normalptr != normals.end(); ++normalptr )
        {
            plane[0] = addTwoVectors(Steinerpoint, *normalptr);               
            plane[1] = subtractsTwoVectors(Steinerpoint, *normalptr);
            intersected = computeRecoveryNormal( xsi, plane, intersectingCutterElements_[faceMarker], false);
            if(intersected)
            {
                storeHigherOrderNode( true, adjacentFacesList[steinerIndex][0], -1, xsi,
                                      intersectingCutterElements_[faceMarker], xfemElement, out);
                break;
            }
        }
        if(!intersected)
        {
            countMissedPoints_++;
            printf("STEINER POINT NOT LIFTED\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"); 
        }  
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    lifts Steiner points lying                       u.may 11/07|
 |          on the edge of a  cutter element                            |   
 *----------------------------------------------------------------------*/
void Intersection::liftSteinerPointOnEdge(
    const int                         steinerIndex,
    int                         lineIndex,     // TODO:??????????????? why is this an argument to this function? 
    const int                         cutterIndex,   
    Epetra_SerialDenseVector&   edgePoint,
    Epetra_SerialDenseVector&   oppositePoint,
    const vector< vector <int> >&     adjacentFacesList,
    const DRT::Element*         xfemElement, 
    tetgenio&                   out)
{

    Epetra_SerialDenseVector    Steinerpoint(3);
    
     // get Steiner point coordinates 
    for(int j=0; j<3; j++)
        Steinerpoint[j] = out.pointlist[adjacentFacesList[steinerIndex][0]*3 + j];
    
    elementToCurrentCoordinates(xfemElement, Steinerpoint);   
    elementToCurrentCoordinates(xfemElement, edgePoint);   
    elementToCurrentCoordinates(xfemElement, oppositePoint);  
     
    const Epetra_SerialDenseVector r1 = subtractsTwoVectors( edgePoint, Steinerpoint);
    const Epetra_SerialDenseVector r2 = subtractsTwoVectors( oppositePoint, Steinerpoint);

    Epetra_SerialDenseVector n1 = computeCrossProduct( r1, r2);
    Epetra_SerialDenseVector n2 = computeCrossProduct( r1, n1);
    
    normalizeVector(n1);
    normalizeVector(n2);

    vector<Epetra_SerialDenseVector> plane(4, Epetra_SerialDenseVector(3));      
    plane[0] = addTwoVectors(Steinerpoint, n1);               
    plane[1] = subtractsTwoVectors(Steinerpoint, n1);
    plane[2] = addTwoVectors(plane[1], n2);
    plane[3] = addTwoVectors(plane[0], n2);
     
    Epetra_SerialDenseVector xsi(3);
    const bool intersected = computeRecoveryPlane( lineIndex, xsi, plane, intersectingCutterElements_[cutterIndex]);
    
    if(intersected)
    {   
        storeHigherOrderNode(   false, adjacentFacesList[steinerIndex][0], lineIndex, xsi,
                                intersectingCutterElements_[cutterIndex], xfemElement, out);
    }
    else
    {
        countMissedPoints_++;
        printf("STEINER POINT NOT LIFTED\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");   
    }
}


/*----------------------------------------------------------------------*
 |  RCI:    lifts Steiner points lying                       u.may 11/07|
 |          on the boundary of the xfem element                         |   
 *----------------------------------------------------------------------*/
void Intersection::liftSteinerPointOnBoundary( 
        const int                         steinerIndex, 
        const vector< vector <int> >&     adjacentFacesList, 
        const vector< vector <int> >&     adjacentFacemarkerList, 
        const DRT::Element*               xfemElement, 
        tetgenio&                         out
        )
{
    int edgeIndex = 0;
    int oppositeIndex = 0;
    int faceIndex = 0;
    int facemarkerIndex = 0;
   
    // find egde on boundary
    bool edgeFound = false;
    for(unsigned int i = 1; i < adjacentFacesList[steinerIndex].size(); i++ )
    {
        edgeIndex = adjacentFacesList[steinerIndex][i];
        if(out.pointmarkerlist[edgeIndex] == 3)
        {
            edgeFound = true;
            facemarkerIndex = (int) (i+1)/2;
            break;
        }
    }   
        
    faceIndex = adjacentFacemarkerList[steinerIndex][facemarkerIndex];
    // locate triangle on boundary
    bool oppositeFound = false;
    for(int i = 0; i < out.numberoftrifaces; i++)
    {   
        if( (out.trifacemarkerlist[i]-facetMarkerOffset_) == -1)
        {
            int countIndex = 0;
            for(int j = 0; j < 3; j++)
            {
                const int index = out.trifacelist[i*3+j];
                if(index == steinerIndex || index == edgeIndex)
                    countIndex++;
            }
            
            if(countIndex == 2)
            {
                for(int j = 0; j < 3; j++)
                {
                    const int index = out.trifacelist[i*3+j];
                    if(index != steinerIndex && index != edgeIndex)
                    {
                        oppositeIndex = index;
                        oppositeFound = true;
                        break;
                    }
                } 
            }
        }
        if(!oppositeFound)
            break;
    }
        
    // compute normal through Steiner point lying in the boundary triangle
    vector <Epetra_SerialDenseVector> plane(5, Epetra_SerialDenseVector(3)); 
    computeIntersectionNormal( adjacentFacesList[steinerIndex][0], edgeIndex, oppositeIndex, plane,
                               xfemElement, out);
               
    // compute intersection normal on boundary
    Epetra_SerialDenseVector xsi(3);
    bool intersected = computeRecoveryNormal(xsi, plane, intersectingCutterElements_[faceIndex], true);
    
    if(intersected)
    {   
        storeHigherOrderNode(   true, adjacentFacesList[steinerIndex][0], -1, xsi,
                                intersectingCutterElements_[faceIndex], xfemElement, out);
    }
    else
    {
        countMissedPoints_++;
        printf("STEINER POINT NOT LIFTED\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");   
    }

}



/*----------------------------------------------------------------------*
 |  RCI:    returns information of the tetrahedra            u.may 08/07|
 *----------------------------------------------------------------------*/  
void Intersection::getTetrahedronInformation( 
    const int           tetIndex,
    const int           faceIndex,  
    vector<int>&        tetraCornerIndices,
    vector<int>&        order,
    const tetgenio&     out) const
{            
    
    // store boundary face node indices
    for(int j=0; j<3; j++)
        tetraCornerIndices[j] = out.trifacelist[faceIndex*3+j];
   
    // store node index opposite to the boundary face of the tetrahedron
    for(int j=0; j<4; j++)
    {
        const int nodeIndex = out.tetrahedronlist[tetIndex*out.numberofcorners+j];
        if(nodeIndex != tetraCornerIndices[0] && nodeIndex != tetraCornerIndices[1] && nodeIndex != tetraCornerIndices[2])
        {
            tetraCornerIndices[3] = out.tetrahedronlist[tetIndex*out.numberofcorners+j];
            break;
        }
    }
    
    for(int j=0; j<4; j++)
    {
        const int nodeIndex = out.tetrahedronlist[tetIndex*out.numberofcorners+j];
        for(int k=0; k<3; k++)
            if(nodeIndex == tetraCornerIndices[k])
            {
                order[k] = j;
                break;
            }
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    collects the tetrahedron corner nodes            u.may 09/07|
 |          transforms them into current coordinates                    |
 |          of the xfem-element                                         |
 *----------------------------------------------------------------------*/  
void Intersection::getTetrahedronNodes(
        vector<Epetra_SerialDenseVector>&       tetraCornerNodes,
        const vector<int>&                      tetraCornerIndices,
        const DRT::Element*                     xfemElement,
        const tetgenio&                         out
        ) const
{

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 3; j++)
            tetraCornerNodes[i][j] = out.pointlist[tetraCornerIndices[i]*3+j];
        
        elementToCurrentCoordinates(xfemElement, tetraCornerNodes[i]);   
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    lifts the higher-order point of an edge of the   u.may 09/07|
 |          linearized interface onto the curved interface              |   
 *----------------------------------------------------------------------*/
void Intersection::computeHigherOrderPoint(    
        const int                                 index1, 
        const int                                 index2, 
        const int                                 faceIndex, 
        const int                                 faceMarker, 
        const int                                 globalHigherOrderIndex, 
        const vector<int>&                        tetraCornerIndices, 
        const vector<Epetra_SerialDenseVector>&   tetraCornerNodes, 
        const DRT::Element*                       xfemElement, 
        tetgenio&                           out
        )
{ 
                            
    bool                                    intersected             = false;   
    bool                                    intersectionNormal      = true; 
    int                                     lineIndex               = -1;
    int                                     adjacentFaceMarker      = -1;
    int                                     adjacentFaceIndex       = -1;  
    Epetra_SerialDenseVector                xsi(3);                
    vector < Epetra_SerialDenseVector >     plane(5, Epetra_SerialDenseVector(3)); 
                         
         
    findAdjacentFace(  tetraCornerIndices[index1], tetraCornerIndices[index2], 
                       faceMarker, adjacentFaceMarker, faceIndex, adjacentFaceIndex, out); 
                           
    
    // edge lies within the xfem element
    if(adjacentFaceMarker  > -1)
    {
        computeIntersectionNormal(  tetraCornerIndices[index1], tetraCornerIndices[index2], faceIndex, 
                                    adjacentFaceIndex, globalHigherOrderIndex, plane, xfemElement, out);   
                                   
        // higher order node lies within the cutter element
        if(adjacentFaceMarker == faceMarker)
        {                     
            intersected = computeRecoveryNormal(xsi, plane,intersectingCutterElements_[faceMarker], false);   
            intersectionNormal = true;
        }
        // higher-order point lies on one of the boundary lines of the cutter element
        else if(adjacentFaceMarker != faceMarker)
        {       
            int cutterIndex        = -1;
            findCommonCutterLine(faceMarker, adjacentFaceMarker, lineIndex, cutterIndex); 
            
            if(lineIndex != -1)
            {
                intersected = computeRecoveryPlane( lineIndex, xsi, plane, intersectingCutterElements_[cutterIndex]);
                intersectionNormal = false;
            }
        }
    }
    // edge lies on the surface of the xfem element
    else if(adjacentFaceMarker  == -1)
    {
        int oppositeIndex = findEdgeOppositeIndex(  tetraCornerIndices[index1], tetraCornerIndices[index2], 
                                                    adjacentFaceIndex, out);  
                                                                                  
        //printf("oppo = %d\n", oppositeIndex);
        
        computeIntersectionNormal(  true,  index1, index2, oppositeIndex, globalHigherOrderIndex,   
                                    tetraCornerIndices, tetraCornerNodes, plane, xfemElement, out);
                               
        intersected = computeRecoveryNormal(    xsi, plane, intersectingCutterElements_[faceMarker], 
                                                true);
        intersectionNormal = true;
                
        if(!intersected)
        {
            printf("REFERNCE DOMAIN");
            lineIndex = findIntersectingSurfaceEdge(    xfemElement, intersectingCutterElements_[faceMarker],
                                                        tetraCornerNodes[index1], tetraCornerNodes[index2]);
            if(lineIndex != -1)
            {                 
                intersected = computeRecoveryPlane( lineIndex, xsi, plane, intersectingCutterElements_[faceMarker]);
                intersectionNormal = false;
            }
        }
    }


    if(intersected)
    {   
        storeHigherOrderNode(   intersectionNormal, globalHigherOrderIndex, lineIndex,
                                xsi, intersectingCutterElements_[faceMarker], xfemElement, out);
    }
    else 
    {
        countMissedPoints_++;
        printf("NO INTERSECTION POINT FOUND!!!!! adjacentFaceMarker = %d\n", adjacentFaceMarker);
    }

}



/*----------------------------------------------------------------------*
 |  RCI:    returns the other two point indices belonging    u.may 09/07|
 |          to a triface that obtains a Steiner point                   |           
 *----------------------------------------------------------------------*/
vector<int> Intersection::getPointIndices(
        const tetgenio&   out, 
        const int         trifaceIndex, 
        const int         steinerPointIndex
        ) const
{
    
    int count = 0;
    vector<int> pointIndices(2);
    
    for(int i = 0; i < 3; i++)
        if(i != steinerPointIndex)
            pointIndices[count++] = out.trifacelist[trifaceIndex*3+i];
           
    return pointIndices;
}



/*----------------------------------------------------------------------*
 |  RCI:    computes the intersection between a              u.may 08/07|
 |          line  and a surface                    RCI                  |           
 *----------------------------------------------------------------------*/
bool Intersection::computeRecoveryNormal( 
    Epetra_SerialDenseVector&                   xsi,
    const vector<Epetra_SerialDenseVector>&     normal,
    const DRT::Element*                         cutterElement,
    const bool                                  onBoundary
    ) const
{
    bool                        intersection = true;
    int                         iter = 0;
    int                         countSingular = 0;
    const int                   maxiter = 50;
    double                      residual = 1.0;
    Epetra_SerialDenseMatrix    A(3,3);
    Epetra_SerialDenseVector    b(3);
    Epetra_SerialDenseVector    dx(3);
    
    xsi.Scale(0.0);
    updateRHSForRCINormal( b, xsi, normal, cutterElement, onBoundary);
                                
    while(residual > TOL14)
    {   
        updateAForRCINormal( A, xsi, normal, cutterElement, onBoundary);
         
        if(!solveLinearSystemWithSVD(A, b, dx, 3))
            countSingular++;
    
        if(countSingular > 5)
        {
            intersection = false;
            break;  
        }
        
        //xsi = addTwoVectors(xsi,dx);
        xsi += dx;
        //printf("dx0 = %20.16f\t, dx1 = %20.16f\t, dx2 = %20.16f\n", dx[0], dx[1], dx[2]);
        if(iter >= maxiter)
        {   
            intersection = false;
            break;
        }       
        
        updateRHSForRCINormal( b, xsi, normal, cutterElement, onBoundary); 
        residual = b.Norm2(); 
        iter++;
        
        //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi[0], xsi[1], xsi[2], residual, TOL14);
    } 
    
    if( (fabs(xsi[0])-1.0) > TOL7  || (fabs(xsi[1])-1.0) > TOL7 )    // line coordinate may be bigger than 1
    {
        //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi[0], xsi[1], xsi[2], residual, TOL14);
        intersection = false;
    }
        
    return intersection;
}
    
    

/*----------------------------------------------------------------------*
 |  RCI:    updates the systemmatrix for the                 u.may 08/07|
 |          computation of a curve-surface intersection                 |
 |          for the recovery of the curved interface                    |
 |          (line-surface intersection)                                 |   
 *----------------------------------------------------------------------*/
void Intersection::updateAForRCINormal(   
    Epetra_SerialDenseMatrix&                   A,
    const Epetra_SerialDenseVector&             xsi,
    const vector<Epetra_SerialDenseVector>&     normal,
    const DRT::Element*                         surfaceElement,
    const bool                                  onBoundary
    ) const
{   
    const int numNodesSurface = surfaceElement->NumNode();
    //Epetra_SerialDenseMatrix surfaceDeriv1(2,numNodesSurface);
    
    A.Scale(0.0);
    //shape_function_2D_deriv1( surfaceDeriv1, xsi[0], xsi[1], surfaceElement->Shape());
    const BlitzMat surfaceDeriv1(shape_function_2D_deriv1( xsi[0], xsi[1], surfaceElement->Shape()));
   
    //printf("num of nodes surfaces = %d\n", numOfNodesSurface);
    if(!onBoundary)
    {
        for(int i=0; i<numNodesSurface; i++)
        {
            const double* pos = surfaceElement->Nodes()[i]->X();
            for(int dim=0; dim<3; dim++)
            {
//                A[dim][0] += node->X()[dim] * surfaceDeriv1[i][0];
//                A[dim][1] += node->X()[dim] * surfaceDeriv1[i][1];
                A[dim][0] += pos[dim] * surfaceDeriv1(0, i);
                A[dim][1] += pos[dim] * surfaceDeriv1(1, i);
            }
        }
        
        for(int dim=0; dim<3; dim++)
            A[dim][2] = (-1.0) * ( normal[0][dim] * (-0.5) + normal[1][dim] * 0.5 );  
    }
    else
    {
//        const int numNodesLine = 3;   
//        Epetra_SerialDenseMatrix lineDeriv1(1,numNodesLine);
//        shape_function_1D_deriv1(lineDeriv1, xsi[2], DRT::Element::line3);
        const BlitzMat lineDeriv1(shape_function_1D_deriv1( xsi[2], DRT::Element::line3));
     
        for(int i=0; i<numNodesSurface; i++)
        {
            const double* pos = surfaceElement->Nodes()[i]->X();
            for(int dim=0; dim<3; dim++)
            {
//                A[dim][0] += node->X()[dim] * surfaceDeriv1[i][0];
//                A[dim][1] += node->X()[dim] * surfaceDeriv1[i][1]; // TODO: was this a bug??? switched indices
                A[dim][0] += pos[dim] * surfaceDeriv1(0,i);
                A[dim][1] += pos[dim] * surfaceDeriv1(1,i);
            }
        }
        
        for(int dim=0; dim<3; dim++)    
            for(int i = 0; i < 3; i++)
            {
                int index = i;
                if(i > 1)   index = 4;
                //A[dim][2] += (-1.0) * normal[index][dim] * lineDeriv1[i][0]; 
                A[dim][2] += (-1.0) * normal[index][dim] * lineDeriv1(0,i);
            }  
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    updates the right-hand-side for the              u.may 08/07|
 |          computation of a curve-surface intersection                 |     
 |          for the recovery of the curved surface (RCI)                |
 |          (line-surface intersection)                                 |   
 *----------------------------------------------------------------------*/
void Intersection::updateRHSForRCINormal( 
    Epetra_SerialDenseVector&                   b,
    const Epetra_SerialDenseVector&             xsi,    
    const vector<Epetra_SerialDenseVector>&     normal,
    const DRT::Element*                         surfaceElement,
    const bool                                  onBoundary
    ) const
{
    const int numNodesSurface = surfaceElement->NumNode();
    Epetra_SerialDenseVector surfaceFunct(numNodesSurface);
  
    b.Scale(0.0);  
    shape_function_2D( surfaceFunct, xsi[0], xsi[1], surfaceElement->Shape());     
  
    if(!onBoundary)
    {
        for(int i=0; i<numNodesSurface; i++)
        {
            const DRT::Node* node = surfaceElement->Nodes()[i];
            for(int dim=0; dim<3; dim++)    
                b[dim] += (-1.0) * node->X()[dim] * surfaceFunct[i];
        }
            
        for(int dim=0; dim<3; dim++)
            b[dim] += normal[0][dim] * 0.5*(1.0 - xsi[2]) + normal[1][dim] * 0.5*(1.0 + xsi[2]);
    }
    else
    {
        const int numNodesLine = 3;
        Epetra_SerialDenseVector lineFunct(numNodesLine);
        shape_function_1D( lineFunct, xsi[2], DRT::Element::line3);
       
        for(int i=0; i<numNodesSurface; i++)
        {
            const DRT::Node* node = surfaceElement->Nodes()[i];
            for(int dim=0; dim<3; dim++)
                b[dim] += (-1.0) * node->X()[dim] * surfaceFunct[i];
        }
        
       
        for(int i=0; i<numNodesLine; i++)
            for(int dim=0; dim<3; dim++)
            {
                int index = i;
                if(i > 1)   index = 4;
                
                b[dim] += normal[index][dim] * lineFunct[i]; 
            } 
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    computes the intersection between a              u.may 08/07|
 |          curve and a plane                      RCI                  |           
 *----------------------------------------------------------------------*/
bool Intersection::computeRecoveryPlane( 
        int&                                        lineIndex,
        Epetra_SerialDenseVector&                   xsi,
        const vector<Epetra_SerialDenseVector>&     plane,
        DRT::Element*                               surfaceElement
        ) const
{
    bool    intersection = true;
    int     begin, end;
    const int     numLines = surfaceElement->NumLine();
   
    // run over all lines (curves)
    if(lineIndex == -1)
    {
        begin = 0;
        end = numLines;
    }
    else
    {
        begin = lineIndex;
        end   = lineIndex + 1;
    }
    
    for(int i = begin; i < end; i++)
    {
        int                         iter = 0;
        const int                   maxiter = 50;
        double                      residual = 1.0;
        const DRT::Element*         lineElement = surfaceElement->Lines()[i];
        Epetra_SerialDenseMatrix    A(3,3);
        Epetra_SerialDenseVector    b(3);
        Epetra_SerialDenseVector    dx(3);
      
        intersection = true;
        xsi.Scale(0.0); 
        
        updateRHSForRCIPlane( b, xsi, plane, lineElement);
                        
        while( residual > TOL14 )
        {   
            updateAForRCIPlane( A, xsi, plane, lineElement, surfaceElement);
            
            if(!gaussElimination(A, b, dx, true, 3, 1))
            {
                intersection = false;
                break;  
            } 
             
            if(iter >= maxiter)
            {   
                intersection = false;
                break;
            }       
        
            //xsi = addTwoVectors(xsi, dx);
            xsi += dx;
  
            updateRHSForRCIPlane( b, xsi, plane, lineElement);
            residual = b.Norm2(); 
            iter++;
        } 
    
        if( (fabs(xsi[2])-1.0) > TOL7 )     // planes coordinate may be bigger than 1
        {   printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi[0], xsi[1], xsi[2], residual, TOL14);
            intersection = false;
        }  
        
        if(intersection)
        {   
            lineIndex = begin;
            break;
        }
    }
           
    return intersection;
}



/*----------------------------------------------------------------------*
 |  RCI:    updates the systemmatrix for the                 u.may 09/07|
 |          computation of a curve-surface intersection                 |
 |          for the recovery of the curved interface                    |
 |          (curve-plane intersection)                                  |   
 *----------------------------------------------------------------------*/
void Intersection::updateAForRCIPlane(   
        Epetra_SerialDenseMatrix&                   A,
        const Epetra_SerialDenseVector&             xsi,
        const vector<Epetra_SerialDenseVector>&     plane,
        const DRT::Element*                         lineElement,
        const DRT::Element*                         surfaceElement
        ) const
{   
    const int numNodesLine = lineElement->NumNode();
    const int numNodesSurface = 4;
    Epetra_SerialDenseMatrix surfaceDeriv(2, numNodesSurface);
    Epetra_SerialDenseMatrix lineDeriv(1,numNodesLine);
   
    A.Scale(0.0);
    shape_function_2D_deriv1(surfaceDeriv,  xsi[0],  xsi[1], DRT::Element::quad4);
    shape_function_1D_deriv1(lineDeriv, xsi[2], lineElement->Shape());
    
    for(int dim=0; dim<3; dim++)
        for(int i=0; i<numNodesSurface; i++)
        {
            A[dim][0] += plane[i][dim] * surfaceDeriv[i][0];
            A[dim][1] += plane[i][dim] * surfaceDeriv[i][1];
        }
        
    for(int i=0; i<numNodesLine; i++)
    {
        const DRT::Node* node = lineElement->Nodes()[i];
        for(int dim=0; dim<3; dim++)
            A[dim][2] += (-1.0) * node->X()[dim] * lineDeriv[i][0];   
    }        
}



/*----------------------------------------------------------------------*
 |  RCI:    updates the right-hand-side for the              u.may 09/07|
 |          computation of a curve-surface intersection                 |     
 |          for the recovery of the curved surface (RCI)                |
 |          (curve-plane intersection)                                  |   
 *----------------------------------------------------------------------*/
void Intersection::updateRHSForRCIPlane( 
    Epetra_SerialDenseVector&                   b,
    const Epetra_SerialDenseVector&             xsi,    
    const vector <Epetra_SerialDenseVector>&    plane,
    const DRT::Element*                         lineElement
    ) const
{
    const int numNodesLine    = lineElement->NumNode();
    const int numNodesSurface = 4;
//    Epetra_SerialDenseVector surfaceFunct(numNodesSurface);
//    Epetra_SerialDenseVector lineFunct(numNodesLine);
   
    const BlitzVec surfaceFunct(shape_function_2D( xsi[0], xsi[1], DRT::Element::quad4 )); 
    const BlitzVec lineFunct(shape_function_1D( xsi[2], lineElement->Shape()));
    
    b.Scale(0.0);
    for(int dim=0; dim<3; dim++)
        for(int i=0; i<numNodesSurface; i++)
        {
            b[dim] -= plane[i][dim] * surfaceFunct(i);
        }
    
    const DRT::Node** lineNodes = lineElement->Nodes();
    for(int i=0; i<numNodesLine; i++)
    { 
        const double* pos = lineNodes[i]->X();
        for(int dim=0; dim<3; dim++)    
        {
            b[dim] +=  pos[dim]  * lineFunct(i);        
        }
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    computes the normal to the interface edge of     u.may 08/07|
 |          the tetrahedron facet lying within this facet               |
 *----------------------------------------------------------------------*/  
void Intersection::computeIntersectionNormal( 
    const bool                              onBoundary,
    const int                               index1,
    const int                               index2,
    const int                               oppositePointIndex,
    const int                               globalHigherOrderIndex, 
    const vector<int>&                      tetraCornerIndices,
    const vector<Epetra_SerialDenseVector>& tetraCornerNodes,
    vector<Epetra_SerialDenseVector>&       plane,
    const DRT::Element*                     xfemElement,
    const tetgenio&                         out) const
{            
    
    Epetra_SerialDenseVector  p1(3);
    Epetra_SerialDenseVector  p2(3);
    Epetra_SerialDenseVector  p3(3);        
    
    
    if(!onBoundary)
    {
        for(int i=0; i<3; i++)
        {    
            p1[i] = tetraCornerNodes[3][i];
            p2[i] = tetraCornerNodes[index1][i];
            p3[i] = tetraCornerNodes[index2][i];             
        }
    }
    else
    {
        for(int i=0; i<3; i++)
        {
            p1[i] = out.pointlist[oppositePointIndex*3          + i];
            p2[i] = out.pointlist[tetraCornerIndices[index1]*3  + i];
            p3[i] = out.pointlist[tetraCornerIndices[index2]*3  + i];                       
        }
    }
                 
    // compute direction vectors of the plane 
    Epetra_SerialDenseVector  r1 = subtractsTwoVectors(p1, p2);
    Epetra_SerialDenseVector  r2 = subtractsTwoVectors(p3, p2);
 
    // normal of the plane
    Epetra_SerialDenseVector  n = computeCrossProduct(r1, r2);
    normalizeVector(n);
 
    // direction vector of the intersection line
    Epetra_SerialDenseVector  r = computeCrossProduct(n, r2);  
    normalizeVector(r);
 
    // computes the start point of the line
    Epetra_SerialDenseVector  m = computeLineMidpoint(p2, p3);
    
    if(!onBoundary)
        m = computeLineMidpoint(p2, p3);
    else
    {
        for(int i = 0; i < 3; i++)
            m[i] = out.pointlist[globalHigherOrderIndex*3+i];
    }
    
    // compute nodes of the normal to the interface edge of the tetrahedron
    plane[0] = addTwoVectors(m, r);               
    plane[1] = subtractsTwoVectors(m, r);
    plane[2] = addTwoVectors(plane[1], n);               
    plane[3] = addTwoVectors(plane[0], n);
    
    
    if(onBoundary)
    {
        for(int i = 0; i < 4; i++)
            elementToCurrentCoordinates(xfemElement, plane[i]);
        
        elementToCurrentCoordinates(xfemElement, m);
          
        plane[4] = m;
      
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    computes the normal to the interface edge of     u.may 11/07|
 |          two adjacent triangular faces                               |
 *----------------------------------------------------------------------*/  
void Intersection::computeIntersectionNormal( 
        const int                                 index1,
        const int                                 index2, 
        const int                                 faceIndex,
        const int                                 adjacentFaceIndex,
        const int                                 globalHigherOrderIndex,
        vector<Epetra_SerialDenseVector>&         plane,
        const DRT::Element*                       xfemElement,
        const tetgenio&                           out
        ) const
{

    int oppositePointIndex = -1;
    int adjacentOppositePointIndex = -1;
    
    for(int i = 0; i < 3; i++)
        if((out.trifacelist[faceIndex*3+i] != index1) && (out.trifacelist[faceIndex*3+i] != index2 ))
        {
            oppositePointIndex = out.trifacelist[faceIndex*3+i];
            break;
        }
        
    for(int i = 0; i < 3; i++)
        if((out.trifacelist[adjacentFaceIndex*3+i] != index1) && (out.trifacelist[adjacentFaceIndex*3+i] != index2 ))
        {
            adjacentOppositePointIndex = out.trifacelist[adjacentFaceIndex*3+i];
            break;
        }
         
     // compute averageNormal of two faces
     Epetra_SerialDenseVector p1(3);
     Epetra_SerialDenseVector p2(3);
     Epetra_SerialDenseVector p3(3);
     Epetra_SerialDenseVector p4(3);
     
    for(int i=0; i<3; i++)
    {
        p1[i] = out.pointlist[index1*3 + i];
        p2[i] = out.pointlist[index2*3 + i];
        p3[i] = out.pointlist[oppositePointIndex*3 + i];              
        p4[i] = out.pointlist[adjacentOppositePointIndex*3 + i];                       
    }
    
    elementToCurrentCoordinates(xfemElement, p1);  
    elementToCurrentCoordinates(xfemElement, p2);  
    elementToCurrentCoordinates(xfemElement, p3);  
    elementToCurrentCoordinates(xfemElement, p4);  
    
    Epetra_SerialDenseVector r1 = subtractsTwoVectors(p1,p2);
    Epetra_SerialDenseVector r2 = subtractsTwoVectors(p3,p2);
    Epetra_SerialDenseVector r3 = subtractsTwoVectors(p4,p2);
    
    Epetra_SerialDenseVector n1 = computeCrossProduct(r2, r1);
    Epetra_SerialDenseVector n2 = computeCrossProduct(r1, r3);
    
    Epetra_SerialDenseVector averageNormal = addTwoVectors(n1, n2);
    Epetra_SerialDenseVector rPlane = computeCrossProduct(n1, r1);
    
    for(int i = 0; i < 3; i++)
        averageNormal[i] = 0.5*averageNormal[i];
        
    normalizeVector(averageNormal);
    normalizeVector(rPlane);
 
    Epetra_SerialDenseVector m(3);
    for(int i = 0; i < 3; i++)
        m[i] = out.pointlist[globalHigherOrderIndex*3+i];

    elementToCurrentCoordinates(xfemElement, m);  
    
    // compute nodes of the normal to the interface edge of the tetrahedron
    plane[0] = addTwoVectors(m, averageNormal);               
    plane[1] = subtractsTwoVectors(m, averageNormal);
    plane[2] = addTwoVectors(plane[1], rPlane);               
    plane[3] = addTwoVectors(plane[0], rPlane);
    
   /* for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 3 ; j++)
            printf("plane = %f\t", plane[i][j]);
    
        printf("\n");
    }
    */
}



/*----------------------------------------------------------------------*
 |  RCI:    computes the normal to the interface edge of     u.may 08/07|
 |          the tetrahedron facet lying within this facet               |
 *----------------------------------------------------------------------*/  
void Intersection::computeIntersectionNormal( 
    const int                               steinerIndex,
    const int                               edgeIndex,
    const int                               oppositeIndex,
    vector<Epetra_SerialDenseVector>&       plane,
    const DRT::Element*                     xfemElement,
    const tetgenio&                         out) const
{            
    
    Epetra_SerialDenseVector  p1(3);
    Epetra_SerialDenseVector  p2(3);
    Epetra_SerialDenseVector  p3(3);        
    
  
    for(int i=0; i<3; i++)
    {
        p1[i] = out.pointlist[oppositeIndex*3  + i];
        p2[i] = out.pointlist[steinerIndex*3   + i];
        p3[i] = out.pointlist[edgeIndex*3      + i];                       
    }
           
    // compute direction vectors of the plane 
    Epetra_SerialDenseVector  r1 = subtractsTwoVectors(p1, p2);
    Epetra_SerialDenseVector  r2 = subtractsTwoVectors(p3, p2);
 
    // normal of the plane
    Epetra_SerialDenseVector  n = computeCrossProduct(r1, r2);
    normalizeVector(n);
 
    // direction vector of the intersection line
    Epetra_SerialDenseVector  r = computeCrossProduct(n, r2);  
    normalizeVector(r);
 
    // compute nodes of the normal to the interface edge of the tetrahedron
    plane[0] = addTwoVectors(p2, r);               
    plane[1] = subtractsTwoVectors(p2, r);
    plane[2] = addTwoVectors(plane[1], n);               
    plane[3] = addTwoVectors(plane[0], n);
    
    
    for(int i = 0; i < 4; i++)
        elementToCurrentCoordinates(xfemElement, plane[i]);
    
    elementToCurrentCoordinates(xfemElement, p2);
      
    plane[4] = p2;
}



/*----------------------------------------------------------------------*
 |  RCI:    computes the midpoint of a line                  u.may 08/07|
 *----------------------------------------------------------------------*/  
Epetra_SerialDenseVector Intersection::computeLineMidpoint( 
    const Epetra_SerialDenseVector& p1, 
    const Epetra_SerialDenseVector& p2) const
{
    Epetra_SerialDenseVector midpoint(3);
    
    for(int i=0; i<3; i++)
        midpoint[i] = (p1[i] + p2[i])*0.5;
        
    return midpoint;
}



/*----------------------------------------------------------------------*
 |  RCI:    searches for the face marker                     u.may 10/07|
 |          of a facet adjacent to a given edge of                      |
 |          of a given facet                                            |
 *----------------------------------------------------------------------*/       
void Intersection::findAdjacentFace(
    const int         edgeIndex1, 
    const int         edgeIndex2, 
    const int         faceMarker,
    int&              adjacentFaceMarker,
    int               faceIndex,
    int&              adjacentFaceIndex,
    tetgenio&   out) const
{

    bool    faceMarkerFound     = false;
    
    for(int i=0; i<out.numberoftrifaces; i++)
    {
        adjacentFaceMarker = out.trifacemarkerlist[i] - facetMarkerOffset_;
        adjacentFaceIndex = i;
        
        if((adjacentFaceMarker > -2) && (faceIndex != adjacentFaceIndex))
        {        
            int countPoints = 0;
            for(int j = 0; j < 3 ; j++)
            {
                int pointIndex = out.trifacelist[ i*3 + j ];
                if(pointIndex == edgeIndex1 || pointIndex == edgeIndex2)
                    countPoints++;
            }
            
            if(countPoints == 2) 
                faceMarkerFound = true;
        }
        if(faceMarkerFound)
            break;
        
    }
    
    if(!faceMarkerFound)
        adjacentFaceMarker = -2;
    
}



/*----------------------------------------------------------------------*
 |  RCI:    finds the global index of the point opposite     u.may 08/07|
 |          to an edge in the adjacent triangular face                  |
 *----------------------------------------------------------------------*/  
int Intersection::findEdgeOppositeIndex( 
    int                                 edgeIndex1,
    int                                 edgeIndex2, 
    int                                 adjacentFaceIndex,
    tetgenio&                           out) const
{
    int oppositePointIndex = -1;
  
    for(int i = 0; i < 3; i++)
        if((out.trifacelist[adjacentFaceIndex*3+i] != edgeIndex1) && (out.trifacelist[adjacentFaceIndex*3+i] != edgeIndex2 ))
        {
            oppositePointIndex = out.trifacelist[adjacentFaceIndex*3+i];
            break;
        }

    return oppositePointIndex;
}



/*----------------------------------------------------------------------*
 |  RCI:    searchs for the common edge of two               u.may 10/07|
 |          adjacent facets                                             |     
 *----------------------------------------------------------------------*/
bool Intersection::findCommonFaceEdge( 
        const int                           faceIndex1, 
        const int                           faceIndex2, 
        const vector<int>&                  adjacentFacesList,
        Epetra_SerialDenseVector&           edgePoint,
        Epetra_SerialDenseVector&           oppositePoint,
        const tetgenio&                     out
        ) const
{
    bool edgeFound = false;
    
    for(int i = 0; i < 2; i++)
    {
        for(int j = 0; j < 2; j++)
        {
            if(adjacentFacesList[faceIndex1*2+i+1] == adjacentFacesList[faceIndex2*2+j+1])
            {
                int index = 0;
                if(i == 0)  index = 1;
                
                for(int k = 0; k < 3; k++)
                {
                    edgePoint[k]        = out.pointlist[adjacentFacesList[faceIndex1*2+i+1]*3 + k];
                    oppositePoint[k]    = out.pointlist[adjacentFacesList[faceIndex1*2+index+1]*3 + k];
                }
                edgeFound = true;
                break;
            }   
        }
        if(edgeFound)
            break; 
    }
    
    return edgeFound;
}
                            
                    
                    
/*----------------------------------------------------------------------*
 |  RCI:    searches for the common edge of two              u.may 10/07|
 |          adjacent cutter elements                                    |
 |          corresponding to the common face edge                       |
 |          of face 1 and facet 2                                       |     
 *----------------------------------------------------------------------*/       
bool Intersection::findCommonCutterLine(  
    const int                                       faceIndex1, 
    const int                                       faceIndex2,
    int&                                            lineIndex,
    int&                                            cutterIndex) const
{
    bool comparison = false;
    // get line arrays, which are computed onthe fly within the cutter elements
    DRT::Element** lines1 = intersectingCutterElements_[faceIndex1]->Lines();
    DRT::Element** lines2 = intersectingCutterElements_[faceIndex2]->Lines();

    const int numLines1 = intersectingCutterElements_[faceIndex1]->NumLine(); 
    const int numLines2 = intersectingCutterElements_[faceIndex2]->NumLine(); 
    const int numNodes  = lines2[0]->NumNode();
    
    for(int i = 0; i < numLines1; i++)
    {
        for(int j = 0; j < numLines2; j++)
        {
            comparison = true;
            for(int k  = 0; k < numNodes; k++)
            {
                const DRT::Node* node1 = lines1[i]->Nodes()[k];
                const DRT::Node* node2 = lines2[j]->Nodes()[k];
                if(!comparePoints(node1, node2))
                {
                    comparison = false;
                    break;   
                }
            }  
             
            if(!comparison)
            {
                comparison = true;
                for(int k  = 0; k < numNodes; k++)
                {
                    DRT::Node* node1;
                    DRT::Node* node2;
                    if(k==2)
                    {
                        node1 = lines1[i]->Nodes()[k];
                        node2 = lines2[j]->Nodes()[k];
                    }
                    else
                    {
                        node1 = lines1[i]->Nodes()[k];
                        node2 = lines2[j]->Nodes()[1-k];
                    }
                    
                    if(!comparePoints(node1, node2))
                        comparison = false;
                }   
            }
            
            if(comparison)
            {
                lineIndex    =  i;
                cutterIndex =  faceIndex1;
                break;
            }
        }
        if(comparison) 
            break;
    }
    return comparison;
}



/*----------------------------------------------------------------------*
 |  RCI:    for the recovery of a higher-order node          u.may 10/07|
 |          by a plane - line element intersection                      |
 |          this method finds the line element of the given cutter      |
 |          element intersecting the plane                              |
 |          checking if the edge nodes of the correspondning            |
 |          facet edge are lying on this line element                   |
 *----------------------------------------------------------------------*/   
int Intersection::findIntersectingSurfaceEdge(
        const DRT::Element*                       xfemElement,
        DRT::Element*                             cutterElement,
        const Epetra_SerialDenseVector&           edgeNode1,
        const Epetra_SerialDenseVector&           edgeNode2) const
{

    int lineIndex = -1;
    Epetra_SerialDenseVector x1(1);
    Epetra_SerialDenseVector x2(1);
    
    Epetra_SerialDenseVector node1(edgeNode1);
    Epetra_SerialDenseVector node2(edgeNode2);
      
    elementToCurrentCoordinates(xfemElement, node1);
    elementToCurrentCoordinates(xfemElement, node2);
    
    x1[0] = node1[0];
    x2[0] = node2[0];
    
    DRT::Element*const* lines = cutterElement->Lines();
    for(int i = 0; i < cutterElement->NumLine(); i++)
    {
        const DRT::Element*  lineElement = lines[i];

        Epetra_SerialDenseVector xsi1(1);
        Epetra_SerialDenseVector xsi2(1);
        currentToElementCoordinates( lineElement, x1, xsi1);
        currentToElementCoordinates( lineElement, x2, xsi2);
        
        const bool check1 = checkPositionWithinElementParameterSpace( xsi1, lineElement->Shape());
        const bool check2 = checkPositionWithinElementParameterSpace( xsi2, lineElement->Shape());
        if( check1  &&  check2 )
        {   
            lineIndex = i;  //countIndex;
            break;
        }
    }
    
    /*
    for(int i = 0; i < lines[lineIndex]->NumNode(); i++ )
    {
        lines[lineIndex]->Nodes()[i]->Print(cout); 
        printf("\n");  
    }
    */
    return lineIndex;
}



/*----------------------------------------------------------------------*
 |  RCI:    stores the higher-order node in the pointlist    u.may 08/07|
 |          at the place of the linear node                             |
 *----------------------------------------------------------------------*/  
void Intersection::storeHigherOrderNode( 
    const bool                                  normal,
    const int                                   globalHigherOrderIndex,
    const int                                   lineIndex, 
    Epetra_SerialDenseVector&                   xsi, 
    DRT::Element*                               surfaceElement, 
    const DRT::Element*                         xfemElement, 
    tetgenio&                                   out
    ) const
{
    
    Epetra_SerialDenseVector xsiLine(3);
    for(int i = 0; i < 3; i++)
        xsiLine[i] = 0.0;
        
    if(normal)     elementToCurrentCoordinates(surfaceElement, xsi);
    else           
    {   
        xsiLine[0] = xsi[2];
        const DRT::Element* lineele = surfaceElement->Lines()[lineIndex];
        elementToCurrentCoordinates(lineele, xsiLine);
        for(int i=0; i<3; i++)
            xsi[i] = xsiLine[i];
    }       
    currentToElementCoordinates(xfemElement, xsi);
    
    //printf("xsiold0 = %20.16f\t, xsiold1 = %20.16f\t, xsiold2 = %20.16f\n", out.pointlist[index*3], out.pointlist[index*3+1], out.pointlist[index*3+2]);
    
    for(int i = 0; i < 3; i++)
        out.pointlist[globalHigherOrderIndex*3+i]   = xsi[i];  
   
    //printf("xsi0    = %20.16f\t, xsi1    = %20.16f\t, xsi2    = %20.16f\n", xsi[0], xsi[1], xsi[2]);
    //printf("\n");
}



/*----------------------------------------------------------------------*
 |  RCI:    stores domain integration cells                  u.may 11/07|
 *----------------------------------------------------------------------*/  
void Intersection::addCellsToDomainIntCellsMap(
        const DRT::Element*             xfemElement,
        map< int, DomainIntCells >&     domainintcells,
        tetgenio&                       out) const
{   
    
    // store domain integration cells
    DomainIntCells   listDomainICPerElement;

    for(int i=0; i<out.numberoftetrahedra; i++ )
    {   
        vector< vector<double> > tetrahedronCoord;
        for(int j = 0; j < out.numberofcorners; j++)
        {
            vector<double> tetnodes(3);
            for(int k = 0; k < 3; k++)
                tetnodes[k] = out.pointlist[out.tetrahedronlist[i*out.numberofcorners+j]*3+k];
              
            tetrahedronCoord.push_back(tetnodes);    
        }
        listDomainICPerElement.push_back(DomainIntCell(DRT::Element::tet10, tetrahedronCoord));                 
    }
    domainintcells.insert(make_pair(xfemElement->Id(),listDomainICPerElement));
}



/*----------------------------------------------------------------------*
 |  RCI:    stores boundary integration cells                u.may 11/07|
 *----------------------------------------------------------------------*/  
void Intersection::addCellsToBoundaryIntCellsMap(
    const int                         		trifaceIndex,
    const int                         		cornerIndex, 
    const int                         		globalHigherOrderIndex, 
    const int                         		faceMarker, 
    vector<vector<double> >&                domainCoord, 
    vector<vector<double> >&                boundaryCoord, 
    const DRT::Element*						xfemElement,
    const tetgenio&               			out) const
{                
	
    //vector<double> trinodes(3);
  
//    printf("i length = %d\n", domainCoord.size());
//    printf("j length = %d\n", domainCoord[0].size());
//    for(int i = 0; i < 6; i++)
//    	for(int j = 0; j < 3; j++)
//    		printf("vector = %f\n ", domainCoord[i][j]);
    
    
    // store corner node
    {
    BlitzVec eleCoordDomainCorner(3);
    for(int k = 0; k < 3; k++)
        eleCoordDomainCorner(k) = out.pointlist[(out.trifacelist[trifaceIndex*3+cornerIndex])*3+k];
  
    domainCoord[cornerIndex][0] = eleCoordDomainCorner(0);
    domainCoord[cornerIndex][1] = eleCoordDomainCorner(1);
    domainCoord[cornerIndex][2] = eleCoordDomainCorner(2);
    
    const BlitzVec physCoordCorner = elementToCurrentCoordinates(xfemElement, eleCoordDomainCorner);
    
    //domainCoord.push_back(trinodes);
    
    const BlitzVec eleCoordBoundaryCorner = CurrentToSurfaceElementCoordinates(intersectingCutterElements_[faceMarker], physCoordCorner);
    
//    cout << "physcood = " << physCoord << "    ";
//    cout << "elecood = " << eleCoord << endl;
//    cout << "cornerIndex = " << cornerIndex << endl;
    boundaryCoord[cornerIndex][0] = eleCoordBoundaryCorner(0);
    boundaryCoord[cornerIndex][1] = eleCoordBoundaryCorner(1);
    boundaryCoord[cornerIndex][2] = 0.0;
    
    //boundaryCoord.push_back(trinodes);
    }
    
    {
    // store higher order node
    BlitzVec eleCoordDomaninHO(3);
    for(int k = 0; k < 3; k++)
        eleCoordDomaninHO(k) = out.pointlist[globalHigherOrderIndex*3+k];
 
    
    
    domainCoord[cornerIndex+3][0] = eleCoordDomaninHO(0);
    domainCoord[cornerIndex+3][1] = eleCoordDomaninHO(1);
    domainCoord[cornerIndex+3][2] = eleCoordDomaninHO(2);
    
    const BlitzVec physCoordHO = elementToCurrentCoordinates(xfemElement, eleCoordDomaninHO);
    
    //domainCoord.push_back(trinodes); 
    
    const BlitzVec eleCoordBoundaryHO = CurrentToSurfaceElementCoordinates(intersectingCutterElements_[faceMarker], physCoordHO);
    
//    cout << "physcood = " << physCoord << "    ";
//    cout << "elecood = " << eleCoord << endl;
    
    boundaryCoord[cornerIndex+3][0] = eleCoordBoundaryHO(0);
    boundaryCoord[cornerIndex+3][1] = eleCoordBoundaryHO(1);
    boundaryCoord[cornerIndex+3][2] = 0.0;
    
    //boundaryCoord.push_back(trinodes);
    }
}



/*----------------------------------------------------------------------*
 |  DB:     Debug only                                       u.may 06/07|
 *----------------------------------------------------------------------*/  
void Intersection::debugXAABBIntersection(
        const Epetra_SerialDenseMatrix cutterXAABB, 
        const Epetra_SerialDenseMatrix xfemXAABB,
        const DRT::Element* cutterElement,
        const DRT::Element* xfemElement,
        const int noC,
        const int noX
        ) const
{
	cout << endl;
    cout << "===============================================================" << endl;
	cout << "Debug Intersection of XAABB's" << endl;
	cout << "===============================================================" << endl;
	cout << endl;
	cout << "CUTTER ELEMENT " << noC << " :" << endl;
	cout << endl;
	for(int jE = 0; jE < cutterElement->NumNode(); jE++)
	{
		cutterElement->Nodes()[jE]->Print(cout);
		cout << endl;
	}
	cout << endl;
    cout << endl;
    cout << "CUTTER XAABB: "<< "                     " << "XFEM XAABB: " << endl;
    cout << endl;
    cout <<  "minX = " << cutterXAABB(0,0) << "      " << "maxX = " << cutterXAABB(0,1) << "      " ;
    cout <<  "minX = " << xfemXAABB(0,0)   << "      " << "maxX = " << xfemXAABB(0,1)   << endl;  
    cout <<  "minY = " << cutterXAABB(1,0) << "      " << "maxY = " << cutterXAABB(1,1) << "      " ;
    cout <<  "minY = " << xfemXAABB(1,0)   << "      " << "maxY = " << xfemXAABB(1,1)   << endl;   
    cout <<  "minZ = " << cutterXAABB(2,0) << "      " << "maxZ = " << cutterXAABB(2,1) << "      " ;
    cout <<  "minZ = " << xfemXAABB(2,0)   << "      " << "maxZ = " << xfemXAABB(2,1) << endl;
    cout << endl;
    cout << endl;
	cout << "XFEM ELEMENT " << noX << " :" << endl;
	cout << endl; 
	for(int jE = 0; jE < xfemElement->NumNode(); jE++)
	{
		xfemElement->Nodes()[jE]->Print(cout);
		cout << endl;
	}
    cout << endl;
    cout << endl;
    cout << "CUTTER XAABB: "<< "                     " << "XFEM XAABB: " << endl;
    cout << endl;
    cout <<  "minX = " << cutterXAABB(0,0) << "      " << "maxX = " << cutterXAABB(0,1) << "      " ;
    cout <<  "minX = " << xfemXAABB(0,0)   << "      " << "maxX = " << xfemXAABB(0,1)   << endl;  
    cout <<  "minY = " << cutterXAABB(1,0) << "      " << "maxY = " << cutterXAABB(1,1) << "      " ;
    cout <<  "minY = " << xfemXAABB(1,0)   << "      " << "maxY = " << xfemXAABB(1,1)   << endl;   
    cout <<  "minZ = " << cutterXAABB(2,0) << "      " << "maxZ = " << cutterXAABB(2,1) << "      " ;
    cout <<  "minZ = " << xfemXAABB(2,0)   << "      " << "maxZ = " << xfemXAABB(2,1) << endl;
	cout << endl;
	cout << endl;    
    cout << "===============================================================" << endl;
    cout << "End Debug Intersection of XAABB's" << endl;
	cout << "===============================================================" << endl;
	cout << endl; cout << endl; cout << endl;


}



/*----------------------------------------------------------------------*
 |  DB:     Debug only                                       u.may 06/07|
 *----------------------------------------------------------------------*/  
void Intersection::debugNodeWithinElement(
        const DRT::Element* element,
        const DRT::Node* node,
        const Epetra_SerialDenseVector& xsi,
        const int noE,
        const int noN,
        const bool within
        ) const
{
    const int numnodes = element->NumNode();
    vector<int> actParams(1,0);
    Epetra_SerialDenseVector funct(numnodes);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    Epetra_SerialDenseVector x(3);
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
      
    params.set("action","calc_Shapefunction");
    actParams[0] = numnodes;   
    
    dserror("we don't use Evaluate anymore, so thius function does not make sence!");
    //element->Evaluate(params, dummyDis, actParams, emptyM, emptyM , funct, xsi, emptyV);        
    
    for(int dim=0; dim<3; dim++)
        x(dim) = 0.0;
    
    for(int dim=0; dim<3; dim++)
        for(int i=0; i<numnodes; i++)
        {
            x(dim) += element->Nodes()[i]->X()[dim] * funct(i);
        }
        
    cout << endl;
    cout << "===============================================================" << endl;
    cout << "Debug Node within element" << endl;
    cout << "===============================================================" << endl;
    cout << endl;
    cout << "ELEMENT " << noE << " :" << endl;
    cout << endl;
/*    for(int jE = 0; jE < element->NumNode(); jE++)
    {
        element->Nodes()[jE]->Print(cout);
        cout << endl;
    }
*/
    cout << endl;
    cout << endl;
    cout << "NODE " << noN << " :" << endl;
    cout << endl; 
        node->Print(cout);
    cout << endl;
    cout << endl;
    cout << "XSI :";
    cout << "   r = " << xsi[0] << "     s = " <<  xsi[1] << "     t = "  << xsi[2] << endl;
    cout << endl; 
    cout << endl;
    cout << "CURRENT COORDINATES :";
    cout << "   x = " << x[0] << "     y = " <<  x[1] << "     z = "  << x[2] << endl;
    cout << endl; 
    cout << endl;
    if(within) cout << "NODE LIES WITHIN ELEMENT" << endl;
    else            cout << "NODE DOES NOT LIE WITHIN ELEMENT" << endl;
    cout << endl;
    cout << endl;
    cout << "===============================================================" << endl;
    cout << "End Debug Node within element" << endl;
    cout << "===============================================================" << endl;
    cout << endl; cout << endl; cout << endl;
}



/*----------------------------------------------------------------------*
 |  DB:     Debug only                                       u.may 06/07|
 *----------------------------------------------------------------------*/  
void Intersection::debugTetgenDataStructure(
        const DRT::Element*               element) const
{    
    cout << endl;
    cout << "===============================================================" << endl;
    cout << "Debug Tetgen Data Structure " << endl;
    cout << "===============================================================" << endl;
    cout << endl;
    cout << "POINT LIST " << " :" << endl;
    cout << endl;
    Epetra_SerialDenseVector xsi(3);
    for(unsigned int i = 0; i < pointList_.size(); i++)
    {
        for(int j = 0; j< 3; j++)
        {
            xsi[j] = pointList_[i].coord[j];
        }
        elementToCurrentCoordinates(element, xsi);
        
        cout << i << ".th point:   ";
        for(int j = 0; j< 3; j++)
        {
            //cout << setprecision(10) << pointList_[i].coord[j] << "\t"; 
             printf("%20.16f\t", pointList_[i].coord[j] );
        }
        cout << endl;
        cout << endl;
        
      /*  for(int j = 0; j< 3; j++)
        {
            cout << xsi[j] << "\t";
        }
        cout << endl;
        cout << endl;*/
    }
    cout << endl;
    cout << endl;
    
    cout << endl;
    cout << "SEGMENT LIST " << " :" << endl;
    cout << endl;
    for(unsigned int i = 0; i < segmentList_.size(); i++)
    {
        cout << i << ".th segment:   ";
        int count = 0;
        for(unsigned int j = 0; j < segmentList_[i].size(); j++)
                cout << segmentList_[i][count++] << "\t";
        
        count = 0;
        for(unsigned int j = 0; j < surfacePointList_[i].size(); j++)
                cout << surfacePointList_[i][count++] << "\t";
                
        cout << endl;
        cout << endl;
    }
    cout << endl;
    cout << endl;
    
    cout << endl;
    cout << "TRIANGLE LIST " << " :" << endl;
    cout << endl;
    for(unsigned int i = 0; i < triangleList_.size(); i++)
    {
        cout << i << ".th triangle:   ";
        for(int j = 0; j< 3; j++)
        {
            cout << triangleList_[i][j] << "\t";
        }
        cout << endl;
        cout << endl;
    }
    cout << endl;
    cout << endl;
  
    cout << "===============================================================" << endl;
    cout << "Debug Tetgen Data Structure" << endl;
    cout << "===============================================================" << endl;
    cout << endl; cout << endl; cout << endl;
}



/*----------------------------------------------------------------------*
 |  Debug only                                               u.may 06/07|
 *----------------------------------------------------------------------*/
void Intersection::debugTetgenOutput( 	tetgenio& in,
										tetgenio& out, 
										const DRT::Element*   element,
    									vector<int>& elementIds) const
{
	char* tetgenIn = "tetgenPLC";
	char* tetgenOut = "tetgenMesh";
	char tetgenInId[30];
	char tetgenOutId[30];
		
	for(unsigned int i = 0; i < elementIds.size(); i++)
	{
		if(element->Id()== elementIds[i])
		{
			// change filename
			sprintf(tetgenInId,"%s%d", tetgenIn, elementIds[i]);
			sprintf(tetgenOutId,"%s%d", tetgenOut, elementIds[i]);
			
			// write piecewise linear complex
			in.save_nodes(tetgenInId);
	    	in.save_poly(tetgenInId);
	      
	    	// write tetrahedron mesh
	    	out.save_elements(tetgenOutId);
	    	out.save_nodes(tetgenOutId);
	    	out.save_faces(tetgenOutId);
	    	
	    	cout << "Saving tetgen output for the " <<  elementIds[i] << ".xfem element" << endl;
            flush(cout);
	    }
    }
}



/*----------------------------------------------------------------------*
 |  DB:     Debug only                                       u.may 09/07|
 *----------------------------------------------------------------------*/  
void Intersection::printTetViewOutput(
    int             index,
    tetgenio&       out) const
{    
    
    FILE *outFile;
    char filename[100];
  

    sprintf(filename, "tetgenMesh%d.node", index);
   
    outFile = fopen(filename, "w");
    fprintf(outFile, "%d  %d  %d  %d\n", out.numberofpoints, out.mesh_dim,
          out.numberofpointattributes, out.pointmarkerlist != NULL ? 1 : 0);
    for (int i = 0; i < out.numberofpoints; i++) 
    {
        fprintf(outFile, "%d  %.16g  %.16g  %.16g", i, out.pointlist[i * 3], out.pointlist[i * 3 + 1], out.pointlist[i * 3 + 2]);
    
        for (int j = 0; j < out.numberofpointattributes; j++) 
        {
            fprintf(outFile, "  %.16g", out.pointattributelist[i * out.numberofpointattributes + j]);
        }
        if (out.pointmarkerlist != NULL) 
        {
            fprintf(outFile, "  %d", out.pointmarkerlist[i]);
        }
        fprintf(outFile, "\n");
    }
    fclose(outFile);    
}



/*----------------------------------------------------------------------*
 |  DB:     Debug only                                       u.may 09/07|
 *----------------------------------------------------------------------*/  
void Intersection::printTetViewOutputPLC(
        const DRT::Element*   xfemElement,
        int                   index,
        tetgenio&             in) const
{    
    
    FILE *outFile;
    char filename[100];
    Epetra_SerialDenseVector xsi(3);

    sprintf(filename, "tetgenPLC%d.node", index);
   
    outFile = fopen(filename, "w");
    fprintf(outFile, "%d  %d  %d  %d\n", in.numberofpoints, in.mesh_dim,
          in.numberofpointattributes, in.pointmarkerlist != NULL ? 1 : 0);
    for (int i = 0; i < in.numberofpoints; i++) 
    {
        
        for(int j = 0; j < 3; j++)
            xsi[j] = in.pointlist[i*3 + j];
        
        elementToCurrentCoordinates(xfemElement, xsi);
        
        fprintf(outFile, "%d  %.16g  %.16g  %.16g", i, xsi[0], xsi[1], xsi[2]);
    
        for (int j = 0; j < in.numberofpointattributes; j++) 
        {
            fprintf(outFile, "  %.16g", in.pointattributelist[i * in.numberofpointattributes + j]);
        }
        if (in.pointmarkerlist != NULL) 
        {
            fprintf(outFile, "  %d", in.pointmarkerlist[i]);
        }
        fprintf(outFile, "\n");
    }
    fclose(outFile);    
}


void Intersection::debugFaceMarker(
		const int 						eleId,
		tetgenio&						out) const
{
	
	ofstream f_system("element_faceMarker.pos");
	f_system << "View \" Face Markers" << " \" {" << endl;
	
	for(int i=0; i<out.numberoftrifaces; i++)
    {
        int trifaceMarker = out.trifacemarkerlist[i] - facetMarkerOffset_;
       
        if(trifaceMarker > -2)
        {
        	vector< vector<double> > triface(3, vector<double>(3));
        	for(int j = 0; j < 3; j++)
        		for(int k = 0; k < 3; k++)
        			triface[j][k] = out.pointlist[out.trifacelist[i*3+j]*3 + k];
        		
        	f_system << IO::GMSH::trifaceToString(trifaceMarker, triface) << endl;
        }
	}
	f_system << "};" << endl;
	//f_system << IO::GMSH::getConfigString(2);	
	f_system.close();
}


void Intersection::debugXFEMConditions(
	const RCP<DRT::Discretization>  cutterdis) const
{
	
	vector< DRT::Condition * >	xfemConditions;
	cutterdis->GetCondition ("XFEMCoupling", xfemConditions);
	
	ofstream f_system("element_xfemconditions.pos");
	f_system << "View \" XFEM conditions" << " \" {" << endl;
	
	for(unsigned int i=0; i<xfemConditions.size(); i++)
    {
        map<int, RCP<DRT::Element > >  geometryMap = xfemConditions[i]->Geometry();
        map<int, RCP<DRT::Element > >::iterator iterGeo;   
        //if(geometryMap.size()==0)   printf("geometry does not obtain elements\n");
      	//printf("size of %d.geometry map = %d\n",i, geometryMap.size());
     	
        for(iterGeo = geometryMap.begin(); iterGeo != geometryMap.end(); iterGeo++ )
        { 
            const DRT::Element*  cutterElement = iterGeo->second.get();
            f_system << IO::GMSH::elementToString(i, cutterElement) << endl;
        }
    }
	
	f_system << "};" << endl;
	f_system.close();
}


void Intersection::debugIntersection(
		const DRT::Element*			xfemElement,
		vector <DRT::Element*> 	cutterElements) const
{
	
	ofstream f_system("intersection.pos");
	f_system << "View \" Intersection" << " \" {" << endl;
	
	f_system << IO::GMSH::elementToString(0, xfemElement) << endl;
	
	for(unsigned int i=0; i<cutterElements.size(); i++)
		f_system << IO::GMSH::elementToString(i+1, cutterElements[i]) << endl;
    
	
	f_system << "};" << endl;
	f_system.close();
}


void Intersection::debugXAABBs(
	const int							id,
	const Epetra_SerialDenseMatrix&     cutterXAABB, 
	const Epetra_SerialDenseMatrix&     xfemXAABB) const
{
	char filename[100];
	sprintf(filename, "element_XAABB%d.pos", id);
	   
	ofstream f_system(filename);
	f_system << "View \" XAABB " << " \" {" << endl;
	std::vector < std::vector <double> > nodes(8, vector<double> (3, 0.0));
	
	//cutterXAABB
	nodes[0][0] = cutterXAABB(0,0);	nodes[0][1] = cutterXAABB(1,0);	nodes[0][2] = cutterXAABB(2,0);	// node 0	
	nodes[1][0] = cutterXAABB(0,1);	nodes[1][1] = cutterXAABB(1,0);	nodes[1][2] = cutterXAABB(2,0);	// node 1
	nodes[2][0] = cutterXAABB(0,1);	nodes[2][1] = cutterXAABB(1,1);	nodes[2][2] = cutterXAABB(2,0);	// node 2
	nodes[3][0] = cutterXAABB(0,0);	nodes[3][1] = cutterXAABB(1,1);	nodes[3][2] = cutterXAABB(2,0);	// node 3
	nodes[4][0] = cutterXAABB(0,0);	nodes[4][1] = cutterXAABB(1,0);	nodes[4][2] = cutterXAABB(2,1);	// node 4
	nodes[5][0] = cutterXAABB(0,1);	nodes[5][1] = cutterXAABB(1,0);	nodes[5][2] = cutterXAABB(2,1);	// node 5
	nodes[6][0] = cutterXAABB(0,1);	nodes[6][1] = cutterXAABB(1,1);	nodes[6][2] = cutterXAABB(2,1);	// node 6
	nodes[7][0] = cutterXAABB(0,0);	nodes[7][1] = cutterXAABB(1,1);	nodes[7][2] = cutterXAABB(2,1);	// node 7
	
	f_system << IO::GMSH::XAABBToString(id+1, nodes) << endl;
	
	//xfemXAABB
	nodes[0][0] = xfemXAABB(0,0);	nodes[0][1] = xfemXAABB(1,0);	nodes[0][2] = xfemXAABB(2,0);	// node 0	
	nodes[1][0] = xfemXAABB(0,1);	nodes[1][1] = xfemXAABB(1,0);	nodes[1][2] = xfemXAABB(2,0);	// node 1
	nodes[2][0] = xfemXAABB(0,1);	nodes[2][1] = xfemXAABB(1,1);	nodes[2][2] = xfemXAABB(2,0);	// node 2
	nodes[3][0] = xfemXAABB(0,0);	nodes[3][1] = xfemXAABB(1,1);	nodes[3][2] = xfemXAABB(2,0);	// node 3
	nodes[4][0] = xfemXAABB(0,0);	nodes[4][1] = xfemXAABB(1,0);	nodes[4][2] = xfemXAABB(2,1);	// node 4
	nodes[5][0] = xfemXAABB(0,1);	nodes[5][1] = xfemXAABB(1,0);	nodes[5][2] = xfemXAABB(2,1);	// node 5
	nodes[6][0] = xfemXAABB(0,1);	nodes[6][1] = xfemXAABB(1,1);	nodes[6][2] = xfemXAABB(2,1);	// node 6
	nodes[7][0] = xfemXAABB(0,0);	nodes[7][1] = xfemXAABB(1,1);	nodes[7][2] = xfemXAABB(2,1);	// node 7
	
	f_system << IO::GMSH::XAABBToString(0, nodes) << endl;
	
	f_system << "};" << endl;
	f_system.close();
}


#endif  // #ifdef CCADISCRET
