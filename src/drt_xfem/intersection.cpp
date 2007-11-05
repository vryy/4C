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
#include "../drt_xfem/intersection_math.H"
#include "../drt_xfem/integrationcell.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_element.H"
//#include "../drt_f3/fluid3.H"

#ifdef PARALLEL
#include "../drt_lib/drt_exporter.H"
#include <mpi.h>
#endif 

using namespace XFEM;
using namespace DRT::Utils;

/*----------------------------------------------------------------------*
 |  MAIN:   computes the interface between the xfem          u.may 06/07|
 |          discretization and the cutter discretization.               |
 |          It returns a list of intersected xfem elements              |
 |          and their integrations cells.                               |
 *----------------------------------------------------------------------*/  
void Intersection::computeIntersection( const RefCountPtr<DRT::Discretization>  xfemdis,
                                        const RefCountPtr<DRT::Discretization>  cutterdis,  
                                        map< int, DomainIntCells >&   domainintcells,
                                        map< int, BoundaryIntCells >&   boundaryintcells)
{
    
	
    bool xfemIntersection; 
    vector< DRT::Condition * >      		xfemConditions;
    map< int, vector< DRT::Element* > >     cutterElementMap;
    DRT::Element*                   		xfemElement;
    vector< DRT::Element* > 				cutterElements;

    countMissedPoints_ = 0;
    
#ifdef PARALLEL
	const int cmyrank = cutterdis->Comm().MyPID();
	const int xnumproc = xfemdis->Comm().NumProc();
	const int cnumproc = cutterdis->Comm().NumProc();
	map< int, RefCountPtr<DRT::Node> >  cutterNodeMap;
	map< int, set<int> > xfemCutterIdMap;
	vector< int > conditionEleCount;
	
	if(cnumproc != xnumproc)
		dserror("the number of processors for xfem and cutter discretizations have to equal each other");		
#endif

	
    // obtain vector of pointers to all xfem conditions of all cutter discretizations
    cutterdis->GetCondition ("XFEMCoupling", xfemConditions);
        
    if(xfemConditions.size()==0)
        dserror("number of fsi xfem conditions = 0");
      
       
#ifdef PARALLEL
		adjustCutterElementNumbering(cutterdis, xfemConditions, conditionEleCount);
#endif         
       
       
    //  k < xfemdis->NumMyColElements()
    for(int k = 0; k < xfemdis->NumMyColElements(); k++)
    {
        xfemIntersection = false;
        xfemElement = xfemdis->lColElement(k);
        initializeXFEM(xfemElement);  
        
        Epetra_SerialDenseMatrix xfemXAABB = computeFastXAABB(xfemElement);
            
        startPointList();
        
        //printf("size of xfem condition = %d\n", xfemConditions.size());
        for(unsigned int i=0; i<xfemConditions.size(); i++)
        {
            map< int, RefCountPtr<DRT::Element > >  geometryMap = xfemConditions[i]->Geometry();
            map<int, RefCountPtr<DRT::Element > >::iterator iterGeo;   
            //if(geometryMap.size()==0)   printf("geometry does not obtain elements\n");
          	//printf("size of %d.geometry map = %d\n",i, geometryMap.size());
         	
            for(iterGeo = geometryMap.begin(); iterGeo != geometryMap.end(); iterGeo++ )
            { 
                DRT::Element*  cutterElement = iterGeo->second.get();
                if(cutterElement == NULL) dserror("geometry does not obtain elements");
            	
                Epetra_SerialDenseMatrix    cutterXAABB = computeFastXAABB(cutterElement);                                 
                bool intersected = intersectionOfXAABB(cutterXAABB, xfemXAABB);    
                //debugXAABBIntersection( cutterXAABB, xfemXAABB, cutterElement, xfemElement, i, k);
                                          
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
					if(i > 0) addToCutterId = conditionEleCount[i];

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
          	}// for-loop over all geometryMap.size()
		 }// for-loop over all xfemConditions.size()           

#ifdef PARALLEL
	} // for loop over all xfem elements} 

	getCutterElementsInParallel(xfemConditions,  conditionEleCount, cutterElementMap,
	 							cutterNodeMap,  xfemCutterIdMap,  xfemdis, cutterdis);
	 										   
    for(int k = 0; k < xfemdis->NumMyColElements(); k++)
    {
        xfemIntersection = false;
        xfemElement = xfemdis->lColElement(k);
        initializeXFEM(xfemElement);    
        startPointList();

#endif         
        
        if(cutterElementMap.find(xfemElement->LID()) != cutterElementMap.end())
            cutterElements = cutterElementMap.find(k)->second;
        else
            cutterElements.resize(0);
            
        for(unsigned int i = 0; i < cutterElements.size(); i++ )
        {
            numInternalPoints_= 0; 
            numBoundaryPoints_ = 0;
            vector< InterfacePoint >  interfacePoints; 
            
            if(cutterElements[i] == NULL) dserror("cutter element is null\n");
           
            initializeCutter(cutterElements[i]); 

            // collect internal points
            for(int m=0; m<cutterElements[i]->NumLine() ; m++)                    
                collectInternalPoints( xfemElement, cutterElements[i], cutterElements[i]->Nodes()[m],
                                        interfacePoints, k, m);
          
            // collect intersection points                                   
            for(int m=0; m<xfemElement->NumLine() ; m++) 
                if(collectIntersectionPoints(   cutterElements[i], xfemElement->Lines()[m],
                                                interfacePoints, 0, m, false, xfemIntersection))                                         
                    storeIntersectedCutterElement(cutterElements[i]); 
            
            
            for(int m=0; m<cutterElements[i]->NumLine() ; m++)                                              
                for(int p=0; p<xfemElement->NumSurface() ; p++) 
                    if(collectIntersectionPoints(   xfemElement->Surfaces()[p], 
                                                    cutterElements[i]->Lines()[m], interfacePoints,
                                                    p, m, true, xfemIntersection)) 
                        storeIntersectedCutterElement(cutterElements[i]); 
            
            // order interface points           
            if(interfacePoints.size()!=0)
            {
#ifdef QHULL
                computeConvexHull( xfemElement, cutterElements[i], interfacePoints);
#else
                dserror("Set QHULL flag to use XFEM intersections!!!");
#endif
            }    
            interfacePoints.clear();     
        }// for-loop over all cutter elements
        
        
        if(xfemIntersection)
        {                                                                          
            //debugTetgenDataStructure(xfemElement);
            computeCDT(xfemElement, cutterElements[0], domainintcells);
        }
       
    }// for-loop over all  actdis->NumMyColElements()
    
    //debugDomainIntCells(domainintcells,2);
    cout << endl;
    if(countMissedPoints_ > 0)
    	cout << "Number of missed points during the recovery copy = " << countMissedPoints_ << endl;
    	
    cout << "Intersection computed sucessfully ";
#ifdef PARALLEL
    flush(cout);
	cout << " rank = " << cmyrank;
    flush(cout);
#endif
    cout << endl;
    cout << endl;
  
}



/*----------------------------------------------------------------------*
 |  INIT:   initializes the private members of the           u.may 08/07|
 |          current xfem element                                        |
 *----------------------------------------------------------------------*/
void Intersection::initializeXFEM(  
    DRT::Element* xfemElement)
{ 
    xfemDistype_ = xfemElement->Shape();
    
    numXFEMSurfaces_ = xfemElement->NumSurface(); 
    numXFEMLines_ = xfemElement->NumLine(); 
   
    numXFEMCornerNodes_  = getNumberOfElementCornerNodes(xfemDistype_);
    
    pointList_.clear();
    segmentList_.clear();                 
    surfacePointList_.clear();                
    triangleList_.clear();
    
    intersectingCutterElements_.clear();
    faceMarker_.clear();
    
    eleLinesSurfaces_.clear();
    eleNodesSurfaces_.clear();
    eleNumberingSurfaces_.clear();
    eleRefCoordinates_.clear(); 
     
    segmentList_.resize(numXFEMSurfaces_);
    surfacePointList_.resize(numXFEMSurfaces_);    
    
    eleLinesSurfaces_ = getEleNodeNumbering_lines_surfaces(xfemDistype_);
    eleNodesSurfaces_ = getEleNodeNumbering_nodes_surfaces(xfemDistype_);
    eleNumberingSurfaces_ = getEleNodeNumberingSurfaces(xfemDistype_);
    eleRefCoordinates_ = getEleNodeNumbering_nodes_reference(xfemDistype_);
    
    higherorder_ = false;
}



/*----------------------------------------------------------------------*
 |  INIT:   initializes the private members of the           u.may 08/07|
 |          current cutter element                                      |
 *----------------------------------------------------------------------*/
void Intersection::initializeCutter(  
    DRT::Element* cutterElement)
{  
    cutterDistype_ = cutterElement->Shape();
}


#ifdef PARALLEL
/*----------------------------------------------------------------------*
 |  PARALLEL:   collects all cutter elements on              u.may 10/07|
 |              each processor   										|
 *----------------------------------------------------------------------*/
void Intersection::adjustCutterElementNumbering( 
	const RefCountPtr<DRT::Discretization>&  	cutterdis,  
	vector< DRT::Condition * >&   				xfemConditions,
	vector<int>& 								conditionEleCount)
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
    const RefCountPtr<DRT::Discretization>      cutterdis,
    vector<int>&                                conditionSend, 
    vector<int>&                                lengthSend, 
    int&                                        nodeSetSizeSend,  
    vector<int>&                                nodeVectorSend,
    vector<char>&                               cutterDataSend )
{
    
    set<int> nodeSet;
    vector< DRT::Condition * > xfemConditions;
    // obtain vector of pointers to all xfem conditions of all cutter discretizations
    cutterdis->GetCondition ("XFEMCoupling", xfemConditions);
    
    lengthSend[0] = 0;
    // pack data on all processors
    for(unsigned int i=0; i<xfemConditions.size(); i++)
    {
        map< int, RefCountPtr<DRT::Element > >  geometryMap = xfemConditions[i]->Geometry();
        map< int, RefCountPtr<DRT::Element > >::iterator iterGeo;   
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
    int                                     index, 
    vector<char>&                           cutterDataRecv, 
    vector<int>&                            nodeVectorRecv, 
    map< int, RefCountPtr<DRT::Node> >&     nodeMap )
{       
    int count = 0;
       
    while (index < (int) cutterDataRecv.size())
    {
        vector<char> data;
        DRT::ParObject::ExtractfromPack(index, cutterDataRecv, data);

        // allocate a node. Fill it with info from extracted element data
        DRT::ParObject* o = DRT::Utils::Factory(data);

        // cast ParObject to node
        DRT::Node* actNode = dynamic_cast< DRT::Node* >(o);
        RefCountPtr<DRT::Node> nodeRCP = rcp(actNode);
    
        // store nodes  
        nodeMap.insert(make_pair(nodeVectorRecv[count++], nodeRCP));
    }
}



/*----------------------------------------------------------------------*
 |  PARALLEL:   collects all cutter elements on              u.may 10/07|
 |              each processor   										|
 *----------------------------------------------------------------------*/
void Intersection::getCutterElementsInParallel(	
	vector< DRT::Condition * >&      		xfemConditions,
	vector<int>&							conditionEleCount, 
	map< int, vector< DRT::Element* > >&   	cutterElementMap,
	map< int, RefCountPtr<DRT::Node> >&  	cutterNodeMap,
	map<int, set<int> >& 					xfemCutterIdMap,
	const RefCountPtr<DRT::Discretization>& xfemdis,
   	const RefCountPtr<DRT::Discretization>& cutterdis  )
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
			map< int, RefCountPtr<DRT::Node> > nodeMap;
			
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
  				DRT::ParObject* o = DRT::Utils::Factory(data);

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
				
                int actCutterId = actCutter->Id() + cutterIdAdd;
  				count++;
  				
  				Epetra_SerialDenseMatrix    cutterXAABB = computeFastXAABB(actCutter);	
                
  				for(int k = 0; k < xfemdis->NumMyColElements(); k++)
    			{      
        			DRT::Element* xfemElement = xfemdis->lColElement(k);
       				initializeXFEM(xfemElement);    
  					Epetra_SerialDenseMatrix    xfemXAABB    = computeFastXAABB(xfemElement);                               
            		bool intersected = intersectionOfXAABB(cutterXAABB, xfemXAABB);    
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
           						int nodeId =  actCutter->Nodes()[inode]->Id();
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
 | GM:      transforms a node in reference coordinates       u.may 07/07|
 |          into current coordinates                                    |
 *----------------------------------------------------------------------*/  
void Intersection::referenceToCurrentCoordinates(   
    DRT::Element* element, 
    Epetra_SerialDenseVector& xsi)
{
    const int numNodes = element->NumNode();
    vector<int> actParams(1,0);
    Epetra_SerialDenseVector funct(numNodes);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;   
    params.set("action","calc_Shapefunction");
           
    actParams[0] = numNodes;            
    element->Evaluate(params, dummyDis, actParams, emptyM, emptyM , funct, xsi, emptyV);           
        
    for(int i=0; i<3; i++)
    {
        xsi[i] = 0.0;
        for(int j=0; j<numNodes; j++)
            xsi[i] += element->Nodes()[j]->X()[i] * funct(j);
    }
}



/*----------------------------------------------------------------------*
 | GM:      transforms a node in reference coordinates       u.may 07/07|
 |          into current coordinates                                    |
 *----------------------------------------------------------------------*/  
void Intersection::referenceToCurrentCoordinates(   
    DRT::Element*                               surfaceElement, 
    Epetra_SerialDenseVector&                   xsi,
    const vector<Epetra_SerialDenseVector>&     surfaceNodes)
{
    const int numNodes = surfaceElement->NumNode();
    vector<int> actParams(1,0);
    Epetra_SerialDenseVector funct(numNodes);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;   
    params.set("action","calc_Shapefunction");
           
    actParams[0] = numNodes;            
    surfaceElement->Evaluate(params, dummyDis, actParams, emptyM, emptyM , funct, xsi, emptyV);           
        
    for(int dim=0; dim<3; dim++)
    {
        xsi[dim] = 0.0;
        for(int i=0; i<numNodes; i++)
            xsi[dim] += surfaceNodes[i][dim] * funct(i);
    }
}



/*----------------------------------------------------------------------*
 | GM:      transforms a node in reference coordinates       u.may 07/07|
 |          into current coordinates                                    |
 *----------------------------------------------------------------------*/  
void Intersection::referenceToCurrentCoordinates(  
    Epetra_SerialDenseVector&                   xsi,
    const vector<Epetra_SerialDenseVector>&     plane)
{
    const int numNodes = 4;
    Epetra_SerialDenseVector funct(numNodes);
    
    shape_function_2D(funct, xsi[0], xsi[1], DRT::Element::quad4);
        
    for(int dim=0; dim<3; dim++)
    {
        xsi[dim] = 0.0;
        for(int i=0; i<numNodes; i++)
            xsi[dim] += plane[i][dim] * funct(i);
    }
}



/*----------------------------------------------------------------------*
 | GM:  transforms a node in current coordinates            u.may 07/07 |
 |      into reference coordinates                                      |
 *----------------------------------------------------------------------*/  
void Intersection::currentToReferenceCoordinates(   
    DRT::Element*               element, 
    Epetra_SerialDenseVector&   xsi)
{
    
    Epetra_SerialDenseVector x(3);
    
    for(int i = 0; i < 3; i++)
    {
        x[i] = xsi[i];
        xsi[i] = 0.0;
    }
    
    checkNodeWithinElement(element, x, xsi);
    
/*    if(!nodeWithinElement)
    {
        printf("xsi = %f \t %f \t %f\n", xsi[0], xsi[1], xsi[2]);
    }
    //    dserror("node not within element");
*/        
        
    // rounding 1 and -1 to be exact for the CDT
    for(int j = 0; j < 3; j++)
    {
        if( fabs((fabs(xsi[j])-1.0)) < TOL7_ &&  xsi[j] < 0)    xsi[j] = -1.0;
        if( fabs((fabs(xsi[j])-1.0)) < TOL7_ &&  xsi[j] > 0)    xsi[j] =  1.0;      
    }  
}     



/*----------------------------------------------------------------------*
 |  GM:     compares two nodes                               u.may 06/07|
 |          overloaded method:                                          |
 |          double*  and  double*                                       |
 *----------------------------------------------------------------------*/  
bool Intersection::comparePoints(    
    const double*     point1,
    const double*     point2,
    const int         length)
{   
    bool equal = true;
             
    for(int i = 0; i < length; i++)
        if(fabs(point1[i] - point2[i]) > TOL7_)
        {
            equal = false;
            break;
        }
  
    return equal;
}



/*----------------------------------------------------------------------*
 |  GM:     compares two nodes                               u.may 06/07|
 |          overloaded method:  vector<double>  and double*             |
 *----------------------------------------------------------------------*/  
bool Intersection::comparePoints(   
    const vector<double>&     point1,
    const double*             point2)
{   
    bool equal = true;
    
    if(point2 == NULL)
        dserror("array is NULL");
        
    for(unsigned int i = 0; i < point1.size() ; i++)
        if(fabs(point1[i] - point2[i]) > TOL7_)
        {
            equal = false;
            break;
        }
    
    return equal;
}



/*----------------------------------------------------------------------*
 |  GM:     compares two nodes                               u.may 06/07|
 |          overloaded method:  vector<double>  and  vector<double>     |
 *----------------------------------------------------------------------*/  
bool Intersection::comparePoints(   
    const vector<double>& point1,
    const vector<double>& point2)
{   
    bool equal = true;
    
    if(point1.size() != point2.size())
        dserror("arrays of nodes need to have the same length");
             
    for(unsigned int i = 0; i < point1.size() ; i++)
        if(fabs(point1[i] - point2[i]) > TOL7_)
        {
            equal = false;
            break;
        }
  
    return equal;
}



/*----------------------------------------------------------------------*
 |  GM:     compares two nodes                               u.may 06/07|
 |          overloaded method:  DRT::Node*  and  DRT::Node*             |
 *----------------------------------------------------------------------*/  
bool Intersection::comparePoints( 
    const DRT::Node*     point1,
    const DRT::Node*     point2)
{
    bool equal = true;
             
    for(unsigned int i = 0; i < 3 ; i++)
        if(fabs(point1->X()[i] - point2->X()[i]) > TOL7_)
        {
            equal = false;
            break;
        }
  
    return equal;
}



/*----------------------------------------------------------------------*
 |  GM:     compares two nodes                               u.may 06/07|
 |          overloaded method:                                          |
 |          Epetra_SerialDenseVector  and  Epetra_SerialDenseVector     |
 *----------------------------------------------------------------------*/  
bool Intersection::comparePoints(    
    const Epetra_SerialDenseVector&     point1,
    const Epetra_SerialDenseVector&     point2)
{   
    bool equal = true;
    
    if(point1.Length() != point2.Length())
        dserror("arrays of nodes need to have the same length");
             
    for(int i = 0; i < point1.Length() ; i++)
        if(fabs(point1[i] - point2[i]) > TOL7_)
        {
            equal = false;
            break;
        }
  
    return equal;
}



/*----------------------------------------------------------------------*
 |  GM:     checks if a certain element is a                 u.may 06/07|
 |          volume  element with help of the discretization type        |
 *----------------------------------------------------------------------*/  
bool Intersection::checkIfVolumeElement(
    DRT::Element* element)
{
    bool isVolume = false;
    DRT::Element::DiscretizationType distype = element->Shape();
    
    if( distype == DRT::Element::hex8  ||
        distype == DRT::Element::hex20 ||
        distype == DRT::Element::hex27 ||
        distype == DRT::Element::tet4  ||
        distype == DRT::Element::tet10  )
    {
        isVolume = true;        
    }
    return isVolume;
}



/*----------------------------------------------------------------------*
 |  GM:     checks if a certain element is a                 u.may 06/07|
 |          with help of the discretization type                        |
 *----------------------------------------------------------------------*/  
bool Intersection::checkIfSurfaceElement(
    DRT::Element* element)
{
    bool isSurface = false;
    DRT::Element::DiscretizationType distype = element->Shape();
    
    if( distype == DRT::Element::quad4 ||
        distype == DRT::Element::quad8 ||
        distype == DRT::Element::quad9 ||
        distype == DRT::Element::tri3  ||
        distype == DRT::Element::tri6  )
    {
        isSurface = true;       
    }
    return isSurface;
}



/*----------------------------------------------------------------------*
 |  GM:     checks if a certain element is a                 u.may 06/07|
 |          line element with help of the discretization type           |
 *----------------------------------------------------------------------*/  
bool Intersection::checkIfLineElement(
    DRT::Element* element)
{
    bool isLine = false;
    DRT::Element::DiscretizationType distype = element->Shape();
    
    if( distype == DRT::Element::line2 ||
        distype == DRT::Element::line3 )
    {
        isLine = true;      
    }
    return isLine;
}



/*----------------------------------------------------------------------*
 |  ICS:    computes an extended axis-aligned bounding box   u.may 06/07|
 |          XAABB for a given element                                   |
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix Intersection::computeFastXAABB( 
    DRT::Element* element)
{
    double  maxDistance; 
    Epetra_SerialDenseMatrix XAABB(3, 2);
    
	for(int dim=0; dim<3; dim++)
	{
		XAABB(dim, 0) = element->Nodes()[0]->X()[dim] - TOL7_;
   		XAABB(dim, 1) = element->Nodes()[0]->X()[dim] + TOL7_;
	}
    
    for(int i=1; i<element->NumNode(); i++)
        for(int dim=0; dim<3; dim++)
		{
            XAABB(dim, 0) = std::min( XAABB(dim, 0), element->Nodes()[i]->X()[dim] - TOL7_);
			XAABB(dim, 1) = std::max( XAABB(dim, 1), element->Nodes()[i]->X()[dim] + TOL7_);
		}
 
    maxDistance = fabs(XAABB(0,1) - XAABB(0,0));
 	for(int dim=1; dim<3; dim++)
	   maxDistance = std::max(maxDistance, fabs(XAABB(dim,1)-XAABB(dim,0)) );
	
    // subtracts half of the maximal distance to minX, minY, minZ
    // adds half of the maximal distance to maxX, maxY, maxZ 
	for(int dim=0; dim<3; dim++)
	{
		XAABB(dim, 0) = XAABB(dim, 0) - 0.5*maxDistance;
		XAABB(dim, 1) = XAABB(dim, 1) + 0.5*maxDistance;
	}	
	
	/*
    printf("\n");
	printf("axis-aligned bounding box:\n minX = %f\n minY = %f\n minZ = %f\n maxX = %f\n maxY = %f\n maxZ = %f\n", 
			  XAABB(0,0), XAABB(1,0), XAABB(2,0), XAABB(0,1), XAABB(1,1), XAABB(2,1));
	printf("\n");
	*/
    
	return XAABB;
}



/*----------------------------------------------------------------------*
 |  ICS:    checks if a node is within an XAABB               u.may 06/07|
 *----------------------------------------------------------------------*/
bool Intersection::isNodeWithinXAABB(    
    const std::vector<double>&         node,
    const Epetra_SerialDenseMatrix&    XAABB)
{
	bool isWithin = true;
    const double tol = TOL7_;
	
	for (int dim=0; dim<3; dim++)
	{
        double diffMin = node[dim] - XAABB(dim,0);
        double diffMax = XAABB(dim,1) - node[dim];
        
   	    if((diffMin < -tol)||(diffMax < -tol)) //check again !!!!!      
            isWithin = false;
    }
		 	
	return isWithin;
}



/*----------------------------------------------------------------------*
 |  ICS:    checks if two XAABB's intersect                  u.may 06/07|
 *----------------------------------------------------------------------*/
bool Intersection::intersectionOfXAABB(  
    const Epetra_SerialDenseMatrix&     cutterXAABB, 
    const Epetra_SerialDenseMatrix&     xfemXAABB)
{
	
  /*====================================================================*/
  /* bounding box topology*/
  /*--------------------------------------------------------------------*/
  /* parameter coordinates (x,y,z) of nodes
   * node 0: (minX, minY, minZ)
   * node 1: (maxX, minY, minZ)
   * node 2: (maxX, maxY, minZ)
   * node 3: (minX, maxY, minZ)
   * node 4: (minX, minY, maxZ)
   * node 5: (maxX, minY, maxZ)
   * node 6: (maxX, maxY, maxZ)
   * node 7: (minX, maxY, maxZ)
   * 
   *                      z
   *                      |           
   *             4========|================7
   *           //|        |               /||
   *          // |        |              //||
   *         //  |        |             // ||
   *        //   |        |            //  ||
   *       //    |        |           //   ||
   *      //     |        |          //    ||
   *     //      |        |         //     ||
   *     5=========================6       ||
   *    ||       |        |        ||      ||
   *    ||       |        o--------||---------y
   *    ||       |       /         ||      ||
   *    ||       0------/----------||------3
   *    ||      /      /           ||     //
   *    ||     /      /            ||    //
   *    ||    /      /             ||   //
   *    ||   /      /              ||  //
   *    ||  /      /               || //
   *    || /      x                ||//
   *    ||/                        ||/
   *     1=========================2
   *
   */
  /*====================================================================*/
	
    bool intersection =  false;
	std::vector < std::vector <double> > nodes(8, vector<double> (3, 0.0));
	
	nodes[0][0] = cutterXAABB(0,0);	nodes[0][1] = cutterXAABB(1,0);	nodes[0][2] = cutterXAABB(2,0);	// node 0	
	nodes[1][0] = cutterXAABB(0,1);	nodes[1][1] = cutterXAABB(1,0);	nodes[1][2] = cutterXAABB(2,0);	// node 1
	nodes[2][0] = cutterXAABB(0,1);	nodes[2][1] = cutterXAABB(1,1);	nodes[2][2] = cutterXAABB(2,0);	// node 2
	nodes[3][0] = cutterXAABB(0,0);	nodes[3][1] = cutterXAABB(1,1);	nodes[3][2] = cutterXAABB(2,0);	// node 3
	nodes[4][0] = cutterXAABB(0,0);	nodes[4][1] = cutterXAABB(1,0);	nodes[4][2] = cutterXAABB(2,1);	// node 4
	nodes[5][0] = cutterXAABB(0,1);	nodes[5][1] = cutterXAABB(1,0);	nodes[5][2] = cutterXAABB(2,1);	// node 5
	nodes[6][0] = cutterXAABB(0,1);	nodes[6][1] = cutterXAABB(1,1);	nodes[6][2] = cutterXAABB(2,1);	// node 6
	nodes[7][0] = cutterXAABB(0,0);	nodes[7][1] = cutterXAABB(1,1);	nodes[7][2] = cutterXAABB(2,1);	// node 7
	
	for (int i = 0; i < 8; i++)
		if(isNodeWithinXAABB(nodes[i], xfemXAABB))
		{
			intersection = true;
			break;
		}
	
		
	if(!intersection)
	{
		nodes[0][0] = xfemXAABB(0,0);	nodes[0][1] = xfemXAABB(1,0);	nodes[0][2] = xfemXAABB(2,0);	// node 0	
		nodes[1][0] = xfemXAABB(0,1);	nodes[1][1] = xfemXAABB(1,0);	nodes[1][2] = xfemXAABB(2,0);	// node 1
		nodes[2][0] = xfemXAABB(0,1);	nodes[2][1] = xfemXAABB(1,1);	nodes[2][2] = xfemXAABB(2,0);	// node 2
		nodes[3][0] = xfemXAABB(0,0);	nodes[3][1] = xfemXAABB(1,1);	nodes[3][2] = xfemXAABB(2,0);	// node 3
		nodes[4][0] = xfemXAABB(0,0);	nodes[4][1] = xfemXAABB(1,0);	nodes[4][2] = xfemXAABB(2,1);	// node 4
		nodes[5][0] = xfemXAABB(0,1);	nodes[5][1] = xfemXAABB(1,0);	nodes[5][2] = xfemXAABB(2,1);	// node 5
		nodes[6][0] = xfemXAABB(0,1);	nodes[6][1] = xfemXAABB(1,1);	nodes[6][2] = xfemXAABB(2,1);	// node 6
		nodes[7][0] = xfemXAABB(0,0);	nodes[7][1] = xfemXAABB(1,1);	nodes[7][2] = xfemXAABB(2,1);	// node 7
	
		for (int i = 0; i < 8; i++)
			if(isNodeWithinXAABB(nodes[i], cutterXAABB))
			{
				intersection = true;
				break;
			}
	}	
	return intersection;
}



/*----------------------------------------------------------------------*
 |  CLI:    collects points that belong to the interface     u.may 06/07|
 |          and lie within an xfem element                              |
 *----------------------------------------------------------------------*/
bool Intersection::collectInternalPoints(   
    DRT::Element*                   xfemElement,
    DRT::Element*                   cutterElement,
    DRT::Node*                      node,
    std::vector< InterfacePoint >&  interfacePoints,
    const int                       elemId,
    const int                       nodeId)
{
    Epetra_SerialDenseVector xsi(3);
    Epetra_SerialDenseVector x(3);
       
    x[0] = node->X()[0];
    x[1] = node->X()[1];
    x[2] = node->X()[2];
    
    bool nodeWithinElement = checkNodeWithinElement(xfemElement, x, xsi);
    // debugNodeWithinElement(xfemElement,node, xsi, elemId ,nodeId, nodeWithinElement);  
    
    if(nodeWithinElement)
    {   
        InterfacePoint ip;
        //debugNodeWithinElement(xfemElement,node,xsi,elemId ,nodeId, nodeWithinElement);  
          
        numInternalPoints_++;
        
        // check if node lies on the boundary of the xfem element
        if(checkIfOnBoundary(xsi, ip))  numBoundaryPoints_++;
                                   
        // intersection coordinates in the surface 
        // element reference coordinate system
        getNodeCoordinates(nodeId, ip.coord, cutterElement->Shape());
                          
        interfacePoints.push_back(ip);
        
        storeIntersectedCutterElement(cutterElement); 
    }
      
    return nodeWithinElement;
}



/*----------------------------------------------------------------------*
 |  CLI:    updates the Jacobi matrix for the computation    u.may 06/07|
 |          if a node is in a given element                             |
 *----------------------------------------------------------------------*/
void Intersection::updateAForNWE(   
    const int                   dim,
    Epetra_SerialDenseMatrix&   A,
    Epetra_SerialDenseVector&   xsi,
    DRT::Element*               element)                                                  
{	
    const int numNodes = element->NumNode();
    vector<int> actParams(1,0);
    Epetra_SerialDenseMatrix deriv1(dim, numNodes);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
   
    
    params.set("action","calc_ShapeDeriv1");
    actParams[0] = numNodes;  
       
    element->Evaluate(params, dummyDis, actParams, deriv1, emptyM , emptyV, xsi, emptyV);       
       
    for(int i=0; i<dim; i++)
        for(int j=0; j<dim; j++)
            A[i][j] = 0.0;
          
    for(int i=0; i<dim; i++)
        for(int k=0; k<dim; k++)
            for(int j=0; j<numNodes; j++)
                A[i][k] += element->Nodes()[j]->X()[i] * deriv1[j][k];
        
}



/*----------------------------------------------------------------------*
 |  CLI:    updates the rhs for the computation if a         u.may 06/07|
 |          node is in a given element                                  |
 *----------------------------------------------------------------------*/
void Intersection::updateRHSForNWE( 
    const int                           dim,
    Epetra_SerialDenseVector&           b,
    Epetra_SerialDenseVector&           xsi,
    const Epetra_SerialDenseVector&     x,
    DRT::Element*                       element)                                                  
{
    const int numNodes = element->NumNode();
    vector<int> actParams(1,0);
    Epetra_SerialDenseVector funct(numNodes);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
      
    params.set("action","calc_Shapefunction");
    actParams[0] = numNodes;       
      
    element->Evaluate(params, dummyDis, actParams, emptyM, emptyM , funct, xsi, emptyV);        
    
    for(int i=0; i<dim; i++)
        b[i] = 0.0;
    
    for(int i=0; i<dim; i++)
        for(int j=0; j<numNodes; j++)
            b[i] += (-1.0) * element->Nodes()[j]->X()[i] * funct[j];
        
     for(int i=0; i<dim; i++)
        b[i] += x[i];
}



/*----------------------------------------------------------------------*
 |  CLI:    checks if a node is within a given element       u.may 06/07|	
 *----------------------------------------------------------------------*/
bool Intersection::checkNodeWithinElement(	
    DRT::Element*                       element,
    const Epetra_SerialDenseVector&     x,
    Epetra_SerialDenseVector&           xsi)
{

	bool nodeWithinElement = true;
    int iter = 0;
    int dim = getDimension(element);
    const int maxiter = 50;
    double residual = 1.0;
    
    Epetra_SerialDenseMatrix A(dim,dim);
    Epetra_SerialDenseVector b(dim);
    
    Epetra_SerialDenseMatrix A_old(dim,dim);
    Epetra_SerialDenseVector b_old(dim);
    Epetra_SerialDenseVector dx(dim);
    
    for(int i = 0; i<dim; i++)
        xsi[i] = 0.0;

    dx = xsi;
            
    updateRHSForNWE( dim, b, xsi, x, element);
   
    while(residual > TOL14_)
    { 
        updateAForNWE( dim, A, xsi, element);
        A_old = A;
        b_old = b;
        //!gaussElimination(A, b, dx, true, dim, 1)
        //!solveLinearSystemWithSVD(A, b, dx, dim)
        if(!gaussElimination(A, b, dx, true, dim, 1))
        {
            printf("MATRIX SINGULAR\n");
            nodeWithinElement = false;
            break;
        }   
        
        Epetra_SerialDenseVector b_new(dim);
        
        /*for(int i = 0 ; i < dim ; i++)
        {
            b_new[i] = 0.0;
            for(int j = 0 ; j < dim ; j++)
                b_new[i] += A_old[i][j]*dx[j];
        }
        printf("x = %20.16f\t, x = %20.16f\t, x = %20.16f\n", x[0],x[1],x[2]);
        printf("dx = %20.16f\t, dx = %20.16f\t, dx = %20.16f\n", dx[0],dx[1],dx[2]);
        printf("b = %20.16f\t, b = %20.16f\t, b = %20.16f\n", b_old[0],b_old[1],b_old[2]);
        printf("b1 = %20.16f\t, b1 = %20.16f\t, b1 = %20.16f\n", b_new[0],b_new[1],b_new[2]);
        printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\n", xsi[0],xsi[1],xsi[2]);
        
        printf("\n");
        */
        xsi = addTwoVectors(xsi,dx);
        
        if(iter >= maxiter)
        {   
            nodeWithinElement = false;
            break;
        }   
        
            
        updateRHSForNWE(dim, b, xsi, x, element);
        residual = b.Norm2();
        iter++; 
    }
    
    //printf("iter = %d\n", iter);
    //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi[0],xsi[1],xsi[2], residual, TOL14_);
    
    for(int i=0; i<dim; i++)
        if( (fabs(xsi[i])-1.0) > TOL7_)     
        {    
            nodeWithinElement = false;
            break;
        }
        
    return nodeWithinElement;
}



/*----------------------------------------------------------------------*
 |  CLI:    checks if a node that lies within an element     u.may 06/07|
 |          lies on one of its surfaces or nodes                        |                                           
 *----------------------------------------------------------------------*/
bool Intersection::checkIfOnBoundary( 
    Epetra_SerialDenseVector&       xsi,    
    InterfacePoint&                 ip)
{
    bool onSurface = false;
    const int count = getSurfaces(xsi, ip.surfaces, xfemDistype_);
        
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
    DRT::Element*                   surfaceElement,
    DRT::Element*                   lineElement,
    std::vector<InterfacePoint>&    interfacePoints, 
    const int                       surfaceId,
    const int                       lineId,
    const bool                      lines,
    bool&                           xfemIntersection)
{
    Epetra_SerialDenseVector xsi(3);
    Epetra_SerialDenseVector upLimit(3);
    Epetra_SerialDenseVector loLimit(3);
   
    if(!checkIfSurfaceElement(surfaceElement))
        dserror("surface element has to be a surface element\n");
    if(!checkIfLineElement(lineElement))
        dserror("line element has to be a line element\n");
  
  
    for(int i = 0; i < 3; i++)
    {
       xsi[i] =  0.0; 
       upLimit[i]  =  1.0;  
       loLimit[i]  = -1.0;
       // extend for triangle
    }
    
    bool intersected = computeCurveSurfaceIntersection(surfaceElement, lineElement, xsi, upLimit, loLimit);
                                        
    if(intersected) 
        addIntersectionPoint( surfaceElement, lineElement,xsi, upLimit, loLimit, 
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
    DRT::Element*               surfaceElement,
    DRT::Element*               lineElement,
    Epetra_SerialDenseVector&   xsi,
    Epetra_SerialDenseVector&   upLimit,
    Epetra_SerialDenseVector&   loLimit)
{
    bool intersection = true;
    int iter = 0;
    const int maxiter = 50;
    double residual = 1.0;
    Epetra_SerialDenseMatrix A(3,3);
    Epetra_SerialDenseVector b(3);
    Epetra_SerialDenseVector dx(3);
   
    Epetra_SerialDenseMatrix A_old(3,3);
    Epetra_SerialDenseVector b_old(3);
    Epetra_SerialDenseVector b_gauss(3);
    dx = xsi;
 
    updateRHSForCSI( b, xsi, surfaceElement, lineElement);
              
    while(residual > TOL14_)
    {   
        updateAForCSI( A, xsi, surfaceElement, lineElement);
   
        A_old = A;
        b_old = b;
        
      /*if(!gaussElimination(A, b, dx, true, 3, 1))
        {
            intersection = false;
            break;
        } 
        */
        if(!solveLinearSystemWithSVD(A, b, dx, 3))
        {
            intersection = false;
            break;
        } 
   
      /*printf("\n");  
        printf("Intersection\n");  
        printf("=============================================================================\n");  
        printf("\n");  
        printf("dx = %20.16f\t, dx = %20.16f\t, dx = %20.16f\n", dx[0], dx[1], dx[2]);
        
       for(int i = 0; i < 3; i++ )
       {
            b_gauss[i] = 0.0;
            for(int  k = 0; k < 3; k++ )
            {  
                b_gauss[i] += A_old[i][k]*dx[k];
            }
            printf("b = %f\t  b_gauss = %f\n", b_old[i], b_gauss[i]);
        }
        printf("\n");
        printf("\n");  
        printf("=============================================================================\n");  
        printf("\n");  
        */
        
        xsi = addTwoVectors(xsi,dx);
       
        if(iter >= maxiter)
        {   
            intersection = false;
            break;
        }       
       
        updateRHSForCSI( b, xsi, surfaceElement, lineElement);
        residual = b.Norm2(); 
        iter++;
    } 
    
    
    if( (xsi[0] > upLimit[0]+TOL7_) || (xsi[1] > upLimit[1]+TOL7_) || (xsi[2] > upLimit[2]+TOL7_)  || 
        (xsi[0] < loLimit[0]-TOL7_) || (xsi[1] < loLimit[1]-TOL7_) || (xsi[2] < loLimit[2]-TOL7_)) 
            intersection = false;
        
        
    return intersection;
}



/*----------------------------------------------------------------------*
 |  CLI:    updates the systemmatrix for the computation     u.may 06/07| 
 |          of a curve-surface intersection (CSI)                       |
 *----------------------------------------------------------------------*/
void Intersection::updateAForCSI(  	Epetra_SerialDenseMatrix&   A,
                                    Epetra_SerialDenseVector&   xsi,
                                    DRT::Element*               surfaceElement,
                                    DRT::Element*               lineElement)        											
{	
	const int numNodesSurface = surfaceElement->NumNode();
   	const int numNodesLine = lineElement->NumNode();
	vector<int> actParams(1,0);
	Epetra_SerialDenseMatrix surfaceDeriv1(2,numNodesSurface);
	Epetra_SerialDenseMatrix lineDeriv1(1,numNodesLine);
	Epetra_SerialDenseMatrix emptyM;
	Epetra_SerialDenseVector emptyV;
    Epetra_SerialDenseVector xsiSurface(2);
    Epetra_SerialDenseVector xsiLine(1);
	DRT::Discretization dummyDis("dummy discretization", null);
	ParameterList params;
	 
    xsiSurface(0) = xsi[0]; // r-coordinate surface
    xsiSurface(1) = xsi[1]; // s-coordinate surface
    xsiLine(0) = xsi[2];    // r-coordinate line
	actParams[0] = numNodesSurface;     
	
    for(unsigned int i=0; i<3; i++)
        for(unsigned int j=0; j<3; j++)
            A[i][j] = 0.0;
	
	params.set("action","calc_ShapeDeriv1");
	surfaceElement->Evaluate(params, dummyDis, actParams, surfaceDeriv1, emptyM , emptyV, xsiSurface, emptyV);			
    //printf("num of nodes surfaces = %d\n", numOfNodesSurface);
       	
    for(int dim=0; dim<3; dim++)
   	    for(int i=0; i<numNodesSurface; i++)
		{
			A[dim][0] += surfaceElement->Nodes()[i]->X()[dim] * surfaceDeriv1(0,i);
			A[dim][1] += surfaceElement->Nodes()[i]->X()[dim] * surfaceDeriv1(1,i);
		}
		
   
    actParams[0] = numNodesLine;
    lineElement->Evaluate(params, dummyDis, actParams, lineDeriv1, emptyM , emptyV, xsiLine, emptyV);

    for(int dim=0; dim<3; dim++)
        for(int i=0; i<numNodesLine; i++)
        {
            A[dim][2] +=  (-1) * lineElement->Nodes()[i]->X()[dim] * lineDeriv1(0,i);
        }
        
}



/*----------------------------------------------------------------------*
 |  CLI:    updates the right-hand-side for the              u.may 06/07|
 |          computation of                                              |
 |          a curve-surface intersection (CSI)                          |
 *----------------------------------------------------------------------*/
void Intersection::updateRHSForCSI( 
    Epetra_SerialDenseVector&   b,
    Epetra_SerialDenseVector&   xsi,
    DRT::Element*               surfaceElement,
    DRT::Element*               lineElement)        											
{
    int numNodesSurface = surfaceElement->NumNode();
    int numNodesLine = lineElement->NumNode();
    vector<int> actParams(1,0);
    Epetra_SerialDenseVector surfaceFunct(numNodesSurface);
    Epetra_SerialDenseVector lineFunct(numNodesLine);
    Epetra_SerialDenseMatrix emptyML(1,numNodesLine);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    Epetra_SerialDenseVector xsiSurface(2);
    Epetra_SerialDenseVector xsiLine(1);
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
    
    
    xsiSurface(0) = xsi[0]; // r-coordinate surface
    xsiSurface(1) = xsi[1]; // s-coordinate surface
    xsiLine(0) = xsi[2];    // r-coordinate line
    params.set("action","calc_Shapefunction");
    		
   	for(unsigned int i=0; i<3; i++)   
		b[i] = 0.0;
    
   	actParams[0] = numNodesSurface;            
	surfaceElement->Evaluate(params, dummyDis, actParams, emptyM, emptyM , surfaceFunct, xsiSurface, emptyV);			
  
	actParams[0] = numNodesLine;     
	lineElement->Evaluate(params, dummyDis, actParams, emptyM, emptyML, lineFunct, xsiLine, emptyV);
		
   	for(int dim=0; dim<3; dim++)
   		for(int i=0; i<numNodesSurface; i++)
		{
			b[dim] += (-1) * surfaceElement->Nodes()[i]->X()[dim] * surfaceFunct(i);
		}
	
	for(int dim=0; dim<3; dim++)
   		for(int i=0; i<numNodesLine; i++)
		{
			b[dim] += lineElement->Nodes()[i]->X()[dim] * lineFunct(i);
		}
}



/*----------------------------------------------------------------------*
 |  CLI:    computes a new starting point for the            u.may 06/07| 
 |          Newton-method in order to find all intersection points      | 
 |          of a curve-surface intersection                             |
 *----------------------------------------------------------------------*/  
int Intersection::computeNewStartingPoint(
    DRT::Element*                surfaceElement,
    DRT::Element*                lineElement,
    const int                    surfaceId,
    const int                    lineId,
    Epetra_SerialDenseVector&    upLimit,
    Epetra_SerialDenseVector&    loLimit,
    std::vector<InterfacePoint>& interfacePoints,
    const bool                   lines)
{	
    bool interval = true;
	int numInterfacePoints = 0;
	Epetra_SerialDenseVector xsi(3);
	
    //printf("xsi2 = %f   %f   %f\n", fabs(xsi[0]), fabs(xsi[1]), fabs(xsi[2]) ); 
    //printf("lolimit = %f   %f   %f\n", fabs(loLimit[0]), fabs(loLimit[1]), fabs(loLimit[2]) ); 
    //printf("uplimit = %f   %f   %f\n", fabs(upLimit[0]), fabs(upLimit[1]), fabs(upLimit[2]) ); 
    
    if(comparePoints(upLimit, loLimit))
        interval = false;
        
	for(int i = 0; i < 3; i++)
		xsi[i] = (double) (( upLimit[i] + loLimit[i] )/2.0);
       
 
	bool intersected = computeCurveSurfaceIntersection(surfaceElement, lineElement, xsi, upLimit, loLimit);
   
    
    if( comparePoints(xsi,upLimit))     intersected = false;
    if( comparePoints(xsi,loLimit))     intersected = false;
       							 	
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
    DRT::Element*                	surfaceElement,
    DRT::Element*                	lineElement,
    Epetra_SerialDenseVector&		xsi,
    Epetra_SerialDenseVector&     	upLimit,
    Epetra_SerialDenseVector&     	loLimit,
    std::vector<InterfacePoint>& 	interfacePoints,
    const int                       surfaceId,
    const int                       lineId,							  
    const bool 					    lines)
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
        //printf("alreadyinlist = false\n");
        interfacePoints.push_back(ip);  
        numInterfacePoints++;
    }   
    
    // recursive call 8 times !!!!
    numInterfacePoints =  numInterfacePoints +
      computeNewStartingPoint(	surfaceElement, lineElement, surfaceId, lineId, 
      							upLimit, xsi, interfacePoints, lines) +
      computeNewStartingPoint(	surfaceElement, lineElement, surfaceId, lineId, 
      							xsi, loLimit, interfacePoints, lines);

	return numInterfacePoints;
}



/*----------------------------------------------------------------------*
 |  ICS:    computes the convex hull of a set of             u.may 06/07|
 |          interface points and stores resulting points,               |
 |          segments and triangles for the use with Tetgen (CDT)        |
 *----------------------------------------------------------------------*/
#ifdef QHULL
void Intersection::computeConvexHull(   
    DRT::Element*           xfemElement,
    DRT::Element*           surfaceElement,
    vector<InterfacePoint>& interfacePoints)
{
    
    vector<int>                 positions;
    vector<double>              searchPoint(3,0);
    vector<double>              vertex(3,0);
    vector< vector<double> >    vertices;  
    Epetra_SerialDenseVector    curCoord(3);  
    InterfacePoint              midpoint;
    
           
    if(!checkIfSurfaceElement(surfaceElement))
   		dserror("surface element has to be a surface element\n");
          
    if(interfacePoints.size() > 2)  
    {
        midpoint = computeMidpoint(interfacePoints);
        // transform it into current coordinates
        for(int j = 0; j < 2; j++)      curCoord[j]  = midpoint.coord[j];            
        referenceToCurrentCoordinates(surfaceElement, curCoord);    
        currentToReferenceCoordinates(xfemElement, curCoord);    
        for(int j = 0; j < 3; j++)      midpoint.coord[j] = curCoord[j]; 
     
        // store coordinates in 
        // points has numInterfacePoints*dim-dimensional components
        // points[0] is the first coordinate of the first point
        // points[1] is the second coordinate of the first point
        // points[dim] is the first coordinate of the second point                     
        coordT* coordinates = (coordT *)malloc((2*interfacePoints.size())*sizeof(coordT));
        int fill = 0;
        for(unsigned int i = 0; i < interfacePoints.size(); i++)
        {
            for(int j = 0; j < 2; j++)
            {
                coordinates[fill++] = interfacePoints[i].coord[j]; 
               // printf("coord = %f\t", interfacePoints[i].coord[j]);
            }
            // printf("\n");
            // transform interface points into current coordinates
            for(int j = 0; j < 2; j++)      
                curCoord[j]  = interfacePoints[i].coord[j];   
                          
            referenceToCurrentCoordinates(surfaceElement, curCoord);  
            currentToReferenceCoordinates(xfemElement, curCoord);
                     
            for(int j = 0; j < 3; j++)         
                interfacePoints[i].coord[j] = curCoord[j];
                
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
                                
                for(int m = 0; m < 2; m++)      curCoord[m]  = vertex[m];           
                // surface reference coordinates to current coordinates       
                referenceToCurrentCoordinates(surfaceElement, curCoord); 
                // current coordinates to xfem element reference coordinates
                currentToReferenceCoordinates(xfemElement, curCoord);    
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
        for(unsigned int i = 0; i < interfacePoints.size(); i++)
        {
            // transform interface points into current coordinates
            for(int j = 0; j < 2; j++)         
                curCoord[j]  = interfacePoints[i].coord[j];   
            
            // surface reference coordinates to current coordinates       
            referenceToCurrentCoordinates(surfaceElement, curCoord); 
            // current coordinates to xfem element reference coordinates
            currentToReferenceCoordinates(xfemElement, curCoord);   
              
            for(int j = 0; j < 3; j++)      
            {   
                interfacePoints[i].coord[j] = curCoord[j];
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
    vector<double>&             searchPoint)
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
 |  for an intersected xfem element in reference configuration          |
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
    DRT::Element*               			element,
    DRT::Element*               			cutterElement,
    map< int, DomainIntCells >&	domainintcells)
{
    int dim = 3;
    int nsegments = 0; 
    int nsurfPoints = 0;
    tetgenio in;
    tetgenio out;
    char switches[] = "pnno2Q";    //o2 Y
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
	bool regions = false;	
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
            in.regionlist[fill++] = 0.9;
            in.regionlist[fill++] = 0.9;
            in.regionlist[fill++] = 0.9;
	    		
	    	// store regional attribute (switch A)   i=0 cutter i=1 fluid
	    	in.regionlist[fill++] = 0;
	    	
	    	// store volume constraint (switch a)   
	    	in.regionlist[fill++] = 0.0;
            
            in.regionlist[fill++] = -1.0;
            in.regionlist[fill++] = -1.0;
            in.regionlist[fill++] = -1.0;
                
            // store regional attribute (switch A)   i=0 cutter i=1 fluid
            in.regionlist[fill++] = 1;
            
            // store volume constraint (switch a)   
            in.regionlist[fill++] = 0.0;
		//} 
	}
	   
    in.save_nodes("tetin");   
    in.save_poly("tetin");   
    //  Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
    //  do quality mesh generation (q) with a specified quality bound
    //  (1.414), and apply a maximum volume constraint (a0.1)
    tetrahedralize(switches, &in, &out); 
 
    
    //Debug
    vector<int> elementIds;
    for(int i = 0; i<numXFEMCornerNodes_; i++)
        elementIds.push_back(i);
    
    debugTetgenOutput(in, out, element, elementIds);
    //printTetViewOutputPLC( element, element->Id(), in);
    
    // recover curved interface for higher order meshes
    bool curvedInterface = true;
    if(curvedInterface)
        recoverCurvedInterface(element, out);
 
    //printTetViewOutput(element->Id(), out);
   
    
   
    // store integrationcells
    DomainIntCells listperElement;
    
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
        listperElement.push_back(DomainIntCell(tetrahedronCoord));                 
    }
    domainintcells.insert(make_pair(element->Id(),listperElement));
}



/*----------------------------------------------------------------------*
 |  CDT:    fills the point list with the corner points      u.may 06/07|
 |          in reference coordinates of the xfem element                |
 *----------------------------------------------------------------------*/  
void Intersection::startPointList(
    )
{
    InterfacePoint ip;
        
    if(!pointList_.empty())
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
    vector<InterfacePoint>::iterator it;
        

    for(unsigned int i = 0; i < interfacePoints.size(); i++ )
    {   
        if(comparePoints(point, interfacePoints[i].coord))
        {         
            alreadyInList = false;
            int count = 0;
            for(it = pointList_.begin(); it != pointList_.end(); it++ )  
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
                pointList_.push_back(interfacePoints[i]);
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
    const vector<InterfacePoint>& interfacePoints)
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
    int pos1 = 0;
    int pos2 = 0;
    
    for(unsigned int i = 0; i < positions.size(); i++ )
    {
        pos1 = positions[i];
        
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
  
    for(unsigned int i = 0; i < intersectingCutterElements_.size(); i++)
        if(intersectingCutterElements_[i] == surfaceElement)
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
    DRT::Element*   xfemElement, 
    tetgenio&       out
    )
{
    vector<int>                                     order(3,0);
    vector<int>                                     tetraCornerIndices(4,0); 
    vector < Epetra_SerialDenseVector >             tetraCornerNodes(4, Epetra_SerialDenseVector(3));
    
    
    // list of point markers , if already visited = 1 , if not = 0
    int* visitedPointIndexList = new int[out.numberofpoints];      
    for(int i = 0; i<out.numberofpoints; i++)
        visitedPointIndexList[i] = 0;
        
    // lifts all corner points into the curved interface
    liftAllSteinerPoints(xfemElement, out);
      
    for(int i=0; i<out.numberoftrifaces; i++)
    {
        // run over all faces not lying in on of the xfem element planes
        int faceMarker = out.trifacemarkerlist[i] - facetMarkerOffset_;
        
        if(faceMarker > -1)
        {        
            int tetIndex = out.adjtetlist[i*2];
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
            }
        }
    }
    
    delete [] visitedPointIndexList;
    visitedPointIndexList = (int *) NULL;
    intersectingCutterElements_.clear();
}



/*----------------------------------------------------------------------*
 |  RCI:    checks if all tetrahedra corner points are lying u.may 09/07|
 |          in a surface element                                        |   
 |          if not corner points is recovered on the surface element    | 
 *----------------------------------------------------------------------*/
void Intersection::liftAllSteinerPoints(
    DRT::Element*                                   xfemElement,
    tetgenio&                                       out)
{
    int lineIndex, cutterIndex;
    Epetra_SerialDenseVector edgePoint(3);
    Epetra_SerialDenseVector oppositePoint(3);
    vector< vector <int> > adjacentFacesList;
    vector< vector <int> > adjacentFacemarkerList;
    
    locateSteinerPoints(adjacentFacesList, adjacentFacemarkerList, out);
    
    if(!adjacentFacesList.empty())
    {
        for(unsigned int i=0; i<adjacentFacesList.size(); i++)
        {        
            bool  normalSteiner = decideNormalOrPlane(  i, lineIndex, cutterIndex, adjacentFacesList, adjacentFacemarkerList,
                                                        edgePoint, oppositePoint, xfemElement, out);
              
            if(normalSteiner)
                liftSteinerPointOnSurface(i, adjacentFacesList, adjacentFacemarkerList, xfemElement, out);
            else
                liftSteinerPointOnEdge( i, lineIndex, cutterIndex, edgePoint, oppositePoint, 
                                        adjacentFacesList, xfemElement, out);
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
    tetgenio&                   out)
{

    for(int i = 0; i < out.numberoftrifaces; i++)
    {   
        if( (out.trifacemarkerlist[i]-facetMarkerOffset_) > -1)
        {  
            for (int j = 0; j < 3; j++)
            {
                int pointIndex = out.trifacelist[i*3+j];
                
                // check if point is a Steiner point
                if(out.pointmarkerlist[pointIndex] != 2 && out.pointmarkerlist[pointIndex] != 3 )
                {
                    bool alreadyInList = false;
                    vector<int> pointIndices = getPointIndices(out, i, j);
                    
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
bool Intersection::decideNormalOrPlane( 
    int                         steinerIndex,
    int&                        lineIndex, 
    int&                        cutterIndex,  
    vector< vector <int> >&     adjacentFacesList,
    vector< vector <int> >&     adjacentFacemarkerList,
    Epetra_SerialDenseVector&   edgePoint,
    Epetra_SerialDenseVector&   oppositePoint,
    DRT::Element*               xfemElement, 
    tetgenio&                   out)
{

    bool                        normalSteiner = true;
    int                         pointIndex = adjacentFacesList[steinerIndex][0];
    Epetra_SerialDenseVector    x(3);
    Epetra_SerialDenseVector    xsi(3);
        
    for(int k=0; k<3; k++)
    {
        x[k]   = out.pointlist[pointIndex*3 + k];
        xsi[k] = 0.0;
    }
   
    InterfacePoint emptyIp;
    checkNodeWithinElement(xfemElement, x, xsi);
    if(checkIfOnBoundary(xsi, emptyIp))     out.pointmarkerlist[pointIndex] = 3;    // on xfem boundary
    else                                    out.pointmarkerlist[pointIndex] = 2;    // not on xfem boundary
    
   
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
    return normalSteiner;
}



/*----------------------------------------------------------------------*
 |  RCI:    lifts Steiner points lying                       u.may 11/07|
 |          within a cutter element                                     |   
 *----------------------------------------------------------------------*/
void Intersection::liftSteinerPointOnSurface(
    int                         steinerIndex,
    vector< vector <int> >&     adjacentFacesList,
    vector< vector <int> >&     adjacentFacemarkerList,
    DRT::Element*               xfemElement, 
    tetgenio&                   out)
{
   
    Epetra_SerialDenseVector    Steinerpoint(3);
    Epetra_SerialDenseVector    averageNormal(3);
    vector<Epetra_SerialDenseVector> normals; 
    
     // get Steiner point coordinates 
    for(int j=0; j<3; j++)
    {
        averageNormal[j] = 0.0;
        Steinerpoint[j] = out.pointlist[adjacentFacesList[steinerIndex][0]*3 + j];
    }
    referenceToCurrentCoordinates(xfemElement, Steinerpoint);   
   
    double length = (int) ( ( (double) (adjacentFacesList[steinerIndex].size()-1))*0.5 );
    
    for(int j = 0; j < length; j++)
    {
        int pointIndex1 = adjacentFacesList[steinerIndex][1 + 2*j];
        int pointIndex2 = adjacentFacesList[steinerIndex][1 + 2*j + 1];
        Epetra_SerialDenseVector    p1(3);
        Epetra_SerialDenseVector    p2(3);
     
        for(int k = 0; k < 3; k++)
        {
            p1[k] =  out.pointlist[pointIndex1*3 + k];
            p2[k] =  out.pointlist[pointIndex2*3 + k];
        }
    
        referenceToCurrentCoordinates(xfemElement, p1);   
        referenceToCurrentCoordinates(xfemElement, p2);   
        Epetra_SerialDenseVector n1 = subtractsTwoVectors(p1, Steinerpoint);
        Epetra_SerialDenseVector n2 = subtractsTwoVectors(p2, Steinerpoint);
    
        Epetra_SerialDenseVector normal = computeCrossProduct( n1, n2);
        normal = normalizeVector(normal);
    
        for(int k=0; k<3; k++)
            averageNormal[k] += normal[k];
            
         normals.push_back(normal);
    }

    // compute average normal
    for(int j=0; j<3; j++)
        averageNormal[j] = averageNormal[j]/( (double)length);
           
        
    vector<Epetra_SerialDenseVector> plane(4);                      
    plane[0] = addTwoVectors(Steinerpoint, averageNormal);               
    plane[1] = subtractsTwoVectors(Steinerpoint, averageNormal);
 
 
    int faceMarker = adjacentFacemarkerList[steinerIndex][0];
    Epetra_SerialDenseVector xsi(3);
    
    bool intersected = computeRecoveryNormal( xsi, plane, intersectingCutterElements_[faceMarker],false);
    if(intersected)
    {
        storeHigherOrderNode(   true, adjacentFacesList[steinerIndex][0], -1, xsi,
                                intersectingCutterElements_[faceMarker], xfemElement, out);
    }
    else
    {
        // loop over all individual normals
        for(unsigned int j = 0; j < normals.size(); j++ )
        {
            plane[0] = addTwoVectors(Steinerpoint, normals[j]);               
            plane[1] = subtractsTwoVectors(Steinerpoint, normals[j]);
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
    int                         steinerIndex,
    int                         lineIndex,
    int                         cutterIndex,   
    Epetra_SerialDenseVector&   edgePoint,
    Epetra_SerialDenseVector&   oppositePoint,
    vector< vector <int> >&     adjacentFacesList,
    DRT::Element*               xfemElement, 
    tetgenio&                   out)
{

    Epetra_SerialDenseVector    Steinerpoint(3);
    
     // get Steiner point coordinates 
    for(int j=0; j<3; j++)
        Steinerpoint[j] = out.pointlist[adjacentFacesList[steinerIndex][0]*3 + j];
    
    referenceToCurrentCoordinates(xfemElement, Steinerpoint);   
    referenceToCurrentCoordinates(xfemElement, edgePoint);   
    referenceToCurrentCoordinates(xfemElement, oppositePoint);  
     
    Epetra_SerialDenseVector r1 = subtractsTwoVectors( edgePoint, Steinerpoint);
    Epetra_SerialDenseVector r2 = subtractsTwoVectors( oppositePoint, Steinerpoint);

    Epetra_SerialDenseVector n1 = computeCrossProduct( r1, r2);
    Epetra_SerialDenseVector n2 = computeCrossProduct( r1, n1);
    
    n1 = normalizeVector(n1);
    n2 = normalizeVector(n2);

    vector<Epetra_SerialDenseVector> plane(4);      
    plane[0] = addTwoVectors(Steinerpoint, n1);               
    plane[1] = subtractsTwoVectors(Steinerpoint, n1);
    plane[2] = addTwoVectors(plane[1], n2);
    plane[3] = addTwoVectors(plane[0], n2);
     
    Epetra_SerialDenseVector xsi(3);
    bool intersected = computeRecoveryPlane( lineIndex, xsi, plane, intersectingCutterElements_[cutterIndex]);
    
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
 |  RCI:    returns information of the tetrahedra            u.may 08/07|
 *----------------------------------------------------------------------*/  
void Intersection::getTetrahedronInformation( 
    const int           tetIndex,
    const int           faceIndex,  
    vector<int>&        tetraCornerIndices,
    vector<int>&        order,
    const tetgenio&     out)
{            
    
    // store boundary face node indices
    for(int j=0; j<3; j++)
        tetraCornerIndices[j] = out.trifacelist[faceIndex*3+j];
   
    // store node index opposite to the boundary face of the tetrahedron
    for(int j=0; j<4; j++)
    {
        int nodeIndex = out.tetrahedronlist[tetIndex*out.numberofcorners+j];
        if(nodeIndex != tetraCornerIndices[0] && nodeIndex != tetraCornerIndices[1] && nodeIndex != tetraCornerIndices[2])
        {
            tetraCornerIndices[3] = out.tetrahedronlist[tetIndex*out.numberofcorners+j];
            break;
        }
    }
    
    for(int j=0; j<4; j++)
    {
        int nodeIndex = out.tetrahedronlist[tetIndex*out.numberofcorners+j];
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
    vector<int>&                            tetraCornerIndices,
    DRT::Element*                           xfemElement,
    const tetgenio&                         out)
{

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 3; j++)
            tetraCornerNodes[i][j] = out.pointlist[tetraCornerIndices[i]*3+j];
        
        referenceToCurrentCoordinates(xfemElement, tetraCornerNodes[i]);   
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    lifts the higher-order point of an edge of the   u.may 09/07|
 |          linearized interface onto the curved interface              |   
 *----------------------------------------------------------------------*/
void Intersection::computeHigherOrderPoint(    
    int                                 index1, 
    int                                 index2, 
    int                                 faceIndex, 
    int                                 faceMarker, 
    int                                 globalHigherOrderIndex, 
    vector<int>&                        tetraCornerIndices, 
    vector<Epetra_SerialDenseVector>&   tetraCornerNodes, 
    DRT::Element*                       xfemElement, 
    tetgenio&                           out)
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
                           
    //printf("adjacentfaceMarker = %d\n",adjacentFaceMarker);
    
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
        printf("NO INTERSECTION POINT FOUND!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    }

}



/*----------------------------------------------------------------------*
 |  RCI:    returns the other two point indices belonging    u.may 09/07|
 |          to a triface that obtains a Steiner point                   |           
 *----------------------------------------------------------------------*/
vector<int> Intersection::getPointIndices(
    tetgenio&   out, 
    int         trifaceIndex, 
    int         steinerPointIndex)
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
    DRT::Element*                               cutterElement,
    const bool                                  onBoundary)
{
    bool                        intersection = true;
    int                         iter = 0;
    int                         countSingular = 0;
    const int                   maxiter = 50;
    double                      residual = 1.0;
    Epetra_SerialDenseMatrix    A(3,3);
    Epetra_SerialDenseVector    b(3);
    Epetra_SerialDenseVector    dx(3);
    
    for(int i = 0; i < 3; i++)
        xsi[i]= 0.0; 
        
    dx = xsi;
    
    updateRHSForRCINormal( b, xsi, normal, cutterElement, onBoundary);
                                
    while(residual > TOL14_)
    {   
        updateAForRCINormal( A, xsi, normal, cutterElement, onBoundary);
         
        if(!solveLinearSystemWithSVD(A, b, dx, 3))  
            countSingular++;
    
        if(countSingular > 5)
        {
            intersection = false;
            break;  
        }
        
        xsi = addTwoVectors(xsi,dx);
        //printf("dx0 = %20.16f\t, dx1 = %20.16f\t, dx2 = %20.16f\n", dx[0], dx[1], dx[2]);
        if(iter >= maxiter)
        {   
            intersection = false;
            break;
        }       
        
        updateRHSForRCINormal( b, xsi, normal, cutterElement, onBoundary); 
        residual = b.Norm2(); 
        iter++;
        
        //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi[0], xsi[1], xsi[2], residual, TOL14_);
    } 
    
    if( (fabs(xsi[0])-1.0) > TOL7_  || (fabs(xsi[1])-1.0) > TOL7_ )    // line coordinate may be bigger than 1
    {
        printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi[0], xsi[1], xsi[2], residual, TOL14_);
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
    DRT::Element*                               surfaceElement,
    const bool                                  onBoundary)                                                 
{   
    const int numNodesSurface = surfaceElement->NumNode();
    const int numNodesLine = 3;
    vector<int> actParams(1,0);
    Epetra_SerialDenseMatrix surfaceDeriv1(2,numNodesSurface);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    Epetra_SerialDenseVector xsiSurface(2);
    Epetra_SerialDenseVector xsiLine(1);
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
     
    xsiSurface(0) = xsi[0]; // r-coordinate surface
    xsiSurface(1) = xsi[1]; // s-coordinate surface
    xsiLine(0)    = xsi[2]; // r-coordinate line
    actParams[0] = numNodesSurface;     
    
    for(unsigned int i=0; i<3; i++)
        for(unsigned int j=0; j<3; j++)
            A[i][j] = 0.0;
    
    params.set("action","calc_ShapeDeriv1");
    surfaceElement->Evaluate(params, dummyDis, actParams, surfaceDeriv1, emptyM , emptyV, xsiSurface, emptyV);          
    //printf("num of nodes surfaces = %d\n", numOfNodesSurface);
        
    if(!onBoundary)
    {
        for(int dim=0; dim<3; dim++)
        {
            for(int i=0; i<numNodesSurface; i++)
            {
                A[dim][0] += surfaceElement->Nodes()[i]->X()[dim] * surfaceDeriv1[i][0];
                A[dim][1] += surfaceElement->Nodes()[i]->X()[dim] * surfaceDeriv1[i][1];
            }
            A[dim][2] = (-1.0) * ( normal[0][dim] * (-0.5) + normal[1][dim] * 0.5 );  
        }
    }
    else
    {
        Epetra_SerialDenseMatrix lineDeriv1(1,numNodesLine);
        shape_function_1D_deriv1(lineDeriv1, xsiLine[0], DRT::Element::line3);
     
        for(int dim=0; dim<3; dim++)
        {
            for(int i=0; i<numNodesSurface; i++)
            {
                A[dim][0] += surfaceElement->Nodes()[i]->X()[dim] * surfaceDeriv1[i][0];
                A[dim][1] += surfaceElement->Nodes()[i]->X()[dim] * surfaceDeriv1[i][1];
            }
            
            for(int i = 0; i < 3; i++)
            {
                int index = i;
                if(i > 1)   index = 4;
                A[dim][2] += (-1.0) * normal[index][dim] * lineDeriv1[i][0]; 
            } 
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
    Epetra_SerialDenseVector&                   xsi,    
    const vector<Epetra_SerialDenseVector>&     normal,
    DRT::Element*                               surfaceElement,
    const bool                                  onBoundary)                                                    
{
    const int numNodesSurface = surfaceElement->NumNode();
    const int numNodesLine = 3;
    vector<int> actParams(1,0);
    Epetra_SerialDenseVector surfaceFunct(numNodesSurface);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    Epetra_SerialDenseVector xsiSurface(2);
    Epetra_SerialDenseVector xsiLine(1);
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
    
    
    xsiSurface(0) = xsi[0]; // r-coordinate surface
    xsiSurface(1) = xsi[1]; // s-coordinate surface
    xsiLine(0)    = xsi[2]; // r-coordinate line
    params.set("action","calc_Shapefunction");
            
    for(unsigned int i=0; i<3; i++)   
        b[i] = 0.0;
    
    actParams[0] = numNodesSurface;            
    surfaceElement->Evaluate(params, dummyDis, actParams, emptyM, emptyM , surfaceFunct, xsiSurface, emptyV);           
  
    if(!onBoundary)
    {
        for(int dim=0; dim<3; dim++)
            for(int i=0; i<numNodesSurface; i++)
            {
                b[dim] += (-1.0) * surfaceElement->Nodes()[i]->X()[dim] * surfaceFunct[i];
            }
            
        for(int dim=0; dim<3; dim++)
            b[dim] += normal[0][dim] * 0.5*(1.0 - xsi[2]) + normal[1][dim] * 0.5*(1.0 + xsi[2]);
    }
    else
    {
        Epetra_SerialDenseVector lineFunct(numNodesLine);
        shape_function_1D( lineFunct, xsiLine[0], DRT::Element::line3);
       
        for(int dim=0; dim<3; dim++)
            for(int i=0; i<numNodesSurface; i++)
            {
                b[dim] += (-1.0) * surfaceElement->Nodes()[i]->X()[dim] * surfaceFunct[i];
            }
        
        for(int dim=0; dim<3; dim++)
            for(int i=0; i<numNodesLine; i++)
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
    DRT::Element*                               surfaceElement)
{
    bool    intersection = true;
    int     begin, end;
    int     numLines = surfaceElement->NumLine();
   
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
        DRT::Element*               lineElement = surfaceElement->Lines()[i];
        Epetra_SerialDenseMatrix    A(3,3);
        Epetra_SerialDenseVector    b(3);
        Epetra_SerialDenseVector    dx(3);
      
        intersection = true;
        // starting value equals (0,0,0)
        for(int j = 0; j < 3; j++)
            xsi[j]= 0.0; 
        
        dx = xsi;
    
        updateRHSForRCIPlane( b, xsi, plane, lineElement);
                        
        while( residual > TOL14_ )
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
        
            xsi = addTwoVectors(xsi, dx);
  
            updateRHSForRCIPlane( b, xsi, plane, lineElement);
            residual = b.Norm2(); 
            iter++;
        } 
    
        if( (fabs(xsi[2])-1.0) > TOL7_ )     // planes coordinate may be bigger than 1
        {   printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi[0], xsi[1], xsi[2], residual, TOL14_);
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
    DRT::Element*                               lineElement,
    DRT::Element*                               surfaceElement)
{   
    const int numNodesLine = lineElement->NumNode();
    const int numNodesSurface = 4;
    vector<int> actParams(1,0);
    Epetra_SerialDenseMatrix surfaceDeriv(2, numNodesSurface);
    Epetra_SerialDenseMatrix lineDeriv(1,numNodesLine);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    Epetra_SerialDenseVector xsiSurface(2);
    Epetra_SerialDenseVector xsiLine(1);
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
     
    xsiSurface[0] = xsi[0]; // r-coordinate surface
    xsiSurface[1] = xsi[1]; // s-coordinate surface
    xsiLine[0]    = xsi[2]; // r-coordinate line
    
    for(unsigned int i=0; i<3; i++)
        for(unsigned int j=0; j<3; j++)
            A[i][j] = 0.0;
    
    
    params.set("action","calc_ShapeDeriv1");
    //actParams[0] = surfaceElement->NumNode();  
    //surfaceElement->Evaluate(params, dummyDis, actParams, surfaceDeriv, emptyM , emptyV, xsiSurface, emptyV);   
    shape_function_2D_deriv1(surfaceDeriv,  xsiSurface[0],  xsiSurface[1], DRT::Element::quad4);
       
    actParams[0] = numNodesLine;    
    lineElement->Evaluate(params, dummyDis, actParams, lineDeriv, emptyM , emptyV, xsiLine, emptyV);               
   
    for(int dim=0; dim<3; dim++)
    {
        for(int i=0; i<numNodesSurface; i++)
        {
            A[dim][0] += plane[i][dim] * surfaceDeriv[i][0];
            A[dim][1] += plane[i][dim] * surfaceDeriv[i][1];
        }
    }
        
    for(int dim=0; dim<3; dim++)
    {
        for(int i=0; i<numNodesLine; i++)
        {
            A[dim][2] += (-1.0) * lineElement->Nodes()[i]->X()[dim] * lineDeriv[i][0];   
        }
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
    Epetra_SerialDenseVector&                   xsi,    
    const vector <Epetra_SerialDenseVector>&    plane,
    DRT::Element*                               lineElement)                                                    
{
    const int numNodesLine    = lineElement->NumNode();
    const int numNodesSurface = 4;
    vector<int> actParams(1,0);
    Epetra_SerialDenseVector surfaceFunct(numNodesSurface);
    Epetra_SerialDenseVector lineFunct(numNodesLine);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    Epetra_SerialDenseVector xsiSurface(2);
    Epetra_SerialDenseVector xsiLine(1);
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
    
    xsiSurface[0] = xsi[0]; // r-coordinate surface
    xsiSurface[1] = xsi[1]; // s-coordinate surface
    xsiLine[0]    = xsi[2]; // r-coordinate line
    params.set("action","calc_Shapefunction");
            
    for(unsigned int i=0; i<3; i++)   
        b[i] = 0.0;
    
    // shape function for the plane 
    shape_function_2D( surfaceFunct, xsiSurface[0], xsiSurface[1], DRT::Element::quad4 ); 
    // shape function for the curve
    actParams[0] = numNodesLine;   
    lineElement->Evaluate(params, dummyDis, actParams, emptyM, emptyV , lineFunct, xsiLine, emptyV);    
  
    for(int dim=0; dim<3; dim++)
        for(int i=0; i<numNodesSurface; i++)
        {
            b[dim] += (-1.0) * plane[i][dim] * surfaceFunct[i];
        }
        
    for(int dim=0; dim<3; dim++)
        for(int i=0; i<numNodesLine; i++)
        {
            b[dim] +=  lineElement->Nodes()[i]->X()[dim]  * lineFunct[i];        
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
    DRT::Element*                           xfemElement,
    const tetgenio&                         out)
{            
    
    Epetra_SerialDenseVector  p1(3);
    Epetra_SerialDenseVector  p2(3);
    Epetra_SerialDenseVector  p3(3);        
    Epetra_SerialDenseVector  m(3);
    Epetra_SerialDenseVector  n(3);
    Epetra_SerialDenseVector  r(3);
    Epetra_SerialDenseVector  r1(3);
    Epetra_SerialDenseVector  r2(3);
    
    
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
    r1 = subtractsTwoVectors(p1, p2);
    r2 = subtractsTwoVectors(p3, p2);
 
    // normal of the plane
    n = computeCrossProduct(r1, r2);
    n = normalizeVector(n);
 
    // direction vector of the intersection line
    r = computeCrossProduct(n, r2);  
    r = normalizeVector(r);
 
    // computes the start point of the line
    m = computeLineMidpoint(p2, p3);
    
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
            referenceToCurrentCoordinates(xfemElement, plane[i]);
        
        referenceToCurrentCoordinates(xfemElement, m);
          
        plane[4] = m;
      
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    computes the normal to the interface edge of     u.may 11/07|
 |          two adjacent triangular faces                               |
 *----------------------------------------------------------------------*/  
void Intersection::computeIntersectionNormal( 
    int                                 index1,
    int                                 index2, 
    int                                 faceIndex,
    int                                 adjacentFaceIndex,
    int                                 globalHigherOrderIndex,
    vector<Epetra_SerialDenseVector>&   plane,
    DRT::Element*                       xfemElement,
    tetgenio&                           out)
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
    
    referenceToCurrentCoordinates(xfemElement, p1);  
    referenceToCurrentCoordinates(xfemElement, p2);  
    referenceToCurrentCoordinates(xfemElement, p3);  
    referenceToCurrentCoordinates(xfemElement, p4);  
    
    Epetra_SerialDenseVector r1 = subtractsTwoVectors(p1,p2);
    Epetra_SerialDenseVector r2 = subtractsTwoVectors(p3,p2);
    Epetra_SerialDenseVector r3 = subtractsTwoVectors(p4,p2);
    
    Epetra_SerialDenseVector n1 = computeCrossProduct(r2, r1);
    Epetra_SerialDenseVector n2 = computeCrossProduct(r1, r3);
    
    Epetra_SerialDenseVector averageNormal = addTwoVectors(n1, n2);
    Epetra_SerialDenseVector rPlane = computeCrossProduct(n1, r1);
    
    for(int i = 0; i < 3; i++)
        averageNormal[i] = 0.5*averageNormal[i];
        
    averageNormal  = normalizeVector(averageNormal);
    rPlane  = normalizeVector(rPlane);
 
    Epetra_SerialDenseVector m(3);
    for(int i = 0; i < 3; i++)
        m[i] = out.pointlist[globalHigherOrderIndex*3+i];

    referenceToCurrentCoordinates(xfemElement, m);  
    
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
 |  RCI:    computes the midpoint of a line                  u.may 08/07|
 *----------------------------------------------------------------------*/  
Epetra_SerialDenseVector Intersection::computeLineMidpoint( 
    const Epetra_SerialDenseVector& p1, 
    const Epetra_SerialDenseVector& p2)
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
    tetgenio&   out)
{

    bool    faceMarkerFound     = false;
    
    for(int i=0; i<out.numberoftrifaces; i++)
    {
        adjacentFaceMarker = out.trifacemarkerlist[i] - facetMarkerOffset_;
        adjacentFaceIndex = i;
        
        if(adjacentFaceMarker >= -1)
        {        
            int countPoints = 0;
            for(int j = 0; j < 3 ; j++)
            {
                int pointIndex = out.trifacelist[ i*3 + j ];
                if(pointIndex == edgeIndex1 || pointIndex == edgeIndex2)
                    countPoints++;
            }
            
            if((countPoints == 2) && (faceIndex != adjacentFaceIndex))
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
    tetgenio&                           out)
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
    vector<int>&                        adjacentFacesList,
    Epetra_SerialDenseVector&           edgePoint,
    Epetra_SerialDenseVector&           oppositePoint,
    tetgenio&                           out)
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
    int&                                            cutterIndex)
{
    bool comparison = false;
    DRT::Node* node1;
    DRT::Node* node2;
    const int numLines1 = intersectingCutterElements_[faceIndex1]->NumLine(); 
    const int numLines2 = intersectingCutterElements_[faceIndex2]->NumLine(); 
    const int numNodes  = intersectingCutterElements_[faceIndex2]->Lines()[0]->NumNode();
    
    
    for(int i = 0; i < numLines1; i++)
    {
        for(int j = 0; j < numLines2; j++)
        {
            comparison = true;
            for(int k  = 0; k < numNodes; k++)
            {
                node1 = intersectingCutterElements_[faceIndex1]->Lines()[i]->Nodes()[k];
                node2 = intersectingCutterElements_[faceIndex2]->Lines()[j]->Nodes()[k];
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
                    if(k==2)
                    {
                        node1 = intersectingCutterElements_[faceIndex1]->Lines()[i]->Nodes()[k];
                        node2 = intersectingCutterElements_[faceIndex2]->Lines()[j]->Nodes()[k];
                    }
                    else
                    {
                        node1 = intersectingCutterElements_[faceIndex1]->Lines()[i]->Nodes()[k];
                        node2 = intersectingCutterElements_[faceIndex2]->Lines()[j]->Nodes()[1-k];
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
    DRT::Element*                       xfemElement,
    DRT::Element*                       cutterElement,
    Epetra_SerialDenseVector&           edgeNode1,
    Epetra_SerialDenseVector&           edgeNode2)
{

    int lineIndex = -1;
    Epetra_SerialDenseVector x1(1);
    Epetra_SerialDenseVector x2(1);
    
    Epetra_SerialDenseVector node1(3);
    Epetra_SerialDenseVector node2(3);
    
    for(int i = 0; i < 3; i++)
    {
        node1[i] = edgeNode1[i];
        node2[i] = edgeNode2[i];
    }
      
    referenceToCurrentCoordinates(xfemElement, node1);
    referenceToCurrentCoordinates(xfemElement, node2);
    
    x1[0] = node1[0];
    x2[0] = node2[0];
    
    for(int i = 0; i < cutterElement->NumLine(); i++)
    {
        Epetra_SerialDenseVector xsi1(1);
        Epetra_SerialDenseVector xsi2(1);
        DRT::Element*  lineElement = cutterElement->Lines()[i];
        
        bool check1 = checkNodeWithinElement( lineElement, x1, xsi1);
        bool check2 = checkNodeWithinElement( lineElement, x2, xsi2);
        if( check1  &&  check2 )
        {   
            lineIndex = i;  //countIndex;
            break;
        }
    }
    
    /*
    for(int i = 0; i < cutterElement->Lines()[lineIndex]->NumNode(); i++ )
    {
        cutterElement->Lines()[lineIndex]->Nodes()[i]->Print(cout); 
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
    DRT::Element*                               xfemElement, 
    tetgenio&                                   out)
{
    
    Epetra_SerialDenseVector xsiLine(3);
    for(int i = 0; i < 3; i++)
        xsiLine[i] = 0.0;
        
    if(normal)     referenceToCurrentCoordinates(surfaceElement, xsi);
    else           
    {   
        xsiLine[0] = xsi[2];
        referenceToCurrentCoordinates(surfaceElement->Lines()[lineIndex], xsiLine);
        for(int i=0; i<3; i++)
            xsi[i] = xsiLine[i];
    }       
    currentToReferenceCoordinates(xfemElement, xsi);
    
    //printf("xsiold0 = %20.16f\t, xsiold1 = %20.16f\t, xsiold2 = %20.16f\n", out.pointlist[index*3], out.pointlist[index*3+1], out.pointlist[index*3+2]);
    
    for(int i = 0; i < 3; i++)
        out.pointlist[globalHigherOrderIndex*3+i]   = xsi[i];  
   
    //printf("xsi0    = %20.16f\t, xsi1    = %20.16f\t, xsi2    = %20.16f\n", xsi[0], xsi[1], xsi[2]);
    //printf("\n");
}



/*----------------------------------------------------------------------*
 |  DB:     Debug only                                       u.may 06/07|
 *----------------------------------------------------------------------*/  
void Intersection::debugXAABBIntersection( 	Epetra_SerialDenseMatrix cutterXAABB, 
											Epetra_SerialDenseMatrix xfemXAABB,
											DRT::Element* cutterElement,
				 							DRT::Element* xfemElement,
				 							int noC,
				 							int noX)
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
void Intersection::debugNodeWithinElement(  DRT::Element* element,
                                            DRT::Node* node,
                                            Epetra_SerialDenseVector& xsi,
                                            int noE,
                                            int noN,
                                            bool within)
{
    int numnodes = element->NumNode();
    vector<int> actParams(1,0);
    Epetra_SerialDenseVector funct(numnodes);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    Epetra_SerialDenseVector x(3);
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
      
    params.set("action","calc_Shapefunction");
    actParams[0] = numnodes;   
      
    element->Evaluate(params, dummyDis, actParams, emptyM, emptyM , funct, xsi, emptyV);        
    
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
void Intersection::debugTetgenDataStructure(    DRT::Element*               element)
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
        referenceToCurrentCoordinates(element, xsi);
        
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
    									DRT::Element * element,
    									vector<int>& elementIds)
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
 |  DB:     Debug only                                       u.may 06/07|
 *----------------------------------------------------------------------*/  
void Intersection::debugDomainIntCells(	map< int, DomainIntCells >&	domainintcells,
											int id)
{    
    cout << endl;
    cout << "===============================================================" << endl;
    cout << "Debug RESULTING INTEGRATION CELLS " << endl;
    cout << "===============================================================" << endl;
    cout << endl;
    
    map<int, DomainIntCells >::iterator iter;
  	for( iter = domainintcells.begin(); iter != domainintcells.end(); ++iter ) 
    {
		cout << "XFEM ELEMENT " << iter->first << " :" << endl;	
		cout << endl;
		
		if(iter->first == id)
		{
			cout << "NUMBER OF INTEGRATIONCELLS : " << iter->second.size() << endl;	
			cout << endl;
			for(unsigned int j = 0; j < iter->second.size(); j++ )
			{
				cout << "IC " << j << ":  " << endl;
				cout << endl;	
				for(unsigned int k = 0; k < iter->second[j].GetCoord().size(); k++)
				{
					//cout << k << "\t";
					for(unsigned int m = 0; m < iter->second[j].GetCoord()[k].size(); m++)
					{
						cout << iter->second[j].GetCoord()[k][m] << "\t";	
					}
					cout << endl;  
				}
				cout << endl;	
			}
			cout << endl;	cout << endl;
		}	
	}		
       
    cout << endl;
    cout << endl;
    cout << endl;
    cout << "===============================================================" << endl;
    cout << "Debug RESULTING INTEGRATION CELLS" << endl;
    cout << "===============================================================" << endl;
    cout << endl; cout << endl; cout << endl;
}



/*----------------------------------------------------------------------*
 |  DB:     Debug only                                       u.may 09/07|
 *----------------------------------------------------------------------*/  
void Intersection::printTetViewOutput(
    int             index,
    tetgenio&       out)
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
    DRT::Element*   xfemElement,
    int             index,
    tetgenio&       in)
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
        
        referenceToCurrentCoordinates(xfemElement, xsi);
        
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




/*----------------------------------------------------------------------*
 |  computes the coordinates of a point within a region for  u.may 07/07|
 |  the Tetgen data structure       (only for visualization)            |
 *----------------------------------------------------------------------*/  
// dependency on Fluid3Surface. This cannot be used for the standard fluid tests
// so switch on, if needed, but don't commit
#if 0
void Intersection::computeRegionCoordinates(    
    DRT::Element*  xfemElement,
    DRT::Element*  cutterElement,
    double* regionCoordinates)
{

    bool inStored = false;
    bool outStored = false;
    bool nodeWithin = false;
    int fill = 0;
    int numCornerPoints = 8;
    DRT::Element*  cutterVolume;
    Epetra_SerialDenseVector xsi(3);
    
    
    const DRT::Element::ElementType etype = cutterElement->Type();
  
    switch(etype)
    {
        case DRT::Element::element_fluid3surface:
        {       
            DRT::Elements::Fluid3Surface* fluid3surface = (DRT::Elements::Fluid3Surface *) (cutterElement);
            
            if(fluid3surface == NULL) 
                dserror("fluid3surface is NULL");
                
            cutterVolume = fluid3surface->GetParent();
            
            if(cutterVolume == NULL) 
                dserror("cuttervolume is NULL");
            break;
        }
        default:
            dserror("elementtype not yet handled");
    }

    for(int i = 0; i < numCornerPoints; i++)
    {
        Epetra_SerialDenseVector x;       
        for(int ii = 0; ii < 3; ii++)
            x[ii] = xfemElement->Nodes()[i]->X()[ii];
            
        nodeWithin = checkNodeWithinElement(cutterVolume, x, xsi);
                                
        if(nodeWithin && !inStored)
        {
            for(int j = 0; j < 3; j++)
                regionCoordinates[fill++] = xfemElement->Nodes()[i]->X()[j];
            
            inStored = true;
        }
        
        if(!nodeWithin && !outStored)
        {
            for(int j = 0; j < 3; j++)
                regionCoordinates[fill++] = xfemElement->Nodes()[i]->X()[j];
    
            outStored = true;
        }
        
        if(inStored && outStored)
        {
            break;
        }
    }
}
#endif




#endif  // #ifdef CCADISCRET
