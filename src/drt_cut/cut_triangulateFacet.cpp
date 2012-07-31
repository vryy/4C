#include "cut_facet.H"
#include "cut_triangulateFacet.H"
#include "cut_kernel.H"

/*-----------------------------------------------------------------------------------------------------------*
              Split the facet into appropriate number of tri and quad                           Sudhakar 04/12
              Work well for both convex and concave facets
*------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::SplitFacet()
{

//  std::cout<<"number of points before deleting = "<<ptlist_.size()<<"\n";

  KERNEL::DeleteInlinePts( ptlist_ );

  /*****************************************************************************************/
  /*std::cout<<"NUmber of points after deleting = "<<ptlist_.size()<<"\n";//blockkk
  std::cout<<"The coordinates of the facet after deleting\n";
  for( unsigned i=0;i<ptlist_.size();i++ )
  {
    Point* ptx = ptlist_[i];
    double coox[3];
    ptx->Coordinates(coox);
    std::cout<<coox[0]<<"\t"<<coox[1]<<"\t"<<coox[2]<<"\n";
  }*/
  /*****************************************************************************************/


  // Deal with zero point facet -- this should never occur, but happens in axel10 cut-test
  if( ptlist_.size()==0 )
    return;

  //std::cout<<"number of points after deleting = "<<ptlist_.size()<<"\n";

  int numpts = ptlist_.size();
  //std::cout<<"numpts = "<<numpts<<"\n";

  if( numpts<3 )
  {
    dserror("A facet must have atleast 3 corners");
  }
  else if ( numpts==3 )
  {
    split_.push_back( ptlist_ );
    return;
  }
  else
  {
    split_.clear();

    std::string geoType;
    std::vector<int> ptConcavity = KERNEL::CheckConvexity(  ptlist_, geoType );

    if( geoType=="convex" || geoType=="1ptConcave" )
      SplitConvex_1ptConcave_Facet( ptConcavity );
    else
      SplitGeneralFacet( ptConcavity );
  }
}

/*-------------------------------------------------------------------------------------*
   If a 4 noded facet is convex, a quad is made. If it is concave,
   it is split into appropriate triangles.
   The Gauss integration rule is available for quad only if it is convex

                                   /\
                                  / .\
                                 /A .B\                                   sudhakar 04/12
                                /  ++  \
                               / +    + \
                               +        +
*--------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::Split4nodeFacet( std::vector<Point*> &poly,
                                                  bool callFromSplitAnyFacet )
{
//  std::cout<<"I am in 4 noded facet program\n";
  if( poly.size()!=4 )
    dserror("This is not a 4 noded facet");

  std::string geoType;
  std::vector<int> ptConcavity = KERNEL::CheckConvexity(  poly, geoType );


  /*************************************************************/
  /*std::cout<<"the 4 noded facet is  "<<geoType<<"\n";//blockkk
  std::cout<<"The coordinates of the facet after deleting\n";
  for( unsigned i=0;i<poly.size();i++ )
  {
    Point* ptx = poly[i];
    double coox[3];
    ptx->Coordinates(coox);
    std::cout<<coox[0]<<"\t"<<coox[1]<<"\t"<<coox[2]<<"\n";
  }*/
/***********************************************************/
  int indStart=0;
  if( geoType=="convex" )
  {
    split_.push_back( poly );
    return;
  }
  else if( geoType=="1ptConcave" )
  {
    indStart = ptConcavity[0];

    std::vector<Point*> temp1(3),temp2(3);
    temp1[0] = poly[indStart];
    temp1[1] = poly[(indStart+1)%4];
    temp1[2] = poly[(indStart+2)%4];

    temp2[0] = poly[indStart];
    temp2[1] = poly[(indStart+2)%4];
    temp2[2] = poly[(indStart+3)%4];

    if( split_.size()==0 )
    {
      split_.resize(2);
      split_[0] = temp1;
      split_[1] = temp2;
    }
    else
    {
      split_.push_back( temp1 );
      split_.push_back( temp2 );
    }
    return;
  }
  else // if there are two concave pts for 4 noded facet --> selfcut
  {
    // this means splitAnyFacet is failed --> call earclipping
    /*if( callFromSplitAnyFacet )
    {
      std::cout<<"WARNING!!!! calling earclipping because splitanyfacet failed\n";
      ptConcavity = KERNEL::CheckConvexity(  ptlist_, geoType );
      split_.clear();
      EarClipping( ptConcavity, false );
      return;
    }*/
    std::cout<<"the points are\n";
    for( unsigned i=0;i<poly.size();i++ )
    {
      Point* p1 = poly[i];
      double x1[3];
      p1->Coordinates(x1);
      std::cout<<x1[0]<<"\t"<<x1[1]<<"\t"<<x1[2]<<"\n";
    }
    dserror( "a 4 noded facet cannot have more than 1 concave point:"
        "This means that the facet is selfcut" );
  }
}

/*---------------------------------------------------------------------------------------------------*
   A facet that is convex or having one concave point is split into 1 Tri and few Quad cells
   splitting starts with concave pt --- eliminate the need to check whether any other pt is
   inside the newCell formed                                                            Sudhakar 07/12
*----------------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::SplitConvex_1ptConcave_Facet( std::vector<int> ptConcavity )
{
  //std::cout<<"Splitting the convex or 1ptConcave facet\n";//blockkk

  /************************************************************************/
  /*std::cout<<"After this loop the no of points = "<<ptlist_.size()<<"\n";
  for( unsigned i=0;i<ptlist_.size();i++ )
  {
    Point* ptx = ptlist_[i];
    double coox[3];
    ptx->Coordinates(coox);
    std::cout<<coox[0]<<"\t"<<coox[1]<<"\t"<<coox[2]<<"\n";
  }*/
  /************************************************************************/

  if( ptConcavity.size() > 1 )
    dserror( "should be called only when the facet has one or no concave points" );

  int num = ptlist_.size();

  if( ptlist_.size()==3 )
  {
    split_.push_back( ptlist_ );
    return;
  }
  else if ( ptlist_.size()==4 )
  {
    Split4nodeFacet( ptlist_, false );
    return;
  }

  bool triDone = false, convex = true;;
  int firstPt=0,secondPt=0,thirdPt=0,lastPt=1,endPt=num-1;
  if( ptConcavity.size()==1 )
  {
    convex = false;
    firstPt = ptConcavity[0];
    lastPt = (firstPt+1)%num; //initialization
    if( firstPt!=0 )
      endPt = firstPt-1;
  }
  std::vector<Point*> newCell;

  while(1)
  {
    newCell.clear();

    secondPt = lastPt;  //properly initialized
    thirdPt = (secondPt+1)%num;

    newCell.push_back(ptlist_[firstPt]);
    newCell.push_back(ptlist_[secondPt]);
    newCell.push_back(ptlist_[thirdPt]);    // tri cell is now formed

    if( thirdPt==endPt )
      triDone = true;
    else if ( !triDone )      // check whether tri can be extended to quad cell
    {
      lastPt = (thirdPt+1)%num;
      newCell.push_back(ptlist_[lastPt]);
      if( lastPt==endPt )
        triDone = true;
    }

    KERNEL::DeleteInlinePts( newCell );

    if( newCell.size()==3 || convex )
      split_.push_back( newCell );
    else if( newCell.size()==4 ) //check the resulting quad is convex, else split into 2 tri
      Split4nodeFacet( newCell, true );
    else
      dserror( "should have either 2 or 3 points" );

    if( triDone )
      break;
  }
}

/*---------------------------------------------------------------------------------------------------*
 * Generalized facet splitting procedure which works for simple facets with any number  sudhakar 08/12
 * of concave points. Involves checking whether a reflex point is inside formed cell
 *---------------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::SplitGeneralFacet( std::vector<int> ptConcavity )
{
  //std::cout<<"splitting general facet\n"; //blockkk
  if( ptConcavity.size() < 2 )
    dserror( "Call TriangulateFacet::SplitConvex_1ptConcave_Facet in such cases" );

  if( ptlist_.size()==3 ) // directly form a Tri cell
  {
    split_.push_back( ptlist_ );
    return;
  }
  else if ( ptlist_.size()==4 ) // can be a convex or concave quad
  {
    Split4nodeFacet( ptlist_, false );
    return;
  }

  int num = ptlist_.size();
  int concsize = ptConcavity.size();

  int firstPt=0,secondPt=1,thirdPt=2,fourthPt=3;
  std::vector<Point*> newCell;

  while(1)
  {
    if( (num-concsize) < 4 ) // this means that no Quad cells can be formed for this geometry
    {
      EarClipping( ptConcavity );
      return;
    }

    newCell.clear();
    int ncross = 0;
    for( int i=0;i<num;++i )
    {
      newCell.clear();
      ncross++;
      int concNo = 0;

      if( i==0 )
      {
        firstPt = ptConcavity[concNo];
        secondPt = (firstPt+1)%num;
      }
      else
      {
        firstPt++;
        secondPt = (firstPt+1)%num;
      }
      if(std::find(ptConcavity.begin(), ptConcavity.end(), secondPt) != ptConcavity.end())
        continue;

      thirdPt = (secondPt+1)%num;

      newCell.push_back(ptlist_[firstPt]);
      newCell.push_back(ptlist_[secondPt]);
      newCell.push_back(ptlist_[thirdPt]);    // tri cell is now formed

      if(std::find(ptConcavity.begin(), ptConcavity.end(), thirdPt) == ptConcavity.end())
      {
        fourthPt = (thirdPt+1)%num;
        newCell.push_back(ptlist_[fourthPt]);
      }

      KERNEL::DeleteInlinePts( newCell );

      /********************************************************/
      /*std::cout<<"the newCell points are\n";
      for( unsigned i=0;i<newCell.size();i++ )
      {
        Point* ptx = newCell[i];
        double coox[3];
        ptx->Coordinates(coox);
        std::cout<<coox[0]<<"\t"<<coox[1]<<"\t"<<coox[2]<<"\n";
      }*/
      /********************************************************/

      bool isEar=true;

      if( newCell.size()==3 )
      {
        for( unsigned j=0;j<ptConcavity.size();j++ )
        {
          unsigned reflInd = ptConcavity[j];
          if(std::find(newCell.begin(), newCell.end(), ptlist_[reflInd]) != newCell.end())
            continue;
          if( KERNEL::PtInsideTriangle( newCell, ptlist_[reflInd]) )
          {
            isEar = false;
            break;
          }
        }
        if( isEar )
        {
          ptlist_.erase( ptlist_.begin()+secondPt ); // erase a pt to form new polygon
          split_.push_back( newCell );
          break;
        }
      }
      else if( newCell.size()==4 ) //check the resulting quad is convex, else split into 2 tri
      {
        for( unsigned j=0;j<ptConcavity.size();j++ )
        {
          unsigned reflInd = ptConcavity[j];
          if(std::find(newCell.begin(), newCell.end(), ptlist_[reflInd]) != newCell.end())
            continue;
          if( KERNEL::PtInsideQuad( newCell, ptlist_[reflInd]) )
          {
            isEar = false;
            break;
          }
        }

        if( isEar )
        {
          Split4nodeFacet( newCell, true );

          // erase internal points of cell to form new polygon
          // when an element is deleted, all other elements are renumbered
          // if-else condition to make sure correct points are deleted
          if( thirdPt==0 )
          {
            ptlist_.erase( ptlist_.begin()+secondPt );
            ptlist_.erase( ptlist_.begin()+thirdPt );
          }
          else
          {
            ptlist_.erase( ptlist_.begin()+thirdPt );
            ptlist_.erase( ptlist_.begin()+secondPt );
          }
          break;
        }
      }
      else
      {
        std::cout<<"number of points in the cell = "<<newCell.size()<<"\n";
        dserror( "neither tri nor quad: Something went wrong in SplitAnyFacet" );
      }

      if( ncross==num )
        dserror("cannot form cell even after making one cycle");
    }

    KERNEL::DeleteInlinePts( ptlist_ );
    num = ptlist_.size();
    if( num==3 )
    {
      split_.push_back( ptlist_ );
      return;
    }
    else if ( num==4 )
    {
      Split4nodeFacet( ptlist_, false );
      return;
    }
    else
    {
      ptConcavity.clear();
      std::string geoType;

      ptConcavity = KERNEL::CheckConvexity(  ptlist_, geoType );

      concsize = ptConcavity.size();
      if( concsize < 2 ) // new ptlist_ forms a convex facet
      {
        SplitConvex_1ptConcave_Facet( ptConcavity );
        return;
      }
    }
    /*****************************************************************************/
    /*std::cout<<"After this loop the no of points = "<<ptlist_.size()<<"\n";
    for( unsigned i=0;i<ptlist_.size();i++ )
    {
      Point* ptx = ptlist_[i];
      double coox[3];
      ptx->Coordinates(coox);
      std::cout<<coox[0]<<"\t"<<coox[1]<<"\t"<<coox[2]<<"\n";
    }
    std::cout<<"the concave points are = ";
    for( unsigned i=0;i<ptConcavity.size();i++ )
      std::cout<<ptConcavity[i]<<"\t";
    std::cout<<"\n";*/
    /*****************************************************************************/
  }
}

/*------------------------------------------------------------------------------------------------------*
            check whether the polygon has two continuous concave points.
            if it is true, then ear clipping should be used to triangulate.              sudhakar 04/12
            split facet will fail in such cases
*-------------------------------------------------------------------------------------------------------*/
bool GEO::CUT::TriangulateFacet::HasTwoContinuousConcavePts( std::vector<int> ptConcavity )
{
  int siz = ptConcavity.size();
  if( siz<2 )
    return false;

  for( int i=0;i<siz;i++ )
  {
    int firstPt = ptConcavity[i];
    int seconPt = ptConcavity[(i+1)%siz];
    if( firstPt!=siz-1 )
    {
      if( (seconPt-firstPt)==1 )
        return true;
    }
    else if( (seconPt-firstPt)==siz-1 )
      return true;
  }
  return false;
}

/*-------------------------------------------------------------------------------------------------*
    Triangulation by ear clipping. Works for all cases, but costly.                 sudhakar 04/12
    Called when facets have two adjacent concave points
    During the process, if facet is free of adjacent concave points, splitanyfacet() is called
*--------------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::EarClipping( std::vector<int> ptConcavity,   // list of concave points
                                              bool triOnly )                  // whether to create triangles only?
{
  //std::cout<<"I am ear clipped\n"; //blockkk
  std::vector<int> convex;

  /*std::cout<<"number of points before deleting = "<<ptlist_.size()<<"\n";
  std::cout<<"the points are\n";
  for( unsigned i=0;i<ptlist_.size();i++ )
  {
    Point* ptt = ptlist_[i];
    double z[3];
    ptt->Coordinates(z);
    std::cout<<z[0]<<"\t"<<z[1]<<"\t"<<z[2]<<"\n";
  }*/

  if( ptlist_.size()==3 )
  {
    split_.push_back( ptlist_ );
    return;
  }

  if( triOnly ) // creates only triangles; when ear clipping is called directly from other functions
  {
    split_.clear();
    KERNEL::DeleteInlinePts( ptlist_ );

    if( ptlist_.size()==3 ) // after deleting the inline points, it may have only 3 points
    {
      split_.push_back( ptlist_ );
      return;
    }

    ptConcavity.clear();
    std::string geoType;

    ptConcavity = KERNEL::CheckConvexity(  ptlist_, geoType );
  }

  /*************************************************************************/
  /*std::cout<<"number of points after deleting = "<<ptlist_.size()<<"\n";

  std::cout<<"the points are\n";
  for( unsigned i=0;i<ptlist_.size();i++ )
  {
    Point* ptt = ptlist_[i];
    double z[3];
    ptt->Coordinates(z);
    std::cout<<z[0]<<"\t"<<z[1]<<"\t"<<z[2]<<"\n";
  }
  std::cout<<"the concave pints are = ";
  for( unsigned i=0;i<ptConcavity.size();i++)
    std::cout<<ptConcavity[i]<<"\t";
  std::cout<<"\n";*/
  /*************************************************************************/

  while(1)
  {
    std::vector<int> reflex( ptConcavity );

    /**********************************************/
    /*std::cout<<"the reflex points are = ";
    for( unsigned i=0;i<reflex.size();i++ )
      std::cout<<reflex[i]<<"\t";
    std::cout<<"\n";*/
    /**********************************************/

    int polPts = ptlist_.size();
    convex.resize( polPts-reflex.size() );

    // Find the convex points
    // They are non-reflex points of the polygon
    int conNum=0;

    for( int i=0;i<polPts;i++ )
    {
      if( ptConcavity.size()>0 && i==ptConcavity[0] )
      {
        ptConcavity.erase( ptConcavity.begin() );
        continue;
      }
      else
      {
        convex[conNum] = i;
        conNum++;
      }
    }

    /**********************************************/
    /*std::cout<<"the convex points are = ";
    for( unsigned i=0;i<convex.size();i++ )
      std::cout<<convex[i]<<"\t";
    std::cout<<"\n";*/
    /**********************************************/

    // if (i) is an ear, the triangle formed by (i-1),i and (i+1) should be completely within the polygon
    // Find first ear point, and make the triangle and remove this ear from the polygon
    std::vector<Point*> tri(3);
    for( int i=0;i<polPts;i++ )
    {
      // a reflex point cannot be a ear
      if(std::find(reflex.begin(), reflex.end(), i) != reflex.end())
        continue;

      unsigned ind0 = i-1;
      if(i==0)
        ind0 = polPts-1;
      unsigned ind2 = (i+1)%polPts;
      tri[0] = ptlist_[ind0];
      tri[1] = ptlist_[i];
      tri[2] = ptlist_[ind2];

      bool isEar = true;

      // only reflex points are to be checked whether they are inside the tri
      for( unsigned j=0;j<reflex.size();j++ )
      {
        unsigned reflInd = reflex[j];
        if( reflInd==ind0 || reflInd==ind2 )
          continue;

        if( KERNEL::PtInsideTriangle( tri, ptlist_[reflInd]) )
        {
          isEar = false;
          break;
        }
      }

      if( !isEar )
        continue;

      /***********************************************/
      /*double ptt1[3],ptt2[3],ptt3[3];
      tri[0]->Coordinates(ptt1);
      tri[1]->Coordinates(ptt2);
      tri[2]->Coordinates(ptt3);
      std::cout<<"the triangle points are\n";
      std::cout<<ptt1[0]<<"\t"<<ptt1[1]<<"\t"<<ptt1[2]<<"\n";
      std::cout<<ptt2[0]<<"\t"<<ptt2[1]<<"\t"<<ptt2[2]<<"\n";
      std::cout<<ptt3[0]<<"\t"<<ptt3[1]<<"\t"<<ptt3[2]<<"\n";*/
      /***********************************************/

      split_.push_back(tri);
      ptlist_.erase( ptlist_.begin()+i ); // the main pt of ear is removed, and new polygon is formed
      break;
    }

    if( ptlist_.size()<3 )
      dserror( "ear clipping produced 2 vertices polygon" );

    KERNEL::DeleteInlinePts( ptlist_ ); // delete inline points in the new polygon

    if( ptlist_.size()==3 )
    {
      split_.push_back( ptlist_ );
      break;
    }

    else if( ptlist_.size()==4 && !triOnly )
    {
      Split4nodeFacet( ptlist_, false );
      break;
    }

    std::string str1;
    ptConcavity.clear();
    ptConcavity = KERNEL::CheckConvexity(  ptlist_, str1 ); // concave points for the new polygon

    if( triOnly==false ) // if possible it shifts to splitGeneralFacet so that no of cells are reduced
    {
      if( ptConcavity.size() < 2 )
      {
        SplitConvex_1ptConcave_Facet( ptConcavity );
        return;
      }
      else if( (ptlist_.size()-ptConcavity.size()) > 3 )
      {
        SplitGeneralFacet( ptConcavity );
        return;
      }
    }

    /*std::cout<<"the reflex points are = ";
    for( unsigned i=0;i<reflex.size();i++ )
      std::cout<<reflex[i]<<"\t";
    std::cout<<"\n";
    std::cout<<"the convex points are = ";
    for( unsigned i=0;i<convex.size();i++ )
      std::cout<<convex[i]<<"\t";
    std::cout<<"\n";*/
  }

  /*std::cout<<"the ear clipped points\n";
  for( unsigned i=0;i<split_.size();i++ )
  {
    std::cout<<"triangle\n";
    std::vector<Point*> tr=split_[i];
    for( unsigned j=0;j<tr.size();j++ )
    {
      Point* pp = tr[j];
      double coo[3];
      pp->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/
  //dserror("ear clipped\n");
}


