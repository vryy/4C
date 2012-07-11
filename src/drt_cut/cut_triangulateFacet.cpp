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
  else
  {
    split_.clear();
    switch( numpts )
    {
    case 3:
    {
      split_.push_back( ptlist_ );
      break;
    }
    case 4:
    {
      Split4nodeFacet( ptlist_ );
      break;
    }
    case 5:
    {
      Split5nodeFacet();
      break;
    }
    case 6:
    {
      Split6nodeFacet();
      break;
    }
    case 7:
    {
      Split7nodeFacet();
      break;
    }
    case 8:
    {
      Split8nodeFacet();
      break;
    }
    case 9:
    {
      Split9nodeFacet();
      break;
    }
    case 10:
    {
      Split10nodeFacet();
      break;
    }
    default:
    {
      std::string geoType;
      std::vector<int> ptConcavity = KERNEL::CheckConvexity(  ptlist_, geoType );

      if( ptConcavity.size()==0 ) //this means a convex facet with more than 10 pts
        ptConcavity.push_back(0);

  //    std::cout<<"size of concavity = "<<ptConcavity.size()<<"\n";

      SplitAnyFacet( ptConcavity );
      break;
    }
    }
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
    if( callFromSplitAnyFacet )
    {
      std::cout<<"WARNING!!!! calling earclipping because splitanyfacet failed\n";//blockkk
      ptConcavity = KERNEL::CheckConvexity(  ptlist_, geoType );
      split_.clear();
      EarClipping( ptConcavity, false );
      return;
    }
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

/*---------------------------------------------------------------------------------------------*
   A 5 noded facet is split into a tri and a quad if it does not have two adjacent concave pts
   For concave volumes the point which produces the concavity is found, and used accordingly

           4  _______  3                   3
              |     .\                     |\       1                            Sudhakar 04/12
              | B  .  \                    | \    /|
              |   .    \ 2                 |  \ 2/ |
              |  .  A  /                   | B \/  |
              | .     /                    |  .  A |
            0 |._____/ 1                 0 |.______| 1
*---------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::Split5nodeFacet()
{
//  std::cout<<"I am in 5 noded facet program\n";
  if( ptlist_.size()!=5 )
    dserror("This is not a 5 noded facet");

  split_.clear();

  std::string geoType;
  std::vector<int> ptConcavity = KERNEL::CheckConvexity(  ptlist_, geoType );

  int indStart=0;
  if( geoType=="convex" || geoType=="1ptConcave" )
  {
    if( geoType=="convex" )
      indStart = 0;
    else if( geoType=="1ptConcave" )
      indStart = ptConcavity[0];

    std::vector<Point*> temp1(3),temp2(4);
    temp1[0] = ptlist_[indStart];
    temp1[1] = ptlist_[(indStart+1)%5];
    temp1[2] = ptlist_[(indStart+2)%5];

    temp2[0] = ptlist_[(indStart+2)%5];
    temp2[1] = ptlist_[(indStart+3)%5];
    temp2[2] = ptlist_[(indStart+4)%5];
    temp2[3] = ptlist_[indStart];

    split_.resize(2);
    split_[0] = temp1;
    split_[1] = temp2;
  }
  else
  {
    bool conti = HasTwoContinuousConcavePts( ptConcavity );
    if( conti )
      EarClipping( ptConcavity );
    else
      SplitAnyFacet( ptConcavity );
  }
}

/*---------------------------------------------------------------------------------------------*
   A 6 noded facet is split into 2 quad if it does not have two adjacent concave pts
   For concave volumes the point which produces the concavity is found, and used accordingly

               4  _______  3                   5 ______ 4
                 /      .\                      |      |        Sudhakar 04/12
                /   B  .  \                     |  B   |3____
               /      .    \ 2                  |     .      | 2
            5  \     .  A  /                    |   .        |
                \   .     /                     | .     A    |
                0\ ._____/ 1                  0 |.___________| 1
*---------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::Split6nodeFacet()
{
//  std::cout<<"I am in 6 noded facet program\n";
  if(ptlist_.size()!=6)
    dserror("This is not a 6 noded facet");

  split_.clear();

  std::string geoType;
  std::vector<int> ptConcavity = KERNEL::CheckConvexity(  ptlist_, geoType );

  int indStart=0;
  if( geoType=="convex" || geoType=="1ptConcave" )
  {
    if( geoType=="convex" )
      indStart = 0;
    else if( geoType=="1ptConcave" )
      indStart = ptConcavity[0];

    std::vector<Point*> temp1(4),temp2(4);
    temp1[0] = ptlist_[indStart];
    temp1[1] = ptlist_[(indStart+1)%6];
    temp1[2] = ptlist_[(indStart+2)%6];
    temp1[3] = ptlist_[(indStart+3)%6];

    temp2[0] = ptlist_[(indStart+3)%6];
    temp2[1] = ptlist_[(indStart+4)%6];
    temp2[2] = ptlist_[(indStart+5)%6];
    temp2[3] = ptlist_[indStart];

    split_.resize(2);
    split_[0] = temp1;
    split_[1] = temp2;
  }
  else
  {
    bool conti = HasTwoContinuousConcavePts( ptConcavity );
    if( conti )
      EarClipping( ptConcavity );
    else
      SplitAnyFacet( ptConcavity );
  }
}

/*---------------------------------------------------------------------------------------------*
   A 7 noded facet is split into 2 quad and a tri if it does not have two adjacent concave pts
   For concave volumes the point which produces the concavity is found, and used accordingly


           4  __________  3                  6  ______________ 5
             /.        .\                      | .            |                       Sudhakar 04/12
            / .       .  \                     |   .     A    |
           /  .  B   .    \ 2                  |     .       / 4
        5 |   .     .      \                   |       .    /
          | C .    .        \                  |         . /
          |   .   .         /                  |  B     .  \ 3
        6 |   .  .     A   /                   |      .     \
           \  . .         /                    |    .        \ 2
            \ ..         /                     |  .      C    |
            0\._________/ 1                  0 |._____________| 1
*---------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::Split7nodeFacet()
{
  //std::cout<<"I am in 7 noded split\n";
  if(ptlist_.size()!=7)
    dserror("This is not a 7 noded facet");

  split_.clear();

  std::string geoType;
  std::vector<int> ptConcavity = KERNEL::CheckConvexity(  ptlist_, geoType );

  int indStart=0;
  if( geoType=="convex" || geoType=="1ptConcave" )
  {
    if( geoType=="convex" )
      indStart = 0;
    else if( geoType=="1ptConcave" )
      indStart = ptConcavity[0];

    std::vector<Point*> temp1(4),temp2(3),temp3(4);

    temp1[0] = ptlist_[indStart];
    temp1[1] = ptlist_[(indStart+1)%7];
    temp1[2] = ptlist_[(indStart+2)%7];
    temp1[3] = ptlist_[(indStart+3)%7];

    temp2[0] = ptlist_[(indStart+3)%7];
    temp2[1] = ptlist_[(indStart+4)%7];
    temp2[2] = ptlist_[indStart];

    temp3[0] = ptlist_[(indStart+4)%7];
    temp3[1] = ptlist_[(indStart+5)%7];
    temp3[2] = ptlist_[(indStart+6)%7];
    temp3[3] = ptlist_[indStart];

    split_.resize(3);
    split_[0] = temp1;
    split_[1] = temp2;
    split_[2] = temp3;
  }
  else
  {
    bool conti = HasTwoContinuousConcavePts( ptConcavity );
    if( conti )
      EarClipping( ptConcavity );
    else
      SplitAnyFacet( ptConcavity );
  }
}

/*---------------------------------------------------------------------------------------------*
   A 8 noded facet is split into 3 quad if it does not have two adjacent concave pts

                 5  __________  4                5  __________  4
                   /.         \                    /          |               Sudhakar 04/12
                  / .          \                  /      A    |
                 /  .  B        \ 3              /............|3
              6 |   .         .  |            6 |            .\
                | C .       .    |              |           .  \
                |   .     .      |              |   B     .     \
              7 |   .    .   A   | 2          7 |        .   C   / 2
                 \  .  .        /                \     .        /
                  \ . .        /                  \   .        /
                  0\._________/ 1                 0\._________/ 1
*---------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::Split8nodeFacet()
{
//  std::cout<<"I am in 8 noded facet program\n";
  if(ptlist_.size()!=8)
    dserror("This is not a 8 noded facet");

  split_.clear();

  std::string geoType;
  std::vector<int> ptConcavity = KERNEL::CheckConvexity(  ptlist_, geoType );

  int indStart=0;
  if( geoType=="convex" || geoType=="1ptConcave" )
  {
    if( geoType=="convex" )
      indStart = 0;
    else if( geoType=="1ptConcave" )
      indStart = ptConcavity[0];

    std::vector<Point*> temp1(4),temp2(4),temp3(4);

    temp1[0] = ptlist_[indStart];
    temp1[1] = ptlist_[(indStart+1)%8];
    temp1[2] = ptlist_[(indStart+2)%8];
    temp1[3] = ptlist_[(indStart+3)%8];

    temp2[0] = ptlist_[(indStart+3)%8];
    temp2[1] = ptlist_[(indStart+4)%8];
    temp2[2] = ptlist_[(indStart+5)%8];
    temp2[3] = ptlist_[indStart];

    temp3[0] = ptlist_[(indStart+5)%8];
    temp3[1] = ptlist_[(indStart+6)%8];
    temp3[2] = ptlist_[(indStart+7)%8];
    temp3[3] = ptlist_[indStart];

    split_.resize(3);
    split_[0] = temp1;
    split_[1] = temp2;
    split_[2] = temp3;
  }
  else
  {
    bool conti = HasTwoContinuousConcavePts( ptConcavity );
    if( conti )
      EarClipping( ptConcavity );
    else
      SplitAnyFacet( ptConcavity );
  }
}

/*---------------------------------------------------------------------------------------------*
   A 9 noded facet is split into 3 quad and a tri if it does not have two adjacent concave pts
*----------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::Split9nodeFacet()
{
  if(ptlist_.size()!=9)
    dserror("This is not a 9 noded facet");

  split_.clear();

  std::string geoType;
  std::vector<int> ptConcavity = KERNEL::CheckConvexity(  ptlist_, geoType );

  int indStart=0;
  if( geoType=="convex" || geoType=="1ptConcave" )
  {
    if( geoType=="convex" )
      indStart = 0;
    else if( geoType=="1ptConcave" )
      indStart = ptConcavity[0];

    std::vector<Point*> temp1(4),temp2(4),temp3(4),temp4(3);

    temp1[0] = ptlist_[indStart];
    temp1[1] = ptlist_[(indStart+1)%9];
    temp1[2] = ptlist_[(indStart+2)%9];
    temp1[3] = ptlist_[(indStart+3)%9];

    temp2[0] = ptlist_[(indStart+3)%9];
    temp2[1] = ptlist_[(indStart+4)%9];
    temp2[2] = ptlist_[(indStart+5)%9];
    temp2[3] = ptlist_[indStart];

    temp3[0] = ptlist_[(indStart+5)%9];
    temp3[1] = ptlist_[(indStart+6)%9];
    temp3[2] = ptlist_[(indStart+7)%9];
    temp3[3] = ptlist_[indStart];

    temp4[0] = ptlist_[(indStart+7)%9];
    temp4[1] = ptlist_[(indStart+8)%9];
    temp4[2] = ptlist_[indStart];

    split_.resize(4);
    split_[0] = temp1;
    split_[1] = temp2;
    split_[2] = temp3;
    split_[3] = temp4;
  }
  else
  {
    bool conti = HasTwoContinuousConcavePts( ptConcavity );
    if( conti )
      EarClipping( ptConcavity );
    else
      SplitAnyFacet( ptConcavity );
  }
}

/*---------------------------------------------------------------------------------------------------*
   A 10 noded facet is split into 4 quad if it does not have two adjacent concave pts
*----------------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::Split10nodeFacet()
{
  if( ptlist_.size()!=10 )
    dserror("This is not a 10 noded facet");

  split_.clear();

  std::string geoType;
  std::vector<int> ptConcavity = KERNEL::CheckConvexity(  ptlist_, geoType );

  int indStart=0;
  if( geoType=="convex" || geoType=="1ptConcave" )
  {
    if( geoType=="convex" )
      indStart = 0;
    else if( geoType=="1ptConcave" )
      indStart = ptConcavity[0];

    std::vector<Point*> temp1(4),temp2(4),temp3(4),temp4(4);

    temp1[0] = ptlist_[indStart];
    temp1[1] = ptlist_[(indStart+1)%10];
    temp1[2] = ptlist_[(indStart+2)%10];
    temp1[3] = ptlist_[(indStart+3)%10];

    temp2[0] = ptlist_[(indStart+3)%10];
    temp2[1] = ptlist_[(indStart+4)%10];
    temp2[2] = ptlist_[(indStart+5)%10];
    temp2[3] = ptlist_[indStart];

    temp3[0] = ptlist_[(indStart+5)%10];
    temp3[1] = ptlist_[(indStart+6)%10];
    temp3[2] = ptlist_[(indStart+7)%10];
    temp3[3] = ptlist_[indStart];

    temp4[0] = ptlist_[(indStart+7)%10];
    temp4[1] = ptlist_[(indStart+8)%10];
    temp4[2] = ptlist_[(indStart+9)%10];
    temp4[3] = ptlist_[indStart];

    split_.resize(4);
    split_[0] = temp1;
    split_[1] = temp2;
    split_[2] = temp3;
    split_[3] = temp4;
  }
  else
  {
    bool conti = HasTwoContinuousConcavePts( ptConcavity );
    if( conti )
      EarClipping( ptConcavity );
    else
      SplitAnyFacet( ptConcavity );
  }
}

/*---------------------------------------------------------------------------------------------------*
     A concave facet which has minimum 2 concavity points are split into appropriate cells
     Also facets those have more than 10 corner points are split
                                                                                       Sudhakar 04/12
*----------------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::SplitAnyFacet( std::vector<int> ptConcavity )
{
  // if the polygon has two adjacent concave points this method fails
  // earclipping is called in such situations
  if( ptConcavity.size()>1 )
  {
    bool conti = HasTwoContinuousConcavePts( ptConcavity );
    if( conti )
    {
      EarClipping( ptConcavity );
      return;
    }
  }

  int num = ptlist_.size();
  int concsize = ptConcavity.size();

  int concNum = 0;
  bool triDone=false,fresh=true;
  std::vector<Point*> newCell;
  int secondPt,lastPt=0,endPt=0;

  // endPt is the last point to encounter in splitting
  // once algorithm reached endPT, splitting is finished
  if( ptConcavity[0]==0 && concsize==1 )
    endPt = ptlist_.size()-1;
  else if( concsize==1 )
    endPt = ptConcavity[0]-1;
  else
    endPt = ptConcavity[0];

  while(1)
  {
    newCell.clear();
    int indStart = ptConcavity[concNum];

    // the location of second point is decided by whether we encounter the next
    // concave point or not
    if( fresh )
      secondPt = (indStart+1)%num;
    else
      secondPt = lastPt;

    newCell.push_back(ptlist_[indStart]);
    newCell.push_back(ptlist_[secondPt]);

    lastPt = (secondPt+1)%num;
    newCell.push_back(ptlist_[lastPt]); // now a Tri is formed, and will be checked below whether a Quad can be formed

    fresh = false;

    if( lastPt==endPt )
    {
      triDone = true;
    }
    if( concNum!=concsize-1 && lastPt==ptConcavity[concNum+1] ) // next concave pt reached,
    {                                                           // fresh=true ---> should change indStart for next newCell
      fresh = true;                                             // forming a Quad is impossible in this case
      concNum++;
    }
    else if( !triDone )   // "newCell" which has Tri now, is modified to have a Quad
    {
      lastPt = (secondPt+2)%num;
      newCell.push_back(ptlist_[lastPt]);

      if( concNum<num-1 && lastPt==ptConcavity[concNum+1] )
      {
        fresh = true;
        concNum++;
      }
      if( lastPt==endPt )
        triDone = true;
    }

    if( newCell.size()==3 )
    {
      split_.push_back( newCell );
    }
    else if( newCell.size()==4 ) //check the resulting quad is convex, else split into 2 tri
    {
      KERNEL::DeleteInlinePts( newCell );
      if( newCell.size()<3 )
        dserror( "not a valid cell" );
      else if( newCell.size()==3 )
        split_.push_back( newCell );
      else
      {
        Split4nodeFacet( newCell, true );
      }
    }
    else
      dserror( "neither tri nor quad: Something went wrong in SplitAnyFacet" );

    if( triDone )
      break;
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
  //std::cout<<"I am ear clipped\n";
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

      //std::cout<<"Ear = "<<ind0<<"\t"<<i<<"\t"<<ind2<<"\n";//blockkk
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

    std::string str1;
    ptConcavity.clear();
    ptConcavity = KERNEL::CheckConvexity(  ptlist_, str1 ); // concave points for the new polygon

    if( triOnly==false ) // if possible it shifts to splitAnyFacet so that no of cells are reduced
    {
      if( ptConcavity.size()<2 )
      {
        if( ptConcavity.size()==0 )
          ptConcavity.push_back(0);
        SplitAnyFacet( ptConcavity );
        return;
      }

      else
      {
        bool conti = HasTwoContinuousConcavePts( ptConcavity );
        if( !conti )
        {
          SplitAnyFacet( ptConcavity );
          return;
        }
      }
    }

    /***************************************************************///blockkk
    /*std::cout<<"the modified points are\n";
    for( unsigned mm=0;mm<ptlist_.size();mm++ )
    {
      Point* pt1 = ptlist_[mm];
      double com[3];
      pt1->Coordinates(com);
      std::cout<<com[0]<<"\t"<<com[1]<<"\t"<<com[2]<<"\n";
    }*/
    /***************************************************************/

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
