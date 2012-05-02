#include "cut_triangulateFacet.H"
#include "cut_kernel.H"
//#include "cut_mesh.H"

/*-----------------------------------------------------------------------------------------------------------*
              Split the facet into appropriate number of tri and quad                           Sudhakar 04/12
              Work well for both convex and concave facets
*------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::SplitFacet()
{
  DeleteInlinePts( ptlist_ );


  //Deal with zero point facet -- this should never occur, but happens in axel10
  if( ptlist_.size()==0 )
    return;

/*  std::cout<<"number of points = "<<ptlist_.size()<<"\n";
  std::cout<<"the point list are\n";
  for( unsigned i=0;i<ptlist_.size();i++ )
  {
    Point* ptt = ptlist_[i];
    double coo[3];
    ptt->Coordinates(coo);
    std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
  }*/

  /*********************************************************************/
/*  ptlist_.clear();
  double x[3];
  x[2]=0.0;

  x[0]=0.0;x[1]=0.0;
  GEO::CUT::Point * p1 = mesh_.NewPoint( x, NULL, NULL );
  ptlist_.push_back(p1);

  x[0]=1.0;x[1]=0.0;
  GEO::CUT::Point * p2 = mesh_.NewPoint( x, NULL, NULL );
  ptlist_.push_back(p2);

  x[0]=1.0;x[1]=1.0;
  GEO::CUT::Point * p3 = mesh_.NewPoint( x, NULL, NULL );
  ptlist_.push_back(p3);

  x[0]=2.0;x[1]=1.0;
  GEO::CUT::Point * p4 = mesh_.NewPoint( x, NULL, NULL );
  ptlist_.push_back(p4);

  x[0]=2.0;x[1]=0.0;
  GEO::CUT::Point * p5 = mesh_.NewPoint( x, NULL, NULL );
  ptlist_.push_back(p5);

  x[0]=3.0;x[1]=0.0;
  GEO::CUT::Point * p6 = mesh_.NewPoint( x, NULL, NULL );
  ptlist_.push_back(p6);

  x[0]=3.0;x[1]=2.0;
  GEO::CUT::Point * p7 = mesh_.NewPoint( x, NULL, NULL );
  ptlist_.push_back(p7);

  x[0]=0.0;x[1]=2.0;
  GEO::CUT::Point * p8 = mesh_.NewPoint( x, NULL, NULL );
  ptlist_.push_back(p8);*/
  /*********************************************************************/

  int numpts = ptlist_.size();
//  std::cout<<"numpts = "<<numpts<<"\n";

  if( numpts<3 )
  {
    dserror("not a valid facet");
  }
  else if( numpts==3 )
  {
    split_.clear();
    split_.push_back( ptlist_ );
    return;
  }
  else if( numpts==4 )
  {
    split_.clear();
    Split4nodeFacet( ptlist_ );
    return;
  }
  else if( numpts==5 )
  {
    Split5nodeFacet();
    return;
  }
  else if( numpts==6 )
  {
    Split6nodeFacet();
    return;
  }
  else if( numpts==7 )
  {
    Split7nodeFacet();
    return;
  }
  else if( numpts==8 )
  {
    Split8nodeFacet();
    return;
  }
  else if( numpts==9 )
  {
    Split9nodeFacet();
    return;
  }
  else if( numpts==10 )
  {
    Split10nodeFacet();
    return;
  }
  else
  {
    std::string geoType;
    std::vector<int> ptConcavity = KERNEL::CheckConvexity(  ptlist_, geoType );

    if( ptConcavity.size()==0 ) //this means a convex facet with more than 10 pts
      ptConcavity.push_back(0);

//    std::cout<<"size of concavity = "<<ptConcavity.size()<<"\n";

    SplitAnyFacet( ptConcavity );
  }

  /*std::cout<<"the "<<ptlist_.size()<<" points are\n";
  for( unsigned i=0;i<ptlist_.size();i++ )
  {
    GEO::CUT::Point * p1 = ptlist_[i];
    double x11[3];
    p1->Coordinates(x11);
    std::cout<<x11[0]<<"\t"<<x11[1]<<"\t"<<x11[2]<<"\n";
  }*/
}

/*--------------------------------------------------------------------------------------*
    If more than two points are on a line, all points except the end points
    are deleted                                                             Sudhakar 04/12
*---------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::DeleteInlinePts( std::vector<Point*>& poly )
{
  bool anyInLine = false;
  unsigned num = poly.size();

  for( unsigned i=0;i<num;i++ )
  {
    Point* pt1 = poly[i];
    Point* pt2 = poly[(i+1)%num];
    unsigned ind = i-1;
    if(i==0)
      ind = num-1;
    Point* pt3 = poly[ind];
    std::vector<Point*>::iterator delPt = poly.begin()+i;
    anyInLine = KERNEL::IsOnLine( pt3, pt1, pt2 );

    if( anyInLine )
    {
      poly.erase(delPt);
      break;
    }
  }
  if( anyInLine )
    DeleteInlinePts( poly );
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
void GEO::CUT::TriangulateFacet::Split4nodeFacet( std::vector<Point*> &poly )
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
  else
  {
/*    std::cout<<"the points are\n";
    for( unsigned i=0;i<poly.size();i++ )
    {
      Point* p1 = poly[i];
      double x1[3];
      p1->Coordinates(x1);
      std::cout<<x1[0]<<"\t"<<x1[1]<<"\t"<<x1[2]<<"\n";
    }*/
    dserror( "a 4 noded facet cannot have more than 1 concave point" );
  }
}

/*---------------------------------------------------------------------------------------------*
   A 5 noded facet is split into a tri and a quad.
   If any 3 points are on the same line, then the middle point is deleted and quad is created
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
   A 6 noded facet is split into 2 quad.
   If any 3 points are on the same line, then the middle point is deleted 5 point function is called
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
   A 7 noded facet is split into 2 quad and a tri.
   If any 3 points are on the same line, then the middle point is deleted 6 point function is called
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
  std::cout<<"I am in 7 noded split\n";
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
   A 8 noded facet is split into 3 quad
   If any 3 points are on the same line, then the middle point is deleted 7 point function is called

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
   A 9 noded facet is split into 3 quad and a tri
   If any 3 points are on the same line, then the middle point is deleted and 8 point function is called
*----------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::Split9nodeFacet()
{
//  std::cout<<"I am in 9 noded\n";
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
   A 10 noded facet is split into 4 quad
   If any 3 points are on the same line, then the middle point is deleted; 9 point function is called
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
     Also facets those have more than 10 corner points are triangulated
                                                                                       Sudhakar 04/12
*----------------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::SplitAnyFacet( std::vector<int> ptConcavity )
{
/*  std::cout<<"I am in concave split program\n";
  std::cout<<"the concave points are = ";
  for( unsigned i=0;i<ptConcavity.size();i++ )
    std::cout<<ptConcavity[i]<<"\t";
  std::cout<<"\n";*/

  // if the polygon two adjacent concave points call earclipping
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

  //split_.clear();
  int concNum = 0;
  bool triDone=false,fresh=true;
  std::vector<Point*> newCell;
  int secondPt,lastPt=0,endPt=0;

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
    newCell.push_back(ptlist_[lastPt]);

    fresh = false;
#if 1
    if( lastPt==endPt )
    {
      triDone = true;
    }
    if( concNum!=concsize-1 && lastPt==ptConcavity[concNum+1] )
    {
      fresh = true;
      concNum++;
    }
    else if( !triDone )
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
      DeleteInlinePts( newCell );
      if( newCell.size()<3 )
        dserror( "not a valid cell" );
      else if( newCell.size()==3 )
        split_.push_back( newCell );
      else
      {
        Split4nodeFacet( newCell );
      }
    }
    else
      dserror( "neither tri nor quad" );
#endif

#if 0
    if( concNum<num-1 )
    {
      if( lastPt==ptConcavity[concNum+1] ) // the third point is a concave point
      {                                    // make tri cell and change the concave point
        concNum++;
        fresh = true;
      }
      if( ((secondPt+2)%num)==ptConcavity[concNum+1] ) //the fourth point is a concave point
      {                                                // make a quad cell and change the concave point
        newCell.push_back( ptlist_[(secondPt+2)%num] );
        concNum++;
        lastPt = (secondPt+2)%num;
        fresh = true;
      }
    }

    if( ptConcavity[0]==0 )
    {
      if( lastPt==ptlist_.size()-1 )
        triDone = true;
    }
    else
    {
      if( lastPt==(ptConcavity[0]-1)%num )
      {
        newCell.push_back( ptlist_[ptConcavity[0]] );
        triDone = true;
      }
    }
    split_.push_back( newCell );

    if( lastPt==ptConcavity[0] )
      triDone = true;
#endif

    if( triDone )
      break;
  }
}


/*------------------------------------------------------------------------------------------------------*
            check whether the polygon has two continuous concave points.
            if it is true, then ear clipping should be used to triangulate              sudhakar 04/12
            split facet will fail in such cases
*-------------------------------------------------------------------------------------------------------*/
bool GEO::CUT::TriangulateFacet::HasTwoContinuousConcavePts( std::vector<int> ptConcavity )
{
  int siz = ptConcavity.size();
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
            Triangulation by ear clipping. Works for all cases, but costly.         sudhakar 04/12
*--------------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::EarClipping( std::vector<int> ptConcavity )
{
  std::cout<<"i am ear clipped\n";
  std::vector<int> convex;

/*  std::cout<<"the points are\n";
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

  while(1)
  {
    std::vector<int> reflex( ptConcavity );

    int polPts = ptlist_.size();
    convex.resize( polPts-reflex.size() );

    // Find the convex points
    // They are non-reflex points of the polygon
    int conNum=0;
    if( ptConcavity.size()<2 )
    {
      if( ptConcavity.size()==0 )
        ptConcavity.push_back(0);
      SplitAnyFacet( ptConcavity );
      return;
    }
#if 0
    if( ptConcavity.size()==0 ) //if there are no concave points, or after removing all concave points
    {
      convex.resize( polPts );
      for( int i=0;i<polPts;i++ )
        convex[i] = i;
    }
#endif
    else  // all the non-concave points are copied to the convex
    {
      for( int i=0;i<polPts;i++ )
      {
        if( i==ptConcavity[0] )
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
    }

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
        if( j==ind0 || j==ind2 )
          continue;
        if( PtInsideTriangle( tri,ptlist_[j]) )
        {
          isEar = false;
          break;
        }
      }

      if( !isEar )
        continue;

      split_.push_back(tri);
      ptlist_.erase( ptlist_.begin()+i );
      break;
    }


    if( ptlist_.size()<3 )
      dserror( "ear clipping produced 2 vertices polygon" );
    if( ptlist_.size()==3 )
    {
      split_.push_back( ptlist_ );
      break;
    }

    DeleteInlinePts( ptlist_ ); //delete inline points in the modified polygon

    std::string str1;
    ptConcavity = KERNEL::CheckConvexity(  ptlist_, str1 );

    /*std::cout<<"the reflex points are = ";
    for( unsigned i=0;i<reflex.size();i++ )
      std::cout<<reflex[i]<<"\t";
    std::cout<<"\n";
    std::cout<<"the convex points are = ";
    for( unsigned i=0;i<convex.size();i++ )
      std::cout<<convex[i]<<"\t";
    std::cout<<"\n";*/
  }

/*  for( unsigned i=0;i<split_.size();i++ )
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
  }
  dserror("ear clipped\n");*/
}

/*------------------------------------------------------------------------------------------------------------*
           check whether the point "check" is inside the triangle formed by tri               sudhakar 04/09
                 uses barycentric coordinates as it is faster
*-------------------------------------------------------------------------------------------------------------*/
bool GEO::CUT::TriangulateFacet::PtInsideTriangle( std::vector<Point*> tri, Point* check )
{
  LINALG::Matrix<3,1> t1,t2,t3,pt, v0(0.0), v1(0.0), v2(0.0);
  tri[0]->Coordinates( t1.A() );
  tri[1]->Coordinates( t2.A() );
  tri[2]->Coordinates( t3.A() );
  check->Coordinates( pt.A() );

  v0.Update(1.0, t3, -1.0, t1);
  v0.Update(1.0, t2, -1.0, t1);
  v0.Update(1.0, pt, -1.0, t1);

  double dot00,dot01,dot02,dot11,dot12;
  dot00 = v0.Dot(v0);
  dot01 = v0.Dot(v1);
  dot02 = v0.Dot(v2);
  dot11 = v1.Dot(v1);
  dot12 = v1.Dot(v2);

  double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
  double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
  double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

  if( (u >= 0) && (v >= 0) && (u + v < 1) )
    return true;
  return false;
}

