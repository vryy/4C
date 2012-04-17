#include "cut_triangulateFacet.H"
#include "cut_kernel.H"

/*-----------------------------------------------------------------------------------------------------------*
              Split the facet into appropriate number of tri and quad                           Sudhakar 04/12
              Work well for both convex and concave facets
*------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::SplitFacet()
{
  DeleteInlinePts();

  int numpts = ptlist_.size();
//  std::cout<<"numpts = "<<numpts<<"\n";

  if( numpts<3 )
    dserror("not a valid facet");
  else if( numpts==3 || numpts==4 )
  {
    split_.push_back( ptlist_ );
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
void GEO::CUT::TriangulateFacet::DeleteInlinePts()
{
  bool anyInLine = false;
  unsigned num = ptlist_.size();

  for( unsigned i=0;i<num;i++ )
  {
    Point* pt1 = ptlist_[i];
    Point* pt2 = ptlist_[(i+1)%num];
    unsigned ind = i-1;
    if(i==0)
      ind = num-1;
    Point* pt3 = ptlist_[ind];
    std::vector<Point*>::iterator delPt = ptlist_.begin()+i;
    anyInLine = KERNEL::IsOnLine( pt3, pt1, pt2 );

    if( anyInLine )
    {
      ptlist_.erase(delPt);
      break;
    }
  }
  if( anyInLine )
    DeleteInlinePts();
}

/*---------------------------------------------------------------------------------------------*
   A 5 noded facet is split into a tri and a quad.
   If any 3 points are on the same line, then the middle point is deleted and quad is created
   For concave volumes the point which produces the concavity is found, and used accordingly

           4  _______  3                   3
              |     .\                     |\       1        Sudhakar 04/12
              | B  .  \                    | \    /|
              |   .    \ 2                 |  \ 2/ |
              |  .  A  /                   | B \/  |
              | .     /                    |  .  A |
            0 |._____/ 1                 0 |.______| 1
*---------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::Split5nodeFacet()
{
  std::cout<<"I am in 5 noded facet program\n";
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
    SplitAnyFacet( ptConcavity );
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
  std::cout<<"I am in 6 noded facet program\n";
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
    SplitAnyFacet( ptConcavity );
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
    SplitAnyFacet( ptConcavity );
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
    SplitAnyFacet( ptConcavity );
}

/*---------------------------------------------------------------------------------------------*
   A 9 noded facet is split into 3 quad and a tri
   If any 3 points are on the same line, then the middle point is deleted and 8 point function is called
*----------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::Split9nodeFacet()
{
  std::cout<<"I am in 9 noded\n";
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
    SplitAnyFacet( ptConcavity );
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
    SplitAnyFacet( ptConcavity );


}

/*---------------------------------------------------------------------------------------------------*
     A concave facet which has minimum 2 concavity points are split into appropriate cells
     Also facets those have more than 10 corner points are triangulated
                                                                                       Sudhakar 04/12
*----------------------------------------------------------------------------------------------------*/
void GEO::CUT::TriangulateFacet::SplitAnyFacet( std::vector<int> ptConcavity )
{
  std::cout<<"I am in concave split program\n";
  std::cout<<"the concave points are = ";
  for( unsigned i=0;i<ptConcavity.size();i++ )
    std::cout<<ptConcavity[i]<<"\t";
  std::cout<<"\n";

  int num = ptlist_.size();

  if( ptConcavity.size()<2 && num<11 ) //specific funtions are available to handle such issues
    dserror("Call specific split program, the facet has %d concavity and %d corners",ptConcavity.size(),num);

  split_.clear();
  int concNum = 0;
  bool triDone=false,fresh=true;
  std::vector<Point*> newCell;
  int secondPt,lastPt=0,endPt=0;

  if( ptConcavity[0]==0 && ptConcavity.size()==1 )
    endPt = ptlist_.size()-1;
  else if( ptConcavity.size()==1 )
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
    newCell.push_back(ptlist_[(secondPt+1)%num]);
    lastPt = (secondPt+1)%num;
    fresh = false;
#if 1
    if( lastPt==endPt )
      triDone = true;
    if( lastPt==ptConcavity[concNum+1] )
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
    }
    split_.push_back( newCell );

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
